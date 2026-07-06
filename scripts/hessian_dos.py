"""Vibrational density of states for 2D budding-cell (mother--bud dumbbell) packings.

This module builds the dynamical matrix of a mechanically arrested budding-cell
packing and returns its normal-mode spectrum, following the model implemented in
``src/jamming_by_growth.f`` / ``src/shear_yeast_linearshear.f``.

------------------------------------------------
Each non-rattler cell ``i`` carries generalized coordinates ``q_i = (x_i, y_i, th_i)``:
the molecule reference point (area-weighted centroid of the two lobes) and the
single rigid-body orientation angle. The bud is rigidly collinear with the cell
axis ``th`` (there is no independent bud-hinge coordinate). With ``dd = alpha - 1``:

    dr1 = (1 + dd) / (1 + dd**2) * dd**2 / 2          (mother offset, +axis)
    dr2 = -(1 + dd) / (1 + dd**2) / 2                 (bud offset, -axis)
    mother_center = (x + dr1*cos th, y + dr1*sin th),   R_mother = 1/2
    bud_center    = (x + dr2*cos th, y + dr2*sin th),   R_bud    = dd/2

Cells interact through purely repulsive harmonic springs between inter-cell lobe
pairs (mother-mother, mother-bud, bud-bud), with unit stiffness ``k = 1``:

    U = sum_{contacts} (k/2) * delta**2,
    delta = (R_a + R_b - r),   r = min-image lobe-center distance,  delta > 0.

The box is square and fully periodic (orthogonal min-image; the zero-strain
Hessian uses no Lees-Edwards shift).

The quasi-zero vibrational modes are the cell-orientation (th) rotations of cells
whose bud is unconstrained: when the bud is small the mother sits near the cell
centroid, so rotating ``th`` mostly swings the contact-free bud at near-zero
energy cost. Completing a bud (adding bud contacts) stiffens that mode into the
bulk. Hence ``N_zero ~ N_u`` and the low-frequency band drains as ``u -> 0``.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh

K_STIFF = 1.0  # harmonic contact stiffness (a.u.); U = (k/2) delta**2
R_MOTHER = 0.5  # mother lobe radius (D0 = 1)


# ----------------------------------------------------------------------------
# Packing I/O
# ----------------------------------------------------------------------------
def _open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


@dataclass
class Packing:
    """A budding-cell packing read from an ``LF_*`` file."""

    x: np.ndarray  # molecule reference x, shape (N,)
    y: np.ndarray  # molecule reference y, shape (N,)
    alpha: np.ndarray  # aspect parameter, shape (N,)
    theta: np.ndarray  # physical cell orientation, shape (N,)
    phi: float
    Lx: float
    Ly: float

    @property
    def N(self):
        return self.x.shape[0]


def read_packing(path, Lx, Ly=None):
    """Read an ``LF_JAMM`` / ``LF_DPHI`` / ``LF_PHI2`` packing file.

    Header line is ``N phi``; each of the ``N`` following rows is
    ``x  y  D0  alpha  theta`` (``D0 == 1``).
    """
    Ly = Lx if Ly is None else Ly
    with _open_text(path) as handle:
        header = handle.readline().split()
        n = int(header[0])
        phi = float(header[1])
        xs = np.empty(n)
        ys = np.empty(n)
        alphas = np.empty(n)
        thetas = np.empty(n)
        for i in range(n):
            cols = handle.readline().split()
            xs[i] = float(cols[0])
            ys[i] = float(cols[1])
            # cols[2] is D0 == 1
            alphas[i] = float(cols[3])
            thetas[i] = float(cols[4])
    return Packing(x=xs, y=ys, alpha=alphas, theta=thetas, phi=phi, Lx=float(Lx), Ly=float(Ly))


# ----------------------------------------------------------------------------
# Geometry: per-cell lobe offsets, radii, lobe centers
# ----------------------------------------------------------------------------
def lobe_offsets(alpha):
    """Return (dr1, dr2, R_bud) for each cell given alpha (bud diameter = alpha-1)."""
    dd = alpha - 1.0
    base = (1.0 + dd) / (1.0 + dd**2)
    dr1 = base * dd**2 / 2.0
    dr2 = -base / 2.0
    r_bud = dd / 2.0
    return dr1, dr2, r_bud


def lobe_centers(x, y, theta, alpha):
    """Return mother and bud lobe centers and radii.

    Returns
    -------
    mother_xy : (N, 2)
    bud_xy    : (N, 2)
    r_mother  : (N,)   (all 0.5)
    r_bud     : (N,)
    """
    dr1, dr2, r_bud = lobe_offsets(alpha)
    c = np.cos(theta)
    s = np.sin(theta)
    mother = np.stack([x + dr1 * c, y + dr1 * s], axis=1)
    bud = np.stack([x + dr2 * c, y + dr2 * s], axis=1)
    r_mother = np.full_like(alpha, R_MOTHER)
    return mother, bud, r_mother, r_bud


# ----------------------------------------------------------------------------
# Contact detection (inter-cell lobe-lobe, periodic min-image)
# ----------------------------------------------------------------------------
@dataclass
class ContactNetwork:
    """A list of inter-cell lobe-lobe contacts in a packing."""

    ci: np.ndarray  # owner cell of lobe a, shape (M,)
    ki: np.ndarray  # lobe index of a (0 mother, 1 bud), shape (M,)
    cj: np.ndarray  # owner cell of lobe b
    kj: np.ndarray  # lobe index of b
    dij: np.ndarray  # contact distance R_a + R_b, shape (M,)

    @property
    def M(self):
        return self.ci.shape[0]


def _min_image(d, L):
    return d - L * np.round(d / L)


def detect_contacts(pack, tol=0.0):
    """Detect inter-cell lobe-lobe overlaps under periodic boundaries.

    A pair (lobe a, lobe b) on different cells is a contact when
    r_ab < R_a + R_b + tol.
    """
    mother, bud, r_mother, r_bud = lobe_centers(pack.x, pack.y, pack.theta, pack.alpha)
    n = pack.N
    # Stack 2 lobes per cell -> 2N lobes; lobe l = 2*cell + k.
    pos = np.empty((2 * n, 2))
    rad = np.empty(2 * n)
    owner = np.empty(2 * n, dtype=int)
    kind = np.empty(2 * n, dtype=int)
    pos[0::2] = mother
    pos[1::2] = bud
    rad[0::2] = r_mother
    rad[1::2] = r_bud
    owner[0::2] = np.arange(n)
    owner[1::2] = np.arange(n)
    kind[0::2] = 0
    kind[1::2] = 1

    nl = 2 * n
    ia, ib = np.triu_indices(nl, k=1)
    dx = _min_image(pos[ia, 0] - pos[ib, 0], pack.Lx)
    dy = _min_image(pos[ia, 1] - pos[ib, 1], pack.Ly)
    r = np.hypot(dx, dy)
    dcontact = rad[ia] + rad[ib]
    mask = (r < dcontact + tol) & (owner[ia] != owner[ib])
    ia, ib = ia[mask], ib[mask]
    return ContactNetwork(
        ci=owner[ia], ki=kind[ia], cj=owner[ib], kj=kind[ib], dij=(rad[ia] + rad[ib]),
    )


def virial_pressure(pack):
    """2D virial pressure over all lobe-lobe contacts.

    P = (1 / (2 Lx Ly)) * sum_contacts delta * r,  with delta = dij - r, k = 1.
    Sums every inter-cell lobe-lobe contact (rattlers included).
    """
    mother, bud, r_mother, r_bud = lobe_centers(pack.x, pack.y, pack.theta, pack.alpha)
    n = pack.N
    pos = np.empty((2 * n, 2))
    rad = np.empty(2 * n)
    owner = np.empty(2 * n, dtype=int)
    pos[0::2], pos[1::2] = mother, bud
    rad[0::2], rad[1::2] = r_mother, r_bud
    owner[0::2] = owner[1::2] = np.arange(n)
    ia, ib = np.triu_indices(2 * n, k=1)
    dx = _min_image(pos[ia, 0] - pos[ib, 0], pack.Lx)
    dy = _min_image(pos[ia, 1] - pos[ib, 1], pack.Ly)
    r = np.hypot(dx, dy)
    dcontact = rad[ia] + rad[ib]
    mask = (r < dcontact) & (owner[ia] != owner[ib])
    delta = dcontact[mask] - r[mask]
    return float(np.sum(delta * r[mask]) / (2.0 * pack.Lx * pack.Ly))


# ----------------------------------------------------------------------------
# Rattler / unconstrained-bud classification
# ----------------------------------------------------------------------------
@dataclass
class Topology:
    is_rattler: np.ndarray  # (N,) bool
    bud_contacts: np.ndarray  # (N,) int, contacts on the bud lobe
    mother_contacts: np.ndarray  # (N,) int, contacts on the mother lobe
    bud_unconstrained: np.ndarray  # (N,) bool: non-rattler cell with 0 bud contacts
    n_rattler: int
    n_nonrattler: int
    n_u: int  # unconstrained lobes among non-rattler cells (Fortran Nu)
    n_bud_unconstrained: int


def classify(pack, contacts):
    """Replicate the ``contacts_yeast`` (Fortran) rattler / Nu logic.

    A cell is a rattler if it has fewer than 3 lobe contacts, or exactly 3 split
    2--1 between its two lobes (one lobe with 2 contacts, the other with 1); a
    3--0 split is not a rattler. Nu counts contact-free lobes among non-rattler
    cells.
    """
    n = pack.N
    mc = np.zeros(n, dtype=int)  # mother-lobe contacts
    bc = np.zeros(n, dtype=int)  # bud-lobe contacts
    for arr_c, arr_k in ((contacts.ci, contacts.ki), (contacts.cj, contacts.kj)):
        for c, k in zip(arr_c, arr_k):
            if k == 0:
                mc[c] += 1
            else:
                bc[c] += 1
    total = mc + bc
    is_rattler = (total < 3) | ((total == 3) & ((mc == 2) | (bc == 2)))
    nonr = ~is_rattler
    # Nu: contact-free lobes belonging to non-rattler cells (mother or bud).
    n_u = int(np.sum(nonr & (mc == 0)) + np.sum(nonr & (bc == 0)))
    bud_unconstrained = nonr & (bc == 0)
    return Topology(
        is_rattler=is_rattler,
        bud_contacts=bc,
        mother_contacts=mc,
        bud_unconstrained=bud_unconstrained,
        n_rattler=int(np.sum(is_rattler)),
        n_nonrattler=int(np.sum(nonr)),
        n_u=n_u,
        n_bud_unconstrained=int(np.sum(bud_unconstrained)),
    )


# ----------------------------------------------------------------------------
# Energy gradient on a fixed contact network (analytic), for finite-diff Hessian
# ----------------------------------------------------------------------------
class FixedNetworkEnergy:
    """Harmonic energy / gradient of a *fixed* inter-cell contact network.

    The contact list is frozen at the reference configuration so the energy is a
    smooth function of ``q`` (no contact opening/closing kinks); finite-differencing
    its gradient yields the full Hessian (stiffness + prestress + geometric terms).
    Coordinates of non-rattler cells only.
    """

    def __init__(self, pack, contacts, topo):
        self.Lx = pack.Lx
        self.Ly = pack.Ly
        self.alpha = pack.alpha
        # Map global cell index -> reduced (non-rattler) DOF index.
        nonr_cells = np.where(~topo.is_rattler)[0]
        self.cells = nonr_cells
        self.cell_to_dof = -np.ones(pack.N, dtype=int)
        self.cell_to_dof[nonr_cells] = np.arange(nonr_cells.shape[0])
        self.ncell = nonr_cells.shape[0]
        self.ndof = 3 * self.ncell

        # Keep only contacts between two non-rattler cells.
        keep = (~topo.is_rattler[contacts.ci]) & (~topo.is_rattler[contacts.cj])
        self.ci = contacts.ci[keep]
        self.ki = contacts.ki[keep]
        self.cj = contacts.cj[keep]
        self.kj = contacts.kj[keep]
        self.dij = contacts.dij[keep]
        self.ncontact = self.ci.shape[0]

        dr1, dr2, r_bud = lobe_offsets(self.alpha)
        self.dr = np.stack([dr1, dr2], axis=1)  # (N,2): offset for lobe 0/1
        self.r_bud = r_bud

        # base coordinates (reference)
        self.q0 = self._pack_q(pack.x, pack.y, pack.theta)

    def _pack_q(self, x, y, theta):
        q = np.empty(self.ndof)
        q[0::3] = x[self.cells]
        q[1::3] = y[self.cells]
        q[2::3] = theta[self.cells]
        return q

    def _lobe_positions(self, q):
        """Lobe centers for the two endpoints of every contact, given reduced q."""
        x = q[0::3]
        y = q[1::3]
        th = q[2::3]
        c = np.cos(th)
        s = np.sin(th)
        # reduced index of each contact's owner cells
        ri = self.cell_to_dof[self.ci]
        rj = self.cell_to_dof[self.cj]
        dr_i = self.dr[self.ci, self.ki]
        dr_j = self.dr[self.cj, self.kj]
        ax = x[ri] + dr_i * c[ri]
        ay = y[ri] + dr_i * s[ri]
        bx = x[rj] + dr_j * c[rj]
        by = y[rj] + dr_j * s[rj]
        return ri, rj, dr_i, dr_j, c, s, ax, ay, bx, by

    def energy(self, q):
        """Harmonic energy U = sum (k/2) delta^2 on the fixed contact network."""
        _, _, _, _, _, _, ax, ay, bx, by = self._lobe_positions(q)
        dx = _min_image(ax - bx, self.Lx)
        dy = _min_image(ay - by, self.Ly)
        r = np.hypot(dx, dy)
        delta = self.dij - r
        return 0.5 * K_STIFF * float(np.sum(delta**2))

    def gradient(self, q):
        ri, rj, dr_i, dr_j, c, s, ax, ay, bx, by = self._lobe_positions(q)
        dx = _min_image(ax - bx, self.Lx)
        dy = _min_image(ay - by, self.Ly)
        r = np.hypot(dx, dy)
        r = np.where(r == 0.0, 1e-15, r)
        delta = self.dij - r
        f = K_STIFF * delta  # spring force magnitude (>0 compressive)
        # dU/d(lobe a position) = -f * n ; n = (dx,dy)/r ; force on a = +f n
        nx = dx / r
        ny = dy / r
        # gradient of U wrt lobe a center = -f * n ; wrt lobe b center = +f * n
        ga_x = -f * nx
        ga_y = -f * ny
        gb_x = f * nx
        gb_y = f * ny

        grad = np.zeros(self.ndof)
        # accumulate into (x,y,theta) of owner cells via chain rule
        # d(lobe)/dx = (1,0); d/dy=(0,1); d/dtheta = dr*(-s, c)
        # lobe a (owner ri):
        np.add.at(grad, 3 * ri + 0, ga_x)
        np.add.at(grad, 3 * ri + 1, ga_y)
        np.add.at(grad, 3 * ri + 2, dr_i * (-s[ri] * ga_x + c[ri] * ga_y))
        # lobe b (owner rj):
        np.add.at(grad, 3 * rj + 0, gb_x)
        np.add.at(grad, 3 * rj + 1, gb_y)
        np.add.at(grad, 3 * rj + 2, dr_j * (-s[rj] * gb_x + c[rj] * gb_y))
        return grad

    def hessian(self, q=None, h=1e-6):
        """Full Hessian via central finite differences of the analytic gradient."""
        if q is None:
            q = self.q0
        n = self.ndof
        H = np.zeros((n, n))
        for a in range(n):
            qp = q.copy()
            qm = q.copy()
            qp[a] += h
            qm[a] -= h
            H[:, a] = (self.gradient(qp) - self.gradient(qm)) / (2.0 * h)
        return 0.5 * (H + H.T)

    def mass_matrix(self, q=None):
        """Consistent mass matrix at constant density (diagonal, per cell).

        The mother lobe carries unit mass; the bud, at the same areal density,
        carries its area fraction, ``m_bud = (r_bud/R_MOTHER)**2``. Each lobe is a
        rigid disk (not a point), so its moment of inertia about ``theta`` is the
        parallel-axis term ``m*dr**2`` plus the disk self-spin ``(1/2) m R**2``:

        M_xx = M_yy = m1 + m2,
        M_thth = (m1*dr1**2 + m2*dr2**2) + (1/2)(m1*R_MOTHER**2 + m2*r_bud**2).

        The reference point is the area-weighted centroid, which (mass ~ area)
        coincides with the centre of mass, so ``m1*dr1 + m2*dr2 = 0`` and the
        translation--rotation coupling vanishes: M is diagonal and orientation
        independent (``q`` is unused, kept for call-signature symmetry with
        ``hessian``).
        """
        dr1 = self.dr[self.cells, 0]
        dr2 = self.dr[self.cells, 1]
        r_bud = self.r_bud[self.cells]
        m1 = np.ones(self.ncell)  # mother lobe, unit mass
        m2 = (r_bud / R_MOTHER) ** 2  # bud, mass ~ area
        mtrans = m1 + m2
        i_theta = (m1 * dr1**2 + m2 * dr2**2) + 0.5 * (m1 * R_MOTHER**2 + m2 * r_bud**2)
        diag = np.empty(self.ndof)
        diag[0::3] = mtrans
        diag[1::3] = mtrans
        diag[2::3] = i_theta
        return np.diag(diag)


# ----------------------------------------------------------------------------
# Spectrum and observables
# ----------------------------------------------------------------------------
@dataclass
class Spectrum:
    omega: np.ndarray  # sorted mode frequencies (global translations removed)
    eigvals: np.ndarray  # corresponding eigenvalues (omega**2)
    n_nonrattler: int
    n_u: int
    n_bud_unconstrained: int
    Z: float  # doubled-contact coordination of the non-rattler network
    u: float
    dZ_mod: float  # Z - Z_iso(u),  Z_iso = 6 - 2u - 2/N_nr
    dZ_naive: float  # Z - 6
    l_c: float  # mean contact length (center-to-center) of the network
    grad_norm: float  # ||dU/dq|| at the reference config (force-balance check)
    n_machine_zero: int  # eigenvalues below 1e-10 before translation removal
    phi: float
    ndof: int


def omega_star(spectrum, lambda_tol=1e-6, n_edge=3):
    """Low-frequency edge of the extended (bulk) band.

    The two global translations are already removed from ``spectrum.eigvals``.
    We further drop the quasi-zero band (the ``N_u`` unconstrained-bud modes,
    ``lambda < lambda_tol``) and estimate the band edge as the root-mean
    frequency of the lowest ``n_edge`` extended eigenvalues,
    ``sqrt(mean(lambda[:n_edge]))``. This is a discrete, reproducible proxy for
    the band edge (not a true inflection-point fit of the cumulative count); it
    vanishes at marginality and grows with the excess coordination.
    """
    lam = np.sort(spectrum.eigvals)
    bulk = lam[lam >= lambda_tol]
    if bulk.size < n_edge + 2:
        return float("nan")
    return float(np.sqrt(np.mean(bulk[:n_edge])))


def analyze_packing(path, Lx, Ly=None, h=1e-6):
    """Full pipeline: read -> contacts -> Hessian -> spectrum + observables."""
    pack = read_packing(path, Lx, Ly)
    contacts = detect_contacts(pack)
    topo = classify(pack, contacts)
    energy = FixedNetworkEnergy(pack, contacts, topo)

    grad0 = energy.gradient(energy.q0)
    grad_norm = float(np.linalg.norm(grad0))

    H = energy.hessian(h=h)
    M = energy.mass_matrix()
    eigvals = eigh(H, M, eigvals_only=True)
    eigvals = np.sort(eigvals)

    # Remove the 2 exact global-translation zero modes (smallest |lambda|).
    # Marginal packings legitimately carry many more quasi-zero modes (the
    # unconstrained-bud rotations, often pushed slightly negative by prestress),
    # so we do NOT require exactly two. The only invariant that must hold is that
    # the two modes we drop are quasi-zero (the uniform x/y translations) and not
    # a real finite-frequency mode -- which fails only for an unrelaxed packing or
    # a disconnected network.
    n_machine_zero = int(np.sum(np.abs(eigvals) < 1e-10))
    assert eigvals[1] < 1e-6, (
        f"2nd-smallest eigenvalue {eigvals[1]:.3e} is not quasi-zero; dropping "
        f"eigvals[2:] would discard a real mode (unrelaxed packing or "
        f"disconnected contact network)"
    )
    nonglobal = eigvals[2:]  # drop the two lowest (uniform x/y translation)
    lam = np.clip(nonglobal, 0.0, None)
    omega = np.sqrt(lam)

    # network coordination (doubled contacts over non-rattlers)
    n_nr = topo.n_nonrattler
    Z = 2.0 * energy.ncontact / n_nr if n_nr else float("nan")
    u = topo.n_u / n_nr if n_nr else float("nan")
    z_iso = 6.0 - 2.0 * u - 2.0 / n_nr
    dZ_mod = Z - z_iso
    dZ_naive = Z - 6.0

    # mean contact length (center-to-center distance of contacting lobes)
    _, _, _, _, _, _, ax, ay, bx, by = energy._lobe_positions(energy.q0)
    dx = _min_image(ax - bx, pack.Lx)
    dy = _min_image(ay - by, pack.Ly)
    l_c = float(np.mean(np.hypot(dx, dy))) if energy.ncontact else float("nan")

    return Spectrum(
        omega=omega,
        eigvals=nonglobal,
        n_nonrattler=n_nr,
        n_u=topo.n_u,
        n_bud_unconstrained=topo.n_bud_unconstrained,
        Z=Z,
        u=u,
        dZ_mod=dZ_mod,
        dZ_naive=dZ_naive,
        l_c=l_c,
        grad_norm=grad_norm,
        n_machine_zero=n_machine_zero,
        phi=pack.phi,
        ndof=energy.ndof,
    )


def count_zero_modes(spectrum, lambda_tol=1e-6):
    """Number of quasi-zero modes (eigenvalue < lambda_tol), translations removed."""
    return int(np.sum(spectrum.eigvals < lambda_tol))


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("usage: python hessian_dos.py <LF_packing[.gz]> [L]")
        raise SystemExit(1)
    path = sys.argv[1]
    L = float(sys.argv[2]) if len(sys.argv) > 2 else 15.0
    spec = analyze_packing(path, L)
    print(f"phi               = {spec.phi:.6f}")
    print(f"N_nonrattler      = {spec.n_nonrattler}")
    print(f"N_u (Fortran)     = {spec.n_u}")
    print(f"N_bud_unconstr.   = {spec.n_bud_unconstrained}")
    print(f"Z (doubled)       = {spec.Z:.4f}")
    print(f"u                 = {spec.u:.4f}")
    print(f"dZ_mod (Z-Ziso)   = {spec.dZ_mod:.4f}")
    print(f"dZ_naive (Z-6)    = {spec.dZ_naive:.4f}")
    print(f"l_c               = {spec.l_c:.4f}")
    print(f"||grad U||        = {spec.grad_norm:.3e}")
    print(f"machine-zero modes= {spec.n_machine_zero} (expect ~2 global translations)")
    print(f"ndof              = {spec.ndof}")
    for tol in (1e-5, 1e-6, 1e-7, 1e-8):
        print(f"N_zero(<{tol:.0e})  = {count_zero_modes(spec, tol)}")
