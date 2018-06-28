#! /usr/bin/env  python

import sys
import numpy as np
from scipy import stats


def G_func(x, G, c):
    return G*x + c


def contact_indicies(data_contcs):
    X0 = data_contcs[0]
    N = len(data_contcs)
    for i in range(N):
        if data_contcs[i] < X0:
            return i

    return N


def get_file(fname):
    data = np.loadtxt(fname)
    data = data.T

    strain = data[0]
    stress = data[2]
    contacts = data[4]
    
    last_idx = contact_indicies(contacts)
    strain = strain[:last_idx]
    stress = stress[:last_idx]

    return strain, stress


def fit(strain_, stress_):
    if len(strain) < 3:
        return -1,-1

    slope, intercept, r_value, p_value, std_err = stats.linregress(strain_, stress_)
    
    return slope, intercept


if __name__ == "__main__":
    fin = sys.argv[1]
    dphi = float(sys.argv[2])
    P0 = float(sys.argv[3])

    strain, stress = get_file(fin)
    G_, c_ = fit(strain, stress)

    foutname = 'GFIT_' + fin
    fout = open(foutname, 'w')
    fout.write(str(dphi) + " " + str(P0) + " "+ str(G_) + " " + str(c_) + "\n")
    fout.close()


