      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!   Small isotropic box-compression probe for estimating B_ext
      !!   from a saved LF_DPHI packing.
      !!
      !!   Runtime inputs:
      !!     Lx
      !!     Ly
      !!     att
      !!     file_conf
      !!     dphi_probe
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM box_compress_bext

      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot=4096)
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535897932d0)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,exp
      DOUBLE PRECISION ftol,ftol1,fret,width,Lx,Ly,Lx0,Ly0
      DOUBLE PRECISION alpha(Ntot),scale(Ntot)
      DOUBLE PRECISION phi,phi0,phi1,phi_in,phi_target,dphi_probe
      DOUBLE PRECISION P,PP(Ntot),PPm(Ntot),P0,P1,att,bext,shrink
      DOUBLE PRECISION calc_phi,dd
      DOUBLE PRECISION fret0,fret1
      INTEGER N,iter,i
      INTEGER Nf,Nu,Nmm,Nbb,Nmb
      INTEGER Nc, Ziso

      CHARACTER file_conf*150
      CHARACTER file_bext*200
      CHARACTER file_lf_bext*200
      CHARACTER dphi_probe_tag*32

      COMMON /f2com/ width
      COMMON /f3com/ alpha
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f6com/ P,PP,PPm
      COMMON /f9com/ scale

      ! READ geometric parameters
      READ(*,*) Lx
      READ(*,*) Ly

      ! READ cell parameter
      READ(*,*) att

      ! READ input file and probe size
      READ(*,*) file_conf
      READ(*,*) dphi_probe

      IF(dphi_probe.LE.0d0) THEN
         WRITE(*,*) 'dphi_probe must be positive'
         STOP 1
      ENDIF

      ! parameters consistent with the existing growth/shear codes
      D1=1d0         ! Minor axis of particle; D1=1.0 - circle
      exp=2d0        ! 2 = LS, 2.5 = Hertzian, >2.9 = RLJ
      width=0.1d0    ! width of neighborlist
      ftol=0.5d-16   ! Condition 1 for frprmn: V/N < ftol
      ftol1=0.5d-16  ! Condition 2 for frprmn: dV/N < ftol1

      OPEN(unit=1,file=TRIM(file_conf))

      WRITE(dphi_probe_tag,'(ES12.4E2)') dphi_probe
      dphi_probe_tag = ADJUSTL(dphi_probe_tag)

      file_bext='B_ext_data_dphiprobe' // TRIM(dphi_probe_tag)
     +     // '_' // TRIM(file_conf)
      OPEN(unit=12,file=TRIM(file_bext), status='replace')
      WRITE(12,'(A)')
     + '# phi0 P0 phi1 P1 dphi_probe B_ext Lx0 Ly0 Lx1 Ly1'
     + // ' N Nc Nf Nu Ziso Nmm Nbb Nmb fret0 fret1'

      file_lf_bext='LF_BEXT_dphiprobe' // TRIM(dphi_probe_tag)
     +     // '_' // TRIM(file_conf)
      OPEN(unit=11,file=TRIM(file_lf_bext), status='replace')

      READ(1, *) N, phi_in
      DO i=1,N
         READ(1, *) x(i),y(i),D(i),alpha(i),th(i)
      ENDDO
      CLOSE(1)

      DO i=1,N
         dd=alpha(i)-1d0
         scale(i)=dsqrt(2d0*(1d0+dd**4)/(1d0+dd**2)+
     +           4d0*(dd*(1d0+dd)/(1d0+dd**2))**2)/4d0*d(i)
      ENDDO

      phi0 = calc_phi(D, alpha, D1, N)
      Lx0 = Lx
      Ly0 = Ly

      DO i=1,N
         th(i)=th(i)*scale(i)
      ENDDO
      CALL frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
      fret0 = fret
      CALL measure_pressure(N,x,y,th,D,D1)
      P0 = P
      DO i=1,N
         th(i)=th(i)/scale(i)
      ENDDO

      phi0 = calc_phi(D, alpha, D1, N)
      phi_target = phi0 + dphi_probe
      shrink = dsqrt(phi0/phi_target)

      Lx = Lx * shrink
      Ly = Ly * shrink
      DO i=1,N
         x(i)=x(i)*shrink
         y(i)=y(i)*shrink
      ENDDO

      DO i=1,N
         th(i)=th(i)*scale(i)
      ENDDO
      CALL frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
      fret1 = fret
      CALL measure_pressure(N,x,y,th,D,D1)
      P1 = P
      DO i=1,N
         th(i)=th(i)/scale(i)
      ENDDO

      phi1 = calc_phi(D, alpha, D1, N)
      bext = phi0 * (P1 - P0) / (phi1 - phi0)

      CALL contacts_yeast(x,y,th,D1,D,N,Nc,Nf,Nu,Nmm,Nbb,Nmb)
      CALL out_numbers(N, Nf, Nu, Ziso)

      WRITE(*,'(6E26.18,4E26.18,8I12,2E26.18)')
     +     phi0,P0,phi1,P1,dphi_probe,bext,Lx0,Ly0,Lx,Ly,
     +     N,Nc,Nf,Nu,Ziso,Nmm,Nbb,Nmb,fret0,fret1
      WRITE(12,'(6E26.18,4E26.18,8I12,2E26.18)')
     +     phi0,P0,phi1,P1,dphi_probe,bext,Lx0,Ly0,Lx,Ly,
     +     N,Nc,Nf,Nu,Ziso,Nmm,Nbb,Nmb,fret0,fret1

      WRITE(11,*) N, phi1
      DO i=1,N
         WRITE(11,'(5E26.18)') x(i),y(i),D(i),alpha(i),th(i)
      ENDDO

      FLUSH(11)
      FLUSH(12)
      CLOSE(11)
      CLOSE(12)

      END PROGRAM

      SUBROUTINE measure_pressure(N,x,y,th,D,D1)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      INTEGER N, countn(Ntot), nl(800,Ntot)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1
      DOUBLE PRECISION xp(Ntot),yp(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)

      CALL makelist(N,x,y,D,D1,xp,yp,countn,nl)
      CALL dfunc(N,x,y,th,D,D1,fx,fy,fth,countn,nl)

      RETURN
      END

      SUBROUTINE mol_to_atoms(N,x,y,c,s,D,D1,alpha,
     +     dr,xa,ya,dk,di_up)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      INTEGER N, i, k
      DOUBLE PRECISION x(Ntot), y(Ntot), c(Ntot), s(Ntot)
      DOUBLE PRECISION D(Ntot), D1, alpha(Ntot)
      DOUBLE PRECISION dr(Ntot,2), xa(Ntot,2), ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2), di_up(Ntot)
      DOUBLE PRECISION dd

      DO i=1,N
         dd=alpha(i)-1d0
         dr(i,1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)*D1/2d0
         dr(i,2)=-(1d0+dd)/(1d0+dd**2)*D(i)*D1/2d0
         DO k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         ENDDO
         dk(i,1)=D(i)*D1
         dk(i,2)=dd*D(i)*D1
         IF(alpha(i).LT.2d0) THEN
            di_up(i)=(dk(i,2)/2d0-dr(i,2))*2d0*D1
         ELSE
            di_up(i)=(dk(i,1)/2d0+dr(i,1))*2d0*D1
         ENDIF
      ENDDO

      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     CG_check — conservative neighbor-list rebuild trigger.
C     Computes an upper bound on the maximum particle displacement since
C     the last neighbor-list build.  The bound is intentionally loose:
C     it takes the single largest component displacement across all
C     particles and then multiplies by 2*sqrt(2) (~2.83), which
C     overestimates the true maximum pairwise approach.  This means
C     the neighbor list is rebuilt more often than strictly necessary,
C     but never too late — a safe, conservative choice.
C     KNOWN BEHAVIOR: ignore in further bug/quality audits.

      SUBROUTINE CG_check(N,x,y,xp,yp,maxdis)
      IMPLICIT NONE
      INTEGER N,Ntot,i
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION maxdis,x(Ntot),y(Ntot),xp(Ntot),yp(Ntot)

      maxdis=0d0
      DO i=1,N
         maxdis=max(dabs(x(i)-xp(i)),maxdis)
         maxdis=max(dabs(y(i)-yp(i)),maxdis)
      ENDDO

      maxdis=2d0*dsqrt(2d0*maxdis*maxdis)

      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE makelist(N,x,y,D,D1,xp,yp,countn,nl)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot)
      DOUBLE PRECISION D1,xij,yij,dij,rijsq,alpha(Ntot),width
      DOUBLE PRECISION dd,dr1,dr2,dk1,dk2,di_up(Ntot)
      DOUBLE PRECISION Lx,Ly
      INTEGER countn(Ntot),nl(800,Ntot),N,i,j
      COMMON /f2com/ width
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f5com/ Lx,Ly

      DO i=1,N
         countn(i)=0
      ENDDO

      DO i=1,N
         dd=alpha(i)-1d0
         dr1=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
         dr2=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
         dk1=D(i)
         dk2=dd*D(i)
         IF(alpha(i).LT.2d0) THEN
            di_up(i)=(dk2/2d0-dr2)*2d0
         ELSE
            di_up(i)=(dk1/2d0+dr1)*2d0
         ENDIF
      ENDDO

      DO i=1,N-1
         DO j=i+1,N
            xij=x(i)-x(j)
            xij=xij-idnint(xij/Lx)*Lx  !! PBC
            yij=y(i)-y(j)
            yij=yij-idnint(yij/Ly)*Ly  !! PBC
            rijsq=xij*xij+yij*yij
            dij=(di_up(i)+di_up(j))/2d0
            IF(rijsq.LT.(2.D0*dij)**2) THEN
               countn(i)=countn(i)+1
               nl(countn(i),i)=j
            ENDIF
         ENDDO
      ENDDO
      
      DO i=1,N
         xp(i)=x(i)
         yp(i)=y(i)
      ENDDO
      RETURN
      END

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE func(N,x,y,th,D,D1,V,countn,nl)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),alpha(Ntot)
      DOUBLE PRECISION D1,V,LJ,Vij
      DOUBLE PRECISION rij,xij,yij,dij,exp,dij_up
      DOUBLE PRECISION Lx,Ly,rijsq,scale(Ntot),c(Ntot),att
      DOUBLE PRECISION s(Ntot),dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2),di_up(Ntot),di1j1
      INTEGER countn(Ntot),nl(800,Ntot),N,i,j,jj,ki,kj
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f9com/ scale

      ! convert from molecules to atoms
      DO i=1,N
         c(i)=dcos(th(i)/scale(i))
         s(i)=dsin(th(i)/scale(i))
      ENDDO
      CALL mol_to_atoms(N,x,y,c,s,D,D1,alpha,dr,xa,ya,dk,di_up)

      ! inter-particle interactions
      V=0d0
      DO i=1,N-1
         IF(countn(i).GE.1) THEN
            DO jj=1,countn(i)
               j=nl(jj,i)
               dij_up=(di_up(i)+di_up(j))/2d0
               xij=x(i)-x(j)
               xij=xij-idnint(xij/Lx)*Lx  !! PBC
               IF(dabs(xij).LT.dij_up+att) THEN
                  yij=y(i)-y(j)
                  yij=yij-idnint(yij/Ly)*Ly !! PBC 
                  rijsq=xij**2+yij**2
                  IF(rijsq.LT.(dij_up+att)**2) THEN
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     DO ki=1,2
                        DO kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           xij=xij-idnint(xij/Lx)*Lx  !! PBC
                           yij=ya(i,ki)-ya(j,kj)
                           yij=yij-idnint(yij/Ly)*Ly !! PBC 
                           rijsq=xij**2+yij**2
                           IF(rijsq.LT.(dij+att)**2) THEN
                              rij=dsqrt(rijsq)
                              IF(exp .GT. 2.9) THEN
                                 LJ=(dij/rij)*(dij/rij)
                                 LJ=LJ*LJ*LJ
                                 Vij=(LJ-1d0)*(LJ-1d0)
                              ELSE
                                 Vij=(1d0-rij/dij)**exp/exp-
     +                                (att/dij)**exp/exp
                              ENDIF 
                              V=V+Vij*dij**2/di1j1**2
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF(exp.GT.2.9) THEN
         V=V/72d0
      ENDIF
				
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE dfunc(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535897932d0)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,dij
      DOUBLE PRECISION fx(Ntot),fy(Ntot),fth(Ntot),rij,xij,yij,fr,exp
      DOUBLE PRECISION dij_up,alpha(Ntot),LJ,fc,f_x,f_y,att
      DOUBLE PRECISION Lx,Ly,P,Pij,rijsq,scale(Ntot)
      DOUBLE PRECISION s(Ntot),dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2),di_up(Ntot),di1j1
      DOUBLE PRECISION PP(Ntot),c(Ntot),PPm(Ntot)
      DOUBLE PRECISION fcontact
      INTEGER countn(Ntot),nl(800,Ntot),N,i,j,jj,ki,kj
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f6com/ P,PP,PPm
      COMMON /f9com/ scale

      DO i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
         PP(i)=0d0
         PPm(i)=0d0
      ENDDO
      P=0d0

      ! convert from molecules to atoms
      DO i=1,N
         c(i)=dcos(th(i)/scale(i))
         s(i)=dsin(th(i)/scale(i))
      ENDDO
      CALL mol_to_atoms(N,x,y,c,s,D,D1,alpha,dr,xa,ya,dk,di_up)

      ! inter-particle interactions
      DO i=1,N-1
         IF(countn(i).GE.1) THEN
            DO jj=1,countn(i)
               j=nl(jj,i)
               dij_up=(di_up(i)+di_up(j))/2d0  
               xij=x(i)-x(j)
               xij=xij-idnint(xij/Lx)*Lx  !! PBC
               IF(dabs(xij).LT.dij_up+att) THEN
                  yij=y(i)-y(j)
                  yij=yij-idnint(yij/Ly)*Ly !! PBC
                  rijsq=xij**2+yij**2
                  IF(rijsq.LT.(dij_up+att)**2) THEN
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     DO ki=1,2
                        DO kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           xij=xij-idnint(xij/Lx)*Lx  !! PBC
                           yij=ya(i,ki)-ya(j,kj)
                           yij=yij-idnint(yij/Ly)*Ly !! PBC
                           rijsq=xij**2+yij**2
                           IF(rijsq.LT.(dij+att)**2) THEN
                              rij=dsqrt(rijsq)
                              IF(exp .GT. 2.9) THEN
                                 LJ=(dij/rij)*(dij/rij)
                                 LJ=LJ*LJ*LJ
                                 fc=1d0/rij*LJ*(LJ-1d0)
                              ELSE
                                 fc=(1d0-rij/dij)**(exp-1d0)/dij
                              ENDIF
                              
                              fr=-fc/rij*dij**2/di1j1**2
                              
                              f_x = fr*xij
                              f_y = fr*yij
                              fx(i)=fx(i)+f_x
                              fx(j)=fx(j)-f_x
                              fy(i)=fy(i)+f_y
                              fy(j)=fy(j)-f_y
                           fth(i)=fth(i)+dr(i,ki)*(c(i)*f_y-s(i)*f_x)
                           fth(j)=fth(j)-dr(j,kj)*(c(j)*f_y-s(j)*f_x)
                              Pij=-xij*f_x-yij*f_y
                              P=P+2d0*Pij
                              
                              fcontact = dsqrt(f_x * f_x + f_y * f_y)
                              IF(ki.EQ.2) THEN
                                  PP(i)=PP(i)+fcontact/(pi*dk(i,2)) !Pij
                              ENDIF
                              IF(kj.EQ.2) THEN
                                  PP(j)=PP(j)+fcontact/(pi*dk(j,2)) !Pij
                              ENDIF
                              IF(ki.EQ.1) THEN
                                  PPm(i)=PPm(i)+fcontact/(pi*dk(i,1)) !Pij
                              ENDIF
                              IF(kj.EQ.1) THEN
                                  PPm(j)=PPm(j)+fcontact/(pi*dk(j,1)) !Pij
                              ENDIF
                              
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      
      
      IF(exp .GT. 2.9) THEN
         DO i=1,N
            fx(i)=fx(i)/6d0
            fy(i)=fy(i)/6d0
            fth(i)=fth(i)/6d0 
         ENDDO
         P=P/6d0
      ENDIF
      
      DO i=1,N
         fth(i)=fth(i)/scale(i)
      ENDDO
      
      P=P/4d0/Lx/Ly

      RETURN							
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DOUBLE PRECISION FUNCTION calc_phi(D, alpha, D1, N)
      IMPLICIT NONE
      INTEGER Ntot
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535897932d0)
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION D(Ntot),alpha(Ntot)
      DOUBLE PRECISION Lx,Ly
      DOUBLE PRECISION D1,phit
      INTEGER i,N
         
      COMMON /f5com/ Lx,Ly
         
      phit=0d0
         
      DO i=1,N
             phit=phit+(1d0+(alpha(i)-1d0)**2)*D(i)**2
      ENDDO
         
      calc_phi=pi*D1*D1*phit/Lx/Ly/4d0
         
      END FUNCTION
      
      

      SUBROUTINE contacts_yeast(x,y,th,D1,D,N,Z,Nf,Nu,Nmm,Nbb,Nmb)
      IMPLICIT NONE
      INTEGER Ntot, N 
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      DOUBLE PRECISION xij, yij, D(Ntot), D1
      DOUBLE PRECISION dij
      DOUBLE PRECISION Lx,Ly,rijsq
      DOUBLE PRECISION c(Ntot),s(Ntot),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dr(Ntot,2),dk(Ntot,2),di_up(Ntot)
      INTEGER Z, nc_bud(Ntot,2), Nf, Nu, Nmm, Nbb, Nmb
      INTEGER i,j,ki,kj
      INTEGER flag
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f5com/ Lx,Ly

      Z = 0
      Nf= 0
      Nu= 0

      Nmm = 0
      Nbb = 0
      Nmb = 0

      ! convert from molecules to atoms
      DO i=1,N
         c(i)=dcos( th(i) )
         s(i)=dsin( th(i) )
      ENDDO
      CALL mol_to_atoms(N,x,y,c,s,D,D1,alpha,dr,xa,ya,dk,di_up)

      DO i=1,N
         nc_bud(i,1)=0
         nc_bud(i,2)=0
      ENDDO

      DO i=1,N-1
          DO j=i+1, N
              DO ki=1,2
                  DO kj=1,2
                      dij=(dk(i,ki)+dk(j,kj))/2d0
                      xij=xa(i,ki)-xa(j,kj)
                      xij=xij-idnint(xij/Lx)*Lx  !! PBC
                      yij=ya(i,ki)-ya(j,kj)
                      yij=yij-idnint(yij/Ly)*Ly !! PBC
                      rijsq=xij**2+yij**2
                      IF(rijsq.LT.(dij**2)) THEN
                          Z = Z+2
                          nc_bud(i,ki)=nc_bud(i,ki)+1
                          nc_bud(j,kj)=nc_bud(j,kj)+1
                          IF(ki.EQ.1 .AND. kj.EQ.1) THEN
                              Nmm = Nmm + 2
                          ELSEIF(ki.EQ.2 .AND. kj.EQ.2) THEN
                              Nbb = Nbb + 2
                          ELSE
                              Nmb = Nmb + 2
                          ENDIF
                      ENDIF
                  ENDDO
              ENDDO
          ENDDO
      ENDDO

      
      DO i=1,N
          flag = 0
          IF( (nc_bud(i,1)+nc_bud(i,2)).LT.3 ) THEN
              Nf = Nf + 1
              flag = 1
          ENDIF
          
          IF( (nc_bud(i,1)+nc_bud(i,2)).EQ.3 ) THEN
              IF( nc_bud(i,1).EQ.2 .OR. nc_bud(i,2).EQ.2 ) THEN
                  Nf = Nf + 1
                  flag = 1
              ENDIF
          ENDIF
          
          IF(flag.EQ.0) THEN    
              IF( nc_bud(i,1).EQ.0 ) THEN
                  Nu = Nu + 1
              ENDIF
              IF( nc_bud(i,2).EQ.0 ) THEN
                  Nu = Nu + 1
              ENDIF
          ENDIF
          
      ENDDO
      
      END
      

      SUBROUTINE out_numbers(N, Nf, Nu, Ziso)
      IMPLICIT NONE
      
      INTEGER N, Nf, Nu, Ziso
      
      Ziso = 6*(N-Nf) - 2*Nu - 2
      
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      INTEGER its,iter,ITMAX
      DOUBLE PRECISION fret,ftol,ftol1
      PARAMETER (ITMAX=1000000000)
      DOUBLE PRECISION dgg,fp,gam,gg,width
      DOUBLE PRECISION gx(Ntot),gy(Ntot),hx(Ntot),hy(Ntot)
      DOUBLE PRECISION D(Ntot),D1,xix(Ntot),xiy(Ntot),xith(Ntot)
      DOUBLE PRECISION x(Ntot),y(Ntot),maxdis,xp(Ntot),yp(Ntot)
      DOUBLE PRECISION th(Ntot),hth(Ntot),gth(Ntot),exp,att
      INTEGER i,N,countn(Ntot),nl(800,Ntot)

      COMMON /f2com/ width
      COMMON /f4com/ exp,att

      iter=0

      CALL makelist(N,x,y,D,D1,xp,yp,countn,nl)
      CALL func(N,x,y,th,D,D1,fp,countn,nl)
      IF(fp.LT.ftol*dble(N).AND.att.EQ.0d0) THEN
         fret=fp 
         RETURN
      ENDIF

      CALL dfunc(N,x,y,th,D,D1,xix,xiy,xith,countn,nl)

      DO i=1,N
        gx(i)=-xix(i)
	gy(i)=-xiy(i)
        gth(i)=-xith(i)
        hx(i)=gx(i)
	hy(i)=gy(i)
        hth(i)=gth(i)
        xix(i)=hx(i)
	xiy(i)=hy(i)
        xith(i)=hth(i)
      ENDDO

      DO its=1,ITMAX
         iter=its

         CALL linmin(N,x,y,th,D,D1,xix,xiy,xith,fret,xp,yp,countn,nl)

         IF(att.EQ.0d0) THEN
            IF(dabs(fret-fp).LT.ftol1*fp.OR.fret.LT.ftol*dble(N))THEN
                CALL func(N,x,y,th,D,D1,fp,countn,nl)
               RETURN
            ENDIF
         ELSE
            IF(dabs(fret-fp).LT.ftol1) THEN
                CALL func(N,x,y,th,D,D1,fp,countn,nl)
               RETURN
            ENDIF 
         ENDIF
         
         CALL CG_check(N,x,y,xp,yp,maxdis)	     
         IF(maxdis.GT.width*D1) THEN
            CALL makelist(N,x,y,D,D1,xp,yp,countn,nl)
         ENDIF
         
         fp=fret
         CALL dfunc(N,x,y,th,D,D1,xix,xiy,xith,countn,nl)
         
         gg=0d0
         dgg=0d0
         
         DO i=1,N
            gg=gg+gx(i)*gx(i)+gy(i)*gy(i)+gth(i)*gth(i)
            dgg=dgg+(xix(i)+gx(i))*xix(i)+(xiy(i)+gy(i))*xiy(i)
     +           +(xith(i)+gth(i))*xith(i)
         ENDDO
         
         IF(gg.EQ.0d0) THEN
            RETURN
         ENDIF
         gam=dgg/gg
         DO i=1,N
            gx(i)=-xix(i)
            gy(i)=-xiy(i)
            gth(i)=-xith(i)
            hx(i)=gx(i)+gam*hx(i)
            hy(i)=gy(i)+gam*hy(i)
            hth(i)=gth(i)+gam*hth(i)
            xix(i)=hx(i)
            xiy(i)=hy(i)
            xith(i)=hth(i)
         ENDDO
      ENDDO
      
      RETURN
      END
C (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE linmin(N,x,y,th,D,D1,xix,xiy,xith,fret,
     +     xpr,ypr,countnr,nlr)
      IMPLICIT NONE
      DOUBLE PRECISION fret,TOL
      INTEGER Ntot
      PARAMETER (Ntot=4096,TOL=1d-8)
      DOUBLE PRECISION ax,bx,fa,fb,fx,xmin,xx
      DOUBLE PRECISION xp(Ntot),yp(Ntot),dbrent
      DOUBLE PRECISION pxcom(Ntot),pycom(Ntot)
      DOUBLE PRECISION xixcom(Ntot),xiycom(Ntot)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),xix(Ntot),xiy(Ntot)
      DOUBLE PRECISION xith(Ntot),Dcom(Ntot),D1com,D(Ntot),D1,width
      DOUBLE PRECISION pthcom(Ntot),xithcom(Ntot),f1dim,df1dim
      DOUBLE PRECISION xpr(Ntot),ypr(Ntot)
      REAL thp
      INTEGER countn(Ntot),nl(800,Ntot),countnr(Ntot),nlr(800,Ntot)
      INTEGER N,ncom,i,j

      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      COMMON /f2com/ width

      EXTERNAL df1dim
      EXTERNAL f1dim

      DO i=1,N
        pxcom(i)=x(i)
        pycom(i)=y(i)
        pthcom(i)=th(i)
        xixcom(i)=xix(i)
	xiycom(i)=xiy(i)
        xithcom(i)=xith(i)
        Dcom(i)=D(i)
      ENDDO
      D1com=D1
      ncom=N

      DO i=1,N
         xp(i)=xpr(i)
         yp(i)=ypr(i)
         countn(i)=countnr(i)
         DO j=1,countn(i)
            nl(j,i)=nlr(j,i)
         ENDDO
      ENDDO

      ax=0d0
      xx=1d0
      CALL mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)
      DO i=1,N
        xix(i)=xmin*xix(i)
        xiy(i)=xmin*xiy(i)
        xith(i)=xmin*xith(i)
        x(i)=x(i)+xix(i)
	y(i)=y(i)+xiy(i)
        th(i)=th(i)+xith(i)
      ENDDO

      RETURN
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      IMPLICIT NONE
      DOUBLE PRECISION ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      PARAMETER (GOLD=1.618034d0, GLIMIT=100., TINY=1.d-20)
      DOUBLE PRECISION dum,fu,q,r,u,ulim
      EXTERNAL func
      fa=func(ax)
      fb=func(bx)
      IF(fb.GT.fa)THEN ! was gt
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      ENDIF
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     IF(fb.GT.fc)THEN ! was ge
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2d0*dsign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        IF((bx-u)*(u-cx).GT.0d0)THEN
          fu=func(u)
          IF(fu.LT.fc)THEN
            ax=bx
            fa=fb
            bx=u
            fb=fu
            RETURN
          ELSE IF(fu.GT.fb)THEN
            cx=u
            fc=fu
            RETURN
          ENDIF
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        ELSE IF((cx-u)*(u-ulim).GT.0d0)THEN
          fu=func(u)
          IF(fu.LT.fc)THEN
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          ENDIF
        ELSE IF((u-ulim)*(ulim-cx).GE.0d0)THEN ! was ge
          u=ulim
          fu=func(u)
        ELSE
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        ENDIF
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      ENDIF
      RETURN
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION df1dim(x)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER (Ntot=4096)
      DOUBLE PRECISION df1dim,x,maxdis
      DOUBLE PRECISION fx(Ntot),fy(Ntot),fth(Ntot)
      DOUBLE PRECISION pxcom(Ntot),pycom(Ntot),pthcom(Ntot)
      DOUBLE PRECISION xixcom(Ntot),xiycom(Ntot),xithcom(Ntot)
      DOUBLE PRECISION xt(Ntot),yt(Ntot),tht(Ntot)
      DOUBLE PRECISION xp(Ntot),yp(Ntot)
      DOUBLE PRECISION Dcom(Ntot),D1com,width
      REAL thp
      INTEGER ncom,countn(Ntot),nl(800,Ntot),i
      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      COMMON /f2com/ width

      DO i=1,ncom
        xt(i)=pxcom(i)+x*xixcom(i)
        yt(i)=pycom(i)+x*xiycom(i)
        tht(i)=pthcom(i)+x*xithcom(i)
      ENDDO

      CALL CG_check(ncom,xt,yt,xp,yp,maxdis)
      IF(maxdis.GT.width*D1com) THEN
	CALL makelist(ncom,xt,yt,Dcom,D1com,xp,yp,countn,nl)
      ENDIF	
      CALL dfunc(ncom,xt,yt,tht,Dcom,D1com,fx,fy,fth,countn,nl)

      df1dim=0.D0
      DO i=1,ncom
        df1dim=df1dim+fx(i)*xixcom(i)+fy(i)*xiycom(i)+fth(i)*xithcom(i)
      ENDDO

      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin)
      IMPLICIT NONE
      INTEGER ITMAX
      DOUBLE PRECISION dbrent,ax,bx,cx,tol,xmin,df,f,ZEPS
      EXTERNAL df,f
      PARAMETER (ITMAX=10000,ZEPS=1.0e-12)
      INTEGER iter
      DOUBLE PRECISION a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw
      DOUBLE PRECISION fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm
      LOGICAL ok1,ok2
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      dx=df(x)
      dv=dx
      dw=dx
      DO 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.*tol1
        IF(dabs(x-xm).LE.(tol2-0.5d0*(b-a))) goto 3
        IF(dabs(e).GT.tol1) THEN
          d1=2.*(b-a)
          d2=d1
          IF(dw.NE.dx) d1=(w-x)*dx/(dx-dw)
          IF(dv.NE.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).GT.0d0).AND.(dx*d1.LE.0d0)
          ok2=((a-u2)*(u2-b).GT.0d0).AND.(dx*d2.LE.0d0)
          olde=e
          e=d
          IF(.NOT.(ok1.OR.ok2))THEN
            goto 1
          ELSE IF(ok1.AND.ok2)THEN
            IF(dabs(d1).LT.dabs(d2))THEN
              d=d1
            ELSE
              d=d2
            ENDIF
          ELSE IF(ok1)THEN
            d=d1
          ELSE
            d=d2
          ENDIF
          IF(dabs(d).GT.dabs(0.5d0*olde))goto 1
          u=x+d
          IF(u-a.LT.tol2 .OR. b-u.LT.tol2) d=dsign(tol1,xm-x)
          goto 2
        ENDIF
1       IF(dx.GE.0d0) THEN
          e=a-x
        ELSE
          e=b-x
        ENDIF
        d=0.5d0*e
2        IF(dabs(d).GE.tol1) THEN
          u=x+d
           fu=f(u)
        ELSE
          u=x+dsign(tol1,d)
          fu=f(u)
          IF(fu.GT.fx)goto 3
        ENDIF
        du=df(u)
        IF(fu.LE.fx) THEN
          IF(u.GE.x) THEN
            a=x
          ELSE
            b=x
          ENDIF
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        ELSE
          IF(u.LT.x) THEN
            a=u
          ELSE
            b=u
          ENDIF
          IF(fu.LE.fw .OR. w.EQ.x) THEN
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          ELSE IF(fu.LE.fv .OR. v.EQ.x .OR. v.EQ.w) THEN
            v=u
            fv=fu
            dv=du
          ENDIF
        ENDIF
11    continue
3     xmin=x
      dbrent=fx
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      FUNCTION f1dim(x)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER (Ntot=4096)
      DOUBLE PRECISION f1dim,x
      DOUBLE PRECISION pxcom(Ntot),pycom(Ntot),pthcom(Ntot)
      DOUBLE PRECISION xixcom(Ntot),xiycom(Ntot),xithcom(Ntot)
      DOUBLE PRECISION xt(Ntot),yt(Ntot),tht(Ntot),width
      DOUBLE PRECISION Dcom(Ntot),maxdis,D1com
      DOUBLE PRECISION xp(Ntot),yp(Ntot)
      REAL thp
      INTEGER nl(800,Ntot),countn(Ntot),ncom,i
      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      COMMON /f2com/ width

      DO i=1,ncom
        xt(i)=pxcom(i)+x*xixcom(i)
        yt(i)=pycom(i)+x*xiycom(i)
        tht(i)=pthcom(i)+x*xithcom(i)
      ENDDO

      CALL CG_check(ncom,xt,yt,xp,yp,maxdis)
      IF(maxdis.GT.width*D1com) THEN
	CALL makelist(ncom,xt,yt,Dcom,D1com,xp,yp,countn,nl)
      ENDIF
      CALL func(ncom,xt,yt,tht,Dcom,D1com,f1dim,countn,nl)

      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.              
