      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!
      !!   Quasi-static growth of cells in a rectangular box
      !!   with Periodic Boundary Conditions and
      !!   via conjugate gradient energy minimization. 
      !! 
      !! 
      !!       Cell type - 1: ellipse (don't use)
      !!                   2: budding
      !!                   3: disk (don't use)
      !!
      !!   Division type - 1: -> ->
      !!                   2: -> <-
      !!                   3: <- ->
      !!                   4: random
      !!
      !!        Feedback - growth rate ~ e^(-P/P0)
      !!                   P0=-1: no feedback
      !!
      !!      Author:    Pawel Gniewek, Carl Schreck
      !!      Email(PG): pawel.gniewek@berkeley.edu
      !!      Email(CS): carl.schreck@berkeley.edu 
      !!      License:   BSD 3
      !!      Reference: "Jamming by growth"; Gniewek, P. and Schreck, C.S. and Hallatschek, O.; 2018
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM jammed_packing

      IMPLICIT NONE
      INTEGER Ntot,Ngen
      PARAMETER(Ntot=4096,Ngen=4096)
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535897932d0)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,exp,ran2
      DOUBLE PRECISION ftol,ftol1,fret,alpha0,width,Lx,Ly,ratei
      DOUBLE PRECISION alpha(Ntot),rate(Ntot),scale(Ntot)
      DOUBLE PRECISION rate0,desync,phi,P,PP(Ntot),D0(Ntot)
      DOUBLE PRECISION dispcm,xa(2),ya(2),PR,PT,P0
      DOUBLE PRECISION cc,ss,dr(2),dd,att,rat
      DOUBLE PRECISION dt, PPm(Ntot), total_growthrate
      INTEGER dtstatus, terminate
      INTEGER N,Nr,seed,iter,i,j,k,kk,m,skip
      INTEGER celltype,divtype
      INTEGER div, Nf, Nu, Nmm,Nbb,Nmb
      CHARACTER file1*80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      CHARACTER file_LF_JAMM*120
      CHARACTER file_LF_DPHI*120
      CHARACTER STATS_file_LF_JAMM*120
      CHARACTER STATS_file_LF_DPHI*120
      CHARACTER file_NC*100
      
      DOUBLE PRECISION x_copy(Ntot),y_copy(Ntot),th_copy(Ntot)
      DOUBLE PRECISION D_copy(Ntot),D0_copy(Ntot)
      DOUBLE PRECISION alpha_copy(Ntot),rate_copy(Ntot)
      DOUBLE PRECISION scale_copy(Ntot)
      DOUBLE PRECISION P_copy,PP_copy(Ntot),PPm_copy(Ntot)
      DOUBLE PRECISION phitemp, calc_phi, wide
      DOUBLE PRECISION phi_j, dphi
      INTEGER N_copy
      INTEGER F(Ntot), Nc, Ziso, F_e(Ntot)
      INTEGER bud_count(Ntot), bud_count_copy(Ntot), num_zero_buds
      INTEGER N_j,Ziso_j,Nc_j,Nf_j,Nu_j,Nmm_j,Nbb_j,Nmb_j
      INTEGER before_jamming, at_jamming, above_jamming
            
      
      COMMON /f2com/ width
      COMMON /f3com/ alpha
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f6com/ P,PP,PT,PR,PPm
      COMMON /f8com/ alpha0
      COMMON /f9com/ scale
      COMMON /f10com/ celltype

      ! READ geometric parameters
      READ(*,*) alpha0
      READ(*,*) Lx
      READ(*,*) Ly

      ! READ cell parameters
      READ(*,*) celltype
      READ(*,*) divtype
      READ(*,*) P0
      READ(*,*) att

      ! READ run parameters
      READ(*,*) rate0
      READ(*,*) desync
      READ(*,*) seed
      READ(*,*) skip

      ! READ pdhi exponent
      READ(*,*) dphi
      
      ! READ output files
      READ(*,*) file1

      IF(celltype.EQ.1 .OR. celltype.EQ.3) THEN
           WRITE(*,*) "ELLIPSE AND DISK CELL TYPE IS NOT SUPPORTED"
           CALL EXIT(0)
      ENDIF
     
      IF(P0.EQ.0d0) THEN
           WRITE(*,*) "P0 = 0 not supported"
           CALL EXIT(0)
      ENDIF
      
      ! parameters
      D1=1d0       ! Minor axis of particle; D1=1.0 - circle
      exp=2d0      ! 2 =  LS, 2.5 = Hertzian, >2.9 = RLJ
      width=0.1d0  ! width of neighborlist 
      ftol=1d-16   ! Condition 1 for frprmn: V/N < ftol 
      ftol1=1d-16  ! Condition 2 for frprmn: dV/N < ftol1
      rat=1.5d0    ! ratio of initial cell 2 to cell 1 volume
      
      phitemp = 0.0d0
      wide = 2.0d0
      total_growthrate = 0.0d0
      
      !dphi = 1d-6

      before_jamming =  0 
      at_jamming = 0
      above_jamming = 0

      ! FILES
 123  OPEN(unit=1,file=TRIM(file1), status='replace') ! CONFIGURATION FILE
      
      ! JAMMED FILES    
      file_LF_JAMM='LF_JAMM_' // TRIM(file1)            ! unit=11
      STATS_file_LF_JAMM='STATS_LF_JAMM_' // TRIM(file1)! unit=21
      
      
      ! DPHI OUTPUT FILES
      file_LF_DPHI='LF_DPHI_' // TRIM(file1)
      OPEN(unit=12,file=TRIM(file_LF_DPHI), status='replace')
      STATS_file_LF_DPHI='STATS_LF_DPHI_' // TRIM(file1)
      OPEN(unit=22,file=TRIM(STATS_file_LF_DPHI), status='replace')
      
      ! SOME STATS
      file_NC='NC_' // TRIM(file1)
      OPEN(unit=13,file=TRIM(file_NC), status='replace')
      
      ! initial size & aspect ratios - total vol = (1+rat)*alpha0
      d(1)=1d0
      d(2)=1d0
      alpha(1)=alpha0
      alpha(2)=dsqrt(alpha0*(1d0+rat)-(alpha0-1d0)**2-2d0)+1d0
     

      ! random initial config
      N=2
      DO i=1,N 
         x(i)=Lx/2 + 0d0
         y(i)=Ly/2 + (dble(i)-1.5d0)*d(i)*D1
         th(i)=(ran2(seed)-0.5d0)*2d0*pi
         D0(i)=D1
         rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0
      ENDDO
      
      dt = 1.0
      dtstatus = 0
      terminate = 0
      
      k=0
      DO WHILE (terminate.NE.1)
          k=k+1
          CALL copy_everything(Ntot,Ngen,
     +     x,y,th,D,alpha,rate,
     +     scale, P, PP, D0,N,
     +     x_copy, y_copy, th_copy, D_copy, alpha_copy, rate_copy,
     +     scale_copy, P_copy, PP_copy, D0_copy, N_copy,
     +     bud_count, bud_count_copy, PPm,PPm_copy)
    
          ! GROW PARTICLES
          DO i=1,N
               ratei=rate(i)*dt
               
               IF(P0.GT.0d0.AND.PP(i).GT.0d0) THEN
                  ratei=ratei*dexp(-PP(i)/P0) !* dt
               ENDif
               
               alpha(i)=1d0+dsqrt((1d0+ratei)*
     +                 (1d0+(alpha(i)-1d0)**2)-1d0)


               IF(alpha(i).GT.2d0*alpha0) THEN
                  dispcm=alpha0/2d0
                  div=1
                  
                  ! divide into 2 - 1st assigned index N+1
                  N=N+1                  
                  D(N)=D(i)
                  x(N)=x(i)+dispcm*dcos(th(i))
                  y(N)=y(i)+dispcm*dsin(th(i))
                  rate(N)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
                  alpha(N)=alpha0
                  
                  ! divide into 2 - 1st assigned index i
                  x(i)=x(i)-dispcm*dcos(th(i))
                  y(i)=y(i)-dispcm*dsin(th(i))
                  rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
                  alpha(i)=alpha0

               ELSE
                   div=0
               ENDif
               
               ! types of division: ->->, <-->, -><-, random
               IF(div.EQ.1) THEN
                  IF(divtype.EQ.1) THEN
                     th(N)=th(i)
                  ELSEIF(divtype.EQ.2) THEN
                     th(N)=th(i)+pi
                  ELSEIF(divtype.EQ.3) THEN
                     th(N)=th(i)
                     th(i)=th(i)+pi
                  ELSEIF(divtype.EQ.4) THEN
                     th(N)=(ran2(seed)-0.5d0)*2d0*pi
                     th(i)=(ran2(seed)-0.5d0)*2d0*pi
                  ENDif

                  ! temp, remove
                  th(i)=th(i) + 1d-4*(ran2(seed)-0.5d0)
                  th(N)=th(N) + 1d-4*(ran2(seed)-0.5d0)
               ENDif
         ENDDO
         
         
         ! convert from angle to length scale = sqrt(I/m) * angle
         DO i=1,N
            dd=alpha(i)-1d0
            scale(i)=dsqrt(2d0*(1d0+dd**4)/(1+dd**2)+
     +               4d0*(dd*(1d0+dd)/(1+dd**2))**2)/4d0*d(i)
     
            th(i)=th(i)*scale(i)
         ENDDO
           
         ! minimize energy
         CALL frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
        
         
         phitemp = calc_phi(D, alpha, D1, N)
         
         WRITE(*,*) k, N, Nc,fret/dble(N), P, dt, phitemp,
     +     before_jamming, at_jamming, above_jamming   

         CALL determine_system_status(N,D,D1,alpha,ftol,wide,dt,
     +     fret, phi_j, dphi,
     +     dtstatus, terminate,
     +     before_jamming, at_jamming, above_jamming)
     
         
         ! REJECT THE MOVE
         IF(dtstatus.EQ.1) THEN
         CALL copy_back_everything(Ntot, Ngen,
     +     x, y, th, D, alpha, rate,
     +     scale, P, PP, D0, N,
     +     x_copy, y_copy, th_copy, D_copy, alpha_copy, rate_copy,
     +     scale_copy, P_copy, PP_copy, D0_copy, N_copy,
     +     bud_count, bud_count_copy,PPm,PPm_copy)       
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! convert back to angles
         DO i=1,N
            th(i)=th(i)/scale(i)
         ENDDO       
         
         IF(mod(k,skip).EQ.0) THEN
            phi = calc_phi(D, alpha, D1, N) 
            WRITE(1,*) 2*N, phi
            DO i=1,N              
               cc=dcos(th(i))
               ss=dsin(th(i))
               dd=alpha(i)-1d0
               dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
               dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
               DO kk=1,2
                  xa(kk)=x(i)+dr(kk)*cc
                  ya(kk)=y(i)+dr(kk)*ss
               ENDDO
               WRITE(1,'(3E26.18,I12)')xa(1),ya(1),d(i),0
               WRITE(1,'(3E26.18,I12)')xa(2),ya(2),d(i)*dd,1
            ENDDO
            FLUSH(1) 
         ENDif
         
         

         IF( at_jamming.EQ.1 ) THEN
            OPEN(unit=11,file=TRIM(file_LF_JAMM), status='replace')
      OPEN(unit=21,file=TRIM(STATS_file_LF_JAMM), status='replace')
            phi = calc_phi(D, alpha, D1, N)
            
            CALL growth_rate(total_growthrate,N,rate,PP,P0,D,alpha)
            
            WRITE(1,*) 2*N, phi
            WRITE(11,*) N, phi
            WRITE(21,*) 2*N, phi, total_growthrate
            DO i=1,N
                WRITE(11,'(5E26.18)') x(i),y(i),D(i),alpha(i),th(i)
                
                cc=dcos(th(i))
                ss=dsin(th(i))
                dd=alpha(i)-1d0
                dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
                dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
                DO kk=1,2
                    xa(kk)=x(i)+dr(kk)*cc
                    ya(kk)=y(i)+dr(kk)*ss
                ENDDO
                WRITE(1,'(3E26.18,I12)')xa(1),ya(1),d(i),0     ! TRA
                WRITE(1,'(3E26.18,I12)')xa(2),ya(2),d(i)*dd,1  ! TRA
                
              WRITE(21,'(3E26.18,I12,4E26.18)')xa(1),ya(1),d(i),0,
     +     PPm(i),rate(i)*dexp(-PP(i)/P0),rate(i),P0
              WRITE(21,'(3E26.18,I12,4E26.18)')xa(2),ya(2),d(i)*dd,1,
     +     PP(i),rate(i)*dexp(-PP(i)/P0),rate(i),P0
            ENDDO
            FLUSH(1)
            FLUSH(11)
            CLOSE(11)
            FLUSH(21)
            CLOSE(21)
      
            CLOSE(13)
            OPEN(unit=13,file=TRIM(file_NC), status='replace')
            CALL contacts_yeast(x,y,th,D1,D,N,Nc,Nf,Nu,Nmm,Nbb,Nmb)
            CALL out_numbers(N, Nf, Nu, Ziso)

           WRITE(13,'(8I8,5E26.18)')N,Ziso,Nc,Nf,Nu,Nmm,Nbb,Nmb,
     +     phi,P,fret,P0,total_growthrate
            FLUSH(13)
            above_jamming = 1
            at_jamming = 0
         ENDif
      ENDDO
      
!!!!!! THE MAIN LOOP ENDS HERE  !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      
      phi = calc_phi(D, alpha, D1, N)
      CALL growth_rate(total_growthrate,N,rate,PP,P0,D,alpha)   
      ! save the last configuration to the file
      WRITE(1,*) 2*N, phi
      WRITE(12,*) N, phi
      WRITE(22,*) 2*N, phi, total_growthrate
      DO i=1,N
          WRITE(12,'(5E26.18)') x(i),y(i),D(i),alpha(i),th(i)
                    
          cc=dcos(th(i))
          ss=dsin(th(i))
          dd=alpha(i)-1d0
          dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
          dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
          DO kk=1,2
              xa(kk)=x(i)+dr(kk)*cc
              ya(kk)=y(i)+dr(kk)*ss
          ENDDO
          
          WRITE(1,'(3E26.18,I12)')xa(1),ya(1),d(i),0     ! TRA
          WRITE(1,'(3E26.18,I12)')xa(2),ya(2),d(i)*dd,1  ! TRA
          

          WRITE(22,'(3E26.18,I12,4E26.18)')xa(1),ya(1),d(i),0,
     +     PPm(i),rate(i)*dexp(-PP(i)/P0),rate(i),P0
          WRITE(22,'(3E26.18,I12,4E26.18)')xa(2),ya(2),d(i)*dd,1,
     +     PP(i),rate(i)*dexp(-PP(i)/P0),rate(i),P0
      ENDDO
      FLUSH(1)
      FLUSH(12)
      FLUSH(22)
      
!!!!!!!!!!!!!

      CALL contacts_yeast(x,y,th,D1,D,N,Nc,Nf,Nu,Nmm,Nbb,Nmb)
      CALL out_numbers(N, Nf, Nu, Ziso)

      WRITE(13,'(8I8,5E26.18)')N,Ziso,Nc,Nf,Nu,Nmm,Nbb,Nmb,
     +     phi,P,fret,P0,total_growthrate
      FLUSH(13)
      
      CLOSE(1)
      CLOSE(12)
      CLOSE(13)
      CLOSE(22)
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      INTEGER Ntot,N
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot),D1
      INTEGER countn(Ntot),nl(800,Ntot)
      
      CALL makelist_dimer(N,x,y,D,D1,xp,yp,countn,nl)

      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE makelist_dimer(N,x,y,D,D1,xp,yp,countn,nl) 
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot)
      DOUBLE PRECISION D1,xij,yij,rij,dij,rijsq,alpha(Ntot),width
      DOUBLE PRECISION dd,dr1,dr2,dk2,di_up(Ntot)
      DOUBLE PRECISION Lx,Ly
      INTEGER countn(Ntot),nl(800,Ntot),N
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
         dk2=dd*D(i)
         di_up(i)=(dk2/2d0-dr2)*2d0
      ENDDO

      DO i=1,N-1
         DO j=i+1,N
            xij=x(i)-x(j)
            xij=xij-idnint(xij/Lx)*Lx  !! PBC
            yij=y(i)-y(j)
            yij=yij-idnint(yij/Ly)*Ly  !! PBC
            rijsq=xij*xij+yij*yij
            dij=(di_up(i)+di_up(j))/2d0
            dij=D1*( D(i) + D(j) )/2d0   
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
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE func(N,x,y,th,D,D1,V,countn,nl)
      IMPLICIT NONE
      INTEGER N,Ntot
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,V
      INTEGER countn(Ntot),nl(800,Ntot)

      CALL func_dimer(N,x,y,th,D,D1,V,countn,nl)

      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE dfunc(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      IMPLICIT NONE
      INTEGER N,Ntot
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1
      DOUBLE PRECISION fx(Ntot),fy(Ntot),fth(Ntot)
      INTEGER countn(Ntot),nl(800,Ntot)

      CALL dfunc_dimer(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE func_dimer(N,x,y,th,D,D1,V,countn,nl)
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535897932d0)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),alpha(Ntot)
      DOUBLE PRECISION D1,V,LJ,Vij
      DOUBLE PRECISION rij,xij,yij,dij,exp,dlnsig,dij_up,sigma
      DOUBLE PRECISION Lx,Ly,rijsq,dijsq_up,scale(Ntot),c(Ntot),att
      DOUBLE PRECISION s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2),di_up(Ntot),di1j1
      INTEGER countn(Ntot),nl(800,Ntot),N
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f9com/ scale

      ! convert from molecules to atoms
      DO i=1,N
         c(i)=dcos(th(i)/scale(i))
         s(i)=dsin(th(i)/scale(i))
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
         ENDif
      ENDDO

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
                              ENDif 
                              V=V+Vij*dij**2/di1j1**2
                           ENDif
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF(exp.GT.2.9) THEN
         V=V/72d0
      ENDif
				
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE dfunc_dimer(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535897932d0)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),sigma,D(Ntot),D1,dij
      DOUBLE PRECISION fx(Ntot),fy(Ntot),fth(Ntot),rij,xij,yij,fr,exp
      DOUBLE PRECISION dij_up,alpha(Ntot),LJ,fc,ft,f_x,f_y,att
      DOUBLE PRECISION fthi,fthj,fth_c,Lx,Ly,P,Pij,rijsq,scale(Ntot)
      DOUBLE PRECISION s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2),di_up(Ntot),di1j1
      DOUBLE PRECISION PP(Ntot),c(Ntot),Vij,PT,PR,PPm(Ntot)
      DOUBLE PRECISION fcontact
      INTEGER countn(Ntot),nl(800,Ntot),N,i
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f6com/ P,PP,PT,PR,PPm
      COMMON /f9com/ scale

      DO i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
         PP(i)=0d0
         PPm(i)=0d0
      ENDDO
      P=0d0
      fcontact=0d0 

      ! convert to from molecules to atoms
      DO i=1,N
         c(i)=dcos(th(i)/scale(i))
         s(i)=dsin(th(i)/scale(i))
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
         ENDif
      ENDDO

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
                              ENDif
                              
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
                              ENDif
                              IF(kj.EQ.2) THEN
                                  PP(j)=PP(j)+fcontact/(pi*dk(j,2)) !Pij
                              ENDif
                              IF(ki.EQ.1) THEN
                                  PPm(i)=PPm(i)+fcontact/(pi*dk(i,1)) !Pij
                              ENDif
                              IF(kj.EQ.1) THEN
                                  PPm(j)=PPm(j)+fcontact/(pi*dk(j,1)) !Pij
                              ENDif
                              
                           ENDif
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
      ENDif
      
      DO i=1,N
         fth(i)=fth(i)/scale(i)
      ENDDO
      
      PT=PT/Lx
      PR=PR/Ly
      P=P/4d0/Lx/Ly

      RETURN							
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
      PARAMETER(Ntot = 4096)
      INTEGER its,iter,ITMAX
      DOUBLE PRECISION fret,ftol,EPS,ftol1
      PARAMETER (EPS=1d-10,ITMAX=1000000000)
      DOUBLE PRECISION dgg,fp,gam,gg,gx(Ntot),gy(Ntot)
      DOUBLE PRECISION hx(Ntot),hy(Ntot)
      DOUBLE PRECISION D(Ntot),D1,xix(Ntot),xiy(Ntot),xith(Ntot)
      DOUBLE PRECISION x(Ntot),y(Ntot),maxdis,xp(Ntot),yp(Ntot)
      DOUBLE PRECISION th(Ntot),hth(Ntot),gth(Ntot),exp,att,V,width
      INTEGER i,N,countn(Ntot),nl(800,Ntot)

      ! not needed
      DOUBLE PRECISION f1,f2,f3,fxe,fye,fthe
      DOUBLE PRECISION xi,yi,thi,max1,max2,max3,del
      DOUBLE PRECISION alpha(Ntot),Lx,Ly,alpha0

      ! not needed
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f5com/ Lx,Ly
      COMMON /f8com/ alpha0

      COMMON /f2com/ width      
      COMMON /f4com/ exp,att

      iter=0

      CALL makelist(N,x,y,D,D1,xp,yp,countn,nl)
      CALL func(N,x,y,th,D,D1,fp,countn,nl)
      IF(fp.LT.ftol*dble(N).AND.att.EQ.0d0) THEN
         fret=fp 
         RETURN
      ENDif

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
         
c         WRITE(*,*) its, fret/dble(N), fret, fp

         IF(att.EQ.0d0) THEN
            IF(dabs(fret-fp).LT.ftol1*fp.OR.fret.LT.ftol*dble(N))THEN
                CALL func(N,x,y,th,D,D1,fp,countn,nl)
               RETURN
            ENDif
         ELSE
            IF(dabs(fret-fp).LT.ftol1) THEN
                CALL func(N,x,y,th,D,D1,fp,countn,nl)
               RETURN
            ENDIF 
         ENDif
         
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
         ENDif
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
      
c     pause 'frprmn maximum iterations exceeded'
      RETURN
      END
C (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE linmin(N,x,y,th,D,D1,xix,xiy,xith,fret,
     +     xpr,ypr,countnr,nlr)
      DOUBLE PRECISION fret,TOL
      PARAMETER (Ntot=4096,TOL=1d-8)
C     USES dbrent,df1dim,mnbrak
      DOUBLE PRECISION ax,bx,fa,fb,fx,xmin,xx,xp(Ntot),yp(Ntot),dbrent
      DOUBLE PRECISION pxcom(Ntot),pycom(Ntot),xixcom(Ntot),xiycom(Ntot)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),xix(Ntot),xiy(Ntot)
      DOUBLE PRECISION xith(Ntot),Dcom(Ntot),D1com,D(Ntot),D1,width
      DOUBLE PRECISION pthcom(Ntot),xithcom(Ntot),f1dim,df1dim
      DOUBLE PRECISION xpr(Ntot),ypr(Ntot)
      INTEGER countn(Ntot),nl(800,Ntot),countnr(Ntot),nlr(800,Ntot)
      INTEGER N,ncom

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
      DOUBLE PRECISION ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
c      EXTERNAL func
      PARAMETER (GOLD=1.618034d0, GLIMIT=100., TINY=1.d-20)
      DOUBLE PRECISION dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      IF(fb.GT.fa)THEN ! was gt
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      ENDif
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
          ENDif
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
          ENDif
        ELSE IF((u-ulim)*(ulim-cx).GE.0d0)THEN ! was ge
          u=ulim
          fu=func(u)
        ELSE
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        ENDif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      ENDif
      RETURN
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION df1dim(x)
      PARAMETER (Ntot=4096)
      DOUBLE PRECISION df1dim,x,maxdis
      DOUBLE PRECISION fx(Ntot),fy(Ntot),fth(Ntot)
      DOUBLE PRECISION pxcom(Ntot),pycom(Ntot),pthcom(Ntot)
      DOUBLE PRECISION xixcom(Ntot),xiycom(Ntot),xithcom(Ntot)
      DOUBLE PRECISION xt(Ntot),yt(Ntot),tht(Ntot)
      DOUBLE PRECISION xp(Ntot),yp(Ntot)
      DOUBLE PRECISION Dcom(Ntot),D1com
      DOUBLE PRECISION alpha(Ntot),width
      INTEGER ncom,countn(Ntot),nl(800,Ntot)
      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      COMMON /f2com/ width
      COMMON /f3com/ alpha ! aspect ratio

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
            ENDif
          ELSE IF(ok1)THEN
            d=d1
          ELSE
            d=d2
          ENDif
          IF(dabs(d).GT.dabs(0.5d0*olde))goto 1
          u=x+d
          IF(u-a.LT.tol2 .OR. b-u.LT.tol2) d=dsign(tol1,xm-x)
          goto 2
        ENDif
1       IF(dx.GE.0d0) THEN
          e=a-x
        ELSE
          e=b-x
        ENDif
        d=0.5d0*e
2        IF(dabs(d).GE.tol1) THEN
          u=x+d
           fu=f(u)
        ELSE
          u=x+dsign(tol1,d)
          fu=f(u)
          IF(fu.GT.fx)goto 3
        ENDif
        du=df(u)
        IF(fu.LE.fx) THEN
          IF(u.GE.x) THEN
            a=x
          ELSE
            b=x
          ENDif
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
          ENDif
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
          ENDif
        ENDif
11    continue
c      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
      RETURN
      END
c  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION f1dim(x)
      PARAMETER (Ntot=4096)
      DOUBLE PRECISION f1dim,x
      DOUBLE PRECISION pxcom(Ntot),pycom(Ntot),pthcom(Ntot)
      DOUBLE PRECISION xixcom(Ntot),xiycom(Ntot),xithcom(Ntot)
      DOUBLE PRECISION xt(Ntot),yt(Ntot),tht(Ntot),width
      DOUBLE PRECISION Dcom(Ntot),maxdis,D1com
      DOUBLE PRECISION xp(Ntot),yp(Ntot),alpha(Ntot)
      INTEGER nl(800,Ntot),countn(Ntot),ncom
      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      
      COMMON /f2com/ width
      COMMON /f3com/ alpha ! aspect ratio

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      DOUBLE PRECISION ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      IF(idum.LE.0) THEN
        idum=max(-idum,1)
        idum2=idum
        DO 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          IF(idum.LT.0) idum=idum+IM1
          IF(j.LE.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      ENDif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      IF(idum.LT.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      IF(idum2.LT.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      IF(iy.LT.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .                 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! BELOW IS A NEW CODE ADDED BY PAWEL GNIEWEK !!!!!!!!
      SUBROUTINE determine_system_status(N,D,D1,alpha,ftol,wide,dt,
     +     fret, phi_j,dphi,
     +     dtstatus, terminate,
     +     before_jamming, at_jamming, above_jamming)
      IMPLICIT NONE
      INTEGER Ntot, N
      PARAMETER(Ntot = 4096)
      
      DOUBLE PRECISION delta, phi, calc_phi
      DOUBLE PRECISION D(Ntot),alpha(Ntot), D1
      DOUBLE PRECISION fret,phitemp,phi_j,ftol,wide,dt,dphi
      INTEGER dtstatus, terminate
      INTEGER before_jamming, at_jamming, above_jamming
      
      delta = 0.5d-8
      
      
      phi = calc_phi(D, alpha, D1, N)
      IF(fret.LT.ftol*N .OR. phi.LT.0.25) THEN
          terminate = 0
          before_jamming = 1
          at_jamming = 0
          above_jamming = 0
          dtstatus = 0
          
          dt = 1d0
          phi_j = 0d0
          RETURN
      ENDif
      
      
      IF(above_jamming.EQ.0 .AND. at_jamming .EQ. 0)THEN
          IF(fret.GT.(ftol * wide * N) )  THEN
             dt = 0.5*dt
             dtstatus = 1
             terminate = 0
             RETURN 
          ELSE IF(fret.GT.(ftol*N) .AND. fret.LT.(ftol*wide*N))THEN
              before_jamming = 0
              at_jamming = 1
              above_jamming = 1
              dt = 1d0
              dtstatus = 0
              terminate = 0
              phi_j = phi
              RETURN
          ! ELSE - not needed. this point should never been reached!
              
        ENDif
      ENDif
      
      
      IF(phi_j .LT. 0.25) THEN
          WRITE(*,*) "ERROR: SHOULD NOT REACH THIS POINT"
      ENDif
      
      IF(at_jamming.NE.0) THEN
          WRITE(*,*) "SOMETHING WRONG: should be at_jamming = 0"
      ENDif

      IF(above_jamming.NE.1) THEN
          WRITE(*,*) "SOMETHING WRONG: should be above_jamming = 1"
      ENDif      
      
      at_jamming = 0
      IF(phi .LT. phi_j+dphi-delta) THEN
          dt = 1
          terminate = 0
          dtstatus = 0
          RETURN
      ELSE IF(phi .GT. phi_j+dphi+delta) THEN
          dt = 0.5 * dt
          dtstatus = 1
          terminate = 0
          RETURN
      ELSE
          dtstatus = 0
          terminate = 1
          RETURN
      ENDif
      
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
      
      
      SUBROUTINE copy_everything(Ntot, Ngen,
     +     x, y, th, D, alpha, rate,      
     +     scale, P, PP, D0,N,
     +     x_copy, y_copy, th_copy, D_copy, alpha_copy, rate_copy,
     +     scale_copy, P_copy, PP_copy, D0_copy, N_copy,
     +     bud_count, bud_count_copy,PPm,PPm_copy)
! ...cc stands for copy
      IMPLICIT NONE
      INTEGER Ntot,Ngen
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot)
      DOUBLE PRECISION alpha(Ntot),rate(Ntot),scale(Ntot)
      DOUBLE PRECISION P,PP(Ntot),D0(Ntot),PPm(Ntot)
      DOUBLE PRECISION x_copy(Ntot),y_copy(Ntot),th_copy(Ntot)
      DOUBLE PRECISION D_copy(Ntot),D0_copy(Ntot)
      DOUBLE PRECISION alpha_copy(Ntot),rate_copy(Ntot)
      DOUBLE PRECISION scale_copy(Ntot)
      DOUBLE PRECISION P_copy,PP_copy(Ntot),PPm_copy(Ntot)
      INTEGER N, N_copy,i,j
      INTEGER bud_count(Ntot), bud_count_copy(Ntot)
      
!      
      P_copy = P
      N_copy = N   
      
! ZERO EVERYTHING !!!
      DO i=1,Ntot
          x_copy(i) = 0d0
          y_copy(i) = 0d0
          th_copy(i) = 0d0
          D_copy(i) = 0d0
          alpha_copy(i) = 0d0
          rate_copy(i) = 0d0
          scale_copy(i) = 0d0
          PP_copy(i) = 0d0
          PPm_copy(i) = 0d0
          D0_copy(i) = 0d0
          bud_count_copy(i) = 0
      ENDDO
     

!!!!!!!!!!!!!!!!!!!!!   MAKE COPY   
      DO i=1,Ntot
          x_copy(i) = x(i)
          y_copy(i) = y(i)
          th_copy(i) = th(i)
          D_copy(i) = D(i)
          alpha_copy(i) = alpha(i)
          rate_copy(i) = rate(i)
          scale_copy(i) = scale(i)
          PP_copy(i) = PP(i)
          PPm_copy(i) = PPm(i)
          D0_copy(i) = D0(i)
          bud_count_copy(i) = bud_count(i)
      ENDDO
      
      
      END
      
      SUBROUTINE copy_back_everything(Ntot, Ngen,  
     +     x, y, th, D, alpha, rate,
     +     scale, P, PP, D0, N,
     +     x_copy, y_copy, th_copy, D_copy, alpha_copy, rate_copy,
     +     scale_copy, P_copy, PP_copy, D0_copy, N_copy,
     +     bud_count, bud_count_copy, PPm, PPm_copy)
! ...cc stands for copy
      IMPLICIT NONE
      INTEGER Ntot,Ngen
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot)
      DOUBLE PRECISION alpha(Ntot),rate(Ntot),scale(Ntot)
      DOUBLE PRECISION P,PP(Ntot),PPm(Ntot)
      DOUBLE PRECISION x_copy(Ntot),y_copy(Ntot),th_copy(Ntot)
      DOUBLE PRECISION alpha_copy(Ntot),rate_copy(Ntot)
      DOUBLE PRECISION scale_copy(Ntot)
      DOUBLE PRECISION D0(Ntot),D0_copy(Ntot),D(Ntot),D_copy(Ntot)
      DOUBLE PRECISION P_copy,PP_copy(Ntot),PPm_copy(Ntot)
      DOUBLE PRECISION dd
      INTEGER N,N_copy,i,j
      INTEGER celltype, bud_count(Ntot), bud_count_copy(Ntot)
      
      COMMON /f10com/ celltype
            
!
      P = P_copy
      N = N_copy
!
      DO i=1,Ntot
          x(i) = x_copy(i)
          y(i) = y_copy(i)
          th(i) = th_copy(i)
          D(i) = D_copy(i)
          alpha(i) = alpha_copy(i)
          rate(i) = rate_copy(i)
          scale(i) = scale_copy(i)
          PP(i) = PP_copy(i)
          PPm(i) = PPm_copy(i)
          D0(i) = D0_copy(i)
          bud_count(i) = bud_count_copy(i)
      ENDDO
      

         DO i=1,N
            IF(celltype.EQ.1) THEN
               scale(i)=dsqrt(1d0+alpha(i)**2)/4d0*d(i)         
            ELSEIF(celltype.EQ.2) THEN
               dd=alpha(i)-1d0
               scale(i)=dsqrt(2d0*(1d0+dd**4)/(1+dd**2)+
     +              4d0*(dd*(1d0+dd)/(1+dd**2))**2)/4d0*d(i)
            ELSEIF(celltype.EQ.3) THEN
               scale(i)=dsqrt(2d0)/4d0*d(i)
            ENDif
            th(i)=th(i)*scale(i)
         ENDDO
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! remove floaters !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE remove_floaters(x,y,th,phi,D1,D,N,F)
      IMPLICIT NONE
      INTEGER Ntot, N 
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION pi
      PARAMETER(pi=3.141592653589793238d0)
      INTEGER R, Rnew, i, j, k, II, JJ, c(Ntot), listP(Ntot), DoF
      DOUBLE PRECISION overlap, aspect_ratio
      INTEGER l, F(Ntot)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),phi,alpha(Ntot)
      DOUBLE PRECISION thIJ, t(3), sign(3), AR2, cont,scale(Ntot)
      DOUBLE PRECISION xij, yij, rij, D(Ntot), D1, z_ave
      DOUBLE PRECISION dij,dtij,dthi,dthj,dtij2,dthi2,dthj2
      DOUBLE PRECISION dtijthi,dtijthj,dthithj, compression
      DOUBLE PRECISION x2(Ntot), y2(Ntot), th2(Ntot), Lx,Ly
      INTEGER count, Z, nn
      
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f5com/ Lx,Ly
      COMMON /f9com/ scale

      ! Remove floaters
      k=0
      phi=0.d0
      DO i=1,N
         ii=i-k
         nn=n-k
         IF(F(ii).EQ.1) THEN
            k=k+1
            DO j=ii,nn
               F(j) = F(j+1)
               x(j) = x(j+1)
               y(j) = y(j+1)
               th(j) = th(j+1)
               D(j) = D(j+1)
               alpha(j) = alpha(j+1)
               scale(j) = scale(j+1)              
            ENDDO
         ENDif
      ENDDO

      N=N-k
      
      DO i=1,N
         phi=phi+D(i)*D(i)
      ENDDO
      phi=pi*D1*D1*phi/Lx/Ly/4d0
      
      END

      
      
!!!!!!!!!!!!!!!! counts contacts !!!!!!!!!!!!!
      SUBROUTINE contacts(x,y,th,D1,D,N,Z,F,Nf,Nu,Nmm,Nbb,Nmb)
      IMPLICIT NONE
      INTEGER Ntot, N        
      PARAMETER(Ntot = 4096)
      INTEGER R, Rnew, i, j, k, II, JJ, c(Ntot), listP(Ntot)
      DOUBLE PRECISION overlap, aspect_ratio
      INTEGER l, F(Ntot), Nf, Nu,Nmm,Nbb,Nmb
      DOUBLE PRECISION x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      DOUBLE PRECISION thIJ, t(3), sign(3),cont
      DOUBLE PRECISION xij, yij, rij, D(Ntot), D1, z_ave
      DOUBLE PRECISION dij,dtij,dthi,dthj,dtij2,dthi2,dthj2
      DOUBLE PRECISION dtijthi,dtijthj,dthithj
      DOUBLE PRECISION x2(Ntot), y2(Ntot), th2(Ntot)
      INTEGER Z
      COMMON /f3com/ alpha ! aspect ratio
      
      Z = 0
      R = 0

      Rnew = 1
      Nf=0
      Nu=0
      Nmm = 0 
      Nbb = 0 
      Nmb = 0
      
      DO i=1, N
         c(i)=0
         listP(i)=i
         F(i) = 0
      ENDDO
      
      DO WHILE(Rnew>0)
         Rnew = 0
         DO i=1, N-R
            II=listP(i)
            c(II)=0
         ENDDO
         
         ! Find Contacts
         DO i=1,N-R
            II=listP(i)
            DO j=1, i-1
               JJ=listP(j)	    
!               IF(overlap(N,x,y,th,D,D1,II,JJ).LT.1.d0) THEN
                  c(II) = c(II)+overlap(N,x,y,th,D,D1,II,JJ) !1
                  c(JJ) = c(JJ)+overlap(N,x,y,th,D,D1,II,JJ) !1
!               ENDif
            ENDDO
         ENDDO

         ! Remove floaters
         i=1
         DO WHILE (i.LE.N-R)
            II=listP(i)
            IF(c(II).LT.3) THEN
               Nf = Nf + 1
               Rnew=Rnew+1
               R=R+1
               DO j=i, N-R
                  listP(j) = listP(j+1)
               ENDDO
               F(II) = 1
            ELSE
               i=i+1
            ENDif
         ENDDO
      ENDDO
      
      
      !Count Contacts
      Z=0
      DO i=1, N-R
         II=listP(i)
         DO j=1, i-1
            JJ=listP(j)	    
!            IF(overlap(N,x,y,th,D,D1,II,JJ).LT.1.d0) THEN
               Z = Z+2*overlap(N,x,y,th,D,D1,II,JJ)
!            ENDif
         ENDDO
      ENDDO

      z_ave = dble(Z)/dble(N-R)

      END
      
      
!     CUSTOM OVERLAP FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION overlap(N,x,y,th,D,D1,i,j)
      IMPLICIT NONE
      INTEGER Ntot, N     
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1
      DOUBLE PRECISION rij,xij,yij,dij,ep,overlap
      DOUBLE PRECISION dtij,dthi,dthj,dtij2,dthi2,dthj2
      DOUBLE PRECISION dtijthi,dtijthj,dthithj, width
      DOUBLE PRECISION alpha(Ntot)
      DOUBLE PRECISION dd,dr1,dr2,dk2,exp,att
      DOUBLE PRECISION Lx,Ly
      !
      DOUBLE PRECISION dijsq_up,scale(Ntot)
      DOUBLE PRECISION s(Ntot),dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2),di_up(Ntot),di1j1      
      INTEGER celltype, ni, i, j
 
      
      COMMON /f2com/ width
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f9com/ scale
      COMMON /f10com/ celltype
      
      IF(D1.NE.1.d0) THEN
          WRITE(*,*) "SOMETHING IS WRONG. D1 SHALL BE EQ. TO 1.0"
      ENDif
      
      overlap = 0.d0

      ! BUDDING CELLS
      IF(celltype.EQ.2) THEN
          WRITE(*,*) "THERE IS SPECIAL CONTACT FUNCTION FOR YEAST"

      ! DISK CELLS      
      ELSEIF(celltype.EQ.3) THEN
                
          xij=x(i)-x(j)
          xij=xij-idnint(xij/Lx)*Lx  !! PBC
          yij=y(i)-y(j)
          yij=yij-idnint(yij/Ly)*Ly  !! PBC
          
          rij=dsqrt(xij**2+yij**2)  
          dij=D1*(d(i)+d(j))/2d0

          IF(rij.LT.dij) THEN
              overlap = 1.d0
          ELSE
              overlap = 0.d0
          ENDIF
      ENDif
          
      END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE contacts_yeast(x,y,th,D1,D,N,Z,Nf,Nu,Nmm,Nbb,Nmb)
      IMPLICIT NONE
      INTEGER Ntot, N 
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      DOUBLE PRECISION xij, yij, D(Ntot), D1
      DOUBLE PRECISION dij
      DOUBLE PRECISION Lx,Ly,rijsq
      DOUBLE PRECISION c(Ntot),s(Ntot),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dd,dr(Ntot,2),dk(Ntot,2)
      INTEGER Z, nc_bud(Ntot,2), Nf, Nu, Nmm, Nbb, Nmb
      INTEGER i,j,ki,kj,k
      INTEGER flag
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f5com/ Lx,Ly

      Z = 0
      Nf= 0
      Nu= 0
      
      Nmm = 0 
      Nbb = 0 
      Nmb = 0
      
      flag = 0
      
      ! convert to from molecules to atoms
      DO i=1,N
         c(i)=dcos( th(i) )
         s(i)=dsin( th(i) )         
         dd=alpha(i)-1d0
         dr(i,1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)*D1/2d0
         dr(i,2)=-(1d0+dd)/(1d0+dd**2)*D(i)*D1/2d0
         DO k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         ENDDO
         dk(i,1)=D(i)*D1
         dk(i,2)=dd*D(i)*D1
      ENDDO
      
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
                          ENDif
                      ENDif
                  ENDDO
              ENDDO
          ENDDO
      ENDDO

      
      DO i=1,N
          flag = 0
          IF( (nc_bud(i,1)+nc_bud(i,2)).LT.3 ) THEN
              Nf = Nf + 1
              flag = 1
          ENDif
          
          IF( (nc_bud(i,1)+nc_bud(i,2)).EQ.3 ) THEN
              IF( nc_bud(i,1).EQ.2 .OR. nc_bud(i,2).EQ.2 ) THEN
                  Nf = Nf + 1
                  flag = 1
              ENDif
          ENDif
          
          IF(flag.EQ.0) THEN    
              IF( nc_bud(i,1).EQ.0 ) THEN
                  nu = nu + 1
              ENDif
              IF( nc_bud(i,2).EQ.0 ) THEN
                  nu = nu + 1
              ENDif
          ENDif
          
      ENDDO
      
      END
      
      SUBROUTINE out_numbers(N, Nf, Nu, Ziso)
      IMPLICIT NONE
      
      INTEGER N, Nf, Nu, Ziso
      
      Ziso = 6*(N-Nf) - 2*Nu - 2
      
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      SUBROUTINE bud_contacts(x,y,th,D1,D,N,bud_count,num_zero_buds)
      IMPLICIT NONE
      INTEGER Ntot, N
      PARAMETER(Ntot = 4096)
      INTEGER i,j,k,ki,kj
      INTEGER nc_bud(Ntot,2),num_zero_buds,bud_count(Ntot)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),alpha(Ntot)
      DOUBLE PRECISION xij,yij,D(Ntot),D1
      DOUBLE PRECISION Lx,Ly,dij,rijsq
      DOUBLE PRECISION c(Ntot),s(Ntot),dd,xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2),dr(Ntot,2)

      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f5com/ Lx,Ly

      num_zero_buds = 0
      
      ! convert to from molecules to atoms
      DO i=1,N
         c(i)=dcos( th(i) )
         s(i)=dsin( th(i) )         
         dd=alpha(i)-1d0
         dr(i,1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)*D1/2d0
         dr(i,2)=-(1d0+dd)/(1d0+dd**2)*D(i)*D1/2d0
         DO k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         ENDDO
         dk(i,1)=D(i)*D1
         dk(i,2)=dd*D(i)*D1
      ENDDO
      
      DO i=1,N
         nc_bud(i,1)=0
         nc_bud(i,2)=0
         bud_count(i) = 0
      ENDDO

      DO i=1,N-1
          DO j=i+1, N
               DO ki=1,2
                  DO kj=1,2
                      dij=(dk(i,ki)+dk(j,kj))/2d0
                      xij=xa(i,ki)-xa(j,kj)
                      xij=xij-idnint(xij/Lx)*Lx !! PBC
                      yij=ya(i,ki)-ya(j,kj)
                      yij=yij-idnint(yij/Ly)*Ly !! PBC
                      rijsq=xij**2+yij**2
                      IF(rijsq.LT.(dij**2)) THEN
                          nc_bud(i,ki)= nc_bud(i,ki)+1
                          nc_bud(j,kj)= nc_bud(j,kj)+1
                      ENDif
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      
      DO i=1,N
          bud_count(i) = nc_bud(i,2)
          IF(bud_count(i).EQ.0)THEN
              num_zero_buds = num_zero_buds + 1
          ENDif
      ENDDO

     
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      SUBROUTINE growth_rate(gr,N,rate,PP,P0,D,alpha)    
      IMPLICIT NONE
      INTEGER Ntot,i,N
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION ratei,P0
      DOUBLE PRECISION rate(Ntot),PP(Ntot),D(Ntot)
      DOUBLE PRECISION alpha(Ntot)
      
      DOUBLE PRECISION Dm, Db, dd, gr
      DOUBLE PRECISION total_area, total_new_mass
      
      gr = 0.0d0
      total_area = 0.0d0
      total_new_mass = 0.0d0
      
      DO i=1,N
          ratei=rate(i)
          IF(P0.GT.0d0.AND.PP(i).GT.0d0) THEN
              ratei=ratei*dexp(-PP(i)/P0)
          ENDif
          
          Dm = D(i)
          Db = alpha(i)-1d0
          total_area = total_area + Dm*Dm + Db*Db
          total_new_mass = total_new_mass + ratei*(Dm*Dm + Db*Db)
      ENDDO
      
      gr = total_new_mass / total_area
      
      END
