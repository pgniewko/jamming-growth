      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!
      !!   Quasi-static growth of cells in a rectangular pbc box
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
      !!      Author:   Pawel Gniewek, Carl Schreck
      !!      Email(PG):pawel.gniewek@berkeley.edu
      !!      Email(CS): ...
      !!      License:  BSD 3
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
      DOUBLE PRECISION alpha(Ntot),rate(Ntot),alphar(Ntot),scale(Ntot)
      DOUBLE PRECISION rate0,desync,phi,flow,tdiv,P,PP(Ntot),D0(Ntot)
      DOUBLE PRECISION aclone(Ntot*Ngen),dispcm,xa(2),ya(2),PR,PT,P0
      DOUBLE PRECISION tbirth(Ntot*Ngen),cc,ss,dr(2),dd,corr,att,rat
      DOUBLE PRECISION xbirth(Ntot*Ngen),ybirth(Ntot*Ngen)
      DOUBLE PRECISION xdiv(999,2),ydiv(999,2),thdiv(999,2)
      DOUBLE PRECISION dt, PPm(Ntot), total_growthrate
      INTEGER dtstatus, terminate
      INTEGER N,Nr,seed,iter,i,j,k,kk,c(Ntot,Ngen),m,skip,Nexist
      INTEGER celltype,divtype,nclone(Ntot*Ngen),nclonebox(Ntot*Ngen)
      INTEGER age(Ntot),agehist(Ngen,2),agetot1,agetot2
      INTEGER div,ndiv,idiv(999,2),kdiv, Nf, Nu, Nmm,Nbb,Nmb
      CHARACTER file1*80
      
! NEW DATA 
      CHARACTER file_LF_JAMM*120
      CHARACTER file_LF_DPHI*120
      CHARACTER STATS_file_LF_JAMM*120
      CHARACTER STATS_file_LF_DPHI*120
      CHARACTER file_NC*100
      
      DOUBLE PRECISION xcc(Ntot),ycc(Ntot),thcc(Ntot),Dcc(Ntot)
      DOUBLE PRECISION alphacc(Ntot),ratecc(Ntot),alpharcc(Ntot)
      DOUBLE PRECISION scalecc(Ntot)
      DOUBLE PRECISION Pcc,PPcc(Ntot),D0cc(Ntot),PPmcc(Ntot)
      DOUBLE PRECISION aclonecc(Ntot*Ngen)
      DOUBLE PRECISION tbirthcc(Ntot*Ngen)
      DOUBLE PRECISION xbirthcc(Ntot*Ngen),ybirthcc(Ntot*Ngen)
      DOUBLE PRECISION xdivcc(999,2),ydivcc(999,2),thdivcc(999,2)
      INTEGER Ncc,ccc(Ntot,Ngen),Nexistcc
      INTEGER nclonecc(Ntot*Ngen),ncloneboxcc(Ntot*Ngen)
      INTEGER agecc(Ntot),agehistcc(Ngen,2)
      INTEGER idivcc(999,2), ndivcc
      INTEGER F(Ntot), Nc, Ziso, F_e(Ntot)
      INTEGER BUDCONT(Ntot), BUDCONTCC(Ntot), ZEROBUDS
      DOUBLE PRECISION phitemp, calc_phi, wide
      INTEGER N_j,Ziso_j,Nc_j,Nf_j,Nu_j,Nmm_j,Nbb_j,Nmb_j
      DOUBLE PRECISION phi_j, dphi
      
      INTEGER before_jamming, at_jamming, above_jamming
      
! END OF NEW DATA      
      
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

      IF(celltype.eq.1 .or. celltype.eq.3) THEN
           WRITE(*,*) "ELLIPSE AND DISK CELL TYPE IS NOT SUPPORTED"
           CALL EXIT(0)
      ENDif
     
      IF(P0.eq.0d0) THEN
           WRITE(*,*) "P0 = 0 not supported"
           CALL EXIT(0)
      ENDif
      
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
      Nexist=2
      DO i=1,N 
         c(i,1)=1
         c(i,2)=i
         x(i)=Lx/2 + 0d0
         y(i)=Ly/2 + (dble(i)-1.5d0)*d(i)*D1
         th(i)=(ran2(seed)-0.5d0)*2d0*pi
         D0(i)=D1
         rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
         age(i)=0
      ENDDO

      ! calculate # steps until division
      tdiv=dlog10(2d0)/dlog10(1d0+rate0)
      
      
      dt = 1.0
      dtstatus = 0
      terminate = 0
      ! loop over division
      k=0
      DO while (terminate.ne.1)
      k=k+1
      
      CALL copy_everything(Ntot,Ngen,ndiv,x,y,th,D,alpha,rate,
     +     alphar,
     +     scale, P, PP, D0, aclone,
     +     tbirth, xbirth, ybirth, xdiv, ydiv, thdiv,
     +     N, c, Nexist, nclone, nclonebox, age, agehist, idiv,
     +     xcc, ycc, thcc, Dcc, alphacc, ratecc, alpharcc,
     +     scalecc, Pcc, PPcc, D0cc, aclonecc,
     +     tbirthcc, xbirthcc, ybirthcc, xdivcc, ydivcc, thdivcc,
     +     Ncc,ccc,Nexistcc,nclonecc,ncloneboxcc,agecc,agehistcc,
     +     idivcc, BUDCONT, BUDCONTCC, PPm,PPmcc)
    
         ! GROW PARTICLES
         ndiv=0
         DO i=1,N
               ratei=rate(i)*dt
               
               IF(P0.gt.0d0.and.PP(i).gt.0d0) THEN
                  ratei=ratei*dexp(-PP(i)/P0) !* dt
               ENDif
               
               alpha(i)=1d0+dsqrt((1d0+ratei)*
     +                 (1d0+(alpha(i)-1d0)**2)-1d0)


               IF(alpha(i).gt.2d0*alpha0) THEN
                  dispcm=alpha0/2d0
                  div=1
                  
                  ! divide into 2 - 1st assigned index N+1
                  N=N+1
                  c(N,1)=c(i,1)+1
                  DO j=2,c(N,1)
                     c(N,j)=c(i,j)
                  ENDDO
                  c(N,c(N,1)+1)=Nexist+1                  
                  D(N)=D(i)
                  x(N)=x(i)+dispcm*dcos(th(i))
                  y(N)=y(i)+dispcm*dsin(th(i))
                  rate(N)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
                  alpha(N)=alpha0
                  age(N)=0
                  
                  ! divide into 2 - 1st assigned index i
                  c(i,1)=c(i,1)+1
                  c(i,c(i,1)+1)=Nexist+2
                  x(i)=x(i)-dispcm*dcos(th(i))
                  y(i)=y(i)-dispcm*dsin(th(i))
                  rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
                  alpha(i)=alpha0
                  age(i)=age(i)+1

                  ! initialize # & area of clones
                  nclone(Nexist+1)=0
                  nclone(Nexist+2)=0
                  aclone(Nexist+1)=0d0
                  aclone(Nexist+2)=0d0

                  ! keep track of clones in system
                  nclonebox(Nexist+1)=1
                  nclonebox(Nexist+2)=1
                  DO j=2,c(i,1)
                     nclonebox(c(i,j))=nclonebox(c(i,j))+1
                  ENDDO

                  ! keep track of cell position at birth
                  tbirth(Nexist+1)=dble(k)/tdiv
                  tbirth(Nexist+2)=dble(k)/tdiv
                  xbirth(Nexist+1)=x(N)
                  xbirth(Nexist+2)=x(i)
                  ybirth(Nexist+1)=y(N)
                  ybirth(Nexist+2)=y(i)

                  Nexist=Nexist+2
               else
                   div=0
               ENDif
               
               ! types of division: ->->, <-->, -><-, random
               IF(div.eq.1) THEN
                  IF(divtype.eq.1) THEN
                     th(N)=th(i)
                  elseIF(divtype.eq.2) THEN
                     th(N)=th(i)+pi
                  elseIF(divtype.eq.3) THEN
                     th(N)=th(i)
                     th(i)=th(i)+pi
                  elseIF(divtype.eq.4) THEN
                     th(N)=(ran2(seed)-0.5d0)*2d0*pi
                     th(i)=(ran2(seed)-0.5d0)*2d0*pi
                  ENDif

                  ! temp, remove
                  th(i)=th(i) + 1d-4*(ran2(seed)-0.5d0)
                  th(N)=th(N) + 1d-4*(ran2(seed)-0.5d0)
               ENDif

               IF(div.eq.1) THEN
                  ndiv=ndiv+1
                  idiv(ndiv,1)=i
                  idiv(ndiv,2)=N
                  xdiv(ndiv,1)=x(i)
                  xdiv(ndiv,2)=x(N)
                  ydiv(ndiv,1)=y(i)
                  ydiv(ndiv,2)=y(N)
                  thdiv(ndiv,1)=th(i)
                  thdiv(ndiv,2)=th(N)
               ENDif
         ENDDO
         
         
         ! convert from angle to length scale = sqrt(I/m) * angle
         DO i=1,N
            dd=alpha(i)-1d0
            scale(i)=dsqrt(2d0*(1d0+dd**4)/(1+dd**2)+
     +              4d0*(dd*(1d0+dd)/(1+dd**2))**2)/4d0*d(i)
     
            th(i)=th(i)*scale(i)
         ENDDO
           
         ! minimize energy
         CALL frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
        
         
         phitemp = calc_phi(D, alpha, D1, N)
         
         WRITE(*,*) k, N,Nc,fret/dble(N), P, dt, phitemp,
     +     before_jamming, at_jamming, above_jamming   

         CALL determine_system_status(N,D,D1,alpha,ftol,wide,dt,
     +     fret, phi_j, dphi,
     +     dtstatus, terminate,
     +     before_jamming, at_jamming, above_jamming)
     
         
         ! REJECT THE MOVE
         IF(dtstatus.eq.1) THEN
         CALL copy_back_everything(Ntot, Ngen,ndiv,x,y,th,D,alpha,
     +     rate, 
     +     alphar,
     +     scale, P, PP, D0, aclone,
     +     tbirth, xbirth, ybirth, xdiv, ydiv, thdiv,
     +     N, c, Nexist, nclone, nclonebox, age, agehist, idiv,
     +     xcc, ycc, thcc, Dcc, alphacc, ratecc, alpharcc,
     +     scalecc, Pcc, PPcc, D0cc, aclonecc,
     +     tbirthcc, xbirthcc, ybirthcc, xdivcc, ydivcc, thdivcc,
     +     Ncc, ccc,Nexistcc,nclonecc,ncloneboxcc,agecc,agehistcc,
     +     idivcc, BUDCONT, BUDCONTCC,PPm,PPmcc)       
         ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

!         phitemp = calc_phi(D, alpha, D1, N)
!         WRITE(*,*) "*",k, N,fret/dble(N), P, dt, phitemp,
!     +     before_jamming, at_jamming, above_jamming
         
         ! convert back to angles
         DO i=1,N
            th(i)=th(i)/scale(i)
         ENDDO       
         
         IF(mod(k,skip).eq.0) THEN
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
         
         

         if ( at_jamming.eq.1 ) THEN
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
            close(11)
            FLUSH(21)
            close(21)
      
            close(13)
            OPEN(unit=13,file=TRIM(file_NC), status='replace')
            CALL contacts_yeast(x,y,th,D1,D,N,Nc,F,Nf,Nu,Nmm,Nbb,Nmb)
            CALL out_numbers(N, Nf, Nu, Ziso)

           WRITE(13,'(8I8,5E26.18)')N,Ziso,Nc,Nf,Nu,Nmm,Nbb,Nmb,
     +     phi,P,fret,P0,total_growthrate
            FLUSH(13)
            above_jamming = 1
            at_jamming = 0
         ENDif

         
      ENDDO
!!!!!! THE MAIN LOOP IS OVER
      
      
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

      CALL contacts_yeast(x,y,th,D1,D,N,Nc,F,Nf,Nu,Nmm,Nbb,Nmb)
      CALL out_numbers(N, Nf, Nu, Ziso)

      WRITE(13,'(8I8,5E26.18)')N,Ziso,Nc,Nf,Nu,Nmm,Nbb,Nmb,
     +     phi,P,fret,P0,total_growthrate
      FLUSH(13)
      
      close(1)
      close(12)
      close(13)
      close(22)
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
      DOUBLE PRECISION dd,dr1,dr2,dk2,di_up(Ntot),exp,att
      DOUBLE PRECISION Lx,Ly
      INTEGER countn(Ntot),nl(800,Ntot),N
      COMMON /f2com/ width
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
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
!            IF(rijsq.lt.(2.D0*dij+(att+width)*D1)**2) THEN
            IF(rijsq.lt.(2.D0*dij)**2) THEN
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
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,V,alpha(Ntot)
      DOUBLE PRECISION rij,xij,yij,dij,exp,dlnsig,dij_up,sigma,LJ
      DOUBLE PRECISION Lx,Ly,rijsq,dijsq_up,scale(Ntot),c(Ntot),att
      DOUBLE PRECISION s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2),Vij
      DOUBLE PRECISION dk(Ntot,2),di_up(Ntot),di1j1,xhit,yhit,yhitout
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
         IF(alpha(i).lt.2d0) THEN 
            di_up(i)=(dk(i,2)/2d0-dr(i,2))*2d0*D1
         else
            di_up(i)=(dk(i,1)/2d0+dr(i,1))*2d0*D1         
         ENDif
      ENDDO

      ! inter-particle interactions
      V=0d0
      DO i=1,N-1
         IF(countn(i).ge.1) THEN
            DO jj=1,countn(i)
               j=nl(jj,i)
               dij_up=(di_up(i)+di_up(j))/2d0
               xij=x(i)-x(j)
               xij=xij-idnint(xij/Lx)*Lx  !! PBC
               IF(dabs(xij).lt.dij_up+att) THEN
                  yij=y(i)-y(j)
                  yij=yij-idnint(yij/Ly)*Ly !! PBC 
                  rijsq=xij**2+yij**2
                  IF(rijsq.lt.(dij_up+att)**2) THEN
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     DO ki=1,2
                        DO kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           xij=xij-idnint(xij/Lx)*Lx  !! PBC
                           yij=ya(i,ki)-ya(j,kj)
                           yij=yij-idnint(yij/Ly)*Ly !! PBC 
                           rijsq=xij**2+yij**2
                           IF(rijsq.lt.(dij+att)**2) THEN
                              rij=dsqrt(rijsq)
                              IF(exp .gt. 2.9) THEN
                                 LJ=(dij/rij)*(dij/rij)
                                 LJ=LJ*LJ*LJ
                                 Vij=(LJ-1d0)*(LJ-1d0)
                              else
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

      IF(exp.gt.2.9) THEN
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
      DOUBLE PRECISION dij_up,alpha(Ntot),LJ,fc,ft,f_x,f_y,scale(Ntot)
      DOUBLE PRECISION fthi,fthj,fth_c,Lx,Ly,P,Pij,rijsq,dijsq_up
      DOUBLE PRECISION s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2),att
      DOUBLE PRECISION dk(Ntot,2),di_up(Ntot),di1j1,xhit,yhit,yhitout
      DOUBLE PRECISION PP(Ntot),c(Ntot),Vij,PT,PR,PPm(Ntot)
      DOUBLE PRECISION fcontact
      INTEGER countn(Ntot),nl(800,Ntot),N !,growth_flag
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f6com/ P,PP,PT,PR,PPm
      COMMON /f9com/ scale
!      COMMON /f11com/ growth_flag

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
         IF(alpha(i).lt.2d0) THEN 
            di_up(i)=(dk(i,2)/2d0-dr(i,2))*2d0*D1
         else
            di_up(i)=(dk(i,1)/2d0+dr(i,1))*2d0*D1         
         ENDif
      ENDDO

      ! inter-particle interactions
      DO i=1,N-1
         IF(countn(i).ge.1) THEN
            DO jj=1,countn(i)
               j=nl(jj,i)
               dij_up=(di_up(i)+di_up(j))/2d0  
               xij=x(i)-x(j)
               xij=xij-idnint(xij/Lx)*Lx  !! PBC
               IF(dabs(xij).lt.dij_up+att) THEN
                  yij=y(i)-y(j)
                  yij=yij-idnint(yij/Ly)*Ly !! PBC
                  rijsq=xij**2+yij**2
                  IF(rijsq.lt.(dij_up+att)**2) THEN
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     DO ki=1,2
                        DO kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           xij=xij-idnint(xij/Lx)*Lx  !! PBC
                           yij=ya(i,ki)-ya(j,kj)
                           yij=yij-idnint(yij/Ly)*Ly !! PBC
                           rijsq=xij**2+yij**2
                           IF(rijsq.lt.(dij+att)**2) THEN
                              rij=dsqrt(rijsq)
                              IF(exp .gt. 2.9) THEN
                                 LJ=(dij/rij)*(dij/rij)
                                 LJ=LJ*LJ*LJ
                                 fc=1d0/rij*LJ*(LJ-1d0)
                              else
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
                              if (ki.eq.2) THEN
                                  PP(i)=PP(i)+fcontact/(pi*dk(i,2)) !Pij
                              ENDif
                              if (kj.eq.2) THEN
                                  PP(j)=PP(j)+fcontact/(pi*dk(j,2)) !Pij
                              ENDif
                              if (ki.eq.1) THEN
                                  PPm(i)=PPm(i)+fcontact/(pi*dk(i,1)) !Pij
                              ENDif
                              if (kj.eq.1) THEN
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

      
      
      IF(exp .gt. 2.9) THEN
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
!      DO i=1,N
!         PP(i)=PP(i)*dble(N)/4d0/Lx/Ly
!         PPm(i)=PPm(i)*dble(N)/4d0/Lx/Ly
!      ENDDO

      RETURN							
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
      PARAMETER(Ntot = 4096)
      INTEGER its,iter,ITMAX
      DOUBLE PRECISION fret,ftol,EPS,ftol1
      PARAMETER (EPS=1d-10,ITMAX=1000000000)
      DOUBLE PRECISION dgg,fp,gam,gg,gx(Ntot),gy(Ntot),hx(Ntot),hy(Ntot)
      DOUBLE PRECISION D(Ntot),D1,xix(Ntot),xiy(Ntot),xith(Ntot),width
      DOUBLE PRECISION x(Ntot),y(Ntot),maxdis,xp(Ntot),yp(Ntot)
      DOUBLE PRECISION th(Ntot),hth(Ntot),gth(Ntot),exp,att,V
      INTEGER N,countn(Ntot),nl(800,Ntot)

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
      IF(fp.lt.ftol*dble(N).and.att.eq.0d0) THEN
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

         IF(att.eq.0d0) THEN
            IF(dabs(fret-fp).lt.ftol1*fp.or.fret.lt.ftol*dble(N))THEN
                CALL func(N,x,y,th,D,D1,fp,countn,nl)
               RETURN
            ENDif
         else
            IF(dabs(fret-fp).lt.ftol1) THEN
                CALL func(N,x,y,th,D,D1,fp,countn,nl)
               RETURN
            ENDIF 
         ENDif
         
         CALL CG_check(N,x,y,xp,yp,maxdis)	     
         IF(maxdis.gt.width*D1) THEN
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
         
         IF(gg.eq.0d0) THEN
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
      IF(fb.gt.fa)THEN ! was gt
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      ENDif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     IF(fb.gt.fc)THEN ! was ge
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2d0*dsign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        IF((bx-u)*(u-cx).gt.0d0)THEN
          fu=func(u)
          IF(fu.lt.fc)THEN
            ax=bx
            fa=fb
            bx=u
            fb=fu
            RETURN
          else IF(fu.gt.fb)THEN
            cx=u
            fc=fu
            RETURN
          ENDif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else IF((cx-u)*(u-ulim).gt.0d0)THEN
          fu=func(u)
          IF(fu.lt.fc)THEN
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          ENDif
        else IF((u-ulim)*(ulim-cx).ge.0d0)THEN ! was ge
          u=ulim
          fu=func(u)
        else
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
      IF(maxdis.gt.width*D1com) THEN
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
        IF(dabs(x-xm).le.(tol2-0.5d0*(b-a))) goto 3
        IF(dabs(e).gt.tol1) THEN
          d1=2.*(b-a)
          d2=d1
          IF(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          IF(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0d0).and.(dx*d1.le.0d0)
          ok2=((a-u2)*(u2-b).gt.0d0).and.(dx*d2.le.0d0)
          olde=e
          e=d
          IF(.not.(ok1.or.ok2))THEN
            goto 1
          else if (ok1.and.ok2)THEN
            IF(dabs(d1).lt.dabs(d2))THEN
              d=d1
            else
              d=d2
            ENDif
          else if (ok1)THEN
            d=d1
          else
            d=d2
          ENDif
          IF(dabs(d).gt.dabs(0.5d0*olde))goto 1
          u=x+d
          IF(u-a.lt.tol2 .or. b-u.lt.tol2) d=dsign(tol1,xm-x)
          goto 2
        ENDif
1       IF(dx.ge.0d0) THEN
          e=a-x
        else
          e=b-x
        ENDif
        d=0.5d0*e
2        IF(dabs(d).ge.tol1) THEN
          u=x+d
           fu=f(u)
        else
          u=x+dsign(tol1,d)
          fu=f(u)
          IF(fu.gt.fx)goto 3
        ENDif
        du=df(u)
        IF(fu.le.fx) THEN
          IF(u.ge.x) THEN
            a=x
          else
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
        else
          IF(u.lt.x) THEN
            a=u
          else
            b=u
          ENDif
          IF(fu.le.fw .or. w.eq.x) THEN
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else IF(fu.le.fv .or. v.eq.x .or. v.eq.w) THEN
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
      IF(maxdis.gt.width*D1com) THEN
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
      if (idum.le.0) THEN
        idum=max(-idum,1)
        idum2=idum
        DO 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      ENDif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      IF(iy.lt.1)iy=iy+IMM1
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
      if (fret.lt.ftol*N .or. phi.lt.0.25) THEN
          terminate = 0
          before_jamming = 1
          at_jamming = 0
          above_jamming = 0
          dtstatus = 0
          
          dt = 1d0
          phi_j = 0d0
          RETURN
      ENDif
      
      
      if (above_jamming.eq.0 .and. at_jamming .eq. 0)THEN
          if (fret.gt.(ftol * wide * N) )  THEN
             dt = 0.5*dt
             dtstatus = 1
             terminate = 0
             RETURN 
          else if (fret.gt.(ftol*N) .and. fret.lt.(ftol*wide*N))THEN
              before_jamming = 0
              at_jamming = 1
              above_jamming = 1
              dt = 1d0
              dtstatus = 0
              terminate = 0
              phi_j = phi
              RETURN
          ! else - not needed. this point should never been reached!
              
        ENDif
      ENDif
      
      
      if (phi_j .lt. 0.25) THEN
          WRITE(*,*) "ERROR: SHOULD NOT REACH THIS POINT"
      ENDif
      
      if (at_jamming.ne.0) THEN
          WRITE(*,*) "SOMETHING WRONG: should be at_jamming = 0"
      ENDif

      if (above_jamming.ne.1) THEN
          WRITE(*,*) "SOMETHING WRONG: should be above_jamming = 1"
      ENDif      
      
      at_jamming = 0
      if (phi .lt. phi_j+dphi-delta) THEN
          dt = 1
          terminate = 0
          dtstatus = 0
          RETURN
      else if (phi .gt. phi_j+dphi+delta) THEN
          dt = 0.5 * dt
          dtstatus = 1
          terminate = 0
          RETURN
      else
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
      
      
      SUBROUTINE copy_everything(Ntot, Ngen, ndiv, x, y, th, D, alpha, 
     +     rate, alphar,
     +     scale, P, PP, D0, aclone,
     +     tbirth, xbirth, ybirth, xdiv, ydiv, thdiv,
     +     N, c, Nexist, nclone, nclonebox, age, agehist, idiv,
     +     xcc, ycc, thcc, Dcc, alphacc, ratecc, alpharcc,
     +     scalecc, Pcc, PPcc, D0cc, aclonecc,
     +     tbirthcc, xbirthcc, ybirthcc, xdivcc, ydivcc, thdivcc,
     +     Ncc, ccc, Nexistcc, nclonecc, ncloneboxcc, agecc, agehistcc,
     +     idivcc, BUDCONT, BUDCONTCC,PPm,PPmcc)
! ...cc stands for copy
      implicit none
      INTEGER Ntot,Ngen
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot)
      DOUBLE PRECISION alpha(Ntot),rate(Ntot),alphar(Ntot),scale(Ntot)
      DOUBLE PRECISION P,PP(Ntot),D0(Ntot),PPm(Ntot)
      DOUBLE PRECISION aclone(Ntot*Ngen)
      DOUBLE PRECISION tbirth(Ntot*Ngen)
      DOUBLE PRECISION xbirth(Ntot*Ngen),ybirth(Ntot*Ngen)
      DOUBLE PRECISION xdiv(999,2),ydiv(999,2),thdiv(999,2)
      INTEGER N,c(Ntot,Ngen),Nexist
      INTEGER nclone(Ntot*Ngen),nclonebox(Ntot*Ngen)
      INTEGER age(Ntot),agehist(Ngen,2)
      INTEGER idiv(999,2), ndiv
      DOUBLE PRECISION xcc(Ntot),ycc(Ntot),thcc(Ntot),Dcc(Ntot)
      DOUBLE PRECISION alphacc(Ntot),ratecc(Ntot),alpharcc(Ntot)
      DOUBLE PRECISION scalecc(Ntot)
      DOUBLE PRECISION Pcc,PPcc(Ntot),D0cc(Ntot),PPmcc(Ntot)
      DOUBLE PRECISION aclonecc(Ntot*Ngen)
      DOUBLE PRECISION tbirthcc(Ntot*Ngen)
      DOUBLE PRECISION xbirthcc(Ntot*Ngen),ybirthcc(Ntot*Ngen)
      DOUBLE PRECISION xdivcc(999,2),ydivcc(999,2),thdivcc(999,2)
      INTEGER Ncc,ccc(Ntot,Ngen),Nexistcc
      INTEGER nclonecc(Ntot*Ngen),ncloneboxcc(Ntot*Ngen)
      INTEGER agecc(Ntot),agehistcc(Ngen,2)
      INTEGER idivcc(999,2), ndivcc
      INTEGER i,j, BUDCONT(Ntot), BUDCONTCC(Ntot)
      
!      
      Pcc = P
      Ncc = N
      Nexistcc = Nexist
      ndivcc = ndiv     
      
! ZERO EVERYTHING !!!
      DO i=1,Ntot
          xcc(i) = 0d0
          ycc(i) = 0d0
          thcc(i) = 0d0
          Dcc(i) = 0d0
          alphacc(i) = 0d0
          ratecc(i) = 0d0
          alpharcc(i) = 0d0
          scalecc(i) = 0d0
          PPcc(i) = 0d0
          PPmcc(i) = 0d0
          D0cc(i) = 0d0
          agecc(i) = 0d0
          BUDCONTCC(i) = 0
      ENDDO
      
      DO i=1,999
          xdivcc(i,1) = 0d0
          ydivcc(i,1) = 0d0
          thdivcc(i,1) = 0d0
          idivcc(i,1) = 0d0
          xdivcc(i,2) = 0d0
          ydivcc(i,2) = 0d0
          thdivcc(i,2) = 0d0
          idivcc(i,2) = 0d0
      ENDDO
      
      DO i=1,Ngen
          agehistcc(i,1) = 0d0
          agehistcc(i,2) = 0d0
      ENDDO
      
      DO i=1,Nexistcc+2!Ntot*Ngen
          tbirthcc(i) = 0d0
          xbirthcc(i) = 0d0
          ybirthcc(i) = 0d0
          aclonecc(i) = 0d0
          nclonecc(i) = 0d0
          ncloneboxcc(i) = 0d0
      ENDDO 
      
      DO i=1,Nexistcc+2!Ntot
          DO j=1,Nexistcc+2!Ngen
              ccc(i,j) = 0
          ENDDO
      ENDDO
     

!!!!!!!!!!!!!!!!!!!!!   MAKE COPY   
      DO i=1,Ntot
          xcc(i) = x(i)
          ycc(i) = y(i)
          thcc(i) = th(i)
          Dcc(i) = D(i)
          alphacc(i) = alpha(i)
          ratecc(i) = rate(i)
          alpharcc(i) = alphar(i)
          scalecc(i) = scale(i)
          PPcc(i) = PP(i)
          PPmcc(i) = PPm(i)
          D0cc(i) = D0(i)
          agecc(i) = age(i)
          BUDCONTCC(i) = BUDCONT(i)
      ENDDO
      
      DO i=1,999
          xdivcc(i,1) = xdiv(i,1)
          ydivcc(i,1) = ydiv(i,1)
          thdivcc(i,1) = thdiv(i,1) 
          idivcc(i,1) = idiv(i,1)
          xdivcc(i,2) = xdiv(i,2)
          ydivcc(i,2) = ydiv(i,2)
          thdivcc(i,2) = thdiv(i,2) 
          idivcc(i,2) = idiv(i,2) 
      ENDDO
      
      DO i=1,Ngen
          agehistcc(i,1) = agehist(i,1)
          agehistcc(i,2) = agehist(i,2)
      ENDDO
      
      DO i=1,Nexistcc+2!Ntot*Ngen
          tbirthcc(i) = tbirth(i)
          xbirthcc(i) = xbirth(i)
          ybirthcc(i) = ybirth(i)
          aclonecc(i) = aclone(i)          
          nclonecc(i) = nclone(i)
          ncloneboxcc(i) = nclonebox(i)
      ENDDO
      
      DO i=1,Nexistcc+2!Ntot
          DO j=1,Nexistcc+2!Ngen
              ccc(i,j) = c(i,j)
          ENDDO
      ENDDO
      END
      
      SUBROUTINE copy_back_everything(Ntot, Ngen, ndiv, x, y, th, D, 
     +     alpha, 
     +     rate, alphar,
     +     scale, P, PP, D0, aclone,
     +     tbirth, xbirth, ybirth, xdiv, ydiv, thdiv,
     +     N, c, Nexist, nclone, nclonebox, age, agehist, idiv,
     +     xcc, ycc, thcc, Dcc, alphacc, ratecc, alpharcc,
     +     scalecc, Pcc, PPcc, D0cc, aclonecc,
     +     tbirthcc, xbirthcc, ybirthcc, xdivcc, ydivcc, thdivcc,
     +     Ncc, ccc, Nexistcc, nclonecc, ncloneboxcc, agecc, agehistcc,
     +     idivcc, BUDCONT, BUDCONTCC,PPm,PPmcc)
! ...cc stands for copy
      implicit none
      INTEGER Ntot,Ngen
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot)
      DOUBLE PRECISION alpha(Ntot),rate(Ntot),alphar(Ntot),scale(Ntot)
      DOUBLE PRECISION P,PP(Ntot),D0(Ntot),PPm(Ntot)
      DOUBLE PRECISION aclone(Ntot*Ngen)
      DOUBLE PRECISION tbirth(Ntot*Ngen)
      DOUBLE PRECISION xbirth(Ntot*Ngen),ybirth(Ntot*Ngen)
      DOUBLE PRECISION xdiv(999,2),ydiv(999,2),thdiv(999,2)
      INTEGER N,c(Ntot,Ngen),Nexist
      INTEGER nclone(Ntot*Ngen),nclonebox(Ntot*Ngen)
      INTEGER age(Ntot),agehist(Ngen,2)
      INTEGER idiv(999,2), ndiv
      DOUBLE PRECISION xcc(Ntot),ycc(Ntot),thcc(Ntot),Dcc(Ntot)
      DOUBLE PRECISION alphacc(Ntot),ratecc(Ntot),alpharcc(Ntot)
      DOUBLE PRECISION scalecc(Ntot)
      DOUBLE PRECISION Pcc,PPcc(Ntot),D0cc(Ntot),PPmcc(Ntot)
      DOUBLE PRECISION aclonecc(Ntot*Ngen)
      DOUBLE PRECISION tbirthcc(Ntot*Ngen)
      DOUBLE PRECISION xbirthcc(Ntot*Ngen),ybirthcc(Ntot*Ngen)
      DOUBLE PRECISION xdivcc(999,2),ydivcc(999,2),thdivcc(999,2)
      DOUBLE PRECISION dd
      INTEGER Ncc,ccc(Ntot,Ngen),Nexistcc
      INTEGER nclonecc(Ntot*Ngen),ncloneboxcc(Ntot*Ngen)
      INTEGER agecc(Ntot),agehistcc(Ngen,2)
      INTEGER idivcc(999,2), ndivcc
      INTEGER i,j, celltype, BUDCONT(Ntot), BUDCONTCC(Ntot)
      COMMON /f10com/ celltype
            
!
      P = Pcc
      N = Ncc
      Nexist = Nexistcc
      ndiv = ndivcc
!
      DO i=1,Ntot
          x(i) = xcc(i)
          y(i) = ycc(i)
          th(i) = thcc(i)
          D(i) = Dcc(i)
          alpha(i) = alphacc(i)
          rate(i) = ratecc(i)
          alphar(i) = alpharcc(i)
          scale(i) = scalecc(i)
          PP(i) = PPcc(i)
          PPm(i) = PPmcc(i)
          D0(i) = D0cc(i)
          age(i) = agecc(i)
          BUDCONT(i) = BUDCONTCC(i)
      ENDDO
      
      DO i=1,999
          xdiv(i,1) = xdivcc(i,1)
          ydiv(i,1) = ydivcc(i,1)
          thdiv(i,1) = thdivcc(i,1) 
          idiv(i,1) = idivcc(i,1)
          xdiv(i,2) = xdivcc(i,2)
          ydiv(i,2) = ydivcc(i,2)
          thdiv(i,2) = thdivcc(i,2) 
          idiv(i,2) = idivcc(i,2) 
      ENDDO
      
      DO i=1,Ngen
          agehist(i,1) = agehistcc(i,1)
          agehist(i,2) = agehistcc(i,2)
      ENDDO
      
      DO i=1,Nexist+2!Ntot*Ngen
          tbirth(i) = tbirthcc(i)
          xbirth(i) = xbirthcc(i)
          ybirth(i) = ybirthcc(i)
          aclone(i) = aclonecc(i)          
          nclone(i) = nclonecc(i)
          nclonebox(i) = ncloneboxcc(i)
      ENDDO
      
      DO i=1,Nexist+2!Ntot
          DO j=1,Nexist+2!Ngen
              c(i,j) = ccc(i,j)
          ENDDO
      ENDDO

         DO i=1,N
            IF(celltype.eq.1) THEN
               scale(i)=dsqrt(1d0+alpha(i)**2)/4d0*d(i)         
            elseIF(celltype.eq.2) THEN
               dd=alpha(i)-1d0
               scale(i)=dsqrt(2d0*(1d0+dd**4)/(1+dd**2)+
     +              4d0*(dd*(1d0+dd)/(1+dd**2))**2)/4d0*d(i)
            elseIF(celltype.eq.3) THEN
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
         IF(F(ii).eq.1) THEN
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
      
      DO while(Rnew>0)
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
!               IF(overlap(N,x,y,th,D,D1,II,JJ).lt.1.d0) THEN
                  c(II) = c(II)+overlap(N,x,y,th,D,D1,II,JJ) !1
                  c(JJ) = c(JJ)+overlap(N,x,y,th,D,D1,II,JJ) !1
!               ENDif
            ENDDO
         ENDDO

         ! Remove floaters
         i=1
         DO while (i.le.N-R)
            II=listP(i)
            IF(c(II).lt.3) THEN
               Nf = Nf + 1
               Rnew=Rnew+1
               R=R+1
               DO j=i, N-R
                  listP(j) = listP(j+1)
               ENDDO
               F(II) = 1
            else
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
!            IF(overlap(N,x,y,th,D,D1,II,JJ).lt.1.d0) THEN
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
      DOUBLE PRECISION dijsq_up,scale(Ntot),c(Ntot)
      DOUBLE PRECISION s(Ntot),dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2),di_up(Ntot),di1j1      
      INTEGER celltype, ni, i, j
 
      
      COMMON /f2com/ width
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f9com/ scale
      COMMON /f10com/ celltype
      
      if (D1.ne.1.d0) THEN
          WRITE(*,*) "SOMETHING IS WRONG. D1 SHALL BE EQ. TO 1.0"
      ENDif
      
      overlap = 0.d0

      ! BUDDING CELLS
      if (celltype.eq.2) THEN
          WRITE(*,*) "THERE IS SPECIAL CONTACT FUNCTION FOR YEAST"

      ! DISK CELLS      
      elseif (celltype.eq.3) THEN
                
          xij=x(i)-x(j)
          xij=xij-idnint(xij/Lx)*Lx  !! PBC
          yij=y(i)-y(j)
          yij=yij-idnint(yij/Ly)*Ly  !! PBC
          
          rij=dsqrt(xij**2+yij**2)  
          dij=D1*(d(i)+d(j))/2d0

          IF(rij.lt.dij) THEN
              overlap = 1.d0
          else
              overlap = 0.d0
          ENDIF
      ENDif
          
      END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! inuse
      SUBROUTINE contacts_yeast(x,y,th,D1,D,N,Z,F,Nf,Nu,Nmm,Nbb,Nmb)
      IMPLICIT NONE
      INTEGER Ntot, N 
      PARAMETER(Ntot = 4096)
      INTEGER i, j, F(Ntot), Z, NCBUD(Ntot,2), Nf, Nu, Nmm, Nbb, Nmb
      INTEGER ki,kj,k
      DOUBLE PRECISION overlap, aspect_ratio
      DOUBLE PRECISION x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      DOUBLE PRECISION xij, yij, D(Ntot), D1
      DOUBLE PRECISION exp,dij_up,dij
      DOUBLE PRECISION Lx,Ly,rijsq,c(Ntot),att
      DOUBLE PRECISION s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2)
      INTEGER flag
      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
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
         NCBUD(i,1)=0
         NCBUD(i,2)=0
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
                      IF(rijsq.lt.(dij**2)) THEN
                          Z = Z+2
                          NCBUD(i,ki)=NCBUD(i,ki)+1
                          NCBUD(j,kj)=NCBUD(j,kj)+1
                          if (ki.eq.1 .and. kj.eq.1) THEN
                              Nmm = Nmm + 2
                          elseif (ki.eq.2 .and. kj.eq.2) THEN
                              Nbb = Nbb + 2
                          else
                              Nmb = Nmb + 2
                          ENDif
                      ENDif
                  ENDDO
              ENDDO
          ENDDO
      ENDDO

      
      DO i=1,N
          flag = 0
          if ( (NCBUD(i,1)+NCBUD(i,2)).lt.3 ) THEN
              nf = nf + 1
              flag = 1
          ENDif
          
          if ( (NCBUD(i,1)+NCBUD(i,2)).eq.3 ) THEN
              if ( NCBUD(i,1).eq.2 .or. NCBUD(i,2).eq.2 ) THEN
                  nf = nf + 1
                  flag = 1
              ENDif
          ENDif
          
          if (flag.eq.0) THEN    
              if ( NCBUD(i,1).eq.0 ) THEN
                  nu = nu + 1
              ENDif
              if ( NCBUD(i,2).eq.0 ) THEN
                  nu = nu + 1
              ENDif
          ENDif
          
      ENDDO
      
      END
      
      SUBROUTINE out_numbers(N, Nf, Nu, Ziso)
      implicit none
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      INTEGER N, Nf, Nu, Ziso
      
      Ziso = 6*(N-Nf) - 2*Nu - 2
      
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      SUBROUTINE bud_contacts(x,y,th,D1,D,N,BUDCONT,ZEROBUDS)
      IMPLICIT NONE
      INTEGER Ntot, N
      PARAMETER(Ntot = 4096)
      INTEGER i, j, NCBUD(Ntot,2), ZEROBUDS, BUDCONT(Ntot),k,ki,kj
      DOUBLE PRECISION x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      DOUBLE PRECISION xij, yij, D(Ntot), D1
      DOUBLE PRECISION exp,dij
      DOUBLE PRECISION Lx,Ly,rijsq,c(Ntot),att 
      DOUBLE PRECISION s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dk(Ntot,2)

      COMMON /f3com/ alpha ! aspect ratio
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly

      ZEROBUDS = 0
      
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
         NCBUD(i,1)=0
         NCBUD(i,2)=0
         BUDCONT(i) = 0
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
                      IF(rijsq.lt.(dij**2)) THEN
                          NCBUD(i,ki)= NCBUD(i,ki)+1
                          NCBUD(j,kj)= NCBUD(j,kj)+1
                      ENDif
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      
      DO i=1,N
          BUDCONT(i) = NCBUD(i,2)
          if (BUDCONT(i).eq.0)THEN
              ZEROBUDS = ZEROBUDS + 1
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
      DOUBLE PRECISION total_area, total_newmass
      
      gr = 0.0d0
      total_area = 0.0d0
      total_newmass = 0.0d0
      
      DO i=1,N
          ratei=rate(i)
          IF(P0.gt.0d0.and.PP(i).gt.0d0) THEN
              ratei=ratei*dexp(-PP(i)/P0)
          ENDif
          
          Dm = D(i)
          Db = alpha(i)-1d0
          total_area = total_area + Dm*Dm + Db*Db
          total_newmass = total_newmass + ratei*(Dm*Dm + Db*Db)
      ENDDO
      
      gr = total_newmass / total_area
      
      END
