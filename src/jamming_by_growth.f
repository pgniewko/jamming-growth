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
      integer Ntot,Ngen
      PARAMETER(Ntot=4096,Ngen=4096)
      double precision pi
      PARAMETER(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,exp,ran2
      double precision ftol,ftol1,fret,alpha0,width,Lx,Ly,ratei
      double precision alpha(Ntot),rate(Ntot),alphar(Ntot),scale(Ntot)
      double precision rate0,desync,phi,flow,tdiv,P,PP(Ntot),D0(Ntot)
      double precision aclone(Ntot*Ngen),dispcm,xa(2),ya(2),PR,PT,P0
      double precision tbirth(Ntot*Ngen),cc,ss,dr(2),dd,corr,att,rat
      double precision xbirth(Ntot*Ngen),ybirth(Ntot*Ngen)
      double precision xdiv(999,2),ydiv(999,2),thdiv(999,2)
      double precision dt, PPm(Ntot), total_growthrate
      integer dtstatus, terminate
      integer N,Nr,seed,iter,steps,i,j,k,kk,c(Ntot,Ngen),m,skip,Nexist
      integer celltype,divtype,nclone(Ntot*Ngen),nclonebox(Ntot*Ngen)
      integer age(Ntot),agehist(Ngen,2),agetot1,agetot2
      integer div,ndiv,idiv(999,2),kdiv, Nf, Nu, Nmm,Nbb,Nmb
      character file1*80
      
! NEW DATA 
      character file_LF_JAMM*120
      character file_LF_DPHI*120
      character STATS_file_LF_JAMM*120
      character STATS_file_LF_DPHI*120
      character file_NC*100
      
      double precision xcc(Ntot),ycc(Ntot),thcc(Ntot),Dcc(Ntot)
      double precision alphacc(Ntot),ratecc(Ntot),alpharcc(Ntot)
      double precision scalecc(Ntot)
      double precision Pcc,PPcc(Ntot),D0cc(Ntot),PPmcc(Ntot)
      double precision aclonecc(Ntot*Ngen)
      double precision tbirthcc(Ntot*Ngen)
      double precision xbirthcc(Ntot*Ngen),ybirthcc(Ntot*Ngen)
      double precision xdivcc(999,2),ydivcc(999,2),thdivcc(999,2)
      integer Ncc,ccc(Ntot,Ngen),Nexistcc
      integer nclonecc(Ntot*Ngen),ncloneboxcc(Ntot*Ngen)
      integer agecc(Ntot),agehistcc(Ngen,2)
      integer idivcc(999,2), ndivcc
      integer F(Ntot), Nc, Ziso, F_e(Ntot)
      integer BUDCONT(Ntot), BUDCONTCC(Ntot), ZEROBUDS
      double precision phitemp, calc_phi, wide
      integer N_j,Ziso_j,Nc_j,Nf_j,Nu_j,Nmm_j,Nbb_j,Nmb_j
      double precision phi_j, dphi
      
      integer before_jamming, at_jamming, above_jamming
      
! END OF NEW DATA      
      
      common /f2com/ width
      common /f3com/ alpha
      common /f4com/ exp,att
      common /f5com/ Lx,Ly
      common /f6com/ P,PP,PT,PR,PPm
      common /f8com/ alpha0
      common /f9com/ scale
      common /f10com/ celltype

      ! read geometric parameters
      read(*,*) alpha0
      read(*,*) Lx
      read(*,*) Ly

      ! read cell parameters
      read(*,*) celltype
      read(*,*) divtype
      read(*,*) P0
      read(*,*) att

      ! read run parameters
      read(*,*) rate0
      read(*,*) desync
      read(*,*) steps
      read(*,*) seed
      read(*,*) skip

      ! read pdhi exponent
      read(*,*) dphi
      
      ! read output files
      read(*,*) file1

      if(celltype.eq.1 .or. celltype.eq.3) then
           write(*,*) "ELLIPSE AND DISK CELL TYPE IS NOT SUPPORTED"
           call EXIT(0)
      endif
     
      if(P0.eq.0d0) then
           write(*,*) "P0 = 0 not supported"
           call EXIT(0)
      endif
      
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
 123  open(unit=1,file=TRIM(file1), status='replace') ! CONFIGURATION FILE
      
      ! JAMMED FILES    
      file_LF_JAMM='LF_JAMM_' // TRIM(file1)            ! unit=11
      STATS_file_LF_JAMM='STATS_LF_JAMM_' // TRIM(file1)! unit=21
      
      
      ! DPHI OUTPUT FILES
      file_LF_DPHI='LF_DPHI_' // TRIM(file1)
      open(unit=12,file=TRIM(file_LF_DPHI), status='replace')
      STATS_file_LF_DPHI='STATS_LF_DPHI_' // TRIM(file1)
      open(unit=22,file=TRIM(STATS_file_LF_DPHI), status='replace')
      
      ! SOME STATS
      file_NC='NC_' // TRIM(file1)
      open(unit=13,file=TRIM(file_NC), status='replace')
      
      ! initial size & aspect ratios - total vol = (1+rat)*alpha0
      d(1)=1d0
      d(2)=1d0
      alpha(1)=alpha0
      alpha(2)=dsqrt(alpha0*(1d0+rat)-(alpha0-1d0)**2-2d0)+1d0
     

      ! random initial config
      N=2
      Nexist=2
      do i=1,N 
         c(i,1)=1
         c(i,2)=i
         x(i)=Lx/2 + 0d0
         y(i)=Ly/2 + (dble(i)-1.5d0)*d(i)*D1
         th(i)=(ran2(seed)-0.5d0)*2d0*pi
         D0(i)=D1
         rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
         age(i)=0
      enddo

      ! calculate # steps until division
      tdiv=dlog10(2d0)/dlog10(1d0+rate0)

      write(*,*) steps, tdiv, 2*int(steps/tdiv)

      ! set histogram of cell ages to 0
      do k=1,2*int(steps/tdiv)
         agehist(k,1)=0
         agehist(k,2)=0
      enddo
      
      dt = 1.0
      dtstatus = 0
      terminate = 0
      ! loop over division
      k=0
      do while (terminate.ne.1)
      k=k+1
      
      call copy_everything(Ntot,Ngen,ndiv,x,y,th,D,alpha,rate,
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
         do i=1,N
               ratei=rate(i)*dt
               
               if(P0.gt.0d0.and.PP(i).gt.0d0) then
                  ratei=ratei*dexp(-PP(i)/P0) !* dt
               endif
               
               alpha(i)=1d0+dsqrt((1d0+ratei)*
     +                 (1d0+(alpha(i)-1d0)**2)-1d0)


               if(alpha(i).gt.2d0*alpha0) then
                  dispcm=alpha0/2d0
                  div=1
                  
                  ! divide into 2 - 1st assigned index N+1
                  N=N+1
                  c(N,1)=c(i,1)+1
                  do j=2,c(N,1)
                     c(N,j)=c(i,j)
                  enddo
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
                  do j=2,c(i,1)
                     nclonebox(c(i,j))=nclonebox(c(i,j))+1
                  enddo

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
               endif
               
               ! types of division: ->->, <-->, -><-, random
               if(div.eq.1) then
                  if(divtype.eq.1) then
                     th(N)=th(i)
                  elseif(divtype.eq.2) then
                     th(N)=th(i)+pi
                  elseif(divtype.eq.3) then
                     th(N)=th(i)
                     th(i)=th(i)+pi
                  elseif(divtype.eq.4) then
                     th(N)=(ran2(seed)-0.5d0)*2d0*pi
                     th(i)=(ran2(seed)-0.5d0)*2d0*pi
                  endif

                  ! temp, remove
                  th(i)=th(i) + 1d-4*(ran2(seed)-0.5d0)
                  th(N)=th(N) + 1d-4*(ran2(seed)-0.5d0)
               endif

               if(div.eq.1) then
                  ndiv=ndiv+1
                  idiv(ndiv,1)=i
                  idiv(ndiv,2)=N
                  xdiv(ndiv,1)=x(i)
                  xdiv(ndiv,2)=x(N)
                  ydiv(ndiv,1)=y(i)
                  ydiv(ndiv,2)=y(N)
                  thdiv(ndiv,1)=th(i)
                  thdiv(ndiv,2)=th(N)
               endif
         enddo
         
         
         ! convert from angle to length scale = sqrt(I/m) * angle
         do i=1,N
            dd=alpha(i)-1d0
            scale(i)=dsqrt(2d0*(1d0+dd**4)/(1+dd**2)+
     +              4d0*(dd*(1d0+dd)/(1+dd**2))**2)/4d0*d(i)
     
            th(i)=th(i)*scale(i)
         enddo
           
         ! minimize energy
         call frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
        
         
         phitemp = calc_phi(D, alpha, D1, N)
         
         write(*,*) k, N,Nc,fret/dble(N), P, dt, phitemp,
     +     before_jamming, at_jamming, above_jamming   

         call determine_system_status(N,D,D1,alpha,ftol,wide,dt,
     +     fret, phi_j, dphi,
     +     dtstatus, terminate,
     +     before_jamming, at_jamming, above_jamming)
     
         
         ! REJECT THE MOVE
         if(dtstatus.eq.1) then
         call copy_back_everything(Ntot, Ngen,ndiv,x,y,th,D,alpha,
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
         end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

!         phitemp = calc_phi(D, alpha, D1, N)
!         write(*,*) "*",k, N,fret/dble(N), P, dt, phitemp,
!     +     before_jamming, at_jamming, above_jamming
         
         ! convert back to angles
         do i=1,N
            th(i)=th(i)/scale(i)
         enddo       
         
         if(mod(k,skip).eq.0) then
            phi = calc_phi(D, alpha, D1, N) 
            write(1,*) 2*N, phi
            do i=1,N              
               cc=dcos(th(i))
               ss=dsin(th(i))
               dd=alpha(i)-1d0
               dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
               dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
               do kk=1,2
                  xa(kk)=x(i)+dr(kk)*cc
                  ya(kk)=y(i)+dr(kk)*ss
               enddo
               write(1,'(3E26.18,I12)')xa(1),ya(1),d(i),0
               write(1,'(3E26.18,I12)')xa(2),ya(2),d(i)*dd,1
            enddo
            flush(1) 
         endif
         
         

         if ( at_jamming.eq.1 ) then
            open(unit=11,file=TRIM(file_LF_JAMM), status='replace')
      open(unit=21,file=TRIM(STATS_file_LF_JAMM), status='replace')
            phi = calc_phi(D, alpha, D1, N)
            
            call growth_rate(total_growthrate,N,rate,PP,P0,D,alpha)
            
            write(1,*) 2*N, phi
            write(11,*) N, phi
            write(21,*) 2*N, phi, total_growthrate
            do i=1,N
                write(11,'(5E26.18)') x(i),y(i),D(i),alpha(i),th(i)
                
                cc=dcos(th(i))
                ss=dsin(th(i))
                dd=alpha(i)-1d0
                dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
                dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
                do kk=1,2
                    xa(kk)=x(i)+dr(kk)*cc
                    ya(kk)=y(i)+dr(kk)*ss
                enddo
                write(1,'(3E26.18,I12)')xa(1),ya(1),d(i),0     ! TRA
                write(1,'(3E26.18,I12)')xa(2),ya(2),d(i)*dd,1  ! TRA
                
                write(21,'(3E26.18,I12,4E26.18)')xa(1),ya(1),d(i),0,
     +     PPm(i),rate(i)*dexp(-PP(i)/P0),rate(i),P0
                write(21,'(3E26.18,I12,4E26.18)')xa(2),ya(2),d(i)*dd,1,
     +     PP(i),rate(i)*dexp(-PP(i)/P0),rate(i),P0
            enddo
            flush(1)
            flush(11)
            close(11)
            flush(21)
            close(21)
      
            close(13)
            open(unit=13,file=TRIM(file_NC), status='replace')
            call contacts_yeast(x,y,th,D1,D,N,Nc,F,Nf,Nu,Nmm,Nbb,Nmb)
            call out_numbers(N, Nf, Nu, Ziso)

           write(13,'(8I8,5E26.18)')N,Ziso,Nc,Nf,Nu,Nmm,Nbb,Nmb,
     +     phi,P,fret,P0,total_growthrate
            flush(13)
            above_jamming = 1
            at_jamming = 0
         endif

         
      enddo
!!!!!! THE MAIN LOOP IS OVER
      
      
      phi = calc_phi(D, alpha, D1, N)
      call growth_rate(total_growthrate,N,rate,PP,P0,D,alpha)   
      ! save the last configuration to the file
      write(1,*) 2*N, phi
      write(12,*) N, phi
      write(22,*) 2*N, phi, total_growthrate
      do i=1,N
          write(12,'(5E26.18)') x(i),y(i),D(i),alpha(i),th(i)
                    
          cc=dcos(th(i))
          ss=dsin(th(i))
          dd=alpha(i)-1d0
          dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
          dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
          do kk=1,2
              xa(kk)=x(i)+dr(kk)*cc
              ya(kk)=y(i)+dr(kk)*ss
          enddo
          
          write(1,'(3E26.18,I12)')xa(1),ya(1),d(i),0     ! TRA
          write(1,'(3E26.18,I12)')xa(2),ya(2),d(i)*dd,1  ! TRA
          

          write(22,'(3E26.18,I12,4E26.18)')xa(1),ya(1),d(i),0,
     +     PPm(i),rate(i)*dexp(-PP(i)/P0),rate(i),P0
          write(22,'(3E26.18,I12,4E26.18)')xa(2),ya(2),d(i)*dd,1,
     +     PP(i),rate(i)*dexp(-PP(i)/P0),rate(i),P0
      enddo
      flush(1)
      flush(12)
      flush(22)
      
!!!!!!!!!!!!!

      call contacts_yeast(x,y,th,D1,D,N,Nc,F,Nf,Nu,Nmm,Nbb,Nmb)
      call out_numbers(N, Nf, Nu, Ziso)

      write(13,'(8I8,5E26.18)')N,Ziso,Nc,Nf,Nu,Nmm,Nbb,Nmb,
     +     phi,P,fret,P0,total_growthrate
      flush(13)
      
      close(1)
      close(12)
      close(13)
      close(22)
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CG_check(N,x,y,xp,yp,maxdis)
      PARAMETER(Ntot = 4096)
      double precision maxdis,x(Ntot),y(Ntot),xp(Ntot),yp(Ntot)
      integer N

      maxdis=0d0
      do i=1,N
         maxdis=max(dabs(x(i)-xp(i)),maxdis)
         maxdis=max(dabs(y(i)-yp(i)),maxdis)
      enddo
      maxdis=2d0*dsqrt(2d0*maxdis*maxdis)

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makelist(N,x,y,D,D1,xp,yp,countn,nl)
      PARAMETER(Ntot = 4096)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot),D1
      double precision Lx,Ly
      integer countn(Ntot),nl(800,Ntot),N,celltype
      common /f10com/ celltype
      common /f5com/ Lx,Ly
      
      call makelist_dimer(N,x,y,D,D1,xp,yp,countn,nl)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makelist_dimer(N,x,y,D,D1,xp,yp,countn,nl) 
      PARAMETER(Ntot = 4096)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot)
      double precision D1,xij,yij,rij,dij,rijsq,alpha(Ntot),width
      double precision dd,dr1,dr2,dk2,di_up(Ntot),exp,att
      double precision Lx,Ly
      integer countn(Ntot),nl(800,Ntot),N
      common /f2com/ width
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly

      do i=1,N
         countn(i)=0
      enddo

      do i=1,N
         dd=alpha(i)-1d0
         dr1=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
         dr2=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
         dk2=dd*D(i)
         di_up(i)=(dk2/2d0-dr2)*2d0
      enddo

      do i=1,N-1
         do j=i+1,N
            xij=x(i)-x(j)
            xij=xij-idnint(xij/Lx)*Lx  !! PBC
            yij=y(i)-y(j)
            yij=yij-idnint(yij/Ly)*Ly  !! PBC
            rijsq=xij*xij+yij*yij
            dij=(di_up(i)+di_up(j))/2d0
            dij=D1*( D(i) + D(j) )/2d0   
!            if(rijsq.lt.(2.D0*dij+(att+width)*D1)**2) then
            if(rijsq.lt.(2.D0*dij)**2) then
               countn(i)=countn(i)+1
               nl(countn(i),i)=j
            end if
         enddo
      enddo
      
      do i=1,N
         xp(i)=x(i)
         yp(i)=y(i)
      enddo
      return
      end
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine func(N,x,y,th,D,D1,V,countn,nl)
      PARAMETER(Ntot = 4096)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,V
      integer countn(Ntot),nl(800,Ntot),N,celltype
      common /f10com/ celltype


      call func_dimer(N,x,y,th,D,D1,V,countn,nl)


      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dfunc(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      PARAMETER(Ntot = 4096)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1
      double precision fx(Ntot),fy(Ntot),fth(Ntot)
      integer countn(Ntot),nl(800,Ntot),N,celltype
      common /f10com/ celltype


      call dfunc_dimer(N,x,y,th,D,D1,fx,fy,fth,countn,nl)



      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine func_dimer(N,x,y,th,D,D1,V,countn,nl)
      PARAMETER(Ntot = 4096)
      double precision pi
      PARAMETER(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,V,alpha(Ntot)
      double precision rij,xij,yij,dij,exp,dlnsig,dij_up,sigma,LJ
      double precision Lx,Ly,rijsq,dijsq_up,scale(Ntot),c(Ntot),att
      double precision s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2),Vij
      double precision dk(Ntot,2),di_up(Ntot),di1j1,xhit,yhit,yhitout
      integer countn(Ntot),nl(800,Ntot),N
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly
      common /f9com/ scale

      ! convert from molecules to atoms
      do i=1,N
         c(i)=dcos(th(i)/scale(i))
         s(i)=dsin(th(i)/scale(i))
         dd=alpha(i)-1d0
         dr(i,1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)*D1/2d0
         dr(i,2)=-(1d0+dd)/(1d0+dd**2)*D(i)*D1/2d0
         do k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         enddo
         dk(i,1)=D(i)*D1
         dk(i,2)=dd*D(i)*D1
         if(alpha(i).lt.2d0) then 
            di_up(i)=(dk(i,2)/2d0-dr(i,2))*2d0*D1
         else
            di_up(i)=(dk(i,1)/2d0+dr(i,1))*2d0*D1         
         endif
      enddo

      ! inter-particle interactions
      V=0d0
      do i=1,N-1
         if(countn(i).ge.1) then
            do jj=1,countn(i)
               j=nl(jj,i)
               dij_up=(di_up(i)+di_up(j))/2d0
               xij=x(i)-x(j)
               xij=xij-idnint(xij/Lx)*Lx  !! PBC
               if(dabs(xij).lt.dij_up+att) then
                  yij=y(i)-y(j)
                  yij=yij-idnint(yij/Ly)*Ly !! PBC 
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij_up+att)**2) then
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     do ki=1,2
                        do kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           xij=xij-idnint(xij/Lx)*Lx  !! PBC
                           yij=ya(i,ki)-ya(j,kj)
                           yij=yij-idnint(yij/Ly)*Ly !! PBC 
                           rijsq=xij**2+yij**2
                           if(rijsq.lt.(dij+att)**2) then
                              rij=dsqrt(rijsq)
                              if(exp .gt. 2.9) then
                                 LJ=(dij/rij)*(dij/rij)
                                 LJ=LJ*LJ*LJ
                                 Vij=(LJ-1d0)*(LJ-1d0)
                              else
                                 Vij=(1d0-rij/dij)**exp/exp-
     +                                (att/dij)**exp/exp
                              endif 
                              V=V+Vij*dij**2/di1j1**2
                           endif
                        enddo
                     enddo
                  end if
               end if
            enddo
         end if
      enddo

      if(exp.gt.2.9) then
         V=V/72d0
      endif

      return				
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dfunc_dimer(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      PARAMETER(Ntot = 4096)
      double precision pi
      PARAMETER(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),sigma,D(Ntot),D1,dij
      double precision fx(Ntot),fy(Ntot),fth(Ntot),rij,xij,yij,fr,exp
      double precision dij_up,alpha(Ntot),LJ,fc,ft,f_x,f_y,scale(Ntot)
      double precision fthi,fthj,fth_c,Lx,Ly,P,Pij,rijsq,dijsq_up
      double precision s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2),att
      double precision dk(Ntot,2),di_up(Ntot),di1j1,xhit,yhit,yhitout
      double precision PP(Ntot),c(Ntot),Vij,PT,PR,PPm(Ntot)
      double precision fcontact
      integer countn(Ntot),nl(800,Ntot),N !,growth_flag
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly
      common /f6com/ P,PP,PT,PR,PPm
      common /f9com/ scale
!      common /f11com/ growth_flag

      do i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
         PP(i)=0d0
         PPm(i)=0d0
      enddo
      P=0d0
      fcontact=0d0 

      ! convert to from molecules to atoms
      do i=1,N
         c(i)=dcos(th(i)/scale(i))
         s(i)=dsin(th(i)/scale(i))
         dd=alpha(i)-1d0
         dr(i,1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)*D1/2d0
         dr(i,2)=-(1d0+dd)/(1d0+dd**2)*D(i)*D1/2d0
         do k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         enddo
         dk(i,1)=D(i)*D1
         dk(i,2)=dd*D(i)*D1
         if(alpha(i).lt.2d0) then 
            di_up(i)=(dk(i,2)/2d0-dr(i,2))*2d0*D1
         else
            di_up(i)=(dk(i,1)/2d0+dr(i,1))*2d0*D1         
         endif
      enddo

      ! inter-particle interactions
      do i=1,N-1
         if(countn(i).ge.1) then
            do jj=1,countn(i)
               j=nl(jj,i)
               dij_up=(di_up(i)+di_up(j))/2d0  
               xij=x(i)-x(j)
               xij=xij-idnint(xij/Lx)*Lx  !! PBC
               if(dabs(xij).lt.dij_up+att) then
                  yij=y(i)-y(j)
                  yij=yij-idnint(yij/Ly)*Ly !! PBC
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij_up+att)**2) then
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     do ki=1,2
                        do kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           xij=xij-idnint(xij/Lx)*Lx  !! PBC
                           yij=ya(i,ki)-ya(j,kj)
                           yij=yij-idnint(yij/Ly)*Ly !! PBC
                           rijsq=xij**2+yij**2
                           if(rijsq.lt.(dij+att)**2) then
                              rij=dsqrt(rijsq)
                              if(exp .gt. 2.9) then
                                 LJ=(dij/rij)*(dij/rij)
                                 LJ=LJ*LJ*LJ
                                 fc=1d0/rij*LJ*(LJ-1d0)
                              else
                                 fc=(1d0-rij/dij)**(exp-1d0)/dij
                              endif
                              
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
                              if (ki.eq.2) then
                                  PP(i)=PP(i)+fcontact/(pi*dk(i,2)) !Pij
                              endif
                              if (kj.eq.2) then
                                  PP(j)=PP(j)+fcontact/(pi*dk(j,2)) !Pij
                              endif
                              if (ki.eq.1) then
                                  PPm(i)=PPm(i)+fcontact/(pi*dk(i,1)) !Pij
                              endif
                              if (kj.eq.1) then
                                  PPm(j)=PPm(j)+fcontact/(pi*dk(j,1)) !Pij
                              endif
                              
                           endif
                        enddo
                     enddo
                  end if
               end if
            enddo
         end if
      enddo

      
      
      if(exp .gt. 2.9) then
         do i=1,N
            fx(i)=fx(i)/6d0
            fy(i)=fy(i)/6d0
            fth(i)=fth(i)/6d0 
         enddo
         P=P/6d0
      endif
      
      do i=1,N
         fth(i)=fth(i)/scale(i)
      enddo
      

      PT=PT/Lx
      PR=PR/Ly
      P=P/4d0/Lx/Ly
!      do i=1,N
!         PP(i)=PP(i)*dble(N)/4d0/Lx/Ly
!         PPm(i)=PPm(i)*dble(N)/4d0/Lx/Ly
!      enddo

      return							
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
      PARAMETER(Ntot = 4096)
      integer its,iter,ITMAX
      double precision fret,ftol,EPS,ftol1
      PARAMETER (EPS=1d-10,ITMAX=1000000000)
      double precision dgg,fp,gam,gg,gx(Ntot),gy(Ntot),hx(Ntot),hy(Ntot)
      double precision D(Ntot),D1,xix(Ntot),xiy(Ntot),xith(Ntot),width
      double precision x(Ntot),y(Ntot),maxdis,xp(Ntot),yp(Ntot)
      double precision th(Ntot),hth(Ntot),gth(Ntot),exp,att,V
      integer N,countn(Ntot),nl(800,Ntot)

      ! not needed
      double precision f1,f2,f3,fxe,fye,fthe
      double precision xi,yi,thi,max1,max2,max3,del
      double precision alpha(Ntot),Lx,Ly,alpha0

      ! not needed
      common /f3com/ alpha ! aspect ratio
      common /f5com/ Lx,Ly
      common /f8com/ alpha0

      common /f2com/ width      
      common /f4com/ exp,att

      iter=0

      call makelist(N,x,y,D,D1,xp,yp,countn,nl)
      call func(N,x,y,th,D,D1,fp,countn,nl)
      if (fp.lt.ftol*dble(N).and.att.eq.0d0) then
         fret=fp 
         return
      endif

      call dfunc(N,x,y,th,D,D1,xix,xiy,xith,countn,nl)

      do i=1,N
        gx(i)=-xix(i)
	gy(i)=-xiy(i)
        gth(i)=-xith(i)
        hx(i)=gx(i)
	hy(i)=gy(i)
        hth(i)=gth(i)
        xix(i)=hx(i)
	xiy(i)=hy(i)
        xith(i)=hth(i)
      enddo

      do its=1,ITMAX
         iter=its

         call linmin(N,x,y,th,D,D1,xix,xiy,xith,fret,xp,yp,countn,nl)
         
c         write(*,*) its, fret/dble(N), fret, fp

         if(att.eq.0d0) then
            if(dabs(fret-fp).lt.ftol1*fp.or.fret.lt.ftol*dble(N))then
                call func(N,x,y,th,D,D1,fp,countn,nl)
               return
            endif
         else
            if(dabs(fret-fp).lt.ftol1) then
                call func(N,x,y,th,D,D1,fp,countn,nl)
               return
            end if 
         endif
         
         call CG_check(N,x,y,xp,yp,maxdis)	     
         if(maxdis.gt.width*D1) then
            call makelist(N,x,y,D,D1,xp,yp,countn,nl)
         end if
         
         fp=fret
         call dfunc(N,x,y,th,D,D1,xix,xiy,xith,countn,nl)
         
         gg=0d0
         dgg=0d0
         
         do i=1,N
            gg=gg+gx(i)*gx(i)+gy(i)*gy(i)+gth(i)*gth(i)
            dgg=dgg+(xix(i)+gx(i))*xix(i)+(xiy(i)+gy(i))*xiy(i)
     +           +(xith(i)+gth(i))*xith(i)
         enddo
         
         if(gg.eq.0d0) then
            return
         endif
         gam=dgg/gg
         do i=1,N
            gx(i)=-xix(i)
            gy(i)=-xiy(i)
            gth(i)=-xith(i)
            hx(i)=gx(i)+gam*hx(i)
            hy(i)=gy(i)+gam*hy(i)
            hth(i)=gth(i)+gam*hth(i)
            xix(i)=hx(i)
            xiy(i)=hy(i)
            xith(i)=hth(i)
         enddo
      enddo
      
c     pause 'frprmn maximum iterations exceeded'
      return
      END
C (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE linmin(N,x,y,th,D,D1,xix,xiy,xith,fret,
     +     xpr,ypr,countnr,nlr)
      double precision fret,TOL
      PARAMETER (Ntot=4096,TOL=1d-8)
C     USES dbrent,df1dim,mnbrak
      double precision ax,bx,fa,fb,fx,xmin,xx,xp(Ntot),yp(Ntot),dbrent
      double precision pxcom(Ntot),pycom(Ntot),xixcom(Ntot),xiycom(Ntot)
      double precision x(Ntot),y(Ntot),th(Ntot),xix(Ntot),xiy(Ntot)
      double precision xith(Ntot),Dcom(Ntot),D1com,D(Ntot),D1,width
      double precision pthcom(Ntot),xithcom(Ntot),f1dim,df1dim
      double precision xpr(Ntot),ypr(Ntot)
      integer countn(Ntot),nl(800,Ntot),countnr(Ntot),nlr(800,Ntot)
      integer N,ncom

      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      common /f2com/ width

      EXTERNAL df1dim
      EXTERNAL f1dim

      do i=1,N
        pxcom(i)=x(i)
        pycom(i)=y(i)
        pthcom(i)=th(i)
        xixcom(i)=xix(i)
	xiycom(i)=xiy(i)
        xithcom(i)=xith(i)
        Dcom(i)=D(i)
      enddo
      D1com=D1
      ncom=N

      do i=1,N
         xp(i)=xpr(i)
         yp(i)=ypr(i)
         countn(i)=countnr(i)
         do j=1,countn(i)
            nl(j,i)=nlr(j,i)
         enddo
      enddo

      ax=0d0
      xx=1d0
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)
      do i=1,N
        xix(i)=xmin*xix(i)
        xiy(i)=xmin*xiy(i)
        xith(i)=xmin*xith(i)
        x(i)=x(i)+xix(i)
	y(i)=y(i)+xiy(i)
        th(i)=th(i)+xith(i)
      enddo

      return
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      double precision ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
c      EXTERNAL func
      PARAMETER (GOLD=1.618034d0, GLIMIT=100., TINY=1.d-20)
      double precision dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then ! was gt
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.gt.fc)then ! was ge
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2d0*dsign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0d0)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0d0)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0d0)then ! was ge
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION df1dim(x)
      PARAMETER (Ntot=4096)
      double precision df1dim,x,maxdis
      double precision fx(Ntot),fy(Ntot),fth(Ntot)
      double precision pxcom(Ntot),pycom(Ntot),pthcom(Ntot)
      double precision xixcom(Ntot),xiycom(Ntot),xithcom(Ntot)
      double precision xt(Ntot),yt(Ntot),tht(Ntot)
      double precision xp(Ntot),yp(Ntot)
      double precision Dcom(Ntot),D1com
      double precision alpha(Ntot),width
      integer ncom,countn(Ntot),nl(800,Ntot)
      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      common /f2com/ width
      common /f3com/ alpha ! aspect ratio

      do i=1,ncom
        xt(i)=pxcom(i)+x*xixcom(i)
        yt(i)=pycom(i)+x*xiycom(i)
        tht(i)=pthcom(i)+x*xithcom(i)
      enddo

      call CG_check(ncom,xt,yt,xp,yp,maxdis)
      if(maxdis.gt.width*D1com) then
	call makelist(ncom,xt,yt,Dcom,D1com,xp,yp,countn,nl)
      end if	
      call dfunc(ncom,xt,yt,tht,Dcom,D1com,fx,fy,fth,countn,nl)

      df1dim=0.D0
      do i=1,ncom
        df1dim=df1dim+fx(i)*xixcom(i)+fy(i)*xiycom(i)+fth(i)*xithcom(i)
      enddo

      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin)
      INTEGER ITMAX
      double precision dbrent,ax,bx,cx,tol,xmin,df,f,ZEPS
      EXTERNAL df,f
      PARAMETER (ITMAX=10000,ZEPS=1.0e-12)
      INTEGER iter
      double precision a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw
      double precision fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm
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
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.*tol1
        if(dabs(x-xm).le.(tol2-0.5d0*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          d1=2.*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0d0).and.(dx*d1.le.0d0)
          ok2=((a-u2)*(u2-b).gt.0d0).and.(dx*d2.le.0d0)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(dabs(d1).lt.dabs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(dabs(d).gt.dabs(0.5d0*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=dsign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0d0) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5d0*e
2        if(dabs(d).ge.tol1) then
          u=x+d
           fu=f(u)
        else
          u=x+dsign(tol1,d)
          fu=f(u)
          if(fu.gt.fx)goto 3
        endif
        du=df(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
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
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
c      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
      return
      END
c  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION f1dim(x)
      PARAMETER (Ntot=4096)
      double precision f1dim,x
      double precision pxcom(Ntot),pycom(Ntot),pthcom(Ntot)
      double precision xixcom(Ntot),xiycom(Ntot),xithcom(Ntot)
      double precision xt(Ntot),yt(Ntot),tht(Ntot),width
      double precision Dcom(Ntot),maxdis,D1com
      double precision xp(Ntot),yp(Ntot),alpha(Ntot)
      integer nl(800,Ntot),countn(Ntot),ncom
      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      
      common /f2com/ width
      common /f3com/ alpha ! aspect ratio

      do i=1,ncom
        xt(i)=pxcom(i)+x*xixcom(i)
        yt(i)=pycom(i)+x*xiycom(i)
        tht(i)=pthcom(i)+x*xithcom(i)
      enddo

      call CG_check(ncom,xt,yt,xp,yp,maxdis)
      if(maxdis.gt.width*D1com) then
	call makelist(ncom,xt,yt,Dcom,D1com,xp,yp,countn,nl)
      end if
      call func(ncom,xt,yt,tht,Dcom,D1com,f1dim,countn,nl)

      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ran2(idum)
      integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      integer idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .                 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! BELOW IS A NEW CODE ADDED BY PAWEL GNIEWEK !!!!!!!!
      SUBROUTINE determine_system_status(N,D,D1,alpha,ftol,wide,dt,
     +     fret, phi_j,dphi,
     +     dtstatus, terminate,
     +     before_jamming, at_jamming, above_jamming)
      IMPLICIT NONE
      integer Ntot, N
      PARAMETER(Ntot = 4096)
      
      double precision delta, phi, calc_phi
      double precision D(Ntot),alpha(Ntot), D1
      double precision fret,phitemp,phi_j,ftol,wide,dt,dphi
      integer dtstatus, terminate
      integer before_jamming, at_jamming, above_jamming
      
      delta = 0.5d-8
      
      
      phi = calc_phi(D, alpha, D1, N)
      if (fret.lt.ftol*N .or. phi.lt.0.25) then
          terminate = 0
          before_jamming = 1
          at_jamming = 0
          above_jamming = 0
          dtstatus = 0
          
          dt = 1d0
          phi_j = 0d0
          return
      endif
      
      
      if (above_jamming.eq.0 .and. at_jamming .eq. 0)then
          if (fret.gt.(ftol * wide * N) )  then
             dt = 0.5*dt
             dtstatus = 1
             terminate = 0
             return 
          else if (fret.gt.(ftol*N) .and. fret.lt.(ftol*wide*N))then
              before_jamming = 0
              at_jamming = 1
              above_jamming = 1
              dt = 1d0
              dtstatus = 0
              terminate = 0
              phi_j = phi
              return
          ! else - not needed. this point should never been reached!
              
        endif
      endif
      
      
      if (phi_j .lt. 0.25) then
          write(*,*) "ERROR: SHOULD NOT REACH THIS POINT"
      endif
      
      if (at_jamming.ne.0) then
          write(*,*) "SOMETHING WRONG: should be at_jamming = 0"
      endif

      if (above_jamming.ne.1) then
          write(*,*) "SOMETHING WRONG: should be above_jamming = 1"
      endif      
      
      at_jamming = 0
      if (phi .lt. phi_j+dphi-delta) then
          dt = 1
          terminate = 0
          dtstatus = 0
          return
      else if (phi .gt. phi_j+dphi+delta) then
          dt = 0.5 * dt
          dtstatus = 1
          terminate = 0
          return
      else
          dtstatus = 0
          terminate = 1
          return
      endif
      
      END
      
      
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DOUBLE PRECISION FUNCTION calc_phi(D, alpha, D1, N)
      IMPLICIT NONE
      integer Ntot
      double precision pi
      PARAMETER(pi=3.1415926535897932d0)
      PARAMETER(Ntot = 4096)
      double precision D(Ntot),alpha(Ntot)
      double precision Lx,Ly
      double precision D1,phit
      integer i,N
         
      common /f5com/ Lx,Ly
         
      phit=0d0
         
      do i=1,N
             phit=phit+(1d0+(alpha(i)-1d0)**2)*D(i)**2
      enddo
         
      calc_phi=pi*D1*D1*phit/Lx/Ly/4d0
         
      END FUNCTION
      
      
      subroutine copy_everything(Ntot, Ngen, ndiv, x, y, th, D, alpha, 
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
      integer Ntot,Ngen
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot)
      double precision alpha(Ntot),rate(Ntot),alphar(Ntot),scale(Ntot)
      double precision P,PP(Ntot),D0(Ntot),PPm(Ntot)
      double precision aclone(Ntot*Ngen)
      double precision tbirth(Ntot*Ngen)
      double precision xbirth(Ntot*Ngen),ybirth(Ntot*Ngen)
      double precision xdiv(999,2),ydiv(999,2),thdiv(999,2)
      integer N,c(Ntot,Ngen),Nexist
      integer nclone(Ntot*Ngen),nclonebox(Ntot*Ngen)
      integer age(Ntot),agehist(Ngen,2)
      integer idiv(999,2), ndiv
      double precision xcc(Ntot),ycc(Ntot),thcc(Ntot),Dcc(Ntot)
      double precision alphacc(Ntot),ratecc(Ntot),alpharcc(Ntot)
      double precision scalecc(Ntot)
      double precision Pcc,PPcc(Ntot),D0cc(Ntot),PPmcc(Ntot)
      double precision aclonecc(Ntot*Ngen)
      double precision tbirthcc(Ntot*Ngen)
      double precision xbirthcc(Ntot*Ngen),ybirthcc(Ntot*Ngen)
      double precision xdivcc(999,2),ydivcc(999,2),thdivcc(999,2)
      integer Ncc,ccc(Ntot,Ngen),Nexistcc
      integer nclonecc(Ntot*Ngen),ncloneboxcc(Ntot*Ngen)
      integer agecc(Ntot),agehistcc(Ngen,2)
      integer idivcc(999,2), ndivcc
      integer i,j, BUDCONT(Ntot), BUDCONTCC(Ntot)
      
!      
      Pcc = P
      Ncc = N
      Nexistcc = Nexist
      ndivcc = ndiv     
      
! ZERO EVERYTHING !!!
      do i=1,Ntot
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
      enddo
      
      do i=1,999
          xdivcc(i,1) = 0d0
          ydivcc(i,1) = 0d0
          thdivcc(i,1) = 0d0
          idivcc(i,1) = 0d0
          xdivcc(i,2) = 0d0
          ydivcc(i,2) = 0d0
          thdivcc(i,2) = 0d0
          idivcc(i,2) = 0d0
      enddo
      
      do i=1,Ngen
          agehistcc(i,1) = 0d0
          agehistcc(i,2) = 0d0
      enddo
      
      do i=1,Nexistcc+2!Ntot*Ngen
          tbirthcc(i) = 0d0
          xbirthcc(i) = 0d0
          ybirthcc(i) = 0d0
          aclonecc(i) = 0d0
          nclonecc(i) = 0d0
          ncloneboxcc(i) = 0d0
      enddo 
      
      do i=1,Nexistcc+2!Ntot
          do j=1,Nexistcc+2!Ngen
              ccc(i,j) = 0
          enddo
      enddo
     

!!!!!!!!!!!!!!!!!!!!!   MAKE COPY   
      do i=1,Ntot
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
      enddo
      
      do i=1,999
          xdivcc(i,1) = xdiv(i,1)
          ydivcc(i,1) = ydiv(i,1)
          thdivcc(i,1) = thdiv(i,1) 
          idivcc(i,1) = idiv(i,1)
          xdivcc(i,2) = xdiv(i,2)
          ydivcc(i,2) = ydiv(i,2)
          thdivcc(i,2) = thdiv(i,2) 
          idivcc(i,2) = idiv(i,2) 
      enddo
      
      do i=1,Ngen
          agehistcc(i,1) = agehist(i,1)
          agehistcc(i,2) = agehist(i,2)
      enddo
      
      do i=1,Nexistcc+2!Ntot*Ngen
          tbirthcc(i) = tbirth(i)
          xbirthcc(i) = xbirth(i)
          ybirthcc(i) = ybirth(i)
          aclonecc(i) = aclone(i)          
          nclonecc(i) = nclone(i)
          ncloneboxcc(i) = nclonebox(i)
      enddo
      
      do i=1,Nexistcc+2!Ntot
          do j=1,Nexistcc+2!Ngen
              ccc(i,j) = c(i,j)
          enddo
      enddo
      end
      
      subroutine copy_back_everything(Ntot, Ngen, ndiv, x, y, th, D, 
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
      integer Ntot,Ngen
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot)
      double precision alpha(Ntot),rate(Ntot),alphar(Ntot),scale(Ntot)
      double precision P,PP(Ntot),D0(Ntot),PPm(Ntot)
      double precision aclone(Ntot*Ngen)
      double precision tbirth(Ntot*Ngen)
      double precision xbirth(Ntot*Ngen),ybirth(Ntot*Ngen)
      double precision xdiv(999,2),ydiv(999,2),thdiv(999,2)
      integer N,c(Ntot,Ngen),Nexist
      integer nclone(Ntot*Ngen),nclonebox(Ntot*Ngen)
      integer age(Ntot),agehist(Ngen,2)
      integer idiv(999,2), ndiv
      double precision xcc(Ntot),ycc(Ntot),thcc(Ntot),Dcc(Ntot)
      double precision alphacc(Ntot),ratecc(Ntot),alpharcc(Ntot)
      double precision scalecc(Ntot)
      double precision Pcc,PPcc(Ntot),D0cc(Ntot),PPmcc(Ntot)
      double precision aclonecc(Ntot*Ngen)
      double precision tbirthcc(Ntot*Ngen)
      double precision xbirthcc(Ntot*Ngen),ybirthcc(Ntot*Ngen)
      double precision xdivcc(999,2),ydivcc(999,2),thdivcc(999,2)
      double precision dd
      integer Ncc,ccc(Ntot,Ngen),Nexistcc
      integer nclonecc(Ntot*Ngen),ncloneboxcc(Ntot*Ngen)
      integer agecc(Ntot),agehistcc(Ngen,2)
      integer idivcc(999,2), ndivcc
      integer i,j, celltype, BUDCONT(Ntot), BUDCONTCC(Ntot)
      common /f10com/ celltype
            
!
      P = Pcc
      N = Ncc
      Nexist = Nexistcc
      ndiv = ndivcc
!
      do i=1,Ntot
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
      enddo
      
      do i=1,999
          xdiv(i,1) = xdivcc(i,1)
          ydiv(i,1) = ydivcc(i,1)
          thdiv(i,1) = thdivcc(i,1) 
          idiv(i,1) = idivcc(i,1)
          xdiv(i,2) = xdivcc(i,2)
          ydiv(i,2) = ydivcc(i,2)
          thdiv(i,2) = thdivcc(i,2) 
          idiv(i,2) = idivcc(i,2) 
      enddo
      
      do i=1,Ngen
          agehist(i,1) = agehistcc(i,1)
          agehist(i,2) = agehistcc(i,2)
      enddo
      
      do i=1,Nexist+2!Ntot*Ngen
          tbirth(i) = tbirthcc(i)
          xbirth(i) = xbirthcc(i)
          ybirth(i) = ybirthcc(i)
          aclone(i) = aclonecc(i)          
          nclone(i) = nclonecc(i)
          nclonebox(i) = ncloneboxcc(i)
      enddo
      
      do i=1,Nexist+2!Ntot
          do j=1,Nexist+2!Ngen
              c(i,j) = ccc(i,j)
          enddo
      enddo

         do i=1,N
            if(celltype.eq.1) then
               scale(i)=dsqrt(1d0+alpha(i)**2)/4d0*d(i)         
            elseif(celltype.eq.2) then
               dd=alpha(i)-1d0
               scale(i)=dsqrt(2d0*(1d0+dd**4)/(1+dd**2)+
     +              4d0*(dd*(1d0+dd)/(1+dd**2))**2)/4d0*d(i)
            elseif(celltype.eq.3) then
               scale(i)=dsqrt(2d0)/4d0*d(i)
            endif
            th(i)=th(i)*scale(i)
         enddo
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! remove floaters !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE remove_floaters(x,y,th,phi,D1,D,N,F)
      IMPLICIT NONE
      integer Ntot, N 
      PARAMETER(Ntot = 4096)
      double precision pi
      PARAMETER(pi=3.141592653589793238d0)
      integer R, Rnew, i, j, k, II, JJ, c(Ntot), listP(Ntot), DoF
      double precision overlap, aspect_ratio
      integer l, F(Ntot)
      double precision x(Ntot),y(Ntot),th(Ntot),phi,alpha(Ntot)
      double precision thIJ, t(3), sign(3), AR2, cont,scale(Ntot)
      double precision xij, yij, rij, D(Ntot), D1, z_ave
      double precision dij,dtij,dthi,dthj,dtij2,dthi2,dthj2
      double precision dtijthi,dtijthj,dthithj, compression
      double precision x2(Ntot), y2(Ntot), th2(Ntot), Lx,Ly
      integer count, Z, nn
      
      common /f3com/ alpha ! aspect ratio
      common /f5com/ Lx,Ly
      common /f9com/ scale

      ! Remove floaters
      k=0
      phi=0.d0
      do i=1,N
         ii=i-k
         nn=n-k
         if(F(ii).eq.1) then
            k=k+1
            do j=ii,nn
               F(j) = F(j+1)
               x(j) = x(j+1)
               y(j) = y(j+1)
               th(j) = th(j+1)
               D(j) = D(j+1)
               alpha(j) = alpha(j+1)
               scale(j) = scale(j+1)              
            enddo
         endif
      enddo

      N=N-k
      
      do i=1,N
         phi=phi+D(i)*D(i)
      enddo
      phi=pi*D1*D1*phi/Lx/Ly/4d0
      
      END

      
      
!!!!!!!!!!!!!!!! counts contacts !!!!!!!!!!!!!
      SUBROUTINE contacts(x,y,th,D1,D,N,Z,F,Nf,Nu,Nmm,Nbb,Nmb)
      IMPLICIT NONE
      integer Ntot, N        
      PARAMETER(Ntot = 4096)
      integer R, Rnew, i, j, k, II, JJ, c(Ntot), listP(Ntot)
      double precision overlap, aspect_ratio
      integer l, F(Ntot), Nf, Nu,Nmm,Nbb,Nmb
      double precision x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      double precision thIJ, t(3), sign(3),cont
      double precision xij, yij, rij, D(Ntot), D1, z_ave
      double precision dij,dtij,dthi,dthj,dtij2,dthi2,dthj2
      double precision dtijthi,dtijthj,dthithj
      double precision x2(Ntot), y2(Ntot), th2(Ntot)
      integer Z
      common /f3com/ alpha ! aspect ratio
      
      Z = 0
      R = 0

      Rnew = 1
      Nf=0
      Nu=0
      Nmm = 0 
      Nbb = 0 
      Nmb = 0
      
      do i=1, N
         c(i)=0
         listP(i)=i
         F(i) = 0
      enddo
      
      do while(Rnew>0)
         Rnew = 0
         do i=1, N-R
            II=listP(i)
            c(II)=0
         enddo
         
         ! Find Contacts
         do i=1,N-R
            II=listP(i)
            do j=1, i-1
               JJ=listP(j)	    
!               if(overlap(N,x,y,th,D,D1,II,JJ).lt.1.d0) then
                  c(II) = c(II)+overlap(N,x,y,th,D,D1,II,JJ) !1
                  c(JJ) = c(JJ)+overlap(N,x,y,th,D,D1,II,JJ) !1
!               endif
            enddo
         enddo

         ! Remove floaters
         i=1
         do while (i.le.N-R)
            II=listP(i)
            if(c(II).lt.3) then
               Nf = Nf + 1
               Rnew=Rnew+1
               R=R+1
               do j=i, N-R
                  listP(j) = listP(j+1)
               enddo
               F(II) = 1
            else
               i=i+1
            endif
         enddo
      enddo
      
      
      !Count Contacts
      Z=0
      do i=1, N-R
         II=listP(i)
         do j=1, i-1
            JJ=listP(j)	    
!            if(overlap(N,x,y,th,D,D1,II,JJ).lt.1.d0) then
               Z = Z+2*overlap(N,x,y,th,D,D1,II,JJ)
!            endif
         enddo
      enddo

      z_ave = dble(Z)/dble(N-R)

      END
      
      
!     CUSTOM OVERLAP FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION overlap(N,x,y,th,D,D1,i,j)
      IMPLICIT NONE
      integer Ntot, N     
      PARAMETER(Ntot = 4096)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1
      double precision rij,xij,yij,dij,ep,overlap
      double precision dtij,dthi,dthj,dtij2,dthi2,dthj2
      double precision dtijthi,dtijthj,dthithj, width
      double precision alpha(Ntot)
      double precision dd,dr1,dr2,dk2,exp,att
      double precision Lx,Ly
      !
      double precision dijsq_up,scale(Ntot),c(Ntot)
      double precision s(Ntot),dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      double precision dk(Ntot,2),di_up(Ntot),di1j1      
      integer celltype, ni, i, j
 
      
      common /f2com/ width
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly
      common /f9com/ scale
      common /f10com/ celltype
      
      if (D1.ne.1.d0) then
          write(*,*) "SOMETHING IS WRONG. D1 SHALL BE EQ. TO 1.0"
      endif
      
      overlap = 0.d0

      ! BUDDING CELLS
      if (celltype.eq.2) then
          write(*,*) "THERE IS SPECIAL CONTACT FUNCTION FOR YEAST"

      ! DISK CELLS      
      elseif (celltype.eq.3) then
                
          xij=x(i)-x(j)
          xij=xij-idnint(xij/Lx)*Lx  !! PBC
          yij=y(i)-y(j)
          yij=yij-idnint(yij/Ly)*Ly  !! PBC
          
          rij=dsqrt(xij**2+yij**2)  
          dij=D1*(d(i)+d(j))/2d0

          if(rij.lt.dij) then
              overlap = 1.d0
          else
              overlap = 0.d0
          end if
      endif
          
      END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! inuse
      SUBROUTINE contacts_yeast(x,y,th,D1,D,N,Z,F,Nf,Nu,Nmm,Nbb,Nmb)
      IMPLICIT NONE
      integer Ntot, N 
      PARAMETER(Ntot = 4096)
      integer i, j, F(Ntot), Z, NCBUD(Ntot,2), Nf, Nu, Nmm, Nbb, Nmb
      integer ki,kj,k
      double precision overlap, aspect_ratio
      double precision x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      double precision xij, yij, D(Ntot), D1
      double precision exp,dij_up,dij
      double precision Lx,Ly,rijsq,c(Ntot),att
      double precision s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      double precision dk(Ntot,2)
      integer flag
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly

      Z = 0
      Nf= 0
      Nu= 0
      
      Nmm = 0 
      Nbb = 0 
      Nmb = 0
      
      flag = 0
      
      ! convert to from molecules to atoms
      do i=1,N
         c(i)=dcos( th(i) )
         s(i)=dsin( th(i) )         
         dd=alpha(i)-1d0
         dr(i,1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)*D1/2d0
         dr(i,2)=-(1d0+dd)/(1d0+dd**2)*D(i)*D1/2d0
         do k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         enddo
         dk(i,1)=D(i)*D1
         dk(i,2)=dd*D(i)*D1
      enddo
      
      do i=1,N
         NCBUD(i,1)=0
         NCBUD(i,2)=0
      enddo

      do i=1,N-1
          do j=i+1, N
              do ki=1,2
                  do kj=1,2
                      dij=(dk(i,ki)+dk(j,kj))/2d0
                      xij=xa(i,ki)-xa(j,kj)
                      xij=xij-idnint(xij/Lx)*Lx  !! PBC
                      yij=ya(i,ki)-ya(j,kj)
                      yij=yij-idnint(yij/Ly)*Ly !! PBC
                      rijsq=xij**2+yij**2
                      if(rijsq.lt.(dij**2)) then
                          Z = Z+2
                          NCBUD(i,ki)=NCBUD(i,ki)+1
                          NCBUD(j,kj)=NCBUD(j,kj)+1
                          if (ki.eq.1 .and. kj.eq.1) then
                              Nmm = Nmm + 2
                          elseif (ki.eq.2 .and. kj.eq.2) then
                              Nbb = Nbb + 2
                          else
                              Nmb = Nmb + 2
                          endif
                      endif
                  enddo
              enddo
          enddo
      enddo

      
      do i=1,N
          flag = 0
          if ( (NCBUD(i,1)+NCBUD(i,2)).lt.3 ) then
              nf = nf + 1
              flag = 1
          endif
          
          if ( (NCBUD(i,1)+NCBUD(i,2)).eq.3 ) then
              if ( NCBUD(i,1).eq.2 .or. NCBUD(i,2).eq.2 ) then
                  nf = nf + 1
                  flag = 1
              endif
          endif
          
          if (flag.eq.0) then    
              if ( NCBUD(i,1).eq.0 ) then
                  nu = nu + 1
              endif
              if ( NCBUD(i,2).eq.0 ) then
                  nu = nu + 1
              endif
          endif
          
      enddo
      
      end
      
      subroutine out_numbers(N, Nf, Nu, Ziso)
      implicit none
      integer Ntot
      PARAMETER(Ntot = 4096)
      integer N, Nf, Nu, Ziso
      
      Ziso = 6*(N-Nf) - 2*Nu - 2
      
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      SUBROUTINE bud_contacts(x,y,th,D1,D,N,BUDCONT,ZEROBUDS)
      IMPLICIT NONE
      integer Ntot, N
      PARAMETER(Ntot = 4096)
      integer i, j, NCBUD(Ntot,2), ZEROBUDS, BUDCONT(Ntot),k,ki,kj
      double precision x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      double precision xij, yij, D(Ntot), D1
      double precision exp,dij
      double precision Lx,Ly,rijsq,c(Ntot),att 
      double precision s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      double precision dk(Ntot,2)

      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly

      ZEROBUDS = 0
      
      ! convert to from molecules to atoms
      do i=1,N
         c(i)=dcos( th(i) )
         s(i)=dsin( th(i) )         
         dd=alpha(i)-1d0
         dr(i,1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)*D1/2d0
         dr(i,2)=-(1d0+dd)/(1d0+dd**2)*D(i)*D1/2d0
         do k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         enddo
         dk(i,1)=D(i)*D1
         dk(i,2)=dd*D(i)*D1
      enddo
      
      do i=1,N
         NCBUD(i,1)=0
         NCBUD(i,2)=0
         BUDCONT(i) = 0
      enddo

      do i=1,N-1
          do j=i+1, N
               do ki=1,2
                  do kj=1,2
                      dij=(dk(i,ki)+dk(j,kj))/2d0
                      xij=xa(i,ki)-xa(j,kj)
                      xij=xij-idnint(xij/Lx)*Lx !! PBC
                      yij=ya(i,ki)-ya(j,kj)
                      yij=yij-idnint(yij/Ly)*Ly !! PBC
                      rijsq=xij**2+yij**2
                      if(rijsq.lt.(dij**2)) then
                          NCBUD(i,ki)= NCBUD(i,ki)+1
                          NCBUD(j,kj)= NCBUD(j,kj)+1
                      endif
                  enddo
              enddo
          enddo
      enddo
      
      do i=1,N
          BUDCONT(i) = NCBUD(i,2)
          if (BUDCONT(i).eq.0)then
              ZEROBUDS = ZEROBUDS + 1
          endif
      enddo

     
      END
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      SUBROUTINE growth_rate(gr,N,rate,PP,P0,D,alpha)    
      IMPLICIT NONE
      integer Ntot,i,N
      PARAMETER(Ntot = 4096)
      double precision ratei,P0
      double precision rate(Ntot),PP(Ntot),D(Ntot)
      double precision alpha(Ntot)
      
      double precision Dm, Db, dd, gr
      double precision total_area, total_newmass
      
      gr = 0.0d0
      total_area = 0.0d0
      total_newmass = 0.0d0
      
      do i=1,N
          ratei=rate(i)
          if(P0.gt.0d0.and.PP(i).gt.0d0) then
              ratei=ratei*dexp(-PP(i)/P0)
          endif
          
          Dm = D(i)
          Db = alpha(i)-1d0
          total_area = total_area + Dm*Dm + Db*Db
          total_newmass = total_newmass + ratei*(Dm*Dm + Db*Db)
      enddo
      
      gr = total_newmass / total_area
      
      END
