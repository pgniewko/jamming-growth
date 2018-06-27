      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!
      !!
      !!      Author: Pawel Gniewek
      !!      Email:
      !!      License:
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PROGRAM shear_packing

      IMPLICIT NONE
      integer Ntot
      parameter(Ntot=4096)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,exp
      double precision ftol,ftol1,fret,width,Lx,Ly
      double precision alpha(Ntot),scale(Ntot)
      double precision phi,flow,P,PP(Ntot),D0(Ntot)
      double precision xa(2),ya(2),PR,PT
      double precision cc,ss,dr(2),dd,att
      integer N,Nr,iter,i,kk
      integer Nf,Nu,Nmm,Nbb,Nmb
      
      character file_conf*150
      character file_traj*150
      character file_G*150
      
      integer F(Ntot), Nc, Ziso
      double precision calc_phi, wide
      double precision phi_init,ddelrx,delrx
      
      integer shear_steps, step_index
      double precision SSTRESS, cory, SSTRESS0

! END OF NEW DATA      
      
      common /f2com/ width
      common /f3com/ alpha
      common /f4com/ exp,att
      common /f5com/ Lx,Ly
      common /f6com/ P,PP,PT,PR
      common /f9com/ scale
      common /f11com/ delrx
      common /f12com/ SSTRESS

      ! read geometric parameters
      read(*,*) Lx
      read(*,*) Ly

      ! read cell parameters
      read(*,*) att
      
      ! read output files
      read(*,*) file_conf
      read(*,*) ddelrx
      read(*,*) shear_steps
      
      ! parameters
      D1=1d0         ! Minor axis of particle; D1=1.0 - circle
      exp=2d0        ! 2 =  LS, 2.5 = Hertzian, >2.9 = RLJ
      width=0.1d0    ! width of neighborlist 
      ftol=0.5d-16   ! Condition 1 for frprmn: V/N < ftol 
      ftol1=0.5d-16  ! Condition 2 for frprmn: dV/N < ftol1
      
      wide = 2.0d0

      SSTRESS = 0d0
      SSTRESS0= 0d0

      ! FILES
 123  open(unit=1,file=TRIM(file_conf)) ! CONFIGURATION FILE
 
      file_traj='prod_shear_' // TRIM(file_conf)
      open(unit=11,file=TRIM(file_traj), status='replace')
      
      file_G ='G_data_' // TRIM(file_conf)
      open(unit=12,file=TRIM(file_G), status='replace')

!     READ PACKING FROM FILE      
      read(1, *) N, phi_init
      do i=1,N
          read(1, *) x(i),y(i),D(i),alpha(i),th(i)
      enddo
      
      do step_index=1,shear_steps
          
         delrx = dble(step_index-1)*ddelrx

         ! convert from angle to length scale = sqrt(I/m) * angle
         do i=1,N
            dd=alpha(i)-1d0
            scale(i)=dsqrt(2d0*(1d0+dd**4)/(1+dd**2)+
     +              4d0*(dd*(1d0+dd)/(1+dd**2))**2)/4d0*d(i)
     
            th(i)=th(i)*scale(i)
         enddo
         
         ! minimize energy
         call frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)

         if (delrx.eq.0d0)then
             SSTRESS0 = SSTRESS
         endif
         
         write(*,*) N,fret/dble(N), 
     +     P,delrx,ddelrx,step_index,SSTRESS,SSTRESS-SSTRESS0
    

         flush(12)
         ! convert back to angles
         do i=1,N
            th(i)=th(i)/scale(i)
         enddo
         
         call contacts_yeast(x,y,th,D1,D,N,Nc,F,Nf,Nu,Nmm,Nbb,Nmb)
         call out_numbers(N, Nf, Nu, Ziso)
         write(12,'(3E26.18,5I12)')delrx,SSTRESS,SSTRESS-SSTRESS0,
     +         N,Nc,Nf,Nu,Ziso   
       
         phi = calc_phi(D, alpha, D1, N) 
         write(11,*) 2*N, phi
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
              write(11,'(3E26.18,I12)')xa(1),ya(1),d(i),0
              write(11,'(3E26.18,I12)')xa(2),ya(2),d(i)*dd,1
         enddo         
         
      enddo
      
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CG_check(N,x,y,xp,yp,maxdis)
      IMPLICIT NONE
      integer Ntot,N,i
      parameter(Ntot = 4096)
      double precision maxdis,x(Ntot),y(Ntot),xp(Ntot),yp(Ntot)

      maxdis=0d0
      do i=1,N
         maxdis=max(dabs(x(i)-xp(i)),maxdis)
         maxdis=max(dabs(y(i)-yp(i)),maxdis)
      enddo
      maxdis=2d0*dsqrt(2d0*maxdis*maxdis)

      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE makelist(N,x,y,D,D1,xp,yp,countn,nl)
      IMPLICIT NONE
      integer Ntot,N
      parameter(Ntot = 4096)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot),D1
      integer countn(Ntot),nl(800,Ntot)
      
      call makelist_dimer(N,x,y,D,D1,xp,yp,countn,nl)

      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE makelist_dimer(N,x,y,D,D1,xp,yp,countn,nl)
      IMPLICIT NONE
      integer Ntot
      parameter(Ntot = 4096)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot)
      double precision D1,xij,yij,rij,dij,rijsq,alpha(Ntot),width
      double precision dd,dr1,dr2,dk2,di_up(Ntot),exp,att
      double precision Lx,Ly,delrx,cory
      integer countn(Ntot),nl(800,Ntot),N,i,j
      common /f2com/ width
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly
      common /f11com/ delrx

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
            yij=y(i)-y(j)
            cory=idnint(yij/Ly)
            xij=xij-cory*delrx*Lx
            xij=xij-idnint(xij/Lx)*Lx  !! PBC
            yij=yij-cory*Ly  !! PBC
            
            rijsq=xij*xij+yij*yij
            dij=(di_up(i)+di_up(j))/2d0
            dij=D1*( D(i) + D(j) )/2d0   
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
      
      END
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE func(N,x,y,th,D,D1,V,countn,nl)
      IMPLICIT NONE
      integer Ntot,N    
      parameter(Ntot = 4096)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,V
      integer countn(Ntot),nl(800,Ntot)
      
      call func_dimer(N,x,y,th,D,D1,V,countn,nl)
      
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE dfunc(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      parameter(Ntot = 4096)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1
      double precision fx(Ntot),fy(Ntot),fth(Ntot)
      integer countn(Ntot),nl(800,Ntot),N

      call dfunc_dimer(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE func_dimer(N,x,y,th,D,D1,V,countn,nl)
      parameter(Ntot = 4096)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,V,alpha(Ntot)
      double precision rij,xij,yij,dij,exp,dlnsig,dij_up,sigma,LJ
      double precision Lx,Ly,rijsq,dijsq_up,scale(Ntot),c(Ntot),att
      double precision s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2),Vij
      double precision dk(Ntot,2),di_up(Ntot),di1j1,xhit,yhit,yhitout
      double precision delrx, cory
      
      integer countn(Ntot),nl(800,Ntot),N
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly
      common /f9com/ scale
      common /f11com/ delrx

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
               yij=y(i)-y(j)
               cory=idnint(yij/Ly)
               xij=xij-cory*delrx*Lx
               xij=xij-idnint(xij/Lx)*Lx  !! PBC
               yij=yij-cory*Ly  !! PBC

               if(dabs(xij).lt.dij_up+att) then
!                  yij=y(i)-y(j)
!                  yij=yij-cory*Ly  !! PBC
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij_up+att)**2) then
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     do ki=1,2
                        do kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           yij=ya(i,ki)-ya(j,kj)
                           cory = idnint(yij/Ly)
                           xij=xij-cory*delrx*Lx
                           xij=xij-idnint(xij/Lx)*Lx  !! PBC
                           yij=yij-cory*Ly  !! PBC
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
			
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE dfunc_dimer(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      parameter(Ntot = 4096)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),sigma,D(Ntot),D1,dij
      double precision fx(Ntot),fy(Ntot),fth(Ntot),rij,xij,yij,fr,exp
      double precision dij_up,alpha(Ntot),LJ,fc,ft,f_x,f_y,scale(Ntot)
      double precision fthi,fthj,fth_c,Lx,Ly,P,Pij,rijsq,dijsq_up
      double precision s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2),att
      double precision dk(Ntot,2),di_up(Ntot),di1j1,xhit,yhit,yhitout
      double precision PP(Ntot),c(Ntot),Vij,PT,PR, delrx, cory
      double precision SSTRESS
      integer countn(Ntot),nl(800,Ntot),N !,growth_flag
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly
      common /f6com/ P,PP,PT,PR
      common /f9com/ scale
      common /f11com/ delrx
      common /f12com/ SSTRESS

      do i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
         PP(i)=0d0
      enddo
      P=0d0
      SSTRESS=0d0

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
               yij=y(i)-y(j)
               cory=idnint(yij/Ly)
               xij=xij-cory*delrx*Lx
               xij=xij-idnint(xij/Lx)*Lx  !! PBC
               yij=yij-cory*Ly  !! PBC
               
               if(dabs(xij).lt.dij_up+att) then
!                  yij=y(i)-y(j)
!                  yij=yij-idnint(yij/Ly)*Ly !! PBC
!                  yij=yij-idnint(xij/Lx)*Lx*delrx
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij_up+att)**2) then
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     do ki=1,2
                        do kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           yij=ya(i,ki)-ya(j,kj)
                           cory = idnint(yij/Ly)
                           xij=xij-cory*delrx*Lx
                           xij=xij-idnint(xij/Lx)*Lx  !! PBC
                           yij=yij-cory*Ly  !! PBC
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
                              
                              SSTRESS=SSTRESS-xij*yij/rij*fc
                              
                              if (ki.eq.2) then
                                  PP(i)=PP(i)+Pij
                              endif
                              if (kj.eq.2) then
                                  PP(j)=PP(j)+Pij
                              endif
                              
                           endif
                        enddo
                     enddo
                  end if
               end if
            enddo
         end if
      enddo

      !write(*,*)  SSTRESS
      
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
      do i=1,N
         PP(i)=PP(i)*dble(N)/4d0/Lx/Ly
      enddo
      
      SSTRESS=SSTRESS/Lx/Ly
							
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DOUBLE PRECISION FUNCTION calc_phi(D, alpha, D1, N)
      IMPLICIT NONE
      integer Ntot
      parameter(Ntot = 4096)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision D(Ntot),alpha(Ntot)
      double precision Lx,Ly,phis,D1,phit
      integer N,i
      common /f5com/ Lx,Ly
         
      phit=0d0
      do i=1,N
         phit=phit+(1d0+(alpha(i)-1d0)**2)*D(i)**2
      enddo
         
      calc_phi=pi*D1*D1*phit/Lx/Ly/4d0
      END FUNCTION
     
      
      ! INUSE
      SUBROUTINE contacts_yeast(x,y,th,D1,D,N,Z,F,Nf,Nu,Nmm,Nbb,Nmb)
      IMPLICIT NONE
      integer Ntot, N 
      parameter(Ntot = 4096)
      integer i, j, F(Ntot), Z, NCBUD(Ntot,2), Nf, Nu, Nmm, Nbb, Nmb
      integer ki,kj,k
      double precision overlap, aspect_ratio
      double precision x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      double precision xij, yij, D(Ntot), D1
      double precision exp,dij_up,dij
      double precision Lx,Ly,rijsq,c(Ntot),att ! ,scale(Ntot)
      double precision s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2)
      double precision dk(Ntot,2)
      double precision delrx,cory
      integer flag
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly
      common /f11com/ delrx
!      common /f9com/ scale

      Z = 0
      Nf= 0
      Nu= 0
      
      Nmm = 0 
      Nbb = 0 
      Nmb = 0
      
      flag = 0
      
      ! convert to from molecules to atoms
      do i=1,N
!         c(i)=dcos(th(i)/scale(i))
!         s(i)=dsin(th(i)/scale(i))
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
                      yij=ya(i,ki)-ya(j,kj)
                      cory=idnint(yij/Ly)
                      xij=xij-cory*delrx*Lx
                      xij=xij-idnint(xij/Lx)*Lx  !! PBC
                      yij=yij-cory*Ly  !! PBC
                      
!                      xij=xa(i,ki)-xa(j,kj)
!                      xij=xij-idnint(xij/Lx)*Lx  !! PBC
!                      yij=ya(i,ki)-ya(j,kj)
!                      yij=yij-idnint(yij/Ly)*Ly !! PBC
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
      
      END
      
      
      SUBROUTINE out_numbers(N, Nf, Nu, Ziso)
      IMPLICIT NONE
      integer Ntot
      parameter(Ntot = 4096)
      integer N, Nf, Nu, Ziso
      
      Ziso = 6*(N-Nf) - 2*Nu - 2
      
      END      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!   BELOW THERE IS STANDARD NUMERICAL CODE   !!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      SUBROUTINE frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
      parameter(Ntot = 4096)
      integer its,iter,ITMAX
      double precision fret,ftol,EPS,ftol1
      parameter (EPS=1d-10,ITMAX=1000000000)
      double precision dgg,fp,gam,gg,gx(Ntot),gy(Ntot),hx(Ntot),hy(Ntot)
      double precision D(Ntot),D1,xix(Ntot),xiy(Ntot),xith(Ntot),width
      double precision x(Ntot),y(Ntot),maxdis,xp(Ntot),yp(Ntot)
      double precision th(Ntot),hth(Ntot),gth(Ntot),exp,att,V
      integer N,countn(Ntot),nl(800,Ntot)

      ! not needed
      double precision f1,f2,f3,fxe,fye,fthe
      double precision xi,yi,thi,max1,max2,max3,del
      double precision alpha(Ntot),Lx,Ly

      ! not needed
      common /f3com/ alpha ! aspect ratio
      common /f5com/ Lx,Ly

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
            if(dabs(fret-fp).lt.ftol1*fp.or.fret.lt.ftol*dble(N)) then
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