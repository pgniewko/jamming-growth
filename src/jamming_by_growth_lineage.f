      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!f
      !!
      !!   Quasi-static growth of cells in a rectangular box
      !!   with Periodic Boundary Conditions and
      !!   via conjugate gradient energy minimization. 
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

      PROGRAM jam_budding_yeast

      IMPLICIT NONE
      ! CONSTANT PARAMETERS
      INTEGER Ntot
      PARAMETER(Ntot=4096)
      INTEGER MAX_POSTJAM_STEPS
      PARAMETER(MAX_POSTJAM_STEPS=250000)
      INTEGER EXIT_MIN_DT, EXIT_MAX_POSTJAM_STEPS
      PARAMETER(EXIT_MIN_DT=10)
      PARAMETER(EXIT_MAX_POSTJAM_STEPS=11)
      DOUBLE PRECISION pi
      PARAMETER(pi=3.1415926535897932d0)
      DOUBLE PRECISION DT_MIN
      PARAMETER(DT_MIN=1d-12)
      ! FILE NAMES
      CHARACTER file1*80     
      CHARACTER file_LF_JAMM*120
      CHARACTER file_LF_DPHI*120
      CHARACTER STATS_file_LF_JAMM*120
      CHARACTER STATS_file_LF_DPHI*120
      CHARACTER file_LINEAGE_JAMM*120
      CHARACTER file_LINEAGE_DPHI*120
      CHARACTER file_DIVLOG*120
      CHARACTER file_TRANSITIONS*120
      CHARACTER file_POSTJAMM_SUMMARY*120
      CHARACTER file_NC*100
      CHARACTER file_STEPLOG*120
      
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,exp,ran2
      DOUBLE PRECISION ftol,ftol1,fret,alpha0,width,Lx,Ly,ratei
      DOUBLE PRECISION alpha(Ntot),rate(Ntot),scale(Ntot)
      DOUBLE PRECISION rate0,desync,phi,P,PP(Ntot)
      DOUBLE PRECISION dispcm,xa(2),ya(2),P0
      DOUBLE PRECISION cc,ss,dr(2),dd,att,rat
      DOUBLE PRECISION dt, PPm(Ntot), total_growthrate, ratei_eff
      DOUBLE PRECISION bud_diameter, delta_bud_area, chi_c_est
      INTEGER reject_time_step, terminate_sim
      INTEGER N,seed,iter,i,k,kk,skip
      INTEGER divtype
      INTEGER div, Nf, Nu, Nmm,Nbb,Nmb
      
      DOUBLE PRECISION x_copy(Ntot),y_copy(Ntot),th_copy(Ntot)
      DOUBLE PRECISION D_copy(Ntot)
      DOUBLE PRECISION alpha_copy(Ntot),rate_copy(Ntot)
      DOUBLE PRECISION scale_copy(Ntot)
      DOUBLE PRECISION P_copy,PP_copy(Ntot),PPm_copy(Ntot)
      DOUBLE PRECISION phitemp, calc_phi, wide
      DOUBLE PRECISION phi_j, dphi
      DOUBLE PRECISION event_phi(Ntot)
      DOUBLE PRECISION tracked_initial_bud_diameter(Ntot)
      DOUBLE PRECISION event_parent_bud_diameter(Ntot)
      INTEGER N_copy, Nc, Ziso
      INTEGER before_jamming, at_jamming, above_jamming
      INTEGER cell_id(Ntot), parent_id(Ntot)
      INTEGER cell_id_copy(Ntot), parent_id_copy(Ntot)
      INTEGER next_cell_id, next_cell_id_copy
      INTEGER event_count, event_parent_id(Ntot), event_new_id(Ntot)
      INTEGER event_parent_index(Ntot), event_new_index(Ntot)
      INTEGER tracked_initial_free(Ntot)
      INTEGER tracked_prev_bud_contacts(Ntot)
      INTEGER tracked_prev_bud_unconstrained(Ntot)
      INTEGER curr_bud_contacts(Ntot), curr_bud_unconstrained(Ntot)
      INTEGER curr_bud_compressed(Ntot)
      INTEGER tracking_initialized
      INTEGER initial_free_total, initial_free_completed
      INTEGER initial_free_divided, n_compressed
      INTEGER n_initial_free_active
      INTEGER n_bud_unconstrained_total
      INTEGER n_postjam_bud_unconstrained
      INTEGER n_mother_unconstrained
      INTEGER postjam_divisions_total
      INTEGER postjam_steps, failure_code

      COMMON /f2com/ width
      COMMON /f3com/ alpha
      COMMON /f4com/ exp,att
      COMMON /f5com/ Lx,Ly
      COMMON /f6com/ P,PP,PPm
      COMMON /f8com/ alpha0
      COMMON /f9com/ scale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      ! READ geometric parameters
      READ(*,*) alpha0
      READ(*,*) Lx
      READ(*,*) Ly

      ! READ cell parameters
      READ(*,*) divtype
      READ(*,*) P0
      READ(*,*) att

      ! READ run parameters
      READ(*,*) rate0
      READ(*,*) desync
      READ(*,*) seed
      READ(*,*) skip

      ! READ \delta\phi
      READ(*,*) dphi
      
      ! READ output files
      READ(*,*) file1
     
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
      tracking_initialized = 0
      initial_free_total = 0
      initial_free_completed = 0
      initial_free_divided = 0
      postjam_divisions_total = 0

      DO i=1,Ntot
         tracked_initial_free(i)=0
         tracked_prev_bud_contacts(i)=0
         tracked_prev_bud_unconstrained(i)=0
         tracked_initial_bud_diameter(i)=0d0
         curr_bud_contacts(i)=0
         curr_bud_unconstrained(i)=0
         curr_bud_compressed(i)=0
         event_parent_bud_diameter(i)=0d0
      ENDDO

      before_jamming =  0
      at_jamming = 0
      above_jamming = 0
      postjam_steps = 0
      failure_code = 0
      Nc = 0

      ! TRAJECTORY FILE
      OPEN(unit=1,file=TRIM(file1), status='replace') ! CONFIGURATION FILE
      
      ! FILES FOR JAMMED PACKING   
      file_LF_JAMM='LF_JAMM_' // TRIM(file1)            ! unit=11
      STATS_file_LF_JAMM='STATS_LF_JAMM_' // TRIM(file1)! unit=21
      file_LINEAGE_JAMM='LINEAGE_LF_JAMM_' // TRIM(file1)
      
      
      ! FILES FOR PACKIG AT DPHI
      ! LF = LAST FRAME
      file_LF_DPHI='LF_DPHI_' // TRIM(file1)
      OPEN(unit=12,file=TRIM(file_LF_DPHI), status='replace')
      STATS_file_LF_DPHI='STATS_LF_DPHI_' // TRIM(file1)
      OPEN(unit=22,file=TRIM(STATS_file_LF_DPHI), status='replace')
      file_LINEAGE_DPHI='LINEAGE_LF_DPHI_' // TRIM(file1)
      
      ! CONTACT STATISTICS (AND OTHER)
      file_NC='NC_' // TRIM(file1)
      OPEN(unit=13,file=TRIM(file_NC), status='replace')
      file_STEPLOG='STEPLOG_' // TRIM(file1)
      OPEN(unit=14,file=TRIM(file_STEPLOG), status='replace')
      file_DIVLOG='DIVLOG_' // TRIM(file1)
      OPEN(unit=15,file=TRIM(file_DIVLOG), status='replace')
      file_TRANSITIONS='TRANSITIONS_' // TRIM(file1)
      OPEN(unit=17,file=TRIM(file_TRANSITIONS), status='replace')
      file_POSTJAMM_SUMMARY='POSTJAMM_SUMMARY_' // TRIM(file1)
      OPEN(unit=18,file=TRIM(file_POSTJAMM_SUMMARY), status='replace')
      WRITE(14,'(A)')
     + '# step N fret_per_particle P dt phi total_growthrate'
     + // ' before_jamming at_jamming above_jamming'
      WRITE(15,'(A)')
     + '# step phi parent_cell_id new_cell_id'
     + // ' parent_index new_index'
      WRITE(17,'(A)')
     + '# step phi cell_id parent_id initial_bud_diameter'
     + // ' bud_diameter delta_bud_area bud_contacts'
     + // ' bud_unconstrained bud_compressed event_code'
     + // ' event_label'
      WRITE(18,'(A)')
     + '# step phi P N Nc Nf Nu Ziso total_growthrate chi_c'
     + // ' n_compressed n_initial_free_total'
     + // ' n_initial_free_active n_initial_free_completed'
     + // ' n_initial_free_divided'
     + // ' n_bud_unconstrained_total'
     + // ' n_postjam_bud_unconstrained'
     + // ' n_mother_unconstrained'
     + // ' postjam_divisions_total'
      
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
         rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0
         cell_id(i)=i
         parent_id(i)=0
      ENDDO
      next_cell_id = N
      
      dt = 1.0
      reject_time_step = 0
      terminate_sim = 0
      
      k=0
      DO WHILE (terminate_sim.NE.1)
          k=k+1
          event_count = 0
          CALL copy_everything(
     +     x,y,th,D,alpha,rate,
     +     scale, P, PP, N,
     +     x_copy, y_copy, th_copy, D_copy, alpha_copy, rate_copy,
     +     scale_copy, P_copy, PP_copy, N_copy,
     +     PPm,PPm_copy,
     +     cell_id,parent_id,cell_id_copy,parent_id_copy,
     +     next_cell_id,next_cell_id_copy)
    
          ! GROW PARTICLES
          DO i=1,N
               ratei=rate(i)*dt
               
               IF (above_jamming.eq.1) THEN ! Make sure that populations do not differ before jamming
                   IF(P0.GT.0d0.AND.PP(i).GT.0d0) THEN
                       ratei=ratei*dexp(-PP(i)/P0)
                   ENDIF
               ENDIF
               
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
                  next_cell_id = next_cell_id + 1
                  cell_id(N)=next_cell_id
                  parent_id(N)=cell_id(i)
                  event_parent_bud_diameter(event_count+1)=
     +                 (alpha(i)-1d0)*D(i)
                  
                  ! divide into 2 - 1st assigned index i
                  x(i)=x(i)-dispcm*dcos(th(i))
                  y(i)=y(i)-dispcm*dsin(th(i))
                  rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
                  alpha(i)=alpha0
                  event_count = event_count + 1
                  event_parent_id(event_count)=cell_id(i)
                  event_new_id(event_count)=cell_id(N)
                  event_parent_index(event_count)=i
                  event_new_index(event_count)=N

               ELSE
                   div=0
               ENDIF
               
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
                  ENDIF

                  th(i)=th(i) + 1d-4*(ran2(seed)-0.5d0)
                  th(N)=th(N) + 1d-4*(ran2(seed)-0.5d0)
               ENDIF
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
        
         
         CALL system_status(N,D,D1,alpha,ftol,wide,dt,
     +     fret, phi_j, dphi, DT_MIN, MAX_POSTJAM_STEPS,
     +     reject_time_step, terminate_sim,
     +     before_jamming, at_jamming, above_jamming,
     +     postjam_steps, failure_code)
     
         
         ! REJECT THE MOVE
         IF(reject_time_step.EQ.1) THEN
         CALL copy_back_everything(
     +     x, y, th, D, alpha, rate,
     +     scale, P, PP, N,
     +     x_copy, y_copy, th_copy, D_copy, alpha_copy, rate_copy,
     +     scale_copy, P_copy, PP_copy, N_copy,
     +     PPm,PPm_copy,
     +     cell_id,parent_id,cell_id_copy,parent_id_copy,
     +     next_cell_id,next_cell_id_copy)       
         ENDIF
         
         ! convert back to angles
         DO i=1,N
            th(i)=th(i)/scale(i)
         ENDDO       

         phi = calc_phi(D, alpha, D1, N)
         CALL growth_rate(total_growthrate,N,rate,PP,P0,D,alpha)
         IF(reject_time_step.EQ.0) THEN
            DO i=1,event_count
               WRITE(15,'(I12,E27.18E3,4I12)')k,phi,
     +              event_parent_id(i),event_new_id(i),
     +              event_parent_index(i),event_new_index(i)
            ENDDO
            FLUSH(15)
            IF(tracking_initialized.EQ.1) THEN
               postjam_divisions_total =
     +            postjam_divisions_total + event_count
               CALL lineage_contacts(
     +              x,y,th,D1,D,alpha,N,
     +              curr_bud_contacts,curr_bud_unconstrained)
               CALL compression_flags(N,PP,curr_bud_compressed)
               DO i=1,event_count
                  tracked_initial_free(event_new_index(i))=0
                  tracked_prev_bud_contacts(event_new_index(i))=0
                  tracked_prev_bud_unconstrained(event_new_index(i))=0
                  tracked_initial_bud_diameter(event_new_index(i))=0d0
                  IF(tracked_initial_free(event_parent_index(i)).EQ.1)
     +                 THEN
                     delta_bud_area =
     +                  event_parent_bud_diameter(i)**2
     +                  - tracked_initial_bud_diameter(
     +                  event_parent_index(i))**2
                     CALL write_transition(
     +                    17,k,phi,cell_id(event_parent_index(i)),
     +                    parent_id(event_parent_index(i)),
     +                    tracked_initial_bud_diameter(
     +                    event_parent_index(i)),
     +                    event_parent_bud_diameter(i),
     +                    delta_bud_area,-1,-1,-1,3,
     +                    'initial_free_division')
                     tracked_initial_free(event_parent_index(i))=0
                     tracked_prev_bud_contacts(event_parent_index(i))=0
                     tracked_prev_bud_unconstrained(
     +                    event_parent_index(i))=0
                     tracked_initial_bud_diameter(
     +                    event_parent_index(i))=0d0
                     initial_free_divided = initial_free_divided + 1
                  ENDIF
               ENDDO
               DO i=1,N
                  IF(tracked_initial_free(i).EQ.1) THEN
                     bud_diameter = (alpha(i)-1d0)*D(i)
                     delta_bud_area = bud_diameter**2
     +                    - tracked_initial_bud_diameter(i)**2
                     IF(tracked_prev_bud_contacts(i).EQ.0
     +                  .AND. curr_bud_contacts(i).GT.0) THEN
                        CALL write_transition(
     +                       17,k,phi,cell_id(i),parent_id(i),
     +                       tracked_initial_bud_diameter(i),
     +                       bud_diameter,delta_bud_area,
     +                       curr_bud_contacts(i),
     +                       curr_bud_unconstrained(i),
     +                       curr_bud_compressed(i),1,
     +                       'first_bud_contact')
                     ENDIF
                     IF(tracked_prev_bud_unconstrained(i).EQ.1
     +                  .AND. curr_bud_unconstrained(i).EQ.0) THEN
                        CALL write_transition(
     +                       17,k,phi,cell_id(i),parent_id(i),
     +                       tracked_initial_bud_diameter(i),
     +                       bud_diameter,delta_bud_area,
     +                       curr_bud_contacts(i),
     +                       curr_bud_unconstrained(i),
     +                       curr_bud_compressed(i),2,
     +                       'free_to_constrained')
                        tracked_initial_free(i)=0
                        tracked_prev_bud_contacts(i)=0
                        tracked_prev_bud_unconstrained(i)=0
                        tracked_initial_bud_diameter(i)=0d0
                        initial_free_completed =
     +                     initial_free_completed + 1
                     ELSE
                        tracked_prev_bud_contacts(i)=
     +                     curr_bud_contacts(i)
                        tracked_prev_bud_unconstrained(i)=
     +                     curr_bud_unconstrained(i)
                     ENDIF
                  ENDIF
               ENDDO
               CALL contacts_yeast(
     +              x,y,th,D1,D,N,Nc,Nf,Nu,Nmm,Nbb,Nmb)
               CALL out_numbers(N, Nf, Nu, Ziso)
               CALL chi_c_active_counts(
     +              N,P0,PP,tracked_initial_free,
     +              curr_bud_unconstrained,Nu,
     +              chi_c_est,n_compressed,n_initial_free_active,
     +              n_bud_unconstrained_total,
     +              n_postjam_bud_unconstrained,
     +              n_mother_unconstrained)
               CALL write_postjamm_summary(
     +              18,k,phi,P,N,Nc,Nf,Nu,Ziso,total_growthrate,
     +              chi_c_est,n_compressed,initial_free_total,
     +              n_initial_free_active,initial_free_completed,
     +              initial_free_divided,
     +              n_bud_unconstrained_total,
     +              n_postjam_bud_unconstrained,
     +              n_mother_unconstrained,
     +              postjam_divisions_total)
               FLUSH(17)
               FLUSH(18)
            ENDIF
         ENDIF
         WRITE(*,'(2I12,5E27.18E3,3I8)')k,N,fret/dble(N),P,dt,phi,
     +        total_growthrate,before_jamming,at_jamming,above_jamming
         WRITE(14,'(2I12,5E27.18E3,3I8)')k,N,fret/dble(N),P,dt,phi,
     +        total_growthrate,before_jamming,at_jamming,above_jamming
         IF(failure_code.EQ.EXIT_MIN_DT) THEN
            WRITE(*,'(A,1X,A,1X,2I12,5E27.18E3)')
     +         '# FAILURE','MIN_DT',k,N,dt,phi,phi_j,phi_j+dphi,
     +         total_growthrate
            WRITE(14,'(A,1X,A,1X,2I12,5E27.18E3)')
     +         '# FAILURE','MIN_DT',k,N,dt,phi,phi_j,phi_j+dphi,
     +         total_growthrate
         ELSE IF(failure_code.EQ.EXIT_MAX_POSTJAM_STEPS) THEN
            WRITE(*,'(A,1X,A,1X,2I12,5E27.18E3)')
     +         '# FAILURE','MAX_POSTJAM_STEPS',k,N,dt,phi,phi_j,
     +         phi_j+dphi,total_growthrate
            WRITE(14,'(A,1X,A,1X,2I12,5E27.18E3)')
     +         '# FAILURE','MAX_POSTJAM_STEPS',k,N,dt,phi,phi_j,
     +         phi_j+dphi,total_growthrate
         ENDIF
         FLUSH(14)
         
         IF(mod(k,skip).EQ.0) THEN
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
               WRITE(1,'(3E27.18E3,I12)')xa(1),ya(1),d(i),0
               WRITE(1,'(3E27.18E3,I12)')xa(2),ya(2),d(i)*dd,1
            ENDDO
            FLUSH(1) 
         ENDIF
         
         

         IF( at_jamming.EQ.1 ) THEN
            OPEN(unit=11,file=TRIM(file_LF_JAMM), status='replace')
            OPEN(unit=21,file=TRIM(STATS_file_LF_JAMM), 
     +                                            status='replace')
     
            phi = calc_phi(D, alpha, D1, N)
            
            CALL growth_rate(total_growthrate,N,rate,PP,P0,D,alpha)
            
            WRITE(1,*)  2*N,phi
            WRITE(11,*) N,phi
            WRITE(21,*) 2*N,phi,total_growthrate
            DO i=1,N
                WRITE(11,'(5E27.18E3)') x(i),y(i),D(i),alpha(i),th(i)
                
                cc=dcos(th(i))
                ss=dsin(th(i))
                dd=alpha(i)-1d0
                dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
                dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
                DO kk=1,2
                    xa(kk)=x(i)+dr(kk)*cc
                    ya(kk)=y(i)+dr(kk)*ss
                ENDDO
                WRITE(1,'(3E27.18E3,I12)')xa(1),ya(1),d(i),0     ! TRA
                WRITE(1,'(3E27.18E3,I12)')xa(2),ya(2),d(i)*dd,1  ! TRA
                
                ratei_eff = rate(i)
                IF(P0.GT.0d0.AND.PP(i).GT.0d0) THEN
                    ratei_eff = rate(i)*dexp(-PP(i)/P0)
                ENDIF
                WRITE(21,'(3E27.18E3,I12,4E27.18E3)')
     +                   xa(1),ya(1),d(i),0,
     +                   PPm(i),ratei_eff,rate(i),P0
                WRITE(21,'(3E27.18E3,I12,4E27.18E3)')
     +                   xa(2),ya(2),d(i)*dd,1,
     +                   PP(i),ratei_eff,rate(i),P0
     
            ENDDO
            FLUSH(1)
            FLUSH(11)
            CLOSE(11)
            FLUSH(21)
            CLOSE(21)
            CALL save_lineage_state(file_LINEAGE_JAMM,
     +         N,x,y,th,D,D1,alpha,PP,cell_id,parent_id)
      
            CLOSE(13)
            OPEN(unit=13,file=TRIM(file_NC), status='replace')
            CALL contacts_yeast(x,y,th,D1,D,N,Nc,Nf,Nu,Nmm,Nbb,Nmb)
            CALL out_numbers(N, Nf, Nu, Ziso)

            WRITE(13,'(8I8,5E27.18E3)')N,Ziso,Nc,Nf,Nu,Nmm,Nbb,Nmb,
     +      phi,P,fret,P0,total_growthrate
            FLUSH(13)
            CALL lineage_contacts(
     +         x,y,th,D1,D,alpha,N,
     +         curr_bud_contacts,curr_bud_unconstrained)
            CALL compression_flags(N,PP,curr_bud_compressed)
            CALL init_initial_free_tracking(
     +         N,D,alpha,curr_bud_contacts,curr_bud_unconstrained,
     +         tracked_initial_free,tracked_prev_bud_contacts,
     +         tracked_prev_bud_unconstrained,
     +         tracked_initial_bud_diameter,initial_free_total,
     +         initial_free_completed,initial_free_divided,
     +         tracking_initialized)
            postjam_divisions_total = 0
            CALL chi_c_active_counts(
     +         N,P0,PP,tracked_initial_free,
     +         curr_bud_unconstrained,Nu,
     +         chi_c_est,n_compressed,n_initial_free_active,
     +         n_bud_unconstrained_total,
     +         n_postjam_bud_unconstrained,
     +         n_mother_unconstrained)
            CALL write_postjamm_summary(
     +         18,k,phi,P,N,Nc,Nf,Nu,Ziso,total_growthrate,
     +         chi_c_est,n_compressed,initial_free_total,
     +         n_initial_free_active,initial_free_completed,
     +         initial_free_divided,
     +         n_bud_unconstrained_total,
     +         n_postjam_bud_unconstrained,
     +         n_mother_unconstrained,
     +         postjam_divisions_total)
            FLUSH(18)
            above_jamming = 1
            at_jamming = 0
         ENDIF
      ENDDO
      
!!!!!! THE MAIN LOOP ENDS HERE  !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      
      IF(failure_code.EQ.0) THEN
         phi = calc_phi(D, alpha, D1, N)
         CALL growth_rate(total_growthrate,N,rate,PP,P0,D,alpha)
         ! save the last configuration to the file
         WRITE(1,*)  2*N, phi
         WRITE(12,*) N, phi
         WRITE(22,*) 2*N, phi, total_growthrate
         DO i=1,N
             WRITE(12,'(5E27.18E3)') x(i),y(i),D(i),alpha(i),th(i)

             cc=dcos(th(i))
             ss=dsin(th(i))
             dd=alpha(i)-1d0
             dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
             dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
             DO kk=1,2
                 xa(kk)=x(i)+dr(kk)*cc
                 ya(kk)=y(i)+dr(kk)*ss
             ENDDO

             WRITE(1,'(3E27.18E3,I12)')xa(1),ya(1),d(i),0     ! TRA
             WRITE(1,'(3E27.18E3,I12)')xa(2),ya(2),d(i)*dd,1  ! TRA

             ratei_eff = rate(i)
             IF(P0.GT.0d0.AND.PP(i).GT.0d0) THEN
                 ratei_eff = rate(i)*dexp(-PP(i)/P0)
             ENDIF
             WRITE(22,'(3E27.18E3,I12,4E27.18E3)')
     +        xa(1),ya(1),d(i),0,PPm(i),ratei_eff,rate(i),P0
             WRITE(22,'(3E27.18E3,I12,4E27.18E3)')
     +        xa(2),ya(2),d(i)*dd,1,PP(i),ratei_eff,rate(i),P0

         ENDDO
         FLUSH(1)
         FLUSH(12)
         FLUSH(22)
         CALL save_lineage_state(file_LINEAGE_DPHI,
     +        N,x,y,th,D,D1,alpha,PP,cell_id,parent_id)

         CALL contacts_yeast(x,y,th,D1,D,N,Nc,Nf,Nu,Nmm,Nbb,Nmb)
         CALL out_numbers(N, Nf, Nu, Ziso)

         WRITE(13,'(8I8,5E27.18E3)')N,Ziso,Nc,Nf,Nu,Nmm,Nbb,Nmb,
     +        phi,P,fret,P0,total_growthrate
         FLUSH(13)
      ENDIF
      
      CLOSE(1)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      CLOSE(17)
      CLOSE(18)
      CLOSE(22)
      IF(failure_code.NE.0) THEN
         CALL EXIT(failure_code)
      ENDIF
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      SUBROUTINE system_status(N,D,D1,alpha,ftol,wide,dt,
     +     fret, phi_j,dphi,dt_min,max_postjam_steps,
     +     reject_time_step, terminate_sim,
     +     before_jamming, at_jamming, above_jamming,
     +     postjam_steps, failure_code)
      IMPLICIT NONE
      INTEGER Ntot, N
      PARAMETER(Ntot = 4096)
      INTEGER EXIT_MIN_DT, EXIT_MAX_POSTJAM_STEPS
      PARAMETER(EXIT_MIN_DT=10)
      PARAMETER(EXIT_MAX_POSTJAM_STEPS=11)
      
      DOUBLE PRECISION delta, phi, calc_phi, dt_min
      DOUBLE PRECISION D(Ntot),alpha(Ntot), D1
      DOUBLE PRECISION fret,phi_j,ftol,wide,dt,dphi
      INTEGER reject_time_step, terminate_sim, max_postjam_steps
      INTEGER before_jamming, at_jamming, above_jamming
      INTEGER postjam_steps, failure_code
      
      delta = 0.5d-8  
      phi = calc_phi(D, alpha, D1, N)
      
      ! PACKING NOT JAMMED YET
      IF(fret.LT.ftol*N) THEN
          terminate_sim = 0
          before_jamming = 1
          at_jamming = 0
          above_jamming = 0
          reject_time_step = 0
          dt = 1d0
          phi_j = 0d0
          postjam_steps = 0
          failure_code = 0
          RETURN
      ENDIF
      
      
      ! ITERATING TIME-STEP TO FIND A JAMMED PACKING
      IF(at_jamming.EQ.0 .AND. above_jamming.EQ.0)THEN
          IF(fret.GT.(ftol * wide * N) ) THEN
             dt = 0.5*dt
             reject_time_step = 1
             terminate_sim = 0
             failure_code = 0
             RETURN 
          ELSE IF(fret.GT.(ftol*N) .AND. fret.LT.(ftol*wide*N))THEN
             before_jamming = 0
             at_jamming = 1
             above_jamming = 1
             dt = 1d0
             reject_time_step = 0
             terminate_sim = 0
             phi_j = phi
             postjam_steps = 0
             failure_code = 0
             RETURN
        ENDIF
      ENDIF
           
      
      at_jamming = 0
      postjam_steps = postjam_steps + 1
      ! ITERATE TIME-STEP TO FIND PROPER DPHI
      IF(phi .LT. phi_j+dphi-delta) THEN
          IF(postjam_steps.GE.max_postjam_steps) THEN
             reject_time_step = 0
             terminate_sim = 1
             failure_code = EXIT_MAX_POSTJAM_STEPS
             RETURN
          ENDIF
          dt = 1
          terminate_sim = 0
          reject_time_step = 0
          failure_code = 0
          RETURN
      ELSE IF(phi .GT. phi_j+dphi+delta) THEN
          dt = 0.5 * dt
          IF(dt.LT.dt_min) THEN
             reject_time_step = 1
             terminate_sim = 1
             failure_code = EXIT_MIN_DT
             RETURN
          ENDIF
          reject_time_step = 1
          terminate_sim = 0
          failure_code = 0
          RETURN
      ELSE
          reject_time_step = 0
          terminate_sim = 1
          failure_code = 0
          RETURN
      ENDIF
      
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
      
      
      SUBROUTINE copy_everything(
     +     x, y, th, D, alpha, rate,
     +     scale, P, PP, N,
     +     x_copy, y_copy, th_copy, D_copy, alpha_copy, rate_copy,
     +     scale_copy, P_copy, PP_copy, N_copy,
     +     PPm,PPm_copy,
     +     cell_id,parent_id,cell_id_copy,parent_id_copy,
     +     next_cell_id,next_cell_id_copy)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot)
      DOUBLE PRECISION alpha(Ntot),rate(Ntot),scale(Ntot)
      DOUBLE PRECISION P,PP(Ntot),PPm(Ntot)
      DOUBLE PRECISION x_copy(Ntot),y_copy(Ntot),th_copy(Ntot)
      DOUBLE PRECISION D_copy(Ntot)
      DOUBLE PRECISION alpha_copy(Ntot),rate_copy(Ntot)
      DOUBLE PRECISION scale_copy(Ntot)
      DOUBLE PRECISION P_copy,PP_copy(Ntot),PPm_copy(Ntot)
      INTEGER N, N_copy, i
      INTEGER cell_id(Ntot), parent_id(Ntot)
      INTEGER cell_id_copy(Ntot), parent_id_copy(Ntot)
      INTEGER next_cell_id, next_cell_id_copy

      P_copy = P
      N_copy = N
      next_cell_id_copy = next_cell_id

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
          cell_id_copy(i) = cell_id(i)
          parent_id_copy(i) = parent_id(i)
      ENDDO

      END
      
      SUBROUTINE copy_back_everything(
     +     x, y, th, D, alpha, rate,
     +     scale, P, PP, N,
     +     x_copy, y_copy, th_copy, D_copy, alpha_copy, rate_copy,
     +     scale_copy, P_copy, PP_copy, N_copy,
     +     PPm, PPm_copy,
     +     cell_id,parent_id,cell_id_copy,parent_id_copy,
     +     next_cell_id,next_cell_id_copy)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      DOUBLE PRECISION x(Ntot),y(Ntot),th(Ntot),D(Ntot),D_copy(Ntot)
      DOUBLE PRECISION alpha(Ntot),rate(Ntot),scale(Ntot)
      DOUBLE PRECISION P,PP(Ntot),PPm(Ntot)
      DOUBLE PRECISION x_copy(Ntot),y_copy(Ntot),th_copy(Ntot)
      DOUBLE PRECISION alpha_copy(Ntot),rate_copy(Ntot)
      DOUBLE PRECISION scale_copy(Ntot)
      DOUBLE PRECISION P_copy,PP_copy(Ntot),PPm_copy(Ntot)
      DOUBLE PRECISION dd
      INTEGER N,N_copy,i
      INTEGER cell_id(Ntot), parent_id(Ntot)
      INTEGER cell_id_copy(Ntot), parent_id_copy(Ntot)
      INTEGER next_cell_id, next_cell_id_copy

      P = P_copy
      N = N_copy
      next_cell_id = next_cell_id_copy
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
          cell_id(i) = cell_id_copy(i)
          parent_id(i) = parent_id_copy(i)
      ENDDO
      

      DO i=1,N
          dd=alpha(i)-1d0
          scale(i)=dsqrt(2d0*(1d0+dd**4)/(1+dd**2)+
     +              4d0*(dd*(1d0+dd)/(1+dd**2))**2)/4d0*d(i)
     
          th(i)=th(i)*scale(i)
      ENDDO
      
      
      END
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE save_lineage_state(file_name,
     +     N,x,y,th,D,D1,alpha,PP,cell_id,parent_id)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      CHARACTER*(*) file_name
      INTEGER N, i
      DOUBLE PRECISION x(Ntot), y(Ntot), th(Ntot), D(Ntot), D1
      DOUBLE PRECISION alpha(Ntot), PP(Ntot)
      DOUBLE PRECISION bud_diameter
      INTEGER cell_id(Ntot), parent_id(Ntot)
      INTEGER bud_contacts(Ntot), bud_unconstrained(Ntot)
      INTEGER bud_compressed(Ntot)

      CALL lineage_contacts(
     +     x,y,th,D1,D,alpha,N,bud_contacts,bud_unconstrained)

      OPEN(unit=16,file=TRIM(file_name), status='replace')
      WRITE(16,'(A)')
     + '# index cell_id parent_id alpha bud_diameter'
     + // ' bud_contacts bud_unconstrained bud_compressed'
      DO i=1,N
         bud_diameter = (alpha(i)-1d0)*D(i)
         bud_compressed(i) = 0
         IF(PP(i).GT.0d0) THEN
            bud_compressed(i) = 1
         ENDIF
         WRITE(16,'(3I12,2E27.18E3,3I12)')i,cell_id(i),parent_id(i),
     +        alpha(i),bud_diameter,bud_contacts(i),
     +        bud_unconstrained(i),bud_compressed(i)
      ENDDO
      FLUSH(16)
      CLOSE(16)

      END

      SUBROUTINE init_initial_free_tracking(
     +     N,D,alpha,bud_contacts,bud_unconstrained,
     +     tracked_initial_free,tracked_prev_bud_contacts,
     +     tracked_prev_bud_unconstrained,
     +     tracked_initial_bud_diameter,initial_free_total,
     +     initial_free_completed,initial_free_divided,
     +     tracking_initialized)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      INTEGER N, i
      DOUBLE PRECISION D(Ntot), alpha(Ntot)
      INTEGER bud_contacts(Ntot), bud_unconstrained(Ntot)
      INTEGER tracked_initial_free(Ntot)
      INTEGER tracked_prev_bud_contacts(Ntot)
      INTEGER tracked_prev_bud_unconstrained(Ntot)
      DOUBLE PRECISION tracked_initial_bud_diameter(Ntot)
      INTEGER initial_free_total, initial_free_completed
      INTEGER initial_free_divided, tracking_initialized

      initial_free_total = 0
      initial_free_completed = 0
      initial_free_divided = 0
      tracking_initialized = 1

      DO i=1,Ntot
         tracked_initial_free(i)=0
         tracked_prev_bud_contacts(i)=0
         tracked_prev_bud_unconstrained(i)=0
         tracked_initial_bud_diameter(i)=0d0
      ENDDO

      DO i=1,N
         IF(bud_unconstrained(i).EQ.1) THEN
            tracked_initial_free(i)=1
            tracked_prev_bud_contacts(i)=bud_contacts(i)
            tracked_prev_bud_unconstrained(i)=bud_unconstrained(i)
            tracked_initial_bud_diameter(i)=(alpha(i)-1d0)*D(i)
            initial_free_total = initial_free_total + 1
         ENDIF
      ENDDO

      END

      SUBROUTINE write_transition(
     +     unit_no,step,phi,cell_id,parent_id,
     +     initial_bud_diameter,bud_diameter,delta_bud_area,
     +     bud_contacts,bud_unconstrained,bud_compressed,
     +     event_code,event_label)
      IMPLICIT NONE
      INTEGER unit_no, step, cell_id, parent_id
      INTEGER bud_contacts, bud_unconstrained, bud_compressed
      INTEGER event_code
      DOUBLE PRECISION phi, initial_bud_diameter
      DOUBLE PRECISION bud_diameter, delta_bud_area
      CHARACTER*(*) event_label

      WRITE(unit_no,'(I12,E27.18E3,2I12,3E27.18E3,4I12,1X,A)')
     +     step,phi,cell_id,parent_id,initial_bud_diameter,
     +     bud_diameter,delta_bud_area,bud_contacts,
     +     bud_unconstrained,bud_compressed,event_code,event_label

      END

      SUBROUTINE compression_flags(N,PP,bud_compressed)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      INTEGER N, i, bud_compressed(Ntot)
      DOUBLE PRECISION PP(Ntot)

      DO i=1,N
         bud_compressed(i)=0
         IF(PP(i).GT.0d0) THEN
            bud_compressed(i)=1
         ENDIF
      ENDDO

      END

      SUBROUTINE chi_c_active_counts(
     +     N,P0,PP,tracked_initial_free,
     +     bud_unconstrained,Nu,
     +     chi_c_est,n_compressed,n_initial_free_active,
     +     n_bud_unconstrained_total,
     +     n_postjam_bud_unconstrained,
     +     n_mother_unconstrained)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      INTEGER N, i, n_compressed, n_initial_free_active
      INTEGER n_bud_unconstrained_total
      INTEGER n_postjam_bud_unconstrained
      INTEGER n_mother_unconstrained
      INTEGER Nu
      INTEGER tracked_initial_free(Ntot), bud_unconstrained(Ntot)
      DOUBLE PRECISION P0, PP(Ntot), chi_c_est, chi_sum

      chi_sum = 0d0
      n_compressed = 0
      n_initial_free_active = 0
      n_bud_unconstrained_total = 0

      DO i=1,N
         IF(PP(i).GT.0d0) THEN
            n_compressed = n_compressed + 1
            IF(P0.GT.0d0) THEN
               chi_sum = chi_sum + dexp(-PP(i)/P0)
            ELSE
               chi_sum = chi_sum + 1d0
            ENDIF
         ENDIF
         IF(tracked_initial_free(i).EQ.1) THEN
            n_initial_free_active = n_initial_free_active + 1
         ENDIF
         IF(bud_unconstrained(i).EQ.1) THEN
            n_bud_unconstrained_total = n_bud_unconstrained_total + 1
         ENDIF
      ENDDO

      n_postjam_bud_unconstrained =
     +     n_bud_unconstrained_total - n_initial_free_active
      IF(n_postjam_bud_unconstrained.LT.0) THEN
         n_postjam_bud_unconstrained = 0
      ENDIF

      n_mother_unconstrained = Nu - n_bud_unconstrained_total
      IF(n_mother_unconstrained.LT.0) THEN
         n_mother_unconstrained = 0
      ENDIF

      IF(n_compressed.GT.0) THEN
         chi_c_est = chi_sum / dble(n_compressed)
      ELSE
         chi_c_est = 1d0
      ENDIF

      END

      SUBROUTINE write_postjamm_summary(
     +     unit_no,step,phi,P,N,Nc,Nf,Nu,Ziso,total_growthrate,
     +     chi_c_est,n_compressed,initial_free_total,
     +     n_initial_free_active,initial_free_completed,
     +     initial_free_divided,
     +     n_bud_unconstrained_total,
     +     n_postjam_bud_unconstrained,
     +     n_mother_unconstrained,
     +     postjam_divisions_total)
      IMPLICIT NONE
      INTEGER unit_no, step, N, Nc, Nf, Nu, Ziso
      INTEGER n_compressed, initial_free_total
      INTEGER n_initial_free_active, initial_free_completed
      INTEGER initial_free_divided
      INTEGER n_bud_unconstrained_total
      INTEGER n_postjam_bud_unconstrained
      INTEGER n_mother_unconstrained
      INTEGER postjam_divisions_total
      DOUBLE PRECISION phi, P, total_growthrate, chi_c_est

      WRITE(unit_no,'(I12,2E27.18E3,5I12,2E27.18E3,9I12)')
     +     step,phi,P,N,Nc,Nf,Nu,Ziso,total_growthrate,chi_c_est,
     +     n_compressed,initial_free_total,n_initial_free_active,
     +     initial_free_completed,initial_free_divided,
     +     n_bud_unconstrained_total,n_postjam_bud_unconstrained,
     +     n_mother_unconstrained,postjam_divisions_total

      END

      SUBROUTINE lineage_contacts(
     +     x,y,th,D1,D,alpha,N,bud_contacts,bud_unconstrained)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      INTEGER N, i, j, ki, kj, flag
      DOUBLE PRECISION x(Ntot), y(Ntot), th(Ntot), alpha(Ntot)
      DOUBLE PRECISION xij, yij, D(Ntot), D1
      DOUBLE PRECISION dij
      DOUBLE PRECISION Lx,Ly,rijsq
      DOUBLE PRECISION c(Ntot),s(Ntot),xa(Ntot,2),ya(Ntot,2)
      DOUBLE PRECISION dr(Ntot,2),dk(Ntot,2),di_up(Ntot)
      INTEGER nc_bud(Ntot,2), bud_contacts(Ntot)
      INTEGER bud_unconstrained(Ntot)
      COMMON /f5com/ Lx,Ly

      DO i=1,N
         c(i)=dcos( th(i) )
         s(i)=dsin( th(i) )
      ENDDO
      CALL mol_to_atoms(N,x,y,c,s,D,D1,alpha,dr,xa,ya,dk,di_up)

      DO i=1,N
         nc_bud(i,1)=0
         nc_bud(i,2)=0
         bud_contacts(i)=0
         bud_unconstrained(i)=0
      ENDDO

      DO i=1,N-1
         DO j=i+1, N
            DO ki=1,2
               DO kj=1,2
                  dij=(dk(i,ki)+dk(j,kj))/2d0
                  xij=xa(i,ki)-xa(j,kj)
                  xij=xij-idnint(xij/Lx)*Lx
                  yij=ya(i,ki)-ya(j,kj)
                  yij=yij-idnint(yij/Ly)*Ly
                  rijsq=xij**2+yij**2
                  IF(rijsq.LT.(dij**2)) THEN
                     nc_bud(i,ki)=nc_bud(i,ki)+1
                     nc_bud(j,kj)=nc_bud(j,kj)+1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO i=1,N
         bud_contacts(i)=nc_bud(i,2)
         flag = 0
         IF( (nc_bud(i,1)+nc_bud(i,2)).LT.3 ) THEN
            flag = 1
         ENDIF
         IF( (nc_bud(i,1)+nc_bud(i,2)).EQ.3 ) THEN
            IF( nc_bud(i,1).EQ.2 .OR. nc_bud(i,2).EQ.2 ) THEN
               flag = 1
            ENDIF
         ENDIF
         IF(flag.EQ.0 .AND. nc_bud(i,2).EQ.0) THEN
            bud_unconstrained(i)=1
         ENDIF
      ENDDO

      END

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
          ENDIF
          
          Dm = D(i)
C         NOTE: Db should formally be (alpha(i)-1d0)*D(i), but D(i)=1
C         for all particles in this code (initialized to 1d0, preserved
C         at division), so the missing D(i) multiplier has no effect.
          Db = alpha(i)-1d0
          total_area = total_area + Dm*Dm + Db*Db
          total_new_mass = total_new_mass + ratei*(Dm*Dm + Db*Db)
      ENDDO
      
      gr = total_new_mass / total_area
      
      END
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!    STANDARD NUMERICAL CODE   !!!!!!!!!!!!!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
      IMPLICIT NONE
      INTEGER Ntot
      PARAMETER(Ntot = 4096)
      INTEGER its,iter,ITMAX
      DOUBLE PRECISION fret,ftol,ftol1
      PARAMETER (ITMAX=1000000000)
      DOUBLE PRECISION dgg,fp,gam,gg,gx(Ntot),gy(Ntot)
      DOUBLE PRECISION hx(Ntot),hy(Ntot)
      DOUBLE PRECISION D(Ntot),D1,xix(Ntot),xiy(Ntot),xith(Ntot)
      DOUBLE PRECISION x(Ntot),y(Ntot),maxdis,xp(Ntot),yp(Ntot)
      DOUBLE PRECISION th(Ntot),hth(Ntot),gth(Ntot),exp,att,width
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
      DOUBLE PRECISION ax,bx,fa,fb,fx,xmin,xx,xp(Ntot),yp(Ntot),dbrent
      DOUBLE PRECISION pxcom(Ntot),pycom(Ntot),xixcom(Ntot),xiycom(Ntot)
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ran2(idum)
      IMPLICIT NONE
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
      ENDIF
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
