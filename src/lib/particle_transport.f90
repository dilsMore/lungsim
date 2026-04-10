! Description:
! This module contains code specific to particle transport problems

module particle_transport

  use other_consts
  use indices
  use arrays
  implicit none
  
  !Module Parameters

  !Module types

  !Module variables

  !Interfaces
  
  public solve_particles_decoupled
  
  public where_inlist, inlist
  public calc_mass_particles
  public write_terminal, write_airway

  !public branch_length


contains

  subroutine solve_particles_decoupled(initial_concentration, inlet_concentration, particle_size)

!!!##########################################################################
!!!##########################################################################
!!!##########################################################################

!!! Assemble matrices for 1D particle transport, and solve. Based on the gas mixing model
!!! with adaptations by Falko Schmidt for particles.
    use arrays
    use exports
    use geometry
    use indices
    use solve,only: pmgmres_ilu_cr
    use species_transport, only: set_elem_volume,initialise_transport,reduce_transport_matrix
    use other_consts, only: MAX_FILENAME_LEN
    use ventilation, only: read_params_evaluate_flow
    use diagnostics, only: enter_exit

    implicit none

    real(dp), intent(in) :: initial_concentration
    real(dp), intent(in) :: inlet_concentration
    real(dp), intent(in) :: particle_size
  
    type(particle_parameters) :: part_param
    type(transport_parameters) :: tp
    integer :: fileid = 10
    real(dp) :: err
    character(len=9) :: problem_type = 'particles'
    character(len=200) :: op_name

    ! Local variables
    real(dp),allocatable :: solution(:)

    integer :: MatrixSize,nonzeros,ncol,nentry, &
         noffset_entry,noffset_row,np,nrow,nrow_BB
    real(dp) :: AA,BB,current_mass,current_volume,inlet_flow, &
         mass_deposit,mass_error,theta,time,volume_error,volume_tree,&
         cpu_dt_start,cpu_dt_end,temp_mass
    
    real(dp) :: ttime, turnovers, time_start, time_end,ccun
    integer :: Gdirn
    real(dp) :: chest_wall_compliance,constrict, COV, FRC, i_to_e_ratio, pmus_step, press_in,&
       refvol, RMaxMean, RMinMean, T_interval, volume_target,undef
    character :: expiration_type*(10)

    integer :: Nexp, nstep, nbreath, j
    logical :: carryon
    character(len=60) :: sub_name
    logical :: last_breath = .false.
    logical :: write_mass = .false.
    character(len=4) :: char_int
    character(len=MAX_FILENAME_LEN) :: file_export
    real(dp) :: Eo, inlet_mouth_concentration, dep_eff_total_prev
    ! TJ -  Check if it's the SINGLE inhalation phase
    logical :: expiration, inspiration
    real(dp) :: extra_mass, last_extra_mass_inspiration, extra_mass_inspiration, dep_eff_total
    real(dp) :: conc_exp, conc_exp_prev
    sub_name = 'solve_particle_decoupled'
    call enter_exit(sub_name,1)
    call initialise_transport(initial_concentration,inlet_concentration,tp)
    
    !Set up parameters, NB any of these parameters could be passed from python in the long term
    part_param%num_brths_gm = 12
    part_param%solve_tolerance = 1.0e-8_dp
    part_param%diffusion_coeff = 22.5_dp ! mm
    part_param%pdia = particle_size * 1.0e-3_dp  ! micron --> mm
    part_param%dt_gm = 0.01_dp
    part_param%VtotTLC = 186.89_dp
    part_param%totacinarLength = 7.73_dp
    part_param%mu = 18.69e-6_dp  ! Pa.s = g/mm/s
    part_param%prho =  1.0e-3_dp             ! density of particles [g/mm^3]
    part_param%lambda = 7.022e-5_dp  ! [mm] mean free path necessary
    part_param%kBoltz = 1.38e-14_dp  ! ! Boltzmann constant [J/K*1d9]=[kg*m^2/s^2/K*1d9]=[g*mm^2/s^2/K]
    part_param%Temperature = 36.0_dp+273.15_dp ! Temperature [K] from rho*R*T
    part_param%out_itr_max = 200      ! max # (outer) iterations using GMRES solver.
    part_param%inr_itr_max = 100      ! max # (inner) iterations using GMRES solver.
    !!! assume a generic 9-generation acinus geometry, based on Haefeli-Bleuer & Weibel, 1988
!!! radius acinar ducts at TLC for 9 generations (HAEFELI-BLEUER, 1988)
    part_param%RacTLC = (/0.25_dp,0.25_dp,0.245_dp,0.2_dp,0.19_dp,0.18_dp,0.17_dp,0.155_dp,0.145_dp,0.125_dp/)
!!! length acinar ducts for terminal bronchiole 9 generations (HAEFELI-BLEUER, 1988)
    part_param%LacTLC = (/0.8_dp,1.33_dp,1.12_dp,0.93_dp,0.83_dp,0.7_dp,0.7_dp,0.7_dp,0.67_dp,0.75_dp/)
!!! acinar volume at TLC for 9 generations (estimated from HAEFELI-BLEUER, 1988)
    part_param%VacTLC = (/0.87_dp,1.70_dp,2.49_dp,4.81_dp,7.56_dp,14.07_dp,24.60_dp,43.93_dp,86.86_dp/) ! mm^3
!!! Cunningham slip correction factor 
    Ccun = 1.0_dp + 2.0_dp * part_param%lambda/part_param%pdia &
         * (1.257_dp+0.4_dp*exp(-0.55_dp*part_param%pdia/part_param%lambda))

    part_param%diffu = part_param%kBoltz*part_param%Temperature*Ccun/&
            (3.0_dp*pi*part_param%mu*part_param%pdia) ! diffusion constant D_B [mm^2/s]


    call read_params_evaluate_flow(Gdirn, chest_wall_compliance, &
       constrict, COV, FRC, i_to_e_ratio, pmus_step, press_in,&
       refvol, RMaxMean, RMinMean, T_interval, volume_target, expiration_type)

    if(GDirn.eq.1) then 
      part_param%gravityx = 9.81_dp*1.0e3_dp
      part_param%gravityy = 0.0_dp
      part_param%gravityz = 0.0_dp
      part_param%grav_factor = 1
    else if(GDirn.eq.2) then
      part_param%gravityx = 0.0_dp
      part_param%gravityy = 9.81_dp*1.0e3_dp
      part_param%gravityz = 0.0_dp
      part_param%grav_factor = 1
    else if (GDirn.eq.3) then !Needs checking, may not always be consistent with mesh
      part_param%gravityx = 0.0_dp
      part_param%gravityy = 0.0_dp
      part_param%gravityz = 9.81_dp*1.0e3_dp ! *MHT change*
      part_param%grav_factor = 1
    else if (GDirn.eq.0) then
        part_param%gravityx = 0.0_dp
        part_param%gravityy = 0.0_dp
        part_param%gravityz = 0.0_dp ! *MHT change*
        part_param%grav_factor = 0
    else
     print *, "ERROR: Posture not recognised (currently only x=1,y=2,z=3))"
     call exit(0)
    endif

    ! TJ - uncomment for zero gravity
    print *, 'Gravity Factor', part_param%grav_factor

    part_param%time_inspiration = 2.0_dp
    part_param%time_breath_hold = 1.0_dp
    part_param%time_expiration = 2.0_dp
    part_param%tidal_volume = volume_target! tidal volume target, mm^3
    part_param%FRC = FRC
    part_param%initial_volume = FRC*1.0e+6_dp !initial volume of air in lungs, mm^3

    call set_elem_volume()

    ! TJ - use this volume_of_mesh() to get deadspace volume
    call volume_of_mesh(part_param%initial_volume,volume_tree) ! to get deadspace volume
    ! TJ - NOW part_param%initial_volume = volume_tree = deadspace volume
    ! TJ - TO scale deadspace volume, go to volume_of_mesh() in geometry.f90

!!! distribute the initial tissue unit volumes along the gravitational axis.
    call set_initial_volume(Gdirn,COV,FRC*1.0e+6_dp,RMaxMean,RMinMean)
    undef = refvol * (FRC*1.0e+6_dp-volume_tree)/DBLE(elem_units_below(1))

!!! calculate the total model volume
    ! TJ - use this volume_of_mesh to reset initial volume to FRC, comment set_initial_volume() in .py
    call volume_of_mesh(part_param%initial_volume, volume_tree)

    write(*,'('' Anatomical deadspace = '',F8.3,'' ml'')') volume_tree/1.0e+3_dp ! in mL
    write(*,'('' Respiratory volume   = '',F8.3,'' L'')') (part_param%initial_volume-volume_tree)/1.0e+6_dp !in L
    write(*,'('' Total lung volume    = '',F8.3,'' L'')') part_param%initial_volume/1.0e+6_dp !in L
    write(*,'('' Inlet flow           = '',f8.3,'' L.s^-1'')') abs(elem_field(ne_Vdot,1))/1.0e+6_dp
    write(*,'('' Inlet concentration (mouth) = '',f8.3,'' g.mm^-3'')') node_field(nj_conc1,1)
    write(*,'('' Inlet concentration (lung) = '',f8.5,'' g.mm^-3'')') tp%inlet_concentration(1)
    write(*,'('' Particle size        = '',f8.3,'' micron m^-3'')') particle_size
    write(*,'('' Particle size        = '',f10.7,'' mm^-3'')') part_param%pdia !(part_param%pdia * 1.0e+3_dp)
    write(*,'('' ccun   = '',f10.7,'' '')') Ccun
    write(*,'('' Diffusion constant   = '',f10.8,'' mm^2.s^-1'')') (part_param%diffu)

    !write(*,'('' Extrathoracic deposition   = '',f8.5,'' '')') part_param%Eo

     op_name = 'file_particle' !ARC TEMP placeholder

     !!! ###########  INITIAL & BOUNDARY CONDITIONS FOR GAS MIXING & EXCHANGE   ###########
     ! sets initial concentration at all nodes, concentration at the inlet node
     ! (for inspiration), and allocates memory for solution matrices.
     turnovers = part_param%num_brths_gm * part_param%tidal_volume / &
                 (part_param%FRC*1.0e+6_dp)

     if(turnovers.lt.6.0_dp)then
        write(*,'(''Warning:'',i3,'' breaths equals'',f5.2,'' turnovers;'// &
           ' set num_brths_gm to'',i3,'' or enter return to continue:'')') &
           part_param%num_brths_gm,turnovers,int(6.0_dp*(part_param%FRC*1.0e+6_dp)/&
           part_param%tidal_volume +1.0_dp)
        read(*,*)
     endif

!!! Eo was here

    if(write_mass)then
       write(*,'(''  Time|   Flow|    dVol|   Vol|  IdlMass|   Total|BrnLumen|BrnDepos|AlvLumen|AlvDepos|'',&
            &''   volEr|  MassEr|    C_inlet|  Mass by gen'')')
       write(*,'(''   (s)|  (L/s)|     (L)|   (L)|    (g)  |    Mass|    Mass|    Mass|    Mass|    Mass|'',&
            &''    %   |    %   |  (g.mm^-3)|       (g)   '')')
    else
       write(*,'(58x, ''||---------------------------------- Mass (g) in:---------------------------------||'')')
       write(*,'(''  Time|   dVol|     Vol ( %Err) |IdlMass|    Mass ( %Err) |'',&
            &'' Extra| Bronch| Bronch| -diff-| -sedi-| -impc-| Alveol| Alveol| '',&
            &''-diff-| -sedi-  | DF lung| DF brn| DF alv|   C_inlet '')')
       write(*,'(''   (s)|    (L)|     (L)         |   (g) |     (g)         |'',&
            &''      | lumen|   wall|       |       |       |  lumen|   wall|'',&
            &''       |         |       |       |       |   |       |'')')
    endif

    time_end = 0.0_dp
    time_start = 0.0_dp
    ttime = 0.0_dp
    time = 0.0_dp


! Initialize extra_mass and last_extra_mass_inspiration
    extra_mass = 0.0_dp
    last_extra_mass_inspiration = 0.0_dp

    ! Initialize nstep to keep track of the number of breaths
    nstep = 0

    !###########################################################
      do nbreath = 1, part_param%num_brths_gm !1!ARC TEMP
          !!! Inhalation
        inlet_flow = abs(elem_field(ne_Vdot, 1))
        call scale_flow_field(inlet_flow)

        time_start = time_end
        time_end = time_start + part_param%time_inspiration

        write(*,'('' Inhaling breath'',i4,'' from'',f8.3,'' s to'',f8.3,'' s'')') nbreath, &
                  time_start, time_end

        !node_field(nj_conc1,1) = tp%inlet_concentration(1) ! need to set here

!       ! Modify inlet concentration only during the first inhalation
        part_param%Eo = calculate_Eo(tp, part_param, .true., particle_size)
        print *, 'Eo (inhalation) = ', part_param%Eo
!
!        ! Apply modification based on Eo for inhalation
!
!        tp%inlet_concentration(:) = tp%inlet_concentration(:) * (1.0_dp - part_param%Eo)
!        print *, 'inhalation concentration', tp%inlet_concentration(1)

!       ! Apply modification based on Eo for the first inhalation
        if (nbreath == 1) then
            tp%inlet_concentration(:) = tp%inlet_concentration(:) * (1.0_dp - part_param%Eo)
            print *, 'inhalation concentration', tp%inlet_concentration(1)
        endif

        node_field(nj_conc1,1) = tp%inlet_concentration(1)! * (1.0_dp - part_param%Eo)! need to set here
        print *, 'inhalation concentration after modification', node_field(nj_conc1,1)

        call solve_particles(fileid, time_end, time_start, .true., last_breath, tp, part_param, &
        write_mass, extra_mass_inspiration, last_extra_mass_inspiration, dep_eff_total)

        ! Apply modification based on Eo for inhalation
        !tp%inlet_concentration(:) = tp%inlet_concentration(:) * (1.0_dp - part_param%Eo)
        !print *, 'inhalation concentration', tp%inlet_concentration(1)

        if (inspiration) then
            ! Accumulate extra_mass during both inhalation and exhalation
            !extra_mass_inspiration = extra_mass
            last_extra_mass_inspiration = extra_mass_inspiration + extra_mass
            if (nbreath > 1) extra_mass_inspiration = last_extra_mass_inspiration
        else
            ! For exhalation, maintain the same extra_mass from the previous inhalation
            extra_mass = last_extra_mass_inspiration
            ! Update extra_mass_inspiration with the value from the previous inhalation
            extra_mass_inspiration = last_extra_mass_inspiration
        endif

        !!!  ! Breath Hold
        ! if (part_param%time_breath_hold.gt.0.0_dp) then
        !     do j = 1, part_param%n_export
        !         time_start = time_end
        !         time_end = time_start + part_param%time_breath_hold
        !
        !         call solve_particles(fileid, time_end, time_start, .true., last_breath, tp, part_param)
        !
        !         Nexp = Nexp + 1
        !         write(char_int, '(i3)') Nexp
        !
        !         file_export = trim(op_name) // '_conc_' // trim(adjustl(char_int))
        !         !ARC TEMP call export_node_field(nj_conc1, file_export, part_param%group_name, 'conc')
        !     enddo
        ! endif

        !!! Exhalation
        inlet_flow = -abs(elem_field(ne_Vdot,1))
        call scale_flow_field(inlet_flow)

        time_start = time_end
        time_end = time_start + part_param%time_expiration

        write(*,'('' Exhaling breath'',i4,'' from'',f8.3,'' s to'',f8.3,'' s'')') nbreath, &
              time_start, time_end

        !inlet_flow = -250.0e+3_dp ! 250 ml/s, for 2 s
        !call scale_flow_field(inlet_flow)

        ! tj - MODIFY CONCENTRATION HERE -
!        if(.not. inspiration) tp%inlet_concentration(1)  = tp%inlet_concentration(1) * (1.0_dp - dep_eff_total)
!        node_field(nj_conc1,1) = tp%inlet_concentration(1)!
        print *, 'exhalation concentration modified', node_field(nj_conc1,1)

        call solve_particles(fileid, time_end, time_start, .false., last_breath, tp, part_param, &
                write_mass, last_extra_mass_inspiration , extra_mass_inspiration, dep_eff_total)

        !!! Calculate Eo for exhalation
        part_param%Eo = calculate_Eo(tp, part_param, .false., particle_size)
        ! Print Eo for exhalation
        print *, 'Eo (exhalation) = ', part_param%Eo
        ! Note: You can optionally modify the concentration during exhalation here...?
        ! Apply modification based on Eo for exhalation
        conc_exp = node_field(nj_conc1,1)
        print *, 'exhalation concentration', conc_exp ! node_field(nj_conc1,1)
        nstep = nbreath

        ! Check if it's not the first breath
        if (nbreath > 1) then
            ! Check if the absolute difference between dep_eff_total values falls below the threshold
            if (abs(dep_eff_total - dep_eff_total_prev) < 0.01) then
                ! If the condition is met, print a message and exit the loop
                print *, "Steady state reached. Simulation stopped."
                exit
            elseif (abs((conc_exp - conc_exp_prev)/conc_exp) < 0.01) then
                ! If the condition is met, print a message and exit the loop
                print *, "Concentration reached Steady state. Simulation stopped.", &
                        abs((conc_exp - conc_exp_prev)/conc_exp)
                exit
            endif
        endif

        print *, "conc change %", abs((conc_exp - conc_exp_prev)/conc_exp)*100
        ! Store the current dep_eff_total for comparison in the next iteration
        dep_eff_total_prev = dep_eff_total

        conc_exp_prev = conc_exp
        ! Increment the breath count
        nstep = nstep + 1

    enddo !Nbreath

    call enter_exit(sub_name,2)

  end subroutine solve_particles_decoupled

    !###################################################################################

function calculate_Eo(tp, part_param, inspiration,  particle_size) result(Eo)
    type(transport_parameters) :: tp
    type(particle_parameters):: part_param
    logical, intent(in) :: inspiration
    real(dp) :: Eo
    real(dp), intent(in) :: particle_size

    if (inspiration) then
        ! M2: Cheng 2003 - extra thoracic airway (Cheng, 2003), D = 0.022 cm^2/s
        ! 1 - exp(-0.000278_dp*(da**2)*Q - 20.4_dp*D**0.66*Q**(-0.31))
        part_param%inlet_flow = abs(elem_field(ne_Vdot,1))/1.0e+3_dp*0.06_dp ! [L/min]
        write(*,'('' inlet_flow   = '',f8.5,'' L.min^-1 '')') part_param%inlet_flow
        Eo = 1 - exp(-0.000278_dp*(particle_size**2)* part_param%inlet_flow - &
            20.4_dp* (part_param%diffu/100)**(0.66_dp) * (part_param%inlet_flow**(-0.31_dp)) )
    else
        ! Skip extra-thoracic calculation for multiple inhalations
        !Eo = 0.0_dp
        ! TJ - MAR 28 2024 - TEST result IF there is new deposition during exhale.
        Eo = 1 - exp(-0.000278_dp*(particle_size**2)* node_field(nj_conc1,1) - &
            20.4_dp* (part_param%diffu/100)**(0.66_dp) * (node_field(nj_conc1,1)**(-0.31_dp)) )
    end if
end function calculate_Eo

!###################################################################################

  subroutine solve_particles(fileid,time_end,time_start,inspiration,last_breath,&
    tp,part_param,write_mass, last_extra_mass_inspiration, extra_mass_inspiration, dep_eff_total)

!!! Assemble matrices for 1D particle transport, and solve. Based on the gas mixing model
!!! with adaptations by Falko Schmidt for particles.

!   use arrays, only: dp,particle_parameters,part_acinus_field,node_field,num_nodes,num_elems
    use arrays
    use indices
    use exports!,only: export_node_field
    use geometry, only: volume_of_mesh
    use solve,only: pmgmres_ilu_cr
    use species_transport, only: assemble_transport_matrix,reduce_transport_matrix,calc_mass
    use other_consts
    use diagnostics, only: enter_exit
    
    implicit none
  
    type(particle_parameters) :: part_param
    type(transport_parameters) :: tp
    integer,intent(in) :: fileid
    real(dp),intent(in) :: time_end,time_start
    logical,intent(in) :: inspiration,last_breath,write_mass

    real(dp) :: err

    ! Local variables
    real(dp),allocatable :: solution(:)

    integer :: i,MatrixSize,nonzeros,ncol,nentry, &
         noffset_entry,noffset_row,np,nrow,nrow_BB
    real(dp) :: AA,BB,current_mass,current_volume,dep_eff_alv,dep_eff_bronch,inlet_flow, &
         lumen_mass,mass_by_gen(50),mass_deposit,mass_error,theta,time,volume_error,volume_tree,&
         cpu_dt_start,cpu_dt_end,temp_mass,unit_mass,unit_wall,last_mass


    real(dp) :: Ccun
    real(dp) :: dt
    logical :: carryon
    character(len=60) :: sub_name, name, groupname, field_name
    integer :: ne_field
    logical :: debug_on = .true., deposition_on = .true., diffusion_on = .true.

    integer :: SOLVER_FLAG

    real(dp) :: extra_mass, dep_eff_extra
    logical :: write_airway_on = .true., write_terminal_on = .true.
    character(len=60) :: export_directory, filename, EXELEMFIELD, EXNODEFIELD
    integer :: ne_depos

    real(dp), intent(in) :: last_extra_mass_inspiration
    real(dp), intent(out) :: extra_mass_inspiration, dep_eff_total

    character(len=5) :: file_id

    character(len=9) :: problem_type = 'particles'

    integer :: num_iter
   real(dp)::          DepFrac(4:6)

    ! #############################################################################

!   set_diagnostics = .true.
    name = 'terminal'
    groupname = 'vent_model'
    field_name = 'bronchial_deposition'
    ne_field = 13

    sub_name = 'solve_particles'
    call enter_exit(sub_name,1)

    ! allocatable array to store the current solution 
    if(.not.allocated(solution)) allocate(solution(num_nodes))
!   if(.not.allocated(M_C_M)) allocate(M_C_M(num_nodes))

    if(.not.allocated(part_acinus_field))then
       allocate(part_acinus_field(20,num_units)) !ARC allocating a 20 long aray?
       part_acinus_field(:,:) = 0.0_dp
    endif
    
    !ARC particle acinus field, should be a 'unit field'

    inlet_flow = elem_field(ne_Vdot,1) ! flow at entry element
    theta = 2.0_dp/3.0_dp
    dt = part_param%dt_gm
    ! get the sparsity arrays for the reduced system. uses compressed row format.
    call reduce_transport_matrix(MatrixSize,nonzeros,noffset_entry,noffset_row,inspiration)
    time = time_start ! initialise the time
    carryon = .true. ! logical for whether solution continues
    current_mass = 0.0_dp
    extra_mass = 0.0_dp
    ! Initialise extra_mass_inspiration
    extra_mass_inspiration = 0.0_dp


    ! main time-stepping loop:  time-stepping continues while 'carryon' is true
    do while (carryon) !
       
       call cpu_time(cpu_dt_start)
       time = time + dt ! increment time
       num_iter = (time + dt - time_start)/dt ! should be renamed
       if(inlet_flow.gt.0.0_dp) node_field(nj_conc1,1) = tp%inlet_concentration(1) ! reset the inlet concentration
       if(abs(inlet_flow).gt.loose_tol)then ! i.e. not for breath hold
          call airway_mesh_deform(dt,part_param%initial_volume,part_param%coupled,'particles',tp,part_param) ! change model size by dV
          call general_track(dt,.true.,tp,part_param)
          call particle_velocity(time,part_param)
       endif

       !tp%ideal_mass = tp%ideal_mass + inlet_flow*dt*tp%inlet_concentration(1)
       tp%ideal_mass = tp%ideal_mass + inlet_flow*dt*node_field(nj_conc1,1)

       ! assemble the element matrices. Element matrix calculation can be done directly 
       ! (based on assumption of interpolation functions) or using Gaussian interpolation.
       call assemble_transport_matrix(part_param%diffusion_coeff) ! also for particles

       ! initialise the values in the solution matrices
       global_AA(1:nonzeros) = 0.0_dp ! equivalent to M in Tawhai thesis
       global_BB(1:num_nodes) = 0.0_dp ! equivalent to K in Tawhai thesis
       global_R(1:num_nodes) = 0.0_dp


       !extra_mass = extra_mass + inlet_flow*dt*node_field(nj_conc1,1) * part_param%Eo
       extra_mass = extra_mass + abs(inlet_flow)*dt*node_field(nj_conc1,1) * part_param%Eo
       extra_mass_inspiration = extra_mass + last_extra_mass_inspiration

       !!! TJ - this is assuming no extra-thoracic deposition for exhale
!       if (inspiration) then
!           ! Calculate extra_mass as before
!           extra_mass_inspiration = extra_mass + last_extra_mass_inspiration
!       else
!           extra_mass = extra_mass_inspiration
!       endif

       ! Assemble the reduced system of matrices
       do np=1,num_nodes ! Loop over rows of unreduced system
          nrow = np ! conveniently true for the way we set up our models
          ! different boundary conditions are applied during inspiration and 
          ! expiration: Dirichlet at model entry during inspiration (concentration
          ! = inlet_concentration), and Neumann at model entry during expiration (dcdt = 0)
          if(.not.inspiration)then
             BB = global_R(nrow)         !get reduced R.H.S.vector
             do nentry = sparsity_row(nrow),sparsity_row(nrow+1)-1  !each row entry
                ncol = sparsity_col(nentry)
                BB = BB-global_K(nentry)*node_field(nj_conc1,ncol)*dt ! -K*c^(n)*dt
                AA = global_M(nentry) + dt*theta*global_K(nentry) ! M+K*dt*theta
                global_AA(nentry) = global_AA(nentry) + AA
             enddo
             global_BB(nrow) =  global_BB(nrow) + BB
          elseif(inspiration.and.np.ne.1)then !not first row
             BB = global_R(nrow)         !get reduced R.H.S.vector
             do nentry = sparsity_row(nrow),sparsity_row(nrow+1)-1  !each row entry
                ncol = sparsity_col(nentry)
                BB = BB-global_K(nentry)*node_field(nj_conc1,ncol)*dt
                AA = global_M(nentry) + dt*theta*global_K(nentry) !M+K*dt*theta
                if(ncol.ne.1)then ! not first column
                   global_AA(nentry-noffset_entry) = &
                        global_AA(nentry-noffset_entry) + AA
                endif
             enddo
             global_BB(nrow-noffset_row) = &
                  global_BB(nrow-noffset_row) + BB
          endif
       enddo

       if(diffusion_on)then
          solution(1:num_nodes) = zero_tol
          ! Call a solver to solve the system of reduced equations. 
          ! Here we use an iterative solver (GMRES == Generalised Minimal 
          ! RESidual method). The solver requires the solution matrices to 
          ! be represented in compressed row format.
          call pmgmres_ilu_cr (MatrixSize,NonZeros,reduced_row,&
               reduced_col,global_AA,solution,global_BB,&
               part_param%out_itr_max,part_param%inr_itr_max,&
               part_param%solve_tolerance,part_param%solve_tolerance,SOLVER_FLAG)
          
          ! transfer the solver solution (in 'Solution') to the node field array          
          do np = 1,num_nodes
             if(.not.inspiration.or.(inspiration.and.np.gt.1))then
                if(inspiration)then
                   nrow_BB=np-1
                else
                   nrow_BB=np
                endif
                node_field(nj_conc1,np) = node_field(nj_conc1,np) &
                     + Solution(nrow_BB) !c^(n+1)=c^(n)+dc
                node_field(nj_conc1,np) = max(0.0_dp,node_field(nj_conc1,np))
                node_field(nj_conc1,np) = min(1.0_dp,node_field(nj_conc1,np))
             endif
          enddo
       endif

       call volume_of_mesh(current_volume,volume_tree)

       if(deposition_on)then
          if(inlet_flow.lt.0.0_dp)then
              call particle_deposition(dt,.false.,part_param, num_iter, DepFrac)
             !node_field(nj_conc1,1) = tp%inlet_concentration(1)
          else
              call particle_deposition(dt,.true.,part_param, num_iter, DepFrac)
              node_field(nj_conc1,1) = tp%inlet_concentration(1)
          endif
       endif
       ! estimate the volume and mass errors
       call calc_mass_particles(nj_conc1,nu_conc1,lumen_mass,mass_deposit,mass_by_gen, &
            part_param,unit_mass,unit_wall)

!!! TJ - add extra-thoracic
    ! Calculate extra_mass_inspiration when necessary
       if (inspiration) then
            ! Accumulate extra_mass during both inhalation and exhalation
            extra_mass_inspiration = extra_mass + last_extra_mass_inspiration
            ! Calculate extra_mass_inspiration for the current inhalation phase
!            print *, "extra_mass", extra_mass, "extra_mass_inspiration", extra_mass_inspiration,&
!                    "last_extra_mass_inspiration", last_extra_mass_inspiration
       else
           extra_mass_inspiration = extra_mass + last_extra_mass_inspiration
           !!! TJ - this is assuming no extra-thoracic deposition for exhale
!            ! For exhalation, maintain the same extra_mass from the previous inhalation
!            extra_mass = last_extra_mass_inspiration
!            ! Update extra_mass_inspiration with the value from the previous inhalation
!            extra_mass_inspiration = last_extra_mass_inspiration

!            print *, "extra_mass", extra_mass, "extra_mass_inspiration", extra_mass_inspiration,&
!                    "last_extra_mass_inspiration", last_extra_mass_inspiration
       endif

       current_mass = lumen_mass + mass_deposit +  unit_mass + unit_wall + extra_mass_inspiration
       
       if(.not.deposition_on) mass_deposit = 0.0_dp

       volume_error = 1.0e+2_dp*(current_volume - (part_param%initial_volume +&
            tp%total_volume_change))/(part_param%initial_volume + tp%total_volume_change)

       if(tp%ideal_mass.gt.0.0_dp)then
          mass_error = 1.0e+2_dp*(current_mass - extra_mass_inspiration - tp%ideal_mass)/tp%ideal_mass
       else
          mass_error = 0.0_dp
       endif
       call cpu_time(cpu_dt_end)

       ! TJ - 25 NOV 2022: here we re-define DE and DF.
       dep_eff_alv = unit_wall/current_mass
       dep_eff_bronch = mass_deposit/current_mass
       dep_eff_extra = extra_mass_inspiration/current_mass

       if(write_mass)then
          write(*,'(f7.3,3(f8.3),f10.2,7(f9.2),f10.5,'' |'',16(f6.1))') &
               time,inlet_flow/1.0e+6_dp,tp%total_volume_change/1.0e+6_dp,&
               current_volume/1.0e+6_dp,tp%ideal_mass,current_mass,lumen_mass, &
               mass_deposit,unit_mass, unit_wall, volume_error, &
               mass_error,node_field(nj_conc1,1),(mass_by_gen(i),i=1,16)
       else
          write(*,'(f7.3,2(f8.3),'' ('',f5.2,'') |'',2(f8.2),'' ('',f5.2,'') |'',10(f8.2),'' |'',3(f8.3), f10.5)')&
               time,tp%total_volume_change/1.0e+6_dp,current_volume/1.0e+6_dp, volume_error, &
               tp%ideal_mass, &
               (current_mass-extra_mass_inspiration) , &  ! the intra-thoracic mass in the entire model
               mass_error, &    ! error in total mass
               extra_mass_inspiration,&
               lumen_mass, &    ! the mass in the bronchial airway lumen
               mass_deposit, &  ! the total mass deposited on bronchial airway walls
               sum(node_field(nj_loss_dif,1:num_nodes)), &  ! bronchial dep mass by diffusion
               sum(node_field(nj_loss_sed,1:num_nodes)), &  ! bronchial dep mass by sedimentation
               sum(node_field(nj_loss_imp,1:num_nodes)), &  ! bronchial dep mass by impaction
               unit_mass, &     ! the mass in the alveolar airway lumen
               unit_wall, &     ! the total mass deposited on alveolar walls
               sum(unit_field(nu_loss_dif,1:num_units)), &  ! alveolar dep mass by diffusion
               sum(unit_field(nu_loss_sed,1:num_units)), &  ! alveolar dep mass by sedimentation
               ! TJ - NEEDS TO CHANGE ABOVE
               (dep_eff_alv+dep_eff_bronch), &  ! total deposition fraction in lung
               dep_eff_bronch, &  ! airway deposition fraction
               dep_eff_alv, &       ! acinar deposition fraction
               node_field(nj_conc1,1)

!               (1-part_param%Eo)*(dep_eff_alv+dep_eff_bronch)+ part_param%Eo, &  ! total deposition fraction
!               (1-part_param%Eo)*dep_eff_bronch, &  ! airway deposition fraction
!               (1-part_param%Eo)*dep_eff_alv        ! acinar deposition fraction
       endif

       err = mass_error
       dep_eff_total = dep_eff_alv+dep_eff_bronch
!      if(last_breath.and.coupled)then
       if(last_breath)then
          if(elem_field(ne_Vdot,1).lt.0.0_dp) then
             carryon = .true.
          else
             carryon = .false.
          endif
       else
          if(time_end - time .le. zero_tol)then
             carryon = .false.
          else
             carryon = .true.
          endif
       endif
    enddo !carryon

    write(*,'(''---------------------------------------'')')
    write(*,'('' End of breath deposition efficiencies:'')')
    write(*,'(5x, ''TOTAL DE ='', f6.3, ''   BRONCHIAL DE ='',f6.3, ''   ALVEOLAR DE='',f6.3, ''    &
            EXTRATHORACIC DE ='',f6.3)' ) &
            (dep_eff_alv+dep_eff_bronch+dep_eff_extra), dep_eff_bronch, dep_eff_alv, dep_eff_extra
    write(*,'(5x, ''Diffusive-    '',5x, f8.3, 15x, f8.3)') &
         sum(node_field(nj_loss_dif,1:num_nodes)),sum(unit_field(nu_loss_dif,1:num_units))
    write(*,'(5x, ''Sedimentation-'',5x, f8.3, 15x, f8.3)') &
         sum(node_field(nj_loss_sed,1:num_nodes)), sum(unit_field(nu_loss_sed,1:num_units))
    write(*,'(5x, ''Impaction-    '',5x, f8.3, 15x, ''   -'')') &
         sum(node_field(nj_loss_imp,1:num_nodes))
    write(*,'(5x, ''Totals--      '',5x, f8.3, 15x, f8.3)') dep_eff_bronch, dep_eff_alv
    write(*,'(''---------------------------------------'')')

    if(write_airway_on) call write_airway(ne_field, 'airway_result', groupname, field_name)
    if(write_terminal_on) call write_terminal('terminal', name)

    call enter_exit(sub_name,2)
  end subroutine solve_particles


!!!###################################################################################
    subroutine write_airway(ne_field, filename, groupname, field_name)

        use arrays
        use exports,only: export_node_field
        use geometry, only: volume_of_mesh
        use species_transport, only: assemble_transport_matrix,reduce_transport_matrix,calc_mass
        use other_consts
        use diagnostics, only: enter_exit

        implicit none

        type(particle_parameters) :: part_param
        type(transport_parameters) :: tp
            integer :: len_end,ne
            integer :: np0, np1, np2, ne_field, nj, ne0, ne_stem
            real(dp) :: midpoint(3)
            logical :: CHANGED
            character(len=300) :: writefile
            character :: sub_name, groupname, field_name
        character :: filename*(*)
        character(LEN=10) :: lobe
        integer :: ifile=10
        character :: readfile*150
        groupname = 'vent_model'
        field_name = 'bronchial_deposition'
        ne_field = 13
        np0 = 1

        if(index(filename, ".txt")>0)then ! full filename is given
           readfile = trim(filename)
        else! append correct extension
           readfile = trim(filename)//'.txt'
        endif
        open(ifile, file=readfile, status='replace')

        do ne = 1,num_elems
          np1 = elem_nodes(1,ne)
          np2 = elem_nodes(2,ne)

          midpoint(:) = (node_xyz(:, np1) + node_xyz(:, np2))/2
          ne0 = element_lobe(ne)

          select case(ne0)
            case(166, 177)
                lobe = 'LUL'
            case(199, 200)
                lobe = 'LLL'
            case(201, 202)
                lobe = 'RUL'
            case(203)
                lobe = 'RLL'
            case(204)
                lobe = 'RML'
            case default
                lobe = 'Trachea' !'copious free time'
          end select

          write(10,'(I6, 4(F10.2),  2(I5), 2x, A5, 2(F7.3), 7(D11.3))') &
                  ne, (node_xyz(:, np1) + node_xyz(:, np2))/2, &
                  sqrt(sum((midpoint(:) - node_xyz(:, np0))**2)), & ! distance
                  elem_ordrs(1,ne), elem_ordrs(2,ne), lobe, elem_field(ne_length,ne), &
                  elem_field(ne_radius,ne), elem_field(ne_Vdot,ne), elem_field(ne_mass,ne), & !change ne_flow to ne_vdot
                  node_field(nj_loss_dif, ne), node_field(nj_loss_imp, ne), &
                  node_field(nj_loss_sed, ne), & ! alveolar dep mass by sed
                  node_field(nj_conc1, ne) , elem_field(ne_part_vel,ne) ! concentration, vp
        end do
        close(ifile)
    end subroutine write_airway

!!!#########################################################################

    subroutine write_terminal(filename, name)

        use arrays
        use exports,only: export_node_field
        use geometry, only: volume_of_mesh,group_elem_parent_term
        use species_transport, only: assemble_transport_matrix,reduce_transport_matrix,calc_mass
        use other_consts
        use diagnostics, only: enter_exit

        implicit none

        type(particle_parameters) :: part_param
        type(transport_parameters) :: tp

        integer :: ifile=10
        character :: readfile*150
        character :: filename*(*), name*(*)
        character(LEN=5) :: lobe
        !!! Local Variables
        integer :: len_end,ne,nj,NOLIST,np,np_last,VALUE_INDEX
        integer :: np0, ne0 !_stem, num_list_total, elem_list_total ! initial node
        logical :: FIRST_NODE

        name = 'terminal'
        np0 = 1
        !ne_parent = 211!(/211, 214, 215, 223, 224/)

        if(index(filename, ".txt")>0)then ! full filename is given
           readfile = trim(filename)
        else! append correct extension
           readfile = trim(filename)//'.txt'
        endif
        open(ifile, file=readfile, status='replace')

        do nolist = 1,num_units
          ne = units(nolist)
          np = elem_nodes(2,ne)

          ne0 = element_lobe(ne)

          select case(ne0)
            case(166, 177)
                lobe = 'LUL'
            case(199, 200)
                lobe = 'LLL'
            case(201, 202)
                lobe = 'RUL'
            case(203)
                lobe = 'RLL'
            case default
                lobe = 'RML' !'copious free time'
          end select

          ! TJ - find the concentration array!

          write(10,'(I6, 3(f10.2), f10.2, A7, f8.3, f10.2, 4(D11.3))') &
                  np, (node_xyz(:,np)), sqrt(sum((node_xyz(:, np) - node_xyz(:, np0))**2)), &  ! distance
                  lobe, unit_field(nu_vol,nolist), unit_field(nu_vdot0,nolist), &
                  unit_field(nu_vol,nolist)*unit_field(nu_conc1,nolist), &
                  unit_field(nu_loss_dif, nolist), unit_field(nu_loss_sed, nolist), elem_field(ne_part_vel,ne)

          !write(10,'(1X,''Lobe: '', A10)') lobe !element_lobe(ne)
        end do!nolist
        close(ifile)
    end subroutine write_terminal

!!!###################################################################################

  subroutine scale_flow_field(inlet_flow)

!   use arrays,only: dp,elem_field,num_elems,num_units,unit_field

    use arrays
    use other_consts
    use indices

!   use indices,only: ne_dvdt,nu_Vdot0

    real(dp),intent(in) :: inlet_flow
    real(dp) :: ratio

    if(abs(elem_field(ne_Vdot,1)).gt.zero_tol)then
       ratio = inlet_flow/elem_field(ne_Vdot,1)

       unit_field(nu_Vdot0,1:num_units) = unit_field(nu_Vdot0,1:num_units)*ratio
       elem_field(ne_Vdot,1:num_elems) = elem_field(ne_Vdot,1:num_elems)*ratio
       elem_field(ne_dvdt,1:num_elems) = elem_field(ne_dvdt,1:num_elems)*ratio

    else
       write(*,'('' Cannot scale to zero flow'')')
    endif

  end subroutine scale_flow_field

!###################################################################################


  subroutine general_track(dt,update,tp,part_param)  
    ! the argument 'update' is a tag to decide whether the concentrations have to updated or not
    ! this should be false in case of mass calculation which is used in the "update_unit_mass" subroutine

!   use arrays,only: dp,node_field,elem_field,num_elems,num_units,unit_field, elem_nodes, elem_cnct,&
!                    units, elems_at_node
    use arrays

!   use indices, only: nj_conc1, nj_conc2, nj_conc3
    use indices

    use geometry
    use species_transport, only: calc_mass
    use other_consts
    use diagnostics, only: enter_exit
    
    implicit none
  
    type(particle_parameters) :: part_param
    type(transport_parameters) :: tp
    real(dp), intent(in) :: dt
    logical, intent(in) :: update

    integer,parameter :: n_gases = 3
    integer :: c_kount,i,i_elem,j,j1,jj2,kount,ne,ne_stem,nj_g(3),np,np1,np2,np1_parent,np2_parent, &
         num_list_total,nunit,nu_g(3),o_elem,sum_inout,sum_dir
    integer,parameter :: elem_f_to_node = 1, f_sign = 2
    integer,allocatable :: TMAT(:,:,:),elem_list_total(:)
    real(dp) :: c_mass(3),concens(3),concentration(3),flow_fraction, fl_temp, &
         inlet_conc(3),local_xi,mass_fraction,mean_c(3),parent_fraction, &
         time_through_element,total_time, unit_mass(3), &
         ratio_mass,total_mass, max_concentration
!!! currently unused but might be again
!!!  ideal_mass_total,mass_at_max,mass_below_max, total_mass,ratio_mass,
    real(dp), allocatable :: concs_at_node(:,:), concent(:,:)
    logical :: go_on, cont, neg_fl
    character(len=60) :: sub_name

    sub_name = 'general_track'
    call enter_exit(sub_name,1)

    nj_g(1) = nj_conc1
    nj_g(2) = nj_conc2
    nj_g(3) = nj_conc3
    nu_g(1) = nu_conc1
    nu_g(2) = nu_conc2
    nu_g(3) = nu_conc3
    inlet_conc(:) = tp%inlet_concentration(:)

    if (.not.allocated(TMAT)) allocate(TMAT(num_nodes,2,3))
    
    !second index: entering/exiting (1/2) flow towards/from a node 
    ! and positive/negative (1/2) flow at each element (positive: downward,
    ! negative:upward), Third index: the elements at the node
    
!!! If the flow is negative in an element, the negative flow always come from its children
!!! however, if the flow is positive the flow might be coming from its parents or its parents' other children
    
    call track_mat_coupled(elem_f_to_node,f_sign,TMAT) !calling the matrix containing flow directions and nodes in/out data
!!! TMAT(np,2,ne) = 1 when flow in element ne is distal (to small airways)
!!! TMAT(np,2,ne) = 2 when flow in element ne is proximal (to large airways)
!!! TMAT(np,1,ne) = 1 when element ne brings flow to node np
!!! TMAT(np,1,ne) = 2 when element ne takes flow from node np
    
    allocate(concent(num_nodes,n_gases))
    
    concent = 0.0_dp
    do nunit = 1,num_units
       ne = units(nunit)
       if(elem_cnct(1,0,ne).eq.0)then ! a terminal
          np = elem_nodes(2,ne)
          if(elem_field(ne_Vdot,ne).lt.0.0_dp)then
             do i = 1,n_gases
                unit_mass(i) = unit_field(nu_vol,nunit)*unit_field(nu_g(i),nunit)
                c_mass(i) = elem_field(ne_Vdot,ne)*node_field(nj_g(i),np)*dt
                unit_mass(i) = unit_mass(i) + c_mass(i)
                !             unit_field(nu_g(i),nunit) = unit_mass(i)/unit_field(nu_vol,nunit)
                node_field(nj_g(i),np) = unit_mass(i)/unit_field(nu_vol,nunit)
                forall (i=1:n_gases) concent(np,i) = node_field(nj_g(i),np)
             enddo
          endif
       endif ! terminal
    enddo
    
    do np = 1,num_nodes
       allocate(concs_at_node(elems_at_node(np,0),n_gases))
       concs_at_node = 0.0_dp
!!!    check that the flows at the node make sense. Can't have flow in 
!!!    all elements directed towards the node, and can't have all flows directed out from the node.
       sum_inout = 0
       sum_dir = 0
       do J = 1,elems_at_node(np,0)
          sum_inout = sum_inout + TMAT(np,elem_f_to_node,J) !sum of flow directions: proximal (2) and distal (1) 
          sum_dir = sum_dir + TMAT(np,f_sign,J)     !sum of flow directions w.r.t. node np
       enddo
       go_on = .true.
       if ((sum_inout == 6) .and. (sum_dir == 4)) go_on = .false. !impossible
       if ((sum_inout == 3) .and. (sum_dir == 5)) go_on = .false. !trapping
       if ((sum_inout == 4) .and. (sum_dir == 3)) go_on = .false. !impossible
       if ((sum_inout == 2) .and. (sum_dir == 3)) go_on = .false. !trapping
       if (go_on) then !meets criteria for sensible flow
          concs_at_node = 0.0_dp
          do J = 1,elems_at_node(np,0)
             ne = elems_at_node(np,J)
             if (TMAT(np,1,J) == 1) then !flow entering the node so contributing to the concentration
                if (TMAT(np,2,J) == 1) then !positive flow
                   c_mass = 0.0_dp  !in case of any appended units!
                   time_through_element = abs(elem_field(ne_vol,ne)/elem_field(ne_Vdot,ne))
                   total_time = time_through_element
                   np1 = elem_nodes(1,ne)
                   np2 = elem_nodes(2,ne)
                   if (total_time.ge.dt) then !location is within this element
                      local_xi = dt/time_through_element
                      do i = 1,n_gases
                         concentration(i) = ((1.0_dp-local_xi)*node_field(nj_g(i),np2)) + &
                              (local_xi*node_field(nj_g(i),np1))
                         mean_c(i) = 0.5_dp*(concentration(i) + node_field(nj_g(i),np2))
                         c_mass(i) = c_mass(i) + mean_c(i) * local_xi &                     ! *MHT change*
                              *abs(elem_field(ne_vol,elems_at_node(np2,1))) !!! for acinus
                      enddo
                   else
                      do i = 1,n_gases
                         mean_c(i) = 0.5_dp*(node_field(nj_g(i),np1)+node_field(nj_g(i),np2))
                         c_mass(i) = c_mass(i) + mean_c(i)*abs(elem_field(ne_vol,elems_at_node(np2,1)))
                      enddo
                      cont = .true.
                      JJ2 = elems_at_node(np,J) ! same value as ne?
                      concentration = 0.0_dp
                      parent_fraction = 1.0_dp
                      J1 = JJ2
                      do while (cont)
                         JJ2 = elem_cnct(-1,1,JJ2) !parent element
                         if (JJ2.gt.0) then !parent is not the trachea
                            np1_parent = elem_nodes(1,JJ2)
                            np2_parent = elem_nodes(2,JJ2)
                            ! Because different parents and children are supplying the node
                            ! they have their own flow fractions. Parents are getting flow from
                            ! their parents for which they might get flow from their children.
                            ! For the children fractions are coming from their nodes but for the
                            ! parents the fraction is from their parents and their parent's children
                            ! going backward to the trachea.
                            ! summing up concentrations for child suppliers. 
                            if (elems_at_node(np2_parent,0) == 3) then
                               ! check whether all of the flow comes from the parent, or some from a sibling
                               neg_fl = .false.
                               do kount = 2,3
                                  ! one of these is the node itself which has positive flow so
                                  ! automatically the other one is selected
                                  if (TMAT(np2_parent,2,kount) == 2) then ! negative flow in sibling
                                     i_elem = elems_at_node(np2_parent,kount) ! sibling element number
                                     neg_fl = .true.
                                     do i = 1,n_gases
                                        call negative_flow(dt,nj_g(i),np2_parent,kount,TMAT,concens(i))
                                     enddo
                                     fl_temp = abs(elem_field(ne_Vdot,i_elem))+ &
                                          abs(elem_field(ne_Vdot,JJ2)) ! sum parent and sibling flows
!                                  else
                                     o_elem = elems_at_node(np2_parent,kount)
                                  endif !TMAT
                               enddo !kount
                               if (neg_fl) then
                                  if (TMAT(np2_parent,2,1) == 2) then ! negative in parent; all comes from sibling
                                     flow_fraction = 1.0_dp ! MHT: haven't checked that this part is correct
                                  else
!                                     fl_temp = abs(elem_field(ne_Vdot,o_elem)+ &
!                                          elem_field(ne_Vdot,JJ2))
                                     flow_fraction = abs(elem_field(ne_Vdot,o_elem)/fl_temp)
                                  endif
                                  !calculate the total fraction up to the supplying parent
                                  parent_fraction = parent_fraction * abs(1.0_dp-flow_fraction)
                                  do i = 1,n_gases
                                     concentration(i) = concentration(i) + (flow_fraction * concens(i))
                                     mean_c(i) = 0.5d0*(concens(i) + node_field(nj_g(i),np2_parent))
                                     c_mass(i) = c_mass(i) + abs(mean_c(i) * elem_field(ne_vol,i_elem)* &
                                          flow_fraction)
                                  enddo
                                  ! the mass the negative flow chidren are adding are just depending
                                  ! on their own flow fraction (not others)
                               endif
                            endif !elems_at_...
                            
                            ! adding the concentration for the tracked parent supplier with respect
                            ! to its fraction down to the node to the concentration of the child
                            ! suppliers calculated above
                            mass_Fraction = elem_field(ne_Vdot,J1)/elem_field(ne_Vdot,JJ2)
                            !if (JJ2.gt.0) then !parent is not the trachea
                            time_through_element = abs(elem_field(ne_vol,JJ2)/&
                                 (elem_field(ne_Vdot,JJ2)))
                            total_time = total_time + time_through_element
                            if (TMAT(np2_parent,2,1) == 2) then
                               forall (i = 1:n_gases) concentration(i) = concentration(i) + (parent_fraction * &
                                    node_field(nj_g(i),np2_parent)) 
                               cont = .false.
                            else if (total_time.ge.dt) then ! Material is within this element
                               local_xi = (total_time-dt)/time_through_element
                               do i = 1,n_gases
                                  concentration(i) = concentration(i) + (parent_fraction * &
                                       (((1.d0-local_xi)*node_field(nj_g(i),np1_parent)) + &
                                       (local_xi*node_field(nj_g(i),np2_parent))))
                                  mean_c(i) = 0.5_dp*(concentration(i) + node_field(nj_g(i),np2_parent))
                                  c_mass(i) = c_mass(i) + mean_c(i) * abs(elem_field(ne_vol,JJ2) * &
                                       (1.d0-local_xi) * mass_fraction) !!!acinus
                               enddo
                               cont = .false.
                            else
                               do i = 1,n_gases
                                  mean_c(i) = 0.5_dp*(node_field(nj_g(i),np1_parent) + node_field(nj_g(i),np2_parent))
                                  c_mass(i) = c_mass(i) + mean_c(i) * abs(elem_field(ne_vol,JJ2) * mass_fraction)
                               enddo
                            endif !total_time.ge.dt
                         else
                            mass_fraction = elem_field(ne_Vdot,J1)/elem_field(ne_Vdot,1)
                            forall (i = 1:n_gases) concentration(i) = concentration(i) + (parent_fraction * inlet_conc(i))
!!! if the fluid is being advected from outside the model, then the 'interpolated mass' should be such that the
!!! concentration is the same as the inspired concentration (inlet concentration), i.e. M = effective_volume*c_inlet
                            !                            interp_mass = interp_mass + (parent_fraction * node_field(nj_mass,1)*2.0_dp)
                            forall (i = 1:n_gases) c_mass(i) = c_mass(i) + concentration(i)*(dt-total_time)*mass_fraction
                            cont = .false.
                         endif !JJ2.gt.0
                      enddo !cont
                   endif !total_time.ge.dt
                   forall (i = 1:n_gases) concs_at_node(J,i) = concentration(i)
                endif !TMAT(np,2,J) == 1
                if (TMAT(np,2,J) == 2) then !negative flow
                   do i = 1,n_gases
                      call negative_flow(dt,nj_g(i),np,J,TMAT,concens(i))
                      concs_at_node(J,i) = concens(i)
                      ! estimate the expired mass:
                      c_mass(i) = 0.5_dp*(node_field(nj_g(i),np)+concs_at_node(J,i))* &
                           elem_field(ne_Vdot,elems_at_node(np,1))*dt
                   enddo
                endif !TMAT(np,2,J) == 2
             endif !if (TMAT(np,1,J) == 1)
          enddo !J=1...
       else
          write(*,'('' WARNING: ignored node:'',i8)') np
          if ((sum_inout == 6) .and. (sum_dir == 4)) write(*,'('' Impossible'')')
          if ((sum_inout == 3) .and. (sum_dir == 5)) write(*,'('' Trapping'')')
          if ((sum_inout == 4) .and. (sum_dir == 3)) write(*,'('' Impossible'')')
          if ((sum_inout == 2) .and. (sum_dir == 3)) write(*,'('' Trapping'')')
          do J = 1,elems_at_node(np,0)
             ne = elems_at_node(np,j)
             write(*,'('' Flow from element'',i8,'' ='',d12.4,'' mm^3/s'')') ne,elem_field(ne_Vdot,ne)
          enddo
          read(*,*)
       endif !if (go_on)
!!! update the temporary array for nodal concentrations (concent)
       c_kount = 0
       fl_temp = 0.0_dp
       do J = 1,elems_at_node(np,0)
          ne = elems_at_node(np,J)
          if (TMAT(np,1,J) == 1) then
             fl_temp = fl_temp + abs(elem_field(ne_Vdot,ne))!-elem_field(ne_dvdt,elems_at_node(np,J)))
             c_kount = c_kount + 1  !check the number of incoming flows to the node
             forall (i = 1:n_gases) concent(np,i) = concs_at_node(J,i) !if there's only one supplier concentration remains the same
           else
             o_elem = elems_at_node(np,J) !the outgoing element
          endif !conc_at_...
       enddo !J=1...
       if (c_kount == 2) then !two elements supplying the node so check the flow fractions
         concent(np,1:n_gases) = 0.0_dp
         do J = 1,elems_at_node(np,0)
           if (elems_at_node(np,J).ne.o_elem) then
             flow_fraction = abs((elem_field(ne_Vdot,elems_at_node(np,J)))/fl_temp)
             forall(i = 1:n_gases) concent(np,i)=concent(np,i)+(flow_fraction*concs_at_node(J,i))
           endif !elems_at...
         enddo !J
       endif ! c_kount

       do i = 1,n_gases
          if (concent(np,i) .le. 1.0e-10_dp) concent(np,i) = 0.0_dp
       enddo
       
       if ((np.eq.1).and.(TMAT(np,2,1) == 1))then
          ! an entrance node inspiring
          forall (i = 1:n_gases) concent(np,i) = inlet_conc(i)
       endif
       if ((np/=1).and.(elems_at_node(np,0)==1).and.(TMAT(np,1,1) == 2)) then
          ! a terminal node expiring
          ne = elems_at_node(np,1)
          if(inlist(ne,units))then ! has a lumped unit attached
             nunit = where_inlist(ne,units)
!             forall (i=1:2) node_field(nj_g(i),np) = unit_field(nu_g(i),nunit)
             forall (i = 1:n_gases) c_mass(i) = elem_field(ne_Vdot,ne)*node_field(nj_g(i),np)*dt
          endif
          forall (i = 1:n_gases) concent(np,i) = node_field(nj_g(i),np)
       endif
       
       deallocate(concs_at_node)

!!! add the accumulated mass to lumped units
       ne = elems_at_node(np,1) ! the first element that the node is in
       if ((np/=1).and.(elems_at_node(np,0) == 1))then ! a terminal node
          if(inlist(ne,units))then ! is terminal with a lumped unit attached
             nunit = where_inlist(ne,units)
             do i = 1,n_gases
                unit_mass(i) = unit_field(nu_vol,nunit)*unit_field(nu_g(i),nunit)
                unit_mass(i) = unit_mass(i) + c_mass(i)
                unit_field(nu_g(i),nunit) = unit_mass(i)/unit_field(nu_vol,nunit)
             enddo
          endif
       else ! is not terminal
          if(inlist(ne,units))then ! this element supplies a branching acinus within a unit
             ! record the cumulative mass entering the acinus + initial mass (dummy storage)
             elem_field(ne_resist,ne) = elem_field(ne_mass,ne) + c_mass(1)
          endif
       endif

    enddo !np

    max_concentration = max(tp%initial_concentration(1),tp%inlet_concentration(1)) ! for primary gas
    if (update) then
       forall(i = 1:n_gases) node_field(nj_g(i),1:num_nodes)= concent(1:num_nodes,i)

!!! adjust the concentrations to conserve mass. Only done for nodes where the concentration
!!! is less than the maximum (set in 'initial_gasmix')
!   ideal_mass_total = elem_field(ne_mass,1) + elem_field(ne_Vdot,1)*dt*node_field(nj_conc1,1)
       call calc_mass(nj_conc1,nu_conc1,total_mass)
       allocate(elem_list_total(num_elems))
       do nunit = 1,num_units
          ne_stem = units(nunit)
          num_list_total = 0
          if(elem_cnct(1,0,ne_stem).ne.0)then
             elem_list_total = 0
!!! get a list of all elements that are within an elastic unit
             call group_elem_by_parent(ne_stem,elem_list_total)
             num_list_total = count(elem_list_total.ne.0)
          endif
          !! the 'ideal' mass for the unit is elem_field(ne_resist,ne_stem)
          !! the 'current' mass for the unit is elem_field(ne_mass,ne_stem)
          if(elem_field(ne_mass,ne_stem).gt.zero_tol)then
             ratio_mass = elem_field(ne_resist,ne_stem)/elem_field(ne_mass,ne_stem)
          else
             ratio_mass = 1.0_dp
          endif
          do i = 1,num_list_total
            ne = elem_list_total(i)
            np2 = elem_nodes(2,ne)
            if(node_field(nj_conc1,np2).lt.max_concentration) then
                node_field(nj_conc1,np2) = node_field(nj_conc1,np2)*ratio_mass
            endif
          enddo
       enddo
       deallocate(elem_list_total)
    endif

    deallocate(concent)
    deallocate(TMAT)

    call enter_exit(sub_name,2)
    
  end subroutine general_track
  

!!!############################################################################################


  subroutine track_mat_coupled(elem_f_to_node,f_sign,TMAT)

!   use arrays, only: dp, elems_at_node, elem_field
    use arrays
    use diagnostics, only: enter_exit

    implicit none
  
    integer,intent(in) :: elem_f_to_node,f_sign
    integer :: TMAT(:,:,:)
    integer :: I, J
    real(dp) :: flo
    character(len=60) :: sub_name
    
    sub_name = 'track_mat_coupled'
    call enter_exit(sub_name,1)
    
    ! put these two below also in the array.f90 as global parameters
    TMAT = 0
    do I = 1,num_nodes
       do J = 1,elems_at_node(I,0)
          flo = elem_field(ne_Vdot,elems_at_node(I,J))
          if (flo .lt. 0.0_dp) then
             if (I.eq.1) then
                TMAT(I,elem_f_to_node,J) = 1
                TMAT(I,f_sign,J) = 2
             elseif (J.eq.1) then !parent element
                TMAT(I,elem_f_to_node,J) = 2
                TMAT(I,f_sign,J) = 2
             else !child elements
                TMAT(I,elem_f_to_node,J) = 1
                TMAT(I,f_sign,J) = 2
             endif
          else
             if (I.eq.1) then
                TMAT(I,elem_f_to_node,J) = 2
                TMAT(I,f_sign,J) = 1
             elseif (J.eq.1) then !parent element
                TMAT(I,elem_f_to_node,J) = 1
                TMAT(I,f_sign,J) = 1
             else !child elements
                TMAT(I,elem_f_to_node,J) = 2
                TMAT(I,f_sign,J) = 1
             endif
          endif
       enddo
    enddo
    
    call enter_exit(sub_name,2)

  end subroutine track_mat_coupled

!!!############################################################################################

  subroutine negative_flow(dt,nj_g,I,J,TMAT,concs)

!   use arrays,only: dp, elem_cnct, elem_symmetry, elems_at_node, node_field, elem_nodes
    use arrays
    use diagnostics, only: enter_exit

    implicit none
  
    integer, intent(in) :: I, J, nj_g,TMAT(:,:,:)
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: concs
    
    integer :: I1, I2, JJ3, kount, child, c_kount, IK
    integer, allocatable :: parents(:), children(:)
    real(dp) :: time_through_element, total_time, concentration, local_xi, fl_temp 
    real(dp) :: flow_fraction, flow_pa
    real(dp), allocatable :: parent_times(:), child_times(:)
    logical :: cont, in_cont
    character(len=60) :: sub_name
    
    sub_name = 'negative_flow'
    call enter_exit(sub_name,1)
    
    ! Ok now we are creating a matrix containing the final nodes from which
    ! the materials are reaching the current node. Then we will calculate
    ! the summation of the concentrations coming from these finals nodes by
    ! going back from them towards the current node and taking into account
    ! the flow fractions where two adjacent elements are feeding their parent
    ! element with negative flow
    
    allocate(parents(1))
    allocate(parent_times(1))
    allocate(children(128))
    allocate(child_times(128))
    
    cont = .true.
    in_cont = .true.
    parents(1) = elem_nodes(2,elems_at_node(I,J)) !second node of the child element
    parent_times(1) = abs(elem_field(ne_vol,elems_at_node(I,J))/ &
         (elem_field(ne_Vdot,elems_at_node(I,J))))
    if (parent_times(1).gt.dt) cont = .false.
    do while (in_cont)
       in_cont = .false.
       c_kount = 0
       do kount = 1,size(parents)
          cont = .false.
          do JJ3 = 2,elems_at_node(parents(kount),0) !find the children for parents (parents are nodes here)
             child = elems_at_node(parents(kount),JJ3) ! Child element for that node
             time_through_element = abs(elem_field(ne_vol,child)/&
                  (elem_field(ne_Vdot,child)))
             total_time = time_through_element + parent_times(kount) !calculate time in each child to see if it's less than dt so the material passes it
             if ((elems_at_node(parents(kount),0).gt.1) .or. (I == 1)) then !parent node is not a terminal node or it is at trachea
                if (TMAT(parents(kount),2,JJ3) == 2) then !the child has negative flow
                   if (parent_times(kount) .lt. dt) then !material is not in this elemen
                      cont = .true. !if the parents are not terminals
                      ! ,their children's flow are negative and the material is not within them
                      ! then we have to go furthur down
                      in_cont = .true. ! still need to check elems thru which material passes; might be some -ve flow in other subtended branches
                      c_kount = c_kount + 1
                      children(c_kount) = elem_nodes(2,child)
                      child_times(c_kount) = total_time !saving the time for the element to be used in calculating the concentrations at its nodes later
                   endif ! time through
                endif !TMAT
             endif ! long if condition above (elem_at_ .....)
          enddo !JJ3
          if (.not.cont) then !either the child elements have positive flow or their parent is a terminal
             !the parent should be preserved
             c_kount = c_kount + 1
             if(c_kount.gt.128)then
                write(*,*) 'Increase size of children array in general_track, or decrease dt'
                read(*,*)
             endif
             children(c_kount) = parents(kount)
             child_times(c_kount) = parent_times(kount)
          endif ! not cont
       enddo !kount
       deallocate(parents)
       allocate(parents(c_kount))
       deallocate(parent_times)
       allocate(parent_times(c_kount))
       parents = children(1:c_kount)
       parent_times = child_times(1:c_kount)
    enddo !cont
    deallocate(children)
    deallocate(child_times)

    concs = 0.0_dp
    do kount = 1,size(parents)
       concentration = 0.0_dp
       I1 = elem_nodes(1,elems_at_node(parents(kount),1)) !first node of the element containing the node
       I2 = elem_nodes(2,elems_at_node(parents(kount),1)) !second node of the element containing the node (actually the node itself!!)
       time_through_element = abs(elem_field(ne_vol,elems_at_node(parents(kount),1))/&
            (elem_field(ne_Vdot,elems_at_node(parents(kount),1))))
       local_xi = abs((parent_times(kount)-dt))/time_through_element
       if ((parent_times(kount) .lt. dt) ) local_xi = 0.0_dp ! might only happen at the final elements(terminals)
       concentration =  ((1.0_dp-local_xi)*node_field(nj_g,I2)) + &
            (local_xi*(node_field(nj_g,I1)))
       JJ3 = elems_at_node(parents(kount),1) !element contatining the node
       do while (JJ3.gt.elems_at_node(I,J))
          c_kount = elems_at_node(I1,0)
          if (c_kount == 3) then ! the node has 3 elements
             if ((TMAT(I1,1,2) == 1) .and. (TMAT(I1,1,3) == 1)) then !both child elements are supplying the node
                
                fl_temp = 0.0_dp
                do IK = 2,c_kount
                   fl_temp = fl_temp + abs(elem_field(ne_Vdot,elems_at_node(I1,IK)))
                enddo
                flow_fraction = abs((elem_field(ne_Vdot,JJ3))/fl_temp)
                !if the other children is also supplying the same node it should be already in the list (parents) and taken care of for its flow fraction here so no need to calculate for both. 
                concentration = flow_fraction * concentration
             endif
          elseif (c_kount == 2) then
             flow_pa = elem_field(ne_Vdot,elem_cnct(-1,1,JJ3)) - elem_field(ne_dvdt,elem_cnct(-1,1,JJ3))
             flow_fraction = abs((elem_field(ne_Vdot,JJ3) * elem_symmetry(JJ3))/flow_pa)
             concentration = flow_fraction * concentration
          endif
          JJ3 = elem_cnct(-1,1,JJ3) !set the element to its parent
          I1 = elem_nodes(1,JJ3)
       enddo ! JJ3
       concs = concs + concentration
    enddo !kount
    
    deallocate(parents)
    deallocate(parent_times)

    call enter_exit(sub_name,2)
    
  end subroutine negative_flow
  
!!!############################################################################################

  subroutine track_mat(TMAT)

    use arrays
    use diagnostics, only: enter_exit

    implicit none
  
    integer :: TMAT(:,:,:)
    integer :: ne,noelem,np
    integer,parameter :: elem_f_to_node = 1, f_sign = 2
    real(dp) :: flo
    character(len=60) :: sub_name

    sub_name = 'track_mat'
    call enter_exit(sub_name,1)

    TMAT = 0
    do np = 1,num_nodes
       do noelem = 1,elems_at_node(np,0)
          ne = elems_at_node(np,noelem)
          flo = elem_field(ne_Vdot,ne)
          if (flo .lt. 0.0_dp) then
!!! flow is directed proximally
             TMAT(np,f_sign,noelem) = 2  !labels flow as in direction 2 (proximal)
             if (np.eq.1) then
                TMAT(np,elem_f_to_node,noelem) = 1 !brings flow to the node
             else if (noelem.eq.1) then !parent element
                TMAT(np,elem_f_to_node,noelem) = 2 !takes flow away from the node
             else !child elements
                TMAT(np,elem_f_to_node,noelem) = 1 !brings flow to the node
             endif
          else
!!! flow is directed distally
             TMAT(np,f_sign,noelem) = 1  !labels flow as in direction 1 (distal)
             if (np.eq.1) then
                TMAT(np,elem_f_to_node,noelem) = 2 !takes from away from the node
             else if (noelem.eq.1) then !parent element
                TMAT(np,elem_f_to_node,noelem) = 1 !brings flow to the node
             else !child elements
                TMAT(np,elem_f_to_node,noelem) = 2 !takes flow away from the node
             endif
          endif
       enddo
    enddo

    call enter_exit(sub_name,2)
    
  end subroutine track_mat

  
!!!############################################################################################

  function effective_volume_at_node(np)

    use arrays, only: dp, elems_at_node, elem_nodes, elem_field, elem_symmetry
    integer,intent(in) :: np

    integer :: i,ne,np1,np2
    real(dp) :: effective_volume_at_node

    effective_volume_at_node = 0.0_dp
    do i = 1,elems_at_node(np,0)
       ne = elems_at_node(np,i)
       np1 = elem_nodes(1,ne)
       np2 = elem_nodes(2,ne)
       effective_volume_at_node = effective_volume_at_node + 0.5_dp*elem_field(ne_vol,ne)*elem_symmetry(ne)
    enddo

  end function effective_volume_at_node

!!!################################################################################
  
  subroutine particle_velocity(time,part_param)
    ! Created by Falko Schmidt 05/2011
    ! calculates the particle velocity of spherical particles
    ! from the fluid flow field and the previous particle velocity due to drag
    
!   use arrays, only: dp, particle_parameters, elem_field,&
!                     node_xyz,elem_nodes,node_xyz,&
!                     num_elems,num_nodes
    use arrays
    use mesh_utilities,only: angle_btwn_vectors
    use diagnostics, only: enter_exit

    implicit none
  
    type(particle_parameters) :: part_param
    real(dp),intent(in) :: time
    integer :: j,ne,np,np2
    real(dp) :: ref_flow,rel_flow,alpha,abbr,vs, Ccun
    real(dp) :: vector1(3),vector2(3),vseff
    character(len=60) :: sub_name

    sub_name = 'particle_velocity'
    call enter_exit(sub_name,1)

!!! abbreviation for velocity equation
    Ccun = 1.0_dp + 2.0_dp * part_param%lambda/part_param%pdia &
         * (1.257_dp+0.4_dp*exp(-0.55_dp*part_param%pdia/part_param%lambda))

    abbr = -(18.0_dp*part_param%mu) * time /(part_param%prho*Ccun * part_param%pdia**2.0_dp)

!!! settlement velocity
    ! TJ - might check with conc?
    vs = (part_param%prho)*part_param%gravityz*part_param%pdia**2.0_dp&
                                             *Ccun/(18.0_dp*part_param%mu)
    
    do ne = 1,num_elems  !loop over all 1D elements
       np = elem_nodes(1,ne)  ! proximal node number
       np2 = elem_nodes(2,ne) ! distal node number
!!! calculate angle between branch axis and gravitational direction (alpha)
       forall(j=1:3) vector1(j) = node_xyz(j,np2)-node_xyz(j,np)

       vector2 = 0.0_dp
!!!! MHT 021221 - this is wrong here: vector2 should reflect the actual gravitational vector. here it assumes Z axis
!       vector2(3) = node_xyz(3,np)+1.0_dp  ! *MHT original*
       vector2(3) = -1.0_dp                   ! *MHT change*
       alpha = angle_btwn_vectors(vector1,vector2)
       !if (alpha .gt.pi/2) print *, "Alpha > pi/2. Elem: ", ne

!!! settlement velocity w.r.t. flow direction; set equal to zero to deactivate gravity effect
       !tj: this is a flow!!!!
       vseff = cos(alpha)* vs * (pi*elem_field(ne_radius,ne)**2.0_dp) * (part_param%grav_factor)
       !vseff = cos(alpha)* vs * part_param%grav_factor ! TJ *change*

!!! note following:
!!! ne_part_vel was stored as nej_flow in cmiss
!!! ne_part_conc was stored as nj_source in cmiss
!!! ne_part_var was stored as nj_loss in cmiss

!!! CALCULATE PARTICLE VELOCITY (stored in elem_field(ne_part_vel,ne))
!!! note that the influence of the Reynolds number is neglected because
!!! for particles <10 micron they are nearly as fast as the fluid
       ref_flow = elem_field(ne_Vdot,ne) ! vf, the actual flow (instead of normalized)
       elem_field(ne_part_vel,ne) = ref_flow - vseff ! tj: this is a flow!!!!
       ! rel_flow = vp0 - vf + cos(alpha) * vs
       rel_flow = elem_field(ne_part_vel,ne) - ref_flow + vseff  ! relative flow particle / fluid

       !elem_field(ne_part_vel,ne) = ref_flow - vseff
       if(abbr.gt.-25.0_dp) elem_field(ne_part_vel,ne) = elem_field(ne_part_vel,ne) &
                                                          + rel_flow*exp(abbr)
!!! calculate spatial derivative of particle/ velocity
!       elem_field(ne_dpart_vel,ne) = 0.0_dp !(XP(nk,nv2,nej_flow,np2)/XP(1,nv2,nj_radius,np2)**2 &
!            -XP(nk,nv,nej_flow,np)/XP(1,nv2,nj_radius,np)**2 )/elem_field(ne_length,ne)/PI    ! spatial derivative of particle velocity
!       XP(2,nv2,nej_flow,np2) = (XP(nk,nv,nej_flow,np)/XP(1,nv2,nj_radius,np)**2 &
!            -XP(nk,nv2,nej_flow,np2)/XP(1,nv2,nj_radius,np2)**2)/elem_field(ne_length,ne)/PI    
    enddo
    !stop
    !print *, "vf", elem_field(ne_Vdot, 40492), "vp", elem_field(ne_part_vel,40492)
    call enter_exit(sub_name,2)
    !print *, 'vseff', vseff
  end subroutine particle_velocity

!!!#######################################################################################
  

  subroutine particle_deposition(dt,inspiration,part_param, num_iter, DepFrac)

!!!    Calculates amount of deposited particles in each element,
!!!    reduces respective concentration field 'nj_source' and
!!!    stores the accumulated deposition quantity in field 'nj_loss'
!!!    ADDIONALLY WRITES current solution vector YP to source
!!!    field 'nj_source'

    use arrays, only: dp, particle_parameters, node_field, elem_field, &
        part_acinus_field,elem_cnct,elem_nodes, elems_at_node, node_xyz,&
        num_nodes
    !use arrays
    use indices
    use other_consts
    use mesh_utilities,only: angle_btwn_points,angle_btwn_vectors,cumulative_branch_length
    use diagnostics, only: enter_exit

    implicit none

    type(particle_parameters) :: part_param
    real(dp),intent(in) :: dt
    logical,intent(in) :: inspiration
    !     Local Variables
    integer :: gen,i,j,maxi,ne,ne0,noelem,np,npend,npstart, &
         npstep,np0,np1,np2,np3,nunit
    integer :: igreat_flag
    real(dp) :: A(3),abbr(9),abbr_t1,abbr_t2,alpha,Atube,B(3),C(3),coeffDiffSph(1:1000,2), &
         crossec(0:9),current_volume,Dalv,deltaV(-1:9), h,j2,lduct,length,radius(0:9),Rin, &
         Vduct(0:9),veloc(0:9),vfluid,volume(9),vol_croot_scale,Vtot,Vtube,Z(3)!Vdep(0:6),
    real(dp), intent(out) :: DepFrac(4:6)
    real(dp) :: Vdep(0:7), volflow
    real(dp) :: Ccun
    ! TJ add sedimentation in alveoli
    real(dp) :: mass_loss_np,mean_conc ,part_acinus_old(9)
    real(dp),allocatable :: part_concentration(:)
    logical :: acinar_deposition = .true., flag
    character(len=60) :: sub_name
    !!! TJ - to enhance particle deposition
    real(dp) :: unit_loss(0:9), unit_loss_sed(0:9), unit_loss_dif(0:9)
    real(dp) :: reduced_alveolar_volume(0:9), prev_reduced_alveolar_volume(0:9)

    integer, intent(in) :: num_iter
    real(dp) :: typical_time !T, kappa

    sub_name = 'particle_deposition'
    call enter_exit(sub_name,1)

    allocate(part_concentration(num_nodes))


    Ccun = 1.0_dp + 2.0_dp * part_param%lambda/part_param%pdia &
         * (1.257_dp+0.4_dp*exp(-0.55_dp*part_param%pdia/part_param%lambda))

!!! Copy particle 'concentration' solution to temporary array
    ! TJ - where to put extrathoracic deposition?
    part_concentration(:) = node_field(nj_conc1,:)!*(1-part_param%Eo)

    if(inspiration)then
       npstart = 1           ! start from node 1 to calculate deposition -- loop in flow direction
       npend = num_nodes    ! REQUIRES correct order of nodes!!!!!!!!!!!!!!!!
       npstep = 1
    else
       npstart = num_nodes   ! start from terminal nodes to calculate deposition -- loop in flow direction
       npend = 1
       npstep = -1
    endif

!!! in below loops, assuming the length and radius changes have already been done
    do np = npstart,npend,npstep ! loop over all nodes in the model
        Vtot = 0.0_dp
        Vdep = 0.0_dp    ! deposition volume
       !!! loop over surrounding elements
       do noelem = 1,elems_at_node(np,0)
          ne = elems_at_node(np,noelem)
          if(np.eq.elem_nodes(1,ne))then
             np2 = elem_nodes(2,ne)
             A(:) = node_xyz(:,np)  ! *MHT change*
             B(:) = node_xyz(:,np2) ! *MHT change*
          else
             np2 = elem_nodes(1,ne)
             A(:) = node_xyz(:,np2) ! *MHT change*
             B(:) = node_xyz(:,np)  ! *MHT change*
          endif
    ! *MHT change* for following: not clear if should be branch length of segment length? currently original
          lduct = cumulative_branch_length(ne) ! this is the length of the whole (unrefined) branch, i.e. between bifurcations
          !element (tube) length - since loop is over every node adjacent element, half length is used
          length = elem_field(ne_length,ne)/2.0_dp
          !lduct = length ! *MHT change*- actually keeps the same as original
          radius(0) = elem_field(ne_radius,ne) ! element (tube) radius
          Atube = pi*radius(0)**2.0_dp ! element (tube) cross section
          Vtube = Atube*length ! element (tube) volume
          Vtot = Vtot+Vtube  ! sum up total volume
          vfluid = elem_field(ne_Vdot,ne)/Atube  ! fluid velocity [mm/s]

!!!! MHT 021221 - this is wrong here: vector2 should reflect the actual gravitational vector. here it assumes Z axis
!!! CALCULATE ANGLE BETWEEN TUBE AXIS AND GRAVITATION DIRECTION
          z = 0.0_dp ! *MHT change*
          z(3) = -1.0_dp  ! *MHT change*
          alpha = abs(angle_btwn_vectors(B-A,Z)) ! calculate angle between tube axis and gravity vector! *MHT change*

!! sedimentation
          ! *MHT  -original ! TJ change back, should be prho (the density of the sphere)
!          Vdep(1) = (Vdep(1) + abs(sin(alpha))*part_param%prho* 9.81e3_dp * &
!                    part_param%pdia**2.0_dp*length*radius(0)*dt/9.0_dp/part_param%mu) &
!                    * part_param%grav_factor
          ! *MHT change*
          mean_conc = 0.5_dp * (node_field(nj_conc1,np) + node_field(nj_conc1,np2))
          Vdep(1) = Vdep(1) + abs(sin(alpha)) * mean_conc * 9.81e3_dp * part_param%pdia**2.0_dp * &
               length*radius(0)*dt/9.0_dp/part_param%mu ! deposition volume due to sedimentation
!          Vdep(1) = Vdep(1) + abs(sin(alpha)) * mean_conc * 9.81e3_dp * part_param%pdia**2.0_dp * &
!               length*dt/18.0_dp/part_param%mu ! deposition volume due to sedimentation ! TJ change

!!! impaction
!!! MHT *change: should be using the angle between an unfrefined branch and its parent, then taking
!!!! just a proportion of this for impaction
          if(inspiration)then ! impact only during inhalation (why?)
              ne0 = element_top_of_branch(ne)  !originally was elem_cnct(-1,1,ne) ! parent element
             if(elem_cnct(-1,1,ne0).gt.0)then ! parent element exists, means not inflow location
! was here but np0 not used
                np1 = elem_nodes(1,ne0)
                np2 = elem_nodes(2,ne0)
                np3 = elem_nodes(1,elem_cnct(-1,1,ne0)) ! parent start node
!!! CALCULATE ANGLE TO DOWNSTREAM BRANCH AS WELL !WHAT ABOUT IMPACT DURING EXHALATION????
                A(:) = node_xyz(:,np1) ! A is the connection node of parent and daughter
                B(:) = node_xyz(:,np2)
                C(:) = node_xyz(:,np3)
                alpha = angle_btwn_points(B,A,C)-pi ! calculate angle between parent and daughter tube axis

                if(abs(sin(alpha)).gt.zero_tol)then
          !*MHT  -original ! TJ change back
!                   abbr(1) = 18.0_dp*part_param%mu/part_param%prho&
!                                                  /part_param%pdia**2.0_dp&
!                                                  /abs(sin(alpha))
          ! *MHT change*
                   abbr(1) = 18.0_dp*part_param%mu/mean_conc &
                                                  /part_param%pdia**2.0_dp&
                                                  /abs(sin(alpha))
                   ! \bar(vp) eqn(38)
                   h = 2.0_dp*abs(sin(alpha))/lduct*abs(vfluid)/ &
                        abbr(1)*abs(elem_field(ne_part_vel,ne))/ &   !!! warning: note that cmiss used nodal value for particle velocity
                        abs(elem_field(ne_radius,ne))**2.0_dp/pi  !!! but only one version for value (??)

                      !h = h * (1.0_dp-exp(-abbr(1)*lduct/abs(vfluid)))

                   if(abbr(1).lt.25.0_dp)then
                      h = h * (1.0_dp-exp(-abbr(1)*lduct/abs(vfluid)))
                   endif
                   Vdep(2) = Vdep(2)+2.0_dp*h*dt*length*2.0_dp*radius(0)
                endif
             endif         !ne0
          endif!inspiration

!!! brownian diffusion
          ! tj - 04 mar 2024 - test for exhale
          if(abs(vfluid).gt.zero_tol)then
              !h = 2.0_dp/3.0_dp*(4.0_dp*part_param%diffu*lduct/vfluid)**0.5_dp/pi
             h = 2.0_dp/3.0_dp*(4.0_dp*part_param%diffu*lduct/abs(vfluid))**0.5_dp/pi ! average travelling distance within duct
          else
             h = 0.0_dp           ! equation not valid for keeping breath (??)
          endif
          Vdep(3) = Vdep(3) + pi * vfluid*h* ((2.0_dp*radius(0)-h)*length/lduct) * dt


!!!!!!!! TJ - Acinar deposition
!!!! calculate deposition in acini
          if(num_units.gt.0.and.elem_cnct(1,0,ne).eq.0.and.acinar_deposition)then
              nunit = where_inlist(ne,units) ! get the unit number
              current_volume = unit_field(nu_vol,nunit) ! current volume of acinus

          !!! scaling factor for either radius or length from Haefeli-Bleuer & Weibel TLC size to current size
          !!! TLC volume - dependent average diameter of an alveolus (Weibel, 1962). original used total lung volume. here we assume 2^15 acini
          ! TJ - 29 Mar 2024 - this is lambda for isotropic expansion.
              vol_croot_scale = (current_volume/part_param%VtotTLC)**(1.0_dp/3.0_dp)
              !Dalv = 0.154_dp * vol_croot_scale
              !Dalv = 1.54e-3_dp * vol_croot_scale ! thesis pg 25
              !ARC should be of acinar volume only?
              Dalv = 1.54e-3_dp * (current_volume*2.0_dp**15/part_param%VtotTLC)**(1.0_dp/3.0_dp) ! [mm] ! *MHT change*

              !Dalv = 0.154_dp  ! *MHT change*

             !!! inner radius for diffusion in alveoli related to inlet diameter (Hansen 1975)
              Rin = Dalv* ((1.0_dp -0.853_dp)*0.853_dp)**0.5_dp !choi & kim (2007)
              !Rin = Dalv*0.325_dp
              abbr(7) = Rin/Dalv*2.0_dp
              abbr(8) = pi*abbr(7)
              ! abbreviation for acinar diffusive deposition
              abbr(2) = -4.0_dp*pi**2.0_dp*part_param%diffu/Dalv**2.0_dp
              flag = .true.
              maxi = 1

! TJ - 28 Nov 2022 - seems this is random branching, is this correct?
              j = 0 ! TJ - j is the number of connecting tube (Koblinger and Hofmann, 1990)
              do while(flag)   ! find loop number for series in acinar diffusion
                  j = j+1
                  j2 = real(j)**2.0_dp
                  coeffDiffSph(j,1) = (-1.0_dp)**(real(j)+1.0_dp)/j2* &
                         (sin(real(j)*abbr(8))/real(j)/pi**2.0_dp - &
                         cos(real(j)*abbr(8))*abbr(7)/pi) ! set up coefficients for sum to speed up
                  coeffDiffSph(j,2) = j2*abbr(2)
                  h = exp(coeffDiffSph(j,2)*dt)/j2
                  if((h.lt.5.0e-7_dp).or.j.eq.500)then
                       maxi = j+1         ! set loop number for series in acinar diffusion
                       flag = .false.
                  endif
              enddo ! while
              maxi = maxi-1 ! not here in the original

              radius(0) = elem_field(ne_radius,ne) !mm
              crossec(0) = radius(0)**2.0_dp*pi ! cross-sectional area of terminal bronchiole in FEM model
            ! difference of FEM-mesh and measurements (HAEFELI-BLEUER,1988)

              !veloc(0) = elem_field(ne_part_vel,ne)/(pi *radius(0)**2.0_dp) ! velocity terminal bronchi

              volflow = abs(elem_field(ne_part_vel,ne))                ! volume flow particle in terminal bronchiole
              if(volflow .lt.0.0_dp) volflow = 0.0_dp ! TJ -change flow to zero
              veloc(0) = volflow/crossec(0) ! unit ????

              deltaV(-1) = dt* elem_field(ne_Vdot,ne) ! total volume change of acinus in one time step
              deltaV(0) = 0.0_dp ! assumption that volume change very small since no alveoli

              abbr(3) = radius(0) - radius(0)*vol_croot_scale
              abbr(4) = 0.0_dp  ! variable to store axial position for radius scaling

              reduced_alveolar_volume(:) = 0.0_dp
              prev_reduced_alveolar_volume(:) = 0.0_dp
              DepFrac(:) =0.0_dp

              do gen = 1,9 ! loop over acinar generations - deposition efficiency the same in all gens
                  abbr(4) = abbr(4)+part_param%LacTLC(gen+1)*vol_croot_scale ! axial position of node (length of acinar duct)
!                  radius(gen) = (part_param%totacinarLength -abbr(4))/(part_param%totacinarLength *abbr(3)) &
!                        + part_param%RacTLC(gen+1)*vol_croot_scale
                  radius(gen) = part_param%RacTLC(gen+1)*vol_croot_scale
                  crossec(gen) = (pi*radius(gen)**2.0_dp) *(2.0_dp**gen) ! accumulated duct cross-sectional area
                  Vduct(gen) = crossec(gen)* (part_param%LacTLC(gen+1)*vol_croot_scale)  ! duct volume in acinar region

                  !volume(gen) =(part_param%VacTLC(gen) * (current_volume/part_param%VtotTLC) - Vduct(gen))

                  ! TJ - mar 19 change - as here we do not have a clear out mechanism.
                  volume(gen) =(part_param%VacTLC(gen) * (current_volume/part_param%VtotTLC) &
                          - Vduct(gen) - reduced_alveolar_volume(gen))*(1-DepFrac(6))


!!!.............diffusion in alveolar tissue
               ! TJ - 28 Nov 2022 - it does not match with Falko's thesis
                  abbr(5) = 0.0_dp
                  abbr(6) = 0.0_dp
                  flag = .true.
                  j = 1
                  do while(flag)
                   !do j = 1,maxi
                       abbr_t1 = abbr(5) ! new param
                       abbr_t2 = abbr(6)
                  !!! TJ - here we treat part_acinus_field(10+gen, nunit) as bar(tau)
                       ! bar(tau)^(t) = tau(t-1) * (conc(t-1) * Va(t-1)/conc(t)/Va(t)) + dt, obtain by acinus_transport
    ! part_acinus_field(10+gen, nunit) = bar(tau), mean residence time
                       abbr(5) = abbr(5)+coeffDiffSph(j,1)*exp(coeffDiffSph(j,2)*part_acinus_field(10+gen,nunit))
                       abbr(6) = abbr(6)+coeffDiffSph(j,1)*exp(coeffDiffSph(j,2)*(part_acinus_field(10+gen,nunit)+dt))
                       if((abs(abbr_t1-abbr(5)).lt.1.0e-8_dp).and.(abs(abbr_t2-abbr(6)).lt.1.0e-8_dp)) flag = .false.
                       if(j.ge.maxi) flag = .false.
                       j = j + 1
                  enddo !j

                  if((abbr(6).gt.0.0_dp).and.(abbr(5).gt.abbr(6))) then !(abs(abbr(6)-abbr(5)).ge.zero_tol))then
                       ! TJ - add (abbr(6).le.abbr(5)) to avoid negative DepFrac(4)
                       ! can become zero for large T (accuracy of real*8)
                       ! diffusion fraction out of a sphere (Diffusion,Jost,1960) (0.853d0 is area correction ChoiKim2007)
                       ! TJ - 29 NOV 2022 - change back to solve negative DepFrac(4) issue
                       ! TJ - 08 Dec 2022 - deprecated 0.853 for further detection
                       DepFrac(4) = 1.0_dp - (abbr(6)/abbr(5))
                  else
                       DepFrac(4) = 0.0_dp
                  endif

                  DepFrac(5) = part_acinus_field(1+gen,nunit) * 9.81e3_dp & !part_param%gravityy &
                                 *part_param%pdia**2.0_dp &
                                 *dt/12.0_dp/part_param%mu/Dalv  !*0.853_dp

               ! sum deposition fractions without mutually eliminating part
                  DepFrac(6) = DepFrac(4)/2.0_dp + DepFrac(5) + DMAX1(DepFrac(4)/2.0_dp - DepFrac(5), 0.0_dp)

                  do j = 4,6
                        if(DepFrac(j).ge.0.2_dp)then
                            if(DepFrac(j).gt.1.0_dp) DepFrac(j) = 1.0_dp ! check if deposition volume extends volume in Acinus
                        endif
                  enddo !j

!!!!............ TJ - 12 mar 2024  ------------------------------------------------------

                   ! volume change each generation within on time step
                  deltaV(gen) = (part_param%VacTLC(gen)/part_param%VtotTLC)*deltaV(-1)
                       ! velocities in acinar generations

                  ! TJ - 30 Apr 2024 - check this with his thesis.
                  !veloc(gen) = (veloc(gen-1)*crossec(gen-1) - deltaV(gen-1)/dt) /crossec(gen)  !* Origin*
                  veloc(gen) = (veloc(gen-1)*crossec(gen-1) - deltaV(gen)/dt) /crossec(gen) !thesis pg37 EQN 90

              !!!............ sedimentation in alveolar tissue
!                  Vdep(5) = (2.0_dp**gen)*2.0_dp/pi *part_acinus_field(1+gen,nunit)&
!                                *9.81e3_dp *part_param%pdia**2.0_dp &
!                                *part_param%LacTLC(gen+1)*vol_croot_scale*radius(gen)*dt &
!                                /(9.0_dp*part_param%mu) * part_param%grav_factor
                  Vdep(7) = 0.0_dp

                  Vdep(5) = (pi * part_acinus_field(1+gen,nunit) * 9.81e3_dp &
                                *part_param%pdia**2.0_dp *dt *Dalv**2 &
                                /(72.0_dp*part_param%mu))* part_param%grav_factor * (2.0_dp**gen)

              !!! ----------- Brownian diffusion
                  if (.not. inspiration) then
                    !!! TJ - MAR 24 2024 - REMEDY FOR BROWNIAN DIFFUSION under exhalation.
                      h = (Rin**3/9/part_param%diffu)* pi/ num_iter
                      typical_time = (Rin)**2 /part_param%diffu/6

                      if (dt > typical_time) then
                          Vdep(4) = 4*pi*(Rin)**2 * h* dt * (2.0_dp**gen)
                      else
                          Vdep(4) = 0.0_dp
                      endif
!                      Vdep(4) = 0.0_dp
                      !Vdep(5) = 0.0_dp
                      Vdep(6) = Vdep(4)/2.0_dp + (Vdep(7)+Vdep(5)) + DMAX1(Vdep(4)/2.0_dp-(Vdep(7)+Vdep(5)),0.0_dp)
                  else
                   !!----------- Brownian diffusion - Falko Schmidt 2011  --------------------------------
                   !! (average travelling distance within alveolar duct)
                      h = (2.0_dp/3.0_dp) * (4.0_dp*part_param%diffu * part_param%LacTLC(gen+1)& !*vol_croot_scale&
                                /abs(veloc(gen)))**0.5_dp /pi
                      if(h.gt.radius(gen))then ! all particles are deposited
                          !TJ - 08 APL -TEST
                           Vdep(4) = deltaV(gen)!Vduct(gen)
                           Vdep(6) = deltaV(gen)!Vduct(gen)
                      else
                           ! deposition with referencing to time step DT
                           Vdep(4) = pi*h* (2.0_dp*radius(gen)-h) *dt*abs(veloc(gen)) * (2.0_dp**gen)
                              !!Total alveolar deposition
                              ! sum deposition volumes without mutually eliminating volume
        !                      Vdep(6) = Vdep(4)/2.0_dp+Vdep(5)+DMAX1(Vdep(4)/2.0_dp-Vdep(5),0.0_dp)
                           Vdep(6) = Vdep(4)/2.0_dp + (Vdep(7)+Vdep(5)) + DMAX1(Vdep(4)/2.0_dp-(Vdep(7)+Vdep(5)),0.0_dp) !TJ change
                      endif
!               !!! --------------- Brownian diffusion - Falko Schmidt 2011 end --------------------------------
                  endif    ! not inhalation

                  !!! TJ - 2024 MAY 03 - There should be a constraint on particle deposition volume or mass,
                  !!! deposition volume cannot exceed either vduct or volume as time increases

                  do j = 4,6 ! modify vdep4 and vdep6 for large diffusion
                      if(Vdep(j).gt.(0.2_dp*Vduct(gen)))then
                             ! check if deposition volume extends volume in Acinus
                          if(Vdep(j).gt.Vduct(gen)) then
                              Vdep(j) = Vduct(gen)!Vduct(gen) TJ MAY 14 LOCAL 17
!                                   write(*,'('' Warning: mechanism '',I2,'' deposits '' , D12.6, &
!                                   '' %  in acinar generation '', I2, '', reduce time step'')') &
!                                   j, Vdep(j)*100.0_dp/deltaV(gen) , gen
!                                   exit
                          endif
                      endif
                     ! if(Vdep(j).ne.Vdep(j)) Vdep(j) = 0.0_dp ! function ISNAN does not work
                  enddo !j


                  reduced_alveolar_volume(gen) = reduced_alveolar_volume(gen) + Vdep(6)

                  unit_loss_dif(gen) = Vdep(4) * part_acinus_field(1 + gen, nunit) &
                          + deltaV(gen)*DepFrac(4)*part_acinus_field(1+gen,nunit)

                  unit_loss_sed(gen) = Vdep(5)*part_acinus_field(1+gen,nunit)&
                          + deltaV(gen)*DepFrac(5)*part_acinus_field(1+gen,nunit)

                  unit_loss(gen) = Vdep(6)*part_acinus_field(1+gen,nunit) &
                          + deltaV(gen)*DepFrac(6)*part_acinus_field(1+gen,nunit)
                  if (volume(gen) .lt. 0.0_dp) print *, 'error: volume(gen)', volume(gen)


                  part_acinus_field(1+gen,nunit) = part_acinus_field(1+gen,nunit) &
                       - part_acinus_field(MIN(gen+2,10),nunit) * volume(gen)*DepFrac(6)/ (volume(gen)+Vduct(gen))  &
                       - part_acinus_field(1+gen,nunit) * Vdep(6)/ (volume(gen)+Vduct(gen)) ! assign new concentration acini


                  if(part_acinus_field(1+gen,nunit).ge.zero_tol)then
                   ! NOTHING because thsi way NaN values are covered too
                  else
                      part_acinus_field(1+gen,nunit) = 0.0_dp
                  endif

                  if(part_acinus_field(1+gen,nunit).gt.zero_tol)then
                      part_acinus_field(10+gen,nunit) = part_acinus_field(10+gen,nunit)+dt ! increase time particles spent in alveoli (in BBM(11..19,ne))
                  else
                      part_acinus_field(10+gen,nunit) = 0.0_dp ! avoid division by zero
                  endif   ! part_acinus_field
              enddo      ! for nine assumed acinar generations

              unit_loss(0) = sum(unit_loss(1:9))
              unit_loss_dif(0) = sum(unit_loss_dif(1:9))
              unit_loss_sed(0) = sum(unit_loss_sed(1:9))

              unit_field(nu_loss,nunit) = unit_field(nu_loss,nunit) + unit_loss(0)
              unit_field(nu_loss_dif,nunit) = unit_field(nu_loss_dif,nunit)+ unit_loss_dif(0)
              unit_field(nu_loss_sed,nunit) = unit_field(nu_loss_sed,nunit) + unit_loss_sed(0)
          endif         ! if a terminal element with acini attached

!          if (nunit == 16384)  then
!              print *, 'DepFrac(4:6)', DepFrac(4:6)
!              print *, 'Vdep(4:6)', Vdep(4:6)
!              print *, 'conc',part_acinus_field(2:10,nunit)
!!           print *, "deltaV(-1) ", deltaV(-1) !"reduced_alveolar_volume(:)", reduced_alveolar_volume(:)
!!           print *, "unit_field(nu_Vdot0)", unit_field(nu_Vdot0, num_units)
!          endif

       enddo            ! noelem

       !!! summation of all deposition effects
       ! sum up diffusion and sedimentation without mutually eliminating volume

        if(.not.inspiration) Vdep(2)=0.0_dp

       Vdep(0) = Vdep(3)/2.0_dp + Vdep(1) + max(Vdep(3)/2.0_dp - Vdep(1),0.0_dp)
       ! add impaction with correction term
       Vdep(0) = Vdep(0) + Vdep(2) - Vdep(0) * Vdep(2)/Vtot

       mass_loss_np = Vdep(0)*part_concentration(np)

!!! store the deposition quantities of each type
       ! *MHT change* : the field names were mislabelled (Vdep(1) and Vdep(3) mixed up
!!! Vdep(3) is diffusion, Vdep(2) is impaction, Vdep(1) is sedimentation
       node_field(nj_loss,np) = node_field(nj_loss,np)+mass_loss_np ! deposition quantity [g]
       node_field(nj_loss_sed,np) = node_field(nj_loss_sed,np)+Vdep(1)*part_concentration(np) ! deposition quantity [g]
       node_field(nj_loss_imp,np) = node_field(nj_loss_imp,np)+Vdep(2)*part_concentration(np) ! deposition quantity [g]
       node_field(nj_loss_dif,np) = node_field(nj_loss_dif,np)+Vdep(3)*part_concentration(np) ! deposition quantity [g]

       !! assign new particle concentration to all nodes
       part_concentration(np) = (part_concentration(np)*Vtot - mass_loss_np)/Vtot
    enddo !np

!!! copy nj_source (concentration - deposition) field to nj_conc1 field
    node_field(nj_conc1,:) = part_concentration(:)

    deallocate(part_concentration)

    call enter_exit(sub_name,2)

  end subroutine particle_deposition



!!! ############################################################
!!! ############################################################



  subroutine acinus_transport(nunit,dt,part_param)

!!! Evaluates particle transport in acini governed by convection diffusion equation
!!! using FDM for 9 'virtual' nodes; using operator splitting (FDM for diffusion and
!!! dV/dt=const for each generation for convection);
!!! geometric data are based on HAEFELI-BLEUER and WEIBEL,1988
!!! Created September, 2011, Falko Schmidt

!   use arrays, only: dp, particle_parameters, elem_field, node_field, &
!                     elem_nodes, part_acinus_field
    use arrays
    use indices
    use other_consts, only: pi
    use diagnostics, only: enter_exit

    implicit none
    type(particle_parameters) :: part_param
    integer :: nunit
    real(dp) :: dt
!!!     Local variables
    integer :: i,ne,np,gen,genBEG,genEND,genSTEP!, igreat_flag
    real(dp) :: current_volume,Lac_TLC_gen_p1,Lac_TLC_gen_p2,radius(-1:9),conc(-1:9), &
         volume(0:9),crossec(-1:9),veloc(0:9),abbr(6), &
         deltaV(-1:9),volflow,vol_croot_scale,courant,effDiffu, &
         part_acinus_old(9),alv_mass_end,alv_mass_start,scale_mass
    character(len=60) :: sub_name

    sub_name = 'acinus_transport'
    call enter_exit(sub_name,1)

    ne = units(nunit)
    np = elem_nodes(2,ne) ! the distal node; adjacent to the acinus
    current_volume = unit_field(nu_vol,nunit)
    vol_croot_scale = (current_volume/part_param%VtotTLC)**(1.0_dp/3.0_dp)
    alv_mass_start = unit_field(nu_vol,nunit)*unit_field(nu_conc1,nunit) - unit_field(nu_loss,nunit)

!!! initialize parameter - CAUTION! FEM model and data above may not match (radius adapted)
    radius(0) = elem_field(ne_radius,ne)
    crossec(0) = radius(0)**2.0_dp * pi                                  ! cross-sectional area terminal bronchiole FEM model
    conc(0) = node_field(nj_conc1,np)                                    ! concentration in terminal bronchiole, defined that mass is concerved
    if(node_field(nj_conc1,np).lt.0.0_dp)then
       !ARC temp write(*,'('' Warning: particle concentration < 0 at unit'',i6)') nunit
       conc(0) = 0.0_dp
       node_field(nj_conc1,np) = 0.0_dp
    endif

    volflow = abs(elem_field(ne_part_vel,ne))                                 ! volume flow particle in terminal bronchiol
    !!! TJ - without this we receive exhalation wrong
    if(volflow .lt.0.0_dp) volflow = 0.0_dp ! TJ -change flow to zero

    veloc(0) = volflow/crossec(0)
    ! total volume change in one time step
    deltaV(-1) = elem_field(ne_Vdot,ne)*dt

    ! assumptions volume change very small since no alveoli
    deltaV(0) = 0.0_dp

    ! difference of FEM-mesh and measurements (HAEFELI-BLEUER,1988)
!!    abbr(1) = radius(0)-scale_radius_acinus(part_param%RacTLC(1),current_volume)
    abbr(1) = radius(0) - radius(0)*vol_croot_scale

    ! variable to store axial position for radius scaling
    abbr(2) = 0.0_dp

    ! ------------------------------------------------
    do gen = 1,9  ! loop over acinar generations
    ! ------------------------------------------------
        !!!! TJ - MAR 20 2024 - SHOULD THIS MATCH WITH subroutine paticle_deposition???
        abbr(2) = abbr(2)+part_param%LacTLC(gen+1)*vol_croot_scale
        radius(gen) = part_param%RacTLC(gen+1)*vol_croot_scale
        crossec(gen) = pi*radius(gen)**2.0_dp*(2.0_dp**gen)               ! accumulated duct cross-sectional area
        !volume(gen) = part_param%VacTLC(gen) * (current_volume/part_param%VtotTLC)
!       radius(gen) = part_param%RacTLC(gen+1)*vol_croot_scale + &     ! scale duct radius wrt. current acinar volume and
!            (part_param%totacinarLength*vol_croot_scale-abbr(2))/(part_param%totacinarLength*vol_croot_scale)*abbr(1) !linear scaling to adapt HBW to FEM mesh
        volume(gen) = DMAX1(part_param%VacTLC(gen)/part_param%VtotTLC*current_volume,0.0_dp)    ! linear scaling of volume
        deltaV(gen) = part_param%VacTLC(gen)/part_param%VtotTLC*deltaV(-1)                      ! volume change each generation within on time step
        conc(gen) = part_acinus_field(1+gen,nunit)                        ! concentration in acinar generations
        veloc(gen) = (veloc(gen-1)*crossec(gen-1)-deltaV(gen-1)/dt) &
                         /crossec(gen)                                   ! velocities corrected by volume change of previous generations
        courant = abs(veloc(gen))*dt/(part_param%LacTLC(gen+1)*vol_croot_scale)                          ! courant number for stability check

       ! warning if CFL-condition is not fulfilled
       if(courant.ge.0.5_dp)then
          write(*,'(''Transport in acinar region may be unstable, Cr='',d12.3, &
               '' decrease DT to at least '',d12.3,''s'')') &
               courant,0.5_dp*(part_param%LacTLC(gen+1)*vol_croot_scale)/abs(veloc(gen))                ! warning if CFL-condition is not fulfilled
          !min_dt = min(0.5_dp*(part_param%LacTLC(gen+1)*vol_croot_scale)/abs(veloc(gen)), min_dt)
           !pause
       endif

       part_acinus_old(gen) = part_acinus_field(1+gen,nunit)                      ! store old concentrations
    enddo !gen

!!! define direction-explicit transport
    if(volflow.ge.0.0_dp)THEN
       genBEG = 1                                                      ! start generation convective transport
       genEND = 9                                                      ! end generation convective transport
       genSTEP = 1                                                     ! direction convective transport
    else                                                               ! else even time step
       genBEG = 8                                                      ! start generation convective transport
       genEND = 0                                                      ! end generation convective transport
       genSTEP = -1                                                    ! direction convective transport
       conc(-1) = node_field(nj_conc1,elem_nodes(1,ne))                ! concentration of proximal node of terminal element
       part_param%LacTLC(1) = elem_field(ne_length,ne)                         ! length of terminal element of 3D tree
       crossec(-1) = crossec(0)                                       ! cross-sectional area of proximal node of terminal element
       volume(0) = crossec(0)*part_param%LacTLC(1)*vol_croot_scale
    endif

!!! convective transport
    do gen = genBEG,genEND,genSTEP
       if(gen.le.8)then                                                ! inside domain use central scheme in space
          conc(gen) = (conc(gen)*(volume(gen)-deltaV(gen))+dt* &
               (conc(MIN0(gen-genSTEP,gen))*veloc(gen)*crossec(gen)- &
               (conc(MAX0(gen-genSTEP,gen))*veloc(gen+1)* &
               crossec(gen+1))))/volume(gen)                           ! new concentration by mass/volume balance
       else
          conc(gen) = (conc(gen)*(volume(gen)-deltaV(gen))+ &
               conc(MIN0(gen-genSTEP,gen))*veloc(gen)* &
               crossec(gen)*dt)/(volume(gen))                          ! new concentration in last generation by mass/volume balance
       endif
    enddo

!!! explicit FDM scheme - for diffusion
    do gen = 1,9            ! loop over acinar generations in alternating directions
        !effDiffu = part_param%diffu + 0.0867_dp* veloc(gen) *part_param%LacTLC(gen+1)*vol_croot_scale
       effDiffu = part_param%diffu + 0.0867_dp*abs(veloc(gen))*part_param%LacTLC(gen+1)*vol_croot_scale            ! effective diffusivity acc. to. Lee(2001)
       Lac_TLC_gen_p1 = part_param%LacTLC(gen+1) * vol_croot_scale
       if(gen.le.8)then                                                  ! inside domain use central scheme in space
          Lac_TLC_gen_p2 = part_param%LacTLC(gen+2) * vol_croot_scale
          abbr(1) = Lac_TLC_gen_p1 * Lac_TLC_gen_p2*&
                 (Lac_TLC_gen_p1 + Lac_TLC_gen_p2)       ! abbreviation (M)
          abbr(2) = Lac_TLC_gen_p1* dt /volume(gen)/abbr(1)
!          abbr(1) = part_param%LacTLC(gen+1)*part_param%LacTLC(gen+2)*(part_param%LacTLC(gen+1)+part_param%LacTLC(gen+2)) &
!               *vol_croot_scale**3.0_dp ! abbreviation to reduce computational effort ! M in EQN 94-95
          !abbr(2) = part_param%LacTLC(gen+1)*dt/volume(gen)/abbr(1) * vol_croot_scale
          part_acinus_field(1+gen,nunit) = conc(gen) + abbr(2) * effDiffu * &
               ((Lac_TLC_gen_p1**2.0_dp*(crossec(gen+1)-crossec(gen)) - Lac_TLC_gen_p2**2.0_dp*(crossec(gen)-crossec(gen-1))) &
               *(Lac_TLC_gen_p1**2.0_dp*(conc(gen+1)-conc(gen)) - Lac_TLC_gen_p2**2.0_dp*(conc(gen)-conc(gen-1)))/abbr(1) &
               +2.0_dp*crossec(gen)* (Lac_TLC_gen_p1*(conc(gen+1)-conc(gen)) - Lac_TLC_gen_p2*(conc(gen)-conc(gen-1))))                  ! FDM with central differences in space
       else                                                              ! domain boundary use different scheme
          abbr(1) = Lac_TLC_gen_p1**2.0_dp                                  ! abbreviation to reduce computational effort
          abbr(2) = Lac_TLC_gen_p1*dt/volume(gen)/abbr(1)
          part_acinus_field(1+gen,nunit) = conc(gen)+abbr(2)*effDiffu* &
               (2.0_dp*crossec(gen)*(conc(gen-1)-conc(gen)))             ! FDM with entral differences in space at boundary (using condition dC/dx=0)
       endif
    enddo !gen

!!! assign new concentrations to 'part_acinus_field' and calculate average time particles spend in acinus
    do gen = 1,9
       if(part_acinus_field(1+gen,nunit).lt.zero_tol) part_acinus_field(1+gen,nunit) = 0.0_dp
       if((volflow.gt.0.0_dp).and.(part_acinus_field(1+gen,nunit).gt.zero_tol))then  ! average time does not change when particles leave alveolus
          part_acinus_field(10+gen,nunit) = part_acinus_old(gen)*part_acinus_field(10+gen,nunit)* &
               (volume(gen)-deltaV(gen))/ &
               part_acinus_field(1+gen,nunit)/volume(gen)                   ! calculate new average time partices spent in alveoli (in part_acinus_field(11..19,ne))
       elseif(part_acinus_field(1+gen,nunit).le.zero_tol)then                 ! concentration is zero
          part_acinus_field(10+gen,nunit)=0.0_dp                            ! reset time
       endif ! volflow

      ! igreat_flag = 0
       !if (nunit.eq.1) then
       !    if ((volume(gen)-deltaV(gen)).gt.volume(gen)) igreat_flag = 1
       !    write(94,*) gen, part_acinus_field(10+gen,nunit), part_acinus_old(gen), igreat_flag
       !end if
    enddo


!!! for exhalation calculate concentration of terminal node of branching tree
    if(volflow.lt.0)then                                                 ! exhalation
       !write(*,*) 'warning: exhalation wrong in acinus_transport'
       gen = 0
       effDiffu = part_param%diffu + 0.0867_dp*abs(veloc(gen))*part_param%LacTLC(gen+1)         ! effective diffusivity acc. to. Lee(2001)

!       Lac_TLC_gen_p1 = part_param%LacTLC(gen+1) * vol_croot_scale
!       Lac_TLC_gen_p2 = part_param%LacTLC(gen+2) * vol_croot_scale
!
!       abbr(1) = Lac_TLC_gen_p1 * Lac_TLC_gen_p2*&
!                 (Lac_TLC_gen_p1 + Lac_TLC_gen_p2)       ! abbreviation (M)
!       abbr(2) = Lac_TLC_gen_p1 * dt /Lac_TLC_gen_p1/crossec(gen)/abbr(1)
!
!       node_field(nj_conc1,np) = conc(gen) + abbr(2) * effDiffu * &
!            ((Lac_TLC_gen_p1**2.0_dp * (crossec(gen+1)-crossec(gen)) &
!            - Lac_TLC_gen_p2**2.0_dp*(crossec(gen)-crossec(gen-1))) &
!            *(Lac_TLC_gen_p1**2.0_dp*(conc(gen+1)-conc(gen)) &
!            -Lac_TLC_gen_p2**2.0_dp*(conc(gen)-conc(gen-1)))/abbr(1) &
!            + 2.0_dp * crossec(gen)* &
!            (Lac_TLC_gen_p1 * (conc(gen+1)-conc(gen)) &
!            -Lac_TLC_gen_p2*(conc(gen)-conc(gen-1))))

       abbr(1) = part_param%LacTLC(gen+1)*part_param%LacTLC(gen+2)*&
                 (part_param%LacTLC(gen+1)+part_param%LacTLC(gen+2))     ! abbreviation
       abbr(2) = part_param%LacTLC(gen+1) * dt /part_param%LacTLC(gen+1)/crossec(gen)/abbr(1)

       node_field(nj_conc1,np) = conc(gen) + abbr(2) * effDiffu * &
            ((part_param%LacTLC(gen+1)**2.0_dp * (crossec(gen+1)-crossec(gen)) &
            -part_param%LacTLC(gen+2)**2.0_dp*(crossec(gen)-crossec(gen-1))) &
            *(part_param%LacTLC(gen+1)**2.0_dp*(conc(gen+1)-conc(gen)) &
            -part_param%LacTLC(gen+2)**2.0_dp*(conc(gen)-conc(gen-1)))/abbr(1) &
            + 2.0_dp * crossec(gen)* &
            (part_param%LacTLC(gen+1)*(conc(gen+1)-conc(gen)) &
            -part_param%LacTLC(gen+2)*(conc(gen)-conc(gen-1))))                     ! FDM with central differences in space

       if(node_field(nj_conc1,np).lt.0.0_dp)then
           write(*,'(''ACIN_TRANSP exhaling concentration < 0 ='',d12.3)') node_field(nj_conc1,np)
           node_field(nj_conc1,np) = 0.0_dp
       endif
    endif

    ! scale the concentrations to maintain acinar mass as per the unit mass
    alv_mass_end = 0.0_dp
    do i = 1,9
       alv_mass_end = alv_mass_end + part_acinus_field(1+i,nunit)* &
            (part_param%VacTLC(i) * unit_field(nu_vol,nunit)/part_param%VtotTLC)
    enddo

    !!! TJ - change **
!    if(alv_mass_end.gt.zero_tol)then
!       scale_mass = alv_mass_start/alv_mass_end
!       do i = 1,9
!          part_acinus_field(1+i,nunit) = part_acinus_field(1+i,nunit) * scale_mass
!       enddo
!    endif

    scale_mass = alv_mass_start/alv_mass_end
    do i = 1,9
        if(alv_mass_end.gt.zero_tol) part_acinus_field(1+i,nunit) = part_acinus_field(1+i,nunit) * scale_mass
        !if(abs(alv_mass_end).gt.zero_tol) part_acinus_field(1+i,nunit) = part_acinus_field(1+i,nunit) * scale_mass
        if(part_acinus_field(1+i,nunit).lt. zero_tol)  part_acinus_field(1+i,nunit) = 0.0_dp
    enddo

!    if (nunit == 16384) then
!        ! Print information for nunit 16384
!        print *, "-----------------------------------------------------------------"
!        print *, "Acinar Tansport Information for nunit 16384:"
!        print *, "loss", unit_field(nu_loss,nunit)
!        print *, "conc", part_acinus_field(1:10,nunit)
!        print *, "current_volume", current_volume, "unit_field(nu_vol,nunit)", unit_field(nu_vol,nunit)
!        print *, "tau", part_acinus_field(11:20,nunit)
!        print *, "volume(gen)", volume(:)
!        print *, "veloc", veloc(:)
!    endif

    call enter_exit(sub_name,2)

  end subroutine acinus_transport

!#######################################################################
  
  function scale_radius_acinus(initial_radius,volume_lung)

!!!   calculates new airway diameter from 
!!!   lung volume, equation is linear interploation of HABIB1994
!!!    REQUIRES reference data stored are from TLC
    use arrays, only: dp
    
    implicit none
  
!!!     Parameter List
    real(dp),intent(in) :: initial_radius,volume_lung
!!!     Local variables
    real(dp) :: alpha,beta,dref,fac
    real(dp),parameter :: TLC = 6.0e+3_dp, FRC = 3.6e+3_dp
    real(dp) :: scale_radius_acinus
      
    alpha = 1.51_dp ! parameter 
    beta = 0.442_dp ! parameter
    dref = 2.0_dp     ! reference diameter [mm]

    fac = alpha/(1.0_dp+(2.0_dp*initial_radius/dref)**(-beta))          ! scaling acc.to HABIB1994
    fac = DMIN1(fac,1.0_dp)                                             ! trachea otherwise would increase
    scale_radius_acinus = initial_radius*(volume_lung-FRC+fac*(TLC-volume_lung))/(TLC-FRC)  ! linear interploation for all lung volumes      

  end function scale_radius_acinus

!#########################################################################

  subroutine calc_mass_particles(nj,nu_field,gas_mass,deposit_mass,mass_by_gen, &
       part_param,unit_mass,unit_wall)
    
!   use arrays, only: dp,num_nodes,elem_nodes,node_field,elem_ordrs,&
!                     elem_field,elem_symmetry,num_elems
    use arrays
    use indices
    use diagnostics, only: enter_exit

    implicit none
  
    integer,intent(in) :: nj,nu_field
    type(particle_parameters) :: part_param
    real(dp) :: gas_mass,deposit_mass,mass_by_gen(:),unit_mass,unit_wall
    !     Local Variables
    integer :: i,ne,ne0,ngen,nmax_gen,np,np1,np2,nunit
    real(dp) :: average_conc,sum_mass,total_mass,vol_by_gen(50)
    real(dp),allocatable :: tree_mass(:),unit_tree_mass(:),wall_mass(:),unit_wall_mass(:)
    character(len=60):: sub_name

    !...........................................................................
    
    sub_name = 'calc_mass_particles'
    call enter_exit(sub_name,1)
    
    if(.not.allocated(tree_mass)) allocate(tree_mass(num_nodes))
    if(.not.allocated(unit_tree_mass)) allocate(unit_tree_mass(num_nodes))
    if(.not.allocated(wall_mass)) allocate(wall_mass(num_nodes))
    if(.not.allocated(unit_wall_mass)) allocate(unit_wall_mass(num_nodes))

    tree_mass = 0.0_dp
    unit_tree_mass = 0.0_dp
    wall_mass = 0.0_dp
    unit_wall_mass = 0.0_dp
    mass_by_gen = 0.0_dp
    vol_by_gen = 0.0_dp

    ! initialise to the mass in each element
    do ne = 1,num_elems
       np1 = elem_nodes(1,ne)
       np2 = elem_nodes(2,ne)
       average_conc = (node_field(nj,np1)+node_field(nj,np2))/2.0_dp
       tree_mass(ne) = average_conc*elem_field(ne_vol,ne) ! mmol/mm^3 * mm^3
       wall_mass(ne) = node_field(nj_loss,np1)/real(elems_at_node(np1,0)) + &
            node_field(nj_loss,np2)/real(elems_at_node(np2,0))
       elem_field(ne_mass,ne) = tree_mass(ne)
       elem_field(ne_loss,ne) = wall_mass(ne)/ &
            (pi*2.0_dp*elem_field(ne_radius,ne)*elem_field(ne_length,ne))
       ngen = elem_ordrs(1,ne)
       mass_by_gen(ngen) = mass_by_gen(ngen) + wall_mass(ne) + tree_mass(ne)
       vol_by_gen(ngen) = vol_by_gen(ngen) + elem_field(ne_vol,ne)
    enddo
    nmax_gen = maxval(elem_ordrs(1,:))

    ! add the mass in each elastic unit to terminal elements
    do nunit = 1,num_units
       ne = units(nunit)
       np1 = elem_nodes(1,ne)
       np2 = elem_nodes(2,ne)

       unit_tree_mass(ne) = unit_field(nu_vol,nunit)*unit_field(nu_field,nunit) &
            - unit_field(nu_loss,nunit) ! units of mass in acinus lumen
       unit_wall_mass(ne) = unit_field(nu_loss,nunit)      ! units of mass deposit
       !if (unit_tree_mass(ne) .lt. 0.0_dp)  unit_tree_mass(ne) = 0.0_dp
    enddo

    ! sum mass recursively up the tree
    do ne = num_elems,2,-1 ! not for the stem branch; parent = 0
       ne0 = elem_cnct(-1,1,ne)
       tree_mass(ne0) = tree_mass(ne0) + dble(elem_symmetry(ne))*tree_mass(ne)
       wall_mass(ne0) = wall_mass(ne0) + dble(elem_symmetry(ne))*wall_mass(ne)
       unit_tree_mass(ne0) = unit_tree_mass(ne0) + dble(elem_symmetry(ne))*unit_tree_mass(ne)
       unit_wall_mass(ne0) = unit_wall_mass(ne0) + dble(elem_symmetry(ne))*unit_wall_mass(ne)

       !if (unit_wall_mass(ne0) .lt. 0.0_dp) print *, 'error:', ne0

       elem_field(ne_mass,ne0) = elem_field(ne_mass,ne0) + dble(elem_symmetry(ne))*elem_field(ne_mass,ne)
    enddo !noelem

    !print *, "unit_wall_mass(ne)", size(unit_wall_mass(:))
    !pause
    gas_mass = tree_mass(1)        ! mass in lumen of bronchi
    deposit_mass = wall_mass(1)    ! mass on walls of bronchi
    unit_mass = unit_tree_mass(1)  ! mass in the lumen of the acini
    unit_wall = unit_wall_mass(1)  ! mass on alveolar walls

    deallocate(tree_mass)
    deallocate(unit_tree_mass)
    deallocate(wall_mass)
    deallocate(unit_wall_mass)

    call enter_exit(sub_name,2)
    
  end subroutine calc_mass_particles
  

!#########################################################################
!!!####################################################################

  subroutine airway_mesh_deform(dt,initial_volume,coupled,problem_type,tp,part_param)

!!! changes the size of elastic units and airways, based on pre-computed
!!! fields for flow into the units (ne_Vdot) and element volume change (ne_dvdt).
!!! The concentration in the units is adjusted to make sure that mass is conserved
!!! when the unit changes volume. If a unit changes in volume, then both the radius
!!! and length scale as the cube root of volume change. For alveolated airways (where
!!! a/A is < 1) the outer radius (indexed as nj_radius) is scaled and a/A unchanged.

!!! note that scaling the length or scaling length and radius gives no difference in mass error
    use indices
    use geometry,only: volume_of_mesh
    use arrays
    use diagnostics, only: enter_exit

    implicit none
    type(transport_parameters) :: tp
    type(particle_parameters) :: part_param
    real(dp),intent(in) :: dt,initial_volume
    logical,intent(in) :: coupled
    character(len=9) :: problem_type

    ! Local parameters
    integer :: i,ne,ne0,np,np1,np2,nunit
    real(dp) :: current_volume,new_volume,previous_volume,ratio, &
         sum_old_vols,tree_volume,unit_mass(2),volume_change
    character(len=60) :: sub_name


    sub_name = 'airway_mesh_deform'
    call enter_exit(sub_name,1)

    call volume_of_mesh(current_volume,tree_volume)

!!! calculate the change in volume of elastic units:
    ! mass(unit) = volume(unit) * concentration(unit)
    ! dv(unit) = flow(unit) * dt
    ! volume(new) = volume(old) + dv
    ! concentration(new) = mass/volume(new)
    do nunit = 1,num_units !for each of the elastic units
       ne = units(nunit) ! local element number
       if(ne.ne.0)then
          np = elem_nodes(2,ne) ! end node, attaches to unit
          volume_change = elem_field(ne_Vdot,ne)*dt
          !volume_change = elem_field(ne_Vdot,ne)/real(elems_at_node(np,0))*dt ! dv = q(unit)*dt
          if(.not.coupled)then ! this is where the unit volume is updated, if not false = if true
             previous_volume = unit_field(nu_vol,nunit)
             new_volume = unit_field(nu_vol,nunit) + volume_change
             unit_field(nu_vol,nunit) = new_volume
          else ! for coupled model, volume has already been updated
             previous_volume = unit_field(nu_vol,nunit) - volume_change
             new_volume = unit_field(nu_vol,nunit)
          endif

          if(elem_cnct(1,0,ne).eq.0)then ! no child branches
             unit_mass(1) = previous_volume*unit_field(nu_conc1,nunit) ! m = v(unit)*c
             unit_mass(2) = previous_volume*unit_field(nu_conc2,nunit) ! m = v(unit)*c
             unit_field(nu_vol,nunit) = new_volume ! i.e. unchanged if coupled
             ! important: adjust the concentration to maintain mass conservation
             if(.not.coupled)then
                unit_field(nu_conc1,nunit) = unit_mass(1)/unit_field(nu_vol,nunit)
                unit_field(nu_conc2,nunit) = unit_mass(2)/unit_field(nu_vol,nunit)
             endif
            select case(model_type)
              case('particle_transport') ! for problem type particles
               call acinus_transport(nunit,dt,part_param)
             end select
          endif
       endif
    enddo

!!! calculate the inflation/deflation of multi-branching acini (or other branches)
    ! dv(element) = dvdt(element) * total_volume_change
    ! ratio_of_volumes = (volume(element)+dv(element))/volume(element)
    ! volume(new) = volume(old) * ratio_of_volumes
    ! length(new) = length(old) * cuberoot(ratio_of_volumes)
    ! radius(new) = radius(old) * cuberoot(ratio_of_volumes)
    if(.not.coupled)then
       do ne = 1,num_elems
          !mht!       volume_change = elem_field(ne_dvdt,ne)*inlet_flow*dt
          np1 = elem_nodes(1,ne)
          np2 = elem_nodes(2,ne)
          volume_change = elem_field(ne_dvdt,ne)*dt
          
          ! adjust the concentrations at both nodes to maintain the mass
          if(abs(volume_change).gt.zero_tol)then
             ne0 = elem_cnct(-1,1,ne) ! parent element
             sum_old_vols = elem_field(ne_vol,ne0)
             do i = 1,elem_cnct(1,0,ne0)
                sum_old_vols = sum_old_vols + elem_field(ne_vol,elem_cnct(1,i,ne0))
             enddo
             node_field(nj_conc1,np1) = node_field(nj_conc1,np1) * sum_old_vols/&
                  (sum_old_vols+volume_change)

             ne0 = ne
             sum_old_vols = elem_field(ne_vol,ne0)
             do i = 1,elem_cnct(1,0,ne0)
                sum_old_vols = sum_old_vols + elem_field(ne_vol,elem_cnct(1,i,ne0))
             enddo
             node_field(nj_conc1,np2) = node_field(nj_conc1,np2) * sum_old_vols/&
                  (sum_old_vols+volume_change)
          endif

          ! adjust the element size
          ratio = (volume_change+elem_field(ne_vol,ne))/elem_field(ne_vol,ne)
          elem_field(ne_vol,ne) = elem_field(ne_vol,ne)*ratio
          elem_field(ne_length,ne) = elem_field(ne_length,ne)*(ratio)**(1.0_dp/3.0_dp)
          elem_field(ne_radius,ne) = elem_field(ne_radius,ne)*(ratio)**(1.0_dp/3.0_dp)

       enddo
    endif

!!! calculate the total model volume
    call volume_of_mesh(current_volume,tree_volume)

    ! TJ - FINAL total_volume_change equals tidal volume (target volume)
    tp%total_volume_change = current_volume-part_param%initial_volume

    call enter_exit(sub_name,2)

  end subroutine airway_mesh_deform
  
  
!#########################################################################

!  function branch_length(ne)
!!!! dummy arguments
!    integer :: ne
!    ! local variables
!    integer :: generation,ne0
!    real(dp) :: length
!    logical :: continue
!    real(dp) branch_length
!
!    generation = elem_ordrs(1,ne)
!    length = elem_field(ne_length,ne)
!    ne0 = elem_cnct(-1,1,ne)
!    continue = .true.
!    if(ne0.eq.0.or.elem_ordrs(1,ne0).ne.generation) continue = .false.
!
!    do while(continue)
!       length = length + elem_field(ne_length,ne0)
!       ne0 = elem_cnct(-1,1,ne0)
!       if(ne0.eq.0.or.elem_ordrs(1,ne0).ne.generation) continue = .false.
!    enddo
!
!    ne0 = elem_cnct(1,1,ne)
!    continue = .true.
!    if(ne0.eq.0.or.elem_ordrs(1,ne0).ne.generation) continue = .false.
!
!    do while(continue)
!       length = length + elem_field(ne_length,ne0)
!       ne0 = elem_cnct(1,1,ne0)
!       if(ne0.eq.0.or.elem_ordrs(1,ne0).ne.generation) continue = .false.
!    enddo
!
!    branch_length = length
!
!  end function branch_length

!#########################################################################

  function inlist(item,ilist)
!!! dummy arguments
    integer :: item,ilist(:)
! local variables
    integer :: n
    logical :: inlist

    inlist = .false.
    do n=1,size(ilist)
       if(item == ilist(n)) inlist = .true.
    enddo

  end function inlist

!!!###############################################################

  function where_inlist(item,ilist)
!!! dummy arguments
    integer :: item,ilist(:)
! local variables
    integer :: n
    integer :: where_inlist

    do n=1,size(ilist)
       if(item == ilist(n)) where_inlist = n
    enddo

  end function where_inlist

! #################################################################

  function element_lobe(ne)
      !!! change to find which lobe belongs to
    integer :: ne
    ! local variables
    integer :: generation,ne0
    logical :: continue
    integer :: element_lobe

    element_lobe = ne
    generation = elem_ordrs(1,ne)
    ne0 = elem_cnct(-1,1,ne)
    continue = .true.

    if(ne0.eq.0)then
       continue = .false.
    else if(elem_ordrs(1,ne0).le.3)then
       continue = .false.
    endif
    !if(elem_ordrs(1,ne0).le.3) continue = .false.

    do while(continue)
       element_lobe = ne0
       ne0 = elem_cnct(-1,1,ne0)
       if(ne0.eq.0)then
           continue = .false.
        else if(elem_ordrs(1,ne0).le.3)then
           continue = .false.
        endif
    enddo
  end function element_lobe


  ! #################################################################

  function element_top_of_branch(ne)
    integer :: ne
    ! local variables
    integer :: generation,ne0
    logical :: continue
    integer :: element_top_of_branch
    
    element_top_of_branch = ne
    generation = elem_ordrs(1,ne)
    ne0 = elem_cnct(-1,1,ne)
    continue = .true.

    if(ne0.eq.0)then
       continue = .false.
    else if(elem_ordrs(1,ne0).ne.generation)then
       continue = .false.
    endif
    !if(ne0.eq.0.or.elem_ordrs(1,ne0).ne.generation) continue = .false.
    
    do while(continue)
       element_top_of_branch = ne0
       ne0 = elem_cnct(-1,1,ne0)
       if(ne0.eq.0)then
          continue = .false.
       else if(elem_ordrs(1,ne0).ne.generation)then
          continue = .false.
       endif
       !if(ne0.eq.0.or.elem_ordrs(1,ne0).ne.generation) continue = .false.
    enddo
    
  end function element_top_of_branch
  
! #################################################################
end module particle_transport
