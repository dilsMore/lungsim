module species_transport
!*Brief Description:* This module contains all the subroutines common
!to species transport models, this includes gas exchange, gas mixing,
!and particle transport models
!*LICENSE:*
!
!
!
!*Full Description:*
!More info on what the module does if necessary
  !
  use arrays
  use diagnostics
  use gas_exchange
  use indices
  use other_consts
  
  implicit none

  !Module parameters

  !Module types

  !Module variables

  !Interfaces
  private 
  public initialise_transport
  public initialise_exchange
  public assemble_transport_matrix
  public reduce_transport_matrix
  public set_elem_volume
  public calc_mass

contains
!
!##############################################################################
!

 !
!##############################################################################
!
! subroutine solve_transport()
! !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT" :: INITIALISE_TRANSPORT
!   use indices
!   use arrays, only: dp
!   use gas_exchange, only: steadystate_gasexchange
!   use particle_transport, only: controller_particlesolve
!   use diagnostics, only: enter_exit
!
!   !local variables
!   real(dp) c_art_o2, c_ven_o2,p_art_co2,p_art_o2, p_ven_co2,p_ven_o2!
!
!   character(len=60) :: sub_name
!
!   sub_name = 'solve_transport'
!   call enter_exit(sub_name,1)
!
!
!   select case (model_type)
!     case ('gas_exchange')
!       print *, 'Nothing implemented'
!       !Note that as V, Q are prerequisites something needs to be added here that checks
!       !these have been read in and if not sets up linear gradient based on some default parameters
!     case ('gas_mix')
!       print *, 'Nothing implemented'
!       !Note that as V is prerequisites something needs to be added here that checks
!       !these have been read in and if not sets up linear gradient based on some default parameters
!     case ('gas_transfer')
!       print *, 'Calling gas transfer model '
!       p_art_co2=40.0_dp
!       p_ven_co2=45.0_dp
!       p_art_o2=100.0_dp
!       p_ven_o2=40.0_dp
!       call steadystate_gasexchange(c_art_o2,c_ven_o2,&
!       p_art_co2,p_art_o2,149.0_dp,p_ven_co2,p_ven_o2,0.03_dp,&
!       0.8_dp*(260.0_dp*1.0e+3_dp/60.0_dp),260.0_dp*1.0e+3_dp/60.0_dp )
!    case('particle_transport')
!       call controller_particlesolve('.')
!       
!     case DEFAULT
!        print*, 'The problem does not exist, exiting'
!        stop
!
!    end select
!   call enter_exit(sub_name,2)
! end subroutine solve_transport


!
!###########################################################################################
!
 subroutine initialise_transport(initial_concentration,inlet_concentration,tp)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT" :: INITIALISE_TRANSPORT
   use indices
   use arrays, only: dp,gasex_field,num_units,node_field,transport_parameters
   use diagnostics, only: enter_exit
   
   real(dp), intent(in) :: initial_concentration,inlet_concentration
   type(transport_parameters) :: tp
   character(len=60) :: sub_name
   sub_name = 'initialise_transport'
   call enter_exit(sub_name,1)
   
   write(*,*) 'Allocating memory and initialising arrays for species transport problems'
   select case (model_type)
     case ('gas_mix')
       print *, 'You are solving a gas mixing model'
       !Note that as V is prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
    case('particle_transport')
        print *, 'You are solving a particle transport model'
        call intial_transport(initial_concentration,inlet_concentration,tp)

    case DEFAULT
        print*, 'The problem does not exist, exiting'
        stop
    end select
   
   call enter_exit(sub_name,2)
 end subroutine initialise_transport
 
 
 !
!###########################################################################################
!
 subroutine initialise_exchange()
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_EXCHANGE" :: INITIALISE_EXCHANGE
   use indices
   use arrays, only: dp,gasex_field,num_units,node_field
   use diagnostics, only: enter_exit

   character(len=60) :: sub_name
   sub_name = 'initialise_exchange'
   call enter_exit(sub_name,1)
   
   write(*,*) 'Allocating memory and initialising arrays for species transport problems'
   select case (model_type)
     case ('gas_exchange')
       print *, 'You are solving a gas exchange model'
       !Note that as V, Q are prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
     case ('gas_transfer')
       print *, 'You are solving a gas transfer model'
       !Note that as V, Q are prerequisites something needs to be added here that checks
       !these have been read in and if not sets up linear gradient based on some default parameters
       !note a linear q gradient should  be set up to scale for shunt fraction automatically
       call initial_gasexchange(149.0_dp)
    case DEFAULT
        print*, 'The problem does not exist, exiting'
        stop
    end select
   
   call enter_exit(sub_name,2)
 end subroutine initialise_exchange
 
 !!!#########################################################################

  subroutine assemble_transport_matrix(diffusion_coeff)

    use indices
    use arrays,only: dp,elem_nodes,num_elems,sparsity_col,reduced_col,sparsity_row,&
      reduced_row, global_K, global_M, global_AA, global_BB, global_R,nonzeros_unreduced
    use geometry, only: volume_of_mesh
    use diagnostics, only: enter_exit
    implicit none

    real(dp),intent(in) :: diffusion_coeff

    integer :: i,j,ncol,ne,nentry,nrow
    real(dp) :: elem_K(2,2),elem_M(2,2),elem_R(2)
    logical :: found
    character(len=60) :: sub_name

!!!................................................................

    sub_name = 'assemble_transport_matrix'
    call enter_exit(sub_name,1)

    global_K(1:nonzeros_unreduced) = 0.0_dp
    global_M(1:nonzeros_unreduced) = 0.0_dp

    do ne=1,num_elems
       select case (model_type)
       case ('gas_mix')
          call element_gasmix(ne,elem_K,elem_M,elem_R,diffusion_coeff)
       case ('particle_transport')
          call element_particles(ne,elem_K,elem_M,elem_R)
        end select
       do i=1,2
          nrow = elem_nodes(i,ne)
          do j=1,2
             ncol = elem_nodes(j,ne)
             found=.false.
             nentry = sparsity_row(nrow) ! start check at start of row
             do while (.not.found)
                if(ncol.eq.sparsity_col(nentry))then
                   found = .true.
                else
                   nentry = nentry+1
                endif
             enddo
             global_K(nentry) = global_K(nentry) + elem_K(i,j)
             global_M(nentry) = global_M(nentry) + elem_M(i,j)
          enddo !j
       enddo !i
    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine assemble_transport_matrix
  
!!!########################################################################

  subroutine intial_transport(initial_concentration,inlet_concentration,tp)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIAL_TRANSPORT" :: INITIAL_TRANSPORT

    use arrays,only: dp,node_field,num_nodes,sparsity_col,reduced_col,sparsity_row,&
      reduced_row, global_K, global_M, global_AA, global_BB, global_R,transport_parameters
    use indices,only: nj_conc1,nu_conc1
    use diagnostics, only: enter_exit
    implicit none
    type(transport_parameters), intent(out) :: tp
    real(dp),intent(in) :: initial_concentration,inlet_concentration
    
    real(dp):: initial_mass

    character(len=60):: sub_name

    ! #########################################################################

    sub_name = 'initial_transport'
    call enter_exit(sub_name,1)
    
    tp%initial_concentration = initial_concentration
    node_field(nj_conc1,1:num_nodes) = tp%initial_concentration(1)

    ! initialise the 'ideal mass' to the mass of gas initially in model
    call calc_mass(nj_conc1,nu_conc1,tp%ideal_mass)
    
    tp%inlet_concentration = inlet_concentration
    node_field(nj_conc1,1) = tp%inlet_concentration(1)
    tp%total_volume_change = 0.0_dp ! records the volume change from FRC

! allocate the arrays for solving
    if(.not.allocated(sparsity_col)) allocate(sparsity_col(1+3*(num_nodes-1)))
    if(.not.allocated(reduced_col))  allocate(reduced_col(1+3*(num_nodes-1)))
    if(.not.allocated(sparsity_row)) allocate(sparsity_row(num_nodes+1))
    if(.not.allocated(reduced_row))  allocate(reduced_row(num_nodes+1))
    if(.not.allocated(global_K))     allocate(global_K(1+3*(num_nodes-1)))
    if(.not.allocated(global_M))     allocate(global_M(1+3*(num_nodes-1)))
    if(.not.allocated(global_AA))    allocate(global_AA(1+3*(num_nodes-1)))
    if(.not.allocated(global_BB))    allocate(global_BB(num_nodes))
    if(.not.allocated(global_R))     allocate(global_R(num_nodes))

    ! calculate the sparsity pattern for unreduced system
    call sparse_transport

    call enter_exit(sub_name,2)

  end subroutine intial_transport
  
  !!!##########################################################################

  subroutine reduce_transport_matrix(MatrixSize,NonZeros,noffset_entry,noffset_row,&
       inspiration)
    use arrays,only: num_nodes,NonZeros_unreduced,reduced_row,reduced_col, &
      sparsity_row,sparsity_col
    implicit none

    integer :: MatrixSize,NonZeros,&
         noffset_entry,noffset_row
    logical,intent(in) :: inspiration

    integer :: i

    if(inspiration)then !remove first row and column (note: also for breath-hold)

       do i=1,num_nodes ! one more than # of rows
          reduced_row(i) = sparsity_row(i+1)-3
       enddo
       NonZeros = NonZeros_unreduced - 3
       do i=1,NonZeros
          reduced_col(i) = sparsity_col(i+3)-1
       enddo
       reduced_row(1)=1
       MatrixSize = num_nodes - 1
       noffset_entry = 3
       noffset_row = 1

    else !expiration

       do i=1,num_nodes+1
          reduced_row(i) = sparsity_row(i)
       enddo
       NonZeros = NonZeros_unreduced
       do i=1,NonZeros
          reduced_col(i) = sparsity_col(i)
       enddo
       reduced_row(1)=1
       MatrixSize = num_nodes
       noffset_entry = 0
       noffset_row = 0

    endif

  end subroutine reduce_transport_matrix
  
  
!
!##############################################################################
!
 subroutine initial_gasexchange(initial_concentration,surface_area,V_cap)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIAL_GASEXCHANGE" :: INITIAL_GASEXCHANGE
   use indices
   use arrays, only: dp,elem_units_below,gasex_field,node_field,num_nodes,&
         num_units,unit_field
   use diagnostics, only: enter_exit
   

   !local variables
   real(dp),intent(in) :: initial_concentration
   real(dp), optional ::  surface_area,V_cap
   real(dp),parameter :: o2molvol = 25.44e+3_dp ! mm^3/mmol (converted from 22.41e3 at STP using V2=T2*V1/T1)


    integer :: nunit
    real(dp) :: Vcap_unit
    real(dp),parameter :: p_water = 47.0_dp
    real(dp),parameter :: press_atm=760.0_dp !atmospheric pressure, mmHg

   character(len=60) :: sub_name

   sub_name = 'initial_gasexchange'
   call enter_exit(sub_name,1)

!!! allocate memory for the gasex_field array, if not already allocated
    if(.not.allocated(gasex_field)) allocate(gasex_field(num_gx,num_units))

!!! initialiase nj_conc2 (for CO2 concentration); currently hardcoded to 40 mmHg
    node_field(nj_conc2,1:num_nodes) = 40.0_dp/(o2molvol*(press_atm-p_water))
    write(*,'('' Initialising Palv_CO2 to 40 mmHg'')')

!!! initialise the gas exchange field for o2 partial pressures
    gasex_field(ng_p_alv_o2,1:num_units) = initial_concentration* &
         o2molvol*(press_atm-p_water)
    gasex_field(ng_p_cap_o2,1:num_units) = initial_concentration*&
         o2molvol*(press_atm-p_water)

    gasex_field(ng_p_alv_co2,1:num_units) = 40.0_dp ! mmHg; should make this user defined
    gasex_field(ng_p_ven_o2,1:num_units) = 40.0_dp ! mmHg; should make this user defined

    unit_field(nu_conc1,1:num_units) = gasex_field(ng_p_alv_o2,1:num_units)/&
         (o2molvol*(press_atm-p_water)) ! from mmHg to mmol/mm^3
    unit_field(nu_conc2,1:num_units) = gasex_field(ng_p_alv_co2,1:num_units)/&
         (o2molvol*(press_atm-p_water)) ! from mmHg to mmol/mm^3

!!! initialise the gas exchange field for co2 partial pressures
    gasex_field(ng_p_alv_co2,1:num_units) = 40.0_dp ! mmHg; should make this user defined
    gasex_field(ng_p_cap_co2,1:num_units) = 40.0_dp ! mmHg; should make this user defined
    gasex_field(ng_p_ven_co2,1:num_units) = 45.0_dp ! mmHg; should make this user defined
    if(present(surface_area))then
      !!! initialise the time blood has been in capillaries
      gasex_field(ng_time,1:num_units) = 0.0_dp
      !!! capillary volume per gas exchange unit = transit time * flow
      ! elem_units_below is the EFFECTIVE number of units, so this is correct
      !Note that these are calculated on a per unit basis in the perfusion model so can be read in for future iterations
      Vcap_unit = V_cap/elem_units_below(1) ! the capillary volume per gas exchange unit
      gasex_field(ng_Vc,1:num_units) = Vcap_unit
      gasex_field(ng_sa,1:num_units) = surface_area/elem_units_below(1)

!!! transit time through the gas exchange unit = capillary volume/flow
      forall (nunit=1:num_units) gasex_field(ng_tt,nunit) = &
           Vcap_unit/unit_field(nu_perf,nunit)
    endif

   call enter_exit(sub_name,2)
 end subroutine initial_gasexchange
 
 !!!################################################################################

  subroutine calc_mass(nj,nu_field,gas_mass)
    use arrays,only: dp,elem_cnct,elem_nodes,elem_symmetry,&
         node_field,num_elems,num_nodes,num_units,elem_field,units,unit_field
    use indices,only: ne_vol,nu_vol
    use diagnostics, only: enter_exit
    implicit none

    integer,intent(in) :: nj,nu_field
    real(dp) :: gas_mass
    !     Local Variables
    integer :: ne,ne0,np1,np2,nunit
    real(dp) :: average_conc
    real(dp),allocatable :: tree_mass(:)
    character(len=60):: sub_name

    !...........................................................................

    sub_name = 'calc_mass'
    call enter_exit(sub_name,1)

    if(.not.allocated(tree_mass)) allocate(tree_mass(num_nodes))

    ! initialise to the mass in each element
    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)
       average_conc = (node_field(nj,np1)+node_field(nj,np2))/2.0_dp
       tree_mass(ne) = average_conc*elem_field(ne_vol,ne)
    enddo

    ! add the mass in each elastic unit to terminal elements
    do nunit=1,num_units
       ne=units(nunit)
       tree_mass(ne) = tree_mass(ne) + &
            unit_field(nu_vol,nunit)*unit_field(nu_field,nunit)
    enddo

    ! sum mass recursively up the tree
    do ne=num_elems,2,-1 ! not for the stem branch; parent = 0
       ne0=elem_cnct(-1,1,ne)
       tree_mass(ne0) = tree_mass(ne0) + dble(elem_symmetry(ne))*tree_mass(ne)
    enddo !noelem

    gas_mass = tree_mass(1)

    deallocate(tree_mass)

    call enter_exit(sub_name,2)

  end subroutine calc_mass
  
  
!!!########################################################################

  subroutine sparse_transport
    use arrays,only: elem_cnct,elem_nodes,num_elems,&
      NonZeros_unreduced,reduced_row,reduced_col, &
      sparsity_row,sparsity_col

    implicit none

    integer :: n_unreduced,i,ncol,ne,ne2,np1,np2,nrow

    sparsity_row(1) = 1
    n_unreduced = 1

    do ne=1,num_elems ! note using local numbering
       if(elem_cnct(-1,0,ne).eq.0)then !at the inlet
          np1=elem_nodes(1,ne) ! start node
          nrow=np1
          do i=1,2
             np2=elem_nodes(i,ne)
             ncol=np2
             sparsity_col(n_unreduced)=ncol
             n_unreduced=n_unreduced+1
          enddo
          sparsity_row(nrow+1)=n_unreduced
       endif

       np1=elem_nodes(2,ne) !end node
       nrow=np1
       do i=1,2
          np2=elem_nodes(i,ne)
          ncol=np2
          sparsity_col(n_unreduced)=ncol
          n_unreduced=n_unreduced+1
       enddo
       do i=1,elem_cnct(1,0,ne) ! for each child branch
          ne2=elem_cnct(1,i,ne)
          np2=elem_nodes(2,ne2)
          ncol=np2
          sparsity_col(n_unreduced)=ncol
          n_unreduced=n_unreduced+1
       enddo
       sparsity_row(nrow+1)=n_unreduced
    enddo !noelem

    NonZeros_unreduced = n_unreduced - 1
    
  end subroutine sparse_transport
  
 !########################## 
  subroutine set_elem_volume
    use arrays,only: dp,elem_field,elem_nodes,node_xyz,num_elems
    use diagnostics, only: enter_exit
    use indices,only: ne_Vdot,ne_length, &
         ne_radius,ne_resist,ne_t_resist,ne_vol,ne_a_A
    use other_consts
    implicit none

    ! Local variables
    integer :: ne,np1,np2

    character(len=60) :: sub_name

    ! ###########################################################################

    sub_name = 'set_elem_volume'
    call enter_exit(sub_name,1)

    do ne=1,num_elems
       np1=elem_nodes(1,ne)
       np2=elem_nodes(2,ne)

       ! element length
       elem_field(ne_length,ne) = DSQRT((node_xyz(1,np2) - &
            node_xyz(1,np1))**2 + (node_xyz(2,np2) - &
            node_xyz(2,np1))**2 + (node_xyz(3,np2) - &
            node_xyz(3,np1))**2)

       ! element volume
       elem_field(ne_vol,ne) = PI * elem_field(ne_radius,ne)**2 * &
            elem_field(ne_length,ne)
            
       elem_field(ne_a_A,ne) = 1.0_dp ! set default for ratio a/A


    enddo !noelem

    call enter_exit(sub_name,2)

  end subroutine set_elem_volume
  
  
  subroutine element_gasmix(ne,elem_K,elem_M,elem_R,diffusion_coeff)
    use arrays,only: dp,elem_field,elem_symmetry
    use indices,only: ne_a_A,ne_length,ne_radius
    use other_consts
    implicit none

    integer,intent(in) :: ne
    real(dp) :: elem_K(2,2),elem_M(2,2),elem_R(2)
    real(dp),intent(in) :: diffusion_coeff

    !Local variables
    real(dp) :: a_A_ratio,inner_area,length,outer_area,radius

    radius = elem_field(ne_radius,ne)
    length = elem_field(ne_length,ne)
    a_A_ratio = elem_field(ne_a_A,ne)
    outer_area=PI*radius**2
    inner_area=outer_area*a_A_ratio

    elem_M(1,1) = outer_area*length/3.0_dp*DBLE(elem_symmetry(ne))
    elem_M(1,2) = outer_area*length/3.0_dp/2.0_dp*DBLE(elem_symmetry(ne))
    elem_M(2,1) = outer_area*length/3.0_dp/2.0_dp
    elem_M(2,2) = outer_area*length/3.0_dp

    elem_K(1,1) = (inner_area*diffusion_coeff/length)*DBLE(elem_symmetry(ne))
    elem_K(1,2) = (-inner_area*diffusion_coeff/length)*DBLE(elem_symmetry(ne))
    elem_K(2,1) = -inner_area*diffusion_coeff/length
    elem_K(2,2) = inner_area*diffusion_coeff/length

    elem_R = 0.0_dp

  end subroutine element_gasmix
  
    !!!########################################################################
!!! original code developed by Falko Schmidt (2011). Adapted by Merryn Tawhai
  subroutine element_particles(ne,elem_K,elem_M,elem_R) 

    use other_consts
    use arrays
    use indices
    use diagnostics

    implicit none

    type(particle_parameters) :: part_param
    integer,intent(in) :: ne
    real(dp) :: elem_K(2,2),elem_M(2,2),elem_R(2)
!!!    real(dp),intent(in) :: diffusion_coeff

    !Local variables
    real(dp) :: a_A_ratio,inner_area,kappa,length,outer_area,radius

    radius = elem_field(ne_radius,ne)
    length = elem_field(ne_length,ne)
    a_A_ratio = elem_field(ne_a_A,ne)
    outer_area=PI*radius**2
    inner_area=outer_area*a_A_ratio

    ! TJ - 2022 NOV 02: change ne_flow to be ne_Vdot
    if(elem_field(ne_Vdot,1).gt.0.0_dp)then ! inhalation
       ! apparent diffusion acc.to. Lee2001 exhalation
       kappa = 0.26_dp*abs(elem_field(ne_Vdot,ne))*2.0_dp/pi/radius
       kappa = kappa/6.0_dp/radius*elem_field(ne_length,ne)
       ! using this kappa is f(l) instead f(d) - see results validation Gomes1993
    else ! exhalation
       kappa = 0.26_dp*abs(elem_field(ne_Vdot,ne))*2.0_dp/pi/radius ! apparent diffusion acc.to. Lee2001 exhalation
    endif

    elem_M(1,1) = outer_area*length/3.0_dp*DBLE(elem_symmetry(ne))
    elem_M(1,2) = outer_area*length/3.0_dp/2.0_dp*DBLE(elem_symmetry(ne))
    elem_M(2,1) = outer_area*length/3.0_dp/2.0_dp
    elem_M(2,2) = outer_area*length/3.0_dp

    elem_K(1,1) = (inner_area*(part_param%diffu+kappa)/length)*DBLE(elem_symmetry(ne))
    elem_K(1,2) = (-inner_area*(part_param%diffu+kappa)/length)*DBLE(elem_symmetry(ne))
    elem_K(2,1) = -inner_area*(part_param%diffu+kappa)/length
    elem_K(2,2) = inner_area*(part_param%diffu+kappa)/length

    elem_R = 0.0_dp

  end subroutine element_particles

end module species_transport
