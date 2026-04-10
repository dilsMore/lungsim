module arrays
!*Brief Description:* This module defines arrays.
!
!*LICENSE:*
!
!
!*Contributor(s):* Merryn Tawhai, Alys Clark
!
!*Full Description:*
!
!This module defines arrays

  use precision
  
  implicit none

  integer :: num_elems,num_elems_2d,num_nodes,num_data,num_nodes_2d,num_units,num_lines_2d,maxgen

  integer,allocatable :: nodes(:) !allocated in define_node_geometry
  integer,allocatable :: nodes_2d(:) !allocated in define_node_geometry_2d
  integer,allocatable :: node_versn_2d(:) !allocated in define_node_geometry_2d
  integer,allocatable :: elems(:) !allocated in define_1d_elements
  integer,allocatable :: lines_2d(:)
  integer,allocatable :: parentlist(:)
  integer,allocatable :: line_versn_2d(:,:,:)
  integer,allocatable :: lines_in_elem(:,:)
  integer,allocatable :: nodes_in_line(:,:,:)
  integer,allocatable :: elems_2d(:) !allocated in define_elem_geometry_2d
  integer,allocatable :: elem_cnct(:,:,:)  !NXI(-ni:ni,1,ne)
  integer,allocatable :: elem_cnct_2d(:,:,:)
  integer,allocatable :: elem_nodes(:,:)
  integer,allocatable :: elem_nodes_2d(:,:)
  integer,allocatable :: elem_versn_2d(:,:)
  integer,allocatable :: elem_lines_2d(:,:)
  integer,allocatable :: elem_ordrs(:,:)
  integer,allocatable :: elem_symmetry(:)
  integer,allocatable :: elem_units_below(:)
  integer,allocatable :: elems_at_node(:,:)
  integer,allocatable :: elems_at_node_2d(:,:)
  integer,allocatable :: units(:)

  real(dp),allocatable :: arclength(:,:)
  real(dp),allocatable :: elem_field(:,:) !properties of elements
  real(dp),allocatable :: elem_direction(:,:)
  real(dp),allocatable :: node_xyz(:,:)
  real(dp),allocatable :: data_xyz(:,:)
  real(dp),allocatable :: data_weight(:,:)
  real(dp),allocatable :: node_xyz_2d(:,:,:,:)
  real(dp),allocatable :: gasex_field(:,:) !gasexchange specific fields
  real(dp),allocatable :: unit_field(:,:) !properties of elastic units
  
  !FEM Matrices for a given problem
  integer,allocatable :: sparsity_col(:),reduced_col(:)
  integer,allocatable :: sparsity_row(:),reduced_row(:)
  real(dp),allocatable :: global_K(:),global_M(:),global_AA(:),global_BB(:)
  real(dp),allocatable :: global_R(:)
  integer :: NonZeros_unreduced
  
  !TEMP: ARC SOME SORT OF ACINUS FIELD FOR PARTICLES, SHOULS be a unit_field
  real(dp),allocatable :: part_acinus_field(:,:) !for particle deposition problems

  real(dp),allocatable :: node_field(:,:)
  real(dp),allocatable :: scale_factors_2d(:,:)

  logical,allocatable :: expansile(:)

  type capillary_bf_parameters
    integer :: num_symm_gen=9 !no units
    real(dp) :: total_cap_area=0.63000e02_dp !m
    real(dp) :: Palv=0.0_dp!Pa
    real(dp) :: H0=0.35000e-05_dp !m
    real(dp) :: K_cap=0.12000e02_dp
    real(dp) :: F_cap=0.18000e01_dp
    real(dp) :: F_sheet=0.10400e00_dp
    real(dp) :: sigma_cap=0.43637e03_dp !Pa
    real(dp) :: mu_c=0.19200e-02_dp !Pa.s
    real(dp) :: alpha_a=2.33e-08_dp !m/Pa
    real(dp) :: alpha_v=2.33e-08_dp !m/Pa
    real(dp) :: F_rec=0.64630e00_dp
    real(dp) :: sigma_rec=0.22300e04_dp
    real(dp) :: L_c=0.11880e-02_dp !m
    real(dp) :: Plb_c=0.0_dp !Pa
    real(dp) :: Pub_c=3138.24_dp !Pa
    real(dp) :: Pub_a_v=3138.24_dp !Pa
    real(dp) :: L_art_terminal=0.13000e-03_dp !m
    real(dp) :: L_vein_terminal=0.13000e-03_dp !m
    real(dp) :: R_art_terminal=0.10000e-04_dp !m
    real(dp) :: R_vein_terminal=0.90000e-05!m
  end type capillary_bf_parameters
  
  type transport_parameters
    real(dp) :: ideal_mass
    real(dp) :: total_volume_change
    real(dp) :: inlet_concentration(3)!currently hardcoded to up to three different materials
    real(dp) :: initial_concentration(3)
    !real(dp) :: inlet_mouth_concentration(3) ! currently hardcoded to up to three different materials
  end type transport_parameters
  
  !TEMP: ARC: particle transport parameters
  type particle_parameters
    integer :: num_brths_gm
    real(dp) :: solve_tolerance, initial_volume, diffusion_coeff, gravityx,&
      gravityy,gravityz,pdia,time_inspiration,time_breath_hold,time_expiration,&
      dt_gm, VtotTLC,totacinarLength
    real(dp) :: inlet_flow, Eo, inlet_mouth_concentration ! TJ - To calculate extrathoracic deposition
    real(dp) :: tidal_volume = 1.e+06_dp! tidal volume target, mm^3
    real(dp) :: FRC = 3.36
    real(dp) :: mu = 18.69e-6_dp
    real(dp) :: prho =  1.0e-3_dp             ! ! density of particles [g/mm^3]
    real(dp) :: lambda = 7.022e-5_dp  ! [mm] mean free path necessary
    real(dp) :: kBoltz = 1.38e-14_dp  ! ! Boltzmann constant [J/K*1d9]=[kg*m^2/s^2/K*1d9]=[g*mm^2/s^2/K]
    real(dp) :: Temperature = 36.0_dp+273.15_dp ! Temperature [K] from rho*R*T
    integer :: out_itr_max = 200      ! max # (outer) iterations using GMRES solver.
    integer :: inr_itr_max = 100      ! max # (inner) iterations using GMRES solver.

    logical :: coupled = .false.
    logical :: last_breath, inspiration
    integer :: n_export
    character(len=200) :: lung_root
    character(len=200) :: results_location
    character(len=20) :: study
    character(len=20) :: subject
    character(len=20) :: protocol
    character(len=100) :: group_name
    real(dp) :: diffu,LacTLC(10),RacTLC(10),VacTLC(9)
    integer :: grav_factor
  end type particle_parameters

  type admittance_param
    character (len=20) :: admittance_type
    character (len=20) :: bc_type
  end type admittance_param
  type, EXTENDS (admittance_param) :: two_parameter
     real(dp) :: admit_P1=1.0_dp
     real(dp) :: admit_P2=1.0_dp
  end type two_parameter
  type, EXTENDS (two_parameter) :: three_parameter
    real(dp) :: admit_P3=1.0_dp
  end type three_parameter
  type, EXTENDS (three_parameter) :: four_parameter
    real(dp) :: admit_P4=1.0_dp
  end type four_parameter
  type,EXTENDS (four_parameter) :: all_admit_param
  end type all_admit_param

  type elasticity_vessels
    character(len=20) ::vessel_type
  end type elasticity_vessels
  type, EXTENDS(elasticity_vessels) :: elasticity_param
    real(dp) :: elasticity_parameters(3)=0.0_dp
  end type elasticity_param

  type fluid_properties
    real(dp) :: blood_viscosity=0.33600e-02_dp !Pa.s
    real(dp) :: blood_density=0.10500e-02_dp !kg/cm3
    real(dp) :: air_viscosity
    real(dp) :: air_density
  end type fluid_properties

! temporary, for debugging:
  real(dp) :: unit_before

  private

  public set_node_field_value, elem_field, num_elems, num_elems_2d, elem_nodes, node_xyz, &
         nodes,nodes_2d, elems, num_nodes, num_nodes_2d, num_data, data_xyz, data_weight, &
         node_xyz_2d, node_versn_2d, units, num_units, unit_field, node_field, dp, &
         elem_cnct, elem_ordrs, elem_direction, elems_at_node, elem_symmetry, expansile, &
         elem_units_below, maxgen,capillary_bf_parameters, zero_tol,loose_tol,gasex_field, &
         num_lines_2d, lines_2d, line_versn_2d, lines_in_elem, nodes_in_line, elems_2d, &
         elem_cnct_2d, elem_nodes_2d, elem_versn_2d, elem_lines_2d, elems_at_node_2d, arclength, &
         scale_factors_2d, parentlist, fluid_properties, elasticity_vessels, admittance_param, &
         elasticity_param, all_admit_param,transport_parameters
  !TEMP ARC particle stuff that is wrong
  public part_acinus_field, particle_parameters
  
  !FEM ARRAYS
  public sparsity_col,reduced_col,sparsity_row,&
      reduced_row, global_K, global_M, global_AA, global_BB, global_R,NonZeros_unreduced

contains
  subroutine set_node_field_value(row, col, value)
  !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_SET_NODE_FIELD_VALUE" :: SET_NODE_FIELD_VALUE
    implicit none

    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

    node_field(row, col) = value

  end subroutine set_node_field_value


end module arrays
