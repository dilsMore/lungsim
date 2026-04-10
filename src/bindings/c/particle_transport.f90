module particle_transport_c
  implicit none
  private

contains

!!!###################################################################################

  subroutine solve_particles_decoupled_c(initial_concentration,inlet_concentration,particle_size) &
    bind(C, name="solve_particles_decoupled_c")
    use arrays, only: dp
    use particle_transport, only: solve_particles_decoupled
    implicit none
    
    real(dp),intent(in) :: initial_concentration, inlet_concentration, particle_size

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_solve_particles_decoupled(initial_concentration,inlet_concentration,particle_size)
#else
    call solve_particles_decoupled(initial_concentration, inlet_concentration,particle_size)
#endif

  end subroutine solve_particles_decoupled_c

!  subroutine write_airway_c(ne_field, EXELEMFILE, filename_len, group_name, group_name_len, field_name, field_name_len) &
!    bind(C, name="write_airway_c")
!
!    use iso_c_binding, only: c_ptr
!    use utils_c, only: strncpy
!    use particle_transport, only: write_airway
!    !use exports, only: export_1d_elem_field
!    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
!    implicit none
!    integer,intent(in) :: ne_field, filename_len, group_name_len, field_name_len
!    type(c_ptr), value, intent(in) :: EXELEMFILE, group_name, field_name
!    character(len=MAX_FILENAME_LEN) :: filename_f
!    character(len=MAX_STRING_LEN) :: group_name_f, field_name_f
!
!    call strncpy(filename_f, EXELEMFILE, filename_len)
!    call strncpy(group_name_f, group_name, group_name_len)
!    call strncpy(field_name_f, field_name, field_name_len)
!
!#if defined _WIN32 && defined __INTEL_COMPILER
!    call so_write_airway(ne_field, filename_f, group_name_f, field_name_f)
!#else
!    call write_airway(ne_field, filename_f, group_name_f, field_name_f)
!#endif
!
!  end subroutine write_airway_c
!
!
!  subroutine write_terminal_c(EXNODEFILE, filename_len, name, name_len) bind(C, name="write_terminal_c")
!
!    use iso_c_binding, only: c_ptr
!    use utils_c, only: strncpy
!    use particle_transport, only: write_terminal
!    use other_consts, only: MAX_STRING_LEN, MAX_FILENAME_LEN
!    implicit none
!    integer,intent(in) :: filename_len, name_len
!    type(c_ptr), value, intent(in) :: EXNODEFILE, name
!    character(len=MAX_FILENAME_LEN) :: filename_f
!    character(len=MAX_STRING_LEN) :: name_f
!
!    call strncpy(filename_f, EXNODEFILE, filename_len)
!    call strncpy(name_f, name, name_len)
!
!#if defined _WIN32 && defined __INTEL_COMPILER
!    call so_write_terminal(filename_f, name_f)
!#else
!    call write_terminal(filename_f, name_f)
!#endif
!
!  end subroutine write_terminal_c
!###################################################################################
end module particle_transport_c
