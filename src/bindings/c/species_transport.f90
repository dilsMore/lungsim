module species_transport_c
  implicit none
  private

contains


  !!!###################################################################################

  subroutine initialise_exchange_c() bind(C, name="initialise_exchange_c")

    use species_transport, only: initialise_exchange
    implicit none

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_initialise_exchange
#else
    call initialise_exchange
#endif

  end subroutine initialise_exchange_c

end module species_transport_c
