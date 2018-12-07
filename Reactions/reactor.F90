module reactor_module

  use amrex_fort_module, only : amrex_real
  use network
  use eos_module
  use react_type_module
  use actual_reactor_module
#ifdef USE_SDC_FORTRAN
  use actual_sdc_module, only : actual_reactor_init_sdc, actual_reactor_close_sdc, actual_react_sdc
#endif

  implicit none

  logical, save, private :: reactor_initialized = .false.

contains

  !Original DVODE version
  subroutine reactor_init(iE) bind(C, name="reactor_init")

    implicit none
    integer(c_int),  intent(in   ) :: iE

    write(*,*) "Using Original dvode react subroutine"

    !$omp parallel
    call actual_reactor_init(iE)
    !$omp end parallel
    
    reactor_initialized = .true.

  end subroutine reactor_init

#ifdef USE_SUNDIALS3x4x
  ! Call to CVODE: only in Fuego 
  subroutine reactor_init_cvode(imethod, iiter, iJac, iE, iDense) bind(C, name="reactor_init_cvode")

    use, intrinsic :: iso_c_binding

    implicit none
    integer(c_int),  intent(in   ) :: imethod, iiter, iJac, iE, iDense
    write(*,*) "Using improved cvode react subroutine"

    !$omp parallel
    call actual_reactor_init_cvode(imethod, iiter, iJac, iE, iDense)
    !$omp end parallel
    
    reactor_initialized = .true.

  end subroutine reactor_init_cvode
#endif

#ifdef USE_SDC_FORTRAN
  ! Call to SDC: only in Fuego 
  subroutine reactor_init_sdc(nLobato, nsdcite) bind(C, name="reactor_init_sdc")

    use, intrinsic :: iso_c_binding

    implicit none
    integer,  intent(in   ) :: nLobato, nsdcite 

    write(*,*) "Using SDC integrator for chemistry"

    !$omp parallel
    call actual_reactor_init_sdc(nLobato, nsdcite)
    !$omp end parallel
    
    reactor_initialized = .true.

  end subroutine reactor_init_sdc
#endif


  subroutine reactor_close() bind(C, name="reactor_close")

    implicit none

    call actual_reactor_close()
    
    reactor_initialized = .false.

  end subroutine reactor_close


#ifdef USE_SDC_FORTRAN
  subroutine reactor_close_sdc() bind(C, name="reactor_close_sdc")

    implicit none

    call actual_reactor_close_sdc()

    reactor_initialized = .false.

  end subroutine reactor_close_sdc
#endif


  function ok_to_react(state)

    implicit none
    type (react_t),intent(in) :: state
    logical                   :: ok_to_react

    ok_to_react = actual_ok_to_react(state)

  end function ok_to_react

  !Original DVODE version
  function react(react_state_in, react_state_out, dt_react, time)

    use amrex_error_module

    type(react_t),  intent(in    ) :: react_state_in
    type(react_t),  intent(inout ) :: react_state_out
    real(c_double), intent(in    ) :: dt_react, time
    type(reaction_stat_t)          :: react

    if (.not. reactor_initialized) then
       call amrex_error('reactor::react called before initialized')
    endif

    !write(*,*) "Using Original dvode react subroutine \n"

    if ( ok_to_react(react_state_in) ) then

       react = actual_react(react_state_in, react_state_out, dt_react, time)

    else

       react = actual_react_null(react_state_in, react_state_out, dt_react, time)

    endif

  end function react

#ifdef USE_SUNDIALS3x4x
  ! Call to CVODE. Only with Fuego
  function react_cvode(react_state_in, react_state_out, dt_react, time)

    use amrex_error_module
    !use, intrinsic :: iso_c_binding

    type(react_t),  intent(in    ) :: react_state_in
    type(react_t),  intent(inout ) :: react_state_out
    real(c_double), intent(in    ) :: dt_react, time
    type(reaction_stat_t)          :: react_cvode

    if (.not. reactor_initialized) then
       call amrex_error('reactor::react_cvode called before initialized')
    endif

    !write(*,*) "Using improved cvode react subroutine \n"

    if ( ok_to_react(react_state_in) ) then

       react_cvode = actual_react_cvode(react_state_in, react_state_out, dt_react, time)

    else

       react_cvode = actual_react_null(react_state_in, react_state_out, dt_react, time)

    endif

  end function react_cvode
#endif

#ifdef USE_SDC_FORTRAN
  ! Call to SDC. Only with Fuego
  function react_sdc(react_state_in, react_state_out, dt_react, time)

    use amrex_error_module

    type(react_t),  intent(in    ) :: react_state_in
    type(react_t),  intent(inout ) :: react_state_out
    real(c_double), intent(in    ) :: dt_react, time
    type(reaction_stat_t)          :: react_sdc

    if (.not. reactor_initialized) then
       call amrex_error('reactor::react_sdc called before initialized')
    endif

    !write(*,*) "Using sdc integration"

    if ( ok_to_react(react_state_in) ) then

       react_sdc = actual_react_sdc(react_state_in, react_state_out, dt_react, time)

    else

       react_sdc = actual_react_null(react_state_in, react_state_out, dt_react, time)

    endif

  end function react_sdc
#endif

end module reactor_module
