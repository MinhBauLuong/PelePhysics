module actual_reactor_module

  use, intrinsic :: iso_c_binding
  use amrex_fort_module, only : amrex_real
  use network, only: nspec, spec_names
  use react_type_module
  use eos_type_module
  use mod_timers    

  implicit none

  real(amrex_real), private, allocatable :: vodeVec(:),cdot(:),rhoydot_ext(:), ydot_ext(:)
  real(amrex_real), private:: rhoedot_ext, rhoe_init, time_init, time_out, rhohdot_ext, &
                              rhoh_init, hdot_ext, h_init, time_old, pressureInit
  integer,private :: iloc, jloc, kloc, iE, iDense
  type (eos_t) :: eos_state
  !$omp threadprivate(vodeVec,cdot,rhoydot_ext,ydot_ext,rhoedot_ext,rhoe_init,time_init,time_out,rhohdot_ext,rhoh_init,hdot_ext,h_init,time_old,pressureInit,iloc,jloc,kloc,eos_state)
  ! CVODE STUFF
  integer                              :: iwrk
  real(amrex_real)                     :: rwrk
  real(c_double), pointer              :: yvec(:)
  type(c_ptr)                          :: sunvec_y
  type(c_ptr)                          :: CVmem
  real(amrex_real)                     :: rhoInit
  real(c_double), allocatable          :: Jmat_KLU(:,:)
  integer(c_long)                      :: NNZ
  integer, allocatable                 :: colPtrs(:), rowVals(:), Jdata(:)

contains

  !DVODE VERSION
  subroutine actual_reactor_init(iE_in)

    use, intrinsic :: iso_c_binding
    use vode_module, only : vode_init

    integer(c_int),  intent(in   ) :: iE_in
    integer :: neq, verbose, itol, order, maxstep
    real(amrex_real) :: rtol, atol
    logical :: use_ajac, save_ajac, always_new_j, stiff

    neq = nspec + 1
    verbose = 0
    itol = 1
    order = 5
    maxstep = 5000
    use_ajac = .false.
    save_ajac = .false.
    always_new_J = .true.
    stiff = .true.
    rtol = 1.d-10
    atol = 1.d-10

    call vode_init(neq,verbose,itol,rtol,atol,order,&
         maxstep,use_ajac,save_ajac,always_new_j,stiff)

    print *,"Using good ol' dvode"
    iE = iE_in
    if (iE == 1) then
        print *," ->with internal energy (UV cst)"
        allocate(rhoydot_ext(nspec))
    else
        print *," ->with enthalpy (HP cst)"
        allocate(ydot_ext(nspec))
    end if 

    allocate(vodeVec(neq))
    allocate(cdot(nspec))

    call build(eos_state)

  end subroutine actual_reactor_init

  !CVODE VERSION
  subroutine actual_reactor_init_cvode(imethod, iiter, iJac, iE_in, iDense_in)

    use, intrinsic :: iso_c_binding
    use cvode_interface
    use fnvector_serial_mod
    use amrex_error_module

    integer(c_int),  intent(in   ) :: imethod, iiter, iJac, iE_in, iDense_in
    !Local
    integer(c_long)               :: neq, i
    integer(c_int)                :: ierr ! CVODE return status
    integer(c_long), parameter    :: mxsteps = 50000
    real(c_double)                :: time
    real(c_double), allocatable   :: atol(:)
    real(c_double)                :: rtol
    type(c_ptr)                   :: atol_cptr
    integer                       ::  verbose
    integer                       :: nJdata(1)

    CALL t_total%init("Total")
    CALL t_init%init("Initialization")
    CALL t_ck%init("Chemkin calls (outside AJac)")
    CALL t_eos%init("EOS calls")
    CALL t_ReInit%init("REInit calls")
    CALL t_CVODE%init("GLOBAL CVODE calls")
    CALL t_AJac%init( " --> Total Analytical Jac call") 
    CALL t_ckJac%init("     --> Chemkin call portion of total Analytical Jac call") 

    CALL t_total%start
    CALL t_init%start

    print *,"Using cvode"
    iE = iE_in
    if (iE == 1) then
        print *," ->with internal energy (UV cst)"
    else
        print *," ->with enthalpy (HP cst)"
    end if 

    verbose = 0
    neq = nspec + 1

    ! Setup CVODE crap
    allocate(yvec(neq))
    sunvec_y = FN_VMake_Serial(neq, yvec)
    if (.not. c_associated(sunvec_y)) call amrex_abort("actual_reactor: failed in N_VMake_Serial()")

    ! Allocate necessary memory
    CVmem = FCVodeCreate(imethod, iiter)
    if (.not. c_associated(CVmem)) call amrex_abort("actual_reactor: failed in FCVodeCreate()")

    time = 0.0d0
    CALL t_ReInit%start
    ierr = FCVodeInit(CVmem, c_funloc(F_RHS_F), time, sunvec_y)
    if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVodeInit()")
    CALL t_ReInit%stop

    ! Set up tolerances
    allocate(atol(neq))
    do i=1,neq
        atol(i) = 1.0d-10
    end do
    rtol = 1.0d-10
    atol_cptr = FN_VMake_Serial(neq, atol)
    ierr = FCVodeSVtolerances(CVmem, rtol, atol_cptr)
    if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVodeSVtolerances()")

    iDense = iDense_in
    if (iDense == 1) then
        ! Tell CVODE to use a dense linear solver.
        print *,"--> DENSE solver "
        ierr = FCVDense(CVmem, neq)
        if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVDense()")
    else if (iDense == 99) then
        ierr = FCVIter(CVmem, neq, 0)
    else
        ! Get some sort of sparsity pattern to fill NNZ...
        call SPARSITY_INFO(nJdata)
        print *,"--> SPARSE solver -- non zero entries ", nJdata(1), " represents ", nJdata(1)/float(neq * neq) *100.0, "% sparsity pattern"
        ! Tell CVODE to use a KLU linear solver.
        NNZ = nJdata(1)
        ierr = FCVSparse(CVmem, neq, neq, NNZ, CSC_MAT)
        ! Allocate and fill additional data
        allocate(Jdata(NNZ))
        allocate(rowVals(NNZ))
        allocate(colPtrs(neq+1))
        call SPARSITY_PREPROC(rowVals,colPtrs)
        !print *, size(rowVals)
        !print *,"rowVals ", rowVals(:)
        !print *,"colPtrs(neq+1) ", colPtrs(neq+1)
        ! iJac = 1
        allocate(Jmat_KLU(neq,neq))
    end if

    ! set Jacobian routine
    if ((iJac == 1).and.(iDense == 1)) then
        print *,"--> With Analytical J"
        if (iE == 1) then
            ierr = FCVDlsSetJacFn(CVmem, c_funloc(f_jac_cvode))
            if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVDlsSetDenseJacFn()")
        else 
            ierr = FCVDlsSetJacFn(CVmem, c_funloc(f_jac_cvode_HP))
            if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVDlsSetDenseJacFn()")
        end if
    else if ((iDense /= 1).and.(iDense /= 99)) then
        print *,"--> SPARSE solver -- always with Analytical J"
        if (iE == 1) then
            ierr = FCVDlsSetJacFn(CVmem, c_funloc(f_jac_cvode_KLU))
            if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVDlsSetDenseJacFn()")
        else 
            ierr = FCVDlsSetJacFn(CVmem, c_funloc(f_jac_cvode_HP_KLU))
            if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVDlsSetDenseJacFn()")
            !print *,"--> SPARSE solver -- HP solve not yet implemented ..."
            !stop
        end if
    !else if (iDense == 99) then
    !    ierr = FCVSpilsSetJacTimes(CVmem, NULL, NULL);
    else
        print *,"--> Without Analytical J"
    end if

    ! increase the defaultmaxstep to 5000.
    ierr = FCVodeSetMaxNumSteps(CVmem, mxsteps)
    if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVodeSetMaxNumSteps()")

    ! Set max order
    !ierr = FCVodeSetMaxOrd(CVmem, 2)
    !if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVodeSetMaxOrd()")

    !ierr = FCVodeSetStabLimDet(CVmem, 0)

    if (iE == 1) then
        allocate(rhoydot_ext(nspec))
    else
        allocate(ydot_ext(nspec))
    end if 

    allocate(cdot(nspec))

    call build(eos_state)

    deallocate(atol)

    CALL t_init%stop

  end subroutine actual_reactor_init_cvode


  subroutine actual_reactor_close()

    if (allocated(vodeVec)) deallocate(vodeVec)
    if (allocated(cdot)) deallocate(cdot)
    if (allocated(rhoydot_ext)) deallocate(rhoydot_ext)
    if (allocated(ydot_ext)) deallocate(ydot_ext)
    if (associated(yvec)) nullify(yvec)
    if (allocated(Jdata)) deallocate(Jdata)
    if (allocated(rowVals)) deallocate(rowVals)
    if (allocated(colPtrs)) deallocate(colPtrs)

    call destroy(eos_state)
   
  end subroutine actual_reactor_close


  function actual_ok_to_react(state)

    use extern_probin_module, only: react_T_min, react_T_max, react_rho_min, react_rho_max

    implicit none

    type (react_t),intent(in) :: state
    logical                   :: actual_ok_to_react
    real(amrex_real)           :: rho

    actual_ok_to_react = .true.

    rho = sum(state % rhoY)
    if (state % T   < react_T_min   .or. state % T   > react_T_max .or. &
        rho         < react_rho_min .or. rho         > react_rho_max) then

       actual_ok_to_react = .false.

    endif

  end function actual_ok_to_react


  function actual_react_null(react_state_in, react_state_out, dt_react, time)
    
    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: actual_react_null

    react_state_out = react_state_in
    actual_react_null % cost_value = 0.d0
    actual_react_null % reactions_succesful = .true.

  end function actual_react_null


  ! Original DVODE version
  function actual_react(react_state_in, react_state_out, dt_react, time)
    
    use amrex_error_module
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, always_new_j, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar
    use chemistry_module, only : molecular_weight
    use eos_module

    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(amrex_real), intent(in   ) :: dt_react, time
    type(reaction_stat_t)          :: actual_react

    external dvode

    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail, neq
    real(amrex_real) :: vodeTime, vodeEndTime, rhoInv, press_recalc 
    real(amrex_real) :: Y_div_W(nspec)


    eos_state % rho               = sum(react_state_in % rhoY(:))
    rhoInv                        = 1.d0 / eos_state % rho
    eos_state % T                 = react_state_in % T
    eos_state % massfrac(1:nspec) = react_state_in % rhoY(1:nspec) * rhoInv
    ! for cst HP
    eos_state % p                 = react_state_in % p

    if (iE == 1) then
        eos_state % e = react_state_in % e
        call eos_re(eos_state)
    else
        eos_state % h = react_state_in % h
        call eos_ph(eos_state)
    end if

    if (always_new_j) call setfirst(.true.)

    MF          = vode_MF
    vodeTime    = time
    vodeEndTime = time + dt_react
    neq         = nspec + 1
    time_old    = time

    vodeVec(neq)     = eos_state % T

    if (iE == 1) then
        rhoe_init            = eos_state % e  *  eos_state % rho
        rhoedot_ext          = react_state_in % rhoedot_ext
        rhoydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec)
        vodeVec(1:nspec) = react_state_in % rhoY(:)
    else
        h_init            = eos_state % h  
        hdot_ext          = react_state_in % rhohdot_ext / eos_state % rho
        ydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec) / eos_state % rho
        vodeVec(1:nspec) = react_state_in % rhoY(:) / eos_state % rho
    end if

    time_init = time
    iloc      = react_state_in % i
    jloc      = react_state_in % j
    kloc      = react_state_in % k

    ! Vode: istate
    ! in:  1: init, 2: continue (no change), 3: continue (w/changes)
    ! out: 1: nothing was done, 2: success
    !      -1: excessive work, -2: too much accuracy, -3: bad input
    !      -4: repeated step failure, -5: repeated conv failure
    !      -6: EWT became 0
    istate = 1

    call dvode(f_rhs, neq, vodeVec(:), vodeTime, vodeEndTime,&
         itol, rtol, atol, itask, istate, iopt, voderwork, lvoderwork, &
         vodeiwork, lvodeiwork, f_jac, MF, voderpar, vodeipar)

    if (verbose .ge. 1) then
       write(6,*) '......dvode done:'
       write(6,*) ' last successful step size = ',voderwork(11)
       write(6,*) '          next step to try = ',voderwork(12)
       write(6,*) '   integrated time reached = ',voderwork(13)
       write(6,*) '      number of time steps = ',vodeiwork(11)
       write(6,*) '    number of fs(RHS EVAL) = ',vodeiwork(12)
       write(6,*) '              number of Js = ',vodeiwork(13)
       write(6,*) '    method order last used = ',vodeiwork(14)
       write(6,*) '   method order to be used = ',vodeiwork(15)
       write(6,*) '            number of LUDs = ',vodeiwork(19)
       write(6,*) ' number of Newton iterations ',vodeiwork(20)
       write(6,*) ' number of Newton failures = ',vodeiwork(21)
       if (istate.eq.-4 .or. istate.eq.-5) then
          ifail = vodeiwork(16)
          if (ifail .eq. nspec+1) then
             write(6,*) '   T has the largest error'
          else
             write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
          end if
       end if
    end if

    if (istate > 0) then

       actual_react % reactions_succesful = .true.
       actual_react % cost_value = DBLE(vodeiwork(12)) ! number of f evaluations

       if (iE == 1) then
           eos_state % rho               = sum(vodeVec(1:nspec))
           rhoInv                        = 1.d0 / eos_state % rho
           eos_state % massfrac(1:nspec) = vodeVec(1:nspec) * rhoInv
           eos_state % T                 = vodeVec(neq)
           eos_state % e                 = (rhoe_init  +  dt_react*rhoedot_ext) * rhoInv
           call eos_re(eos_state)
           react_state_out % rhoY(:)     = vodeVec(1:nspec) 
           react_state_out % rho         = sum(vodeVec(1:nspec))
       else
           eos_state % massfrac(1:nspec) = vodeVec(1:nspec)
           eos_state % T                 = vodeVec(neq)
           eos_state % h                 = (h_init  +  dt_react*hdot_ext)
           call eos_ph(eos_state)
           react_state_out % rhoY(:)     = vodeVec(1:nspec)*eos_state % rho
           react_state_out % rho         = eos_state % rho 
       end if

       react_state_out % T = eos_state % T
       react_state_out % e = eos_state % e
       react_state_out % h = eos_state % h
       react_state_out % p = eos_state % p

       !Y_div_W(:)   = eos_state % massfrac(:) / molecular_weight(:)
       !press_recalc = eos_state % rho * eos_state % T * 8.31451e+07 * sum(Y_div_W(:))
       !write(*,*) "e,h,p,rho,prec ? ", eos_state % e, eos_state % h, eos_state % p, react_state_out % rho, press_recalc


    else

       actual_react % reactions_succesful = .false.

       print *,'vode failed at',react_state_in % i,react_state_in % j,react_state_in % k
       print *,'input state:'
       print *,'T',react_state_in%T
       if (iE == 1) then
           print *,'e',react_state_in%e
       else
           print *,'h',react_state_in%h
       end if
       print *,'rho',eos_state%rho
       print *,'rhoY',react_state_in%rhoY
       if (iE == 1) then
           print *,'rhoe forcing',react_state_in%rhoedot_ext
       else
           print *,'rhoh forcing',react_state_in%rhohdot_ext
       end if
       print *,'rhoY forcing',react_state_in%rhoydot_ext(1:nspec)

       write(6,*) '......dvode data:'
       write(6,*) ' last successful step size = ',voderwork(11)
       write(6,*) '          next step to try = ',voderwork(12)
       write(6,*) '   integrated time reached = ',voderwork(13)
       write(6,*) '      number of time steps = ',vodeiwork(11)
       write(6,*) '              number of fs = ',vodeiwork(12)
       write(6,*) '              number of Js = ',vodeiwork(13)
       write(6,*) '    method order last used = ',vodeiwork(14)
       write(6,*) '   method order to be used = ',vodeiwork(15)
       write(6,*) '            number of LUDs = ',vodeiwork(19)
       write(6,*) ' number of Newton iterations ',vodeiwork(20)
       write(6,*) ' number of Newton failures = ',vodeiwork(21)
       if (istate.eq.-4 .or. istate.eq.-5) then
          ifail = vodeiwork(16)
          if (ifail .eq. nspec+1) then
             write(6,*) '   T has the largest error'
          else
             write(6,*) '   spec with largest error = ', trim(spec_names(ifail))
          end if
       end if

       print *,'Final T',vodeVec(neq)
       print *,'Final rhoY',vodeVec(1:nspec)

       call amrex_error('vode failed')

    end if

  end function actual_react


  subroutine f_rhs(neq, time, y, ydot, rpar, ipar)

    use chemistry_module, only : molecular_weight
    use eos_module

    integer,         intent(in)   :: neq, ipar(*)
    real(amrex_real), intent(in)  :: y(neq), time, rpar(*)
    real(amrex_real), intent(out) :: ydot(neq)
    integer          :: n
    real(amrex_real) :: rhoInv, dPdt


    if (iE == 1) then
        eos_state % rho               = sum(y(1:nspec))
        rhoInv                        = 1.d0 / eos_state % rho
        eos_state % massfrac(1:nspec) = y(1:nspec) * rhoInv
        eos_state % T                 = y(neq) ! guess
        eos_state % e = (rhoe_init + (time - time_init) * rhoedot_ext) * rhoInv
        call eos_re(eos_state)
        call eos_get_activity(eos_state)
    else
        eos_state % massfrac(1:nspec) = y(1:nspec) 
        eos_state % T                 = y(neq) ! guess
        eos_state % h = h_init + (time - time_init) * hdot_ext
        call eos_ph(eos_state)
        call eos_get_activity_h(eos_state)
    end if

    call ckwc(eos_state % T, eos_state % Acti, iwrk, rwrk, cdot)

    if (iE == 1) then
        ydot(neq)    = rhoedot_ext 
        do n=1,nspec
           ydot(n)   = cdot(n) * molecular_weight(n) + rhoYdot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%ei(n)*ydot(n)
        end do
        ydot(neq)    = ydot(neq)/(eos_state%rho * eos_state%cv)
    else
        ydot(neq)    = hdot_ext
        do n=1,nspec
           ydot(n)   = cdot(n) * molecular_weight(n) / eos_state%rho + ydot_ext(n)
           ydot(neq) = ydot(neq) - eos_state%hi(n)*ydot(n) 
        end do
        ydot(neq)    = ydot(neq)/eos_state%cp
    end if

  end subroutine f_rhs


  subroutine f_jac(neq, npt, y, t, pd)
    use amrex_error_module

    integer,        intent(in)  :: neq, npt
    real(amrex_real),intent(in)  :: y(neq,npt), t
    real(amrex_real),intent(out) :: pd(neq,neq)

    call amrex_error('DVODE version: Analytic Jacobian not yet implemented')

  end subroutine f_jac


  !CVODE VERSION
  function actual_react_cvode(react_state_in, react_state_out, dt_react, time)
    
    use amrex_error_module
    use fnvector_serial_mod 
    use cvode_interface
    use eos_module
    use, intrinsic :: iso_c_binding

    type(react_t),   intent(in   ) :: react_state_in
    type(react_t),   intent(inout) :: react_state_out
    real(c_double),  intent(in   ) :: time, dt_react
    type(reaction_stat_t)          :: actual_react_cvode

    integer(c_long) :: neq
    integer(c_int)  :: i
    integer(c_int)  :: ierr ! CVODE return status
    real(amrex_real)              :: rhoInv

    integer            :: verbose
    integer(c_long)    :: nfevals,nfevals_jac,nlinsetups,njevals,nstp
    integer(c_int)     :: qlast,qcur 
    real(c_double)     :: hlast,hcur,tcur


    ! Fill the state to advance
    CALL t_eos%start          
    eos_state % rho               = sum(react_state_in % rhoY(:))
    rhoInv                        = 1.d0 / eos_state % rho
    eos_state % T                 = react_state_in % T
    eos_state % massfrac(1:nspec) = react_state_in % rhoY(1:nspec) * rhoInv
    ! for cst HP
    eos_state % p                 = react_state_in % p
    pressureInit                  = react_state_in % p

    if (iE == 1) then
        eos_state % e = react_state_in % e
        call eos_re(eos_state)
    else
        eos_state % h = react_state_in % h
        call eos_ph(eos_state)
    end if
    CALL t_eos%stop          

    time_init = time
    time_out  = time + dt_react
    neq       = nspec + 1

    yvec(neq) = eos_state % T

    if (iE == 1) then
        rhoe_init            = eos_state % e  *  eos_state % rho
        rhoedot_ext          = react_state_in % rhoedot_ext
        rhoydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec)
        yvec(1:nspec)        = react_state_in % rhoY(1:nspec)
    else
        h_init               = eos_state % h  
        hdot_ext             = react_state_in % rhohdot_ext / eos_state % rho
        ydot_ext(1:nspec)    = react_state_in % rhoydot_ext(1:nspec) / eos_state % rho
        yvec(1:nspec)        = react_state_in % rhoY(:) / eos_state % rho
    end if

    iloc = react_state_in % i
    jloc = react_state_in % j
    kloc = react_state_in % k

    ! Do we really need this ?
    CALL t_ReInit%start
    ierr = FCVodeReInit(CVmem, time_init, sunvec_y)
    if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVodeReInit()")  
    CALL t_ReInit%stop

    ! cvode call
    CALL t_CVODE%start
    ierr = FCVode(CVmem, time_out, sunvec_y, time_init, CV_NORMAL)
    if (ierr /= 0) call amrex_abort("actual_reactor: failed in FCVode()")
    CALL t_CVODE%stop

    CALL t_eos%start          
    if (iE == 1) then
        eos_state % rho = sum(yvec(1:nspec))
        rhoInv          = 1.d0 / eos_state % rho
        eos_state % T   = yvec(neq)
        eos_state % massfrac(1:nspec) = yvec(1:nspec) * rhoInv
        eos_state % e                 = (rhoe_init  +  dt_react*rhoedot_ext) * rhoInv
        call eos_re(eos_state)
        react_state_out % rhoY(:)     = yvec(1:nspec)
        react_state_out % rho         = sum(yvec(1:nspec))
    else
        eos_state % massfrac(1:nspec) = yvec(1:nspec)
        eos_state % T                 = yvec(neq)
        eos_state % h                 = (h_init  +  dt_react*hdot_ext)
        call eos_ph(eos_state)
        react_state_out % rhoY(:)     = yvec(1:nspec)*eos_state % rho
        react_state_out % rho         = eos_state % rho 
    end if
    CALL t_eos%stop        

    react_state_out % T = eos_state % T
    react_state_out % e = eos_state % e
    react_state_out % h = eos_state % h
    react_state_out % p = eos_state % p

    actual_react_cvode % reactions_succesful = .true.

    verbose = 0
    if (verbose .ge. 1) then
       write(6,*) '......cvode done:'
       ierr = FCVodeGetLastStep(CVmem, hlast)
       write(6,*) ' last successful step size = ',hlast
       ierr = FCVodeGetCurrentStep(CVmem, hcur)
       write(6,*) '          next step to try = ',hcur
       ierr = FCVodeGetCurrentTime(CVmem, tcur)
       write(6,*) '   integrated time reached = ',tcur
       ierr = FCVodeGetNumSteps(CVmem, nstp)
       write(6,*) '      number of time steps = ',nstp
       ierr = FCVodeGetNumRhsEvals(CVmem,nfevals)
       ierr = FCVDlsGetNumRhsEvals(CVmem,nfevals_jac)
       nfevals = nfevals + nfevals_jac
       write(6,*) '    number of fs(RHS EVAL) = ',nfevals
       ierr = FCVDlsGetNumJacEvals(CVmem, njevals)
       write(6,*) '              number of Js = ',njevals
       ierr = FCVodeGetLastOrder(CVmem,qlast)
       write(6,*) '    method order last used = ',qlast
       ierr = FCVodeGetCurrentOrder(CVmem,qcur) 
       write(6,*) '   method order to be used = ',qcur
       !write(6,*) '            number of LUDs = ',vodeiwork(19)
       ierr = FCVodeGetNumNonlinSolvIters(CVmem, nlinsetups)
       write(6,*) ' number of Newton iterations ',nlinsetups
       !write(6,*) ' number of Newton failures = ',vodeiwork(21)
    end if

    actual_react_cvode % cost_value = nfevals ! number of f evaluations

  end function actual_react_cvode


  integer(c_int)  function F_RHS_F(time, sunvec_y_in, sunvec_f_in, userdata) &
                  result(ierr) bind(C,name='F_RHS_F')

    use, intrinsic :: iso_c_binding
    use chemistry_module, only : molecular_weight
    use network, only: nspec
    use eos_module
    use cvode_interface
    use fnvector_serial_mod

    real(c_double), value :: time
    type(c_ptr), value    :: sunvec_y_in
    type(c_ptr), value    :: sunvec_f_in
    type(c_ptr), value    :: userdata

    ! pointers to data in SUNDAILS vectors
    real(c_double), pointer :: yvec_wk(:)
    real(c_double), pointer :: fvec(:)

    ! Local
    integer(c_int)                :: n
    integer(c_long)               :: neq
    real(amrex_real)              :: rhoInv

    neq = nspec + 1

    ! get data arrays from SUNDIALS vectors
    call FN_VGetData_Serial(sunvec_f_in, fvec)
    call FN_VGetData_Serial(sunvec_y_in, yvec_wk)

    CALL t_eos%start          
    if (iE == 1) then
        ! rhoY and T
        eos_state % rho = sum(yvec_wk(1:nspec))
        rhoInv          = 1.d0 / eos_state % rho
        eos_state % massfrac(1:nspec) = yvec_wk(1:nspec) * rhoInv
        eos_state % T   = yvec_wk(neq) ! guess
        eos_state % e   = (rhoe_init + (time - time_init) * rhoedot_ext) * rhoInv
        call eos_re(eos_state)
        call eos_get_activity(eos_state)
    else
        ! Y and T: probably not entirely correct
        eos_state % massfrac(1:nspec) = yvec_wk(1:nspec) 
        eos_state % T                 = yvec_wk(neq) ! guess
        eos_state % h   = h_init + (time - time_init) * hdot_ext
        call eos_ph(eos_state)
        call eos_get_activity_h(eos_state)
    end if
    CALL t_eos%stop        

    CALL t_ck%start          
    call ckwc(eos_state % T, eos_state % Acti, iwrk, rwrk, cdot)
    CALL t_ck%stop       

    if (iE == 1) then
        ! rhoY and T
        fvec(neq) = rhoedot_ext
        do n=1,nspec
           fvec(n) = cdot(n) * molecular_weight(n) + rhoydot_ext(n)
           fvec(neq) = fvec(neq) - eos_state%ei(n)*fvec(n)
        end do
        fvec(neq) = fvec(neq)/(eos_state%rho * eos_state%cv)
    else
        ! Y and T: probably not entirely correct
        fvec(neq) = hdot_ext
        do n=1,nspec
           fvec(n) = cdot(n) * molecular_weight(n)/eos_state%rho + ydot_ext(n)
           fvec(neq) = fvec(neq) - eos_state%hi(n)*fvec(n) 
        end do
        fvec(neq) = fvec(neq)/eos_state%cp
    end if

    ierr = 0
    return
  end function F_RHS_F


  integer(c_int) function f_jac_cvode(tn, sunvec_y_in, sunvec_f_in, sunMat_J, &
           user_data, tmp1, tmp2, tmp3) result(ierr) bind(C,name='f_jac_cvode')

        use, intrinsic :: iso_c_binding
        use amrex_error_module
        use chemistry_module, only : molecular_weight, inv_mwt
        use fnvector_serial_mod
        use fsunmat_dense_mod
        !use ode_params

        implicit none
        real(c_double),  value :: tn
        type(c_ptr),     value :: sunvec_y_in
        type(c_ptr),     value :: sunvec_f_in
        type(c_ptr),     value :: sunmat_J
        type(c_ptr),     value :: user_data
        type(c_ptr),     value :: tmp1, tmp2, tmp3

        ! pointers to data in SUNDAILS vector and matrix
        real(c_double), pointer :: yvec_wk(:)
        real(c_double), pointer :: Jmat(:,:)

        ! local variables
        integer            :: i, j
        integer(c_long)    :: neq
        real(amrex_real)   :: rho0, rhoinv, Temp 
        real(amrex_real)   :: YT0(nspec), C(nspec)
        integer, parameter :: consP = 0

        CALL t_AJac%start

        neq = nspec + 1

        ! get data array from SUNDIALS vector
        call FN_VGetData_Serial(sunvec_y_in, yvec_wk)
    
        ! get data array from SUNDIALS matrix
        call FSUNMatGetData_Dense(sunmat_J, Jmat)

        rho0 = sum(yvec_wk(1:nspec))
        YT0(:) = yvec_wk(1:nspec)/rho0
        Temp = yvec_wk(neq)

        call ckytcr(rho0, Temp, YT0, iwrk, rwrk, C)
        ! C in mol/cm3
        CALL t_ckJac%start          
        call DWDOT(Jmat, C, Temp, consP)
        CALL t_ckJac%stop          
        ! J(specs, specs) in 1/s

        !plot J
        !!!write(*,*) "m = ["
        !!!write(*,*) abs(Jmat(neq,neq)), abs(Jmat(1:neq-1,neq))
        !!!do i=1,neq-1
        !!!    write(*,*) abs(Jmat(neq,i)), abs(Jmat(1:neq-1,i))
        !!!write(*,*) "]"
        !!!end do

        write(34,*) abs(Jmat(neq,neq)), abs(Jmat(neq,1:neq-1))
        do i=1,neq-1
            write(34,*) abs(Jmat(i,neq)), abs(Jmat(i,1:neq-1))
        end do
        stop

        do j=1,nspec
           do i=1,nspec
              Jmat(i,j) = Jmat(i,j) * molecular_weight(i) * inv_mwt(j)
           end do
           i=neq
           Jmat(i,j) = Jmat(i,j) * inv_mwt(j) * rho0
        end do
    
        j = neq
        rhoinv = 1.d0/rho0    
        do i=1,nspec
           Jmat(i,j) = Jmat(i,j) * molecular_weight(i) * rhoinv
        enddo

        CALL t_AJac%stop

  end function f_jac_cvode


  integer(c_int) function f_jac_cvode_HP(tn, sunvec_y_in, sunvec_f_in, sunMat_J, &
           user_data, tmp1, tmp2, tmp3) result(ierr) bind(C,name='f_jac_cvode_HP')

        use, intrinsic :: iso_c_binding
        use amrex_error_module
        use chemistry_module, only : molecular_weight, inv_mwt
        use fnvector_serial_mod
        use fsunmat_dense_mod

        implicit none
        real(c_double),  value :: tn
        type(c_ptr),     value :: sunvec_y_in
        type(c_ptr),     value :: sunvec_f_in
        type(c_ptr),     value :: sunmat_J
        type(c_ptr),     value :: user_data
        type(c_ptr),     value :: tmp1, tmp2, tmp3

        ! pointers to data in SUNDAILS vector and matrix
        real(c_double), pointer :: yvec_wk(:)
        real(c_double), pointer :: Jmat(:,:)

        ! local variables
        integer            :: i, j
        integer(c_long)    :: neq
        real(amrex_real)   :: rho0, rhoinv, Temp 
        real(amrex_real)   :: YT0(nspec), C(nspec)
        integer, parameter :: consP = 0

        CALL t_AJac%start

        neq = nspec + 1

        ! get data array from SUNDIALS vector
        call FN_VGetData_Serial(sunvec_y_in, yvec_wk)
    
        ! get data array from SUNDIALS matrix
        call FSUNMatGetData_Dense(sunmat_J, Jmat)

        !rho0 = sum(yvec_wk(1:nspec))
        YT0(:) = yvec_wk(1:nspec)
        Temp = yvec_wk(neq)
        call ckrhoy(pressureInit, Temp, YT0, iwrk, rwrk, rho0)
        call ckytcr(rho0, Temp, YT0, iwrk, rwrk, C)
        ! C in mol/cm3
        CALL t_ckJac%start          
        call DWDOT(Jmat, C, Temp, consP)
        CALL t_ckJac%stop   
        ! J(specs, specs) in 1/s for specs, K is there for othres

    
        do j=1,nspec
           do i=1,nspec
              Jmat(i,j) = Jmat(i,j) * molecular_weight(i) * inv_mwt(j)
           end do
           i=neq
           Jmat(i,j) = Jmat(i,j) * inv_mwt(j) * rho0
        end do
    
        j = neq
        rhoinv = 1.d0/rho0    
        do i=1,nspec
           Jmat(i,j) = Jmat(i,j) * molecular_weight(i) * rhoinv
        enddo

        CALL t_AJac%stop

  end function f_jac_cvode_HP


  integer(c_int) function f_jac_cvode_KLU(tn, sunvec_y_in, sunvec_f_in, sunmat_J, &
           user_data, tmp1, tmp2, tmp3) result(ierr) bind(C,name='f_jac_cvode_KLU')

        use, intrinsic :: iso_c_binding
        use amrex_error_module
        use chemistry_module, only : molecular_weight, inv_mwt
        use fnvector_serial_mod
        use fsunmat_dense_mod
        use fsunmat_sparse_mod

        implicit none
        real(c_double),  value :: tn
        type(c_ptr),     value :: sunvec_y_in
        type(c_ptr),     value :: sunvec_f_in
        type(c_ptr),     value :: sunmat_J
        type(c_ptr),     value :: user_data
        type(c_ptr),     value :: tmp1, tmp2, tmp3

        ! pointers to data in SUNDAILS vector and matrix
        real(c_double), pointer :: yvec_wk(:)
        !type(c_ptr)             :: colptrs,rowvals,matdata
        real(c_double),  pointer :: f_data(:)
        integer(c_long), pointer :: f_idxval(:)
        integer(c_long), pointer :: f_idxptr(:)


        ! local variables
        integer            :: i, j, nbVals
        integer(c_long)    :: neq
        real(amrex_real)   :: rho0, rhoinv, Temp 
        real(amrex_real)   :: YT0(nspec), C(nspec)
        integer, parameter :: consP = 0

        CALL t_AJac%start

        neq = nspec + 1
        Jmat_KLU(:,:) = 0.0d0

        ! get data array from SUNDIALS vector
        call FN_VGetData_Serial(sunvec_y_in, yvec_wk)
    
        ! get data array from SUNDIALS matrix
        call FSUNMatGetData_Sparse(sunmat_J, f_data, f_idxval, f_idxptr)

        ! Stuff to calculate Jacobian
        rho0 = sum(yvec_wk(1:nspec))
        YT0(:) = yvec_wk(1:nspec)/rho0
        Temp = yvec_wk(neq)
        ! C in mol/cm3
        call ckytcr(rho0, Temp, YT0, iwrk, rwrk, C)
        ! J(specs, specs) in 1/s
        CALL t_ckJac%start          
        call DWDOT(Jmat_KLU, C, Temp, consP)
        CALL t_ckJac%stop   
        ! Renormalizations
        do j=1,nspec
           do i=1,nspec
              Jmat_KLU(i,j) = Jmat_KLU(i,j) * molecular_weight(i) * inv_mwt(j)
           end do
           i=neq
           Jmat_KLU(i,j) = Jmat_KLU(i,j) * inv_mwt(j) * rho0
        end do
        j = neq
        rhoinv = 1.d0/rho0    
        do i=1,nspec
           Jmat_KLU(i,j) = Jmat_KLU(i,j) * molecular_weight(i) * rhoinv
        enddo

        ! Fill sparse attributes
        !do i=0,neq
        !    f_idxptr(i+1) = (i)*neq  
        !    print *, (i)*neq, colPtrs(i+1)
        !end do
        !stop
        f_idxptr(:) = colPtrs(:)

        !do i=0,neq-1
        !    do j=1,neq
        !        f_idxval(i*neq+j) = j-1
        !        f_data(i*neq+j) = Jmat_KLU(i+1,j) 
        !    end do
        !end do
        f_idxval(:) = rowVals(:) - 1

        do i=2,neq+1
            nbVals = colPtrs(i)-colPtrs(i-1)
            !print *, "For column ", i-1, " nbVals ", nbVals
            do j=1,nbVals
                !print *, "   idx ", rowVals(colPtrs(i-1) + j),i-1 
                f_data(colPtrs(i-1)+j) = Jmat_KLU(i-1, rowVals(colPtrs(i-1) + j))
                !print *, "     -> seek data in Jmat_KLU in row ", rowVals(colPtrs(i-1)+j), Jmat_KLU(rowVals(colPtrs(i-1) + j),i-1)
            end do
        end do

        CALL t_AJac%stop

  end function f_jac_cvode_KLU

  integer(c_int) function f_jac_cvode_HP_KLU(tn, sunvec_y_in, sunvec_f_in, sunmat_J, &
           user_data, tmp1, tmp2, tmp3) result(ierr) bind(C,name='f_jac_cvode_HP_KLU')

        use, intrinsic :: iso_c_binding
        use amrex_error_module
        use chemistry_module, only : molecular_weight, inv_mwt
        use fnvector_serial_mod
        use fsunmat_dense_mod
        use fsunmat_sparse_mod

        implicit none
        real(c_double),  value :: tn
        type(c_ptr),     value :: sunvec_y_in
        type(c_ptr),     value :: sunvec_f_in
        type(c_ptr),     value :: sunmat_J
        type(c_ptr),     value :: user_data
        type(c_ptr),     value :: tmp1, tmp2, tmp3

        ! pointers to data in SUNDAILS vector and matrix
        real(c_double), pointer :: yvec_wk(:)
        real(c_double),  pointer :: f_data(:)
        integer(c_long), pointer :: f_idxval(:)
        integer(c_long), pointer :: f_idxptr(:)


        ! local variables
        integer            :: i, j, nbVals
        integer(c_long)    :: neq
        real(amrex_real)   :: rho0, rhoinv, Temp 
        real(amrex_real)   :: YT0(nspec), C(nspec)
        integer, parameter :: consP = 0

        CALL t_AJac%start

        neq = nspec + 1
        Jmat_KLU(:,:) = 0.0d0

        ! get data array from SUNDIALS vector
        call FN_VGetData_Serial(sunvec_y_in, yvec_wk)
    
        ! get data array from SUNDIALS matrix
        call FSUNMatGetData_Sparse(sunmat_J, f_data, f_idxval, f_idxptr)

        ! Stuff to calculate Jacobian
        YT0(:) = yvec_wk(1:nspec)
        Temp = yvec_wk(neq)
        call ckrhoy(pressureInit, Temp, YT0, iwrk, rwrk, rho0)
        ! C in mol/cm3
        call ckytcr(rho0, Temp, YT0, iwrk, rwrk, C)
        ! J(specs, specs) in 1/s
        CALL t_ckJac%start          
        call DWDOT(Jmat_KLU, C, Temp, consP)
        CALL t_ckJac%stop   

        ! Renormalizations
        do j=1,nspec
           do i=1,nspec
              Jmat_KLU(i,j) = Jmat_KLU(i,j) * molecular_weight(i) * inv_mwt(j)
           end do
           i=neq
           Jmat_KLU(i,j) = Jmat_KLU(i,j) * inv_mwt(j) * rho0
        end do
        j = neq
        rhoinv = 1.d0/rho0    
        do i=1,nspec
           Jmat_KLU(i,j) = Jmat_KLU(i,j) * molecular_weight(i) * rhoinv
        enddo

        ! Fill sparse attributes
        f_idxptr(:) = colPtrs(:)
        f_idxval(:) = rowVals(:) - 1

        do i=2,neq+1
            nbVals = colPtrs(i)-colPtrs(i-1)
            do j=1,nbVals
                f_data(colPtrs(i-1)+j) = Jmat_KLU(i-1, rowVals(colPtrs(i-1) + j))
            end do
        end do

        CALL t_AJac%stop

  end function f_jac_cvode_HP_KLU

end module actual_reactor_module
