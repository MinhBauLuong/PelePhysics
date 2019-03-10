module actual_sdc_module

  use, intrinsic :: iso_c_binding
  !use amrex_fort_module, only : amrex_real
  use mod_sdc_defs
  use network, only: nspec, spec_names
  use react_type_module
  use eos_type_module
  use mod_timers    

  implicit none

  real(amrex_real), private, allocatable :: cdot_k(:,:), cdot_kp1(:,:), I_k(:,:)
  real(amrex_real), private, allocatable :: delta_SDC(:,:)
  real(amrex_real), private, allocatable :: delta_SDC_max(:)
  type (eos_t), private, allocatable     :: eos_state_k(:), eos_state_kp1(:)
  real(amrex_real), private, allocatable :: dtLobato(:)
  integer                                :: nLobato, nsdcite

! Parameters of high order Gauss-Lobatto Lagrange weights       
  real(amrex_real), parameter :: sqrt5  = SQRT(5.0d0)
  real(amrex_real), parameter :: sqrt21 = SQRT(21.0d0)
  real(amrex_real), parameter :: tol    = 1.0d-04

contains

  !--------------------------!
  subroutine actual_reactor_init_sdc(iE_in,nLobato_in,nsdcite_in, iverbose_in)

      use bechem_module, only: init_bechem

      implicit none

      integer, intent(in)  :: iE_in, iverbose_in
      integer, intent(in)  :: nLobato_in, nsdcite_in
      integer              :: i

      !CALL t_total%init("Total")
      !CALL t_init%init("Initialization")
      !CALL t_bechem%init("Backward-Euler")
      !CALL t_ck%init("Chemkin calls")
      !CALL t_eos%init("EOS calls")

      !CALL t_total%start
      !CALL t_init%start

      iE      = iE_in
      nLobato = nLobato_in
      nsdcite = nsdcite_in
      verbose = iverbose_in

      allocate(eos_state_k(nLobato))
      allocate(eos_state_kp1(nLobato))
      if (iE == 1) then
          print *," ->with internal energy (UV cst)"
          allocate(rhoe_init(nLobato))
      else
          print *," ->with enthalpy (sort of HP cst)"
          allocate(rhoh_init(nLobato))
      end if
      allocate(rhoydot_ext(nspec))
      allocate(cdot_k(nLobato,nspec+1))
      allocate(cdot_kp1(nLobato,nspec+1))
      allocate(delta_SDC(nLobato,nspec+1))
      allocate(delta_SDC_max(nLobato))
      allocate(I_k(nLobato-1,nspec+1))
      allocate(dtLobato(nLobato-1))
      allocate(rhs(nspec+1))

      do i=1,nLobato
          call build(eos_state_k(i))
          call build(eos_state_kp1(i))
      end do

      ! Need molecular weights
      ! Redundant with use chemistry_module, only : molecular_weight
      allocate(mwt(nspec))
      allocate(invmwt(nspec))

      call CKWT(iwrk, rwrk, mwt)
      invmwt(:)   = 1.0/mwt(:)

      ! Allocate stuff for BE
      call init_bechem

      !CALL t_init%stop

  end subroutine actual_reactor_init_sdc
  !--------------------------!


  !--------------------------!
  function actual_react_sdc(react_state_in, react_state_out, dt_react, time)

      use amrex_error_module  
      use eos_module
      use chemistry_module, only : molecular_weight
      use, intrinsic :: iso_c_binding

      implicit none

      type(react_t),   intent(in   ) :: react_state_in
      type(react_t),   intent(inout) :: react_state_out
      real(amrex_real), intent(in  ) :: dt_react, time
      type(reaction_stat_t)          :: actual_react_sdc
      integer                        :: n, i, j, sdc
      real(amrex_real)               :: rhoInv

      if (verbose .ge. 1) then
          print *, "----------------------------" 
          print *, "    STARTING CHEMISTRY SOLVE" 
          print *, "----------------------------" 
      end if

!     Initialize dt intervals
!     Compute the Gauss-Lobatto intervals between 0 and dt
      call compute_dt_intervals(dt_react)

!     Setup the initial state of each iterations (tn)  
      !CALL t_eos%start          
      eos_state_k(1) % rho = sum(react_state_in % rhoY(:))  
      rhoInv               = 1.d0 / eos_state_k(1) % rho
      eos_state_k(1) % T   = react_state_in % T  
      eos_state_k(1) % massfrac(1:nspec) = react_state_in % rhoY(1:nspec) * rhoInv
      ! Deal with energy and compute T and conc + other useful eos var
      if (iE == 1) then
          eos_state_k(1) % e   = react_state_in % e
          rhoe_init(1)         = eos_state_k(1) % e  *  eos_state_k(1) % rho
          rhoedot_ext          = react_state_in % rhoedot_ext
          call eos_re(eos_state_k(1))
          call eos_get_activity(eos_state_k(1))
      else
          eos_state_k(1) % h   = react_state_in % h
          rhoh_init(1)         = eos_state_k(1) % h  *  eos_state_k(1) % rho
          rhohdot_ext          = react_state_in % rhohdot_ext
          call eos_rh(eos_state_k(1))
          call eos_get_activity_h(eos_state_k(1))
      end if

!     Determine if the system to solve is stiff 
      iStiff = 1
      if (eos_state_k(1) % T < 500.0) then
              iStiff = 1
      else if (eos_state_k(1) % T > 3000.0) then
              iStiff = 1
      end if

      !Compute species prod rates and correct with ext source term
      !CALL t_ck%start          
      call ckwc(eos_state_k(1) % T, eos_state_k(1) % Acti, iwrk, rwrk, cdot_k(1,1:nspec))
      !CALL t_ck%stop         
      !External forces for species
      rhoydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec)
      cdot_k(1,1:nspec) = cdot_k(1,1:nspec)* molecular_weight(1:nspec) + rhoydot_ext(:)
      !Compute T source term
      if (iE == 1) then
          cdot_k(1,nspec+1) = rhoedot_ext
          do n=1,nspec
              cdot_k(1,nspec+1) = cdot_k(1,nspec+1) - eos_state_k(1) % ei(n)*cdot_k(1,n)
          end do
          cdot_k(1,nspec+1) = cdot_k(1,nspec+1)/(eos_state_k(1) % rho * eos_state_k(1) % cv)
      else
          cdot_k(1,nspec+1) = rhohdot_ext
          do n=1,nspec
              cdot_k(1,nspec+1) = cdot_k(1,nspec+1) - eos_state_k(1) % hi(n)*cdot_k(1,n)
          end do
          cdot_k(1,nspec+1) = cdot_k(1,nspec+1)/(eos_state_k(1) % rho * eos_state_k(1) % cp)
      end if

!     Initialize the SDC loop by setting all Gauss-Lobatto points with the initial state (1)       
      do j=2,nLobato
          eos_state_k(j)   = eos_state_k(1)
          ! this is useless will need to remove
          if (iE == 1) then
              rhoe_init(j)     = rhoe_init(1)
          else
              rhoh_init(j)     = rhoh_init(1)
          end if
          cdot_k(j,:)      = cdot_k(1,:)
      end do

!     Start the loop
      if (verbose .ge. 1) then
          print *, "----------------------------" 
          print *, "    STARTING SDC LOOPS" 
          print *, "----------------------------" 
      end if
      do sdc = 1, nsdcite
           if (verbose .ge. 1) then
               write(*,'(A,I2,A)') " SDC ite ", sdc ,"  -------------------- "

               write(*,'(A,3(ES16.8))') "    T:",eos_state_k(1)%T, eos_state_k(2)%T, eos_state_k(3)%T 
               write(*,'(A,3(ES16.8))') "  rho:",eos_state_k(1)%rho, eos_state_k(2)%rho, eos_state_k(3)%rho 
           end if

!         First Lobatto point of each sdc iteration is always the same (previous dt)
          eos_state_kp1(1)     = eos_state_k(1)
          cdot_kp1(1,:)        = cdot_k(1,:)
          delta_SDC(1,:)       = 0.0d0
          delta_SDC_max(1)     = 0.0d0

!         Compute the quadrature for the state k     
          call compute_quadrature(I_k, cdot_k, dt_react)

!         Update the state with the correction integrals & the quadratures   
          do j=1,nLobato-1
             if (verbose .ge. 2) then
                 print *," - Working on the ", j+1, " Lobatto node"
             end if
             call sdc_advance_chem(eos_state_kp1(j+1), eos_state_kp1(j), eos_state_k(j+1), &
                      cdot_k(j,:), cdot_kp1(j,:), cdot_k(j+1,:), cdot_kp1(j+1,:), &
                      I_k(j,:), j)
          end do

          if (verbose .ge. 5) then
          !   Only works for ndodecane wang !
              write(*,'(A,3(ES16.8))') "   I_k(O2):",I_k(:,8)
              write(*,'(A,3(ES16.8))') "   c.k(O2):",cdot_k(:,8)
              write(*,'(A,3(ES16.8))') " c.kp1(O2):",cdot_kp1(:,8)
          end if

!         Copy state k+1 into k for the next iteration      
          if (verbose .ge. 2) then
              print *," - Update state and wdot for next iteration"
          end if

          do j=2,nLobato
              delta_SDC(j,nspec+1) = eos_state_k(j)%T - eos_state_kp1(j)%T
              delta_SDC(j,1:nspec) = eos_state_k(j)%massfrac(1:nspec) - eos_state_kp1(j)%massfrac(1:nspec)
              delta_SDC_max(j)     = maxval(abs(delta_SDC(j,:)))
              eos_state_k(j) = eos_state_kp1(j) 
              cdot_k(j,:)    = cdot_kp1(j,:)
          end do
          if (verbose .ge. 3) then
              write(*,'(A,3(ES16.8))') "    deltaYO2 (on each nLobatto) :",delta_SDC(:,8)
              write(*,'(A,3(ES16.8))') "    deltaM (on each nLobatto) :",delta_SDC_max(:) 
          end if


          if (verbose .ge. 1) then
              write(*,'(A,3(ES16.8))') "    T (end):",eos_state_kp1(1)%T, eos_state_kp1(2)%T, eos_state_kp1(3)%T 
              write(*,'(A,3(ES16.8))') "        rho:",eos_state_kp1(1)%rho, eos_state_kp1(2)%rho, eos_state_kp1(3)%rho 
              print *," "
              print *," "
          end if

          if ((maxval(delta_SDC_max(:)) < tol) .and. (sdc > 2*nLobato - 2)) then
              exit
          end if
      end do

!     Update out state      
      react_state_out % rho     = eos_state_kp1(nLobato)%rho
      react_state_out % rhoY(:) = react_state_out % rho * eos_state_kp1(nLobato)%massfrac(:)
      react_state_out % T = eos_state_kp1(nLobato) % T
      react_state_out % e = eos_state_kp1(nLobato) % e 
      react_state_out % h = eos_state_kp1(nLobato) % h
      if (iE == 1) then
          !react_state_out % e = eos_state_kp1(nLobato) % e
          react_state_out % rhoedot_ext = react_state_in % rhoedot_ext
      else
          !react_state_out % h = eos_state_kp1(nLobato) % h
          react_state_out % rhohdot_ext = react_state_in % rhohdot_ext
      end if
      react_state_out % rhoydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec)


!     Get a sense of cell difficulty 
      actual_react_sdc%reactions_succesful = .true.
      actual_react_sdc%cost_value = 1.0 

!     Timers 
      !CALL t_total%stop
      !CALL t_total%norm(t_total)
      !CALL t_bechem%norm(t_total)
      !CALL t_eos%norm(t_total)
      !CALL t_ck%norm(t_total)
      !!PRINTS
      !CALL t_total%print
      !CALL t_bechem%print
      !CALL t_eos%print
      !CALL t_ck%print
      !!RESETS
      !CALL t_total%reset
      !CALL t_bechem%reset
      !CALL t_eos%reset
      !CALL t_ck%reset

  end function actual_react_sdc

  subroutine actual_reactor_close_sdc()

      use bechem_module, only: close_bechem

      implicit none

      integer :: j

      if (allocated(cdot_k))      deallocate(cdot_k)
      if (allocated(cdot_kp1))    deallocate(cdot_kp1)
      if (allocated(rhoydot_ext)) deallocate(rhoydot_ext)
      if (allocated(I_k))         deallocate(I_k)
      if (allocated(delta_SDC))   deallocate(delta_SDC)
      if (allocated(delta_SDC_max))   deallocate(delta_SDC_max)
      if (allocated(dtLobato))    deallocate(dtLobato)
      if (allocated(rhoe_init))   deallocate(rhoe_init)
      if (allocated(rhoh_init))   deallocate(rhoh_init)

      do j=1,nLobato
          call destroy(eos_state_k(j))
          call destroy(eos_state_kp1(j))
      end do

      deallocate(eos_state_k)
      deallocate(eos_state_kp1)

      if (allocated(mwt))    deallocate(mwt)
      if (allocated(invmwt)) deallocate(invmwt)
      
      call close_bechem

  end subroutine actual_reactor_close_sdc


  subroutine compute_dt_intervals(dt)

      double precision, intent(in ) :: dt

      if (verbose .ge. 1) then
          print *,"Computing dt intervals using ", nLobato, " Lobatto nodes"
      end if

      if (nLobato .eq. 2) then
          dtLobato(1) = dt
      else if (nLobato .eq. 3) then
          dtLobato(1) = 0.5d0*dt
          dtLobato(2) = 0.5d0*dt
      else if (nLobato .eq. 4) then
          dtLobato(1) = (0.5d0-0.1d0*sqrt5)*dt
          dtLobato(2) = (0.2d0*sqrt5)*dt
          dtLobato(3) = (0.5d0-0.1d0*sqrt5)*dt
      else if (nLobato .eq. 5) then
          dtLobato(1) = (0.5d0-sqrt21/14.0)*dt
          dtLobato(2) = (sqrt21/14.0)*dt
          dtLobato(3) = (sqrt21/14.0)*dt
          dtLobato(4) = (0.5d0-sqrt21/14.0)*dt
      else
          write(*,*) "ERROR : nLobato cannot exceed 5 right now"
          stop
      end if

  end subroutine compute_dt_intervals


  subroutine compute_quadrature(I, f, dt)

      double precision, intent(out) :: I(nLobato-1,nspec+1)
      double precision, intent(in ) :: f(nLobato,nspec+1) 
      double precision, intent(in ) :: dt 

      if (verbose .ge. 1) then
          print *," - Computing quadrature on ", nLobato, " Lobatto nodes"
      end if

      ! We get the interpolating legendre poly bet 0 and dt
      if (nLobato .eq. 2) then
          I(1,:) = 0.5d0*dt*(f(1,:)+f(2,:)) 
      else if (nLobato .eq. 3) then
          I(1,:) = (5.0d0*f(1,:) + 8.0d0*f(2,:) - f(3,:))*dt/24.d0
          I(2,:) = (- f(1,:) + 8.0d0*f(2,:) + 5.0d0*f(3,:))*dt/24.d0
      else if (nLobato .eq. 4) then
          I(1,:) = ((11.0d0+sqrt5)*f(1,:) + &
                   (25.0d0-sqrt5)*f(2,:) + &
                   (25.0d0-13.0d0*sqrt5)*f(3,:) + &
                   (-1.0d0+sqrt5)*f(4,:))*dt/120.d0
          I(2,:) = ((-2.0d0*sqrt5)*f(1,:) + &
                   (14.0d0*sqrt5)*f(2,:) + &
                   (14.0d0*sqrt5)*f(3,:) + &
                   (-2.0d0*sqrt5)*f(4,:))*dt/120.d0
          I(3,:) = ((-1.0d0+sqrt5)*f(1,:) + &
                   (25.0d0-13.0d0*sqrt5)*f(2,:) + &
                   (25.0d0-sqrt5)*f(3,:) + &
                   (11.0d0+sqrt5)*f(4,:))*dt/120.d0
      else if (nLobato .eq. 5) then
          I(1,:) = ((119.0d0+3.0d0*sqrt21)*f(1,:)/1960.0d0 + &    
                   (343.0d0-9.0d0*sqrt21)*f(2,:)/2520.0d0 + & 
                   (392.0d0-96.0d0*sqrt21)*f(3,:)/2205.0d0 + & 
                   (343.0d0-69.0d0*sqrt21)*f(4,:)/2520.0d0 + & 
                   (-21.0d0+3.0d0*sqrt21)*f(5,:)/1960.0d0)*dt
          I(2,:) = ((-315.0d0-24.0d0*sqrt21)*f(1,:)/15680.0d0 + &    
                   (2421.0d0*sqrt21)*f(2,:)/6048.0d0 + & 
                   (96.0d0*sqrt21)*f(3,:)/2205.0d0 + & 
                   (-549.0d0*sqrt21)*f(4,:)/6048.0d0 + & 
                   (315.0d0-24.0d0*sqrt21)*f(5,:)/15680.0d0)*dt
          I(3,:) = ((315.0d0-24.0d0*sqrt21)*f(1,:)/15680.0d0 + &    
                   (-549.0d0*sqrt21)*f(2,:)/6048.0d0 + & 
                   (96.0d0*sqrt21)*f(3,:)/2205.0d0 + & 
                   (2421.0d0*sqrt21)*f(4,:)/6048.0d0 + & 
                   (-315.0d0-24.0d0*sqrt21)*f(5,:)/15680.0d0)*dt
          I(4,:) = ((-21.0d0+3.0d0*sqrt21)*f(1,:)/1960.0d0 + &    
                   (343.0d0-69.0d0*sqrt21)*f(2,:)/2520.0d0 + & 
                   (392.0d0-96.0d0*sqrt21)*f(3,:)/2205.0d0 + & 
                   (343.0d0-9.0d0*sqrt21)*f(4,:)/2520.0d0 + & 
                   (119.0d0+3.0d0*sqrt21)*f(5,:)/1960.0d0)*dt
      else
          write(*,*) "ERROR : nLobato cannot exceed 5 right now"
          stop
      end if

  end subroutine compute_quadrature


  subroutine sdc_advance_chem(state_kp1_jp1, state_kp1_j, state_k_jp1, cdot_k_j, cdot_kp1_j, &
                              cdot_k_jp1, cdot_kp1_jp1, I_k_lcl, nlobato)

      use bechem_module, only: bechem
      use chemistry_module, only : molecular_weight
      use eos_module

      implicit none

      type(eos_t),   intent(inout  ) :: state_kp1_jp1
      type(eos_t),   intent(in     ) :: state_kp1_j
      type(eos_t),   intent(in     ) :: state_k_jp1
      double precision, intent(in  ) :: cdot_k_j(nspec+1)
      double precision, intent(in  ) :: cdot_kp1_j(nspec+1)
      double precision, intent(in  ) :: cdot_k_jp1(nspec+1)
      double precision, intent(out ) :: cdot_kp1_jp1(nspec+1)
      double precision, intent(in  ) :: I_k_lcl(nspec+1)
      integer, intent(in  )          :: nlobato

      integer           :: n, i, ierr
      double precision  :: rho, Temp
      double precision  :: dtLobato_lcl
      !double precision  :: rYguess(nspec)
      double precision  :: rY(nspec)
      double precision  :: hguess, eguess

      dtLobato_lcl = dtLobato(nlobato)

!     SDC rhs 
      if (verbose .ge. 3) then
          print *,"   - Updating the RHS "
      end if
      if (iStiff == 1) then
          do i=1,nspec
              rhs(i) = state_kp1_j % rho * state_kp1_j % massfrac(i) - (dtLobato_lcl * cdot_k_jp1(i) - I_k_lcl(i))
          end do
          rhs(nspec+1) = state_kp1_j % T - (dtLobato_lcl * cdot_k_jp1(nspec+1) -  I_k_lcl(nspec+1))
      else
          do i=1,nspec
              rhs(i) = state_kp1_j % rho * state_kp1_j % massfrac(i)  + dtLobato_lcl * (cdot_kp1_j(i) - cdot_k_j(i)) + I_k_lcl(i)
          end do
          rhs(nspec+1) = state_kp1_j % T + dtLobato_lcl * (cdot_kp1_j(nspec+1) - cdot_k_j(nspec+1)) +  I_k_lcl(nspec+1)
      end if

!     Define initial state
      rho        = state_k_jp1 % rho
      Temp       = state_k_jp1 % T

!     ... and solve with newton iterations
      if (iStiff == 1) then
            if (verbose .ge. 3) then
                print *,"   - Call the BE solver (iStiff == 1) "
            end if
            !rYguess(:) = state_k_jp1 % massfrac(:) * rho
            rY(:) = state_k_jp1 % massfrac(:) * rho
            if (iE == 1) then
                call bechem(rY, rho, Temp, dtLobato_lcl)
            else
                call bechem(rY, rho, Temp, dtLobato_lcl)
            end if
!     ... or simple explicit BE scheme
      else
            if (verbose .ge. 3) then
                print *,"   - Simple BE solve (iStiff == 0) "
            end if
            do i=1,nspec
                rY(i) = rhs(i)
            end do
            rho = sum(rY(:))
            Temp = rhs(1+nspec)
      end if

!     Update output state (kp1, jp1)
      state_kp1_jp1%massfrac(1:nspec) = rY(:) / rho
      state_kp1_jp1%T                 = Temp
      state_kp1_jp1%rho               = rho
      !print *, Temp, state_kp1_jp1%T
      if (iE == 1) then
          state_kp1_jp1%e             = (rhoe_init(1) + rhoedot_ext * sum(dtLobato(1:nlobato))) / rho
          !print *, state_kp1_j % rho, rho
          call eos_re(state_kp1_jp1)
          if (iStiff == 0) then
              call eos_get_activity(state_kp1_jp1)
          end if
      else
          state_kp1_jp1%h             = (rhoh_init(1) + rhohdot_ext * sum(dtLobato(1:nlobato))) / rho
          call eos_rh(state_kp1_jp1)
          if (iStiff == 0) then
              call eos_get_activity_h(state_kp1_jp1)
          end if
      end if
      !print *, Temp, state_kp1_jp1%T


!     Output new cdot (kp1, jp1)
      if (iStiff == 1) then
          do i=1,nspec
              cdot_kp1_jp1(i) = (rY(i) - rhs(i)) / dtLobato_lcl
          end do
      else
          call ckwc(state_kp1_jp1 % T, state_kp1_jp1 % Acti, iwrk, rwrk, cdot_kp1_jp1(1:nspec))
          cdot_kp1_jp1(1:nspec) = cdot_kp1_jp1(1:nspec)* molecular_weight(1:nspec) + rhoydot_ext(:)
      end if

      !Compute T source term
      if (iE == 1) then
          cdot_kp1_jp1(nspec+1) = rhoedot_ext
          do n=1,nspec
              cdot_kp1_jp1(nspec+1) = cdot_kp1_jp1(nspec+1) - state_kp1_jp1 % ei(n)*cdot_kp1_jp1(n)
          end do
          cdot_kp1_jp1(nspec+1) = cdot_kp1_jp1(nspec+1)/(state_kp1_jp1 % rho * state_kp1_jp1 % cv)
      else
          cdot_kp1_jp1(nspec+1) = rhohdot_ext
          do n=1,nspec
              cdot_kp1_jp1(nspec+1) = cdot_kp1_jp1(nspec+1) - state_kp1_jp1 % hi(n)*cdot_kp1_jp1(n)
          end do
          cdot_kp1_jp1(nspec+1) = cdot_kp1_jp1(nspec+1)/(state_kp1_jp1 % rho * state_kp1_jp1 % cp)
      end if

  end subroutine sdc_advance_chem

end module actual_sdc_module
