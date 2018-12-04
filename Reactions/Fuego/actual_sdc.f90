module actual_sdc_module

!#include "petsc/finclude/petscsys.h"
!#include "petsc/finclude/petscvec.h"
!#include "petsc/finclude/petscmat.h"
!#include "petsc/finclude/petscpc.h"
!#include "petsc/finclude/petscksp.h"

  !USE petscvec
  !USE petscsys
  !USE petscmat
  !USE petscpc
  !USE petscksp

  use amrex_fort_module, only : amrex_real
  use network, only: nspec, spec_names
  use react_type_module
  use eos_type_module
  use mod_timers    

  implicit none

  real(amrex_real), private, allocatable :: cdot_k(:,:), cdot_kp1(:,:), I_k(:,:)
  real(amrex_real), private, allocatable :: rhoydot_ext(:), rhoh_init(:)
  type (eos_t), private, allocatable     :: eos_state_k(:), eos_state_kp1(:)
  real(amrex_real), private, allocatable :: dtLobato(:)
  integer                   :: nLobato, nsdcite
  real(amrex_real)          :: rhohdot_ext

contains

  subroutine actual_reactor_init_sdc(nLobato_in,nsdcite_in)

      implicit none

      integer, intent(in)  :: nLobato_in, nsdcite_in
      integer              :: i
      !PetscErrorCode       :: ierr

      !CALL t_total%init("Total")
      !CALL t_init%init("Initialization")
      !CALL t_bechem%init("Backward-Euler")
      !CALL t_ck%init("Chemkin calls")
      !CALL t_eos%init("EOS calls")

      !CALL t_total%start
      !CALL t_init%start

      nLobato = nLobato_in
      nsdcite = nsdcite_in

      allocate(eos_state_k(nLobato))
      allocate(eos_state_kp1(nLobato))
      allocate(rhoh_init(nspec))
      allocate(rhoydot_ext(nspec))
      allocate(cdot_k(nLobato,nspec+1))
      allocate(cdot_kp1(nLobato,nspec+1))
      allocate(I_k(nLobato-1,nspec+1))
      allocate(dtLobato(nLobato-1))

      !call PetscInitialize(PETSC_NULL_CHARACTER,ierr)

      do i=1,nLobato
          call build(eos_state_k(i))
          call build(eos_state_kp1(i))
      end do

      !CALL t_init%stop

  end subroutine actual_reactor_init_sdc

  function actual_react_sdc(react_state_in, react_state_out, dt_react, time)

      use amrex_error_module  
      use eos_module
      use chemistry_module, only : molecular_weight
      use, intrinsic :: iso_c_binding

      implicit none

      type(react_t),   intent(in   ) :: react_state_in
      type(react_t),   intent(inout) :: react_state_out
      real(amrex_real), intent(in   ) :: dt_react, time
      type(reaction_stat_t)          :: actual_react_sdc
      integer          :: i, j, sdc
      real(amrex_real) :: rhoInv

!     Initialize dt intervals
      call compute_dt_intervals(dt_react)

!     Setup the initial state of each iterations (tn)  
      !CALL t_eos%start          
      eos_state_k(1) % rho = sum(react_state_in % rhoY(:))  
      rhoInv               = 1.d0 / eos_state_k(1) % rho
      eos_state_k(1) % T   = react_state_in % T  
      eos_state_k(1) % massfrac(1:nspec) = react_state_in % rhoY(1:nspec) * rhoInv
      eos_state_k(1) % h   = react_state_in % h
      rhoh_init(1)         = eos_state_k(1) % h  *  eos_state_k(1) % rho
      eos_state_k(1) % p   = react_state_in % p

!     External forces
      rhohdot_ext          = react_state_in % rhohdot_ext
      rhoydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec)

!     Compute wdot of initial state
      call eos_get_activity_h(eos_state_k(1))
      !CALL t_eos%stop         

      !CALL t_ck%start          
      call ckwc(eos_state_k(1) % T, eos_state_k(1) % Acti, iwrk, rwrk, cdot_k(1,1:nspec))
      !CALL t_ck%stop         
      cdot_k(1,1:nspec) = cdot_k(1,1:nspec)* molecular_weight(1:nspec) + rhoydot_ext(:)
      cdot_k(1,nspec+1) = rhohdot_ext

!     Initialize the SDC loop by setting all Gauss-Lobato point with the initial state        
      do j=2,nLobato
          eos_state_k(j)   = eos_state_k(1)
          rhoh_init(j)     = rhoh_init(1)
          cdot_k(j,:)      = cdot_k(1,:)
      end do

!     Start the loop
      do sdc = 1, nsdcite

!          write(*,'(A,I2,A)') " SDC ite ", sdc ,"  -------------------- "
!          write(*,'(A,3(ES12.4))') "   T:",eos_state_k(1)%T, eos_state_k(2)%T, eos_state_k(3)%T      
!          write(*,'(A,3(ES12.4))') " rho:",eos_state_k(1)%rho, eos_state_k(2)%rho, eos_state_k(3)%rho      
!          write(*,'(A,3(ES12.4))') "  Yf:",eos_state_k(1)%massfrac(1), eos_state_k(2)%massfrac(1),eos_state_k(3)%massfrac(1)      

!         First Lobato point of each sdc iteration is always the same
          eos_state_kp1(1)     = eos_state_k(1)
          cdot_kp1(1,:)        = cdot_k(1,:)

!         Compute the quadrature for the state k     
          call compute_quadrature(I_k, cdot_k, dt_react)

!         Update the state with the correction integrals & the quadratures   
          do j=1,nLobato-1
             ! eos_state_kp1(j+1) = eos_state_k(j+1) 
              call sdc_advance_chem(eos_state_kp1(j+1), eos_state_kp1(j), eos_state_k(j+1), &
                      cdot_k(j+1,:), cdot_kp1(j+1,:), I_k(j,:), dtLobato(j), dt_react)
          end do
!          write(*,'(A,3(ES12.4))') "   I_k:",I_k(:,1)
!          write(*,'(A,3(ES12.4))') "   c.k:",cdot_k(:,1)
!          write(*,'(A,3(ES12.4))') " c.kp1:",cdot_kp1(:,1)

!         Copy state k+1 into k for the next iteration      
          do j=2,nLobato
             eos_state_k(j) = eos_state_kp1(j) 
             cdot_k(j,:)    = cdot_kp1(j,:)
          end do

      end do

!     Update out state      
      react_state_out % rho     = eos_state_kp1(nLobato)%rho
      !write(*,*) " eos_state_out % rho", react_state_out % rho
      react_state_out % rhoY(:) = react_state_out % rho * eos_state_kp1(nLobato)%massfrac(:)
      react_state_out % T = eos_state_kp1(nLobato) % T
      react_state_out % h = eos_state_kp1(nLobato) % h

      react_state_out % rhohdot_ext = react_state_in % rhohdot_ext
      react_state_out % rhoydot_ext(1:nspec) = react_state_in % rhoydot_ext(1:nspec)


      actual_react_sdc%reactions_succesful = .true.
      actual_react_sdc%cost_value = 1.0 

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

      implicit none
      integer :: j

      if (allocated(cdot_k)) deallocate(cdot_k)
      if (allocated(cdot_kp1)) deallocate(cdot_kp1)

      do j=1,nLobato
          call destroy(eos_state_k(j))
          call destroy(eos_state_kp1(j))
      end do

      deallocate(eos_state_k)
      deallocate(eos_state_kp1)
   
  end subroutine actual_reactor_close_sdc


  subroutine compute_dt_intervals(dt)

      double precision, intent(in ) :: dt

      if (nLobato .eq. 2) then
          dtLobato(1) = dt
      else if (nLobato .eq. 3) then
          dtLobato(:) = 0.5d0*dt
      else
          write(*,*) "ERROR : nLobato cannot exceed 3 right now"
          stop
      end if

  end subroutine compute_dt_intervals


  subroutine compute_quadrature(I, f, dt)

      double precision, intent(out) :: I(nLobato-1,nspec+1)
      double precision, intent(in ) :: f(nLobato,nspec+1) 
      double precision, intent(in ) :: dt 

      ! We get the interpolating legendre poly bet 0 and dt
      if (nLobato .eq. 2) then
          I(1,:) = 0.5d0*dt*(f(1,:)+f(2,:)) 
      else if (nLobato .eq. 3) then
          I(1,:) = (5.0d0*f(1,:) + 8.0d0*f(2,:) - f(3,:))*dt/24.d0
          I(2,:) = (- f(1,:) + 8.0d0*f(2,:) + 5.0d0*f(3,:))*dt/24.d0
      else
          write(*,*) "ERROR : nLobato cannot exceed 3 right now"
          stop
      end if

  end subroutine compute_quadrature


  subroutine sdc_advance_chem(state_kp1_jp1, state_kp1_j, state_k_jp1, cdot_k, cdot_kp1, I_k_lcl, dtLobato, dt)

      use bechem_module, only: bechem
      use eos_module
      !include 'spec.h'

      implicit none

      type(eos_t),   intent(inout  ) :: state_kp1_jp1
      type(eos_t),   intent(in     ) :: state_kp1_j
      type(eos_t),   intent(in     ) :: state_k_jp1
      double precision, intent(in   )  :: cdot_k(nspec+1)
      double precision, intent(out  )  :: cdot_kp1(nspec+1)
      double precision, intent(in   )  :: I_k_lcl(nspec+1)
      double precision, intent(in   )  :: dtLobato 
      double precision, intent(in   )  :: dt 

      integer :: i, ierr
      double precision  :: rho_init
      double precision  :: T_init
      double precision  :: Yguess(nspec)
      double precision  :: Y(nspec)
      double precision  :: hguess
      double precision  :: rhs(nspec+1)

!     SDC rhs 
      do i=1,nspec
          rhs(i) = state_kp1_j % massfrac(i) - (dtLobato * cdot_k(i) - I_k_lcl(i)) / state_k_jp1 % rho
      end do
      rhs(nspec+1) = state_kp1_j % h - (dtLobato * cdot_k(nspec+1) -  I_k_lcl(nspec+1)) / state_k_jp1 % rho

!     Define initial state
      rho_init  = state_k_jp1 % rho
      T_init    = state_k_jp1 % T
      Yguess(:) = state_k_jp1 % massfrac(:) 
      hguess    = state_k_jp1 % h
!     ... and solve with newton iterations
      call bechem(Y, Yguess, hguess, rho_init, T_init, rhs, rhohdot_ext, dtLobato)

!     Output new cdot (kp1, jp1)
      do i=1,nspec
          cdot_kp1(i) = rho_init * (Y(i) - rhs(i)) / dtLobato 
      end do
      cdot_kp1(1+nspec) = cdot_k(nspec+1)

      !state_kp1_jp1%massfrac(:) = Y(:)
      state_kp1_jp1%massfrac(1:nspec) = state_kp1_j % massfrac(1:nspec) + (dtLobato * (cdot_kp1(1:nspec) - cdot_k(1:nspec)) + I_k_lcl(1:nspec))/ rho_init 
      !state_kp1_jp1%h           = hguess
      state_kp1_jp1%h           = state_kp1_j % h + I_k_lcl(nspec+1)/rho_init 
      !CALL t_eos%start
      call get_T_given_hY(hguess,Y,iwrk,rwrk,T_init,ierr)
      !CALL t_eos%stop
      state_kp1_jp1%T           = T_init
      state_kp1_jp1%p           = state_k_jp1 % p
      !CALL t_eos%start
      call ckrhoy(state_kp1_jp1%p,T_init,state_kp1_jp1%massfrac(:),iwrk,rwrk,rho_init)
      !CALL t_eos%stop
      state_kp1_jp1%rho =  rho_init

  end subroutine sdc_advance_chem

end module actual_sdc_module
