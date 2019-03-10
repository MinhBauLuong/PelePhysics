module bechem_module

   use network, only: nspec
   use mod_timers
   use mod_sdc_defs
   
   implicit none

   private

   !     Jacobian matrix and factorization
   double precision, allocatable, save :: Jac(:,:), A(:,:)
   double precision, allocatable, save :: AT(:,:)
   !     Pivot vector returned from LINPACK
   integer, allocatable, save :: ipvt(:)

   double precision :: rcond
   integer          :: consP

   integer, parameter          :: max_iter = 100
   double precision, parameter :: tol = 1.d-06

   logical                     :: recompute_jac

   public :: bechem, init_bechem, close_bechem

contains

   !-----------------------!
   subroutine init_bechem

     ! Allocate space for Jac, factorization and pivot
     !if (.not. allocated(A)) then
        allocate(Jac(Nspec+1,Nspec+1))
        allocate(A(Nspec+1,Nspec+1))
        allocate(AT(Nspec+1,Nspec+1))
        allocate(ipvt(Nspec+1))
     !end if

     ! Initially assume Jac is wrong
     recompute_jac = .true.

     if (iE == 1) then
          consP = 0
     else
          consP = 1
     end if
     
   end subroutine init_bechem
   !-----------------------!


   !-----------------------!
   subroutine close_bechem

        if (allocated(Jac)) deallocate(Jac)
        if (allocated(A)) deallocate(A)
        if (allocated(AT)) deallocate(AT)
        if (allocated(ipvt)) deallocate(ipvt)

   end subroutine close_bechem
   !-----------------------!


   !-----------------------!
   subroutine bechem(rY, rho, T, dt)
     ! Do a Backward Euler solve for the chemistry using Newton's method
     !        rY : solution (output)
     !       rY0 : initial guess
     !       rho : initial density
     !         T : inout temperature
     !        dt : timestep

     double precision, intent(inout)  :: rY(Nspec)
     double precision, intent(inout ) :: rho
     double precision, intent(inout ) :: T
     double precision, intent(in )    :: dt
     
     integer          :: iter, n, ierr
     integer          :: flag
     double precision :: rmax, rho_inv, Tdot, res_nl_norm
     double precision, dimension(Nspec+1) :: res_nl
     double precision, dimension(Nspec+1) :: dx_nl
     double precision, dimension(Nspec)   :: Y, wdot
     double precision, dimension(Nspec+1) :: val_ref

     !CALL t_bechem%start

     if (verbose .ge. 4) then
         print *,"   +++++++++++++++++++++++++++++++++" 
         print *,"     STARTING NEWTON ITERATIONS     "
     end if

     !val_ref(:) = 1.0/(res_nl_init(:) + 1.0d-10) + 1.0d-10*res_nl_init(:)
     val_ref(:) = 1.0

     ! maximum number of iterations is max_iter
     do iter = 0, max_iter
        ! Newton's method: iteratively solve M(x_{n+1} - x_n) = -F(x_n)
        ! M is given by (I - dt*wdot) 
        ! F(x) is given by (x_n - dt*wdot - rhs)

        !if (iter .eq. 0) then
            rho_inv = 1.d0 / rho
            Y  = rY(:) * rho_inv

            if (verbose .ge. 5) then
                print *,"     Working on the ", iter, " iteration (T, Y(O2)) ", T, Y(8)
            end if

            ! compute Tdot and wdot
            call compute_source_terms(Y,T,rho,wdot,Tdot)

            ! Compute initial residuals
            call compute_residuals(rY,T,wdot,Tdot,dt,res_nl)

        !end if

        call check_max_residuals(res_nl,val_ref,iter,rho,rY,T,rmax,res_nl_norm,flag)
        if (flag == 1) then
            return
        end if

        !     compute the Jacobian, and then compute its LU factorization
        if (recompute_jac) then
            call compute_jac(rho,Y,T)
        end if

        call LUA(dt)
        !     call LINPACK to solve the linear system M (dx) = res
        !     using the output from the factorization given by dgefa

        dx_nl = res_nl
        call dgesl(A, Nspec+1, Nspec+1, ipvt, dx_nl, 0)
        !     in dx_nl now is delta_x and res_nl is the Newton residual

        !     solve for the difference x = x_{n+1} - x_n
        !     so the next iterate is given by x_{n+1} = x_n + x
        call linesearch(rY,T,rho,res_nl,dt,dx_nl,flag)
        if (flag == 1) then
            return
        end if
        !call update_vals(rY,T,rho,dx_nl)

     end do

     if (iter .ge. max_iter) then
        print *,"     Newton solve failed to converge !! " 
        print *,'     bechem: iter=',iter
        print *,'     bechem: rmax=',rmax
        print *,'     bechem: rL2=',res_nl_norm
        stop
     endif

     !CALL t_bechem%stop

   end subroutine bechem
   !-----------------------!


   !-----------------------!
   subroutine compute_source_terms(Y,T,rho,wdot,Tdot)

     double precision, intent(in)     :: Y(Nspec)
     double precision, intent(in )    :: T, rho
     double precision, intent(out)    :: wdot(Nspec)
     double precision, intent(out)    :: Tdot

     double precision                     :: cp, cp_inv, cv, cv_inv
     double precision                     :: rho_inv
     double precision, dimension(Nspec)   :: hi, ei
     integer                              :: n

     ! multiply by molecular weight to get the right units
     ! and add ext source term
     call CKWYR(rho, T, Y, iwrk, rwrk, wdot)
     wdot(:) = wdot(:) * mwt(:) + rhoydot_ext(:)

     ! compute C_p/v 
     if (iE == 1) then
         call CKCVBS(T, Y, iwrk, rwrk, cv)
         cv_inv = 1.d0/cv
     else 
         call CKCPBS(T, Y, iwrk, rwrk, cp)
         cp_inv = 1.d0/cp
     end if
        
     ! compute T&spec src terms
     rho_inv = 1.0d0 / rho
     if (iE == 1) then
         call CKUMS (T, iwrk, rwrk, ei)
         Tdot = rhoedot_ext 
         do n=1,Nspec
             Tdot = Tdot - ei(n) * wdot(n) 
         end do
         Tdot = Tdot * cv_inv * rho_inv
     else 
         call CKHMS (T, iwrk, rwrk, hi)
         Tdot = rhohdot_ext 
         do n=1,Nspec
             Tdot = Tdot - hi(n) * wdot(n) 
         end do
         Tdot = Tdot * cp_inv * rho_inv
     end if

   end subroutine compute_source_terms
   !-----------------------!


   !-----------------------!
   subroutine compute_residuals(rY,T,wdot,Tdot,dt,res_nl)

     double precision, intent(in)     :: rY(Nspec),wdot(Nspec)
     double precision, intent(in )    :: T,Tdot,dt
     double precision, intent(out)    :: res_nl(Nspec+1)

        res_nl(1:Nspec) = -(rY(:) - dt*wdot(:) - rhs(1:Nspec))
        res_nl(Nspec+1) = -(T - dt * Tdot - rhs(Nspec+1))

   end subroutine compute_residuals
   !-----------------------!


   !-----------------------!
   subroutine check_max_residuals(res_nl,val_ref,iter,rho,rY,T,rmax,res_nl_norm,flag)

     double precision, intent(in)     :: res_nl(Nspec+1),val_ref(Nspec+1)
     double precision, intent(in )    :: rho,T
     double precision, intent(in)     :: rY(Nspec)
     integer, intent(in)              :: iter
     double precision, intent(out)    :: rmax, res_nl_norm
     integer, intent(out)             :: flag

     flag = 0

     rmax = maxval(abs(res_nl))
     res_nl_norm = 0.5*NORM2(res_nl(:)*val_ref(:))**2.0
     if (verbose .ge. 5) then
         print *,"     L2 residual, Max residual: ", res_nl_norm, rmax
     end if

     if (isnan(rmax)) then
           print *," "
           print *,"     BE solve returned NaN !! " 
           print *," "
           print *,'     bechem: iteration: ', iter
           print *,'     bechem: reciprocal condition number of J = ', rcond
           print *,'     bechem: density = ', rho
           print *,'     bechem: temperature = ', T
           print *,'     bechem: rY = ',rY
           print *,"   +++++++++++++++++++++++++++++++++" 
           stop
     endif

     ! if we have reached the desired tolerance then we are done
     if (rmax .le. tol) then
        if (verbose .ge. 4) then
            print *,"       --> Newton has converged !! <--  " 
            print *,"   +++++++++++++++++++++++++++++++++" 
        end if
        flag = 1 
        return
     endif

   end subroutine check_max_residuals
   !-----------------------!


   !-----------------------!
   subroutine compute_jac(rho,Y,T)

     double precision, intent(inout)  :: Y(Nspec)
     double precision, intent(in )    :: rho,T

     double precision, dimension(Nspec) :: C
     integer                            :: i,j

     !    compute the molar concentrations
     call CKYTCR(rho, T, Y, iwrk, rwrk, C)
     !     use the concentrarions to compute the reaction Jacobian
     call DWDOT(Jac, C, T, consP)

     !     convert to the appropriate units
     do j=1,Nspec
        do i=1,Nspec
           Jac(i,j) = Jac(i,j) * mwt(i) * invmwt(j)
        end do
        i=Nspec+1
        Jac(i,j) = Jac(i,j) * invmwt(j) 
     end do

     j = Nspec+1
     do i=1,Nspec
        Jac(i,j) = Jac(i,j) * mwt(i) 
     end do

   end subroutine compute_jac
   !-----------------------!


   !-----------------------!
   subroutine LUA(dt)
     !     compute the LU factorization of the Newton system's Jacobian

     double precision, intent(in )    :: dt

     integer          :: i,j
     double precision :: z(Nspec+1)

     !     we are computing the Jacobian of the function (I - dt w)
     !     so the Jacobian is given by I - dt * chem Jac
     do j=1,Nspec+1
       do i=1,Nspec+1
          A(i,j) = -dt*Jac(i,j)
       end do
     end do

     do i=1,Nspec+1
        A(i,i) = (1.d0 + A(i,i))
     end do

     ! Compute transpose of NL system Jacobian
     AT = TRANSPOSE(A)

     !     call LINPACK to get the LU factorization of the matrix A
     !     call dgefa(A, Nspec+1, Nspec+1, ipvt, info)
     call DGECO(A, Nspec+1, Nspec+1, ipvt, rcond, z)

   end subroutine LUA
   !-----------------------!


   !-----------------------!
   subroutine update_vals(rY,T,rho,dx_nl)

       implicit none

       double precision, intent(inout  ) :: rY(nspec)
       double precision, intent(inout  ) :: T,rho
       double precision, intent(in  )    :: dx_nl(nspec+1)

       rY(:) = rY(:) + dx_nl(1:Nspec)        
       T     = T     + dx_nl(1+Nspec)
       rho   = sum(rY(:))

   end subroutine update_vals
   !-----------------------!


   !-----------------------!
   subroutine linesearch(rY,T,rho,res_nl,dt,dx_nl,flag)

       implicit none

       double precision, intent(inout  ) :: rY(nspec)
       double precision, intent(inout  ) :: T,rho
       double precision, intent(in  )    :: dx_nl(nspec+1)
       double precision, intent(in )     :: dt
       double precision, intent(in )     :: res_nl(nspec+1)
       integer, intent(out)              :: flag

       integer          :: count_linesearch
       double precision :: rho_inv
       double precision :: lambda, alpha
       double precision :: T_tmp, Tdot_tmp
       double precision :: Y_tmp(Nspec), rY_tmp(Nspec)
       double precision :: wdot_tmp(Nspec)
       double precision :: res_nl_tmp_norm, rmax_tmp
       double precision :: res_nl_tmp(Nspec+1), bound_norm(Nspec+1)

       logical          :: satisfied

       integer, parameter :: max_linesearch = 5

       flag = 0

       satisfied = .false.
       lambda = 1.0d0
       alpha = 1.0d-04
       count_linesearch = 0
       do while (.not. satisfied)
          rY_tmp(:) = rY(:) + lambda * dx_nl(1:Nspec)        
          T_tmp     = T     + lambda * dx_nl(1+Nspec)

          ! Compute rho and Y
          rho     = sum(rY_tmp(:))

          rho_inv = 1.0d0 / rho
          Y_tmp   = rY_tmp(:) * rho_inv

          ! compute wdot
          call compute_source_terms(Y_tmp,T_tmp,rho,wdot_tmp,Tdot_tmp)

          ! Compute residuals (a)
          call compute_residuals(rY_tmp,T_tmp,wdot_tmp,Tdot_tmp,dt,res_nl_tmp)

          ! compute norm of r + alpha * lambda * (AT:dx) (b)
          CALL compute_bound(alpha, lambda, dx_nl, res_nl_tmp, bound_norm)
          !CALL compute_bound(alpha, lambda, dx_nl, res_nl, bound_norm)
          
          ! min of (b) - (a)
          res_nl_tmp_norm = minval(bound_norm(:) - res_nl_tmp(:)) 

          ! check if a < b and if yes then satisfied if no update lambda
          if ( 0.5*NORM2(res_nl_tmp)**2 < 0.5*NORM2(bound_norm)**2) then
                  satisfied = .true.
                  rY(:)     = rY_tmp(:)
                  T         = T_tmp

                  if (verbose .ge. 6) then
                      print *,"       *INFO: number of linesearch steps is ", count_linesearch
                      print *,"       Armijo condition: ", lambda, res_nl_tmp_norm, 0.5*NORM2(res_nl_tmp)**2, 0.5*NORM2(bound_norm)**2
                      print *," "
                  end if

                  ! if we have reached the desired tolerance then we are done
                  rmax_tmp = maxval(abs(res_nl_tmp))
                  if (rmax_tmp .le. tol) then
                     if (verbose .ge. 4) then
                         print *,"       --> Newton has converged !! <--  " 
                         print *,"   +++++++++++++++++++++++++++++++++" 
                     end if
                     flag = 1 
                     return
                  endif
          else
                  lambda =  lambda*0.5
                  count_linesearch = count_linesearch + 1
                  if (verbose .ge. 6) then
                      print *,"       Armijo condition: ", lambda, res_nl_tmp_norm, 0.5*NORM2(res_nl_tmp)**2, 0.5*NORM2(bound_norm)**2
                  end if
          end if

          if (count_linesearch > max_linesearch) then
              if (verbose .ge. 6) then
                  print *,"       *Max linesearch reached !! " 
                  print *,'       linesearch: itmax =',max_linesearch
                  print *,'       linesearch: lambda =', lambda*2.0d0
              end if
              satisfied = .true.
              rY(:)     = rY_tmp(:)
              T         = T_tmp
          end if

       end do

   end subroutine linesearch
     
   subroutine compute_bound(a,l,dx_nl,res_nl_lcl,b)

         double precision, intent(in  ) :: a,l
         double precision, intent(in  ) :: res_nl_lcl(Nspec+1)
         double precision, intent(in  ) :: dx_nl(nspec+1)
         double precision, intent(out ) :: b(Nspec+1)

         double precision :: al
         double precision :: vect_tmp1(Nspec+1)

         al = a*l
         vect_tmp1 = MATMUL(AT,dx_nl)
         b(:) = res_nl_lcl(:) + al * vect_tmp1(:)

   end subroutine compute_bound


end module bechem_module
