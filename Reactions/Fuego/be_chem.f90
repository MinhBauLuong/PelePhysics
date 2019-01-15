module bechem_module

   use network, only: nspec
   use mod_timers
   !use actual_sdc_module, only: rhohdot_ext
   
   implicit none
   private

   !     Jacobian matrix and factorization
   double precision, allocatable, save :: Jac(:,:), A(:,:)
   double precision, allocatable, save :: JacT(:,:), AT(:,:)
   !     Pivot vector returned from LINPACK
   integer, allocatable, save :: ipvt(:)

   public :: bechem

contains

   !       do a Backward Euler solve for the chemistry using Newton's method
   !         Y : solution (output)
   !        Y0 : initial guess
   !      hmix : enthalpy (guess) // OR internal energy (guess)
   !       rho : initial density
   !         T : initial temperature
   !       rhs : input right-hand side
   !rhohdot_ext: ext term source of energy (h OR e)
   !rhoydot_ext: ext term source of species
   ! NB: we should be able to get rid of these ! 
   !        dt : timestep
   !        iE : 1 (UV) // 2 (rhoH)
   subroutine bechem(Y, Y0, hmix, rho, T_init, rhs, rhohdot_ext, rhoydot_ext, dt, iE)

     double precision, intent(out) :: Y(Nspec)
     double precision, intent(in ) :: Y0(Nspec)
     double precision, intent(inout ) :: hmix
     double precision, intent(in ) :: rho
     double precision, intent(in ) :: T_init
     double precision, intent(in ) :: rhs(Nspec+1)
     double precision, intent(in ) :: rhohdot_ext
     double precision, intent(in ) :: rhoydot_ext(Nspec)
     double precision, intent(in ) :: dt
     integer, intent(in )          :: iE
     
     integer          :: iwrk, iter, n, ierr
     double precision :: rwrk, rmax, rho_inv, T
     double precision :: cp, cp_inv, cv, cv_inv, Tdot
     double precision :: rcond, dum
     double precision :: res_nl_norm
     double precision, dimension(Nspec+1) :: res_nl,res_nl_init
     double precision, dimension(Nspec+1) :: res_nl_divided_val_ref
     double precision, dimension(Nspec+1) :: dx_nl
     double precision, dimension(Nspec)   :: wdot, mwt, invmwt
     double precision, dimension(Nspec)   :: hi, ei
     integer, parameter :: max_iter = 1000, NiterMAX = 400
     double precision, parameter :: tol = 1.d-10
     double precision, dimension(Nspec+1) :: val_ref
     logical          :: recompute_jac

     !CALL t_bechem%start
     ! calc val ref

     ! Need molecular weights
     call CKWT(iwrk, rwrk, mwt)
     invmwt(:) = 1.0/mwt(:)

     ! assume rho does not change
     rho_inv = 1.d0/rho
   
     ! Allocate space for Jac, factorization and pivot
     if (.not. allocated(A)) then
        allocate(Jac(Nspec+1,Nspec+1))
        !allocate(JacT(Nspec+1,Nspec+1))
        allocate(A(Nspec+1,Nspec+1))
        allocate(AT(Nspec+1,Nspec+1))
        allocate(ipvt(Nspec+1))
     end if
     
     ! start with the initial guess
     Y = Y0
     T = T_init

     ! Initially assume Jac is wrong
     recompute_jac = .true.
    
     print *,"   +++++++++++++++++++++++++++++++++" 
     print *,"     STARTING NEWTON ITERATIONS     "

     ! maximum number of iterations is max_iter
     do iter = 0, max_iter
        ! Newton's method: iteratively solve J(x_{n+1} - x_n) = -F(x_n)
        ! F(x) is given by (I - dt*wdot/rho - rhs/rho)
        
        print *,"     Working on the ", iter, " iteration (T, Y(O2)) ", T, Y(8)

        ! get the temperature
        if (iter .eq. 0) then

            ! recompute T to be sure
            ! Not sure we need this though
            !if (iE == 1) then
            !    call get_T_given_eY(hmix,Y,iwrk,rwrk,T,ierr)
            !else
            !    call get_T_given_hY(hmix,Y,iwrk,rwrk,T,ierr)
            !end if
            !if (ierr/=0) then
            !  print *,'bechem: H/E to T solve failed'
            !  stop
            !end if

            ! compute wdot
            call CKWYR(rho, T, Y, iwrk, rwrk, wdot)

            !call CKPY(rho, T, Y, iwrk, rwrk, dum)
            !print *, "P is ? ", dum
            !stop

            ! compute C_p/v 
            if (iE == 1) then
                call CKCVBS(T, Y, iwrk, rwrk, cv)
                cv_inv = 1.d0/cv
            else 
                call CKCPBS(T, Y, iwrk, rwrk, cp)
                cp_inv = 1.d0/cp
            end if
        
            ! compute T&spec src terms
            if (iE == 1) then
                call CKUMS (T, iwrk, rwrk, ei)
                ! it is called rhoh but it is e
                Tdot = rhohdot_ext 
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
            ! multiply by molecular weight to get the right units
            wdot(:) = wdot(:) * mwt(:) + rhoydot_ext(:)

            ! Compute initial residuals
            res_nl_init(1:Nspec) = -(Y(:) - dt*wdot(:)*rho_inv - rhs(1:Nspec))
            res_nl_init(Nspec+1) = -(T - dt * Tdot - rhs(Nspec+1))

            res_nl = res_nl_init
            !print *, " -- > res_nl_init ", res_nl_init
            !print *, " "

!            val_ref(:) = 1.0/(res_nl_init(:) + 1.0d-10) + 1.0d-10*res_nl_init(:)
            val_ref(:) = 1.0
        end if

        rmax = maxval(abs(res_nl))
        res_nl_norm = 0.5*NORM2(res_nl(:)*val_ref(:))**2.0
        print *,"     L2 residual, Max residual: ", res_nl_norm, rmax

        if (isnan(rmax)) then
           print *," "
           print *,"     BE solve returned NaN !! " 
           print *," "
           print *,'     bechem: iteration: ', iter
           print *,'     bechem: reciprocal condition number of J = ', rcond
           print *,'     bechem: Cp = ', cp
           print *,'     bechem: density = ', rho
           print *,'     bechem: temperature = ', T
           print *,'     bechem: Y0:',Y0
           print *,"   +++++++++++++++++++++++++++++++++" 
           stop
        endif

        ! if we have reached the desired tolerance then we are done
        !if (res_nl_norm .le. tol) then
        if (rmax .le. tol) then
           print *,"       --> Newton has converged !! <--  " 
           print *,"   +++++++++++++++++++++++++++++++++" 
           exit
        endif

        !     compute the Jacobian, and then compute its LU factorization
        call LUA
        !     call LINPACK to solve the linear system Ax = res
        !     using the output from the factorization given by dgefa
        dx_nl = res_nl
        call dgesl(A, Nspec+1, Nspec+1, ipvt, dx_nl, 0)
        !     in dx_nl now is delta_x and res_nl is the Newton residual

        ! if (maxval(abs(r)) > max_change) then
        !    m = -1
        !    do n = 1,Nspec
        !       if (abs(r(n)) > max_change) then
        !          m = n
        !       end if
        !    end do
        !    print *, r(m)
        !    r = max_change * r / abs(r(m))
        ! end if
        
        !     solve for the difference x = x_{n+1} - x_n
        !     so the next iterate is given by x_{n+1} = x_n + x
        print *,"     Calling a linesearch"
        call linesearch

        !Y(:) = Y(:) + dx_nl(1:Nspec)        
        !T = T + dx_nl(1+Nspec)
     end do

     if (iter .ge. max_iter) then
        print *,"     Newton solve failed to converge !! " 
        print *,'     bechem: iter=',iter
        print *,'     bechem: rmax=',rmax
        print *,'     bechem: rL2=',res_nl_norm
        stop
     endif

     !CALL t_bechem%stop

   contains

     !     compute the LU factorization of the Jacobian
     subroutine LUA
       integer :: i, j, iwrk     !, info
       double precision :: rwrk
       double precision, dimension(Nspec) :: C
       integer, parameter :: consP = 0
       double precision :: z(Nspec+1)

       if (recompute_jac) then
           !    compute the molar concentrations
           call CKYTCR(rho, T, Y, iwrk, rwrk, C)
           !     use the concentrarions to compute the reaction Jacobian
           call DWDOT(Jac, C, T, consP)

           !     convert to the appropriate units
           do j=1,Nspec
              do i=1,Nspec
                 Jac(i,j) = Jac(i,j) * mwt(i) * invmwt(j)
                ! JacT(j,i) = Jac(i,j) 
              end do
              i=Nspec+1
              Jac(i,j) = Jac(i,j) * invmwt(j) * rho 
              !JacT(j,i) = Jac(i,j) 
           end do

           j = Nspec+1
           do i=1,Nspec
              Jac(i,j) = Jac(i,j) * mwt(i) * rho_inv
              !JacT(j,i) = Jac(i,j) 
           end do

           !JacT(j,j) = Jac(j,j)

       end if

       !     we are computing the Jacobian of the function (I - dt w)
       !     so the Jacobian is given by I - dt Jac
       do j=1,Nspec+1
         do i=1,Nspec+1
            A(i,j) = -dt*Jac(i,j)
            !AT(i,j) = -dt*JacT(i,j)
         end do
       end do

       do i=1,Nspec+1
          A(i,i) = (1.d0 + A(i,i))
          !AT(i,i) = (1.d0 + AT(i,i))
       end do

       ! Compute transpose of NL system Jacobian
       AT = TRANSPOSE(A)

       !     call LINPACK to get the LU factorization of the matrix A
       !     call dgefa(A, Nspec+1, Nspec+1, ipvt, info)
       call DGECO(A, Nspec+1, Nspec+1, ipvt, rcond, z)
       !write(*,*) 'condit ? ', rcond

     end subroutine LUA


     subroutine linesearch 

         double precision :: lambda, alpha
         double precision :: T_tmp, Tdot_tmp
         double precision :: res_nl_tmp_norm
         double precision :: Y_tmp(Nspec)
         double precision :: wdot_tmp(Nspec)
         double precision :: res_nl_tmp(Nspec+1), bound_norm(Nspec+1)
         logical          :: satisfied
         integer, parameter :: max_linesearch = 5
         integer            :: count_linesearch


         satisfied = .false.
         lambda = 1.0d0
         ! alpha should be small... like 1e-4 ??
         alpha = 1.0d-04
         count_linesearch = 0
         do while (.not. satisfied)
             Y_tmp(:) = Y(:) + lambda * dx_nl(1:Nspec)        
             T_tmp = T + lambda * dx_nl(1+Nspec)

             ! compute wdot
             call CKWYR(rho, T_tmp, Y_tmp, iwrk, rwrk, wdot_tmp)

             ! compute C_p/v 
             if (iE == 1) then
                 call CKCVBS(T_tmp, Y_tmp, iwrk, rwrk, cv)
                 cv_inv = 1.d0/cv
             else 
                 call CKCPBS(T_tmp, Y_tmp, iwrk, rwrk, cp)
                 cp_inv = 1.d0/cp
             end if
        
             ! Compute src terms 
             if (iE == 1) then
                 ! it is called rhoh but it is e
                 call CKUMS (T_tmp, iwrk, rwrk, ei)
                 Tdot_tmp = rhohdot_ext 
                 do n=1,Nspec
                     Tdot_tmp = Tdot_tmp - ei(n) * wdot_tmp(n) 
                 end do
                 Tdot_tmp = Tdot_tmp * cv_inv * rho_inv
             else 
                 call CKHMS (T_tmp, iwrk, rwrk, hi)
                 Tdot_tmp = rhohdot_ext 
                 do n=1,Nspec
                     Tdot_tmp = Tdot_tmp - hi(n) * wdot_tmp(n) 
                 end do
                 Tdot_tmp = Tdot_tmp * cp_inv * rho_inv
             end if
             ! multiply by molecular weight to get the right units
             wdot_tmp(:) = wdot_tmp(:) * mwt(:) + rhoydot_ext(:)

             ! Compute residuals (a)
             res_nl_tmp(1:Nspec) = -(Y_tmp(:) - dt * wdot_tmp(:) * rho_inv - rhs(1:Nspec))
             res_nl_tmp(Nspec+1) = -(T_tmp - dt * Tdot_tmp - rhs(Nspec+1))
             ! compute norm of res_nl_tmp 
             !res_nl_tmp_norm = 0.5*NORM2(res_nl_tmp(:)*val_ref(:))**2.0

             ! compute norm of r + alpha * lambda * (AT:dx) (b)
             CALL compute_bound(alpha, lambda, res_nl, bound_norm)
             !print *,"       New residual: ", res_nl_tmp(:)
             !print *,"       Norm of new residu: ", 0.5*NORM2(res_nl_tmp)**2
             
             ! min of (b) - (a)
             res_nl_tmp_norm = minval(bound_norm(:) - res_nl_tmp(:)) 
             !res_nl_tmp_norm = maxval(bound_norm(:) - res_nl_tmp(:)) 
             !print *, "bound_norm(:) - res_nl_tmp(:) ", (bound_norm(:) - res_nl_tmp(:))
             !print *, "Out condition ? ", res_nl_tmp_norm

             ! check if a < b and if yes then satisfied if no update lambda
             !if ( res_nl_tmp_norm > -0.000000001) then
             if ( 0.5*NORM2(res_nl_tmp)**2 < 0.5*NORM2(bound_norm)**2) then
                     satisfied = .true.
                     res_nl(:) = res_nl_tmp(:)
                     !dx_nl(:)  = lambda * dx_nl(:)
                     Y(:)      = Y_tmp(:)
                     T         = T_tmp
                     print *,"       *INFO: number of linesearch steps is ", count_linesearch
                     print *,"       Armijo condition: ", lambda, res_nl_tmp_norm, 0.5*NORM2(res_nl_tmp)**2, 0.5*NORM2(bound_norm)**2
                     print *," "
             else
                     lambda =  lambda*0.5
                     count_linesearch = count_linesearch + 1
                     print *,"       Armijo condition: ", lambda, res_nl_tmp_norm, 0.5*NORM2(res_nl_tmp)**2, 0.5*NORM2(bound_norm)**2
             end if

             if (count_linesearch > max_linesearch) then
                 print *,"       *Max linesearch reached !! " 
                 print *,'       linesearch: itmax =',max_linesearch
                 print *,'       linesearch: lambda =', lambda*2.0d0
                 satisfied = .true.
                 res_nl(:) = res_nl_tmp(:)
                 Y(:)      = Y_tmp(:)
                 T         = T_tmp
             end if

         end do

     end subroutine linesearch
     
     subroutine compute_bound(a,l,res_nl_lcl,b)

         double precision :: a,l
         double precision :: res_nl_lcl(Nspec+1), b(Nspec+1)
         double precision :: al
         double precision :: vect_tmp1(Nspec+1)
         integer          :: i

         al = a*l
         vect_tmp1 = MATMUL(AT,dx_nl)
         b(:) = res_nl_lcl(:) + al * vect_tmp1(:)
         !CHECKS
         !do i=1,Nspec+1
         !    dum = AT(Nspec+1,i) * dx_nl(i)
         !end do
         !print *," dum, vect_tmp1(Nspec+1) ", dum, vect_tmp1(Nspec+1) 
         !print *," al*dum, al*vect_tmp1(Nspec+1) ", al*dum, al*vect_tmp1(Nspec+1) 
         !print *,"  res_nl_lcl(Nspec+1) ", res_nl_lcl(Nspec+1) 
         !print *," res_nl_lcl(Nspec+1) + al*dum, b(Nspec+1) ", res_nl_lcl(Nspec+1) + al*dum, b(Nspec+1) 
         !stop
         !print *,"       delta_x: ", dx_nl(:) 
         !print *," "
         !print *,"       Obj func slope : ", vect_tmp1
         !print *,"       Bound: ",  b(:) 
         !print *,"       Norm of bound: ", 0.5*NORM2(b)**2
         !print *," "

     end subroutine compute_bound

   end subroutine bechem

end module bechem_module
