module bechem_module

   use network, only: nspec
   use mod_timers
   !use actual_sdc_module, only: rhohdot_ext
   
   implicit none
   private

   !     Jacobian matrix and factorization
   double precision, allocatable, save :: Jac(:,:), A(:,:)
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
     double precision :: rcond
     double precision, dimension(Nspec+1) :: r
     double precision, dimension(Nspec)   :: wdot, mwt, invmwt
     double precision, dimension(Nspec)   :: hi(Nspec), ei(Nspec)
     integer, parameter :: max_iter = 1000, NiterMAX = 400
     double precision, parameter :: tol = 1.d-07

     !CALL t_bechem%start

     rho_inv = 1.d0/rho
   
     if (.not. allocated(A)) then
        allocate(Jac(Nspec+1,Nspec+1))
        allocate(A(Nspec+1,Nspec+1))
        allocate(ipvt(Nspec+1))
     end if
     
     ! start with the initial guess
     Y = Y0
    
     print *,"   +++++++++++++++++++++++++++++++++" 
     print *,"     STARTING NEWTON ITERATIONS     "
     !print *,"   +++++++++++++++++++++++++++++++++" 
     ! maximum number of iterations is max_iter
     do iter = 0, max_iter
        ! Newton's method: iteratively solve J(x_{n+1} - x_n) = -F(x_n)
        ! F(x) is given by (I - dt*wdot/rho - rhs/rho)
        
        ! get the temperature
        if (iter .eq. 0) then
            T = T_init
            if (iE == 1) then
                call get_T_given_eY(hmix,Y,iwrk,rwrk,T,ierr)
            else
                call get_T_given_hY(hmix,Y,iwrk,rwrk,T,ierr)
            end if
            if (ierr/=0) then
              print *,'bechem: H/E to T solve failed'
              stop
            end if
            !write (*,*) "T ? ", T
        end if
        print *,"     Working on the ", iter, " iteration (T, Y(O2)) ", T, Y(8)

        
        ! compute wdot
        call CKWYR(rho, T, Y, iwrk, rwrk, wdot)

        ! compute spec energy
        ! compute C_p/v 
        if (iE == 1) then
            call CKUMS (T, iwrk, rwrk, ei)
            call CKCVBS(T, Y, iwrk, rwrk, cv)
            cv_inv = 1.d0/cv
        else 
            call CKHMS (T, iwrk, rwrk, hi)
            call CKCPBS(T, Y, iwrk, rwrk, cp)
            cp_inv = 1.d0/cp
        end if
        
        ! multiply by molecular weight to get the right units
        call CKWT(iwrk, rwrk, mwt)
        if (iE == 1) then
            ! it is called rhoh but it is e
            Tdot = rhohdot_ext 
            do n=1,Nspec
                Tdot = Tdot - ei(n) * wdot(n) 
            end do
            Tdot = Tdot * cv_inv * rho_inv
        else 
            Tdot = rhohdot_ext 
            do n=1,Nspec
                Tdot = Tdot - hi(n) * wdot(n) 
            end do
            Tdot = Tdot * cp_inv * rho_inv
        end if
        wdot(:) = wdot(:) * mwt(:) + rhoydot_ext(:)
        r(1:Nspec) = -(Y(:) - dt*wdot(:)*rho_inv - rhs(1:Nspec))
        r(Nspec+1) = -(T - dt * Tdot - rhs(Nspec+1))
        !write (*,*) "res sur h? ", r(Nspec+1), rhs(Nspec+1)

        rmax = maxval(abs(r))
        print *,"     Max residual, T residual: ", rmax, r(Nspec+1)

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
        if (rmax .le. tol) then
           !print *,'iters=',iter
           !print *,'I shall exit now '
           print *,"       --> Newton has converged !! <--  " 
           print *,"   +++++++++++++++++++++++++++++++++" 
           exit
        endif

        !     compute the Jacobian, and then compute its LU factorization
        call LUA
        !     call LINPACK to solve the linear system Ax = r
        !     using the output from the factorization given by dgefa
        call dgesl(A, Nspec+1, Nspec+1, ipvt, r, 0)
        
        !     solve for the difference x = x_{n+1} - x_n
        !     so the next iterate is given by x_{n+1} = x_n + x

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
        
        Y(:) = Y(:) + r(1:Nspec)        
        T = T + r(1+Nspec)


     end do

     if (iter .ge. max_iter) then
        print *,'wchem: Newton solve failed to converge'
        print *,'wchem: iter=',iter
        print *,'wchem: rmax=',rmax
        stop
     endif

     !CALL t_bechem%stop

   contains

     !     compute the LU factorization of the Jacobian
     subroutine LUA
       integer :: i, j, iwrk     !, info
       double precision :: rwrk
       double precision, dimension(Nspec) :: C
       integer, parameter :: consP = 1
       double precision :: z(Nspec+1)

       !    compute the molar concentrations
       call CKYTCR(rho, T, Y, iwrk, rwrk, C)
       !     use the concentrarions to compute the reaction Jacobian
       call DWDOT(Jac, C, T, consP)

       !     convert to the appropriate units
       call CKWT(iwrk, rwrk, mwt)
       invmwt(:) = 1.0/mwt(:)
       do j=1,Nspec
          do i=1,Nspec
             Jac(i,j) = Jac(i,j) * mwt(i) * invmwt(j)
          end do
          i=Nspec+1
          Jac(i,j) = Jac(i,j) * invmwt(j) !* rho 
       end do

       j = Nspec+1
       rho_inv = 1/rho
       do i=1,Nspec
          Jac(i,j) = Jac(i,j) * mwt(i) !* rho_inv
       end do

       !     we are computing the Jacobian of the function (I - dt w)
       !     so the Jacobian is given by I - dt Jac
       do j=1,Nspec+1
         do i=1,Nspec+1
            A(i,j) = -dt*Jac(i,j)
         end do
       end do

       do i=1,Nspec+1
          A(i,i) = 1.d0 + A(i,i)
       end do

       !     call LINPACK to get the LU factorization of the matrix A
       !     call dgefa(A, Nspec+1, Nspec+1, ipvt, info)
       call DGECO(A, Nspec+1, Nspec+1, ipvt, rcond, z)
       !write(*,*) 'condit ? ', rcond

     end subroutine LUA

   end subroutine bechem

end module bechem_module
