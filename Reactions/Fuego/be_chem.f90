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

   !     do a Backward Euler solve for the chemistry using Newton's method
   !      Y : solution (output)
   !     Y0 : initial guess
   !   hmix : enthalpy (guess)
   !    rho : initial density
   !      T : initial temperature
   !    rhs : input right-hand side
   !     dt : timestep
   subroutine bechem(Y, Y0, hmix, rho, T_init, rhs, rhohdot_ext, dt)

     double precision, intent(out) :: Y(Nspec)
     double precision, intent(in ) :: Y0(Nspec)
     double precision, intent(inout ) :: hmix
     double precision, intent(in ) :: rho
     double precision, intent(in ) :: T_init
     double precision, intent(in ) :: rhs(Nspec+1)
     double precision, intent(in ) :: dt
     double precision, intent(in ) :: rhohdot_ext
     
     integer          :: iwrk, iter, n, ierr
     double precision :: rwrk, rmax, rho_inv, cp, cp_inv, T
     double precision :: rcond
     double precision, dimension(Nspec+1) :: r
     double precision, dimension(Nspec)   :: wdot, mwt, invmwt
     integer, parameter :: max_iter = 1000, NiterMAX = 400
     double precision, parameter :: tol = 1.d-16

     CALL t_bechem%start

     rho_inv = 1.d0/rho
   
     if (.not. allocated(A)) then
        allocate(Jac(Nspec+1,Nspec+1))
        allocate(A(Nspec+1,Nspec+1))
        allocate(ipvt(Nspec+1))
     end if
     
     ! start with the initial guess
     Y = Y0
     
     ! maximum number of iterations is max_iter
     do iter = 0, max_iter
        ! Newton's method: iteratively solve J(x_{n+1} - x_n) = -F(x_n)
        ! F(x) is given by (I - dt*wdot/rho - rhs/rho)
        
        ! get the temperature
        if (iter .eq. 0) then
            T = T_init
        end if
        call get_T_given_hY(hmix,Y,iwrk,rwrk,T,ierr)
        if (ierr/=0) then
          print *,'bechem: H to T solve failed'
          stop
        end if
        !write (*,*) "T ? ", T
        
        ! compute wdot
        call CKWYR(rho, T, Y, iwrk, rwrk, wdot)
        !write (*,*) "wdot ? ", wdot(:)
        ! compute C_p and 1/C_p
        call CKCPBS(T, Y, iwrk, rwrk, cp)
        cp_inv = 1.d0/cp
        !write (*,*) " cp ? ", cp
        
        ! multiply by molecular weight to get the right units
        call CKWT(iwrk, rwrk, mwt)
        do n=1,Nspec
           wdot(n) = wdot(n) * mwt(n)
        end do
        r(1:Nspec) = -(Y(:) - dt*wdot(:)*rho_inv - rhs(1:Nspec))
        r(Nspec+1) = -(hmix - dt*rhohdot_ext*rho_inv - rhs(Nspec+1))
        !write (*,*) "res sur h? ", r(Nspec+1), rhs(Nspec+1)

        rmax = maxval(abs(r))
        !write (*,*) "res sur Y ?", rmax

        if (isnan(rmax)) then
           print *,'wchem: backward Euler solve returned NaN'
           print *,'wchem: iteration: ', iter
           print *,'wchem: reciprocal condition number of J = ', rcond
           print *,'wchem: Cp = ', cp
           print *,'wchem: density = ', rho
           print *,'wchem: temperature = ', T
           
           print *,'wchem: Y0:',Y0
           stop
        endif

        ! if we have reached the desired tolerance then we are done
        if (rmax .le. tol) then
           !print *,'iters=',iter
           !print *,'I shall exit now '
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
        hmix = hmix + r(1+Nspec)


     end do

     if (iter .ge. max_iter) then
        print *,'wchem: Newton solve failed to converge'
        print *,'wchem: iter=',iter
        print *,'wchem: rmax=',rmax
        stop
     endif

     CALL t_bechem%stop

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
          Jac(i,j) = Jac(i,j) * invmwt(j) * rho 
       end do

       j = Nspec+1
       rho_inv = 1/rho
       do i=1,Nspec
          Jac(i,j) = Jac(i,j) * mwt(i) * rho_inv
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
