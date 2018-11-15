module main_module

  use amrex_fort_module, only : amrex_real
  use mod_timers   

  implicit none

  real(amrex_real), dimension(:,:), allocatable :: Y_in, Y_forcing_in
  real(amrex_real), dimension(:), allocatable :: temp
  real(amrex_real)                            :: pressure
  integer :: nlin

contains

    subroutine extern_init(name,namlen) &
                    bind(C, name="extern_init")

    use, intrinsic :: iso_c_binding
    use network
    use eos_module
    use transport_module
    use reactor_module

    integer :: namlen
    integer :: name(namlen)

    real (kind=dp_t) :: small_temp = 1.d-200
    real (kind=dp_t) :: small_dens = 1.d-200
    integer :: nLobato, nsdcite 

    ! initialize the external runtime parameters in
    ! extern_probin_module
    call runtime_init(name,namlen)

    call network_init()

    call eos_init(small_temp, small_dens)

    call transport_init()

    nLobato = 3
    nsdcite = 3
    call reactor_init_sdc(nLobato, nsdcite)

  end subroutine extern_init


  subroutine extern_close() bind(C, name="extern_close")

    use transport_module
    use reactor_module

    call transport_close()

    call reactor_close_sdc()

  end subroutine extern_close


  subroutine get_num_spec(nspec_out) bind(C, name="get_num_spec")

    use network, only : nspec

    implicit none

    integer, intent(out) :: nspec_out

    nspec_out = nspec

  end subroutine get_num_spec


  subroutine read_data_from_txt(name, namlen) &
          bind(C, name="read_data_from_txt")

    use network, only: nspec

    implicit none
    integer :: namlen
    integer :: name(namlen)
    ! Local var
    integer            :: i
    integer, parameter :: maxlen = 256
    character (len=maxlen) :: probin
    character(len=6)   :: a
    real(amrex_real)   :: y,x_velocity,y_velocity,density,dum

    ! create the filename
    if (namlen > maxlen) then
       print *, 'Initialization file name too long'
       stop
    endif
    do i = 1, namlen
       probin(i:i) = char(name(i))
    end do
    print *, "Initializing from output txt file ",probin(1:namlen)

    ! read in the file
    open (unit=49, file=probin(1:namlen), form='formatted', status='old')
    ! Start with number of useful lines
    read(49,*) nlin
    write(*,*) "  --> txt file has ", nlin, " lines"
    ! allocate stuff
    allocate(Y_in(nlin,nspec))
    allocate(Y_forcing_in(nlin,nspec+1))
    allocate(temp(nlin))
    ! ignore headers
    read(49,*) a
    read(49,*) a
    read(49,*) a
    ! read useful lines
    !y  x_velocity  y_velocity  density  rhoh  tracer  temp  RhoRT  divu  dsdt  FuncCount  Y(H2)  Y(H)  Y(O)  Y(O2)  Y(OH)  Y(H2O)  Y(HO2)  Y(CH2)  Y(CH2(S))  Y(CH3)  Y(CH4)  Y(CO)  Y(CO2)  Y(HCO)  Y(CH2O)  Y(CH3O)  Y(C2H4)  Y(C2H5)  Y(C2H6)  Y(N2)  Y(AR)  CH4_ConsumptionRate  YH2_ForcingTerm  YH_ForcingTerm  YO_ForcingTerm  YO2_ForcingTerm  YOH_ForcingTerm  YH2O_ForcingTerm  YHO2_ForcingTerm  YCH2_ForcingTerm  YCH2(S)_ForcingTerm  YCH3_ForcingTerm  YCH4_ForcingTerm  YCO_ForcingTerm  YCO2_ForcingTerm  YHCO_ForcingTerm  YCH2O_ForcingTerm  YCH3O_ForcingTerm  YC2H4_ForcingTerm  YC2H5_ForcingTerm  YC2H6_ForcingTerm  YN2_ForcingTerm  YAR_ForcingTerm  Temperature_ForcingTerm  HeatRelease
    DO i = 1, nlin
      read(49,*) y,x_velocity,y_velocity,density,dum,dum,temp(i),dum,dum,dum,dum,Y_in(i,:),dum,Y_forcing_in(i,:)
      !read(49,*) y,x_velocity,y_velocity,density,dum,dum,temp(i),dum,dum,dum,dum,Y_in(i,:),dum,dum,dum,dum,dum,dum,Y_forcing_in(i,:)
      print *, sum(Y_in(i,:))
    END DO

    close (unit=49)

    ! Conversion MKS to CGS for PelePhys
    !! nspec + 1 is energy which here is enthalpy !!
    DO i = 1, nlin
      Y_forcing_in(i,1:nspec) = Y_forcing_in(i,1:nspec)*1.d-3
      Y_forcing_in(i,nspec+1) = Y_forcing_in(i,nspec+1)*10.0
    END DO

  end subroutine read_data_from_txt


  subroutine initialize_data( &
       lo,           hi, &
       rhoY,         rY_lo, rY_hi, &
       rhoY_src,     rY_src_lo, rY_src_hi) &
       bind(C, name="initialize_data")

    use amrex_constants_module, only: M_PI, HALF, ONE, TWO, ZERO
    use network, only: nspec
    use eos_type_module
    use eos_module

    implicit none
    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  rY_lo(3), rY_hi(3)
    integer         , intent(in   ) ::  rY_src_lo(3), rY_src_hi(3)
    real(amrex_real), intent(inout) ::  rhoY(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspec+2)
    real(amrex_real), intent(inout) ::  rhoY_src(rY_src_lo(1):rY_src_hi(1),rY_src_lo(2):rY_src_hi(2),rY_src_lo(3):rY_src_hi(3),nspec+1)

    ! local variables
    integer          :: i, j, k
    !real(amrex_real) :: pressure
    type(eos_t) :: eos_state

    call build(eos_state)

    ! CGS UNITS
    pressure = 1013250.d0

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             eos_state % p          = pressure
             eos_state % T          = temp(i+1)
             eos_state % massfrac(:)     = Y_in(i+1,:)
             eos_state % massfrac(nspec) = ONE - sum(Y_in(i+1,1:nspec-1))

             call eos_tp(eos_state)

             ! rhoY(:nspec) = rhoY, rhoY(nspec+1) = nrg, rhoY(nspec+2) = T
             rhoY(i,j,k,1:nspec) = eos_state % massfrac * eos_state % rho
             ! full enthalpy mode
             rhoY(i,j,k,nspec+1) = eos_state % h

             rhoY(i,j,k,nspec+2) = eos_state % T

             ! rhoY_src(:nspec) = rhoForcingSpecs, rhoY_src(nspec+1) = rhoForcingNRG
             rhoY_src(i,j,k,1:nspec) = Y_forcing_in(i+1,1:nspec)
             rhoY_src(i,j,k,nspec+1) = Y_forcing_in(i+1,nspec+1)

          end do
       end do
    end do

    call destroy(eos_state)

  end subroutine initialize_data


  subroutine initialize_data_byhand( &
       lo,hi, &
       rhoY,         rY_lo, rY_hi, &
       dx, plo, phi) &
       bind(C, name="initialize_data_byhand")

    use amrex_constants_module, only: M_PI, HALF, ONE, TWO, ZERO
    use network, only: nspec
    use eos_type_module
    use eos_module

    implicit none

    integer         , intent(in   ) ::     lo(3),    hi(3)
    integer         , intent(in   ) ::  rY_lo(3), rY_hi(3)
    real(amrex_real), intent(in   ) ::     dx(3)
    real(amrex_real), intent(in   ) ::     plo(3),  phi(3)
    real(amrex_real), intent(inout) ::  rhoY(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3),nspec+2)

    ! local variables
    integer            :: i, j, k, ii
    real(amrex_real)   :: y
    type(eos_t)        :: eos_state
    real(amrex_real)   :: dum
    character(len=6)   :: a

    CALL t_eos%start
    call build(eos_state)
    CALL t_eos%stop

    CALL t_readData%init("Read Data")   
    CALL t_readData%start

    ! read in the file
    open (unit=49, file="datafromSC.dat", form='formatted', status='old')
    ! allocate stuff
    allocate(Y_in(1,nspec))
    allocate(temp(1))
    read(49,*) a
    read(49,*) pressure,temp(1),dum,Y_in(1,:)
    eos_state % molefrac(:) = Y_in(1,:)
    pressure = pressure*10.d0
    print *, "data read from datafromSC.dat ", pressure,temp(1),eos_state % molefrac(:)
    print *, "sum mole frac ", sum(eos_state % molefrac(:))
    close (unit=49)


    call eos_xty(eos_state)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             eos_state % p               = pressure
             eos_state % T               = temp(1) 
             eos_state % massfrac(nspec) = ONE - sum(eos_state % massfrac(1:nspec-1))

             call eos_tp(eos_state)

             rhoY(i,j,k,1:nspec) = eos_state % massfrac * eos_state % rho
             ! full enthalpy mode
             !stop 'Full H mode not allowed when data init by hand'
             rhoY(i,j,k,nspec+1) = eos_state % h

             rhoY(i,j,k,nspec+2) = eos_state % T

          end do
       end do
    end do

    call destroy(eos_state)

    CALL t_readData%stop

  end subroutine initialize_data_byhand


  subroutine react_state(lo,hi, &
                         mold,mo_lo,mo_hi, &
                         mnew,mn_lo,mn_hi, &
                         ysrc,ys_lo,ys_hi, &
                         mask,m_lo,m_hi, &
                         cost,c_lo,c_hi, &
                         time,dt_react ) bind(C, name="react_state")

    use network           , only : nspec
    use react_type_module
    use reactor_module
    use react_type_module
    use, intrinsic :: iso_c_binding

    implicit none

    integer         , intent(in   ) ::    lo(3),    hi(3)
    integer         , intent(in   ) :: mo_lo(3), mo_hi(3)
    integer         , intent(in   ) :: mn_lo(3), mn_hi(3)
    integer         , intent(in   ) :: ys_lo(3), ys_hi(3)
    integer         , intent(in   ) ::  m_lo(3),  m_hi(3)
    integer         , intent(in   ) ::  c_lo(3),  c_hi(3)
    real(amrex_real), intent(inout)  :: mold(mo_lo(1):mo_hi(1),mo_lo(2):mo_hi(2),mo_lo(3):mo_hi(3),nspec+2)
    real(amrex_real), intent(inout)  :: mnew(mn_lo(1):mn_hi(1),mn_lo(2):mn_hi(2),mn_lo(3):mn_hi(3),nspec+2)
    real(amrex_real), intent(inout)  :: ysrc(ys_lo(1):ys_hi(1),ys_lo(2):ys_hi(2),ys_lo(3):ys_hi(3),nspec+1)
    integer, intent(inout)           :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    real(amrex_real), intent(inout)  :: cost(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(c_double), intent(in)       :: time, dt_react

    ! local variables
    integer          :: i, j, k, ii
    real(amrex_real) :: rho
    integer          :: istop 
    real(c_double)   :: dt_react_tmp, dt_react_incr, time_tmp

    type (react_t)          :: react_state_in !, react_state_out
    type (reaction_stat_t)  :: stat

    call build(react_state_in)

    istop = 0
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

                write(*,*) ""
                write(*,*) "Dealing with cell ", i,j,k

                react_state_in %              h = mold(i,j,k,nspec+1)
                react_state_in %    rhohdot_ext = ysrc(i,j,k,nspec+1)
                react_state_in %              T = mold(i,j,k,nspec+2)
                react_state_in %        rhoY(:) = mold(i,j,k,1:nspec)
                react_state_in %            rho = sum(react_state_in % rhoY(:))
                react_state_in % rhoYdot_ext(:) = ysrc(i,j,k,1:nspec)
                react_state_in % i = i
                react_state_in % j = j
                react_state_in % k = k
                react_state_in % p = pressure

                if (dt_react < 0.0) then
                        time_tmp = time
                        dt_react_incr =  - dt_react / 3000.0
                        !if (react_state_in % T > 1500.0) then
                        !    istop = 1
                            write(12,*) "#dt_react_incr ", dt_react_incr
                            write(12,*) "#time, T, h, rho, Yks in cell "
                        !end if
                        do ii= 1, 3000
                            dt_react_tmp   = ii* dt_react_incr
                            stat           = react_sdc(react_state_in, react_state_in, dt_react_incr, time)
                            time_tmp       = dt_react_tmp
                            rho            = sum(react_state_in % rhoY(1:nspec))
                            !if (istop.eq.1) then
                                write(12,*) time_tmp, react_state_in % T, react_state_in % h, rho, react_state_in % rhoY(1:nspec)/rho
                            !end if
                        end do
                else
                        stat = react_sdc(react_state_in, react_state_in, dt_react, time)
                end if
                !cost(i,j,k) = stat % cost_value

                ! Export e whenever
                mnew(i,j,k,nspec+1)             = react_state_in % e
                mnew(i,j,k,nspec+2)             = react_state_in % T
                mnew(i,j,k,1:nspec)             = react_state_in % rhoY(1:nspec)
                rho                             = sum(mnew(i,j,k,1:nspec))
                mnew(i,j,k,1:nspec)             = mnew(i,j,k,1:nspec)/rho
                write(*,*) "T in cell ", react_state_in % T
                write(*,*) "Y in cell ", mnew(i,j,k,1:nspec)
                !if (istop.eq.1) then
                !        stop
                !end if

                CALL t_total%stop
                CALL t_total%norm(t_total)
                CALL t_bechem%norm(t_total)
                CALL t_eos%norm(t_total)
                CALL t_ck%norm(t_total)
                if ((i == 0).and.(j==0).and.(k==0)) then
                    CALL t_readData%norm(t_total)
                end if
                !PRINTS
                CALL t_total%print
                CALL t_bechem%print
                CALL t_eos%print
                CALL t_ck%print
                if ((i == 0).and.(j==0).and.(k==0))  then
                    CALL t_readData%print
                end if
                !RESETS
                CALL t_total%reset
                CALL t_bechem%reset
                CALL t_eos%reset
                CALL t_ck%reset

          end do
       enddo
    enddo

    call destroy(react_state_in)

  end subroutine react_state



end module main_module
