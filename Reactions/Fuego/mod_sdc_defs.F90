module mod_sdc_defs

      use amrex_fort_module, only : amrex_real

      implicit none

      integer                       :: iE, verbose, iStiff
      real(amrex_real)              :: rhohdot_ext, rhoedot_ext
      real(amrex_real), allocatable :: rhoh_init(:), rhoe_init(:)
      real(amrex_real), allocatable :: rhoydot_ext(:)
      double precision, allocatable :: mwt(:), invmwt(:)
      double precision, allocatable :: rhs(:)
      integer                       :: iwrk
      double precision              :: rwrk

end module mod_sdc_defs
