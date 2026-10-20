      module eirmod_hecr
cdr: comments....?
cdr: notational cleanup.
cdr: separation of atomic structure from rate-coefficients
cdr:  form I --> form I condensed to form II --> form. II
cdr:   .._f1 -->                      .._f1c -->  .._f2
cdr: added: rr12, rr13, rr02, rr03,  for f1c sanity checks f1 vs. f2

      implicit none
     
      private

      public :: atomic_structure, rate_coefficient,
     .          population_coefficient,
     .          cr_rate_coefficient, cooling_rate_coefficient,
     .          eirene_init_rate_arrays, eirene_init_constants, 
     .          eirene_init_popcoe,
     .          eirene_read_ralchenko_collision_data,
     .          eirene_read_ralchenko_ionization_data,
     .          ral_col_liste, ral_ion_liste
      public :: EV2J, J2EV, CM2J, J2CM, CM2EV, EV2CM, J2K, EV2K, CM2K

cdr   integer, parameter, public :: ULH = 65  for testing: remove LTE fringe
      integer, parameter, public :: ULH = 59
      integer, parameter, public :: ULL = 59

      type atomic_structure
cdr needs to be set only once, not per cell
        real*8 :: e(ULH)      ! atomic structure, 
        real*8 :: en(ULH)     ! atomic structure, 
        real*8 :: g(ULH)      ! atomic structure, 
        real*8 :: f(ULH,ULH)  ! atomic structure. Oszillator strength
        real*8 :: a(ULH,ULH)  ! atomic structure, A_ik
      end type atomic_structure

      type rate_coefficient
cdr single step rate coefficients, and Saha-Boltzmann coeff. vs. Te
        real*8 :: s(ULH)      ! rate coeff. ionization 
        real*8 :: alpha(ULH)  ! rate coeff. threebody rec.
        real*8 :: beta(ULH)   ! rate coeff. rad. rec
        real*8 :: betad(ULH)  ! rate coeff. rec. di-electronic
        real*8 :: c(ULH,ULH)  ! rate coeff. excit, de-exit
        real*8 :: z(ULH)      ! Saha-Boltzmann population at Te, 
      end type rate_coefficient

      type population_coefficient
cdr  form. I, ms resolved
	real*8 :: r0(ULH)
	real*8 :: r1(ULH)
	real*8 :: r2(ULH)
	real*8 :: r3(ULH)
	real*8 :: r_ext(ULH) ! R_EXT for metastable resolved
cdr  metastable equil. populations, for sanity checks
        real*8 :: rr12,rr13,rr02,rr03
cdr  form II, ms unresolved
	real*8 :: rr0(ULH)   ! R_0 
	real*8 :: rr1(ULH)   ! R_1 
	real*8 :: rr_ext(ULH)! R_EXT for metastable unresolved 
      end type population_coefficient

      type cr_rate_coefficient
cdr  form. I,  MS resolved
	real*8 :: k01  
	real*8 :: k1   
	real*8 :: k21
	real*8 :: k31
	real*8 :: k02
	real*8 :: k12
	real*8 :: k2
	real*8 :: k32
	real*8 :: k03
	real*8 :: k13
	real*8 :: k23
	real*8 :: k3
!  form. I condensed to form. II (for testing):
	real*8 :: scr_f1c   
	real*8 :: alpcr_f1c
	real*8 :: scr_ext_f1c
cdr scr_f2, alpcr_f2, form II
	real*8 :: scr_f2
	real*8 :: alpcr_f2
	real*8 :: scr_ext_f2
      end type cr_rate_coefficient

      type cooling_rate_coefficient
!  form. I condensed to form. II (for testing):
	real*8 :: ec1_f1c
        real*8 :: ec0_f1c
        real*8 :: ec_ext_f1c
!  form. I:
        real*8 :: ec1_f1
	real*8 :: ec2_f1
	real*8 :: ec3_f1
	real*8 :: ec0_f1

!  form. I condensed to form. II (for testing):
	real*8 :: rad1_f1c
        real*8 :: rad0_f1c
        real*8 :: rad_ext_f1c
!  form. I         
	real*8 :: rad1_f1
	real*8 :: rad2_f1
	real*8 :: rad3_f1
	real*8 :: rad0_f1
        real*8 :: rad_ext_f1
!  form. II:
	real*8 :: ec1_f2
	real*8 :: ec0_f2
	real*8 :: ec_ext_f2
	real*8 :: rad1_f2
	real*8 :: rad0_f2
        real*8 :: rad_ext_f2
      end type cooling_rate_coefficient

      type ral_col_liste
        integer :: isc, fsc 
        real*8 :: ftype, x0
        real*8 :: fitpara(0:5)
        type(ral_col_liste), pointer :: next
      end type ral_col_liste

      type ral_ion_liste
        integer :: isc 
        real*8 :: fitpara(0:5)
        type(ral_ion_liste), pointer :: next
      end type ral_ion_liste

      type(atomic_structure), public, save :: at_strc
      type(rate_coefficient), public, save :: rate
      type(population_coefficient), public, save :: popcoe
      type(cr_rate_coefficient), public, save :: crrate
      type(cooling_rate_coefficient), public, save :: coolrate

      type(ral_col_liste), pointer, public, save :: ral_col
      type(ral_ion_liste), pointer, public, save :: ral_ion

      real*8, public, save :: PI, BOHR, EMAS, PLANCK, CSPEED, BOLTZMANN,
     .                        ELECTRON_VOLT, RYDBERG, RYD, IPG, E_IPG

      real*8, public, save :: a_init(ULH,ULH)
      real*8, public, save :: cl_q(ULH,ULH), cl_b(ULH,ULH)
      real*8, public, save :: cl_beta(ULH,ULH), cl_gamma(ULH,ULH), 
     .                        cl_delta(ULH,ULH)
      real*8, public, save :: gt(ULH), ga(ULH), gb(ULH)

      contains

      subroutine eirene_init_constants
      
        PI = 4.d0 * ATAN(1.d0)

!+ * Updated Jan 2006 according to NIST
!-#define GSL_CONST_MKSA_BOHR_RADIUS (5.291772083e-11) /* m */
!+#define GSL_CONST_MKSA_BOHR_RADIUS (5.291772108e-11) /* m */
        BOHR = 5.291772083d-11  ! Bohr radius in m  ! from C

!-#define GSL_CONST_MKSA_MASS_ELECTRON (9.10938188e-31) /* kg */
!+#define GSL_CONST_MKSA_MASS_ELECTRON (9.1093826e-31) /* kg */
        EMAS = 9.10938188d-31  ! electron mass in kg ! from C

!-#define GSL_CONST_MKSA_PLANCKS_CONSTANT_H (6.62606876e-34) /* kg m^2 / s */
!+#define GSL_CONST_MKSA_PLANCKS_CONSTANT_H (6.6260693e-34) /* kg m^2 / s */
        PLANCK = 6.62606896d-34  ! Planck's H

! #define GSL_CONST_MKSA_SPEED_OF_LIGHT (2.99792458e8) /* m / s */
        CSPEED = 2.99792458d8  ! speed of light

!-#define GSL_CONST_MKSA_BOLTZMANN (1.3806503e-23) /* kg m^2 / K s^2 */
!+#define GSL_CONST_MKSA_BOLTZMANN (1.3806505e-23) /* kg m^2 / K s^2 */
        BOLTZMANN = 1.3806504d-23  ! Boltzmann constant

!-#define GSL_CONST_MKSA_ELECTRON_VOLT (1.602176462e-19) /* kg m^2 / s^2 */
!+#define GSL_CONST_MKSA_ELECTRON_VOLT (1.60217653e-19) /* kg m^2 / s^2 */
      electron_volt = 1.602176487d-19

      IPG = 198305.0    ! ionization potential of ground state in cm^{-1}
      E_IPG=CM2EV(IPG)  ! ionization potential of ground state in eV

!-#define GSL_CONST_MKSA_RYDBERG (2.17987190389e-18) /* kg m^2 / s^2 */
!+#define GSL_CONST_MKSA_RYDBERG (2.17987209e-18) /* kg m^2 / s^2 */
      RYDBERG = 2.17987196968d-18      ! taken from C program
      RYD = J2CM(RYDBERG)  ! Rydberg constant in cm^{-1}

      end subroutine eirene_init_constants

      function EV2J (x)
!     eV -> J 
      implicit none
      real*8, intent(in) :: x
      real*8 EV2J
      EV2J = electron_volt * x
      return
      end function EV2J

      function J2EV (x)
!     J -> eV
      implicit none
      real*8, intent(in) :: x
      real*8 J2EV
      J2EV = x / electron_volt 
      return
      end function J2EV

      function CM2J (x)
!     cm^{-1} -> J
      implicit none
      real*8, intent(in) :: x
      real*8 CM2J
      CM2J = PLANCK * CSPEED * (x * 1D2)
      return
      end function CM2J

      function J2CM (x)
!     J -> cm^{-1} 
      implicit none
      real*8, intent(in) :: x
      real*8 J2CM
      J2CM = x / (PLANCK * CSPEED) * 1D-2
      return
      end function J2CM

      function CM2EV (x)
!     cm^{-1} -> eV
      implicit none
      real*8, intent(in) :: x
      real*8 CM2EV
      CM2EV = J2EV(CM2J(x))
      return
      end function CM2EV

      function EV2CM (x)
!     eV -> cm^{-1}
      real*8, intent(in) :: x
      real*8 EV2CM
      EV2CM = J2CM(EV2J(x))
      return
      end function EV2CM

      function J2K (x)
!     J -> K
      real*8, intent(in) :: x
      real*8 J2K
      J2K = x / BOLTZMANN
      return
      end function J2K

      function EV2K (x)
!     eV -> K
      real*8, intent(in) :: x
      real*8 EV2K
      EV2K = J2K(EV2J(x))
      return
      end function EV2K

      function CM2K (x)
!     cm^{-1} -> K
      real*8, intent(in) :: x
      real*8 CM2K
      CM2K = J2K(CM2J(x))
      return
      end function CM2K

      subroutine eirene_init_rate_arrays
cdr  atomic structure
      at_strc%e  = 0.d0
      at_strc%en = 0.d0
      at_strc%g  = 0.d0
      at_strc%f  = 0.d0
      at_strc%a  = 0.d0
cdr single step rates, vs. Te
      rate%s = 0.d0
      rate%alpha = 0.d0
      rate%beta = 0.d0
      rate%betad = 0.d0
      rate%c = 0.d0
      rate%z  = 0.d0 ! Saha-Boltzmann coeff. vs. Te

      end subroutine eirene_init_rate_arrays

******************************************************************

      subroutine eirene_init_popcoe(iform)

      integer, intent(in) :: iform

      select case (iform)
      case(1)
        popcoe%r0 = 0.d0
        popcoe%r1 = 0.d0
        popcoe%r2 = 0.d0
        popcoe%r3 = 0.d0
        popcoe%r_ext = 0.d0
cdr derived from form 1: f1c (metastable condensed) populations
        popcoe%rr12 = 0.d0
        popcoe%rr13 = 0.d0
        popcoe%rr02 = 0.d0
        popcoe%rr03 = 0.d0
      case(2)
        popcoe%rr0 = 0.d0
        popcoe%rr1 = 0.d0
        popcoe%rr_ext = 0.d0
      end select

      end subroutine eirene_init_popcoe

******************************************************************

      subroutine eirene_read_ralchenko_collision_data

      implicit none
      integer :: iun, io, ncol, i
      character(500) :: zeile
      character(20) :: token(100)     
      type(ral_col_liste), pointer :: element, last

      nullify(ral_col)

      open (newunit=iun, file='HECRDATA/CRS_FIT_DATA_RALCHENKO_PAPER',
     .      status='OLD',iostat=io)
      if (io /= 0) then
        WRITE(*,*)  'failed to open: ',
     .              'HECRDATA/CRS_FIT_DATA_RALCHENKO_PAPER'
        CALL EIRENE_EXIT_OWN(1)
      END IF

      do 
        read (iun,'(A500)',iostat=io) zeile 
        if (io /= 0) exit
        if (zeile(1:1) == '#') cycle
        zeile = adjustl(zeile)

        allocate(element)
        element%fitpara = 0.d0
        nullify(element%next)

!        token=repeat(' ',20)
!        call eirene_split_line(zeile,token,ncol)
!        read (token(1),*) element%isc
!        read (token(2),*) element%fsc
!        read (token(3),*) element%ftype
!        read (token(4),*) element%x0
!        do i=5, ncol
!          read(token(i),*) element%fitpara(i-5) 
!        end do
        read (zeile,*) element%isc, element%fsc, element%ftype,
     .                 element%x0, element%fitpara(0:5) 

        if (.not.associated(ral_col)) then
          ral_col => element
          last => ral_col
        else
          last%next => element
          last => last%next
        end if
      end do

      close(iun)

      end subroutine eirene_read_ralchenko_collision_data

******************************************************************

      subroutine eirene_read_ralchenko_ionization_data

      implicit none
      integer :: iun, io, ncol, i
      character(500) :: zeile
      character(20) :: token(100)     
      type(ral_ion_liste), pointer :: element, last

      nullify(ral_ion)

      open (newunit=iun, file='HECRDATA/IONIZE_FIT_DATA_RALCHENKO',
     .      status='OLD',iostat=io)
      if (io /= 0) then
        WRITE(*,*)  'failed to open: ',
     .              'HECRDATA/IONIZ_FIT_DATA_RALCHENKO'
        CALL EIRENE_EXIT_OWN(1)
      END IF

      do 
        read (iun,'(A500)',iostat=io) zeile 
        if (io /= 0) exit
        if (zeile(1:1) == '#') cycle
        zeile = adjustl(zeile)

        allocate(element)
        element%fitpara = 0.d0
        nullify(element%next)

!        token=repeat(' ',20)
!        call eirene_split_line(zeile,token,ncol)
!        read (token(1),*) element%isc
!        do i=2, ncol
!          read(token(i),*) element%fitpara(i-2) 
!        end do
        read (zeile,*) element%isc, element%fitpara(0:5)

        if (.not.associated(ral_ion)) then
          ral_ion => element
          last => ral_ion
        else
          last%next => element
          last => last%next
        end if
      end do

      close (iun)

      end subroutine eirene_read_ralchenko_ionization_data

******************************************************************

      subroutine eirene_split_line (line, token, no)
cdr  comments missing
      character(*) :: line
      character(20) :: token(*)
      integer, intent(out) :: no
      integer :: ll, ic

      no = 0
      ll = len_trim(line)
!      token=repeat(' ',20)

      do while(ll > 0)
        line=adjustl(line)
        ic = scan(line,' ')
        no = no + 1
        token(no)(1:ic-1)= line(1:ic-1)
        line=line(ic+1:ll)
        ll=len_trim(line)
      end do

      end subroutine eirene_split_line

      end module eirmod_hecr
