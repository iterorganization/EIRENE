!pb  22.11.06: flag ip2shft for shift of second parameter to rate_coeff introduced
!pb  24.11.06: get extrapolation parameters for polynomial fit only

!pb  07.12.06: double declaration of dsub removed
!dr  19.02.14: COMMENTS
!dr  29.10.15: bug fix: lexp option for ifit=5 was missing.
!dr            no consequences for any earlier runs, except CRM option with H.2 (corona) rates.
cdr  nov. 15:  indicators ip1, ip2 for extrapolation or interpolation added,
cdr            in intp_tab1d and intp_tab2d
cdr            rename q1,q2 to pp1,pp2: modified input parameters p1, p2.
cdr sept. 22:  redundancy and ambiguity in flags: reacdat and modclf.
cdr            Consistency tests added.
cdr            Added: MODC, and distinguish: MODC=2 and MODC=3

      function EIRENE_rate_coeff (ir, ic, p1, p2, lexp, ip2shft)
     .                     result (rate)

!  evaluate reaction rate coefficient (cm^3/s), or log_e thereof,
!  and return this as "rate"

!  currently 5 different options controlled by 'reacdat(ir)%rtc%ifit'
!  ifit=1:   single polynomial fit, use P1, (e.g. HYDHEL, AMJUEL, H.2)
!  ifit=2:   double polynomial fit, use P1, P2, (e.g. HYDHEL, H.3, AMJUEL, H.4,...)
!  ifit=3:   interpolation in 2 parameter table (e.g. ADAS)
!  ifit=4:   interpolation in single parameter table (e.g. open ADAS, ...)
!  ifit=5:   use internal eirene collision radiative code. To be generalized

!   input:
!   ir:        Reaction number, as stored in eirene input arrays, e.g. reacdat or modclf.
!              Negative values of IR (-1 to -11): minimal internal eirene A&M models
!              are only available for REACDAT and IFTFLG
!   ic:        cell number
!   p1:        first parameter (usually: log_e temperature,...)
!   p2:        second parameter (if any, e.g. log_e (density),...,log_e(test particle energy),...)
!   lexp:      return rate=rate coefficient in cm**3/sec
!   not lexp:  return rate=log_e(rate coefficient) with rate coefficient in cm**3/sec
!   ip2shft:   >0: carry out shift in parameter p2 for fit expression evaluation,
!              currently hard-wired: 1e-8.
!             (currently : only for ifit=2, polynomial fits vs. ne, T, ne in units 1e8 *cm**-3)

! to be done:
!              ip2shft option: currently hard-wired only for ifit=2 and shift = 1e-8
!              what happens if later call with other shift ?  coding to be reconsidered !

!              remove ifirst and ifsub conditions and set the data once, and save.

      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_ctrcei, only: trcamd
      use EIRMOD_comprt, only: iunout
      use EIRMOD_COLRAD, ONLY: EIRENE_COLRAD

      implicit none

      integer, intent(in) :: ir, ip2shft, ic
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp
      real(dp) :: res, rate, EIRENE_sngl_poly, dum(9),
     .            pp1, rc1min,  rc1max, fp1(6),
     .            pp2, rc2min,  rc2max, fp2(6),
     .                 earrh0,
     .                 rrc2min, rrc2max
c  factors to convert from ln to log10 and back
      real(dp), save :: xlog10e = 4.34294482d-01, ! 1./ln(10) = log10(e)
     .                  xln10   = 2.30258509299_dp, ! ln(10)  = loge(10)
c  transformation of parameters p1 and p2:
     .                  dsub    = 18.420680744_dp ! ln(1e8), hard-wired.
                                ! But should come from database

      integer :: jfex1mn, jfex1mx,jfex2mn, jfex2mx
      integer :: ip1, ip2, icrm, ivar, iform
      integer :: modc
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      EXTERNAL :: EIRENE_SNGL_POLY, EIRENE_DBL_POLY, EIRENE_EXIT_OWN

      interface
        function EIRENE_intp_tab2d (ad,p1,p2,ip1,ip2) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: adas_data
          type(adas_data), pointer :: ad
          real(dp), intent(in) :: p1, p2
          integer, intent(out) :: ip1,ip2
          real(dp) :: res
        end function EIRENE_intp_tab2d

        function EIRENE_intp_tab1d (tb,p1,ip1) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: tab1d_data
          type(tab1d_data), pointer :: tb
          real(dp), intent(in) :: p1
          integer, intent(out) :: ip1
          real(dp) :: res
        end function EIRENE_intp_tab1d
      end interface

      if (.not.reacdat(ir)%lrtc) then
        write (iunout,*) ' no data for rate',
     .                   ' coefficient available for reaction ',ir
        call EIRENE_exit_own(1)
      end if

      rate = 0._dp

      modc = 1
!pb modclf is only defined for non-default reactions
      if (ir > 0) then
        modc=eirene_idez(modclf(ir),3,5)
      end if

c.............................................................

      if (mod(iftflg(ir,2),100) == 10) then          ! iftflg(-11:nreac)

!  SET A CONSTANT RATE
cdr potential hidden link: this iftflg option only works with IFIT=1 ?
        rate = reacdat(ir)%rtc%poly%dblpol(1,1)

cdr   missing: iftflg < 100: multiply density,  else: not

cdr   lexp missing

c.............................................................

      elseif (reacdat(ir)%rtc%ifit == 1) then

        if (modc.ne.1) then
          write (iunout,*) 'inconsistency in rate_coeff.f '
          write (iunout,*) 'modc, reacdat ',modc, reacdat(ir)%rtc%ifit
        endif

!  SINGLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE), FOR LN(RATE)

c  extrapolation data: for 1d polynomial fits
        rc1min  = reacdat(ir)%rtc%rc1min
        rc1max  = reacdat(ir)%rtc%rc1max
        fp1(1:3)= reacdat(ir)%rtc%fp1l
        fp1(4:6)= reacdat(ir)%rtc%fp1r
        jfex1mn = reacdat(ir)%rtc%jfex1mn
        jfex1mx = reacdat(ir)%rtc%jfex1mx
        earrh0  = reacdat(ir)%earrh0

        rate = eirene_sngl_poly(reacdat(ir)%rtc%poly%dblpol(1:9,1),
     .                   p1, rc1min, rc1max, fp1, jfex1mn, jfex1mx,
     .                   earrh0, trcamd, lexp)

C lexp condition: already done in sngl_poly


c..............................................................


      else if (reacdat(ir)%rtc%ifit == 2) then

!  DOUBLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE) AND P2=LN(...), FOR LN(RATE)

c  extrapolation data: for 2d polynomial fits
        rc1min  = reacdat(ir)%rtc%rc1min
        rc1max  = reacdat(ir)%rtc%rc1max
        rc2min  = reacdat(ir)%rtc%rc2min ! IF H.4 HERE 1E8,
                                 ! ALREADY SET IN CALLING PROGRAM
        rc2max  = reacdat(ir)%rtc%rc2max ! if H.4 HERE 1E16,
                                         ! NOT YET SET.
        fp1(1:3)= reacdat(ir)%rtc%fp1l
        fp1(4:6)= reacdat(ir)%rtc%fp1r
        fp2(1:3)= reacdat(ir)%rtc%fp2b
        fp2(4:6)= reacdat(ir)%rtc%fp2t
        jfex1mn = reacdat(ir)%rtc%jfex1mn
        jfex1mx = reacdat(ir)%rtc%jfex1mx
        jfex2mn = reacdat(ir)%rtc%jfex2mn
        jfex2mx = reacdat(ir)%rtc%jfex2mx


c  ip2shft: rescale parameter p2 (currently only by factor 1e-8 for electron density): pp2
c  In this case: p2 = ln(ne), and density ne in cm**-3, pp2= ln(ne/1e8) as in amjuel fits.

        pp2 = p2
        rrc2min=rc2min
        rrc2max=rc2max
        if (ip2shft > 0) then
          pp2 = pp2 - dsub
          rrc2min=rc2min - dsub
          rrc2max=rc2max - dsub
          if (modc.ne.3) then
cdr density parameter shift 1e-8
            write (iunout,*) 'inconsistency in rate_coeff.f '
            write (iunout,*) 'modc, reacdat, ip2shft ',
     .                        modc, reacdat(ir)%rtc%ifit,ip2shft
          endif
        else
          if (modc.ne.2) then
cdr energy parameter, no shift
            write (iunout,*) 'inconsistency in rate_coeff.f '
            write (iunout,*) 'modc, reacdat, ip2shft ',
     .                        modc, reacdat(ir)%rtc%ifit,ip2shft
          endif
        endif
cdr     write (iunout,*) 'particle rate '

        call EIRENE_dbl_poly
     .       (reacdat(ir)%rtc%poly%dblpol,p1,pp2,rate,dum,
     .        rc1min,  rc1max,  fp1, jfex1mn, jfex1mx,
     .        rrc2min, rrc2max, fp2, jfex2mn, jfex2mx,
     .        trcamd)

C       if (.not. lexp)  rate=rate
        if (lexp)        rate = exp(max(-100._dp,rate))

c..............................................................

      else if (reacdat(ir)%rtc%ifit == 3) then

! 2D TABULAR INPUT, FOR LOG10 OF RATE, cm^3/s
! E.G.: ADAS adf11 ACD and SCD FILES
cdr  extrapolation data: for 2d tabulated data, option not ready
cdr  to be added here

!  currently hard-wired: TAB2D input parameters pp1, pp2 and table coefficients are log10

c  convert parameters p1,p2 from ln to log10: pp1,pp2
        pp1 = xlog10e*p1
        pp2 = xlog10e*p2
C  assume here: tabulated data are log10 (to be generalized)
        res = eirene_intp_tab2d(reacdat(ir)%rtc%adas,pp1,pp2,ip1,ip2)
c  Now: res is log10 of a rate

        if (lexp) then
          rate=10._dp**res
        else
          rate = xln10*res ! convert from log10(rate) to ln(rate)
        end if
cdr this unit conversion must be wrong in case lexp !!


c..............................................................

      else if (reacdat(ir)%rtc%ifit == 4) then

! SINGLE PARAMETER TABLE
cdr  extrapolation data: for 1d tabulated data: option not ready (only CxHy data ?)
cdr  to be added here

! currently hard-wired: input parameters q1 and table coefficients are neither ln nor log10

        pp1 = exp(p1)
C  assume here: tabulated data are neither ln nor log10 (to be generalized)
        rate = eirene_intp_tab1d(reacdat(ir)%rtc%tab1d,pp1,ip1)

!  lexp option not connected here !

c..............................................................

      else if (reacdat(ir)%rtc%ifit == 5) then

! INTERNAL COLLISION RADIATIVE CODE
c  returns rates, rate coefs, etc, not LN or LOG of rates

c  convert (log) parameters p1, p2 to exp(p1), exp(p2): PP1,PP2
        PP1 = EXP(P1)
        PP2 = EXP(P2)

        icrm = reacdat(ir)%rtc%crm%iflav
        ivar = reacdat(ir)%rtc%crm%ivarst
        iform = reacdat(ir)%rtc%crm%iformul

        CALL EIRENE_COLRAD(IR, ICRM, RES,
     .                     IFORM, IVAR, IC, PP1, PP2)

! lexp option was not connected here, but used in xstei.f
! corrected, Oct. 28th 2015

        if (lexp) then
          rate = res
        elseif (res.gt.0.0) then
          rate = log(res)
        elseif (res.eq.0.0) then
          rate =-50.
        else
          write (iunout,*) 'prate_coef, KK ',IR
          write (iunout,*) 'wrong sign from cr model'
          write (iunout,*) 'p1,p2,rate ',pp1,pp2,res
          write (iunout,*) 'return exp(-50)'
          rate =-50.
        endif

      end if

      return

      end function EIRENE_rate_coeff
