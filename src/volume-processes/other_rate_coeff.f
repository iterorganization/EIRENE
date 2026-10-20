cdr  added: cell number ic, for optimization in called CR routines colrad.f
cdr  call driver routine colrad for CR models. new variable: iflavor
!pb  01.06.2017  copied from rate_coeff.f

c  to be done: h_colrad called twice per cell ??
c              re-use erate from previous call to rate_coeff

cdr  19.02.14: COMMENTS
cdr sept.15:  lexp not fully written, in case of adas 2d tables
cdr           ifit=4 option was missing (1D tables). added, but not checked.





      function EIRENE_other_rate_coeff (ir, ic, p1, p2, lexp, ip2shft)
     .                     result (orate)

!  evaluate other atomic data: mostly: population coefficients, density ratios
!                 e.g. from AMJUEL H.11, H.12 sections
!  and return this as "orate"

!  currently 5 different options controlled by 'reacdat(ir)%rtc%ifit'
!  Only ifit=2 and ifit=3 tested so far. Caution!
!  ifit=1: single polynomial fit, use P1, (e.g. HYDHEL, H.2)
!  ifit=2: double polynomial fit, use P1, P2, (e.g. HYDHEL, H.3, AMJUEL, H.4, H.12,...)
!  ifit=3: interpolation in 2-parameter table (e.g. ADAS)
!  ifit=4: interpolation in single parameter table (e.g. open ADAS,...)
!  ifit=5: use internal eirene collision radiative code. To be generalized
!          (currently here also other rates, orate, for this particular option).
!          More logical if the latter are moved
!          to routine "eirene_energy_rate_coeff"

!   input:
!   ir:    reaction number, as stored in eirene arrays.
!   ic:    cell number
!   p1:    first parameter (usually: log_e temperature,...)
!   p2:    second parameter  (if any, e.g. log_e (density),...,log_e(test particle energy),...)
!   lexp:  return orate=rate coefficient in ... units
!   not lexp: return orate=log_e(rate coefficient) with rate coefficient in ...units
!   ip2shft: >0: carry out shift in parameter p2 for fit expression evaluation,
!             currently hard-wired: factor 1e-8.  p2 --> p2*factor
!             Currently : only for ifit=2, polynomial fits vs. ne, T, ne in units 1e8 *cm**-3.
!             emissivity.f relies on the current use of ip2shft in the tested cases!

! to be done: lexp option for ifit=4
!             ip2shft option: currently hard-wired only for ifit=2 and shift = 1e-8
!             what happens if later call with other shift ?  coding to be reconsidered !

!             remove ifirst and ifsub conditions and set the data once, and save. DONE (Nov. 15)

      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_ccona
      use EIRMOD_ctrcei, only: trcamd
      use EIRMOD_comprt, only: iunout
      use EIRMOD_COLRAD, ONLY: EIRENE_COLRAD

      implicit none

      integer, intent(in) :: ir, ip2shft, ic
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp

      real(dp) :: res, orate, EIRENE_sngl_poly, dum(9),
     .            pp1, rc1min,  rc1max, fp1(6),
     .            pp2, rc2min,  rc2max, fp2(6),
     .                 earrh0,
     .                 rrc2min, rrc2max,
     .            O_SCR
      real(dp), save :: xlog10e =  4.34294482d-01, !1./ln(10) = log10(e)
     .                  xln10   =  2.30258509299_dp, !ln(10)
c  transformation of parameters p1 and p2:
     .                  dsub    = 18.420680744_dp  !ln(1e8), hard-wired.
                                        !But should come from database

      integer :: jfex1mn, jfex1mx,jfex2mn, jfex2mx
      integer :: ip1, ip2, icrm, ivar, iform
      external :: eirene_sngl_poly, eirene_dbl_poly, eirene_exit_own

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

      if (.not.reacdat(ir)%loth) then
        write (iunout,*) ' no data for other reaction available',
     .                   ' for reaction ',ir
        call EIRENE_exit_own(1)
      end if

      orate = 0._dp

c.............................................................


      if (mod(iftflg(ir,2),100) == 10) then

!  SET A CONSTANT RATE
        orate = reacdat(ir)%oth%poly%dblpol(1,1)

cdr   missing: iftflg < 100:  multiply density,  else: not

cdr   lexp missing

c.............................................................

      elseif (reacdat(ir)%oth%ifit == 1) then

!  SINGLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE), FOR LN(OTHER RATE)

c  extrapolation data: for 1d polynomial fits
        rc1min  = reacdat(ir)%oth%rc1min
        rc1max  = reacdat(ir)%oth%rc1max
        fp1(1:3)= reacdat(ir)%oth%fp1l
        fp1(4:6)= reacdat(ir)%oth%fp1r
        jfex1mn = reacdat(ir)%oth%jfex1mn
        jfex1mx = reacdat(ir)%oth%jfex1mx
        earrh0  = reacdat(ir)%earrh0
cdr  careful: Arrhenius factor for orate?
        earrh0=0._DP

        orate = eirene_sngl_poly(reacdat(ir)%oth%poly%dblpol(1:9,1),
     .                   p1, rc1min, rc1max, fp1, jfex1mn, jfex1mx,
     .                   earrh0, trcamd, lexp)

c..............................................................


      else if (reacdat(ir)%oth%ifit == 2) then

!  DOUBLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE) AND P2, FOR LN(OTHER RATE)

c  extrapolation data: for 2d polynomial fits
        rc1min  = reacdat(ir)%oth%rc1min
        rc1max  = reacdat(ir)%oth%rc1max
        rc2min  = reacdat(ir)%oth%rc2min ! IF H.4 HERE 1E8,
                                  ! ALREADY SET IN CALLING PROGRAM
        rc2max  = reacdat(ir)%oth%rc2max ! if H.4 HERE 1E16,
                                         ! NOT YET SET.
        fp1(1:3)= reacdat(ir)%oth%fp1l
        fp1(4:6)= reacdat(ir)%oth%fp1r
        fp2(1:3)= reacdat(ir)%oth%fp2b
        fp2(4:6)= reacdat(ir)%oth%fp2t
        jfex1mn = reacdat(ir)%oth%jfex1mn
        jfex1mx = reacdat(ir)%oth%jfex1mx
        jfex2mn = reacdat(ir)%oth%jfex2mn
        jfex2mx = reacdat(ir)%oth%jfex2mx


c  rescale parameter p2 (currently only by 1e-8 for density): pp2
        pp2 = p2
        rrc2min=rc2min
        rrc2max=rc2max
        if (ip2shft > 0) then
cdr hidden link?  hard coded MODC=3 option: density dependence
          pp2 = pp2 - dsub
          rrc2min=rc2min - dsub
          rrc2max=rc2max - dsub
        endif
cdr     write (iunout,*) 'other rate '

        call EIRENE_dbl_poly
     .       (reacdat(ir)%oth%poly%dblpol,p1,pp2,orate,dum,
     .        rc1min,  rc1max,  fp1, jfex1mn, jfex1mx,
     .        rrc2min, rrc2max, fp2, jfex2mn, jfex2mx,
     .        trcamd)

C       if (.not. lexp)  orate=orate
        if (lexp)        orate = exp(max(-100._dp,orate))

c..............................................................

      else if (reacdat(ir)%oth%ifit == 3) then

! 2D TABULAR INPUT, FOR LOG10 OF RATE,  cm^3/s
! E.G.: ADAS FILES
cdr  extrapolation data: for 2d tabulated data, option not ready
cdr  to be added here

!  currently hard-wired: input parameters pp1, pp2 and table coefficients are log10

c  convert parameters p1 and p2 from ln to log10: pp1,pp2
        pp1 = xlog10e*p1
        pp2 = xlog10e*p2
C  assume here: tabulated data are log10 (to be generalized)
        orate = eirene_intp_tab2d(reacdat(ir)%oth%adas,pp1,pp2,ip1,ip2)

        if (lexp) then
          orate=10._dp**orate
        else
          orate = xln10*orate ! convert from log10(orate) to ln(orate)

        end if
cdr this unit conversion must be wrong in case lexp !!


c..............................................................

      else if (reacdat(ir)%oth%ifit == 4) then

! SINGLE PARAMETER TABLE
cdr  extrapolation data: for 1d tabulated data: option not ready (only CxHy data ?)
cdr  to be added here

! currently hard-wired: input parameters q1 and table coefficients are neither ln nor log10

        pp1 = exp(p1)
C  assume here: tabulated data are neither ln nor log10 (to be generalized)
        orate = eirene_intp_tab1d(reacdat(ir)%oth%tab1d,pp1,ip1)

!  lexp option not connected here !

c..............................................................

      else if (reacdat(ir)%oth%ifit == 5) then

! INTERNAL COLLISION RADIATIVE CODE
c  returns rates, rate coefs, etc, not LN or LOG of rates

c  convert (log) parameters p1, p2 to exp(p1), exp(p2): PP1,PP2
        PP1 = EXP(P1)
        PP2 = EXP(P2)

        icrm = reacdat(ir)%oth%crm%iflav
        ivar = reacdat(ir)%oth%crm%ivarst
        iform = reacdat(ir)%oth%crm%iformul

        CALL EIRENE_COLRAD(IR, ICRM, RES,
     .                     IFORM, IVAR, IC, PP1, PP2)

!  lexp option was not connected here, but used in xstei.f ! corrected, Oct. 28th 2015

        if (lexp) then
          orate = res
c  not exp: return log of rate
        elseif (res.gt.0.0) then 
          orate = log(res)
        elseif (res.eq.0.0) then
          orate =-50.
        else
          write (iunout,*) 'orate_coef, KK ',IR
          write (iunout,*) 'wrong sign from cr model'
          write (iunout,*) 'p1,p2,orate ',pp1,pp2,res
          write (iunout,*) 'return orate=exp(-50)'
          orate =-50.
        endif

      end if


      return

      end function EIRENE_other_rate_coeff
