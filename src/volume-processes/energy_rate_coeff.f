!pb  22.11.06: flag ip2shft for shift of second parameter to rate_coeff introduced
!pb  24.11.06: get extrapolation parameters for polynomial fit only
!pb  30.11.06: divide energy-weighted rate coefficient by ELCHA to get correct units

c  to be done: h_colrad called twice per cell ??
c              re-use erate from previous call to rate_coeff

cdr  19.02.14: COMMENTS
cdr sept.15:  lexp not fully written, in case of adas 2d tables
cdr           also unit conversion incorrect in that case. --> tbd
cdr           ifit=4 option was missing (1D tables). added, but not checked.

cdr  16.11.15: bug fix: error in arguments in call to H_colrad
cdr            added: ifit=4 option (1d table interpolation)
cdr            additional parameters in calls to 1d and 2d table interpolation
cdr  26.11.15: additional parameter IC in call to H_colrad,
cdr            for later use to identify "visited cells"
cdr  sept. 16: started to add extrapolation options. not ready....
cdr  jan.  18: call driver routine for CR models: colrad.f
cdr            tbd: electron energy loss rates can change sign.
cdr            Be careful with log(e_src). Routine should only be called with
cdr            LEXP=.true.

      function EIRENE_energy_rate_coeff (ir, ic, p1, p2, lexp, ip2shft)
     .                            result (erate)

!  evaluate energy-weighted rate coefficient, eV/s per incident particle,
!  and return this as "erate"

! erate is a rate coefficient weighted with an energy, e.g. an energy cost,
!       or energy gain.
! It must be positive, because also ln(erate) or log(erate) is used,
!    in certain data formats.
! If it is a loss, rather than a gain, sign change to be done in calling routine,
!   as well as a shift (if any) by potential energy loss rate
!   (e.g. conversion from electron cooling rate to radiation loss rate)

!  currently 5 different options controlled by 'reacdat(ir)%rtcew%ifit'
!  ifit=1: single polynomial fit, use P1, (e.g. HYDHEL, AMJUEL, H.8)
!  ifit=2: double polynomial fit, use P1, P2, (e.g. HYDHEL, H.9, AMJUEL, H.10,...)
!  ifit=3: interpolation in 2 parameter table (e.g. ADAS)
!  ifit=4: interpolation in single parameter table (e.g. open ADAS, ...)
!  ifit=5: use internal eirene collision radiative code. To be generalized

!   input:
!   ir:    reaction number, as stored in eirene arrays.
!   ic:    cell number (e.g. for internal CR models).
!   p1:    first parameter (usually: log_e temperature,...)
!   p2:    second parameter (if any, e.g. log_e (density),...,log_e(test particle energy),...)
!   lexp:  return erate=energy-weighted rate coefficient in eV*cm**3/sec
!   not lexp: return erate=log_e(erate coefficient) with rate coefficient in cm**3/sec
!   ip2shft:  >0: carry out shift in parameter p2 for fit expression evaluation,
!                 currently hard-wired: 1e-8.
!                (currently : only for ifit=2, polynomial fits vs. ne, T, ne in units 1e8 *cm**-3)

! to be done:
!
!             ip2shft option: currently hard-wired only for ifit=2 and shift = 1e-8
!             what happens if later call with other shift ?  coding to be reconsidered !

!             remove ifirst and ifsub conditions and set the data once, and save.  DONE (Nov. 15)

      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_ccona
      use EIRMOD_ctrcei, only: trcamd
      use EIRMOD_comprt, only: iunout
      use EIRMOD_COLRAD, ONLY: EIRENE_COLRAD

      implicit none

      integer, intent(in) :: ir, ic, ip2shft
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp

      real(dp) :: res, erate, EIRENE_sngl_poly, dum(9),
     .            pp1, rc1min,  rc1max, fp1(6),
     .            pp2, rc2min,  rc2max, fp2(6),
     .                 earrh0,
     .                 rrc2min, rrc2max
      real(dp), save :: xlog10e =  4.34294482d-01, !1./ln(10) = log10(e)
     .                  xln10   =  2.30258509299_dp, !ln(10)
     .                  dsub    = 18.420680744_dp, !ln(1e8), hard-wired.
                                    ! But should come from database
     .                  xlnelch =-43.2777390821    !ln(elcha)
      integer :: jfex1mn, jfex1mx, jfex2mn, jfex2mx
      integer :: ip1, ip2, icrm, ivar, iform
      external :: EIRENE_sngl_poly

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

      EXTERNAL :: EIRENE_DBL_POLY, EIRENE_EXIT_OWN

      if (.not.reacdat(ir)%lrtcew) then
        write (iunout,*) ' no data for energy-weighted rate',
     .                   ' coefficient available for reaction ',ir
        call EIRENE_exit_own(1)
      end if

      erate = 0._dp

c.............................................................


      if (mod(iftflg(ir,4),100) == 10) then

!  SET A CONSTANT ENERGY-WEIGHTED RATE COEFFICIENT.
!  re-use format "poly":  LN of energy-weighted rate.
        res = reacdat(ir)%rtcew%poly%dblpol(1,1)

        if (lexp) then
          erate = exp(max(-100._dp,res))
        else
          erate=res
        endif

cdr missing: iftflg < 100: multiply density, else: not


c.............................................................

      elseif (reacdat(ir)%rtcew%ifit == 1) then

!  SINGLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE), FOR LN(ENERGY-WEIGHTED RATE)

c  extrapolation data: for 1d polynomial fits
        rc1min  = reacdat(ir)%rtcew%rc1min
        rc1max  = reacdat(ir)%rtcew%rc1max
        fp1(1:3)= reacdat(ir)%rtcew%fp1l
        fp1(4:6)= reacdat(ir)%rtcew%fp1r
        jfex1mn = reacdat(ir)%rtcew%jfex1mn
        jfex1mx = reacdat(ir)%rtcew%jfex1mx
cdr  careful:  arrhenius factor for energy rate?
        earrh0  = reacdat(ir)%earrh0
        earrh0  = 0._DP

        erate = eirene_sngl_poly(reacdat(ir)%rtcew%poly%dblpol(1:9,1),
     .                           p1,rc1min,rc1max,fp1,jfex1mn,jfex1mx,
     .                           earrh0,trcamd,lexp)

! lexp=false: erate is ln(energy rate), with energy rate >0.
! lexp=true : erate is the energy rate
! If it is loss, rather than a gain, sign change to be done in calling routine,
! as well as shift (if any) by potential energy loss rate


c..............................................................


      else if (reacdat(ir)%rtcew%ifit == 2) then

!  DOUBLE POLYNOMIAL FIT VS. P1 =LN(TEMPERATURE) AND P2,  FOR LN(ENERGY-WEIGHTED RATE)

c  extrapolation data: for 2d polynomial fits
        rc1min  = reacdat(ir)%rtcew%rc1min
        rc1max  = reacdat(ir)%rtcew%rc1max
        rc2min  = reacdat(ir)%rtcew%rc2min
        rc2max  = reacdat(ir)%rtcew%rc2max
        fp1(1:3)= reacdat(ir)%rtcew%fp1l
        fp1(4:6)= reacdat(ir)%rtcew%fp1r
        fp2(1:3)= reacdat(ir)%rtcew%fp2b
        fp2(4:6)= reacdat(ir)%rtcew%fp2t
        jfex1mn = reacdat(ir)%rtcew%jfex1mn
        jfex1mx = reacdat(ir)%rtcew%jfex1mx
        jfex2mn = reacdat(ir)%rtcew%jfex2mn
        jfex2mx = reacdat(ir)%rtcew%jfex2mx


c  rescale parameter p2 (currently only by 1e-8 for density): pp2
        pp2 = p2
        rrc2min=rc2min
        rrc2max=rc2max
        if (ip2shft > 0) then
cdr  hidden link?  hard coded MODC=3 option (density dependence)
          pp2 = pp2 - dsub
          rrc2min=rc2min - dsub
          rrc2max=rc2max - dsub
        endif

        call EIRENE_dbl_poly
     .       (reacdat(ir)%rtcew%poly%dblpol,p1,pp2,res,dum,
     .        rc1min,  rc1max,  fp1, jfex1mn, jfex1mx,
     .        rrc2min, rrc2max, fp2, jfex2mn, jfex2mx,
     .        trcamd)

! RES is ln(energy rate), with energy rate >0.
! If it is loss, rather than a gain, sign change to be done in calling routine,
! as well as shift (if any) by potential energy loss rate

        if (lexp) then
          erate = exp(max(-100._dp,res))
        else
          erate=res
        endif

c..............................................................

      else if (reacdat(ir)%rtcew%ifit == 3) then

! 2D TABULAR INPUT, FOR LOG10 OF ENERGY-WEIGHTED RATE, Joule*cm^3/s
! E.G.: ADAS adf11 PLT and PRB FILES
cdr  extrapolation data: for 2d tabulated data, option not ready
cdr  to be added here

!  currently hard-wired: input parameters pp1, pp2 and table coefficients are log10

c  convert parameters p1 and p2 from ln to log10: pp1,pp2
        pp1 = xlog10e*p1
        pp2 = xlog10e*p2
C  assume here: tabulated data are log10 (to be generalized)
c  and in joule cm^3/s
        res=eirene_intp_tab2d(reacdat(ir)%rtcew%adas,pp1,pp2,ip1,ip2)

        if (lexp) then
          erate=10._dp**res/elcha
        else
c  in this database model erate is strictly positive, and log10(erate) is returned
          erate = xln10*res       !  convert from log10(erate)
                                  !  to ln(erate)
          erate = erate - xlnelch !  convert ln(erate) from
                                  !  ln[joule cm^3/s] to ln[eV cm^3/s]
        end if

c  deal with bremsstrahlung. currently we assume that bremsstrahlung (free-free)
c  electron energy cost is included in the rates, if they come from ADAS.
c  (PRB coefficients contain bremsstrahlung, but PLT do not.)

c  tbd: So subtract this part in read_tab2d already, not here nor in calling routine.




c..............................................................

      else if (reacdat(ir)%rtcew%ifit == 4) then

!  proprietary option: not ready
        if (.true.) goto 990
! SINGLE PARAMETER 1D TABLE
cdr  extrapolation data: for 1D tabulated data: option not ready (only CxHy data ?)
cdr  to be added here

! currently hard-wired: input parameters q1 and table coefficients are neither ln nor log10

        pp1 = exp(p1)
C  assume here: tabulated data are neither ln nor log10 (to be generalized)
        res = eirene_intp_tab1d(reacdat(ir)%rtcew%tab1d,pp1,ip1)

!  lexp option not connected here !

c..............................................................

      else if (reacdat(ir)%rtcew%ifit == 5) then

! INTERNAL COLLISION RADIATIVE CODE
c  returns rates, rate coefs, etc, not LN or LOG of rates

c  convert (log) parameters p1, p2 to exp(p1), exp(p2): PP1,PP2
        PP1 = EXP(P1)
        PP2 = EXP(P2)

        icrm = reacdat(ir)%rtcew%crm%iflav
        ivar = reacdat(ir)%rtcew%crm%ivarst
        iform = reacdat(ir)%rtcew%crm%iformul

        CALL EIRENE_COLRAD(IR, ICRM, RES,
     .                     IFORM, IVAR, IC, PP1, PP2)

!  electron energy-weighted loss rates are taken positive in CRM COLRAD, and negative if
!  it is a gain. For negative (i.e. gain) rates, the log(e-rate) return is not possible.
!  energy rate coefficient should always only be called with LEXP=.TRUE. for such processes

        IF (LEXP) then
          erate = res
c  not exp: return log of rate
        elseif (res.gt.0.0) then
          erate = log(res)
        elseif (res.eq.0.0) then
          erate =-50.
        else
          write (iunout,*) 'erate_coef, KK ',IR
          write (iunout,*) 'wrong sign from cr model'
          write (iunout,*) 'p1,p2,erate ',pp1,pp2,res
          write (iunout,*) 'return exp(-50)'
          erate =-50.
        endif

      end if


      return

  990 continue
      write (iunout,*) 'Proprietary (unfinished) option ifit=4 '
      write (iunout,*) 'encountered in routine: energy_rate_coeff.f '
      call eirene_exit_own(1)

      end function EIRENE_energy_rate_coeff
