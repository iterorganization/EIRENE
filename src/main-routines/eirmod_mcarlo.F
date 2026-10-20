      MODULE EIRMOD_MCARLO

cdr  Called from:...

cdr  Perform a (linear) Monte Carlo simulation, loop over strata,
cdr  for one single iteration or time step.
cdr  Within each stratum: loop over Monte Carlo histories.

cdr  Scale and save tallies, including additional post processed tallies.
cdr  Non-linearities may be accounted for in calling program,
cdr  via iterative loops or stepping.

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CAI, ONLY: XMCT
      USE EIRMOD_COMUSR
      USE EIRMOD_CTSURF
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CLOGAU
      USE EIRMOD_CRAND
      USE EIRMOD_CGRID
      USE EIRMOD_CSPEZ
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_COMSPL
      USE EIRMOD_COMSIG
      USE EIRMOD_CLGIN
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEI
      USE EIRMOD_CUPD
      USE EIRMOD_PHOTON
      USE EIRMOD_MPI
      USE EIRMOD_SECOND_OWN, ONLY: EIRENE_SECOND_OWN
      USE EIRMOD_SAMVOL, ONLY: EIRENE_SAMVL0
      USE EIRMOD_RANF, ONLY: RANF_EIRENE, RANSET_EIRENE
      USE EIRMOD_SWITCH_PARTINFO, ONLY: EIRENE_SWITCH_PARTINFO,
     .                                  EIRENE_OUTPUT_PARTINFO
      USE EIRMOD_LOCATE, ONLY: EIRENE_LOCAT0, EIRENE_LOCAT1
      USE EIRMOD_SAMSRF, ONLY: EIRENE_SAMSF0
      USE EIRMOD_SPUTER, ONLY: EIRENE_SPUTR0
#ifdef USE_OPENMP
#ifndef USE_EXT_OPENMP
cym these variables need copyin
      USE EIRMOD_COMXS, ONLY: XSTOR, XSTORV
      USE EIRMOD_SPUTER, ONLY: ETH, Q, M2M1, ES, ETF
      USE EIRMOD_SAMVOL, ONLY: FREC, VSOURC, VSMXI, RQ21, PS21, ASIMP,
     .                         ISOURC, ICMX, IFREC, LPLSSR, LATSSR
#endif
#endif
      USE EIRMOD_STATIS, ONLY: EIRENE_STATS0, EIRENE_STATS1,
     .                         EIRENE_STATS2
      USE EIRMOD_UPDLIN
      USE EIRMOD_REFLEC, ONLY: EIRENE_REFLC0,
cym variables that need to be allocated/associated for worker threads
cym will disappear when parallel zone will encompass the whole code
     .                         IREDUC, FREDUC, EREDUC
#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
      USE EIRMOD_REFLEC, ONLY: ZIDED, ZENGY, ZRANGE,
     .                         ZDE, ZDEL, ZR, ZIDE,
     .                         E0AV, QUOTE, QUOTR
#endif
      USE EIRMOD_PLT2D, ONLY: EIRENE_CHCTRC
      USE EIRMOD_OPENMP
      USE EIRMOD_INFCOP, ONLY: EIRENE_IF3COP, EIRENE_INFCOP_PRE_STRATA,
     .                         EIRENE_INFCOP_POST_STRATUM,
     .                         EIRENE_IF3COP_SUM
      USE EIRMOD_CALSTR_BUFFERED
      USE EIRMOD_BALANCED_STRATEGY
      USE EIRMOD_IMPROVED_STRATEGY
      USE EIRMOD_PRESSURELOOP

      IMPLICIT NONE
      PRIVATE
      PUBLIC :: EIRENE_MCARLO, EIRENE_MCARLO2

      REAL(DP), ALLOCATABLE, SAVE :: DUMMY(:),
     .                               ZVOLIN(:), ZVOLIW(:), SCLTAL(:,:)

      CONTAINS

c  nov.16th 2005: npts_save = npts always, not only for nlmovie option
c                 because otherwise in iterative mode a stratum cannot be
c                 re-activated, once it was deactivated in a particlar iteration.
c                 v.kotov
c  19.12.05:  bug: no printout of surface tally std. dev., for sum over strata
c             bug fix: here in mcarlo.f: sigmaw = stvw and sgmws=stvws added
c             also needed for this bug fix: clear_sumostra, stat_sumostra

!PB 02.03.06: storing of trajectories
!pb 08.11.06: definition of splitting arrays changed
!             RSPLST(NLEVEL,1:NPARTC) --> RSPLST(1:NPARTC,NLEVEL)
!             ISPLST(NLEVEL,1:MPARTC) --> ISPLST(1:MPARTC,NLEVEL)
!pb 01.12.06: open and close of fort.10 moved to WRSTRT
!pb 05.12.06: COLLECT_CENSUS introduced to allow for time-dependent mode in
!             parallel calculation
!pb 15.12.06: COLSUM replaced by COLLECT_COUTAU
!pb 18.12.06: call to PEDIST only done by processor 0 to avoid trouble because
!             of inaccuracies, call to broad_pedist needed to distribute
!             information calculated in pedist
!pb 05.02.07: copy ALGV to help array before integrating to avoid array bound violation
!pb 22.05.07: option introduced to set the random number seed after a specified
!             number of particles (used to check parallelization)
!pb 15.11.07: Meaning of NTCPU changed:
!             now NTCPU is the amount of cpu time used for particle tracing
!             times used for initialization and integration of result is not
!             taken into account
!   21.07.09: Meaning of XTIM changed: now XTIM is the time allocated for each stratum
!             no longer the end time
!dr 10.05.10: LOCAT0 might also turn off a stratum. Then: skip this is MCARLO, added after call to LOCAT0
cdr 22.09.14: upfcop only to be called in coupled mode: nmode.gt.0
cdr 22.09.14: cpu time output removed. To be collected and printout made conditional
cdr dec. 15 : 'upfcop.f' now 'updlin.f', moved from couple-specific part to main eirene code,
cdr           under scoring/updlin.f.
!pb 18.01.16: for totally random particle trajectories avoid usage of same random numbers in consecutive calls to mcarlo
cdr april 16: use: nstrai rather than nstra in do-loops. Bug fix from ITER-IO
!pb may 16  : for NPRS<NSTRAI call to IF3COP moved out of strata loop
cdr Nov. 20 : Do RANMAR random generator test only once, not per stratum.
cdr           It uses its own (Marsagia-) seed for initialization, independent
cdr           of eirene seed per stratum.
cdr fix ph6:  remove obsolete option LINE_CUTOFF here and from SAMVOL
c
      SUBROUTINE EIRENE_MCARLO
C
C  MONTE CARLO CALCULATION
C
      IMPLICIT NONE
C
      CHARACTER(6), SAVE :: CIS
      CHARACTER(10), SAVE :: CDATE, CTIME

      REAL(DP), ALLOCATABLE :: OUTAU(:)

      REAL(DP), ALLOCATABLE :: XTIM(:)

      REAL(DP) :: XFL1,
     .          XPRNLS, XFACT, XPRNLI,
     .          TIMI, XPT, XX1, XPT1, XFL, SECND, XX,
     .          ZW, ZWW, ZVOLWT, ZVOLNT, FSIG, ZFLUX,
     .          OVER, WTT, timan, timen,
     .          tim1, tim2, timst,
     .          rn1, ran, xmax
      REAL(DP), SAVE :: SECND1, SECND2, SECDEL
!     REAL(DP) :: OVER_ACC

cym IC ?
      INTEGER :: NPTS_SAVE(NSTRA), NINITL_SAVE(NSTRA), IPTSI
      INTEGER, SAVE :: ISDV, IALS, ISTRAA, ISTRAE, ICELL, ISTRAI,
     .           IENTRY, IGFFT, IALV, IDV, I, IRC, NMX,
     .           NINIST, IPANU, ISEED_ISTRA, ISEED_IPTSI, IDUMRAN,
     .           NDIGITS, NINIMAX, NAVAIL,
     .           ISTR, NPTTOT, NREC11,
     .           IADD, ICLV, IADV,
     .           INODES, J, IT, IMCP,
     .           ISUM, NPX, IS, NEW_ITER, ISPC, IN,
     .           JATM, JMOL, JION, JPHOT, JPLS,
     .           IERR, IFIRST
#ifdef CHECKBIN
      INTEGER :: ipb, jpb
#endif

      LOGICAL, SAVE :: LGSTOP
      LOGICAL :: LGABORT
cpg
      INTEGER, SAVE :: MY_ID

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP THREADPRIVATE(IPANU,LGSTOP,
!$OMP& MY_ID,NINIST,NINIMAX,CDATE,CTIME,
!$OMP& ISEED_ISTRA)
!$OMP THREADPRIVATE(CIS)
#endif
C     INTEGER :: N2
cdr 4.11.20
      INTEGER, SAVE :: ICO_CALL=0, IRAN_TEST=0

      real(dp) :: timused
      integer :: itimstart, itimend, itimrate

C     used only for IMPROVED strategy
      integer :: mapping(nstrai), order(nstrai), istr1, istr2, priority
C
      LOGICAL :: LOGHELP(NSTRA)
      LOGICAL :: LSTP
      EXTERNAL :: EIRENE_ALGTAL, EIRENE_BROAD_IESTR, EIRENE_CALSTR,
     .            EIRENE_CLEAR_STRATUM, EIRENE_CLEAR_SUMOSTRA,
     .            EIRENE_COLLECT_CENSUS, EIRENE_COLLECT_COUTAU,
     .            EIRENE_COLLECT_DATA_USR, EIRENE_EMISSIVITY,
     .            EIRENE_FOLION, EIRENE_FOLNEUT, EIRENE_FTCRI,
     .            EIRENE_INTEGRATE_TALLIES, EIRENE_INTTAL,
     .            EIRENE_LEER, EIRENE_MASAGE, EIRENE_MASBOX,
     .            EIRENE_MASJ1, EIRENE_MASJ1R, EIRENE_MASJ2,
     .            EIRENE_MASJ2R, EIRENE_MASR1, EIRENE_MASR5,
     .            EIRENE_NCELLN, EIRENE_OUTLST, EIRENE_REFLC0_PHOTON,
     .            EIRENE_SCAL_SURF_TALLIES, EIRENE_SCAL_VOLAV_TALLIES,
     .            EIRENE_SCALE_DEVIATION, EIRENE_SCALE_TALLIES,
     .            EIRENE_SET_SCAL_CONST, EIRENE_STAT_SUMOSTRA,
     .            EIRENE_STATS0_SPC, EIRENE_STATS1_SPC,
     .            EIRENE_STATS2_SPC, EIRENE_SUM_AVERAGE,
     .            EIRENE_SUMOSTRA, EIRENE_SYMET, EIRENE_WRSTRT,
     .            IOFLUSH_USR, RANMAR_TEST, EIRENE_CHECK_EXIT
C  OVERHEAD FOR POSTPROCESSING (SECONDS)
C      DATA N2/2/
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP BARRIER
#endif
      ALLOCATE(XTIM(0:NSTRA))
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP MASTER
#endif
      TIMI=EIRENE_SECOND_OWN()
      timan= timi
      tim1 = timi
      timst= timi
C
      IF (NFILEN.NE.0) THEN
        NREC11=NOUTAU
        OPEN (UNIT=11+ifoff,ACCESS='DIRECT',FORM='UNFORMATTED',
     .        RECL=8*NREC11)
      ENDIF

      IF (.NOT.ALLOCATED(DUMMY)) THEN
        ALLOCATE ( DUMMY(NRTAL),
     .             ZVOLIN(NRTAL),
     .             ZVOLIW(NRTAL),
     .             SCLTAL(N1MX,NTALV))
      END IF

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C-------------------------------------------------------------------
C
C** INITIALIZE SOME DATA AND SUBROUTINES (ONCE FOR ALL STRATA);
c   NEEDED FOR TRACING, STATISTICS AND SCALING. ONCE PER ITERATION OR STEP
C
C  SCLTAL: FLAG FOR SCALING OF VOLUME-AVERAGED TALLY
C  SCLTAL =0  1.
C         =1  ZVOLIN(ICELL)
C         =2  ZW
C         =3  ZVOLIW(ICELL)
C         =4  ZWW
C
      SCLTAL=0.D0
      SCLTAL(1,1:8)=1
      SCLTAL(1,9:56)=3
C
C  INITIALIZE RANDOM NUMBER ARRAYS
      INIV1=0
      INIV2=0
      INIV3=0
      INIV4=0
C  DETERMINE MAXIMAL INTEGER (DEPENDING ON MACHINE)
      IF (NLCRR) THEN
        INTMAX=HUGE(1)
      ENDIF
C
C  IRNDVC MUST BE EVEN AND NOT LARGER THAN 64 (COMMON CRAND)
C  IRNDVC IS THE NUMBER OF RANDOM VECTORS PRODUCED IN ONE CALL TO
C  TO RANDOM SAMPLING ROUTINES
      IF (NLCRR) THEN
        IRNDVC=2
      ELSE
        IRNDVC=64
      ENDIF
      IRNDVH=IRNDVC/2

      TIMen=EIRENE_SECOND_OWN()
      IF (TRCTIM)
     > write (iunout,*) 'CPU time for init of mcarlo ', timen-tim1
      tim1 = timen
C
C  INITIALIZE SUBR. STATIS
C
      CALL EIRENE_LEER(1)
      CALL EIRENE_STATS0
      TIM2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) 'CPU time for stats0 ', tim2-tim1
      tim1 = tim2
      CALL EIRENE_STATS0_SPC
      TIM2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) 'CPU time for stats0_spc ', tim2-tim1
      tim1 = tim2
C  INITIALIZE SUBR. REFLEC AND SPUTER
      CALL EIRENE_REFLC0
      TIM2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) 'CPU time for reflc0 ', tim2-tim1
      tim1 = tim2
      IF (NPHOT > 0) THEN
        CALL EIRENE_REFLC0_PHOTON
        TIM2=EIRENE_SECOND_OWN()
        IF (TRCTIM)
     >   write (iunout,*) 'CPU time for reflc0_photon ', tim2-tim1
        tim1 = tim2
      END IF
      CALL EIRENE_SPUTR0
      TIM2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) 'CPU time for sputr0 ', tim2-tim1
      tim1 = tim2
C  INITIALIZE SUBR. SAMVOL
      CALL EIRENE_SAMVL0
      TIM2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) 'CPU time for samvl0 ', tim2-tim1
      tim1 = tim2
C  INITIALIZE SUBR. SAMSRF
      CALL EIRENE_SAMSF0
      TIM2=EIRENE_SECOND_OWN()
      IF (TRCTIM) write (iunout,*) 'CPU time for samsf0 ', tim2-tim1
      tim1 = tim2
C  INITIALIZE SUBR. UPDLIN
      CALL EIRENE_PREPARE_UPDLIN
      TIM2=EIRENE_SECOND_OWN()
      IF (TRCTIM)
     >  write (iunout,*) 'CPU time for prepare_updlin ', tim2-tim1
      tim1 = tim2
!pb moved EIRENE_INFCOP_PRE_MCARLO called from EIRENE_EIRENE
csw 28jan2011 extra B25, produce the files for surface properties visualisation if required
!pb   call eirene_extrab25_srfprvsl
csw
cxpb zero out the output arrays to B2.5
!pb   call eirene_wneutrals_clean(.false.)
cxpb
C
C
      IESTR=-1
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP END MASTER
#endif
      IF (NFILEN.EQ.2.OR.NFILEN.EQ.7) GOTO 2000
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP MASTER
#endif
C
C**** CLEAR WORK AREA FOR SUM OVER STRATA ****************************
C
      CALL EIRENE_CLEAR_SUMOSTRA

c now initialized in eirene_init_cspez
!pb      LOGATM=.FALSE.
!pb      LOGION=.FALSE.
!pb      LOGMOL=.FALSE.
!pb      LOGPLS=.FALSE.
!pb      LOGPHOT=.FALSE.
C
C
C   MAXIMAL CALCULATION TIME ALLOWED FOR EACH STRATUM,
C   PROPORTIONAL TO NPTS(ISTRA), OR FLUX(ISTRA) (INPUT)
C   OR LINEAR COMBINATION THEREOF
C   THEREFORE NUMBER OF TEST PARTICLES MAY BE LESS THAN NPTS
C   BUT DO AT LEAST 2 PARTICLES, IN CASE NPTS(ISTRA).GE.2
C
      XX = NTCPU
CVKMPI      XTIM(0)=EIRENE_SECOND_OWN()
CVKMPI      SECND=XTIM(0)
      SECND=EIRENE_SECOND_OWN()

      NPTS_SAVE=NPTS
      NINITL_SAVE = NINITL
!pb 28012016
!   count number of times MCARLO has been called
      ICO_CALL = ICO_CALL + 1

      timan=secnd   ! wall clock time, after initial overhead
C
C  REMAINING CPU TIME, SUBTRACT N2 SECONDS FOR PRINTOUT AND PLOTS
!pb   XX1=XX-N2

C  CHANGED: use XX=NTCPU seconds of cpu time for calculation of trajectories
      XX1 = XX

      XPT=0.
      XFL=0.
      DO 7 ISTRAI=1,NSTRAI
        ISTRA=ISTRAI
        NLSRON(ISTRA) = .TRUE.
        IF (NPTS(ISTRA).LE.0.AND.FLUX(ISTRA).GT.0.D0) THEN
          FLUX(ISTRA)=0.D0
          NLSRON(ISTRA) = .FALSE.
          WRITE (iunout,'(1x,a,i3,a)') 'STRATUM ISTRA= ',ISTRA,
     .                     ' TURNED OFF, BECAUSE NPTS=0'
          CALL EIRENE_LEER(1)
        ELSEIF (NPTS(ISTRA).GT.0 .AND. FLUX(ISTRA).LE.0.D0) THEN
          NPTS(ISTRA)=0
          NLSRON(ISTRA) = .FALSE.
          WRITE (iunout,'(1x,a,i3,a)') 'STRATUM ISTRA= ',ISTRA,
     .                     ' TURNED OFF, BECAUSE FLUX=0.0'
          CALL EIRENE_LEER(1)
        ELSEIF (NPTS(ISTRA).LE.0 .AND. FLUX(ISTRA).LE.0.D0) THEN
          NPTS(ISTRA)=0
          NLSRON(ISTRA) = .FALSE.
          WRITE (iunout,'(1x,a,i3,a)') 'STRATUM ISTRA= ',ISTRA,
     .                     ' TURNED OFF, BECAUSE NPTS=0 AND FLUX=0.0'
          CALL EIRENE_LEER(1)
        ENDIF
        IF (SUM(SORWGT(1:NSRFSI(ISTRA),ISTRA)).LE.0.D0) THEN
          NPTS(ISTRA)=0
          NLSRON(ISTRA)=.FALSE.
          WRITE (iunout,*) 'STRATUM ISTRA= ',ISTRA,
     .                     ' TURNED OFF, BECAUSE',
     .                     ' THE SUM OF THE FLUXES'
          WRITE (iunout,*) 'FROM THE SUBSTRATA DEFINED BY',
     .                     ' SORWGT(SUBSTRATUM,STRATUM) IS .LE. ZERO'
          CALL EIRENE_LEER(1)
        ENDIF
        IF (NPTS(ISTRA).LE.0.OR.FLUX(ISTRA).LE.0.D0)
     .     NLSRON(ISTRA) = .FALSE.
        XPT=XPT+FLOAT(NPTS(ISTRA))
        XFL=XFL+FLUX(ISTRA)
    7 CONTINUE
      XPT1=0.
      XFL1=0.

      xtim = 0._dp
      DO 8 ISTRAI=1,NSTRAI
        IF (NLSRON(ISTRAI)) THEN
          XPT1=NPTS(ISTRAI) !VKMPI
          XFL1=FLUX(ISTRAI) !VKMPI
          XTIM(ISTRAI)=XX1*((1.-ALLOC)*XPT1/(XPT+EPS60)+
     +                      (   ALLOC)*XFL1/(XFL+EPS60)) !VKMPI
        ELSE
          xtim(istrai)=0.0 !VKMPI
        END IF
    8 CONTINUE

C  REDISTRIBUTE XTIM IN CASE THAT SOURCES ARE TEMPORARILY SWITCHED OFF (SHORT CYCLE)
CVKMPI      DO ISTRAI=1,NSTRAI
CVKMPI        DXTIM(ISTRAI)=XTIM(ISTRAI)-XTIM(ISTRAI-1)
CVKMPI        IF (.NOT.NLSRON(ISTRAI)) DXTIM(ISTRAI)=0._DP
CVKMPI      END DO

CVKMPI      DO ISTRAI=1,NSTRAI
CVKMPI        XTIM(ISTRAI)=XTIM(ISTRAI-1)+DXTIM(ISTRAI)
CVKMPI      END DO

      XTIM(0)=SUM(XTIM(1:NSTRAI))
C
      TIMen=EIRENE_SECOND_OWN()

      CALL EIRENE_LEER(2)
      CALL EIRENE_MASAGE('LOOP OVER STRATA STARTS AT CPU TIME (SEC) :')
      CALL EIRENE_MASR1 ('STARTTIM',SECND)
      CALL EIRENE_MASAGE('CPU TIME ASSIGNED TO STRATA (SEC) :')
      IF (ALLOC.EQ.0.D0) THEN
        CALL EIRENE_MASAGE('PROPORTIONAL NPTS(ISTRA)')
      ELSEIF (ALLOC.EQ.1.) THEN
        CALL EIRENE_MASAGE('PROPORTIONAL FLUX(ISTRA)')
      ELSE
        CALL EIRENE_MASAGE('WEIGHTED ALLOCATION BETWEEN NPTS AND FLUX')
      ENDIF
      DO 9 ISTRAI=1,NSTRAI
CVKMPI        DELT=XTIM(ISTRAI)-XTIM(ISTRAI-1)
CVKMPI        CALL EIRENE_MASJ1R ('STRATUM, TIME   ',ISTRAI,DELT)
        CALL EIRENE_MASJ1R ('STRATUM, TIME   ',ISTRAI,XTIM(ISTRAI))
CVKMPI
    9 CONTINUE
      CALL EIRENE_LEER(2)
C
C  ASSIGN NUMBER OF PARTICLES TO BE STORED ON CENSUS, PROPORTIONAL
C  TO CPU TIME ASSIGNED TO EACH STRATUM
C
      IF (NPRNLI.GT.0) THEN
        WRITE(iunout,'(1x,a,a)')
     .   'MAXIMUM NUMBER OF PARTICLES THAT WILL BE SAVED ',
     .   'ON CENSUS (TIME DEP MODE):'
        WRITE(iunout,*) 'PROP. TO CPU TIME ALLOCATED FOR EACH STRATUM'
        DO ISTRAI=1,NSTRAI
          XFACT=XTIM(ISTRAI)/XX1 !VKMPI
          XPRNLS        =NPRNLI*XFACT+0.5
          NPRNLS(ISTRAI)=INT(XPRNLS)
        ENDDO
   10   ISUM=SUM(NPRNLS(1:NSTRAI))
        IF (ISUM.NE.NPRNLI) THEN
C  ROUND-OFF ERRORS
          NMX=0
          NPX=-1
          DO ISTRAI=1,NSTRAI
            IF (NPRNLS(ISTRAI).GT.NPX) THEN
              NMX=ISTRAI
              NPX=NPRNLS(ISTRAI)
            ENDIF
          ENDDO
          IS=ISIGN(1,ISUM-NPRNLI)
          NPRNLS(NMX)=NPRNLS(NMX)-IS
          GOTO 10
        ENDIF
        WRITE (iunout,*) 'TOTAL CENSUS NPRNLI ',NPRNLI
        CALL EIRENE_LEER(1)
        DO ISTRAI=1,NSTRAI
          CALL EIRENE_MASJ2 ('STRATUM, NUMBER ',ISTRAI,NPRNLS(ISTRAI))
        ENDDO
      ENDIF

      IF (NPTST.LT.0.OR.NLMOVIE) THEN
!tf If the number of particles has to change, then we should do it here,
!tf before the parallel initialization
c
C  revert sequence of strata, so that census stratum is dealt with first
c  to ensure: ALL particles from census are re-launched, one by one (not: by sampling)
!tf  We are outside the stata loop, so we do not need the following if condition
!tf          ISTRA=NSTRAI-ISTR+1
!tf          IF (ISTRA.EQ.NSTRAI-1) THEN
C  TOTAL STORAGE STILL AVAILABLE ON NEW CENSUS
C  AFTER ONE TO ONE RESTART FROM OLD CENSUS IS COMPLETED
#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
!$OMP ATOMIC WRITE
#endif
        NPTTOT=NPRNLI-IPRNLI
C  REDEFINE NPTS ACCORDING TO XTIM(ISTRA)
        CALL EIRENE_LEER(2)
        WRITE (iunout,*)
     .   'REDEFINE NPTS FOR ONE-BY-ONE RELAUNCH FROM CENSUS'
        ISUM=0
        DO IS=1,NSTRAI-1
          XFACT=XTIM(IS)/XTIM(0)  !PB  XTIM(IS):
                                  !    CPU TIME ASSIGNED TO STRATUM IS
          XPRNLI=NPTTOT*XFACT+0.5
          NPTS(IS)=INT(XPRNLI)
          ISUM=ISUM+NPTS(IS)
          WRITE(iunout,*) 'ISTRA, NPTS = ',IS,NPTS(IS)
        ENDDO
!tf          ENDIF
      ENDIF
C
C  ASSIGN PEs TO STRATA
C
      call init_work_distribution(xtim,xx1,input_distribution_strategy)
      if (nmode > 0) call eirene_infcop_pre_strata

      if (nprs > 1) then
c nlident: jeder proc. von einer quelle istra bekommt gleichen seed gem. ninitl(istra).
c         erzeugt bei zwei gleichen quellen (istra) identische ergebnisse.
c not nlident: ninitl wird auf dem processor geaendert, add my_pe*10000

CHJL/YM NLIDENT code was moved from here to within the threaded region
        NPRNLS=NPRNLI
      ENDIF
C
C**** INITIALIZE COMMONS COUTAU AND CSPEZ
C
csw 19mar2013 moved to here after call to pedist (xmct/xmcp)
cdr  presumably because pedist uses xmct,xmcp from previous cycle with
cdr  external code.
cdr  In pedist.f we currently hope that COUTAU has not been deallocated
cdr  between the present and the previous cycle.
!pb copy NLSRON to LOGHELP to avoid warnings from Intel compiler
!pb      CALL EIRENE_INIT_COUTAU(NLSRON)

      LOGHELP(1:NSTRA) = NLSRON(1:NSTRA)
      CALL EIRENE_INIT_COUTAU(LOGHELP)
      XMCT(0)=0.
      FASCL(:,0)=1.
      FMSCL(:,0)=1.
      FISCL(:,0)=1.
      FPHSCL(:,0)=1.

cdr feb 2020
      if (any(ninitl(1:nstrai) .lt. 0)) then
        write (iunout,*) 'negative random seeds found'
        if (nlcrr)
     .    write (iunout,*) 'correlated sampling turned off'
        nlcrr=.false.
      else
        ninimax=0
        do istr=1,nstrai
          ninimax = max(ninimax,ninitl(istr)+my_pe*10000)
        end do
        write (iunout,*) 'largest preselected random seed ',ninimax
cdr  0 (zero) is now a perfectly legal seed
        xmax=max(ninimax-1.0_dp,1._dp)
        ndigits=int(log10(xmax))
cdr  Do not touch the last NDIGITS digits for seeding
cdr  produced by code.
cdr  With RANMAR, we can use seeds up to 1e9 (precisely: 950.000.000) seeds
cdr  We block the last NDIGITS digits of the seed,
cdr  so for each NINITL, we have NAVAIL further save seeds,
cdr  e.g. for correlated sampling or for multi-processor runs.
        NAVAIL=950000000/10**ndigits  ! available seeds per NINITL,
cdr  for seeding multiple processors per NINITL, or correlated sampling
      endif
cdr
C
C
C**** STRATA LOOP ****************************************************
C
CHJL Need to check this removal of over_acc
!      OVER_ACC=0.D0
      NEW_ITER=0
C
C     ordering of the strata for IMPROVED strategy
      if (work_distribution_strategy == STRATEGY_IMPROVED) then
        do  istr1 = 1, nstrai
          if (.not. calc_stratum(istr1)) then
            order(istr1) = 1
          elseif (nparts_loc(istr1) == npts(istr1)) then
            order(istr1) = 2
          elseif (.not. i_am_leader(istr1)) then
            order(istr1) = 3
          else
            order(istr1) = 4
          endif
        enddo
        istr2 = 1
        do priority = 1, 4
          do istr1 = 1, nstrai
            if (order(istr1) == priority) then
              mapping(istr2) = istr1
              istr2 = istr2 + 1
            endif
          enddo
        enddo
      endif
C
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP END MASTER
!$OMP BARRIER
#endif
! Initialise openMP including allocating THREADPRIVATE arrays
      CALL EIRENE_PREPARE_OPENMP()
      NPANU=0
#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
!$OMP  PARALLEL DEFAULT(SHARED)
!$OMP& COPYIN(ETH,Q,M2M1,ES,ETF,ISDVI,IPSTD,RPST,
!$OMP& RCMSPL,ICMSPL,LCMSOU,XSTOR,XSTORV,
!$OMP& ZIDED,ZENGY,ZRANGE,ZDE,ZDEL,ZR,ZIDE,E0AV,QUOTE,QUOTR,
!$OMP& FREC,VSOURC,VSMXI,RQ21,PS21,ASIMP,
!$OMP& ISOURC,ICMX,IFREC,LPLSSR,LATSSR)
!$OMP& PRIVATE(ISTR)
#endif
      DO ISTR=1,NSTRAI   ! main loop over strata

        timan=EIRENE_second_own()

C  SPECIAL TREATMENT FOR MOVIE OPTION, OR FOR ONE-BY ONE RELAUNCH FROM CENSUS ARRAY
C  IN TIME DEP. MODE
        IF (NPTST.LT.0.OR.NLMOVIE) THEN
ctf Changing the number of particles moved before the parallel initialization
c
C  revert sequence of strata, so that census stratum is dealt with first
c  to ensure: ALL particles from census are re-launched, one by one (not: by sampling)
          ISTRA=NSTRAI-ISTR+1
cdr  now istr=1, istra=nstrai, is the stratum from old census
          IF (ISTRA.EQ.NSTRAI-1) THEN
cdr  something wrong here: this was originally stratum ISTR=2.
C  TOTAL STORAGE STILL AVAILABLE ON NEW CENSUS
C  AFTER ONE TO ONE RESTART FROM OLD CENSUS IS COMPLETED
            NPTTOT=NPRNLI-IPRNLI
C  REDEFINE NPTS ACCORDING TO XTIM(ISTRA)
            CALL EIRENE_LEER(2)
            WRITE (iunout,*)
     .        'REDEFINE NPTS FOR OF ONE-BY-ONE RELAUNCH FROM CENSUS'
            ISUM=0
            DO IS=1,NSTRAI-1
              XFACT=XTIM(IS)/XTIM(0)  !PB  XTIM(IS):
                                      ! CPU TIME ASSIGNED TO STRATUM IS
              XPRNLI=NPTTOT*XFACT+0.5
              NPTS(IS)=INT(XPRNLI)
              ISUM=ISUM+NPTS(IS)
              WRITE(iunout,*) 'ISTRA, NPTS = ',IS,NPTS(IS)
            ENDDO
          ENDIF
        ELSE
          if (work_distribution_strategy == STRATEGY_IMPROVED) then
            istra = mapping(istr)
          else
            istra = istr
          endif
        ENDIF

C  MOVIE OPTION (NLMOVIE): DONE,
C    if     nlmovie: sequence of strata is reversed, census stratum istra=nstrai comes first!
C                    one by one re-launch of ALL particles from census
c    if not nlmovie: census stratum istra=nstrai comes last.
#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
        CALL EIRENE_CHECK_EXIT
!$OMP BARRIER
#endif
        IF (.NOT.NLSRON(ISTRA)) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP MASTER
#endif
          CALL EIRENE_LEER(2)
          WRITE (iunout,*) 'STRATUM NO. ',ISTRA,' ABANDONED, NLSRON=F'
          CALL EIRENE_LEER(2)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END MASTER
#endif
          CYCLE
        ENDIF

        IF (CALC_STRATUM(ISTRA)) THEN
          CALL EIRENE_LEER(2)
          WRITE (iunout,*) 'BEGIN TO WORK ON STRATUM NO. ',ISTRA
          CALL EIRENE_LEER(2)
          XMCP(ISTRA)=0.
c??
          XMCT(ISTRA)=0.
          IPANU=0
C
C  INITIALIZE RANDOM NUMBER GENERATOR FOR STRATUM ISTRA
cym test for multiple strata / restores the status after call to REFLC0
cym ideally would need initializing values on the master thread ...
          EREDUC=0._dp
          FREDUC=0._dp
          IREDUC=0

C  find random number generator seed, from input flag NINITL(ISTRA)
          IF (NINITL(ISTRA).GT.0) THEN
            IF (.NOT.NLIDENT) THEN
              MY_ID= EIRENE_ITHREAD + MY_PE*EIRENE_NTHREADS
              NINIST=NINITL(ISTRA)+MY_ID*10000 !one seed per thread
                                     !to generate independent numbers
            ELSE

              CALL EIRENE_LEER(1)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP SINGLE
#endif
          WRITE(IUNOUT,*)'.........................................'
          WRITE(IUNOUT,*)'NLIDENT: '
          WRITE(IUNOUT,*)'DEBUG MODE FOR PARALLELIZATION IS ACTIVE'
          WRITE(IUNOUT,*)'IF MULTIPLE CORES PER STRATUM, THEN ALL'
          WRITE(IUNOUT,*)'ASSIGNED CORES KEEP IDENTICAL RANDOM SEED.'
          WRITE(IUNOUT,*)'FOR ANY GIVEN STRATUM ISTRA, ALL NCIS CORES'
          WRITE(IUNOUT,*)'ASSIGNED TO ISTRA MUST PRODUCE IDENTICAL'
          WRITE(IUNOUT,*)'OUTPUT. ALSO VARIANCES PER STRATUM MUST'
          WRITE(IUNOUT,*)'SCALE EXACTLY WITH 1/NCIS(ISTRA),'
          WRITE(IUNOUT,*)'NOT ONLY ON STATISTICAL AVERAGE'
          WRITE(IUNOUT,*)'.........................................'

              CALL EIRENE_LEER(1)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END SINGLE
#endif
cym same seed for everybody (all threads of all processes -> debug)
              NINIST=NINITL(ISTRA)
            ENDIF

c  initialize random number generator with chosen input seed NINIST
c  ranset checks, if this is a legal seed for a particular generator,
c  and otherwise enforces that or stops the run.
            ISEED_ISTRA=RANSET_EIRENE(NINIST)
cdr
cdr  RANMAR is a fully portable generator.
cdr  Here we test if that is the case indeed in this run, on this machine.
cdr  Testing once should suffice.
            iadd=0
            ierr=0
            if (nlranmar .and. iran_test.eq.0) then
             call ranmar_test(iadd,ierr,ran)
             iran_test=iran_test+1
             if (ierr.eq.0) then
              write (iunout,*) 'RANMAR test passed.'
              if (IADD.eq.0) then
                write (iunout,*) 'Next random number will be: ',RAN
              else
                write (iunout,*) 'Skip next IADD rand. numbers: ',IADD
                write (iunout,*) 'Last random number used in test: ',RAN
              endif
             else
              write (iunout,*) 'RANMAR test failed. Error code:'
              write (iunout,'(I7)') IERR
              write (iunout,*) 'Last random number used in test: ',RAN
             endif
            endif
            call eirene_leer(1)

c  find random number seed from truly random procedure from wall clock time (use date and time)
cdr  This is quite arbitrary. And probably even wrong:
cdr  High chances of unintended correlations (the "birthday paradox"...)
          ELSEIF (NINITL(ISTRA).LT.0) THEN
cdr  format of CDATE: hhmmss.xxx
            CALL DATE_AND_TIME(CDATE,CTIME)
cdr we only use: hhmmss: i.e. a number between 0 and 235959,
cdr              in total: only about different 80000 seeds left
            READ(CTIME(1:6),*) NINITL(ISTRA)
!pb 28012016
!  add number of calls to MCARLO in order to avoid same random seeds in very short
!  cycles (within the same CDATE value still) with an external code
            NINITL(ISTRA) = NINITL(ISTRA) + ICO_CALL

            IF (.NOT.NLIDENT) THEN
              MY_ID= EIRENE_ITHREAD + MY_PE*EIRENE_NTHREADS
              NINIST=NINITL(ISTRA)+MY_ID*10000
            ELSE
              NINIST=NINITL(ISTRA)
            ENDIF
            WRITE (iunout,*) 'NINITL(ISTRA) SET TO ',NINIST
            ISEED_ISTRA=RANSET_EIRENE(NINIST)

          ELSEIF (NINITL(ISTRA).EQ.0) THEN
C  DO NOT RE-INITIALIZE RANDOM GENERATOR FOR THIS STRATUM, NOTHING TO BE DONE HERE
C  INTERNAL DEFAULT FIRST SEED IS TAKEN FOR FIRST STRATUM. FROM THEN ON: NO FURTHER SEEDING.
            ISEED_ISTRA=RANSET_EIRENE(0)

          ENDIF

          IF (TRCRNF) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
            WRITE (iunout,*) 'INITIALIZE RANDOM NUMBERS FOR STRATUM ',
     .                       'ISTRA= ',ISTRA,'thread# = ',EIRENE_ITHREAD
            WRITE (iunout,*) 'NINITL(ISTRA) SET TO ',NINITL(ISTRA)
            WRITE (iunout,*) 'ISEED_ISTRA (LEGAL SEED, AS USED) ',
     .                        ISEED_ISTRA
            CALL EIRENE_LEER(1)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
          ENDIF

c  remove remaining old generated random number vectors from earlier strata
          INIV1=0
          INIV2=0
          INIV3=0
          INIV4=0

C  ISEED_ISTRA IS SET NOW
C
          FASCL(:,ISTRA)=1.
          FMSCL(:,ISTRA)=1.
          FISCL(:,ISTRA)=1.
          FPHSCL(:,ISTRA)=1.

          LOGATM(:,ISTRA)=.FALSE.
          LOGION(:,ISTRA)=.FALSE.
          LOGMOL(:,ISTRA)=.FALSE.
          LOGPLS(:,ISTRA)=.FALSE.
          LOGPHOT(:,ISTRA)=.FALSE.

          timen=EIRENE_second_own()
C
C  CLEAR WORK AREA FOR THIS STRATUM
C
cym check whether that is needed on all threads ?
cym previously called by master thread only
          CALL EIRENE_CLEAR_STRATUM
C
C  ENFORCE TOROIDAL OR POLOIDAL SYMMETRY FOR THIS STRATUM
C
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP MASTER
#endif
          IF (NLAVRP(ISTRA)) THEN
            NLPOL=.FALSE.
          ENDIF
C
          IF (NLAVRT(ISTRA)) THEN
            NLTOR=.FALSE.
          ENDIF
C
C Should this not go into EIRENE_CLEAR_STRATUM as it resets a quantity
C for this stratum?
          IPRNLS=0
C
          LVNMXTR=NVNMXTR.GT.0.AND.KVNMXTR(ISTRA).LT.NVNMXTR
          KCFF=0
          XVNMX=0._DP
          XVNMX2=0._DP
C
C This will never happen as NLSRON status have not changed...
cym fortunate that it 'never happens' because not compatible with OpenMP
cym if NLSRON status is 'fixed' would require another way of implementing
cym          IF (.NOT.NLSRON(ISTRA)) CYCLE
cym ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C  INITIALIZE SUBR. LOCATE
C
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END MASTER
#endif
          CALL EIRENE_LOCAT0
C
C  LOCATE AND FOLLOW MC-PARTICLES
C
          CALL EIRENE_FTCRI(ISTRA,CIS)
          CALL EIRENE_MASBOX
     .     ('LAUNCH PARTICLES FOR STRATUM NUMBER ISTRA='//CIS)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP MASTER
#endif
          OVER=EIRENE_SECOND_OWN()-SECND
C  ACCUMULATED OVERHEAD BETWEEN STRATA
CHJL Need to check this removal of over_acc
!          OVER_ACC=OVER_ACC+OVER
CVKMPI        CALL EIRENE_MASR1 ('OVERHEAD',OVER)
CVKMPI        XTIM(ISTRA)=XTIM(ISTRA)+OVER_ACC
          WRITE (iunout,*) 'XTIM(ISTRA)= ',XTIM(ISTRA)

          TIMI=EIRENE_SECOND_OWN()     !VKMPI
          XTIM(ISTRA)=XTIM(ISTRA)+TIMI !VKMPI
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END MASTER
#endif
C
          LGABORT=.FALSE.

cdr  LGLAST = T: LAST TRAJECTORY OF PRESENT STRATUM ISTRA
          LGLAST=.FALSE.
cdr  LGSTOP = T: same as lglast, but only due to npts or cpu time criterion
cdr  LGLAST      may also have been set during particle tracking, for other reasons.
cdr              presently: e.g. if census array is full, set in TIMCOL
          LGSTOP=.FALSE.
C
C
csw 19feb2013
!pb 03122013        timstart=mpi_wtime()
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP MASTER
#endif
          call system_clock (itimstart, itimrate)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END MASTER
#endif
csw
#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
          call eirene_check_exit
!$OMP BARRIER
#endif
C  PARTICLE LOOP WITHIN STRATUM ISTRA
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP  DO REDUCTION(.OR.:LOGATM,LOGMOL,LOGION,LOGPLS,
!$OMP&                   LOGPHOT,LMETSPW,LMETSP)
!$OMP& REDUCTION(+:WTOTA,WTOTM,WTOTI,WTOTPH,WTOTP,WTOTE,
!$OMP&             ETOTA,ETOTM,ETOTI,ETOTPH,ETOTP,EELFI,XMCP,
!$OMP&             PPATI,PPMLI,PPIOI,PPPHTI,
!$OMP&             EPATI,EPMLI,EPIOI,EPPHTI,EPPLI,
!$OMP&             PAATI,PAMLI,PAIOI,PRFAAI,PRFAMI,PRFAII,
!$OMP&             ESTIMV,ESTIMS)
!$OMP& PRIVATE(I,J,IN,ISPC,INODES,ISEED_IPTSI,IDUMRAN)
!$OMP& PRIVATE(SECND1,SECND2,SECDEL)
#endif
          DO 100 IPTSI=1,NPARTS_LOC(ISTRA)

! Check that another thread has not aborted the particle loop
!pb#ifdef USE_OPENMP
            IF(LGABORT) cycle
!pb#endif

            call calstr_progress_message(iptsi, istra, nparts_loc)
C  SOME PREPARATORY WORK, ONCE FOR EACH NEW PARTICLE HISTORY
C
C  RE-INITIALIZE INDEX ARRAYS: VISITED CELLS, VISITED WALL SEGMENTS
            NCLMT = 0
            DO I=1,NCLMTS
              IN=ICLMT(I)
              IMETCL(IN) = 0
            END DO
c  LMETSP: array for 1st ("species") index of volume-averaged tallies,
c  which is scored along a trajectory
            LMETSP=.FALSE.
            NCLMTS = 0

            DO I=1,NWLMT
              IN=IWLMT(I)
              IMETWL(IN) = 0
            END DO
c  LMETSPW: array for 1st ("species") index of surface-averaged tallies,
c  which is scored along a trajectory
            LMETSPW=.FALSE.
            NWLMT = 0
            NWLMTS = 0

            IF (NADSPC > 0) THEN
              DO ISPC=1,NADSPC
                ESTIML(ISPC)%IMETSP = 0
              END DO
            END IF
C...........................................................................

C LGSTOP is always equal LGLAST, see line 681
C Should not this be (.NOT.LGLAST.AND.LGSTOP)?

cym cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cym this is NOT yet working with OpenMP - goto 101 removed
cym correct solution may involve an IF enclosing the particle loop
cym cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            IF (LGLAST.AND.LGSTOP) THEN
              CALL EIRENE_LEER(1)
              WRITE (iunout,*)
     .          'NO FURTHER COMP. TIME AVAIL. FOR THIS STRATUM'
              WRITE (iunout,'(1x,a,a)')
     .          'M.C. HISTORIES FOLLOWED UNTIL THAT TIME FOR ',
     .          'THIS STRATUM'
              call system_clock (itimend, itimrate)
              timused=real(itimend-itimstart,DP)/REAL(itimrate,DP)
              CALL EIRENE_MASJ2R('ISTRA,IPANU,TIMUSED     ',
     .                            ISTRA,IPANU,TIMUSED)
              IF (NPRNLI.GT.0) THEN
                WRITE (iunout,*) 'M.C. HISTORIES THAT SCORED AT CENSUS'
                CALL EIRENE_MASJ1 ('IPRNLS= ',IPRNLS)
              ENDIF
              IF (TRCLST) CALL EIRENE_OUTLST
cym cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cym exit from the loop not allowed with OpenMP
cym currently time limit not enforced
cym cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!pb#ifdef USE_OPENMP
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC WRITE
#endif
              LGABORT = .TRUE.   ! Abort this stratum,
!pb#else
!pb!pb17.11.2021              GOTO 101
!pb              LGABORT = .TRUE.   ! Abort this stratum,
!pb#endif
            ELSEIF (LGLAST.AND..NOT.LGSTOP) THEN
              CALL EIRENE_LEER(1)
              WRITE (iunout,*) 'CENSUS ARRAYS FILLED FOR THIS STRATUM'
              WRITE (iunout,'(1x,a,a)')
     .          'M.C. HISTORIES FOLLOWED UNTIL THAT TIME FOR ',
     .          'THIS STRATUM'
              call system_clock (itimend, itimrate)
              timused=real(itimend-itimstart,DP)/REAL(itimrate,DP)
              CALL EIRENE_MASJ2R('ISTRA,IPANU,TIMUSED     ',
     .                            ISTRA,IPANU,TIMUSED)
              WRITE (iunout,*) 'M.C. HISTORIES THAT SCORED AT CENSUS'
              CALL EIRENE_MASJ1 ('IPRNLS= ',IPRNLS)
              IF (TRCLST) CALL EIRENE_OUTLST
cym cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cym see comment above
cym cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!pb#ifdef USE_OPENMP
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP ATOMIC WRITE
#endif
              LGABORT = .TRUE.   ! Abort this stratum,
!pb#else
!pb!pb17.11.2021              GOTO 101
!pb              LGABORT = .TRUE.   ! Abort this stratum
!pb#endif
            ENDIF

C  WALL CLOCK TIME AT START OF NEXT MONTE CARLO HISTORY
            SECND1=EIRENE_SECOND_OWN()
C
C  LAST HISTORY FOR PRESENT STRATUM ?
! all requested particles done
            LGLAST = IPTSI.EQ.NPTS(ISTRA)
! cpu limit reached and minimum no. of part. done
            LGLAST = LGLAST.OR.(SECND1.GT.XTIM(ISTRA).AND.
     .                          IPTSI.GE.NMINPTS(ISTRA).AND.
     .                          .NOT.NLMOVIE)
CDR         LGLAST = LGLAST.OR.(CENSUS FILLED ?)  CURRENTLY DONE IN TIMCOL

            LGSTOP = LGLAST

C.......................................................................
C  CORRELATED SAMPLING: CREATE A RANDOM NUMBER GENERATOR SEED FOR NEXT PARTICLE
C  FROM THE SEED USED FOR THE CURRENT PARTICLE
            IF (NLCRR) THEN
C
C  RE-INITIALIZE RANDOM NUMBERS FOR PARTICLE IPTSI, TO GENERATE CORRELATION
C
c  current seed within current stratum is iseed_istra.

c  NLCRR: get new tentative seed. Build a fixed sequence
c         of seeds starting from NINITL,
c         then re-initialize each particle seed from that sequence.
cdr  Compose a seed  iseed=AAABBB, where BBB is fixed seed NINITL per stratum
cdr  and AAA a fixed sequence for all the particles in this stratum.
              iseed_iptsi=iseed_istra+(10**ndigits+iptsi)
              idumran=ranset_eirene(iseed_iptsi)
! now we have the seed iseed_iptsi to start the history.

C  FOR TEST ONLY: PRINT FIRST RANDOM NUMBER PER TRAJECTORY
              IF (TRCRNF) THEN
                call eirene_leer(1)
                write (iunout,*) 'new particle ',iptsi
                RN1=RANF_EIRENE()  ! sacrifice one random number
                                   ! for testing random sequence
                write (iunout,*) 'iseed,iseed_next,rn',
     .                            iseed_iptsi,iseed_istra,RN1
              ENDIF

c  derive one more seed, for reflec.f. cdr: unfinished....
cdr           ISEEDR=INT(ISEED_ISTRA*0.3D0)

              INIV1=0
              INIV2=0
              INIV3=0
              INIV4=0
C......................................................................


            ELSE IF (
     .             (MY_PE == 0) .AND.
     .             (NPTSDEL(ISTRA).GT.0) .AND.
     .             (NINITL(ISTRA).GT.0) .AND.
     .             (.NOT.NLCRR)
     .                                  ) THEN
c  Proprietary option for internal testing of MPI parallelization only:
c  NPTSDEL(ISTRA)  (input block 7)

c  IDENTICALLY MATCHING MULTIPROC. RUNS ON A SINGLE PROC:
c
c  only in runs with no correlated sampling, and with NINITL(ISTRA) ge 0:
c  simulate seeds of a multi-processor run, if run on a single processor.
c  Set a new seed after nptsdel particles to exactly match the seeds used in
c  a corresponding multi-processor run.

c  The multiprocessor run must have been set up such that it completed
c  exactly nptsdel(istra) trajectories on each of the iproc(istra) processors
c  which ran on stratum ISTRA.
c  The corresponding single processor run must complete exactly
c  nptsdel(istra)*iproc(istra) trajectories for stratum ISTRA

              IF (MOD(IPTSI-1,NPTSDEL(ISTRA)) == 0) THEN
                NINIST=NINITL(ISTRA)+IPTSI/NPTSDEL(ISTRA)*10000
c  initialize random number generator with a "legal" seed,
                idumran=ranset_eirene(ninist)

                INIV1=0
                INIV2=0
                INIV3=0
                INIV4=0
              END IF
            ENDIF
C...................................................................
C  NEXT MONTE CARLO HISTORY
c
c   LAUNCH A NEW PARTICLE NOW
c...................................................................

            XMCP(ISTRA)=XMCP(ISTRA)+1.
            NPANU=NPANU+1
            IPANU=IPANU+1
            ITRJ = NCHORI + MOD(IPANU,NTRJ) + 1
            NLEVEL=0
            CALL EIRENE_LOCAT1(IPANU)

C  IS BIRTH PROCESS SURVIVED?
            IF (.NOT.LGPART) GOTO 110
C
  102       CONTINUE
C  FOLLOW NEUTRAL PARTICLE
            IF (ITYP.EQ.0.OR.ITYP.EQ.1.OR.ITYP.EQ.2) THEN
              CALL EIRENE_FOLNEUT
C  FOLLOW TEST ION
            ELSEIF (ITYP.EQ.3) THEN
              CALL EIRENE_FOLION
            ENDIF
C  NEXT GENERATION ?
            IF (LGPART) GOTO 102
C
  110       CONTINUE

            IF (NLRAY(ISTRA)) THEN
              CALL EIRENE_CLEAR_TRAJECTORY (ITRJ)
            END IF

C  NUMBER OF REMAINING NODES AND NUMBER OF LEVELS AT NEXT NODE
            IF (NLEVEL.GT.0) THEN
  104         INODES=NODES(NLEVEL)-1
              NODES(NLEVEL)=INODES
              IF(INODES.LE.0) GO TO 103
C  RESTORE VARIABLES AND START NEW TRACK
              DO 105 J=1,NPARTC
                RPST(J)=RSPLST(J,NLEVEL)
  105         CONTINUE
              DO 106 J=1,MPARTC
                IPST(J)=ISPLST(J,NLEVEL)
  106         CONTINUE
              ITYP=ISPEZI(ISPZ,-1)
              IPHOT=ISPEZI(ISPZ,0)
              IATM=ISPEZI(ISPZ,1)
              IMOL=ISPEZI(ISPZ,2)
              IION=ISPEZI(ISPZ,3)
              IPLS=ISPEZI(ISPZ,4)
              CALL EIRENE_NCELLN(NCELL,NRCELL,NPCELL,NTCELL,NACELL,
     .         NBLOCK,NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
              NBLCKA=NSTRD*(NBLOCK-1)+NACELL
              NLSRFX=MRSURF.GT.0
              NLSRFY=MPSURF.GT.0
              NLSRFZ=MTSURF.GT.0
              NLSRFA=MASURF.GT.0
              IF (NLTRC) THEN
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP CRITICAL
#endif
                CALL EIRENE_CHCTRC(X0,Y0,Z0,0,12)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END CRITICAL
#endif
              ENDIF

!  PARTICLE TYPE AND SPECIES MAY HAVE CHANGED
!  PREPARE POINTER FOR UNIFIED SUBROUTINES FPATH, UPDATE, ETC.
              CALL EIRENE_SWITCH_PARTINFO

              IC_NEUT=0
              IC_ION=0
              GOTO 102
C  RETURN TO PREVIOUS LEVEL
  103         CONTINUE
              NLEVEL=NLEVEL-1
              IF(NLEVEL.GT.0) GOTO 104
            ENDIF
C  HISTORY HAS ENDED
C
C  IN CASE NLERG: EITHER LOGATM(1,ISTRA) OR LOGMOL(1,ISTRA)
C  ACTIVATE CORRESPONDING STANDARD DEVIATION ESTIMATOR
C
            IF (NLERG.AND.IPTSI.EQ.1) THEN
              IF (NMOLI >= 1) THEN
                IF (LOGMOL(1,ISTRA)) THEN
                  IIH(1)=2
                  CALL EIRENE_STATS0
                ENDIF
              ENDIF
            ENDIF

cdr  dec. 15:
cdr  update linear algebraic combinations of tallies after finishing trajectory
cdr  this enables also statistical variances for those tallies, avoiding covariance estimators.
            if (nmode.gt.0) call eirene_updlin  !cdr

C   MEAN SQUARE
            IF (NSIGI.GT.0) CALL EIRENE_STATS1
     .                                    (NSBOX_TAL,NR1TAL,NP2TAL,
     .                                     NT3TAL,NLIMPS,
     .                                     NLSYMP(ISTRA),NLSYMT(ISTRA))
            IF (NSIGI_SPC.GT.0) CALL EIRENE_STATS1_SPC
     .                                    (NSBOX_TAL,NR1TAL,NP2TAL,
     .                                     NT3TAL,NLIMPS,
     .                                     NLSYMP(ISTRA),NLSYMT(ISTRA))
C
            IF (TRCTIM) THEN
              SECND2=EIRENE_SECOND_OWN( )
              SECDEL=SECND2-SECND1
              CALL EIRENE_MASJ1R('PART., CPU TIME ',NPANU,SECDEL)
              if (secdel.gt.0.1)
     .          write (iunout,*) 'LONG HISTORY ',npanu,secdel
            ENDIF
  100     CONTINUE    !  nparts_loc(istra)

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END DO
#endif

          IF (.NOT.LGABORT) THEN
            CALL EIRENE_LEER(1)

            WRITE (iunout,*) 'ALL REQUESTED TRAJECTORIES COMPLETED'
            WRITE (iunout,'(1x,a,a)')
     .       'M.C. HISTORIES FOLLOWED UNTIL THAT TIME FOR ',
     .       'THIS STRATUM'
            CALL EIRENE_MASJ2('ISTRA,IPANU     ',ISTRA,IPANU)

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP MASTER
#endif
            call system_clock (itimend, itimrate)
            timused=real(itimend-itimstart,DP)/REAL(itimrate,DP)
            CALL EIRENE_MASR1('TIMUSED ',TIMUSED)

            IF (NPRNLI.GT.0) THEN
              WRITE (iunout,*) 'M.C. HISTORIES THAT SCORED AT CENSUS'
              CALL EIRENE_MASJ1 ('IPRNLS= ',IPRNLS)
            ENDIF
            IF (TRCLST) CALL EIRENE_OUTLST
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END MASTER
#endif
C           GOTO 101
          END IF
cym ccccccccccccccccccccccccccccccccccccccccccccc
cym see comment above
cym ccccccccccccccccccccccccccccccccccccccccccccc
#ifndef USE_OPENMP
C 101     CONTINUE
#endif
C
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP BARRIER
#endif
C
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP MASTER
#endif
C
          IF (KCFF.GT.0) THEN
            CALL EIRENE_LEER(1)
            WRITE(iunout,*) KCFF,' HISTORIES CUT OFF(',
     .        NEVNTMAX,' EVENTS)'
            DO I=1,3
              XVNMX2(I)=SQRT((XVNMX2(I)-XVNMX(I)**2./KCFF)/KCFF)
              XVNMX(I)=XVNMX(I)/KCFF
            ENDDO
            WRITE(iunout,'(t12,3(6x,a1,5x))') 'X','Y','Z'
            WRITE(iunout,'(a,1p,t12,3e12.4)') 'AVERAGE',XVNMX,
     .        'SPREAD',XVNMX2
            KVNMXTR(ISTRA)=KVNMXTR(ISTRA)+1
          END IF
C
          XMCT(istra)=timused
          call time_particles(istra, timused, ipanu)
          call time_particles_improved(istra, timused, ipanu)
csw
          SECND=EIRENE_SECOND_OWN()
C
C**** PARTICLE TRACING FOR THIS STRATUM FINISHED **********************
C
c
c    collect data for one stratum ISTRA from all PEs performing calculations
c    for this stratum
c
          if (need_calstr(istra)) then
            call system_clock (itimstart)
            call EIRENE_calstr
            call system_clock(itimend, itimrate)
            call time_calstr(istra, real(itimend-itimstart,DP)
     &                     /REAL(itimrate,DP))
            call time_calstr_improved(istra,
     &           real(itimend-itimstart,DP) / REAL(itimrate,DP))
          else
            call time_calstr(istra, 0.0_dp)
            call time_calstr_improved(istra, 0.0_dp)
          endif

          if (I_AM_LEADER(istra)) then
C
C  UPDATE AND CHECK LOGICALS FOR TALLIES
C
            call system_clock (itimstart)
            DO 120  JMOL=1,NMOLI
              LOGMOL(0,ISTRA)=LOGMOL(0,ISTRA).OR.LOGMOL(JMOL,ISTRA)
              LOGMOL(JMOL,0) =LOGMOL(JMOL,0) .OR.LOGMOL(JMOL,ISTRA)
  120       CONTINUE
            DO 130  JATM=1,NATMI
              LOGATM(JATM,0) =LOGATM(JATM,0) .OR.LOGATM(JATM,ISTRA)
              LOGATM(0,ISTRA)=LOGATM(0,ISTRA).OR.LOGATM(JATM,ISTRA)
  130       CONTINUE
            DO 133  JION=1,NIONI
              LOGION(JION,0) =LOGION(JION,0) .OR.LOGION(JION,ISTRA)
              LOGION(0,ISTRA)=LOGION(0,ISTRA).OR.LOGION(JION,ISTRA)
  133       CONTINUE
            DO 135  JPLS=1,NPLSI
              LOGPLS(JPLS,0) =LOGPLS(JPLS,0) .OR.LOGPLS(JPLS,ISTRA)
              LOGPLS(0,ISTRA)=LOGPLS(0,ISTRA).OR.LOGPLS(JPLS,ISTRA)
  135       CONTINUE
            DO JPHOT=1,NPHOTI
              LOGPHOT(JPHOT,0)=LOGPHOT(JPHOT,0).OR.LOGPHOT(JPHOT,ISTRA)
              LOGPHOT(0,ISTRA)=LOGPHOT(0,ISTRA).OR.LOGPHOT(JPHOT,ISTRA)
            END DO
            LOGMOL(0,0)=LOGMOL(0,0).OR.LOGMOL(0,ISTRA)
            LOGION(0,0)=LOGION(0,0).OR.LOGION(0,ISTRA)
            LOGATM(0,0)=LOGATM(0,0).OR.LOGATM(0,ISTRA)
            LOGPLS(0,0)=LOGPLS(0,0).OR.LOGPLS(0,ISTRA)
            LOGPHOT(0,0)=LOGPHOT(0,0).OR.LOGPHOT(0,ISTRA)
C
C  NUMBER OF LOCATED M.C. HISTORIES FOR THIS STRATUM: XMCP(ISTRA)
C
            IF(XMCP(ISTRA).LT.1.) GOTO 1111
C
            WTT=0._DP
            DO JPHOT=1,NPHOTI
              WTOTPH(0,ISTRA)=WTOTPH(0,ISTRA)+WTOTPH(JPHOT,ISTRA)
              WTT=WTT+WTOTPH(JPHOT,ISTRA)*NPRT(JPHOT)
            END DO
            DO 200 JATM=1,NATMI
              WTOTA(0,ISTRA)=WTOTA(0,ISTRA)+WTOTA(JATM,ISTRA)
              WTT=WTT+WTOTA(JATM,ISTRA)*NPRT(NSPH+JATM)
  200       CONTINUE
            DO 201 JMOL=1,NMOLI
              WTOTM(0,ISTRA)=WTOTM(0,ISTRA)+WTOTM(JMOL,ISTRA)
              WTT=WTT+WTOTM(JMOL,ISTRA)*NPRT(NSPA+JMOL)
  201       CONTINUE
            DO 202 JION=1,NIONI
              WTOTI(0,ISTRA)=WTOTI(0,ISTRA)+WTOTI(JION,ISTRA)
              WTT=WTT+WTOTI(JION,ISTRA)*NPRT(NSPAM+JION)
  202       CONTINUE
            WTOTE(ISTRA)=0._DP
            DO 203 JPLS=1,NPLSI
              WTOTP(0,ISTRA)=WTOTP(0,ISTRA)+WTOTP(JPLS,ISTRA)
              WTOTE(ISTRA)=WTOTE(ISTRA)+WTOTP(JPLS,ISTRA)*NCHRGP(JPLS)
              WTT=WTT-WTOTP(JPLS,ISTRA)*NPRT(NSPAMI+JPLS)
  203       CONTINUE
            CALL EIRENE_LEER(2)
            WRITE (iunout,*) 'TOTAL WEIGHT OF PRIMARY SOURCE PARTICLES'
            WRITE (iunout,*) 'BULK IONS, ATOMS, MOLECULES, TEST IONS'
            CALL EIRENE_MASR5
     .         ('WTPLS,WTATM,WTMOL,WTION,WTPHOT          ',
     .           WTOTP(0,ISTRA),WTOTA(0,ISTRA),WTOTM(0,ISTRA),
     .           WTOTI(0,ISTRA),WTOTPH(0,ISTRA))
            WRITE (iunout,*) 'TOTAL NUMBER OF MONTE CARLO HISTORIES'
            IMCP=INT(XMCP(ISTRA))
            CALL EIRENE_MASJ1 ('NPART   ',IMCP)
C
C
C  SET SOME SCALING CONSTANTS
C
            CALL EIRENE_SET_SCAL_CONST (ISTRA, WTT,
     .       ZWW, ZW, ZVOLNT, ZVOLWT, ZVOLIN, ZVOLIW, SCLTAL, N1MX)
C
C   STATISTICS, IF REQUESTED
C
            IF (XMCP(ISTRA).LE.1.) GOTO 219
C
C  FACTORS FOR STANDARD DEVIATION
            ZFLUX=FLXFAC(ISTRA)*XMCP(ISTRA)
            FSIG=SQRT(XMCP(ISTRA)/(XMCP(ISTRA)-1.))
C
            IF (NSIGI.GT.0) THEN
              CALL EIRENE_STATS2(XMCP(ISTRA),FSIG,ZFLUX)
C  CONVERT TO %
              SDVI1=MAX(0._DP,SDVI1-EPS6)*100.D0
              SDVI2=MAX(0._DP,SDVI2-EPS6)*100.D0
            ENDIF
            IF (NSIGI_SPC.GT.0) THEN
              CALL EIRENE_STATS2_SPC(XMCP(ISTRA),FSIG,ZFLUX)
C  CONVERT TO %
              DO ISPC=1,NADSPC
                ESTIML(ISPC)%SGMS=
     .           MAX(0._DP,ESTIML(ISPC)%SGMS-EPS6)*100.D0
                DO J=0,ESTIML(ISPC)%NSPC+1
                  ESTIML(ISPC)%SGM(J)=
     .             MAX(0._DP,ESTIML(ISPC)%SGM(J)-EPS6)*100.D0
                END DO
              END DO
            ENDIF
C
  219       CONTINUE

            CALL EIRENE_SCAL_VOLAV_TALLIES (ISTRA, ZWW, ZW,
     .           ZVOLIN, ZVOLIW, SCLTAL, N1MX)
C
C   REPLACE DEFAULT TALLIES BY USER-SUPPLIED
C   COLLISION-ESTIMATED TALLIES
C   THIS IS DONE BEFORE VOLUME INTEGRATION THUS THERE IS THE RISK TO
C   DESTROY TERMS NEEDED FOR GLOBAL BALANCES
C
            IF (LCOLV) THEN
              DO 285 ICLV=1,NCLVI
                IS=ICLVS(ICLV)
                IT=ICLVT(ICLV)
                IF (IT.LE.0.OR.IT.GE.NTALA) GOTO 285
                IGFFT=NFSTVI(IT)
                IF (IS.LE.0.OR.IS.GT.IGFFT) GOTO 285
                IADD=NADDV(IT)
                DO 286 J=1,NSBOX_TAL
                  ESTIMV(IADD+IS,J)=COLV(ICLV,J)
  286           CONTINUE
  285         CONTINUE
            END IF
C
C   REPLACE DEFAULT TALLIES BY USER-SUPPLIED
C   TRACKLENGTH-ESTIMATED TALLIES
C   THIS IS DONE BEFORE VOLUME INTEGRATION THUS THERE IS THE RISK TO
C   DESTROY TERMS NEEDED FOR GLOBAL BALANCES
C
            IF (LADDV) THEN
              DO 290 IADV=1,NADVI
                IS=IADVS(IADV)
                IT=IADVT(IADV)
                IF (IT.LE.0.OR.IT.GE.NTALA) GOTO 290
                IGFFT=NFSTVI(IT)
                IF (IS.LE.0.OR.IS.GT.IGFFT) GOTO 290
                IADD=NADDV(IT)
                DO 295 J=1,NSBOX_TAL
                  ESTIMV(IADD+IS,J)=ADDV(IADV,J)
  295           CONTINUE
  290         CONTINUE
            END IF
C
C
C   INTEGRATE VOLUME-AVERAGED PROFILES   450 --- 459
C
            CALL EIRENE_INTEGRATE_TALLIES (ISTRA)
C
C   SYMMETRISE VOLUME-AVERAGED TALLIES?
            IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
!pb        CALL EIRENE_SYMET(ESTIMV,NTALV,NRTAL,NR1TAL,NP2TAL,NT3TAL,
!pb     .             NLSYMP(ISTRA),NLSYMT(ISTRA))
              CALL EIRENE_SYMET(ESTIMV,NVOLTL,NRTAL,
     .             NR1TAL,NP2TAL,NT3TAL,
     .             NLSYMP(ISTRA),NLSYMT(ISTRA))
            ENDIF
C
C  WORK WITH VOLUME-AVERAGED TALLIES FOR THIS STRATUM FINISHED
C
C  SCALE SURFACE-AVERAGED ESTIMATORS AND OTHER FLUXES 600 - 630
C
            CALL EIRENE_SCAL_SURF_TALLIES (ISTRA,ZWW,ZW)
C
C
C   SUM OVER SURFACE INDEX
C   IN THE SURFACE-AVERAGED ESTIMATORS
C
C
C  SUM OVER SPECIES INDEX FOR INTEGRATED VOLUME-AVERAGED TALLIES
C                         AND INTEGRATED SURFACE-AVERAGED TALLIES
C
            CALL EIRENE_SUM_AVERAGE (ISTRA)
C
            CALL EIRENE_SCALE_TALLIES (ISTRA)
C
C   ALGEBRAIC EXPRESSION IN TALLIES 801--900
C
            IF (NALVI.GT.0.OR.NALSI.GT.0) THEN
C
              CALL EIRENE_ALGTAL
C
              IF (LALGV) THEN
                DO 830 IALV=1,NALVI
                  DUMMY(1:NSBOX_TAL) = ALGV(IALV,1:NSBOX_TAL)
                  CALL EIRENE_INTTAL (DUMMY,VOLTAL,1,1,
     .                 NSBOX_TAL,ALGVI(IALV,ISTRA),
     .                 NR1TAL,NP2TAL,NT3TAL,NBMLT)
                  ALGV(IALV,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)
  830           CONTINUE
              END IF
C
              IF (LALGS) THEN
                DO 832 IALS=1,NALSI
                  ALGSI(IALS,ISTRA)=0.
                  DO 831 J=1,NLIMPS
                    ALGSI(IALS,ISTRA)=ALGSI(IALS,ISTRA)+ALGS(IALS,J)
  831             CONTINUE
  832           CONTINUE
              END IF
C
            ENDIF
C
C  CALCULATE VOLUMETRIC LINE EMISSIVITIES FOR SELECTED SPECIES AND LINES
C
            IF (NLEMIS) THEN
              CALL EIRENE_EMISSIVITY (ISTRA,1,NUM_LINES,0)
            END IF
C
C  SCALE STANDARD DEVIATIONS, WHICH ARE NOT GIVEN IN % REL.ERROR
C  1/XMCP IS INCLUDED IN ZVOLIN,ZW,ZWW,... FOR TALLY AVERAGING
C  THEREFORE IT MUST BE MULTIPLIED HERE BECAUSE ONLY FLUX SCALING
C
            IF (XMCP(ISTRA).LE.1.D0) GOTO 950
C
            CALL EIRENE_SCALE_DEVIATION(ISTRA, ZWW, ZW, ZVOLNT, ZVOLWT,
     .                                  ZVOLIN, ZVOLIW, SCLTAL, N1MX)
C
  950       CONTINUE

            IESTR=ISTRA
C
C  CALL INTERFACE TO OTHER CODES TO RETURN DATA. STRATUM ISTRA
csw 13mar2013 ONLY WHEN RUN IN NON-PARALLEL MODE
C  OR WHEN RUN WITH EQUAL NUMBER OR MORE STRATA THAN PROCESSES
C
            IF (NMODE.GT.0) THEN
              CALL EIRENE_INFCOP_POST_STRATUM(ISTRA)
!tf Every stratum leader calls if3cop and the results are summed up
!tf by calling if3cop_sum after the strata loop. The SOLPS-ITER coupling
!tf is implemented this way. Coupling with other codes should either be
!tf modified accordingly, or the coupling routine should have some flag
!tf which tells whether to call if3cop from here (within strata loop),
!tf or after it.
!pb Therefore logical LIF3COP_FROM_LOOP has been added. It is set in
!pb EIRENE_DEFAULTS_USR. For the SOLPS-ITER coupling LIF3COP_FROM_LOOP=.T.
              IF (LIF3COP_FROM_LOOP) THEN
                IESTR=ISTRA
                ISTRAA=ISTRA
                ISTRAE=ISTRA
                IENTRY=0
                CALL EIRENE_IF3COP(IENTRY,LSTP,
     .                             IFIRST,ISTRAA,ISTRAE,NEW_ITER)
                NEW_ITER=1
              ENDIF
            ENDIF
C
C  WRITE RESULTS FOR THIS STRATUM ON TEMP. FILE
C

            IF (NFILEN.EQ.1) THEN
csw 18jul2011
csw 08mar2013 added check nprs < nstrai
cdr npesta is the master processor for stratum no ISTRA
              if (nprs==1.or.I_am_leader(istra)) then
                CALL EIRENE_WRSTRT(ISTRA,NSTRAI,
     .              NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .              TRCFLE)
              endif
            ENDIF
C
C  UPDATE TALLIES FOR "SUM OVER STRATA"
C
            IF (NSTRAI.EQ.1) GOTO 1111
C
            CALL EIRENE_SUMOSTRA (ISTRA) ! Integrals only
C
C
            IF (NSMSTRA == 1) THEN
C  VOLUME TALLIES
              do idv=1,nidv
                smestv(idv,1:nrtal) = smestv(idv,1:nrtal) +
     .                                estimv(idv,1:nrtal)
              end do
C  SURFACE TALLIES
              SMESTS = SMESTS + ESTIMS
C  SPECTRA
              DO ISPC=1,NADSPC
                SMESTL(ISPC)%SPC = SMESTL(ISPC)%SPC +
     .                             ESTIML(ISPC)%SPC
                SMESTL(ISPC)%SPCS = SMESTL(ISPC)%SPCS +
     .                              ESTIML(ISPC)%SPCS
                IF (ESTIML(ISPC)%ISPCOPT==2) THEN
                  SMESTL(ISPC)%SPCAN = SMESTL(ISPC)%SPCAN +
     .                   ESTIML(ISPC)%SPCAN
                END IF
              END DO
            END IF

C  covariances
            DO 1170 ISDV=1,NSIGCI
              DO 1172 ICELL=1,NSBOX_TAL
                STVC(0,ISDV,ICELL)=
     .           STVC(0,ISDV,ICELL)+SIGMAC(0,ISDV,ICELL)
                STVC(1,ISDV,ICELL)=
     .           STVC(1,ISDV,ICELL)+SIGMAC(1,ISDV,ICELL)**2
                STVC(2,ISDV,ICELL)=
     .           STVC(2,ISDV,ICELL)+SIGMAC(2,ISDV,ICELL)**2
 1172         CONTINUE
              STVCS(0,ISDV)=STVCS(0,ISDV)+SGMCS(0,ISDV)
              STVCS(1,ISDV)=STVCS(1,ISDV)+SGMCS(1,ISDV)**2
              STVCS(2,ISDV)=STVCS(2,ISDV)+SGMCS(2,ISDV)**2
 1170       CONTINUE
C
 1111       CONTINUE
            WRITE(iunout,*)
     .           'CUMULATED CPU TIME USED UNTIL END OF STRATUM ISTRA'
            WRITE(iunout,*) 'ISTRA, CPU(S) ',ISTRA,
     .           EIRENE_SECOND_OWN()-TIMST
            CALL EIRENE_LEER(2)
            call system_clock(itimend, itimrate)
            call time_postproc(istra, real(itimend-itimstart,DP)/
     &                         REAL(itimrate,DP))
            call time_postproc_improved(istra,
     &           real(itimend-itimstart,DP) / REAL(itimrate,DP))
          else
            IF (NMODE.GT.0) THEN
!  This routine has to be called by all processors
              CALL EIRENE_INFCOP_POST_STRATUM(ISTRA)
            END IF
          endif  ! I_AM_LEADER(ISTRA)
#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP END MASTER
#endif

        ELSE
         ! We update calstr timings even if we skipped the stratum
          call time_calstr(istra, 0.0_dp)
          call time_calstr_improved(istra, 0.0_dp)
        END IF                  ! CALC_STRATUM(ISTRA)
cdr  do another sub-stratum, for missing rare events stored on census
cdr     If...
cdr  not ready. currently we do this in next time-cycle

        call ioflush_usr
#if( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
        CALL EIRENE_CHECK_EXIT
!$OMP BARRIER
#endif

      END DO                    ! ISTR

#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP MASTER
#else
!pb 17.11.2021
#if ( defined(USE_OPENMP) && !defined(USE_EXT_OPENMP) )
!$OMP END PARALLEL
#endif
#endif
C
C*** STRATA LOOP FINISHED *******************************************
C
      call calstr_buffered_finish
      NPTS=NPTS_SAVE
      NINITL = NINITL_SAVE
C
      IF (NPRS > 1) THEN
cxpb  The (NPESTA/=0 .AND. NLSRON) mask breaks the output
cxpb        IF (ANY(NPESTA(1:NSTRAI) /= 0 .AND. NLSRON(1:NSTRAI)))
cxpb     &   CALL EIRENE_COLLECT_COUTAU
        CALL EIRENE_COLLECT_COUTAU
        IF (NPRNLI > 0) CALL EIRENE_COLLECT_CENSUS
      END IF

      IF (TRCHKTIM) CALL EIRENE_OUTPUT_PARTINFO

C
C  CALL INTERFACE TO OTHER CODES TO RETURN DATA. STRATUM ISTRA
C
csw 08mar2013 shifted behind STRATA LOOP, do all strata in one go
csw 13mar2013 do it here if in parallel mode
C   AND MORE PROCESSES THAN STRATA
      IF (NMODE.GT.0) THEN
        IF (LIF3COP_FROM_LOOP) THEN
!tf For SOLPS-ITER IF3COP is already called in the strata loop.
!tf we only need to sum up the results
          CALL EIRENE_IF3COP_SUM
        ELSE
          IF (NPRS > 1) THEN
C This is very case-specific and may be different for each plasma code.
C Introducing another interfacing subroutine within the strata loop
C solves this issue much more flexibly.
C This if block needs to go into the if3cop, if relevant for the
C plasma code.
            ISTRAA=1
            ISTRAE=NSTRAI
            CALL EIRENE_IF3COP(IENTRY,LSTP,
     .                         IFIRST,ISTRAA,ISTRAE,NEW_ITER)
            NEW_ITER=1
          ENDIF  ! nprs  > 1
        END IF
      ENDIF    ! nmode > 0
csw
      IF (NPRS > 1) THEN
        CALL EIRENE_COLLECT_DATA_USR
      END IF
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP END MASTER
#endif
      IF ((MY_PE .EQ. 0) .AND. (NSTRAI.EQ.1)) THEN
C
C  WRITE RESULTS FOR SUM OVER STRATA ON TEMP. FILE
C  USE THE DATA FOR STRATUM NO. 1 FOR THIS, RATHER THAN DOING
C  A USELESS SUMMATION
C
C  INDICATE: DATA FOR ISTRA=1 ARE ON CESTIM, BUT WRITE AS SUM OVER
C  STRATA
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP MASTER
#endif
        IESTR=1
        IF (NFILEN.EQ.1.OR.NFILEN.EQ.6) THEN
          CALL EIRENE_WRSTRT(0,NSTRAI,
     .              NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .              TRCFLE)
        ENDIF
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP END MASTER
#endif
        GOTO 2000
      ENDIF
      IF (XMCP(0).LE.1) GOTO 2000
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP MASTER
#endif
C SEQUENTIAL REGION

      IF(MY_PE .EQ. 0) THEN

C
C    STATISTICS, SUM OVER STRATA
C
      CALL EIRENE_STAT_SUMOSTRA
C
C  PUT SUM OVER STRATA BACK ONTO CESTIM, CSDVI, ESTIML...
C
      IF (NSMSTRA == 1) THEN
C  VOLUME-AVERAGED TALLIES
        ESTIMV(1:NIDV,1:NRTAL) = SMESTV(1:NIDV,1:NRTAL)
C  SURFACE-AVERAGED TALLIES
        ESTIMS = SMESTS
C  SPECTRA TALLIES
        DO ISPC=1,NADSPC
          ESTIML(ISPC)%SPC = SMESTL(ISPC)%SPC
          ESTIML(ISPC)%SPCS = SMESTL(ISPC)%SPCS
          IF (ESTIML(ISPC)%ISPCOPT==2) THEN
            ESTIML(ISPC)%SPCAN = SMESTL(ISPC)%SPCAN
          END IF
        END DO
C
C  NOW PUT VARIANCES FOR SUM OVER STRATA BACK ONTO VARIANCE TALLIES
C
C  SPECTRA TALLY VARIANCES
        DO ISPC=1,NADSPC
          IF (NSIGI_SPC > 0) THEN
            ESTIML(ISPC)%SGM = SMESTL(ISPC)%STV
            ESTIML(ISPC)%SGMS = SMESTL(ISPC)%STVS
          END IF
        END DO
C  CELL- AND SURFACE-AVERAGED DEFAULT TALLY VARIANCES
        SIGMA  = STV
        SGMS   = STVS
        SIGMAW = STVW
        SGMWS  = STVWS
C
C   ALGEBRAIC EXPRESSION IN TALLIES, SUM OVER STRATA  1571--1579
C
        IF (NALVI.GT.0.OR.NALSI.GT.0) THEN
C
          CALL EIRENE_ALGTAL
C
          DO 1571 IALV=1,NALVI
            DUMMY(1:NSBOX_TAL) = ALGV(IALV,1:NSBOX_TAL)
            CALL EIRENE_INTTAL (DUMMY,VOLTAL,1,1,
     .                   NSBOX_TAL,ALGVI(IALV,0),
     .                   NR1TAL,NP2TAL,NT3TAL,NBMLT)
            ALGV(IALV,1:NSBOX_TAL) = DUMMY(1:NSBOX_TAL)
 1571     CONTINUE
C
          DO 1572 IALS=1,NALSI
            ALGSI(IALS,0)=0.
            DO 1573 J=1,NLIMPS
              ALGSI(IALS,0)=ALGSI(IALS,0)+ALGS(IALS,J)
 1573       CONTINUE
 1572     CONTINUE
C
        ENDIF

!     Update pressure feedback control loop:
      DO I=1, NLIMPS
        CALL updatePressureFeedback(RPRESSFED(I))
      END DO

C
C  CALCULATE VOLUMETRIC LINE EMISSIVITIES, SUM OVER STRATA
C
        IF (NLEMIS) THEN
          CALL EIRENE_EMISSIVITY (0,1,NUM_LINES,0)
        ENDIF
C
C  WRITE RESULTS FOR SUM OVER STRATA ON TEMP. FILE
C
cdr IESTR = 0 tells the code that all tallies now contain
cdr           "sum over strata"
        IESTR=0

        IF (NFILEN.EQ.1.OR.NFILEN.EQ.6) THEN
cdr spectrum tally variances are already in ESTIML
          CALL EIRENE_WRSTRT(0,NSTRAI,
     .                NESTM1,NESTM2,NADSPC,
     .                ESTIMV,ESTIMS,ESTIML,
     .                NSDVI1,SDVI1,NSDVI2,SDVI2,
     .                NSDVC1,SIGMAC,NSDVC2,SGMCS,NSIGI_SPC,
     .                TRCFLE)
        ENDIF
      ENDIF ! NSMSTRA == 1
C
      ENDIF ! MY_PE .EQ. 0
C
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP END MASTER
#endif
 2000 CONTINUE
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP MASTER
#endif
      CALL EIRENE_BROAD_IESTR(IESTR)

      IF(MY_PE .EQ. 0) THEN
C
C  SAVE OR RESTORE SOME DATA FOR "EIRENE RECALL OPTION NFILE.NE.0"
C  FROM FILE "FT11" (all data in module COUTAU)
C  NOTE: RECORD IRC=3 MAY BE USED IN INTERFACING ROUTINE INFCOP
C
      IF (NFILEN.EQ.1.OR.NFILEN.EQ.6) THEN
        IF (TRCFLE) WRITE (iunout,*) 'WRITE DATA FOR RECALL OPTION '
        IRC=1
        WRITE (11+ifoff,REC=IRC) LOGATM,LOGION,LOGMOL,LOGPLS,LOGPHOT
#ifdef CHECKBIN
        WRITE (111+ifoff,*) 'LOGATM,LOGION,LOGMOL,LOGPLS,LOGPHOT ',
     .                 LOGATM,LOGION,LOGMOL,LOGPLS,LOGPHOT
#endif
        IF (TRCFLE) WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL EIRENE_WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11+ifoff,REC=IRC) OUTAU
#ifdef CHECKBIN
        WRITE (111+ifoff,*) 'OUTAU '
        do ipb=1,noutau,5
          WRITE (111+ifoff,*) ipb,(OUTAU(jpb),jpb=ipb,min(ipb+4,noutau))
        end do
#endif
        DEALLOCATE (OUTAU)
        IF (TRCFLE) WRITE (iunout,*) 'WRITE 11  IRC= ',IRC

      ELSEIF (NFILEN.EQ.2.OR.NFILEN.EQ.7) THEN
cdr  in this case the entire MC calculation has been skipped ("recall option" only)
cdr  Nothing has been recalculated in present cycle.
cdr  Both Monte Carlo loops:
cdr     DO ISTRA=1,NSTRAI              (strata)
cdr       DO 100 IPTSI=1,NPTS(ISTRA)   (histories)
cdr  are bypassed.
        IF (TRCFLE) WRITE (iunout,*) 'READ DATA FOR RECALL OPTION'
        IRC=1
        READ (11+ifoff,REC=IRC) LOGATM,LOGION,LOGMOL,LOGPLS,LOGPHOT
        IF (TRCFLE) WRITE (iunout,*) 'READ 11  IRC= ',IRC
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        READ (11+ifoff,REC=IRC) OUTAU
        CALL EIRENE_READ_COUTAU (OUTAU, IUNOUT)
        DEALLOCATE (OUTAU)
        IF (TRCFLE) WRITE (iunout,*) 'READ 11  IRC= ',IRC
      ENDIF

C END SEQUENTIAL REGION
      ENDIF

cdr  dec. 15
cdr  see above. Routine UPDLIN.f contains linear combination of tallies
      if (nmode.gt.0) call eirene_reset_updlin

      CALL EIRENE_CHECK_EXIT
#if ( defined(USE_EXT_OPENMP) && !defined(USE_OPENMP) )
!$OMP END MASTER
#endif
C
      RETURN
      END SUBROUTINE EIRENE_MCARLO

      SUBROUTINE EIRENE_MCARLO2
      IMPLICIT NONE

      IF (ALLOCATED(DUMMY)) THEN
         DEALLOCATE (DUMMY,ZVOLIN,ZVOLIW,SCLTAL)
      END IF

      RETURN
      END SUBROUTINE EIRENE_MCARLO2

      END MODULE EIRMOD_MCARLO
