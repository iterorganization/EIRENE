cdr   mar   20: remove routines LAX, GALPD (linear alg. eqs. solver)
cdr             and put them to "mathematics", for common use by all
cdr             internal CR models.
cdr   mar   19: change sequence in parameter list (pop_esc),
cdr             to better sync with logic in He_Colrad.
cdr   feb   18: sync with h_colrad, H,H2 CRM, Sawada-Fujimoto-Reiter
cdr             cleaned up, more comments.
cdr             added lopaque, ebeta, e_alpcr
cdr      tbd:   avoid 3rd right  hand side in popcof if Q_ext==0: done. l_ext
cdr      tbd:   better do: gauss laguerre integration of gaunt3 and gaunt4 ?
cdr      tbd:   avoid execution of develop-tests: e_alpcr_T, e-scr_T,....: done. ctt
cdr   Jan   18: parameter ICELL removed, icell is now controlled by
cdr             new calling CRM-driver routine COLRAD
cdr             popcof and matrix solver LAX, GALPD extended to handle
c               up to three right hand sides (parent states) simultaneously
c               rather than inverting the matrix three times.
cdr   July  17: bug fix in function mmdei (exp. integr.)
cdr             A typo during synchronisation with solps-iter.
cdr             correct: z=0.25*y, rather than z=0.25+0*y
cdr:  April 17: synchronized with version from solps-iter: spelling errors in comments,
cdr             use EIRMOD_PRECISION instead of real*8
cdr             (this may complicate standalone use, outside eirene)
cdr             call "exit_own" rather than "stop", further cleanup...
cdr             remaining differences:
cdr                 use eirmod_ccrm (Vlad Kotov) in solps-iter (commented out)
cdr                 lima=34 or lima=40, lima undefined in solps-iter?
CDR  INSIDE EIRENE: THIS H_COLRAD IS CALLED FROM EIRENE_COLRAD.
CDR  done:       there: introducedan array 'visited(icell)' and store e-rate, etc..., further possible data
cdr              for next call to H_colrad.
cdr  EIRENE_COLRAD is called from rate_coef.f, from energy_rate_coef.f,
cdr                 and from other_rate_coef.f
cdr           to provide ionization, radiation and electron cooling rates, 
cdr           for given Te, ne., as well as to provide reduced CR population coefficients
c****************************************************************************************************
C*
C*     COLLISIONAL-RADIATIVE MODEL OF
C*
C*     ATOMIC HYDROGEN
C*
C*
C   ASSUME: SLOWLY EVOLVING SPECIES: H,H+
C   ASSUME: QUASI STEADY STATE OF H*(N) WITH H, H+
C
C   INPUT:

C   TEMP      : ELECTRON TEMPERATURE
C   DENSEL    : ELECTRON DENSITY
C
C   L_EXT     : 3RD (EXTERNAL) SOURCE OF EXCITED STATES (E.G. PHOTO-EXCITATION)
C   L_EXT, Q_EXT(N): ??? ->  H*(N) external source rate, e.g. molecules,
C                                  or photo-excitation

C   LOPAQUE: make ALL Lyman transitions entirely black (older option)
C   POP_ESC: population escape factor matrix for individual transitions

C
C
C   OUTPUT:
C   R0(..)    : TRAIN OF H* TRAVELING WITH H+
C   R1(..)    : TRAIN OF H* TRAVELING WITH H
C
C   R_EXT(..) : TRAIN OF H* DUE TO EXTERNAL SOURCE FOR H*(N)
C
c   C: elec. impact excitation processes
c   F: elec. impact de-excitation processes (inverse to C, detailed balance)
c   A: spontaneous radiative decay
c   S: ionization
c
c   ALPHA(N): H+    ->  H*(N)  three-body recombination from H+
c                              (inverse to S: elect. impact ionization)
c   BETA(N) : H+    ->  H*(N)  radiative rec. from H+
C   EBETA      :                  CORRESPONDING ENERGY-WEIGHTED RATE
c   C(1,N)  : H(1)  ->  H*(N)  excitation from ground state
C   Q_EXT(N): ???   ->  H*(N)  external source


c   reduced pop coeff r0,r1,r_ext are in [units-S]/[units-M].
c   eff rates scr, alp,  are in [units-S]
c  Here new: 
c   units matrix M:  rate  (Aik) 1/s) or ne * rate-coeff (also: 1/s)
c   units source S:  S, beta, ne*alpha: rate coef (cm**3/s)
c   hence: for pop0,pop1 taken times "densel", to make them dimensionless
c   
c   pop_ext (=pop2) - arrays of reduced population coefficients....depends (linearly) on units of Q_ext
c   scr_ext, alp_ext. same units as Q_ext.
C*
C***********************************************************************
      SUBROUTINE EIRENE_H_COLRAD (TEMP, DENSEL, IFORM,                        ! in
     .                            POP_ESC,                                    ! in,  selected pop esc factors for some lines
     .                            Q_EXT, L_EXT,                               ! in,  rate coef. cm**3/s
     .                            AIKEINST, ENERGLEV,                         ! out, in first call
     .                            STATWGHT,                                   ! out, in first call
     .                            POP0, POP1, POP_EXT,                        ! out
     .                            ALPCR, SCR, SCR_EXT, ALP_EXT,               ! out
     .                            E_ALPCR, E_SCR, E_SCR_EXT, E_ALP_EXT)       ! out
ctt  .                           ,E_ALPCR_T, E_SCR_T, E_SCR_EXT_T  ! for testing

      USE EIRMOD_PRECISION
C     USE EIRMOD_CCRM
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE

C--------- ATOMIC PARAMETER ------------------------------------------
      INTEGER, INTENT(IN) :: IFORM
      REAL(DP), INTENT(IN) :: TEMP, DENSEL
      REAL(DP), INTENT(OUT) :: AIKEINST(40,40), ENERGLEV(40),
     .                         STATWGHT(40)
      REAL(DP), INTENT(IN) :: Q_EXT(40), POP_ESC(40,40)
      logical lopaque,l_ext

      REAL(DP), INTENT(OUT) ::   ALPCR,    SCR,     SCR_EXT, ALP_EXT
      REAL(DP), INTENT(OUT) :: E_ALPCR,  E_SCR,   E_SCR_EXT, E_ALP_EXT
ctt   REAL(DP), INTENT(OUT) :: E_ALPCR_T,E_SCR_T, E_SCR_EXT_T
      REAL(DP), INTENT(OUT) :: POP0(40), POP1(40), POP_EXT(40)

      REAL(DP), SAVE :: A(40,40), E_AT(40), OSC(40,40)
c     REAL(DP), SAVE :: POP_ESC(40,40)

      REAL(DP) :: C(40,40),S(40),F(40,40)
     &           ,SAHA(40),BETA(40),ALPHA(40),EBETA(40)
      REAL(DP) :: R1(40),R0(40),R_EXT(40)
C
cdr   integer, parameter :: lupa=34, lima=40
      integer, parameter :: lupa=34, lima=34
      INTEGER, SAVE :: IFRST=0
      INTEGER :: IP
      REAL(DP) :: P
      EXTERNAL :: EIRENE_CLSAHA, EIRENE_EINSTN,
     .            EIRENE_E_IONREC, EIRENE_IONREC,
     .            EIRENE_POPCOF_M, EIRENE_RATCOF,
     .            EIRENE_EXIT_OWN
c
cdr atomic structure: to be done only once,not for each cell
      if (ifrst == 0) then
        aikeinst=0.0_dp
        energlev=0.0_dp
        statwght=0.0_dp
      endif

      IF (LIMA.GT.40.OR.LUPA.GT.LIMA) THEN
        WRITE (iunout,*) 'LIMA, LUPA ??? ',LIMA,LUPA
        CALL EIRENE_EXIT_OWN(1)
      ENDIF
C
C ATOM
c  pop_esc   : population escape factor
c  pop_esc= 1: opt. thin
c  pop_esc= 0: opt. thick
c
c  only once and for all for a given CR model dataset on REACDAT
c
      lopaque=.false.  ! .true.: fully Lyman opaque, all transitions to ground state blocked


      IF (IFRST == 0) THEN  ! must be redone, if lopaque or pop_esc change
cdr  better: move lopaque, pop_esc outside this routine. And check always
        CALL EIRENE_EINSTN(OSC,A,E_AT,40,lopaque,POP_ESC)
        aikeinst=A
        energlev=E_AT
        IFRST = 1
      END IF
c
c
C ATOM
      CALL EIRENE_CLSAHA(TEMP,SAHA)
      do ip=1,40
        P=ip
        statwght(ip)=P*P      
      enddo

      CALL EIRENE_RATCOF(TEMP,OSC,SAHA,C,F,S,ALPHA,BETA,EBETA,lopaque)

C
C***********************************************************************
C
C ATOMIC HYDROGEN
C

      CALL EIRENE_POPCOF_M(DENSEL,SAHA,C,F,S,A,ALPHA,BETA,LUPA,LIMA,
     &             R0,R1,R_EXT,
     &                   Q_EXT,L_EXT)
! ebeta is not needed here

C TRAINS OF ELECTRONICALLY EXCITED H
      DO IP=2,LIMA
C
CDR COUPLING TO H+ IONS
        POP0(IP)=R0(IP)*DENSEL
CDR COUPLING TO H ATOMS
        POP1(IP)=R1(IP)*DENSEL
CDR COUPLING TO EXTERNAL SOURCE Q_ext  FOR H*(N)
cdr  q_ext:  source per second:  1/s:  no density factor densel, since CR matrix: 1/s. also R_ext [1]
cdr  aber Vorsicht: dann ist auch src_ext und alp_ext in 1/s, nicht in cm**3/s
        POP_EXT(IP)=R_EXT(IP)

      END DO
C  EFFECTIVE COLLISION RATE COEFFICIENTS, ATOMS
      CALL EIRENE_IONREC(C,S,SAHA,A,ALPHA,BETA,
     &            R0,R1,DENSEL,LUPA,LIMA,
     &            F,R_EXT,Q_EXT,L_EXT,
     &            ALPCR,SCR,
     &            ALP_EXT,SCR_EXT)
C***********************************************************************
C  EFFECTIVE ELECTRON COOLING RATE COEFFICIENTS, ATOMS
      CALL EIRENE_E_IONREC(C,S,SAHA,A,ALPHA,BETA,EBETA,
     &              R0,R1,DENSEL,LUPA,LIMA,
     &              F,R_EXT,Q_EXT,L_EXT,E_AT,
     &              ALPCR,      SCR,    SCR_EXT,
     &              E_ALPCR,  E_SCR,  E_SCR_EXT, E_ALP_EXT
ctt  &             ,E_ALPCR_T,E_SCR_T,E_SCR_EXT_T
     &              )
C***********************************************************************
C

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      RETURN
      END SUBROUTINE EIRENE_H_COLRAD

C**********************************************************************
      SUBROUTINE EIRENE_EINSTN(F,A,E_AT,LIM,LOPAQUE,POP_ESC)
C
C     CALCULATION OF OSCILLATOR STRENGTH AND EINSTEIN COEFFICIENT
C     FOR ATOMIC HYDROGEN
C     E_AT are the energy levels. Ground state is E_AT(1)=0.
C
C     L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972).
C
      USE EIRMOD_PRECISION
      USE EIRMOD_CCONA, ONLY : PIA
      IMPLICIT NONE
      INTEGER LIM
      REAL(DP) F(40,40)
      REAL(DP) A(40,40)
      REAL(DP) E_AT(40)
      REAL(DP) POP_ESC(40,40)
      INTEGER I, J
      REAL(DP) AI, AJ, G, P, UH, X
      logical lopaque

      UH=13.595_DP
      DO 100 I=1,LIM
        P=I
        E_AT(I)=UH*(1.0_DP-1.0_DP/P**2)
  100 CONTINUE

      DO 101 I=1,LIM-1
       DO 102 J=I+1,LIM
        AI=REAL(I,DP)
        AJ=REAL(J,DP)
        X=1.0_DP-(AI/AJ)**2

cdr: Johnson Gaunt factor approx.
cdr  gaunt=g(i,x)=G(I,J)
        IF(I.GE.3) THEN

          G=(.9935_DP+.2328_DP/AI-.1296_DP/AI**2)
     *     -(.6282_DP-.5598_DP/AI+.5299_DP/AI**2)/(AI*X)
C Feb. 2008: TYPO CORRECTED: .3387 --> .3887
     *     +(.3887_DP-1.181_DP/AI+1.470_DP/AI**2)/(AI*X)**2
        ELSEIF(I.EQ.2) THEN
          G=1.0785_DP-.2319_DP/X+.02947_DP/X**2
        ELSEIF(I.EQ.1) THEN
          G=1.1330_DP-.4059_DP/X+.07014_DP/X**2
        ENDIF

cdr  gaunt(x) done. gaunt=g(i,j)

        F(I,J)=2.0_DP**6/(3.0_DP*SQRT(3.0_DP)*PIA)*(AI/AJ)**3/
     .        (2.0_DP*AI**2)*G/X**3
        A(J,I)=8.03E9_DP*AI**2/AJ**2*(AI**(-2)-AJ**(-2))**2*F(I,J)

cdr  Fully Ly opaque: all a(j-->1) transitions are removed.
        if (lopaque.and.i.eq.1) a(j,i)=0.0_DP
cdr apply population escape factor to transition J-->I
        A(J,I)=A(J,I) * POP_ESC(J,I)
  102  CONTINUE
  101 CONTINUE

      RETURN
      END SUBROUTINE EIRENE_EINSTN

C***********************************************************************
      SUBROUTINE EIRENE_CLSAHA(TEMP,SAHA)
C
C     SAHA-BOLTZMANN COEFFICIENT FOR ATOMIC HYDROGEN
C
C
C
      USE EIRMOD_PRECISION
      USE EIRMOD_CCONA, ONLY : EVKEL
      IMPLICIT NONE
      REAL(DP) TEMP,SAHA(40)
      INTEGER I
      REAL(DP) P, TE, UION

      TE=TEMP/EVKEL

      DO 101 I=1,40
        P=REAL(I,DP)
        UION=13.595_DP/TEMP/P**2
c     if (uion.gt.65_DP) then
c       saha(i)=1.e30_DP
c     else
        SAHA(I)=P**2*EXP(UION)/2.414E15_DP/SQRT(TE**3)
c     endif
  101 CONTINUE

      RETURN
      END SUBROUTINE EIRENE_CLSAHA

C***********************************************************************
      SUBROUTINE EIRENE_RATCOF(TEMP,OSC,SAHA,C,F,S,ALPHA,BETA,EBETA
     .                  ,lopaque)
C
C     RATE COEFFICIENTS FOR ATOMIC HYDROGEN
CDR   ALSO: ERATE, ENERGY-WEIGHTED RATE FOR RAD. RECOMB.
C
C
C
      USE EIRMOD_PRECISION
      USE EIRMOD_CCONA, ONLY : EVKEL
      IMPLICIT NONE
      REAL(DP) TEMP
      REAL(DP) OSC(40,40),C(40,40),F(40,40),U(40,40)
      REAL(DP) SAHA(40),S(40),ALPHA(40),BETA(40),EBETA(40),UION(40)
      INTEGER I, II, III, J, JJ
      REAL(DP) DCCX, ELN, EP, EXS, P, TE, TEL10, UH, XP, XS
      EXTERNAL EIRENE_CLBETA, EIRENE_EXCOFF

cdr  above some critical Te0 value the radiative rate coefficients beta become unphysical
c    perhaps due to numerical integration, or due to fit expression for integrand.
c    Extrapolate energy rate coeff ebeta beyond Te0 =500 eV

      logical lopaque
      REAL(DP) RECMB1(9)
      DATA RECMB1/
     .-2.955397218588D+01,-5.415952519274D-01, -1.898544362547D-02,
     .-6.296638424302D-03,-6.590547592879D-04,  1.979096752346D-04,
     .-3.384091511722D-06,-1.312178532481D-06,  6.794867069246D-08/
      REAL(DP) RECMB2(9)
      DATA RECMB2/
     .-3.028258203128D+01,-6.489344359811D-01,-5.044516146247D-02,
     .-6.834257068606D-03, 7.621603311294D-04, 1.820369471925D-04,
     .-4.172173992553D-05, 3.151129832749D-06,-8.557790319804D-08/
      REAL(DP) RECMB3(9)
      DATA RECMB3/
     .-3.081308111854D+01,-7.506832188514D-01,-6.392987382349D-02,
     .-3.136753737114D-03, 1.130486114398D-03,-5.891101787416D-06,
     .-1.871544818151D-05, 1.924155011543D-06,-6.113117867247D-08/
      REAL(DP) RECMB4(9)
      DATA RECMB4/
     .-3.124178248661D+01,-8.324721778000D-01,-6.749686238962D-02,
     .-1.819027715538D-04, 1.029662235699D-03,-9.947563495630D-05,
     .-1.157212109508D-06, 6.419046650810D-07,-2.632212378124D-08/
      REAL(DP) RECMB5(9)
      DATA RECMB5/
     .-3.160248102907D+01,-8.970582941706D-01,-6.672021078686D-02,
     . 1.735998035356D-03, 8.102846437404D-04,-1.293548645460D-04,
     . 7.653597399125D-06,-1.007290294963D-07,-4.396689144443D-09/
      REAL(DP) RECMB6(9)
      DATA RECMB6/
     .-3.191200354649D+01,-9.487943192650D-01,-6.416290359056D-02,
     . 2.915417324842D-03, 5.904290511293D-04,-1.293159094005D-04,
     . 1.121762622770D-05,-4.647730339162D-07, 7.248953722425D-09/
      REAL(DP) RECMB7(9)
      DATA RECMB7/
     .-3.218096339788D+01,-9.910879396074D-01,-6.095219508517D-02,
     . 3.617704413680D-03, 4.023923659235D-04,-1.170991503299D-04,
     . 1.208395104539D-05,-6.153747484991D-07, 1.276002457208D-08/
      REAL(DP) RECMB8(9)
      DATA RECMB8/
     .-3.241655828927D+01,-1.026367760090D+00,-5.758688104786D-02,
     . 4.015706361086D-03, 2.503898666107D-04,-1.009126065615D-04,
     . 1.164859556536D-05,-6.527211813278D-07, 1.485788926575D-08/
      REAL(DP) RECMB9(9)
      DATA RECMB9/
     .-3.262387280670D+01,-1.056350328756D+00, -5.428573116785D-02,
     . 4.218200737528D-03, 1.304699903027D-04, -8.440144760676D-05,
     . 1.063623086507D-05,-6.324515306256D-07,  1.511971184207D-08/
      REAL(DP) RECMB10(9)
      DATA RECMB10/
     .-3.280658322283D+01,-1.082264892917D+00,-5.113789272221D-02,
     . 4.293438883183D-03,3.692139459358D-05 ,-6.905886517440D-05,
     . 9.413397820660D-06,-5.846761030723D-07, 1.443147593158D-08/
      REAL(DP) RECMB11(9)
      DATA RECMB11/
     .-3.296736752352D+01, -1.105008301698D+00, -4.817218718568D-02,
     . 4.284573570871D-03, -3.554562646018D-05, -5.539693364836D-05,
     . 8.160180100312D-06, -5.253448069028D-07,  1.327912727320D-08/
      REAL(DP) RECMB12(9)
      DATA RECMB12/
     .-3.310815524446D+01,-1.125249327580D+00,-4.538927682514D-02,
     . 4.219407043180D-03,-9.129558854381D-05,-4.350496815893D-05,
     . 6.963937973383D-06,-4.630004759641D-07, 1.193307256583D-08/
      REAL(DP) RECMB13(9)
      DATA RECMB13/
     .-3.323028186699D+01,-1.143497947426D+00,-4.277662562220D-02,
     . 4.116029910827D-03,-1.338072127899D-04,-3.328224832213D-05,
     . 5.863111240449D-06,-4.021451443280D-07, 1.054515552989D-08/
      REAL(DP) RECMB14(9)
      DATA RECMB14/
     .-3.333457692604D+01,-1.160151559549D+00,-4.031606772897D-02,
     . 3.986144192297D-03,-1.657642895871D-04,-2.454533247143D-05,
     . 4.866420567932D-06,-3.445615518446D-07, 9.181805185443D-09/
      REAL(DP) RECMB15(9)
      DATA RECMB15/
     .-3.342140353146D+01, -1.175528330345D+00, -3.798784446956D-02,
     . 3.837630497852D-03, -1.892109056864D-04, -1.713312193944D-05,
     . 3.978313590670D-06, -2.915668119657D-07,  7.895773130286D-09/
      REAL(DP) RECMB16(9)
      DATA RECMB16/
     .-3.349065709684D+01,-1.189890025875D+00,-3.577113355659D-02,
     . 3.675431118626D-03,-2.057472771812D-04,-1.086107122550D-05,
     . 3.190471670743D-06,-2.432219763675D-07, 6.698796511582D-09/
      REAL(DP) RECMB17(9)
      DATA RECMB17/
     .-3.354171350627D+01, -1.203459634476D+00,-3.364574539433D-02,
     . 3.502705068758D-03, -2.165263655319D-04,-5.593693576382D-06,
     . 2.496628300465D-06,-1.995492586410D-07 , 5.598476502405D-09/
      REAL(DP) RECMB18(9)
      DATA RECMB18/
     .-3.357331947394D+01,-1.216435628777D+00, -3.159069352529D-02,
     . 3.321134001198D-03,-2.224762687832D-04, -1.195846543818D-06,
     . 1.887944138620D-06,-1.603648021903D-07, 4.597768509898D-09 /
      REAL(DP) RECMB19(9)
      DATA RECMB19/
     .-3.358338560247D+01,-1.229004753335D+00,-2.958476231981D-02,
     . 3.131349929429D-03,-2.242129098657D-04, 2.426519089379D-06,
     . 1.357374516638D-06,-1.253982655621D-07, 3.692268883687D-09/
      REAL(DP) RECMB20(9)
      DATA RECMB20/
     .-3.356862240333D+01, -1.241355838486D+00, -2.760425715542D-02,
     . 2.932832208231D-03, -2.221882606505D-04,  5.370072522835D-06,
     . 8.951422982914D-07, -9.414351171950D-08,  2.871350578727D-09/
      REAL(DP) RECMB21(9)
      DATA RECMB21/
     .-3.352387025175D+01, -1.253698383107D+00, -2.562222079503D-02,
     . 2.724279920091D-03, -2.166083216770D-04, 7.671932594855D-06 ,
     . 4.995725943925D-07, -6.664545124985D-08,  2.139078683720D-09/
      REAL(DP) RECMB22(9)
      DATA RECMB22/
     .-3.344080404025D+01,-1.266290224370D+00, -2.360370763255D-02,
     . 2.502743536910D-03,-2.074812094740D-04,  9.376815010545D-06,
     . 1.640662776104D-07,-4.245674204110D-08,  1.483365500922D-09/
      REAL(DP) RECMB23(9)
      DATA RECMB23/
     .-3.330515902115D+01,-1.279492872276D+00, -2.149943684787D-02,
     . 2.263297751642D-03,-1.945289210035D-04,  1.048130070914D-05,
     .-1.111440892860D-07,-2.168045908152D-08,  9.097953739845D-10/
      REAL(DP) RECMB24(9)
      DATA RECMB24/
     .-3.308989756561D+01,-1.293895127861D+00,-1.922974855128D-02,
     . 1.996974143238D-03,-1.769459921577D-04, 1.093203266053D-05,
     .-3.242765500448D-07,-4.272006247928D-09, 4.140158369421D-10/
      REAL(DP) RECMB25(9)
      DATA RECMB25/
     .-3.273412834255D+01,-1.310670020686D+00,-1.663863373590D-02,
     . 1.685306322247D-03,-1.528534971941D-04, 1.057088450235D-05,
     .-4.649188077982D-07, 9.294915640077D-09, 6.752259753282D-12/
      REAL(DP) RECMB26(9)
      DATA RECMB26/
     .-3.204416240346D+01,-1.333113656592D+00,-1.330417762787D-02,
     . 1.278449892734D-03,-1.168945644350D-04, 8.891790020206D-06,
     .-4.953005276433D-07, 1.718610439449D-08,-2.710780782107D-10/
      REAL(DP) RECOMB(9,26)
      EQUIVALENCE (RECOMB(1,1),RECMB1(1))
      EQUIVALENCE (RECOMB(1,2),RECMB2(1))
      EQUIVALENCE (RECOMB(1,3),RECMB3(1))
      EQUIVALENCE (RECOMB(1,4),RECMB4(1))
      EQUIVALENCE (RECOMB(1,5),RECMB5(1))
      EQUIVALENCE (RECOMB(1,6),RECMB6(1))
      EQUIVALENCE (RECOMB(1,7),RECMB7(1))
      EQUIVALENCE (RECOMB(1,8),RECMB8(1))
      EQUIVALENCE (RECOMB(1,9),RECMB9(1))
      EQUIVALENCE (RECOMB(1,10),RECMB10(1))
      EQUIVALENCE (RECOMB(1,11),RECMB11(1))
      EQUIVALENCE (RECOMB(1,12),RECMB12(1))
      EQUIVALENCE (RECOMB(1,13),RECMB13(1))
      EQUIVALENCE (RECOMB(1,14),RECMB14(1))
      EQUIVALENCE (RECOMB(1,15),RECMB15(1))
      EQUIVALENCE (RECOMB(1,16),RECMB16(1))
      EQUIVALENCE (RECOMB(1,17),RECMB17(1))
      EQUIVALENCE (RECOMB(1,18),RECMB18(1))
      EQUIVALENCE (RECOMB(1,19),RECMB19(1))
      EQUIVALENCE (RECOMB(1,20),RECMB20(1))
      EQUIVALENCE (RECOMB(1,21),RECMB21(1))
      EQUIVALENCE (RECOMB(1,22),RECMB22(1))
      EQUIVALENCE (RECOMB(1,23),RECMB23(1))
      EQUIVALENCE (RECOMB(1,24),RECMB24(1))
      EQUIVALENCE (RECOMB(1,25),RECMB25(1))
      EQUIVALENCE (RECOMB(1,26),RECMB26(1))

C     INITIALIZATION
      S=0.0_DP
      ALPHA=0.0_DP
      BETA=0.0_DP
      EBETA=0.0_DP
      UION=0.0_DP
      C=0.0_DP
      F=0.0_DP
      U=0.0_DP

      TE=TEMP/EVKEL    !  Te in Kelvin, TEMP in eV
      TEL10=log10(temp)
      UH=13.595_DP

      DO 101 I=1,40
        P=REAL(I,DP)
        UION(I)=13.595_DP/TEMP/P**2
  101 CONTINUE
      DO 102 I=1,40
       DO J=1,40
        U(I,J)=UION(I)-UION(J)
       END DO
  102 CONTINUE

c  excitation   :            C
c  de-excitation:            F  (by detailed balance)
c  ionization:               S
c  three body recombination: ALPHA (inverse to ionization S)
      IF(TE.GT.5.0E3) THEN
! Te gt than 5000 Kelvin: calculate S, and derive alpha
        CALL EIRENE_EXCOFF(U,OSC,TEMP,C,F,S,ALPHA)

        DO 105 I=1,40
          ALPHA(I)=S(I)*SAHA(I)
  105   CONTINUE

      ELSE
! Te lt than 5000 Kelvin: calculate ALPHA, and derive S for transition to n=1
!                         calculate S, and derive ALPHA for n=2,...40
        CALL EIRENE_EXCOFF(U,OSC,TEMP,C,F,S,ALPHA)
        S(1)=ALPHA(1)/SAHA(1)
        DO 19 I=2,40
          ALPHA(I)=S(I)*SAHA(I)
   19   CONTINUE
      END IF

c   next: radiative recombination: BETA

      DO 602 I=1,40

        P=REAL(I,DP)
        XP=UH/TEMP/P**2
        EP=UH/P**2
        CALL EIRENE_CLBETA(XP,P,XS,EXS) ! return XS,EXS for
c                  rad. rec. rate and electron energy-weighted rate, both: into P state at T= temp,


        BETA(I) = 5.197D-14*(UH/TEMP)**.5/P   *XS
        EBETA(I)= 5.197D-14*(UH/TEMP)**.5/P*EP*EXS

c  the high temp behaviour (above about 500 - 1000 eV) is spurious, and gets very wrong for ebeta
c         if (i.eq.1) then
c           ratio=ebeta(i)/beta(i)
c           write (iunout,*) 'calc: temp, i, B, EB, RATIO ',
c    .                              temp, I, BETA(i),EBETA(I), RATIO
c         endif

C  for high Te use the d ln(K)/d ln(T) from the old Claudine (Mahn-Welge) code.
          IF (TEMP.GT.500.0_DP) THEN  ! tried also higher temp,
                                      ! from 800 eV: spurious
            iii=min(i,26)  ! assume the same correction term
                           ! for all I > 26
            ELN=LOG(TEMP)
C  EKIN=  <SIGMA * V * EKIN(ELEC)>   (EV*CM**3/S)
            DCCX=RECOMB(9,III)*8
            DO 57 II=1,7
              JJ=8-II
              DCCX=DCCX*ELN+RECOMB(JJ+1,III)*JJ
   57       CONTINUE
            EBETA(I)=TEMP*BETA(I)*(DCCX+1.5)
          ENDIF

  602 CONTINUE
c
      if (lopaque) beta(1)=0._DP
c

      RETURN
      END SUBROUTINE EIRENE_RATCOF

C***********************************************************************
      SUBROUTINE EIRENE_EXCOFF(U,OSC,TEMP,C,F,S,ALPHA)
C
C     EXCITATION RATE COEFFICIENT
C
C
C
      USE EIRMOD_PRECISION
      USE EIRMOD_CCONA, ONLY : EVKEL
      IMPLICIT NONE
      REAL(DP) TEMP
      REAL(DP) U(40,40),OSC(40,40),C(40,40),F(40,40)
      REAL(DP) S(40),ALPHA(40)
      INTEGER I, J
      REAL(DP) AL, CJ, CV, F1, GG, GGG, P, Q, S1, SJ, SV, TE
      EXTERNAL EIRENE_COF1N, EIRENE_COFJO, EIRENE_COFJS, EIRENE_COFJS2,
     .         EIRENE_COFVR, EIRENE_COFVS
C
      S=0.0_DP
      ALPHA=0.0_DP
      C=0.0_DP
      F=0.0_DP

      TE=TEMP/EVKEL

C*********  1 -> J
      I=1
      P=REAL(I,DP)
      DO 100 J=2,40
        Q=REAL(J,DP)
        CALL EIRENE_COF1N(U(I,J),OSC(I,J),TE,F1,I,J)
        F(J,1)=F1
        C(1,J)=Q**2/EXP(U(1,J))*F(J,1)  !/P**2, but P=1 here
  100 CONTINUE

C*********  2-10 -> J
      DO 110 I=2,10
        P=REAL(I,DP)
        DO J=I+1,40
          Q=REAL(J,DP)
          CALL EIRENE_COFVR(U(I,J),OSC(I,J),TEMP,CV,I,J)
          CALL EIRENE_COFJO(U(I,J),OSC(I,J),TE,CJ,I,J)

          GG=((P-2.0_DP)/8.0_DP)**0.25
          C(I,J)=(1.0_DP-GG)*CJ+GG*CV
          F(J,I)=P**2/Q**2*EXP(U(I,J))*C(I,J)
        END DO
  110 CONTINUE


C*********  I(>11) -> J

      DO 120 I=11,39
        P=REAL(I,DP)
        DO J=I+1,40
          Q=REAL(J,DP)
          CALL EIRENE_COFVR(U(I,J),OSC(I,J),TEMP,CV,I,J)
          C(I,J)=CV
          F(J,I)=P**2/Q**2*EXP(U(I,J))*C(I,J)
        END DO
  120 CONTINUE

C*********  S  1 ->  ionization
      I=1

      IF(TE.GT.5.0E3) THEN
        CALL EIRENE_COFJS(TE,S1,I)
        S(1)=S1

      ELSE      !  Te  <= 5000
        CALL EIRENE_COFJS2(TE,AL,I)
        ALPHA(1)=AL
      END IF
C
C*********  S  2-10 ->
      DO 210 I=2,10
        P=REAL(I,DP)

        CALL EIRENE_COFJS(TE,SJ,I)
        CALL EIRENE_COFVS(TEMP,SV,I)

        GGG=((P-2.0_DP)/8.0_DP)**0.25
        S(I)=(1.0_DP-GGG)*SJ+GGG*SV
  210 CONTINUE

C*********  S  I(>11) ->
      DO 220 I=11,40
        CALL EIRENE_COFVS(TEMP,SV,I)
        S(I)=SV
  220 CONTINUE
      RETURN
      END SUBROUTINE EIRENE_EXCOFF

C***********************************************************************
      SUBROUTINE EIRENE_COF1N(U,OSC,TE,F,I,J)
C
C     C(1,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN 1 -> J
C
C
C     K.SAWADA, K.ERIGUCHI, T.FUJIMOTO
C     J. APPL. PHYS. 73, 8122 (1993).
C
C
C
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE

      REAL(DP) U, OSC, TE, F
      INTEGER I, J, ICON
      REAl(DP) A, B, C1, E1, E1Y, E1Z, E2, E2Y, E2Z,
     .         R, RX, P, Q, X, Y, Y1, Z, Z1
      REAL(DP) EIRENE_GINT
      EXTERNAL EIRENE_GINT, EIRENE_EXPI

      P=REAL(I,DP)
      Q=REAL(J,DP)
      X=1.0_DP-P**2/Q**2
      Y=U
C
      R=0.45_DP*X
      IF(J.EQ.2) THEN
        RX=-0.95_DP
      ELSE IF(J.EQ.3) THEN
        RX=-0.95_DP
      ELSE IF(J.EQ.4) THEN
        RX=-0.95_DP
      ELSE IF(J.EQ.5) THEN
        RX=-0.95_DP
      ELSE IF(J.EQ.6) THEN
        RX=-0.94_DP
      ELSE IF(J.EQ.7) THEN
        RX=-0.93_DP
      ELSE IF(J.EQ.8) THEN
        RX=-0.92_DP
      ELSE IF(J.EQ.9) THEN
        RX=-0.91_DP
      ELSE
        RX=-0.9_DP
      END IF
C
      Z=R+Y
      A=2.0_DP*P**2*OSC/X
      B=4.0_DP*P**4/Q**3/X**2*(1.0_DP+1.3333_DP/X-0.603_DP/X**2)
      C1=2.0_DP*P**2/X
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z

      IF(Y.GT.50.0_DP) THEN
        E1Y=EXP(-Y)*EIRENE_GINT(Y)/Y
        E1Z=EXP(-Z)*EIRENE_GINT(Z)/Z
        E2Y=EXP(-Y)*(1.0-EIRENE_GINT(Y))
        E2Z=EXP(-Z)*(1.0-EIRENE_GINT(Z))
        E1=(1.0_DP/Y+0.5_DP)*E1Y+RX*(1.0_DP/Z+0.5_DP)*E1Z
        E2=E2Y/Y+RX*E2Z/Z
        F=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)*P**2/Q**2*EXP(Y)
      ELSE

cdr  use exponential integral here

        CALL EIRENE_EXPI(Y1,E1Y,ICON)
        CALL EIRENE_EXPI(Z1,E1Z,ICON)

        E1Y=-E1Y
        E1Z=-E1Z
        E2Y=EXP(-Y)-Y*E1Y
        E2Z=EXP(-Z)-Z*E1Z
        E1=(1.0_DP/Y+0.5_DP)*E1Y+RX*(1.0_DP/Z+0.5_DP)*E1Z
        E2=E2Y/Y+RX*E2Z/Z

        F=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)*P**2/Q**2*EXP(Y)
      END IF
      RETURN
      END SUBROUTINE EIRENE_COF1N

C***********************************************************************
      SUBROUTINE EIRENE_COFJO(U,OSC,TE,C,I,J)
C
C     C(I,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN I -> J
C
C     L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972).
C
C
C
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
      REAL(DP) U, OSC, TE, C
      INTEGER I, J, ICON
      REAL(DP) A, B, C1, BN, E1, E1Y, E1Z, E2, E2Y, E2Z, P, Q, R,
     .         X, Y, Y1, Z, Z1
      EXTERNAL EIRENE_EXPI, EIRENE_EXIT_OWN

      P=REAL(I,DP)
      BN=(4.0_DP-18.63_DP/P+36.24_DP/P**2-28.09_DP/P**3)/P
      Q=REAL(J,DP)
      X=1.0_DP-P**2/Q**2
      Y=U
      R=1.94_DP/P**1.57_DP*X
      Z=R+Y
      A=2.0_DP*P**2*OSC/X
      B=4.0_DP*P**4/Q**3/X**2*(1.0_DP+1.3333_DP/X+BN/X**2)
      C1=2.0_DP*P**2/X
      B=B-A*LOG(C1)
      Y1=-Y
      Z1=-Z
      CALL EIRENE_EXPI(Y1,E1Y,ICON)
      CALL EIRENE_EXPI(Z1,E1Z,ICON)
      IF (ICON.NE.0) GOTO 1000
      E1Y=-E1Y
      E1Z=-E1Z
      E2Y=EXP(-Y)-Y*E1Y
      E2Z=EXP(-Z)-Z*E1Z
      E1=(1.0_DP/Y+0.5_DP)*E1Y-(1.0_DP/Z+0.5_DP)*E1Z
      E2=E2Y/Y-E2Z/Z

      C=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)

      RETURN
 1000 WRITE(iunout,*) 'ERROR IN COFJO        ICON = ',ICON
      CALL EIRENE_EXIT_OWN(1)
      END SUBROUTINE EIRENE_COFJO

C***********************************************************************
      SUBROUTINE EIRENE_COFVR(U,OSC,TEMP,C,I,J)
C
C     C(I,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN I -> J
C
C     L.VRIENS, A.H.M.SMEETS, PHYS. REV. A 22, 940 (1980).
C
C
C
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP) U, OSC, TEMP, C
      INTEGER I, J
      REAL(DP) A, B, C1, C2, BN, DELTA,
     .         G1, G2, G3, GAMMA, P, Q, S1, UH, X

      UH=13.595_DP
      P=REAL(I,DP)
      BN=1.4_DP/P*LOG(P)-0.7_DP/P-0.51_DP/P**2+1.16_DP/P**3-0.55_DP/P**4
      Q=REAL(J,DP)
      S1=Q-P
      X=1.0_DP-P**2/Q**2
      A=2.0_DP*P**2*OSC/X
      B=4.0_DP*P**4/Q**3/X**2*(1.0_DP+1.3333_DP/X+BN/X**2)
      DELTA=EXP(-B/A)+0.06_DP*S1**2/P**2/Q
      G1=1.0_DP+TEMP/UH*P**3
      G2=3.0_DP+11.0_DP*S1**2/P**2
      G3=6.0_DP+1.6_DP*S1*Q+0.3_DP/S1**2+
     .   0.8_DP*Q**1.5/S1**0.5*(S1-0.6_DP)
      GAMMA=UH*LOG(G1)*G2/G3
      C1=1.6D-7*TEMP**0.5/(TEMP+GAMMA)*EXP(-U)
      C2=0.3_DP*TEMP/UH+DELTA

      C=C1*(A*LOG(C2)+B)

      RETURN
      END SUBROUTINE EIRENE_COFVR

C***********************************************************************
      SUBROUTINE EIRENE_COFJS(TE,S,I)
C
C     S(I); IONIZATION RATE COEFFICIENT FOR ATOMIC HYDROGEN  I ->
C
C
C
C     K.SAWADA, K.ERIGUCHI, T.FUJIMOTO
C     J. APPL. PHYS. 73, 8122 (1993).
C
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP) TE, S
      INTEGER I, N, K, ICON
      REAL(DP) G(0:2,40)
      REAL(DP) A, B, BN, C1,
     .         E0Y, E0Z, E1, E1Y, E1Z, E2, E2Y, E2Z, EGY, EGZ,
     .         P, R, RX, Y, Y1, Z, Z1
      EXTERNAL EIRENE_EXPI

      G(0,1)=1.1330_DP
      G(1,1)=-0.4059_DP
      G(2,1)=0.07014_DP
      G(0,2)=1.0785_DP
      G(1,2)=-0.2319_DP
      G(2,2)=0.02947
      DO 350 N=3,40
        G(0,N)= 0.9935_DP     +0.2328_DP/N   -0.1296_DP/N**2
        G(1,N)=-0.6282_DP/N   +0.5598_DP/N**2-0.5299_DP/N**3
        G(2,N)= 0.3887_DP/N**2-1.1810_DP/N**3+1.4700_DP/N**4
  350 CONTINUE

      IF (I.EQ.1) THEN

        P=1.0_DP

        Y=1.57770E5_DP/TE
        R=0.45_DP
        Z=R+Y
        RX=-0.59_DP
        A=0.0_DP
        DO 223 K=0,2
          A=A+G(K,1)/(K+3)
  223   CONTINUE
        A=A*1.9603_DP*P

        B=0.66667_DP*P**2*(5.0_DP-0.603_DP)
        C1=2.0_DP*P**2
        B=B-A*LOG(C1)
        Y1=-Y
        Z1=-Z
        CALL EIRENE_EXPI(Y1,E1Y,ICON)
        CALL EIRENE_EXPI(Z1,E1Z,ICON)

        E1Y=-E1Y
        E1Z=-E1Z
        E2Y=EXP(-Y)-Y*E1Y
        E2Z=EXP(-Z)-Z*E1Z
        E1=1/Y*E1Y+RX*1/Z*E1Z
        E0Y=EXP(-Y)/Y
        E0Z=EXP(-Z)/Z
        EGY=E0Y-2*E1Y+E2Y
        EGZ=E0Z-2*E1Z+E2Z
        E2=EGY+RX*EGZ

        S=1.093D-10*SQRT(TE)*P**2*Y**2*(A*E1+B*E2)

      ELSE
        P=REAL(I,DP)
        BN=(4.0_DP-18.63_DP/P+36.24_DP/P**2-28.09_DP/P**3)/P
        Y=1.57770E5_DP/TE/P**2
        R=1.94_DP/P**1.57_DP
        Z=R+Y
        A=0.0_DP
        DO 222 K=0,2
          A=A+G(K,I)/(K+3)
  222   CONTINUE
        A=A*1.9603_DP*P
        B=0.66667_DP*P**2*(5.0_DP+BN)
        C1=2.0_DP*P**2
        B=B-A*LOG(C1)
        Y1=-Y
        Z1=-Z
        CALL EIRENE_EXPI(Y1,E1Y,ICON)
        CALL EIRENE_EXPI(Z1,E1Z,ICON)
        E1Y=-E1Y
        E1Z=-E1Z
        E2Y=EXP(-Y)-Y*E1Y
        E2Z=EXP(-Z)-Z*E1Z
        E1=1/Y*E1Y-1/Z*E1Z
        E0Y=EXP(-Y)/Y
        E0Z=EXP(-Z)/Z
        EGY=E0Y-2*E1Y+E2Y
        EGZ=E0Z-2*E1Z+E2Z
        E2=EGY-EGZ

        S=1.093D-10*SQRT(TE)*P**2*Y**2*(A*E1+B*E2)
      END IF
      RETURN
      END SUBROUTINE EIRENE_COFJS

C***********************************************************************
      SUBROUTINE EIRENE_COFJS2(TE,W,I)
C
C     S(I); IONIZATION RATE COEFFICIENT FOR ATOMIC HYDROGEN  I ->
C     FOR LOW TEMPERATURE
C
C
C     K.SAWADA, K.ERIGUCHI, T.FUJIMOTO
C     J. APPL. PHYS. 73, 8122 (1993).
C
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP) TE, W
      INTEGER I, N, K
      REAL(DP) G(0:2,40)
      REAL(DP) A, B, C1, E0Y, E0Z, E1, E1Y, E1Z, E2, E2Y, E2Z, EGY, EGZ,
     .         P, PP, R, RX, Y, Y1, Z, Z1
      REAL(DP) EIRENE_GINT
      EXTERNAL EIRENE_GINT

      G(0,1)= 1.1330_DP
      G(1,1)=-0.4059_DP
      G(2,1)= 0.07014_DP
      G(0,2)= 1.0785_DP
      G(1,2)=-0.2319_DP
      G(2,2)= 0.02947_DP
      DO 350 N=3,40
        G(0,N)= 0.9935_DP     +0.2328_DP/N   -0.1296_DP/N**2
        G(1,N)=-0.6282_DP/N   +0.5598_DP/N**2-0.5299_DP/N**3
        G(2,N)= 0.3887_DP/N**2-1.1810_DP/N**3+1.4700_DP/N**4
  350 CONTINUE

C     IF (I.EQ.1) THEN
        PP=REAL(I,DP)
        P=1.0_DP
        Y=1.57770E5_DP/TE
        R=0.45_DP
        Z=R+Y
        RX=-0.59_DP

        A=0.0_DP
        DO 223 K=0,2
          A=A+G(K,1)/(K+3)
  223   CONTINUE
        A=A*1.9603_DP*P

        B=0.66667_DP*P**2*(5.0_DP-0.603_DP)
        C1=2.0_DP*P**2
        B=B-A*LOG(C1)
        Y1=-Y
        Z1=-Z

        E1Y=EXP(-Y)*EIRENE_GINT(Y)/Y
        E1Z=EXP(-Z)*EIRENE_GINT(Z)/Z
        E2Y=EXP(-Y)*(1.0_DP-EIRENE_GINT(Y))
        E2Z=EXP(-Z)*(1.0_DP-EIRENE_GINT(Z))
        E1=1.0_DP/Y*E1Y+RX*1.0_DP/Z*E1Z
        E0Y=EXP(-Y)/Y
        E0Z=EXP(-Z)/Z
        EGY=E0Y-2*E1Y+E2Y
        EGZ=E0Z-2*E1Z+E2Z
        E2=EGY+RX*EGZ

        W=1.093D-10*SQRT(TE)*P**2*Y**2*(A*E1+B*E2)*
     *    EXP(13.595_DP/TE*1.1605E4_DP)/2.414D15/SQRT(TE**3)

      RETURN
      END SUBROUTINE EIRENE_COFJS2

C***********************************************************************
      SUBROUTINE EIRENE_COFVS(TEMP,S,I)
C
C     S(I); IONIZATION RATE COEFFICIENT FOR ATOMIC HYDROGEN  I ->
C
C
C     L.VRIENS, A.H.M.SMEETS, PHYS. REV. A 22, 940 (1980).
C
C
C
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP) TEMP, S
      INTEGER I
      REAL(DP) P, UI, UIZ

      P=REAL(I,DP)
      UI=13.595_DP/TEMP/P**2
      UIZ=UI**2.33_DP+4.38_DP*UI**1.72_DP+1.32_DP*UI
      S=9.56D-6/TEMP**1.5*EXP(-UI)/UIZ

      RETURN
      END SUBROUTINE EIRENE_COFVS

cdr
      FUNCTION EIRENE_GINT(XX)
      USE EIRMOD_PRECISION
      REAL(DP) EIRENE_GINT
      REAL(DP) XX
      REAL(DP) X,GG(8)
cdr
      DATA GG/0.2677737343_DP, 8.6347608925_DP,
     *       18.0590169730_DP, 8.5733287401_DP,
     *        3.9584960228_DP,21.0996530827_DP,
     *       25.6329561486_DP, 9.5733223454_DP/
cdr
      x=xx
cdr
      EIRENE_GINT=(GG(1)+GG(2)*X+GG(3)*X**2+GG(4)*X**3+X**4)/
     *            (GG(5)+GG(6)*X+GG(7)*X**2+GG(8)*X**3+X**4)
      RETURN
      END FUNCTION EIRENE_GINT


C***********************************************************************
      SUBROUTINE EIRENE_CLBETA(XP,P,S,ES)
C
C     RADIATIVE RECOMBINATION RATE COEFFICIENT, TO ATOMIC STATE P
C   input :
c      p=  bound state I, into which rad rec goes    --> common
C     xp=  uh/temp / p^2                             --> common    = ep/temp
C     note: ep = binding energy in p state. ep(1)=UH, ep(inf)=0.0
c   output:
C      S=  normalized rate coeff for rad rec to state P,
C          then in calling routine ratcof:
c          BETA(P) =5.197D-14*(UH/TEMP)**.5/P * S      = 5.197D-14*(EP/TEMP)**.5 * S
C      ES= normalized electron energy-weighted rate coefficient
C          then in calling routine ratcof:
c          EBETA(P)=5.197D-14*(UH/TEMP)**.5/P * EP* ES = 5.197D-14*(EP/TEMP)**.5 * EP *ES
C
C   integrate functions gaunt3(x), for S and gaunt4, for ES, from 0 to 20
C
      USE EIRMOD_PRECISION

      IMPLICIT NONE
      REAL(DP) XP,P,S,ES
      REAL(DP) EIRENE_GAUNT3,EIRENE_GAUNT4,PP,XPP,A,B,EPSR
      INTEGER II, NMIN, NMAX
      COMMON PP,XPP
      EXTERNAL EIRENE_GAUNT3, EIRENE_GAUNT4, EIRENE_AQC8

      II=INT(P)
      PP=P
      XPP=XP   !  uh/temp/p**2

      A=0.0_DP
      B=20.0_DP
      EPSR=1.0D-4
cdr   EPSR=1.0D-5  slowed down code by factor of 100 !!
      NMIN=15
      NMAX=511

c  integrate functions GAUNT from x=0 to 20:
      CALL EIRENE_AQC8(A,B,EIRENE_GAUNT3,EPSR,NMIN,NMAX, S)
      CALL EIRENE_AQC8(A,B,EIRENE_GAUNT4,EPSR,NMIN,NMAX,ES)
C
      RETURN
      END SUBROUTINE EIRENE_CLBETA

C***********************************************************************

      FUNCTION EIRENE_GAUNT3(X)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP) PP,XPP,U,B,X,EIRENE_GAUNT3
      COMMON PP,XPP

      U=X/XPP
      B=PP
      EIRENE_GAUNT3=(1.0_DP/(U+1.0_DP)+
     .               0.1728_DP*(U-1.0_DP)/B**(2.0_DP/3.0_DP)/
     .                         (U+1.0_DP)**(5.0_DP/3.0_DP)-
     .               0.0496_DP*(U**2+4.0_DP/3.0_DP*U+1.0_DP)/
     .                                    B**(4.0_DP/3.0_DP)/
     .               (U+1.0_DP)**(7.0_DP/3.0_DP))*EXP(-X)
      RETURN
      END FUNCTION EIRENE_GAUNT3

      FUNCTION EIRENE_GAUNT4(X)
c  same as gaunt3, but for energy-weighted rate coeff (juel rep,3858 (2001) crmol manual, , eq. 9a, 9b)
c  gaunt4= u * gaunt3,  and additional pre-factor Ep = UH/(p^2), additional to (EP/T)^...

cdr careful: integration of gaunt4 fails above Te gt 4500 eV
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      REAL(DP) PP,XPP,U,B,X,EIRENE_GAUNT4
      COMMON PP,XPP
      U=X/XPP
      B=PP
      EIRENE_GAUNT4=U*
     .    (1.0_DP/(U+1.0_DP)+0.1728_DP*(U-1.0_DP)/B**(2.0_DP/3.0_DP)/
     .    (U+1.0_DP)**(5.0_DP/3.0_DP)-0.0496_DP*
     .    (U**2+4.0_DP/3.0_DP*U+1.0_DP)/B**(4.0_DP/3.0_DP)/
     .    (U+1.0_DP)**(7.0_DP/3.0_DP))*EXP(-X)
      RETURN
      END FUNCTION EIRENE_GAUNT4

C***********************************************************************
      SUBROUTINE EIRENE_POPCOF_M(DENSEL,SAHA,C,F,S,A,ALPHA,BETA,LUP,LIM,
     &      R0,R1,R_EXT,
     &            Q_EXT,L_EXT)
C
C     SOLUTION OF RATE EQUATION FOR ATOMIC HYDROGEN
C
C  COPY OF ORIGINAL ROUTINE EIRENE_POPCOF FOR SIMULTANEOUS SOLUTION
c  WITH OF UP TO 3 RIGHT HAND SIDES. (DEFAULT: 2: IONISATION, RECOMBINATION)
c  l_ext: indicate that 3rd right hand side term is requested
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER LUP, LIM
      REAL(DP) DENSEL
      REAL(DP) C(40,40),F(40,40),A(40,40),W(40,40),
     &         SAHA(40),S(40),ALPHA(40),BETA(40),R0(40),R1(40),
     &         Q_EXT(40),R_EXT(40),
     &         VW(40),WA(40,40)
      REAL(DP) :: BLAX(3,40)
      INTEGER I,J,K,L,IE,ICON,IS
      LOGICAL :: L_EXT
      integer ip(40)
      REAL(DP) SUMA, SUMAS, SUMC, SUMF, SUMFS
      EXTERNAL :: EIRENE_LAX_M

      DO 201 K=2,LUP-1

        DO 202 L=2,K
cdr stoss bevoelkerung von k von unten
          W(K,L)=C(L,K)*DENSEL
  202   CONTINUE

cc diagonale
cc entvoelkerung durch stoesse nach unten
        SUMF=0.0_DP
        DO 301 I=1,K-1
          SUMF=SUMF+F(K,I)
  301   CONTINUE
cc entvoelkerung durch stoesse nach oben
        SUMC=0.0_DP
        DO 302 I=K+1,LIM
          SUMC=SUMC+C(K,I)
  302   CONTINUE
cc  spontan nach unten
        SUMA=0.0_DP
        DO 303 I=1,K-1
          SUMA=SUMA+A(K,I)
  303   CONTINUE
cdr entvoelkerung von k: stoesse nach unten, nach oben, ionis, spontan
cdr                      nach unten
        W(K,K)=-(DENSEL*(SUMF+SUMC+S(K))+SUMA)

cc diagonale fertig

        DO 203 L=K+1,LUP
cdr bevoelkerung durch: stoesse von oben, spontan von oben
          W(K,L)=DENSEL*F(L,K)+A(L,K)
  203   CONTINUE

  201 CONTINUE
cdr k loop finished, k=2, lup-1 (i.e. without last state lup)


c  special treatment of last state lup: 2-->lup, 3-->lup,..., only for excitation to it,
c                                  but not for de-exit to it, also not for rad rec to it.
      DO 211 L=2,LUP-1
        W(LUP,L)=C(L,LUP)*DENSEL
  211 CONTINUE
c  beitrag des letzten zustandes lup zu diagonal
      SUMF=0.0_DP
      DO 311 I=1,LUP-1
        SUMF=SUMF+F(LUP,I)
  311 CONTINUE
      SUMC=0.0_DP
      DO 313 I=LUP+1,LIM   !Boltzmann LTE contribution for LIM gt. LUP
        SUMC=SUMC+C(LUP,I)
  313 CONTINUE
      SUMA=0.0_DP
      DO 312 I=1,LUP-1
        SUMA=SUMA+A(LUP,I)
  312 CONTINUE

      W(LUP,LUP)=-(DENSEL*(SUMF+SUMC+S(LUP))+SUMA)

C  RECHTE SEITEN:
      DO 550 K=2,LUP

c  vorbereiten fuer recombination, Saha correction
        SUMFS=0.0_DP
        DO 500 I=LUP+1,LIM
          SUMFS=SUMFS+F(I,K)*SAHA(I)
  500   CONTINUE
        SUMAS=0.0_DP
        DO 501 I=LUP+1,LIM
          SUMAS=SUMAS+SAHA(I)*A(I,K)
  501   CONTINUE
c
c  matrixelemente: 1/s  (densel*rate coeff. )
c  rechte seiten : cm**3/s, nicht: 1/s, also fuer elektronendichte=1
c                                      (bzw: stosspartnerdichte =1)
c  geht wg. linearitaet.
c  recombination e + H+ --> H*
        W(K,LUP+1)=-(DENSEL*SUMFS+SUMAS+(DENSEL*ALPHA(K)+BETA(K)))
c  ionisation e + H --> H*
        W(K,LUP+2)=-C(1,K)
c  external source: Q_EXT
        W(K,LUP+3)=-Q_EXT(K)
  550 CONTINUE

cdr w besetzt fuer w(i,j) i=2,lup,j=2,lup+3
cdr OK, AS LONG AS lup<38


cdr reduziere w indices um 1: auf wa: i=1,lup-1,j=1,(lup-1)+3

      DO 402 I=1,LUP-1
        DO J=1,LUP-1+3
          WA(I,J)=W(I+1,J+1)
        END DO
  402 CONTINUE

c  two or three right linearly additive hand side terms?
      ie=2
      if (L_EXT) ie=3
      DO 3000 J=1,LUP-1
        BLAX(1:ie,J)=WA(J,LUP:LUP+ie-1)
 3000 CONTINUE

      CALL EIRENE_LAX_M(WA,40,LUP-1,BLAX,3,ie,0.0_DP,1,IS,VW,IP,ICON)
c

        DO J=1,LUP-1
C coupling to H+
            R0(J+1)   =BLAX(1,J)
C coupling to H-groundstate
            R1(J+1)   =BLAX(2,J)
        ENDDO
C coupling to Q_EXT: external sources for H*(n) population
        IF (L_EXT) THEN
          DO J=1,LUP-1
            R_EXT(J+1)=BLAX(3,J)
          ENDDO
        ELSE
          R_EXT=0.
        ENDIF

      RETURN
      END SUBROUTINE EIRENE_POPCOF_M

C***********************************************************************
      SUBROUTINE EIRENE_IONREC
     &                 (C,S,SAHA,A,ALPHA,BETA,R0,R1,DENSEL,LUP,LIM,
     &                  F,R_EXT,Q_EXT,L_EXT,
     &                  ALPCR,SCR,
     &                  ALP_EXT,SCR_EXT)
C
C     EFFECTIVE IONIZATION AND RECOMBINATION RATE COEFFICIENTS
C     FOR ATOMIC HYDROGEN
C  1) FROM GROUND STATE TO CONTINUUM: SCR
C  2) FROM CONTINUUM TO GROUND STATE: ALPCR
C  3) FROM EXTERNAL TO GROUND OR CONTINUUM: SCR_EXT, ALP_EXT
C
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER LUP, LIM
      REAL(DP) ALPCR, DENSEL, SCR, SCR_EXT
      REAL(DP) ALP_EXT, ALP1_EXT, ALP2_EXT
      REAL(DP) C(40,40),S(40),SAHA(40),A(40,40),ALPHA(40),BETA(40),
     &         R0(40),R1(40),R_EXT(40),Q_EXT(40),F(40,40)
      INTEGER I
      REAL(DP) ALPCR1, ALPCR2, SUSCR, SUSRAD
      LOGICAL :: L_EXT
C
      DO 5000 I=LUP+1,LIM
        R0(I)=1.0_DP*SAHA(I)
        R1(I)=0.0_DP
        R_EXT(I)=0.0_DP
 5000 CONTINUE
C
C  IONIS. TRANSITION TO CONTINUUM
      SCR=S(1)
      DO 5001 I=2,LUP
        SUSCR=C(1,I)-R1(I)*(F(I,1)*DENSEL+A(I,1))
        SCR=SCR+SUSCR
 5001 CONTINUE

C  RECOMB. TRANSITION TO (1)
      ALPCR1=DENSEL*ALPHA(1)+BETA(1)

      ALPCR2=0.0_DP

      DO 5003 I=2,LIM
        ALPCR2=ALPCR2+R0(I)*(DENSEL*F(I,1)+A(I,1))
 5003 CONTINUE

      ALPCR=ALPCR1+ALPCR2

C  FROM EXTERNAL, TRANSITION TO CONTINUUM
      SCR_EXT=0.0_DP

      IF (.NOT.L_EXT) RETURN

      DO 5002 I=2,LUP
        SUSRAD=Q_EXT(I)-R_EXT(I)*(F(I,1)*DENSEL+A(I,1))
        SCR_EXT=SCR_EXT+SUSRAD
 5002 CONTINUE

C  FROM EXTERNAL, TRANSITION TO (1)
      ALP1_EXT=Q_EXT(1)

      ALP2_EXT=0.0

      DO 5004 I=2,LIM
        ALP2_EXT=ALP2_EXT+R_EXT(I)*(DENSEL*F(I,1)+A(I,1))
 5004 CONTINUE

      ALP_EXT=ALP1_EXT+ALP2_EXT

      RETURN
      END SUBROUTINE EIRENE_IONREC

C***********************************************************************
      SUBROUTINE EIRENE_E_IONREC
     &                   (C,S,SAHA,A,ALPHA,BETA,EBETA,
     &                    R0,R1,DENSEL,LUP,LIM,
     &                    F,R_EXT,Q_EXT,L_EXT,E_AT,
     &                    ALPCR,      SCR,    SCR_EXT,            !  ordinary rate coeffcients
     &                    E_ALPCR,  E_SCR,  E_SCR_EXT, E_ALP_EXT  !  electron energy-weighted rate coefficients
ctt  &                   ,E_ALPCR_T,E_SCR_T,E_SCR_EXT_T           !  radiation energy losses only, for consistency testing
     &                    )
C
C     EFFECTIVE ELECTRON ENERGY LOSS IONIZATION AND RECOMBINATION
C     RATE COEFFICIENTS FOR ATOMIC HYDROGEN AT DENSITY DENSEL, TEMPERATURE TEMP

C   energies are counted from the point of view of electrons: energy loss: negative, energy gain: positive
C   RADIATED ENERGIES ARE TAKEN NEGATIVE (LOSS).. IF THEY ARE POSITIVE (E.G. IN TEST CASES) Q_EXT,..ETC.
C   THEN ENERGY MUST BE ABSORBED FROM THE RADIATION FIELD TO ENABLE A PARTICULAR TRANSIION.
C   A WARNING IS THEN PRINTED

C     E_SCR,  E_ALPCR,....  IS TOTAL ELECTRON COOLING RATE. NEGATIVE IF LOSS FOR ELECTRONS, POSITIVE ELSE.
C     E_SCR_T, ....._T,.... IS RADIATION ENERGY RATE ONLY, SHOULD BE NEGATIVE ALWAYS.

C     ENERGY LEVELS: E_AT(1)=0, E_AT(H+)=13.595=UH
C
C
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER LUP, LIM
      REAL(DP) ALPCR, DENSEL, E_ALPCR, E_SCR, E_SCR_EXT, SCR
      REAL(DP) SCR_EXT, ALP_EXT, E_ALP_EXT
      REAL(DP) C(40,40),S(40),SAHA(40),A(40,40),F(40,40),
     &         ALPHA(40),BETA(40),EBETA(40),
     &         R0(40),R1(40),R_EXT(40),Q_EXT(40),E_AT(40)
      INTEGER I, I1, I2
      REAL(DP) DE, SUSCR, UH
      LOGICAL :: L_EXT

      REAL(DP) EMEAN_REC(41)
C
      DO 5000 I=LUP+1,LIM
        R0(I)=1.0_DP*SAHA(I)
        R1(I)=0.0_DP
cdr  no Saha-fringe for external contributions
        R_EXT(I)=0.0_DP
 5000 CONTINUE
C
      E_SCR=0.0_DP
      E_SCR_EXT=0.0_DP
      E_ALPCR=0.0_DP
      E_ALP_EXT=0.0_DP !rc exists but not used
C  FOR TESTING: ONLY CUMULATE RADIATION LOSSES HERE.
ctt   E_SCR_T=0.0_DP
ctt   E_SCR_EXT_T=0.0_DP
ctt   E_ALPCR_T=0.0_DP
C
C  EFFECTIVE ELECTRON COOLING CORRESPONDING TO
C  "ORDINARY" COUPLING TO GROUND STATE S(I)
C
C  PART I
C  1(EXTERN)--> inf.  R1(1)=1
      UH=13.595_DP
      DE=(E_AT(1)-UH)
      E_SCR  =S(1)*DE  !  *POPULATION OF H(1) = R1(1) =1
C  1(EXTERN) --> 1  !  here no contribution

C  PART II
C  i1--> inf.  R1(i1)=...
      do i1=2,lup
        DE=E_AT(i1)-UH
        E_SCR  =E_SCR+r1(i1)*densel*S(I1)*DE
      enddo

C  PART III
C  1--> I2.  R1(1)=1
      DO I2=2,LUP
        DE=(E_AT(1)-E_AT(I2))
        SUSCR=C(1,I2)*DE
        E_SCR=E_SCR+SUSCR
      ENDDO
C  PART IV
C  i1--> i2. i1<i2,  R1(i1)=...
      DO  I1=2,LUP-1
        DO  I2=I1+1,LUP
          DE=(E_AT(I1)-E_AT(I2))
          SUSCR=r1(i1)*densel*C(I1,I2)*DE
          E_SCR=E_SCR+SUSCR
        ENDDO
      ENDDO

C  PART V
C  i2--> i1.  i2>i1, R1(i2)=...
C            (includes i1=1, i.e., inverse to PART III)
      DO I1=1,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I2)-E_AT(I1))  !  POSITIVE, GAIN FOR ELECTRON ENERGY
          SUSCR  =r1(i2)*densel*F(I2,I1)*DE
          E_SCR=E_SCR+SUSCR
C  separate treatment of radiation losses alone, only for testing
ctt       SUSCR_T=r1(i2)*A(I2,I1)*(-1.)*DE   !  NEGATIVE, RADIATIVE LOSS
ctt       E_SCR_T=E_SCR_T+SUSCR_T
        ENDDO
      ENDDO

C  for test only.  evaluate second formula for E_SCR,
C                  using radiation loss E_SCR_T and effective rate SCR
ctt   UH=13.595_DP
ctt   DE=(E_AT(1)-UH)  !  POTENTIAL DIFFERENCE IS NEGATIVE, ELECTRONS LOSE ENERGY IN IONISATION
ctt   E_SCR_T=E_SCR_T+SCR *DE

C  contribution for "ordinary" coupling to ground state done



c.....................................H+ RECOMBINATION, EXTERNAL: H(INF) = CONTIN = H+.........................................
c  next: contribution for coupling to H+

C  EFFECTIV ELECTRON COOLING CORRESPONDING TO
C  "ORDINARY" COUPLING TO CONTINUUM: E_ALPCR

C  WE NEED MEAN ELECTRON KINETIC ENERGY LOST IN RECOMBINATION TO STATE I: EMEAN_REC(I)  (TAKEN NEGATIVE)
      do I=1,lup
        emean_rec(i)=-EBETA(I)/BETA(I)
      enddo

C
C  PART I

C  CONTIN(EXTERN)--> CONTIN

C  CONTIN(EXTERN)--> 1 , POPULATION R0 OF H+: = 1
C
      UH=13.595_DP
      E_ALPCR  =DENSEL*ALPHA(1)*(UH-E_AT(1))  !  * POPULATION OF H+
      ! ENERGY GAINED BY ELECTRON GAS, THREE-BODY, SPECTATOR ELECTRON
      E_ALPCR  =E_ALPCR +BETA(1)*EMEAN_REC(1) !  * POPULATION OF H+
      ! ENERGY LOST FROM ELECTRON GAS, RAD REC, NEGATIVE
C  for test only: SEPARATE TREATMENT OF RADIATION LOSS ALONE
ctt   DE=0.
ctt   E_ALPCR_T  =DENSEL*ALPHA(1)*DE    ! no radiation loss involved in three-body rec. ALPHA
ctt   DE=EMEAN_REC(1)-(UH-E_AT(1))      ! TOTAL ENERGY CARRIED AWAY IN RADIATION IN RAD REC, SHOULD BE NEGATIVE
ctt   E_ALPCR_T  =E_ALPCR_T+BETA(1)*DE  !  * POPULATION OF H+,   radiation loss (continuum) involved in rad. rec BETA


C  PART II
C  i1--> inf.  R0(i1)=...  !negative, loss for electron energy
      do i1=2,lup
        DE=(E_AT(i1)-UH)
        E_ALPCR  =E_ALPCR+r0(i1)*densel*S(I1)*DE
      enddo

C  PART III
C  cont--> I2.  POPULATION R0 OF H+: = 1  !  gain from three body rec ALPHA, loss of electron energy in rad. rec BETA
      DO I2=2,LUP

        SUSCR=DENSEL*ALPHA(I2)*(UH-E_AT(I2))  ! THREE BODY
        SUSCR=SUSCR+BETA(I2)*EMEAN_REC(I2)        ! RAD REC
        E_ALPCR=E_ALPCR+SUSCR
C  FOR TEST ONLY: SEPARATE TREATMENT OF RADIATION LOSS ALONE
ctt     DE=0.
ctt     E_ALPCR_T  =DENSEL*ALPHA(I2)*DE    ! no radiation loss involved in three-body rec. ALPHA
ctt     DE=EMEAN_REC(I2)-(UH-E_AT(I2))     ! TOTAL ENERGY CARRIED AWAY IN RADIATION, ALWAYS NEGATIVE (LOSS)
ctt     E_ALPCR_T  =E_ALPCR_T+BETA(I2)*DE  !  * POPULATION OF H+,   no radiation loss involved in three-body rec. ALPHA
      ENDDO

C  PART IV
C  i1--> i2. i1<i2,  R0(i1)=...
      DO I1=2,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I1)-E_AT(I2))
          SUSCR=r0(i1)*densel*C(I1,I2)*DE
          E_ALPCR=E_ALPCR+SUSCR
        ENDDO
      ENDDO

C  PART V
C  i2--> i1.  i2>i1, R0(i2)=...
C            (includes i1=1)
      DO I1=1,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I2)-E_AT(I1))
          SUSCR  =r0(i2)*densel*F(I2,I1)*DE
          E_ALPCR=E_ALPCR+SUSCR
C  separate treatment of radiation losses alone, only for testing
ctt       SUSCR_T=r0(i2)*A(I2,I1)*(-1.)*DE   !  NEGATIVE (= LOSS) BY ANSATZ
ctt       E_ALPCR_T=E_ALPCR_T+SUSCR_T
        ENDDO
      ENDDO

C  for test only for lupa=lima=32.  evaluate second formula for E_ALPCR,
C                  using radiation loss E_ALPCR_T, effective rate ALPCR AND DELPOT
C  TEST DONE, OK.
ctt   UH=13.595
ctt   DEALP=(UH-E_AT(1))  !  ELECTRONS GAIN POTENTIAL ENERGY PER RECOMBINATION EVENT, DEALP IS THEREFORE TAKEN POSITIVE
ctt   E_ALPCR_TT=E_ALPCR_T+ALPCR *DEALP
c  e_alpcr_TT: test quantity, should be equal to E_alpcr
c  e_alpcr_T : this would be the shifted (radiation alone part) e.g. fitted in amjuel format

c     WRITE (IUNOUT,*) 'DE,TE, E_ALPCR, E_ALPCR_TT, E_ALPCR_T , alpcr'
c     WRITE (IUNOUT,*) log10(DENSEL),TEMP,E_ALPCR,E_ALPCR_TT,E_alpcr_T,alpcr

C  contribution for "ordinary" coupling to H+ state done

c..................................................................
c  next: contribution for coupling to external source Q_EXT
      IF (.NOT.L_EXT) RETURN
c        use same codes as for E_SCR, but with Q_EXT(K) <-- C(1,K)
c                                     and S(1)=0
c        and, of course, R_EXT(k) instead of R1(k)
C        furthermore: no electron energy losses associated with Q_EXT: E_Q_EXT=0
C
C
C  PART I
C  zero, because S(1)=0
C  PART II
C  i1--> contin  R_EXT(i1)=...
      do i1=2,lup
        DE=(E_AT(i1)-UH)
        E_SCR_EXT  =E_SCR_EXT+R_EXT(i1)*densel*S(I1)*DE
      enddo
C  PART III
C  1--> I2.  R_EXT(1)=1
      DO I2=2,LUP
C  NO ELECTRON ENERGY LOSS IS ASSOCIATED WITH Q_EXT
C  STRICTLY: TOGETHER WITH Q_EXT THERE SHOULD ALSO BE AN E_Q_EXT
C  FOR EXTERNAL SOURCES.
        DE=0
        SUSCR=Q_EXT(I2)*DE
        E_SCR_EXT=E_SCR_EXT+SUSCR
      ENDDO

C  PART IV
C  i1--> i2. i1<i2,  R_EXT(i1)=...
      DO I1=2,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I1)-E_AT(I2))
          SUSCR=R_EXT(i1)*densel*C(I1,I2)*DE
          E_SCR_EXT=E_SCR_EXT+SUSCR
        ENDDO
      ENDDO

C  PART V
C  i2--> i1.  i2>i1, POPULATION R_EXT(i2)=...
C            (includes i1=1, i.e., inverse to PART III)
      DO I1=1,LUP-1
        DO I2=I1+1,LUP
          DE=(E_AT(I2)-E_AT(I1))
          SUSCR  =R_EXT(i2)*densel*F(I2,I1)*DE
          E_SCR_EXT=E_SCR_EXT+SUSCR
C  separate treatment of radiation losses alone, only for testing
ctt       SUSCR_T=R_EXT(i2)*A(I2,I1)*(-1.)*DE
ctt       E_SCR_EXT_T=E_SCR_EXT_T+SUSCR_T
        ENDDO
      ENDDO

C  for test only.  evaluate second formula for E_SCR_EXT,
C                  using radiation loss E_SCR_EXT_T and effective rate SCR_EXT
ctt   UH=13.595
ctt   DE=(E_AT(1)-UH)
ctt   E_SCR_EXT_T=E_SCR_EXT_T+SCR_EXT *DE

c  this test only works if in PART III one would set:
C       DE=(E_AT(1)-E_AT(I2))
C  But for external sources this is not necessarily the case
C  test done, 11.01.05, o.k., then test switched off.

      RETURN
      END SUBROUTINE EIRENE_E_IONREC

C********************************************************************
C            C(I,J)  FROM  JOHNSON
C *******************************************************************
      SUBROUTINE EIRENE_JOHN(TEMP,CJ,OSC)
C
C
C     C(I,J); EXCITATION RATE COEFFICIENT FOR ATOMIC
C     HYDROGEN I -> J
C
C     L.C.JOHNSON, ASTROPHYS. J. 174, 227 (1972).
C
C
C
C
C
      USE EIRMOD_PRECISION
      USE EIRMOD_CCONA, ONLY : EVKEL
      IMPLICIT NONE
      REAL(DP) TEMP
      REAL(DP) OSC(40,40),CJ(40,40)
      INTEGER I, J, ICON
      REAL(DP) A, B, BN, C1, E1, E1Y, E1Z, E2, E2Y, E2Z, EP, EQ,
     .         P, Q, R, TE, U, X, Y, Y1, Z, Z1
      EXTERNAL EIRENE_EXPI

      TE=TEMP/EVKEL
C
      DO 1 I=1,40
       DO 2 J=1,40
        P=REAL(I,DP)
        Q=REAL(J,DP)
        IF (I.EQ.1) THEN
          BN=-0.603_DP
        ELSE
          BN=(4.0_DP-18.63_DP/P+36.24_DP/P**2-28.09_DP/P**3)/P
        END IF
        EP=13.6_DP/P**2
        EQ=13.6_DP/Q**2
        U=(EP-EQ)/TEMP
C
        IF (U.GT.0.0_DP) THEN
          X=1.0_DP-P**2/Q**2
          Y=U
          R=1.94_DP/P**1.57_DP*X
          Z=R+Y
          A=2.0_DP*P**2*OSC(I,J)/X
          B=4.0_DP*P**4/Q**3/X**2*(1.0_DP+1.3333_DP/X+BN/X**2)
          C1=2.0_DP*P**2/X
          B=B-A*LOG(C1)
          Y1=-Y
          Z1=-Z
          CALL EIRENE_EXPI(Y1,E1Y,ICON)
          CALL EIRENE_EXPI(Z1,E1Z,ICON)
          E1Y=-E1Y

          E1Z=-E1Z

          E2Y=EXP(-Y)-Y*E1Y
          E2Z=EXP(-Z)-Z*E1Z
          E1=(1.0_DP/Y+0.5_DP)*E1Y-(1.0_DP/Z+0.5_DP)*E1Z
          E2=E2Y/Y-E2Z/Z

          CJ(I,J)=1.093D-10*SQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)
        ELSE
          CJ(I,J)=0.0_DP

        END IF
C
    2  CONTINUE
    1 CONTINUE
C
C
      RETURN
      END SUBROUTINE EIRENE_JOHN
c
      subroutine EIRENE_aqc8(a,b,f,epsr,nmin,nmax,S)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
c  S=integral from a to b, of the function f(x) (external).
c  epsr : relative errors, input
c   nmin,nmax:  min and max number of function calls
c   (nmax<511, nmin>15)
c
c output
c S
c err: estim absolut error
c n  anzahl der functionsaufrufe
c icon: error code
      integer nmin, nmax
      real(dp) a, b, f, epsr, s
      external f
      external EIRENE_midpnt, EIRENE_qromo
cdr  perhaps Romberg integration?
      call EIRENE_qromo(f,a,b,s,EIRENE_midpnt,epsr)
      return
      end subroutine EIRENE_aqc8


      SUBROUTINE EIRENE_QROMO(FUNC,A,B,SS,CHOOSE,epsr)
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
      INTEGER JMAX, JMAXP, KM, K
      PARAMETER (JMAX=14,JMAXP=JMAX+1,KM=4,K=KM+1)
      REAL(DP) A, B, DSS, SS, EPSR
      REAL(DP) S(JMAXP),H(JMAXP)
      INTEGER J
      REAL(DP) FUNC
      external choose,func,eirene_polint

      H(1)=1.0_DP
      DO 11 J=1,JMAX

        CALL CHOOSE(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL EIRENE_POLINT(H(J-KM),S(J-KM),K,0.0_DP,SS,DSS)
          IF (ABS(DSS).LT.EPSR*ABS(SS)) then
            RETURN
          endif
        ENDIF
        S(J+1)=S(J)
        H(J+1)=H(J)/9.0_DP
   11 CONTINUE
      write(iunout,*) '(W) Too many steps.'
      END SUBROUTINE EIRENE_QROMO

      SUBROUTINE EIRENE_POLINT(XA,YA,N,X,Y,DY)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER NMAX, N
      PARAMETER (NMAX=10)
      REAL(DP) DY, X, Y
      REAL(DP) XA(N),YA(N),C(NMAX),D(NMAX)
      INTEGER I, NS, M
      REAL(DP) DEN, DIF, DIFT, HO, HP, W

      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N

        DIFT=ABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
          NS=I

          DIF=DIFT

        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
   11 CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO=XA(I)-X
          HP=XA(I+M)-X
          W=C(I+1)-D(I)
          DEN=HO-HP
          DEN=W/DEN
          D(I)=HP*DEN
          C(I)=HO*DEN
   12   CONTINUE
        IF (2*NS.LT.N-M) THEN
          DY=C(NS+1)
        ELSE
          DY=D(NS)
          NS=NS-1
        ENDIF
        Y=Y+DY
   13 CONTINUE
      RETURN
      END SUBROUTINE EIRENE_POLINT
C
      SUBROUTINE EIRENE_MIDPNT(FUNC,A,B,S,N)
      USE EIRMOD_PRECISION
      IMPLICIT NONE
      INTEGER N, IT, J
      REAL(DP) A, B, S
      REAl(DP) DEL, DDEL, SUM, TNM, X
      REAL(DP) FUNC
      external func
      save

      IF (N.EQ.1) THEN
        S=(B-A)*FUNC(0.5_DP*(A+B))
        IT=1
      ELSE
        TNM=REAL(IT,DP)
        DEL=(B-A)/(3.0_DP*TNM)
        DDEL=DEL+DEL
        X=A+0.5_DP*DEL
        SUM=0.0_DP
        DO 11 J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DDEL
          SUM=SUM+FUNC(X)
          X=X+DEL
   11   CONTINUE
        S=(S+(B-A)*SUM/TNM)/3.0_DP
        IT=3*IT
      ENDIF
      RETURN
      END SUBROUTINE EIRENE_MIDPNT
