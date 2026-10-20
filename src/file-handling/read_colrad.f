cdr aug 19: remove argument iz1 (unused)
cdr sept 18:
cdr alternative reaction component pop. coeff: revised.
cdr tbd.: revise: form 2 vs. form 1, IFLAV and ICRM may be different now.


cdr  1) character(*) "filnam" is interpreted to identify
cdr               a particular "internal CR-model" ICRM.
cdr     character(*) "filnam" --> iflav (=1,2, or 4), iformul (=1,2)
cdr     e.g.  filnam='CR' (default, = 'CR_H_2'),  or filnam='CR_He_2',...
cdr  2) recognized data from internal CR-Models: crl(..)%ihsw
cdr                                              crl(..)%hsrt  (character(10)


       subroutine EIRENE_read_colrad (ir,filnam,reac,isw,
     .                                m_popesc, m_upper, m_lower,
     .                                irow_esc,icol_esc,pop_esc,
     .                                m_qext,
     .                                irow_ext,ital_ext)

cdr  Called from SLREAC.f in initialization phase

cdr  Purpose: prepare A&M data structure REACDAT(IR) for an internal, built-in,
cdr            collisional-radiative code:
cdr  1)  H_colrad,     (iflav=1)  (code: K. Sawada, adapted for eirene: D. Reiter)
cdr  2)  H2_colrad.... (iflav=2),  not ready
cdr  3)  Li_colrad.... (iflav=3),  not ready
cdr  4)  He_colrad,    (iflav=4)  (code: M. Goto, translated to fortran: P. Boerner
cdr                                      adapted for eirene: D. Reiter)
cdr
cdr  input:
c           ir:           internal reaction number on eirene structure REACDAT
c           filnam:       key for distinguishing internal collisional radiative model
c                         filnam=CR_A_I, A= H, HE H2, I=1,2
c                         if only filnam=CR    --> CH_H_2
c                         if only filnam=CR_A  --> CH_A_2
c           reac(10):         2.1.5... or 2.1.8... or 2.3.9.... or 2.3.13... (character(10))

c           isw:   =0     data for interaction potential                 (not in use)
c                  =1     data for collision cross-section               (not in use)
c                  =2-4   data for reaction rate coefficient             (only = 4  in use)
c                  =5-7   data for momentum-weighted rate coefficient    (not in use)
c                  =8-10  data for energy-weighted rate coefficient      (only = 10 in use)
c                  =11,12 other data, such as red. pop. coefficients     (only = 12 in use)
c
cdr:  for H: currently used only
cdr                       H.4, 2.1.5  and H.10, 2.1.5,  EI,  ionisation
cdr                       H.4, 2.1.PHs and H.10, 2.1.PHs, EI,  ionisation
cdr                       H.4. 2.1.8  and H.10, 2.1.8, RC   recombination
cdr                       H.4. 2.1.PHa and H.10, 2.1.PHa, RC   recombination
cdr                       and  H.11, H.12: selected population coefficients
cdr:  for He: currently used only
cdr                       H.4, 2.2 and H.10, 2.2, EI,  ionisation
cdr                       H.4. 2.3.9 and H.10, 2.3.9, RC   recombination
cdr                       and  H.11, H.12: selected population coefficients

cdr unfinished, but partly available:
cdr available :  iform=1 (MS resolved), l_ext=t  (external source of excited states,
cdr              other than ionising or recombining component
c
c  to be done: units, log-lin, scaling, asymptotics

      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout

      implicit none

      integer, intent(in) :: ir, isw
      character(len=*), intent(in) :: filnam, reac
c  optional input parameters
      integer, intent(in), optional :: m_popesc, m_upper, m_lower
      integer, intent(in), optional :: irow_esc(:), icol_esc(:)
      real(dp) , intent(in), optional :: pop_esc(:)
      integer, intent(in), optional :: m_qext
      integer, intent(in), optional :: irow_ext(:), ital_ext(:,:)

      integer, save :: ifirst = 0
      integer :: ivar, i, istor,
     .           iflav, iformul, icrm,
     .           indu, indu2, ival,
     .           nitems, nform
      character(:), allocatable :: strng
c
      type cr_liste
c  list of  data available from internal cr-code models ICRM
        integer :: items
        integer, allocatable :: ihsw(:)
        character(10), allocatable :: hstr(:)
      end type cr_liste
CDR  MAXCRM=4, currently. Should be set after reading block 4, as needed.
      type(cr_liste), save :: crl(maxcrm)

      external :: eirene_exit_own, eirene_leer

      close (29+ifoff)  ! nothing further to be read, currently

cdr  error exit for unfinished options
      if (isw.ne.4 .and. isw.ne.10 .and. isw.ne.12)  goto 1000

cdr  tbd: also exit unless "REAC" contains 2.1.5, OR 2.1.8,
cdr                                     OR 2.3.9, OR 2.3.13
cdr                                     OR 2.2,   OR 2.3.2
cdr
cdr  other reactions are not yet programmed, neither in xsectp,
cdr  nor in RATE_COEFF, nor in ENERGY_RATE_COEFF nor in OTHER_RATE_COEFF.

      if (ifirst == 0) then
        ifirst = 1

        crl%items = 0

c  list of identifiers for data available from intrinsic CR codes
c....................................................................
c
c     crl(1): H_COLRAD (atomic), formulation II only
c....................................................................

        call eirene_init_identifiers_H_colrad

c....................................................................
c
c  crl(4): HE_COLRAD (atomic), formulation II only
c.....................................................................

        call eirene_init_identifiers_HE_colrad


      end if

!pb check optional arguments
c     irow_esc = 0
c     if (present(ir_esc)) irow_esc = ir_esc
c     icol_esc = 0
c     if (present(ic_esc)) icol_esc = ic_esc
c     pop_esc = 1._dp
c     if (present(p_esc)) pop_esc = p_esc

!pb  identify, which collisional radiative model is to be used:
cdr  set variables: iflav, iformul, parsing from FILNAM = CR....

      strng = trim(filnam)
      indu = index(strng,'_')
      if (indu == 0) then
!  only key 'CR' found, default behaviour, use H coll. rad. model
!                       in MS unresolved mode F2
        IFLAV = 1
        iformul = 0
      else
        indu2 = index(strng(indu+1:),'_')
        if (indu2 == 0) then
!  no specification of MS condensation type (F1 or F2) available, use default: F2
          iformul = 0
        else
!  specification of MS condensation type (MS resolution) is available
          write (IUNOUT,*) 'READ_COLRAD ', STRNG, INDU,INDU2
          read (strng(indu+indu2+1:),*) ival
          if ((ival == 1) .or. (ival == 2)) then
            iformul = ival
          else
            write (iunout,*) ' ERROR READING REACTION IR = ', ir
            write (iunout,*) ' WRONG MS RESOLUTION TYPE FOUND '
            write (iunout,*) ' FORMTYPE = ', ival
            write (iunout,*) ' USE DEFAULT FORMULATION FOR THIS CRM '
            iformul = 0
          end if
        end if
!
        if (index(strng(indu:),'_HE') > 0) then
cdr  He_Colrad code:
          iflav = 4
          if (iformul == 0) iformul = 2  ! default, unless otherwise specified
          IF (iformul == 2) NFORM=1
          IF (iformul == 1) NFORM=3  ! (one ground and two MS states) 
        elseif (index(strng(indu:),'_H2') > 0) then
cdr  H2 cr-code: not ready
          iflav = 2
          if (iformul == 0) iformul = 2  ! default, unless otherwise specified
          IF (iformul == 2) NFORM=1
          IF (iformul == 1) NFORM=15  ! (one ground and 14 MS states H2(v)) 
        elseif (index(strng(indu:),'_H') > 0) then
cdr  H_colrad code:
          iflav = 1
          if (iformul == 0) iformul = 2  ! default, unless otherwise specified
            IF (iformul == 2) NFORM=1
            IF (iformul == 1) NFORM=2  ! (one ground and n=2 or 2s state)
        elseif (index(strng(indu:),'_LI') > 0) then
cdr  Li cr-code: not ready
          iflav = 3
          if (iformul == 0) iformul = 2  ! default, unless otherwise specified
          NFORM=1
        else
          write (iunout,*) ' ERROR READING REACTION IR = ', ir
          write (iunout,*) ' UNKNOW TYPE OF COLL. RAD. MODEL'
          write (iunout,*) ' KEY STRING : ',strng
          CALL EIRENE_EXIT_OWN(1)
        end if
      end if

cdr  now we have: iflav and iformul
cdr  find corresponding CR model: ICRM
cdr  (Petras code was : icrm=iflav, to be generalized for different IFORM
      ICRM=IFLAV  !...for the time being. Must be changed.
CDR  Next: interpret character string 'REAC', for CR model ICRM:  IFLAV/IFORMUL


cdr  IDENTIFY THE INDEX IVAR OF THE PRE-DEFINED VARIABLE HSTR(IVAR)
cdr  TO BE STORED ON M_HCOL(1:NHCOL_STORE).
      IVAR = 0
      DO I = 1, CRL(ICRM)%ITEMS
        IF (ISW /= CRL(ICRM)%IHSW(I)) CYCLE
        IF (REAC(1:10) == CRL(ICRM)%HSTR(I)(1:10)) THEN
             IVAR = I
             EXIT
        END IF
      END DO
      IF (IVAR == 0) GOTO 1000

!  CHECK IF VARIABLE HAS ALREADY BEEN MARKED FOR STORING EARLIER FOR THIS ICRM

      ISTOR = PCRM(ICRM)%P%NHCOL_STORE + 1
      DO I = 1, PCRM(ICRM)%P%NHCOL_STORE
        IF (IVAR == PCRM(ICRM)%P%M_HCOL(I)) THEN
          ISTOR = I
          EXIT
        END IF
      END DO
!  VARIABLE NOT YET MARKED FOR STORING --> MARK NOW
      IF (ISTOR > PCRM(ICRM)%P%NHCOL_STORE) THEN
        PCRM(ICRM)%P%NHCOL_STORE = ISTOR
        PCRM(ICRM)%P%M_HCOL(PCRM(ICRM)%P%NHCOL_STORE) = IVAR
      END IF

      SELECT CASE (ISW)
        CASE (2:4)
          IF (REACDAT(IR)%LRTC) THEN
            WRITE (IUNOUT,*) ' RATE COEFFICIENT ALREADY SPECIFIED',
     .                       ' FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTC)

          REACDAT(IR)%LRTC = .TRUE.
          REACDAT(IR)%RTC%IFIT = 5

          CALL EIRENE_STORE_CRM_PARAMETERS(REACDAT(IR)%RTC)

        CASE (5:7)
          IF (REACDAT(IR)%LRTCMW) THEN
            WRITE (IUNOUT,*) ' MOMENTUM-WEIGHTED RATE COEFFICIENT',
     .                       ' ALREADY SPECIFIED FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCMW)

          REACDAT(IR)%LRTCMW = .TRUE.
          REACDAT(IR)%RTCMW%IFIT = 5

          CALL EIRENE_STORE_CRM_PARAMETERS(REACDAT(IR)%RTCMW)

        CASE (8:10)
          IF (REACDAT(IR)%LRTCEW) THEN
            WRITE (IUNOUT,*) ' ENERGY-WEIGHTED RATE COEFFICIENT',
     .                       ' ALREADY SPECIFIED FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%RTCEW)

          REACDAT(IR)%LRTCEW = .TRUE.
          REACDAT(IR)%RTCEW%IFIT = 5

          CALL EIRENE_STORE_CRM_PARAMETERS(REACDAT(IR)%RTCEW)

        CASE (11:12)
          IF (REACDAT(IR)%LOTH) THEN
            WRITE (IUNOUT,*) ' OTHER RATE COEFFICIENT',
     .                       ' ALREADY SPECIFIED FOR REACTION', IR
            WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
            CALL EIRENE_EXIT_OWN(1)
          END IF

          CALL EIRENE_ALLOC_FIT_FORM (REACDAT(IR)%OTH)

          REACDAT(IR)%LOTH = .TRUE.
          REACDAT(IR)%OTH%IFIT = 5

          CALL EIRENE_STORE_CRM_PARAMETERS(REACDAT(IR)%OTH)

        CASE DEFAULT
          GOTO 1000
        END SELECT
      RETURN

 1000 continue
      CALL EIRENE_LEER(1)
      WRITE (IUNOUT,*) ' ERROR IN "READ_COLRAD" : '
      WRITE (IUNOUT,*) ' WRONG DATA TYPE FOR INTERNAL COLRAD OPTION'
      WRITE (IUNOUT,*) ' REACTION NO. ', IR
      WRITE (IUNOUT,'(1X,A,I0)') ' DATA TYPE H.', ISW
      WRITE (IUNOUT,'(1X,A,A10)') ' DATA NR.    ', REAC(1:10)
      CALL EIRENE_EXIT_OWN(1)
      RETURN

      CONTAINS


!***********************************************************************


      SUBROUTINE EIRENE_STORE_CRM_PARAMETERS (RP)
cdr  RP  is: REACDAT(KK)%A, with A = RTC, RTCMW, RTCEW or OTH
cdr      for reaction card no. KK 

      TYPE(FIT_FORMS),POINTER :: RP   


      ALLOCATE (RP%CRM)

      RP%CRM%ICRM = ICRM      
      RP%CRM%IFLAV = IFLAV
      RP%CRM%IFORMUL = IFORMUL

      RP%CRM%IVARST = ISTOR  ! Number of stored variables (per cell) for ICRM
      RP%CRM%M_POPESC = M_POPESC
cdr  tbd.: next two parameters: if both are .ne. 0:
      RP%CRM%M_UPPER = M_UPPER !dr turn pop. coeff into emissivity
      RP%CRM%M_LOWER = M_LOWER !dr turn pop. coeff into emissivity
      RP%CRM%M_QEXT   = M_QEXT   ! external sources (photo excitation) for model ICRM
      
      IF (M_POPESC > 0) THEN
cdr optional CRM parameters: population escape factors
         ALLOCATE(RP%CRM%IROW_ESC(M_POPESC))
         ALLOCATE(RP%CRM%ICOL_ESC(M_POPESC))
         ALLOCATE(RP%CRM%POP_ESC(M_POPESC))
         RP%CRM%IROW_ESC = IROW_ESC(1:M_POPESC)
         RP%CRM%ICOL_ESC = ICOL_ESC(1:M_POPESC)
         RP%CRM%POP_ESC  = POP_ESC(1:M_POPESC)
      END IF

      IF (M_QEXT > 0) THEN
cdr optional CRM parameters: population escape factors
         ALLOCATE(RP%CRM%IROW_EXT(M_QEXT))       !  position of external source rate in q_ext array
         ALLOCATE(RP%CRM%ITAL_EXT(3,1:M_QEXT))   !  ispz ,ital, and ispza (absorber) for external source rate
         RP%CRM%IROW_EXT    = IROW_EXT(1:M_QEXT) 
         RP%CRM%ITAL_EXT    = ITAL_EXT(1:3,1:M_QEXT)
      END IF


      END SUBROUTINE EIRENE_STORE_CRM_PARAMETERS


!***********************************************************************


      SUBROUTINE EIRENE_INIT_IDENTIFIERS_H_COLRAD
cdr  this coding is very inflexible. It is essentially impossible to
cdr  add further parameters in a logical way, e.g. MS resolved data (form-I)
cdr  For MS resolved data, hstr(1) --> hstr(1,nms,nms)
cdr                        hstr(3) --> hstr(3,nms)
cdr mar 2023: added alp_ext, ealp_ext,with label 109 and 110, while 
cdr           logically it should come after labels 5,6

c....................................................................
c
cdr   crl(1): H_COLRAD (atomic), formulation II only
c.....................................................................

      nitems=110
      crl(1)%items = nitems     !  6 + 5 + 5 + 5 + 29 + 29 + 29

      allocate (crl(1)%ihsw(nitems)) ! values: 4, 10, 12
      allocate (crl(1)%hstr(nitems)) ! character(10)

c   cr rate, H.4, coupling to ground state H(1), ionisation
      crl(1)%ihsw(1) = 4
      crl(1)%hstr(1) = '2.1.5     '
c  electron cooling rate coeff. H.10, e_scr is negative from h-colrad
c  note: with delpot=-13.6 (input): this becomes the radiation loss rate coeff. alone
      crl(1)%ihsw(2) = 10
      crl(1)%hstr(2) = '2.1.5     '

c   cr rate, H.4, coupling to continuum H+, recombination
      crl(1)%ihsw(3) = 4
      crl(1)%hstr(3) = '2.1.8     '
c  electron cooling/heating rate coeff. H.10, (both signs possible. loss: negative e_alpcr))
c  note: with delpot=+13.6 (input): this becomes the radiation loss rate coeff. alone
      crl(1)%ihsw(4) = 10
      crl(1)%hstr(4) = '2.1.8     '
cdr
c   cr rate and e_rate, ionising coupling to external sources, e.g. molecules,
cdr                     radiation field, photo-excitation. Unfinished
      crl(1)%ihsw(5) = 4
      crl(1)%hstr(5) = '2.1.PHs   '     ! also '2.1.EXs1', '2.1.EXs2',...
      crl(1)%ihsw(6) = 10
      crl(1)%hstr(6) = '2.1.PHs   '

cdr  missing:  7,8:
c   cr rate and e_rate, recomb. coupling to external sources, e.g. molecules,
cdr                     radiation field, photo-excitation. Unfinished
cdr mar. 2023: Now added. 

c  reduced population coefficient, H(n=3,2,4,5,6) states,
c  AMJUEL notation
c  no.  7-11:  H*(n) states, n=2,...6
c  COMPONENT: coupling to ground state H(1), ionisation, numbering as in AMJUEL file.
      crl(1)%ihsw(7) = 12
      crl(1)%hstr(7) = '2.1.5a    ' ! n=3
      crl(1)%ihsw(8) = 12
      crl(1)%hstr(8) = '2.1.5b    ' ! n=2
      crl(1)%ihsw(9) = 12
      crl(1)%hstr(9) = '2.1.5c    ' ! n=4
      crl(1)%ihsw(10) = 12
      crl(1)%hstr(10) = '2.1.5d    ' ! n=5
      crl(1)%ihsw(11) = 12
      crl(1)%hstr(11) = '2.1.5e    ' ! n=6

c  reduced population coefficient, H(n=3,2,4,5,6) states,
c  COMPONENT: coupling to H+, recombination, numbering as in AMJUEL file.
      crl(1)%ihsw(12) = 12
      crl(1)%hstr(12) = '2.1.8a    ' ! n=3
      crl(1)%ihsw(13) = 12
      crl(1)%hstr(13) = '2.1.8b    ' ! n=2
      crl(1)%ihsw(14) = 12
      crl(1)%hstr(14) = '2.1.8c    ' ! n=4
      crl(1)%ihsw(15) = 12
      crl(1)%hstr(15) = '2.1.8d    ' ! n=5
      crl(1)%ihsw(16) = 12
      crl(1)%hstr(16) = '2.1.8e    ' ! n=6

c  reduced population coefficient, H(n=3,2,4,5,6) states,
c  COMPONENT: coupling to radiation field (photo-excitation)
      crl(1)%ihsw(17) = 12
      crl(1)%hstr(17) = '2.1.PHa  ' ! n=3
      crl(1)%ihsw(18) = 12
      crl(1)%hstr(18) = '2.1.PHb  ' ! n=2
      crl(1)%ihsw(19) = 12
      crl(1)%hstr(19) = '2.1.PHc  ' ! n=4
      crl(1)%ihsw(20) = 12
      crl(1)%hstr(20) = '2.1.PHd  ' ! n=5
      crl(1)%ihsw(21) = 12
      crl(1)%hstr(21) = '2.1.PHe  ' ! n=6


cdr: here to be done: further components: coupling to H2, H2+, H+, H3+

!  Coupling to H ground state
!  alternative naming convention for excited states, n=2,...30
cdr  duplicating storage needs for i=2,..6?
cdr   H.12  2.1.5.i,   for i=2,30
      crl(1)%ihsw(22:50) = 12
      do i=2, 30
        crl(1)%hstr(21+i-1) = '2.1.5.    '
        write(crl(1)%hstr(21+i-1)(7:8),'(i0)') i
      end do

!  Coupling to H+ ion state
!  alternative naming convention for excited states, n=2,...30
cdr   H.12  2.1.8.i,   for i=2,30
      crl(1)%ihsw(51:79) = 12
      do i=2, 30
        crl(1)%hstr(50+i-1) = '2.1.8.    '
        write(crl(1)%hstr(50+i-1)(7:8),'(i0)') i
      end do

!  Coupling to ext. source of extited states, e.g. photo-excitation
!  alternative naming convention for excited states, n=2,...30
cdr   H.12  2.1.5PH.i,   for i=2,30
      crl(1)%ihsw(80:108) = 12
      do i=2, 30
        crl(1)%hstr(79+i-1) = '2.1.PH    '
        write(crl(1)%hstr(79+i-1)(7:8),'(i0)') i
      end do
cdr added in march 2023
      crl(1)%ihsw(109) = 4
      crl(1)%hstr(109) = '2.1.PHa   '
      crl(1)%ihsw(110) = 10
      crl(1)%hstr(110) = '2.1.PHa   '
      
      END SUBROUTINE EIRENE_INIT_IDENTIFIERS_H_COLRAD


 !***********************************************************************
     
      
      SUBROUTINE EIRENE_INIT_IDENTIFIERS_HE_COLRAD

c....................................................................
c
cdr   crl(4): HE_COLRAD (atomic), formulation II only
c.....................................................................

      crl(4)%items = 195
      allocate (crl(4)%ihsw(crl(4)%items))
      allocate (crl(4)%hstr(crl(4)%items))

c  effective ionisation rate, formulation II, MS condensed
      crl(4)%ihsw(1) = 4
      crl(4)%hstr(1) = '2.3.9a  '
c  electron cooling rate coeff. e_scr is negative from h-colrad
c  note: with delpot=-24.588 (input): this becomes the radiation loss rate coeff. alone
      crl(4)%ihsw(2) = 10
      crl(4)%hstr(2) = '2.3.9a  '
c  effective recombination rate
      crl(4)%ihsw(3) = 4
      crl(4)%hstr(3) = '2.3.13a '
c  electron cooling/heating rate coeff. (both signs possible. loss: negative e_alpcr))
c  note: with delpot=+24.588 (input): this becomes the radiation loss rate coeff. alone
      crl(4)%ihsw(4) = 10
      crl(4)%hstr(4) = '2.3.13a '
c  external source driven ionisation rate, e.g. photo-excitation driven ionisation
      crl(4)%ihsw(5) = 4
      crl(4)%hstr(5) = '2.3.9PH '
      crl(4)%ihsw(6) = 10
      crl(4)%hstr(6) = '2.3.9PH '
C
c  population coefficients, FII, coupling to ground state He(1) atom
c  AMJUEL notation
c  Strangely, here the AMJUEL reaction notation is that from B (Strahl, Behringer)
c  and C (Colrad, McWhirter):  2.2 and 2.3.2
      crl(4)%ihsw(7) = 12
      crl(4)%hstr(7) = '2.2a    ' ! I=6
      crl(4)%ihsw(8) = 12
      crl(4)%hstr(8) = '2.2b    ' ! I=7
      crl(4)%ihsw(9) = 12
      crl(4)%hstr(9) = '2.2c    ' ! I=8
      crl(4)%ihsw(10) = 12
      crl(4)%hstr(10) = '2.2d    ' ! I=10
      crl(4)%ihsw(11) = 12
      crl(4)%hstr(11) = '2.2e    ' ! I=16


c  population coefficients, FII, coupling to He+ ion
      crl(4)%ihsw(12) = 12
      crl(4)%hstr(12) = '2.3.2a  ' ! I=6
      crl(4)%ihsw(13) = 12
      crl(4)%hstr(13) = '2.3.2b  ' ! I=7
      crl(4)%ihsw(14) = 12
      crl(4)%hstr(14) = '2.3.2c  ' ! I=8
      crl(4)%ihsw(15) = 12
      crl(4)%hstr(15) = '2.3.2d  ' ! I=10
      crl(4)%ihsw(16) = 12
      crl(4)%hstr(16) = '2.3.2e  ' ! I=16



c  population coefficients, coupling to external source of excitation (e.g. photons)
c  only available if L_EXT=.TRUE. in call to HE_COLRAD
      crl(4)%ihsw(17) = 12
      crl(4)%hstr(17) = '2.2PHa  ' ! I=6
      crl(4)%ihsw(18) = 12
      crl(4)%hstr(18) = '2.2PHb  ' ! I=7
      crl(4)%ihsw(19) = 12
      crl(4)%hstr(19) = '2.2PHc  ' ! I=8
      crl(4)%ihsw(20) = 12
      crl(4)%hstr(20) = '2.2PHd  ' ! I=10
      crl(4)%ihsw(21) = 12
      crl(4)%hstr(21) = '2.2PHe  ' ! I=16
      

!  alternative naming convention,
!  He data numbering from HYDHEL: reactions 2.3.9.i and 2.3.13.i, i=2,59
      crl(4)%ihsw(22:79) = 12
      do i=2, 59
        crl(4)%hstr(21+i-1) = '2.3.9.    '
        write(crl(4)%hstr(21+i-1)(7:8),'(i0)') i
      end do

!  alternative naming convention
      crl(4)%ihsw(80:137) = 12
      do i=2, 59
        crl(4)%hstr(79+i-1) = '2.3.13.   '
        write(crl(4)%hstr(79+i-1)(8:9),'(i0)') i
      end do

!  alternative naming convention
      crl(4)%ihsw(138:195) = 12
      do i=2, 59
        crl(4)%hstr(137+i-1) = '2.3.9PH.  '
        write(crl(4)%hstr(137+i-1)(9:10),'(i0)') i
      end do
      
      END SUBROUTINE EIRENE_INIT_IDENTIFIERS_HE_COLRAD

      end subroutine EIRENE_read_colrad
