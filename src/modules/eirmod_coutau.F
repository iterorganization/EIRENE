c nov. 2015:  species index ipls added for energy-pl tallies:
c             eapli,empli,eipli,ephpli,eppli
cdr dec. 15:  comments added. missing tallies ppeli, epeli, etc..??
cpb jan. 18:  array NFSTPI moved to module EIRMOD_COMUSR
cdr jan. 22:  vectorial MXPL_VEC tallies
cdr mar. 23:  fetch_outau --> fetch_outaui renamed, returns selected integrated tallies
cdr           fetch_outau added. returns selected tallies

      MODULE EIRMOD_COUTAU
cdr  mostly global tallies integrated from more deeply
cdr  resolved tallies,
cdr  integration done during post MC-loop processing.
cdr  Few tallies are directly scored: ptrash, etrash, .... xmcp

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CLOGAU, ONLY: NLSPCSCL_ATM, NLSPCSCL_MOL,
     .                         NLSPCSCL_ION, NLSPCSCL_PHOT

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_COUTAU, EIRENE_DEALLOC_COUTAU,
     .          EIRENE_INIT_COUTAU, 
     .          EIRENE_FETCH_OUTAU, EIRENE_FETCH_OUTAUI,
     .          EIRENE_WRITE_COUTAU, EIRENE_READ_COUTAU,
     .          EIRENE_INIT_COUTAU_REINIT,
     .          EIRENE_BROADCAST_COUTAU

! INTEGRALS OF VOLUME-AVERAGED TALLIES
      REAL(DP), PUBLIC, ALLOCATABLE, TARGET, SAVE ::
     R PDENAI(:,:), PDENMI(:,:), PDENII(:,:), PDENPHI(:,:),
     R EDENAI(:,:), EDENMI(:,:), EDENII(:,:), EDENPHI(:,:),
     R PAELI(:),    PAATI(:,:),  PAMLI(:,:),  PAIOI(:,:),  PAPHTI(:,:),
     R PAPLI(:,:),
     R PMELI(:),    PMATI(:,:),  PMMLI(:,:),  PMIOI(:,:),  PMPHTI(:,:),
     R PMPLI(:,:),
     R PIELI(:),    PIATI(:,:),  PIMLI(:,:),  PIIOI(:,:),  PIPHTI(:,:),
     R PIPLI(:,:),
     R PPHELI(:),   PPHATI(:,:), PPHMLI(:,:), PPHIOI(:,:), PPHPHTI(:,:),
     R PPHPLI(:,:),
     R EAELI(:),    EAATI(:),    EAMLI(:),  EAIOI(:),  EAPHTI(:),
     R EAPLI(:,:),
     R EMELI(:),    EMATI(:),    EMMLI(:),  EMIOI(:),  EMPHTI(:),
     R EMPLI(:,:),
     R EIELI(:),    EIATI(:),    EIMLI(:),  EIIOI(:),  EIPHTI(:),
     R EIPLI(:,:),
     R EPHELI(:),   EPHATI(:),   EPHMLI(:), EPHIOI(:), EPHPHTI(:),
     R EPHPLI(:,:),
     R RAELI(:,:),  RMELI(:,:),  RIELI(:,:),
     R ADDVI(:,:),
     R COLVI(:,:),  SNAPVI(:,:), COPVI(:,:),  BGKVI(:,:), ALGVI(:,:),
     R PGENAI(:,:), PGENMI(:,:), PGENII(:,:), PGENPHI(:,:),
     R EGENAI(:,:), EGENMI(:,:), EGENII(:,:), EGENPHI(:,:),
     R VGENAI(:,:), VGENMI(:,:), VGENII(:,:), VGENPHI(:,:),
     R PPATI(:,:),  PPMLI(:,:),  PPIOI(:,:),  PPPHTI(:,:), PPPLI(:,:),
     R EPATI(:),    EPMLI(:),    EPIOI(:),    EPPHTI(:),   EPPLI(:,:),
     R VXDENAI(:,:), VXDENMI(:,:), VXDENII(:,:), VXDENPHI(:,:),
     R VYDENAI(:,:), VYDENMI(:,:), VYDENII(:,:), VYDENPHI(:,:),
     R VZDENAI(:,:), VZDENMI(:,:), VZDENII(:,:), VZDENPHI(:,:),
     R MAPLI_VEC(:,:),  MMPLI_VEC(:,:),  MIPLI_VEC(:,:),
     R MPHPLI_VEC(:,:)

! INTEGRALS OF SURFACE TALLIES: PARTICLE FLUXES
      REAL(DP), PUBLIC, ALLOCATABLE, TARGET, SAVE ::
     R POTATI(:,:), PRFAAI(:,:), PRFMAI(:,:), PRFIAI(:,:), PRFPHAI(:,:),
     R PRFPAI(:,:),
     R POTMLI(:,:), PRFAMI(:,:), PRFMMI(:,:), PRFIMI(:,:), PRFPHMI(:,:),
     R PRFPMI(:,:),
     R POTIOI(:,:), PRFAII(:,:), PRFMII(:,:), PRFIII(:,:), PRFPHII(:,:),
     R PRFPII(:,:),
     R POTPHTI(:,:), PRFAPHTI(:,:), PRFMPHTI(:,:), PRFIPHTI(:,:),
     R PRFPHPHTI(:,:), PRFPPHTI(:,:),
     R POTPLI(:,:)
! INTEGRALS OF SURFACE TALLIES: ENERGY FLUXES
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R EOTATI(:,:), ERFAAI(:,:), ERFMAI(:,:), ERFIAI(:,:), ERFPHAI(:,:),
     R ERFPAI(:,:),
     R EOTMLI(:,:), ERFAMI(:,:), ERFMMI(:,:), ERFIMI(:,:), ERFPHMI(:,:),
     R ERFPMI(:,:),
     R EOTIOI(:,:), ERFAII(:,:), ERFMII(:,:), ERFIII(:,:), ERFPHII(:,:),
     R ERFPII(:,:),
     R EOTPHTI(:,:), ERFAPHTI(:,:), ERFMPHTI(:,:), ERFIPHTI(:,:),
     R ERFPHPHTI(:,:), ERFPPHTI(:,:),
     R EOTPLI(:,:)
! INTEGRALS OF SURFACE TALLIES: SPUTTERING
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R SPTAATI(:,:), SPTMATI(:,:), SPTIATI(:,:), SPTPHATI(:,:),
     R SPTPATI(:,:),
     R SPTAMLI(:,:), SPTMMLI(:,:), SPTIMLI(:,:), SPTPHMLI(:,:),
     R SPTPMLI(:,:),
     R SPTAIOI(:,:), SPTMIOI(:,:), SPTIIOI(:,:), SPTPHIOI(:,:),
     R SPTPIOI(:,:),
     R SPTAPHTI(:,:), SPTMPHTI(:,:), SPTIPHTI(:,:), SPTPHPHTI(:,:),
     R SPTPPHTI(:,:),
     R SPTAPLI(:,:), SPTMPLI(:,:), SPTIPLI(:,:), SPTPHPLI(:,:),
     R SPTPPLI(:,:),
     R sptatti(:), sptmtti(:), sptitti(:), sptphtti(:), sptpltti(:),
     R SPTTTI(:),
! INTEGRALS OF ADDITIONAL TALLIES
     R ADDSI(:,:),  ALGSI(:,:),
     R SPUMPI(:,:)

! INTEGRAL VALUES, global balances, scaling, from primary sources
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R WTOTA(:,:),  WTOTM(:,:),  WTOTI(:,:),  WTOTP(:,:),  WTOTPH(:,:),
     R WTOTE(:),
     R ETOTA(:),    ETOTM(:),    ETOTI(:),    ETOTP(:),    ETOTPH(:),
cdr scored along the flights
     R XMCP(:),     EELFI(:,:),
     R PTRASH(:),   ETRASH(:),

cdr for scaling
     R FLUXT(:),    FLXFAC(:),
     R FASCL(:,:),  FMSCL(:,:),  FISCL(:,:),  FPHSCL(:,:)

! helper pointers for species resolved integrals
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R PAATI2(:,:), PAMLI2(:,:), PAIOI2(:,:), PAPHTI2(:,:), PAPLI2(:,:),
     R PMATI2(:,:), PMMLI2(:,:), PMIOI2(:,:), PMPHTI2(:,:), PMPLI2(:,:),
     R PIATI2(:,:), PIMLI2(:,:), PIIOI2(:,:), PIPHTI2(:,:), PIPLI2(:,:),
     R PPHATI2(:,:), PPHMLI2(:,:), PPHIOI2(:,:), PPHPHTI2(:,:),
     R PPHPLI2(:,:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R PRFAAI2(:,:), PRFMAI2(:,:), PRFIAI2(:,:), PRFPHAI2(:,:),
     R PRFAMI2(:,:), PRFMMI2(:,:), PRFIMI2(:,:), PRFPHMI2(:,:),
     R PRFAII2(:,:), PRFMII2(:,:), PRFIII2(:,:), PRFPHII2(:,:),
     R PRFAPHTI2(:,:), PRFMPHTI2(:,:), PRFIPHTI2(:,:),
     R PRFPHPHTI2(:,:)

#if ( defined(USE_OPENMP) || defined(USE_EXT_OPENMP) )
!$OMP  THREADPRIVATE(PAATI2, PAMLI2, PAIOI2, PAPHTI2, PAPLI2,
!$OMP& PMATI2, PMMLI2, PMIOI2, PMPHTI2, PMPLI2,
!$OMP& PIATI2, PIMLI2, PIIOI2, PIPHTI2, PIPLI2,
!$OMP& PPHATI2, PPHMLI2, PPHIOI2, PPHPHTI2, PPHPLI2,
!$OMP& PRFAAI2, PRFMAI2, PRFIAI2, PRFPHAI2,
!$OMP& PRFAMI2, PRFMMI2, PRFIMI2, PRFPHMI2,
!$OMP& PRFAII2, PRFMII2, PRFIII2, PRFPHII2,
!$OMP& PRFAPHTI2, PRFMPHTI2, PRFIPHTI2, PRFPHPHTI2)
#endif

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NADDI(:),  NFRSTI(:), NDDWI(:),  NFRTWI(:),
     I NFSTVI(:), NFSTWI(:),
     I NEXTVI(:), NEXTWI(:)

      INTEGER, PUBLIC, SAVE ::
     I NOUTA1, NOUTA2, NOUTAS, NOUTAU, NOUTTL
      INTEGER, SAVE :: IFRST=0

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_COUTAU
      LOGICAL :: LOGHELP(NSTRA)

      IF (ALLOCATED(PDENAI)) RETURN
cdr integrated volume tallies
      NOUTA1 = NVLTLP*NSTRAP
cdr integrated surface tallies
      NOUTA2 = NSFTLP*NSTRAP
      NOUTAS = (2*NPHOTP+3*NATMP+3*NMOLP+1*NPLSP+4*NIONP+11)*NSTRAP
      NOUTAU = NOUTA1+NOUTA2+NOUTAS
      NOUTTL = 4*(NTALV+NTALS)

      ALLOCATE (PDENAI(0:NATM,0:NSTRA))
      ALLOCATE (PDENMI(0:NMOL,0:NSTRA))
      ALLOCATE (PDENII(0:NION,0:NSTRA))
      ALLOCATE (PDENPHI(0:NPHOT,0:NSTRA))

      ALLOCATE (EDENAI(0:NATM,0:NSTRA))
      ALLOCATE (EDENMI(0:NMOL,0:NSTRA))
      ALLOCATE (EDENII(0:NION,0:NSTRA))
      ALLOCATE (EDENPHI(0:NPHOT,0:NSTRA))
cdr  volumetric particles source tallies, from incident atoms, sources for e,a,m,i,ph,pl
      ALLOCATE (PAELI(0:NSTRA))
      IF (NLSPCSCL_ATM) THEN
        ALLOCATE (PAATI(0:NATMP*NATMP-1,0:NSTRA))
        ALLOCATE (PAMLI(0:NMOLP*NATMP-1,0:NSTRA))
        ALLOCATE (PAIOI(0:NIONP*NATMP-1,0:NSTRA))
        ALLOCATE (PAPHTI(0:NPHOTP*NATMP-1,0:NSTRA))
        ALLOCATE (PAPLI(0:NPLSP*NATMP-1,0:NSTRA))
      ELSE
        ALLOCATE (PAATI(0:NATM,0:NSTRA))
        ALLOCATE (PAMLI(0:NMOL,0:NSTRA))
        ALLOCATE (PAIOI(0:NION,0:NSTRA))
        ALLOCATE (PAPHTI(0:NPHOT,0:NSTRA))
        ALLOCATE (PAPLI(0:NPLS,0:NSTRA))
      END IF
cdr  volumetric particles source tallies, from incident molec., sources for e,a,m,i,ph,pl
      ALLOCATE (PMELI(0:NSTRA))
      IF (NLSPCSCL_MOL) THEN
        ALLOCATE (PMATI(0:NATMP*NMOLP-1,0:NSTRA))
        ALLOCATE (PMMLI(0:NMOLP*NMOLP-1,0:NSTRA))
        ALLOCATE (PMIOI(0:NIONP*NMOLP-1,0:NSTRA))
        ALLOCATE (PMPHTI(0:NPHOTP*NMOLP-1,0:NSTRA))
        ALLOCATE (PMPLI(0:NPLSP*NMOLP-1,0:NSTRA))
      ELSE
        ALLOCATE (PMATI(0:NATM,0:NSTRA))
        ALLOCATE (PMMLI(0:NMOL,0:NSTRA))
        ALLOCATE (PMIOI(0:NION,0:NSTRA))
        ALLOCATE (PMPHTI(0:NPHOT,0:NSTRA))
        ALLOCATE (PMPLI(0:NPLS,0:NSTRA))
      END IF
cdr  volumetric particles source tallies, from incident test ions, sources for e,a,m,i,ph,pl
      ALLOCATE (PIELI(0:NSTRA))
      IF (NLSPCSCL_ION) THEN
        ALLOCATE (PIATI(0:NATMP*NIONP-1,0:NSTRA))
        ALLOCATE (PIMLI(0:NMOLP*NIONP-1,0:NSTRA))
        ALLOCATE (PIIOI(0:NIONP*NIONP-1,0:NSTRA))
        ALLOCATE (PIPHTI(0:NPHOTP*NIONP-1,0:NSTRA))
        ALLOCATE (PIPLI(0:NPLSP*NIONP-1,0:NSTRA))
      ELSE
        ALLOCATE (PIATI(0:NATM,0:NSTRA))
        ALLOCATE (PIMLI(0:NMOL,0:NSTRA))
        ALLOCATE (PIIOI(0:NION,0:NSTRA))
        ALLOCATE (PIPHTI(0:NPHOT,0:NSTRA))
        ALLOCATE (PIPLI(0:NPLS,0:NSTRA))
      END IF
cdr  volumetric particles source tallies, from incident photons, sources for e,a,m,i,ph,pl
      ALLOCATE (PPHELI(0:NSTRA))
      IF (NLSPCSCL_PHOT) THEN
        ALLOCATE (PPHATI(0:NATMP*NPHOTP-1,0:NSTRA))
        ALLOCATE (PPHMLI(0:NMOLP*NPHOTP-1,0:NSTRA))
        ALLOCATE (PPHIOI(0:NIONP*NPHOTP-1,0:NSTRA))
        ALLOCATE (PPHPHTI(0:NPHOTP*NPHOTP-1,0:NSTRA))
        ALLOCATE (PPHPLI(0:NPLSP*NPHOTP-1,0:NSTRA))
      ELSE
        ALLOCATE (PPHATI(0:NATM,0:NSTRA))
        ALLOCATE (PPHMLI(0:NMOL,0:NSTRA))
        ALLOCATE (PPHIOI(0:NION,0:NSTRA))
        ALLOCATE (PPHPHTI(0:NPHOT,0:NSTRA))
        ALLOCATE (PPHPLI(0:NPLS,0:NSTRA))
      END IF

      ALLOCATE (EAELI(0:NSTRA))
      ALLOCATE (EAATI(0:NSTRA))
      ALLOCATE (EAMLI(0:NSTRA))
      ALLOCATE (EAIOI(0:NSTRA))
      ALLOCATE (EAPHTI(0:NSTRA))
      ALLOCATE (EAPLI(0:NPLS,0:NSTRA))

      ALLOCATE (EMELI(0:NSTRA))
      ALLOCATE (EMATI(0:NSTRA))
      ALLOCATE (EMMLI(0:NSTRA))
      ALLOCATE (EMIOI(0:NSTRA))
      ALLOCATE (EMPHTI(0:NSTRA))
      ALLOCATE (EMPLI(0:NPLS,0:NSTRA))

      ALLOCATE (EIELI(0:NSTRA))
      ALLOCATE (EIATI(0:NSTRA))
      ALLOCATE (EIMLI(0:NSTRA))
      ALLOCATE (EIIOI(0:NSTRA))
      ALLOCATE (EIPHTI(0:NSTRA))
      ALLOCATE (EIPLI(0:NPLS,0:NSTRA))

      ALLOCATE (EPHELI(0:NSTRA))
      ALLOCATE (EPHATI(0:NSTRA))
      ALLOCATE (EPHMLI(0:NSTRA))
      ALLOCATE (EPHIOI(0:NSTRA))
      ALLOCATE (EPHPHTI(0:NSTRA))
      ALLOCATE (EPHPLI(0:NPLS,0:NSTRA))

      ALLOCATE (ADDVI(0:NADV,0:NSTRA))
      ALLOCATE (COLVI(0:NCLV,0:NSTRA))
      ALLOCATE (SNAPVI(0:NSNV,0:NSTRA))
      ALLOCATE (COPVI(0:NCPV,0:NSTRA))
      ALLOCATE (BGKVI(0:NBGV,0:NSTRA))
      ALLOCATE (ALGVI(0:NALV,0:NSTRA))

      ALLOCATE (PGENAI(0:NATM,0:NSTRA))
      ALLOCATE (PGENMI(0:NMOL,0:NSTRA))
      ALLOCATE (PGENII(0:NION,0:NSTRA))
      ALLOCATE (PGENPHI(0:NPHOT,0:NSTRA))

      ALLOCATE (EGENAI(0:NATM,0:NSTRA))
      ALLOCATE (EGENMI(0:NMOL,0:NSTRA))
      ALLOCATE (EGENII(0:NION,0:NSTRA))
      ALLOCATE (EGENPHI(0:NPHOT,0:NSTRA))

      ALLOCATE (VGENAI(0:NATM,0:NSTRA))
      ALLOCATE (VGENMI(0:NMOL,0:NSTRA))
      ALLOCATE (VGENII(0:NION,0:NSTRA))
      ALLOCATE (VGENPHI(0:NPHOT,0:NSTRA))

cdr  volumetric particles source tallies, from bulk ions, sources for e,a,m,i,ph,pl
cdr             tally for electron source missing ??
      ALLOCATE (PPATI(0:NATM,0:NSTRA))
      ALLOCATE (PPMLI(0:NMOL,0:NSTRA))
      ALLOCATE (PPIOI(0:NION,0:NSTRA))
      ALLOCATE (PPPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (PPPLI(0:NPLS,0:NSTRA))

cdr  volumetric energy source tallies, from bulk ions, sources for e,a,m,i,ph,pl
cdr             tally for electron source missing ??
      ALLOCATE (EPATI(0:NSTRA))
      ALLOCATE (EPMLI(0:NSTRA))
      ALLOCATE (EPIOI(0:NSTRA))
      ALLOCATE (EPPHTI(0:NSTRA))
      ALLOCATE (EPPLI(0:NPLS,0:NSTRA))

      ALLOCATE (VXDENAI(0:NATM,0:NSTRA))
      ALLOCATE (VXDENMI(0:NMOL,0:NSTRA))
      ALLOCATE (VXDENII(0:NION,0:NSTRA))
      ALLOCATE (VXDENPHI(0:NPHOT,0:NSTRA))
      ALLOCATE (VYDENAI(0:NATM,0:NSTRA))
      ALLOCATE (VYDENMI(0:NMOL,0:NSTRA))
      ALLOCATE (VYDENII(0:NION,0:NSTRA))
      ALLOCATE (VYDENPHI(0:NPHOT,0:NSTRA))
      ALLOCATE (VZDENAI(0:NATM,0:NSTRA))
      ALLOCATE (VZDENMI(0:NMOL,0:NSTRA))
      ALLOCATE (VZDENII(0:NION,0:NSTRA))
      ALLOCATE (VZDENPHI(0:NPHOT,0:NSTRA))
cdr  vectorial tallies, component-wise species summed up on -1,-2,-3
cdr                     sum over components and species (vector norm?): (0,0:nstra)
      ALLOCATE (MAPLI_VEC(-3:3*NPLS,0:NSTRA))
      ALLOCATE (MMPLI_VEC(-3:3*NPLS,0:NSTRA))
      ALLOCATE (MIPLI_VEC(-3:3*NPLS,0:NSTRA))
      ALLOCATE (MPHPLI_VEC(-3:3*NPLS,0:NSTRA))

cdr surface-averaged tallies
cdr  note: distinct from volumetric sources, these are ordered here
cdr  by receiving species
cdr  surface rates, for particle balance atoms
      ALLOCATE (POTATI(0:NATM,0:NSTRA))
      IF (NLSPCSCL_ATM) THEN
        ALLOCATE (PRFAAI(0:NATMP*NATMP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFAAI(0:NATM,0:NSTRA))
      END IF
      IF (NLSPCSCL_MOL) THEN
        ALLOCATE (PRFMAI(0:NATMP*NMOLP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFMAI(0:NATM,0:NSTRA))
      END IF
      IF (NLSPCSCL_ION) THEN
        ALLOCATE (PRFIAI(0:NATMP*NIONP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFIAI(0:NATM,0:NSTRA))
      END IF
      IF (NLSPCSCL_PHOT) THEN
        ALLOCATE (PRFPHAI(0:NATMP*NPHOTP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFPHAI(0:NATM,0:NSTRA))
      END IF
      ALLOCATE (PRFPAI(0:NATM,0:NSTRA))
cdr  surface rates, for particle balance molecules
      ALLOCATE (POTMLI(0:NMOL,0:NSTRA))
      IF (NLSPCSCL_ATM) THEN
        ALLOCATE (PRFAMI(0:NMOLP*NATMP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFAMI(0:NMOL,0:NSTRA))
      END IF
      IF (NLSPCSCL_MOL) THEN
        ALLOCATE (PRFMMI(0:NMOLP*NMOLP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFMMI(0:NMOL,0:NSTRA))
      END IF
      IF (NLSPCSCL_ION) THEN
        ALLOCATE (PRFIMI(0:NMOLP*NIONP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFIMI(0:NMOL,0:NSTRA))
      END IF
      IF (NLSPCSCL_PHOT) THEN
        ALLOCATE (PRFPHMI(0:NMOLP*NPHOTP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFPHMI(0:NMOL,0:NSTRA))
      END IF
      ALLOCATE (PRFPMI(0:NMOL,0:NSTRA))
cdr  surface rates, for particle balance test ions
      ALLOCATE (POTIOI(0:NION,0:NSTRA))
      IF (NLSPCSCL_ATM) THEN
        ALLOCATE (PRFAII(0:NIONP*NATMP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFAII(0:NION,0:NSTRA))
      END IF
      IF (NLSPCSCL_MOL) THEN
        ALLOCATE (PRFMII(0:NIONP*NMOLP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFMII(0:NION,0:NSTRA))
      END IF
      IF (NLSPCSCL_ION) THEN
        ALLOCATE (PRFIII(0:NIONP*NIONP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFIII(0:NION,0:NSTRA))
      END IF
      IF (NLSPCSCL_PHOT) THEN
        ALLOCATE (PRFPHII(0:NIONP*NPHOTP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFPHII(0:NION,0:NSTRA))
      END IF
      ALLOCATE (PRFPII(0:NION,0:NSTRA))
cdr  surface rates, for particle balance photons
      ALLOCATE (POTPHTI(0:NPHOT,0:NSTRA))
      IF (NLSPCSCL_ATM) THEN
        ALLOCATE (PRFAPHTI(0:NPHOTP*NATMP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFAPHTI(0:NPHOT,0:NSTRA))
      END IF
      IF (NLSPCSCL_MOL) THEN
        ALLOCATE (PRFMPHTI(0:NPHOTP*NMOLP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFMPHTI(0:NPHOT,0:NSTRA))
      END IF
      IF (NLSPCSCL_ION) THEN
        ALLOCATE (PRFIPHTI(0:NPHOTP*NIONP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFIPHTI(0:NPHOT,0:NSTRA))
      END IF
      IF (NLSPCSCL_PHOT) THEN
        ALLOCATE (PRFPHPHTI(0:NPHOTP*NPHOTP-1,0:NSTRA))
      ELSE
        ALLOCATE (PRFPHPHTI(0:NPHOT,0:NSTRA))
      END IF
      ALLOCATE (PRFPPHTI(0:NPHOT,0:NSTRA))

      ALLOCATE (EOTATI(0:NATM,0:NSTRA))
      ALLOCATE (ERFAAI(0:NATM,0:NSTRA))
      ALLOCATE (ERFMAI(0:NATM,0:NSTRA))
      ALLOCATE (ERFIAI(0:NATM,0:NSTRA))
      ALLOCATE (ERFPHAI(0:NATM,0:NSTRA))
      ALLOCATE (ERFPAI(0:NATM,0:NSTRA))

      ALLOCATE (EOTMLI(0:NMOL,0:NSTRA))
      ALLOCATE (ERFAMI(0:NMOL,0:NSTRA))
      ALLOCATE (ERFMMI(0:NMOL,0:NSTRA))
      ALLOCATE (ERFIMI(0:NMOL,0:NSTRA))
      ALLOCATE (ERFPHMI(0:NMOL,0:NSTRA))
      ALLOCATE (ERFPMI(0:NMOL,0:NSTRA))

      ALLOCATE (EOTIOI(0:NION,0:NSTRA))
      ALLOCATE (ERFAII(0:NION,0:NSTRA))
      ALLOCATE (ERFMII(0:NION,0:NSTRA))
      ALLOCATE (ERFIII(0:NION,0:NSTRA))
      ALLOCATE (ERFPHII(0:NION,0:NSTRA))
      ALLOCATE (ERFPII(0:NION,0:NSTRA))

      ALLOCATE (EOTPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (ERFAPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (ERFMPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (ERFIPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (ERFPHPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (ERFPPHTI(0:NPHOT,0:NSTRA))

cdr  surface incident fluxes of bulk ions (no fluxes emitted for bulk ions from surfaces)
      ALLOCATE (POTPLI(0:NPLS,0:NSTRA))
      ALLOCATE (EOTPLI(0:NPLS,0:NSTRA))

      ALLOCATE (SPTAATI(0:NATM,0:NSTRA))
      ALLOCATE (SPTMATI(0:NATM,0:NSTRA))
      ALLOCATE (SPTIATI(0:NATM,0:NSTRA))
      ALLOCATE (SPTPHATI(0:NATM,0:NSTRA))
      ALLOCATE (SPTPATI(0:NATM,0:NSTRA))
      ALLOCATE (SPTAMLI(0:NMOL,0:NSTRA))
      ALLOCATE (SPTMMLI(0:NMOL,0:NSTRA))
      ALLOCATE (SPTIMLI(0:NMOL,0:NSTRA))
      ALLOCATE (SPTPHMLI(0:NMOL,0:NSTRA))
      ALLOCATE (SPTPMLI(0:NMOL,0:NSTRA))
      ALLOCATE (SPTAIOI(0:NION,0:NSTRA))
      ALLOCATE (SPTMIOI(0:NION,0:NSTRA))
      ALLOCATE (SPTIIOI(0:NION,0:NSTRA))
      ALLOCATE (SPTPHIOI(0:NION,0:NSTRA))
      ALLOCATE (SPTPIOI(0:NION,0:NSTRA))
      ALLOCATE (SPTAPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (SPTMPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (SPTIPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (SPTPHPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (SPTPPHTI(0:NPHOT,0:NSTRA))
      ALLOCATE (SPTAPLI(0:NPLS,0:NSTRA))
      ALLOCATE (SPTMPLI(0:NPLS,0:NSTRA))
      ALLOCATE (SPTIPLI(0:NPLS,0:NSTRA))
      ALLOCATE (SPTPHPLI(0:NPLS,0:NSTRA))
      ALLOCATE (SPTPPLI(0:NPLS,0:NSTRA))
      ALLOCATE (SPTATTI(0:NSTRA))
      ALLOCATE (SPTMTTI(0:NSTRA))
      ALLOCATE (SPTITTI(0:NSTRA))
      ALLOCATE (SPTPHTTI(0:NSTRA))
      ALLOCATE (SPTPLTTI(0:NSTRA))
      ALLOCATE (SPTTTI(0:NSTRA))

      ALLOCATE (ADDSI(0:NADS,0:NSTRA))
      ALLOCATE (ALGSI(0:NALS,0:NSTRA))
      ALLOCATE (SPUMPI(0:NSPZ,0:NSTRA))

      ALLOCATE (RAELI(0:NATM,0:NSTRA))
      ALLOCATE (RMELI(0:NMOL,0:NSTRA))
      ALLOCATE (RIELI(0:NION,0:NSTRA))

      ALLOCATE (WTOTA(0:NATM,0:NSTRA))
      ALLOCATE (WTOTM(0:NMOL,0:NSTRA))
      ALLOCATE (WTOTI(0:NION,0:NSTRA))
      ALLOCATE (WTOTPH(0:NPHOT,0:NSTRA))
      ALLOCATE (WTOTP(0:NPLS,0:NSTRA))

      ALLOCATE (WTOTE(0:NSTRA))

cdr  etote still missing ??
      ALLOCATE (ETOTA(0:NSTRA))
      ALLOCATE (ETOTM(0:NSTRA))
      ALLOCATE (ETOTI(0:NSTRA))
      ALLOCATE (ETOTPH(0:NSTRA))
      ALLOCATE (ETOTP(0:NSTRA))

cdr  scored along the flight
      ALLOCATE (XMCP(0:NSTRA))
      ALLOCATE (EELFI(0:NION,0:NSTRA))
      ALLOCATE (PTRASH(0:NSTRA))
      ALLOCATE (ETRASH(0:NSTRA))
cdr  for scaling
      ALLOCATE (FLUXT(0:NSTRA))
      ALLOCATE (FLXFAC(0:NSTRA))
      ALLOCATE (FASCL(0:NATM,0:NSTRA))
      ALLOCATE (FMSCL(0:NMOL,0:NSTRA))
      ALLOCATE (FISCL(0:NION,0:NSTRA))
      ALLOCATE (FPHSCL(0:NPHOT,0:NSTRA))

      ALLOCATE (NADDI(NTALV))
      ALLOCATE (NFRSTI(NTALV))
      ALLOCATE (NDDWI(NTALS))
      ALLOCATE (NFRTWI(NTALS))
      ALLOCATE (NFSTVI(NTALV))
      ALLOCATE (NFSTWI(NTALS))
      ALLOCATE (NEXTVI(NTALV))
      ALLOCATE (NEXTWI(NTALS))

      WRITE (IUNMEM,'(A,T25,I15)')
     .       ' COUTAU ',NOUTAU*8 + NOUTTL*4

cdr this should go to another place?
      NADDI  = 0
      NFRSTI = 0
      NDDWI  = 0
      NFRTWI = 0
      NFSTVI = 0
      NFSTWI = 0
      NEXTVI = 0
      NEXTWI = 0

csw 19mar2013
      xmcp=0  !dr  this is done in init_coutau. So can go out here?
      LOGHELP = .TRUE.
      CALL EIRENE_INIT_COUTAU(LOGHELP)

      RETURN
      END SUBROUTINE EIRENE_ALLOC_COUTAU


      SUBROUTINE EIRENE_DEALLOC_COUTAU

      IF (.NOT.ALLOCATED(PDENAI)) RETURN

c tallies obtained by integration from volumetric tallies
      DEALLOCATE (PDENAI)
      DEALLOCATE (PDENMI)
      DEALLOCATE (PDENII)
      DEALLOCATE (PDENPHI)
      DEALLOCATE (EDENAI)
      DEALLOCATE (EDENMI)
      DEALLOCATE (EDENII)
      DEALLOCATE (EDENPHI)
      DEALLOCATE (PAELI)
      DEALLOCATE (PAATI)
      DEALLOCATE (PAMLI)
      DEALLOCATE (PAIOI)
      DEALLOCATE (PAPHTI)
      DEALLOCATE (PAPLI)
      DEALLOCATE (PMELI)
      DEALLOCATE (PMATI)
      DEALLOCATE (PMMLI)
      DEALLOCATE (PMIOI)
      DEALLOCATE (PMPHTI)
      DEALLOCATE (PMPLI)
      DEALLOCATE (PIELI)
      DEALLOCATE (PIATI)
      DEALLOCATE (PIMLI)
      DEALLOCATE (PIIOI)
      DEALLOCATE (PIPHTI)
      DEALLOCATE (PIPLI)
      DEALLOCATE (PPHELI)
      DEALLOCATE (PPHATI)
      DEALLOCATE (PPHMLI)
      DEALLOCATE (PPHIOI)
      DEALLOCATE (PPHPHTI)
      DEALLOCATE (PPHPLI)

      DEALLOCATE (EAELI)
      DEALLOCATE (EAATI)
      DEALLOCATE (EAMLI)
      DEALLOCATE (EAIOI)
      DEALLOCATE (EAPHTI)
      DEALLOCATE (EAPLI)
      DEALLOCATE (EMELI)
      DEALLOCATE (EMATI)
      DEALLOCATE (EMMLI)
      DEALLOCATE (EMIOI)
      DEALLOCATE (EMPHTI)
      DEALLOCATE (EMPLI)
      DEALLOCATE (EIELI)
      DEALLOCATE (EIATI)
      DEALLOCATE (EIMLI)
      DEALLOCATE (EIPHTI)
      DEALLOCATE (EIIOI)
      DEALLOCATE (EIPLI)
      DEALLOCATE (EPHELI)
      DEALLOCATE (EPHATI)
      DEALLOCATE (EPHMLI)
      DEALLOCATE (EPHPHTI)
      DEALLOCATE (EPHIOI)
      DEALLOCATE (EPHPLI)
      DEALLOCATE (RAELI)
      DEALLOCATE (RMELI)
      DEALLOCATE (RIELI)
      DEALLOCATE (ADDVI)
      DEALLOCATE (COLVI)
      DEALLOCATE (SNAPVI)
      DEALLOCATE (COPVI)
      DEALLOCATE (BGKVI)
      DEALLOCATE (ALGVI)
      DEALLOCATE (PGENAI)
      DEALLOCATE (PGENMI)
      DEALLOCATE (PGENII)
      DEALLOCATE (PGENPHI)
      DEALLOCATE (EGENAI)
      DEALLOCATE (EGENMI)
      DEALLOCATE (EGENII)
      DEALLOCATE (EGENPHI)
      DEALLOCATE (VGENAI)
      DEALLOCATE (VGENMI)
      DEALLOCATE (VGENII)
      DEALLOCATE (VGENPHI)
      DEALLOCATE (PPATI)
      DEALLOCATE (PPMLI)
      DEALLOCATE (PPIOI)
      DEALLOCATE (PPPHTI)
      DEALLOCATE (PPPLI)
      DEALLOCATE (EPATI)
      DEALLOCATE (EPMLI)
      DEALLOCATE (EPIOI)
      DEALLOCATE (EPPHTI)
      DEALLOCATE (EPPLI)
      DEALLOCATE (VXDENAI)
      DEALLOCATE (VXDENMI)
      DEALLOCATE (VXDENII)
      DEALLOCATE (VXDENPHI)
      DEALLOCATE (VYDENAI)
      DEALLOCATE (VYDENMI)
      DEALLOCATE (VYDENII)
      DEALLOCATE (VYDENPHI)
      DEALLOCATE (VZDENAI)
      DEALLOCATE (VZDENMI)
      DEALLOCATE (VZDENII)
      DEALLOCATE (VZDENPHI)
      DEALLOCATE (MAPLI_VEC)
      DEALLOCATE (MMPLI_VEC)
      DEALLOCATE (MIPLI_VEC)
      DEALLOCATE (MPHPLI_VEC)

c  integrated surface-averaged tallies
      DEALLOCATE (POTATI)
      DEALLOCATE (PRFAAI)
      DEALLOCATE (PRFMAI)
      DEALLOCATE (PRFIAI)
      DEALLOCATE (PRFPHAI)
      DEALLOCATE (PRFPAI)
      DEALLOCATE (POTMLI)
      DEALLOCATE (PRFAMI)
      DEALLOCATE (PRFMMI)
      DEALLOCATE (PRFIMI)
      DEALLOCATE (PRFPHMI)
      DEALLOCATE (PRFPMI)
      DEALLOCATE (POTIOI)
      DEALLOCATE (PRFAII)
      DEALLOCATE (PRFMII)
      DEALLOCATE (PRFIII)
      DEALLOCATE (PRFPHII)
      DEALLOCATE (PRFPII)
      DEALLOCATE (POTPHTI)
      DEALLOCATE (PRFAPHTI)
      DEALLOCATE (PRFMPHTI)
      DEALLOCATE (PRFIPHTI)
      DEALLOCATE (PRFPHPHTI)
      DEALLOCATE (PRFPPHTI)
      DEALLOCATE (POTPLI)
c
      DEALLOCATE (EOTATI)
      DEALLOCATE (ERFAAI)
      DEALLOCATE (ERFMAI)
      DEALLOCATE (ERFIAI)
      DEALLOCATE (ERFPHAI)
      DEALLOCATE (ERFPAI)
      DEALLOCATE (EOTMLI)
      DEALLOCATE (ERFAMI)
      DEALLOCATE (ERFMMI)
      DEALLOCATE (ERFIMI)
      DEALLOCATE (ERFPHMI)
      DEALLOCATE (ERFPMI)
      DEALLOCATE (EOTIOI)
      DEALLOCATE (ERFAII)
      DEALLOCATE (ERFMII)
      DEALLOCATE (ERFIII)
      DEALLOCATE (ERFPHII)
      DEALLOCATE (ERFPII)
      DEALLOCATE (EOTPHTI)
      DEALLOCATE (ERFAPHTI)
      DEALLOCATE (ERFMPHTI)
      DEALLOCATE (ERFIPHTI)
      DEALLOCATE (ERFPHPHTI)
      DEALLOCATE (ERFPPHTI)
      DEALLOCATE (EOTPLI)
c
      DEALLOCATE (SPTAATI)
      DEALLOCATE (SPTAMLI)
      DEALLOCATE (SPTAIOI)
      DEALLOCATE (SPTAPHTI)
      DEALLOCATE (SPTAPLI)
      DEALLOCATE (SPTMATI)
      DEALLOCATE (SPTMMLI)
      DEALLOCATE (SPTMIOI)
      DEALLOCATE (SPTMPHTI)
      DEALLOCATE (SPTMPLI)
      DEALLOCATE (SPTIATI)
      DEALLOCATE (SPTIMLI)
      DEALLOCATE (SPTIIOI)
      DEALLOCATE (SPTIPHTI)
      DEALLOCATE (SPTIPLI)
      DEALLOCATE (SPTPHATI)
      DEALLOCATE (SPTPHMLI)
      DEALLOCATE (SPTPHIOI)
      DEALLOCATE (SPTPHPHTI)
      DEALLOCATE (SPTPHPLI)
      DEALLOCATE (SPTPATI)
      DEALLOCATE (SPTPMLI)
      DEALLOCATE (SPTPIOI)
      DEALLOCATE (SPTPPHTI)
      DEALLOCATE (SPTPPLI)
      DEALLOCATE (SPTATTI)
      DEALLOCATE (SPTMTTI)
      DEALLOCATE (SPTITTI)
      DEALLOCATE (SPTPHTTI)
      DEALLOCATE (SPTPLTTI)
      DEALLOCATE (SPTTTI)
      DEALLOCATE (ADDSI)
      DEALLOCATE (ALGSI)
      DEALLOCATE (SPUMPI)

c  global fluxes from primary source
      DEALLOCATE (WTOTA)
      DEALLOCATE (WTOTM)
      DEALLOCATE (WTOTI)
      DEALLOCATE (WTOTP)
      DEALLOCATE (WTOTPH)
      DEALLOCATE (WTOTE)
      DEALLOCATE (ETOTA)
      DEALLOCATE (ETOTM)
      DEALLOCATE (ETOTI)
      DEALLOCATE (ETOTP)
      DEALLOCATE (ETOTPH)

c scored "on the fly"
      DEALLOCATE (XMCP)
cdr At surface events of test ions:
cdr from energy gained from (sheath) electric field
      DEALLOCATE (EELFI)
cdr trash: particles killed during tracing, "fail-safe"
      DEALLOCATE (PTRASH)
      DEALLOCATE (ETRASH)

c stuff for scaling tallies
      DEALLOCATE (FLUXT)
      DEALLOCATE (FLXFAC)
      DEALLOCATE (FASCL)
      DEALLOCATE (FMSCL)
      DEALLOCATE (FISCL)
      DEALLOCATE (FPHSCL)

c size of tallies
      DEALLOCATE (NADDI)
      DEALLOCATE (NFRSTI)
      DEALLOCATE (NDDWI)
      DEALLOCATE (NFRTWI)
      DEALLOCATE (NFSTVI)
      DEALLOCATE (NFSTWI)
      DEALLOCATE (NEXTVI)
      DEALLOCATE (NEXTWI)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_COUTAU


      SUBROUTINE EIRENE_INIT_COUTAU(LSTR)
cdr input: lstr(istra)=true : this stratum ISTRA exists
cdr                    false: stratum ISTRA turned off in this run

      LOGICAL, INTENT(IN) :: LSTR(NSTRA)
      INTEGER :: JSTRA, ISTRA

      DO JSTRA=0,NSTRA
        ISTRA = JSTRA

cdr  what is this ?
cpb  for short cycle:
cpb  keep results of strata which are not calculated in this timestep/iteration
        IF ((ISTRA >= 1) .AND. (IFRST > 0)) THEN
          IF (.NOT. LSTR(ISTRA)) CYCLE
        END IF

cdr particle densities
        PDENAI(:,ISTRA)  = 0._DP
        PDENMI(:,ISTRA)  = 0._DP
        PDENII(:,ISTRA)  = 0._DP
        PDENPHI(:,ISTRA) = 0._DP
cdr  energy densities
        EDENAI(:,ISTRA)  = 0._DP
        EDENMI(:,ISTRA)  = 0._DP
        EDENII(:,ISTRA)  = 0._DP
        EDENPHI(:,ISTRA) = 0._DP

cdr volumetric tallies for particle balance, sources from atoms for el, a,m,i,ph,pl
        PAELI(ISTRA)     = 0._DP
        PAATI(:,ISTRA)   = 0._DP
        PAMLI(:,ISTRA)   = 0._DP
        PAIOI(:,ISTRA)   = 0._DP
        PAPHTI(:,ISTRA)  = 0._DP
        PAPLI(:,ISTRA)   = 0._DP

        PMELI(ISTRA)     = 0._DP
        PMATI(:,ISTRA)   = 0._DP
        PMMLI(:,ISTRA)   = 0._DP
        PMIOI(:,ISTRA)   = 0._DP
        PMPHTI(:,ISTRA)  = 0._DP
        PMPLI(:,ISTRA)   = 0._DP

        PIELI(ISTRA)     = 0._DP
        PIATI(:,ISTRA)   = 0._DP
        PIMLI(:,ISTRA)   = 0._DP
        PIIOI(:,ISTRA)   = 0._DP
        PIPHTI(:,ISTRA)  = 0._DP
        PIPLI(:,ISTRA)   = 0._DP
cdr
        PPHELI(ISTRA)    = 0._DP
        PPHATI(:,ISTRA)  = 0._DP
        PPHMLI(:,ISTRA)  = 0._DP
        PPHIOI(:,ISTRA)  = 0._DP
        PPHPHTI(:,ISTRA) = 0._DP
        PPHPLI(:,ISTRA)  = 0._DP

cdr volumetric tallies for energy balance, sources from atoms for el, a,m,i,ph,pl
        EAELI(ISTRA)   = 0._DP
        EAATI(ISTRA)   = 0._DP
        EAMLI(ISTRA)   = 0._DP
        EAIOI(ISTRA)   = 0._DP
        EAPHTI(ISTRA)  = 0._DP
        EAPLI(:,ISTRA) = 0._DP

        EMELI(ISTRA)   = 0._DP
        EMATI(ISTRA)   = 0._DP
        EMMLI(ISTRA)   = 0._DP
        EMIOI(ISTRA)   = 0._DP
        EMPHTI(ISTRA)  = 0._DP
        EMPLI(:,ISTRA) = 0._DP

        EIELI(ISTRA)   = 0._DP
        EIATI(ISTRA)   = 0._DP
        EIMLI(ISTRA)   = 0._DP
        EIIOI(ISTRA)   = 0._DP
        EIPHTI(ISTRA)  = 0._DP
        EIPLI(:,ISTRA) = 0._DP

        EPHELI(ISTRA)  = 0._DP
        EPHATI(ISTRA)  = 0._DP
        EPHMLI(ISTRA)  = 0._DP
        EPHIOI(ISTRA)  = 0._DP
        EPHPHTI(ISTRA) = 0._DP
        EPHPLI(:,ISTRA)= 0._DP

        RAELI(:,ISTRA) = 0._DP
        RMELI(:,ISTRA) = 0._DP
        RIELI(:,ISTRA) = 0._DP

        ADDVI(:,ISTRA)   = 0._DP
        COLVI(:,ISTRA)   = 0._DP
        SNAPVI(:,ISTRA)  = 0._DP
        COPVI(:,ISTRA)   = 0._DP
        BGKVI(:,ISTRA)   = 0._DP
        ALGVI(:,ISTRA)   = 0._DP

        PGENAI(:,ISTRA)  = 0._DP
        PGENMI(:,ISTRA)  = 0._DP
        PGENII(:,ISTRA)  = 0._DP
        PGENPHI(:,ISTRA) = 0._DP
        EGENAI(:,ISTRA)  = 0._DP
        EGENMI(:,ISTRA)  = 0._DP
        EGENII(:,ISTRA)  = 0._DP
        EGENPHI(:,ISTRA) = 0._DP
        VGENAI(:,ISTRA)  = 0._DP
        VGENMI(:,ISTRA)  = 0._DP
        VGENII(:,ISTRA)  = 0._DP
        VGENPHI(:,ISTRA) = 0._DP

cdr  particle sources from pl, for electrons: tally ppeli missing ??
        PPATI(:,ISTRA)  = 0._DP
        PPMLI(:,ISTRA)  = 0._DP
        PPIOI(:,ISTRA)  = 0._DP
        PPPHTI(:,ISTRA) = 0._DP
        PPPLI(:,ISTRA)  = 0._DP

cdr  energy sources from pl, for electrons: tally epeli missing ??
        EPATI(ISTRA)  = 0._DP
        EPMLI(ISTRA)  = 0._DP
        EPIOI(ISTRA)  = 0._DP
        EPPHTI(ISTRA) = 0._DP
        EPPLI(:,ISTRA)= 0._DP

        VXDENAI(:,ISTRA)  = 0._DP
        VXDENMI(:,ISTRA)  = 0._DP
        VXDENII(:,ISTRA)  = 0._DP
        VXDENPHI(:,ISTRA) = 0._DP
        VYDENAI(:,ISTRA)  = 0._DP
        VYDENMI(:,ISTRA)  = 0._DP
        VYDENII(:,ISTRA)  = 0._DP
        VYDENPHI(:,ISTRA) = 0._DP
        VZDENAI(:,ISTRA)  = 0._DP
        VZDENMI(:,ISTRA)  = 0._DP
        VZDENII(:,ISTRA)  = 0._DP
        VZDENPHI(:,ISTRA) = 0._DP
        MAPLI_VEC(:,ISTRA)    = 0._DP
        MMPLI_VEC(:,ISTRA)    = 0._DP
        MIPLI_VEC(:,ISTRA)    = 0._DP
        MPHPLI_VEC(:,ISTRA)   = 0._DP

cdr  surface tallies
        POTATI(:,ISTRA) = 0._DP
        PRFAAI(:,ISTRA) = 0._DP
        PRFMAI(:,ISTRA) = 0._DP
        PRFIAI(:,ISTRA) = 0._DP
        PRFPHAI(:,ISTRA)= 0._DP
        PRFPAI(:,ISTRA) = 0._DP
        POTMLI(:,ISTRA) = 0._DP
        PRFAMI(:,ISTRA) = 0._DP
        PRFMMI(:,ISTRA) = 0._DP
        PRFIMI(:,ISTRA) = 0._DP
        PRFPHMI(:,ISTRA)= 0._DP
        PRFPMI(:,ISTRA) = 0._DP
        POTIOI(:,ISTRA) = 0._DP
        PRFAII(:,ISTRA) = 0._DP
        PRFMII(:,ISTRA) = 0._DP
        PRFIII(:,ISTRA) = 0._DP
        PRFPHII(:,ISTRA)= 0._DP
        PRFPII(:,ISTRA) = 0._DP
        POTPHTI(:,ISTRA)= 0._DP
        PRFAPHTI(:,ISTRA) = 0._DP
        PRFMPHTI(:,ISTRA) = 0._DP
        PRFIPHTI(:,ISTRA) = 0._DP
        PRFPHPHTI(:,ISTRA)= 0._DP
        PRFPPHTI(:,ISTRA) = 0._DP
        POTPLI(:,ISTRA) = 0._DP
        EOTATI(:,ISTRA) = 0._DP
        ERFAAI(:,ISTRA) = 0._DP
        ERFMAI(:,ISTRA) = 0._DP
        ERFIAI(:,ISTRA) = 0._DP
        ERFPHAI(:,ISTRA)= 0._DP
        ERFPAI(:,ISTRA) = 0._DP
        EOTMLI(:,ISTRA) = 0._DP
        ERFAMI(:,ISTRA) = 0._DP
        ERFMMI(:,ISTRA) = 0._DP
        ERFIMI(:,ISTRA) = 0._DP
        ERFPHMI(:,ISTRA)= 0._DP
        ERFPMI(:,ISTRA) = 0._DP
        EOTIOI(:,ISTRA) = 0._DP
        ERFAII(:,ISTRA) = 0._DP
        ERFMII(:,ISTRA) = 0._DP
        ERFIII(:,ISTRA) = 0._DP
        ERFPHII(:,ISTRA)= 0._DP
        ERFPII(:,ISTRA) = 0._DP
        EOTPHTI(:,ISTRA)= 0._DP
        ERFAPHTI(:,ISTRA) = 0._DP
        ERFMPHTI(:,ISTRA) = 0._DP
        ERFIPHTI(:,ISTRA) = 0._DP
        ERFPHPHTI(:,ISTRA)= 0._DP
        ERFPPHTI(:,ISTRA) = 0._DP
        EOTPLI(:,ISTRA)   = 0._DP
        SPTAATI(:,ISTRA) = 0._DP
        SPTMATI(:,ISTRA) = 0._DP
        SPTIATI(:,ISTRA) = 0._DP
        SPTPHATI(:,ISTRA)= 0._DP
        SPTPATI(:,ISTRA) = 0._DP
        SPTAMLI(:,ISTRA) = 0._DP
        SPTMMLI(:,ISTRA) = 0._DP
        SPTIMLI(:,ISTRA) = 0._DP
        SPTPHMLI(:,ISTRA)= 0._DP
        SPTPMLI(:,ISTRA) = 0._DP
        SPTAIOI(:,ISTRA) = 0._DP
        SPTMIOI(:,ISTRA) = 0._DP
        SPTIIOI(:,ISTRA) = 0._DP
        SPTPHIOI(:,ISTRA)= 0._DP
        SPTPIOI(:,ISTRA) = 0._DP
        SPTAPHTI(:,ISTRA) = 0._DP
        SPTMPHTI(:,ISTRA) = 0._DP
        SPTIPHTI(:,ISTRA) = 0._DP
        SPTPHPHTI(:,ISTRA)= 0._DP
        SPTPPHTI(:,ISTRA) = 0._DP
        SPTAPLI(:,ISTRA) = 0._DP
        SPTMPLI(:,ISTRA) = 0._DP
        SPTIPLI(:,ISTRA) = 0._DP
        SPTPHPLI(:,ISTRA)= 0._DP
        SPTPPLI(:,ISTRA) = 0._DP
        SPTATTI(ISTRA) = 0._DP
        SPTMTTI(ISTRA) = 0._DP
        SPTITTI(ISTRA) = 0._DP
        SPTPHTTI(ISTRA)= 0._DP
        SPTPLTTI(ISTRA)= 0._DP
        SPTTTI(ISTRA)  = 0._DP
        ADDSI(:,ISTRA) = 0._DP
        ALGSI(:,ISTRA) = 0._DP
        SPUMPI(:,ISTRA)= 0._DP

        WTOTA(:,ISTRA)  = 0._DP
        WTOTM(:,ISTRA)  = 0._DP
        WTOTI(:,ISTRA)  = 0._DP
        WTOTP(:,ISTRA)  = 0._DP
        WTOTPH(:,ISTRA) = 0._DP
        WTOTE(ISTRA)  = 0._DP
        ETOTA(ISTRA)  = 0._DP
        ETOTM(ISTRA)  = 0._DP
        ETOTI(ISTRA)  = 0._DP
        ETOTP(ISTRA)  = 0._DP
        ETOTPH(ISTRA) = 0._DP
cdr  scored along the flight
        EELFI(:,ISTRA)= 0._DP
        PTRASH(ISTRA) = 0._DP
        ETRASH(ISTRA) = 0._DP
        XMCP(ISTRA)   = 0._DP

cdr  for scaling
        FLUXT(ISTRA)  = 0._DP
        FLXFAC(ISTRA) = 0._DP
        FASCL(:,ISTRA) = 0._DP
        FMSCL(:,ISTRA) = 0._DP
        FISCL(:,ISTRA) = 0._DP
        FPHSCL(:,ISTRA)= 0._DP

      END DO
      IFRST = 1

      RETURN
      END SUBROUTINE EIRENE_INIT_COUTAU

C     The following SUBROUTINE is for reinitialization of EIRENE

      SUBROUTINE EIRENE_INIT_COUTAU_REINIT
      IMPLICIT NONE
      IFRST = 0
      return
      END SUBROUTINE EIRENE_INIT_COUTAU_REINIT


      SUBROUTINE EIRENE_WRITE_COUTAU (OUTAU,IUNOUT)
      REAL(DP) :: OUTAU(:)
      INTEGER, INTENT(IN) :: IUNOUT
      INTEGER :: IA, IE
      EXTERNAL :: EIRENE_EXIT_OWN

      IF (SIZE(OUTAU) < NOUTAU) THEN
        WRITE (iunout,*)
     .  ' SIZE OF ARRAY OUTAU TOO SMALL IN WRITE_COUTAU '
        WRITE (iunout,*) ' OUTAU(',NOUTAU,') REQUIRED '
        CALL EIRENE_EXIT_OWN(1)
      END IF

      IA = 1
      IE = IA - 1 + SIZE(PDENAI)
      OUTAU(IA:IE) = PACK(PDENAI,.TRUE.)
c     write (iunout,*) '1, NATMP ',(IE-IA+1)/NSTRAP, NATMP

      IA = IE + 1
      IE = IA - 1 + SIZE(PDENMI)
      OUTAU(IA:IE) = PACK(PDENMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PDENII)
      OUTAU(IA:IE) = PACK(PDENII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PDENPHI)
      OUTAU(IA:IE) = PACK(PDENPHI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EDENAI)
      OUTAU(IA:IE) = PACK(EDENAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EDENMI)
      OUTAU(IA:IE) = PACK(EDENMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EDENII)
      OUTAU(IA:IE) = PACK(EDENII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EDENPHI)
      OUTAU(IA:IE) = PACK(EDENPHI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PAELI)
      OUTAU(IA:IE) = PACK(PAELI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PAATI)
      OUTAU(IA:IE) = PACK(PAATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PAMLI)
      OUTAU(IA:IE) = PACK(PAMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PAIOI)
      OUTAU(IA:IE) = PACK(PAIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PAPHTI)
      OUTAU(IA:IE) = PACK(PAPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PAPLI)
      OUTAU(IA:IE) = PACK(PAPLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PMELI)
      OUTAU(IA:IE) = PACK(PMELI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PMATI)
      OUTAU(IA:IE) = PACK(PMATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PMMLI)
      OUTAU(IA:IE) = PACK(PMMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PMIOI)
      OUTAU(IA:IE) = PACK(PMIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PMPHTI)
      OUTAU(IA:IE) = PACK(PMPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PMPLI)
      OUTAU(IA:IE) = PACK(PMPLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PIELI)
      OUTAU(IA:IE) = PACK(PIELI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PIATI)
      OUTAU(IA:IE) = PACK(PIATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PIMLI)
      OUTAU(IA:IE) = PACK(PIMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PIIOI)
      OUTAU(IA:IE) = PACK(PIIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PIPHTI)
      OUTAU(IA:IE) = PACK(PIPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PIPLI)
      OUTAU(IA:IE) = PACK(PIPLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHELI)
      OUTAU(IA:IE) = PACK(PPHELI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHATI)
      OUTAU(IA:IE) = PACK(PPHATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHMLI)
      OUTAU(IA:IE) = PACK(PPHMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHIOI)
      OUTAU(IA:IE) = PACK(PPHIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHPHTI)
      OUTAU(IA:IE) = PACK(PPHPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHPLI)
      OUTAU(IA:IE) = PACK(PPHPLI ,.TRUE.)

cdr  particle sources done.


      IA = IE + 1
      IE = IA - 1 + SIZE(EAELI)
      OUTAU(IA:IE) = PACK(EAELI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EAATI)
      OUTAU(IA:IE) = PACK(EAATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EAMLI)
      OUTAU(IA:IE) = PACK(EAMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EAIOI)
      OUTAU(IA:IE) = PACK(EAIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EAPHTI)
      OUTAU(IA:IE) = PACK(EAPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EAPLI)
      OUTAU(IA:IE) = PACK(EAPLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EMELI)
      OUTAU(IA:IE) = PACK(EMELI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EMATI)
      OUTAU(IA:IE) = PACK(EMATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EMMLI)
      OUTAU(IA:IE) = PACK(EMMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EMIOI)
      OUTAU(IA:IE) = PACK(EMIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EMPHTI)
      OUTAU(IA:IE) = PACK(EMPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EMPLI)
      OUTAU(IA:IE) = PACK(EMPLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EIELI)
      OUTAU(IA:IE) = PACK(EIELI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EIATI)
      OUTAU(IA:IE) = PACK(EIATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EIMLI)
      OUTAU(IA:IE) = PACK(EIMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EIIOI)
      OUTAU(IA:IE) = PACK(EIIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EIPHTI)
      OUTAU(IA:IE) = PACK(EIPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EIPLI)
      OUTAU(IA:IE) = PACK(EIPLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHELI)
      OUTAU(IA:IE) = PACK(EPHELI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHATI)
      OUTAU(IA:IE) = PACK(EPHATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHMLI)
      OUTAU(IA:IE) = PACK(EPHMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHIOI)
      OUTAU(IA:IE) = PACK(EPHIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHPHTI)
      OUTAU(IA:IE) = PACK(EPHPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHPLI)
      OUTAU(IA:IE) = PACK(EPHPLI ,.TRUE.)

cdr energy sources done.Next: additional tallies

      IA = IE + 1
      IE = IA - 1 + SIZE(ADDVI)
      OUTAU(IA:IE) = PACK(ADDVI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(COLVI)
      OUTAU(IA:IE) = PACK(COLVI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SNAPVI)
      OUTAU(IA:IE) = PACK(SNAPVI,.TRUE.)
c     write (iunout,*) 'snap, NSNVP ',(IE-IA+1)/NSTRAP, NSNVP
      IA = IE + 1
      IE = IA - 1 + SIZE(COPVI)
      OUTAU(IA:IE) = PACK(COPVI ,.TRUE.)
c     write (iunout,*) 'copv, NCPVP ',(IE-IA+1)/NSTRAP, NCPVP
      IA = IE + 1
      IE = IA - 1 + SIZE(BGKVI)
      OUTAU(IA:IE) = PACK(BGKVI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ALGVI)
      OUTAU(IA:IE) = PACK(ALGVI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PGENAI)
      OUTAU(IA:IE) = PACK(PGENAI,.TRUE.)
c     write (iunout,*) 'pgena, NATMP ',(IE-IA+1)/NSTRAP, NATMP

      IA = IE + 1
      IE = IA - 1 + SIZE(PGENMI)
      OUTAU(IA:IE) = PACK(PGENMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PGENII)
      OUTAU(IA:IE) = PACK(PGENII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PGENPHI)
      OUTAU(IA:IE) = PACK(PGENPHI,.TRUE.)
c     write (iunout,*) 'pgenph, NPHOTP ',(IE-IA+1)/NSTRAP, NPHOTP

      IA = IE + 1
      IE = IA - 1 + SIZE(EGENAI)
      OUTAU(IA:IE) = PACK(EGENAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EGENMI)
      OUTAU(IA:IE) = PACK(EGENMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EGENII)
      OUTAU(IA:IE) = PACK(EGENII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EGENPHI)
      OUTAU(IA:IE) = PACK(EGENPHI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VGENAI)
      OUTAU(IA:IE) = PACK(VGENAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VGENMI)
      OUTAU(IA:IE) = PACK(VGENMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VGENII)
      OUTAU(IA:IE) = PACK(VGENII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VGENPHI)
      OUTAU(IA:IE) = PACK(VGENPHI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPATI)
      OUTAU(IA:IE) = PACK(PPATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPMLI)
      OUTAU(IA:IE) = PACK(PPMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPIOI)
      OUTAU(IA:IE) = PACK(PPIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPPHTI)
      OUTAU(IA:IE) = PACK(PPPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PPPLI)
      OUTAU(IA:IE) = PACK(PPPLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPATI)
      OUTAU(IA:IE) = PACK(EPATI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPMLI)
      OUTAU(IA:IE) = PACK(EPMLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPIOI)
      OUTAU(IA:IE) = PACK(EPIOI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPPHTI)
      OUTAU(IA:IE) = PACK(EPPHTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EPPLI)
      OUTAU(IA:IE) = PACK(EPPLI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VXDENAI)
      OUTAU(IA:IE) = PACK(VXDENAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VXDENMI)
      OUTAU(IA:IE) = PACK(VXDENMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VXDENII)
      OUTAU(IA:IE) = PACK(VXDENII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VXDENPHI)
      OUTAU(IA:IE) = PACK(VXDENPHI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VYDENAI)
      OUTAU(IA:IE) = PACK(VYDENAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VYDENMI)
      OUTAU(IA:IE) = PACK(VYDENMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VYDENII)
      OUTAU(IA:IE) = PACK(VYDENII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VYDENPHI)
      OUTAU(IA:IE) = PACK(VYDENPHI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VZDENAI)
      OUTAU(IA:IE) = PACK(VZDENAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VZDENMI)
      OUTAU(IA:IE) = PACK(VZDENMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VZDENII)
      OUTAU(IA:IE) = PACK(VZDENII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(VZDENPHI)
      OUTAU(IA:IE) = PACK(VZDENPHI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(MAPLI_VEC)
      OUTAU(IA:IE) = PACK(MAPLI_VEC,.TRUE.)
c     write (iunout,*) 'mapl, NPLSP ',(IE-IA+1)/NSTRAP, NPLSP
      IA = IE + 1
      IE = IA - 1 + SIZE(MMPLI_VEC)
      OUTAU(IA:IE) = PACK(MMPLI_VEC,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(MIPLI_VEC)
      OUTAU(IA:IE) = PACK(MIPLI_VEC,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(MPHPLI_VEC)
      OUTAU(IA:IE) = PACK(MPHPLI_VEC,.TRUE.)

      IF (IE /= NOUTA1) THEN
        WRITE (iunout,*) ' ERROR IN WRITE_COUTAU '
        WRITE (iunout,*) ' NOUTA1 = ',NOUTA1,' IE = ',IE
        CALL EIRENE_EXIT_OWN(1)
      END IF

cdr  surface tallies
      IA = IE + 1
      IE = IA - 1 + SIZE(POTATI)
      OUTAU(IA:IE) = PACK(POTATI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFAAI)
      OUTAU(IA:IE) = PACK(PRFAAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFMAI)
      OUTAU(IA:IE) = PACK(PRFMAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFIAI)
      OUTAU(IA:IE) = PACK(PRFIAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPHAI)
      OUTAU(IA:IE) = PACK(PRFPHAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPAI)
      OUTAU(IA:IE) = PACK(PRFPAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(POTMLI)
      OUTAU(IA:IE) = PACK(POTMLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFAMI)
      OUTAU(IA:IE) = PACK(PRFAMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFMMI)
      OUTAU(IA:IE) = PACK(PRFMMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFIMI)
      OUTAU(IA:IE) = PACK(PRFIMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPHMI)
      OUTAU(IA:IE) = PACK(PRFPHMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPMI)
      OUTAU(IA:IE) = PACK(PRFPMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(POTIOI)
      OUTAU(IA:IE) = PACK(POTIOI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFAII)
      OUTAU(IA:IE) = PACK(PRFAII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFMII)
      OUTAU(IA:IE) = PACK(PRFMII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFIII)
      OUTAU(IA:IE) = PACK(PRFIII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPHII)
      OUTAU(IA:IE) = PACK(PRFPHII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPII)
      OUTAU(IA:IE) = PACK(PRFPII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(POTPHTI)
      OUTAU(IA:IE) = PACK(POTPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFAPHTI)
      OUTAU(IA:IE) = PACK(PRFAPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFMPHTI)
      OUTAU(IA:IE) = PACK(PRFMPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFIPHTI)
      OUTAU(IA:IE) = PACK(PRFIPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPHPHTI)
      OUTAU(IA:IE) = PACK(PRFPHPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPPHTI)
      OUTAU(IA:IE) = PACK(PRFPPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(POTPLI)
      OUTAU(IA:IE) = PACK(POTPLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTATI)
      OUTAU(IA:IE) = PACK(EOTATI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFAAI)
      OUTAU(IA:IE) = PACK(ERFAAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFMAI)
      OUTAU(IA:IE) = PACK(ERFMAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFIAI)
      OUTAU(IA:IE) = PACK(ERFIAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPHAI)
      OUTAU(IA:IE) = PACK(ERFPHAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPAI)
      OUTAU(IA:IE) = PACK(ERFPAI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTMLI)
      OUTAU(IA:IE) = PACK(EOTMLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFAMI)
      OUTAU(IA:IE) = PACK(ERFAMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFMMI)
      OUTAU(IA:IE) = PACK(ERFMMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFIMI)
      OUTAU(IA:IE) = PACK(ERFIMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPHMI)
      OUTAU(IA:IE) = PACK(ERFPHMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPMI)
      OUTAU(IA:IE) = PACK(ERFPMI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTIOI)
      OUTAU(IA:IE) = PACK(EOTIOI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFAII)
      OUTAU(IA:IE) = PACK(ERFAII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFMII)
      OUTAU(IA:IE) = PACK(ERFMII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFIII)
      OUTAU(IA:IE) = PACK(ERFIII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPHII)
      OUTAU(IA:IE) = PACK(ERFPHII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPII)
      OUTAU(IA:IE) = PACK(ERFPII,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTPHTI)
      OUTAU(IA:IE) = PACK(EOTPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFAPHTI)
      OUTAU(IA:IE) = PACK(ERFAPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFMPHTI)
      OUTAU(IA:IE) = PACK(ERFMPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFIPHTI)
      OUTAU(IA:IE) = PACK(ERFIPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPHPHTI)
      OUTAU(IA:IE) = PACK(ERFPHPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPPHTI)
      OUTAU(IA:IE) = PACK(ERFPPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTPLI)
      OUTAU(IA:IE) = PACK(EOTPLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAATI)
      OUTAU(IA:IE) = PACK(SPTAATI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMATI)
      OUTAU(IA:IE) = PACK(SPTMATI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIATI)
      OUTAU(IA:IE) = PACK(SPTIATI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHATI)
      OUTAU(IA:IE) = PACK(SPTPHATI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPATI)
      OUTAU(IA:IE) = PACK(SPTPATI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAMLI)
      OUTAU(IA:IE) = PACK(SPTAMLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMMLI)
      OUTAU(IA:IE) = PACK(SPTMMLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIMLI)
      OUTAU(IA:IE) = PACK(SPTIMLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHMLI)
      OUTAU(IA:IE) = PACK(SPTPHMLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPMLI)
      OUTAU(IA:IE) = PACK(SPTPMLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAIOI)
      OUTAU(IA:IE) = PACK(SPTAIOI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMIOI)
      OUTAU(IA:IE) = PACK(SPTMIOI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIIOI)
      OUTAU(IA:IE) = PACK(SPTIIOI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHIOI)
      OUTAU(IA:IE) = PACK(SPTPHIOI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPIOI)
      OUTAU(IA:IE) = PACK(SPTPIOI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAPHTI)
      OUTAU(IA:IE) = PACK(SPTAPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMPHTI)
      OUTAU(IA:IE) = PACK(SPTMPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIPHTI)
      OUTAU(IA:IE) = PACK(SPTIPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHPHTI)
      OUTAU(IA:IE) = PACK(SPTPHPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPPHTI)
      OUTAU(IA:IE) = PACK(SPTPPHTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAPLI)
      OUTAU(IA:IE) = PACK(SPTAPLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMPLI)
      OUTAU(IA:IE) = PACK(SPTMPLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIPLI)
      OUTAU(IA:IE) = PACK(SPTIPLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHPLI)
      OUTAU(IA:IE) = PACK(SPTPHPLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPPLI)
      OUTAU(IA:IE) = PACK(SPTPPLI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTATTI)
      OUTAU(IA:IE) = PACK(SPTATTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMTTI)
      OUTAU(IA:IE) = PACK(SPTMTTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTITTI)
      OUTAU(IA:IE) = PACK(SPTITTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHTTI)
      OUTAU(IA:IE) = PACK(SPTPHTTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPLTTI)
      OUTAU(IA:IE) = PACK(SPTPLTTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTTTI)
      OUTAU(IA:IE) = PACK(SPTTTI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ADDSI)
      OUTAU(IA:IE) = PACK(ADDSI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ALGSI)
      OUTAU(IA:IE) = PACK(ALGSI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(SPUMPI)
      OUTAU(IA:IE) = PACK(SPUMPI,.TRUE.)


      IA = IE + 1
      IE = IA - 1 + SIZE(RAELI)
      OUTAU(IA:IE) = PACK(RAELI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(RMELI)
      OUTAU(IA:IE) = PACK(RMELI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(RIELI)
      OUTAU(IA:IE) = PACK(RIELI,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTA)
      OUTAU(IA:IE) = PACK(WTOTA ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTM)
      OUTAU(IA:IE) = PACK(WTOTM ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTI)
      OUTAU(IA:IE) = PACK(WTOTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTPH)
      OUTAU(IA:IE) = PACK(WTOTPH ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTP)
      OUTAU(IA:IE) = PACK(WTOTP ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTE)
      OUTAU(IA:IE) = PACK(WTOTE ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTA)
      OUTAU(IA:IE) = PACK(ETOTA ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTM)
      OUTAU(IA:IE) = PACK(ETOTM ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTI)
      OUTAU(IA:IE) = PACK(ETOTI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTPH)
      OUTAU(IA:IE) = PACK(ETOTPH ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTP)
      OUTAU(IA:IE) = PACK(ETOTP ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(XMCP)
      OUTAU(IA:IE) = PACK(XMCP  ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(FLUXT)
      OUTAU(IA:IE) = PACK(FLUXT ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(FLXFAC)
      OUTAU(IA:IE) = PACK(FLXFAC,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(EELFI)
      OUTAU(IA:IE) = PACK(EELFI ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(PTRASH)
      OUTAU(IA:IE) = PACK(PTRASH,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(ETRASH)
      OUTAU(IA:IE) = PACK(ETRASH,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(FASCL)
      OUTAU(IA:IE) = PACK(FASCL ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(FMSCL)
      OUTAU(IA:IE) = PACK(FMSCL ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(FISCL)
      OUTAU(IA:IE) = PACK(FISCL ,.TRUE.)

      IA = IE + 1
      IE = IA - 1 + SIZE(FPHSCL)
      OUTAU(IA:IE) = PACK(FPHSCL ,.TRUE.)

      IF (IE /= NOUTAU) THEN
        WRITE (iunout,*) ' ERROR IN WRITE_COUTAU '
        WRITE (iunout,*) ' NOUTAU = ',NOUTAU,' IE = ',IE
        CALL EIRENE_EXIT_OWN(1)
      END IF

      RETURN
      END SUBROUTINE EIRENE_WRITE_COUTAU


      SUBROUTINE EIRENE_READ_COUTAU (OUTAU, IUNOUT)
      REAL(DP) :: OUTAU(:)
      INTEGER, INTENT(IN) :: IUNOUT
      INTEGER :: IA, IE
      EXTERNAL :: EIRENE_EXIT_OWN

      IF (SIZE(OUTAU) < NOUTAU) THEN
        WRITE (iunout,*)
     .  ' SIZE OF ARRAY OUTAU TOO SMALL IN READ_COUTAU '
        WRITE (iunout,*) ' OUTAU(',NOUTAU,') REQUIRED '
        CALL EIRENE_EXIT_OWN(1)
      END IF

      IA = 1
      IE = IA - 1 + SIZE(PDENAI)
      PDENAI = RESHAPE(OUTAU(IA:IE),SHAPE(PDENAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PDENMI)
      PDENMI = RESHAPE(OUTAU(IA:IE),SHAPE(PDENMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PDENII)
      PDENII = RESHAPE(OUTAU(IA:IE),SHAPE(PDENII))

      IA = IE + 1
      IE = IA - 1 + SIZE(PDENPHI)
      PDENPHI = RESHAPE(OUTAU(IA:IE),SHAPE(PDENPHI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EDENAI)
      EDENAI = RESHAPE(OUTAU(IA:IE),SHAPE(EDENAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EDENMI)
      EDENMI = RESHAPE(OUTAU(IA:IE),SHAPE(EDENMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EDENII)
      EDENII = RESHAPE(OUTAU(IA:IE),SHAPE(EDENII))

      IA = IE + 1
      IE = IA - 1 + SIZE(EDENPHI)
      EDENPHI = RESHAPE(OUTAU(IA:IE),SHAPE(EDENPHI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PAELI)
      PAELI  = RESHAPE(OUTAU(IA:IE),SHAPE(PAELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PAATI)
      PAATI  = RESHAPE(OUTAU(IA:IE),SHAPE(PAATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PAMLI)
      PAMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PAMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PAIOI)
      PAIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(PAIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PAPHTI)
      PAPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(PAPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PAPLI)
      PAPLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PAPLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PMELI)
      PMELI  = RESHAPE(OUTAU(IA:IE),SHAPE(PMELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PMATI)
      PMATI  = RESHAPE(OUTAU(IA:IE),SHAPE(PMATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PMMLI)
      PMMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PMMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PMIOI)
      PMIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(PMIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PMPHTI)
      PMPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(PMPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PMPLI)
      PMPLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PMPLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PIELI)
      PIELI  = RESHAPE(OUTAU(IA:IE),SHAPE(PIELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PIATI)
      PIATI  = RESHAPE(OUTAU(IA:IE),SHAPE(PIATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PIMLI)
      PIMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PIMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PIIOI)
      PIIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(PIIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PIPHTI)
      PIPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(PIPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PIPLI)
      PIPLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PIPLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHELI)
      PPHELI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPHELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHATI)
      PPHATI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPHATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHMLI)
      PPHMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPHMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHIOI)
      PPHIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPHIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHPHTI)
      PPHPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(PPHPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPHPLI)
      PPHPLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPHPLI ))




      IA = IE + 1
      IE = IA - 1 + SIZE(EAELI)
      EAELI  = RESHAPE(OUTAU(IA:IE),SHAPE(EAELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EAATI)
      EAATI  = RESHAPE(OUTAU(IA:IE),SHAPE(EAATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EAMLI)
      EAMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(EAMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EAIOI)
      EAIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(EAIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EAPHTI)
      EAPHTI  = RESHAPE(OUTAU(IA:IE),SHAPE(EAPHTI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EAPLI)
      EAPLI  = RESHAPE(OUTAU(IA:IE),SHAPE(EAPLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EMELI)
      EMELI  = RESHAPE(OUTAU(IA:IE),SHAPE(EMELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EMATI)
      EMATI  = RESHAPE(OUTAU(IA:IE),SHAPE(EMATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EMMLI)
      EMMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(EMMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EMIOI)
      EMIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(EMIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EMPHTI)
      EMPHTI  = RESHAPE(OUTAU(IA:IE),SHAPE(EMPHTI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EMPLI)
      EMPLI  = RESHAPE(OUTAU(IA:IE),SHAPE(EMPLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EIELI)
      EIELI  = RESHAPE(OUTAU(IA:IE),SHAPE(EIELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EIATI)
      EIATI  = RESHAPE(OUTAU(IA:IE),SHAPE(EIATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EIMLI)
      EIMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(EIMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EIIOI)
      EIIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(EIIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EIPHTI)
      EIPHTI  = RESHAPE(OUTAU(IA:IE),SHAPE(EIPHTI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EIPLI)
      EIPLI  = RESHAPE(OUTAU(IA:IE),SHAPE(EIPLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHELI)
      EPHELI  = RESHAPE(OUTAU(IA:IE),SHAPE(EPHELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHATI)
      EPHATI  = RESHAPE(OUTAU(IA:IE),SHAPE(EPHATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHMLI)
      EPHMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(EPHMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHIOI)
      EPHIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(EPHIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHPHTI)
      EPHPHTI  = RESHAPE(OUTAU(IA:IE),SHAPE(EPHPHTI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPHPLI)
      EPHPLI  = RESHAPE(OUTAU(IA:IE),SHAPE(EPHPLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(ADDVI)
      ADDVI  = RESHAPE(OUTAU(IA:IE),SHAPE(ADDVI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(COLVI)
      COLVI  = RESHAPE(OUTAU(IA:IE),SHAPE(COLVI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(SNAPVI)
      SNAPVI = RESHAPE(OUTAU(IA:IE),SHAPE(SNAPVI))

      IA = IE + 1
      IE = IA - 1 + SIZE(COPVI)
      COPVI  = RESHAPE(OUTAU(IA:IE),SHAPE(COPVI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(BGKVI)
      BGKVI  = RESHAPE(OUTAU(IA:IE),SHAPE(BGKVI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(ALGVI)
      ALGVI  = RESHAPE(OUTAU(IA:IE),SHAPE(ALGVI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PGENAI)
      PGENAI = RESHAPE(OUTAU(IA:IE),SHAPE(PGENAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PGENMI)
      PGENMI = RESHAPE(OUTAU(IA:IE),SHAPE(PGENMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PGENII)
      PGENII = RESHAPE(OUTAU(IA:IE),SHAPE(PGENII))

      IA = IE + 1
      IE = IA - 1 + SIZE(PGENPHI)
      PGENPHI = RESHAPE(OUTAU(IA:IE),SHAPE(PGENPHI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EGENAI)
      EGENAI = RESHAPE(OUTAU(IA:IE),SHAPE(EGENAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EGENMI)
      EGENMI = RESHAPE(OUTAU(IA:IE),SHAPE(EGENMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EGENII)
      EGENII = RESHAPE(OUTAU(IA:IE),SHAPE(EGENII))

      IA = IE + 1
      IE = IA - 1 + SIZE(EGENPHI)
      EGENPHI = RESHAPE(OUTAU(IA:IE),SHAPE(EGENPHI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VGENAI)
      VGENAI = RESHAPE(OUTAU(IA:IE),SHAPE(VGENAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VGENMI)
      VGENMI = RESHAPE(OUTAU(IA:IE),SHAPE(VGENMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VGENII)
      VGENII = RESHAPE(OUTAU(IA:IE),SHAPE(VGENII))

      IA = IE + 1
      IE = IA - 1 + SIZE(VGENPHI)
      VGENPHI = RESHAPE(OUTAU(IA:IE),SHAPE(VGENPHI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPATI)
      PPATI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPMLI)
      PPMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPIOI)
      PPIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPPHTI)
      PPPHTI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPPHTI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PPPLI)
      PPPLI  = RESHAPE(OUTAU(IA:IE),SHAPE(PPPLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPATI)
      EPATI  = RESHAPE(OUTAU(IA:IE),SHAPE(EPATI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPMLI)
      EPMLI  = RESHAPE(OUTAU(IA:IE),SHAPE(EPMLI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPIOI)
      EPIOI  = RESHAPE(OUTAU(IA:IE),SHAPE(EPIOI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPPHTI)
      EPPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(EPPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EPPLI)
      EPPLI = RESHAPE(OUTAU(IA:IE),SHAPE(EPPLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VXDENAI)
      VXDENAI = RESHAPE(OUTAU(IA:IE),SHAPE(VXDENAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VXDENMI)
      VXDENMI = RESHAPE(OUTAU(IA:IE),SHAPE(VXDENMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VXDENII)
      VXDENII = RESHAPE(OUTAU(IA:IE),SHAPE(VXDENII))

      IA = IE + 1
      IE = IA - 1 + SIZE(VXDENPHI)
      VXDENPHI = RESHAPE(OUTAU(IA:IE),SHAPE(VXDENPHI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VYDENAI)
      VYDENAI = RESHAPE(OUTAU(IA:IE),SHAPE(VYDENAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VYDENMI)
      VYDENMI = RESHAPE(OUTAU(IA:IE),SHAPE(VYDENMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VYDENII)
      VYDENII = RESHAPE(OUTAU(IA:IE),SHAPE(VYDENII))

      IA = IE + 1
      IE = IA - 1 + SIZE(VYDENPHI)
      VYDENPHI = RESHAPE(OUTAU(IA:IE),SHAPE(VYDENPHI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VZDENAI)
      VZDENAI = RESHAPE(OUTAU(IA:IE),SHAPE(VZDENAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VZDENMI)
      VZDENMI = RESHAPE(OUTAU(IA:IE),SHAPE(VZDENMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(VZDENII)
      VZDENII = RESHAPE(OUTAU(IA:IE),SHAPE(VZDENII))

      IA = IE + 1
      IE = IA - 1 + SIZE(VZDENPHI)
      VZDENPHI = RESHAPE(OUTAU(IA:IE),SHAPE(VZDENPHI))

      IA = IE + 1
      IE = IA - 1 + SIZE(MAPLI_VEC)
      MAPLI_VEC = RESHAPE(OUTAU(IA:IE),SHAPE(MAPLI_VEC))

      IA = IE + 1
      IE = IA - 1 + SIZE(MMPLI_VEC)
      MMPLI_VEC = RESHAPE(OUTAU(IA:IE),SHAPE(MMPLI_VEC))

      IA = IE + 1
      IE = IA - 1 + SIZE(MIPLI_VEC)
      MIPLI_VEC = RESHAPE(OUTAU(IA:IE),SHAPE(MIPLI_VEC))

      IA = IE + 1
      IE = IA - 1 + SIZE(MPHPLI_VEC)
      MPHPLI_VEC = RESHAPE(OUTAU(IA:IE),SHAPE(MPHPLI_VEC))



      IA = IE + 1
      IE = IA - 1 + SIZE(POTATI)
      POTATI = RESHAPE(OUTAU(IA:IE),SHAPE(POTATI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFAAI)
      PRFAAI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFAAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFMAI)
      PRFMAI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFMAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFIAI)
      PRFIAI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFIAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPHAI)
      PRFPHAI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFPHAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPAI)
      PRFPAI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFPAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(POTMLI)
      POTMLI = RESHAPE(OUTAU(IA:IE),SHAPE(POTMLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFAMI)
      PRFAMI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFAMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFMMI)
      PRFMMI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFMMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFIMI)
      PRFIMI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFIMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPHMI)
      PRFPHMI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFPHMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPMI)
      PRFPMI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFPMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(POTIOI)
      POTIOI = RESHAPE(OUTAU(IA:IE),SHAPE(POTIOI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFAII)
      PRFAII = RESHAPE(OUTAU(IA:IE),SHAPE(PRFAII))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFMII)
      PRFMII = RESHAPE(OUTAU(IA:IE),SHAPE(PRFMII))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFIII)
      PRFIII = RESHAPE(OUTAU(IA:IE),SHAPE(PRFIII))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPHII)
      PRFPHII = RESHAPE(OUTAU(IA:IE),SHAPE(PRFPHII))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPII)
      PRFPII = RESHAPE(OUTAU(IA:IE),SHAPE(PRFPII))

      IA = IE + 1
      IE = IA - 1 + SIZE(POTPHTI)
      POTPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(POTPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFAPHTI)
      PRFAPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFAPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFMPHTI)
      PRFMPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFMPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFIPHTI)
      PRFIPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFIPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPHPHTI)
      PRFPHPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFPHPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(PRFPPHTI)
      PRFPPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(PRFPPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(POTPLI)
      POTPLI = RESHAPE(OUTAU(IA:IE),SHAPE(POTPLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTATI)
      EOTATI = RESHAPE(OUTAU(IA:IE),SHAPE(EOTATI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFAAI)
      ERFAAI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFAAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFMAI)
      ERFMAI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFMAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFIAI)
      ERFIAI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFIAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPHAI)
      ERFPHAI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFPHAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPAI)
      ERFPAI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFPAI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTMLI)
      EOTMLI = RESHAPE(OUTAU(IA:IE),SHAPE(EOTMLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFAMI)
      ERFAMI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFAMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFMMI)
      ERFMMI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFMMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFIMI)
      ERFIMI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFIMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPHMI)
      ERFPHMI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFPHMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPMI)
      ERFPMI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFPMI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTIOI)
      EOTIOI = RESHAPE(OUTAU(IA:IE),SHAPE(EOTIOI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFAII)
      ERFAII = RESHAPE(OUTAU(IA:IE),SHAPE(ERFAII))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFMII)
      ERFMII = RESHAPE(OUTAU(IA:IE),SHAPE(ERFMII))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFIII)
      ERFIII = RESHAPE(OUTAU(IA:IE),SHAPE(ERFIII))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPHII)
      ERFPHII = RESHAPE(OUTAU(IA:IE),SHAPE(ERFPHII))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPII)
      ERFPII = RESHAPE(OUTAU(IA:IE),SHAPE(ERFPII))

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTPHTI)
      EOTPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(EOTPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFAPHTI)
      ERFAPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFAPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFMPHTI)
      ERFMPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFMPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFIPHTI)
      ERFIPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFIPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPHPHTI)
      ERFPHPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFPHPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ERFPPHTI)
      ERFPPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(ERFPPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(EOTPLI)
      EOTPLI = RESHAPE(OUTAU(IA:IE),SHAPE(EOTPLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAATI)
      SPTAATI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTAATI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMATI)
      SPTMATI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTMATI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIATI)
      SPTIATI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTIATI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHATI)
      SPTPHATI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPHATI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPATI)
      SPTPATI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPATI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAMLI)
      SPTAMLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTAMLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMMLI)
      SPTMMLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTMMLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIMLI)
      SPTIMLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTIMLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHMLI)
      SPTPHMLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPHMLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPMLI)
      SPTPMLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPMLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAIOI)
      SPTAIOI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTAIOI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMIOI)
      SPTMIOI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTMIOI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIIOI)
      SPTIIOI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTIIOI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHIOI)
      SPTPHIOI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPHIOI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPIOI)
      SPTPIOI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPIOI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAPHTI)
      SPTAPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTAPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMPHTI)
      SPTMPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTMPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIPHTI)
      SPTIPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTIPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHPHTI)
      SPTPHPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPHPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPPHTI)
      SPTPPHTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPPHTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTAPLI)
      SPTAPLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTAPLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMPLI)
      SPTMPLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTMPLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTIPLI)
      SPTIPLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTIPLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHPLI)
      SPTPHPLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPHPLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPPLI)
      SPTPPLI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPPLI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTATTI)
      SPTATTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTATTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTMTTI)
      SPTMTTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTMTTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTITTI)
      SPTITTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTITTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPHTTI)
      SPTPHTTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPHTTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTPLTTI)
      SPTPLTTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTPLTTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPTTTI)
      SPTTTI = RESHAPE(OUTAU(IA:IE),SHAPE(SPTTTI))

      IA = IE + 1
      IE = IA - 1 + SIZE(ADDSI)
      ADDSI  = RESHAPE(OUTAU(IA:IE),SHAPE(ADDSI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(ALGSI)
      ALGSI  = RESHAPE(OUTAU(IA:IE),SHAPE(ALGSI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(SPUMPI)
      SPUMPI = RESHAPE(OUTAU(IA:IE),SHAPE(SPUMPI))

      IA = IE + 1
      IE = IA - 1 + SIZE(RAELI)
      RAELI  = RESHAPE(OUTAU(IA:IE),SHAPE(RAELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(RMELI)
      RMELI  = RESHAPE(OUTAU(IA:IE),SHAPE(RMELI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(RIELI)
      RIELI  = RESHAPE(OUTAU(IA:IE),SHAPE(RIELI ))


      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTA)
      WTOTA  = RESHAPE(OUTAU(IA:IE),SHAPE(WTOTA ))

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTM)
      WTOTM  = RESHAPE(OUTAU(IA:IE),SHAPE(WTOTM ))

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTI)
      WTOTI  = RESHAPE(OUTAU(IA:IE),SHAPE(WTOTI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTPH)
      WTOTPH = RESHAPE(OUTAU(IA:IE),SHAPE(WTOTPH))

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTP)
      WTOTP  = RESHAPE(OUTAU(IA:IE),SHAPE(WTOTP ))

      IA = IE + 1
      IE = IA - 1 + SIZE(WTOTE)
      WTOTE  = RESHAPE(OUTAU(IA:IE),SHAPE(WTOTE ))

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTA)
      ETOTA  = RESHAPE(OUTAU(IA:IE),SHAPE(ETOTA ))

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTM)
      ETOTM  = RESHAPE(OUTAU(IA:IE),SHAPE(ETOTM ))

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTI)
      ETOTI  = RESHAPE(OUTAU(IA:IE),SHAPE(ETOTI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTPH)
      ETOTPH = RESHAPE(OUTAU(IA:IE),SHAPE(ETOTPH))

      IA = IE + 1
      IE = IA - 1 + SIZE(ETOTP)
      ETOTP  = RESHAPE(OUTAU(IA:IE),SHAPE(ETOTP ))

      IA = IE + 1
      IE = IA - 1 + SIZE(XMCP)
      XMCP   = RESHAPE(OUTAU(IA:IE),SHAPE(XMCP  ))

      IA = IE + 1
      IE = IA - 1 + SIZE(FLUXT)
      FLUXT  = RESHAPE(OUTAU(IA:IE),SHAPE(FLUXT ))

      IA = IE + 1
      IE = IA - 1 + SIZE(FLXFAC)
      FLXFAC = RESHAPE(OUTAU(IA:IE),SHAPE(FLXFAC))

      IA = IE + 1
      IE = IA - 1 + SIZE(EELFI)
      EELFI  = RESHAPE(OUTAU(IA:IE),SHAPE(EELFI ))

      IA = IE + 1
      IE = IA - 1 + SIZE(PTRASH)
      PTRASH = RESHAPE(OUTAU(IA:IE),SHAPE(PTRASH))

      IA = IE + 1
      IE = IA - 1 + SIZE(ETRASH)
      ETRASH = RESHAPE(OUTAU(IA:IE),SHAPE(ETRASH))

      IA = IE + 1
      IE = IA - 1 + SIZE(FASCL)
      FASCL  = RESHAPE(OUTAU(IA:IE),SHAPE(FASCL ))

      IA = IE + 1
      IE = IA - 1 + SIZE(FMSCL)
      FMSCL  = RESHAPE(OUTAU(IA:IE),SHAPE(FMSCL ))

      IA = IE + 1
      IE = IA - 1 + SIZE(FISCL)
      FISCL  = RESHAPE(OUTAU(IA:IE),SHAPE(FISCL ))

      IA = IE + 1
      IE = IA - 1 + SIZE(FPHSCL)
      FPHSCL  = RESHAPE(OUTAU(IA:IE),SHAPE(FPHSCL ))

      IF (IE /= NOUTAU) THEN
        WRITE (iunout,*) ' ERROR IN READ_COUTAU '
        WRITE (iunout,*) ' NOUTAU = ',NOUTAU,' IE = ',IE
        CALL EIRENE_EXIT_OWN(1)
      END IF

      RETURN
      END SUBROUTINE EIRENE_READ_COUTAU


      SUBROUTINE EIRENE_FETCH_OUTAUI (OUTAUI,ITAL,ISP,ISTRA,IUNOUT)
cdr  return OUTAU = total (volume integrated tally),
cdr    for tally ITAL, species ISP, stratum ISTRA
      REAL(DP) :: OUTAU
      INTEGER, INTENT(IN) :: ITAL, ISP, ISTRA, IUNOUT
      REAL(DP), INTENT(OUT) :: OUTAUI

      outaui=0._dp
      IF (.not. ALLOCATED(PDENAI)) RETURN

      SELECT CASE (ITAL)
      CASE (1)
        OUTAU = PDENAI(ISP,ISTRA)
      CASE (2)
        OUTAU = PDENMI(ISP,ISTRA)
      CASE (3)
        OUTAU = PDENII(ISP,ISTRA)
      CASE (4)
        OUTAU = PDENPHI(ISP,ISTRA)
      CASE (5)
        OUTAU = EDENAI(ISP,ISTRA)
      CASE (6)
        OUTAU = EDENMI(ISP,ISTRA)
      CASE (7)
        OUTAU = EDENII(ISP,ISTRA)
      CASE (8)
        OUTAU = EDENPHI(ISP,ISTRA)
      CASE (9)
        OUTAU = PAELI(ISTRA)
      CASE (10)
        OUTAU = PAATI(ISP,ISTRA)
      CASE (11)
        OUTAU = PAMLI(ISP,ISTRA)
      CASE (12)
        OUTAU = PAIOI(ISP,ISTRA)
      CASE (13)
        OUTAU = PAPHTI(ISP,ISTRA)
      CASE (14)
        OUTAU = PAPLI(ISP,ISTRA)
      CASE (15)
        OUTAU = PMELI(ISTRA)
      CASE (16)
        OUTAU = PMATI(ISP,ISTRA)
      CASE (17)
        OUTAU = PMMLI(ISP,ISTRA)
      CASE (18)
        OUTAU = PMIOI(ISP,ISTRA)
      CASE (19)
        OUTAU = PMPHTI(ISP,ISTRA)
      CASE (20)
        OUTAU = PMPLI(ISP,ISTRA)
      CASE (21)
        OUTAU = PIELI(ISTRA)
      CASE (22)
        OUTAU = PIATI(ISP,ISTRA)
      CASE (23)
        OUTAU = PIMLI(ISP,ISTRA)
      CASE (24)
        OUTAU = PIIOI(ISP,ISTRA)
      CASE (25)
        OUTAU = PIPHTI(ISP,ISTRA)
      CASE (26)
        OUTAU = PIPLI(ISP,ISTRA)
      CASE (27)
        OUTAU = PPHELI(ISTRA)
      CASE (28)
        OUTAU = PPHATI(ISP,ISTRA)
      CASE (29)
        OUTAU = PPHMLI(ISP,ISTRA)
      CASE (30)
        OUTAU = PPHIOI(ISP,ISTRA)
      CASE (31)
        OUTAU = PPHPHTI(ISP,ISTRA)
      CASE (32)
        OUTAU = PPHPLI(ISP,ISTRA)
      CASE (33)
        OUTAU = EAELI(ISTRA)
      CASE (34)
        OUTAU = EAATI(ISTRA)
      CASE (35)
        OUTAU = EAMLI(ISTRA)
      CASE (36)
        OUTAU = EAIOI(ISTRA)
      CASE (37)
        OUTAU = EAPHTI(ISTRA)
      CASE (38)
        OUTAU = EAPLI(ISP,ISTRA)
      CASE (39)
        OUTAU = EMELI(ISTRA)
      CASE (40)
        OUTAU = EMATI(ISTRA)
      CASE (41)
        OUTAU = EMMLI(ISTRA)
      CASE (42)
        OUTAU = EMIOI(ISTRA)
      CASE (43)
        OUTAU = EMPHTI(ISTRA)
      CASE (44)
        OUTAU = EMPLI(ISP,ISTRA)
      CASE (45)
        OUTAU = EIELI(ISTRA)
      CASE (46)
        OUTAU = EIATI(ISTRA)
      CASE (47)
        OUTAU = EIMLI(ISTRA)
      CASE (48)
        OUTAU = EIIOI(ISTRA)
      CASE (49)
        OUTAU = EIPHTI(ISTRA)
      CASE (50)
        OUTAU = EIPLI(ISP,ISTRA)
      CASE (51)
        OUTAU = EPHELI(ISTRA)
      CASE (52)
        OUTAU = EPHATI(ISTRA)
      CASE (53)
        OUTAU = EPHMLI(ISTRA)
      CASE (54)
        OUTAU = EPHIOI(ISTRA)
      CASE (55)
        OUTAU = EPHPHTI(ISTRA)
      CASE (56)
        OUTAU = EPHPLI(ISP,ISTRA)
      CASE (57)
        OUTAU = ADDVI(ISP,ISTRA)
      CASE (58)
        OUTAU = COLVI(ISP,ISTRA)
      CASE (59)
        OUTAU = SNAPVI(ISP,ISTRA)
      CASE (60)
        OUTAU = COPVI(ISP,ISTRA)
      CASE (61)
        OUTAU = BGKVI(ISP,ISTRA)
      CASE (62)
        OUTAU = ALGVI(ISP,ISTRA)
      CASE (63)
        OUTAU = PGENAI(ISP,ISTRA)
      CASE (64)
        OUTAU = PGENMI(ISP,ISTRA)
      CASE (65)
        OUTAU = PGENII(ISP,ISTRA)
      CASE (66)
        OUTAU = PGENPHI(ISP,ISTRA)
      CASE (67)
        OUTAU = EGENAI(ISP,ISTRA)
      CASE (68)
        OUTAU = EGENMI(ISP,ISTRA)
      CASE (69)
        OUTAU = EGENII(ISP,ISTRA)
      CASE (70)
        OUTAU = EGENPHI(ISP,ISTRA)
      CASE (71)
        OUTAU = VGENAI(ISP,ISTRA)
      CASE (72)
        OUTAU = VGENMI(ISP,ISTRA)
      CASE (73)
        OUTAU = VGENII(ISP,ISTRA)
      CASE (74)
        OUTAU = VGENPHI(ISP,ISTRA)
      CASE (75)
        OUTAU = PPATI(ISP,ISTRA)
      CASE (76)
        OUTAU = PPMLI(ISP,ISTRA)
      CASE (77)
        OUTAU = PPIOI(ISP,ISTRA)
      CASE (78)
        OUTAU = PPPHTI(ISP,ISTRA)
      CASE (79)
        OUTAU = PPPLI(ISP,ISTRA)
      CASE (80)
        OUTAU = EPATI(ISTRA)
      CASE (81)
        OUTAU = EPMLI(ISTRA)
      CASE (82)
        OUTAU = EPIOI(ISTRA)
      CASE (83)
        OUTAU = EPPHTI(ISTRA)
      CASE (84)
        OUTAU = EPPLI(ISP,ISTRA)
      CASE (85)
        OUTAU = VXDENAI(ISP,ISTRA)
      CASE (86)
        OUTAU = VXDENMI(ISP,ISTRA)
      CASE (87)
        OUTAU = VXDENII(ISP,ISTRA)
      CASE (88)
        OUTAU = VXDENPHI(ISP,ISTRA)
      CASE (89)
        OUTAU = VYDENAI(ISP,ISTRA)
      CASE (90)
        OUTAU = VYDENMI(ISP,ISTRA)
      CASE (91)
        OUTAU = VYDENII(ISP,ISTRA)
      CASE (92)
        OUTAU = VYDENPHI(ISP,ISTRA)
      CASE (93)
        OUTAU = VZDENAI(ISP,ISTRA)
      CASE (94)
        OUTAU = VZDENMI(ISP,ISTRA)
      CASE (95)
        OUTAU = VZDENII(ISP,ISTRA)
      CASE (96)
        OUTAU = VZDENPHI(ISP,ISTRA)
cdr Next:
cdr  vectorial tallies: sum over IPLS species for each vector component
cdr                     is on ISP=-1,-2,-3.
cdr                     ISP=0 is for sum over species AND vector components
      CASE (97)
        OUTAU = MAPLI_VEC(ISP,ISTRA)
      CASE (98)
        OUTAU = MMPLI_VEC(ISP,ISTRA)
      CASE (99)
        OUTAU = MIPLI_VEC(ISP,ISTRA)
      CASE (100)
        OUTAU = MPHPLI_VEC(ISP,ISTRA)
      CASE (101)
        OUTAU = RAELI(ISP,ISTRA)
      CASE (102)
        OUTAU = RMELI(ISP,ISTRA)
      CASE (103)
        OUTAU = RIELI(ISP,ISTRA)
      CASE DEFAULT
        WRITE (iunout,*) ' WRONG TALLY NUMBER IN FETCH_OUTAUI '
        WRITE (iunout,*) ' 0 RETURNED '
        OUTAU = 0
      END SELECT
      OUTAUI=OUTAU

      RETURN
      END SUBROUTINE EIRENE_FETCH_OUTAUI


      SUBROUTINE EIRENE_FETCH_OUTAU (OUTAU,ITAL,ISP,ISTR,IUNOUT)
cdr return OUTAU = eirene tally,
cdr for tally ITAL, species ISP, stratum ISTRA=IESTR?
cdr The entire tally is returned with ic1=1, ic2=NSBOX_TAL 
      USE EIRMOD_CESTIM 
      USE EIRMOD_CGRID, ONLY : NSBOX, NSBOX_TAL
      USE EIRMOD_CSPEI, ONLY : IESTR
      USE EIRMOD_CGEOM, ONLY : NCLTAL

      INTEGER, INTENT(IN) :: ITAL, ISP, ISTR, IUNOUT
      REAL(DP), INTENT(OUT) :: OUTAU(NRAD)
      INTEGER :: I, K, NF, IINDEX, ICELL

      OUTAU=0._DP

CDR  CHECK: IESTR=ISTR ??
      IF (IESTR .NE. ISTR) THEN
        WRITE (IUNOUT,*) 'STRATUM NO. CONFLICT IN FETCH_OUTAU '
        WRITE (IUNOUT,*) 'RETURN VALUE = 0.0'
        RETURN
      ENDIF
      
      IF (.NOT.ALLOCATED (ESTIMV)) RETURN 

      NF= NFIRST(ITAL)
      K = ISP


      IF (NSBOX_TAL /= NSBOX) THEN
        DO I=1,NSBOX
          ICELL = NCLTAL(I)
          IF (ICELL > 0) THEN
            OUTAU(I)=ESTIMV(NADDV(ITAL)+K,ICELL)
          ELSE
            OUTAU(I)=0._DP
          END IF
        END DO
      ELSE
        DO I=1,NSBOX_TAL
          OUTAU(I)=ESTIMV(NADDV(ITAL)+K,I)
        END DO
      END IF

      RETURN
      END SUBROUTINE EIRENE_FETCH_OUTAU


      SUBROUTINE EIRENE_BROADCAST_COUTAU(ME)
      USE EIRMOD_MPI
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

      IF (ME /= 0) CALL EIRENE_ALLOC_COUTAU

      CALL MPI_BCAST (NADDI,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFRSTI,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFSTVI,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NEXTVI,NTALV,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NDDWI,NTALS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFRTWI,NTALS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NFSTWI,NTALS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      CALL MPI_BCAST (NEXTWI,NTALS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROADCAST_COUTAU

      END MODULE EIRMOD_COUTAU
