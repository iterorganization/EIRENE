      MODULE EIRMOD_UPTUSR

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CUPD
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEZ
      USE EIRMOD_CGRID
      USE EIRMOD_CLOGAU
      USE EIRMOD_CCONA
      USE EIRMOD_CPOLYG
      USE EIRMOD_CZT1

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_UPTUSR, EIRENE_uptusr_reinit

c      REAL(DP), ALLOCATABLE, SAVE :: CNDYNA(:),CNDYNP(:)
CDR
      REAL(DP), ALLOCATABLE, SAVE :: VPX(:),VPY(:),VRX(:),VRY(:)

      INTEGER, SAVE :: IFIRST, IA0, IA1, IA2, IA3, NA4, INDEXM, INDEXF
      DATA IFIRST/0/

      integer, save :: num=0,eirene_nbirth,eirene_njetto
      character(len=256), save :: eirene_fbirth,eirene_ftransfer,
     &     eirene_fstoreneutflux

      CONTAINS

      SUBROUTINE EIRENE_UPTUSR(XSTOR2,XSTORV2,WV,IFLAG)
C
C  USER-SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME-AVERAGED
C
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 23/08/2006                VPX, VPY, VRX, VRY changed to       |
C|                           ALLOCATABLE, SAVE to speed up       |
C|                           subroutine call (save time in       |
C|                           storage allocation)                 |
C| 16/07/2010   D.Harting    Added two variables to eirene_user  |
C|                           namelist for use of flux dependency |
C|                           in chemical sputtering.             |
C+---------------------------------------------------------------+
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                         XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
CDR
      INTEGER :: IAT, IPL, I, IR, IP, IRD
csw
      real(dp) :: dist,wtr
      integer :: iaei,irei,iacx,ircx
      real(dp) :: eirene_phi_offsets(9)
      integer :: eirene_wallFluxModel ! calculation of wall fluxes for chemical sputtering
c                = 0: no wall fluxes are used (old edge2d model)
c                = 1: only ion fluxes are used
c                = 2: ion fluxes and neutral fluxes from last eirene iteration are used
c                = 3: ion and neutral fluxes are used, and EIRENE is iterated to give
c                     converged neutral fluxes.
      logical :: eirene_use_elstepdat_bug
      namelist /eirene_user/eirene_nbirth,eirene_njetto,
     .                      eirene_fbirth,eirene_ftransfer,
     .                      eirene_phi_offsets,
     .                      eirene_fstoreneutflux,
     .                      eirene_wallFluxModel,
     .                      eirene_use_elstepdat_bug
csw

      IF (IFIRST.EQ.0) THEN
        IFIRST=1
c        ALLOCATE (CNDYNA(NATM))
c        ALLOCATE (CNDYNP(NPLS))
c        DO IAT=1,NATMI
c          CNDYNA(IAT)=1.D3*AMUA*RMASSA(IAT)
c        END DO
c        DO IPL=1,NPLSI
c          CNDYNP(IPL)=1.D3*AMUA*RMASSP(IPL)
c        END DO
C
CDR
CDR  PROVIDE A RADIAL UNIT VECTOR PER CELL
CDR  VPX,VPY,  NEEDED FOR PROJECTING PARTICLE VELOCITIES
CDR  SAME FOR POLOIDAL UNIT VECTOR VRX,VRY
C
        ALLOCATE (VPX(NRAD))
        ALLOCATE (VPY(NRAD))
        ALLOCATE (VRX(NRAD))
        ALLOCATE (VRY(NRAD))
        DO I=1,NRAD
          VPX(I)=0.
          VPY(I)=0.
          VRX(I)=0.
          VRY(I)=0.
        END DO
        DO IR=1,NR1STM
          DO IP=1,NP2NDM
            IRD=IR+(IP-1)*NR1P2
            VPX(IRD)=PLNX(IR,IP)
            VPY(IRD)=PLNY(IR,IP)
            VRX(IRD)=PPLNX(IR,IP)
            VRY(IRD)=PPLNY(IR,IP)
          END DO
        END DO
!pb  MOD_ADDV is no incremental value. It is a flag indicating whether all the
!pb  rates used for emissivity lines are to be stored or whether storage saving
!pb  mode ist to be used, only storing the rates for the latest used line
!pb     IA0=MOD_ADDV
        IA0=0
        IA1=IA0+NATMI+NMOLI
        IA2=IA1+NATMI+NMOLI
        IA3=IA2+NATMI+NMOLI
        NA4=IA3+NATMI+NMOLI
        INDEXM=NPLSI
        INDEXF=2*NPLSI
csw
        eirene_njetto=0
cdmh
        eirene_fstoreneutflux = 'eirene.chemFluxDep'
        eirene_wallFluxModel = 1
cdmh
        open(unit=9998,file='eirene_user.namelist')
        read(9998,eirene_user)
        close(9998)
csw
      ENDIF
csw
csw 24jan08 additional tallies for JETTO cold neutral coupling
csw
      if(       eirene_njetto .gt. 0 .and.
     .          nadv .eq. 8 .and. ityp.eq.1 .and. iatm.eq.1) then
        do i=1,ncou
          dist = clpd(i)
          ird = nrcell+nupc(i)*NR1P2+NBLCKA
          wtr = wv*dist
c         sources/sinks due ionisation:
          iaei = 1
          irei=lgaei(iatm,iaei)
          addv(1,ird) = addv(1,ird)+wtr*sigvei(irei)
          addv(4,ird) = addv(4,ird)+wtr*sigvei(irei)*esigei(irei,5)
          addv(5,ird) = addv(5,ird)+wtr*sigvei(irei)*esigei(irei,4)

c         sources/sinks due recombination:
c           not included

c         net sources due CX:
          iacx = 1
          ircx=lgacx(iatm,iacx,0)
          addv(8,ird) = addv(8,ird)+wtr*sigvcx(ircx)*
     .               (E0 - esigcx(ircx,1))
        enddo
      endif
      RETURN
      END SUBROUTINE EIRENE_UPTUSR


      SUBROUTINE EIRENE_uptusr_reinit
      if(ifirst .ne. 0) then
        ifirst=0
c        if(allocated(cndyna)) deallocate(cndyna)
c        if(allocated(cndynp)) deallocate(cndynp)
        if(allocated(vpx)) deallocate(vpx)
        if(allocated(vpy)) deallocate(vpy)
        if(allocated(vrx)) deallocate(vrx)
        if(allocated(vry)) deallocate(vry)
      endif
      return
      END SUBROUTINE EIRENE_UPTUSR_reinit

      END MODULE EIRMOD_UPTUSR
