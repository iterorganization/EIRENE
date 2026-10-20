      MODULE EIRMOD_UPCUSR

      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      USE EIRMOD_CGRID
      USE EIRMOD_CCONA
      use EIRMOD_CUPD

      IMPLICIT NONE
      PRIVATE

      PUBLIC :: EIRENE_UPCUSR, EIRENE_upcusr_reinit

cswx 24sep07
      logical, save :: lfirst=.false.
      integer, save :: num=0,eirene_nbirth,eirene_njetto
      real*8, save :: delang
      character(len=256), save :: eirene_fbirth,eirene_ftransfer,
     &     eirene_fstoreneutflux

      real*8, allocatable, save :: rdata(:,:)
      integer, allocatable, save :: idata(:,:)

      CONTAINS

      SUBROUTINE EIRENE_UPCUSR(WS,IND,KK)
C
C  USER-SUPPLIED COLLISION ESTIMATOR, VOLUME-AVERAGED
C
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 16/07/2010   D.Harting    Added two variables to eirene_user  |
C|                           namelist for use of flux dependency |
C|                           in chemical sputtering.             |
C+---------------------------------------------------------------+
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: IND, KK
      real*8 :: xx,yy,zz,d,pphi,alph,v0v
      real(dp) :: eirene_phi_offsets(9)
      integer :: eirene_wallFluxModel ! calculation of wall fluxes for chemical sputtering
c                = 0: no wall fluxes are used (old edge2d model)
c                = 1: only ion fluxes are used
c                = 2: ion fluxes and neutral fluxes from last eirene iteration are used
c                = 3: ion and neutral fluxes are used, and EIRENE is iterated to give
c                     converged neutral fluxes.
      logical :: eirene_use_elstepdat_bug
      integer :: i,j
      real*8 :: cosrot, sinrot, my_velx, my_vely, my_velz

c     added rotsav_torcol from old eirene version
c     was added by sven wiesen to calculate diagnostic file eirene_fbirth
c     was defined in COMPRT and set in particle-tracing.f
c     current eirene version does not initialze rotsav_torcol, so the calculation
c     here is not correct. Needs to be reinserted to particle tracing if diagnostic
c     is used again
      real(dp) :: ROTSAV_TORCOL

      namelist /eirene_user/eirene_nbirth,eirene_njetto,
     .                      eirene_fbirth,eirene_ftransfer,
     .                      eirene_phi_offsets,
     .                      eirene_fstoreneutflux,
     .                      eirene_wallFluxModel,
     .                      eirene_use_elstepdat_bug
cswx
C
C     WS=WEIGHT/SIGTOT=WEIGHT/(VEL*ZMFPI)=WEIGHT/(VEL*SIGMA,MACR.)
C
C  FOR PARTICLE DENSITY IN CELL NO. NCELL
C     COLV(1,NCELL)=COLV(1,NCELL)+WS
C
cswx 24sep07
      if(.not.lfirst) then
        lfirst=.true.
        num=0
        delang = pi2a/dble(nttram)

c default maximum number
        eirene_nbirth=2000
        eirene_fbirth = ' '
cdmh
      eirene_fstoreneutflux = 'eirene.chemFluxDep'
      eirene_wallFluxModel = 1
cdmh

        open(unit=9998,file='eirene_user.namelist')
        read(9998,eirene_user)
        close(9998)
        if(eirene_nbirth .gt. 0) then
          allocate(idata(eirene_nbirth,2))
          allocate(rdata(eirene_nbirth,9))
        endif
      endif

      ROTSAV_TORCOL=0.d0
c     skip calculation as it is not correct with ROTSAV_TORCOL=0.d0
      return

      if(         eirene_nbirth .gt. 0
     .      .and. ind .eq. 1
     .      .and. num .lt. eirene_nbirth
     .      .and. weight .gt. 1.d-2
     .      .and. e0.gt.2.e4) then
        num=num+1

        d =sqrt(x0*x0 + z0*z0)
        alph = pi2a/dble(nt3rd-1)
        pphi = mod(phi+pi2a-(alph/2.d0),pi2a)
        !pphi = mod(phi+pi2a,pi2a)
        xx = cos(pphi) * d
        yy = y0
        zz = sin(pphi) * d

        alph=0.
        alph = alph - rotsav_torcol
        ! back-rotate from torcol, periodic boundary
        alph = alph +(360.d0/dble(nt3rd-1))/360d0  *pi2a
        ! add offsets from beamboxes
        alph = alph + eirene_phi_offsets(istra)/360.d0 * pi2a

        cosrot=cos(alph)
        sinrot=sin(alph)
        my_velx = velx*cosrot + velz*sinrot
        my_vely = vely
        my_velz =-velx*sinrot + velz*cosrot
c        write(598, '(a,20(1x,e13.6))') 'UPCUSR:',
c     .        rotsav_torcol, my_velx, my_vely, my_velz
        v0v=my_velx*bxin(ncell)+my_vely*byin(ncell)+my_velz*bzin(ncell)
        idata(num,1) = num
        idata(num,2) = 2
        rdata(num,1) = sqrt(xx*xx + zz*zz)/100.d0
        rdata(num,2) = pphi/degrad
        rdata(num,3) = yy/100.d0
        rdata(num,4) = -v0v ! not used in ascot anymore
        rdata(num,5) = e0
        rdata(num,6) = weight
        rdata(num,7) = my_velx*vel/100.d0
csw change ordering
        rdata(num,8) = -my_velz*vel/100.d0
        rdata(num,9) = my_vely*vel/100.d0
      endif
cswx
      RETURN
      END SUBROUTINE EIRENE_UPCUSR

      SUBROUTINE EIRENE_upcusr_reinit
      IMPLICIT NONE
      INTEGER :: I,J

      if(lfirst) then
        lfirst=.false.

        if(eirene_nbirth .gt. 0) then
          open(unit=9998,file=trim(eirene_fbirth),
     .       form='formatted',access='sequential',status='replace')
          rewind(9998)
          write(9998,'(i7)') num
          write(9998,'(a1,a4,a4,9a14)') '#','NP','NZ','R[m]',
     .                                'PHI[deg]',
     .                                'Z[m]',
     .                                'V0/V','E[eV]','W',
     .                                'VX[m/s]','VY[m/s]','VZ[m/s]'
          do i=1,num
            write(9998,'(i5,i4,20e14.5)')
     .         (idata(i,j), j=1,2),
     .         (rdata(i,j), j=1,9)
          enddo
          close(9998)
          deallocate(idata)
          deallocate(rdata)
        endif
      endif
      return
      END SUBROUTINE EIRENE_UPCUSR_reinit

      END MODULE EIRMOD_UPCUSR
