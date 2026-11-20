      MODULE EIRMOD_WNEUTRALS

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_CCOUPL
      use eirmod_CESTIM
      use eirmod_COMUSR
      use eirmod_CCONA
      use eirmod_CGEOM
      use eirmod_CGRID
      use eirmod_CLGIN
      use eirmod_CTRIG
      use eirmod_COMSIG
      use eirmod_COMSOU
      use eirmod_CPOLYG
      use eirmod_CTRCEI
     , , only: TRCINT
      use eirmod_COMPRT
     , , only: IUNOUT, ISTRA
      use eirmod_CINIT
     , , only: FORT, FORT_LC
      use eirmod_braeir
     , , only: pux, puy, pvx, pvy
      use eirmod_eirbra
      IMPLICIT NONE
      PRIVATE

      public :: b2_cell
      public :: wneutral_fluxes
      public :: eirene_wneutrals_init
      public :: eirene_wneutrals_fill
      public :: eirene_wneutrals_recv
      public :: eirene_wneutrals_save
      public :: eirene_wneutrals_send
      public :: eirene_wneutrals_clean
      public :: eirene_wneutrals_reduce
      public :: eirene_wneutrals_dealloc
      public :: eirene_wneutrals_alloc_arrays

! eirdiag.h/eirdiag.f
!c*** Volume data:
!c***    srcml   :  molecule dissociation rate
!c***    edissml :  power loss due to molecule dissociation
!c***    eneutrad:  power radiated due to neutral atoms
!c***    emolrad :  power radiated due to molecules
!c***    eionrad :  power radiated due to molecular ions
!c*** Surface data:
!c***    wldnek  :  net kinetic energy deposited by neutrals
!c***    wldnep  :  potential energy released by neutrals
!c***    wldna   :  flux of atoms impinging onto the surface
!c***    ewlda   :  their average energy
!c***    wldnm   :  flux of molecules impinging onto the surface
!c***    ewldm   :  their average energy
!c***    wldra   :  flux of reflected atoms
!c***    wldrm   :  flux of reflected molecules
!c***    wldpp   :  flux of plasma ions impinging onto the surface
!c***    wldpa   :  flux of resulting atoms
!c***    wldpm   :  flux of resulting molecules
!c***    wldpeb  :  power carried away by these atoms and molecules
!c***    wldspt  :  flux of sputtered wall material
!c***    isrftype:  surface type (iliin in Eirene)
!c***    wlarea  :  areas of the surface segments from Eirene
!c***    wlabsrp :  absorption at the surfaces
!c***               (1-recyct from Eirene)
!c***    wlpump  :  pumped flux at the surfaces
      real(DP), save, allocatable, dimension(:,:,:), public ::
     , dab2,dmb2,dib2,tab2,tmb2,tib2,rfluxa,rfluxm,refluxa,refluxm,
     , pfluxa,pfluxm,pefluxa,pefluxm,emiss,emissmol,srcml,edissml,
     , tfluxa,tfluxm,tefluxa,tefluxm   ! toroidal fluxes, IYS 20.01.2017
      real(DP), save, allocatable, dimension(:,:), public ::
     , wldnek,wldnep
      real(DP), save, allocatable, dimension (:,:,:), public ::
     , wldna,ewlda,wldnm,ewldm,wldra,wldrm,wldpp,wldpa,wldpm,wldspta,
     , wldsptm
      real(DP), save, allocatable, dimension (:,:), public ::
     , wldpeb,wldspt,wlabsrp,wlpump
      real(DP), save, allocatable, dimension (:,:,:), public ::
     , eneutrad, emolrad, eionrad
      real(DP), save, allocatable, public :: eirpump(:),   ! pumped flux
     ,                                       eirspta(:), eirsptm(:)      ! sputtered flux  !iyv 07.03.18
      real(DP), save, public :: aver_frac, aver_frac46
      real(DP), save, allocatable, dimension(:,:), public ::
     , wldnek_aver,wldnep_aver ! average fluxes, som 02.04.2019
      real(DP), save, allocatable, dimension (:,:,:), public ::
     , wldna_aver,ewlda_aver,wldnm_aver,ewldm_aver,wldra_aver,
     , wldrm_aver,wldpp_aver,wldpa_aver,wldpm_aver,wldspta_aver,
     , wldsptm_aver ! average fluxes, som 02.04.2019
      real(DP), save, allocatable, dimension (:,:), public ::
     , wldpeb_aver,wldspt_aver,wlpump_aver ! average fluxes
      real(DP), save, allocatable, dimension (:,:), public ::
     , PDENA_aver,PDENM_aver,EDENA_aver,EDENM_aver,
     , VXDENA_aver,VXDENM_aver,VYDENA_aver,VYDENM_aver,
     , VZDENA_aver,VZDENM_aver ! average fluxes, som 02.04.2018

      integer, save, allocatable, public ::
     ,  eirdiag_nds_ind(:)

      !c*** NEUTRAL FLUXES,
      !c*** SPATIALLY RESOLVED ON NON-DEFAULT STANDARD SURFACES (NDS)
      !c*** First dimension : index of atom or molecule.
      !c*** Second dimension: index of surface element,
      !                       controlled by 'eirdiag_nds_ind'
      real(DP), save, allocatable, dimension(:), public ::
     ,  wlarea, sarea_res      !area of the surface elements
      real(DP), save, allocatable, public ::
     ,  wldna_res(:,:),    !incident atom flux
     ,  wldnm_res(:,:),    !incident molecule flux
     ,  ewlda_res(:,:),    !incident energy flux of atoms
     ,  ewldm_res(:,:),    !incident energy flux of molecules
     ,  ewldt_res(:),      !total energy wall load from Eirene particles
     ,  ewldea_res(:,:),   !emitted energy flux of atoms
     ,  ewldem_res(:,:),   !net incident energy flux due of molecules
     ,  ewldrp_res(:),     !kinetic energy of reflected neutrals
                           !originated from ions
     ,  ewldmr_res(:,:),   !energy due to recombination of atoms and
                           !atomic ions into molecules
     ,  wldspt_res(:),     !flux of sputtered wall material
     ,  wldspta_res(:,:),  !sputtered flux
                           !for each type of emitted atom
     ,  wldsptm_res(:,:),  !sputtered flux
                           !for each type of emitted molecule
     ,  wlpump_res(:,:)    !pumped flux

      !c*** dissociation energy of the hydrogen molecule
      real(DP), parameter, public :: diss_pot_H2=4.48_DP

      ! Integrals over volume.
      ! _INT: WHOLE GRID; _INT_B2: B2 grid.
      real(DP), save, allocatable, dimension (:,:), public ::
     &  PDENA_INT, PDENA_INT_B2, !Total number of particles
     &  PDENM_INT, PDENM_INT_B2,
     &  PDENI_INT, PDENI_INT_B2,
     &  EDENA_INT, EDENA_INT_B2, !Total energy, Joules
     &  EDENM_INT, EDENM_INT_B2,
     &  EDENI_INT, EDENI_INT_B2

      integer, save, allocatable, public :: isrftype(:)
      logical, save, allocatable, public :: amark(:,:)
      character(8), save, allocatable, public :: eirtxt(:)

      ! wneutral globals
      integer, save :: ifirst_wneutral=0
      integer, save :: ia0,ia1,ia2,ia3
      real(DP), save :: hlp_cnv

      !pb volumes of B2.5 cells
      real(DP), save, allocatable, public :: volcel(:)

      INTEGER, SAVE, PUBLIC ::
     &  nr1tal_save, np2tal_save, nt3tal_save,
     &  nsbox_tal_save, nsurf_tal_save, nradd_tal_save

      REAL(DP), ALLOCATABLE, SAVE, PUBLIC ::
     &  PUXE(:), PUYE(:), PUXN(:), PUYN(:),
     &  PVXE(:), PVYE(:), PVXN(:), PVYN(:)

      CONTAINS

      !
      ! DETERMINING IF CELL (IX,IY) BELONGS TO B2 GRID
      !
      FUNCTION B2_CELL(IX,IY)

      IMPLICIT NONE

      LOGICAL :: B2_CELL
      INTEGER, INTENT(IN) :: IX, IY

      IF(IX.GE.0.AND.IX.LE.NDXP.AND.IY.GE.0.AND.IY.LE.NDYP) THEN
        B2_CELL=.TRUE.
      ELSE
        B2_CELL=.FALSE.
      END IF

      RETURN
      END FUNCTION B2_CELL

      subroutine eirene_wneutrals_alloc_arrays(ncv)
      implicit none
      logical :: exist
      integer, intent(in) :: ncv
      EXTERNAL :: EIRENE_EXIT_OWN

      if(allocated(dab2)) return
      allocate(dab2(ncv,natm,1))
      allocate(dmb2(ncv,nmol,1))
      allocate(dib2(ncv,nion,1))
      allocate(tab2(ncv,natm,1))
      allocate(tmb2(ncv,nmol,1))
      allocate(tib2(ncv,nion,1))
      allocate(rfluxa(ncv,natm,1))
      allocate(rfluxm(ncv,nmol,1))
      allocate(refluxa(ncv,natm,1))
      allocate(refluxm(ncv,nmol,1))
      allocate(pfluxa(ncv,natm,1))
      allocate(pfluxm(ncv,nmol,1))
      allocate(pefluxa(ncv,natm,1))
      allocate(pefluxm(ncv,nmol,1))
      allocate(tfluxa(ncv,natm,1))     ! IYS 20.01.2017
      allocate(tfluxm(ncv,nmol,1))     ! IYS 20.01.2017
      allocate(tefluxa(ncv,natm,1))    ! IYS 20.01.2017
      allocate(tefluxm(ncv,nmol,1))    ! IYS 20.01.2017
      allocate(emiss(ncv,1,1))
      allocate(emissmol(ncv,1,1))
      allocate(srcml(ncv,nmol,1))
      allocate(edissml(ncv,nmol,0:nstra+1))
      allocate(wldnek(nlim+nsts,0:nstra+1))
      allocate(wldnep(nlim+nsts,0:nstra+1))
      allocate(wldna(nlimps,natm,0:nstra+1))
      allocate(ewlda(nlimps,natm,0:nstra+1))
      allocate(wldnm(nlimps,nmol,0:nstra+1))
      allocate(ewldm(nlimps,nmol,0:nstra+1))
      allocate(wldra(nlimps,natm,0:nstra+1))
      allocate(wldrm(nlimps,nmol,0:nstra+1))
      allocate(wldpp(nlimps,npls,0:nstra+1))
      allocate(wldpa(nlimps,natm,0:nstra+1))
      allocate(wldpm(nlimps,nmol,0:nstra+1))
      allocate(wldpeb(nlimps,0:nstra+1))
      allocate(wldspt(nlimps,0:nstra+1))
      allocate(wldspta(nlimps,natm,0:nstra+1))
      allocate(wldsptm(nlimps,nmol,0:nstra+1))
      allocate(eneutrad(ncv,natm,0:nstra+1))
      allocate(emolrad(ncv,nmol,0:nstra+1))
      allocate(eionrad(ncv,nion,0:nstra+1))
      allocate(isrftype(nlimps))
      allocate(wlarea(nlim+nsts))
      allocate(wlabsrp(nspz,nlim+nsts))
      allocate(wlpump(nspz,nlim+nsts))
      allocate(eirpump(nspz))                                           !iyv 07.03.18
      allocate(eirspta(natmi), eirsptm(nmoli))                          !iyv 07.03.18
      allocate(eirtxt(nspz))
      allocate(eirdiag_nds_ind(1:nsts+1))
      allocate(pdena_int(natm,0:nstra),pdena_int_b2(natm,0:nstra),
     &         pdenm_int(nmol,0:nstra),pdenm_int_b2(nmol,0:nstra),
     &         pdeni_int(nion,0:nstra),pdeni_int_b2(nion,0:nstra),
     &         edena_int(natm,0:nstra),edena_int_b2(natm,0:nstra),
     &         edenm_int(nmol,0:nstra),edenm_int_b2(nmol,0:nstra),
     &         edeni_int(nion,0:nstra),edeni_int_b2(nion,0:nstra))

      allocate(volcel(ncv))

      allocate(wldnek_aver(nlim+nsts,0:nstra+1))
      allocate(wldnep_aver(nlim+nsts,0:nstra+1))
      allocate(wldna_aver(nlimps,natm,0:nstra+1))
      allocate(ewlda_aver(nlimps,natm,0:nstra+1))
      allocate(wldnm_aver(nlimps,nmol,0:nstra+1))
      allocate(ewldm_aver(nlimps,nmol,0:nstra+1))
      allocate(wldra_aver(nlimps,natm,0:nstra+1))
      allocate(wldrm_aver(nlimps,nmol,0:nstra+1))
      allocate(wldpp_aver(nlimps,npls,0:nstra+1))
      allocate(wldpa_aver(nlimps,natm,0:nstra+1))
      allocate(wldpm_aver(nlimps,nmol,0:nstra+1))
      allocate(wldpeb_aver(nlimps,0:nstra+1))
      allocate(wldspt_aver(nlimps,0:nstra+1))
      allocate(wldspta_aver(nlimps,natm,0:nstra+1))
      allocate(wldsptm_aver(nlimps,nmol,0:nstra+1))
      allocate(wlpump_aver(nspz,nlim+nsts))

C     READ IN THE NUMBER OF TRIANGLES
      INQUIRE (FILE=FORT_LC//'34',EXIST=EXIST)
      IF (.NOT.EXIST) GOTO 984
      OPEN (UNIT=34,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=984)
      REWIND 34
      READ(34,*,ERR=984) NTRII
      CLOSE(34)

      allocate(PDENA_aver(natmi,ntrii))
      allocate(PDENM_aver(nmoli,ntrii))
      allocate(EDENA_aver(natmi,ntrii))
      allocate(EDENM_aver(nmoli,ntrii))
      allocate(VXDENA_aver(natmi,ntrii))
      allocate(VXDENM_aver(nmoli,ntrii))
      allocate(VYDENA_aver(natmi,ntrii))
      allocate(VYDENM_aver(nmoli,ntrii))
      allocate(VZDENA_aver(natmi,ntrii))
      allocate(VZDENM_aver(nmoli,ntrii))

      !tamas zero init
      wldnep = 0
      wldna  = 0
      ewlda  = 0
      wldnm  = 0
      ewldm  = 0
      wldra  = 0
      wldrm  = 0
      wldpp  = 0
      wldpa  = 0
      wldpm  = 0
      wlarea = 0
      wldspta = 0
      wldsptm = 0

      aver_frac = 0
      wldnek_aver = 0
      wldnep_aver = 0
      wldna_aver  = 0
      ewlda_aver  = 0
      wldnm_aver  = 0
      ewldm_aver  = 0
      wldra_aver  = 0
      wldrm_aver  = 0
      wldpp_aver  = 0
      wldpa_aver  = 0
      wldpm_aver = 0
      wldpeb_aver = 0
      wldspt_aver = 0
      wldspta_aver = 0
      wldsptm_aver = 0
      wlpump_aver = 0

      aver_frac46 = 0
      PDENA_aver = 0
      PDENM_aver = 0
      EDENA_aver = 0
      EDENM_aver = 0
      VXDENA_aver = 0
      VXDENM_aver = 0
      VYDENA_aver = 0
      VYDENM_aver = 0
      VZDENA_aver = 0
      VZDENM_aver = 0

      return
  984 WRITE(*,*)
     w     "ERROR IN WNEUTRALS: CANNOT READ "//FORT//"34 ",
     .     "(TABLE OF CELLS OF TRIANGULAR GRID)"
      CALL EIRENE_EXIT_OWN(1)
      end subroutine eirene_wneutrals_alloc_arrays

      subroutine eirene_wneutrals_init(broadcast)
      use eirmod_mpi
      use eirmod_cpes
     , ,only: new_leader
      implicit none

      logical :: do_broadcast
      logical, optional :: broadcast

!     !c================================================================
      !c---------------------------------------------------------------<
      !c      write(iunout,*) '%%% wneutrals_init'
      !c--------------------------------------------------------------->

      IF (IFIRST_wneutral.EQ.0) then
      ! if it is not the first call, then we skip over initialization
        ifirst_wneutral=1
        call eirene_wneutrals_alloc_arrays(ncvp)

        call eirene_wneutrals_clean(.true.)

        eirdiag_nds_ind=0

        hlp_cnv=1./elcha
!pb  MOD_ADDV is a flag not a counter or index
!pb     ia0=mod_addv
        ia0=0
        ia1=ia0+natmi+nmoli
        ia2=ia1+natmi+nmoli
        ia3=ia2+natmi+nmoli

      endif ! IFIRST_wneutral
      !c---------------------------------------------------------------<
      if (new_leader) then
! there are new PEs among the leaders, we have to distribute the data
! to them
        if (present(broadcast)) then
          do_broadcast = broadcast
        else
          do_broadcast = .false.
        endif
        if (do_broadcast) call eirene_wneutrals_init_broadcast
      endif
      return
      end subroutine eirene_wneutrals_init

      subroutine eirene_wneutrals_init_broadcast
      use eirmod_mpi
      use eirmod_cpes
     , ,only: get_leader_comm
      implicit none
      integer :: comm
      comm = get_leader_comm()
      return
      end subroutine eirene_wneutrals_init_broadcast

      !c
      !c================================================================
      subroutine eirene_wneutrals_fill(istra_in)
      implicit none
      integer, intent(in) :: istra_in
      integer :: i,j,ic,in,istra_save,icell,jatm,jmol,jion
      real(DP) :: hlp
      character*36 hlp_frm
      !c
      !c write (iunout,*),'%%% wneutrals_fill: istra,istra_in = ',
      !c                                       istra,istra_in
      istra_save=istra
      istra=istra_in
      !c--------------------------------------------------------------->

!csw
!csw 21feb2012 corrected radiation from neutrals (atoms only),
!csw taken from SOLPS4.3 (V.Kotov)
!csw 04mar2013 shifted from wneutrals_save to here (wneutrals_fill)
!csw
      write(hlp_frm,'(a,i3,a)') '(a,i6,1p,',natmi+1,'(1x,e13.6))'
      eneutrad(:,:,istra) = 0.0_dp
      emolrad(:,:,istra) = 0.0_dp
      eionrad(:,:,istra) = 0.0_dp
      do icell=1,nCip
        do jatm=1,natmi
          if (lrael) eneutrad(icell,jatm,istra)=
     &      eneutrad(icell,jatm,istra)+rael(jatm,icell)*voltal(icell)
        enddo
        do jmol=1,nmoli
          if (lrmel) emolrad(icell,jmol,istra)=
     &      emolrad(icell,jmol,istra)+rmel(jmol,icell)*voltal(icell)
        enddo
        do jion=1,nioni
          if (lriel) eionrad(icell,jion,istra)=
     &      eionrad(icell,jion,istra)+riel(jion,icell)*voltal(icell)
        enddo
      enddo

      volcel = 0.0_dp
      do in=1,ntrii
        ic=ixtri(in)
        if(ic.gt.0) volcel(ic) = volcel(ic) + vol(in)
      end do

      !c
      !C map 1d EIRENE neutral densities and temperatures on 2d arrays
      !c and change the units to SI for plotting in B2
      !c
      !write(iunout,*) 'istra ',istra
      !write(iunout,*) 'natmi, nmoli, nioni ',natmi,nmoli,nioni
      !write(iunout,*) 'ndxp, ndyp ',ndxp,ndyp
      !crfs     IF (WTOTP(0,ISTRA).EQ.0.) GOTO 60
      do jatm=1,natmi
        do in=1,ncip
            if (lpdena) dab2(in,jatm,1)=dab2(in,jatm,1)+
     &                         pdena(jatm,in)*1.0e6_dp
!   IYS 20.01.2017
! rfluxa, pfluxa and tfluxa are cell-centered particle flux densities
! The Eirene momentum densities are converted from g*(cm/s)/cm^3 to kg*(m/s)/m^3
            if (LVXDENA.AND.LVYDENA) then
              rfluxa(in,jatm,1)=rfluxa(in,jatm,1)+
     &             (VXDENA(jatm,in)*PVX(in)+
     &              VYDENA(jatm,in)*PVY(in))*
     &          10.0_dp/RMASSA(jatm)/AMUAKG
              pfluxa(in,jatm,1)=pfluxa(in,jatm,1)+
     &             (VXDENA(jatm,in)*PUX(in)+
     &              VYDENA(jatm,in)*PUY(in))*
     &          10.0_dp/RMASSA(jatm)/AMUAKG
            else if (laddv) then
             if(ia0+jatm.le.nadv)
     >        rfluxa(in,jatm,1)=rfluxa(in,jatm,1)+
     &                                   addv(ia0+jatm,in)*1.0e4_dp
             if(ia2+jatm.le.nadv)
     >        pfluxa(in,jatm,1)=pfluxa(in,jatm,1)+
     &                                   addv(ia2+jatm,in)*1.0e4_dp
            end if
            if (LVZDENA) tfluxa(in,jatm,1)=tfluxa(in,jatm,1)+
     &             VZDENA(jatm,in)*10.0_dp/RMASSA(jatm)/AMUAKG*
     &             sign(1.0_dp,PUX(in)*PVY(in)-PUY(in)*PVX(in))

!   IYS 20.01.2017
!   refluxa and pefluxa - energy flux density components -
!   cannot be expressed through Eirene standard tallies,
!   and their expression through additional tallies should be double-checked
!   so do not use these arrays in until such a check would be performed

            if (laddv) then
              if(ia1+jatm.le.nadv)
     >         refluxa(in,jatm,1)=refluxa(in,jatm,1)+
     &                            addv(ia1+jatm,in)*1.0e4_dp*elcha
              if(ia3+jatm.le.nadv)
     >         pefluxa(in,jatm,1)=pefluxa(in,jatm,1)+
     &                            addv(ia3+jatm,in)*1.0e4_dp*elcha
            end if
            if (ledena) tab2(in,jatm,1)=tab2(in,jatm,1)+
     &                            edena(jatm,in)
        end do
      end do
      do jmol=1,nmoli
        edissml(:,jmol,istra) = 0.0_dp
        do in=1,ncip
            if (lpdenm) dmb2(in,jmol,1)=dmb2(in,jmol,1)+
     &                         pdenm(jmol,in)*1.0e6_dp
!   IYS 20.01.2017
! rfluxm, pfluxm and tfluxm are cell-centered particle flux densities
! The Eirene momentum densities are converted from g*(cm/s)/cm^3 to kg*(m/s)/m^3
            if (LVXDENM.AND.LVYDENM) then
              rfluxm(in,jmol,1)=rfluxm(in,jmol,1)+
     &             (VXDENM(jmol,in)*PVX(in)+
     &              VYDENM(jmol,in)*PVY(in))*
     &           10.0_dp/RMASSM(jmol)/AMUAKG

              pfluxm(in,jmol,1)=pfluxm(in,jmol,1)+
     &             (VXDENM(jmol,in)*PUX(in)+
     &              VYDENM(jmol,in)*PUY(in))*
     &           10.0_dp/RMASSM(jmol)/AMUAKG
            else if (laddv) then
             if(ia0+natmi+jmol.le.nadv)
     >         rfluxm(in,jmol,1)=rfluxm(in,jmol,1)+
     &                                 addv(ia0+natmi+jmol,in)*1.0e4_dp
             if(ia2+natmi+jmol.le.nadv)
     >         pfluxm(in,jmol,1)=pfluxm(in,jmol,1)+
     &                                 addv(ia2+natmi+jmol,in)*1.0e4_dp
            end if
            if (LVZDENM) tfluxm(in,jmol,1)=tfluxm(in,jmol,1)+
     &             VZDENM(jmol,in)*10.0_dp/RMASSM(jmol)/AMUAKG*
     &             sign(1.0_dp,PUX(in)*PVY(in)-PUY(in)*PVX(in))

!   IYS 20.01.2017
!   refluxm and pefluxm - energy flux density components -
!   cannot be expressed through Eirene standard tallies,
!   and their expression through additional tallies should be double-checked
!   so do not use these arrays in until such a check would be performed

            if(ia1+natmi+jmol.le.nadv .and. laddv)
     >       refluxm(in,jmol,1)=refluxm(in,jmol,1)+
     &                      addv(ia1+natmi+jmol,in)*1.0e4_dp*elcha
            if(ia3+natmi+jmol.le.nadv .and. laddv)
     >       pefluxm(in,jmol,1)=pefluxm(in,jmol,1)+
     &                      addv(ia3+natmi+jmol,in)*1.0e4_dp*elcha
            if (ledenm) tmb2(in,jmol,1)=tmb2(in,jmol,1)+
     &                      edenm(jmol,in)
            if (lpmml) srcml(in,jmol,1)=srcml(in,jmol,1)+
     &                      pmml(jmol,in)*vol(in)*hlp_cnv
            !c*** Potential energy source related to hydrogen molecules
            !c*** for determining the radiation (W per cell)
            if (ncharm(jmol).eq.2 .and. lpmml)
     >       edissml(in,jmol,istra)=edissml(in,jmol,istra)+
     &                      pmml(jmol,in)*vol(in)*diss_pot_H2
        end do
      end do
      do jion=1,nioni
        do in=1,ncip
            if (lpdeni) dib2(in,jion,1)=dib2(in,jion,1)+
     &                         pdeni(jion,in)*1.0e6_dp
            if (ledeni) tib2(in,jion,1)=tib2(in,jion,1)+
     &                         edeni(jion,in)
        end do
      end do

!c
!c*** Rescale the surface data from A to 1/sec and average the energy
!c
      do i=1,nlimps
        wldnek(i,istra)=0.0_dp
        wldnep(i,istra)=0.0_dp
        wldpeb(i,istra)=0.0_dp
!c*** hlp accumulates the power taken away with re-emitted particles
        hlp=0.
        do j=1,natmi
          if (leotat)  wldnek(i,istra)=wldnek(i,istra)+eotat(j,i)
          if (lerfpat) wldpeb(i,istra)=wldpeb(i,istra)+erfpat(j,i)
          if (lerfaat) hlp=hlp+erfaat(j,i)
          if (lerfmat) hlp=hlp+erfmat(j,i)
          if (lerfiat) hlp=hlp+erfiat(j,i)
          if (lerfphat) hlp=hlp+erfphat(j,i)
          if (lpotat.and.leotat) then
            if(potat(j,i).gt.0.) then
              ewlda(i,j,istra)=eotat(j,i)/potat(j,i)
            else
              ewlda(i,j,istra)=0.0_dp
            end if
          else
            ewlda(i,j,istra)=0.0_dp
          end if
          if (lpotat) wldna(i,j,istra)=hlp_cnv*potat(j,i)
          wldra(i,j,istra)=0.
          if (lprfaat) wldra(i,j,istra)=wldra(i,j,istra)+prfaat(j,i)
          if (lprfmat) wldra(i,j,istra)=wldra(i,j,istra)+prfmat(j,i)
          if (lprfiat) wldra(i,j,istra)=wldra(i,j,istra)+prfiat(j,i)
          wldra(i,j,istra)=hlp_cnv*wldra(i,j,istra)
          if (lprfpat) wldpa(i,j,istra)=hlp_cnv*prfpat(j,i)
        end do
        do j=1,nmoli
          if (leotml)  wldnek(i,istra)=wldnek(i,istra)+eotml(j,i)
          if (lerfpml) wldpeb(i,istra)=wldpeb(i,istra)+erfpml(j,i)
          if (lerfaml) hlp=hlp+erfaml(j,i)
          if (lerfmml) hlp=hlp+erfmml(j,i)
          if (lerfiml) hlp=hlp+erfiml(j,i)
          if (lerfphml) hlp=hlp+erfphml(j,i)
          if (lpotml.and.leotml) then
            if(potml(j,i).gt.0.0_dp) then
              ewldm(i,j,istra)=eotml(j,i)/potml(j,i)
            else
              ewldm(i,j,istra)=0.0_dp
            end if
          else
            ewldm(i,j,istra)=0.0_dp
          end if
          if (lpotml) wldnm(i,j,istra)=hlp_cnv*potml(j,i)
          wldrm(i,j,istra)=0.0_dp
          if (lprfaml) wldrm(i,j,istra)=wldrm(i,j,istra)+prfaml(j,i)
          if (lprfmml) wldrm(i,j,istra)=wldrm(i,j,istra)+prfmml(j,i)
          if (lprfiml) wldrm(i,j,istra)=wldrm(i,j,istra)+prfiml(j,i)
          wldrm(i,j,istra)=hlp_cnv*wldrm(i,j,istra)
          if (lprfpml) wldpm(i,j,istra)=hlp_cnv*prfpml(j,i)
        end do
        do j=1,nioni
          if (leotio)  wldnek(i,istra)=wldnek(i,istra)+eotio(j,i)
          if (lerfpio) wldpeb(i,istra)=wldpeb(i,istra)+erfpio(j,i)
          if (lerfaio) hlp=hlp+erfaio(j,i)
          if (lerfmio) hlp=hlp+erfmio(j,i)
          if (lerfiio) hlp=hlp+erfiio(j,i)
          if (lerfphio) hlp=hlp+erfphio(j,i)
        end do
        do j=1,nfla
          if (lpotpl) wldpp(i,j,istra)=hlp_cnv*potpl(j,i)
        end do
        !c*** recombination energy of hydrogen molecules
        do j=1,nmoli
          if(ncharm(j).eq.2) then
            if(lprfaml) wldnep(i,istra)=
     .                  wldnep(i,istra)+diss_pot_H2*prfaml(j,i)
            if(lprfpml) wldnep(i,istra)=
     .                  wldnep(i,istra)+diss_pot_H2*prfpml(j,i)
          end if
        end do
        !c*** subtract the outcoming power
        wldnek(i,istra)=wldnek(i,istra)-hlp

        !c*** calculate the flux of sputtered particles
        if (lspttot) wldspt(i,istra)=hlp_cnv*spttot(i)
        do j=1,natmi
          wldspta(i,j,istra)=0.
          if (lsptaat) wldspta(i,j,istra)=wldspta(i,j,istra)+sptaat(j,i)
          if (lsptmat) wldspta(i,j,istra)=wldspta(i,j,istra)+sptmat(j,i)
          if (lsptiat) wldspta(i,j,istra)=wldspta(i,j,istra)+sptiat(j,i)
          if (lsptpat) wldspta(i,j,istra)=wldspta(i,j,istra)+sptpat(j,i)
          if (lsptphat) wldspta(i,j,istra)=
     .                  wldspta(i,j,istra)+sptphat(j,i)
          wldspta(i,j,istra)=wldspta(i,j,istra)*hlp_cnv
        end do
        do j=1,nmoli
          wldsptm(i,j,istra)=0.
          if (lsptaml) wldsptm(i,j,istra)=wldsptm(i,j,istra)+sptaml(j,i)
          if (lsptmml) wldsptm(i,j,istra)=wldsptm(i,j,istra)+sptmml(j,i)
          if (lsptiml) wldsptm(i,j,istra)=wldsptm(i,j,istra)+sptiml(j,i)
          if (lsptpml) wldsptm(i,j,istra)=wldsptm(i,j,istra)+sptpml(j,i)
          if (lsptphml) wldsptm(i,j,istra)=
     .                  wldsptm(i,j,istra)+sptphml(j,i)
          wldsptm(i,j,istra)=wldsptm(i,j,istra)*hlp_cnv
        end do
      end do
      !c
      call update_integrals(istra)

      !c---------------------------------------------------------------<
      istra=istra_save
      return
      end subroutine eirene_wneutrals_fill

      !c
      subroutine eirene_wneutrals_save
      implicit none
      integer :: i,j,k,l,ic,jatm,jmol,jion
      real(DP), save :: pumpsum, sptsum
      character*36 hlp_frm
      external :: eirene_leer
      !c
      !c      write (iunout,*) '%%% wneutrals_save: istra = ',istra
      !c*** Calculate the totals (stratum 0)
      !c

      ! Zero out arrays for inactive strata
      do k=1,nstrai
         if (.not. nlsron(k)) then
            wldnek(:,k) = 0.0_DP
            wldnep(:,k) = 0.0_DP
            wldpeb(:,k) = 0.0_DP
            wldspt(:,k) = 0.0_DP
            ewlda(:,:,k) = 0.0_DP
            wldna(:,:,k) = 0.0_DP
            wldra(:,:,k) = 0.0_DP
            wldpa(:,:,k) = 0.0_DP
            wldspta(:,:,k) = 0.0_DP
            ewldm(:,:,k) = 0.0_DP
            wldnm(:,:,k) = 0.0_DP
            wldrm(:,:,k) = 0.0_DP
            wldpm(:,:,k) = 0.0_DP
            wldsptm(:,:,k) = 0.0_DP
            wldpp(:,:,k) = 0.0_DP
            eneutrad(:,:,k) = 0.0_DP
            emolrad(:,:,k) = 0.0_DP
            eionrad(:,:,k) = 0.0_DP
            edissml(:,:,k) = 0.0_DP
            wldnek_aver(:,k) = 0.0_DP
            wldnep_aver(:,k) = 0.0_DP
            wldna_aver(:,:,k) = 0.0_DP
            ewlda_aver(:,:,k) = 0.0_DP
            wldnm_aver(:,:,k) = 0.0_DP
            ewldm_aver(:,:,k) = 0.0_DP
            wldra_aver(:,:,k) = 0.0_DP
            wldrm_aver(:,:,k) = 0.0_DP
            wldpp_aver(:,:,k) = 0.0_DP
            wldpa_aver(:,:,k) = 0.0_DP
            wldpm_aver(:,:,k) = 0.0_DP
            wldpeb_aver(:,k) = 0.0_DP
            wldspt_aver(:,k) = 0.0_DP
            wldspta_aver(:,:,k) = 0.0_DP
            wldsptm_aver(:,:,k) = 0.0_DP
         end if
      end do

      sptsum=0.0_DP
      do i=1,nlimps
        wldnek(i,0)=0.0_DP
        wldnep(i,0)=0.0_DP
        wldpeb(i,0)=0.0_DP
        wldspt(i,0)=0.0_DP
        do j=1,natmi
          ewlda(i,j,0)=0.0_DP
          wldna(i,j,0)=0.0_DP
          wldra(i,j,0)=0.0_DP
          wldpa(i,j,0)=0.0_DP
          wldspta(i,j,0)=0.0_DP
        end do
        do j=1,nmoli
          ewldm(i,j,0)=0.0_DP
          wldnm(i,j,0)=0.0_DP
          wldrm(i,j,0)=0.0_DP
          wldpm(i,j,0)=0.0_DP
          wldsptm(i,j,0)=0.0_DP
        end do
        do j=1,nfla
          wldpp(i,j,0)=0.0_DP
        end do
        !c
        do k=1,nstrai
          wldnek(i,0)=wldnek(i,0)+wldnek(i,k)
          wldnep(i,0)=wldnep(i,0)+wldnep(i,k)
          wldpeb(i,0)=wldpeb(i,0)+wldpeb(i,k)
          wldspt(i,0)=wldspt(i,0)+wldspt(i,k)
          do j=1,natmi
            wldna(i,j,0)=wldna(i,j,0)+wldna(i,j,k)
            wldra(i,j,0)=wldra(i,j,0)+wldra(i,j,k)
            wldpa(i,j,0)=wldpa(i,j,0)+wldpa(i,j,k)
            ewlda(i,j,0)=ewlda(i,j,0)+ewlda(i,j,k)*wldna(i,j,k)
            wldspta(i,j,0)=wldspta(i,j,0)+wldspta(i,j,k)
          end do
          do j=1,nmoli
            wldnm(i,j,0)=wldnm(i,j,0)+wldnm(i,j,k)
            wldrm(i,j,0)=wldrm(i,j,0)+wldrm(i,j,k)
            wldpm(i,j,0)=wldpm(i,j,0)+wldpm(i,j,k)
            ewldm(i,j,0)=ewldm(i,j,0)+ewldm(i,j,k)*wldnm(i,j,k)
            wldsptm(i,j,0)=wldsptm(i,j,0)+wldsptm(i,j,k)
          end do
          do j=1,nfla
            wldpp(i,j,0)=wldpp(i,j,0)+wldpp(i,j,k)
          end do
        end do
        do j=1,natmi
          if(wldna(i,j,0).gt.0.0_DP) then
            ewlda(i,j,0)=ewlda(i,j,0)/wldna(i,j,0)
          else
            ewlda(i,j,0)=0.0_DP
          end if
        end do
        do j=1,nmoli
          if(wldnm(i,j,0).gt.0.0_DP) then
            ewldm(i,j,0)=ewldm(i,j,0)/wldnm(i,j,0)
          else
            ewldm(i,j,0)=0.0_DP
          end if
        end do
        sptsum=sptsum+wldspt(i,0)
      end do
!pb
      eneutrad(:,:,0) = 0.0_DP
      emolrad(:,:,0) = 0.0_DP
      eionrad(:,:,0) = 0.0_DP
      edissml(:,:,0) = 0.0_DP
      do ic = 1, ncip
        do jatm=1,natmi
          eneutrad(ic,jatm,0) = sum(eneutrad(ic,jatm,1:nstrai))
        end do
        do jmol=1,nmoli
          edissml(ic,jmol,0) = sum(edissml(ic,jmol,1:nstrai))
          emolrad(ic,jmol,0) = sum(emolrad(ic,jmol,1:nstrai))
        end do
        do jion=1,nioni
          eionrad(ic,jion,0) = sum(eionrad(ic,jion,1:nstrai))
        end do
      end do
      !c--------------------------------------------------------------->


      !c*** Surface type and properties
      wlabsrp = 0.0_DP
      pumpsum = 0.0_DP
      eirpump = 0.0_DP
      eirspta = 0.0_DP
      eirsptm = 0.0_DP
      do i=1,nlim+nsts
        isrftype(i)=iliin(i)
        wlarea(i)=1.0e-4_DP*sarea(i)
        if(iliin(i).eq.1) then !{
          do j=1,nspz !{
            wlabsrp(j,i)=1.0_DP-recyct(j,i)
          end do !}
        else if(iliin(i).eq.2) then !}{
          wlabsrp(1:nspz,i)=1.0_DP
        end if !}
      !c** calculate the pumped flux for each nspz particle
        do j=1,nspz
          IF (LSPUMP) THEN
            wlpump(j,i)=SPUMP(j,i)*hlp_cnv
            pumpsum=pumpsum+wlpump(j,i)
            eirpump(j) = eirpump(j) + wlpump(j,i)
          END IF
        end do
        do j = 1, natmi
          eirspta(j) = eirspta(j) + wldspta(i,j,0)
        enddo
        do j = 1, nmoli
          eirsptm(j) = eirsptm(j) + wldsptm(i,j,0)
        enddo
      end do

      do j=1,nspz !{
        eirtxt(j)=texts(j)
      end do !}

      if (aver_frac.ne.0.0_DP) then ! average fluxes, som 02.04.2019
        wldnek_aver(:,0)=(1._DP-aver_frac)*wldnek(:,0)+
     .   aver_frac*wldnek_aver(:,0)
        wldnep_aver(:,0)=(1._DP-aver_frac)*wldnep(:,0)+
     .   aver_frac*wldnep_aver(:,0)
        wldpeb_aver(:,0)=(1._DP-aver_frac)*wldpeb(:,0)+
     .   aver_frac*wldpeb_aver(:,0)
        wldspt_aver(:,0)=(1._DP-aver_frac)*wldspt(:,0)+
     .   aver_frac*wldspt_aver(:,0)
        wldna_aver(:,:,0)=(1._DP-aver_frac)*wldna(:,:,0)+
     .   aver_frac*wldna_aver(:,:,0)
        wldra_aver(:,:,0)=(1._DP-aver_frac)*wldra(:,:,0)+
     .   aver_frac*wldra_aver(:,:,0)
        wldpa_aver(:,:,0)=(1._DP-aver_frac)*wldpa(:,:,0)+
     .   aver_frac*wldpa_aver(:,:,0)
        ewlda_aver(:,:,0)=(1._DP-aver_frac)*ewlda(:,:,0)+
     .   aver_frac*ewlda_aver(:,:,0)
        wldspta_aver(:,:,0)=(1._DP-aver_frac)*wldspta(:,:,0)+
     .   aver_frac*wldspta_aver(:,:,0)
        wldnm_aver(:,:,0)=(1._DP-aver_frac)*wldnm(:,:,0)+
     .   aver_frac*wldnm_aver(:,:,0)
        wldrm_aver(:,:,0)=(1._DP-aver_frac)*wldrm(:,:,0)+
     .   aver_frac*wldrm_aver(:,:,0)
        wldpm_aver(:,:,0)=(1._DP-aver_frac)*wldpm(:,:,0)+
     .   aver_frac*wldpm_aver(:,:,0)
        ewldm_aver(:,:,0)=(1._DP-aver_frac)*ewldm(:,:,0)+
     .   aver_frac*ewldm_aver(:,:,0)
        wldsptm_aver(:,:,0)=(1._DP-aver_frac)*wldsptm(:,:,0)+
     .   aver_frac*wldsptm_aver(:,:,0)
        wldpp_aver(:,:,0)=(1._DP-aver_frac)*wldpp(:,:,0)+
     .   aver_frac*wldpp_aver(:,:,0)
        wlpump_aver(:,:)=(1._DP-aver_frac)*wlpump(:,:)+
     .   aver_frac*wlpump_aver(:,:)
      endif
      if (aver_frac46.ne.0.0_DP) then ! average fluxes, som 02.04.2019
        PDENA_aver(:,1:ntrii)=(1._DP-aver_frac46)*PDENA(:,1:ntrii)+
     .   aver_frac46*PDENA_aver(:,1:ntrii)
        PDENM_aver(:,1:ntrii)=(1._DP-aver_frac46)*PDENM(:,1:ntrii)+
     .   aver_frac46*PDENM_aver(:,1:ntrii)
        EDENA_aver(:,1:ntrii)=(1._DP-aver_frac46)*EDENA(:,1:ntrii)+
     .   aver_frac46*EDENA_aver(:,1:ntrii)
        EDENM_aver(:,1:ntrii)=(1._DP-aver_frac46)*EDENM(:,1:ntrii)+
     .   aver_frac46*EDENM_aver(:,1:ntrii)
        VXDENA_aver(:,1:ntrii)=(1._DP-aver_frac46)*VXDENA(:,1:ntrii)+
     .   aver_frac46*VXDENA_aver(:,1:ntrii)
        VXDENM_aver(:,1:ntrii)=(1._DP-aver_frac46)*VXDENM(:,1:ntrii)+
     .   aver_frac46*VXDENM_aver(:,1:ntrii)
        VYDENA_aver(:,1:ntrii)=(1._DP-aver_frac46)*VYDENA(:,1:ntrii)+
     .   aver_frac46*VYDENA_aver(:,1:ntrii)
        VYDENM_aver(:,1:ntrii)=(1._DP-aver_frac46)*VYDENM(:,1:ntrii)+
     .   aver_frac46*VYDENM_aver(:,1:ntrii)
        VZDENA_aver(:,1:ntrii)=(1._DP-aver_frac46)*VZDENA(:,1:ntrii)+
     .   aver_frac46*VZDENA_aver(:,1:ntrii)
        VZDENM_aver(:,1:ntrii)=(1._DP-aver_frac46)*VZDENM(:,1:ntrii)+
     .   aver_frac46*VZDENM_aver(:,1:ntrii)
      endif

      !c
      !c*** Calculate the temperatures
      !c
      do ic=1,ncip
        do jatm=1,natmi
          if (dab2(ic,jatm,1).gt.0.0_DP) then
            tab2(ic,jatm,1)=tab2(ic,jatm,1)/dab2(ic,jatm,1)*
     &                            elcha*2.0_DP/3.0_DP*1.0e6_DP
          else
            tab2(ic,jatm,1)=1.0e-6_DP*elcha
          end if
        end do
        do jmol=1,nmoli
          if (dmb2(ic,jmol,1).gt.0.0_DP) then
            tmb2(ic,jmol,1)=tmb2(ic,jmol,1)/dmb2(ic,jmol,1)*
     &                            elcha*2.0_DP/3.0_DP*1.0e6_DP
          else
            tmb2(ic,jmol,1)=1.0e-6_DP*elcha
          end if
        end do
        do jion=1,nioni
          if (dib2(ic,jion,1).gt.0.0_DP) then
            tib2(ic,jion,1)=tib2(ic,jion,1)/dib2(ic,jion,1)*
     &                            elcha*2.0_DP/3.0_DP*1.0e6_DP
          else
            tib2(ic,jion,1)=1.0e-6_DP*elcha
          end if
        end do
      end do

      !c*** print some neutral fluxes across the "non-default" surfaces

      write(hlp_frm,'(a,i3,a)')
     & '(/16x,2(2x,a8,1x),',3*(natmi+nmoli)+nfla,'(3x,a6,i2))'
      write(iunout,hlp_frm) '  area','  power',
     &            ('atflx',i,i=1,natmi),('mlflx',i,i=1,nmoli),
     &            ('atflxr',i,i=1,natmi),('mlflxr',i,i=1,nmoli),
     &            ('atflxp',i,i=1,natmi),('mlflxp',i,i=1,nmoli),
     &            ('plsflx',i,i=1,nfla)
      do i=1,nstsi
        j=nlim+i
        do k=1,nfla
          wldpp(j,k,0)=0.0_DP
          do l=1,nstrai
            wldpp(j,k,0)=wldpp(j,k,0)+wldpp(j,k,l)
          end do
        end do
        do k=1,natmi
          wldpa(j,k,0)=0.0_DP
          do l=1,nstrai
            wldpa(j,k,0)=wldpa(j,k,0)+wldpa(j,k,l)
          end do
        end do
        do k=1,nmoli
          wldpm(j,k,0)=0.0_DP
          do l=1,nstrai
            wldpm(j,k,0)=wldpm(j,k,0)+wldpm(j,k,l)
          end do
        end do

        write(hlp_frm,'(a,i3,a)')
     .   '(a,i4,1p,',2+3*(natmi+nmoli)+nfla,'e11.3)'
        if (nmoli.gt.0) then
          if(wldnek(j,0).ne.0.0_DP .or.
     &       wldna(j,1,0).ne.0.0_DP .or. wldnm(j,1,0).ne.0.0_DP .or.
     &       wldra(j,1,0).ne.0.0_DP .or. wldrm(j,1,0).ne.0.0_DP .or.
     &       wldpa(j,1,0).ne.0.0_DP .or. wldpm(j,1,0).ne.0.0_DP .or.
     &       wldpp(j,1,0).ne.0.0_DP)
     >        write(iunout,hlp_frm)
     &          'non-def-surf ',i, 1.0e-4_DP*sarea(j),
     &               1.0e-6_DP*wldnek(j,0),
     &              (wldna(j,k,0),k=1,natmi),(wldnm(j,k,0),k=1,nmoli),
     &              (wldra(j,k,0),k=1,natmi),(wldrm(j,k,0),k=1,nmoli),
     &              (wldpa(j,k,0),k=1,natmi),(wldpm(j,k,0),k=1,nmoli),
     &              (wldpp(j,k,0),k=1,nfla)
        else
          if(wldnek(j,0).ne.0.0_DP .or.
     &       wldna(j,1,0).ne.0.0_DP .or. wldra(j,1,0).ne.0.0_DP .or.
     &       wldpa(j,1,0).ne.0.0_DP .or. wldpp(j,1,0).ne.0.0_DP)
     >        write(iunout,hlp_frm)
     &          'non-def-surf ',i, 1.0e-4_DP*sarea(j),
     &               1.0e-6_DP*wldnek(j,0),
     &              (wldna(j,k,0),k=1,natmi),(wldra(j,k,0),k=1,natmi),
     &              (wldpa(j,k,0),k=1,natmi),(wldpp(j,k,0),k=1,nfla)
        end if
      end do
      call eirene_leer(1)

      CALL WNEUTRAL_FLUXES !VK
!carlis printout of pumped and sputtered flux for energy balance script
      WRITE(IUNOUT,*) "PUMPED FLUX:",pumpsum
      WRITE(IUNOUT,*) "SPUTTERED FLUX:",sptsum
      CALL EIRENE_LEER(1)
      CALL UPDATE_INTEGRALS(0) !VKI
      return
      end subroutine eirene_wneutrals_save

      subroutine eirene_wneutrals_clean(fullclean)
      implicit none
      logical, intent(in) :: fullclean
      if(allocated(dab2)) then
        dab2=0.0_DP
        dmb2=0.0_DP
        dib2=0.0_DP
        tab2=0.0_DP
        tmb2=0.0_DP
        tib2=0.0_DP
        rfluxa=0.0_DP
        rfluxm=0.0_DP
        refluxa=0.0_DP
        refluxm=0.0_DP
        pfluxa=0.0_DP
        pfluxm=0.0_DP
        pefluxa=0.0_DP
        pefluxm=0.0_DP
        tfluxa=0.0_DP
        tfluxm=0.0_DP
        tefluxa=0.0_DP
        tefluxm=0.0_DP
        emiss=0.0_DP
        emissmol=0.0_DP
        srcml=0.0_DP
        edissml=0.0_DP
        if (fullclean) then
          wldnek=0.0_DP
          wldnep=0.0_DP
          wldna=0.0_DP
          ewlda=0.0_DP
          wldnm=0.0_DP
          ewldm=0.0_DP
          wldra=0.0_DP
          wldrm=0.0_DP
          wldpp=0.0_DP
          wldpa=0.0_DP
          wldpm=0.0_DP
          wldpeb=0.0_DP
          wldspt=0.0_DP
          wldspta=0.0_DP
          wldsptm=0.0_DP
          eneutrad=0.0_DP
          emolrad=0.0_DP
          eionrad=0.0_DP
          pdena_int=0.0_DP
          pdena_int_b2=0.0_DP
          pdenm_int=0.0_DP
          pdenm_int_b2=0.0_DP
          pdeni_int=0.0_DP
          pdeni_int_b2=0.0_DP
          edena_int=0.0_DP
          edena_int_b2=0.0_DP
          edenm_int=0.0_DP
          edenm_int_b2=0.0_DP
          edeni_int=0.0_DP
          edeni_int_b2=0.0_DP
        endif
      endif
      return
      end subroutine eirene_wneutrals_clean

      subroutine eirene_wneutrals_dealloc
      implicit none

      if(allocated(dab2)) then
        deallocate(dab2)
        deallocate(dmb2)
        deallocate(dib2)
        deallocate(tab2)
        deallocate(tmb2)
        deallocate(tib2)
        deallocate(rfluxa)
        deallocate(rfluxm)
        deallocate(refluxa)
        deallocate(refluxm)
        deallocate(pfluxa)
        deallocate(pfluxm)
        deallocate(pefluxa)
        deallocate(pefluxm)
        deallocate(emiss)
        deallocate(emissmol)
        deallocate(srcml)
        deallocate(edissml)
        deallocate(wldnek)
        deallocate(wldnep)
        deallocate(wldna)
        deallocate(ewlda)
        deallocate(wldnm)
        deallocate(ewldm)
        deallocate(wldra)
        deallocate(wldrm)
        deallocate(wldpp)
        deallocate(wldpa)
        deallocate(wldpm)
        deallocate(wldpeb)
        deallocate(wldspt)
        deallocate(wldspta)
        deallocate(wldsptm)
        deallocate(isrftype)
        deallocate(eneutrad,emolrad,eionrad)
      endif

      if(allocated(wlarea)) then
        deallocate(wlarea,wlabsrp,eirtxt,wlpump)
        deallocate(eirdiag_nds_ind)
      endif

      if(allocated(wldnek_aver)) then
        deallocate(wldnek_aver)
        deallocate(wldnep_aver)
        deallocate(wldna_aver)
        deallocate(ewlda_aver)
        deallocate(wldnm_aver)
        deallocate(ewldm_aver)
        deallocate(wldra_aver)
        deallocate(wldrm_aver)
        deallocate(wldpp_aver)
        deallocate(wldpa_aver)
        deallocate(wldpm_aver)
        deallocate(wldpeb_aver)
        deallocate(wldspt_aver)
        deallocate(wldspta_aver)
        deallocate(wldsptm_aver)
        deallocate(wlpump_aver)
      endif

      if(allocated(eirpump)) then                                        !iyv 07.03.18 {
        deallocate(eirpump, eirspta, eirsptm)
      endif                                                              !iyv 07.03.18 }

      ifirst_wneutral=0

      if(allocated(wldna_res)) then
        deallocate(wldna_res, wldnm_res,
     &             ewldrp_res, ewldmr_res,
     &             ewlda_res, ewldm_res, ewldt_res,
     &             ewldea_res, ewldem_res,
     &             wldspt_res, wldspta_res,
     &             wlpump_res, sarea_res)
      end if

      if(allocated(pdena_int))
     > deallocate(pdena_int,pdena_int_b2,pdenm_int,pdenm_int_b2,
     &            pdeni_int,pdeni_int_b2,edena_int,edena_int_b2,
     &            edenm_int,edenm_int_b2,edeni_int,edeni_int_b2)

      if(allocated(pdena_aver))
     > deallocate(pdena_aver,pdenm_aver,
     &            edena_aver,edenm_aver,
     &            vxdena_aver,vxdenm_aver,
     &            vydena_aver,vydenm_aver,
     &            vzdena_aver,vzdenm_aver)

      return
      end subroutine eirene_wneutrals_dealloc

      subroutine eirene_wneutrals_send(istra, rank)
      use eirmod_mpi
      implicit none
      integer, intent(in) :: istra
      integer, intent(in) :: rank
      integer :: ierr
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
      external :: mpi_send
#endif

      call mpi_send(edissml(1,1,istra) , size(edissml(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(eneutrad(1,1,istra) , size(eneutrad(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(emolrad(1,1,istra) , size(emolrad(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(eionrad(1,1,istra) , size(eionrad(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(ewlda(1,1,istra) , size(ewlda(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(ewldm(1,1,istra) , size(ewldm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldna(1,1,istra) , size(wldna(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldnek(1,istra) , size(wldnek(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldnep(1,istra) , size(wldnep(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldnm(1,1,istra) , size(wldnm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldpa(1,1,istra) , size(wldpa(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldpeb(1,istra) , size(wldpeb(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldpm(1,1,istra) , size(wldpm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldpp(1,1,istra) , size(wldpp(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldra(1,1,istra) , size(wldra(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldrm(1,1,istra) , size(wldrm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldspt(1,istra) , size(wldspt(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldspta(1,1,istra) , size(wldspta(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(wldsptm(1,1,istra) , size(wldsptm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(edena_int(:,istra) , size(edena_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(edena_int_b2(:,istra) , size(edena_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(edeni_int(:,istra) , size(edeni_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(edeni_int_b2(:,istra) , size(edeni_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(edenm_int(:,istra) , size(edenm_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(edenm_int_b2(:,istra) , size(edenm_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(pdena_int(:,istra) , size(pdena_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(pdena_int_b2(:,istra) , size(pdena_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(pdeni_int(:,istra) , size(pdeni_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(pdeni_int_b2(:,istra) , size(pdeni_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(pdenm_int(:,istra) , size(pdenm_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      call mpi_send(pdenm_int_b2(:,istra) , size(pdenm_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, 0, istra, MPI_COMM_WORLD,
     &              ierr)

      return
      end subroutine eirene_wneutrals_send

      subroutine eirene_wneutrals_recv(istra, rank)
      use eirmod_mpi
      implicit none
      integer, intent(in) :: istra
      integer, intent(in) :: rank
      integer :: ierr
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) )
      external :: mpi_recv
#endif

      call mpi_recv(edissml(1,1,istra) , size(edissml(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(eneutrad(1,1,istra) , size(eneutrad(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(emolrad(1,1,istra) , size(emolrad(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(eionrad(1,1,istra) , size(eionrad(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(ewlda(1,1,istra) , size(ewlda(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(ewldm(1,1,istra) , size(ewldm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldna(1,1,istra) , size(wldna(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldnek(1,istra) , size(wldnek(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldnep(1,istra) , size(wldnep(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldnm(1,1,istra) , size(wldnm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldpa(1,1,istra) , size(wldpa(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldpeb(1,istra) , size(wldpeb(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldpm(1,1,istra) , size(wldpm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldpp(1,1,istra) , size(wldpp(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldra(1,1,istra) , size(wldra(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldrm(1,1,istra) , size(wldrm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldspt(1,istra) , size(wldspt(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldspta(1,1,istra) , size(wldspta(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(wldsptm(1,1,istra) , size(wldsptm(:,:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(edena_int(:,istra) , size(edena_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(edena_int_b2(:,istra) , size(edena_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(edeni_int(:,istra) , size(edeni_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(edeni_int_b2(:,istra) , size(edeni_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(edenm_int(:,istra) , size(edenm_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(edenm_int_b2(:,istra) , size(edenm_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(pdena_int(:,istra) , size(pdena_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(pdena_int_b2(:,istra) , size(pdena_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(pdeni_int(:,istra) , size(pdeni_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(pdeni_int_b2(:,istra) , size(pdeni_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(pdenm_int(:,istra) , size(pdenm_int(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      call mpi_recv(pdenm_int_b2(:,istra) , size(pdenm_int_b2(:,istra)),
     &              MPI_DOUBLE_PRECISION, rank, istra,
     &              MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

      return
      end subroutine eirene_wneutrals_recv

      subroutine eirene_wneutrals_reduce
      use eirmod_mpi
      use eirmod_cpes
     , ,only: my_pe, get_leader_comm
      implicit none
      integer :: comm, ierr
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      external :: mpi_reduce
#endif

      comm = get_leader_comm()

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, dab2, size(dab2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(dab2, dab2, size(dab2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, dib2, size(dib2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(dib2, dib2, size(dib2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, dmb2, size(dmb2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(dmb2, dmb2, size(dmb2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, pefluxa, size(pefluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(pefluxa, pefluxa, size(pefluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, pefluxm, size(pefluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(pefluxm, pefluxm, size(pefluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, pfluxa, size(pfluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(pfluxa, pfluxa, size(pfluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, pfluxm, size(pfluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(pfluxm, pfluxm, size(pfluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, refluxa, size(refluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(refluxa, refluxa, size(refluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, refluxm, size(refluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(refluxm, refluxm, size(refluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, rfluxa, size(rfluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(rfluxa, rfluxa, size(rfluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, rfluxm, size(rfluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(rfluxm, rfluxm, size(rfluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, tfluxa, size(tfluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(tfluxa, tfluxa, size(tfluxa),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, tfluxm, size(tfluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(tfluxm, tfluxm, size(tfluxm),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, srcml, size(srcml),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(srcml, srcml, size(srcml),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, tab2, size(tab2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(tab2, tab2, size(tab2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, tib2, size(tib2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(tib2, tib2, size(tib2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, tmb2, size(tmb2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(tmb2, tmb2, size(tmb2),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      if (my_pe==0) then
        call mpi_reduce(MPI_IN_PLACE, volcel, size(volcel),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      else
        call mpi_reduce(volcel, volcel, size(volcel),
     &        MPI_DOUBLE_PRECISION, MPI_SUM, 0, comm, ierr)
      endif

      return
      end subroutine eirene_wneutrals_reduce

      SUBROUTINE WNEUTRAL_FLUXES

      use eirmod_braeir
      use eirmod_CGEOM

      IMPLICIT NONE

      INTEGER, SAVE :: NCL
      DATA NCL /-1/

      IF(NLIM+NSTS.GE.NLMPGS) THEN
       WRITE(IUNOUT,*) "ERROR IN WNEUTRAL_FLUXES: ",
     w                 "NO SPATIALLY RESOLVED SURFACE TALLIES"
       WRITE(IUNOUT,*) "THE SUBROUTINE WILL BE SKIPPED"
       RETURN
      END IF

      IF(NCL.LT.0) THEN

C PREPARING TABLES IF START AND END INDICES AND CALCULATING THE
C TOTAL NUMBER OF SURFACE ELEMENTS NCL, S
       CALL PREPARE_TABLES(NCL)

C CHECKING AND PRINTING (IF NECESSARY) TABLES OF INDICES
       CALL CHECK_TABLES(NCL)

C ALLOCATES ARRAYS FOR FLUXES
       IF(.NOT.ALLOCATED(wldna_res))
     >  ALLOCATE(wldna_res(NATM,NCL), wldnm_res(NMOL,NCL),
     .           ewlda_res(NATM,NCL), ewldm_res(NMOL,NCL),
     .           ewldt_res(NCL),
     .           ewldea_res(NATM,NCL), ewldem_res(NMOL,NCL),
     .           ewldrp_res(NCL), ewldmr_res(NMOL,NCL),
     .           wldspt_res(NCL), sarea_res(NCL),
     .           wldspta_res(NATM,NCL), wldsptm_res(NMOL,NCL),
     .           wlpump_res(NSPZ,NCL) )

C CALCULATING AREAS OF SURFACE ELEMENTS
       CALL SAREA_INIT
       CALL CALCULATE_MAPPED(NCL,
     c                       WNEUTRAL_FLUXES_SAREA_UPDATE)
       CALL SAREA_POSTPROC
      END IF

      CALL FLUXES_INIT
      CALL CALCULATE_MAPPED(NCL,
     c                      WNEUTRAL_FLUXES_FLUXES_UPDATE)
C UPDATING INTEGRAL FLUXES FOR EACH SURFACE, RECOMBINATION ENERGY OF MOLECULES
      CALL FLUXES_POSTPROC !UPDATING INTEGRAL FLUXES ON EACH SURFACE

C COMPARING SUMS OVER CALCULATED SPATIALLY RESOLVED DATA WITH THOSE
C ESTIMATED DIRECTLY ON NDSS
      CALL CONTROL_SUMS(NCL)

C UNIT CONVERSION, CALCULATING OF FLUX DENSITIES ...
      CALL POSTPROCESS_FLUXES(NCL)

      RETURN
      END SUBROUTINE WNEUTRAL_FLUXES


C****************************************************************************************
C TEMPLATE SUBROUTINE FOR CALCULATING THE SURFACE DATA MAPPED TO B2 GRID
C****************************************************************************************
      SUBROUTINE CALCULATE_MAPPED(NCL,FLUX_UPDATE)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCL
      INTEGER :: IT,IS,ISS,MS,IR,I

      INTERFACE
       SUBROUTINE FLUX_UPDATE(I,J,K,L,M,IEIR)
       INTEGER,INTENT(IN) :: I,J,K,L,M,IEIR
       END SUBROUTINE FLUX_UPDATE
      END INTERFACE

      ISS = 0
      DO I = 1, NNNMTI_EIR
        IT = TNMTI_EIR(I)
        IS = SNMTI_EIR(I)
        IF (IT < 1) CYCLE
        MS = INSPAT(IS,IT)
        IF(MS.LT.1.OR.MS.GT.NNNMTI_EIR) GOTO 100
        IR = INMTI_EIR(I)
        MS=NLIM+NSTS+INSPAT(IS,IT)
        IF(MS.LT.1.OR.MS.GT.NLMPGS) GOTO 100
        IF(IR.LT.0.OR.IR.GT.NNNMTI_B2) GOTO 100

        CALL FLUX_UPDATE(IR,MS,ISS,IT,IS,I)

        CYCLE

  100   CONTINUE
        WRITE(IUNOUT,*) "ERROR IN WNEUTRAL_FLUXES"
        WRITE(IUNOUT,*) "LOST ELEMENT ON SURFACE ",ISS
        WRITE(IUNOUT,*) "IT, IS, IR, MS: ",
     w                   IT, IS, IR, MS
        WRITE(IUNOUT,*) "NNNMTI_EIR, NNNMTI_B2, NLMPGS: ",
     w                   NNNMTI_EIR, NNNMTI_B2, NLMPGS
      END DO
      END SUBROUTINE CALCULATE_MAPPED

C****************************************************************************************
C CHECKING AND PRINTING (IF NECESSARY) TABLES OF INDICES
C****************************************************************************************
      SUBROUTINE CHECK_TABLES(NCL)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCL
      INTEGER :: ISS,NT

C CHECK TABLES
      IF(TRCINT) THEN
        WRITE(IUNOUT,*) "NCL ", NCL
        DO ISS=1,NSTS
          IF(eirdiag_nds_ind(ISS).GT.0) THEN
            write (iunout,'(a,i6,a)') 'SURFACE',ISS
            write (iunout,'(a,i6)')
     .           'eirdiag_nds_ind',eirdiag_nds_ind(iss)
            nt = eirdiag_nds_ind(iss+1) - eirdiag_nds_ind(iss) - 1
            write (iunout,'(a,i6)') 'NO. OF FACES ',nt
          END IF
        END DO
      END IF
      END SUBROUTINE CHECK_TABLES

C**********************************************************************************************
C SUBROUTINES FOR CALCULATING THE FLUXES MAPPED TO B2 GRID
C**********************************************************************************************
C
C INITIALIZE FLUXES
C
      SUBROUTINE FLUXES_INIT

      wldna_res=0._DP
      wldnm_res=0._DP
      ewlda_res=0._DP
      ewldm_res=0._DP
      ewldt_res=0._DP
      ewldea_res=0._DP
      ewldem_res=0._DP
      ewldrp_res=0._DP
      ewldmr_res=0._DP
      wldspt_res=0._DP
      wldspta_res=0._DP
      wldsptm_res=0._DP
      wlpump_res=0._DP

      END SUBROUTINE FLUXES_INIT

C
C UPDATING INTEGRAL FLUXES ON EACH SURFACE,
C ENERGY FLUX DUE TO DISSOCIATION INTO MOLECULES
C
      SUBROUTINE FLUXES_POSTPROC

      IMPLICIT NONE

      INTEGER :: IS,ISS,I,IND

      DO IS=1,NSTS
       IND=eirdiag_nds_ind(IS)
       IF(IND.LT.0) CYCLE
       ISS=NLIM+IS
       ewldt_res(IND)=0._DP
       ewldrp_res(IND)=0._DP
       wldspt_res(IND)=0._DP
       DO I=1,NATM
        IF (LPOTAT) wldna_res(I,IND)=POTAT(I,ISS)
        IF (LEOTAT) then
          ewlda_res(I,IND)=EOTAT(I,ISS)
          ewldt_res(IND)  =ewldt_res(IND)+EOTAT(I,ISS)
        END IF
        ewldea_res(I,IND)= 0._DP
        IF (LERFAAT) THEN
          ewldea_res(I,IND)=ewldea_res(I,IND)+ERFAAT(I,ISS)
          ewldt_res(IND)   =ewldt_res(IND)-ERFAAT(I,ISS)
        END IF
        IF (LERFMAT) THEN
          ewldea_res(I,IND)=ewldea_res(I,IND)+ERFMAT(I,ISS)
          ewldt_res(IND)   =ewldt_res(IND)-ERFMAT(I,ISS)
        END IF
        IF (LERFIAT) THEN
          ewldea_res(I,IND)=ewldea_res(I,IND)+ERFIAT(I,ISS)
          ewldt_res(IND)   =ewldt_res(IND)-ERFIAT(I,ISS)
        END IF
        IF (LERFPHAT) THEN
          ewldea_res(I,IND)=ewldea_res(I,IND)+ERFPHAT(I,ISS)
          ewldt_res(IND)   =ewldt_res(IND)-ERFPHAT(I,ISS)
        END IF
        IF (LERFPAT) ewldrp_res(IND)=
     .               ewldrp_res(IND)+ERFPAT(I,ISS)
        wldspta_res(I,IND)= 0._DP
        IF (LSPTAAT) wldspta_res(I,IND)=
     .               wldspta_res(I,IND)+SPTAAT(I,ISS)
        IF (LSPTMAT) wldspta_res(I,IND)=
     .               wldspta_res(I,IND)+SPTMAT(I,ISS)
        IF (LSPTIAT) wldspta_res(I,IND)=
     .               wldspta_res(I,IND)+SPTIAT(I,ISS)
        IF (LSPTPAT) wldspta_res(I,IND)=
     .               wldspta_res(I,IND)+SPTPAT(I,ISS)
        IF (LSPTPHAT) wldspta_res(I,IND)=
     .               wldspta_res(I,IND)+SPTPHAT(I,ISS)
       END DO
       DO I=1,NMOL
        IF (LPOTML) wldnm_res(I,IND)= POTML(I,ISS)
        IF (LEOTML) THEN
          ewldm_res(I,IND)= EOTML(I,ISS)
          ewldt_res(IND)  = ewldt_res(IND)+EOTML(I,ISS)
        END IF
        ewldem_res(I,IND)= 0._DP
        IF (LERFAML) THEN
          ewldem_res(I,IND)=ewldem_res(I,IND)+ERFAML(I,ISS)
          ewldt_res(IND)   =ewldt_res(IND)-ERFAML(I,ISS)
        END IF
        IF (LERFMML) THEN
          ewldem_res(I,IND)=ewldem_res(I,IND)+ERFMML(I,ISS)
          ewldt_res(IND)   =ewldt_res(IND)-ERFMML(I,ISS)
        END IF
        IF (LERFIML) THEN
          ewldem_res(I,IND)=ewldem_res(I,IND)+ERFIML(I,ISS)
          ewldt_res(IND)   =ewldt_res(IND)-ERFIML(I,ISS)
        END IF
        IF (LERFPHML) THEN
          ewldem_res(I,IND)=ewldem_res(I,IND)+ERFPHML(I,ISS)
          ewldt_res(IND)   =ewldt_res(IND)-ERFPHML(I,ISS)
        END IF
        IF (LERFPML) ewldrp_res(IND)=
     .               ewldrp_res(IND)+ERFPML(I,ISS)
        ewldmr_res(I,IND)= 0._DP
        IF (LPRFAML) ewldmr_res(I,IND)=
     .               ewldmr_res(I,IND)+PRFAML(I,ISS)
        IF (LPRFPML) THEN
          ewldmr_res(I,IND)=ewldmr_res(I,IND)+PRFPML(I,ISS)
          IF (NCHARM(I).EQ.2)
     .     ewldt_res(IND)=ewldt_res(IND)+DISS_POT_H2*PRFPML(I,ISS)
        END IF
        wldsptm_res(I,IND)= 0._DP
        IF (LSPTAML) wldsptm_res(I,IND)=
     .               wldsptm_res(I,IND)+SPTAML(I,ISS)
        IF (LSPTMML) wldsptm_res(I,IND)=
     .               wldsptm_res(I,IND)+SPTMML(I,ISS)
        IF (LSPTIML) wldsptm_res(I,IND)=
     .               wldsptm_res(I,IND)+SPTIML(I,ISS)
        IF (LSPTPML) wldsptm_res(I,IND)=
     .               wldsptm_res(I,IND)+SPTPML(I,ISS)
        IF (LSPTPHML) wldsptm_res(I,IND)=
     .               wldsptm_res(I,IND)+SPTPHML(I,ISS)
        IF (LPRFAML .and. NCHARM(I).EQ.2)
     .   ewldt_res(IND)=ewldt_res(IND)+DISS_POT_H2*PRFAML(I,ISS)
       END DO
       DO I=1,NION
        IF (LEOTIO) ewldt_res(IND)=ewldt_res(IND)+EOTIO(I,ISS)
        IF (LERFAIO) ewldt_res(IND)=ewldt_res(IND)-ERFAIO(I,ISS)
        IF (LERFMIO) ewldt_res(IND)=ewldt_res(IND)-ERFMIO(I,ISS)
        IF (LERFIIO) ewldt_res(IND)=ewldt_res(IND)-ERFIIO(I,ISS)
        IF (LERFPHIO) ewldt_res(IND)=ewldt_res(IND)-ERFPHIO(I,ISS)
       END DO
       DO I=1,NSPZ
        IF (LSPUMP) wlpump_res(I,IND)=SPUMP(I,ISS)
       END DO
       IF (LSPTTOT) wldspt_res(IND)= wldspt_res(IND)+SPTTOT(ISS)
      END DO

C DISSOCIATION ENERGY OF MOLECULES: ENERGY FLUX DUE TO
C RECOMBINATION OF ATOMS INTO MOLECULES
      DO I=1,NMOL
        IF(NCHARM(I).EQ.2) THEN !HYDROGEN MOLECULE
          ewldmr_res(I,:)=ewldmr_res(I,:)*DISS_POT_H2
        ELSE
          ewldmr_res(I,:)=0._DP
        END IF
      END DO

      END SUBROUTINE FLUXES_POSTPROC

C
C
C
      SUBROUTINE SAREA_INIT

      IMPLICIT NONE

      IF(.NOT.ALLOCATED(AMARK)) ALLOCATE(AMARK(NLMPGS,NSTS))
      SAREA_RES=0._DP
      AMARK=.TRUE.

      END SUBROUTINE SAREA_INIT

C
C
C
      SUBROUTINE SAREA_POSTPROC

      IMPLICIT NONE

      INTEGER :: IS,IND,ISS

      IF(ALLOCATED(AMARK)) DEALLOCATE(AMARK)
C TOTAL AREAS OF THE SURFACES
      DO IS=1,NSTS
       IND=eirdiag_nds_ind(IS)
       IF(IND.LT.0) CYCLE
       ISS=NLIM+IS
       SAREA_RES(IND)=SAREA(ISS)
      END DO
      SAREA_RES=SAREA_RES*1.0e-4_DP !CONVERTING FROM CM**2 TO M**2

      END SUBROUTINE SAREA_POSTPROC

C****************************************************************************************
C COMPARING SUMS OVER CALCULATED SPATIALLY RESOLVED DATA WITH THOSE
C ESTIMATED DIRECTLY ON NDSS
C****************************************************************************************
      SUBROUTINE CONTROL_SUMS(NCL)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCL
      INTEGER :: NCL2,IS,ISS,IA,IE,I,IND
      REAL(DP) :: WLDNA_SUM(NATM,NSTS),WLDNM_SUM(NMOL,NSTS),
     r            EWLDA_SUM(NATM,NSTS),EWLDM_SUM(NMOL,NSTS),
     r            EWLDEA_SUM(NATM,NSTS),EWLDEM_SUM(NMOL,NSTS),
     r            EWLDRP_SUM(NSTS),EWLDMR_SUM(NMOL,NSTS),
     r            EWLDT_SUM(NSTS),WLDSPT_SUM(NSTS),SAREA_SUM(NSTS),
     r            WLDSPTA_SUM(NATM,NSTS),WLDSPTM_SUM(NMOL,NSTS),
     r            WLPUMP_SUM(NSPZ,NSTS)
      REAL(DP) :: FSUM
      CHARACTER(6) :: VNAME

      NCL2=0
      DO IS=1,NSTS
        IF(eirdiag_nds_ind(IS)+1.eq.eirdiag_nds_ind(IS+1)) THEN
          NCL2=NCL2+1
          CYCLE !THIS SURFACE IS SKIPPED
        END IF
        ISS=NLIM+IS
        IND=eirdiag_nds_ind(IS)
        IA=IND+1
        IE=eirdiag_nds_ind(IS+1)-1
        SAREA_SUM(IS)=sum(sarea_res(IA:IE))
        NCL2=NCL2+IE-IA+2
        VNAME=' SAREA'
        CALL CHECK_SUM(SAREA_SUM(IS),SAREA_RES(IND),VNAME,IS,1)
        DO I=1,NATM
          WLDNA_SUM(I,IS)=sum(wldna_res(I,IA:IE))
          VNAME=' WLDNA'
          CALL CHECK_SUM(WLDNA_SUM(I,IS),WLDNA_RES(I,IND),VNAME,IS,I)
          EWLDA_SUM(I,IS)=sum(ewlda_res(I,IA:IE))
          VNAME=' EWLDA'
          CALL CHECK_SUM(EWLDA_SUM(I,IS),EWLDA_RES(I,IND),VNAME,IS,I)
          EWLDEA_SUM(I,IS)=sum(ewldea_res(I,IA:IE))
          VNAME='EWLDEA'
          CALL CHECK_SUM(EWLDEA_SUM(I,IS),EWLDEA_RES(I,IND),VNAME,IS,I)
          WLDSPTA_SUM(I,IS)=sum(wldspta_res(I,IA:IE))
          VNAME='WLDSPA'
          CALL CHECK_SUM(WLDSPTA_SUM(I,IS),WLDSPTA_RES(I,IND),
     .                   VNAME,IS,I)
        END DO
        DO I=1,NMOL
          WLDNM_SUM(I,IS)=sum(wldnm_res(I,IA:IE))
          VNAME=' WLDNM'
          CALL CHECK_SUM(WLDNM_SUM(I,IS),WLDNM_RES(I,IND),VNAME,IS,I)
          EWLDM_SUM(I,IS)=sum(ewldm_res(I,IA:IE))
          VNAME=' EWLDM'
          CALL CHECK_SUM(EWLDM_SUM(I,IS),EWLDM_RES(I,IND),VNAME,IS,I)
          EWLDEM_SUM(I,IS)=sum(ewldem_res(I,IA:IE))
          VNAME='EWLDEM'
          CALL CHECK_SUM(EWLDEM_SUM(I,IS),EWLDEM_RES(I,IND),VNAME,IS,I)
          EWLDMR_SUM(I,IS)=sum(ewldmr_res(I,IA:IE))
          VNAME='EWLDMR'
          CALL CHECK_SUM(EWLDMR_SUM(I,IS),EWLDMR_RES(I,IND),VNAME,IS,I)
          WLDSPTM_SUM(I,IS)=sum(wldsptm_res(I,IA:IE))
          VNAME='WLDSPM'
          CALL CHECK_SUM(WLDSPTM_SUM(I,IS),WLDSPTM_RES(I,IND),
     .                   VNAME,IS,I)
        END DO
        DO I=1,NSPZ
          WLPUMP_SUM(I,IS)=sum(wlpump_res(I,IA:IE))
          VNAME='WLPUMP'
          CALL CHECK_SUM(WLPUMP_SUM(I,IS),WLPUMP_RES(I,IND),VNAME,IS,I)
        END DO
        EWLDRP_SUM(IS)=sum(ewldrp_res(IA:IE))
        VNAME='EWLDRP'
        CALL CHECK_SUM(EWLDRP_SUM(IS),EWLDRP_RES(IND),VNAME,IS,1)
        WLDSPT_SUM(IS)=sum(wldspt_res(IA:IE))
        VNAME='WLDSPT'
        CALL CHECK_SUM(WLDSPT_SUM(IS),WLDSPT_RES(IND),VNAME,IS,1)
        EWLDT_SUM(IS)=sum(ewldt_res(IA:IE))
        VNAME=' EWLDT'
        CALL CHECK_SUM(EWLDT_SUM(IS),EWLDT_RES(IND),VNAME,IS,1)
      END DO
      IF(NCL2.NE.NCL) THEN
        WRITE(IUNOUT,*) "ERROR IN WNEUTRAL_FLUXES: ",
     w                  "NOT ALL ELEMENTS WERE CHECKED"
        WRITE(IUNOUT,*) "TOTAL ", NCL, "CHECKED ", NCL2
      END IF

      IF(TRCINT) THEN
        WRITE(IUNOUT,*)
        WRITE(IUNOUT,*) "WNEUTRAL_FLUX CONTROL SUMS"
        DO IS=1,NSTS
          IND=eirdiag_nds_ind(IS)
          IF(IND.LT.0) CYCLE
          IF(eirdiag_nds_ind(IS)+1.eq.eirdiag_nds_ind(IS+1)) CYCLE
          WRITE(IUNOUT,*)
          WRITE(IUNOUT,*) "SURFACE ", IS

          WRITE(IUNOUT,601) "VARIABLE",
     w                      "SPECIES", "STANDARD",
     w                      "RESOLVED", " DIFFERENCE"
          WRITE(IUNOUT,600) "SAREA", 0, SAREA_SUM(IS),SAREA_RES(IND),
     w                              ABS(SAREA_SUM(IS)-SAREA_RES(IND))
          DO I=1,NATM
            WRITE(IUNOUT,600) "WLDNA",I,WLDNA_RES(I,IND),WLDNA_SUM(I,IS)
     w                        ,ABS(WLDNA_RES(I,IND)-WLDNA_SUM(I,IS))
            WRITE(IUNOUT,600) "EWLDA",I,EWLDA_RES(I,IND),EWLDA_SUM(I,IS)
     w                        ,ABS(EWLDA_RES(I,IND)-EWLDA_SUM(I,IS))
            WRITE(IUNOUT,600) "EWLDEA",
     w                       I,EWLDEA_RES(I,IND),EWLDEA_SUM(I,IS),
     w                     ABS(EWLDEA_RES(I,IND)-EWLDEA_SUM(I,IS))
            WRITE(IUNOUT,600) "WLDSPTA",
     w                       I,WLDSPTA_RES(I,IND),WLDSPTA_SUM(I,IS),
     w                     ABS(WLDSPTA_RES(I,IND)-WLDSPTA_SUM(I,IS))
          END DO
          DO I=1,NMOL
            WRITE(IUNOUT,600) "WLDNM",I,WLDNM_RES(I,IND),WLDNM_SUM(I,IS)
     w                        ,ABS(WLDNM_RES(I,IND)-WLDNM_SUM(I,IS))
            WRITE(IUNOUT,600) "EWLDM",I,EWLDM_RES(I,IND),EWLDM_SUM(I,IS)
     w                        ,ABS(EWLDM_RES(I,IND)-EWLDM_SUM(I,IS))
            WRITE(IUNOUT,600) "EWLDEM",
     w                       I,EWLDEM_RES(I,IND),EWLDEM_SUM(I,IS),
     w                     ABS(EWLDEM_RES(I,IND)-EWLDEM_SUM(I,IS))
            WRITE(IUNOUT,600) "EWLDMR",I,
     w                         EWLDMR_RES(I,IND),EWLDMR_SUM(I,IS),
     w                     ABS(EWLDMR_RES(I,IND)-EWLDMR_SUM(I,IS))
            WRITE(IUNOUT,600) "WLDSPTM",
     w                       I,WLDSPTM_RES(I,IND),WLDSPTM_SUM(I,IS),
     w                     ABS(WLDSPTM_RES(I,IND)-WLDSPTM_SUM(I,IS))
          END DO
          DO I=1,NSPZ
            WRITE(IUNOUT,600) "WLPUMP",
     w                       I,WLPUMP_RES(I,IND),WLPUMP_SUM(I,IS),
     w                     ABS(WLPUMP_RES(I,IND)-WLPUMP_SUM(I,IS))
          END DO
          WRITE(IUNOUT,600) "EWLDRP",0,
     w                       EWLDRP_RES(IND),EWLDRP_SUM(IS),
     w                   ABS(EWLDRP_RES(IND)-EWLDRP_SUM(IS))
          WRITE(IUNOUT,600) "WLDSPT",0,
     w                       WLDSPT_RES(IND),WLDSPT_SUM(IS),
     w                   ABS(WLDSPT_RES(IND)-WLDSPT_SUM(IS))
          WRITE(IUNOUT,600) "EWLDT",0,
     w                       EWLDT_RES(IND),EWLDT_SUM(IS),
     w                   ABS(EWLDT_RES(IND)-EWLDT_SUM(IS))
        END DO !DO IS=1,NSTS

        WRITE(IUNOUT,*)
        WRITE(IUNOUT,*) "WNEUTRAL_FLUX: COMPARING TO OLD DIAGNOSTIC"
        DO IS=1,NSTS
          IF(eirdiag_nds_ind(IS).LT.0) CYCLE
          IF(eirdiag_nds_ind(IS)+1.eq.eirdiag_nds_ind(IS+1)) CYCLE
          ISS=NLIM+IS
          WRITE(IUNOUT,*)
          WRITE(IUNOUT,*) "SURFACE ", IS

          WRITE(IUNOUT,601) "VARIABLE",
     w                      "SPECIES", "NEW",
     w                      "OLD", " DIFFERENCE"
          DO I=1,NATM
            FSUM=WLDNA(ISS,I,0)*ELCHA
            WRITE(IUNOUT,600) "WLDNA",I,WLDNA_SUM(I,IS),FSUM,
     w                              ABS(WLDNA_SUM(I,IS)-FSUM)
          END DO
          DO I=1,NMOL
            FSUM=WLDNM(ISS,I,0)*ELCHA
            WRITE(IUNOUT,600) "WLDNM",I,WLDNM_SUM(I,IS),FSUM,
     w                              ABS(WLDNM_SUM(I,IS)-FSUM)
          END DO
          FSUM=SUM(EWLDA_SUM(:,IS))+SUM(EWLDM_SUM(:,IS))-
     w        SUM(EWLDEA_SUM(:,IS))-SUM(EWLDEM_SUM(:,IS))
          WRITE(IUNOUT,600) "WLDNEK",0,FSUM,wldnek(ISS,0),
     w                             ABS(FSUM-wldnek(ISS,0))
          WRITE(IUNOUT,600) "WLDPEB",0,EWLDRP_SUM(IS),wldpeb(ISS,0),
     w                             ABS(EWLDRP_SUM(IS)-wldpeb(ISS,0))
          FSUM=SUM(EWLDMR_SUM(:,IS))
          WRITE(IUNOUT,600) "WLDNEP",0,FSUM,wldnep(ISS,0),
     w                             ABS(FSUM-wldnep(ISS,0))
        END DO !DO IS=1,NSTS
      END IF !IF(TRCINT)
      RETURN
  600 FORMAT(1X,A12,I8,1P,3E16.7)
  601 FORMAT(1X,A12,A8,3A16)

      END SUBROUTINE CONTROL_SUMS

C******************************************************************************
C CHECKING IF TWO VARIABLES ARE EQUAL AND PRINTING ERROR MESSAGE IF NOT
C******************************************************************************
      SUBROUTINE CHECK_SUM(FSUM1,FSUM2,NAME,IS,ISPC)

      IMPLICIT NONE

      REAL(DP),INTENT(IN) :: FSUM1,FSUM2
      CHARACTER(6),INTENT(IN) :: NAME
      INTEGER,INTENT(IN) :: IS,ISPC
      REAL(DP) :: EPS=1.0e-8_DP

      IF( ABS(FSUM1-FSUM2).GT.EPS*ABS(FSUM1+FSUM2).AND.
     w   (ABS(FSUM1).GT.EPS.OR.ABS(FSUM2).GT.EPS) ) THEN
        IF(ISPC.GT.0) THEN
          WRITE(IUNOUT,'(a,a6,a,i4,a,i3,a)')
     .                    "ERROR IN WNEUTRAL_FLUXES: ",NAME,
     w                    " FOR SURFACE ", IS, " SPECIES ", ISPC,
     w                    " IS INCORRECT"
          WRITE(IUNOUT,*) "FSUM_RES=",FSUM1, " FSUM_STAND=", FSUM2
        ELSE
          WRITE(IUNOUT,'(a,a6,a,i4,a)')
     .                    "ERROR IN WNEUTRAL_FLUXES: ",NAME,
     w                    " FOR SURFACE ", IS, " IS INCORRECT"
          WRITE(IUNOUT,*) "FSUM_RES=",FSUM1, " FSUM_STAND=", FSUM2
        END IF
      END IF

      RETURN
      END SUBROUTINE CHECK_SUM

C****************************************************************************************
C PREPARING THE TABLES OF START AND END INDICES AND CALCULATING THE
C TOTAL NUMBER OF SURFACE ELEMENTS NCL
C****************************************************************************************
      SUBROUTINE PREPARE_TABLES(NCL)

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: NCL
      INTEGER :: ISS, JSS, MXNUM, NEINS, IP1, J,
     .           I, NT, ILA, nt_red, i1, i2, icl
      integer, allocatable :: fc_long(:), fc_short(:), fc_mult(:)
      INTEGER, ALLOCATABLE :: ISTART(:), IVX(:,:), KNOTS(:)

C FIND OUT THE NUMBER OF TRIANGLES IN EACH NDSS
C NDSS=Non-Default Standard Surface
      NCL = 1
      MXNUM = 0
      write (*,*) 'in prepare_tables'
      DO ISS = 1, NSTS
        JSS = ISS + NLIM
        write (*,*) ' jss, numtr ',jss,SURF_TRIAN(JSS)%NUMTR
        IF (SURF_TRIAN(JSS)%NUMTR > 0) THEN
          NCL = NCL + SURF_TRIAN(JSS)%NUMTR
          MXNUM = MAX(MXNUM,SURF_TRIAN(JSS)%NUMTR)
        END IF
      END DO

      IF (ALLOCATED(TNMTI_EIR)) THEN
       IF (SIZE(TNMTI_EIR) < NCL) THEN
        DEALLOCATE(TNMTI_EIR)
        DEALLOCATE(SNMTI_EIR)
        DEALLOCATE(FNMTI_EIR)
        DEALLOCATE(ANMTI_EIR)
        DEALLOCATE(INMTI_EIR)
        DEALLOCATE(FNMTI_B2)
       END IF
      END IF
      IF (.NOT.ALLOCATED(TNMTI_EIR)) THEN
#ifndef LEGACYCOMP
        ALLOCATE(TNMTI_EIR(NCL+NSTS),SOURCE=0)
        ALLOCATE(SNMTI_EIR(NCL+NSTS),SOURCE=0)
        ALLOCATE(FNMTI_EIR(NCL+NSTS),SOURCE=0)
        ALLOCATE(ANMTI_EIR(NCL+NSTS),SOURCE=0)
        ALLOCATE(INMTI_EIR(NCL+NSTS),SOURCE=0)
#else
        ALLOCATE(TNMTI_EIR(NCL+NSTS))
        ALLOCATE(SNMTI_EIR(NCL+NSTS))
        ALLOCATE(FNMTI_EIR(NCL+NSTS))
        ALLOCATE(ANMTI_EIR(NCL+NSTS))
        ALLOCATE(INMTI_EIR(NCL+NSTS))
        TNMTI_EIR=0
        SNMTI_EIR=0
        FNMTI_EIR=0
        ANMTI_EIR=0
        INMTI_EIR=0
#endif
      END IF
      IF (.NOT.ALLOCATED(FNMTI_B2)) THEN
#ifndef LEGACYCOMP
        ALLOCATE(FNMTI_B2(NCL+NSTS),SOURCE=0)
#else
        ALLOCATE(FNMTI_B2(NCL+NSTS))
        FNMTI_B2=0
#endif
      END IF

      ALLOCATE (FC_LONG(MXNUM))
      ALLOCATE (FC_SHORT(MXNUM))
      ALLOCATE (FC_MULT(MXNUM))
      ALLOCATE (IVX(MXNUM,3))
      ALLOCATE (KNOTS(NKNOT))
      ALLOCATE (ISTART(MXNUM))

      NNNMTI_EIR = NCL
      EIRDIAG_NDS_IND(1) = 1
      ICL = 0

C  SORT TRIANGLES ALONG NONDEFAULT SURFACES

      do iss = 1, nsts
        j = iss + nlim
        knots = 0
        nt = surf_trian(j)%numtr
        if (nt <= 1) then
          ! for integral
          eirdiag_nds_ind(iss+1) = eirdiag_nds_ind(iss) + 1
          i1 = eirdiag_nds_ind(iss)
          fnmti_b2(i1) = -1
          if (iss.lt.nsts) fnmti_b2(i1+1) = -1
          nnnmti_b2 = eirdiag_nds_ind(min(iss+1,nsts))
          cycle
        end if
C  find faces of the corresponding B2 surface
        call find_surface_parts (j, mxnum, nt, nt_red,
     .                           fc_long, fc_short, fc_mult, ivx)
!  identify vertices along the B2 surface
        do i = 1, nt_red
          i1 = ivx(i,1)
          i2 = ivx(i,2)
          knots(i1) = knots(i1) + 1  ! count the number of occurences
                                     ! of node i1
          knots(i2) = knots(i2) + 1
        end do
!  identify start or end points
!  for a noncircular surface the endpoints occur only once, inner points
!  occur twice
!  for a circular surface such as the core in single null configuration
!  each vertex occurs twice
        neins = count(knots == 1)
        if (neins /= 2) then
          write (iunout,*) ' surface ', j
          write (iunout,*) ' number of endpoint ',neins
        end if
!

        if (neins > 0) then
!  take one endpoint as starting point
          istart(1:neins) = pack([(i,i=1,nknot)],knots == 1)
          ip1 = istart(1)
          ila = 1
        else
!  use the first face as start of the surface
          ip1 = ivx(1,2)
          ila = 2
        end if

!  sort the B2 faces my means of their vertices
        call sort_faces (mxnum, ivx, nt_red, ip1, ila)

!  store data for this surface

        tnmti_eir(icl+1:icl+nt) = surf_trian(j)%itrias(1:nt)
        snmti_eir(icl+1:icl+nt) = surf_trian(j)%itrisi(1:nt)
        fnmti_eir(icl+1:icl+nt) = fc_long(1:nt)
        anmti_eir(icl+1:icl+nt) = fc_mult(1:nt)

        ! for integral
        eirdiag_nds_ind(iss+1) = eirdiag_nds_ind(iss) + nt_red + 1
        i1 = eirdiag_nds_ind(iss)+1
        i2 = eirdiag_nds_ind(iss)+nt_red
        fnmti_b2(i1-1) = -1
        fnmti_b2(i1:i2) = ivx(1:nt_red,3)
        do i = 1, nt_red
          where (fnmti_eir(icl+1:icl+nt) == ivx(i,3))
            inmti_eir(icl+1:icl+nt) = i1 + i - 1
          end where
        end do
        icl = icl + nt
        nnnmti_b2 = i2
      end do

      ncl = nnnmti_b2

      deallocate (ivx)
      deallocate (knots)
      deallocate (istart)

      END SUBROUTINE PREPARE_TABLES

C****************************************************************************************

      subroutine find_surface_parts (j, nd, ntin, ntout,
     .                               fc_long, fc_short, fc_mult, ivx)

      implicit none
      integer, intent(in) :: j, nd, ntin
      integer, intent(out) :: ntout
      integer, intent(inout) :: fc_long(nd), fc_short(nd), fc_mult(nd),
     .                          IVX(nd,3)
      INTEGER, ALLOCATABLE :: KNOTS(:), IPRTS(:)
      integer :: i, it, is, ifc, nprts, minp, jp,
     .           ip, is1, ic
      logical :: visited(3,ntri)

      ntout = 0
      fc_long = 0
      fc_short = 0

      do i = 1, ntin
        it = surf_trian(j)%itrias(i)
        is = surf_trian(j)%itrisi(i)
        fc_long(i) = ifctri(is,it)
        fc_short(i) = ifctri(is,it)
        fc_mult(i) = 1
      end do

!  remove doubles
      do i = 1, ntin
        if (fc_short(i) == 0) cycle
        ifc = fc_short(i)
        where (fc_short(i+1:ntin) == ifc)
          fc_short(i+1:ntin) = 0
        end where
      end do

      ntout = count(fc_short /= 0)

      ALLOCATE (KNOTS(NKNOT))
      ALLOCATE (IPRTS(NTOUT))

      ivx = 0
      if (ntout < ntin) then
!  remove zeros
        call remove_zeros (ntin,fc_short) !  idenify vertices
!  find vertices of faces
        ic = 0
        visited = .false.
        do i = 1, ntout
          nprts = count(fc_long==fc_short(i))
          iprts(1:nprts) = pack([(i,i=1,ntin)],
     .                          fc_long(1:ntin)==fc_short(i))
          knots = 0
          do ip = 1, nprts
            jp = iprts(ip)
            IT = SURF_TRIAN(J)%ITRIAS(JP)
            IS = SURF_TRIAN(J)%ITRISI(JP)
            IS1 = IS + 1
            IF (IS1 > 3) IS1 = 1
            knots(necke(is,it)) = knots(necke(is,it)) + 1
            knots(necke(is1,it)) = knots(necke(is1,it)) + 1
            visited(is,it) = .true.
            if (nchbar(is,it) > 0) then
              if (visited(nseite(is,it),nchbar(is,it))) then
                fc_mult(jp) = 0
              end if
            end if
          end do
          minp = minval(knots,1,knots > 0)
          ic = ic + 1
          ivx(ic,1:2) = pack([(i,i=1,nknot)],knots == minp)
          ivx(ic,3) = fc_short(i)
        end do

      else
        do i = 1, ntin
          IT = SURF_TRIAN(J)%ITRIAS(I)
          IS = SURF_TRIAN(J)%ITRISI(I)
          IS1 = IS + 1
          IF (IS1 > 3) IS1 = 1
          ivx(i,1) = necke(is,it)
          ivx(i,2) = necke(is1,it)
          ivx(i,3) = fc_long(i)
        end do
      end if

      DEALLOCATE (KNOTS)
      DEALLOCATE (IPRTS)

      end subroutine find_surface_parts

C****************************************************************************************
      subroutine remove_zeros(n,a)
      implicit none
      integer, intent(in) :: n
      integer, intent(inout) :: a(n)
      integer :: i, ian, ien, istore

!  remove zeros
      istore = 0
      ian = 1
      ien = n
      do
        do i = ian, ien
          if ((istore == 0) .and. (a(i) > 0)) cycle
          if (istore == 0) then
            istore = i
          else if (a(i) > 0) then
            a(istore) = a(i)
            a(i) = 0
            ian = istore + 1
            istore = 0
            exit
          end if
        end do
        if (i > ien) exit
      end do

      end subroutine remove_zeros

C****************************************************************************************

      subroutine sort_faces (nd, ivx, nt, ip1_in, ila)

      implicit none
      integer, intent(in) :: nd, nt, ip1_in, ila
      integer, intent(inout) :: ivx(nd,3)
      integer :: il, i, ip1, ihelp

      ip1 = ip1_in

!  sort the B2 faces my means of their vertices
      do il = ila, nt
        do i = il, nt
          if (ivx(i,1) == ip1) then
            if (i > il) then
              ihelp = ivx(il,+1)
              ivx(il,1) = ivx(i,1)
              ivx(i,1) = ihelp
              ihelp = ivx(il,2)
              ivx(il,2) = ivx(i,2)
              ivx(i,2) = ihelp
              ihelp = ivx(il,3)
              ivx(il,3) = ivx(i,3)
              ivx(i,3) = ihelp
            end if
            ip1 = ivx(il,2)
            exit
          elseif (ivx(i,2) == ip1) then
            if (i > il) then
              ihelp = ivx(il,1)
              ivx(il,1) = ivx(i,2)
              ivx(i,2) = ihelp
              ihelp = ivx(il,2)
              ivx(il,2) = ivx(i,1)
              ivx(i,1) = ihelp
              ihelp = ivx(il,3)
              ivx(il,3) = ivx(i,3)
              ivx(i,3) = ihelp
            else if (i == il) then
              ihelp = ivx(il,1)
              ivx(il,1) = ivx(il,2)
              ivx(il,2) = ihelp
            end if
            ip1 = ivx(il,2)
            exit
          end if
        end do
      end do

      end subroutine sort_faces


C******************************************************************************
C POST-PROCESSING OF THE CALCULATED FLUXES:
C UNIT CONVERSION, CALCULATING OF FLUX DENSITIES ...
C******************************************************************************
      SUBROUTINE POSTPROCESS_FLUXES(NCL)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NCL
      INTEGER :: IS,I,IND,IA,IE,ISS
      REAL(DP) :: F_CNV

C UNIT CONVERSION: AMPERE TO S-1
      F_CNV=1./ELCHA
      wldna_res=F_CNV*wldna_res
      wldnm_res=F_CNV*wldnm_res

C CALCULATING THE FLUX DENSITIES
      DO ISS=1,NSTS
       IND=eirdiag_nds_ind(ISS)
       IF(IND.LT.0) CYCLE
       IA=IND+1
       IE=eirdiag_nds_ind(ISS+1)-1
       DO IS=IA,IE
        IF(.NOT.sarea_res(IS).GT.0._DP) CYCLE
        DO I=1,NATM
         wldna_res(I,IS)=wldna_res(I,IS)/sarea_res(IS)
         ewlda_res(I,IS)=ewlda_res(I,IS)/sarea_res(IS)
         ewldea_res(I,IS)=ewldea_res(I,IS)/sarea_res(IS)
         wldspta_res(I,IS)=wldspta_res(I,IS)/sarea_res(IS)
        END DO
        DO I=1,NMOL
         wldnm_res(I,IS)=wldnm_res(I,IS)/sarea_res(IS)
         ewldm_res(I,IS)=ewldm_res(I,IS)/sarea_res(IS)
         ewldem_res(I,IS)=ewldem_res(I,IS)/sarea_res(IS)
         ewldmr_res(I,IS)=ewldmr_res(I,IS)/sarea_res(IS)
         wldsptm_res(I,IS)=wldsptm_res(I,IS)/sarea_res(IS)
        END DO
        DO I=1,NSPZ
         wlpump_res(I,IS)=wlpump_res(I,IS)/sarea_res(IS)
        END DO
        ewldt_res(IS)=ewldt_res(IS)/sarea_res(IS)
        ewldrp_res(IS)=ewldrp_res(IS)/sarea_res(IS)
        wldspt_res(IS)=wldspt_res(IS)/sarea_res(IS)
       END DO
      END DO
      RETURN
      END SUBROUTINE POSTPROCESS_FLUXES

      !
      ! CALCULATING VOLUME INTEGRALS OVER EIRENE TALLIES:
      ! TOTAL NUMBER OF PARTICLES AND ENERGY
      !
      SUBROUTINE UPDATE_INTEGRALS(ISTRA)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA
      INTEGER :: I, IN, JATM, JMOL, JION
      REAL(DP) :: PP,EE
      LOGICAL :: LHELP
      EXTERNAL :: EIRENE_LEER

      LHELP = ISTRA.EQ.0
      IF (.NOT.LHELP) LHELP = NLSRON(ISTRA)

      ! ATOMS
      PDENA_INT(:,ISTRA)=0.0
      PDENA_INT_B2(:,ISTRA)=0.0
      EDENA_INT(:,ISTRA)=0.0
      EDENA_INT_B2(:,ISTRA)=0.0
      IF (LHELP) THEN
       DO JATM=1,NATM
        DO IN=1,NCIP
         IF (LPDENA) PP=PDENA(JATM,IN)*VOL(IN)
         IF (LEDENA) EE=EDENA(JATM,IN)*VOL(IN)
         IF (LPDENA) PDENA_INT(JATM,ISTRA)=PDENA_INT(JATM,ISTRA)+PP
         IF (LEDENA) EDENA_INT(JATM,ISTRA)=EDENA_INT(JATM,ISTRA)+EE
         IF (LPDENA)
     >     PDENA_INT_B2(JATM,ISTRA)=PDENA_INT_B2(JATM,ISTRA)+PP
         IF (LEDENA)
     >     EDENA_INT_B2(JATM,ISTRA)=EDENA_INT_B2(JATM,ISTRA)+EE
        END DO
        EDENA_INT(JATM,ISTRA)=EDENA_INT(JATM,ISTRA)*ELCHA
        EDENA_INT_B2(JATM,ISTRA)=EDENA_INT_B2(JATM,ISTRA)*ELCHA
       END DO
      END IF

      ! MOLECULES
      PDENM_INT(:,ISTRA)=0.0
      PDENM_INT_B2(:,ISTRA)=0.0
      EDENM_INT(:,ISTRA)=0.0
      EDENM_INT_B2(:,ISTRA)=0.0
      IF (LHELP) THEN
       DO JMOL=1,NMOL
        DO IN=1,NCIP
         IF (LPDENM) PP=PDENM(JMOL,IN)*VOL(IN)
         IF (LEDENM) EE=EDENM(JMOL,IN)*VOL(IN)
         IF (LPDENM) PDENM_INT(JMOL,ISTRA)=PDENM_INT(JMOL,ISTRA)+PP
         IF (LEDENM) EDENM_INT(JMOL,ISTRA)=EDENM_INT(JMOL,ISTRA)+EE
         IF (LPDENM)
     >     PDENM_INT_B2(JMOL,ISTRA)=PDENM_INT_B2(JMOL,ISTRA)+PP
         IF (LEDENM)
     >     EDENM_INT_B2(JMOL,ISTRA)=EDENM_INT_B2(JMOL,ISTRA)+EE
        END DO
        EDENM_INT(JMOL,ISTRA)=EDENM_INT(JMOL,ISTRA)*ELCHA
        EDENM_INT_B2(JMOL,ISTRA)=EDENM_INT_B2(JMOL,ISTRA)*ELCHA
       END DO
      END IF

      ! TEST IONS
      PDENI_INT(:,ISTRA)=0.0
      PDENI_INT_B2(:,ISTRA)=0.0
      EDENI_INT(:,ISTRA)=0.0
      EDENI_INT_B2(:,ISTRA)=0.0
      IF (LHELP) THEN
       DO JION=1,NION
        DO IN=1,NCIP
         IF (LPDENI) PP=PDENI(JION,IN)*VOL(IN)
         IF (LEDENI) EE=EDENI(JION,IN)*VOL(IN)
         IF (LPDENI) PDENI_INT(JION,ISTRA)=PDENI_INT(JION,ISTRA)+PP
         IF (LEDENI) EDENI_INT(JION,ISTRA)=EDENI_INT(JION,ISTRA)+EE
         IF (LPDENI)
     >     PDENI_INT_B2(JION,ISTRA)=PDENI_INT_B2(JION,ISTRA)+PP
         IF (LEDENI)
     >     EDENI_INT_B2(JION,ISTRA)=EDENI_INT_B2(JION,ISTRA)+EE
        END DO
        EDENI_INT(JION,ISTRA)=EDENI_INT(JION,ISTRA)*ELCHA
        EDENI_INT_B2(JION,ISTRA)=EDENI_INT_B2(JION,ISTRA)*ELCHA
       END DO
      END IF

      IF(ISTRA.EQ.0) THEN
      ! CALCULATE (CHECK) SUM OVER ALL STRATA

       ! ATOMS
        DO JATM=1,NATM
         PP=0.0
         EE=0.0
         DO I=1,NSTRA
           IF (NLSRON(I)) THEN
             PP=PP+PDENA_INT(JATM,I)
             EE=EE+EDENA_INT(JATM,I)
           END IF
         END DO
         CALL UPDATE_INTEGRALS_ERROR('PDENA',JATM,PP,PDENA_INT(JATM,0))
         CALL UPDATE_INTEGRALS_ERROR('EDENA',JATM,EE,EDENA_INT(JATM,0))
         PP=0.0
         EE=0.0
         DO I=1,NSTRA
           IF (NLSRON(I)) THEN
             PP=PP+PDENA_INT_B2(JATM,I)
             EE=EE+EDENA_INT_B2(JATM,I)
           END IF
         END DO
         CALL UPDATE_INTEGRALS_ERROR('PDENA_B2',JATM,PP,
     .                                PDENA_INT_B2(JATM,0))
         CALL UPDATE_INTEGRALS_ERROR('EDENA_B2',JATM,EE,
     .                                EDENA_INT_B2(JATM,0))
        END DO

        ! MOLECULES
        DO JMOL=1,NMOL
         PP=0.0
         EE=0.0
         DO I=1,NSTRA
           IF (NLSRON(I)) THEN
             PP=PP+PDENM_INT(JMOL,I)
             EE=EE+EDENM_INT(JMOL,I)
           END IF
         END DO
         CALL UPDATE_INTEGRALS_ERROR('PDENM',JMOL,PP,PDENM_INT(JMOL,0))
         CALL UPDATE_INTEGRALS_ERROR('EDENM',JMOL,EE,EDENM_INT(JMOL,0))
         PP=0.0
         EE=0.0
         DO I=1,NSTRA
           IF (NLSRON(I)) THEN
             PP=PP+PDENM_INT_B2(JMOL,I)
             EE=EE+EDENM_INT_B2(JMOL,I)
           END IF
         END DO
         CALL UPDATE_INTEGRALS_ERROR('PDENM_B2',JMOL,PP,
     .                                PDENM_INT_B2(JMOL,0))
         CALL UPDATE_INTEGRALS_ERROR('EDENM_B2',JMOL,EE,
     .                                EDENM_INT_B2(JMOL,0))
        END DO

        ! TEST IONS
        DO JION=1,NION
         PP=0.0
         EE=0.0
         DO I=1,NSTRA
           IF (NLSRON(I)) THEN
             PP=PP+PDENI_INT(JION,I)
             EE=EE+EDENI_INT(JION,I)
           END IF
         END DO
         CALL UPDATE_INTEGRALS_ERROR('PDENI',JION,PP,PDENI_INT(JION,0))
         CALL UPDATE_INTEGRALS_ERROR('EDENI',JION,PP,PDENI_INT(JION,0))
         PP=0.0
         EE=0.0
         DO I=1,NSTRA
           IF (NLSRON(I)) THEN
             PP=PP+PDENI_INT_B2(JION,I)
             EE=EE+EDENI_INT_B2(JION,I)
           END IF
         END DO
         CALL UPDATE_INTEGRALS_ERROR('PDENI_B2',JION,PP,
     .                                PDENI_INT_B2(JION,0))
         CALL UPDATE_INTEGRALS_ERROR('EDENI_B2',JION,EE,
     .                                EDENI_INT_B2(JION,0))
        END DO

        WRITE(IUNOUT,*) "WNEUTRALS: UPDATE_INTEGRALS FINISHED"
        CALL EIRENE_LEER(1)

      END IF !IF(ISTRA.EQ.0) THEN
      RETURN
      END SUBROUTINE UPDATE_INTEGRALS

      !
      ! PRINTING ERROR MESSAGE FOR UPDATE INTEGRALS
      !
      SUBROUTINE UPDATE_INTEGRALS_ERROR(TALLY,ISPZ,SUM1,SUM2)

      IMPLICIT NONE

      CHARACTER(*),INTENT(IN) :: TALLY
      INTEGER,INTENT(IN) :: ISPZ
      REAL(DP),INTENT(IN) :: SUM1,SUM2
      REAL(DP) :: S

      S=SUM1+SUM2+EPS12
      IF(ABS(SUM1-SUM2)/S.LT.EPS12) RETURN !NO ERROR MESSAGE

      WRITE(IUNOUT,*) "ERROR IN WNEUTRAL, UPDATE_INTEGRALS:",
     &                " INCONSISTENT INTEGRALS FOR ",
     &                 TALLY(1:LEN_TRIM(TALLY))
      WRITE(IUNOUT,*) "  ISPZ, SUM1, SUM2 ",ISPZ,SUM1,SUM2

      RETURN
      END SUBROUTINE UPDATE_INTEGRALS_ERROR

C******************************************************************************
C CALCULATING THE AREAS OF SURFACE ELEMENTS MAPPED TO B2 GRID
C******************************************************************************
      SUBROUTINE WNEUTRAL_FLUXES_SAREA_UPDATE(IR,MS,ISS,IT,IS,IEIR)

      use eirmod_COMUSR

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IR,MS,ISS,IT,IS,IEIR

C UPDATING SURFACE AREA, EXCLUDING DOUBLE COUNTING
      sarea_res(IR)=sarea_res(IR)+SAREA(MS)*ANMTI_EIR(IEIR)
      RETURN

      END SUBROUTINE WNEUTRAL_FLUXES_SAREA_UPDATE

C
C UPDATING FLUXES, TO BE CALLED FROM CALCULATE_MAPPED
C
      SUBROUTINE WNEUTRAL_FLUXES_FLUXES_UPDATE(IR,MS,ISS,IT,IS,IEIR)

      use eirmod_COMUSR

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IR,MS,ISS,IT,IS,IEIR
      INTEGER :: I

      DO I=1,NATM
        IF (LPOTAT) wldna_res(I,IR)=wldna_res(I,IR)+POTAT(I,MS)
C INCIDENT ENERGY MINUS REFLECTED ENERGY
        IF (LEOTAT) THEN
          ewlda_res(I,IR)=ewlda_res(I,IR)+EOTAT(I,MS)
          ewldt_res(  IR)=ewldt_res(  IR)+EOTAT(I,MS)
        END IF
        IF (LERFAAT) THEN
          ewldea_res(I,IR)=ewldea_res(I,IR)+ERFAAT(I,MS)
          ewldt_res(IR)=ewldt_res(IR)-ERFAAT(I,MS)
        END IF
        IF (LERFMAT) THEN
          ewldea_res(I,IR)=ewldea_res(I,IR)+ERFMAT(I,MS)
          ewldt_res(IR)=ewldt_res(IR)-ERFMAT(I,MS)
        END IF
        IF (LERFIAT) THEN
          ewldea_res(I,IR)=ewldea_res(I,IR)+ERFIAT(I,MS)
          ewldt_res(IR)=ewldt_res(IR)-ERFIAT(I,MS)
        END IF
        IF (LERFPHAT) THEN
          ewldea_res(I,IR)=ewldea_res(I,IR)+ERFPHAT(I,MS)
          ewldt_res(IR)=ewldt_res(IR)-ERFPHAT(I,MS)
        END IF
        IF (LERFPAT) ewldrp_res(IR)=ewldrp_res(IR)+ERFPAT(I,MS)
        IF (LSPTAAT) wldspta_res(I,IR)=
     .               wldspta_res(I,IR)+SPTAAT(I,MS)
        IF (LSPTMAT) wldspta_res(I,IR)=
     .               wldspta_res(I,IR)+SPTMAT(I,MS)
        IF (LSPTIAT) wldspta_res(I,IR)=
     .               wldspta_res(I,IR)+SPTIAT(I,MS)
        IF (LSPTPAT) wldspta_res(I,IR)=
     .               wldspta_res(I,IR)+SPTPAT(I,MS)
        IF (LSPTPHAT) wldspta_res(I,IR)=
     .               wldspta_res(I,IR)+SPTPHAT(I,MS)
      END DO
      DO I=1,NMOL
        IF (LPOTML) wldnm_res(I,IR)=wldnm_res(I,IR)+POTML(I,MS)
        IF (LEOTML) THEN
          ewldm_res(I,IR)=ewldm_res(I,IR)+EOTML(I,MS)
          ewldt_res(  IR)=ewldt_res(  IR)+EOTML(I,MS)
        END IF
        IF (LERFAML) THEN
          ewldem_res(I,IR)=ewldem_res(I,IR)+ERFAML(I,MS)
          ewldt_res(IR)=ewldt_res(IR)-ERFAML(I,MS)
        END IF
        IF (LERFMML) THEN
          ewldem_res(I,IR)=ewldem_res(I,IR)+ERFMML(I,MS)
          ewldt_res(IR)=ewldt_res(IR)-ERFMML(I,MS)
        END IF
        IF (LERFIML) THEN
          ewldem_res(I,IR)=ewldem_res(I,IR)+ERFIML(I,MS)
          ewldt_res(IR)=ewldt_res(IR)-ERFIML(I,MS)
        END IF
        IF (LERFPHML) THEN
          ewldem_res(I,IR)=ewldem_res(I,IR)+ERFPHML(I,MS)
          ewldt_res(IR)=ewldt_res(IR)-ERFPHML(I,MS)
        END IF
        IF (LERFPML) ewldrp_res(IR)=ewldrp_res(IR)+ERFPML(I,MS)
        IF (LPRFAML) THEN
          ewldmr_res(I,IR)=ewldmr_res(I,IR)+PRFAML(I,MS)
          IF (NCHARM(I).EQ.2)
     .     ewldt_res(IR)=ewldt_res(IR)+DISS_POT_H2*PRFAML(I,MS)
        END IF
        IF (LPRFPML) THEN
          ewldmr_res(I,IR)=ewldmr_res(I,IR)+PRFPML(I,MS)
          IF (NCHARM(I).EQ.2)
     .     ewldt_res(IR)=ewldt_res(IR)+DISS_POT_H2*PRFPML(I,MS)
        END IF
        IF (LSPTAML) wldsptm_res(I,IR)=
     .               wldsptm_res(I,IR)+SPTAML(I,MS)
        IF (LSPTMML) wldsptm_res(I,IR)=
     .               wldsptm_res(I,IR)+SPTMML(I,MS)
        IF (LSPTIML) wldsptm_res(I,IR)=
     .               wldsptm_res(I,IR)+SPTIML(I,MS)
        IF (LSPTPML) wldsptm_res(I,IR)=
     .               wldsptm_res(I,IR)+SPTPML(I,MS)
        IF (LSPTPHML) wldsptm_res(I,IR)=
     .               wldsptm_res(I,IR)+SPTPHML(I,MS)
      END DO
      DO I=1,NION
        IF (LEOTIO) ewldt_res(IR)=ewldt_res(IR)+EOTIO(I,MS)
        IF (LERFAIO) ewldt_res(IR)=ewldt_res(IR)-ERFAIO(I,MS)
        IF (LERFMIO) ewldt_res(IR)=ewldt_res(IR)-ERFMIO(I,MS)
        IF (LERFIIO) ewldt_res(IR)=ewldt_res(IR)-ERFIIO(I,MS)
        IF (LERFPHIO) ewldt_res(IR)=ewldt_res(IR)-ERFPHIO(I,MS)
      END DO
      DO I=1,NSPZ
        IF (LSPUMP) wlpump_res(I,IR)=
     .              wlpump_res(I,IR)+SPUMP(I,MS)
      END DO
      IF (LSPTTOT) wldspt_res(IR)=wldspt_res(IR)+SPTTOT(MS)

      RETURN
      END SUBROUTINE WNEUTRAL_FLUXES_FLUXES_UPDATE

      END MODULE EIRMOD_WNEUTRALS
