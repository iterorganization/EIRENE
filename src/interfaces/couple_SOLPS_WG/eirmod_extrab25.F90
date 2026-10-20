

      module eirmod_extrab25

      use eirmod_precision
      use eirmod_parmmod
      use eirmod_CESTIM
      use eirmod_COMPRT
      use eirmod_COMUSR
      use eirmod_COMXS
      use eirmod_CTEXT
      use eirmod_CTRCEI
      use eirmod_CCONA
      use eirmod_COUTAU
      use eirmod_COMSOU
      use eirmod_CSDVI
      use eirmod_CSPEI
      use eirmod_CADGEO
      use eirmod_CLOGAU
      use eirmod_CPOLYG
      use eirmod_CTRIG
      use eirmod_CZT1
      use eirmod_CCOUPL
      use eirmod_COMSPL
      use eirmod_COMSIG
      use eirmod_CINIT
      use eirmod_CGRID
      use eirmod_CGEOM
      use eirmod_solps
      use eirmod_wneutrals
      use eirmod_braeir &
          , only : pux, puy, pvx, pvy
      use eirmod_eirbra &
          , only : fnmti_b2, ncvtot
#ifdef DIMENSIONS_MODULE
      use b2mod_dimensions  ! IGNORE
#endif
      implicit none
      private

      public :: eirene_extrab25_cleanup
      public :: eirene_extrab25_eirpbls,  eirene_extrab25_eirpbls_init
      public :: eirene_extrab25_srfprvsl, eirene_extrab25_srfprvsl_init
      public :: eirene_extrab25_alloc_mods
      public :: eirene_extrab25_braeir_init
      public :: eirene_extrab25_iniusr_init
      public :: eirene_extrab25_emissivity
      public :: write_f44, read_title, write_title

      integer, save, public :: nnlimi,nnstsi,nnatmi,nnmoli,nnioni
      integer, save, public :: nnplsi,nns

!PB  moved to EIRMOD_COMUSR as the arrays are read within the main input
!pb      ! eirpbls globals
!pb      integer, save, allocatable, public :: lkindp(:), lkindm(:), lkindi(:)

      ! srfprvls
      integer, save, public :: nsrfcls
      integer, save, public, allocatable :: msrfcls(:), lsrfcls(:,:)

      ! B2.5 neutrals parameters modifications
      integer, save, public :: bn_spcsrf
      integer, save, public, allocatable :: bl_spcsrf(:), bi_spcsrf(:)
      integer, save, public, allocatable :: bj_spcsrf(:), bsps_sgrp(:)
      real(DP), save, public, allocatable :: bsps_absr(:), bsps_trno(:)
      real(DP), save, public, allocatable :: bsps_mtri(:), bsps_tmpr(:)
      real(DP), save, public, allocatable :: bsps_trni(:), bsps_spph(:)
      real(DP), save, public, allocatable :: bsps_spch(:)
      character*8, save, public, allocatable :: bsps_mtrl(:), bsps_id(:)

#ifndef NO_B2_CHEM_SPUT
      ! remaining bits and pieces from braeir common
      real(DP), save, public, allocatable :: bchemical_sputter_yield(:)
      real(DP), save, public :: fchar_chemical
      integer, save, public :: igass_chemical, itsput_chemical, &
                               issput_chemical
#endif

      ! normals of B2.5 cell edges per triangle
      real(DP), save, public, allocatable :: plnxtri(:), plnytri(:)
      real(DP), save, public, allocatable :: pplnxtri(:), pplnytri(:)

      !flux_save
      real(DP), save, public, allocatable :: flux_save(:)

!pb 27012016
! flag indicating if subroutine iniusr is called from B2.5
      integer, public, save :: ini_iniusr=0

      !format flags for output files
      integer, public, save :: jvft44=20260128, jvft46=20170930

      !names for atoms
      character*8, save, public, allocatable :: texta(:)

#ifdef B25_EIRENE
#ifndef DIMENSIONS_MODULE
#include <DIMENSIONS.F>
#endif
      !The following should match the definitions found in KOPPLDIM.F
      !dimensions for B2.5 array allocation
      integer, public :: nxdd, nydd
      parameter (nxdd=DEF_NXD+DEF_NCUT*5, nydd=DEF_NYD)
      integer, public :: nsgmx
      parameter (nsgmx=2*DEF_NLIM+4*DEF_NYD)

      !Data for checking strata correspondence
      integer, public, allocatable :: rcpos_eir(:), rcbeg_eir(:), rcend_eir(:), &
                                      rcprt_eir(:), rcspi_eir(:), rcspe_eir(:)
      character*1, public, allocatable :: rcchr_eir(:)
#endif

      contains

      subroutine eirene_extrab25_alloc_mods
      use eirmod_parmmod
      use eirmod_eirbra
      use eirmod_braeir
      implicit none
      if (allocated(flux_save)) return
      allocate (flux_save(nstra))
      flux_save=0.d0
      end subroutine eirene_extrab25_alloc_mods


      subroutine write_f44(edition)
!c     this subroutine produces an output file ft44 for plotting in B2
!c     with the neutral densities, temperatures, and fluxes
!c     (summed up for all strata)
      use eirmod_COMPRT
      implicit none
      character*4, intent(in) :: edition
      !c*** label for fort.44 file
      character*40 :: get_Eir_hash
      external get_Eir_hash
      external eirene_neutr

      character*12 :: filename
      character*3 :: stnum
      integer :: iistra,n2,ncl,i,j,k,ix,jatm,jmol,jion
      !c
      !c*** writing on file ft44
      !c
      filename = fort_lc//'44'
      if (edition.ne.'    ') filename = trim(filename)//'.'//edition
      write (iunout,*) 'Writing ',trim(filename)
      n2 = nCvp
      OPEN (UNIT=44,FILE=trim(FILENAME),ACCESS='SEQUENTIAL',FORM='FORMATTED') ! added 19980603 dpc
      rewind (44)
      WRITE(44,'(i6,2x,i8,2x,a40)') ncvp,jvft44,get_Eir_hash()
      if (jvft44.ge.20260128) then
        write(44,'(i4,2x,i4,2x,i4,2x,i4)') natmi,nmoli,nioni,nfla
      else
        write(44,'(i4,2x,i4,2x,i4)') natmi,nmoli,nioni
      end if
      !cank
      do jatm=1,natmi
        write(44,*) texts(jatm+nsph)
      end do
      do jmol=1,nmoli
        write(44,*) texts(jmol+nspa)
      end do
      do jion=1,nioni
        write(44,*) texts(jion+nspam)
      end do
      call write_title(44,'dab2',n2*natmi)
      call eirene_neutr(44,ncvp,natmi,dab2,ncvp,natm,1,1)
      call write_title(44,'tab2',n2*natmi)
      call eirene_neutr(44,ncvp,natmi,tab2,ncvp,natm,1,1)
      if (nmoli.gt.0) then
        call write_title(44,'dmb2',n2*nmoli)
        call eirene_neutr(44,ncvp,nmoli,dmb2,ncvp,nmol,1,1)
        call write_title(44,'tmb2',n2*nmoli)
        call eirene_neutr(44,ncvp,nmoli,tmb2,ncvp,nmol,1,1)
      end if
      if (nioni.gt.0) then
        call write_title(44,'dib2',n2*nioni)
        call eirene_neutr(44,ncvp,nioni,dib2,ncvp,nion,1,1)
        call write_title(44,'tib2',n2*nioni)
        call eirene_neutr(44,ncvp,nioni,tib2,ncvp,nion,1,1)
      end if
      call write_title(44,'rfluxa',n2*natmi)
      call eirene_neutr(44,ncvp,natmi,rfluxa,ncvp,natm,1,1)
      if (nmoli.gt.0) then
        call write_title(44,'rfluxm',n2*nmoli)
        call eirene_neutr(44,ncvp,nmoli,rfluxm,ncvp,nmol,1,1)
      end if
      call write_title(44,'pfluxa',n2*natmi)
      call eirene_neutr(44,ncvp,natmi,pfluxa,ncvp,natm,1,1)
      if (nmoli.gt.0) then
        call write_title(44,'pfluxm',n2*nmoli)
        call eirene_neutr(44,ncvp,nmoli,pfluxm,ncvp,nmol,1,1)
      end if
      call write_title(44,'refluxa',n2*natmi)
      call eirene_neutr(44,ncvp,natmi,refluxa,ncvp,natm,1,1)
      if (nmoli.gt.0) then
        call write_title(44,'refluxm',n2*nmoli)
        call eirene_neutr(44,ncvp,nmoli,refluxm,ncvp,nmol,1,1)
      end if
      call write_title(44,'pefluxa',n2*natmi)
      call eirene_neutr(44,ncvp,natmi,pefluxa,ncvp,natm,1,1)
      if (nmoli.gt.0) then
        call write_title(44,'pefluxm',n2*nmoli)
        call eirene_neutr(44,ncvp,nmoli,pefluxm,ncvp,nmol,1,1)
      end if
      call write_title(44,'emiss',n2)
      call eirene_neutr(44,ncvp,1,emiss,ncvp,1,1,1)
      call write_title(44,'emissmol',n2)
      call eirene_neutr(44,ncvp,1,emissmol,ncvp,1,1,1)
      !cank 960511
      !c*** save the molecule-related sources...
      if (nmoli.gt.0) then
        call write_title(44,'srcml',n2*nmoli)
        call eirene_neutr(44,ncvp,nmoli,srcml,ncvp,nmol,1,1)
        call write_title(44,'edissml',n2*nmoli)
        call eirene_neutr(44,ncvp,nmoli,edissml,ncvp,nmol,1,1)
      end if
      !c*** and the data on wall loading
      write(44,'(3i6)') nnlimi, nnstsi, nstrai
      call write_title(44,'wldnek(0)',nnlimi+nnstsi)
      call neutrs(44,wldnek,1)
      call write_title(44,'wldnep(0)',nnlimi+nnstsi)
      call neutrs(44,wldnep,1)
      call write_title(44,'wldna(0)',(nnlimi+nnstsi)*natmi)
      call neutrs(44,wldna,natmi)
      call write_title(44,'ewlda(0)',(nnlimi+nnstsi)*natmi)
      call neutrs(44,ewlda,natmi)
      call write_title(44,'wldnm(0)',(nnlimi+nnstsi)*nmoli)
      call neutrs(44,wldnm,nmoli)
      call write_title(44,'ewldm(0)',(nnlimi+nnstsi)*nmoli)
      call neutrs(44,ewldm,nmoli)
      !c*** write down the wall geometry (only valid for polygons!)
      call write_title(44,'wall_geometry',4*nnlimi)
      write (44,'(8f10.4)') (0.01*p1(1,ix), 0.01*p1(2,ix),&
                     0.01*p2(1,ix), 0.01*p2(2,ix), ix=1,nnlimi)
      !c*** from 960623 on:
      call write_title(44,'wldra(0)',(nnlimi+nnstsi)*natmi)
      call neutrs(44,wldra,natmi)
      call write_title(44,'wldrm(0)',(nnlimi+nnstsi)*nmoli)
      call neutrs(44,wldrm,nmoli)
      !c*** from 960727 on:
      if(nstrai.gt.1) then
        do iistra=1,nstrai
          write (stnum,'(i3)') iistra
          call write_title(44,'wldnek('//stnum//')',nnlimi+nnstsi)
          call neutrs(44,wldnek(1,iistra),1)
          call write_title(44,'wldnep('//stnum//')',nnlimi+nnstsi)
          call neutrs(44,wldnep(1,iistra),1)
          call write_title(44,'wldna('//stnum//')',natmi*(nnlimi+nnstsi))
          call neutrs(44,wldna(1,1,iistra),natmi)
          call write_title(44,'ewlda('//stnum//')',natmi*(nnlimi+nnstsi))
          call neutrs(44,ewlda(1,1,iistra),natmi)
          call write_title(44,'wldnm('//stnum//')',nmoli*(nnlimi+nnstsi))
          if (nmoli.gt.0) call neutrs(44,wldnm(1,1,iistra),nmoli)
          call write_title(44,'ewldm('//stnum//')',nmoli*(nnlimi+nnstsi))
          if (nmoli.gt.0) call neutrs(44,ewldm(1,1,iistra),nmoli)
          call write_title(44,'wldra('//stnum//')',natmi*(nnlimi+nnstsi))
          call neutrs(44,wldra(1,1,iistra),natmi)
          call write_title(44,'wldrm('//stnum//')',nmoli*(nnlimi+nnstsi))
          if (nmoli.gt.0) call neutrs(44,wldrm(1,1,iistra),nmoli)
        end do
      end if
      !c*** from 961228 on:
      call write_title(44,'wldpp(0)',(nnlimi+nnstsi)*nfla)
      call neutrs(44,wldpp,nfla)
      call write_title(44,'wldpa(0)',(nnlimi+nnstsi)*natmi)
      call neutrs(44,wldpa,natmi)
      call write_title(44,'wldpm(0)',(nnlimi+nnstsi)*nmoli)
      call neutrs(44,wldpm,nmoli)
      call write_title(44,'wldpeb(0)',nnlimi+nnstsi)
      call neutrs(44,wldpeb,1)
      call write_title(44,'wldspt(0)',nnlimi+nnstsi)
      call neutrs(44,wldspt,1)
      if (jvft44.ge.20170328) then
        call write_title(44,'wldspta(0)',(nnlimi+nnstsi)*natmi)
        call neutrs(44,wldspta,natmi)
        call write_title(44,'wldsptm(0)',(nnlimi+nnstsi)*nmoli)
        call neutrs(44,wldsptm,nmoli)
      end if
      if(nstrai.gt.1) then
        do iistra=1,nstrai
          write (stnum,'(i3)') iistra
          call write_title(44,'wldpp('//stnum//')',(nnlimi+nnstsi)*nfla)
          call neutrs(44,wldpp(1,1,iistra),nfla)
          call write_title(44,'wldpa('//stnum//')',(nnlimi+nnstsi)*natmi)
          call neutrs(44,wldpa(1,1,iistra),natmi)
          call write_title(44,'wldpm('//stnum//')',(nnlimi+nnstsi)*nmoli)
          if (nmoli.gt.0) call neutrs(44,wldpm(1,1,iistra),nmoli)
          call write_title(44,'wldpeb('//stnum//')',nnlimi+nnstsi)
          call neutrs(44,wldpeb(1,iistra),1)
          call write_title(44,'wldspt('//stnum//')',nnlimi+nnstsi)
          call neutrs(44,wldspt(1,iistra),1)
          if (jvft44.ge.20170328) then
            call write_title(44,'wldspta('//stnum//')',(nnlimi+nnstsi)*natmi)
            call neutrs(44,wldspta(1,1,iistra),natmi)
            call write_title(44,'wldsptm('//stnum//')',(nnlimi+nnstsi)*nmoli)
            if (nmoli.gt.0) call neutrs(44,wldsptm(1,1,iistra),nmoli)
          end if
        end do
      end if
      !c*** from 20000727 on:
      call write_title(44,'isrftype',nnlimi+nnstsi)
      write(44,'(18i4)') (isrftype(i),i=1,nnlimi)
      write(44,'(18i4)') (isrftype(nlim+i),i=1,nnstsi)
      if (NLWRMSH) THEN
      !c*** from 20051115 on:
        call write_title(44,'wlarea',nnlimi+nnstsi)
        write(44,'(1p,6e13.5)') (wlarea(i),i=1,nnlimi)
        write(44,'(1p,6e13.5)') (wlarea(nlim+i),i=1,nnstsi)
      !c*** from 20060206 on:
        k=0
        call write_title(44,'wlabsrp(A)',(nnlimi+nnstsi)*nnatmi)
        write(44,'(6a13)') (eirtxt(j+k),j=1,nnatmi)
        write(44,'(1p,6e13.5)') ((wlabsrp(j+k,i),j=1,nnatmi),i=1,nnlimi)
        write(44,'(1p,6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nnatmi),i=1,nnstsi)
        k=k+nnatmi
        call write_title(44,'wlabsrp(M)',(nnlimi+nnstsi)*nnmoli)
        write(44,'(6a13)') (eirtxt(j+k),j=1,nnmoli)
        write(44,'(1p,6e13.5)') ((wlabsrp(j+k,i),j=1,nnmoli),i=1,nnlimi)
        write(44,'(1p,6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nnmoli),i=1,nnstsi)
        k=k+nnmoli
        call write_title(44,'wlabsrp(I)',(nnlimi+nnstsi)*nnioni)
        write(44,'(6a13)') (eirtxt(j+k),j=1,nnioni)
        write(44,'(1p,6e13.5)') ((wlabsrp(j+k,i),j=1,nnioni),i=1,nnlimi)
        write(44,'(1p,6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nnioni),i=1,nnstsi)
        k=k+nnioni
        call write_title(44,'wlabsrp(P)',(nnlimi+nnstsi)*nfla)
        write(44,'(6a13)') (eirtxt(j+k),j=1,nfla)
        write(44,'(1p,6e13.5)') ((wlabsrp(j+k,i),j=1,nfla),i=1,nnlimi)
        write(44,'(1p,6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nfla),i=1,nnstsi)
        if (jvft44.ge.20170328) then
          k=0
          call write_title(44,'wlpump(A)',(nnlimi+nnstsi)*nnatmi)
          write(44,'(6a13)') (eirtxt(j+k),j=1,nnatmi)
          write(44,'(1p,6e13.5)') ((wlpump(j+k,i),j=1,nnatmi),i=1,nnlimi)
          write(44,'(1p,6e13.5)') ((wlpump(j+k,nlim+i),j=1,nnatmi),i=1,nnstsi)
          k=k+nnatmi
          call write_title(44,'wlpump(M)',(nnlimi+nnstsi)*nnmoli)
          write(44,'(6a13)') (eirtxt(j+k),j=1,nnmoli)
          write(44,'(1p,6e13.5)') ((wlpump(j+k,i),j=1,nnmoli),i=1,nnlimi)
          write(44,'(1p,6e13.5)') ((wlpump(j+k,nlim+i),j=1,nnmoli),i=1,nnstsi)
          k=k+nnmoli
          call write_title(44,'wlpump(I)',(nnlimi+nnstsi)*nnioni)
          write(44,'(6a13)') (eirtxt(j+k),j=1,nnioni)
          write(44,'(1p,6e13.5)') ((wlpump(j+k,i),j=1,nnioni),i=1,nnlimi)
          write(44,'(1p,6e13.5)') ((wlpump(j+k,nlim+i),j=1,nnioni),i=1,nnstsi)
          k=k+nnioni
          call write_title(44,'wlpump(P)',(nnlimi+nnstsi)*nfla)
          write(44,'(6a13)') (eirtxt(j+k),j=1,nfla)
          write(44,'(1p,6e13.5)') ((wlpump(j+k,i),j=1,nfla),i=1,nnlimi)
          write(44,'(1p,6e13.5)') ((wlpump(j+k,nlim+i),j=1,nfla),i=1,nnstsi)
        end if
      !c*** from 20071209 on:
        call write_title(44,'eneutrad',n2*natmi)
        call eirene_neutr(44,ncvp,natmi,eneutrad,ncvp,natm,1,1) !ak+vk
      !c*** from 20180323 on:
        if (nmoli.gt.0) then
          call write_title(44,'emolrad',n2*nmoli)
          call eirene_neutr(44,ncvp,nmoli,emolrad,ncvp,nmol,1,1)
        end if
        if (nioni.gt.0) then
          call write_title(44,'eionrad',n2*nioni)
          call eirene_neutr(44,ncvp,nioni,eionrad,ncvp,nion,1,1)
        end if

      !c*** from 20080706 on:

        call write_title(44,'eirdiag',nsts+1)
        WRITE(44,'(12i6)') (eirdiag_nds_ind(i),i=1,NSTS+1)

        NCL=eirdiag_nds_ind(NSTS+1) - 1

        if (allocated(fnmti_b2)) then
          call write_title(44,'fnmti_b2',ncl)
          WRITE(44,'(12i6)') (fnmti_b2(i),i=1,NCL)
        end if

        if(allocated(wldna_res)) then
!tf If there is an "ERROR IN WNEUTRAL_FLUXES", then it returns without
!tf allocating these arrays. I have added the if condition here to avoid
!tf segmentation fault. It might not be the best solution, since if someone
!tf reads the fort.44 file, he might not know that these values were not
!tf printed. Alternative solution would be to always allocate and initialize
!tf the variables in WNEUTRAL_FLUXES before it returns.
          call write_title(44,'sarea_res',ncl)
          WRITE(44,'(1p,6e15.7)') (sarea_res(i),i=1,NCL)
          call write_title(44,'wldna_res',ncl*natm)
          WRITE(44,'(1p,6e15.7)') ((wldna_res(j,i),j=1,natm),i=1,NCL)
          call write_title(44,'wldnm_res',ncl*nmol)
          WRITE(44,'(1p,6e15.7)') ((wldnm_res(j,i),j=1,nmol),i=1,NCL)
          call write_title(44,'ewlda_res',ncl*natm)
          WRITE(44,'(1p,6e15.7)') ((ewlda_res(j,i),j=1,natm),i=1,NCL)
          call write_title(44,'ewldm_res',ncl*nmol)
          WRITE(44,'(1p,6e15.7)') ((ewldm_res(j,i),j=1,nmol),i=1,NCL)
          call write_title(44,'ewldea_res',ncl*natm)
          WRITE(44,'(1p,6e15.7)') ((ewldea_res(j,i),j=1,natm),i=1,NCL)
          call write_title(44,'ewldem_res',ncl*nmol)
          WRITE(44,'(1p,6e15.7)') ((ewldem_res(j,i),j=1,nmol),i=1,NCL)
          call write_title(44,'ewldrp_res',ncl)
          WRITE(44,'(1p,6e15.7)') (ewldrp_res(i),i=1,NCL)
          call write_title(44,'ewldmr_res',ncl*nmol)
          WRITE(44,'(1p,6e15.7)') ((ewldmr_res(j,i),j=1,nmol),i=1,NCL)
          if (jvft44.ge.20170228) then
            call write_title(44,'wldspt_res',ncl)
            WRITE(44,'(1p,6e15.7)') (wldspt_res(i),i=1,NCL)
          end if
          if (jvft44.ge.20170328) then
            call write_title(44,'wldspta_res',ncl*natm)
            WRITE(44,'(1p,6e15.7)') ((wldspta_res(j,i),j=1,natm),i=1,NCL)
            call write_title(44,'wldsptm_res',ncl*nmol)
            WRITE(44,'(1p,6e15.7)') ((wldsptm_res(j,i),j=1,nmol),i=1,NCL)
            k=0
            call write_title(44,'wlpump_res(A)',ncl*nnatmi)
            write(44,'(6a13)') (eirtxt(j+k),j=1,nnatmi)
            write(44,'(1p,6e13.5)') ((wlpump_res(j+k,i),j=1,natm),i=1,NCL)
            k=k+nnatmi
            call write_title(44,'wlpump_res(M)',ncl*nnmoli)
            write(44,'(6a13)') (eirtxt(j+k),j=1,nnmoli)
            write(44,'(1p,6e13.5)') ((wlpump_res(j+k,i),j=1,nnmoli),i=1,NCL)
            k=k+nnmoli
            call write_title(44,'wlpump_res(I)',ncl*nnioni)
            write(44,'(6a13)') (eirtxt(j+k),j=1,nnioni)
            write(44,'(1p,6e13.5)') ((wlpump_res(j+k,i),j=1,nnioni),i=1,NCL)
            k=k+nnioni
            call write_title(44,'wlpump_res(P)',ncl*nfla)
            write(44,'(6a13)') (eirtxt(j+k),j=1,nfla)
            write(44,'(1p,6e13.5)') ((wlpump_res(j+k,i),j=1,nfla),i=1,NCL)
          end if
          if (jvft44.ge.20201006) then
            call write_title(44,'ewldt_res',ncl)
            WRITE(44,'(1p,6e15.7)') (ewldt_res(i),i=1,NCL)
          end if
        end if

      !c*** from 20081111 on. Integrals over EIRENE tallies
        call write_title(44,'pdena_int',(nstrai+1)*natm)
        WRITE(44,'(1p,6e15.7)') ((PDENA_INT(I,J),I=1,NATM),J=0,NSTRAI)
        call write_title(44,'pdenm_int',(nstrai+1)*nmol)
        WRITE(44,'(1p,6e15.7)') ((PDENM_INT(I,J),I=1,NMOL),J=0,NSTRAI)
        call write_title(44,'pdeni_int',(nstrai+1)*nion)
        WRITE(44,'(1p,6e15.7)') ((PDENI_INT(I,J),I=1,NION),J=0,NSTRAI)
        call write_title(44,'pdena_int_b2',(nstrai+1)*natm)
        WRITE(44,'(1p,6e15.7)') ((PDENA_INT_B2(I,J),I=1,NATM),J=0,NSTRAI)
        call write_title(44,'pdenm_int_b2',(nstrai+1)*nmol)
        WRITE(44,'(1p,6e15.7)') ((PDENM_INT_B2(I,J),I=1,NMOL),J=0,NSTRAI)
        call write_title(44,'pdeni_int_b2',(nstrai+1)*nion)
        WRITE(44,'(1p,6e15.7)') ((PDENI_INT_B2(I,J),I=1,NION),J=0,NSTRAI)
        call write_title(44,'edena_int',(nstrai+1)*natm)
        WRITE(44,'(1p,6e15.7)') ((EDENA_INT(I,J),I=1,NATM),J=0,NSTRAI)
        call write_title(44,'edenm_int',(nstrai+1)*nmol)
        WRITE(44,'(1p,6e15.7)') ((EDENM_INT(I,J),I=1,NMOL),J=0,NSTRAI)
        call write_title(44,'edeni_int',(nstrai+1)*nion)
        WRITE(44,'(1p,6e15.7)') ((EDENI_INT(I,J),I=1,NION),J=0,NSTRAI)
        call write_title(44,'edena_int_b2',(nstrai+1)*natm)
        WRITE(44,'(1p,6e15.7)') ((EDENA_INT_B2(I,J),I=1,NATM),J=0,NSTRAI)
        call write_title(44,'edenm_int_b2',(nstrai+1)*nmol)
        WRITE(44,'(1p,6e15.7)') ((EDENM_INT_B2(I,J),I=1,NMOL),J=0,NSTRAI)
        call write_title(44,'edeni_int_b2',(nstrai+1)*nion)
        WRITE(44,'(1p,6e15.7)') ((EDENI_INT_B2(I,J),I=1,NION),J=0,NSTRAI)

      else  ! Run with no triangular mesh
      !c*** from 20071209 on:
        call write_title(44,'eneutrad',n2*natmi)
        call eirene_neutr(44,ncvp,natmi,eneutrad,ncvp,natm,1,1) !ak+vk
      !c*** from 20180323 on:
        if (nmoli.gt.0) then
          call write_title(44,'emolrad',n2*nmoli)
          call eirene_neutr(44,ncvp,nmoli,emolrad,ncvp,nmol,1,1)
        end if
        if (nioni.gt.0) then
          call write_title(44,'eionrad',n2*nioni)
          call eirene_neutr(44,ncvp,nioni,eionrad,ncvp,nion,1,1)
        end if
      endif ! NLWRMSH
      rewind (44)
      close (44)

      filename = fort_lc//'44_aver'
      if (edition.ne.'    ') filename = trim(filename)//'.'//edition
      if (aver_frac.ne.0._DP) then                                             !som 02.04.19
        write (iunout,*) 'Writing ',trim(filename)
        OPEN (UNIT=44,FILE=trim(filename),ACCESS='SEQUENTIAL',FORM='FORMATTED')
        rewind (44)
        WRITE(44,'(i6,2x,i8,2x,a40)') ncvp,jvft44,get_Eir_hash()
        write(44,'(e16.8)') aver_frac
        write(44,'(i4)') nfla
        call write_title(44,'wldnek_aver(0)',nnlimi+nnstsi)
        call remove_small_values(nlimps, wldnek_aver, EPS60)
        call neutrs(44,wldnek_aver,1)
        call write_title(44,'wldnep_aver(0)',nnlimi+nnstsi)
        call remove_small_values(nlimps, wldnep_aver, EPS60)
        call neutrs(44,wldnep_aver,1)
        call write_title(44,'wldna_aver(0)',(nnlimi+nnstsi)*natmi)
        call remove_small_values(nlimps*natmi, wldna_aver, EPS60)
        call neutrs(44,wldna_aver,natmi)
        call write_title(44,'ewlda_aver(0)',(nnlimi+nnstsi)*natmi)
        call remove_small_values(nlimps*natmi, ewlda_aver, EPS60)
        call neutrs(44,ewlda_aver,natmi)
        call write_title(44,'wldnm_aver(0)',(nnlimi+nnstsi)*nmoli)
        call remove_small_values(nlimps*nmoli, wldnm_aver, EPS60)
        call neutrs(44,wldnm_aver,nmoli)
        call write_title(44,'ewldm_aver(0)',(nnlimi+nnstsi)*nmoli)
        call remove_small_values(nlimps*nmoli, ewldm_aver, EPS60)
        call neutrs(44,ewldm_aver,nmoli)
        call write_title(44,'wldra_aver(0)',(nnlimi+nnstsi)*natmi)
        call remove_small_values(nlimps*natmi, wldra_aver, EPS60)
        call neutrs(44,wldra_aver,natmi)
        call write_title(44,'wldrm_aver(0)',(nnlimi+nnstsi)*nmoli)
        call remove_small_values(nlimps*nmoli, wldrm_aver, EPS60)
        call neutrs(44,wldrm_aver,nmoli)
        call write_title(44,'wldpp_aver(0)',(nnlimi+nnstsi)*nfla)
        call remove_small_values(nlimps*nfla, wldpp_aver, EPS60)
        call neutrs(44,wldpp_aver,nfla)
        call write_title(44,'wldpa_aver(0)',(nnlimi+nnstsi)*natmi)
        call remove_small_values(nlimps*natmi, wldpa_aver, EPS60)
        call neutrs(44,wldpa_aver,natmi)
        call write_title(44,'wldpm_aver(0)',(nnlimi+nnstsi)*nmoli)
        call remove_small_values(nlimps*nmoli, wldpm_aver, EPS60)
        call neutrs(44,wldpm_aver,nmoli)
        call write_title(44,'wldpeb_aver(0)',nnlimi+nnstsi)
        call remove_small_values(nlimps, wldpeb_aver, EPS60)
        call neutrs(44,wldpeb_aver,1)
        call write_title(44,'wldspt_aver(0)',nnlimi+nnstsi)
        call remove_small_values(nlimps, wldspt_aver, EPS60)
        call neutrs(44,wldspt_aver,1)
        call write_title(44,'wldspta_aver(0)',(nnlimi+nnstsi)*natmi)
        call remove_small_values(nlimps*natmi, wldspta_aver, EPS60)
        call neutrs(44,wldspta_aver,natmi)
        call write_title(44,'wldsptm_aver(0)',(nnlimi+nnstsi)*nmoli)
        call remove_small_values(nlimps*nmoli, wldsptm_aver, EPS60)
        call neutrs(44,wldsptm_aver,nmoli)
        k=0
        call remove_small_values((nlim+nnstsi)*(nnatmi+nnmoli+nnioni+nfla), &
            &  wlpump_aver, EPS60)
        call write_title(44,'wlpump_aver(A)',(nnlimi+nnstsi)*nnatmi)
        write(44,'(1p,6e13.5)') ((wlpump_aver(j+k,i),j=1,nnatmi),i=1,nnlimi)
        write(44,'(1p,6e13.5)') ((wlpump_aver(j+k,nlim+i),j=1,nnatmi),i=1,nnstsi)
        k=k+nnatmi
        call write_title(44,'wlpump_aver(M)',(nnlimi+nnstsi)*nnmoli)
        write(44,'(1p,6e13.5)') ((wlpump_aver(j+k,i),j=1,nnmoli),i=1,nnlimi)
        write(44,'(1p,6e13.5)') ((wlpump_aver(j+k,nlim+i),j=1,nnmoli),i=1,nnstsi)
        k=k+nnmoli
        call write_title(44,'wlpump_aver(I)',(nnlimi+nnstsi)*nnioni)
        write(44,'(1p,6e13.5)') ((wlpump_aver(j+k,i),j=1,nnioni),i=1,nnlimi)
        write(44,'(1p,6e13.5)') ((wlpump_aver(j+k,nlim+i),j=1,nnioni),i=1,nnstsi)
        k=k+nnioni
        call write_title(44,'wlpump_aver(P)',(nnlimi+nnstsi)*nfla)
        write(44,'(1p,6e13.5)') ((wlpump_aver(j+k,i),j=1,nfla),i=1,nnlimi)
        write(44,'(1p,6e13.5)') ((wlpump_aver(j+k,nlim+i),j=1,nfla),i=1,nnstsi)
        rewind (44)
        close (44)
      endif

      IF (NLWRMSH) THEN
!c
!c*** writing ft46 file
!c    tallies on the triangle mesh
!c    (raw data, EIRENE units)
!c
        filename = fort_lc//'46'
        if (edition.ne.'    ') filename = trim(filename)//'.'//edition
        write (iunout,*) 'Writing ',trim(filename)
        OPEN (UNIT=46,FILE=trim(FILENAME),ACCESS='SEQUENTIAL',FORM='FORMATTED')
        rewind (46)

        write(46,'(i6,2x,i8,2x,a40)') ncvtot, jvft46, get_Eir_hash()
        write(46,'(i4,2x,i4,2x,i4)') natmi, nmoli, nioni
        do jatm=1,natmi
          write(46,*) texts(jatm+nsph)
        end do
        do jmol=1,nmoli
          write(46,*) texts(jmol+nspa)
        end do
        do jion=1,nioni
          write(46,*) texts(jion+nspam)
        end do

        if (lpdena) then
          call write_title(46,'pdena',ncvtot*natmi)
          write(46,'(1p,6e15.7)') ((PDENA(i,j),j=1,ncvtot),i=1,natmi)
        else
          call write_title(46,'pdena (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (PDENA(1,j),j=1,ncvtot)
        end if
        if (lpdenm) then
          call write_title(46,'pdenm',ncvtot*nmoli)
          write(46,'(1p,6e15.7)') ((PDENM(i,j),j=1,ncvtot),i=1,nmoli)
        else
          call write_title(46,'pdenm (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (PDENM(1,j),j=1,ncvtot)
        end if
        if (lpdeni) then
          call write_title(46,'pdeni',ncvtot*nioni)
          write(46,'(1p,6e15.7)') ((PDENI(i,j),j=1,ncvtot),i=1,nioni)
        else
          call write_title(46,'pdeni (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (PDENI(1,j),j=1,ncvtot)
        end if

        if (ledena) then
          call write_title(46,'edena',ncvtot*natmi)
          write(46,'(1p,6e15.7)') ((EDENA(i,j),j=1,ncvtot),i=1,natmi)
        else
          call write_title(46,'edena (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (EDENA(1,j),j=1,ncvtot)
        end if
        if (ledenm) then
          call write_title(46,'edenm',ncvtot*nmoli)
          write(46,'(1p,6e15.7)') ((EDENM(i,j),j=1,ncvtot),i=1,nmoli)
        else
          call write_title(46,'edenm (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (EDENM(1,j),j=1,ncvtot)
        end if
        if (ledeni) then
          call write_title(46,'edeni',ncvtot*nioni)
          write(46,'(1p,6e15.7)') ((EDENI(i,j),j=1,ncvtot),i=1,nioni)
        else
          call write_title(46,'edeni (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (EDENI(1,j),j=1,ncvtot)
        end if

        if (lvxdena) then
          call write_title(46,'vxdena',ncvtot*natmi)
          write(46,'(1p,6e15.7)') ((VXDENA(i,j),j=1,ncvtot),i=1,natmi)
        else
          call write_title(46,'vxdena (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VXDENA(1,j),j=1,ncvtot)
        end if
        if (lvxdenm) then
          call write_title(46,'vxdenm',ncvtot*nmoli)
          write(46,'(1p,6e15.7)') ((VXDENM(i,j),j=1,ncvtot),i=1,nmoli)
        else
          call write_title(46,'vxdenm (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VXDENM(1,j),j=1,ncvtot)
        end if
        if (lvxdeni) then
          call write_title(46,'vxdeni',ncvtot*nioni)
          write(46,'(1p,6e15.7)') ((VXDENI(i,j),j=1,ncvtot),i=1,nioni)
        else
          call write_title(46,'vxdeni (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VXDENI(1,j),j=1,ncvtot)
        end if

        if (lvydena) then
          call write_title(46,'vydena',ncvtot*natmi)
          write(46,'(1p,6e15.7)') ((VYDENA(i,j),j=1,ncvtot),i=1,natmi)
        else
          call write_title(46,'vydena (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VYDENA(1,j),j=1,ncvtot)
        end if
        if (lvydenm) then
          call write_title(46,'vydenm',ncvtot*nmoli)
          write(46,'(1p,6e15.7)') ((VYDENM(i,j),j=1,ncvtot),i=1,nmoli)
        else
          call write_title(46,'vydenm (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VYDENM(1,j),j=1,ncvtot)
        end if
        if (lvydeni) then
          call write_title(46,'vydeni',ncvtot*nioni)
          write(46,'(1p,6e15.7)') ((VYDENI(i,j),j=1,ncvtot),i=1,nioni)
        else
          call write_title(46,'vydeni (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VYDENI(1,j),j=1,ncvtot)
        end if

        if (lvzdena) then
          call write_title(46,'vzdena',ncvtot*natmi)
          write(46,'(1p,6e15.7)') ((VZDENA(i,j),j=1,ncvtot),i=1,natmi)
        else
          call write_title(46,'vzdena (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VZDENA(1,j),j=1,ncvtot)
        end if
        if (lvzdenm) then
          call write_title(46,'vzdenm',ncvtot*nmoli)
          write(46,'(1p,6e15.7)') ((VZDENM(i,j),j=1,ncvtot),i=1,nmoli)
        else
          call write_title(46,'vzdenm (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VZDENM(1,j),j=1,ncvtot)
        end if
        if (lvzdeni) then
          call write_title(46,'vzdeni',ncvtot*nioni)
          write(46,'(1p,6e15.7)') ((VZDENI(i,j),j=1,ncvtot),i=1,nioni)
        else
          call write_title(46,'vzdeni (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VZDENI(1,j),j=1,ncvtot)
        endif

        call write_title(46,'volumes',ntrii)
        write(46,'(1p,6e15.7)') (VOL(j),j=1,ntrii)
        call write_title(46,'pux',ncvp)
        write(46,'(1p,6e15.7)') (PUX(j),j=1,ncvp)
        call write_title(46,'puy',ncvp)
        write(46,'(1p,6e15.7)') (PUY(j),j=1,ncvp)
        call write_title(46,'pvx',ncvp)
        write(46,'(1p,6e15.7)') (PVX(j),j=1,ncvp)
        call write_title(46,'pvy',ncvp)
        write(46,'(1p,6e15.7)') (PVY(j),j=1,ncvp)

!! atom condensation
        if (lpgena) then
          call write_title(46,'pgena',ncvtot*natmi)
          write(46,'(1p,6e15.7)') ((PGENA(i,j),j=1,ncvtot),i=1,natmi)
        else
          call write_title(46,'pgena (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (PGENA(1,j),j=1,ncvtot)
        end if
        if (lvgena) then
          call write_title(46,'vgena',ncvtot*natmi)
          write(46,'(1p,6e15.7)') ((VGENA(i,j),j=1,ncvtot),i=1,natmi)
        else
          call write_title(46,'vgena (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (VGENA(1,j),j=1,ncvtot)
        end if
        if (legena) then
          call write_title(46,'egena',ncvtot*natmi)
          write(46,'(1p,6e15.7)') ((EGENA(i,j),j=1,ncvtot),i=1,natmi)
        else
          call write_title(46,'egena (not computed)',ncvtot)
          write(46,'(1p,6e15.7)') (EGENA(1,j),j=1,ncvtot)
        end if

        rewind(46)
        close (46)

        if (aver_frac46.ne.0._DP) then ! average fluxes, som 02.04.2019
          filename = fort_lc//'46_aver'
          if (edition.ne.'    ') filename = trim(filename)//'.'//edition
          write (iunout,*) 'Writing ',trim(filename)
          OPEN (UNIT=46,FILE=trim(filename),ACCESS='SEQUENTIAL',FORM='FORMATTED')
          rewind (46)
          write(46,'(2x,i8,2x,a40)') jvft46, get_Eir_hash()
          write(46,'(e16.8)') aver_frac46
          write(46,'(i6)') ncvtot
          if (lpdena) then
            call write_title(46,'pdena_aver',ncvtot*natmi)
            call remove_small_values(ncvtot*natmi, pdena_aver, EPS60)
            write(46,'(1p,6e15.7)') ((PDENA_aver(i,j),j=1,ncvtot),i=1,natmi)
          else
            call write_title(46,'pdena (not computed)',ncvtot)
            call remove_small_values(ncvtot, pdena_aver, EPS60)
            write(46,'(1p,6e15.7)') (PDENA_aver(1,j),j=1,ncvtot)
          end if
          if (lpdenm) then
            call write_title(46,'pdenm_aver',ncvtot*nmoli)
            call remove_small_values(ncvtot*nmoli, pdenm_aver, EPS60)
            write(46,'(1p,6e15.7)') ((PDENM_aver(i,j),j=1,ncvtot),i=1,nmoli)
          else
            call write_title(46,'pdenm (not computed)',ncvtot)
            call remove_small_values(ncvtot, pdenm_aver, EPS60)
            write(46,'(1p,6e15.7)') (PDENM_aver(1,j),j=1,ncvtot)
          end if

          if (ledena) then
            call write_title(46,'edena_aver',ncvtot*natmi)
            call remove_small_values(ncvtot*natmi, edena_aver, EPS60)
            write(46,'(1p,6e15.7)') ((EDENA_aver(i,j),j=1,ncvtot),i=1,natmi)
          else
            call write_title(46,'edena (not computed)',ncvtot)
            call remove_small_values(ncvtot, edena_aver, EPS60)
            write(46,'(1p,6e15.7)') (EDENA_aver(1,j),j=1,ncvtot)
          end if
          if (ledenm) then
            call write_title(46,'edenm_aver',ncvtot*nmoli)
            call remove_small_values(ncvtot*nmoli, edenm_aver, EPS60)
            write(46,'(1p,6e15.7)') ((EDENM_aver(i,j),j=1,ncvtot),i=1,nmoli)
          else
            call write_title(46,'edenm (not computed)',ncvtot)
            call remove_small_values(ncvtot, edenm_aver, EPS60)
            write(46,'(1p,6e15.7)') (EDENM_aver(1,j),j=1,ncvtot)
          end if

          if (lvxdena) then
            call write_title(46,'vxdena_aver',ncvtot*natmi)
            call remove_small_values(ncvtot*natmi, vxdena_aver, EPS60)
            write(46,'(1p,6e15.7)') ((VXDENA_aver(i,j),j=1,ncvtot),i=1,natmi)
          else
            call write_title(46,'vxdena (not computed)',ncvtot)
            call remove_small_values(ncvtot, vxdena_aver, EPS60)
            write(46,'(1p,6e15.7)') (VXDENA_aver(1,j),j=1,ncvtot)
          end if
          if (lvxdenm) then
            call write_title(46,'vxdenm_aver',ncvtot*nmoli)
            call remove_small_values(ncvtot*nmoli, vxdenm_aver, EPS60)
            write(46,'(1p,6e15.7)') ((VXDENM_aver(i,j),j=1,ncvtot),i=1,nmoli)
          else
            call write_title(46,'vxdenm (not computed)',ncvtot)
            call remove_small_values(ncvtot, vxdenm_aver, EPS60)
            write(46,'(1p,6e15.7)') (VXDENM_aver(1,j),j=1,ncvtot)
          end if

          if (lvydena) then
            call write_title(46,'vydena_aver',ncvtot*natmi)
            call remove_small_values(ncvtot*natmi, vydena_aver, EPS60)
            write(46,'(1p,6e15.7)') ((VYDENA_aver(i,j),j=1,ncvtot),i=1,natmi)
          else
            call write_title(46,'vydena (not computed)',ncvtot)
            call remove_small_values(ncvtot, vydena_aver, EPS60)
            write(46,'(1p,6e15.7)') (VYDENA_aver(1,j),j=1,ncvtot)
          end if
          if (lvydenm) then
            call write_title(46,'vydenm_aver',ncvtot*nmoli)
            call remove_small_values(ncvtot*nmoli, vydenm_aver, EPS60)
            write(46,'(1p,6e15.7)') ((VYDENM_aver(i,j),j=1,ncvtot),i=1,nmoli)
          else
            call write_title(46,'vydenm (not computed)',ncvtot)
            call remove_small_values(ncvtot, vydenm_aver, EPS60)
            write(46,'(1p,6e15.7)') (VYDENM_aver(1,j),j=1,ncvtot)
          end if

          if (lvzdena) then
            call write_title(46,'vzdena_aver',ncvtot*natmi)
            call remove_small_values(ncvtot*natmi, vzdena_aver, EPS60)
            write(46,'(1p,6e15.7)') ((VZDENA_aver(i,j),j=1,ncvtot),i=1,natmi)
          else
            call write_title(46,'vzdena (not computed)',ncvtot)
            call remove_small_values(ncvtot, vzdena_aver, EPS60)
            write(46,'(1p,6e15.7)') (VZDENA_aver(1,j),j=1,ncvtot)
          end if
          if (lvzdenm) then
            call write_title(46,'vzdenm_aver',ncvtot*nmoli)
            call remove_small_values(ncvtot*nmoli, vzdenm_aver, EPS60)
            write(46,'(1p,6e15.7)') ((VZDENM_aver(i,j),j=1,ncvtot),i=1,nmoli)
          else
            call write_title(46,'vzdenm (not computed)',ncvtot)
            call remove_small_values(ncvtot, vzdenm_aver, EPS60)
            write(46,'(1p,6e15.7)') (VZDENM_aver(1,j),j=1,ncvtot)
          end if
          rewind(46)
          close (46)
        endif

      ENDIF

!cc<<<
!c      write (iunout,'()')
!c      write (iunout,*) '%%% Eirene data in Eirene %%%',nlimps
!c      write (iunout,*) 'nnatmi,nnmoli,nnstsi,nnlimi = ', &
!c                        nnatmi,nnmoli,nnstsi,nnlimi
!c      do k=0,nstrai !{
!c       write (iunout,'()')
!c       if (nmoli.gt.0) then
!c         write (iunout,*) 'wldna, wldnm, wldra, wldrm:',k
!c         write (iunout,'(1x,20(6x,a2,i3.2))') &
!c               ('na',i,i=1,natmi),('nm',i,i=1,nmoli), &
!c               ('ra',i,i=1,natmi),('rm',i,i=1,nmoli)
!c       else
!c         write (iunout,*) 'wldna, wldra:',k
!c         write (iunout,'(1x,20(6x,a2,i3.2))') &
!c               ('na',i,i=1,natmi),('ra',i,i=1,natmi)
!c       end if
!c       do j=1,nlim+nstsi !{
!c         if(j.le.nlimi .or. j.gt.nlim) then !{
!c           hlp_pr=.false.
!c           do i=1,natmi !{
!c             hlp_pr= hlp_pr .or. wldna(j,i,k).ne.0. &
!c                                              .or. wldra(j,i,k).ne.0.
!c           end do !}
!c           do i=1,nmoli !{
!c             hlp_pr= hlp_pr .or. wldnm(j,i,k).ne.0. &
!c                                              .or. wldrm(j,i,k).ne.0.
!c           end do !}
!c           if(hlp_pr) then !{
!c             write (iunout,'(1p,i4,20e11.4)') j, &
!c                   (wldna(j,i,k),i=1,natmi),(wldnm(j,i,k),i=1,nmoli), &
!c                   (wldra(j,i,k),i=1,natmi),(wldrm(j,i,k),i=1,nmoli)
!c           end if !}
!c         end if !}
!c       end do !}
!c
!c       write (iunout,'()')
!c       write (iunout,*) 'wldpp:',k
!c       do j=nlim+1,nlim+nstsi !{
!c         write (iunout,'(1p,i4,20e9.2)') j,(wldpp(j,i,k),i=1,nfla)
!c       end do !}
!c       write (iunout,'()')
!c       if (nmoli.gt.0) then
!c         write (iunout,*) 'wldpeb, wldpa, wldpm:',k
!c       else
!c         write (iunout,*) 'wldpeb, wldpa:',k
!c       end if
!c       do j=nlim+1,nlim+nstsi !{
!c         write (iunout,'(1p,i4,20e9.2)') j,wldpeb(j,k), &
!c                    (wldpa(j,i,k),i=1,natmi),(wldpm(j,i,k),i=1,nmoli)
!c       end do !}
!c      end do !}
!cc>>>
      return

      contains

      subroutine remove_small_values(n, field, small)
      implicit none
      integer i, n
      real(DP) :: field(n), small

      do i = 1, n
        field(i) = min(-small,field(i))+max(small,field(i))
      end do
      return
      end subroutine remove_small_values

      end subroutine write_f44
      !c======================================================================

      subroutine neutrs(kard,dummy,ldmf)
      use eirmod_parmmod
      use eirmod_cadgeo
      use eirmod_clgin
      implicit none
      integer, intent(in) :: kard,ldmf
      integer is,iif
      real(DP), intent(in) :: dummy(nlimps,*)

      do iif=1,ldmf
        if(nlimi.gt.0) write(kard,'(1p,5(1x,1e15.8))') &
                       (dummy(is,iif),is=1,nlimi)
        if(nstsi.gt.0) write(kard,'(1p,5(1x,1e15.8))') &
                       (dummy(is+nlim,iif),is=1,nstsi)
      enddo
      return
      end subroutine neutrs

      subroutine eirene_extrab25_eirpbls(istr)
      !c
      !c old version : 20.02.97 15:22
      !c new version : 26.01.2011 (s.wiesen@fz-juelich.de) - f90 module
      !c
      !c======================================================================
      !c*** This subroutine provides correction of the particle sources in an
      !c*** attempt to reduce the statistical noise in the particle balance.
      !c***
      !c*** All the particle sources are classified by the nucleus kind and
      !c*** the stratum, and treated are those selected "true" in the array
      !c*** "nlpbls". For each class, the actual flux of the plasma ions to
      !c*** the target is calculated and compared with the expected flux which
      !c*** appears as particle sink in B2. All the particle sources for this
      !c*** class are then multiplied by the ratio of "B2" to "Eirene" fluxes
      !c*** in order have the sources accurately corresponding to the sinks -
      !c*** at least, on the integral measure.
      !c----------------------------------------------------------------------
      use eirmod_PARMMOD
      use eirmod_CESTIM
      use eirmod_COUTAU
      use eirmod_COMUSR
      use eirmod_EIRBRA
      use eirmod_CCOUPL
      use eirmod_COMPRT
      use eirmod_COMSPL
      implicit none
      integer :: k,istr,l,i,j
      real(dp) :: w,hlp_flb2, hlp_fleir

      !write (iunout,*) '%% eirpbls: istr= ',istr
      !write (iunout,*) 'nlpbls(istr)= ',(nlpbls(i,istr),i=1,natmi)
      !write (iunout,'(a,20i5)') 'lkindp= ',(lkindp(i),i=1,nfla)
      !write (iunout,'(a,20i5)') 'lkindm= ',(lkindm(i),i=1,nmoli)
      !write (iunout,'(a,20i5)') 'lkindi= ',(lkindi(i),i=1,nioni)
      !write (iunout,'(a,1p/(9e11.2))') 'flxspci(istr):', &
      !                                   (flxspci(l,istr),l=1,nfla)

      do k=1,natmi
        !write (iunout,*) '%% k= ',k
        srccrfc(k,istr)=1.
        if (nlpbls(k,istr)) then
          !c
          !c*** Calculate the correction factor
          !c
          srccrfc(k,istr)=0.
          hlp_flb2=0.
          hlp_fleir=0.
          do l=1,nfla
            if(lkindp(l).eq.k) then
              hlp_flb2=hlp_flb2+flxspci(l,istr)
              hlp_fleir=hlp_fleir-potpli(l,istr)
            end if
          end do  ! do l
          !write (iunout,*) 'hlp_flb2,hlp_fleir=',hlp_flb2,hlp_fleir
          if(hlp_fleir.gt.0.) then
            w=hlp_flb2/hlp_fleir
            srccrfc(k,istr)=w
            !write (iunout,*) 'Correction factor: ',w,k
            !c
            !c*** Correct the B2 sources
            !c
            do l=1,nfla
              if(lkindp(l).eq.k) then
                !write (iunout,*) 'sni,smo: k,l,istr= ',k,l,istr
                do i=1,ncip
                  sni(i,l,istr)=w*sni(i,l,istr)
                  smo(i,l,istr)=w*smo(i,l,istr)
                end do
              end if
            end do  ! do l
            !c
            !c*** Correct the Eirene tallies for this stratum
            !c
            do i=1,nrad
              if (lpdena) pdena(k,i)=w*pdena(k,i)
              if (ledena) edena(k,i)=w*edena(k,i)
              if (lpaat)  paat(k,i)=w*paat(k,i)
              if (lpmat)  pmat(k,i)=w*pmat(k,i)
              if (lpiat)  piat(k,i)=w*piat(k,i)
            end do  ! do i
            do l=1,nmoli
              if(lkindm(l).eq.k) then
                !write (iunout,*) 'vol. mol.: k,l,istr= ',k,l,istr
                do i=1,nrad
                  if (lpdenm) pdenm(l,i)=w*pdenm(l,i)
                  if (ledenm) edenm(l,i)=w*edenm(l,i)
                  if (lpaml)  paml(l,i)=w*paml(l,i)
                  if (lpmml)  pmml(l,i)=w*pmml(l,i)
                  if (lpiml)  piml(l,i)=w*piml(l,i)
                end do
              end if
            end do  ! do l
            do l=1,nioni
              if(lkindi(l).eq.k) then
                !write (iunout,*) 'vol. ion.: k,l,istr= ',k,l,istr
                do i=1,nrad
                  if (lpdeni) pdeni(l,i)=w*pdeni(l,i)
                  if (ledeni) edeni(l,i)=w*edeni(l,i)
                  if (lpaio)  paio(l,i)=w*paio(l,i)
                  if (lpmio)  pmio(l,i)=w*pmio(l,i)
                  if (lpiio)  piio(l,i)=w*piio(l,i)
                end do
              end if
            end do  ! do l
            do i=1,nlmpgs
              if (lpotat)  potat(k,i)=w*potat(k,i)
              if (lprfaat) prfaat(k,i)=w*prfaat(k,i)
              if (lprfmat) prfmat(k,i)=w*prfmat(k,i)
              if (lprfiat) prfiat(k,i)=w*prfiat(k,i)
              if (leotat)  eotat(k,i)=w*eotat(k,i)
              if (lerfaat) erfaat(k,i)=w*erfaat(k,i)
              if (lerfmat) erfmat(k,i)=w*erfmat(k,i)
              if (lerfiat) erfiat(k,i)=w*erfiat(k,i)
              if (lprfpat) prfpat(k,i)=w*prfpat(k,i)
              if (lerfpat) erfpat(k,i)=w*erfpat(k,i)
            end do  ! do i
            do l=1,nmoli
              if(lkindm(l).eq.k) then
                !write (iunout,*) 'surf. mol.: k,l,istr= ',k,l,istr
                do i=1,nlmpgs
                  if (lpotml)  potml(l,i)=w*potml(l,i)
                  if (lprfaml) prfaml(l,i)=w*prfaml(l,i)
                  if (lprfmml) prfmml(l,i)=w*prfmml(l,i)
                  if (lprfiml) prfiml(l,i)=w*prfiml(l,i)
                  if (leotml)  eotml(l,i)=w*eotml(l,i)
                  if (lerfaml) erfaml(l,i)=w*erfaml(l,i)
                  if (lerfmml) erfmml(l,i)=w*erfmml(l,i)
                  if (lerfiml) erfiml(l,i)=w*erfiml(l,i)
                  if (lprfpml) prfpml(l,i)=w*prfpml(l,i)
                  if (lerfpml) erfpml(l,i)=w*erfpml(l,i)
                end do
              end if
            end do  ! do l
            do l=1,nioni
              if(lkindi(l).eq.k) then
                !write (iunout,*) 'surf. ion.: k,l,istr= ',k,l,istr
                do i=1,nlmpgs
                  if (lpotio)  potio(l,i)=w*potio(l,i)
                  if (lprfaio) prfaio(l,i)=w*prfaio(l,i)
                  if (lprfmio) prfmio(l,i)=w*prfmio(l,i)
                  if (lprfiio) prfiio(l,i)=w*prfiio(l,i)
                  if (leotio)  eotio(l,i)=w*eotio(l,i)
                  if (lerfaio) erfaio(l,i)=w*erfaio(l,i)
                  if (lerfmio) erfmio(l,i)=w*erfmio(l,i)
                  if (lerfiio) erfiio(l,i)=w*erfiio(l,i)
                  if (lprfpio) prfpio(l,i)=w*prfpio(l,i)
                  if (lerfpio) erfpio(l,i)=w*erfpio(l,i)
                end do
              end if
            end do  ! do l
            do i=1,nlmpgs
              if (lprfpat) prfpat(k,i)=w*prfpat(k,i)
              if (lerfpat) erfpat(k,i)=w*erfpat(k,i)
            end do  ! do i
            do l=1,nfla
              if(lkindp(l).eq.k) then
                j=iflb(l)
                do i=1,nlmpgs
                  if (lpotpl) potpl(j,i)=w*potpl(j,i)
                  if (leotpl) eotpl(j,i)=w*eotpl(j,i)
                end do
              end if
            end do  ! do l
            do l=1,nmoli
              if(lkindm(l).eq.k) then
                do i=1,nlmpgs
                  if (lprfpml) prfpml(l,i)=w*prfpml(l,i)
                  if (lerfpml) erfpml(l,i)=w*erfpml(l,i)
                end do
              end if
            end do  ! do l
            do l=1,nioni
              if(lkindi(l).eq.k) then
                do i=1,nlmpgs
                  if (lprfpio) prfpio(l,i)=w*prfpio(l,i)
                  if (lerfpio) erfpio(l,i)=w*erfpio(l,i)
                end do
              end if
            end do  ! do l
          end if  ! hlp_fleir.gt.0.
        end if
      end do
      return
      end subroutine eirene_extrab25_eirpbls

      subroutine eirene_extrab25_eirpbls_init(nm,ni,np)
      implicit none
      integer, intent(in) :: nm, ni, np

      if (allocated(lkindm)) return
      allocate(lkindm(nm))
      if(.not.allocated(lkindp)) then
        allocate(lkindp(np))
        lkindp=0
      end if
      if(.not.allocated(lkindi)) then
        allocate(lkindi(ni))
        lkindi=0
      end if
      lkindm=0

      return
      end subroutine eirene_extrab25_eirpbls_init

      subroutine eirene_extrab25_braeir_init(chemical_sputter_yield)
      implicit none
#ifdef B25_EIRENE
#ifndef DIMENSIONS_MODULE
#include <DIMENSIONS.F>
#endif
      real(dp), intent(in) :: chemical_sputter_yield(0:DEF_NLIM+DEF_NSTS)
#else
      real(dp), intent(in) :: chemical_sputter_yield(0:nlim+nsts)
#endif
#ifndef NO_B2_CHEM_SPUT
      integer i
#ifdef B25_EIRENE
      external :: eirene_masage, eirene_masj2, eirene_exit_own

      if (nlim.gt.DEF_NLIM) then
#ifndef DIMENSIONS_MODULE
        call eirene_masage('Increase size of DEF_NLIM in DIMENSIONS.F')
#else
        call eirene_masage('Increase size of DEF_NLIM in b2mod_dimensions')
#endif
        call eirene_masj2 ('DEF_NLIM, NLIM :',DEF_NLIM,nlim)
        call eirene_exit_own(1)
      end if
      if (nsts.gt.DEF_NSTS) then
#ifndef DIMENSIONS_MODULE
        call eirene_masage('Increase size of DEF_NSTS in DIMENSIONS.F')
#else
        call eirene_masage('Increase size of DEF_NSTS in b2mod_dimensions')
#endif
        call eirene_masj2 ('DEF_NSTS, NSTS :',DEF_NSTS,nsts)
        call eirene_exit_own(1)
      end if
#endif
      allocate(bchemical_sputter_yield(0:nlim+nsts))
      do i = 0,nlim+nsts
        bchemical_sputter_yield(i) = chemical_sputter_yield(i)
      end do
#endif
      return
      end subroutine eirene_extrab25_braeir_init

      subroutine eirene_extrab25_srfprvsl_init
      implicit none
      integer :: n
      n = abs(nsrfcls)
      if(n > 0) then
        allocate(msrfcls(n))
        allocate(lsrfcls(msrfclsx,n))
      endif
      end subroutine

      subroutine eirene_extrab25_srfprvsl
      !c old version : 31.01.2005 23:41
      !c new version : 26.01.2011 (s.wiesen@fz-juelich.de) - f90 module
      !c=======================================================================
      !c*** Produces data files for visualisation of surface properties
      !c*** (a technological diagnostics to check the input data).
      !c*** There may be up to NSRFCLSX different sets of selection criteria
      !c*** (different "views"). Every combination of the values of the
      !c*** selected variables found in the input data is reflected in one or
      !c*** two files, separately for the walls ("additional surfaces") and
      !c*** grid edges ("non-default standard surfaces" with recycling).
      !c*** Each file contains a header explaining the selection criteria used,
      !c*** and (R,Z) columns for every found segment of the wall or grid edge,
      !c*** which satisfies the selection criteria (e.g., wall material).
      !c*** The data are arranged to be used directly by gnuplot.
      !c*** If a directory "Eirvsp<view>" (e.g., Eirvsp2) is present, the files
      !c*** are created there, else in a directory "Eirvsp" if it exists,
      !c*** otherwise in the current directory.
      !c***
      !c** Selection criteria (mind PARMUSR and get_crt_value by modification):
      !c***
      !c***  ind   description                 variable        header    letter
      !c***
      !c***   1  wall material               : znml(xmlim)   : xmlim     :  w
      !c***   2  surface group               : isrc(isrsg)   : isrsg     :  g
      !c***   3  sputtering model            : ilspt         : ilspt     :  s
      !c***   4  reflection model            : ilref         : ilref     :  r
      !c***   5  chem. sputtering factor     : recycc        : recycc    :  c
      !c***   6  phys. sputtering factor     : recycs        : recycs    :  p
      !c***   7  surface type                : iliin         : iliin     :  y
      !c***   8  surface side                : ilside        : ilside    :  d
      !c***   9  surface absorption          : 1-recyct      : absrp     :  a
      !c***  10  transparence out            : transp(1,)    : trnspo    :  o
      !c***  11  transparence in             : transp(2,)    : trnspi    :  i
      !c***  12  surface temperature         : ewall         : ewall     :  t

      !c*** Output filenames
      !c***
      !c***  <l1><i1>-<l2><i2>...-<ln><in>.vs[wgp]    (e.g., w1g2c1.vsw)
      !c***
      !c***  Here l1, l2, ..., ln stand for the reference letter for the
      !c***  corresponding criterion, and i1, i2, ..., in stand for the index
      !c***  of the value in the corresponding value list.
      !c***  "w" in the file extension means "wall", "g" means "grid", and "p"
      !c***  means the common picture.
      !c=======================================================================
      use eirmod_PARMMOD
      use eirmod_CLGIN
      use eirmod_CADGEO
      use eirmod_CGEOM
      use eirmod_CTRCEI
      use eirmod_COMSOU
      use eirmod_CCOUPL
      implicit none
      !c*** Input data:
      !c***  nsrfcls       : actual number of views
      !c***  msrfcls(j)    : number of criteria in view j
      !c***  lsrfcls(i,j)  : i-th criterion in view j

      !c*** Internal data:
      !c***  n_spvr        : max. number of real criteria
      !c***  n_spvi        : max. number of integer criteria
      !c***  n_spvl        : max. number of logical criteria
      !c***  l_spvt(k)     : type of the criterion k (1-real, 2-int, 3-log)
      !c***  c_spva(k)     : reference letter for criterion k in filenames
      !c***  c_spvh(k)     : reference name for criterion k in headers
      !c***  j_spvr(k)     : index of criterion k in the real arrays
      !c***  j_spvi(k)     : index of criterion k in the integer arrays
      !c***  j_spvl(k)     : index of criterion k in the logical arrays
      !c***  m_spvv(k)     : actual number of values for criterion k
      !c***  n_spvv(k)     : number of decimal digits in m_spvv(k)
      !c***  r_spvv(i,l)   : list of values for real criterion l
      !c***  i_spvv(i,l)   : list of values for integer criterion l
      !c***  l_spvv(i,l)   : list of values for logical criterion l
      !c***  t_spvr        : relative tolerance for distinguishing real values

      integer, parameter :: n_spvr=7, n_spvi=5, n_spvl=0
      real(dp), parameter :: t_spvr=1.e-7_DP
      integer l_spvt(msrfclsx),j_spvr(msrfclsx),j_spvi(msrfclsx)
      integer j_spvl(msrfclsx),m_spvv(msrfclsx),i_spvv(nlimps,n_spvi)
      integer n_spvv(msrfclsx)
      real(dp) r_spvv(nlimps,n_spvr)
      !c* No logical criteria foreseen at present, hence 0:n_spvl
      logical l_spvv(nlimps,0:n_spvl)
      character c_spva(msrfclsx),c_spvh(msrfclsx)*6
      data l_spvt /1,2,2,2,1,1,2,2,1,1,1,1/
      data c_spvh / &
                'xmlim ','isrsg ','ilspt ','ilref ','recycc','recycs', &
                'iliin ','ilside','absrp ','trnspo','trnspi','ewall '/
      data c_spva / &
                  'w'   ,  'g'   ,  's'   ,  'r'   ,  'c'   ,  'p'   , &
                  'y'   ,  'd'   ,  'a'   ,  'o'   ,  'i'   ,  't'   /
      integer, parameter :: n_luex=2, lu00=34, lu1=99

      integer luex(n_luex)
      data luex / 37, 44 /
      character hlp_fn*(msrfclsx*5+4),hlp_fr*72,hlp_ff*72
      character  hlp_cr(msrfclsx),hlp_hd(msrfclsx)*6,hlp_c,hlp_dir*6,hlp_frp*16
      character  hlp_vl(msrfclsx)*12,hlp_fp*(msrfclsx*5+14)
      data hlp_dir/'Eirvsp'/

      real(dp) :: hlp_ra(nlimps)
      integer :: hlp_ia(nlimps),ind(nlimps),hlp_k(msrfclsx)
      logical :: ex,hlp_l,hlp_d,hlp_i,hlp_v,hlp_do,l1,hlp_exit
      logical :: hlp_la(0:n_spvl)
      integer :: j,k,l,i,ic,kr,ki,iv,is,mc,ie,iw,i1,lun1,lun,ivv,istr
      integer :: j1,ii,i2,j2
      real(dp) :: x,r1
      external :: eirene_exit_own
      !c=======================================================================

      if(nsrfcls.eq.0) return
      hlp_exit=nsrfcls.lt.0
      nsrfcls=abs(nsrfcls)

      !c*** 1. Consistency checks and initialisation

      ex=.false.
      if(msrfclsx.gt.nlimps) then
        write(iunout,*) 'srfprvsl error: nlimps < msrfclsx',nlimps,msrfclsx
        ex=.true.
      end if
      if(n_spvr+n_spvi+n_spvl.ne.msrfclsx) then
        write(iunout,*) 'srfprvsl error: inconsistent dimensions'
        ex=.true.
      end if
      j=0
      k=0
      l=0
      do i=1,msrfclsx
        if(l_spvt(i).eq.1) then
          j=j+1
        else if(l_spvt(i).eq.2) then
          k=k+1
        else if(l_spvt(i).eq.3) then
          l=l+1
        else if(l_spvt(i).ne.0) then
          write(iunout,*) 'srfprvsl error: unexpected value in l_spvt. ', &
                                              'i,l_spvt(i)=',i,l_spvt(i)
          ex=.true.
        end if
      end do
      if(j.gt.n_spvr) then
        write(iunout,*) 'srfprvsl error: too many real criteria ',j,n_spvr
        ex=.true.
      end if
      if(k.gt.n_spvi) then
        write(iunout,*) 'srfprvsl error: too many integer criteria ',k,n_spvi
        ex=.true.
      end if
      if(l.gt.n_spvl) then
        write(iunout,*) 'srfprvsl error: too many logical criteria ',l,n_spvl
        ex=.true.
      end if

      do ic=1,msrfclsx
        j_spvr(ic)=0
        j_spvi(ic)=0
        j_spvl(ic)=0
      end do
      kr=0
      ki=0
      do ic=1,msrfclsx
        if(l_spvt(ic).eq.1) then
          kr=kr+1
          if(kr.gt.n_spvr) then
            write(iunout,*) 'srfprvsl: not enough space for real criteria'
            ex=.true.
          end if
          if(.not.ex) j_spvr(ic)=kr
        else if(l_spvt(ic).eq.2) then
          ki=ki+1
          if(ki.gt.n_spvi) then
            write(iunout,*) 'srfprvsl: not enough space for integer criteria'
            ex=.true.
          end if
          if(.not.ex) j_spvi(ic)=ki
        else if(l_spvt(ic).eq.3) then
          ki=ki+1
          if(ki.gt.n_spvi) then
            write(iunout,*) 'srfprvsl: not enough space for logical criteria'
            ex=.true.
          end if
          if(.not.ex) j_spvi(ic)=ki
        end if
      end do

      if(ex) then
        write(iunout,*) 'srfprvsl: inconsistent internal parameters found'
        return
      end if

      do iv=2,nsrfcls
        ex=.true.
        do j=1,iv-1
          if(ex) then
            if(msrfcls(j).eq.msrfcls(iv)) then
              hlp_l=.true.
              do i=1,msrfcls(iv)
                hlp_i=.false.
                do k=1,msrfcls(j)
                  hlp_i=hlp_i .or. lsrfcls(k,j).eq.lsrfcls(i,iv)
                end do
                hlp_l=hlp_l .and. hlp_i
              end do
              if(hlp_l) then
                write(iunout,*) 'srfprvsl: identical views ',j,iv
                msrfcls(iv)=0
                ex=.false.
              end if
            end if
          end if
        end do
      end do

      !c*** 2. Producing the values lists

      do ic=1,msrfclsx
        m_spvv(ic)=0
      end do
      do iv=1,nsrfcls
        do l=1,msrfcls(iv)
          ic=lsrfcls(l,iv)
          do is=1,nlimi+ntargi
            mc=1
            ii=is-nlimi
            if(ii.gt.0) then
              mc=nsrfsi(ii)
            else if(igjum0(is).ne.0) then
              mc=0
            end if
            do ie=1,mc
              hlp_do=.true.
              iw=0
              do while (hlp_do)
                j=is
                iw=iw+1
                hlp_do=ii.gt.0
                if(hlp_do) j=ind_srf_str(ii,ie,iw)
                if(j.gt.0) then
                  ex=get_crt_value(ic,j,i1,r1,l1,1).eq.l_spvt(ic)
                  if(.not.ex) then
                    write(iunout,*) 'srfprvsl: type inconsistency with ', &
                                                 'get_crt_value. ic=',ic
                    return
                  end if
                  do i=1,m_spvv(ic)
                    if(ex) then
                      if(l_spvt(ic).eq.1) then
                        if(t_spvr*(abs(r_spvv(i,j_spvr(ic)))+abs(r1)) &
                           .ge. abs(r_spvv(i,j_spvr(ic))-r1)) ex=.false.
                        else if(l_spvt(ic).eq.2) then
                          if(i_spvv(i,j_spvi(ic)) .eq. i1) ex=.false.
                        else if(l_spvt(ic).eq.3) then
                          ex= ex .and. (l_spvv(i,j_spvl(ic)).eqv.l1)
                        end if
                    end if
                  end do
                  if(ex .or. m_spvv(ic).eq.0) then
                    m_spvv(ic)=m_spvv(ic)+1
                    if(l_spvt(ic).eq.1) then
                      r_spvv(m_spvv(ic),j_spvr(ic))=r1
                    else if(l_spvt(ic).eq.2) then
                      i_spvv(m_spvv(ic),j_spvi(ic))=i1
                    else if(l_spvt(ic).eq.3) then
                      l_spvv(m_spvv(ic),j_spvl(ic))=l1
                    end if
                  end if
                else
                  hlp_do=.false.
                end if
              end do
            end do
          end do
        end do
      end do

      !c*** 3. Sorting the lists

      do ic=1,msrfclsx
        if(m_spvv(ic).gt.0) then
          if(l_spvt(ic).eq.1) then
            call sortrd(r_spvv(1,j_spvr(ic)),ind,m_spvv(ic))
            do i=1,m_spvv(ic)
              hlp_ra(i)=r_spvv(ind(i),j_spvr(ic))
            end do
            do i=1,m_spvv(ic)
              r_spvv(i,j_spvr(ic))=hlp_ra(i)
            end do
          else if(l_spvt(ic).eq.2) then
            call sortia(i_spvv(1,j_spvi(ic)),ind,m_spvv(ic))
            do i=1,m_spvv(ic)
              hlp_ia(i)=i_spvv(ind(i),j_spvi(ic))
            end do
            do i=1,m_spvv(ic)
              i_spvv(i,j_spvi(ic))=hlp_ia(i)
            end do
          else if(l_spvt(ic).eq.3 .and. m_spvv(ic).gt.1) then
            l_spvv(1,j_spvl(ic))=.true.
            l_spvv(2,j_spvl(ic))=.false.
          end if
        end if
      end do

      !c*** 4. Preparing parameters for the filename generation

      do i=1,msrfclsx
        n_spvv(i)=0
      end do
      do i=1,msrfclsx
        if(m_spvv(i).gt.0) then
          x=m_spvv(i)
          x=log10(x+0.1)
          n_spvv(i)=int(x)+1
        end if
      end do

      !c*** Get free unit numbers to connect the output files and
      !c*** prepare the format specification to generate the file names

      lun1=-1
      lun=igivelun(lu00,lu1, luex,n_luex)
      if(lun.gt.0) lun1=igivelun(lun+1,lu1, luex,n_luex)
      if(lun.lt.0 .or. lun1.lt.0) then
        write(iunout,*) 'srfprvsl: unable to find a free LUN for the output', &
                                                                lun,lun1
        return
      end if

      !c*** 5. Producing the output files

      hlp_i=.true.  ! means "try to open files in the subdirectory"

      if(nsrfcls.lt.10) then
        ivv=1
      else if(nsrfcls.lt.100) then
        ivv=2
      else if(nsrfcls.lt.1000) then
        ivv=3
      else
        write(iunout,*) 'srfprvsl: too many views - no separate directories'
        ivv=0
      end if
      do iv=1,nsrfcls
        if(msrfcls(iv).gt.0) then
          hlp_v=ivv.gt.0 ! means "try a separate directory for the view"

          !c*** Prepare the filenames and headers

          do l=1,msrfcls(iv)
            hlp_cr(l)=c_spva(lsrfcls(l,iv))
            hlp_hd(l)=c_spvh(lsrfcls(l,iv))
            hlp_k(l)=n_spvv(lsrfcls(l,iv))
          end do

          hlp_fr='(a,i2,a)'
          write(hlp_ff,hlp_fr) '(''(''',msrfcls(iv), &
                               '(''1h-,a,i'',i1.1,''.'',i1.1,'',''),''2a''))'

          write(hlp_fr,hlp_ff) (hlp_k(i),hlp_k(i),i=1,msrfcls(iv))

          write(hlp_ff,'(a,i2,a)') &
                       '(''##'',',msrfcls(iv),'(2x,a,'' = '',a))'

          !c*** Check for the value combinations

          do i = 1, nlimps
            hlp_ra(i)=0
            hlp_ia(i)=0
          end do
          do i = 0, n_spvl
            hlp_la(i)=.false.
          end do
          ex=set_crt_values(0,ind,hlp_ra,hlp_ia,hlp_la,r_spvv, &
                                          i_spvv,l_spvv,j_spvr,j_spvi, &
                                               j_spvl,l_spvt,m_spvv)
          do while (set_crt_values(iv,ind,hlp_ra,hlp_ia,hlp_la,r_spvv, &
                                          i_spvv,l_spvv,j_spvr,j_spvi, &
                                               j_spvl,l_spvt,m_spvv))
            hlp_l=.true.  ! means "[gw] file is not yet opened"
            hlp_d=.true.  ! means "p file is not yet opened"
            hlp_c='w'
            do is=1,nlimi+ntargi
              mc=1
              ii=is-nlimi
              if(ii.gt.0) then
                mc=nsrfsi(ii)
              else if(igjum0(is).ne.0) then
                mc=0
              end if
              do ie=1,mc
                hlp_do=.true.
                iw=0
                do while (hlp_do)
                  ex=.true.
                  j=is
                  istr=0
                  iw=iw+1
                  hlp_do=ii.gt.0
                  if(hlp_do) then
                    j=ind_srf_str(ii,ie,iw)
                    istr=j-nlim
                  end if
                  if(j.gt.0) then
                    do l=1,msrfcls(iv)
                      ic=lsrfcls(l,iv)
                      if(get_crt_value(ic,j,i1,r1,l1,1).ne. &
                                                     l_spvt(ic)) then
                        write(iunout,*) 'srfprvsl: type inconsistency with ', &
                                                 'get_crt_value. ic=',ic
                        return
                      end if
                      if(ex) then
                        if(l_spvt(ic).eq.1) then
                          if(t_spvr*(abs(hlp_ra(j_spvr(ic)))+abs(r1)) &
                             .lt. abs(hlp_ra(j_spvr(ic))-r1)) ex=.false.
                        else if(l_spvt(ic).eq.2) then
                          if(hlp_ia(j_spvi(ic)) .ne. i1) ex=.false.
                        else if(l_spvt(ic).eq.3) then
                          !if(hlp_la(j_spvl(ic)) .ne. l1) ex=.false.
                          ex= ex .and. (hlp_la(j_spvl(ic)).neqv.l1)
                        end if
                      end if
                    end do
                    if(ex) then

                      !c*** Combination found

                      if(hlp_l) then
                        !c***    Open the [gw] file
                        inquire(lun,opened=ex)
                        if(ex) close(lun)
                        write(hlp_fn,hlp_fr) &
                          (hlp_cr(i),ind(i),i=1,msrfcls(iv)),'.vs',hlp_c
                        hlp_fn=hlp_fn(2:len(hlp_fn))

  100                   continue
                        if(hlp_v) then
                          write(hlp_frp,'(a4,i1,a1,i1,a7)') &
                                        '(a,i',ivv,'.',ivv,',a)'
                          write(hlp_fp,hlp_frp) hlp_dir,iv,'/'//hlp_fn
                          open(lun,file=hlp_fp,err=120)
                        else if(hlp_i) then
                          hlp_fp=hlp_dir//'/'//hlp_fn
                          open(lun,file=hlp_fp,err=110)
                        else
                          hlp_fp=hlp_fn
                          open(lun,file=hlp_fp,err=9000)
                        end if
                        go to 150
  110                   hlp_i=.false.
                        go to 100
  120                   hlp_v=.false.
                        go to 100
  150                   continue
                        inquire(lun,opened=ex)
                        write (iunout,*) '== file: ',lun,ex,'  ',hlp_fp

                        !c***    Write the header

                        do l=1,msrfcls(iv)
                          ic=lsrfcls(l,iv)
                          if(l_spvt(ic).eq.1) then
                            write(hlp_vl(l),'(1p,g12.5)') &
                                               r_spvv(ind(l),j_spvr(ic))
                          else if(l_spvt(ic).eq.2) then
                            write(hlp_vl(l),'(i6)') &
                                               i_spvv(ind(l),j_spvi(ic))
                          else if(l_spvt(ic).eq.3) then
                            write(hlp_vl(l),'(l1)') &
                                               l_spvv(ind(l),j_spvl(ic))
                          end if
                        end do
                        write(lun,hlp_ff) &
                                   (hlp_hd(i),hlp_vl(i),i=1,msrfcls(iv))
                        if(hlp_c.eq.'w') then
                          write(lun,'(3a)') '# Additional surfaces'
                          write(lun,'(a10,8x,a1,14x,a1)') '# ieir idg','R','Z'
                        else
                          write(lun,'(3a)') '# Non-default ', &
                                                     'standard surfaces'
                          write(lun,'(a,6x,a1,14x,a1)') &
                                    '# stratum sub-stratum non-def-surf','R','Z'
                        end if
                        hlp_l=.false.
                      end if

                      if(hlp_d) then
                        !c***    Open the p file
                        inquire(lun1,opened=ex)
                        if(ex) close(lun1)
                        write(hlp_fn,hlp_fr)  &
                            (hlp_cr(i),ind(i),i=1,msrfcls(iv)),'.vs','p'
                        hlp_fn=hlp_fn(2:len(hlp_fn))

  200                   continue
                        if(hlp_v) then
                          write(hlp_frp,'(a4,i1,a1,i1,a7)') &
                                        '(a,i',ivv,'.',ivv,',a)'
                          write(hlp_fp,hlp_frp) hlp_dir,iv,'/'//hlp_fn
                          open(lun1,file=hlp_fp,err=220)
                        else if(hlp_i) then
                          hlp_fp=hlp_dir//'/'//hlp_fn
                          open(lun1,file=hlp_fp,err=210)
                        else
                          hlp_fp=hlp_fn
                          open(lun1,file=hlp_fp,err=9000)
                        end if
                        go to 250
  210                   hlp_i=.false.
                        go to 200
  220                   hlp_v=.false.
                        go to 200
  250                   continue
                        inquire(lun1,opened=ex)
                        write (iunout,*) '== file: ',lun1,ex,'  ',hlp_fp

                        !c***    Write the header

                        do l=1,msrfcls(iv)
                          ic=lsrfcls(l,iv)
                          if(l_spvt(ic).eq.1) then
                            write(hlp_vl(l),'(1p,g12.5)') &
                                               r_spvv(ind(l),j_spvr(ic))
                          else if(l_spvt(ic).eq.2) then
                            write(hlp_vl(l),'(i6)') &
                                               i_spvv(ind(l),j_spvi(ic))
                          else if(l_spvt(ic).eq.3) then
                            write(hlp_vl(l),'(l1)') &
                                               l_spvv(ind(l),j_spvl(ic))
                          end if
                        end do
                        write(lun1,hlp_ff) &
                                   (hlp_hd(i),hlp_vl(i),i=1,msrfcls(iv))
                        write(lun1,'(3a)') '# All surfaces together'
                        write(lun1,'(a,6x,a1,14x,a1)') &
                                   '#   ieir       idg                ','R','Z'
                        write(lun1,'(a,6x,a1,14x,a1)') &
                                   '# stratum sub-stratum non-def-surf','R','Z'
                        hlp_d=.false.
                      end if

                      !c***      Write the values
                      if(ii.le.0) then
                        !c***        Wall
                        i=index_in_dg(is)
                        write(lun,'(2i5,1p,2e15.5)') is,i, &
                                             0.01*p1(1,is),0.01*p1(2,is)
                        write(lun,'(2i5,1p,2e15.5)') is,i, &
                                             0.01*p2(1,is),0.01*p2(2,is)
                        write(lun,*)
                        write(lun1,'(i6,2i12,4x,1p,2e15.5)') is,i,0, &
                                             0.01*p1(1,is),0.01*p1(2,is)
                        write(lun1,'(i6,2i12,4x,1p,2e15.5)') is,i,0, &
                                             0.01*p2(1,is),0.01*p2(2,is)
                        write(lun1,*)
                      else
                        !c***        Grid

                        i1=0
                        if(indim(ie,ii).eq.1) then             !c radial:
                          i1=insor(ie,ii)
                          j1=ingrda(ie,ii,2)+iw-1
                          i2=i1
                          j2=j1+1
                        else if(indim(ie,ii).eq.2) then        !c poloidal:
                          i1=ingrda(ie,ii,1)+iw-1
                          j1=insor(ie,ii)
                          i2=i1+1
                          j2=j1
                        end if
                        write(lun,'(i6,2i12,4x,1p,2e15.5)') ii,ie,istr, &
                                       0.01*xpol(i1,j1),0.01*ypol(i1,j1)
                        write(lun,'(i6,2i12,4x,1p,2e15.5)') ii,ie,istr, &
                                       0.01*xpol(i2,j2),0.01*ypol(i2,j2)
                        write(lun,*)
                        write(lun1,'(i6,2i12,4x,1p,2e15.5)') ii,ie,istr, &
                                       0.01*xpol(i1,j1),0.01*ypol(i1,j1)
                        write(lun1,'(i6,2i12,4x,1p,2e15.5)') ii,ie,istr, &
                                       0.01*xpol(i2,j2),0.01*ypol(i2,j2)
                        write(lun1,*)
                      end if
                    end if
                  else
                    hlp_do=.false.
                  end if
                end do
              end do
              if(ii.eq.0) then
                hlp_l=.true.
                hlp_c='g'
              end if
            end do
            inquire(lun,opened=ex)
            if(ex) close(lun)
            hlp_l=.true.
            inquire(lun1,opened=ex)
            if(ex) close(lun1)
            hlp_d=.true.
          end do
        end if
      end do

      if(hlp_exit) call eirene_exit_own(1)
      return
      !c=======================================================================
 9000 write(iunout,*) 'srfprvsl: error opening the file ',hlp_fn
      return
      !c=======================================================================
      return

      contains

      subroutine sortrd(a,ind,n)
      !c*** sorting double precision array in ascending order
      implicit none
      real(dp) :: a(*)
      integer :: ind(*)
      integer :: n,i,j,m,ll,k
      real(dp) :: r,u
      do i=1,n
        ind(i)=i
      enddo
      if(n.le.1) return
      do k=2,n
        m=k-1
        ll=ind(k-1)
        r=a(ll)
        do j=k,n
          l=ind(j)
          u=a(l)
          if(u.lt.r) then
            ll=l
            r=u
            m=j
          end if
        enddo
        if(m.ne.k-1) then
          i=ind(k-1)
          ind(k-1)=ll
          ind(m)=i
        end if
      enddo
      end subroutine

      subroutine sortia(a,ind,n)
      !c*** sorting integer array in ascending order
      implicit none
      integer :: a(*),ind(*)
      integer :: n
      integer :: r,u,i,m,ll,j,k
      !c
      do i=1,n
        ind(i)=i
      enddo
      if(n.le.1) return
      do k=2,n
        m=k-1
        ll=ind(k-1)
        r=a(ll)
        do j=k,n
          l=ind(j)
          u=a(l)
          if(u.lt.r) then
            ll=l
            r=u
            m=j
          end if
        enddo
        if(m.ne.k-1) then
          i=ind(k-1)
          ind(k-1)=ll
          ind(m)=i
        end if
      enddo
      end subroutine

      integer function igivelun(l0,l1,lex,nl)
      implicit none
      !c*** This function returns a free LUN for I/O in the range from l0 to l1
      !c*** lex contains reserved numbers (nl) to be excluded from search
      integer :: lex(*),nl,l0,l1
      logical :: op
      integer :: i,j
      do i=l0,l1
        inquire(unit=i,opened=op)
        if(.not.op) then
          do j=1,nl
            op= op .or. i.eq.lex(j)
          end do
          if(.not.op) then
            igivelun=i
            return
          end if
        end if
      end do
      igivelun=-1
      end function

      integer function index_in_dg(is)
      !c*** Returns "element number" from the DG model, corresponding to
      !c*** the additional surface is, or 0 if the correspondence in not found
      use eirmod_PARMMOD
      use eirmod_CADGEO
      use eirmod_CTEXT
      implicit none
      character*72 :: t
      integer :: i,is,l

      index_in_dg=0
      if(is.gt.nlimi .or. is.le.0) return
      l=index(txtsfl(is),':')
      if(l.eq.0) return
      t=txtsfl(is)(l+1:)//'/'
      i=0
      read(t,*,err=10) i
      index_in_dg=i
 10   return
      !c=======================================================================
      end function

      logical function set_crt_values(iw,ind,ra,ia,la,rv,iv,lv, &
                                      jr,ji,jl,lt,nv)
      !c*** Sets the next set of values for a view iw (iw=0 - reset counters)
      !c*** Input:
      !c***  nv  number of values for each criterion
      !c***  lt  criterion type (1-real, 2-int, 3-log)
      !c***  jr,ji,jl  criterion indices (real, integer, logical)
      !c***  rv,iv,lv value sets for each criterion (real, integer, logical)
      use eirmod_PARMMOD
      use eirmod_CTRCEI
      use eirmod_CLGIN
      implicit none
      integer :: iw
      real(dp) :: ra(nlimps),rv(nlimps,*)
      integer :: ind(nlimps),ia(nlimps),iv(nlimps,*),lt(*), &
                 nv(*),jr(*),ji(*),jl(*)
      logical :: la(0:*),lv(nlimps,0:*)
      integer :: i,j
      logical :: ex

      set_crt_values=.true.
      if(iw.eq.0) then

        !c*** Reset the counters

        do i=2,msrfclsx
          ind(i)=1
        end do
        ind(1)=0
        return
      end if
      ex=.true.
      do i=1,msrfcls(iw)
        j=lsrfcls(i,iw)
        if(ex) then
          ind(i)=ind(i)+1
          ex=ind(i).gt.nv(j)
          if(ex) ind(i)=1
        end if
      end do
      if(ex) then
        set_crt_values=.false.
        return
      end if
      do j=1,msrfcls(iw)
        i=lsrfcls(j,iw)
        if(lt(i).eq.1) then
          ra(jr(i))=rv(ind(j),jr(i))
        else if(lt(i).eq.2) then
          ia(ji(i))=iv(ind(j),ji(i))
        else if(lt(i).eq.3) then
          la(jl(i))=lv(ind(j),jl(i))
        else
          set_crt_values=.false.
        end if
      end do
      !c=======================================================================
      end function

      integer function ind_srf_str(is,ic,il)
      !c*** This function returns index of a "non-default standard surface"
      !c*** corresponding to a segment IL of a surface IC of a recycling
      !c*** stratum IS, or 0 if no correspondence found (e.g., segment outside
      !c*** the stratum)
      !c*** Only 2D geometry (no toroidal surfaces) is considered
      use eirmod_PARMMOD
      use eirmod_CLGIN
      use eirmod_COMSOU
      use eirmod_CCOUPL
      implicit none
      integer :: j,is,il,k,ic
      !c=======================================================================
      j=0
      ind_srf_str=0
      if(il.gt.0 .and. is.gt.0 .and. is.le.ntargi .and. &
                                 ic.gt.0 .and. ic.le.nsrfsi(is)) then
        if(indim(ic,is).eq.1) then !                  radial surface:

          k=ingrda(ic,is,2)+il-1
          if(k.lt.ingrde(ic,is,2)) &
                                j=inmp1i(insor(ic,is),k,ingrda(ic,is,3))

        else if(indim(ic,is).eq.2) then !          poloidal surface:

          k=ingrda(ic,is,1)+il-1
          if(k.lt.ingrde(ic,is,1)) &
                                j=inmp2i(k,insor(ic,is),ingrda(ic,is,3))

        end if
      end if
      if(j.gt.0) ind_srf_str=nlim+j
      !c=======================================================================
      end function

      integer function get_crt_value(ic,j,i1,r1,l1,isp)
      !c*** This function returns the type of the criterion and sets i1, r1,
      !c*** or l1 to the parameter value for integer, real, or logical
      !c*** criterion, respectively.
      !c*** The input parameters are criterion index ic and surface index j.
      !csw added isp for ISPEZ and ISRC (CHECK)
      use eirmod_PARMMOD
      use eirmod_CLGIN
      use eirmod_CTRCEI
      implicit none
      integer :: ic,j,i1,isp
      real(dp) :: r1
      logical :: l1
      !c=======================================================================
      if(ic.gt.mgwtiesx .or. ic.le.0) then
        get_crt_value=-1
        return
      end if
      get_crt_value=0
      i1=0
      r1=0.
      l1=.false.
      !c-----------------------------------------------------------------------
      select case (ic)
      case (1)
        r1=100.*znml(j)+zncl(j)
        get_crt_value=1
      case (2)
        i1=isrc(isp,j)
        get_crt_value=2
      case (3)
        i1=ilspt(j)
        get_crt_value=2
      case (4)
        i1=ilref(isp,j)
        get_crt_value=2
      case (5)
        r1=recycc(isp,j)
        get_crt_value=1
      case (6)
        r1=recycs(isp,j)
        get_crt_value=1
      case (7)
        i1=iliin(j)
        get_crt_value=2
      case (8)
        i1=ilside(j)
        get_crt_value=2
      case (9)
        r1=1.-recyct(isp,j)
        get_crt_value=1
      case (10)
        r1=transp(isp,1,j)
        get_crt_value=1
      case (11)
        r1=transp(isp,2,j)
        get_crt_value=1
      case (12)
        r1=ewall(j)
        get_crt_value=1
      end select
      return
      !c=======================================================================
      end function

      end subroutine

      subroutine eirene_extrab25_iniusr_init(n_spcsrf,l_spcsrf,&
                i_spcsrf,&
                j_spcsrf,sps_sgrp,sps_absr,sps_trno,sps_trni,&
                sps_mtri,sps_tmpr,sps_spph,sps_spch,&
                sps_mtrl,sps_id)
      implicit none
      integer, intent(in) :: n_spcsrf
#ifdef B25_EIRENE
#ifndef DIMENSIONS_MODULE
#include <DIMENSIONS.F>
#endif
      integer, intent(in) :: l_spcsrf(DEF_NLIM+DEF_NSTS), &
            i_spcsrf(DEF_NSPCSRFG), j_spcsrf(DEF_NSPCSRFG), &
            sps_sgrp(DEF_NSPCSRFG)
      real(dp), intent(in) :: sps_absr(DEF_NSPCSRFG), &
            sps_trno(DEF_NSPCSRFG), sps_trni(DEF_NSPCSRFG), &
            sps_mtri(DEF_NSPCSRFG), sps_tmpr(DEF_NSPCSRFG), &
            sps_spph(DEF_NSPCSRFG), sps_spch(DEF_NSPCSRFG)
      character*8, intent(in) :: sps_mtrl(DEF_NSPCSRFG), &
                                 sps_id(DEF_NSPCSRFG)
#else
      integer, intent(in) :: l_spcsrf(nlim+nsts), &
            i_spcsrf(*), j_spcsrf(*), sps_sgrp(*)
      real(dp), intent(in) :: sps_absr(*), sps_trno(*), sps_trni(*), &
               sps_mtri(*), sps_tmpr(*), sps_spph(*), sps_spch(*)
      character*8, intent(in) :: sps_mtrl(*), sps_id(*)
#endif

      bn_spcsrf=n_spcsrf
      allocate(bl_spcsrf(nlimps))
      allocate(bi_spcsrf(n_spcsrf))
      allocate(bj_spcsrf(n_spcsrf))
      allocate(bsps_sgrp(n_spcsrf))
      allocate(bsps_absr(n_spcsrf))
      allocate(bsps_trno(n_spcsrf))
      allocate(bsps_trni(n_spcsrf))
      allocate(bsps_mtri(n_spcsrf))
      allocate(bsps_tmpr(n_spcsrf))
      allocate(bsps_spph(n_spcsrf))
      allocate(bsps_spch(n_spcsrf))
      allocate(bsps_mtrl(n_spcsrf))
      allocate(bsps_id(n_spcsrf))

      bl_spcsrf(1:nlimps) = l_spcsrf(1:nlimps)
      bi_spcsrf(1:n_spcsrf) = i_spcsrf(1:n_spcsrf)
      bj_spcsrf(1:n_spcsrf) = j_spcsrf(1:n_spcsrf)
      bsps_sgrp(1:n_spcsrf) = sps_sgrp(1:n_spcsrf)
      bsps_absr(1:n_spcsrf) = sps_absr(1:n_spcsrf)
      bsps_trno(1:n_spcsrf) = sps_trno(1:n_spcsrf)
      bsps_trni(1:n_spcsrf) = sps_trni(1:n_spcsrf)
      bsps_mtri(1:n_spcsrf) = sps_mtri(1:n_spcsrf)
      bsps_tmpr(1:n_spcsrf) = sps_tmpr(1:n_spcsrf)
      bsps_spph(1:n_spcsrf) = sps_spph(1:n_spcsrf)
      bsps_spch(1:n_spcsrf) = sps_spch(1:n_spcsrf)
      bsps_mtrl(1:n_spcsrf) = sps_mtrl(1:n_spcsrf)
      bsps_id(1:n_spcsrf) = sps_id(1:n_spcsrf)

!pb 27012016
! flag indicating if subroutine iniusr is called from B2.5
      ini_iniusr = 1
      end subroutine eirene_extrab25_iniusr_init


      subroutine eirene_extrab25_cleanup
      implicit none

      if(allocated(msrfcls)) then
        deallocate(msrfcls)
        deallocate(lsrfcls)
      endif

      if(allocated(bl_spcsrf)) then
        deallocate(bl_spcsrf)
        deallocate(bi_spcsrf)
        deallocate(bj_spcsrf)
        deallocate(bsps_sgrp)
        deallocate(bsps_absr)
        deallocate(bsps_trno)
        deallocate(bsps_trni)
        deallocate(bsps_mtri)
        deallocate(bsps_tmpr)
        deallocate(bsps_spph)
        deallocate(bsps_spch)
        deallocate(bsps_mtrl)
        deallocate(bsps_id)
      endif

      if(allocated(flux_save)) then
        deallocate(flux_save)
      endif

      if(allocated(plnxtri)) then
        deallocate(plnxtri, plnytri, pplnxtri, pplnytri)
      endif

      if(allocated(bchemical_sputter_yield)) then
        deallocate(bchemical_sputter_yield)
      endif

#ifdef B25_EIRENE
      if(allocated(rcpos_eir)) then
        deallocate(rcpos_eir)
        deallocate(rcbeg_eir)
        deallocate(rcend_eir)
        deallocate(rcprt_eir)
        deallocate(rcspi_eir)
        deallocate(rcspe_eir)
        deallocate(rcchr_eir)
      endif
#endif
      return
      end subroutine eirene_extrab25_cleanup


      subroutine read_title(kard)
      integer kard
      character*80 zeile

!xpb This is an optional title line to an eirene data field
      read(kard,'(a80)',end=12) zeile
      if (index(zeile,'*eirene data field').eq.0) then
        backspace(kard)
      else
        if (trcint) write(iunout,'(a)') 'Reading '//trim(zeile)
      end if
   12 continue
      return
      end subroutine read_title


      subroutine write_title(kard,titel,nsize)
      integer kard, nsize
      character*(*) titel

      write(kard,'(3a,i6)') &
       & '*eirene data field ',trim(titel),' with size ',nsize
      return
      end subroutine write_title


      subroutine eirene_extrab25_emissivity
      implicit none
      integer :: istr, i, j, iadv, icell
      external :: eirene_emissivity, eirene_rstrt

      istr = 0

      IF (IESTR.EQ.ISTR) THEN
!  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=ISTR
        CALL EIRENE_RSTRT(ISTR,NSTRAI,NESTM1,NESTM2,NADSPC, &
                   ESTIMV,ESTIMS,ESTIML, &
                   NSDVI1,SDVI1,NSDVI2,SDVI2, &
                   NSDVC1,SIGMAC,NSDVC2,SGMCS, &
                   NSIGI_SPC,TRCFLE)
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTR.EQ.0) THEN
        IESTR=ISTR
        CALL EIRENE_RSTRT(ISTR,NSTRAI,NESTM1,NESTM2,NADSPC, &
                   ESTIMV,ESTIMS,ESTIML, &
                   NSDVI1,SDVI1,NSDVI2,SDVI2, &
                   NSDVC1,SIGMAC,NSDVC2,SGMCS, &
                   NSIGI_SPC,TRCFLE)
      ELSE
        WRITE (IUNOUT,*) 'ERROR IN EXTRAB25_EMISSIVITY: ' // &
                         'DATA FOR STRATUM ISTRA= ', ISTR
        WRITE (IUNOUT,*) &
          'ARE NOT AVAILABLE. EXTRAB25_EMISSIVITY ABANDONED'
        RETURN
      ENDIF

      emiss(:,1,1) = 0._dp
      emissmol(:,1,1) = 0._dp
      do i = 1, num_lines
        write (iunout,*) 'in EXTRAB25_EMISSIVITY, line no. = ',i
        if (mod_addv == 0) then
! storage saving mode:
! ADDV is overwritten when a new line comes, within a run.
! Thus re-calculate the new emissivity profile on ADDV now
          call eirene_emissivity(istr, i, i, 0)
!       else
! Sufficiently large storage on ADDV additional tally array,
! for all lines and components. No need to reset ADDV tallies.
        end if

!dr run over components
        do j = 1, emis_lines(i)%num_compo
!dr  iadv: tally number on ADDV
          iadv = emis_lines(i)%compo(j)%iadv

!pb  This is dangerous! It is implicitly assumed that the default
!pb  emissivity model is used.
!pb  In case of a user-specific model specified in the Eirene input
!pb  this might produce rubbish.
          if (i.eq.1.and.istr.eq.0) then ! Ba-alpha emissivity for fort.44
            do icell = 1, nsbox_tal
              if (icell.gt.ncvp) cycle   !nh Cell in void region
              if (j.ge.1 .and. j.le.2) then ! Atomic components
                emiss(icell,1,1)=emiss(icell,1,1)+addv(iadv,icell)*1.0d6
              else if (j.ge.3 .and. j.le.6) then ! Molecular components
                emissmol(icell,1,1)=emissmol(icell,1,1)+addv(iadv,icell)*1.0d6
              end if
            end do
          end if
        end do
      end do

      end subroutine eirene_extrab25_emissivity

      end module eirmod_extrab25

!!!Local Variables:
!!! mode: f90
!!! End:
