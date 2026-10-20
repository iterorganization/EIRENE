      MODULE EIRMOD_EIRBRA

C  NEUTRAL SOURCE TERMS: SNI,SMO,SEE,SEI (EIRENE ---> BRAAMS)

      USE EIRMOD_PRECISION
      USE EIRMOD_MPI
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: EIRENE_ALLOC_EIRBRA, EIRENE_DEALLOC_EIRBRA,
     P          EIRENE_INIT_EIRBRA
csw mpi
      public :: eirene_broadcast_eirbra
      public :: eirene_mpisend_eirbra, eirene_mpirecv_eirbra
csw

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R SNI(:,:,:,:), SMO(:,:,:,:),
     R SEE(:,:,:),   SEI(:,:,:)
cdjm Jan2017
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R SNE_PAEL(:,:,:), SNE_PMEL(:,:,:), SNE_PIEL(:,:,:),
     R SNA_PAAT(:,:,:,:), SNA_PMAT(:,:,:,:),
     R SNA_PIAT(:,:,:,:), SNA_PAML(:,:,:,:),
     R SNA_PMML(:,:,:,:), SNA_PIML(:,:,:,:),
     R SNA_PAIO(:,:,:,:), SNA_PMIO(:,:,:,:),
     R SNA_PIIO(:,:,:,:),
     R SNI_PAPL(:,:,:,:), SNI_PMPL(:,:,:,:),
     R SNI_PIPL(:,:,:,:), SNI_PPPL(:,:,:,:),
     R SMO_MAPL(:,:,:,:), SMO_MMPL(:,:,:,:),
     R SMO_MIPL(:,:,:,:), SMO_MPPL(:,:,:,:),
     R SEE_EAEL(:,:,:), SEE_EMEL(:,:,:),
     R SEE_EIEL(:,:,:), SEE_EPEL(:,:,:),
     R SEI_EAPL(:,:,:), SEI_EMPL(:,:,:),
     R SEI_EIPL(:,:,:), SEI_EPPL(:,:,:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R VOLSUMN(:), VOLSUMM(:,:), VOLSUMEI(:), VOLSUMEE(:)

      INTEGER, PUBLIC, SAVE ::
     I NDXD, NDYD, NFLD, NATD, NMLD, NIOD, NSTRAD

C AK
      REAL(DP),PUBLIC,ALLOCATABLE,SAVE ::
     .                srcstrn(:)
     .               ,flxspci(:,:)
     .               ,srccrfc(:,:)
c*** srcstrn: intensities of different neutral sources (fluxt)
c*** flxspci: fluxes of plasma ions to the recycling surfaces
c*** srccrfc: source correction factors (from infcop)

C AK END

      INTEGER, PARAMETER :: LONG = SELECTED_INT_KIND(10)

      CONTAINS


      SUBROUTINE EIRENE_ALLOC_EIRBRA(NDXP, NDYP, NFL,
     .                               NATM, NMOL, NION, NSTRA)

      USE EIRMOD_PARMMOD, ONLY: IUNMEM

      INTEGER, INTENT(IN) :: NDXP, NDYP, NFL, NATM, NMOL, NION, NSTRA
      INTEGER(LONG) :: IMEM

      IF (ALLOCATED(SNI)) RETURN

      NDXD = NDXP+1
      NDYD = NDYP+1
      NATD = NATM
      NMLD = NMOL
      NIOD = NION
      NFLD = NFL
!pb  obviously SOLPS always assumes there is a time stratum
!pb  therefore add additional space in arrays
csw   NSTRAD = NSTRA
      NSTRAD = NSTRA + 1

      ALLOCATE (SNI(0:NDXD,0:NDYD,NFLD,NSTRAD))
cdjm Jan2017
      ALLOCATE (SNE_PAEL(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SNE_PMEL(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SNE_PIEL(0:NDXD,0:NDYD,NSTRAD))
cdjm Jan2017
      ALLOCATE (SNA_PAAT(0:NDXD,0:NDYD,NATD,NSTRAD))
      ALLOCATE (SNA_PMAT(0:NDXD,0:NDYD,NATD,NSTRAD))
      ALLOCATE (SNA_PIAT(0:NDXD,0:NDYD,NATD,NSTRAD))
cdjm Jan2017
      ALLOCATE (SNA_PAML(0:NDXD,0:NDYD,NMLD,NSTRAD))
      ALLOCATE (SNA_PMML(0:NDXD,0:NDYD,NMLD,NSTRAD))
      ALLOCATE (SNA_PIML(0:NDXD,0:NDYD,NMLD,NSTRAD))
cdjm Jan2017
      ALLOCATE (SNA_PAIO(0:NDXD,0:NDYD,NIOD,NSTRAD))
      ALLOCATE (SNA_PMIO(0:NDXD,0:NDYD,NIOD,NSTRAD))
      ALLOCATE (SNA_PIIO(0:NDXD,0:NDYD,NIOD,NSTRAD))
cdjm Jan2017
      ALLOCATE (SNI_PAPL(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SNI_PMPL(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SNI_PIPL(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SNI_PPPL(0:NDXD,0:NDYD,NFLD,NSTRAD))

!dr  may 2022: vectorial momentum sources
      ALLOCATE (SMO(0:NDXD,0:NDYD,3*NFLD,NSTRAD))
cdjm Jan2017
      ALLOCATE (SMO_MAPL(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SMO_MMPL(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SMO_MIPL(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SMO_MPPL(0:NDXD,0:NDYD,NFLD,NSTRAD))

      ALLOCATE (SEE(0:NDXD,0:NDYD,NSTRAD))
cdjm Jan2017
      ALLOCATE (SEE_EAEL(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SEE_EMEL(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SEE_EIEL(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SEE_EPEL(0:NDXD,0:NDYD,NSTRAD))

      ALLOCATE (SEI(0:NDXD,0:NDYD,NSTRAD))
cdjm Jan2017
      ALLOCATE (SEI_EAPL(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SEI_EMPL(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SEI_EIPL(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SEI_EPPL(0:NDXD,0:NDYD,NSTRAD))

      ALLOCATE (VOLSUMN(NSTRAD))
      ALLOCATE (VOLSUMM(0:2,NSTRAD))
      ALLOCATE (VOLSUMEI(NSTRAD))
      ALLOCATE (VOLSUMEE(NSTRAD))

      ALLOCATE (SRCSTRN(NSTRAD))
      ALLOCATE (SRCCRFC(NATD,NSTRAD))
      ALLOCATE (FLXSPCI(NFLD,NSTRAD))

      IMEM = INT((INT(NDXD+1,LONG)*INT(NDYD+1,LONG)*INT(NSTRAD,LONG)*
     .            INT(12*NFLD+3*(NATD+NMLD+NIOD)+13,LONG)+
     .            INT(NATD+NFLD+7,LONG)*INT(NSTRAD,LONG))*8,LONG)
      WRITE (IUNMEM,'(A,T25,I15)')
     .             ' EIRBRA ', IMEM

      CALL EIRENE_INIT_EIRBRA

      RETURN
      END SUBROUTINE EIRENE_ALLOC_EIRBRA


      SUBROUTINE EIRENE_DEALLOC_EIRBRA

      IF (.NOT.ALLOCATED(SNI)) RETURN

      DEALLOCATE (SNI)
cdjm Jan2017
      DEALLOCATE (SNE_PAEL)
      DEALLOCATE (SNE_PMEL)
      DEALLOCATE (SNE_PIEL)
cdjm Jan2017
      DEALLOCATE (SNA_PAAT)
      DEALLOCATE (SNA_PMAT)
      DEALLOCATE (SNA_PIAT)
cdjm Jan2017
      DEALLOCATE (SNA_PAML)
      DEALLOCATE (SNA_PMML)
      DEALLOCATE (SNA_PIML)
cdjm Jan2017
      DEALLOCATE (SNA_PAIO)
      DEALLOCATE (SNA_PMIO)
      DEALLOCATE (SNA_PIIO)
cdjm Jan2017
      DEALLOCATE (SNI_PAPL)
      DEALLOCATE (SNI_PMPL)
      DEALLOCATE (SNI_PIPL)
      DEALLOCATE (SNI_PPPL)

      DEALLOCATE (SMO)
cdjm Jan2017
      DEALLOCATE (SMO_MAPL)
      DEALLOCATE (SMO_MMPL)
      DEALLOCATE (SMO_MIPL)
      DEALLOCATE (SMO_MPPL)

      DEALLOCATE (SEE)
cdjm Jan2017
      DEALLOCATE (SEE_EAEL)
      DEALLOCATE (SEE_EMEL)
      DEALLOCATE (SEE_EIEL)
      DEALLOCATE (SEE_EPEL)

      DEALLOCATE (SEI)
cdjm Jan2017
      DEALLOCATE (SEI_EAPL)
      DEALLOCATE (SEI_EMPL)
      DEALLOCATE (SEI_EIPL)
      DEALLOCATE (SEI_EPPL)

      DEALLOCATE (VOLSUMN)
      DEALLOCATE (VOLSUMM)
      DEALLOCATE (VOLSUMEI)
      DEALLOCATE (VOLSUMEE)

      DEALLOCATE (SRCSTRN)
      DEALLOCATE (SRCCRFC)
      DEALLOCATE (FLXSPCI)

      RETURN
      END SUBROUTINE EIRENE_DEALLOC_EIRBRA


      SUBROUTINE EIRENE_INIT_EIRBRA

      SNI = 0._DP
cdjm Jan2017
      SNE_PAEL = 0._DP
      SNE_PMEL = 0._DP
      SNE_PIEL = 0._DP
cdjm Jan2017
      SNA_PAAT = 0._DP
      SNA_PMAT = 0._DP
      SNA_PIAT = 0._DP
cdjm Jan2017
      SNA_PAML = 0._DP
      SNA_PMML = 0._DP
      SNA_PIML = 0._DP
cdjm Jan2017
      SNA_PAIO = 0._DP
      SNA_PMIO = 0._DP
      SNA_PIIO = 0._DP
cdjm Jan2017
      SNI_PAPL = 0._DP
      SNI_PMPL = 0._DP
      SNI_PIPL = 0._DP
      SNI_PPPL = 0._DP

      SMO = 0._DP
cdjm Jan2017
      SMO_MAPL = 0._DP
      SMO_MMPL = 0._DP
      SMO_MIPL = 0._DP
      SMO_MPPL = 0._DP

      SEE = 0._DP
cdjm Jan2017
      SEE_EAEL = 0._DP
      SEE_EMEL = 0._DP
      SEE_EIEL = 0._DP
      SEE_EPEL = 0._DP

      SEI = 0._DP
cdjm Jan2017
      SEI_EAPL = 0._DP
      SEI_EMPL = 0._DP
      SEI_EIPL = 0._DP
      SEI_EPPL = 0._DP

      VOLSUMN  = 0._DP
      VOLSUMM  = 0._DP
      VOLSUMEI = 0._DP
      VOLSUMEE = 0._DP

      SRCSTRN = 0._DP
      SRCCRFC = 0._DP
      FLXSPCI = 0._DP

      RETURN
      END SUBROUTINE EIRENE_INIT_EIRBRA


csw mpi
      subroutine eirene_broadcast_eirbra

      integer :: inum, ierr
#if( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      external :: mpi_bcast
#endif

      inum = (ndxd+1)*(ndyd+1)*nfld*nstrad
      call mpi_bcast (sni,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
cdjm Jan2017
      call mpi_bcast (sni_papl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sni_pmpl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sni_pipl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sni_pppl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

      inum = (ndxd+1)*(ndyd+1)*3*nfld*nstrad
      call mpi_bcast (smo,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
cdjm Jan2017
      inum = (ndxd+1)*(ndyd+1)*nfld*nstrad
      call mpi_bcast (smo_mapl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (smo_mmpl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (smo_mipl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (smo_mppl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

cdjm Jan2017
      inum = (ndxd+1)*(ndyd+1)*natd*nstrad
      call mpi_bcast (sna_paat,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sna_pmat,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sna_piat,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

cdjm Jan2017
      inum = (ndxd+1)*(ndyd+1)*nmld*nstrad
      call mpi_bcast (sna_paml,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sna_pmml,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sna_piml,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

cdjm Jan2017
      inum = (ndxd+1)*(ndyd+1)*niod*nstrad
      call mpi_bcast (sna_paio,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sna_pmio,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sna_piio,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

      inum = (ndxd+1)*(ndyd+1)*nstrad
      call mpi_bcast (see,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
cdjm Jan2017
      call mpi_bcast (sne_pael,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sne_pmel,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sne_piel,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
cdjm Jan2017
      call mpi_bcast (see_eael,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (see_emel,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (see_eiel,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (see_epel,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

      call mpi_bcast (sei,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
cdjm Jan2017
      call mpi_bcast (sei_eapl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sei_empl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sei_eipl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (sei_eppl,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

      inum = nstrad
      call mpi_bcast (volsumn,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      inum = 3*nstrad
      call mpi_bcast (volsumm,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      inum = nstrad
      call mpi_bcast (volsumee,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)
      call mpi_bcast (volsumei,inum,MPI_DOUBLE_PRECISION,
     .                0,MPI_COMM_WORLD,ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)

      return
      end subroutine eirene_broadcast_eirbra

      subroutine eirene_mpisend_eirbra(istr,irank)
        integer, intent(in) :: istr  !< stratum index
        integer, intent(in) :: irank !< where to send
        integer :: ierr
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
        external :: mpi_send
#endif

        call mpi_send(volsumn(istr), 1, MPI_DOUBLE_PRECISION,
     .                irank, istr, MPI_COMM_WORLD, ierr)
        call mpi_send(volsumm(:,istr), 1, MPI_DOUBLE_PRECISION,
     .                irank, istr, MPI_COMM_WORLD, ierr)
        call mpi_send(volsumee(istr), 1, MPI_DOUBLE_PRECISION,
     .                irank, istr, MPI_COMM_WORLD, ierr)
        call mpi_send(volsumei(istr), 1, MPI_DOUBLE_PRECISION,
     .                irank, istr, MPI_COMM_WORLD, ierr)
      return
      end subroutine eirene_mpisend_eirbra

      subroutine eirene_mpirecv_eirbra(istr,irank)
        integer, intent(in) :: istr !< stratum index
!> rank from where the message was sent
        integer, intent(in) :: irank
        integer :: ierr
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
        external :: mpi_recv
#endif

        call mpi_recv(volsumn(istr), 1, MPI_DOUBLE_PRECISION,
     .                irank, istr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ierr)
        call mpi_recv(volsumm(:,istr), 1, MPI_DOUBLE_PRECISION,
     .                irank, istr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ierr)
        call mpi_recv(volsumee(istr), 1, MPI_DOUBLE_PRECISION,
     .                irank, istr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ierr)
        call mpi_recv(volsumei(istr), 1, MPI_DOUBLE_PRECISION,
     .                irank, istr, MPI_COMM_WORLD, MPI_STATUS_IGNORE,
     .                ierr)
      return
      end subroutine eirene_mpirecv_eirbra
csw
      END MODULE EIRMOD_EIRBRA
