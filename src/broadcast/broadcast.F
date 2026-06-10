!  30.08.06:  broadcast of density models and spectra introduced
!             broadcast of logicals for living and missing tallies added
!             broadcast of reaction data adopted to new structure
!  10.01.07:  increase number of entries in LLOGAU
!  02.03.07:  remove escd2* arrays
!  20.06.07:  MAXLEVEL used for broadcast of RSPLST and ISPLST
!  20.06.07:  NCONA used for broadcast of RCONA
!  20.06.07:  NSPEZV_DIM used for broadcast of NSPEZV
!  20.06.07:  NLOGAU used for broadcast of LLOGAU
!  01.07.09:  broadcast of HFTR0-3 arrays removed. this is done in broadref.f
!  01.03.11:  NSPEZV_DIM removed
cdr 15.10.14:  renaming of arrays for variances for sum over strata 'smestl' spectrum tallies
cdr  JAN  16:  additional species index for eplds-->eplei, eplpi

cdr  unification of naming conventions for electron impact collisions

cdr  sept 16:  ETH (collision threshold energy) added to reaction data
cdr            RTMAX and ERTMAX added to reaction data: max. of "rate" sigma(v_rel)*v_rel
cdr            broadcast data for extrapolation from tables or fits, independent of IFIT
cdr  Nov  16:  nmds --> nmei, nids --> niei.
cdr  Nov  16:  mxcolls --> mstor0
cdr  Jan  17:  only comments
cdr  July 17:  bug fix: dimensioning of LCUT(0:N2NDPLGS) corrected
cdr            remove NCHORD (is: NCHOR)
c    Aug. 17:  NMODE, LSMOPRO: exception wrt. MPI. Why necessary?
c              broadcasting of CHRTLS was done twice. Removed once.
cpb  Dec. 17:  remove type SPECT_ARRAY, not needed in Fortran 2003
cpb  Jan. 18:  remove unused arrays TEDTEDX, TEDTEDY, TEDTEDZ
c    Jan. 18:  new submodule alloc_fit_form used to allocate, and initialize REACDAT(IR)
cdr  May 18 :  broadcast new variables for internal CR code (currently H_COLRAD):
cdr            nhcol_store
cdr            m_hcol(nreac)
cdr  Sept 18:  redundant arrays: JEREARC, JEREAEI removed
cdr            NHVREI (formerly: NREAHV)
cdr            NHVRPI (formerly: NRHVPI)
cdr Oct.  18:  remove CHELP, towards more rational (and unified)
cdr            reading of AM data in input.f
cdr            and setting of "density models" in block 5.
cdr            Line emissivity rates are now just ordinary A&M reactions: no. IRC.
cdr            Population escape factors POP_ESC are
cdr            now only via IRC reaction card, block 4. Not in block 12 any more.
cdr            "Density Models" are now charcterized by a reaction IRC (from block 4)
cdr            and a donor state ISPZ, ITYP, rather than by a full reaction deck CHELP.
cdr  Oct. 18:  input tallies on PLSTLS(NINPTL). Includes 18 old input tallies but
cdr            now also derived tallies: EDRIFT, BV_VEC, PMOM_VEC
cdr  tbd:      broadcast: livtali etc. move to correct position
cdr  Jan. 19:  separate routine for broadcast of CCOUPL
cdr  ???       apparently also COMNNL removed here from broadcasting
cdr  Nov. 19:  bugfix re %poly% dimensioning. ND --> ND1. Now: separate dimensioning
cdr            of poly%... and tab1d%... or tab2d%... data.

      SUBROUTINE EIRENE_BROADCAST
cdr
c     tbd: some text here, about logic of this code ??
c
cdr
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
     >    , ONLY : EIRENE_BROADCAST_PARMMOD
      USE EIRMOD_CESTIM
     >    , ONLY : EIRENE_BROADCAST_CESTIM
      USE EIRMOD_COMUSR
     >    , ONLY : EIRENE_BROADCAST_COMUSR
      USE EIRMOD_CADGEO
     >    , ONLY : EIRENE_BROADCAST_CADGEO
      USE EIRMOD_CAI
     >    , ONLY : EIRENE_BROADCAST_CAI
      USE EIRMOD_CCONA
     >    , ONLY : EIRENE_BROADCAST_CCONA
      USE EIRMOD_CLOGAU
     >    , ONLY : EIRENE_BROADCAST_CLOGAU
      USE EIRMOD_CINIT
     >    , ONLY : EIRENE_BROADCAST_CINIT
      USE EIRMOD_COMSIG
     >    , ONLY : EIRENE_BROADCAST_COMSIG
      USE EIRMOD_CPOLYG
     >    , ONLY : EIRENE_BROADCAST_CPOLYG
      USE EIRMOD_CGRID
     >    , ONLY : EIRENE_BROADCAST_CGRID
      USE EIRMOD_CZT1
     >    , ONLY : EIRENE_BROADCAST_CZT1
      USE EIRMOD_CTRCEI
     >    , ONLY : EIRENE_BROADCAST_CTRCEI
      USE EIRMOD_CCOUPL
     >    , ONLY : EIRENE_BROADCAST_CCOUPL
      USE EIRMOD_CGEOM
     >    , ONLY : EIRENE_BROADCAST_CGEOM
      USE EIRMOD_CSDVI
     >    , ONLY : EIRENE_BROADCAST_CSDVI
      USE EIRMOD_CTETRA
     >    , ONLY : EIRENE_BROADCAST_CTETRA
      USE EIRMOD_COMSOU
     >    , ONLY : EIRENE_BROADCAST_COMSOU
      USE EIRMOD_CSTEP
     >    , ONLY : EIRENE_BROADCAST_CSTEP
      USE EIRMOD_COMSPL
     >    , ONLY : EIRENE_BROADCAST_COMSPL
      USE EIRMOD_CTEXT
     >    , ONLY : EIRENE_BROADCAST_CTEXT
      USE EIRMOD_CLGIN
     >    , ONLY : EIRENE_BROADCAST_CLGIN
      USE EIRMOD_COUTAU
     >    , ONLY : EIRENE_BROADCAST_COUTAU
      USE EIRMOD_COMXS
     >    , ONLY : EIRENE_BROADCAST_COMXS
      USE EIRMOD_CTRIG
     >    , ONLY : EIRENE_BROADCAST_CTRIG
      USE EIRMOD_COMNNL
     >    , ONLY : EIRENE_BROADCAST_COMNNL
      USE EIRMOD_PHOTON
     >    , ONLY : EIRENE_BROADCAST_PHOTON
      USE EIRMOD_CREF
     >    , ONLY : EIRENE_BROADCAST_CREF
      USE EIRMOD_TIMEA, ONLY : EIRENE_TIMEA0_OC
      USE EIRMOD_CPES
     >    , ONLY : MY_PE, NLIDENT,
     >             input_distribution_strategy
      USE EIRMOD_MPI
      IMPLICIT NONE
      INTEGER :: IER
      EXTERNAL :: EIRENE_ALLOCATE_MODULES, EIRENE_BROADCAST_USR,
     .            EIRENE_CHECK_EXIT
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

      CALL EIRENE_CHECK_EXIT
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

!pb  in order to avoid cyclic dependencies in compilation
!pb  hand over processor number via argument list

      CALL EIRENE_BROADCAST_PARMMOD(MY_PE)
      CALL EIRENE_BROADCAST_COMSOU(MY_PE)

      CALL EIRENE_ALLOCATE_MODULES

      CALL EIRENE_BROADCAST_CCONA(MY_PE)

      CALL EIRENE_BROADCAST_CGEOM(MY_PE)
      CALL EIRENE_BROADCAST_CGRID(MY_PE)
      CALL EIRENE_BROADCAST_CPOLYG(MY_PE)
      CALL EIRENE_BROADCAST_CTRIG(MY_PE)
      CALL EIRENE_BROADCAST_CTETRA(MY_PE)
      CALL EIRENE_BROADCAST_CADGEO(MY_PE)
      CALL EIRENE_BROADCAST_CLGIN(MY_PE)

      CALL EIRENE_BROADCAST_COMUSR(MY_PE)
      CALL EIRENE_BROADCAST_CINIT(MY_PE)
      CALL EIRENE_BROADCAST_CZT1(MY_PE)
      CALL EIRENE_BROADCAST_CAI(MY_PE)
      CALL EIRENE_BROADCAST_CESTIM(MY_PE)
      CALL EIRENE_BROADCAST_CLOGAU(MY_PE)
      CALL EIRENE_BROADCAST_COMSIG(MY_PE)
      CALL EIRENE_BROADCAST_COMSPL(MY_PE)
      CALL EIRENE_BROADCAST_COMXS(MY_PE)
      CALL EIRENE_BROADCAST_COUTAU(MY_PE)
      CALL EIRENE_BROADCAST_CSDVI(MY_PE)
      CALL EIRENE_BROADCAST_CSTEP(MY_PE)
      CALL EIRENE_BROADCAST_CTEXT(MY_PE)
      CALL EIRENE_BROADCAST_CTRCEI(MY_PE)

      CALL MPI_BCAST (NLIDENT,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

      CALL MPI_BCAST (input_distribution_strategy,1,MPI_INTEGER,0,
     .                MPI_COMM_WORLD,ier)

cdr  broadcast arrays for interfacing with B2/B2.5 codes family
      CALL EIRENE_BROADCAST_CCOUPL(MY_PE)

      CALL EIRENE_BROADCAST_COMNNL(MY_PE)

      CALL EIRENE_BROADCAST_USR(MY_PE)

      CALL EIRENE_BROADCAST_PHOTON(MY_PE)
      CALL EIRENE_BROADCAST_CREF(MY_PE)

c     ------------------------------------------------------------     c
c     for the trace ion module
!pb   CALL EIRENE_BROAD_TIM_PARMS
c     ------------------------------------------------------------     c
c
      CALL EIRENE_CHECK_EXIT
      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

cOS   now call the octree-builder to build an octree on every node except
c     on the "root" node, where this is already done via timea0 after input
      IF (MY_PE .NE. 0) THEN
        CALL EIRENE_TIMEA0_OC()
      END IF

      RETURN
      END SUBROUTINE EIRENE_BROADCAST
