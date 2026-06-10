!pb 300806  reduce commands for WTOTE and EELFI added
!           reduction of spectrum data corrected
!pb 181206  group management changed
!pb 110707  calls to mpi_reduce corrected
!pb 060309  mpi_real8 --> mpi_double_precision
!pb 090309  rewritten to use automatic arrays as output buffer in mpi_reduce
!pb 090309  loops reorganized
!pb 270309  typos corrected

cdr Nov. 15:  comments needed. copv tallies: variances for coupling ??
cdr                            to be checked again after changes in 2013
cdr dec. 15:  eppli: now resolved wrt. species index ipls, added
cdr july 17:  comments re. call to user routine: calstr_usr.
cpb Dec. 17:  remove type SPECT_ARRAY, not needed in Fortran 2003
cdr Mar. 19:  further comments,
cdr           tbd: sync coding for estiml variances with those of other variance tallies
cdr           tbd: Make allocatable: helps, dummys
cdr           Dimensioning of covariance tallies is likely incorrect:
cdr           Can cause problems in case of few cells (0-d "box" cases)
cdr           with many covariance tallies.

C> \brief Collect results from worker processes onto master process of
C> stratum
C>
C> This subroutine is called from MCARLO.f, from within strata loop, at
C> the end of each stratum, if there are more than one processor
C> working on any strata.
C>
C> It collects data from processors belonging to one particular stratum
C> istra (COMPRT) and stores merged data for output tallies for stratum
C> istra on the master process for this stratum npesta(istra).
C>
C> Input via modules:
C> - istra  (stratum number, from common COMPRT)
C> - npesta(istra)  : number of master processor for stratum istra
C> - npestr(istra)  : total no. of processors working on stratum istra
C> - tallies
      SUBROUTINE EIRENE_CALSTR

      USE EIRMOD_PRECISION, ONLY: DP
      USE EIRMOD_PARMMOD, ONLY: NATM, NION, NMOL, NPHOT, NPLS,
     .                          NLMPGS, NCV, NLIMPS,
     .                          NRTAL, NADSPC, NVOLTL, NSDW, NSD, NSRFTL
      USE EIRMOD_CAI, ONLY: XMCT
      USE EIRMOD_COMUSR, ONLY: NATMI, NIONI, NMOLI, NPHOTI, NPLSI
      USE EIRMOD_CESTIM, ONLY: ESTIML, ESTIMS, ESTIMV
      USE EIRMOD_CSPEZ, ONLY: LOGATM, LOGION, LOGMOL, LOGPHOT, LOGPLS
      USE EIRMOD_COMPRT, ONLY: ISTRA
      USE EIRMOD_CPES, ONLY: NEED_CALSTR, CALC_STRATUM,
     .                       GET_STRATUM_COMM, I_AM_LEADER,
     ,                       STRATEGY_BALANCED, STRATEGY_IMPROVED,
     ,                       STRATEGY_SIMPLE,
     ,                       WORK_DISTRIBUTION_STRATEGY
      USE EIRMOD_CSDVI, ONLY: NSIGI_SPC, SDVI1, SDVI2, SIGMAC, SGMCS
      USE EIRMOD_COUTAU
      USE EIRMOD_CALSTR_BUFFERED
      USE EIRMOD_MPI

      IMPLICIT NONE

      real(dp), allocatable :: helpest(:), helpv(:), dummyv(:)
      real(dp) :: helpa(0:natm), helpm(0:nmol), helpi(0:nion),
     .            helpp(0:npls), helpph(0:nphot),
     .            helps(nlmpgs+1), helpc(1)
      real(dp) :: dummys(nlmpgs+1)
      integer :: calstr_comm
      integer :: ier1, ier, ir, i, ispc, my_pe_gr,
     .           mxdim, ns, j, istr
      logical, allocatable :: lhelp(:)
      logical :: lhelpa(0:natm), lhelpm(0:nmol), lhelpi(0:nion),
     .           lhelpp(0:npls), lhelpph(0:nphot)
      external :: eirene_calstr_usr, eirene_masage, eirene_exit_own
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      external :: mpi_reduce
#endif

C This subroutine is only called if PROCFORSTR(ISTRA,MY_PE), check is duplication.
C Could not check whether need_calstr(istra) can be moved outside of this subroutine.
      if( need_calstr(istra) .and. calc_stratum(istra)) then
CDR  more than one single processor was active on this stratum ISTRA,
CDR  and my_pe is one of them

!pb     istra is a pointer, type check failure with Intel compiler under Windows
        istr = istra
        calstr_comm = get_stratum_comm(istr)

c  my_pe_gr=0 indicates: my_pe is the master processor for istra
C Would it make more sense to turn my_pe_gr into a logical?
C Need to clarify what eirene_calstr_usr does with my_pe_gr.
        call mpi_comm_rank(calstr_comm, my_pe_gr, ier)
        if (ier /= mpi_success) then
          call eirene_masage
     .     ('ERROR IN SUBROUTINE EIRENE_CALSTR at mpi_comm_rank.')
          call eirene_exit_own(1)
        end if

        if (work_distribution_strategy == STRATEGY_BALANCED .or.
     &      work_distribution_strategy == STRATEGY_SIMPLE   .or.
     &      work_distribution_strategy == STRATEGY_IMPROVED) then
          call eirene_calstr_buffered(calstr_comm, my_pe_gr,
     &                                I_am_leader(istra))
          ! Note that if eirene_calstr_usr has any variables to reduce,
          ! then they must be added to eirmod_buffered_calstr. Otherwise
          ! we lose the advantage of non-blocking calls, and we probably
          ! run into a deadlock.
          ! call eirene_calstr_usr (my_pe_gr, calstr_comm)
          return
        end if

cdr missing: wtotph ??

        call mpi_reduce(WTOTM(0:nmoli,istra),helpm,nmoli+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) WTOTM(0:nmoli,istra) = helpm(0:nmoli)

        call mpi_reduce(WTOTA(0:natmi,istra),helpa,natmi+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) WTOTA(0:natmi,istra) = helpa(0:natmi)

        call mpi_reduce(WTOTI(0:nioni,istra),helpi,nioni+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) WTOTI(0:nioni,istra) = helpi(00:nioni)

        call mpi_reduce(WTOTP(0:nplsi,istra),helpp,nplsi+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) WTOTP(0:nplsi,istra) = helpp(0:nplsi)

        call mpi_reduce(WTOTE(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) WTOTE(istra) = helpc(1)

        call mpi_reduce(XMCP(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) XMCP(istra) = helpc(1)

        call mpi_reduce(XMCT(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) XMCT(istra) = helpc(1)

        call mpi_reduce(PTRASH(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) PTRASH(istra) = helpc(1)

        call mpi_reduce(ETRASH(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) ETRASH(istra) = helpc(1)

        call mpi_reduce(ETOTA(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) ETOTA(istra) = helpc(1)

        call mpi_reduce(ETOTM(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) ETOTM(istra) = helpc(1)

        call mpi_reduce(ETOTI(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) ETOTI(istra) = helpc(1)

        call mpi_reduce(ETOTP(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) ETOTP(istra) = helpc(1)

        call mpi_reduce(EELFI(0:nioni,istra),helpi,nioni+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) EELFI(0:nioni,istra) = helpi(0:nioni)

c  particle balance tallies: from bulk (ipls) to species a,m,i,ph,pl
        call mpi_reduce(PPATI(0:natmi,istra),helpa,natmi+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) PPATI(0:natmi,istra) = helpa(0:natmi)

        call mpi_reduce(PPMLI(0:nmoli,istra),helpm,nmoli+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) PPMLI(0:nmoli,istra) = helpm(0:nmoli)

        call mpi_reduce(PPIOI(0:nioni,istra),helpi,nioni+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) PPIOI(0:nioni,istra) = helpi(0:nioni)

        if (nphoti > 0) then
          call mpi_reduce(PPPHTI(0:nphoti,istra),helpph,nphoti+1,
     .         mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
          if (my_pe_gr==0) PPPHTI(0:nphoti,istra) = helpph(0:nphoti)
        end if

        call mpi_reduce(PPPLI(0:nplsi,istra),helpp,nplsi+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) PPPLI(0:nplsi,istra) = helpp(0:nplsi)

c  energy balance tallies: from bulk (ipls) to species a,m,i,ph,pl
        call mpi_reduce(EPATI(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) EPATI(istra) = helpc(1)

        call mpi_reduce(EPMLI(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) EPMLI(istra) = helpc(1)

        call mpi_reduce(EPIOI(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) EPIOI(istra) = helpc(1)

        call mpi_reduce(EPPHTI(istra),helpc(1),1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) EPPHTI(istra) = helpc(1)

        call mpi_reduce(EPPLI(0:nplsi,istra),helpp,nplsi+1,
     .       mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        if (my_pe_gr==0) EPPLI(0:nplsi,istra) = helpp(0:nplsi)

c  all volume-averaged tallies: estimv
        if (my_pe_gr==0) then
          call mpi_reduce(MPI_IN_PLACE,estimv,nvoltl*nrtal,
     .         mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        else
          call mpi_reduce(estimv,estimv,nvoltl*nrtal,
     .         mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        end if

c  all surface-averaged tallies: estims
        if (my_pe_gr==0) then
          call mpi_reduce(MPI_IN_PLACE,estims,nsrftl*nlmpgs,
     .         mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        else
          call mpi_reduce(estims,estims,nsrftl*nlmpgs,
     .         mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
        end if

c  all energy-resolved ("spectra") tallies: estiml%spc
cdr  different treatment because tallies and their variances are mixed
cdr  into a single data structure, distinct from all other tallies?
        do ispc=1,nadspc
          ns = estiml(ispc)%nspc
          allocate (helpest(ns+2))
          call mpi_reduce(estiml(ispc)%spc(0:ns+1),helpest,
     .                    estiml(ispc)%nspc+2,
     .         mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
          if (my_pe_gr==0) estiml(ispc)%spc(0:ns+1)=helpest(1:ns+2)

c  standard deviation of energy-resolved "spectra": estiml%sdv,...?
          if (nsigi_spc > 0) then
            call mpi_reduce(estiml(ispc)%sdv,helpest,
     .                      estiml(ispc)%nspc+2,
     .           mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
            if (my_pe_gr==0)
     .        estiml(ispc)%sdv(0:ns+1) = helpest(1:ns+2)

            call mpi_reduce(estiml(ispc)%sgm,helpest,
     .                      estiml(ispc)%nspc+2,
     .           mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
            if (my_pe_gr==0)
     .        estiml(ispc)%sgm(0:ns+1) = helpest(1:ns+2)

            call mpi_reduce(estiml(ispc)%sgms,helpc(1),1,
     .           mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
            if (my_pe_gr==0) estiml(ispc)%sgms = helpc(1)
          end if

          deallocate (helpest)
        end do   !nadspc

        mxdim = max(nrtal+1, ncv)
        allocate (helpv(mxdim), dummyv(mxdim))

C  standard deviation of volume-averaged tallies. sdvi1: sigma, sgms
        if (nsd > 0) then
          do ir=1,nsd
            dummyv(1:nrtal+1) = sdvi1(ir,1:nrtal+1)
            call mpi_reduce(dummyv,helpv,nrtal+1,
     .           mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
            if (my_pe_gr==0) sdvi1(ir,1:nrtal+1) = helpv(1:nrtal+1)
          end do
        end if

C  standard deviation of surface-averaged tallies.  sdvi2: sigmaw,sgmws
        if (nsdw > 0) then
          do ir=1,nsdw
            dummys(1:nlimps+1) = sdvi2(ir,1:nlimps+1)
            call mpi_reduce(dummys,helps,nlimps+1,
     .           mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
            if (my_pe_gr==0) sdvi2(ir,1:nlimps+1) = helps(1:nlimps+1)
          end do
        end if

C  covariances between two volume-averaged tallies.  sigmac,sgmcs
        if (ncv > 0) then
          do i=0,2
            do j=1,ncv
              dummyv(1:nrtal) = sigmac(i,j,1:nrtal)
              call mpi_reduce(dummyv,helpv,nrtal,
     .             mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
              if (my_pe_gr==0) sigmac(i,j,1:nrtal) = helpv(1:nrtal)
            end do
          end do

          dummyv(1:ncv) = sgmcs(0,1:ncv)
cdr in case of nrtal < ncv: crash ?
          call mpi_reduce(dummyv,helpv,ncv,
     .         mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
          if (my_pe_gr==0) sgmcs(0,1:ncv) = helpv(1:ncv)

          dummyv(1:ncv) = sgmcs(1,1:ncv)
          call mpi_reduce(dummyv,helpv,ncv,
     .         mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
          if (my_pe_gr==0) sgmcs(1,1:ncv) = helpv(1:ncv)

          dummyv(1:ncv) = sgmcs(2,1:ncv)
          call mpi_reduce(dummyv,helpv,ncv,
     .         mpi_double_precision,mpi_sum,0,calstr_comm,ier1)
          if (my_pe_gr==0) sgmcs(2,1:ncv) = helpv(1:ncv)
        end if

        mxdim = max (nmoli+1,natmi+1,nioni+1,nphoti+1,nplsi+1)
        allocate (lhelp(mxdim))

        call mpi_reduce(LOGMOL(0:nmoli,ISTRA),lhelpm,NMOLI+1,
     .       mpi_logical,mpi_LOR,0,calstr_comm,ier1)
        if (my_pe_gr==0) LOGMOL(0:nmoli,ISTRA) = lhelpm(0:nmoli)

        call mpi_reduce(LOGATM(0:natmi,ISTRA),lhelpa,NATMI+1,
     .       mpi_logical,mpi_LOR,0,calstr_comm,ier1)
        if (my_pe_gr==0) LOGATM(0:natmi,ISTRA) = lhelpa(0:natmi)

        call mpi_reduce(LOGION(0:nioni,ISTRA),lhelpi,NIONI+1,
     .       mpi_logical,mpi_LOR,0,calstr_comm,ier1)
        if (my_pe_gr==0) LOGION(0:nioni,ISTRA) = lhelpi(0:nioni)

        if (nphoti > 0) then
          call mpi_reduce(LOGPHOT(0:nphoti,ISTRA),lhelpph,NPHOTI+1,
     .         mpi_logical,mpi_LOR,0,calstr_comm,ier1)
          if (my_pe_gr==0) LOGPHOT(0:nphoti,ISTRA) = lhelpph(0:nphoti)
        end if

        call mpi_reduce(LOGPLS(0:nplsi,ISTRA),lhelpp,NPLSI+1,
     .       mpi_logical,mpi_LOR,0,calstr_comm,ier1)
        if (my_pe_gr==0) LOGPLS(0:nplsi,ISTRA) = lhelpp(0:nplsi)

        deallocate(lhelp)
c
c  collect user or case-specific information from all PEs that worked on
c  stratum no. ISTRA.  Depends on ...usr.f  or ...cop.f routines.
c  Strictly there should also be an analogue call to eirene_calstr_cop.f

        call eirene_calstr_usr (my_pe_gr, calstr_comm)

      endif

      if (allocated(helpv)) deallocate (helpv)
      if (allocated(dummyv)) deallocate (dummyv)
      RETURN
      END SUBROUTINE EIRENE_CALSTR
