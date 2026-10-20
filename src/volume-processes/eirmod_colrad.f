cdr nov. 21:  in he: pop excludes ne factor (p2),
cdr           in h : pop already includes ne factor
cdr dec. 21:  done : now pop always with ne factor.
cdr           done : separate atomic structure data from CR-rates data structure
cdr           done : rename variables, f1, f2, f1c (f1-condensed to f2)
cdr dec- 21:  he:  missing still: Form 1  (the current form 1 code
cdr                actually provides form 2 data (namely: f1c data). 
cdr           sanity checks: compare f2 and f1c, as well as two formulas for e-cool

cdr  indirect adressing istor_crm  where set? meaning?
cdr  iflav and icrm should be made independent, because
cdr  for one iflav (e.g. He) there can be two (iform=1,2) different CR models

      MODULE EIRMOD_COLRAD
      use EIRMOD_precision
      use EIRMOD_comxs
      IMPLICIT NONE
      PRIVATE

c  population coefficients output from H_colrad(....) via parameter list
      real(dp), allocatable, save :: pop0(:), pop1(:), pop_ext(:)

      real(dp), allocatable, save ::  aikeinst(:,:),energlev(:),
     .                                statwght(:)

      PUBLIC :: eirene_colrad, eirene_colrad_reinit,
     .          eirene_dealloc_colrad

cdr jan 18:  distinct from solps4.3 version: e_alpcr is correct now.
cdr          (electron cooling/heating terms associated with recombination)
cdr feb 18:  l_ext, q_ext, lopaque, pop_esc: must not change,
c            after first call to CR model ICRM, otherwise: reset LVIS
c            so far: q_ext is not connected (l_ext=.false.)
cdr may 18:  add population escape factors pop_esc(40,40), for hydrogen atom.
cdr          default: optically thin: pop_esc=1
cdr dec.18:  additional flag: IFORM: MS condensed vs. MS resolved. Not fully available.
c            I.e.: so far: iform=1 is not connected (iform=2 always).


      type cls
c  stored output from CR code model no. ICRM: h_stor(1:nhcol_stor,1:nrad),
c                               and atomic structure: energy-levels, stat weights, Aik,...
c  input for   CR codes: population escape factors, external sources for QSS states
cdr storage footprint: these next tallies are usually far too large
        real(dp), allocatable :: h_stor(:,:), pop_esc(:,:), q_ext(:,:),
     .                           aikeinst(:,:), energlev(:), statwght(:)

c  Try to avoid repeated calls to same CR model in same plasma grid cell
c  for the current run/iteration/time-cycle
        logical, allocatable :: lvisit(:)
c  do we have external components for CR matrix ICRM
        logical :: l_ext
      end type cls
      type(cls), save :: cr(maxcrm)

      integer, save :: icr_dim(maxcrm) = (/ 40, 1, 1, 65 /)

      real(dp) :: ALPCR, SCR, SCR_EXT,  ALP_EXT, 
     .            E_ALPCR, E_SCR, E_SCR_EXT, E_ALP_EXT,
     .            E_ALPCR_T, E_SCR_T, E_SCR_EXT_T

      CONTAINS

      subroutine eirene_colrad (ir, icrm, res, 
cdr optional, not needed for setting up CRM in initialization
     .                          iform, ivar, icell, p1, p2)

!   driver routine for collisional-radiative models.
!     Calls internal CR code no. icrm, for cell no. ICELL
!     Keeps all results from this call and
!     marks the cells already visited (for icrm: cr(icrm)%lvisit),
!     to avoid double calls to CR code for one and the same cell,
!     or for two or more differenct CR code output quantities

!   input:
!   ir <= 0    only initialize CR model ICRM, and reset COLRAD density-model field particles IPLS
!              by calls to plasma_deriv(ipls). iabs(ir)=iteration no. if previous iteration 
!   ir:        additionally: reaction number KK, as stored in eirene REACDAT input data structures.
!   icrm:      choice of internal CR model.
!              Currently icrm=1: H-colrad
!                        icrm=4: He-colrad
!              Tbd:      icrm=2: H2-colrad
!                        icrm=3: Li-colrad,  etc.
!   iform:     flag for CR-condensed vs. CR-uncondensed choices
cdr next parameters: only when called with ir gt 0
!              (formulation I, II,....), i.e. "meta-stable" resolved or not
!   ivar:      this call to colrad picks one particular CR variable no. IVAR for cell ICELL,
cdr            e.g. RES= src,  pop0, e_alpd, etc..
cdr            In the first call of CR model ICRM/IFORM in cell ICELL,
cdr            also all further NSTOR_CRM parameters from this CR code will be saved,
cdr            on list: CR(...
cdr            NSTOR_CRM was already set in first pass SETAMD ?
!
!
!   icell: cell for which collisional-radiative model should be calculated
!   p1:    first parameter (usually: log_e temperature,...)
!   p2:    second parameter (if any, e.g. log_e (density),...,log_e(test particle energy),...)

!   output:
!   res:   result, for cell no. icell.

      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_cinit, only: cdenmodel
      use EIRMOD_comprt, only: iunout
cdr singlet-triplet mixing in He  CR code requires local B-field strength (T)
      use EIRMOD_comusr, only: bfin, diin
      use EIRMOD_hecr

      implicit none

      integer, intent(in) :: ir, icrm
      integer, intent(in), optional :: iform, icell, ivar
      real(dp), intent(in), optional  :: p1, p2
      real(dp), intent(out) :: res

      real(dp) :: E_ALPCR_T, E_SCR_T, E_SCR_EXT_T   these arrays are for testing only
      real(dp), allocatable :: q_ext(:)
      real(dp) :: mfield
      real(dp), save :: sum2=0.0, sum3=0.0, sumsrc=0.0, sumalp=0.0
      logical  :: l_ext

      integer :: i, nesc, irow_esc, icol_esc, irc, ipop, jcrm,
     .           ipe, nstor_crm, next, irow_ext, jpls
      integer, save :: icount=0, icount_cr(4)=0

c  corresponding output from He_colrad: partially via EIRMOD_HECR
c                                       partially via parameter list
c ....
      integer :: ndim, jdim

      external :: eirene_h_colrad, eirene_exit_own

      if (ir.le.0) then 
cdr     call eirene_colrad_reinit
        call eirene_dealloc_colrad
        call eirene_leer(1)
        write (iunout,*) 'from previous iteration: ',iabs(ir)
        write (iunout,'(A15,1X,3(es12.4),1x)') 'sum2, sum3     ',sum2,
     .                                     sum3,sum2+sum3
        write (iunout,'(A15,1X,3(es12.4,1x))') 'sumsrc, sumalp ',sumsrc,
     .                                     sumalp,sumsrc+sumalp
        write (iunout,*) 'colrad: reset test sums ',ir
        sum2=0.0
        sum3=0.0
        sumsrc=0.0
        sumalp=0.0
      endif

c  ndim: size (no. of states) of CR model ICRM
c  nstor_crm: no. of output parameters stored per cell for CR model ICRM/IFORM
      ndim = icr_dim(icrm)
      if (.not. allocated(q_ext)) then
        allocate(q_ext(ndim))        
      endif
      nstor_crm = pcrm(icrm)%p%nhcol_store
      res = 0.0
      q_ext = 0.0_DP 

cdr  dwell on data type cls   (=cs(icrm)%...)
      if (.not. allocated(cr(icrm)%lvisit)) then
c  initialize data for internal CR code no. ICRM:
cdr  IFORM must stay the same for all calls with ICRM ? 
cdr  Seems not guaranteed
        call eirene_leer(1)
        write (iunout,*) 'COLRAD: Set internal CR model ICRM ',ICRM

        allocate (cr(icrm)%lvisit(nrad))
        allocate (cr(icrm)%h_stor(nstor_crm,nrad))

        allocate(cr(icrm)%pop_esc(ndim,ndim)) !   line population escape factor (default:==1)
        allocate(cr(icrm)%q_ext(ndim,nrad))   !   e.g. photo excitation rate for H*(n), He*(n)
        allocate(cr(icrm)%energlev(ndim))
        allocate(cr(icrm)%statwght(ndim))
        allocate(cr(icrm)%aikeinst(ndim,ndim))       
c  initialize
        cr(icrm)%lvisit   = .false.
        cr(icrm)%energlev = 0.0_DP
        cr(icrm)%statwght = 0.0_DP
        cr(icrm)%aikeinst = 0.0_DP
        cr(icrm)%l_ext    = .false.
        cr(icrm)%q_ext    = 0.0_DP
        cr(icrm)%pop_esc  = 1.0_DP    ! optically thin CR model

cdr  accumulate all population escape factors for internal CR model no. ICRM.
cdr  Scan all reaction cards (block 4, data type ifit=5).
cdr  Population escape factors may have been specified via H.4 (rate coeff)
cdr  or H.10  (energy loss rates) or H.12  (population coefficients, ratios).
cdr
 60     format(2i6,2x,a4,2x,i6,2x,a3,es12.4)

        write (iunout,*) 'scan all nreac reaction cards ',nreaci

        do irc = 1, nreaci
cdr  scan over all reaction decks (from input block 4).
cdr  for H.4, H.10 and H.12 reaction cards, resp.

          if (reacdat(irc)%lrtc) then  ! rate coeff
cdr  loop over H.4 rate coefficients from CR model
            call eirene_build_cr_matrix(icrm,ndim,irc,
     .                                       reacdat(irc)%rtc,'H.4 ')        
          elseif (reacdat(irc)%lrtcew) then ! energy rate coeff
cdr  loop over H.10 energy rate coefficients from CR model
            call eirene_build_cr_matrix(icrm,ndim,irc,
     .                                       reacdat(irc)%rtcew,'H.10')  
          elseif (reacdat(irc)%loth) then ! other rate coeff, pop_coef
cdr  loop over H.12 data, popul. coeffs, from CR model
            call eirene_build_cr_matrix(icrm,ndim,irc,
     .                                      reacdat(irc)%oth,'H.12')
          end if
        end do  !  loop over IREAC

        if (ir .le. 0) then
cdr  re-initializatio of ICRM, 
cdr  not called from plasma_deriv, nor from rate_coef.f etc.
          do JPLS=1,NPLS
            if (cdenmodel(jpls) == 'COLRAD    ') then
cdr  tbd.: check whether cr-model ICRM is associated with JPLS
              call eirene_plasma_deriv(1,JPLS)
              write (iunout,*) 'field species ',jpls,
     .                         ' reset in eirene_COLRAD)' 
            endif 
          enddo
cdr  when in iterative mode:
cdr  strictly: we only need to dwell on atomic data involving a modified IPLS field particle
          call eirene_setamd(1)

          call eirene_leer(1)
        endif

      end if    ! initialize ICRM
   
cdr
cdr  Now we have defined population escape factor matrix POP_ESC for CR model ICRM, 
cdr  as well as the external source rates (e.g. due to photo excitation).
cdr  Set for all cells, and for all NSTOR_CRM quantities derived from this CR model ICRM.      

      if (ir .le. 0) return

cdr  dwell on individual cells ICELL, and reactions IR,
cdr  called from "rate-coeff.f" routines

      select case(icrm)

      case(1)
! COLLISIONAL-RADIATIVE MODEL OF ATOMIC HYDROGEN, Form 2 only (MS unresolved)
        if (.not.allocated(pop0)) then
          allocate(pop0(40))
          allocate(pop1(40))
          allocate(pop_ext(40))
        end if

        if (.not.cr(icrm)%lvisit(icell))  then

! cell number ICELL has not yet been visited so far in this run for this ICRM
! CR model needs to be calculated.
! In later calls, for this ICELL,
!    we assume IFORM, POP_ESC to be unchanged, for ICRM, and for all ICELL 
!    while Q_EXT and the flag L_EXT may differ from cell to cell.

cdr here set one or more external contributions for parameters p1,p2

          L_EXT = CR(ICRM)%L_EXT

c          write (iunout,*) 'COLRAD icrm,icell ',icrm,icell
          if (l_ext) then
            q_ext(:) = CR(ICRM)%Q_EXT(:,icell)
            if (any(q_ext(:) .ne. 0._dp)) then
c           write (iunout,*) 'icrm incell icell ',icell,q_ext(2),
c     .                                              q_ext(3)
              sum2=sum2+q_ext(2)
              sum3=sum3+q_ext(3)
           else
c   no external source set for this particular cell           
              l_ext=.false.
            endif 
          endif
          
          CALL EIRENE_H_COLRAD(p1, p2, iform,                           ! in: te,ne,
     .                         CR(ICRM)%POP_ESC,                        ! in: population escape factor matrix
     .                         Q_EXT, L_EXT,                            ! in: external source: 1/s
c
     .                         cr(icrm)%aikeinst,cr(icrm)%energlev,     ! out: atomic constants:
     .                         cr(icrm)%statwght,                       !      energies, weights, Aik,...
cdr  in case of He_colrad: these next output parameters are handled via module EIRMOD_HECR.
     .                         POP0, POP1, POP_EXT,                     ! out: reduced population coefficients
     .                         ALPCR,    SCR,    SCR_EXT, ALP_EXT,      ! out: effective rates
     .                         E_ALPCR,  E_SCR,  E_SCR_EXT, E_ALP_EXT,  ! out: effective energy rates
     .                         E_ALPCR_T,E_SCR_T,E_SCR_EXT_T)           ! for testing only

cdr keep atomic structure data: einstein coeff, energy levels, stat weights
cdr should be needed from CR-code ICRM only once, from first call.
          icount_cr(icrm)=icount_cr(icrm)+1
          if (icount_cr(icrm).eq.1) then
            write (iunout,*) 'H CR code: '//
     .                       'i, level energy, stat. weight, A_i1'
            do jdim=1,ndim
              write (iunout,'(I3,1x,3(1PE12.4))')
     .            jdim,cr(icrm)%energlev(jdim),
     .                 cr(icrm)%statwght(jdim),
     .                 cr(icrm)%aikeinst(jdim,1)
            enddo
            call eirene_leer(1)      
          endif
          if (l_ext) then
            sumsrc=sumsrc+scr_ext
            sumalp=sumalp+alp_ext
c          write (iunout,'(A,1X,2(1PE12.4))') ' test ',e_scr,e_scr_t
c          write (iunout,'(A,1X,2(1PE12.4))') ' test ph',e_scr_ext
c     .     ,e_scr_ext_t
          endif 
          call eirene_store_H_colrad_results(icrm,icell,iform)

          cr(icrm)%lvisit(icell) = .true.
        end if  ! visited(icell)?

cdr  Result RES from h_colrad (ICRM=1) in cell ICELL is available.
cdr  it may be a rate coefficient, a rate, an energy loss rate or a reduced population coefficient
        res = cr(icrm)%h_stor(ivar,icell)
        return


      case (4)

! COLLISIONAL-RADIATIVE MODEL OF FOR Helium, available: IFORM = 1 and IFORM=2

        if (.not.cr(icrm)%lvisit(icell))  then

          cr(icrm)%l_ext = .false.
          mfield = bfin(icell)
          call eirene_he_colrad (p1, p2, mfield, iform,              ! in: Te, ne, B-field, FI or FII,...
     .                           cr(icrm)%pop_esc,                   ! in: population escape factors
     .                           cr(icrm)%q_ext, cr(icrm)%l_ext,     ! in: external sources
c
     .                           cr(icrm)%aikeinst,cr(icrm)%energlev,! out: atomic constants: energies, weights, Aik,...
     .                           cr(icrm)%statwght)
! further out: via eirmod_hecr: pop, eff. rates, etc.
c
          icount_cr(icrm)=icount_cr(icrm)+1

cdr keep atomic structure data: einstein coeff, energy levels, stat weights
cdr should be needed only once, from first call only.
          if (icount_cr(icrm).eq.1) then
            write (iunout,*) 'i, level energy, stat weight, A_i1'
            do jdim=1,ndim
              write (iunout,'(I3,1x,3(1PE12.4))')
     .            jdim,cr(icrm)%energlev(jdim),
     .                 cr(icrm)%statwght(jdim),
     .                 cr(icrm)%aikeinst(jdim,1)
            enddo
          endif
          if (icount_cr(icrm).eq.1) then
            write (iunout,*) 'iform, ip, p1,p2,s(1)'
          endif                                                       
          write (iunout,'(I3,1x,4(1PE12.4))') iform,e_ipg,p1,p2,
     .                                       rate%s(1)

          call eirene_store_HE_colrad_results(icrm,icell,iform)

          cr(icrm)%lvisit(icell) = .true.
        end if  ! visited(icell) ?

cdr  Result RES from he_colrad, ICRM, IFORM, is now calculated.
cdr  It may be a rate, an energy loss rate or a reduced population coefficient
        res = cr(icrm)%h_stor(ivar,icell)
        return

      case default   ! currently: icrm .ne. 1 and icrm .ne. 4
         write (iunout,*) ' REQUESTED COLLISIONAL-RADIATIVE MODEL' //
     .                    ' NOT AVAILABLE'
         WRITE (iunout,*) 'icrm ',icrm
         WRITE (iunout,*) 'kk   ',ir
         call eirene_exit_own(1)
      end select

      call eirene_colrad_reinit

      return
      
      END subroutine eirene_colrad

      
!********************************************************************
      
      
      subroutine eirene_build_cr_matrix(icrm, ndim, irc, rp, csw)
cdr  add population escape factors, external sources, etc... 
cdr  from rp=reacdat(irc)%...  to CR matrix system
cdr  CR%(ICRM)%...

cdr  IC_esc: counter for pop-esc factors for CR model ICRM
cdr  CSW: the pop-esc factor has been set in a reaction of type CSW
cdr       from the reaction card no. IRC
cdr  IC_EXT: counter for non-vanishing external sources on grid.
      use EIRMOD_precision
      use EIRMOD_PARMMOD
      use EIRMOD_comusr, only: diin, dein
      use EIRMOD_coutau, only: eirene_fetch_outau, eirene_fetch_outaui
      use EIRMOD_comprt, only: iunout
      use EIRMOD_comxs, only: iftflg
      USE EIRMOD_CCONA, ONLY : ELCHA,EPS60
      USE eirmod_cgrid, ONLY: NSBOX_TAL
      use eirmod_cgeom, only: NCLTAL

      integer, intent(in) :: icrm, ndim, irc
      TYPE(FIT_FORMS),POINTER :: RP   !  --> reacdat(irc)%xxx%   with xxx=rtc, rtcew, oth
      character(4) :: csw
      integer :: nesc, ipe, irow_esc, icol_esc, ic_esc, icell, i
      integer :: next, irow, ic_ext, isp, itl, ist, ispa
      real(dp) :: tallyi, tally(nrad) !  setting one cell at a time 
      
      if (rp%ifit .ne. 5) return
      if (rp%crm%iflav .ne. icrm) return ! currently: icrm=1(H), =4(He)

!   use internal CR code ICRM

cdr set cr(icrm)....%pop_esc(:,:)
      nesc=rp%crm%m_popesc
      if (nesc .gt.0) then
      ic_esc=0
      write (iunout,*) 'COLRAD: build popesc,icrm,nesc ',icrm,nesc
      do ipe=1,nesc
        irow_esc = rp%crm%irow_esc(ipe)
        icol_esc = rp%crm%icol_esc(ipe)
        if ((irow_esc > 0) .and. (icol_esc > 0)) then
           ic_esc=ic_esc+1
           cr(icrm)%pop_esc(irow_esc,icol_esc) = rp%crm%pop_esc(ipe)
           write (iunout,*) 'cr: pop_esc factor set for transition'
           write (iunout,60) irc,csw,irow_esc,'-->',icol_esc,
     .          'to ',cr(icrm)%pop_esc(irow_esc,icol_esc)
        end if
      enddo  ! ipe
 60   format(i6,1x,a4,i6,a3,2x,i6,2x,a3,es12.4)
      endif
cdr  set cr(icrm)...%q_ext from photonic excitation
      next=rp%crm%m_qext  
      
      if (next .gt. 0) then
      write (iunout,*) 'COLRAD: build q_ext,icrm,next ',icrm,next
      do ipe=1,next
        irow = rp%crm%irow_ext(ipe)
        if (irow > 0 .and. irow .le. ndim) then
          isp = rp%crm%ital_ext(1,ipe)
          itl = rp%crm%ital_ext(2,ipe)
          ispa = rp%crm%ital_ext(3,ipe)
          ist = 0
c  sum over strata: iestr=0 ?
          call eirene_fetch_outaui(TALLYI,itl,isp,ist,
     .                             iunout)
          if (tallyi.eq.0) then
            write (iunout,*) 'tally ',isp,itl,'found empty'
            cycle
          else 
            call eirene_fetch_outau (TALLY, itl,isp,ist,iunout)
            write (iunout,*) 'cr: q_ext(:,:) set for state ',irow
            write (iunout,70) irc,csw,isp,itl,ispa,'total: ',tallyi
            cr(icrm)%l_ext=.true.
cdr  
cdr  tally is a rate/cm**3, i.e., contains the absorber density
cdr  convert to rate (1/s) per absorber, to be used in CR code as external source
            do icell=1,NSBOX_TAL
c            write (iunout,*) 'diin ',ispa,icell,diin(ispa,icell),
c     .                  tally(icell)
crc  checks for diin to avoid arbitrary large values of tally
            if (diin(ispa,icell) .ge. 1.0E2_DP) then
              tally(icell)=tally(icell)/diin(ispa,icell)/elcha
            else
crc  will prompt lext to be false above
              tally(icell)=0.0_DP
            end if
            enddo
            cr(icrm)%q_ext(irow,:)    = tally(:)
cdr  try to tell CR code that this is a rate, not a rate coefficient
            iftflg(irc,2)=iftflg(irc,2)+100
          endif 
        end if
      enddo  ! ipe
 70   format(i6,1x,A4,1x,3(i3,1x),A7,1x,es12.4)
      endif
      
      end subroutine eirene_build_cr_matrix
      
      
!********************************************************************

      
      SUBROUTINE eirene_colrad_reinit
cdr this must be done after each internal iteration or time-cycle
      integer :: jcrm

      do jcrm = 1, maxcrm
        if (allocated(cr(jcrm)%lvisit)) then
          cr(jcrm)%lvisit = .false.          
        end if

        if (allocated(cr(jcrm)%h_stor)) then
          cr(jcrm)%h_stor = 0._dp
        end if
      end do

      RETURN
      END SUBROUTINE eirene_colrad_reinit

      SUBROUTINE eirene_dealloc_colrad
      implicit none

      integer :: jcrm

      do jcrm = 1, maxcrm
        if (allocated(cr(jcrm)%lvisit)) deallocate (cr(jcrm)%lvisit)
        if (allocated(cr(jcrm)%h_stor)) deallocate (cr(jcrm)%h_stor)
        if (allocated(cr(jcrm)%pop_esc)) deallocate (cr(jcrm)%pop_esc)
        if (allocated(cr(jcrm)%q_ext)) deallocate (cr(jcrm)%q_ext)
        if (allocated(cr(jcrm)%energlev)) deallocate(cr(jcrm)%energlev)
        if (allocated(cr(jcrm)%statwght)) deallocate(cr(jcrm)%statwght)
        if (allocated(cr(jcrm)%aikeinst)) deallocate(cr(jcrm)%aikeinst)

        if (allocated(pop0)) then
          deallocate (pop0)
          deallocate (pop1)
          deallocate (pop_ext)
        end if
      end do

      return

      end SUBROUTINE eirene_dealloc_colrad


!********************************************************************


      subroutine eirene_store_H_colrad_results(icrm,icell,iform)
cdr store many cr data from cr model icrm in current cell: here: H_CR code
cdr  NMS: number of metastable (p-space) species
cdr  NEXT: number of further external sources

cdr  istor_crm =1:  src(nms.nms)
cdr            =2:  ec1(nms)
cdr            =3:  alpcr(nms)
cdr            =4:  ec0(1)
cdr            =5:  scr_ext
cdr            =6:  ec_ext
c
cdr            =7:  pop1      not ready
cdr            =8:  pop0      not ready
cdr            =9:  pop_ext   not ready

      use EIRMOD_precision
      use EIRMOD_comprt, only: iunout

      integer, intent(in) :: icrm, icell, iform
      integer :: i, ipop, nstor_crm
c
c  up to nstor_crm=...%nhcol_store parameters from the CR model are stored in cell ICELL
C  TBD: if .NOT.L_EXT: only case(1) to case(4) and case(7) to case(16) are available

cdr currently: for H_colrad:

      nstor_crm = pcrm(icrm)%p%nhcol_store

      do i = 1, nstor_crm

        select case(pcrm(icrm)%p%m_hcol(i))
c  effective ionisation rate
          case (1)                  ! H.4  2.1.5
            cr(icrm)%h_stor(i,icell) = scr
          case (2)                      ! H.10 2.1.5
c  electron cooling rate coeff. e_scr is negative from h-colrad
c  note: with delpot=-13.6 (input): this becomes the radiation loss rate coeff. alone
            cr(icrm)%h_stor(i,icell) = -e_scr
c  effective recombination rate
          case (3)                      ! H.4  2.1.8
            cr(icrm)%h_stor(i,icell) = alpcr
          case (4)                      ! H.10 2.1.8
c  electron cooling/heating rate coeff. (both signs possible. loss: negative e_alpcr))
c  note: with delpot=+13.6 (input): this becomes the radiation loss rate coeff. alone
            cr(icrm)%h_stor(i,icell) = -e_alpcr
c  external source driven ionisation rate, e.g. photo-excitation driven ionisation
          case (5)                      ! H.4  2.1PHs
            cr(icrm)%h_stor(i,icell) = scr_ext
          case (6)                      ! H.10 2.1PHs
            cr(icrm)%h_stor(i,icell) = -e_scr_ext

c  population coefficients, coupling to ground state H(1) atom, n=2,30
cdr two different lables are available, for backward compatibility with AMJUEL database

          case (8,22)                   ! H.12  2.1.5b, or 2.1.5.2
            cr(icrm)%h_stor(i,icell) = pop1(2)
          case (7,23)                   ! H.12  2.1.5a, or 2.1.5.3
            cr(icrm)%h_stor(i,icell) = pop1(3)
          case (9,24)                   ! H.12  2.1.5c, or 2.1.5.4
            cr(icrm)%h_stor(i,icell) = pop1(4)
          case (10,25)                  ! H.12  2.1.5d, or 2.1.5.5
            cr(icrm)%h_stor(i,icell) = pop1(5)
          case (11,26)                  ! H.12  2.1.5e, or 2.1.5.6
            cr(icrm)%h_stor(i,icell) = pop1(6)
c  only alternative notation, higher states
          case (27:50)                  ! H.12  2.1.5.7 - 2.1.5.30
            ipop = pcrm(icrm)%p%m_hcol(i) - 20
            cr(icrm)%h_stor(i,icell) = pop1(ipop)

c  population coefficients, coupling to H+ ion

          case (13,51)                  ! H.12  2.1.8b, or 2.1.8.2
            cr(icrm)%h_stor(i,icell) = pop0(2)
          case (12,52)                  ! H.12  2.1.8a, or 2.1.8.3
            cr(icrm)%h_stor(i,icell) = pop0(3)
          case (14,53)                  ! H.12  2.1.8c, or 2.1.8.4
            cr(icrm)%h_stor(i,icell) = pop0(4)
          case (15,54)                  ! H.12  2.1.8d, or 2.1.8.5
            cr(icrm)%h_stor(i,icell) = pop0(5)
          case (16,55)                  ! H.12  2.1.8e, or 2.1.8.6
            cr(icrm)%h_stor(i,icell) = pop0(6)
c  alternative notation, higher states
          case (56:79)                  ! H.12  2.1.8.7 - 2.1.8.30
            ipop = pcrm(icrm)%p%m_hcol(i) - 49
            cr(icrm)%h_stor(i,icell) = pop0(ipop)

c  population coefficients, coupling to external source of excitation (e.g. photons)
c  only available if L_EXT=.TRUE. in call to H_COLRAD

          case (18,80)                  ! H.12  2.1PHb, or 2.1PH.2
            cr(icrm)%h_stor(i,icell) = pop_ext(2)
          case (17,81)                  ! H.12  2.1PHa, or 2.1PH.3
            cr(icrm)%h_stor(i,icell) = pop_ext(3)
          case (19,82)                  ! H.12  2.1PHc, or 2.1PH.4
            cr(icrm)%h_stor(i,icell) = pop_ext(4)
          case (20,83)                  ! H.12  2.1PHd, or 2.1PH.5
            cr(icrm)%h_stor(i,icell) = pop_ext(5)
          case (21,84)                  ! H.12  2.1PHe, or 2.1PH.6
            cr(icrm)%h_stor(i,icell) = pop_ext(6)
c  alternative notation, higher states
          case (85:108)                 ! H.12  2.1PH.7 - 2.1PH.30
            ipop = pcrm(icrm)%p%m_hcol(i) - 78
            cr(icrm)%h_stor(i,icell) = pop0(ipop)
c  external source driven recombination rate, e.g. photo-excitation driven recombnation
          case (109)                      ! H.4  2.1PHa
            cr(icrm)%h_stor(i,icell) = alp_ext
          case (110)                      ! H.10 2.1PHa
            cr(icrm)%h_stor(i,icell) = -e_alp_ext
          case default
            write (iunout,*) ' ERROR IN COLRAD, M_HCOL(I) '
            write (iunout,*) ' REQUESTED RATE FROM H_COLRAD ?? '
            call eirene_exit_own(1)
          end select

        end do

      end subroutine eirene_store_H_colrad_results


!********************************************************************



      subroutine eirene_store_HE_colrad_results(icrm, icell, iform)
cdr store many cr data from cr model icrm in current cell: here: He_CR code

cdr  istor_crm =1:  src
cdr            =2:  ec1
cdr            =3:  alpcr
cdr            =4:  ec0
cdr            =5:  scr_ext
cdr            =6:  ec_ext
c
cdr            =7:  pop1      not ready, also: need a factor p2.
cdr            =8:  pop0      not ready, also: need a factor p2
cdr            =9:  pop_ext   not ready, also: need a factor p2


cdr strictly: iform should not be needed here, 
cdr We already have the indirect addressing istor_crm for CR model ICRM ?

      use EIRMOD_precision
      use EIRMOD_comprt, only: iunout
      use EIRMOD_hecr

      integer, intent(in) :: icrm, icell, iform
      
      integer :: i, ipop, nstor_crm, istor_crm
      integer, save :: icount=0
      real(dp) :: test, radrate_0, radrate_1, radrate_ext
c
cdr
c  up to nstor_crm parameters from the cr-model ICRM will be stored in cell ICELL.
c        istor_crm is the flag for a particular type of data for model ICRM:
c  currently: 
c        istor_crm=1: eff. ionisation rate
c        istor_crm=2: eff. elect. cooling rate, from ionis. component
c                     (taken positive if cooling)
C  TBD: if .NOT.L_EXT: only case(1) to case(16) are available
cdr  

      nstor_crm = pcrm(icrm)%p%nhcol_store

      do i = 1, nstor_crm
        istor_crm = pcrm(icrm)%p%m_hcol(i) !dr indirect addressing. why?

        select case(istor_crm)
c  effective ionisation rate
        case (1)
          select case (iform)
          case (1)                    ! H.4  2.3.9b,c,d,e,f,g,h,i,j
            cr(icrm)%h_stor(i,icell) = crrate%scr_f1c   ! form 1 condensed to form 2
cdr         write (iunout,*) ' problem in colrad with He crm'
cdr         write (iunout,*)
cdr  .              ' scr not ready in formulation I '
cdr         call eirene_exit_own(1)
          case (2)                    ! H.4  2.3.9a
            cr(icrm)%h_stor(i,icell) = crrate%scr_f2
          end select
        case (2)
c  electron cooling rate coeff. e_scr is negative from he-colrad
c  note: with delpot=-24.6 (input): this becomes the radiation loss rate coeff. alone
          select case (iform)
          case (1)                   ! H.10  2.3.9b,c,d,e,f,g,h,i,j
            cr(icrm)%h_stor(i,icell) = coolrate%ec1_f1c   ! form 1 condensed to form 2
            radrate_1                = coolrate%rad1_f1c 
cdr         write (iunout,*) ' problem in colrad with He crm'
cdr         write (iunout,*)
cdr  .              ' e_scr not ready in formulation I '
cdr         call eirene_exit_own(1)
          case (2)                   ! H.10 2.3.9a
            cr(icrm)%h_stor(i,icell) = coolrate%ec1_f2
            radrate_1                = coolrate%rad1_f2 
          end select
c  effective recombination rate
        case (3)
          select case (iform)
          case (1)                   ! H.4  2.3.13b,c,d
            cr(icrm)%h_stor(i,icell) = crrate%alpcr_f1c  ! form 1 condensed to form 2
            write (iunout,*) ' problem in colrad with He crm'
            write (iunout,*)
     .              ' alpcr not ready in formulation I '
cdr         call eirene_exit_own(1)
          case (2)                   ! H.4  2.3.13a
            cr(icrm)%h_stor(i,icell) = crrate%alpcr_f2
          end select
        case (4)
c  electron cooling/heating rate coeff. (both signs possible. loss: negative e_alpcr))
c  note: with delpot=+24.6 (input): this becomes the radiation loss rate coeff. alone
          select case (iform)
          case (1)                   ! H.10  2.3.13b,c,d
            cr(icrm)%h_stor(i,icell) = coolrate%ec0_f1c ! form 1 condensed to form 2
            radrate_0                = coolrate%rad0_f1c ! for sanity check (not saved)
cdr         write (iunout,*) ' problem in colrad with He crm'
cdr         write (iunout,*)
cdr  .            ' term e_alpcr not ready in formulation I '
cdr         call eirene_exit_own(1)
          case (2)                   ! H.10 2.3.13a
            cr(icrm)%h_stor(i,icell) = coolrate%ec0_f2
            radrate_0                = coolrate%rad0_f2 ! for sanity check (not saved)
          end select
c  external source driven ionisation rate, e.g. photo-excitation driven ionisation
        case (5)
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = crrate%scr_ext_f1c
          case (2)                    ! H.4  2.3.9PH
            cr(icrm)%h_stor(i,icell) = crrate%scr_ext_f2
          end select
        case (6)
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = coolrate%ec_ext_f1c
            radrate_ext              = coolrate%rad_ext_f1c 
cdr         write (iunout,*) ' problem in colrad with He crm'
cdr         write (iunout,*)
cdr  .            ' term e_scr_ext not ready in formulation I '
cdr         call eirene_exit_own(1)
          case (2)                    ! H.10 2.3.9PH
            cr(icrm)%h_stor(i,icell) = coolrate%ec_ext_f2
            radrate_ext              = coolrate%rad_ext_f2 
          end select

c  population coefficients, coupling to ground state He(1) atom
cdr AMJUEL H.12 notation.  IFORM=1: unfinished.
        case (7)
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r1(6)
          case (2)                    ! H.12  2.2a  or 2.2.6
            cr(icrm)%h_stor(i,icell) = popcoe%rr1(6)
          end select
        case (8)
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r1(7)
          case (2)                    ! H.12  2.2b  or 2.2.7
            cr(icrm)%h_stor(i,icell) = popcoe%rr1(7)
          end select
        case (9)
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r1(8)
          case (2)                    ! H.12  2.2c  or 2.2.8
            cr(icrm)%h_stor(i,icell) = popcoe%rr1(8)
          end select
        case (10)
          select case (iform)
          case (1)   !  unfinished
            cr(icrm)%h_stor(i,icell) = popcoe%r1(10)
          case (2)                    ! H.12  2.2d  or 2.2.10
            cr(icrm)%h_stor(i,icell) = popcoe%rr1(10)
          end select
        case (11)
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r1(16)
          case (2)                    ! H.12  2.2e
            cr(icrm)%h_stor(i,icell) = popcoe%rr1(16)
          end select

cdr alternative notation: 2.3.9.i,  i=2,59, rather than ...2.3.9a,b,c,d,e
        case (22:79)
          ipop = pcrm(icrm)%p%m_hcol(i) - 20
          select case (iform)
          case (1)
cdr  this should go from i=4,59, but for 3 "trains" of pop1
            cr(icrm)%h_stor(i,icell) = popcoe%r1(ipop)
          case (2)                   ! H.12  2.3.9.2 - 2.3.9.59
            cr(icrm)%h_stor(i,icell) = popcoe%rr1(ipop)
          end select

c  population coefficients, coupling to He+ ion
cdr  In this case, perhaps, IFORM=1 works ? Because only one train of excited states with He+
cdr  Except: population of MS states 2 and 3. These are solved for by transport eqs.
        case (12)
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r0(6)
          case (2)                    ! H.12  2.3.2a
            cr(icrm)%h_stor(i,icell) = popcoe%rr0(6)
          end select
        case (13)                     ! H.12  2.3.2b
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r0(7)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr0(7)
          end select
        case (14)                     ! H.12  2.3.2c
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r0(8)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr0(8)
          end select
        case (15)                     ! H.12  2.3.2d
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r0(10)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr0(10)
          end select
        case (16)                     ! H.12  2.3.2e
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r0(16)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr0(16)
          end select

c  alternative state notation: 2...59, rather than a,b,c,d,e
        case (80:137)                  ! H.12  2.3.13.2 - 2.3.13.59
          ipop = pcrm(icrm)%p%m_hcol(i) - 78
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r1(ipop)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr1(ipop)
          end select

c  population coefficients, coupling to external source of excitation (e.g. photons)
c  only availabel if L_EXT=.TRUE. in call to He_COLRAD
        case (17)                     ! H.12  2.2PHa
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r_ext(6)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr_ext(6)
          end select
        case (18)                     ! H.12  2.2PHb
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r_ext(7)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr_ext(7)
          end select
        case (19)                     ! H.12  2.2PHc
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r_ext(8)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr_ext(8)
          end select
        case (20)                     ! H.12  2.2PHd
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r_ext(10)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr_ext(10)
          end select
        case (21)                     ! H.12  2.2PHe
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r_ext(16)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr_ext(16)
        end select

c  alternative state notation: 2...59, rather than a,b,c,d,e
        case (138:195)                 ! H.12  2.3.9PH.2 - 2.3.9PH.59
          ipop = pcrm(icrm)%p%m_hcol(i) - 136
          select case (iform)
          case (1)
            cr(icrm)%h_stor(i,icell) = popcoe%r1(ipop)
          case (2)
            cr(icrm)%h_stor(i,icell) = popcoe%rr1(ipop)
          end select
        case default
          write (iunout,*) ' ERROR IN COLRAD, M_HCOL(I) '
          write (iunout,*) ' REQUESTED DATA FROM HE_COLRAD ?? '
          write (iunout,*) ' istor_crm = ',istor_crm
          call eirene_exit_own(1)

        end select ! istor_crm

      end do  ! all NSTOR_CRM parameters done, for model ICRM in ICELL


cdr sanity check: consistency between e_cool and rad-rate, coupling to ground state
      if (icount <= 10) then
              test = cr(icrm)%h_stor(1,icell)*E_IPG 
     .                - radrate_1
              write (iunout,*) ' He-colrad, icrm, icell ',icrm,icell
              write (iunout,*) ' rad_scr ',radrate_1
              write (iunout,*) ' e_scr ',cr(icrm)%h_stor(2,icell)
              write (iunout,*) ' test  ',test
              icount = icount + 1
      end if

      end subroutine eirene_store_HE_colrad_results

c*********************************************************
      END MODULE EIRMOD_COLRAD
