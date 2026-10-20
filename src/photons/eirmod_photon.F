cdr  "photon-module" of eirene code:
cdr
cdr  [1] D.Reiter, et al. PPCF 44 (2002) p1723, and JNM 313-316 (2003) p 845
cdr
cdr  [2] S.Wiesen, PhD Thesis, 2005 (Ruhr Univ. Bochum), unpublished
cdr
cdr  [3] V.Kotov, D.Reiter, A.K.Kukushkin et al, Contrib. Plasma Phys.,
cdr      46, No.7-9,(2006)
cdr      V.Kotov, D.Reiter, A.K.Kukushkin, PhD Thesis (Ruhr Univ. Bochum),
cdr      Report JUEL-4257 (Nov 2007) and
cdr      V.Kotov, D.Reiter, Plasma Phys. Control. Fusion 54(2012) 082003
cdr
cdr This version: Jan. 2018:
cdr    re-work (and document) the so-called photon module.
cdr    (contains earlier revisions, cleanup, etc....Detlev Reiter, 2005, 2006)

cdr    The original version was developed in 2000 -- 2003 within a project (BMBF)
cdr    in collaboration with Philips research Lab. Aachen, for applications to
cdr    to high pressure gas discharge lamps.
cdr    It was largely re-written later, re-structured, several times,
cdr    for use in (nonlinear) standalone eirene
cdr    (applications to high density divertor plasmas with hydrogen resonance line re-absorption)
cdr    Further re-writing (2006-2007): V.Kotov, to support iterations for self-consistently
cdr    coupled 2D B2-EIRENE plasma solutions.
cdr    transport code.
cdr
cdr STARTING ONCE AGAIN:
cdr new data structure REACDAT. Try to re-connect photonic reactions to
cdr rest of code, unify notation....
cdr
cdr jan 18: note: ph_xsectp is still there.
cdr         But corresponds to what would be called XSTPH?,
cdr         (what is XSTRC?)
cdr nov 18: notational cleanup: nrot -> nrph. nnrot  is now obsolete.
cdr nov 20: remove alloc_xsectph, all variables now in comxs
crc oct 24: remove eirene_line_cutoff
!................................................................................

!  photon.f  this modules contains routines to sample, evaluate
!            photon line profiles and photonic rates (absorption, emission, etc.)
!            Only to be used for photon tracking in an atomic background.
!            All references to test atom collisions with photon background removed.
!            Still under construction: ph_post_energy, to sample photon
!            parameters after a collision (needed, e.g., to render photon-atom
!            iterations more implicit than presently).
! Hence: currently photon.f works for purely absorbing media for photons.


! 24.2.05: ph_xsectp removed partially, is now xstrc in 'volume-processes', cleaned up
cdr jan 18:  immer noch da. entspricht aber wohl dem xstot, (xstrc) was es nicht gibt.
!  7.3.05: ph_energy exchanged:  comments, cleaned up
!  7.3.05: ph_sam_lorentz replaced by sam_lorentz: bug fix, was wrong
!                         rescaling from Cauchy to Lorentz, alph-->alphh
!  7.3.05: ph_lorentz replaced by lorentz: comments, speed-up
!  7.3.05: ph_homprof replaced by naturalprof: comments, speed-up
!          still wrong: Sum_Aik missing, over all upper and lower levels.
!  8.3.05: doppler, dopplerprof and sam_doppler: cleaned up, commented
! 10.3.05: voigtprof and sam_voigt corrected (wrong Lorentzian FWHM), cleaned up
!          and commented
!  6.4.05: getcoeff: commented, redundant input removed (idsc), ipl out, option 7 out,
!          net abs. rates (wg. stim emis) to be done in calling program
!          new in: option 6:  spont. emission rate for photons (used in sigrad from diagno)
!          old option 6 (spont. emission, from atoms point of view) is now option 1.
!          meaning of iipl: bulk species associated with reaction kk, from calling program.
! 20.4.05: zeeman profiles newly written. speed-up and clean-up,
!          avoid all tr_...routines, probably buggy.
! mai  05: bulk particle drift in normal zeeman triplet (nldrft). also locate corrected
!          for doppler broadening from an-isotropic bulk distributions
! 27.6.05  phv_nrota, phv_nrotph  removed
! june 05: lorvdw corrected and cleaned up. Now both sampling and evaluation of profile
!          in energy and in wavelength units is correct.
! july 05: lorvdw further rewritten to speed up. Subr. sam_lorvdw split into sam_lorentz
!          and sam_vdwqs (quasistatic vdWaals, red wing).
! 18.8.05: iunout in write statements, index for number of foreign gases in
!          pressure broadening introduced in LORVDWPROF: reaction%ifremd
c  4.1.06: hplnk_bar = hplck/2Pi introduced in ccona, and used here
c          some more speed ups in lorvdwprof. still much more to be done
! 08.2.06: phv_lgprc removed from declaration, no longer needed

! 08.5.06: zm_stark_profile and zm_stark_doppler added, for Lyman_alpha
! 08.5.06: sam_zm_stark added
! 19.12.06: sam_zm_stark rewritten. old version --> sam_zm_stark1
!
!    do be done:
!        replace very unefficient programming of constants
!        unify use of fadeeva function

      MODULE EIRMOD_PHOTON
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_CESTIM
      USE EIRMOD_CGEOM
      USE EIRMOD_CGRID
      USE EIRMOD_CLOGAU
      USE EIRMOD_COMPRT
      USE EIRMOD_COMSOU
      USE EIRMOD_COMUSR
      USE EIRMOD_COMXS
      USE EIRMOD_COUTAU
      USE EIRMOD_CRAND
      USE EIRMOD_CSDVI
      USE EIRMOD_CSPEI
      USE EIRMOD_CSPEZ
      USE EIRMOD_CTEXT
      USE EIRMOD_CTRCEI
      USE EIRMOD_CUPD
      USE EIRMOD_CZT1
      USE EIRMOD_RANF, ONLY: RANF_EIRENE

      IMPLICIT NONE

      PRIVATE


      PUBLIC :: eirene_ph_init, eirene_ph_energy, eirene_ph_getcoeff,
     .    eirene_ph_post_energy,
     .    eirene_ph_lorvdw, eirene_planck,
     .    eirene_ph_b12, eirene_ph_b21

      PUBLIC :: EIRENE_BROADCAST_PHOTON

      integer, public, save  :: phv_muldens
      integer, public, save, allocatable ::
     .   PHV_LGAOT(:,:,:),
cdr  .   LGPHPH(:,:,:),NPHPHI(:),   now: comxs
     .   PHV_NAOTI(:),
cdr  .   IESTPH(:,:,:),  now: comxs
     .   PHV_IESTOTAT(:,:,:),
cdr  .   N1STPH(:,:,:),N2NDPH(:,:,:), now: comxs
     .   PHV_N1STOTAT(:,:,:),PHV_N2NDOTAT(:,:,:),
     .   phv_xistra(:)
csw constants
      real(dp), public, save :: STEFBCON
      real(dp), public, save :: hwvdw

csw external
      integer, external :: eirene_idez

!pb black body removal begin (try to eliminate saturated line core part)
!     unfinished

      real(dp), allocatable, public, save  ::
     .          ecutleft(:,:), ecutright(:,:),
     .          phicut(:,:), phi_rj_left(:,:,:),
     .          phi_rj_right(:,:,:), phi_zero(:,:),phi_inf(:,:),
     .          eintleft(:,:), eintright(:,:), eint_inf(:,:),
     .          xintleft(:,:), xintright(:,:), xint_inf(:,:),
     .          xint_cut(:,:)
      logical, allocatable, public, save :: lsrcpls(:)

!pb black body removal end

      CONTAINS

c I/O & MISC-ROUTINES
      SUBROUTINE EIRENE_PH_INIT(ICAL)
      integer, intent(in) :: ical
      external :: eirene_exit_own

      select case(ical)
csw ICAL == 0
      CASE(0)
c constants
      STEFBCON=8.*PIA**5 / (15. * CLIGHT**3 * HPLNK**3)
c
      return

      case (1)
c ICAL=1
c
      return

      case(2)
c ICAL=2
c
      return

      case(3)
c ICAL=3
c
      return

      case default
         write(iunout,*)
     .  'PHOTON MODULE EIRMOD_(PH_INIT): ICAL=',ical,'N/A'
         call EIRENE_exit_own(1)
      end select
      RETURN
      END SUBROUTINE EIRENE_PH_INIT


c CROSS-SECTIONS, RATES, RATE COEFFICIENTS

      SUBROUTINE EIRENE_PH_GETCOEFF(kkin,isp,ity,icell,iipl,fac,res)
c  evaluate absorption, emission and stim. emission rate coeff.  "res" (cm**3/s)
c  for a photon with energy E=E0, in cell ICELL.
C  This requires evaluation of the absorption line shape profiles "iptype"
C  current version: iptype=0,1,2,3,4,5,6,7,8,9,10,11

c  input:
c          kkin: nrearc(irrc), nreaph(irph), reaction number from input block 4
c          isp :           = iphot, iatm (redundant?)
c          ity : (=ityp),  = 0: test photons point of view
c          ity : (=ityp),  = 1: test atoms point of view (out)
c          iipl: species index of background species for reaction kk
c     derived from kkin:
c          iid = 4  photon absorption
c          iid = 5  stim. emiss. by photon
c          iid = 6  spont. photon emission
c          iid = 7  photon net absorption (absorption - stim. emission)

c  output:
c          fac: value of profile shape function Phi(E)dE, (Phi(E) is normalized to 1)
c          res: value of rate, or rate coeff., at E
c          iid = 4  res= rate coeff. = B12 * E00 * c/4 Pi * Phi(E)
c          iid = 5  res= rate coeff. = B21 * E00 * c/4 Pi * Phi(E)
c          iid = 6  res= rate        = A12 * Phi(E) to be done
c          iid = 7  res= rate coeff. = out
cdr  special treatment of delta function line shape: Phi(E)=1 und E00 weglassen ?
cdr  als: z.b. iid=4  res=B12 *c/4 Pi  --> Planck  bei E (=line center ?)
c
      IMPLICIT NONE
      integer, intent(in) :: kkin,isp,ity,icell,iipl
      real(dp), intent(out) :: fac,res
      integer :: iid,kk,iptype
      real(dp):: gam,e00,
     .           v,dnd,xx,yy,pnue,pnue0,
     .           fwhm,shift,dvdw,drft
c      integer :: ipl2
c      real(dp) :: d1, d2, e1, e2, g1, g2, l0, l00, omega_max,
c     .            omega_min
      real(dp):: ctheta,ctheta2,dbz
      real(dp):: t_e,t_p,t_g
      real(dp):: bx, by, bz, bf
      external :: eirene_exit_own

c  reaction data for process kkin are now loaded into "reaction"
c  includes also reaction%b12.
c
c      if (kkin /= idreac) call EIRENE_get_reaction(kkin)
      kk=kkin
      iid = reacdat(kk)%phr%line%ircart
      fac=0._dp
      hwvdw=0._dp
      ctheta =0._dp
      ctheta2=0._dp

cdr  find $fac(E0)$: the value of the line profile $fac(E)$
cdr  (emission profile or absorption profile) at E0:
cdr  the current photon energy.
cdr  units of fac: 1/eV, because: integral dE fac(E) = 1.0

      select case(iid)
c  case 1,2,3  : atoms point of view in radiation field: out
c  case 4,5,6  : photons point of view in neutral gas field
c  to be done: remove case iid=7 from this routine. And add
c              rates for stim. emission in calling program, e.g.
c              all absorb. and stim emiss rates.

      case(4,5,6)   !iid 4,5,6 only
c     P.2 PH_ABS OT, P.2 PH_STIM OT
         if(lgvac(icell,iipl)) then
            res=0.
            return
         endif

         e00=reacdat(kk)%phr%line%e0  ! line center
         pnue0 = e00*EV2HZ
         pnue  = e0 *EV2HZ

         iptype = reacdat(kk)%phr%line%iprofiletype
         select case(iptype)

         case(0)
c  delta distribution, all mass at e0=e00
cdr         FAC=1./E00  ! cancel e00 factor in rate coeff.
cdr or, old version ??
            FAC=0._DP
            if (abs(E0-e00)/e00.lt.eps12) fac=1._DP
         case(1)
            call EIRENE_dopplerprof(iipl,kk,icell,dnd,drft)
            xx=e0-(e00+drft)
            fac = EIRENE_DOPPLER(xx,dnd)
         case(2)
            call eirene_naturalprof(kk,gam)
            xx = e0-e00
            fac = EIRENE_LORENTZ(xx, gam)
         case(3)
            call EIRENE_voigtprof(iipl,kk,icell,dnd,drft,gam)
            xx=(e0-(e00+drft))/dnd
            yy=gam*0.5_dp/dnd
            fac = REAL(EIRENE_PH_FADDEEVA(xx,yy,dnd), DP)
         case(4)
c  use wavelength scale
c           call  EIRENE_lorvdwprof(icell,fwhm,shift,dvdw,.false.)
c           l00=hpcl/e00
c           l0=hpcl/e0
c           xx = l0-l00
c           hwvdw=fwhm
c           fac = EIRENE_ph_lorvdw(xx, fwhm, shift, dvdw, icell)
c           fac=fac*hpcl/(e0*e0)
c  use energy scale
            call  EIRENE_lorvdwprof(kk,icell,fwhm,shift,dvdw,.true.)
            xx = E0-E00
            hwvdw=fwhm
            fac = EIRENE_ph_lorvdw(-xx, fwhm, -shift, dvdw, icell)
         case(5)
            write (iunout,*) ' ph_getcoeff for photons,'
            write (iunout,*) ' iptype == 5 not implemented yet'
            call EIRENE_exit_own(1)
         case(6,7)
!normal zeeman: delta or doppler components
            gam=0._dp
            call EIRENE_zeeman_normalprof(icell,ctheta2,dbz)
            if(iptype == 6) then
               dnd=0._dp
               drft=0._dp
               fac = EIRENE_zm_profile(e0,ctheta2,dbz,gam,dnd,drft,
     .                                 e00,6)
            elseif (iptype == 7) then
               call EIRENE_dopplerprof(iipl,kk,icell,dnd,drft)
               fac = EIRENE_zm_profile(e0,ctheta2,dbz,gam,dnd,drft,
     .                                 e00,7)
            endif
         case(8,9)
!normal zeeman: lorentz or voigt components
           call EIRENE_zeeman_normalprof(icell,ctheta2,dbz)
           call eirene_naturalprof(kk,gam)
           if(iptype == 8) then
             dnd=0._dp
             drft=0._dp
             fac = EIRENE_zm_profile(e0,ctheta2,dbz,gam,dnd,drft,e00,8)
           elseif (iptype == 9) then
             call EIRENE_dopplerprof(iipl,kk,icell,dnd,drft)
             fac = EIRENE_zm_profile(e0,ctheta2,dbz,gam,dnd,drft,e00,9)
           endif

         case(10,11)
!zeeman-stark, or zeeman-stark-doppler, at photon energy E0

c  Stark  broadening by electrons
           T_e=Tein(icell)
c  Stark  broadening by ions: protons, deuterons, tritons
           T_p=T_e
c  Zeeman splitting
           bx = 0._dp
           by = 0._dp
           bz = 1._dp
           bf = 1._dp
           if (lbxin) bx = bxin(icell)
           if (lbyin) by = byin(icell)
           if (lbzin) bz = bzin(icell)
           if (lbfin) Bf = Bfin(icell)
           ctheta=velx*bx+vely*by+velz*bz
           ctheta2=ctheta*ctheta

           if(iptype == 10) then
c  No Doppler broadening by absorbing gas
c            T_g=0
             v=0.0
c  only one evaluation at E0-E00
c            omega_min=E0-E00
c            omega_max=E0-E00
c            npt=1

             fac = EIRENE_zm_stark_profile(dein(icell),T_e,T_p,T_g,Bf,
     .             ctheta2,v,e00,e0)

c    .             npt,omega_min,omega_max)

           elseif (iptype == 11) then
c  With Doppler broadening by absorbing gas
             call EIRENE_dopplerprof(iipl,kk,icell,dnd,drft)
c            T_g=TIIN(mplsti(iipl),icell)
c            v=0.0
c  only one evaluation at E0-E00
c            omega_min=E0-E00
c            omega_max=E0-E00
c            npt=1

             fac = EIRENE_zm_stark_doppler_profile(dein(icell),
     .             T_e,T_p,T_g,Bf,
     .             ctheta2,v,dnd,drft,e00,e0)

c    .             npt,omega_min,omega_max)
           endif
         case default
           write (iunout,*) 'invalid line shape in ph_getcoeff'
           call eirene_exit_own(1)
        end select ! iptype

      case(1,2,3)  
cdr  atom in radiation field. Currently out.

      end select ! iid

      select case(iid)

         case(4)    ! absorption
            res=e00*reacdat(kk)%phr%line%b12
c  result should have units cm**3/s
c  --> multiply with absorber density to make it a rate!
            phv_muldens=1
            res=res*fac*clight/(4._dp*PIA)

         case(5)    ! stimulated emission
            res=e00*reacdat(kk)%phr%line%b21
c  result should have units cm**3/s
c  --> multiply with upper state density to make it a rate!
            phv_muldens=1
            res=res*fac*clight/(4._dp*PIA)

         case(6)     ! spontaneous emission
c     P.1 AT_AIK OT
           res=reacdat(kk)%phr%line%aik
c  result should have units 1/s,
c  --> do not multiply with density in calling program!
c  --> divide by upper state density to make it a rate coefficient!
           phv_muldens=0
           res=res*fac

        case default   ! unknown line profile
          write(iunout,*)
     .  'PHOTON MODULE EIRMOD_(PH_GETCOEFF): iid=',iid,
     .                    'not in use'
          res=0.
        end select ! iid

c     case(7)  reduce absorption by stimulated emission
c     P.2 PH_CABS OT
c        if(lgvac(icell,iipl)) then
c           res=0.
c           return
c        endif
c
c        g1 = reaction%g1
c        g2 = reaction%g2
c        e00 = reaction%e0
c        e1 = reaction%e1
c        e2 = e00 + e1
c        pnue0 = e00*EV2HZ
c        pnue  = e0 *EV2HZ
c
c        fac = ..... from line profile
c
c        d1=DIIN(iipl,icell)
c        if ((N1STPH(isp,idsc,1) == 4) .and.
c    .       (N1STPH(isp,idsc,3) == 1)) then
c           ipl2=N1STPH(isp,idsc,2)
c           d2=DIIN(ipl2,icell)
c        else
c           d2=0._dp
c           write (iunout,*) ' WARNING !! '
c           write (iunout,*) ' WRONG TYPE OF SECONDARY PARTICLE IN',
c    .                  ' PH_GETCOEFF'
c           write (iunout,*) ' RATE CHANGED ACCORDINGLY '
c        END IF
c
c        res=e00*reaction%B12*(1._dp- g1*d2/g2/d1)
c result should have units cm**3/s --> multiply backgr. density!
c        phv_muldens=1
c        res=res*fac*clight/(4._dp*PIA)
c

c     case default
c        write(iunout,*)'PHOTON MODULE (PH_GETCOEFF): iid=',iid,'not in use'
c        res=0.
c     end select

      return
      END SUBROUTINE EIRENE_PH_GETCOEFF



      FUNCTION EIRENE_ZM_PROFILE(X, CTHETA2, DBZ, GAM,DND,
     .                           DRFT,E00,IPROF)
     .         RESULT(RES)
c  zeeman - profile - splitting
c   input:  x:   photon energy [eV]
c  output:  res: value of zeeman splitted absorption profile
c   iprof:  6:   delta
c           7:   doppler
c           8:   lorentz
c           9:   lorentz+doppler, i.e., voigt
      implicit none
      real(dp), intent(in) :: x, ctheta2,dbz,gam,dnd,drft,e00
      integer, intent(in) :: iprof
      real(dp) :: xx,yy,val,del(-1:1),e00d,res
      integer :: ipol
      external :: eirene_exit_own

      res = 0.
      del(-1)=-dbz
      del(0) = 0._dp
      del(1) = dbz
      do ipol = -1,1
        e00d = e00 + del(ipol)

        select case(iprof)
        case(6)
! delta
          if( abs(x-e00d)/e00 < eps12) then
             val =1.
          else
             val = 0.
          endif
        case(7)
! doppler
          val = EIRENE_doppler(x-(e00d+drft),dnd)
        case(8)
! lorentz
           val = EIRENE_lorentz(x-e00d,gam)
        case(9)
! voigt
          xx = (x-(e00d+drft))/dnd
          yy=gam*0.5_dp/dnd
          val = REAL(EIRENE_PH_FADDEEVA(xx,yy,dnd), DP)
c
c  next cases: perhaps from atoms point of view?
c              all taken out.

        case default
           write(iunout,*) 'zm_profile: error(1)'
           call eirene_exit_own(1)
       end select

       select case(ipol)
       case(1,-1)
          val = val * (1.d0+ctheta2)/4.d0
        case(0)
          val = val * (1.d0-ctheta2)/2.d0
        case default
           write(iunout,*) 'zm_profile: error(2)'
           call eirene_exit_own(1)
        end select

       res=res+val
      end do  !  loop over 3 normal zeeman components done
      return
      end function EIRENE_zm_profile



      function EIRENE_zm_stark_profile(N,Te,Ti,T_g,B,
     .                          ctheta2,v,e00,e0)
     .                          result(res)
c    .                          npt,omega_min,omega_max)

cdr   npt: option to evaluate function at many energies omega
cdr        in the range omega_min,....,omega_max: removed

!********** DEUTERIUM LYMAN ALPHA LINE SHAPE CALCULATION **********
!03-07-2006
!
!  based upon a routine provided by: J. Rosato, Univ. Marseille
!  see: J. Nucl. Mat., PSI 2006, to appear
!  version: v2 , added: natural broadening
!                added: recommended energy range for full line
!
!******************************************************************
!Calculation and recording of the line shape
!Retained effects :
! -> fine structure
! -> Zeeman effect
! -> Stark broadening (impact approximation for ions and electrons)
! -> Natural broadening (NIST data for hydrogen)
! -> Doppler shift
!On label 123 : setting of the estimated interval of omega
! dr:           use T_g instead of T (=Te=Ti) for Doppler width

c  on input:
c   N: Plasma Density (cm-3), ne=ni
c   Te: Plasma Temperature, Electrons (eV)
c   Ti: Plasma Temperature, Ions (eV)
c   T_g: Neutral Gas Temperature (eV) (only for estimating line width)
c   B: Magnetic field (T)
c   ctheta2: cos**2 of: Observation angle with magnetic field
c   v: Emitter/Absorber velocity (m/s)
c   omega_min: lower bound of interval (eV)
c   omega_max: upper bound of interval (eV)

c   npt: Number of points on line, not in use here
c   if omega_min.eq.omega_max: npt is reset to 1
c                              evaluate line shape at single energy: omega
c   if npt=1 : omega_min is used and omega_max is ignored
c
c  on output:
c   if omega_min.gt.omega_max:
c     omega_min: estimated lower bound of interval (eV)
c     omega_max: estimated upper bound of interval (eV)

      implicit none

!Physical and mathematical constants
      real(dp),parameter::e=1.6022e-19
      real(dp),parameter::hbar=1.0546e-34
      real(dp),parameter::me=9.1094e-31
      real(dp),parameter::alpha=7.2974e-3
      real(dp),parameter::EI=13.606
      real(dp),parameter::c=2.9979e8
      real(dp),parameter::A=6.265e+08  ! Natural broadening added in v2

      real(dp),intent(in)::N,Te,Ti,T_g,B,ctheta2,v
c     real(dp),intent(inout)::omega_min,omega_max
c     integer,intent(inout)::npt
      real(dp)::omega_SF,omega_Z,gamma,gam,epsilon,
     .          omega_plus,omega_minus,
     .          omega1,omega2,omega_D
      real(dp)::C1,C2,C3,C4,C5,C6,C7,C8,omega,line_shape,
     .          res,
     .          e00,e0
c      real(dp) :: delta_omega, interval_omega, omega_D_th
c      integer::i

      omega_SF=alpha*alpha*EI/24.
      omega_Z=hbar*B/(2.*me)
      omega_D=.75*EI*v/c
      gamma=EIRENE_coll(N,Te,Ti,epsilon)+(hbar*A)/e
      omega_plus=.25*omega_SF+.5*omega_Z
      omega_minus=.25*omega_SF-.5*omega_Z
      omega1=.25*sqrt(4.*omega_Z*omega_Z+4.*omega_Z*omega_SF+
     .       9.*omega_SF*omega_SF)
      omega2=.25*sqrt(4.*omega_Z*omega_Z-4.*omega_Z*omega_SF+
     .       9.*omega_SF*omega_SF)
      C1=.5-.5*(omega_Z+omega_minus)/omega1
      C2=.5+.5*(omega_Z+omega_minus)/omega1
      C3=.5+.5*(omega_Z-omega_plus)/omega2
      C4=.5-.5*(omega_Z-omega_plus)/omega2
      C5=.5+.5*omega_plus/omega1
      C6=.5-.5*omega_plus/omega1
      C7=.5+.5*omega_minus/omega2
      C8=.5-.5*omega_minus/omega2

c  estimate energy range
c     if (omega_max.lt.omega_min) then
c       omega_D_th=.75*EI*sqrt(2.*e*T_g/m_D)/c
c       interval_omega=4.*(omega_Z+omega_plus+omega2)+20.*gamma+
c    .                 8.*omega_D_th
c       omega_min=-0.5*interval_omega
c       omega_max=+0.5*interval_omega
c     endif

c123  interval_omega=omega_max-omega_min
c     if (interval_omega.eq.0.0) npt=1
c     delta_omega=0.0
c     if (npt.gt.1) delta_omega=interval_omega/real(npt-1,dp)
c     omega=omega_min

      omega=E0-E00

c     do i=1,npt

c  function Lorentz needs FWHM
        gam=2._dp*gamma
        line_shape=(.5+.5*ctheta2)*
     .  (EIRENE_Lorentz(omega-omega_D-2.*omega_plus,gam)
     .  +EIRENE_Lorentz(omega-omega_D-2.*omega_minus,gam)
     .  +C1*EIRENE_Lorentz(omega-omega_D-omega_Z+omega_minus-omega1,gam)
     .  +C2*EIRENE_Lorentz(omega-omega_D-omega_Z+omega_minus+omega1,gam)
     .  +C3*EIRENE_Lorentz(omega-omega_D+omega_Z+omega_plus-omega2,gam)
     .  +C4*EIRENE_Lorentz(omega-omega_D+omega_Z+omega_plus+omega2,gam))
     .  +(1.-ctheta2)*
     .  (C5*EIRENE_Lorentz(omega-omega_D+omega_plus-omega1,gam)
     .  +C6*EIRENE_Lorentz(omega-omega_D+omega_plus+omega1,gam)
     .  +C7*EIRENE_Lorentz(omega-omega_D+omega_minus-omega2,gam)
     .  +C8*EIRENE_Lorentz(omega-omega_D+omega_minus+omega2,gam))

c       omega=omega+delta_omega
c  factor 0.25 included, because original line shape
c  was normalized to 4 for any fixed ctheta2.
        res=0.25* line_shape
c     end do   ! npt
      return
      end function EIRENE_zm_stark_profile

      function EIRENE_zm_stark_doppler_profile(N,Te,Ti,T_g,B,
     .               ctheta2,v,dnd,drft,e00,e0)
     .               result(res)
c    .               npt,omega_min,omega_max

cdr   npt: option to evaluate function at many energies omega
cdr        in the range omega_min,....,omega_max: removed

!********** DEUTERIUM LYMAN ALPHA LINE SHAPE CALCULATION **********
!03-07-2006
!
!  based upon a routine provided by: J. Rosato, Univ. Marseille
!  see: J. Nucl. Mat., PSI 2006, to appear
!  version: v2 , added: natural broadening
!                added: evaluation of recommended energy range for full line
!
!******************************************************************
!Calculation and recording of the line shape
!Retained effects :
! -> fine structure
! -> Zeeman effect
! -> Stark broadening (impact approximation for ions and electrons)
! -> Natural broadening (NIST data for hydrogen)
! -> Doppler shift
!On label 123 : setting of the estimated interval of omega
! dr:           use T_g instead of T (=Te=Ti) for Doppler width

c  on input:
c   N: Plasma Density (cm-3), ne=ni
c   Te: Plasma Temperature, Electrons (eV)
c   Ti: Plasma Temperature, Ions (eV)
c   T_g: Neutral Gas Temperature (eV) (only for estimating line width)
c   B: Magnetic field (T)
c   ctheta2: cos**2 of: Observation angle with magnetic field
c   v: Emitter/Absorber velocity (m/s)
c   omega_min: lower bound of interval (eV)
c   omega_max: upper bound of interval (eV)
c   if omega_min.eq.omega_max: evaluate line shape at single energy: omega
c  on output:
c   if omega_min.gt.omega_max:
c     omega_min: estimated lower bound of interval (eV)
c     omega_max: estimated upper bound of interval (eV)

      implicit none

!Physical and mathematical constants
      real(dp),parameter::e=1.6022e-19
      real(dp),parameter::hbar=1.0546e-34
      real(dp),parameter::me=9.1094e-31
      real(dp),parameter::alpha=7.2974e-3
      real(dp),parameter::EI=13.606
      real(dp),parameter::A=6.265e+08  ! Natural broadening added in v2

      real(dp),intent(in)::N,Te,Ti,T_g,B,ctheta2
      real(dp),intent(in)::dnd,drft
c     real(dp),intent(inout)::omega_min,omega_max
c     integer,intent(inout)::npt
      real(dp)::omega_SF,omega_Z,gamma,gam,epsilon,
     .          omega_plus,omega_minus,omega1,omega2
      real(dp)::C1,C2,C3,C4,C5,C6,C7,C8,omega,line_shape,res,
     .          v,xx,yy,val,ssum,ssum1,e00,e0,x(-1:8),ci(-1:8)
      integer::i


      omega_SF=alpha*alpha*EI/24.
      omega_Z=hbar*B/(2.*me)
      gamma=EIRENE_coll(N,Te,Ti,epsilon)+(hbar*A)/e
      omega_plus=.25*omega_SF+.5*omega_Z
      omega_minus=.25*omega_SF-.5*omega_Z
      omega1=.25*sqrt(4.*omega_Z*omega_Z+4.*omega_Z*omega_SF+
     .       9.*omega_SF*omega_SF)
      omega2=.25*sqrt(4.*omega_Z*omega_Z-4.*omega_Z*omega_SF+
     .       9.*omega_SF*omega_SF)
      C1=.5-.5*(omega_Z+omega_minus)/omega1
      C2=.5+.5*(omega_Z+omega_minus)/omega1
      C3=.5+.5*(omega_Z-omega_plus)/omega2
      C4=.5-.5*(omega_Z-omega_plus)/omega2
      C5=.5+.5*omega_plus/omega1
      C6=.5-.5*omega_plus/omega1
      C7=.5+.5*omega_minus/omega2
      C8=.5-.5*omega_minus/omega2

c     omega=omega_min
      omega=E0-E00

c  function Lorentz needs FWHM
        gam=2._dp*gamma



c        x=e0-e00
! lorentz
c        val = EIRENE_lorentz(x,gam)
! voigt
c        xx = (x-drft)/dnd
c        yy=gam*0.5_dp/dnd
c        val = REAL(EIRENE_PH_FADDEEVA(xx,yy,dnd), DP)


         x(-1)=omega-2.*omega_plus
         Ci(-1)=1.
         x(0) =omega-2.*omega_minus
         Ci(0)=1.
         x(1) =omega-omega_Z+omega_minus-omega1
         Ci(1)=C1
         x(2) =omega-omega_Z+omega_minus+omega1
         Ci(2)=C2
         x(3) =omega+omega_Z+omega_plus-omega2
         Ci(3)=C3
         x(4) =omega+omega_Z+omega_plus+omega2
         Ci(4)=C4
         x(5) =omega+omega_plus-omega1
         Ci(5)=C5
         x(6) =omega+omega_plus+omega1
         Ci(6)=C6
         x(7) =omega+omega_minus-omega2
         Ci(7)=C7
         x(8) =omega+omega_minus+omega2
         Ci(8)=C8
        ssum=0.
        do i=-1,4
! voigt
          xx = (x(i)-drft)/dnd
          yy=gam*0.5_dp/dnd
          val = REAL(EIRENE_PH_FADDEEVA(xx,yy,dnd), DP)
          ssum=ssum+val*ci(i)
        enddo
        ssum1=0.
        do i=5,8
! voigt
          xx = (x(i)-drft)/dnd
          yy=gam*0.5_dp/dnd
          val = REAL(EIRENE_PH_FADDEEVA(xx,yy,dnd), DP)
          ssum1=ssum1+val*ci(i)
        enddo
        ssum=ssum*(.5+.5*ctheta2)
        ssum1=ssum1*(1.-ctheta2)
        line_shape=ssum+ssum1

c  factor 0.25 included, because original line shape
c  was normalized to 4 for any fixed ctheta2.
        res=0.25*line_shape
c     end do
      return
      end function EIRENE_zm_stark_doppler_profile


!******************************************************************
      function EIRENE_coll(N,Te,Ti,epsilon)
!Returns collision operator matrix element for quantum numbers n=2,l=1
!Ion and electron broadening are retained
      implicit none

!Physical and mathematical constants
      real(dp),parameter::e=1.6022e-19
      real(dp),parameter::m_D=3.3445e-27
      real(dp),parameter::hbar=1.0546e-34
      real(dp),parameter::me=9.1094e-31
      real(dp),parameter::epsilon0=8.8542e-12
      real(dp),parameter::pi=3.1416

      real(dp)::EIRENE_coll,N,Te,Ti,epsilon,v0,ve,rhoWi,rhoWe,
     .          lambda_Di,lambda_De,
     .          ymin_i,ymin_e,phi_i,phi_e
      real(dp)::EIRENE_expint
      external::EIRENE_expint

      v0=sqrt(2.*e*Ti/m_D)
      ve=sqrt(2.*e*Te/me)
      rhoWi=sqrt(6.)*hbar/(me*v0)
      rhoWe=sqrt(6.)*hbar/(me*ve)
      lambda_Di=sqrt(epsilon0*Ti/(N*e*1.e6))
      lambda_De=sqrt(epsilon0*Te/(N*e*1.e6))
      ymin_i=(16./9.)*(rhoWi/lambda_Di)*(rhoWi/lambda_Di)
      ymin_e=(16./9.)*(rhoWe/lambda_De)*(rhoWe/lambda_De)
      phi_i=(hbar/e)*(12.*N*(1.e6)*sqrt(pi)/v0)*(hbar/me)*(hbar/me)*
     .      (3.+EIRENE_expint(ymin_i))
      phi_e=(hbar/e)*(12.*N*(1.e6)*sqrt(pi)/ve)*(hbar/me)*(hbar/me)*
     .      (3.+EIRENE_expint(ymin_e))
      EIRENE_coll=phi_i+phi_e
c  parameter for checking validity of approximation
c     epsilon=(e/hbar)*phi_i*((N*1.e6)**(-1./3.))/v0
      return
      end function EIRENE_coll



      REAL(DP) FUNCTION EIRENE_PH_B12(kk) result(res)
c  photon absorption B12 coefficient, obtained from stim. em. B21 coefficient by detailed balancing
      IMPLICIT NONE
c calculates B12 Einstein coefficient in units: cm**2
! change G1, G2 to REAL, better precision, avoid integer division
      real(dp) :: g1,g2
      integer, intent(in) :: kk

      res=EIRENE_PH_B21(kk)  ! units:

      g1=reacdat(kk)%phr%line%g1
      g2=reacdat(kk)%phr%line%g2

      res=res*g2/g1
      return
      END FUNCTION EIRENE_PH_B12


      REAL(DP) FUNCTION EIRENE_PH_B21(kk) result(res)
c stim. photon emission B21 coefficient, obtained from A21, by detailed balancing
      IMPLICIT NONE
c calculates B21 Einstein coefficient in units: cm**2
      integer, intent(in) :: kk
      real(DP) :: e00

      e00=reacdat(kk)%phr%line%e0

      res=reacdat(kk)%phr%line%aik
      res=res * (hplnk*clight)**2 / (2.*e00**3)
c convert [cm^2 / (eV*s) ] --> [cm^2]
      res=res*hplnk

      return
      END FUNCTION EIRENE_PH_B21



!
      REAL(dp) FUNCTION EIRENE_PH_ENERGY(icell,kk,ipl2,VN,nldoppl)
     .  result(res)
      IMPLICIT NONE
cdr  sample frequency (here: energy) from emission profile
!
!  icell: cell number (needed for parameters in sampling distributions)
!  kk   : process number (for data provided by call get_reaction(kk))
!  ipl2 : species index of emitting atom (ityp=4, bulk)
!
cdr  iprofiletype =0  only line centre (delta),  stationary atoms
cdr  iprofiletype =1  ditto plus doppler in calling program
cdr  iprofiletype =2  only lorentz,  stationary atoms
cdr  iprofiletype =3  ditto plus doppler in calling program
cdr  iprofiletype =4  lorentz+vdw,  stationary atoms
cdr  iprofiletype =5  ditto plus doppler in calling program (new, july 04)
cdr
cdr
cdr
cdr
cdr  july 04: nldoppl introduced
cdr             now ph_energy only samples in the rest frame of the
cdr             emitting atom.
cdr             it returns nldoppl=true, if doppler shift is to be
cdr                                      added in calling program
cdr
cdr  march 05: normal zeeman triplet added.
cdr
cdr  iprofiletype =6  only delta plus zeeman shift,  stationary atoms
cdr  iprofiletype =7  ditto plus doppler in calling program
cdr  iprofiletype =8  lorentz plus zeeman shift,  stationary atoms
cdr  iprofiletype =9  ditto plus doppler in calling program
cdr
cpb  iprofiletype =104 =4 but using line cutoffs
cdr
cdr  iprofiletype =10 zeeman-stark,  stationary atoms
cdr  iprofiletype =11 zeeman-stark, doppler (not in calling program)
cdr
cdr  may 2006: argument VN added, so that doppler (and motional Stark)
cdr            can be included in line-shape 10,11 (zm_stark_profile)
cdr
      integer, intent(in) :: icell,kk,ipl2
      real(dp), intent(in) :: vn
      logical, intent(out) :: nldoppl
      real(dp) :: e00,dnd,drft,gam,fwhm,shift,dvdw,de
c      real(dp) :: dl, l0, l00
      real(dp) :: ctheta,ctheta2,dbz,zep1,v,T_e,T_p,T_g
      real(dp) :: bx, by, bz, bf
      integer :: ictoff
      external :: eirene_exit_own

      ctheta =0._dp
      ctheta2=0._dp

!  idreac: kk from last call to get_reaction
      e00=reacdat(kk)%phr%line%e0

      select case(reacdat(kk)%phr%line%iprofiletype)

!  unless otherwise stated, sampling is in eV

      case(0)
c  delta, no doppler
         res=e00
         nldoppl=.false.
      case(1)
c  delta plus doppler
cdr      call dopplerprof(ipl2,icell,dnd,drft)
cdr      res = sam_doppler(dnd,drft,e00)
         res=e00
         nldoppl=.true.
      case(2)
c  lorentz, no doppler
         call eirene_naturalprof(kk,gam)
         res = EIRENE_sam_lorentz(gam,e00)
         nldoppl=.false.
      case(3)
c  lorentz plus doppler
cdr      call voigtprof(ipl2,icell,dnd,drft,gam)
cdr      res = sam_voigt(gam,dnd,drft,e00)
         call eirene_naturalprof(kk,gam)
         res = EIRENE_sam_lorentz(gam,e00)
         nldoppl=.true.
      case(4)
         ictoff = nreact(kk)
         if ( (ictoff == 0) .or.
     .        (ecutleft(ictoff,icell) > ecutright(ictoff,icell))) then
c  lorentz plus quasistatic vanderWaals, no doppler
c
c  sampling l0 and dl=l0-l00 in wavelength units
c  ie. fwhm, shift and dvdw are given in wavelength units
c
c        call lorvdwprof(icell,fwhm,shift,dvdw,.false.)
c        l00=hpcl/E00
c        l0 = EIRENE_sam_lorentz(fwhm,l00+shift)
c        dl = EIRENE_sam_vdwqs(dvdw,1.E30_dp)
c        l0=l0+dl
c  convert sampled wavelength to energy (eV)
c        res = hpcl/l0
c        nldoppl=.false.
c
c  default: sampling e0 and de=e00-e0 in energy (frequency) units
         call EIRENE_lorvdwprof(kk,icell,fwhm,shift,dvdw,.true.)
         e0 = EIRENE_sam_lorentz(fwhm,e00+shift)
         de = EIRENE_sam_vdwqs(dvdw,e0)
         e0=e0-de
         res=e0
         nldoppl=.false.
         else
c  lorentz plus quasistatic vanderWaals, no doppler, line cutoff
c
c  default: sampling e0 and de=e00-e0 in energy (frequency) units
           e0 = EIRENE_sam_cutoff(kk,ictoff,icell)
           res=e0
           nldoppl=.false.
         end if

      case(5)
c  lorentz plus quasistatic vanderWaals plus doppler
c
c  sampling l0 and dl=l0-l00 in wavelength units
c  ie. fwhm, shift and dvdw are given in wavelength units
c
c        call EIRENE_lorvdwprof(icell,fwhm,shift,dvdw,.false.)
c        l00=hpcl/E00
c        l0 = EIRENE_sam_lorentz(fwhm,l00+shift)
c        dl = EIRENE_sam_vdwqs(dvdw,1.e30_dp)
c        l0=l0+dl
c  convert sampled wavelength to energy (eV)
c        res = hpcl/l0
c        nldoppl=.true.

c  default: sampling e0 and de=e00-e0 in energy (frequency) units
         call EIRENE_lorvdwprof(kk,icell,fwhm,shift,dvdw,.true.)
         e0 = EIRENE_sam_lorentz(fwhm,e00+shift)
         de = EIRENE_sam_vdwqs(dvdw,e0)
         e0=e0-de
         nldoppl=.true.
      case(6)
c  delta + normal zeeman triplet, no doppler
         call EIRENE_zeeman_normalprof(icell,ctheta2,dbz)
         gam=0._dp
         dnd=0._dp
         drft=0._dp
         zep1 = EIRENE_sam_zeeman_normal(ctheta2,dbz,gam,dnd,drft,e00,0)
         res = zep1
         nldoppl=.false.
      case(7)
c  delta + normal zeeman triplet, plus doppler
         call EIRENE_zeeman_normalprof(icell,ctheta2,dbz)
cdr      call dopplerprof(ipl2,icell,dnd,drft)
         gam=0._dp
         dnd=0._dp
         drft=0._dp
         zep1 = EIRENE_sam_zeeman_normal(ctheta2,dbz,gam,dnd,drft,e00,0)
         res = zep1
         nldoppl=.true.
      case(8)
c  lorentz + normal zeeman triplet, no doppler
         call EIRENE_zeeman_normalprof(icell,ctheta2,dbz)
         call eirene_naturalprof(kk,gam)
         dnd=0._dp
         drft=0._dp
         zep1 = EIRENE_sam_zeeman_normal(ctheta2,dbz,gam,dnd,drft,e00,2)
         res = zep1
         nldoppl=.false.
      case(9)
c  lorentz + normal zeeman triplet, plus doppler
         call EIRENE_zeeman_normalprof(icell,ctheta2,dbz)
cdr      call EIRENE_voigtprof(ipl2,icell,dnd,drft,gam)
         call eirene_naturalprof(kk,gam)
         dnd=0._dp
         drft=0._dp
         zep1 = EIRENE_sam_zeeman_normal(ctheta2,dbz,gam,dnd,drft,e00,2)
         res = zep1
         nldoppl=.true.
      case(10)
c  zeeman stark, no doppler, no motional stark
c  stark broadening by electrons:
         T_e=Tein(icell)
c  stark broadening by ions: protons, deuterons, tritons
         T_p=T_e
c  doppler broadening by emitting gas: hydrogen, deuterium, tritium
         T_g=0.0
         v=0.0
         bx = 0._dp
         by = 0._dp
         bz = 1._dp
         bf = 1._dp
         if (lbxin) bx = bxin(icell)
         if (lbyin) by = byin(icell)
         if (lbzin) bz = bzin(icell)
         if (lbfin) Bf = Bfin(icell)
         ctheta=velx*bx+vely*by+velz*bz
         ctheta2=ctheta*ctheta
         zep1 = EIRENE_sam_zm_stark(dein(icell),T_e,T_p,T_g,Bf,
     .                              ctheta2,E00,v)
         res = zep1
         nldoppl=.false.
      case(11)
c  zeeman stark, doppler and motional stark is included in sampling
         T_e=Tein(icell)
c  stark broadening by ions: protons, deuterons, tritons
         T_p=T_e
c  doppler broadening by emitting gas: hydrogen, deuterium, tritium
         T_g=TIIN(mplsti(ipl2),icell)
         v=vn/100._dp
         bx = 0._dp
         by = 0._dp
         bz = 1._dp
         bf = 1._dp
         if (lbxin) bx = bxin(icell)
         if (lbyin) by = byin(icell)
         if (lbzin) bz = bzin(icell)
         if (lbfin) Bf = Bfin(icell)
         ctheta=velx*bx+vely*by+velz*bz
         ctheta2=ctheta*ctheta
         zep1 = EIRENE_sam_zm_stark(dein(icell),T_e,T_p,T_g,Bf,
     .                              ctheta2,E00,v)
         res = zep1
         nldoppl=.false.

      case default
        write (iunout,*)
     .  'profiletype in ph_energy ? exit called EIRENE_'
        call EIRENE_exit_own(1)
      end select

      return
      END FUNCTION EIRENE_PH_ENERGY
c


c  this next routine ph_post_energy is currently not in use. It is currently
c  developed for photon re-emission (scattering) during photon tracing.

c POST-COLLISION
      SUBROUTINE EIRENE_PH_POST_ENERGY(icell,kk,iflg,il,
     .                          iold,itypold,vxo,vyo,vzo,vlo,e0o,
     .                          itypnew)
c sample post-collision energy.
c
c incident particle: (iold,itypold,vxo,....e0o)
c already decided: new test particle has type itypnew
c
c called from colatm (itypold=1), or called from colphot (itypold=0)
!  6.7.05: itypold=1 option removed from eirene. now moved to photon_sven.ff
c
c iflg = 0: sp.emission, no bulk pre-collision partner
c        1: absorption
c        2: stim.em
      IMPLICIT NONE
      integer, intent(in) :: icell,kk,iflg,iold,il,itypold,itypnew
      real(dp),intent(in) :: vxo,vyo,vzo,vlo,e0o
      integer :: nrc,ipln,itypn,
     .    ityp0,ityp1,ityp2,ipl0,ipl1,ipl2,
     .    ipl0v
      real(dp) :: vx,vy,vz,vxn,vyn,vzn,cvrss1,velq,e1,e2,e00,l00,gam,
     .    zep1,cangl,een,dl,
     .    fwhm,shift,dvdw,ee,cvel,l0,e00s,
     .    velx_b, vely_b, velz_b, velparm, vel_b
      external :: eirene_fisotr, eirene_fgauss, eirene_exit_own

c      if (idreac /= kk) call EIRENE_get_reaction(kk)

      select case(itypold)
      case(0)
         do nrc=1,nrcph(iold)
            if(kk == ireacph(iold,nrc)) exit
         enddo
         IPL0 =eirene_IDEZ(IBULKPH(IOLD,nrc),3,3)
         IPL1 =eirene_IDEZ(ISCD1PH(IOLD,nrc),3,3)
         IPL2 =eirene_IDEZ(ISCD2PH(IOLD,nrc),3,3)
         ITYP0=eirene_IDEZ(IBULKPH(IOLD,nrc),1,3)
         ITYP1=eirene_IDEZ(ISCD1PH(IOLD,nrc),1,3)
         ITYP2=eirene_IDEZ(ISCD2PH(IOLD,nrc),1,3)
      case default
         nrc=-1
      end select

      ipl0v = mplsv(ipl0)
c
c
      select case(iflg)
      case(0)
c     spont. emission, no bulk particle as collision partner
c
c  this part: itypold=1 : removed
c
c        if (itypold==1) then  .....
c        elseif (itypold==0) then
c     to be written
c        endif

         if(ityp1 == itypnew) then
            itypn = ityp1
            ipln  = ipl1
         elseif(ityp2 == itypnew) then
            itypn = ityp2
            ipln  = ipl2
         else
            write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
            write(iunout,*) '   ityp of secondary wrong'
            write(iunout,*) '   (sp.emission)'
            call EIRENE_exit_own(1)
         endif

         secs: select case(itypn)
         case(0)
c        photon followed:
c
C  SAMPLE ISOTROPIC EMISSION OF PHOTON IN REST FRAME OF EMITTING PARTICLE
            VEL=CLIGHT
            IF (INIV3.EQ.0) CALL EIRENE_FISOTR
            VELX=FI1(INIV3)
            VELY=FI2(INIV3)
            VELZ=FI3(INIV3)


            IF (ITYPOLD == 0) THEN
c  sample velocity from ipl0 bulk, shifted maxw
c  this case has been added (cdr, aug.23, 2004)
c  "complete re-distribution",
c  for doppler here first sample a background velocity of bulk atom
c  (similar to subr. locate)
               IF(INIV2.EQ.0) CALL EIRENE_FGAUSS
               VXN=FG1(INIV2)
               VYN=FG2(INIV2)
               VZN=FG3(INIV2)
               INIV2=INIV2-1

               VELPARM=ZRG(ipl0,icell)

               vx=VELPARM*VXN+VXIN(ipl0v,icell)
               vy=VELPARM*VYN+VYIN(ipl0v,icell)
               vz=VELPARM*VZN+VZIN(ipl0v,icell)
               velq=vx*vx+vy*vy+vz*vz
               vel_b=dsqrt(velq)

               VELX_B=VX/vel_b
               VELY_B=VY/vel_b
               VELZ_B=VZ/vel_b
               cangl = (velx*velx_b + vely*vely_b + velz*velz_b)
               cvel  = vel_b*cangl
            ELSE
c  this was originally for itypold=1,  now removed
            END IF

c        sample new freq, doppler shift due n=2 velocity
            e00=reacdat(kk)%phr%line%e0
            e1=reacdat(kk)%phr%line%e1
            e2=e00+e1

            labelp1: select case(reacdat(kk)%phr%line%iprofiletype)
            case(0)
               ee = e00
            case(1)
               ee = e00*(1.-cvel/clight)
            case(2)
               call eirene_naturalprof(kk,gam)
               zep1=EIRENE_sam_lorentz(gam,e00)
               ee = zep1
            case(3)
               call eirene_naturalprof(kk,gam)
               een = e00*(1.-cvel/clight)
               zep1=EIRENE_sam_lorentz(gam,een)
               ee = zep1
            case(4)
c constants in wavelengths
               call EIRENE_lorvdwprof(kk,icell, fwhm,shift,dvdw,.false.)
               l00 = hpcl/e00
               l0  = EIRENE_sam_lorentz(fwhm,shift+l00)
               dl  = EIRENE_sam_vdwqs(dvdw,1.e30_dp)
               l0=l0+dl
               ee = hpcl/l0
            case(5)
c constants in wavelengths
               call EIRENE_lorvdwprof(kk,icell, fwhm,shift,dvdw,.false.)
               l00 = hpcl/e00
               l0  = EIRENE_sam_lorentz(fwhm,shift+l00)
               dl  = EIRENE_sam_vdwqs(dvdw,1.e30_dp)
               l0=l0+dl
               ee  = hpcl/l0
               ee  = ee*(1.-cvel/clight)
            case(6)
               call EIRENE_strkprof(icell,fwhm,shift)
               call eirene_naturalprof(kk,gam)
               gam=gam+fwhm

               e00s = e00+shift
               een = e00s
               zep1=EIRENE_sam_lorentz(gam,een)
               ee = zep1
            case(7)
               call EIRENE_strkprof(icell,fwhm,shift)
               call eirene_naturalprof(kk,gam)
               gam=gam+fwhm

               e00s = e00+shift
               een = e00s*(1.-cvel/clight)
               zep1=EIRENE_sam_lorentz(gam,een)
               ee = zep1
            end select labelp1
            e0=ee
            return

         case(1)
c     n=1 followed
c     use old values from n=2 testpart.
            velx=vxo
            vely=vyo
            velz=vzo
            vel=vlo
            e0=e0o
            return
         case(4)
c     bulk particle (photon or atom)
            vel=0.
            e0=0.
            return
         case default
            write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
            write(iunout,*) '   default case?!, check input!'
            call EIRENE_exit_own(1)
         end select secs
c
c  spontane emission , for incident atom, finished
c
      case(1)
c     absorption process
         if(ipl1 /= 0 .and. ipl2 /= 0) then
            write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
            write(iunout,*) '       only 2 secondaries allowed in abs.'
            write(iunout,*) '       process, check input file'
            write(iunout,*) '       ph + n=1 -> n=2'
            call EIRENE_exit_own(1)
         endif
         if(ityp1 == itypnew) then
            itypn = ityp1
            ipln  = ipl1
         elseif(ityp2 == itypnew) then
            itypn = ityp2
            ipln  = ipl2
         else
            write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
            write(iunout,*) '   ityp of secondary wrong'
            write(iunout,*) '   (absorption)'
            call EIRENE_exit_own(1)
         endif

         abcs: select case(itypold)
         case(0)
c        photon on n=1 background, n=2 followed
            abcs2: select case(itypn)
            case(1)
c           followed n=2 is atom
c           sample new energy from ipl0 bulk, shifted maxw
               IF(INIV2.EQ.0) CALL EIRENE_FGAUSS
               VXN=FG1(INIV2)
               VYN=FG2(INIV2)
               VZN=FG3(INIV2)
               INIV2=INIV2-1

               cvrss1=cvrssa(ipln)
               VEL=ZRG(ipl0,icell)

c               vx=VEL*VXN+VXIN(ipl0v,icell)
c               vy=VEL*VYN+VYIN(ipl0v,icell)
c               vz=VEL*VZN+VZIN(ipl0v,icell)
               vx=VEL*VXN
               vy=VEL*VYN
               vz=VEL*VZN
               velq=vx*vx+vy*vy+vz*vz
               vel=dsqrt(velq)

               VELX=VX/vel
               VELY=VY/vel
               VELZ=VZ/vel
               E0=CVRSS1*VELQ
               return
            case(4)
c               write(iunout,*) 'PHOTON MODULE (PH_POST_ENERGY):'
c               write(iunout,*) '       secondary in abs.process is bulk'
c               write(iunout,*) '       somethings wrong in input file'
c               write(iunout,*) '       ph + n=1(bulk) --> n=2(bulk)'
c               call EIRENE_exit_own(1)
               e0=0.
               vel=0.
               return
            end select abcs2

         case(1)
c        n=1 on photon background, n=2 followed,  removed in july 05
         case default
            write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
            write(iunout,*) '   default case?!, check input!'
            call EIRENE_exit_own(1)
         end select abcs

      case(2)
c     stimulated process
         if(ipl1 == 0 .or. ipl2 == 0) then
            write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
            write(iunout,*) '       two secondaries are needed in ',
     .                              'stim.em.'
            write(iunout,*) '       check input file'
            write(iunout,*) '       (ph + n=2) -> n=1 + 2ph'
            call EIRENE_exit_own(1)
         endif
         if(ityp1 == itypnew) then
            itypn = ityp1
            ipln  = ipl1
         elseif(ityp2 == itypnew) then
            itypn = ityp2
            ipln  = ipl2
         else
            write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
            write(iunout,*) '   ityp of secondary wrong'
            write(iunout,*) '   (stim.emission)'
            call EIRENE_exit_own(1)
         endif

         stcs: select case(itypold)
         case(0)
c        photon on n=2 background
            stcs2: select case(itypn)

            case(0)
c           two photons followed
c           use old energy from photon(ipl0)
            e0=e0o
            velx=vxo
            vely=vyo
            velz=vzo
            vel =vlo
            velq=vlo*vlo
            return

            case(1)
c           n=1 followed
c           sample new energy from ipl0 bulk, shifted maxwellian
               IF(INIV2.EQ.0) CALL EIRENE_FGAUSS
               VXN=FG1(INIV2)
               VYN=FG2(INIV2)
               VZN=FG3(INIV2)
               INIV2=INIV2-1

               cvrss1=cvrssa(ipln)
               VEL=ZRG(ipl0,icell)

c               vx=VEL*VXN+VXIN(ipl0v,icell)
c               vy=VEL*VYN+VYIN(ipl0v,icell)
c               vz=VEL*VZN+VZIN(ipl0v,icell)
               vx=VEL*VXN
               vy=VEL*VYN
               vz=VEL*VZN
               velq=vx*vx+vy*vy+vz*vz
               vel=dsqrt(velq)

               VELX=VX/vel
               VELY=VY/vel
               VELZ=VZ/vel
               E0=CVRSS1*VELQ
               return
c     bulk particle
            case(4)
               e0=0.
               vel=0.
               return
            case default
               write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
               write(iunout,*) '   default case?!, check input!'
               call EIRENE_exit_own(1)
            end select stcs2

         case(1)
c        n=2 on photon background , removed in july 05
         case default
            write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
            write(iunout,*) '   default case?!, check input!'
            call EIRENE_exit_own(1)
         end select stcs
      end select

      write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_POST_ENERGY):'
      write(iunout,*) '   end of subroutine reached'
      write(iunout,*) '   there has to be an error, check input!'
      call EIRENE_exit_own(1)
      return
      END SUBROUTINE EIRENE_PH_POST_ENERGY


c ROUTINES FOR EVALUATING SOME PARAMETERS OF LINE SHAPE FUNCTIONS:
C  LORENTZ * VD WAALS      = eirene_LORVDWprof
C  LORENTZ * DOPPLER       = eirene_VOIGTprof
C  DOPPLER                 = eirene_Dopplerprof
C  NATURAL LINE BROADENING = NATURALprof
C  STARK FOR LYMAN ALPHA   = STRKprof
C  NORMAL-ZEEMANN

      SUBROUTINE EIRENE_LORVDWPROF(kk,icell,fwhm,shift,dvdw,lscale)
!
!  set parameters for convoluted Lorentz / vdWaals profile
!
!  input: icell:  cell number (for background data for broadening parameters)
!         lscale: true : output parameters in energy units, eV
!         lscale: false: output parameters in wavelength units, cm
!
      implicit none
      integer, intent(in) :: kk,icell
      logical, intent(in) :: lscale
      real(dp),intent(out) :: fwhm,shift,dvdw
      integer :: ipl, ip, iplti
      real(dp) :: nipl,ne,e00,dlst,dlqs,dlvw,te,tipl,dlres,l00,mipl
      real(dp) :: cvel2a_m
      real(dp) :: dlvws, dlqss
      fwhm=0.
      shift=0.
      dvdw=0.

      ipl=reacdat(kk)%phr%line%ignd
      e00=reacdat(kk)%phr%line%e0
      l00=hpcl/e00

c densities [m^-3], temperatures [eV]
c nipl  and tipl for vdWaals broadening with other ipl will be set later
      iplti=mplsti(ipl)
      nipl=DIIN(ipl,icell)*1.e6_dp
      tipl=tiin(iplti,icell)
      ne  =DEIN(icell)*1.e6_dp
      te  =tein(icell)

c mass ipl
      mipl=rmassp(ipl)

!  convert cvel2a from cm/s  to m/s
      cvel2a_m=cvel2a*1.E-2_dp

!  FWHM, resonance broadening  [cm]
      dlres= l00*l00*PIA/clight*reacdat(kk)%phr%line%c3*nipl
!  FWHM, quadratic stark broadening [cm]
      dlst=(sqrt(8./PIA*te/pmasse)*cvel2a_m)**(1./3.)
      dlst=dlst*l00*l00/PI2A/clight*11.37*
     .     reacdat(kk)%phr%line%c4**(2./3.)
      dlst=dlst*ne
!  for vdW, impact and quasistatic approx.: sum over up to 10 foreign gas species
      dlvws = 0._dp
      dlqss = 0._dp
      do ip=1,reacdat(kk)%phr%line%ifremd
        ipl=reacdat(kk)%phr%line%iplsc6(ip)
        if (ipl > 0) then
          iplti=mplsti(ipl)
          nipl=DIIN(ipl,icell)*1.e6_dp
          tipl=tiin(iplti,icell)
          mipl=rmassp(ipl)
!  FWHM  vdWaals, impact approximation [cm]
          dlvw=(sqrt(16._dp/PIA*tipl/mipl)*
     .          cvel2a_m)**(0.6_dp)
          dlvw=dlvw*l00*l00/PI2A/CLIGHT*8.08_dp
     .         *reacdat(kk)%phr%line%c6a(ip)**(0.4_dp)
          dlvw=dlvw*nipl
          dlvws = dlvws + dlvw
!  Parameter for vdWaals, quasi-static approximation [cm]
          dlqs = l00*l00/PI2A/clight*
     .            (4._dp/3._dp*PIA)**2._dp*reacdat(kk)%phr%line%c6a(ip)
          dlqs = dlqs * nipl*nipl
          dlqss = dlqss + dlqs
        end if
      end do

!  Shift [cm]
!  factor 0.5, because this formula for shift applies to HWHM, not FWHM
      shift=0.5*(-1.73*dlst -0.7*dlvws)

!  FWHM for Lorentzian
      fwhm   =dlres + dlst + dlvws
!  Parameter for Exponential
      dvdw =dlqss

c transform wavelength [cm] --> energy [eV]
c dE / dlambda = -hc*lambda^-2 = -E^2/(hc)
      if (lscale) then
        shift  =-shift*e00*e00/hpcl
        fwhm   = fwhm *e00*e00/hpcl
        dvdw   = dvdw *e00*e00/hpcl
      end if

      return
      end subroutine EIRENE_lorvdwprof

      subroutine EIRENE_voigtprof(ipl,kk,icell,dnd,drft,gam)
!  return gam, the FWHM of the Lorentz profile (natural broadening)
!  return dnd, the doppler width (half width at 1/e maximum)
!              dnd/sqrt(2) is the std. deviation of a gaussian
!  return drft, the line shift due to the drift motion in the maxwellian
      implicit none
      integer, intent(in) :: ipl,kk,icell
      real(dp), intent(out) :: dnd,drft,gam

      call EIRENE_dopplerprof(ipl,kk,icell,dnd,drft)
      call eirene_naturalprof(kk,gam)
      return
      end subroutine EIRENE_voigtprof

      subroutine EIRENE_DOPPLERPROF(iipl,kk,icell,dnd,drft)
!  returns the doppler width dnd of a doppler broadened line profile
!     and the shift drft due to a drift in the maxwellian
!  use energy units, eV
!  This profile is a gaussian profile with sig=dnd/sqrt(2) as standard deviation.
!  dnd is defined as Doppler half width at 1/e of maximum.
!  dnd*sqrt(4*ln(2)) is the FWHM  of the doppler  profile.
!
!  dnd = E00/clight*sqrt(2*Temp/Mass)
!                    with sqrt(2*Temp/Mass) in the same units as clight (cm/s)
!  drft= E00* (VEC_V_PHOTON * VEC_V_DRIFT)/CLIGHT^2
!
      implicit none
      integer, intent(in) :: iipl,kk,icell
      real(dp), intent(out) :: dnd,drft
      real(dp) :: e00,t
      integer :: ipl,iplti,iplv

      ipl=iipl
      e00=reacdat(kk)%phr%line%e0
      iplti=mplsti(ipl)
      t=tiin(iplti,icell)
      if (t .le. 0.0) then
        write (iunout,*) 'zero temp. encountered in dopplerprof.'
        write (iunout,*) 'ipl, t ',ipl,t
        call eirene_exit_own(1)
      endif
      dnd=e00*sqrt(t)*rsqdvp(ipl)/clight
      drft=0._dp
      if (nldrft) then
        iplv =mplsv(ipl)
        drft=velx*vxin(iplv,icell)+vely*vyin(iplv,icell)+
     .       velz*vzin(iplv,icell)
        drft=drft*e00/clight
      endif
      return
      end subroutine EIRENE_DOPPLERPROF

      subroutine EIRENE_NATURALPROF(kk,gam)
!  this routines returnes the FWHM parameter (gam)
!  for the Lorentzian Natural line broadening profile
!  still to be done: add other Aik's from upper and lower level of line
!  still to be done: precomputation of scaling factors
      IMPLICIT NONE
      integer, intent(in) :: kk
      real(dp), intent(out) :: gam

      gam=0.
c  NATURAL LINE BROADENING, Radians/s
      gam=gam+reacdat(kk)%phr%line%aik
c  convert to frequency (Hz), and then to energy, eV
      gam=gam*hplnk_bar
      return
      END subroutine EIRENE_NATURALPROF


      subroutine EIRENE_STRKPROF(icell, fwhm,shift)
!
! return FWHM for linear electron stark broadening (energy units, eV)
!  output:      fwhm is the FWHM,
!               shift is the line-shift.
!  currently formula from Sobelman, Vainshtein, for Lyman alpha only.
!
      IMPLICIT NONE
      integer, intent(in) :: icell
      real(dp), intent(out) :: fwhm,shift
      real(dp) :: de,te,nn,
     .            lc,ve,inn,n1
c      real(dp) :: ec, eh, fac, nn2, nn3, nn4, part1, part2, part3, sqeh,
c     .            w, z
      real(dp) :: rd,rw,hw0

c to be written:
c stark-effect: for now only hydrogen treated
c               linear stark effect, and electron contribution only
c               Lyman alpha only, see fixed n1 and nn below.

      de = dein(icell)
      te = tein(icell)
c     lower state
      n1  = 1._dp
c     upper state
      nn  = 2._dp
c compton length h^bar/m_e/c [cm]
      lc = 3.8616e-11_dp

c
c   FWHM from Sobelman, Vainshtein, Yukov, eq. 7.3.35

c  eq.7.3.36
      Inn=n1**4+nn**4

c  debye radius (cm), mean electron velocity
      rd = 7.43e2_dp*dsqrt(te/de)
      ve = dsqrt(8._dp/PIA)*4.19e7_dp*dsqrt(te)
c  weisskopf radius, (cm) eq. 7.3.33,  use for Bohr radius: a_0=h^bar^2/(m_e*e^2)
      rw = dsqrt(2._dp/3._dp*Inn)/ve*lc*clight

      hw0= 32._dp/3._dp*de/ve*lc*lc*clight*clight
     .                 * (dlog(rd/rw)+.215_dp)*Inn
c  hw0 (=gamma) is FWHM, see definition in Sobelman, Lorentzian, eq. 7.1.18
c  convert w (rad/s) to frequency (Hz), and then to energy (eV)
      hw0=hw0*hplnk_bar
csw
csw stark effect (Griem), not in use.
csw
c     eh = 13.606_dp
c     z  = 1._dp
c     nn2 = nn*nn
c     nn3 = nn*nn*nn
c     nn4 = nn2*nn2
c electron plasma frequency [eV]
c     ec =8.98e3_dp*dsqrt(de)*hplnk
c     sqeh=dsqrt(eh)

c     part1 = te**1.5/((z-1._dp)*sqeh*ec + te*nn2*ec/z/sqeh)
c     part1 = dlog(part1) *(3._dp*nn4-9._dp*nn2)

c     part2 = nn4 + nn4*2._dp*te/eh/(1._dp+2._dp*te/eh+z*z/nn4)

c     part3 = nn3/((z-1._dp)*z*z+te*nn2*z/eh)*(te/eh)**1.5+1.4_dp
c     part3 = dlog(part3) *(nn4/3._dp + 17._dp*nn2/3._dp)

c mean velocity of electron sqrt(kT/me) [cm/s]
c     ve = 4.19e7_dp*dsqrt(te)

c     fac = dsqrt(2._dp*PIA)/z/z*lc*lc*clight*clight
c     fac = fac / ve * de

c     w  = hplnk_bar * fac * (part1 + part2 + part3)
c brd.const, shift [eV]
c      hw = w
c  FWHM-Griem, done,  not in use.
!
      fwhm = hw0
      shift = fwhm*dsqrt(3._dp)/2._dp

      return
      end subroutine EIRENE_strkprof


      subroutine EIRENE_zeeman_normalprof(icell,ctheta2,dbz)
!
!  return the angle theta of the photon relative to B
!  and the zeeman splitting paramater dbz for
!  the normal zeeman triplet
!  return only cos(theta)^2, rather than theta itself
!
      IMPLICIT NONE
      integer, intent(in) :: icell
      real(dp), intent(out) :: ctheta2,dbz
      real(dp)  :: ctheta
      real(dp)  :: bx, by, bz, bf

      bx = 0._dp
      by = 0._dp
      bz = 1._dp
      bf = 1._dp
      if (lbxin) bx = bxin(icell)
      if (lbyin) by = byin(icell)
      if (lbzin) bz = bzin(icell)
      if (lbfin) Bf = Bfin(icell)
      ctheta=velx*bx+vely*by+velz*bz
      ctheta2=ctheta*ctheta
c
      dbz=bf*mub

      return
      END subroutine EIRENE_zeeman_normalprof

C.........................................................................................
c  evaluate line shape functions at fixed frequency (or wavelength or energy)
c  1) lorentz                          --> evaluates lorentz profile
c  2) doppler                          --> evaluates doppler profile
c  3) faddeeva  (via humlick function) --> evaluates lorentz-doppler convol. (Voigt profile)
c  faddeeva2
c  lorvdw                              --> evaluates lorentz-vdWalls convolution


      REAL(dp) FUNCTION EIRENE_LORENTZ(xx,gm) result(res)
!  evaluate lorentzian profile at x, with xx=x-shift,
!  gm  is FWHM, xx already includes the shift
!  gmh is HWHM
!  special case: Cauchy distribution, for gm=2, shift=0

      IMPLICIT NONE
      real(dp), intent(in) :: xx,gm
      real(dp) :: gmh

      gmh=gm*0.5
      res = gmh/PIA/(xx*xx + gmh*gmh)
      return
      END FUNCTION EIRENE_LORENTZ

      REAL(dp) FUNCTION EIRENE_DOPPLER(xx,dnd) result(res)
!
!  return value of a doppler profile at xx
!  dnd is the doppler width
!  i.e. this value is from a gaussian, with std. dev. sig= dnd/sqrt(2)
!
      IMPLICIT NONE
      real(dp), intent(in) :: xx,dnd
      res=exp(-(xx/dnd)**2)/(dnd*PISQ)
      return
      END FUNCTION EIRENE_DOPPLER

      real(dp) function EIRENE_planck(E, T, b_nu, imode) result(res)
c  for testing purposes:
c  a) evaluate planck function B_nu(T) for radiation intensity at freq. nu
c  b) evaluate planck function B_E (T) for radiation intensity at energ. E
c
c  i.e.: use energy scale instead of frequency scale
c  B_E = 1/h_planck B_nu(T), with E = h_planck * nu
c  input : E and T in eV
c  output: B_E(T) in 1/cm**3/eV/sterad * cm/s * eV
c          i.e.   in 1/cm**2/s/sterad

cdr  imode: =1 Planck
cdr         =2 Wien limit E >> kT
cdr         =3 Rayleigh-Jeans limit E << kT

      implicit none
      real(dp), intent(in) :: t,e
      integer , intent(in) :: imode
      real(dp), intent(out) :: b_nu

      res=0._dp
      select case(imode)
      case(1)
! planck
        res = 2._dp*e**3
        res=res / (hplnk**2*clight**2)
        res=res / (dexp(e/t)-1._dp)
        b_nu=res
        res=res/hplnk
      case(2)
! wien, e >> kt
        res = 2._dp*e**3
        res = res / (hplnk**3*clight**2)
        res=  res*dexp(-e/t)
      case(3)
! rayleigh-jeans, e << kt
        res = 2*e**2*t/(hplnk**3*clight**2)
      end select

      return
      end function EIRENE_planck

      complex(dp) function EIRENE_ph_faddeeva(x,y,dnd) result(cres)
      implicit none
      real(dp), intent(in) :: x,y,dnd

      cres=EIRENE_ph_humlik(x,y)/(dnd*sqrt(pia))
      return
      end function EIRENE_ph_faddeeva

      real(dp) function EIRENE_ph_lorvdw(dx,fwhm,shift,dvdw,
     .                            icell) result(res)
!  evaluate convolution integral of Lorentzian and Exponential (Stormberg)
!  dx=l0-l00, if wavelength units are used
!  dx=e00-e0, if energy (frequency) units are used
!  fwhm is the FWHM of the Lorentzian
!  shift is the shift of the Lorentzian
!  dvdw is the parameter in the Exponential
      implicit none
      real(dp), intent(in) :: dx,fwhm,shift,dvdw
      real(dp) :: a,b,cpi,rvdw,aa,hw
      complex(dp) :: z1,zz,sz1,fad,z
      integer, intent(in) :: icell

c  this routine works with the hwhm, rather than with the fwhm
      hw=0.5_dp*fwhm

c  lorentzian part
      a=(dx-shift)/hw
      aa  = 1._dp/(1._dp+a*a)
      res = aa/(hw*PIA)

c  next: part with imag. fadeeva function
      b   = PIQU*dvdw/hw
      cpi = dsqrt(dvdw/hw) / fwhm

! kind parameter added for better precision
      z1  = cmplx(-a, -1._dp,dp) * aa
!pb      zz  = cdsqrt(b*z1) *(0._dp,1._dp)
      zz  = sqrt(b*z1) *(0._dp,1._dp)
#ifdef LEGACYCOMP
      fad = EIRENE_ph_faddeeva2(abs(zz), aimag(zz), icell, 0)
#else
      fad = EIRENE_ph_faddeeva2(zz%RE, zz%IM, icell, 0)
#endif
!pb      sz1 = z1*cdsqrt(z1)
      sz1 = z1*sqrt(z1)
      z=sz1*fad
      rvdw = real(AIMAG(z),DP)*cpi

      res=res+rvdw

      return

      CONTAINS

      complex(dp) function EIRENE_ph_faddeeva2(x,y,icell,ipl)
     .  result(cres)
      implicit none
      real(dp), intent(in) :: x,y
      integer, intent(in) :: icell,ipl
      real(dp) :: u,v
      logical :: flag
      external :: eirene_exit_own

      call EIRENE_ph_wofz(x,y,u,v,flag)
      if(flag) then
         write(iunout,*) 'PHOTON MODULE EIRMOD_(PH_FADDEEVA2):'
         write(iunout,*) '   flag=',flag
         call EIRENE_exit_own(1)
      endif
! kind parameter added for better precision
      cres=cmplx(u,v,dp)
      return
      end function EIRENE_ph_faddeeva2

      end function EIRENE_ph_lorvdw

c...............................................................................

c RANDOM SAMPLING x  (frequency, energy or wavelength) from line shape routines
! 020205  comments, cleaned up, otherwise identical to photon.f in eirene_04
c  0) delta          (nothing to be done)
c  1)         lorentz           (Reiter)
c  2)         doppler           (Reiter)
c  6,7,8,9)   normal zeemann    (Reiter)
c  zeemann-stark1               (Rosato)
c  10,11)     zeemann-stark     (Rosato)
c  3) voigt                     (Reiter)
c  vdwqs                        (Reiter)

      REAL(dp) FUNCTION EIRENE_SAM_LORENTZ(alph,shift) result(res)
!  sample from a Lorentzian line profile.
!  alph is the FWHM, shift is the shift.
!  alphh = alph/2 is the HWHM (half width half maximum)
!  the distribution is g(x)=1/pi*[alphh/{(x-shift)**2+alphh**2}]
!  sampling is first from a standardized Lorentzian (alph=2, shift=0)
!  i.e., from a Cauchy-distribution (x),
!  then transform: xx=x*alphh+shift

      IMPLICIT NONE
      real(dp), intent(in) :: alph, shift
      real(dp) :: rr, x, xx, alphh

      alphh=alph*0.5
      do
        rr=PIHA*(RANF_EIRENE()*2._DP-1._DP)
        x=tan(rr)
!  x is sampled from Cauchy. Now transform to Lorentz(gam,shift)
        xx=alphh*x+shift
        res=xx
        if (res > 0._dp) exit
      end do
      return
      END FUNCTION EIRENE_SAM_LORENTZ

      REAL(dp) FUNCTION EIRENE_SAM_DOPPLER(dnd,drft,e00) result(res)
!
!  sample from a central doppler profile
!  dnd is the doppler width. drft is drift contribution in case of drifting maxw.
!  dnd*sqrt(4.*log(2.)) is the FWHM
!  i.e., sample from a Gaussian with standard deviation sig=dnd/sqrt(2)
!
      IMPLICIT NONE
      real(dp), intent(in) :: dnd,drft,e00
      real(dp) :: v1,v2,s,ar,f1,sig

      sig=dnd/sqrt(2._dp)
c  now sample from a gaussian with standard deviation sig
cdr  1st step: uniform on a unit circle, by rejection
      do
         v1=2.*RANF_EIRENE()-1.
         v2=2.*RANF_EIRENE()-1.
         s=v1*v1+v2*v2
         if(s < 1.) exit
      enddo
cdr  2nd step: convert to standard gaussian
      ar=log(s)
      f1=v1*sqrt(-(ar+ar)/s)*sig
cdr  we could save this v2, ar, and s, for next call 
c     f2=v2*sqrt(-(ar+ar)/s)*sig
      res=f1+drft+e00
      return
      END FUNCTION EIRENE_SAM_DOPPLER

      function EIRENE_sam_zeeman_normal
     . (ctheta2,dbz,gam,dnd,drft,e00,iprof) result(res)
!
!  sample energy (eV) from a normal zeeman triplet, each component either
!     iprof=0: delta
!     iprof=1: doppler
!     iprof=2: lorentz
!     iprof=3: voigt
!
      implicit none
      real(dp), intent(in) :: ctheta2,dbz,gam,dnd,drft,e00
      integer, intent(in) :: iprof
      real(dp) :: dsum,e00d,del,r0,ssum,res
      real(dp) :: strength(-1:1)
      integer :: ipol
      external :: eirene_exit_own

      res = 0._dp
c  3 normal zeeman components:
c   build sampling distribution for sampling the component, given theta
      strength(-1) = (1._dp+ctheta2)*0.5_dp
      strength(0)  =  1._dp-ctheta2  ! = sin(theta)^2
      strength(1)  =  strength(-1)
      dsum=2._dp

      r0 = ranf_eirene()*dsum

      ssum = 0.
      do ipol = -1,1
        ssum = ssum + strength(ipol)
        if(r0 <= ssum) then
          del = dbz*real(ipol,DP)
          e00d = e00 + del
          select case(iprof)
          case(0)
             res = e00d
          case(1)
             res = EIRENE_sam_doppler(dnd,drft,e00d)
          case(2)
             res = EIRENE_sam_lorentz(gam, e00d)
          case(3)
             res = EIRENE_sam_voigt(gam,dnd,drft,e00d)
          case default
            write (iunout,*) 'exit from sam_zeeman_normal, iprof ',iprof
            call EIRENE_exit_own(1)
          end select
          return
       endif
      enddo
      return
      end function EIRENE_sam_zeeman_normal


      FUNCTION EIRENE_sam_zm_stark1(N,Te,Ti,T_g,B,
     .                       ctheta2,E00,v)
     .                       result(res)


c  this is basically "function zm_stark_profile", modified for
c  sampling from rather than evaluating of the line profile.
c  v is the velocity of the emitting atom.
c  npt set internally to 100.

c  This routine has been replaced by new (faster) version of sam_zm_stark
c  (see function sam_zm_stark, below).
c  This first version (this routine sam_zm_stark1)
c  is based on a numerical integration of the pdf into
c  a cumulative distribution and then "sampling by inversion"

c  the original line_shape (sum of weighted Lorentzians)
c  is normalized to 4.0, for each given incident angle theta.
c  here: normalized to 1.
c  to be done: speed up of evaluation of sum of Lorentzians

!********** DEUTERIUM LYMAN ALPHA LINE SHAPE CALCULATION **********
!******************************************************************
!Calculation and recording of the line shape
!Retained effects :
! -> fine structure
! -> Zeeman effect
! -> Stark broadening (impact approximation for ions and electrons)
! -> Natural broadening (NIST data for hydrogen)
! -> Doppler shift
!On label 123 : setting of the estimated interval of omega
! dr:           use T_g instead of T (=Te=Ti) for Doppler width

c  on input:
c   N: Plasma Density (cm-3), ne=ni
c   Te: Plasma Temperature, Electrons (eV)
c   Ti: Plasma Temperature, Ions (eV)
c   T_g: Neutral Gas Temperature (eV) (only for estimating line width)
c   B: Magnetic field (T)
c   ctheta2: cos**2 of: Observation angle with magnetic field
c   v: Emitter/Absorber velocity (m/s)
c   npt: Number of points on line
c   e_min: lower bound of interval (eV)
c   e_max: upper bound of interval (eV)
c  on output:
c     e_min: estimated lower bound of interval (eV)
c     e_max: estimated upper bound of interval (eV)

      implicit none

!Physical and mathematical constants
      real(dp),parameter::e=1.6022e-19
      real(dp),parameter::m_D=3.3445e-27  !DEUTERONS
      real(dp),parameter::hbar=1.0546e-34
      real(dp),parameter::me=9.1094e-31
      real(dp),parameter::alpha=7.2974e-3
      real(dp),parameter::EI=13.606
      real(dp),parameter::c=2.9979e8
      real(dp),parameter::A=6.265e+08  ! Natural broadening added in v2

      real(dp),intent(in)::N,Te,Ti,T_g,B,ctheta2,v
      real(dp),intent(inout)::E00
      real(dp)::e_min,e_max
      integer::npt
      real(dp)::omega_SF,omega_Z,gamma,gam,epsilon,
     .          omega_plus,omega_minus,
     .          omega1,omega2,omega_D,omega_D_th
      real(dp)::C1,C2,C3,C4,C5,C6,C7,C8,omega,delta_omega,line_shape,
     .          interval_omega,res
      real(dp)::r0,rr0,ssum,ssum1,del,ls_old,ls_new,
     .          omega_old,slope,ph,q,det
      integer::i

      omega_SF=alpha*alpha*EI/24.
      omega_Z=hbar*B/(2.*me)
      omega_D=.75*EI*v/c
      gamma=EIRENE_coll(N,Te,Ti,epsilon)+(hbar*A)/e
      omega_plus=.25*omega_SF+.5*omega_Z
      omega_minus=.25*omega_SF-.5*omega_Z
      omega1=.25*sqrt(4.*omega_Z*omega_Z+4.*omega_Z*omega_SF+
     .       9.*omega_SF*omega_SF)
      omega2=.25*sqrt(4.*omega_Z*omega_Z-4.*omega_Z*omega_SF+
     .       9.*omega_SF*omega_SF)
      C1=.5-.5*(omega_Z+omega_minus)/omega1
      C2=.5+.5*(omega_Z+omega_minus)/omega1
      C3=.5+.5*(omega_Z-omega_plus)/omega2
      C4=.5-.5*(omega_Z-omega_plus)/omega2
      C5=.5+.5*omega_plus/omega1
      C6=.5-.5*omega_plus/omega1
      C7=.5+.5*omega_minus/omega2
      C8=.5-.5*omega_minus/omega2

c  for sampling: always find omega-range

      npt=100
      omega_D_th=.75*EI*sqrt(2.*e*T_g/m_D)/c
      interval_omega=4.*(omega_Z+omega_plus+omega2)+20.*gamma+
     .               8.*omega_D_th
      e_min=-0.5*interval_omega
      e_max=+0.5*interval_omega

      delta_omega=interval_omega/real(npt-1,dp)
      omega=e_min

      r0 = ranf_eirene()
      ssum=0.0
      ls_old=0.
      del=delta_omega
      omega_old=omega-del
      do i=1,npt
c  EIRENE function "Lorentz" needs FWHM, gamma is HWHM.
        gam=2._dp*gamma
        line_shape=(.5+.5*ctheta2)*
     .  (EIRENE_Lorentz(omega-omega_D-2.*omega_plus,gam)
     .  +EIRENE_Lorentz(omega-omega_D-2.*omega_minus,gam)
     .  +C1*EIRENE_Lorentz(omega-omega_D-omega_Z+omega_minus-omega1,gam)
     .  +C2*EIRENE_Lorentz(omega-omega_D-omega_Z+omega_minus+omega1,gam)
     .  +C3*EIRENE_Lorentz(omega-omega_D+omega_Z+omega_plus-omega2,gam)
     .  +C4*EIRENE_Lorentz(omega-omega_D+omega_Z+omega_plus+omega2,gam))
     .  +(1.-ctheta2)*
     .  (C5*EIRENE_Lorentz(omega-omega_D+omega_plus-omega1,gam)
     .  +C6*EIRENE_Lorentz(omega-omega_D+omega_plus+omega1,gam)
     .  +C7*EIRENE_Lorentz(omega-omega_D+omega_minus-omega2,gam)
     .  +C8*EIRENE_Lorentz(omega-omega_D+omega_minus+omega2,gam))
c  original line_shape (sum of weighted Lorentzians)
c  is normalized to 4.0, for each given incident angle theta.
        ls_new=line_shape*0.25
        slope=(ls_new-ls_old)/del
c  ssum  is the integral up to omega-del
c  ssum1 is the integral up to omega, Trapez rule
        ssum1=ssum+0.5*(ls_new+ls_old)*del
        if (r0.lt.ssum1) then
c  at this point: ssum.le.r0..lt.ssum1
          rr0=r0-ssum
          if (slope.eq.0.0) then
            res=E00+omega_old+rr0
          elseif (slope.gt.0.0) then
            ph=(omega_old-ls_old/slope)
            q=omega_old*omega_old-2./slope*(ls_old*omega_old+rr0)
            det=ph*ph-q
            if (det.le.0.0) then
              write (iunout,*) 'negat. det. in sam_zm_stark'
              det=0.0
            endif
            res=E00+ph+sqrt(det)
          elseif (slope.lt.0.0) then
            ph=(omega_old-ls_old/slope)
            q=omega_old*omega_old-2./slope*(ls_old*omega_old+rr0)
            det=ph*ph-q
            if (det.le.0.0) then
              write (iunout,*) 'negat. det. in sam_zm_stark'
              det=0.0
            endif
            res=E00+ph-sqrt(det)
          endif
          exit
        endif
c  prepare next interval
        ssum=ssum1
        ls_old=ls_new
        omega_old=omega
        omega=omega+del
      end do
      return
      end function EIRENE_sam_zm_stark1


      FUNCTION EIRENE_sam_zm_stark(N,Te,Ti,T_g,B,
     .                      ctheta2,E00,v)
     .                      result(res)

c  this is basically "function zm_stark_profile", modified for
c  sampling from rather than evaluating of the line profile.
c  v is the velocity of the emitting atom.


c  original line_shape is a weighted sum of Lorentzians.
c  It is normalized to 4.0, for each given incident angle theta.
c  here: randomly chose Lorentzian, then sample from this Lorentzian.
c  to be done: speed up of evaluation of sum of Lorentzians

!********** DEUTERIUM LYMAN ALPHA LINE SHAPE CALCULATION **********
!******************************************************************
!Sampling of the line shape
!Retained effects :
! -> fine structure
! -> Zeeman effect
! -> Stark broadening (impact approximation for ions and electrons)
! -> Natural broadening (NIST data for hydrogen)
! -> Doppler shift
!On label 123 : setting of the estimated interval of omega
! dr:           use T_g instead of T (=Te=Ti) for Doppler width

c  on input:
c   N: Plasma Density (cm-3), ne=ni
c   Te: Plasma Temperature, Electrons (eV)
c   Ti: Plasma Temperature, Ions (eV)
c   T_g: Neutral Gas Temperature (eV) (only for estimating line width)
c   B: Magnetic field (T)
c   ctheta2: cos**2 of: Observation angle with magnetic field
c   v: Emitter/Absorber velocity (m/s)
c  on output:
c   res: random number sampled from zeemann-stark-profile
      implicit none

!Physical and mathematical constants
      real(dp),parameter::e=1.6022e-19
      real(dp),parameter::hbar=1.0546e-34
      real(dp),parameter::me=9.1094e-31
      real(dp),parameter::alpha=7.2974e-3
      real(dp),parameter::EI=13.606
      real(dp),parameter::c=2.9979e8
      real(dp),parameter::A=6.265e+08  ! Natural broadening added in v2

      real(dp),intent(in)::N,Te,Ti,T_g,B,ctheta2,v
      real(dp),intent(inout)::E00
      real(dp)::omega_SF,omega_Z,gamma,gam,epsilon,
     .          omega_plus,omega_minus,
     .          omega1,omega2,omega_D
      real(dp)::C1,C2,C3,C4,C5,C6,C7,C8,omega,
     .          res,Ci(-1:8),x(-1:8)
      real(dp)::r0,tha,thb,shift
      integer::i

      omega_SF=alpha*alpha*EI/24.
      omega_Z=hbar*B/(2.*me)
      omega_D=.75*EI*v/c
      gamma=EIRENE_coll(N,Te,Ti,epsilon)+(hbar*A)/e
      omega_plus=.25*omega_SF+.5*omega_Z
      omega_minus=.25*omega_SF-.5*omega_Z
      omega1=.25*sqrt(4.*omega_Z*omega_Z+4.*omega_Z*omega_SF+
     .       9.*omega_SF*omega_SF)
      omega2=.25*sqrt(4.*omega_Z*omega_Z-4.*omega_Z*omega_SF+
     .       9.*omega_SF*omega_SF)
      C1=.5-.5*(omega_Z+omega_minus)/omega1
      C2=.5+.5*(omega_Z+omega_minus)/omega1
      C3=.5+.5*(omega_Z-omega_plus)/omega2
      C4=.5-.5*(omega_Z-omega_plus)/omega2
      C5=.5+.5*omega_plus/omega1
      C6=.5-.5*omega_plus/omega1
      C7=.5+.5*omega_minus/omega2
      C8=.5-.5*omega_minus/omega2

         tha=0.5+.5*ctheta2
         thb=1.0   -ctheta2
         omega=e00
         x(-1)=omega-omega_D-2.*omega_plus
         Ci(-1)=1.*tha
         x(0) =omega-omega_D-2.*omega_minus
         Ci(0)=1.*tha+Ci(-1)
         x(1) =omega-omega_D-omega_Z+omega_minus-omega1
         Ci(1)=C1*tha+Ci(0)
         x(2) =omega-omega_D-omega_Z+omega_minus+omega1
         Ci(2)=C2*tha+Ci(1)
         x(3) =omega-omega_D+omega_Z+omega_plus-omega2
         Ci(3)=C3*tha+Ci(2)
         x(4) =omega-omega_D+omega_Z+omega_plus+omega2
         Ci(4)=C4*tha+Ci(3)
         x(5) =omega-omega_D+omega_plus-omega1
         Ci(5)=C5*thb+Ci(4)
         x(6) =omega-omega_D+omega_plus+omega1
         Ci(6)=C6*thb+Ci(5)
         x(7) =omega-omega_D+omega_minus-omega2
         Ci(7)=C7*thb+Ci(6)
         x(8) =omega-omega_D+omega_minus+omega2
         Ci(8)=C8*thb+Ci(7)


      r0 = ranf_eirene()*Ci(8)
      do i=-1,7
        if (r0.lt.ci(i)) goto 4711
      enddo
      i=8
 4711 continue

c  sample from lorentzian no. i
c  EIRENE function "Lorentz" needs FWHM, gamma is HWHM.
        gam=2._dp*gamma
        shift=x(i)
        res=EIRENE_sam_lorentz(gam,shift)

      return
      end function EIRENE_sam_zm_stark



      REAL(dp) FUNCTION EIRENE_SAM_VOIGT(alph,dnd,drft,e00) result(res)
!
!  sample from voigt profile. convolution of
!  lorentz (FWHM=alph) and doppler (Doppler-width=dnd, drift: drft)
!
      IMPLICIT NONE
      real(dp),intent(in) :: alph,dnd,drft,e00
      real(dp) :: xx,f1,v1,v2,s,f2,ar,sig,alphh

c lorentz fwhm=alph
      alphh=alph*0.5_dp
      xx = PIHA*(RANF_EIRENE()*2.-1.)
      f1=alphh*tan(xx)

c doppler = gauss with sig=dnd/sqrt(2)
      sig=dnd/sqrt(2.)
      do
         v1=2.*RANF_EIRENE()-1.
         v2=2.*RANF_EIRENE()-1.
         s=v1*v1+v2*v2
         if(s <= 1.) exit
      enddo
      ar=log(s)
      f2=v1*sqrt(-(ar+ar)/s)*sig

c gauss*lorentz:
      res=f1+f2+drft+e00
      return
      END FUNCTION EIRENE_SAM_VOIGT








      REAL(dp) FUNCTION EIRENE_SAM_VDWQS(dvdw,x0)
     .         result(res)
c  sample exponential (quasistatic van der Waals, red wing)
c  in wavelength units: return l0-l00, l0>l00
c  in energy (frequency) units: return e00-e0, e0<e00
c  x0 is needed in case of energy sampling to ensure positive energies
c  from convoluted profile.
      IMPLICIT NONE
      real(dp), intent(in) :: dvdw,x0
      real(dp) :: xx,sig,beta,xxq

      beta=piqu*dvdw
      sig=1._dp/sqrt(beta)
      do
        xx=EIRENE_sam_doppler(sig,0._dp,0._dp)
        xxq=xx*xx
        if (xxq > 1._dp/x0) exit
      end do
      res=1._dp / xxq

      return
      end function EIRENE_sam_vdwqs
c
c
c
      SUBROUTINE EIRENE_PH_WOFZ (XI, YI, U, V, FLAG)
C      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 16, NO. 1, PP. 47.
C
C  GIVEN A COMPLEX NUMBER Z = (XI,YI), THIS SUBROUTINE COMPUTES
C  THE VALUE OF THE FADDEEVA FUNCTION W(Z) = EXP(-Z**2)*ERFC(-I*Z),
C  WHERE ERFC IS THE COMPLEX COMPLEMENTARY ERROR FUNCTION AND I
C  MEANS SQRT(-1).
C  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT
C  IS 14 SIGNIFICANT DIGITS; IN THE 3RD AND 4TH IT IS 13 SIGNIFICANT
C  DIGITS OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO
C  OF THE FUNCTION.
C  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION.
C
C
C  THE CODE CONTAINS A FEW COMPILER-DEPENDENT PARAMETERS :
C     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
C                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
C                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
C                FLOATING-POINT ARITHMETIC
C     RMAXEXP  = LN(RMAX) - LN(2)
C     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
C                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
C  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY AR/POE DEFINED WILL
C  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS
C
C
C  PARAMETER LIST
C     XI     = REAL      PART OF Z
C     YI     = IMAGINARY PART OF Z
C     U      = REAL      PART OF W(Z)
C     V      = IMAGINARY PART OF W(Z)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  XI, YI      ARE THE INPUT-PARAMETERS
C  U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C
C  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
C
C  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE
C  PUT TO 0 UPON UNDERFLOW;
C
C  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT COMPUTATION OF
C  THE COMPLEX ERROR FUNCTION, ACM TRANS. MATH. SOFTWARE.
C
*
*
*
*
!      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      LOGICAL :: A, B
      LOGICAL, intent(inout) :: FLAG
      real(dp),PARAMETER   :: FACTOR   = 1.12837916709551257388D0,
     *                        RMAXREAL = 1.34078079299426D+154,
     *                        RMAXEXP  = 709.089565710919D0,
     *                        RMAXGONI = 3.53711887601422D+15
c     *                        RMAXREAL = 0.5D+154,
c     *                        RMAXEXP  = 708.503061461606D0,
c     *                        RMAXGONI = 3.53711887601422D+15
      real(dp), intent(in) :: xi,yi
      real(dp), intent(out):: u,v
      real(dp) :: xabs,yabs,x,y,qrho,xabsq,xquad,yquad,xsum,ysum,
     .            xaux,u1,v1,daux,u2,v2,h,rx,ry,sx,sy,tx,ty,c,qlambda,
     .            w1,h2
      integer :: n,j,i,kapn,nu,np1
*
      FLAG = .FALSE.
*
      XABS = DABS(XI)
      YABS = DABS(YI)
      X    = XABS/6.3
      Y    = YABS/4.4
*
C
C     THE FOLLOWING IF STATEMENT PROTECTS
C     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
C
      IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
*
      QRHO = X**2 + Y**2
*
      XABSQ = XABS**2
      XQUAD = XABSQ - YABS**2
      YQUAD = 2*XABS*YABS
*
      A     = QRHO.LT.0.085264D0
*
      IF (A) THEN
C
C  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA FUNCTION IS EVALUATED
C  USING A POWER SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
C  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
        QRHO  = (1-0.85*Y)*DSQRT(QRHO)
        N     = IDNINT(6 + 72*QRHO)
        J     = 2*N+1
        XSUM  = 1.0/J
        YSUM  = 0.0D0
        DO 10 I=N, 1, -1
          J    = J - 2
          XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
          YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
          XSUM = XAUX + 1.0/J
   10   CONTINUE
        U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
        V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
        DAUX =  DEXP(-XQUAD)
        U2   =  DAUX*DCOS(YQUAD)
        V2   = -DAUX*DSIN(YQUAD)
*
        U    = U1*U2 - V1*V2
        V    = U1*V2 + V1*U2
*
      ELSE
C
C  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
C  CONTINUED FRACTION
C  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
C  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
C  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
C  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
C  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
C  TO OBTAIN THE REQUIRED ACCURACY
C  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
C  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
C
*
        IF (QRHO.GT.1.0) THEN
          H    = 0.0D0
          KAPN = 0
          QRHO = DSQRT(QRHO)
          NU   = INT(3 + (1442/(26*QRHO+77)))
        ELSE
          QRHO = (1-Y)*DSQRT(1-QRHO)
          H    = 1.88*QRHO
          H2   = 2*H
          KAPN = IDNINT(7  + 34*QRHO)
          NU   = IDNINT(16 + 26*QRHO)
        ENDIF
*
        B = (H.GT.0.0)
*
        IF (B) QLAMBDA = H2**KAPN
*
        RX = 0.0
        RY = 0.0
        SX = 0.0
        SY = 0.0
*
        DO 11 N=NU, 0, -1
          NP1 = N + 1
          TX  = YABS + H + NP1*RX
          TY  = XABS - NP1*RY
          C   = 0.5/(TX**2 + TY**2)
          RX  = C*TX
          RY  = C*TY
          IF ((B).AND.(N.LE.KAPN)) THEN
            TX = QLAMBDA + SX
            SX = RX*TX - RY*SY
            SY = RY*TX + RX*SY
            QLAMBDA = QLAMBDA/H2
          ENDIF
   11   CONTINUE
*
        IF (H.EQ.0.0) THEN
          U = FACTOR*RX
          V = FACTOR*RY
        ELSE
          U = FACTOR*SX
          V = FACTOR*SY
        END IF
*
        IF (YABS.EQ.0.0) U = DEXP(-XABS**2)
*
      END IF
*
*
C
C  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
C
*
      IF (YI.LT.0.0) THEN
*
        IF (A) THEN
          U2    = 2*U2
          V2    = 2*V2
        ELSE
          XQUAD =  -XQUAD
*
C
C         THE FOLLOWING IF STATEMENT PROTECTS 2*EXP(-Z**2)
C         AGAINST OVERFLOW
C
          IF ((YQUAD.GT.RMAXGONI).OR.
     *        (XQUAD.GT.RMAXEXP)) GOTO 100
*
          W1 =  2*DEXP(XQUAD)
          U2  =  W1*DCOS(YQUAD)
          V2  = -W1*DSIN(YQUAD)
        END IF
*
        U = U2 - U
        V = V2 - V
        IF (XI.GT.0.0) V = -V
      ELSE
        IF (XI.LT.0.0) V = -V
      END IF
*
      RETURN
*
  100 FLAG = .TRUE.
      RETURN
*
      END SUBROUTINE EIRENE_PH_WOFZ


c
c
c
      complex(dp) FUNCTION EIRENE_PH_HUMLIK(x,y) result(cres)
      IMPLICIT NONE
c arguments
      real(dp), intent(in) :: x,y
c To calculate the FADDEEVA function with relative error less than 10^(-R).
c R0=1.51*EXP(1.144*R) and R1=1.60*EXP(0.554*R) can be set by the the user
c subject to the constraints 14.88<R0<460.4 and 4.85<R1<25.5
      REAL(dp) :: K,L
      real(dp), PARAMETER :: R0 = 146.7_DP, R1 = 14.67_DP ! for R=4,
                                                    ! region boundaries

c Constants
      real(dp), PARAMETER :: RRTPI = 0.56418958_DP  ! 1/sqrt(pi)
      real(dp), PARAMETER :: Y0 = 1.5_DP, Y0PY0 = Y0+Y0, Y0Q = Y0*Y0
                                                    ! for cpf12 algor.
      REAL(dp), save :: C(0:5), S(0:5), T(0:5)
c SAVE preserves values of C, S and T (static) arrays between procedure calls

      DATA C / 1.0117281_DP, -0.75197147_DP,     0.012557727_DP,
     .       0.010022008_DP, -0.00024206814_DP,  0.00000050084806_DP /
      DATA S / 1.393237_DP,   0.23115241_DP,    -0.15535147_DP,
     .     0.0062183662_DP,   0.000091908299_DP,-0.00000062752596_DP /
      DATA T / 0.31424038_DP, 0.94778839_DP,     1.5976826_DP,
     .        2.2795071_DP,   3.0206370_DP,      3.8897249_DP /

c Local variables
      INTEGER :: J                 ! Loop variable
      INTEGER :: RG1, RG2, RG3     ! y polynomial flags
      REAL(dp) :: ABX, XQ, YQ, YRRTPI ! |x|, x^2, y^2, y/SQRT(pi)
      REAL(dp) :: XLIM0, XLIM1, XLIM2, XLIM3, XLIM4 ! |x| on
                                                    ! region boundaries
      REAL(dp) :: A0, D0, D2, E0, E2, E4, H0, H2, H4, H6
                                   ! W4 temporary variables
      REAL(dp) :: P0, P2, P4, P6, P8, Z0, Z2, Z4, Z6, Z8
      real(dp) :: b1,f1,f3,f5,q1,q3,q5,q7
      REAL(dp) :: XP(0:5), XM(0:5), YP(0:5), YM(0:5)
                                   ! CPF12 temporary values
      REAL(dp) :: MQ(0:5), PQ(0:5), MF(0:5), PF(0:5)
      REAL(dp) :: D, YF, YPY0, YPY0Q

c**** Start of executable code *****************************************

      RG1 = 1                   ! Set flags
      RG2 = 1
      RG3 = 1
      YQ  = Y*Y                 ! y^2
      YRRTPI = Y*RRTPI          ! y/SQRT(pi)

c Region boundaries when both K and L are required or when R<>4
      XLIM0 = R0 - Y
      XLIM1 = R1 - Y
      XLIM3 = 3.097_DP*Y - 0.45_DP

      XLIM2 = 6.8_DP - Y
      XLIM4 = 18.1_DP*Y + 1.65_DP
      IF ( Y .LE. 0.000001_DP ) THEN ! When y<10^-6
         XLIM1 = XLIM0          ! avoid W4 algorithm
         XLIM2 = XLIM0
      ENDIF
c.....
      ABX = ABS ( X )           ! |x|
      XQ  = ABX*ABX             ! x^2
      IF ( ABX .GT. XLIM0 ) THEN ! Region 0 algorithm
         K = YRRTPI / (XQ + YQ)
         L = k*x/y
      ELSEIF ( ABX .GT. XLIM1 ) THEN ! Humlicek W4 Region 1
         IF ( RG1 .NE. 0 ) THEN ! First point in Region 1
            RG1 = 0
            A0 = YQ + 0.5_DP    ! Region 1 y-dependence
            D0 = A0*A0
            D2 = YQ + YQ - 1.0_DP
            b1 = yq - 0.5_DP
         ENDIF
         D = RRTPI / (D0 + XQ*(D2 + XQ))
         K = D*Y * (A0 + XQ)
         L = d*x * (b1 + xq)
      ELSEIF ( ABX .GT. XLIM2 ) THEN ! Humlicek W4 Region 2
         IF ( RG2 .NE. 0 ) THEN ! First point in Region 2
            RG2 = 0
            H0 =  0.5625_DP + YQ*(4.5_DP +
     &                        YQ*(10.5_DP + YQ*(6.0_DP + YQ)))
                                ! Region 2 y-dependence
            H2 = -4.5_DP   + YQ*(9.0_DP + YQ*( 6.0_DP + YQ* 4.0_DP))
            H4 = 10.5_DP   - YQ*(6.0_DP - YQ*  6.0_DP)
            H6 = -6.0_DP   + YQ* 4.0_DP
            E0 =  1.875_DP + YQ*(8.25_DP + YQ*(5.5_DP + YQ))
            E2 =  5.25_DP  + YQ*(1.0_DP  + YQ* 3.0_DP)
            E4 =  0.75_DP*H6
            f1 = -1.875_DP + yq*(5.25_DP + yq*(4.5_DP + yq))
            f3 =  8.25_DP  - yq*(1.0_DP  - yq* 3.0_DP)
            f5 = -5.5_DP   + yq* 3.0_DP
         ENDIF
         D = RRTPI / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))))
         K = D*Y   *(E0 + XQ*(E2 + XQ*(E4 + XQ)))
         L = d*x   *(f1 + xq*(f3 + xq*(f5 + xq)))
      ELSEIF ( ABX .LT. XLIM3 ) THEN ! Humlicek W4 Region 3
         IF ( RG3 .NE. 0 ) THEN ! First point in Region 3
            RG3 = 0 ! Region 3 y-dependence
            Z0 = 272.1014_DP  + Y*(1280.829_DP + Y*(2802.870_DP
     &           + Y*(3764.966_DP + Y*(3447.629_DP
     &           + Y*(2256.981_DP + Y*(1074.409_DP
     &           + Y*(369.1989_DP + Y*(88.26741_DP
     &           + Y*(13.39880_DP + Y)))))))))
            Z2 = 211.678_DP   + Y*(902.3066_DP + Y*(1758.336_DP
     &           + Y*(2037.310_DP + Y*(1549.675_DP
     &           + Y*(793.4273_DP + Y*(266.2987_DP
     &           + Y*(53.59518_DP + Y*5.0_DP)))))))
            Z4 = 78.86585_DP  + Y*(308.1852_DP + Y*(497.3014_DP
     &           + Y*(479.2576_DP + Y*(269.2916_DP
     &           + Y*(80.39278_DP + Y*10.0_DP)))))
            Z6 = 22.03523_DP  + Y*(55.02933_DP + Y*(92.75679_DP
     &           + Y*(53.59518_DP + Y*10.0_DP)))
            Z8 = 1.496460_DP  + Y*(13.39880_DP + Y*5.0_DP)
            P0 = 153.5168_DP  + Y*(549.3954_DP + Y*(919.4955_DP
     &           + Y*(946.8970_DP + Y*(662.8097_DP
     &           + Y*(328.2151_DP + Y*(115.3772_DP
     &           + Y*(27.93941_DP + Y*(4.264678_DP
     &           + Y*0.3183291_DP))))))))
            P2 = -34.16955_DP + Y*(-1.322256_DP + Y*(124.5975_DP
     &           + Y*(189.7730_DP + Y*(139.4665_DP
     &           + Y*(56.81652_DP + Y*(12.79458_DP
     &           + Y*1.2733163_DP))))))
            P4 = 2.584042_DP  + Y*(10.46332_DP + Y*(24.01655_DP
     &           + Y*(29.81482_DP + Y*(12.79568_DP
     &           + Y*1.9099744_DP))))
            P6= -0.07272979_DP + Y*(0.9377051_DP + Y*(4.266322_DP
     &           + Y*1.273316_DP))
            P8 = 0.0005480304_DP + Y*0.3183291_DP
            q1 = 173.2355_DP  + y*(508.2585_DP + y*(685.8378_DP
     .           + y*(557.5178_DP + y*(301.3208_DP
     .           + y*(111.0528_DP + y*(27.62940_DP
     .           + y*(4.264130_DP + y*0.3183291_DP)))))))
            q3 = 18.97431_DP  + y*(100.7375_DP + y*(160.4013_DP
     .           + y*(130.8905_DP + y*(55.88650_DP
     .           + y*(12.79239_DP + y*1.273316_DP)))))
            q5 = 7.985877_DP  + y*(19.83766_DP + y*(28.88480_DP
     .           + y*(12.79239_DP + y*1.909974_DP)))
            q7 = 0.6276985_DP + y*(4.264130_DP + y*1.273316_DP)
         ENDIF
         D = 1.7724538_DP /
     .       (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))))
         K = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))))
         L = d*x*(q1+xq*(q3+xq*(q5+xq*(q7 + xq*0.3183291_DP))))
      ELSE                      ! Humlicek CPF12 algorithm
         YPY0 = Y + Y0
         YPY0Q = YPY0*YPY0
         K = 0.0_DP
         L = 0.0_DP
         DO J = 0, 5
            D = X - T(J)
            MQ(J) = D*D
            MF(J) = 1.0 / (MQ(J) + YPY0Q)
            XM(J) = MF(J)*D
            YM(J) = MF(J)*YPY0
            D = X + T(J)
            PQ(J) = D*D
            PF(J) = 1.0 / (PQ(J) + YPY0Q)
            XP(J) = PF(J)*D
            YP(J) = PF(J)*YPY0
            L=L+c(j)*(xm(j)+xp(j)) + s(j)*(ym(j)-yp(j))
         ENDDO

         IF ( ABX .LE. XLIM4 ) THEN ! Humlicek CPF12 Region I
            DO J = 0, 5
               K = K + C(J)*(YM(J)+YP(J)) - S(J)*(XM(J)-XP(J))
            ENDDO
         ELSE                   ! Humlicek CPF12 Region II
            YF   = Y + Y0PY0
            DO J = 0, 5
               K = K
     &     + (C(J)*(MQ(J)*MF(J)-Y0*YM(J)) + S(J)*YF*XM(J)) / (MQ(J)+Y0Q)
     &     + (C(J)*(PQ(J)*PF(J)-Y0*YP(J)) - S(J)*YF*XP(J)) / (PQ(J)+Y0Q)
            ENDDO
            K = Y*K + EXP ( -XQ )
         ENDIF
      ENDIF
*.....
! kind parameter added for better precision
      cres = CMPLX(K,L,dp)
      RETURN
      END FUNCTION EIRENE_PH_HUMLIK

c EIRENE utilities:

c
c  experimental routine for black body removal
c


      function EIRENE_sam_cutoff (kk,ictoff,icell) result (res)
      integer, intent(in) :: kk,ictoff,icell
      real(dp) :: xintsum, xi1, xi2, res, e00, xx, phi_xi1, fwhm,
     .            shift, dvdw, l00, l0cut, e_xi1, a, b, c, rd,
     .            phi_xx, fcut, e1, e2, es, sint, lright, l1, l2, ls,
     .            r, phlamcut, phidif
      real(dp), allocatable, save :: ali(:,:), bli(:,:), xintli(:,:),
     .                               eli(:,:), alnorm(:), phili(:)
      real(dp), allocatable, save :: ari(:,:), bri(:,:), xintri(:,:),
     .                               lri(:,:), arnorm(:), phiri(:)
      logical, allocatable, save :: visitl(:), visitr(:)
      integer, save :: ictsave = -1
      integer :: icount, nrjp, irj
      external :: eirene_exit_own

      nrjp = reacdat(kk)%phr%line%nrjprt

      if ((ictsave /= ictoff) .and. (nrjp > 1)) then
        ictsave = ictoff
        if (allocated(ali)) then
          deallocate (ali)
          deallocate (bli)
          deallocate (xintli)
          deallocate (eli)
          deallocate (alnorm)
          deallocate (phili)
          deallocate (ari)
          deallocate (bri)
          deallocate (xintri)
          deallocate (lri)
          deallocate (arnorm)
          deallocate (phiri)
        end if
        allocate (ali(nrjp,nrad))
        allocate (bli(nrjp,nrad))
        allocate (eli(0:nrjp,nrad))
        allocate (xintli(nrjp,nrad))
        allocate (alnorm(nrad))
        allocate (phili(0:nrjp))
        allocate (ari(nrjp,nrad))
        allocate (bri(nrjp,nrad))
        allocate (lri(0:nrjp,nrad))
        allocate (xintri(nrjp,nrad))
        allocate (arnorm(nrad))
        allocate (phiri(0:nrjp))

        if (.not.allocated(visitl)) then
           allocate (visitl(nrad))
           allocate (visitr(nrad))
        end if
        visitl = .false.
        visitr = .false.
      end if


      xintsum = xintleft(ictoff,icell) +
     .         (xint_inf(ictoff,icell)-xintright(ictoff,icell))
      e00 = reacdat(kk)%phr%line%e0
      l00=hpcl/e00

      icount=0
      if (xintsum*ranf_eirene() < xintleft(ictoff,icell)) then

!  sample energy in left wing
c  use energy scale

        if (nrjp == 0) then

c  use rectangle for rejection

          call  EIRENE_lorvdwprof(kk,icell,fwhm,shift,dvdw,.true.)

          do
            icount = icount + 1
            xi2 = ranf_eirene() * phicut(ictoff,icell)
            xi1 = ranf_eirene() * ecutleft(ictoff,icell)
            xx = xi1-E00
            phi_xi1 = EIRENE_ph_lorvdw(-xx, fwhm, -shift, dvdw, icell)
            if (xi2 <= phi_xi1) then
              res = xi1
c             write (iunout,*) ' left interval in sam_cutoff ',icount,
c    .                           res
              exit
            end if
          end do

        elseif (nrjp == 1) then

!  use triangle for rejection

          call  EIRENE_lorvdwprof(kk,icell,fwhm,shift,dvdw,.true.)
          phidif = phicut(ictoff,icell) - phi_zero(ictoff,icell)

          do
            icount = icount + 1
            r = sqrt(ranf_eirene())
            xi1 = ecutleft(ictoff,icell) * r
!  f(x) = a*x + phi_0, a=phidif/ecutleft
!  f(xi1) = a*xi1 = phidif/ecutleft*ecutleft*r+phi_0 = phidif*r+phi_0
            fcut = phidif * r + phi_zero(ictoff,icell)
            xi2 = fcut * ranf_eirene()
            xx = xi1-E00
            phi_xi1 = EIRENE_ph_lorvdw(-xx, fwhm, -shift, dvdw, icell)
            if (xi2 <= phi_xi1) then
              res = xi1
c             write (iunout,*) ' left interval in sam_cutoff ',icount,
c    .                           res
              exit
            end if
          end do

        elseif (nrjp > 1) then

!  use polygon for rejection

          if (.not.visitl(icell)) then

            eli(:,icell) = 0
            phili(0) = 0
            phili(1:nrjp-1) = phi_rj_left(ictoff,1:nrjp-1,icell)
            phili(nrjp) = phicut(ictoff,icell)

            do irj = 1, nrjp
              eli(irj,icell) = irj*ecutleft(ictoff,icell)/nrjp
            end do

            alnorm(icell) = 0._dp
            do irj = 1, nrjp
              ali(irj,icell) = (phili(irj-1) - phili(irj)) /
     .                         (eli(irj-1,icell) - eli(irj,icell))
              bli(irj,icell) = phili(irj)-ali(irj,icell)*eli(irj,icell)
              xintli(irj,icell) = 0.5_dp * ali(irj,icell) *
     .                           (eli(irj,icell)**2-eli(irj-1,icell)**2)
     .                          + bli(irj,icell) *
     .                           (eli(irj,icell) - eli(irj-1,icell))
              alnorm(icell) = alnorm(icell) + xintli(irj,icell)
            end do

            visitl(icell) = .true.

          end if

          call  EIRENE_lorvdwprof(kk,icell,fwhm,shift,dvdw,.true.)

          do

            icount = icount + 1
            xi1 = ranf_eirene() * alnorm(icell)

            sint = 0._dp
            do irj = 1, nrjp
              if (xi1 <= sint+xintli(irj,icell)) exit
              sint = sint + xintli(irj,icell)
            end do

            if (irj > nrjp) then
              write (iunout,*) ' nonsense in sam_cutoff '
              call EIRENE_exit_own(1)
            end if

            a = 0.5_dp * ali(irj,icell)
            b = bli(irj,icell) / a
            c = (-0.5_dp*ali(irj,icell)*eli(irj-1,icell)**2 -
     .           bli(irj,icell)*eli(irj-1,icell) -
     .           xi1 + sint) / a

            rd = 0.25_dp * b*b - c
            if (rd < 0) then
              write (iunout,*) ' negative wurzel in sam_cutoff '
              call EIRENE_exit_own(1)
            end if

            e1 = -0.5_dp * b + sqrt(rd)
            e2 = -0.5_dp * b - sqrt(rd)
            if ((eli(irj-1,icell)<e1).and.(e1 <= eli(irj,icell))) es=e1
            if ((eli(irj-1,icell)<e2).and.(e2 <= eli(irj,icell))) es=e2

            fcut = ali(irj,icell)*es + bli(irj,icell)
            xi2 = ranf_eirene()*fcut

            xx = es-E00
            phi_xx = EIRENE_ph_lorvdw(-xx, fwhm, -shift, dvdw, icell)
            if (xi2 <= phi_xx) then
              res = es
c             write (iunout,*) ' left interval in sam_cutoff ',icount,
c    .                           res
              exit
            end if
          end do

        end if

      else

!  sample energy in right wing
c  use wavelength scale

        if (nrjp == 0) then

!  use rectangle for rejection

          l0cut = hpcl/ecutright(ictoff,icell)
          phlamcut = phicut(ictoff,icell)*hpcl/(l0cut*l0cut)

          call  EIRENE_lorvdwprof(kk,icell,fwhm,shift,dvdw,.false.)

          do
            icount = icount + 1
            xi2 = ranf_eirene() * phlamcut
            xi1 = ranf_eirene() * l0cut
c  convert sampled wavelength to energy (eV)
            e_xi1 = hpcl/xi1
            xx = xi1-l00
            phi_xi1 = EIRENE_ph_lorvdw(xx, fwhm, shift, dvdw, icell)
!pb            phi_xi1 = phi_xi1*hpcl/(e_xi1*e_xi1)
            if (xi2 <= phi_xi1) then
              res = e_xi1
c             write (iunout,*) ' right interval in sam_cutoff ',icount,
c    .                         res
              exit
            end if
          end do

        elseif (nrjp == 1) then

!  use triangle for rejection

          l0cut = hpcl/ecutright(ictoff,icell)
          phlamcut = phicut(ictoff,icell)*hpcl/(l0cut*l0cut)

          call  EIRENE_lorvdwprof(kk,icell,fwhm,shift,dvdw,.false.)
          phidif = phlamcut - phi_inf(ictoff,icell)

          do
            icount = icount + 1
            r = sqrt(ranf_eirene())
            xi1 = l0cut * r
            e_xi1 = hpcl/xi1
!  f(x) = a*x+phi_inf, a=(phicut-phi_inf)/l0cut
!  f(xi1) = a*xi1+phi_inf = phidif/l0cut*l0cut*r+phi_inf = phidif*r+phi_inf
!pb            fcut = phlamcut * r
            fcut = phidif * r + phi_inf(ictoff,icell)
            xi2 = fcut * ranf_eirene()
            xx =xi1-l00
            phi_xi1 = EIRENE_ph_lorvdw(xx, fwhm, shift, dvdw, icell)
!pb            phi_xi1 = phi_xi1*hpcl/(e_xi1*e_xi1)
            if (xi2 <= phi_xi1) then
              res = e_xi1
c             write (iunout,*) ' right interval in sam_cutoff ',icount,
c    .                           res
              exit
            end if
          end do


        elseif (nrjp > 1) then

!  use polygon for rejection

          if (.not.visitr(icell)) then

            lri(:,icell) = 0
            phiri(0) = 0
            phiri(1:nrjp-1) = phi_rj_right(ictoff,1:nrjp-1,icell)
            phiri(nrjp) = phicut(ictoff,icell)

            lright = hpcl/ecutright(ictoff,icell)
            do irj = 1, nrjp
              lri(irj,icell) = irj*lright/nrjp
            end do

            arnorm(icell) = 0._dp
            do irj = 1, nrjp

              ari(irj,icell) = (phiri(irj-1) - phiri(irj)) /
     .                         (lri(irj-1,icell) - lri(irj,icell))
              bri(irj,icell) = phiri(irj)-ari(irj,icell)*lri(irj,icell)
              xintri(irj,icell) = 0.5_dp * ari(irj,icell) *
     .                           (lri(irj,icell)**2-lri(irj-1,icell)**2)
     .                          + bri(irj,icell) *
     .                           (lri(irj,icell)-lri(irj-1,icell))
              arnorm(icell) = arnorm(icell) + xintri(irj,icell)
            end do

            visitr(icell) = .true.

          end if

          call EIRENE_lorvdwprof(kk,icell,fwhm,shift,dvdw,.true.)

          do
            icount = icount + 1

            xi1 = ranf_eirene() * arnorm(icell)
            sint = 0._dp
            do irj = 1, nrjp
              if (xi1 <= sint+xintri(irj,icell)) exit
              sint = sint + xintri(irj,icell)
            end do

            if (irj > nrjp) then
              write (iunout,*) ' nonsense in sam_cutoff '
              call EIRENE_exit_own(1)
            end if

            a = 0.5_dp * ari(irj,icell)
            b = bri(irj,icell) / a
            c = (-0.5_dp*ari(irj,icell)*lri(irj-1,icell)**2 -
     .           bri(irj,icell)*lri(irj-1,icell) -
     .           xi1 + sint) / a

            rd = 0.25_dp * b*b - c
            if (rd < 0) then
              write (iunout,*) ' negative wurzel in sam_cutoff '
              call EIRENE_exit_own(1)
            end if

            l1 = -0.5_dp * b + sqrt(rd)
            l2 = -0.5_dp * b - sqrt(rd)
            if ((lri(irj-1,icell)<l1).and.(l1 <= lri(irj,icell))) ls=l1
            if ((lri(irj-1,icell)<l2).and.(l2 <= lri(irj,icell))) ls=l2

            fcut = ari(irj,icell)*ls + bri(irj,icell)
            xi2 = ranf_eirene()*fcut

            es = hpcl/ls
            xx = es-E00
            phi_xx = EIRENE_ph_lorvdw(-xx, fwhm, -shift, dvdw, icell)
            if (xi2 <= phi_xx) then
              res = es
c             write (iunout,*) ' right interval in sam_cutoff ',icount,
c    .                           res
              exit
            end if
          end do

        end if
      end if

      return
      end function EIRENE_sam_cutoff

cdr nov 20:  this next routine is now entirely obsolete 

      SUBROUTINE EIRENE_BROADCAST_PHOTON(ME)
      USE EIRMOD_MPI
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ME
      INTEGER :: IER, NNROT, NRC, KK, IATM, IPHOT
#if ( defined(USE_MPI) && !defined(OPEN_MPI) && !defined(GFORTRAN) && !defined(CRAY) )
      EXTERNAL :: MPI_BCAST
#endif

      IF (NPHOTI > 0) THEN
        nnrot=0
        do iatm=1,natmi
          if(nrca(iatm) > 0) then
            do nrc=1,nrca(iatm)
              kk=ireaca(iatm,nrc)
              if(iswr(kk) == 7) then
                NNROT=NNROT+1
              endif
            enddo
          endif
        enddo
!pb  out! not needed for the time being
!pb        IF (ME .NE. 0) call PH_ALLOC_XSECTA(nnrot)
        IF (NNROT > 0) THEN
          CALL MPI_BCAST (PHV_LGAOT,(NATM+1)*(NNROT+1)*6,
     .                    MPI_INTEGER,0,MPI_COMM_WORLD,ier)
          CALL MPI_BCAST (PHV_NAOTI,NATM,
     .                    MPI_INTEGER,0,MPI_COMM_WORLD,ier)
          CALL MPI_BCAST (PHV_IESTOTAT,(NATM+1)*NNROT*3,
     .                    MPI_INTEGER,0,MPI_COMM_WORLD,ier)
          CALL MPI_BCAST (PHV_N1STOTAT,(NATM+1)*NNROT*3,
     .                    MPI_INTEGER,0,MPI_COMM_WORLD,ier)
          CALL MPI_BCAST (PHV_N2NDOTAT,(NATM+1)*NNROT*3,
     .                    MPI_INTEGER,0,MPI_COMM_WORLD,ier)
        END IF

      END IF

      CALL MPI_BARRIER(MPI_COMM_WORLD,ier)

      RETURN
      END SUBROUTINE EIRENE_BROADCAST_PHOTON

      END MODULE EIRMOD_PHOTON
