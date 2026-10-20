      subroutine eirene_write_json_file(filename)

      use eirmod_precision
      use eirmod_parmmod
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_CCONA
      USE EIRMOD_CGRPTL
      USE EIRMOD_CLOGAU
      USE EIRMOD_CPL3D
      USE EIRMOD_CPLOT
      USE EIRMOD_CINIT
      USE EIRMOD_COMSIG
      USE EIRMOD_CREF
      USE EIRMOD_CPOLYG
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CCOUPL
      USE EIRMOD_CGEOM
      USE EIRMOD_CSDVI
      USE EIRMOD_COMPRT
      USE EIRMOD_CPES
      USE EIRMOD_COMNNL
      USE EIRMOD_COMSOU
      USE EIRMOD_CSTEP
      USE EIRMOD_COMSPL
      USE EIRMOD_CTEXT
      USE EIRMOD_CLGIN
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
      USE EIRMOD_CTRIG
      USE EIRMOD_CTETRA
      USE EIRMOD_CESTIM
      USE EIRMOD_CUPD
      USE EIRMOD_PHOTON
      USE EIRMOD_JSON
      USE EIRMOD_PRESSURELOOP
#ifndef NO_JSON
      USE EIRMOD_IOUSR, ONLY: eirene_write_block_11_usr
      use json_module           ! IGNORE
#endif

      implicit none
      character(*) :: filename
#ifndef NO_JSON
      type(json_core) :: json
      type(json_value),pointer :: p, block_1, block_2, block_3a,
     .                            block_3b, block_45, block_6, block_7,
     .                            block_8, block_9, block_10, block_11,
     .                            block_12, block_13, block_14, block_0,
     .                            block_15
      type(json_value),pointer :: cmlines, crs_lines
      type(s_stack), pointer :: cur

! initialize the class
!      call json%initialize(.false.,.true.,.false.,'ES')
!      call json%initialize(compact_reals=.true.,real_format='ES',
!     .          case_sensitive_keys=.false.,compress_vectors=.true.)
      call json%initialize(compact_reals=.true.,real_format='ES',
     .          case_sensitive_keys=.false.,compress_vectors=.false.)

! initialize the structure:
      call json%create_object(p,'')

! add an "inputs" object to the structure:
      call json%create_object(block_0,'HEADER')
      call json%create_object(block_1,'GENERAL_DATA')
      call json%create_object(block_2,'STANDARD_MESH')
      call json%create_object(block_3a,'NONDEF_STD_SURFACES')
      call json%create_object(block_3b,'ADDITIONAL_SURFACES')
      call json%create_object(block_45,'PHYSICS_MODEL')
      call json%create_object(block_6,'REFLECTION_MODELS')
      call json%create_object(block_7,'SOURCES')
      call json%create_object(block_8,'SPECIFIC_ZONES')
      call json%create_object(block_9,'STATISTICS')
      call json%create_object(block_10,'ADDITIONAL_TALLIES')
      call json%create_object(block_11,'OUTPUT')
      call json%create_object(block_12,'DIAGNOSTICS')
      call json%create_object(block_13,'TIMEDEPENDENT_MODE')
      call json%create_object(block_14,'INTERFACING')
      call json%create_object(block_15,'MPI_INFORMATION')

      call json%add(p, block_0)     !add it to the root
      call json%add(p, block_1)     !add it to the root
      call json%add(p, block_2)     !add it to the root
      call json%add(p, block_3a)    !add it to the root
      call json%add(p, block_3b)    !add it to the root
      call json%add(p, block_45)    !add it to the root
      call json%add(p, block_6)     !add it to the root
      call json%add(p, block_7)     !add it to the root
      call json%add(p, block_8)     !add it to the root
      call json%add(p, block_9)     !add it to the root
      call json%add(p, block_10)    !add it to the root
      call json%add(p, block_11)    !add it to the root
      call json%add(p, block_12)    !add it to the root
      call json%add(p, block_13)    !add it to the root
      call json%add(p, block_14)    !add it to the root
      call json%add(p, block_15)    !add it to the root

      call json%add(block_0,'EIRENE_VERSION',eirene_version_string)
      call json%add(block_0,'TXTRUN',txtrun)
! COMMENT lines
      if (associated(cm_stack%head)) then
        call json%create_array(cmlines,'COMMENTS') !an array
        cur => cm_stack%head
        do while (associated(cur))
          call json%add(cmlines,'CM',trim(cur%string))
          cur => cur%next
        end do
        call json%add(block_0,cmlines)
      end if

      call eirene_write_block_1(block_1)
      call eirene_write_block_2(block_2)
      call eirene_write_block_3a(block_3a)
      call eirene_write_block_3b(block_3b)
      call eirene_write_block_45(block_45)
      call eirene_write_block_6(block_6)
      call eirene_write_block_7(block_7)
      call eirene_write_block_8(block_8)
      call eirene_write_block_9(block_9)
      call eirene_write_block_10(block_10)
      call eirene_write_block_11(block_11)
      call eirene_write_block_12(block_12)
      call eirene_write_block_13(block_13)
      call eirene_write_block_14(block_14)
      call eirene_write_block_15(block_15)

! write the file:
      call json%print(p,trim(filename))

!cleanup:
      call json%destroy(p)
      if (json%failed()) stop 1

      call eirene_dealloc_json_arrays
#endif
      return
#ifndef NO_JSON
      contains

!***************************************************************************************

      subroutine eirene_write_block_1(me)

      type(json_value),pointer :: me
      type(json_value),pointer :: cfile
      type(json_value),pointer :: var
      integer :: ifile, ntcpu_orig

      ntcpu_orig = ntcpu
      if (abs(mpts_comsou-1._dp) > eps10)
     .   ntcpu_orig = nint( real(ntcpu_orig,dp) / mpts_comsou )

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.1")
      call json%add(me,'NPRLL',nprll)
      call json%add(me,'NMODE',nmode)
      call json%add(me,'NTCPU',ntcpu_orig)
      call json%add(me,'NFILE',nfile)
      call json%add(me,'NITER0',iiter)
      call json%add(me,'NITER',niter)
      call json%add(me,'NTIME0',itimv)
      call json%add(me,'NTIME',ntime_in)
      call json%add(me,'NEVNTMAX',nevntmax/1000)
      call json%add(me,'NVNMXTR',nvnmxtr)

      call json%add(me,'NOPTIM',noptim_in)
      call json%add(me,'NOPTM1',noptm1)
      call json%add(me,'NGEOM_USR',ngeom_usr)
      call json%add(me,'NCOUP_INPUT',ncoup_input)
      call json%add(me,'NSMSTRA',nsmstra_in)
      call json%add(me,'NSTORAM',nstoram)
      call json%add(me,'NGSTAL',ngstal)
      call json%add(me,'NRTAL',nrtal_in)

      call json%add(me,'NLSCL',nlscl)
      call json%add(me,'NLTEST',nltest)
      call json%add(me,'NLANA',nlana)
      call json%add(me,'NLDRFT',nldrft)
      call json%add(me,'NLCRR',nlcrr)
      call json%add(me,'NLERG',nlerg)
      call json%add(me,'NLIDENT',nlident)
      call json%add(me,'NLONE',nlone)
      call json%add(me,'NLMOVIE',nlmovie)
      call json%add(me,'NLDFST',nldfst)
      call json%add(me,'NLRANMAR',nlranmar)
      call json%add(me,'NLCASCAD',nlcascad)
      call json%add(me,'NLOCTREE',nloctree)
      call json%add(me,'NLWRMSH',nlwrmsh)
      call json%add(me,'NEXVS',nexvs)
      call json%add(me,'NLTRIMESH',nltrimesh)
      call json%add(me,'NLSPCSCL',nlspcscl)
      call json%add(me,'NLSPCSCL_ON',nlspcscl_on)
      call json%add(me,'NLSOLEDGE',nlsoledge)
      call json%add(me,'NLSPCSPT',nlspcspt)

      call json%create_array(cfile,'CFILE') !an array
      do ifile=1, ndbnames
        if (.not.ldbread(ifile)) cycle
        call json%create_object(var,'')    !name does not matter
        call json%add(var, 'FILE',trim(dbhandle(ifile)))
        call json%add(var, 'PATH',trim(dbfname_in(ifile)))
        call json%add(cfile,var)
        nullify(var)
      end do
      call json%add(me,cfile)

      return
      end subroutine eirene_write_block_1

!******************************************************************************

      subroutine eirene_write_block_2(me)

      type(json_value),pointer :: me
      type(json_value),pointer :: rad, pol, tor, mlt, add, pnts, pnt,
     .                            plgs, plg
      integer :: k,i,j,nm

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.2")

      call json%add(me,'INDGRD',indgrd(1:3))
      call json%add(me,'NLRAD',nlrad)
      call json%add(me,'NLPOL',nlpol_in)
      call json%add(me,'NLTOR',nltor)
      call json%add(me,'NLMLT',nlmlt)
      call json%add(me,'NLADD',nladd)

! radial grid
      if (nlrad) then

        call json%create_object(rad,'RADIAL_GRID')
        call json%add(me,rad)

        call json%add(rad,'NLSLB',nlslb)
        call json%add(rad,'NLCRC',nlcrc)
        call json%add(rad,'NLELL',nlell)
        call json%add(rad,'NLTRI',nltri)
        call json%add(rad,'NLPLG',nlplg_in)
        call json%add(rad,'NLFEM',nlfem_in)
        call json%add(rad,'NLTET',nltet)
        call json%add(rad,'NLGEN',nlgen)

        if ((indgrd(1) <= 5) .and. nlplg .and.(plrefl > 0._dp)) then
          call json%add(rad,'NR1ST',nr1st_in)
        else
          call json%add(rad,'NR1ST',nr1st_in)
        end if
        call json%add(rad,'NRSEP',nrsep)
        call json%add(rad,'NRPLG',nrplg)
        call json%add(rad,'NPPLG',npplg)
        call json%add(rad,'NRKNOT',nrknot)
        call json%add(rad,'NCOOR',ncoor)

        if (indgrd(1) .le. 5) then
          if (nlslb .or. nlcrc .or. nlell .or. nltri) then
            call json%add(rad,'RIA',ria)
            call json%add(rad,'RGA',rga)
            call json%add(rad,'RAA',raa)
            call json%add(rad,'RRA',rra)
            if (nlell .or. nltri) then
              call json%add(rad,'EP1IN',ep1in)
              call json%add(rad,'EP1OT',ep1ot)
              call json%add(rad,'EP1CH',ep1ch)
              call json%add(rad,'EXEP1',exep1)
              call json%add(rad,'ELLIN',ellin)
              call json%add(rad,'ELLOT',ellot)
              call json%add(rad,'ELLCH',ellch)
              call json%add(rad,'EXELL',exell)
              if (nltri) then
                call json%add(rad,'TRIIN',triin)
                call json%add(rad,'TRIOT',triot)
                call json%add(rad,'TRICH',trich)
                call json%add(rad,'EXTRI',extri)
              end if
            end if
          end if
          if (nlplg) then
            call json%add(rad,'XPCOR',xpcor)
            call json%add(rad,'YPCOR',ypcor)
            call json%add(rad,'ZPCOR',zpcor)
            call json%add(rad,'PLREFL',plrefl)

            call json%create_array(pnts,'POLYGON_PARTS') !an array
            call json%add(rad,pnts)
            do k = 1, npplg
              call json%create_object(pnt,'') !name does not matter
              call json%add(pnt,'IPART',k)
              call json%add(pnt,'NPOINT',[npoint(1,k),npoint(2,k)])
              call json%add(pnts,pnt)
            end do

            nm = nr1st
            if (plrefl > 0._dp) nm = nm-1
            call json%create_array(plgs,'POLYGONS') !an array
            call json%add(rad,plgs)
            do i=1, nm
              call json%create_object(plg,'') !name does not matter
              call json%add(plg,'IPOLY',i)
              call json%add(plg,'XPOL_YPOL',
     .                          [(xpol(i,j),ypol(i,j),j=1,nrplg)])
              call json%add(plgs,plg)
              nullify(plg)
            end do
          end if
          if (nlfem_in .or. nltet) then
            call json%add(rad,'CASE',casename)
            call json%add(rad,'XPCOR',xpcor)
            call json%add(rad,'YPCOR',ypcor)
            call json%add(rad,'ZPCOR',zpcor)
          end if
        elseif (indgrd(1) .eq. 6) then
          if (nlslb .or. nlcrc .or. nlell .or. nltri) then
            call json%add(rad,'RIA',ria)
            call json%add(rad,'RGA',rga)
            call json%add(rad,'RAA',raa)
          elseif (nlplg .or. nlfem .or. nltet) then
            call json%add(rad,'XPCOR',xpcor)
            call json%add(rad,'YPCOR',ypcor)
            call json%add(rad,'ZPCOR',zpcor)
          end if
        end if  ! indgrd(1)
      end if  ! nlrad

! poloidal grid
      call json%create_object(pol,'POLOIDAL_GRID')
      call json%add(me,pol)

      call json%add(pol,'NLPLY',nlply)
      call json%add(pol,'NLPLA',nlpla)
      call json%add(pol,'NLPLP',nlplp)

      call json%add(pol,'NP2ND',np2nd_in)
      call json%add(pol,'NPSEP',npsep)
      call json%add(pol,'NPPLA',nppla)
      call json%add(pol,'NPPER',npper)

      if (indgrd(2) .le. 5) then
        call json%add(pol,'YIA',yia)
        call json%add(pol,'YGA',yga)
        call json%add(pol,'YAA',yaa)
        call json%add(pol,'YYA',yya)
      end if

! toroidal grid
      call json%create_object(tor,'TOROIDAL_GRID')
      call json%add(me,tor)

      call json%add(tor,'NLTRZ',nltrz)
      call json%add(tor,'NLTRA',nltra)
      call json%add(tor,'NLTRT',nltrt)

      call json%add(tor,'NT3RD',nt3rd_in)
      call json%add(tor,'NTSEP',ntsep)
      call json%add(tor,'NTTRA',nttra)
      call json%add(tor,'NTPER',ntper)

      if (indgrd(3) .le. 5) then
        call json%add(tor,'ZIA',zia)
        call json%add(tor,'ZGA',zga)
        call json%add(tor,'ZAA',zaa)
        call json%add(tor,'ZZA',zza)
        call json%add(tor,'ROA',roa)
      end if

! block multiplication
      if (nlmlt) then
        call json%create_object(mlt,'MESH_MULTIPLICATION')
        call json%add(me,mlt)
        call json%add(mlt,'NBMLT',nbmlt)
        call json%add(mlt,'VOLCOR',volcor(1:nbmlt))
      end if

! additional cells
      if (nladd) then
        call json%create_object(add,'ADD_CELLS')
        call json%add(me,add)
        call json%add(add,'NRADD',nradd)
        call json%add(add,'VOLADD',voladd(1:nradd))
      end if

      return
      end subroutine eirene_write_block_2

!******************************************************************************

      subroutine eirene_write_block_3a(me)

      type(json_value),pointer :: me
      type(json_value),pointer :: nondefs, srf
      integer :: ists, idimp, nlj, nst

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#subsection.2.3.1")

      if (((ntime >= 1) .and. (nprnli > 0)) .or. nlerg) then
        call json%add(me,'NSTSI',nstsi-1)
        nst = nstsi-1
      else
        call json%add(me,'NSTSI',nstsi)
        nst = nstsi
      end if
      call json%create_array(nondefs,'SURFACES') !an array
      do ists=1, nst
        nlj = nlim + ists
        idimp = maxloc(inump(ists,1:3),dim=1)
        call json%create_object(srf,'')    !name does not matter
        call json%add(srf,'TXTSFL',trim(txtsfl(nlj)))
        call json%add(srf,'IDIMP',idimp)
        call json%add(srf,'INUMP',inump(ists,idimp))
        call json%add(srf,'IRPTA1',irpta(ists,1))
        call json%add(srf,'IRPTE1',irpte(ists,1))
        call json%add(srf,'IRPTA2',irpta(ists,2))
        call json%add(srf,'IRPTE2',irpte(ists,2))
        call json%add(srf,'IRPTA3',irpta(ists,3))
        call json%add(srf,'IRPTE3',irpte(ists,3))

        call eirene_write_surf_switches (srf, nlj)
        call eirene_write_ref_model (srf, nlj)

        call json%add(nondefs,srf)
        nullify(srf)
      end do
      call json%add(me,nondefs)

      return
      end subroutine eirene_write_block_3a

!******************************************************************************

      subroutine eirene_write_surf_switches (srf,nlj)

      type(json_value),pointer :: srf
      integer, intent(in) :: nlj

      call json%add(srf,'ILIIN',iliin(nlj))
      call json%add(srf,'ILSIDE',ilside(nlj))
      call json%add(srf,'ILSWCH',ilswch(nlj))
      call json%add(srf,'ILEQUI',ilequi(nlj))
      call json%add(srf,'ILTOR',iltor(nlj))
      call json%add(srf,'ILCOL',ilcol_in(nlj))
      call json%add(srf,'ILFIT',ilfit(nlj))
      call json%add(srf,'ILCELL',ilcell(nlj))
      call json%add(srf,'ILBOX',ilbox(nlj))
      call json%add(srf,'ILPLG',ilplg(nlj))

      return
      end subroutine eirene_write_surf_switches

!******************************************************************************

      subroutine eirene_write_ref_model (srf,nlj)

      type(json_value),pointer :: srf
      integer, intent(in) :: nlj
      real(dp) :: znmlc

      if (len_trim(smod_name(nlj)) > 0) then

         call json%add(srf,'SURFMOD','SURFMOD_'//trim(smod_name(nlj)))

      else

         znmlc = 100*znml(nlj) + zncl(nlj)
         call json%add(srf,'ILREF',ilref(1,nlj))
         call json%add(srf,'ILSPT',ilspt(nlj))
         call json%add(srf,'ISRS',isrs(1,nlj))
         call json%add(srf,'ISRC',isrc(1,nlj))
         call json%add(srf,'LCHSPNWL',lchspnwl_in(1,nlj))
         call json%add(srf,'ZNML',znmlc)
         call json%add(srf,'EWALL',ewall(nlj))
         call json%add(srf,'EWBIN',ewbin(nlj))
         call json%add(srf,'TRANSP',transp(1,1:2,nlj))
         call json%add(srf,'FSHEAT',fsheat(nlj))
         call json%add(srf,'RECYCF',recycf(1,nlj))
         call json%add(srf,'RECYCT',recyct(1,nlj))
         call json%add(srf,'RECPRM',recprm(1,nlj))
         call json%add(srf,'EXPPL',exppl(1,nlj))
         call json%add(srf,'EXPEL',expel(1,nlj))
         call json%add(srf,'EXPIL',expil(1,nlj))
         call json%add(srf,'RECYCS',recycs(1,nlj))
         call json%add(srf,'RECYCC',recycc(1,nlj))
         call json%add(srf,'SPTPRM',sptprm(1,nlj))
         call json%add(srf,'ESPUTS',esputs(1,nlj))
         call json%add(srf,'ESPUTC',esputc(1,nlj))
      end if

      return
      end subroutine eirene_write_ref_model

!******************************************************************************

      subroutine eirene_write_block_3b(me)

      type(json_value),pointer :: me
      type(json_value),pointer :: ext, blk
      type(json_value),pointer :: adds, srf, eqlin, eqscn, trns,
     .                            lin, scn, trn
      type(json_value),pointer :: ch0lines, ch1lines, ch2lines
      type(s_stack), pointer :: cur
      type(transform), pointer :: trf
      integer :: i, j

      if (nlimi == 0) then
        call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#subsection.2.3.2")
        call json%add(me,'NLIMI',nlimi)
        return
      end if

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#subsection.2.3.2")
      call json%add(me,'INCLUDE',
     .     "eirene_add_surfaces.input.json")

! write input for additional surfaces onto separate file
! initialize the structure:
      call json%create_object(ext,'')
      call json%create_object(blk,'ADDITIONAL_SURFACES')
      call json%add(ext,blk)

      call json%add(blk,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#subsection.2.3.2")

      call json%add(blk,'NLIMI',nlimi)

! CH0 lines
      if (associated(ch0_stack%head)) then
        call json%create_array(ch0lines,'CH0-LINES') !an array
        cur => ch0_stack%head
        do while (associated(cur))
          call json%add(ch0lines,'CH0',cur%string)
          cur => cur%next
        end do
        call json%add(blk,ch0lines)
      end if

! loop over additional surfaces
      call json%create_array(adds,'SURFACES') !an array

      do i=1, nlimi

        call json%create_object(srf,'') !name does not matter
        call json%add(srf,'ILIM',i)
        call json%add(srf,'TXTSFL',trim(txtsfl(i)))

! CH1 lines
        if (associated(ch1_stack(i)%head)) then
          call json%create_array(ch1lines,'CH1-LINES') !an array
          cur => ch1_stack(i)%head
          do while (associated(cur))
            call json%add(ch1lines,'CH1',cur%string)
            cur => cur%next
          end do
          call json%add(srf,ch1lines)
        end if

! CH2 lines
        if (associated(ch2_stack(i)%head)) then
          call json%create_array(ch2lines,'CH2-LINES') !an array
          cur => ch2_stack(i)%head
          do while (associated(cur))
            call json%add(ch2lines,'CH2',cur%string)
            cur => cur%next
          end do
          call json%add(srf,ch2lines)
        end if

        call json%add(srf,'RLB',rlb_in(i))
        call json%add(srf,'SAREA',sarea(i))
        call json%add(srf,'RLWMN',rlwmn(i))
        call json%add(srf,'RLWMX',rlwmx(i))

        call eirene_write_surf_switches (srf, i)

! geometrical information
! surface definition
        if (rlb_in(i) < 2.) then
          call json%add(srf,'A0LM',a0lm_in(i))
          call json%add(srf,'A1LM',a1lm_in(i))
          call json%add(srf,'A2LM',a2lm_in(i))
          call json%add(srf,'A3LM',a3lm_in(i))
          call json%add(srf,'A4LM',a4lm_in(i))
          call json%add(srf,'A5LM',a5lm_in(i))
          call json%add(srf,'A6LM',a6lm_in(i))
          call json%add(srf,'A7LM',a7lm_in(i))
          call json%add(srf,'A8LM',a8lm_in(i))
          call json%add(srf,'A9LM',a9lm_in(i))
        else if (rlb_in(i) < 7.) then
! 2 point surface: straight line
          if (rlb_in(i) >= 2.) then
            call json%add(srf,'P1',p1_in(1:3,i))
            call json%add(srf,'P2',p2_in(1:3,i))
          end if
! 3 point surface: triangle
          if (rlb_in(i) >= 3.) then
            call json%add(srf,'P3',p3_in(1:3,i))
          end if
! 4 point surface: flat quadrangle
          if (rlb_in(i) >= 4.) then
            call json%add(srf,'P4',p4_in(1:3,i))
          end if
! 5 point surface
          if (rlb_in(i) >= 5.) then
            call json%add(srf,'P5',p5_in(1:3,i))
          end if
! 6 point surface
          if (rlb_in(i) >= 6.) then
            call json%add(srf,'P6',p6_in(1:3,i))
          end if
        end if

! boundary specification
        if ((rlb_in(i) > 0.) .and. (rlb_in(i) < 2.)) then
          call json%add(srf,'XLIMS1',xlims1_in(1,i))
          call json%add(srf,'YLIMS1',ylims1_in(1,i))
          call json%add(srf,'ZLIMS1',zlims1_in(1,i))
          call json%add(srf,'XLIMS2',xlims2_in(1,i))
          call json%add(srf,'YLIMS2',ylims2_in(1,i))
          call json%add(srf,'ZLIMS2',zlims2_in(1,i))
        else if (rlb_in(i) <= 0.) then

          if (ilin(i) > 0) then
            call json%create_array(eqlin,'LIN_BOUNDS') !an array
            do j=1, ilin(i)
              call json%create_object(lin,'') !name does not matter
              call json%add(lin,'ALIMS',alims_in(j,i))
              call json%add(lin,'XLIMS',xlims_in(j,i))
              call json%add(lin,'YLIMS',ylims_in(j,i))
              call json%add(lin,'ZLIMS',zlims_in(j,i))
              call json%add(eqlin,lin)
              nullify(lin)
            end do
            call json%add(srf,eqlin)
          end if

          if (iscn(i) > 0) then
            call json%create_array(eqscn,'SCN_BOUNDS') !an array
            do j=1, iscn(i)
              call json%create_object(scn,'') !name does not matter
              call json%add(scn,'ALIMS0',alims0_in(j,i))
              call json%add(scn,'XLIMS1',xlims1_in(j,i))
              call json%add(scn,'YLIMS1',ylims1_in(j,i))
              call json%add(scn,'ZLIMS1',zlims1_in(j,i))
              call json%add(scn,'XLIMS2',xlims2_in(j,i))
              call json%add(scn,'YLIMS2',ylims2_in(j,i))
              call json%add(scn,'ZLIMS2',zlims2_in(j,i))
              call json%add(scn,'XLIMS3',xlims3_in(j,i))
              call json%add(scn,'YLIMS3',ylims3_in(j,i))
              call json%add(scn,'ZLIMS3',zlims3_in(j,i))
              call json%add(eqscn,scn)
              nullify(scn)
            end do
            call json%add(srf,eqscn)
          end if

        end if

        call eirene_write_ref_model (srf, i)

! transform
        if (associated(trnsfrm(i)%head)) then
          call json%create_array(trns,'TRANSFORMATION') !an array
          trf => trnsfrm(i)%head
          do while (associated(trf))
            call json%create_object(trn,'') !name does not matter
            call json%add(trn,'ITINI',trf%itini)
            call json%add(trn,'ITEND',trf%itend)
            call json%add(trn,'XLCOR',trf%xlcor)
            call json%add(trn,'YLCOR',trf%ylcor)
            call json%add(trn,'ZLCOR',trf%zlcor)
            call json%add(trn,'XLREF',trf%xlref)
            call json%add(trn,'YLREF',trf%ylref)
            call json%add(trn,'ZLREF',trf%zlref)
            call json%add(trn,'XLROT',trf%xlrot)
            call json%add(trn,'YLROT',trf%ylrot)
            call json%add(trn,'ZLROT',trf%zlrot)
            call json%add(trn,'ALROT',trf%alrot)
            call json%add(trns,trn)
            nullify(trn)
            trf => trf%next
          end do
          call json%add(srf,trns)
        end if

        call json%add(adds,srf)
        nullify(srf)
        nullify(trns)
      end do
      call json%add(blk,adds)

! write the file:
      call json%print(ext,'eirene_add_surfaces.input.json')

!cleanup:
      call json%destroy(ext)
      if (json%failed()) stop 1

      return
      end subroutine eirene_write_block_3b

!******************************************************************************

      subroutine eirene_write_block_45(me)

      type(json_value),pointer :: me
      type(json_value),pointer :: ext, atoms, mols, ions, phots
      type(json_value),pointer :: reacs, rea, blk4, blk5, ssapm
      type(json_value),pointer :: blk4a, blk4b, blk4c, blk4d
      type(json_value),pointer :: bulks, blk5a, blk5b, blk
      type(json_value),pointer :: te, tis, ti, dis, di, vls, vl, bf, v,
     .                            op, ops, zi, zis
      integer :: i, ndum1(1), ndum2(1), ndum3(1), ndum4(1), nti, nv

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#subsection.2.4.1")
      call json%add(me,'INCLUDE',
     .     "eirene_physics_model.input.json")

! initialize the structure:
      call json%create_object(ext,'')
      call json%create_object(blk,'PHYSICS_MODEL')
      call json%add(ext,blk)

      call json%add(blk,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#subsection.2.4.1")

      call json%create_object(blk4,'SPECIES_SPEC')
      call json%create_object(blk5,'BACKGROUND')
      call json%add(blk5,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#subsection.2.5.1")

      call json%add(blk,blk4)
      call json%add(blk,blk5)

! start block 4
      call json%create_object(ssapm,'REACTIONS')
      call json%add(blk4,ssapm)

      call json%add(ssapm,'NREACI',nreaci)

! write classical reaction lines for easy comparison with fixed input
      if (associated(crs_stack%head)) then ! an array
        call json%create_array(crs_lines,'ClassicReactionStrings')
        cur => crs_stack%head
        do while (associated(cur))
          call json%add(crs_lines,'CRS',trim(cur%string))
          cur => cur%next
        end do
        call json%add(ssapm,crs_lines)
      end if

! loop over reaction lines
      call json%create_array(reacs,'REAC_SPECS') !an array

      do i=1, irlines

        call json%create_object(rea,'') !name does not matter

        call json%add(rea,'IR',reaclines(i)%no)
        call json%add(rea,'FILNAM',reaclines(i)%FILE)
        call json%add(rea,'H123',reaclines(i)%h_select)
        call json%add(rea,'REAC',reaclines(i)%reac_string)
        call json%add(rea,'CRC',reaclines(i)%reactyp)
        call json%add(rea,'MASSP',reaclines(i)%mp)
        call json%add(rea,'MASST',reaclines(i)%mt)
        call json%add(rea,'DP',reaclines(i)%dpp)
        call json%add(rea,'R1MN',reaclines(i)%r1mn)
        call json%add(rea,'R1MX',reaclines(i)%r1mx)
        call json%add(rea,'R2MN',reaclines(i)%r2mn)
        call json%add(rea,'R2MX',reaclines(i)%r2mx)
        call json%add(rea,'JFEX1MN',reaclines(i)%jfex1mn)
        call json%add(rea,'JFEX1MX',reaclines(i)%jfex1mx)
        call json%add(rea,'JFEX2MN',reaclines(i)%jfex2mn)
        call json%add(rea,'JFEX2MX',reaclines(i)%jfex2mx)
        call json%add(rea,'FP1',reaclines(i)%fp1)
        call json%add(rea,'FP2',reaclines(i)%fp2)

        if (reaclines(i)%FILE == 'CONST') then
          call json%add(rea,'IFTFLG',reaclines(i)%iftflg)
          call json%add(rea,'NCOEF',reaclines(i)%ncoef)
          call json%add(rea,'COEF',
     .                  reaclines(i)%coef(1:reaclines(i)%ncoef))
        end if

        if ((reaclines(i)%FILE == 'ADAS') .OR.
     .      (reaclines(i)%FILE == 'TAB2D')) then
          call json%add(rea,'ELNAME',reaclines(i)%element)
          call json%add(rea,'IZ',reaclines(i)%iz)
          call json%add(rea,'BUNDLING',trim(reaclines(i)%bundling))
        end if

        if (index(reaclines(i)%FILE,'CR') /= 0) then
          call json%add(rea,'M_POPESC',reaclines(i)%m_popesc)
          if (reaclines(i)%m_popesc > 0) then
           call json%add(rea,'IROW_ESC',
     .      reaclines(i)%irow_esc(1:reaclines(i)%m_popesc))
           call json%add(rea,'ICOL_ESC',
     .      reaclines(i)%icol_esc(1:reaclines(i)%m_popesc))
           call json%add(rea,'POP_ESC',
     .      reaclines(i)%pop_esc(1:reaclines(i)%m_popesc))
          end if
          call json%add(rea,'M_QEXT',reaclines(i)%m_qext)
           if (reaclines(i)%m_qext > 0) then
             call json%add(rea,'IROW_EXT',
     .      reaclines(i)%irow_ext(1:reaclines(i)%m_qext))
             call json%add(rea,'ITAL_EXT_1',
     .      reaclines(i)%ital_ext(1,1:reaclines(i)%m_qext))
             call json%add(rea,'ITAL_EXT_2',
     .      reaclines(i)%ital_ext(2,1:reaclines(i)%m_qext))
             call json%add(rea,'ITAL_EXT_3',
     .      reaclines(i)%ital_ext(3,1:reaclines(i)%m_qext))
           end if
        end if

        call json%add(reacs,rea)

      end do

      call json%add(ssapm,reacs)

! atoms
      ityp = 1
      call json%create_object(atoms,'ATOMS')
      call json%add(atoms,'NATMI',natmi)
      call json%create_array(blk4a,'SPECIES') !an array
      call eirene_write_block_4abcd (blk4a, ityp, nsph, natmi, natm,
     .     texts, 'IATM', 'A', nmassa, nchara, ndum1, ndum2, isrf, isrt,
     .     nrca, nfola, ngena, nhsts, ireaca, ibulka,
     .     iscd1a, iscd2a, iscd3a, iscd4a, iscdea, iestma, ibgka,
     .     eeleca, ebulka, escd1a, freaca, edpota)
      call json%add(atoms,blk4a)
      call json%add(blk4,atoms)

! molecules
      ityp = 2
      call json%create_object(mols,'MOLECULES') !an array
      call json%add(mols,'NMOLI',nmoli)
      call json%create_array(blk4b,'SPECIES') !an array
      call eirene_write_block_4abcd (blk4b, ityp, nspa, nmoli, nmol,
     .     texts, 'IMOL', 'M', nmassm, ncharm, nprt, ndum2, isrf, isrt,
     .     nrcm, nfolm, ngenm, nhsts, ireacm, ibulkm,
     .     iscd1m, iscd2m, iscd3m, iscd4m, iscdem, iestmm, ibgkm,
     .     eelecm, ebulkm, escd1m, freacm, edpotm, LKIND=lkindm)
      call json%add(mols,blk4b)
      call json%add(blk4,mols)

! test ions
      ityp = 3
      call json%create_object(ions,'TEST_IONS') !an array
      call json%add(ions,'NIONI',nioni)
      call json%create_array(blk4c,'SPECIES') !an array
      call eirene_write_block_4abcd (blk4c, ityp, nspam, nioni, nion,
     .     texts, 'IION', 'I', nmassi, nchari, nprt, nchrgi, isrf, isrt,
     .     nrci, nfoli, ngeni, nhsts, ireaci, ibulki,
     .     iscd1i, iscd2i, iscd3i, iscd4i, iscdei, iestmi, ibgki,
     .     eeleci, ebulki, escd1i, freaci, edpoti, LKIND=lkindi)
      call json%add(ions,blk4c)
      call json%add(blk4,ions)

! photons
      ityp = 0
      call json%create_object(phots,'PHOTONS') !an array
      call json%add(phots,'NPHOTI',nphoti)
      call json%create_array(blk4d,'SPECIES') !an array
      call eirene_write_block_4abcd (blk4d, ityp, 0, nphoti, nphot,
     .     texts, 'IPHOT', 'PH', ndum1, ndum2, ndum3, ndum4, isrf, isrt,
     .     nrcph, nfolph, ngenph, nhsts, ireacph, ibulkph,
     .     iscd1ph, iscd2ph, iscd3ph, iscd4ph, iscdeph, iestmph, ibgkph,
     .     eelecph, ebulkph, escd1ph, freacph, edpotph)
      call json%add(phots,blk4d)
      call json%add(blk4,phots)

! start block 5
! bulk ions
      ityp = 4
      call json%create_object(bulks,'BULK_IONS')
      call json%add(bulks,'NPLSI',nplsi)
      call json%create_array(blk5a,'SPECIES') !an array
      call eirene_write_block_4abcd (blk5a, ityp, nspami, nplsi, npls,
     .     texts, 'IPLS', 'P', nmassp, ncharp, nprt, nchrgp, isrf, isrt,
     .     nrcp, ndum1, ndum2, nhsts, ireacp, ibulkp,
     .     iscd1p, iscd2p, iscd3p, iscd4p, iscdep, ndum3, ndum4,
     .     eelecp, ebulkp, escd1p, freacp, edpotp,
     .     DENSLIM=denslim, CDENMODEL=cdenmodel)
      call json%add(bulks,blk5a)

      call json%add(blk5,bulks)

! plasma background
      call json%create_object(blk5b,'PLASMA')

      call json%add(blk5b,'INDPRO',indpro_in)

! Te
      if (indpro(1) <= 5) then
        call json%create_object(te,'TE')
        call json%add(te,'TE0',te0)
        call json%add(te,'TE1',te1)
        call json%add(te,'TE2',te2)
        call json%add(te,'TE3',te3)
        call json%add(te,'TE4',te4)
        call json%add(te,'TE5',te5)
        call json%add(blk5b,te)
      end if

! Ti
      if (indpro(2) <= 5) then
        call json%create_array(tis,'TI') !an array

        nti = 1
        if (mod(iabs(indpro(2)),100) > 9) nti=nplsI

        do i = 1, nti
          call json%create_object(ti,'')

          call json%add(ti,'TI0',ti0(i))
          call json%add(ti,'TI1',ti1(i))
          call json%add(ti,'TI2',ti2(i))
          call json%add(ti,'TI3',ti3(i))
          call json%add(ti,'TI4',ti4(i))
          call json%add(ti,'TI5',ti5(i))

          call json%add(tis,ti)
        end do
        call json%add(blk5b,tis)
      end if
! Ni
      if (indpro(3) <= 5) then
        call json%create_array(dis,'DI') !an array

        do i = 1, nplsi
          call json%create_object(di,'')

          call json%add(di,'DI0',di0(i))
          call json%add(di,'DI1',di1(i))
          call json%add(di,'DI2',di2(i))
          call json%add(di,'DI3',di3(i))
          call json%add(di,'DI4',di4(i))
          call json%add(di,'DI5',di5(i))

          call json%add(dis,di)
        end do
        call json%add(blk5b,dis)
      end if

! VX, VY, VZ
      if (indpro(4) <= 5) then
        call json%create_array(vls,'VEL') !an array

        nv = 1
        if (nlmlv) nv=nplsI

        do i = 1, nv
          call json%create_object(vl,'')

          call json%add(vl,'VX0',vx0(i))
          call json%add(vl,'VX1',vx1(i))
          call json%add(vl,'VX2',vx2(i))
          call json%add(vl,'VX3',vx3(i))
          call json%add(vl,'VX4',vx4(i))
          call json%add(vl,'VX5',vx5(i))

          call json%add(vl,'VY0',vy0(i))
          call json%add(vl,'VY1',vy1(i))
          call json%add(vl,'VY2',vy2(i))
          call json%add(vl,'VY3',vy3(i))
          call json%add(vl,'VY4',vy4(i))
          call json%add(vl,'VY5',vy5(i))

          call json%add(vl,'VZ0',vz0(i))
          call json%add(vl,'VZ1',vz1(i))
          call json%add(vl,'VZ2',vz2(i))
          call json%add(vl,'VZ3',vz3(i))
          call json%add(vl,'VZ4',vz4(i))
          call json%add(vl,'VZ5',vz5(i))

          call json%add(vls,vl)
        end do
        call json%add(blk5b,vls)
      end if

! B field
      if (indpro(5) <= 5) then
        call json%create_object(bf,'B_FIELD')
        call json%add(bf,'B0',b0)
        call json%add(bf,'B1',b1)
        call json%add(bf,'B2',b2)
        call json%add(bf,'B3',b3)
        call json%add(bf,'B4',b4)
        call json%add(bf,'B5',b5)
        call json%add(blk5b,bf)
      end if

! Zi
      if ((indpro(11) > 0) .and. (indpro(11) <= 5)) then
        call json%create_array(zis,'ZI') !an array

        do i = 1, nplsi
          call json%create_object(zi,'')

          call json%add(zi,'ZI0',zi0(i))
          call json%add(zi,'ZI1',zi1(i))
          call json%add(zi,'ZI2',zi2(i))
          call json%add(zi,'ZI3',zi3(i))
          call json%add(zi,'ZI4',zi4(i))
          call json%add(zi,'ZI5',zi5(i))

          call json%add(zis,zi)
        end do

        call json%add(blk5b,zis)
      end if

! VOL
      if (indpro(12) <= 5) then
        call json%create_object(v,'VOL')
        call json%add(v,'VL0',vl0)
        call json%add(v,'VL1',vl1)
        call json%add(v,'VL2',vl2)
        call json%add(v,'VL3',vl3)
        call json%add(v,'VL4',vl4)
        call json%add(v,'VL5',vl5)
        call json%add(blk5b,v)
      end if

! OPTIONAL INPUT FOR INPUT TALLIES
      if (any(intlopts /= 0)) then
        call json%create_array(ops,'INPUT_TALLY_OPTIONS') !an array
        do i = 1, ntali
          if (intlopts(i) == 0) cycle
          call json%create_object(op,'')
          call json%add(op,'ITAL',-i)
          call json%add(op,'IOPT',intlopts(i))
          call json%add(ops,op)
        end do
        call json%add(blk5b,ops)
      end if

      call json%add(blk5,blk5b)
! write the file:
      call json%print(ext,'eirene_physics_model.input.json')

!cleanup:
      call json%destroy(ext)
      if (json%failed()) stop 1

      return
      end subroutine eirene_write_block_45

!******************************************************************************

      subroutine eirene_write_block_4abcd (me, ityp, nbas, nloop, ndim,
     .     texts, cndx, cext, nmass, nchar, nprt, nchrg, isrf, isrt,
     .     nrc, nfol, ngen, nhsts, ireac, ibulk,
     .     iscd1, iscd2, iscd3, iscd4, iscde, iestm, ibgk,
     .     eelec, ebulk, escd1, freac, edpot, lkind, denslim, cdenmodel)

      type(json_value),pointer :: me
      type(json_value),pointer :: elem, reas, re, model
      integer, intent(in) :: ityp, nbas, nloop, ndim
      integer, intent(in) :: nmass(*), nchar(*),
     .         isrf(nspz,*), isrt(nspz,*), nprt(*), nchrg(*),
     .         nrc(*), nfol(*), ngen(*), nhsts(*),
     .         ireac(ndim,*), ibulk(ndim,*), iscd1(ndim,*),
     .         iscd2(ndim,*), iscd3(ndim,*), iscd4(ndim,*),
     .         iscde(ndim,*), iestm(ndim,*), ibgk(ndim,*)
      integer, optional :: lkind(ndim)
      real(dp), optional :: denslim(ndim)
      real(dp), intent(in) :: eelec(ndim,*), ebulk(ndim,*),
     .         escd1(ndim,*), freac(ndim,*), edpot(ndim,*)
      character(*) :: texts(*), cndx, cext
      character(*), optional :: cdenmodel(ndim)
      integer :: i, ispz, numsec, k, j, jdens
      logical :: lden

      lden = present(cdenmodel)

      do i = 1, nloop

        ispz = nbas + i
        numsec = 0
        if (any(iscd1(i,1:nrc(i)) > 0)) numsec = numsec + 1
        if (any(iscd2(i,1:nrc(i)) > 0)) numsec = numsec + 1
        if (any(iscd3(i,1:nrc(i)) > 0)) numsec = numsec + 1
        if (any(iscd4(i,1:nrc(i)) > 0)) numsec = numsec + 1

        call json%create_object(elem,'') !name does not matter

        call json%add(elem,cndx,i)
        call json%add(elem,'SPECIES',texts(ispz))
        call json%add(elem,'NMASS'//cext,nmass(i))
        call json%add(elem,'NCHAR'//cext,nchar(i))

        if (ityp >= 2) then
          call json%add(elem,'NPRT',nprt(ispz))
        end if

        if (ityp >= 3) then
          call json%add(elem,'NCHRG'//cext,nchrg(i))
        end if

        call json%add(elem,'ISRF',isrf(ispz,1))
        call json%add(elem,'ISRT',isrt(ispz,1))
        call json%add(elem,'NUMSEC',numsec)
        call json%add(elem,'NRC'//cext,nrc(i))
        if (ityp < 4) then
          call json%add(elem,'NFOL'//cext,nfol(i))
          call json%add(elem,'NGEN'//cext,ngen(i))
        end if
        call json%add(elem,'NHSTS',nhsts(ispz))

        if (present(lkind)) then
          call json%add(elem,'LKIND'//cext,lkind(i))
        end if

        if (present(denslim)) then
          if (denslim(i) > 9.e29) then
            jdens = 0
          else
            jdens = nint(log10(denslim(i))/log10(2._dp))
          end if
          call json%add(elem,'DENSLIM',jdens)
        end if

        call json%create_array(reas,'REACTIONS') !an array

        do k = 1, nrc(i)

          call json%create_object(re,'') !name does not matter

          call json%add(re,'IREAC'//cext,ireac(i,k))
          call json%add(re,'IBULK'//cext,ibulk(i,k))
          call json%add(re,'ISCD1'//cext,iscd1(i,k))
          call json%add(re,'ISCD2'//cext,iscd2(i,k))
          if (numsec > 2)
     .      call json%add(re,'ISCD3'//cext,iscd3(i,k))
          if (numsec > 3)
     .      call json%add(re,'ISCD4'//cext,iscd4(i,k))
          call json%add(re,'ISCDE'//cext,iscde(i,k))

          if (ityp < 4) then
            call json%add(re,'IESTM'//cext,iestm(i,k))
            call json%add(re,'IBGK'//cext,ibgk(i,k))
          end if

          call json%add(re,'EELEC'//cext,eelec(i,k))
          call json%add(re,'EBULK'//cext,ebulk(i,k))
          call json%add(re,'ESCD1'//cext,escd1(i,k))
          call json%add(re,'FREAC'//cext,freac(i,k))

          if (ityp == 1)
     .      call json%add(re,'EDPOT'//cext,edpot(i,k))

          call json%add(reas,re)

         end do

         if ((ityp == 4) .and. lden) then
           if (len_trim(cdenmodel(i)) > 0) then

             call json%create_object(model,'DENS_MODEL')
             call json%add(model,'NRE',TDMPAR(I)%TDM%NRE)
             SELECT CASE (CDENMODEL(I))
             CASE (FORT//'13')
               call json%add(model,'MODEL',FORT//'13')
               call json%add(model,'ISP',TDMPAR(I)%TDM%ISP(1))
             CASE (FORT//'10')
               call json%add(model,'MODEL',FORT//'10')
               call json%add(model,'ISP',TDMPAR(I)%TDM%ISP(1))
               call json%add(model,'ITP',TDMPAR(I)%TDM%ITP(1))
               call json%add(model,'ISTR',TDMPAR(I)%TDM%ISTR(1))
             CASE ('CONSTANT  ')
               call json%add(model,'MODEL','CONSTANT')
               call json%add(model,'TVAL',TDMPAR(I)%TDM%TVAL)
               call json%add(model,'DVAL',TDMPAR(I)%TDM%DVAL)
               call json%add(model,'VXVAL',TDMPAR(I)%TDM%VXVAL)
               call json%add(model,'VYVAL',TDMPAR(I)%TDM%VYVAL)
               call json%add(model,'VZVAL',TDMPAR(I)%TDM%VZVAL)
             CASE ('MULTIPLY  ')
               call json%add(model,'MODEL','MULTIPLY')
               call json%add(model,'ISP',TDMPAR(I)%TDM%ISP(1))
               call json%add(model,'ITP',TDMPAR(I)%TDM%ITP(1))
               call json%add(model,'ISTR',TDMPAR(I)%TDM%ISTR(1))
               call json%add(model,'DFACTOR',TDMPAR(I)%TDM%DFACTOR)
               call json%add(model,'TFACTOR',TDMPAR(I)%TDM%TFACTOR)
               call json%add(model,'VFACTOR',TDMPAR(I)%TDM%VFACTOR)
             CASE ('SAHA      ')
!PB   TO BE WRITTEN
               call json%add(model,'MODEL','SAHA')
             CASE ('BOLTZMANN ')
               call json%add(model,'MODEL','MULTIPLY')
               call json%add(model,'ISP',TDMPAR(I)%TDM%ISP(1))
               call json%add(model,'ITP',TDMPAR(I)%TDM%ITP(1))
               call json%add(model,'ISTR',TDMPAR(I)%TDM%ISTR(1))
               call json%add(model,'G_BOLTZ',TDMPAR(I)%TDM%G_BOLTZ)
               call json%add(model,'DELTAE',TDMPAR(I)%TDM%DELTAE)
             CASE ('CORONA    ')
               call json%add(model,'MODEL','CORONA')
               call json%add(model,'ISP',TDMPAR(I)%TDM%ISP(1))
               call json%add(model,'ITP',TDMPAR(I)%TDM%ITP(1))
               call json%add(model,'ISTR',TDMPAR(I)%TDM%ISTR(1))
               call json%add(model,'IRC',TDMPAR(I)%TDM%IRC(1))
               call json%add(model,'A_CORONA',TDMPAR(I)%TDM%A_CORONA)
             CASE ('COLRAD    ')
               call json%add(model,'MODEL','COLRAD')
               DO J=1, TDMPAR(I)%TDM%NRE
                 call json%add(model,'ISP',TDMPAR(I)%TDM%ISP(J))
                 call json%add(model,'ITP',TDMPAR(I)%TDM%ITP(J))
                 call json%add(model,'ISTR',TDMPAR(I)%TDM%ISTR(J))
                 call json%add(model,'IRC',TDMPAR(I)%TDM%IRC(J))
               END DO
             CASE DEFAULT
!PB  NOTHING TO BE DONE
             END SELECT
             call json%add(elem, model)
           end if
         end if

         call json%add(elem, reas)

         call json%add(me,elem)

      end do

      return
      end subroutine eirene_write_block_4abcd

!******************************************************************************


      subroutine eirene_write_block_6(me)

      type(json_value),pointer :: me
      type(json_value),pointer :: vrs, vr, smod, sm
      TYPE(REFMODEL), POINTER :: refcur
      type(spec_ref), pointer :: spr
      character(len=:), allocatable :: tps(:)
      character(1) :: slash
      integer :: i1, i, max_len, j

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#subsection.2.6.1")

      call json%add(me,'NLTRIM',nltrim)

      if (nltrim) then
        if (len_trim(reffil(1)) > 0) then
          slash = char(47)
          i1 = scan(reffil(1),slash,.true.)
          if (i1 > 0) then

            call json%add(me,'CPATH',reffil(1)(1:i1))
            max_len = 0
            do i=1,nhd6
              max_len=max(max_len,len_trim(reffil(i)(i1+1:)))
            end do

            allocate (character(len=max_len) :: tps(nhd6))
            tps = repeat(' ',max_len)
!           do i = 1, nhd6
            j = 0
            do i = nhd6, 1, -1
              j = j + 1
              tps(j) = trim(reffil(i)(i1+1:))
            end do

            call json%add(me,'PROJ_ON_MAT',tps)
            deallocate (tps)

          end if
        end if
      end if

      call json%add(me,'DATD',datd_in(1:natmi))
      call json%add(me,'DMLD',dmld_in(1:nmoli))
      call json%add(me,'DIOD',diod_in(1:nioni))
      call json%add(me,'DPLD',dpld_in(1:nplsi))
      if (nphoti > 0)
     .  call json%add(me,'DPHD',dphd_in(1:nphoti))

      call json%add(me,'ERMIN',ermin)
      call json%add(me,'ERCUT',ercut)
      call json%add(me,'RPROB0',rprob0)
      call json%add(me,'RINTEG',rinteg(1))
      call json%add(me,'EINTEG',einteg(1))
      call json%add(me,'AINTEG',ainteg(1))

      refcur => reflist
      call json%create_array(smod,'SURFMODS')
      call json%add(me,smod)

      do while (associated(refcur))
        call json%create_object(sm,'')
        call json%add(sm,'NAME','SURFMOD_'//REFCUR%REFNAME)
        call json%add(sm,'ILREF',REFCUR%JLREF(1))
        call json%add(sm,'ILSPT',REFCUR%JLSPT)
        call json%add(sm,'ISRS',REFCUR%JSRS(1))
        call json%add(sm,'ISRC',REFCUR%JSRC(1))
        call json%add(sm,'LCHSPNWL',REFCUR%JLCHSPNWL(1))
        call json%add(sm,'REFCELL',REFCUR%REFCELL)
        call json%add(sm,'ZNML',REFCUR%ZNMLR)
        call json%add(sm,'EWALL',REFCUR%EWALLR)
        call json%add(sm,'EWBIN',REFCUR%EWBINR)
        call json%add(sm,'TRANSP',REFCUR%TRANSPR(1,1:2))
        call json%add(sm,'FSHEAT',REFCUR%FSHEATR)
        call json%add(sm,'RECYCF',REFCUR%RCYCFR(1))
        call json%add(sm,'RECYCT',REFCUR%RCYCTR(1))
        call json%add(sm,'RECPRM',REFCUR%RCPRMR(1))
        call json%add(sm,'EXPPL',REFCUR%EXPPLR(1))
        call json%add(sm,'EXPEL',REFCUR%EXPELR(1))
        call json%add(sm,'EXPIL',REFCUR%EXPILR(1))
        call json%add(sm,'RECYCS',REFCUR%RCYCSR(1))
        call json%add(sm,'RECYCC',REFCUR%RCYCCR(1))
        call json%add(sm,'SPTPRM',REFCUR%STPRMR(1))
        call json%add(sm,'ESPUTS',REFCUR%ESPTSR(1))
        call json%add(sm,'ESPUTC',REFCUR%ESPTCR(1))
        call json%add(sm,'REFPRESS',REFCUR%REFPRESS)

        call json%create_array(vrs,'CHANGES')
        spr => REFCUR%SPEC_LINES
        do while(associated(spr))
          call json%create_object(vr,' ')
          call json%add(vr,'VARIABLE',spr%varname)
          call json%add(vr,'SPECIES',spr%spcname)
          if ((spr%varname == 'ISRS') .or.
     .        (spr%varname == 'ISRC')) then
            call json%add(vr,'VALUE',spr%ival)
          else
            call json%add(vr,'VALUE',spr%rval)
          end if
          call json%add(vrs,vr)
          spr => spr%next
        end do

        refcur => refcur%next

        call json%add(sm,vrs)
        call json%add(smod,sm)
      end do

      return
      end subroutine eirene_write_block_6


!******************************************************************************

      subroutine eirene_write_block_7(me)

      type(json_value),pointer :: me
      type(json_value),pointer :: srcs, src, subs, sub
      real(dp) :: ampts
      integer :: istra, j, nst, npts_orig, nminpts_orig

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.7")

      nst = nstrai_in

      ampts = mpts_comsou

      call json%add(me,'NSTRAI',nst)
      call json%add(me,'INDSRC',indsrc(1:nst))
      call json%add(me,'ALLOC',alloc)
      call json%add(me,'AMPTS',mpts_comsou)

      call json%create_array(srcs,'STRATA')
      call json%add(me, srcs)

      do istra = 1, nst
        if (indsrc(istra) == 6) cycle
        call json%create_object(src,'')

        call json%add(src,'ISTRA',istra)
        call json%add(src,'TXTSOU',txtsou(istra))

        call json%add(src,'NLAVRP',nlavrp(istra))
        call json%add(src,'NLAVRT',nlavrt(istra))
        call json%add(src,'NLSYMP',nlsymp(istra))
        call json%add(src,'NLSYMT',nlsymt(istra))
C       call json%add(src,'NLRAY',nlray(istra))

!  reverse enhancement of number of histories
        npts_orig = npts(istra)
        nminpts_orig = nminpts(istra)
        if (abs(ampts-1._dp) > eps10) then
          npts_orig = nint( real(npts_orig,dp) / ampts )
          nminpts_orig = nint( real(nminpts_orig,dp) / ampts )
        end if

!  999999 is the largest number of particles in input
!  thus npts was negative originally
        if (npts_orig > 999999) npts_orig = -1

        call json%add(src,'NPTS',npts_orig)
        call json%add(src,'NINITL',ninitl(istra))
        call json%add(src,'NEMODS',nemods(istra))
        call json%add(src,'NAMODS',namods(istra))
        call json%add(src,'NMINPTS',nminpts_orig)
        call json%add(src,'NPTSDEL',nptsdel(istra))
C       call json%add(src,'NRAYEN',nrayen(istra))

        call json%add(src,'FLUX',flux(istra))
        call json%add(src,'SCALV',scalv(istra))
        call json%add(src,'IVLSF',ivlsf(istra))
        call json%add(src,'ISCLS',iscls(istra))
        call json%add(src,'ISCLT',isclt(istra))
        call json%add(src,'ISCL1',iscl1(istra))
        call json%add(src,'ISCL2',iscl2(istra))
        call json%add(src,'ISCL3',iscl3(istra))
        call json%add(src,'ISCLB',isclb(istra))
        call json%add(src,'ISCLA',iscla(istra))

        call json%add(src,'NLATM',nlatm(istra))
        call json%add(src,'NLMOL',nlmol(istra))
        call json%add(src,'NLION',nlion(istra))
        call json%add(src,'NLPLS',nlpls(istra))
        call json%add(src,'NLPHOT',nlphot(istra))

        call json%add(src,'NSPEZ',nspez(istra))

        call json%add(src,'NLPNT',nlpnt(istra))
        call json%add(src,'NLLNE',nllne(istra))
        call json%add(src,'NLSRF',nlsrf(istra))
        call json%add(src,'NLVOL',nlvol(istra))
        call json%add(src,'NLCNS',nlcns(istra))

        call json%add(src,'NSRFSI',nsrfsi(istra))

        call json%create_array(subs,'SUBSTRATA')
        do j=1, nsrfsi(istra)
          call json%create_object(sub,'')

          call json%add(sub,'INUM',j)
          call json%add(sub,'INDIM',indim(j,istra))
          call json%add(sub,'INSOR',insor(j,istra))
          call json%add(sub,'INGRDA',ingrda(j,istra,1:3))
          call json%add(sub,'INGRDE',ingrde(j,istra,1:3))

          call json%add(sub,'SORWGT',sorwgt(j,istra))
          call json%add(sub,'SORLIM',sorlim(j,istra))
          call json%add(sub,'SORIND',sorind(j,istra))
          call json%add(sub,'SOREXP',sorexp(j,istra))
          call json%add(sub,'SORIFL',sorifl(j,istra))

          call json%add(sub,'NRSOR',nrsor(j,istra))
          call json%add(sub,'NPSOR',npsor(j,istra))
          call json%add(sub,'NTSOR',ntsor(j,istra))
          call json%add(sub,'NBSOR',nbsor(j,istra))
          call json%add(sub,'NASOR',nasor(j,istra))
          call json%add(sub,'NISOR',nisor(j,istra))
          call json%add(sub,'ISTOR',istor(j,istra))

          call json%add(sub,'SORAD1',sorad1(j,istra))
          call json%add(sub,'SORAD2',sorad2(j,istra))
          call json%add(sub,'SORAD3',sorad3(j,istra))
          call json%add(sub,'SORAD4',sorad4(j,istra))
          call json%add(sub,'SORAD5',sorad5(j,istra))
          call json%add(sub,'SORAD6',sorad6(j,istra))

          call json%add(subs,sub)
        end do
        call json%add(src, subs)

        call json%add(src,'SORENI',soreni(istra))
        call json%add(src,'SORENE',sorene(istra))
        call json%add(src,'SORVDX',sorvdx(istra))
        call json%add(src,'SORVDY',sorvdy(istra))
        call json%add(src,'SORVDZ',sorvdz(istra))

        call json%add(src,'SORCOS',sorcos_in(istra))
        call json%add(src,'SORMAX',sormax_in(istra))
        call json%add(src,'SORCTX',sorctx(istra))
        call json%add(src,'SORCTY',sorcty(istra))
        call json%add(src,'SORCTZ',sorctz(istra))
        call json%add(src,'RAYFRAC',rayfrac(istra))

        call json%add(srcs,src)
      end do

      return
      end subroutine eirene_write_block_7

!******************************************************************************

      subroutine eirene_write_block_8(me)
      type(json_value),pointer :: me, zones, zone, ch3lines
      type(json_value),pointer :: temp, dens, vels, vols
      type(s_stack), pointer :: cur
      TYPE(TEMPERATURE),POINTER :: TEMPCUR
      TYPE(DENSITY),POINTER :: DENCUR
      TYPE(VELOCITY),POINTER :: VELCUR
      TYPE(VOLUMEP),POINTER :: VOLCUR
      integer :: i

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.8")

      call json%add(me,'NZADD',nzadd)

      call json%create_array(zones,'ZONE')
      call json%add(me,zones)

      do i = 1, nzadd
        call json%create_object(zone,'')

        call json%add(zone,'INI',ini_zone(i))
        call json%add(zone,'INE',ine_zone(i))
! CH3 lines
        if (associated(ch3_stack(i)%head)) then
          call json%create_array(ch3lines,'CH3-LINES') !an array
          cur => ch3_stack(i)%head
          do while (associated(cur))
            call json%add(ch3lines,'CH3',cur%string)
            cur => cur%next
          end do
          call json%add(zone,ch3lines)
        end if

        tempcur => templist
        do while (associated(tempcur))
          if (tempcur%ii == ini_zone(i)) then
            call json%create_object(temp,'T')
            call json%add(temp,'IDION',tempcur%idion)
            call json%add(temp,'TE',tempcur%te)
            call json%add(temp,'TI',tempcur%ti)
            call json%add(zone,temp)
            exit
          end if
          tempcur => tempcur%next
        end do

        dencur => denlist
        do while (associated(dencur))
          if (dencur%ii == ini_zone(i)) then
            call json%create_object(dens,'D')
            call json%add(dens,'IDION',dencur%idion)
            call json%add(dens,'DI',dencur%di)
            call json%add(zone,dens)
            exit
          end if
          dencur => dencur%next
        end do

        velcur => vellist
        do while (associated(velcur))
          if (velcur%ii == ini_zone(i)) then
            if (velcur%iz == 1) then
              call json%create_object(vels,'M')
            else
              call json%create_object(vels,'V')
            endif
            call json%add(vels,'IDION',velcur%idion)
            call json%add(vels,'VX',velcur%vx)
            call json%add(vels,'VY',velcur%vy)
            call json%add(vels,'VZ',velcur%vz)
            call json%add(zone,vels)
            exit
          end if
          velcur => velcur%next
        end do

        volcur => vollist
        do while (associated(volcur))
          if (volcur%ii == ini_zone(i)) then
            call json%create_object(vols,'VL')
            call json%add(vols,'VOL',volcur%vol)
            call json%add(zone,vols)
            exit
          end if
          volcur => volcur%next
        end do

        call json%add(zones,zone)
      end do

      return
      end subroutine eirene_write_block_8

!******************************************************************************

      subroutine eirene_write_block_9(me)
      type(json_value),pointer :: me, svi, sv, ssi, ss, sci, sc
      integer :: iprsf(nlimps+1), i, num, ia

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.9")

      call json%add(me,'NLPRCA',nlprca(1:natmi))
      call json%add(me,'NLPRCM',nlprcm(1:nmoli))
      call json%add(me,'NLPRCI',nlprci(1:nioni))
      call json%add(me,'NLPRCPH',nlprcph(1:nphoti))

      iprsf = pack( [(i, i=0,nlimps)], nlprcs, [(-1,i=0,nlimps)] )
      num = count(iprsf > -1)
      call json%add(me,'NPRCSF',num)
      if (num > 0) call json%add(me,'IPRSF',iprsf(1:num))

      call json%add(me,'MAXLEV',maxlev)
      call json%add(me,'MAXRAD',maxrad)
      call json%add(me,'MAXPOL',maxpol)
      call json%add(me,'MAXTOR',maxtor)
      call json%add(me,'MAXADD',maxadd)

      if (maxrad > 0) then
        call json%add(me,'NSSPL_RAD',NSSPL(1:maxrad))
        call json%add(me,'PRMSPL_RAD',PRMSPL(1:maxrad))
      end if

      if (maxpol > 0) then
        call json%add(me,'NSSPL_POL',NSSPL(n1st+1:n1st+maxpol))
        call json%add(me,'PRMSPL_POL',PRMSPL(n1st+1:n1st+maxpol))
      end if

      if (maxtor > 0) then
        ia = n1st+n2nd
        call json%add(me,'NSSPL_TOR',NSSPL(ia+1:ia+maxtor))
        call json%add(me,'PRMSPL_TOR',PRMSPL(ia+1:ia+maxpol))
      end if

      if (maxadd > 0) then
        ia = n1st+n2nd+n3rd
        call json%add(me,'NSSPL_ADD',NSSPL(ia+1:ia+maxadd))
        call json%add(me,'PRMSPL_ADD',PRMSPL(ia+1:ia+maxadd))
      end if

      call json%add(me,'WMINV',wminv)
      call json%add(me,'WMINS',wmins)
      call json%add(me,'WMINC',wminc)
      call json%add(me,'WMINL',wminl)

      call json%add(me,'SPLPAR',splpar)

      call json%add(me,'NSIGVI',nsigvi)
      call json%add(me,'NSIGSI',nsigsi)
      call json%add(me,'NSIGCI',nsigci)
      call json%add(me,'NSIGI_SPC',nsigi_spc)

      if (nsigvi > 0) then
        call json%create_array(svi,'DEV_VOL_TAL')
        do i = 1, nsigvi
          call json%create_object(sv,'')
          call json%add(sv,'IGH',igh(i))
          call json%add(sv,'IIH',iih(i))
          call json%add(svi,sv)
        end do
        call json%add(me,svi)
      end if

      if (nsigsi > 0) then
        call json%create_array(ssi,'DEV_SRF_TAL')
        do i = 1, nsigsi
          call json%create_object(ss,'')
          call json%add(ss,'IGHW',ighw(i))
          call json%add(ss,'IIHW',iihw(i))
          call json%add(ssi,ss)
        end do
        call json%add(me,ssi)
      end if

      if (nsigci > 0) then
        call json%create_array(sci,'CORR_COEFFS')
        do i = 1, nsigci
          call json%create_object(sc,'')
          call json%add(sc,'IGHC',ighc(1:2,i))
          call json%add(sc,'IIHC',iihc(1:2,i))
          call json%add(sci,sc)
        end do
        call json%add(me,sci)
      end if

      return
      end subroutine eirene_write_block_9

!******************************************************************************

      subroutine eirene_write_block_10(me)
      type(json_value),pointer :: me
      type(json_value),pointer :: blk_10a, blk_10b, blk_10c,
     .                            blk_10d, blk_10e, blk_10f
      type(json_value),pointer :: avt, col, alg, sav, als, spc
      TYPE(EIRENE_SPECTRUM), POINTER :: ESPEC
      integer :: j, isf

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.10")

      call json%add(me,'NADVI',nadvi)
      call json%add(me,'NCLVI',nclvi)
      call json%add(me,'NALVI',nalvi)
      call json%add(me,'NADSI',nadsi)
      call json%add(me,'NALSI',nalsi)
      call json%add(me,'NADSPC',nadspc)

      call json%create_array(blk_10a,'ADD_VOL_TAL')
      call json%create_array(blk_10b,'COLL_EST')
      call json%create_array(blk_10c,'ALG_VOL_TAL')
      call json%create_array(blk_10d,'ADD_SRF_TAL')
      call json%create_array(blk_10e,'ALG_SRF_TAL')
      call json%create_array(blk_10f,'SPECTRA')

      call json%add(me,blk_10a)
      call json%add(me,blk_10b)
      call json%add(me,blk_10c)
      call json%add(me,blk_10d)
      call json%add(me,blk_10e)
      call json%add(me,blk_10f)

! additional volume averaged tallies
      do j = 1, NADVI
        call json%create_object(avt,'')
        call json%add(avt,'IADVE',iadve(j))
        call json%add(avt,'IADVS',iadvs(j))
        call json%add(avt,'IADVT',iadvt(j))
        call json%add(avt,'IADRC',iadrc(j))
        call json%add(avt,'TXTTLA',txttal(j,ntala))
        call json%add(avt,'TXTSCA',txtspc(j,ntala))
        call json%add(avt,'TXTUTA',txtunt(j,ntala))
        call json%add(blk_10a,avt)
      end do

! collision estimator
      do j = 1, NCLVI
        call json%create_object(col,'')
        call json%add(col,'ICLVE',iclve(j))
        call json%add(col,'ICLVS',iclvs(j))
        call json%add(col,'ICLVT',iclvt(j))
        call json%add(col,'ICLRC',iclrc(j))
        call json%add(col,'TXTTLC',txttal(j,ntalc))
        call json%add(col,'TXTSCC',txtspc(j,ntalc))
        call json%add(col,'TXTUTC',txtunt(j,ntalc))
        call json%add(blk_10b,col)
      end do

! algebraic volume tallies
      do j = 1, NALVI
        call json%create_object(alg,'')
        call json%add(alg,'CHRTAL',chrtal(j))
        call json%add(alg,'TXTTLR',txttal(j,ntalr))
        call json%add(alg,'TXTSCR',txtspc(j,ntalr))
        call json%add(alg,'TXTUTR',txtunt(j,ntalr))
        call json%add(blk_10c,alg)
      end do

! additional surface tallies
      do j = 1, NADSI
        call json%create_object(sav,'')
        call json%add(sav,'IADSE',iadse(j))
        call json%add(sav,'IADSS',iadss(j))
        call json%add(sav,'IADST',iadst(j))
        call json%add(sav,'TXTTLW',txttlw(j,ntlsa))
        call json%add(sav,'TXTSPW',txtspw(j,ntlsa))
        call json%add(sav,'TXTUNW',txtunw(j,ntlsa))
        call json%add(blk_10d,sav)
      end do

! algebraic surface tallies
      do j = 1, NALSI
        call json%create_object(als,'')
        call json%add(alg,'CHRTLS',chrtls(j))
        call json%add(alg,'TXTTLW',txttal(j,ntlsr))
        call json%add(alg,'TXTSPW',txtspc(j,ntlsr))
        call json%add(alg,'TXTUNW',txtunt(j,ntlsr))
        call json%add(blk_10e,alg)
      end do

! spectra
      do j = 1, nadspc
        call json%create_object(spc,'')
        espec => estiml(j)
        isf = espec%ispcsrf
        if (isf > nlim) isf = -(isf-nlim)
        call json%add(spc,'ISPSRF',isf)
        call json%add(spc,'IPTYP',espec%iprtyp)
        call json%add(spc,'IPSPZ',espec%iprsp)
        call json%add(spc,'ISPTYP',espec%ispctyp)
        call json%add(spc,'NSPS',espec%nspc)
        call json%add(spc,'ISRFCLL',espec%isrfcll)
        call json%add(spc,'IDIREC',espec%idirec)
        call json%add(spc,'ISPCOPT',espec%ispcopt)
        if (espec%ispcopt == 2) then
          call json%add(spc,'ISPLDEG',espec%ispldeg)
        end if

        call json%add(spc,'SPCMN',espec%spcmin)
        call json%add(spc,'SPCMX',espec%spcmax)
        call json%add(spc,'SPC_SHIFT',espec%esp_00)
        call json%add(spc,'SPCPLT_X',espec%spc_xplt)
        call json%add(spc,'SPCPLT_Y',espec%spc_yplt)
        call json%add(spc,'SPCPLT_SAME',espec%spc_same)

        if (espec%idirec /= 0) then
          call json%add(spc,'SPCVX',espec%spcvx)
          call json%add(spc,'SPCVY',espec%spcvy)
          call json%add(spc,'SPCVZ',espec%spcvz)
        end if

        call json%add(blk_10f,spc)
      end do

      return
      end subroutine eirene_write_block_10

!******************************************************************************

      subroutine eirene_write_block_11(me)
      type(json_value),pointer :: me, tals, tal, srfs, srf
      type(json_value),pointer :: plads, plad, plstds, plstd,
     .                            pltls, pltl, spcs, spc
      integer :: j, nsrf, i

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.11")

! block 11a, print output
      call json%add(me,'TRCPLT',trcplt)
      call json%add(me,'TRCHST',trchst)
      call json%add(me,'TRCNAL',trcnal)
      call json%add(me,'TRCMOD',trcmod)
      call json%add(me,'TRCSIG',trcsig)

      call json%add(me,'TRCGRD',trcgrd)
      call json%add(me,'TRCSUR',trcsur)
      call json%add(me,'TRCREF',trcref)
      call json%add(me,'TRCFLE',trcfle)
      call json%add(me,'TRCAMD',trcamd)

      call json%add(me,'TRCINT',trcint)
      call json%add(me,'TRCLST',trclst)
      call json%add(me,'TRCSOU',trcsou)
      call json%add(me,'TRCREC',trcrec)
      call json%add(me,'TRCTIM',trctim)

      call json%add(me,'TRCBLA',trcbla)
      call json%add(me,'TRCBLM',trcblm)
      call json%add(me,'TRCBLI',trcbli)
      call json%add(me,'TRCBLP',trcblp)
      call json%add(me,'TRCBLE',trcble)

      call json%add(me,'TRCBLPH',trcblph)
      call json%add(me,'TRCTAL',trctal)
      call json%add(me,'TRCOCT',trcoct)
      call json%add(me,'TRCCEN',trccen)
      call json%add(me,'TRCRNF',trcrnf)

!      call json%add(me,'TRCDBG2',trcdbg2)
!      call json%add(me,'TRCDBGE',trcdbge)
!      call json%add(me,'TRCDBGM',trcdbgm)
!      call json%add(me,'TRCDBGF',trcdbgf)
!      call json%add(me,'TRCDBGL',trcdbgf)

!      call json%add(me,'TRCDBGS',trcdbgs)
!      call json%add(me,'TRCDBGG',trcdbgg)
!      call json%add(me,'TRCDBGMPI',trcdbgmpi)
!      call json%add(me,'TRCDBGC',trcdbgc)
      call json%add(me,'TRCHKTIM',trchktim)
      call json%add(me,'TRCSCL',trcscl)

      call json%add(me,'TRCSRC',trcsrc)

! volume averaged tallies
      call json%add(me,'NVOLPR',nvolpr)
      call json%add(me,'NSPCPR',nspcpr)

      call json%create_array(tals,'OUT_VOL_TAL')
      do j = 1, nvolpr
         call json%create_object(tal,'')
         call json%add(tal,'NTLV',nprtlv(j))
         call json%add(tal,'NFLGV',nflagv(j))
         call json%add(tal,'NSPEZV',nspezv(j,1:2))
         call json%add(tal,'NTLVF',ntlvfl(j))
         call json%add(tals,tal)
      end do
      call json%add(me,tals)

! output of surfaces
      call json%add(me,'NSURPR',nsurpr)

      call json%create_array(srfs,'SRF_OUTPUT')
      do j = 1, nsurpr
         call json%create_object(srf,'')
         nsrf = nprsrf(j)
         if (nsrf > nlim) nsrf = -(nsrf-nlim)
         if (nsrf == nstsi+1) nsrf = 0
         call json%add(srf,'NSRF',nsrf)
         call json%add(srf,'NTLS',nprtls(j))
         call json%add(srf,'NFLGS',nflags(j))
         call json%add(srf,'NSPEZS',nspezs(j,1:2))
         call json%add(srf,'NTLSF',ntlsfl(j))
         call json%add(srfs,srf)
      end do
      call json%add(me,srfs)

! switch off tallies
      if (allocated(itlvout)) then
        call json%add(me,'NTLVOUT',ntlvout)
        call json%add(me,'NUMTAL_V',itlvout(1:ntlvout))
        deallocate(itlvout)
      end if

      if (allocated(itlsout)) then
        call json%add(me,'NTLSOUT',ntlsout)
        call json%add(me,'NUMTAL_S',itlsout(1:ntlsout))
        deallocate(itlsout)
      end if

! block 11 b, geometrie plot

      call json%add(me,'PL1ST',pl1st)
      call json%add(me,'PL2ND',pl2nd)
      call json%add(me,'PL3RD',pl3rd)
      call json%add(me,'PLADD',pladd)
      call json%add(me,'PLHST',plhst)

      call json%add(me,'PLCUT',plcut)
      call json%add(me,'PLBOX',plbox)
      call json%add(me,'PLSTOR',plstor)

      call json%add(me,'PLNUMV',plnumv)
      call json%add(me,'PLNUMS',plnums)
      call json%add(me,'PLARR',plarr)
      call json%add(me,'LRPSCUT',lrpscut)
      call json%add(me,'PLIDL',plidl)
      call json%add(me,'PLVTK',plvtk)

      call json%add(me,'NPLINR',nplinr)
      call json%add(me,'NPLOTR',nplotr)
      call json%add(me,'NPLDLR',npldlr)
      call json%add(me,'NPLINP',nplinp)
      call json%add(me,'NPLOTP',nplotp)
      call json%add(me,'NPLDLP',npldlp)
      call json%add(me,'NPLINT',nplint)
      call json%add(me,'NPLOTT',nplott)
      call json%add(me,'NPLDLT',npldlt)

      call eirene_write_block_11_usr(json,me)

      call json%create_array(plads,'3D_ADD_SRF')
      do j = 1, 5
        call json%create_object(plad,'')
        call json%add(plad,'PL3A',pl3a(j))
        call json%add(plad,'TEXTLA',textla(j))
        call json%add(plad,'IPLTA',iplta(j))
        call json%add(plad,'IPLAA',iplaa(j,1:iplta(j)))
        call json%add(plad,'IPLEA',iplea(j,1:iplta(j)))
        call json%add(plads,plad)
      end do
      call json%add(me,plads)

      call json%create_array(plstds,'3D_STD_SRF')
      do j = 1, 3
        call json%create_object(plstd,'')
        call json%add(plstd,'PL3S',pl3s(j))
        call json%add(plstd,'TEXTLS',textls(j))
        call json%add(plstd,'IPLTS',iplts(j))
        call json%add(plstd,'IPLAS',iplas(j,1:iplta(j)))
        call json%add(plstd,'IPLES',iples(j,1:iplta(j)))
        call json%add(plstds,plstd)
      end do
      call json%add(me,plstds)

      call json%add(me,'CH2MX',ch2mx)
      call json%add(me,'CH2MY',ch2my)
      call json%add(me,'CH2X0',ch2x0)
      call json%add(me,'CH2Y0',ch2y0)
      call json%add(me,'CH2Z0',ch2z0)

      call json%add(me,'CH3MX',ch3mx)
      call json%add(me,'CH3MY',ch3my)
      call json%add(me,'CH3MZ',ch3mz)
      call json%add(me,'CH3X0',ch3x0)
      call json%add(me,'CH3Y0',ch3y0)
      call json%add(me,'CH3Z0',ch3z0)

      call json%add(me,'ANGLE1',angle1)
      call json%add(me,'ANGLE2',angle2)
      call json%add(me,'ANGLE3',angle3)

      call json%add(me,'I1TRC',i1trc)
      call json%add(me,'I2TRC',i2trc)
      call json%add(me,'ISYPLT',isyplt)
      call json%add(me,'ILINIE',ilinie)

! block 11 c, plots of tallies

      call json%add(me,'NVOLPL',nvolpl)

      call json%add(me,'PLTSRC',pltsrc)

      if (lrpscut) call json%add(me,'CUTPLANE',cutplane)

      if (nvolpl > 0) then
        call json%create_array(pltls,'PLOT_TALLIES')
        do j = 1, nvolpl
          call json%create_object(pltl,'')

          call json%add(pltl,'NSP',nsptal(j))

          call json%add(pltl,'PLTL2D',pltl2d(j))
          call json%add(pltl,'PLTL3D',pltl3d(j))
          call json%add(pltl,'PLTLLG',pltllg(j))
          call json%add(pltl,'PLTLER',pltler(j))

          call json%add(pltl,'TALZMI',talzmi(j))
          call json%add(pltl,'TALZMA',talzma(j))
          call json%add(pltl,'TALXMI',talxmi(j))
          call json%add(pltl,'TALXMA',talxma(j))
          call json%add(pltl,'TALYMI',talymi(j))
          call json%add(pltl,'TALYMA',talyma(j))

          if (pltl2d(j)) then
            call json%add(pltl,'LHIST2',lhist2(j))
            call json%add(pltl,'LSMOT2',lsmot2(j))
            call json%create_array(spcs,'SPECIES')
            do i = 1, nsptal(j)
              call json%create_object(spc,'')
              call json%add(spc,'ISPTAL',isptal(j,i))
              call json%add(spc,'NTL',nptali(j,i))
              call json%add(spc,'NPLIN2',nplin2(j,i))
              call json%add(spc,'NPLOT2',nplot2(j,i))
              call json%add(spc,'NPLDL2',npldl2(j,i))
              call json%add(spcs,spc)
            end do
            call json%add(pltl,spcs)
          end if

          if (pltl3d(j)) then
            call json%add(pltl,'LHIST3',lhist3(j))
            call json%add(pltl,'LCNTR3',lcntr3(j))
            call json%add(pltl,'LSMOT3',lsmot3(j))
            call json%add(pltl,'LRAPS3',lraps3(j))
            call json%add(pltl,'LVECT3',lvect3(j))
            call json%add(pltl,'LRPVC3',lrpvc3(j))
            call json%add(pltl,'LRPS3D',lraps3d)
            call json%add(pltl,'LRPSCN',lr3dcon)

            call json%add(pltl,'LPRAD3',lprad3(j))
            call json%add(pltl,'LPPOL3',lppol3(j))
            call json%add(pltl,'LPTOR3',lptor3(j))

            call json%create_array(spcs,'SPECIES')
            do i = 1, nsptal(j)
              call json%create_object(spc,'')
              call json%add(spc,'ISPTAL',isptal(j,i))
              call json%add(spc,'NTL',nptali(j,i))
              call json%add(spc,'IPROJ3',iproj3(j,i))
              call json%add(spc,'NPLI13',npli13(j,i))
              call json%add(spc,'NPLO13',nplo13(j,i))
              call json%add(spc,'NPLI23',npli23(j,i))
              call json%add(spc,'NPLO23',nplo23(j,i))
              call json%add(spc,'IPLN',iplane)
              call json%add(spcs,spc)
            end do
            call json%add(pltl,spcs)
          end if

          call json%add(pltl,'TALW1',talw1(j))
          call json%add(pltl,'TALW2',talw2(j))
          call json%add(pltl,'FCABS1',fcabs1(j))
          call json%add(pltl,'FCABS2',fcabs2(j))
          call json%add(pltl,'RPSDL',rapsdel)

          call json%add(pltls,pltl)
        end do
        call json%add(me,pltls)
      end if

      return
      end subroutine eirene_write_block_11

!******************************************************************************

      subroutine eirene_write_block_12(me)
      type(json_value),pointer :: me, chrds, chrd, line, compo, contr,
     .                            rat, lines, comps, conts, rats
      type(tcontrib) :: cnt
      integer :: i, istchr, j, k

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.12")

!pb   if (num_lines > 0) then
      if (ldef_lines) then
        call json%add(me,'KEYWORD', "DEFINE_LINES")
        call json%add(me,'NUM_LINES',num_lines)
        call json%add(me,'MOD_ADDV',mod_addv)

        call json%create_array(lines,'LINES')

        do i = 1, num_lines

          call json%create_object(line,'')
          call json%add(line,'LINE_NAME',emis_lines(i)%line_name)
          call json%add(line,'NUM_COMPO',emis_lines(i)%num_compo)
          call json%add(line,'EINSTEIN',emis_lines(i)%einstein)
          call json%add(line,'TRANS_EN',emis_lines(i)%trans_en)

          call json%create_array(comps,'COMPONENTS')

          do j = 1, emis_lines(i)%num_compo

            call json%create_object(compo,'')
            call json%add(compo,'COMPO_NAME',
     .                    emis_lines(i)%compo(j)%compo_name)
            call json%add(compo,'NUM_CONTRIB',
     .                    emis_lines(i)%compo(j)%num_contrib)
            call json%add(compo,'COMPO_IRC',
     .                    emis_lines(i)%compo(j)%irc)

            call json%create_array(conts,'CONTRIBUTIONS')

            do k = 1, emis_lines(i)%compo(j)%num_contrib
              cnt = emis_lines(i)%compo(j)%contrib(k)
              call json%create_object(contr,'')
              call json%add(contr,'ISP',cnt%isp)
              call json%add(contr,'ITP',cnt%itp)
              call json%add(contr,'IRATIO',cnt%iratio)
              if (cnt%iratio > 0) then
                call json%create_array(rats,'RATIOS')
                call json%create_object(rat,'')
                call json%add(rat,'IRC_RAT',cnt%irc_rat(1))
                call json%add(rats,rat)
                if (cnt%iratio == 2) then
                  call json%create_object(rat,'')
                  call json%add(rat,'IRC_RAT',cnt%irc_rat(2))
                  call json%add(rat,'ISP_RAT',cnt%isp_rat(1:2))
                  call json%add(rat,'ITP_RAT',cnt%itp_rat(1:2))
                  call json%add(rats,rat)
                end if ! iratio == 2
                call json%add(contr,rats)
              end if ! iratio > 0
              call json%add(conts,contr)
            end do ! contributions
            call json%add(compo,conts)
            call json%add(comps,compo)
          end do ! components

          call json%add(line,comps)
          call json%add(lines,line)
        end do  ! lines
        call json%add(me,lines)

      else if (nlemis) then

        call json%add(me,'KEYWORD', "DEFAULT_LINES")

      end if

!  lines of sight

      call json%add(me,'NCHORI',nchori)
      call json%add(me,'NCHENI',ncheni)
      if (nchori > 0) then
        call json%create_array(chrds,'CHORDS')

        do i = 1, nchori
          call json%create_object(chrd,'')
          call json%add(chrd,'TXTSIG',txtsig(i))
          call json%add(chrd,'NCHTAL',nchtal(i))
          call json%add(chrd,'NSPSCL',nspscl(i))
          call json%add(chrd,'NSPNEW',nspnew(i))
          if (nlstchr(i)) then
            istchr = 1
          else
            istchr = 0
          end if
          call json%add(chrd,'ISTCHR',istchr)

          if (len_trim(ch_line_name(i)) > 0)
     .        call json%add(chrd,'USE_LINE',ch_line_name(i))

          call json%add(chrd,'NSPSTR',nspstr(i))
          call json%add(chrd,'NSPSPZ',nspspz(i))
          call json%add(chrd,'NSPINI',nspini(i))
          call json%add(chrd,'NSPEND',nspend(i))
          call json%add(chrd,'NSPBLC',nspblc(i))
          call json%add(chrd,'NSPADD',nspadd(i))

          call json%add(chrd,'EMIN1',emin1(i))
          call json%add(chrd,'EMAX1',emax1(i))
          call json%add(chrd,'ESHIFT',eshift(i))

          call json%add(chrd,'IPIVOT',ipivot(i))
          call json%add(chrd,'XPIVOT',xpivot(i))
          call json%add(chrd,'YPIVOT',ypivot(i))
          call json%add(chrd,'ZPIVOT',zpivot(i))

          call json%add(chrd,'ICHORD',ichord(i))
          call json%add(chrd,'XCHORD',xchord(i))
          call json%add(chrd,'YCHORD',ychord(i))
          call json%add(chrd,'ZCHORD',zchord(i))

          call json%add(chrds,chrd)
        end do

        call json%add(me,chrds)
        call json%add(me,'PLCHOR',plchor)
        call json%add(me,'PLSPEC',plspec)
        call json%add(me,'PRSPEC',prspec)
        call json%add(me,'PLARGL',plargl)
        call json%add(me,'PRARGL',prargl)
      end if

      return
      end subroutine eirene_write_block_12

!******************************************************************************

      subroutine eirene_write_block_13(me)
      type(json_value),pointer :: me, snps, snp
      integer :: i

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.13")

      call json%add(me,'NPRNLI',nprnli_in)
      call json%add(me,'DEFAULT_TIME_HORIZON',ldef_time_horizon)

      if (nprnli_in > 0) then
        call json%add(me,'NINITL_READ',ninitl_read)
        call json%add(me,'NPRMUL',nprmul)

        if (.not.ldef_time_horizon) then
          call json%add(me,'NPTST',nptst)
          call json%add(me,'NTMSTP',ntmstp)

          call json%add(me,'DTIMV',dtimv)
          call json%add(me,'TIME0',time0)

          call json%add(me,'NSNVI',nsnvi)

          if (nsnvi > 0) then
            call json%create_array(snps,'SNAPSHOTS')

            do i = 1, nsnvi
              call json%create_object(snp,'')
              call json%add(snp,'ISNVE',isnve(i))
              call json%add(snp,'ISNVS',isnvs(i))
              call json%add(snp,'ISNVT',isnvt(i))
              call json%add(snp,'ISNRC',isnrc(i))
              call json%add(snp,'TXTTLT',TXTTAL(i,ntalt))
              call json%add(snp,'TXTSCT',TXTSPC(i,ntalt))
              call json%add(snp,'TXTUTT',TXTUNT(i,ntalt))
              call json%add(snps,snp)
            end do
            call json%add(me,snps)
          end if
        end if
      end if

      return
      end subroutine eirene_write_block_13

!******************************************************************************

      subroutine eirene_write_block_14(me)
      type(json_value),pointer :: me, intls, intl
      integer :: i
      interface
        subroutine  eirene_wr0_json (json,this)
        use json_module           !IGNORE
        class(json_core),intent(inout) :: json
        type(json_value),pointer, intent(inout) :: this
        end subroutine  eirene_wr0_json

!        subroutine eirene_wrjson_usr(me)
!        use json_module           !IGNORE
!        type(json_value),pointer :: me
!        end subroutine eirene_wrjson_usr
      end interface

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.14")

      if (nmode == 0) then
        call json%add(me,'NAINI',naini)
!pb     call json%add(me,'NCOPII',ncopii)
        call json%add(me,'NCOPIE',ncpvi)
        if (naini > 0) then
           call json%create_array(intls,'ADD_IN_TAL')
           do i = 1, naini
             call json%create_object(intl,'')
             call json%add(intl,'NAINS',nains(i))
             call json%add(intl,'NAINT',naint(i))
             call json%add(intl,'TXTPLS',txtpls(i,ntaln))
             call json%add(intl,'TXTPSP',txtpsp(i,ntaln))
             call json%add(intl,'TXTPUN',txtpun(i,ntaln))
             call json%add(intls,intl)
           end do
           call json%add(me,intls)
        end if
      else
        call eirene_wr0_json(json,me)
      end if

!      call eirene_wrjson_usr(me)

      return
      end subroutine eirene_write_block_14

!******************************************************************************

      subroutine eirene_write_block_15(me)
      type(json_value),pointer :: me
      character(80) :: strategy

      call json%add(me,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#section.2.15")

      select case(NPRLL)
      case (-1)
        strategy = 'USER_DEFINED'
      case (0)
        strategy = 'EMBARRASS'
      case (1)
        strategy = 'ORIGINAL'
      case (2)
        strategy = 'APCAS'
      case (3)
        strategy = 'BALANCED'
      case (4)
        strategy = 'IMPROVED'
      case default
        strategy = 'AUTOMATIC'
      end select

      call json%add(me,'STRATEGY',strategy)

      return
      end subroutine eirene_write_block_15
#endif
!******************************************************************************

      end subroutine eirene_write_json_file
