      subroutine eirene_write_json_amdata(filename)

      use eirmod_precision
      use eirmod_parmmod
      USE EIRMOD_COMUSR
      USE EIRMOD_JSON
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT
#ifndef NO_JSON
      use json_module           !IGNORE
#endif

      implicit none
      character(*) :: filename

#ifndef NO_JSON
      type(json_core) :: json
      type(json_value),pointer :: pAMD
      type(json_value),pointer :: crs_lines
      type(s_stack), pointer   :: cur

! initialize the class
      call json%initialize(compact_reals=.true.,real_format='ES',
     .          case_sensitive_keys=.false.,compress_vectors=.false.)

! initialize the structure:
      call json%create_object(pAMD,'')

      call eirene_write_block_45(pAMD)

! write the file:
      call json%print(pAMD,trim(filename))

!cleanup:
      call json%destroy(pAMD)
      if (json%failed()) then
         write (iunout, *)'ERROR by writing json AMData <',
     >                    filename,'>'
        stop 1
      endif

!     call eirene_dealloc_json_arrays
#endif
      return

#ifndef NO_JSON
      contains

!***************************************************************************************

      subroutine eirene_write_block_45(pAMD)

      type(json_value),pointer :: pAMD
      type(json_value),pointer :: PART, REACT

      integer :: i
#ifdef NEVERCOMEHERE
      integer :: ndum1(1), ndum2(1), ndum3(1), ndum4(1), nti, nv
#endif

      character(300) :: sCReactStr

      call json%add(pAMD,'MANUAL',
     .     "http://www.eirene.de/eirene.pdf#subsection.2.4.1")
      call json%add(pAMD,'INCLUDE',
     .     "eirene_physics_model.input.json")

! initialize the structure:
      call json%create_object(REACT,'Reactions')
      call json%create_object(PART,'Particles')

      call json%add(pAMD, REACT)     !add it to the root
      call json%add(pAMD, PART)      !add it to the root

! start block 4

      call json%add(REACT,'NREACI',nreaci)

! write classical reaction lines for easy comparison with fixed input
      if (associated(crs_stack%head)) then ! an array
        call json%create_array(crs_lines,'ClassicReactionStrings')
        cur => crs_stack%head
        do while (associated(cur))
          call json%add(crs_lines,'CRS',trim(cur%string))
          cur => cur%next
        end do
        call json%add(REACT,crs_lines)
      end if

      do i=1, IRLINES
 100     format(i3,1x,A8,1x,A4,1x,A,1x,A3,1x,i3,1x,i3,5es12.4)
         write (sCReactStr, 100) reaclines(i)%no, reaclines(i)%FILE,
     >         reaclines(i)%h_select,reaclines(i)%reac_string,
     >         reaclines(i)%reactyp,reaclines(i)%mp,
     >         reaclines(i)%mt,reaclines(i)%dpp,
     >         reaclines(i)%r1mn,reaclines(i)%r1mx,
     >         reaclines(i)%r2mn,reaclines(i)%r2mx

          if (reaclines(i)%r1mn > 0)
     >    write (sCReactStr, "(A/,i3,1x,3es12.4)")  sCReactStr,
     >           reaclines(i)%jfex1mn, reaclines(i)%fp1(1:3)
          if (reaclines(i)%r1mx.gt.0)
     >    write (sCReactStr, "(A,i3,1x,3es12.4)")  sCReactStr,
     >           reaclines(i)%jfex1mx, reaclines(i)%fp1(4:6)
          if (reaclines(i)%r2mn.gt.0)
     >    write (sCReactStr, "(A,i3,1x,3es12.4)")  sCReactStr,
     >           reaclines(i)%jfex2mn, reaclines(i)%fp2(1:3)
          if (reaclines(i)%r2mx.gt.0)
     >    write (sCReactStr, "(A,i3,1x,3es12.4)")  sCReactStr,
     >           reaclines(i)%jfex2mx, reaclines(i)%fp2(4:6)

          call json%add(REACT, 'CStr' ,sCReactStr)
      enddo


#ifdef NEVERCOMEHERE
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
        call json%add(rea,'R1MX',reaclines(i)%r1mn)
        call json%add(rea,'R2MN',reaclines(i)%r2mn)
        call json%add(rea,'R2MX',reaclines(i)%r2mn)
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
        end if

        if (index(reaclines(i)%FILE,'CR') /= 0) then
           call json%add(rea,'M_POPESC',reaclines(i)%m_popesc)
           if (reaclines(i)%m_popesc > 0) then
            call json%add(rea,'IROW_ESC',
     .       reaclines(i)%irow_esc(1:reaclines(i)%m_popesc))
            call json%add(rea,'ICOL_ESC',
     .       reaclines(i)%icol_esc(1:reaclines(i)%m_popesc))
            call json%add(rea,'POP_ESC',
     .       reaclines(i)%pop_esc(1:reaclines(i)%m_popesc))
           end if
          call json%add(rea,'M_QEXT',reaclines(i)%m_qext)
          if (reaclines(i)%m_qext > 0) then
            call json%add(rea,'IROW_EXT',
     .          reaclines(i)%irow_ext(1:reaclines(i)%m_qext))
            call json%add(rea,'ITAL_EXT_1',
     .          reaclines(i)%ital_ext(1,1:reaclines(i)%m_qext))
            call json%add(rea,'ITAL_EXT_2',
     .          reaclines(i)%ital_ext(2,1:reaclines(i)%m_qext))
            call json%add(rea,'ITAL_EXT_3',
     .          reaclines(i)%ital_ext(3,1:reaclines(i)%m_qext))
            end if
         end if

        call json%add(reacs,rea)

      end do

      call json%add(REACT,reacs)

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
     .     eelecm, ebulkm, escd1m, freacm, edpotm)
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
     .     eeleci, ebulki, escd1i, freaci, edpoti)
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
     .     eelecp, ebulkp, escd1p, freacp, edpotp, cdenmodel)
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

#endif
      return
      end subroutine eirene_write_block_45

!******************************************************************************
#ifdef NEVERCOMEHERE

      subroutine eirene_write_block_4abcd (me, ityp, nbas, nloop, ndim,
     .     texts, cndx, cext, nmass, nchar, nprt, nchrg, isrf, isrt,
     .     nrc, nfol, ngen, nhsts, ireac, ibulk,
     .     iscd1, iscd2, iscd3, iscd4, iscde, iestm, ibgk,
     .     eelec, ebulk, escd1, freac, edpot, cdenmodel)

      type(json_value),pointer :: me
      type(json_value),pointer :: elem, reas, re, model
      integer, intent(in) :: ityp, nbas, nloop, ndim
      integer, intent(in) :: nmass(*), nchar(*),
     .         isrf(nspz,*), isrt(nspz,*), nprt(*), nchrg(*),
     .         nrc(*), nfol(*), ngen(*), nhsts(*),
     .         ireac(ndim,*), ibulk(ndim,*), iscd1(ndim,*),
     .         iscd2(ndim,*), iscd3(ndim,*), iscd4(ndim,*),
     .         iscde(ndim,*), iestm(ndim,*), ibgk(ndim,*)
      real(dp), intent(in) :: eelec(ndim,*), ebulk(ndim,*),
     .         escd1(ndim,*), freac(ndim,*), edpot(ndim,*)
      character(*) :: texts(*), cndx, cext
      character(*), optional :: cdenmodel(ndim)
      integer :: i, ispz, numsec, k, j, nre
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

#endif
!******************************************************************************
#endif
      end subroutine eirene_write_json_amdata
