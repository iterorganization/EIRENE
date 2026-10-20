C+---------------------------------------------------------------+
C| Purpose:                                                      |
C| --------                                                      |
C| Write the eirene.tranfer file to give the EIRENE results back |
C| to EDGE2D.                                                    |
C| Write also eirene.chemFluxDep file, to store the neutral      |
C| particle fluxes on to the walls for the next EIRENE call      |
C+---------------------------------------------------------------+
C| Modifications:                                                |
C| --------------                                                |
C| 16/07/2010   D.Harting    Added writing of neutral fluxes to  |
C|                           file eirene.chemFluxDep. Added also |
C|                           two variables to eirene_user        |
C|                           namelist.                           |
C| 23/11/2010   D.Harting    If EDGE2D is used with density      |
C|                           control by puff+recycling, the      |
C|                           number of puffeing surfaces and thus|
C|                           the number of additional surfaces   |
C|                           (NLIM) may vary. Before, the neutral|
C|                           flux file eirene.chemFluxDep from   |
C|                           previous run was checked to have the|
C|                           same number of add. surfaces as the |
C|                           actual run. This forced a stopping  |
C|                           of the code. Now the actual triangle|
C|                           number and its side is checked, to  |
C|                           asure that the right neutral flux   |
C|                           is used.                            |
C| 24/11/2010   D.Harting    Added for backward compatibility    |
C|                           a version number to the neutral flux|
C|                           file eirene.chemFluxDep. If the     |
C|                           Version number in the file and in   |
C|                           code are not matching, the file is  |
C|                           not read and zero neutral flux is   |
C|                           assumed in the actual run. At the   |
C|                           end of the eirene run, a new neutral|
C|                           flux file with the current version  |
C|                           number is generated.                |
C| 14/09/2020   N. Horsten   Writing out if EIRENE is run in     |
C|                           time-dependent mode (NTIME), as     |
C|                           NSTRA and NSTS are only increased by|
C|                           1 if EIRENE is run in time-dependent|
C|                           mode.                               |
C+---------------------------------------------------------------+
      subroutine EIRENE_outusr
      use eirmod_precision
      use eirmod_parmmod
      use eirmod_comusr
      use eirmod_ccona
      use eirmod_cestim
      use eirmod_ctrcei
      use eirmod_ctext
      use eirmod_comxs
      use eirmod_cgeom
      use eirmod_cgrid
      use eirmod_cinit
      use eirmod_ctrig
      use eirmod_csdvi
      use eirmod_clogau
      use eirmod_comsou
      use eirmod_comprt
      use eirmod_czt1
      use eirmod_coutau
      use eirmod_learc1, only : EIRENE_LEARC1
      use eirmod_usrdata

      implicit none
      integer :: fp,mtri,ir,is

      real(dp) :: dummy(1)
      integer :: ier,idum
      integer :: eirene_learc2

      real(dp) :: pa,pi,pm,pph,vvol,val,sumvol,sumval,vdenpara,bb
      integer :: j,i
      real(dp) :: c1,c2,c3,c4,precom,erecom

csw ratecoeff.dat
      integer :: iaei, iacx, imei, imcx, iiei, iicx, iirc
      integer :: irei, ircx, irrc
      real(dp) :: de,di,ti
      integer :: idsc,ipl,kk,nrc,ireac,k,np,nr,msg
csw 25oct07
      real(dp) :: x1,x2,y1,y2,ar,xc
      real(dp),allocatable :: sumpotpl(:)
csw
      logical :: lcxsigma, lexists

      integer, external :: EIRENE_idez

      logical, save :: ldebug

      logical  :: lfound
      real(dp) :: neutralFluxFileVersion, neutralTransferFileVersion
      integer  :: ll

Ccks The variable to hold the name of the eirene.surface_transfer file (="eirene.surface_transfer")
      character(len=64), save :: eirene_surface_transfer_filename

Ccks loop index and eirene.elemente number of triangles
      integer :: IATMsurf, ELEM_TRI_TOTAL
Ccks  dummy variables for reading eirene.elemente
      integer :: dummy1, dummy2, dummy3, dummy4, dummy5, dummy6
      integer :: dummy7, dummy8, dummy9
Ccks  counter for reporting the total number of eirene.sutface_transfer cells
      integer :: counter
Ccks  arrays holding the surface, ii, jj indices from eirene.elemente
      integer, allocatable :: SURF_INDX(:), II_INDX(:), JJ_INDX(:)

c     replicate old sputtered flux arrays sptpl, sptat, sptml, sptio, sptpht
c     for the moment these are filled with values from sptpltot, sptatot,sptmtot,sptitot,sptphtot
c     in the future it is better to pass particle-resolved sputtered fluxes
c     from sptXY X=PH,I,A,M,P Y=PHT,IO,AT,ML,PL
      real(dp),dimension(npls_fix,nlmpgs)  :: sptpl
      real(dp),dimension(natm,nlmpgs)  :: sptat
      real(dp),dimension(nmol,nlmpgs)  :: sptml
      real(dp),dimension(nion,nlmpgs)  :: sptio
      real(dp),dimension(nphot,nlmpgs) :: sptpht

      NeutralFluxFileVersion = 1.2
      NeutralTransferFileVersion = 2.0
      ldebug=.false.
csw
      ll=len_trim(casename)
      eirene_ftransfer = 'eirene.transfer'
      eirene_njetto=0
cdmh
      eirene_fstoreneutflux = 'eirene.chemFluxDep'
      eirene_wallFluxModel = 1
cdmh
      open(unit=9998,file='eirene_user.namelist')
      read(9998,eirene_user)
      close(9998)
csw

!cks  ------------------------------------------------------------------------
!cks  --------------------eirene.surface_transfer creation--------------------

!cks  Creation of file eirene.surface_transfer to output the reflected power from the
!cks  target at higher resolution than just a total (cell-level resolution)

      eirene_surface_transfer_filename = casename(1:ll)
     &     // '.surface_transfer'

!cks  Obtain the II,JJ EDGE2D to surface index correspondence from eirene.elemente
      eirene_felemente = casename(1:ll)
     &     // '.elemente'
      ! check if this is and EDGE2D case
      IF (casename(1:ll) .eq. 'eirene') THEN
         inquire(file=trim(eirene_felemente),exist=lexists)
         IF (.NOT.LEXISTS) THEN
            write(*,*) "OUTUSR: File: '"//trim(eirene_felemente)
     &           //"' does not exist"
            write(*,*) "STOPPING...."
            call EIRENE_exit_own(1)
         ENDIF
        open(9988,file=trim(eirene_felemente))
        read(9988, '(I30)')  ELEM_TRI_TOTAL
        allocate(SURF_INDX(ELEM_TRI_TOTAL))
        allocate(II_INDX(ELEM_TRI_TOTAL))
        allocate(JJ_INDX(ELEM_TRI_TOTAL))
        do k=1,ELEM_TRI_TOTAL
!         read(9988,'(12I10)') SURF_INDX(k),
          read(9988,'(1I10,3I6,1I10,2I6,1I10,2I6,2I10)') SURF_INDX(k),
     &                      dummy1, dummy2, dummy3,
     &                      dummy4, dummy5, dummy6,
     &                      dummy7, dummy8, dummy9,
     &                      II_INDX(k), JJ_INDX(k)
        ENDDO
        close(9988)
!cks  Write the headings for the file (so that it is human-readable)
        counter = 0
        open(unit=9997,file=eirene_surface_transfer_filename)
        write(9997,'(a)',advance='no')
     &'     ISURF     ISIDE     ITRIA        II        JJ SURF_INDX'
!cks  Corresponds to EOTAT
        DO IATMsurf=1,NATM
          write(9997,'(a,1I7)',advance='no') '   Inc.Power Atom ',
     &                                       IATMsurf
        ENDDO
!cks  Corresponds to EOTML
        DO IATMsurf=1,NMOL
          write(9997,'(a,1I7)',advance='no') '    Inc.Power Mol.',
     &                                       IATMsurf
        ENDDO
!cks  Corresponds to ERFPAT
        DO IATMsurf=1,NATM
          write(9997,'(a,1I7)',advance='no') '  Refl.Power BI->A',
     &                                       IATMsurf
        ENDDO
!cks  Corresponds to ERFPML
        DO IATMsurf=1,NMOL
          write(9997,'(a,1I7)',advance='no') '  Refl.Power BI->M',
     &                                       IATMsurf
        ENDDO
!cks  Corresponds to ERFPPHT
        DO IATMsurf=1,NPHOT
          write(9997,'(a,1I7)',advance='no') ' Refl.Power BI->Ph',
     &                                       IATMsurf
        ENDDO
!cks  Corresponds to ERFPIO
        DO IATMsurf=1,NION
          write(9997,'(a,1I7)',advance='no') ' Refl.Power BI->TI',
     &                                       IATMsurf
        ENDDO
!cks  Corresponds to ERFAAT
        DO IATMsurf=1,NATM
          write(9997,'(a,1I7)',advance='no') '   Refl.Power A->A',
     &                                       IATMsurf
        ENDDO
!cks  Corresponds to ERFAML
        DO IATMsurf=1,NMOL
          write(9997,'(a,1I7)',advance='no') '   Refl.Power A->M',
     &                                       IATMsurf
        ENDDO
!cks  Corresponds to ERFMAT
        DO IATMsurf=1,NATM
          write(9997,'(a,1I7)',advance='no') '   Refl.Power M->A',
     &                                       IATMsurf
        ENDDO
!cks  Corresponds to ERFMML
        DO IATMsurf=1,NMOL
          write(9997,'(a,1I7)',advance='no') '   Refl.Power M->M',
     &                                       IATMsurf
        ENDDO

        write(9997,'(a)') " ENDOFLINE"
!cks  Headings are now written.
!cks  Order is EOTAT(NATM), EOTML(NMOL), ERFPAT(NATM), ERFPML(NMOL),
!cks  ERFPPHT(NPHOT), ERFPIO(NION), ERFAAT(NATM), ERFAML(NMOL), ERFMAT(NATM), ERFMML(NMOL)
!cks  (read in EDGE2D: pf2ds/linkeirene.f in "read eirene.surface_transfer")

        DO is = 1, NSTS
          if(surftype(is)/=TARG) cycle
          k = NLIM + is
!cks  loop over all triangle edges
          DO np=1,3
!cks  loop over all triangles
            DO nr=1,NTRII
!cks  if the surface index matches the surface for the triangle-edge
!cks  combination in the inmti array, then we have the right triangle
!cks  edge indices
              IF (k .eq. INMTI(np,nr)) THEN
                 write(9997,'(6I10)',advance='no')  INMTI(np,nr),
     &              np,nr,
     &              II_INDX(nr),
     &              JJ_INDX(nr),
     &              NLIM+NSTS+INSPAT(np,nr)
!cks     .              ESTIMS(NADDW(31)+1,NLIM+NSTS+INSPAT(np,nr))
                 DO IATMsurf=1,NATM
                    write(9997,'(D25.16)',advance='no')
     &                  EOTAT(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 DO IATMsurf=1,NMOL
                    write(9997,'(D25.16)',advance='no')
     &                  EOTML(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 DO IATMsurf=1,NATM
                    write(9997,'(D25.16)',advance='no')
     &                  ERFPAT(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 DO IATMsurf=1,NMOL
                    write(9997,'(D25.16)',advance='no')
     &                  ERFPML(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 DO IATMsurf=1,NPHOT
                    write(9997,'(D25.16)',advance='no')
     &                  ERFPPHT(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 DO IATMsurf=1,NION
                    write(9997,'(D25.16)',advance='no')
     &                  ERFPIO(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 DO IATMsurf=1,NATM
                    write(9997,'(D25.16)',advance='no')
     &                  ERFAAT(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 DO IATMsurf=1,NMOL
                    write(9997,'(D25.16)',advance='no')
     &                  ERFAML(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 DO IATMsurf=1,NATM
                    write(9997,'(D25.16)',advance='no')
     &                  ERFMAT(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 DO IATMsurf=1,NMOL
                    write(9997,'(D25.16)',advance='no')
     &                  ERFMML(IATMsurf,NLIM+NSTS+INSPAT(np,nr))
                 ENDDO
                 write(9997,'(a)') " ENDOFLINE"
                 counter = counter + 1
               ENDIF
             ENDDO
          ENDDO
        ENDDO
        write(9997,'(a)') "Total data rows and columns: "
        write(9997,'(2I10)') counter, 6+4*NATM+4*NMOL+NPHOT
        deallocate(JJ_INDX)
        deallocate(II_INDX)
        deallocate(SURF_INDX)
        close(9997)
      ENDIF

!cks  --------------------eirene.surface_transfer creation--------------------
!cks  ------------------------------------------------------------------------


C     fill replicated sputtered flux arrays
      sptpl(:,:) = 0.d0
      sptpl(1,1:nlmpgs) = sptpltot(1:nlmpgs)
      sptat(:,:) = 0.d0
      sptat(1,1:nlmpgs) = sptatot(1:nlmpgs)
      sptml(:,:) = 0.d0
      sptml(1,1:nlmpgs) = sptmtot(1:nlmpgs)
      sptio(:,:) = 0.d0
      sptio(1,1:nlmpgs) = sptitot(1:nlmpgs)
      if (nphot > 0) then
        sptpht(:,:) = 0.d0
        sptpht(1,1:nlmpgs) = sptphtot(1:nlmpgs)
      end if

csw 25oct07
      allocate(sumpotpl(npls_fix))
csw
      do k=1,npls_fix
         sumpotpl(k) = 0.
         !if(nsurf==0) cycle
         do np=1,3
            do nr=1,ntrii
              if(inmti(np,nr)==0) cycle
              is = inmti(np,nr) - NLIM
              if(surftype(is)==TARG) then
                   msg = nlim+nsts+inspat(np,nr)

                   x1 = xtrian(necke(np,nr))
                   y1 = ytrian(necke(np,nr))
                   if(np == 3) then
                      x2 = xtrian(necke(1,nr))
                      y2 = ytrian(necke(1,nr))
                   else
                      x2 = xtrian(necke(np+1,nr))
                      y2 = ytrian(necke(np+1,nr))
                   endif

                   ar = sqrt( (x1-x2)**2 + (y1-y2)**2)
                   xc = (x1+x2)/2.d0
c     potpl in 1/s/cm
                   sumpotpl(k) = sumpotpl(k) +
     .                  estims(naddw(25)+k,msg)
     .                  /1.6022e-19/2.d0/pia/xc
              endif
            enddo
         enddo
      enddo
csw


      fp = 4999
      open(unit=fp,file=trim(eirene_ftransfer),access='sequential',
     .     status='replace')

      write(fp,'(a32,f14.6)') "* Neutral transfer file version:",
     &     NeutralTransferFileVersion

cnh   Writing out if EIRENE is running in time-dependent mode
      write(fp,'(a)') '* ntime  :'
      write(fp,'(1x,i6)') ntime

      write(fp,'(a)') '* ntrii,nrad  :'
      write(fp,'(3(1x,i6))') ntrii,nrad

      write(fp,'(a)') '* nstordr:'
      write(fp,'(1x,i6)') nstordr

      write(fp,'(a)') '* natm,nmol,nion,nphot  :'
      write(fp,'(4(1x,i6))') natm,nmol,nion,nphot

      write(fp,'(a)') '* npls:'
      write(fp,'(3(1x,i6))') npls_fix

      write(fp,'(a)') '* nlimps,nlim,nsts,nlmpgs:'
      write(fp,'(4(1x,i6))') nlimps,nlim,nsts,nlmpgs

c---------------------------------------
      write(fp,'(a,i6)') '* BULK SPECIES NPLS = ',npls_fix
      do ipls=1,npls_fix
         ityp=4
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IPLS = ',ipls
         do ir=1,ntrii
            c1 = mapl_vec(ipls,ir)*SIGN(1._DP,bv_vec(ipls,ir))
            c2 = mmpl_vec(ipls,ir)*SIGN(1._DP,bv_vec(ipls,ir))
            c3 = mipl_vec(ipls,ir)*SIGN(1._DP,bv_vec(ipls,ir))
c            c4 = mphpl(ipls,ir)
            c4 = 0.
C           Recombination contribution
            precom = 0.d0
            erecom = 0.d0
            DO IIRC=1,NPRCI(IPLS)
               IRRC=LGPRC(IPLS,IIRC)
               precom = precom - TABRC1(IRRC,IR)*DIIN(IPLS,IR)*ELCHA
               erecom = erecom + EELRC1(IRRC,IR)*DIIN(IPLS,IR)*ELCHA
            ENDDO
            ti = 0d0
            if(ipls==1) ti = tiin(ipls,ir)

            write(fp,'(i6,32(1x,e14.6))') ir,
     .           papl(ipls,ir),
     .           pmpl(ipls,ir),
     .           pipl(ipls,ir),
c     .           pphpl(ipls,ir),
     .           0.,
c pspl:
     .           papl(ipls,ir)+pmpl(ipls,ir)+pipl(ipls,ir)
c     .                        +pphpl(ipls,ir),
     .                        +0.,
c mapl:
     .           c1,
c mmpl:
     .           c2,
c mipl:
     .           c3,
c mphpl:
     .           c4,
c mspl:
     .           c1+c2+c3+c4,

c not used anymore... (changed background profiles in extra file?)
     .           diin(ipls,ir),
     .           vxin(ipls,ir),
     .           vyin(ipls,ir),
     .           vzin(ipls,ir),
     .           bv_vec(ipls,ir),
     .           ti,
     .           edrift(ipls,ir),

c     .           eapl(ipls,ir),empl(ipls,ir),eipl(ipls,ir),ephpl(ir),
     .           eapl(ipls,ir),empl(ipls,ir),eipl(ipls,ir),0.,
c espl:
c     .           eapl(ipls,ir)+empl(ipls,ir)+eipl(ipls,ir)+ephpl(ir),
     .           eapl(ipls,ir)+empl(ipls,ir)+eipl(ipls,ir)+0.,

c     .           eael(ir),emel(ir),eiel(ir),ephel(ir),
     .           eael(ir),emel(ir),eiel(ir),0.,
c esel:
c     .           eael(ir)+emel(ir)+eiel(ir)+ephel(ir)
     .           eael(ir)+emel(ir)+eiel(ir)+0.,
c particle recombination:
     .           precom,
c electron energy recombination:
     .           erecom

         enddo
         write(fp,'(a,i6)') ' SRF.AV. IPLS = ',ipls
         do is=1,nlmpgs
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potpl(ipls,is),
     .           eotpl(ipls,is),
     .           sptpl(ipls,is),
     .           spump(ispz,is)
         enddo
         write(fp,'(20(1x,e14.6))') sumpotpl(ipls),srec(ipls,0)
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* ATOMIC SPECIES, NATM = ',natm
      do iatm=1,natm
         ityp=1
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IATM = ',iatm
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdena(iatm,ir)*bxin(ir) + vydena(iatm,ir)
     .              *byin(ir) + vzdena(iatm,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir,
     .           pdena(iatm,ir),
     .           vxdena(iatm,ir),
     .           vydena(iatm,ir),
     .           vzdena(iatm,ir),
     .           vdenpara,
     .           edena(iatm,ir),
     .           edena(iatm,ir)/(pdena(iatm,ir)+1.d-60),
     .           0.,
     .           sigma(1,ir),
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IATM = ',iatm
         do is=1,nlmpgs
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potat(iatm,is),
     .           prfaat(iatm,is),
     .           prfmat(iatm,is),
     .           prfiat(iatm,is),
c     .           prfphat(iatm,is),
     .           0.,
     .           prfpat(iatm,is),
     .           eotat(iatm,is),
     .           erfaat(iatm,is),
     .           erfmat(iatm,is),
     .           erfiat(iatm,is),
c     .           erfphat(iatm,is),
     .           0.,
     .           erfpat(iatm,is),
     .           sptat(iatm,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* MOLECULAR SPECIES, NMOL = ',nmol
      do imol=1,nmol
         ityp=2
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IMOL = ',imol
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdenm(imol,ir)*bxin(ir) + vydenm(imol,ir)
     .              *byin(ir) + vzdenm(imol,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir,
     .           pdenm(imol,ir),
     .           vxdenm(imol,ir),
     .           vydenm(imol,ir),
     .           vzdenm(imol,ir),
     .           vdenpara,
     .           edenm(imol,ir),
     .           edenm(imol,ir)/(pdenm(imol,ir)+1.d-60),
     .           0.,
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IMOL = ',imol
         do is=1,nlmpgs
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potml(imol,is),
     .           prfaml(imol,is),
     .           prfmml(imol,is),
     .           prfiml(imol,is),
c     .           prfphml(imol,is),
     .           0.,
     .           prfpml(imol,is),
     .           eotml(imol,is),
     .           erfaml(imol,is),
     .           erfmml(imol,is),
     .           erfiml(imol,is),
c     .           erfphml(imol,is),
     .           0.,
     .           erfpml(imol,is),
     .           sptml(imol,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* IONIC SPECIES, NION = ',nion
      do iion=1,nion
         ityp=3
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IION = ',iion
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdeni(iion,ir)*bxin(ir) + vydeni(iion,ir)
     .              *byin(ir) + vzdeni(iion,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir,
     .           pdeni(iion,ir),
     .           vxdeni(iion,ir),
     .           vydeni(iion,ir),
     .           vzdeni(iion,ir),
     .           vdenpara,
     .           edeni(iion,ir),
     .           edeni(iion,ir)/(pdeni(iion,ir)+1.d-60),
     .           0.,
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IION = ',iion
         do is=1,nlmpgs
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potio(iion,is),
     .           prfaio(iion,is),
     .           prfmio(iion,is),
     .           prfiio(iion,is),
c     .           prfphio(iion,is),
     .           0.,
     .           prfpio(iion,is),
     .           eotio(iion,is),
     .           erfaio(iion,is),
     .           erfmio(iion,is),
     .           erfiio(iion,is),
c     .           erfphio(iion,is),
     .           0.,
     .           erfpio(iion,is),
     .           sptio(iion,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* PHOTONIC SPECIES, NPHOT = ',nphot
      do iphot=1,nphot
         ityp=0
         ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
         write(fp,'(a,i6)') ' VOL.AV. IPHOT = ',iphot
         do ir=1,ntrii
            bb = sqrt( bxin(ir)**2 + byin(ir)**2 + bzin(ir)**2 )
            if( bb > 0.) then
               vdenpara = vxdenph(iphot,ir)*bxin(ir) + vydenph(iphot,ir)
     .              *byin(ir) + vzdenph(iphot,ir)*bzin(ir)
               vdenpara = vdenpara / bb
            else
               vdenpara=0.
            endif
            write(fp,'(i6,20(1x,e14.6))') ir,
     .           pdenph(iphot,ir),
     .           vxdenph(iphot,ir),
     .           vydenph(iphot,ir),
     .           vzdenph(iphot,ir),
     .           vdenpara,
     .           edenph(iphot,ir),
     .           edenph(iphot,ir)/(pdenph(iphot,ir)+1.d-60),
     .           0.,
     .           0.,
     .           0.
         enddo
         write(fp,'(a,i6)') ' SRF.AV. IPHOT = ',iphot
         do is=1,nlmpgs
            write(fp,'(i6,20(1x,e14.6))') is,
     .           potpht(iphot,is),
     .           prfapht(iphot,is),
     .           prfmpht(iphot,is),
     .           prfipht(iphot,is),
     .           prfphpht(iphot,is),
     .           prfppht(iphot,is),
     .           eotpht(iphot,is),
     .           erfapht(iphot,is),
     .           erfmpht(iphot,is),
     .           erfipht(iphot,is),
     .           erfphpht(iphot,is),
     .           erfppht(iphot,is),
     .           sptpht(iphot,is),
     .           spump(ispz,is)
         enddo
      enddo

c---------------------------------------
      write(fp,'(a,i6)') '* MISC DATA'
      do ir=1,ntrii
         write(fp,'(i6,20(1x,e14.6))') ir,
     .        dble(ncltal(ir)),0.,0.,
     .        vol(ir),voltal(ir),
     .        dein(ir),tein(ir),bxin(ir),byin(ir),bzin(ir),bfin(ir)
      enddo

csw 20dec07-----------------------------
      write(fp,'(a)') '* STRATUM DATA'
      write(fp,'(i6)') nstra
      do istra=1,nstra
        do iatm=1,natm
          write(fp,'(2i6,20(1x,e14.6))') istra,iatm,
     .              wtota(iatm,istra),0.
        enddo
        do imol=1,nmol
          write(fp,'(2i6,20(1x,e14.6))') istra,imol,
     .              wtotm(imol,istra),0.
        enddo
        do iion=1,nion
          write(fp,'(2i6,20(1x,e14.6))') istra,iion,
     .              wtoti(iion,istra),0.
        enddo
        do iphot=1,nphot
          write(fp,'(2i6,20(1x,e14.6))') istra,iphot,
     .              wtotph(iphot,istra),0.
        enddo
        do ipls=1,npls_fix
          write(fp,'(2i6,20(1x,e14.6))') istra,ipls,
     .              wtotp(ipls,istra),0.
        enddo
      enddo

csw 24jan08
c---------------------------------------
      if(eirene_njetto .gt. 0) then
        write(fp,'(a)') '* ADDITIONAL TALLY DATA FOR JETTO NEUTRALS'
        write(fp,'(6i6)') nadv
        do ir=1,ntrii
          write(fp,'(i6,20(1x,e14.6))') ir,( addv(i,ir), i=1,nadv ),
     .                                     ( ppat(i,ir), i=1,natm )
        enddo
      endif

      if(eirene_nbirth .gt. 0) then
        write(fp,'(a)') '* ADDITIONAL TALLY DATA FOR JETTO NBI'
        write(fp,'(6i6)') nstra
        do is=1,nstra
          write(fp,'(i6,20(1x,e14.6))') is, etota(is)
        enddo
      endif

      close(fp)
csw 25oct07
      deallocate(sumpotpl)
csw

cdmh 15jun10
c----------------------------------------
c     store neutral particle fluxes [A] on wall
      fp = 4999
      open(unit=fp,file=trim(eirene_fstoreneutflux),access='sequential',
     .     status='replace')
      write(fp,'(a28,f14.6)') "* Neutral flux file version:",
     &     NeutralFluxFileVersion
      write(fp,'(a,a)') '*  NLIM,   NSTS,  NGITT, NGSTAL,',
     &     '   NATM,  NMOL, NLIMPS,  NTRII'
      write(fp,'(8i8)') NLIM, NSTS, NGITT, NGSTAL, NATM, NMOL, NLIMPS,
     &     NTRII

      do iatm=1,NATM
         write(fp,'(a,i0)') '* neutral fluxes from atom species ',iatm
         write(fp,'(a,a)') "*  idx,  neutral_flux, ",
     &        "SAREA, ITRIA, ISIDE, ISURF"
         do is=1,NLIMPS
            do j= 1, surf_trian(is)%numtr
              nr = surf_trian(is)%itrias(j)
              np = surf_trian(is)%itrisi(j)
              MSURFG=NLIM+NSTS+INSPAT(np,nr)
              write(fp,'(i6,1x,2(e14.6,1x),i6,1x,i6,1x,i6)')
     &        INSPAT(np,nr),POTAT(iatm,MSURFG),SAREA(MSURFG),
     &        nr,np,INMTI(np,nr)
            enddo
         enddo                  !is
      enddo                     !iatm

      do imol=1,NMOL
         write(fp,'(a,i0)') '* neutral fluxes from molecular species ',
     &        imol
         write(fp,'(a,a)') "*  idx,  neutral_flux, ",
     &        "SAREA, ITRIA, ISIDE, ISURF"
         do is=1,NLIMPS
            do j= 1, surf_trian(is)%numtr
              nr = surf_trian(is)%itrias(j)
              np = surf_trian(is)%itrisi(j)
              MSURFG=NLIM+NSTS+INSPAT(np,nr)
              write(fp,'(i6,1x,2(e14.6,1x),i6,1x,i6,1x,i6)')
     &        INSPAT(np,nr),POTML(imol,MSURFG),
     &        SAREA(MSURFG),nr,np,INMTI(np,nr)
            enddo
         enddo                  !is
      enddo                     !imol

      close(fp)

      call eirene_dealloc_usrdata()
      call eirene_prousr(dummy, -1, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0)
      ! for some saved LFIRST flags

      end subroutine eirene_outusr
