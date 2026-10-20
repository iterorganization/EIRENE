cdr  re-activated: Jan 2018
cdr  fix ph4:  cleanup: remove pressure broadening, polari,...
cdr            rename pointer to type(line_data): phline%...  to line_ir%...
cdr  Apr.22: The iprftype card is now already read in read_reaclines,
cdr          as were all other cards. So this exception is removed now. 

      subroutine EIRENE_read_photdbk (ir, reac, isw, iprftype)
c   read parameters relevant "reaction no IR" for line transport (photon gas transport)
c   from photonic database, into EIRENE data structure REACDAT(IR).
c   A photon (IPHOT) "line" is a sharp or broadened "line" photon
c   from a bound-bound light emission source.
c   A continuum emission from a given spectral source distribution is also
c   a "line" (=photon species), by abuse of language.
c

      use EIRMOD_precision
      use EIRMOD_parmmod
      USE EIRMOD_COMXS
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA
      USE EIRMOD_CINIT
      USE EIRMOD_PHOTON

      implicit none

      integer, intent(in) :: ir, isw, iprftype
      CHARACTER(50), INTENT(IN) :: REAC

      real(dp) :: wl, aik, ei, ej, b12, b21
      real(dp) :: rdata(9,1)

      integer :: gi, gj, inep, knep
      integer :: ianf, iend, iblnk, lr, ic,
     .           i1, lel
      character(1000) :: zeile
      character(20) :: elementname
      character(1) :: cha
      type(line_data), pointer :: phline
      external :: eirene_exit_own

      IF (REACDAT(IR)%LPHR) THEN
          WRITE (IUNOUT,*) ' PARAMETER FOR PHOTONIC REACTION ALREADY',
     .                     ' SPECIFIED FOR REACTION', IR
          WRITE (IUNOUT,*) ' CHECK SPECIFICATION OF REACTIONS'
          CALL EIRENE_EXIT_OWN(1)
      END IF

      lr=len_trim(reac)

      read (29+ifoff,*)
      read (29+ifoff,*)
      read (29+ifoff,*)

      do
        read (29+ifoff,'(A1000)',end=990) zeile
        if (zeile(1:2) == '--') cycle
        call EIRENE_chr_subcomma(zeile)

!  read element name
        ianf = 3
        iend = ianf + scan(zeile(ianf:),'|') - 1

        call EIRENE_chr_delete_blanks(zeile(ianf:iend))
        iblnk = scan(zeile(ianf:iend),'|') - 1
        if (iblnk < 0) iblnk = iend-ianf+1

        if (iblnk == 0) then
           write (iunout,*) ' ERROR IN DATABASE PHOTON'
           write (iunout,*) ' NO ELEMENT NAME FOUND'
           call eirene_exit_own(1)
        end if

        elementname = repeat(' ',20)
        elementname(1:iblnk) = zeile(ianf:ianf+iblnk-1)

!  read wavelength WL
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        read (zeile(ianf:iend-1),*) wl

!  read Zaehler
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        ic = ianf + verify(zeile(ianf:iend-1),' ') - 1
        cha = zeile(ic:ic)

        write (elementname(iblnk+1:),'(f10.4,a1)') wl,cha
        lel = len_trim(elementname)
        i1 = index(elementname,' ')
        do while (i1 < lel)
           elementname(i1:lel-1) = elementname(i1+1:lel)
           elementname(lel:lel) = ' '
           lel = lel - 1
           i1 = index(elementname,' ')
        end do

        if (elementname(1:iblnk+10) /= reac(1:iblnk+10)) cycle

!  skip reading 'transition'
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1

!  read aik
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        read (zeile(ianf:iend-1),*) aik

!  skip reading oscillator strength fij
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1

!  read gj
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        read (zeile(ianf:iend-1),*) gj

!  read gi
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        read (zeile(ianf:iend-1),*) gi

!  skip ll
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1


!  skip lu
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1

!  read ej
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          ej = 0._dp
        else
          read (zeile(ianf:iend-1),*) ej
        end if

!  read ei
        ianf = iend + 2
        iend = ianf + scan(zeile(ianf:),'|') - 1
        if (verify(zeile(ianf:iend-1),' ') == 0) then
          ei = 0._dp
        else
          read (zeile(ianf:iend-1),*) ei
        end if


        if (.true.) exit

      end do

      close (unit=29+ifoff)

c  reaction no IR is a "photonic" reaction
      reacdat(ir)%lphr = .true.

      allocate (reacdat(ir)%phr)
      allocate (reacdat(ir)%phr%line)
      nullify (reacdat(ir)%phr%adas)
      nullify (reacdat(ir)%phr%poly)
      nullify (reacdat(ir)%phr%tab1d)
      nullify (reacdat(ir)%phr%crm)

      allocate (phline)

      phline => reacdat(ir)%phr%line

      phline%aik = aik
!  line center wavelength
!     wl is in nm
!  line center energy in eV
      phline%e0 = hpcl / wl *1.E7_DP
!  stat. weights, upper, lower
      phline%g1 = gj
      phline%g2 = gi
!     Ej is in [1/cm] i.e. in 1.E-7 [1/nm]
cdr convert to eV units
      phline%e1 = ej * clight*hplanck*erg_to_ev

      select case (isw)
        case (1)    ! absorption
           phline%ircart = 4
        case (2)    ! emission
           phline%ircart = 4
        case (3)    ! stimulated emission
cdr  ??
        case (4)    ! effective absorption, corrected for stim. emiss.
           phline%ircart = 7
        case default
           phline%ircart = 0
      end select

      phline%iprofiletype = iprftype


cdr  jan 18: try to reconnect photonic data to reacdat structure.
cdr          not finished

      reacdat(ir)%phr%ifit = -1

cdr  fetch data for bound-bound transition line
c      call EIRENE_get_reaction(ir)
      b21=EIRENE_ph_b21(ir)
      b12=b21*gi/gj
      phline%b21 = b21
      phline%b12 = b12

      phline%reacname = reac(1:len_trim(reac))

cdr
c  rest of data: use reacdat(ir)%phr%poly, e.g. for Aik, and volumetric
c                                             source of photons. (RC process)
c  This is done by call to set_reaction_data (better name would be: "set_poly")
c  i.e. a single reaction IR can consist of OT and of RC processes.
c

      modclf(ir) = 100
c  constant rate (1/s)
      iftflg(ir,2) = 110

      rdata = 0._dp
      rdata(1,1) = aik

cdr
c  So far photonic cross-sections, rate coeff. and rates are constant.
c  i.e. special (trivial, 0th-order) cases of polygonial fits.
c  Use REACDAT type "poly" also for photonic data
      inep=1
      knep=1
      call EIRENE_set_reaction_data
     .  (ir,isw,iftflg(ir,2),rdata,inep,knep,
     .   iunout,.false.)

      return

  990 continue
      write (iunout,*) 'REACTION ',reac,' NOT FOUND IN FILE PHOTON'
      call EIRENE_exit_own(1)


      contains


      subroutine EIRENE_chr_subcomma (str)
cdr  replace comma "," with point "." in character string STR
      implicit none
      character(len=*), intent(in out) :: str
      integer :: i

      do
        i=scan(str,',')
        if (i == 0) exit
        str(i:i)='.'
      end do
      return
      end subroutine EIRENE_chr_subcomma


      subroutine EIRENE_chr_delete_blanks (str)
cdr  remove blanks from character string STR
      implicit none
      character(len=*), intent (in out) :: str
      character, allocatable :: compact(:)
      integer :: i, ic, l

      l=len(str)
      allocate (compact(l))
      compact=' '

      i=1
      ic=1
      do while (i<=l)
        if (str(i:i) /= ' ') then
          compact(ic) = str(i:i)
          ic=ic+1
          i=i+1
        else
          i=i+1
        end if
      end do

      do i=1,l
        str(i:i) = compact(i)
      end do

      return
      end subroutine EIRENE_chr_delete_blanks


      end subroutine EIRENE_read_photdbk
