      subroutine eirene_lookup_adasdir_usr(DSN, FOUND,
     .                                     REAC, ELNAME, BUNDLING)
      use eirmod_precision
      use eirmod_parmmod, only: ifoff
      use eirmod_comprt, only: iunout
      use eirmod_cinit, only: master_path
      implicit none
      character(*),intent(inout) :: DSN
      character(*), intent(in) :: REAC
      logical, intent(inout) :: found
      character(*), intent(in), optional :: ELNAME, BUNDLING
      integer :: ierr
      logical :: found0, found1, found2
      EXTERNAL :: EIRENE_EXIT_OWN

!     If open fails, try again interpreting the path as relative
!     to SOLPSTOP/modules/Eirene
      inquire
     .  (FILE=trim(MASTER_PATH)//'/modules/Eirene/Database/'//TRIM(DSN),
     .   exist=found0)
      if (found0) then
        write (iunout,'(2a)') 'TAB1D OR TAB2D OR ADAS: ',
     >  trim(MASTER_PATH)//'/modules/Eirene/Database/'//trim(DSN)

        OPEN (UNIT=29+ifoff,
     .       FILE=trim(MASTER_PATH)//'/modules/Eirene/Database/'//
     .       TRIM(DSN),IOSTAT=ierr)
      else
        ierr = 1
      end if

      if ( ierr /= 0 .or. .not.found0) then
!     Next chance, try again with path relative to SOLPSTOP only
        inquire (FILE=trim(MASTER_PATH)//'/'//TRIM(DSN),
     .       exist=found1)
        if (found1) then
          write (iunout,'(2a)') 'TAB1D OR TAB2D OR ADAS: ',
     >        trim(MASTER_PATH)//'/'//trim(DSN)
          OPEN (UNIT=29+ifoff,
     .          FILE=trim(MASTER_PATH)//'/'//TRIM(DSN),IOSTAT=ierr)
        end if

        if ( ierr /= 0 .or. .not.found1) then
!       Try with default ADAS directory location in SOLPS-ITER
          IF (PRESENT(BUNDLING) .AND. VERIFY(BUNDLING,' ').NE.0) THEN
            DSN = ADJUSTL(TRIM(REAC)) // '/' //
     .            ADJUSTL(TRIM(REAC)) // '_' //
     .            TRIM(ELNAME) // '_' // TRIM(BUNDLING) // '.dat'
          ELSE
            DSN = ADJUSTL(TRIM(REAC)) // '/' //
     .            ADJUSTL(TRIM(REAC)) // '_' //
     .            TRIM(ELNAME) // '.dat'
          ENDIF
          inquire
     .      (FILE=trim(MASTER_PATH)//'/modules/adas/adf11/'//
     .            trim(DSN),exist=found2)
          if (found2) then
            write (iunout,'(2a)') 'TAB1D OR TAB2D OR ADAS: ',
     >           trim(MASTER_PATH)//'/modules/adas/adf11/'//trim(DSN)
            OPEN (UNIT=29+ifoff,
     .            FILE=trim(MASTER_PATH)//'/modules/adas/adf11/'//
     .                 trim(DSN),IOSTAT=ierr)
            if ( ierr /= 0 .or. .not.found2) then
              WRITE (iunout,'(a,a)')
     .             'SLREAC.F: UNABLE TO FIND REACTION FILE '//TRIM(DSN)
              CALL EIRENE_EXIT_OWN(1)
            endif
          endif
        endif
      endif

      found = found .or. found0 .or. found1 .or. found2

      return
      end subroutine eirene_lookup_adasdir_usr
