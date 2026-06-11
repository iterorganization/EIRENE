      SUBROUTINE EIRENE_SHOW_GIT_INFO(LOUT)
      
      use eirene_provenance
      implicit none
      integer, intent(in) :: LOUT

C................. Write out GIT specific configuration ................
C

      WRITE (LOUT,*) ""
      WRITE (LOUT,*) "EIRENE REPOSITORY STATUS DURING COMPILATION :-"
      WRITE (LOUT,*) "----------------------------------------------"

      WRITE(LOUT,"(a,a)") " Current GIT repository:  ",
     &trim(git_repository)
      WRITE(LOUT,"(a,a)") " Current GIT release tag: ",
     &git_release_tag
      WRITE(LOUT,"(a,a)") " Current GIT branch:      ",
     &git_branch
      WRITE(LOUT,"(a,a)") " Last commit SHA1-key:    ",
     &trim(git_head_sha1)
      WRITE(LOUT,"(a,a)") " Hostname at compilation: ",
     &trim(hostname)
      if (git_is_dirty.eq."true") then
         WRITE (LOUT,*) ""
         WRITE (LOUT,*) "********************************************"
         WRITE (LOUT,*) "***  WARNING:                            ***"
         WRITE (LOUT,*) "***  There were uncommitted changes in   ***"
         WRITE (LOUT,*) "***  the repository during compilation.  ***"
         WRITE (LOUT,*) "***  This executable may not be          ***"
         WRITE (LOUT,*) "***  reproducible by its SHA1-key!       ***"
         WRITE (LOUT,*) "********************************************"
      endif

      END SUBROUTINE EIRENE_SHOW_GIT_INFO
C
      SUBROUTINE EIRENE_SHOW_GIT_INFO_SHORT(LOUT)
C***********************************************************************
C   AUTHOR  :  Derek Harting (d.harting@fz-juelich.de)
C
C   PURPOSE :  Print out short format GIT repository status like repository
C              name, last commit SHA1-key and if the repsository was up to date.
C              The information is passed during compilation time via
C              the preprocessor in the __GIT variables.
C
C   INPUT   :   LOUT    - Unit number for output
C
C***********************************************************************
            
      use eirene_provenance
      use cmgutil_util, only: lenstr

      implicit none
      INTEGER, intent(in) :: LOUT
C
      CHARACTER :: CSTR*70

C
C................. Write out GIT specific configuration ................
C
C
      WRITE(LOUT,"(a,a)") " EIRENE GIT repository : ",
     &  trim(git_repository)
      CSTR=git_head_sha1
      if (git_is_dirty.eq."true") then
        CSTR=CSTR(1:LENSTR(CSTR))//" ( + uncommitted changes !! )"
      ENDIF
      WRITE(LOUT,"(a,a)") " EIRENE SHA1-key       : ",CSTR

      END  SUBROUTINE EIRENE_SHOW_GIT_INFO_SHORT
