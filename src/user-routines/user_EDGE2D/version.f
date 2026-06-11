C> \brief Print EIRENEs git version hash.
C>
C> Prints/sends the date, branch and version hash of EIRENE to the
C> usual output unit.
      SUBROUTINE EIRENE_VERSION
      USE EIRENE_PROVENANCE
      USE EIRMOD_PARMMOD, ONLY: EIRENE_VERSION_STRING
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      EXTERNAL :: EIRENE_HEADNG, EIRENE_LEER

      WRITE(IUNOUT,'(2A)') ' EIRENE VERSION ',EIRENE_VERSION_STRING
      CALL EIRENE_LEER(2)

      CALL EIRENE_HEADNG('GIT REVISION USED IN THIS RUN',29)

      WRITE(IUNOUT,*) 'DESCRIPTION:',
     . git_describe
      WRITE(IUNOUT,*) 'DATE:       ',
     . git_commit_date
      WRITE(IUNOUT,*) 'HASH:       ',
     . git_head_sha1
      WRITE(IUNOUT,*) 'BRANCH:     ',
     . git_branch

      RETURN
      END SUBROUTINE EIRENE_VERSION
