

      SUBROUTINE EIRENE_CALSTR_USR(ipe,ic)
c called from broadcast/calstr
c can be used from user or case-specific routines (...usr.f, ...cop.f)
c to collect user/case-specific data from all
c PEs sharing calculations for a particular stratum
      IMPLICIT NONE
      integer, intent(in) :: ipe, ic
      RETURN
      END SUBROUTINE EIRENE_CALSTR_USR
