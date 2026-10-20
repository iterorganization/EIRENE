C  28.6.05
c Counting reactions involving photons (ityp=0 test particles)
C returns: only NRPH, for subr. find_param, setamd.
c currently: no EI, CX, PI, EL reaction for photons so far.
c Further NRPH reactions may be counted in XSECTP_PARAM, e.g. spont. photon emission
c         from bulk particles, iswr(kk)=6 reactions.

C
      SUBROUTINE EIRENE_XSECTPH_PARAM

      USE EIRMOD_PARMMOD
      USE EIRMOD_COMXS
      USE EIRMOD_COMUSR
      USE EIRMOD_PHOTON
      IMPLICIT NONE
csw
csw  PHOTON COLLISION, PH - type
csw
      integer :: iphot,nrc,kk

      do iphot=1,nphoti
         if(nrcph(iphot) > 0) then
            do nrc=1,nrcph(iphot)
               kk=ireacph(iphot,nrc)
               if(iswr(kk) == 7) then
                 NRPH=NRPH+1
               endif
            enddo
         endif
      enddo
cdr
      RETURN

      END SUBROUTINE EIRENE_XSECTPH_PARAM
