!pb  30.08.06: data structure for reaction data redefined
!pb  12.10.06: modcol revised
cdr  Jan. 2014:  minor bug in printout corrected
cdr  currently this routine is not called for RC processes, only for PH processes
cdr  Notation often just the other way round?

      SUBROUTINE EIRENE_XSTRC(ipls,nrc,idsc,irrc)
cdr
cdr  to replace       ISWR==6, non-default RC part, as called from XSECTP.F
cdr  already done for ISWR==7 (photonic reaction part, "lines").
cdr  Prepare volume recombination processes: bulk (+ bulk)--> test (+ bulk)
cdr  e.g.   ISWR==6  (RC)                    H+    +  e   --> H    (+ rad.)
cdr  e.g.   ISWR==7  (PH)                    H(n=2)       --> Ly-alpha (+H(n=1))
cdr
c    ipls: incident bulk
c    nrc : index of reaction in list of all reactions for field particle IPLS
c    idsc: index of RC/PH reaction for species ipls
c    irrr: index of RC/PH reaction in NRRC/NRPH arrays
c
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_COMXS
      USE EIRMOD_PHOTON

      IMPLICIT NONE
      integer, intent(in) :: ipls,nrc,idsc,irrc
      integer :: kk,ipl0,ipl1,ipl2,ityp0,ityp1,ityp2,
     .           j
      integer, external :: EIRENE_idez
      real(dp) :: factkk,aik
      EXTERNAL :: EIRENE_EXIT_OWN

c  fetch data for process nrc of ipls

      kk = IREACP(ipls,nrc)
c      if (kk /= idreac) call EIRENE_get_reaction(kk)

      FACTKK=FREACP(IPLS,NRC)
      IF (FACTKK.EQ.0.D0) FACTKK=1.
      aik=reacdat(kk)%phr%line%aik


      IPL0 =EIRENE_IDEZ(IBULKP(ipls,nrc),3,3)
      IPL1 =EIRENE_IDEZ(ISCD1P(ipls,nrc),3,3)
      IPL2 =EIRENE_IDEZ(ISCD2P(ipls,nrc),3,3)
      ITYP0=EIRENE_IDEZ(IBULKP(ipls,nrc),1,3)
      ITYP1=EIRENE_IDEZ(ISCD1P(ipls,nrc),1,3)
      ITYP2=EIRENE_IDEZ(ISCD2P(ipls,nrc),1,3)

      IF ((IPL0 < 1) .OR. (IPL0 > MAXSPC(ITYP0))) GOTO 993
      IF ((IPL1 < 1) .OR. (IPL1 > MAXSPC(ITYP1))) GOTO 994
      IF ((IPL2 < 1) .OR. (IPL2 > MAXSPC(ITYP2))) GOTO 994

      LGPRC(IPLS,IDSC)=IRRC

      facrrc(irrc,1) = factkk
      facrrc(irrc,2) = log(factkk)
      NREARC(irrc) = kk
      do j=1,nrad
         tabrc1(irrc,j)=aik*factkk
      enddo
      modcol(6,2,irrc)=1

      select case(ityp1)
      case(0)
         NPHPRC(IRRC)=IPL1
      case(1)
         NATPRC(IRRC)=IPL1
      case(2)
         NMLPRC(IRRC)=IPL1
      case(3)
         NIOPRC(IRRC)=IPL1
      case(4)
         NPLPRC(IRRC)=IPL1
      case default
         write(iunout,*) 'volume-processes: xstrc.f:'
         write(iunout,*) '   ityp1=',ityp1,' not allowed'
         call EIRENE_exit_own(1)
      end select
      select case(ityp2)
      case(0)
         NPHPRC_2(IRRC)=IPL2
      case(1)
         NATPRC_2(IRRC)=IPL2
      case(2)
         NMLPRC_2(IRRC)=IPL2
      case(3)
         NIOPRC_2(IRRC)=IPL2
      case(4)
         NPLPRC_2(IRRC)=IPL2
      case default
         write(iunout,*) 'volume-processes: xstrc.f:'
         write(iunout,*) '   ityp2=',ityp2,' not allowed'
         call EIRENE_exit_own(1)
      end select
      return

  993 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTRC: EXIT CALLED'
      WRITE (iunout,*)
     .  'SPECIES INDEX OF INCIDENT BULK PARTICLE OUT OF RANGE'
      WRITE (iunout,*) 'KK, ityp0, IPL0 ',
     .                  KK, ityp0, IPL0
      CALL EIRENE_EXIT_OWN(1)

  994 CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTRC: EXIT CALLED'
      WRITE (iunout,*)
     .  'SPECIES INDEX OF SECONDARY PARTICLE OUT OF RANGE'
      WRITE (iunout,*) 'KK, ityp0, IPL0 ',
     .                  KK, ityp0, IPL0
      WRITE (iunout,*) 'ityp1, IPL1, itpy2, IPL2 ',
     .                  ityp1, IPL1, ityp2, IPL2
      CALL EIRENE_EXIT_OWN(1)

      END SUBROUTINE EIRENE_XSTRC
