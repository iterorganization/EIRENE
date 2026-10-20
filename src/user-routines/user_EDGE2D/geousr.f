      SUBROUTINE EIRENE_GEOUSR
C
C   PREPARE DATA FOR LIMITER SURFACES
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CADGEO
      USE EIRMOD_COMPRT
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTRIG
      USE EIRMOD_CCONA
      USE EIRMOD_CGEOM
      USE EIRMOD_CGRID
      USE EIRMOD_CLGIN
      USE EIRMOD_CINIT
      USE EIRMOD_CPOLYG
      USE EIRMOD_USRDATA
      IMPLICIT NONE
      CHARACTER(80) :: ZEILE
      REAL(DP) :: XCOOR, YCOOR, ZCOOR, xan, xen, yan, yen, zan, zen
      INTEGER :: NADMOD, NASMOD, NORMOD, NRS, IPUNKT, NSSIR, NSSIP,
     .           IDIR, IR, IP, IT, IC, IN, NAS
      logical :: lx1,lx2,lx3,lx4,ly1,ly2,ly3,ly4
      integer :: ists, i, j, is
C
C MODIFY GEOMETRY
C
C
cswx 24sep07
      ! coarse mesh
      do j=1,nrpla
        ncltal(j) = int(real(j-1)/icoarselev)+1
      enddo
      nsbox_tal = int(nrpla/icoarselev)
      do j=1,nr1st-nrpla
        ncltal(nrpla+j)=nsbox_tal+j
      enddo
      write(iunout,('(a,10i8)')),'geousr : coarse mesh',
     &                          nsbox,nsbox_tal,nrpla,nrad,nr1st
      nsbox_tal = nsbox_tal + nr1st - nrpla
      nr1tal    = nsbox_tal

!rf   fix INMTI and INMTINSS as should be done in grid.f
!     needed for diagnostic surfaces
      do ists = 1, nlimps
        if(iliin(ists)==0) cycle
        do j= 1, surf_trian(ists)%numtr
          i = surf_trian(ists)%itrias(j)
          is = surf_trian(ists)%itrisi(j)
          if (nchbar(is,i) > 0) then
            ! surface is seen by neighbor
            inmti(   nseite(is,i),nchbar(is,i)) = inmti(is,i)
            inmtinss(nseite(is,i),nchbar(is,i)) = -1
          end if
        end do
      end do

      return
C
C
C  ABSCHALTEN NICHT ERREICHBARER ODER DOPPELT VORHANDENER FLAECHEN
C
C   HIERHER: LGJUM1, LGJUM2 SETZEN ZUR BESCHLEUNIGUNG (NICHT UNBEDINGT
C   NOETIG)
C
C   LGJUM1(J,I)=.TRUE. :
C   ABSCHALTEN DER FLAECHE I, FALLS TEILCHEN AUF J SITZT
C
C   LGJUM2(J,I)=.TRUE. :
C   ABSCHALTEN DES ERSTEN SCHNITTPUNKTES MIT FLAECHE I, FALLS
C   TEILCHEN AUF J SITZT (FALLS I EINE FLAECHE ZWEITER ORDNUNG IST)
C
C   DEFAULTS: LGJUM1(J,J)=.TRUE. FUER EBENE FLAECHEN,
C             LGJUM2(J,J)=.TRUE. FUER FLAECHEN ZWEITER ORDNUNG
C
C

      RETURN
      END SUBROUTINE EIRENE_GEOUSR
