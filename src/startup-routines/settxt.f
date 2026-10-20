c    2005 bug fix : text(ntalv-13) --> text(71)
c    17.03.06: txttal and txttlw added for additional tallies
cdr  29.09.14: TXTUNT corrected for generation limits, momentum sources
c    oct.14  : input tally 22 (potential) connnected to text arrays
cdr  dec. 15 : energy source tallies for bulk ions: additional species index ipls
cdr            tallies 38,44,50,56 and 84
cdr  dec.17:   pumped flux tally SPUMP: range 1--N5=NSPZ, rather than N7+1--N8
cdr            size of array LMETSPW decreased accordingly
cdr  june 18:  nlemis used to condition some storage setting (FOR REVISED BLOCK 12)
cdr  oct 18 :  setting text and storage range for input tallies: moved to own routines:
cdr            settxt_intal, and setprm_intal, to accommodate also the new input gradient tallies.

      SUBROUTINE EIRENE_SETTXT
c  Set default texts  (volume tallies: name, species, units),
C    ditto: surface and input tallies.
C  Main call: SETTXT
C  Set first (leading) dimension of tally arrays: nfstvi, nfstwi.
C  Subroutine: STTXT1
C  Set 1st index range per tally: nspan(itl), nspen(itl), for vol and surf. tallies,
c                                 for pointers to large tally arrays
c
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
      USE EIRMOD_CLOGAU

      IMPLICIT NONE

      INTEGER :: I, J, KK, ICO
      CHARACTER(24) :: TEXT24
      CHARACTER(72) :: TEXT72

      TXTTAL(1,1)='PARTICLE DENSITY (ATOMS)                         '
      TXTTAL(1,2)='PARTICLE DENSITY (MOLECULES)                     '
      TXTTAL(1,3)='PARTICLE DENSITY (TEST IONS)                     '
      TXTTAL(1,4)='PARTICLE DENSITY (PHOTONS)                       '
      TXTTAL(1,5)='ENERGY DENSITY (ATOMS)                           '
      TXTTAL(1,6)='ENERGY DENSITY (MOLECULES)                       '
      TXTTAL(1,7)='ENERGY DENSITY (TEST IONS)                       '
      TXTTAL(1,8)='ENERGY DENSITY (PHOTONS)                         '
      TXTTAL(1,9)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,10)=
     . 'PARTICLE SOURCE (ATOMS) FROM ATOM-PLASMA INTERACTION        '
      TXTTAL(1,11)=
     . 'PARTICLE SOURCE (MOLECULES) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,12)=
     . 'PARTICLE SOURCE (TEST IONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,13)=
     . 'PARTICLE SOURCE (PHOTONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,14)=
     . 'PARTICLE SOURCE (BULK IONS) FROM ATOM-PLASMA INTERACTION    '
      TXTTAL(1,15)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,16)=
     . 'PARTICLE SOURCE (ATOMS) FROM MOLECULE-PLASMA INTERACTION    '
      TXTTAL(1,17)=
     . 'PARTICLE SOURCE (MOLECULES) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,18)=
     . 'PARTICLE SOURCE (TEST IONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,19)=
     . 'PARTICLE SOURCE (PHOTONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,20)=
     . 'PARTICLE SOURCE (BULK IONS) FROM MOLECULE-PLASMA INTERACTION'
      TXTTAL(1,21)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,22)=
     . 'PARTICLE SOURCE (ATOMS) FROM TEST ION-PLASMA INTERACTION    '
      TXTTAL(1,23)=
     . 'PARTICLE SOURCE (MOLECULES) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,24)=
     . 'PARTICLE SOURCE (TEST IONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,25)=
     . 'PARTICLE SOURCE (PHOTONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,26)=
     . 'PARTICLE SOURCE (BULK IONS) FROM TEST ION-PLASMA INTERACTION'
      TXTTAL(1,27)=
     . 'PARTICLE SOURCE (ELECTRONS) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,28)=
     . 'PARTICLE SOURCE (ATOMS) FROM PHOTON-PLASMA INTERACTION      '
      TXTTAL(1,29)=
     . 'PARTICLE SOURCE (MOLECULES) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,30)=
     . 'PARTICLE SOURCE (TEST IONS) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,31)=
     . 'PARTICLE SOURCE (PHOTONS) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,32)=
     . 'PARTICLE SOURCE (BULK IONS) FROM PHOTON-PLASMA INTERACTION  '
      TXTTAL(1,33)=
     . 'ENERGY SOURCE (ELECTRONS) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,34)=
     . 'ENERGY SOURCE (ATOMS) FROM ATOM-PLASMA INTERACTION          '
      TXTTAL(1,35)=
     . 'ENERGY SOURCE (MOLECULES) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,36)=
     . 'ENERGY SOURCE (TEST IONS) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,37)=
     . 'ENERGY SOURCE (PHOTONS) FROM ATOM-PLASMA INTERACTION        '
      TXTTAL(1,38)=
     . 'ENERGY SOURCE (BULK IONS) FROM ATOM-PLASMA INTERACTION      '
      TXTTAL(1,39)=
     . 'ENERGY SOURCE (ELECTRONS) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,40)=
     . 'ENERGY SOURCE (ATOMS) FROM MOLECULE-PLASMA INTERACTION      '
      TXTTAL(1,41)=
     . 'ENERGY SOURCE (MOLECULES) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,42)=
     . 'ENERGY SOURCE (TEST IONS) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,43)=
     . 'ENERGY SOURCE (PHOTONS) FROM MOLECULE-PLASMA INTERACTION    '
      TXTTAL(1,44)=
     . 'ENERGY SOURCE (BULK IONS) FROM MOLECULE-PLASMA INTERACTION  '
      TXTTAL(1,45)=
     . 'ENERGY SOURCE (ELECTRONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,46)=
     . 'ENERGY SOURCE (ATOMS) FROM TEST ION-PLASMA INTERACTION      '
      TXTTAL(1,47)=
     . 'ENERGY SOURCE (MOLECULES) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,48)=
     . 'ENERGY SOURCE (TEST IONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,49)=
     . 'ENERGY SOURCE (PHOTONS) FROM TEST ION-PLASMA INTERACTION    '
      TXTTAL(1,50)=
     . 'ENERGY SOURCE (BULK IONS) FROM TEST ION-PLASMA INTERACTION  '
      TXTTAL(1,51)=
     . 'ENERGY SOURCE (ELECTRONS) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,52)=
     . 'ENERGY SOURCE (ATOMS) FROM PHOTON-PLASMA INTERACTION        '
      TXTTAL(1,53)=
     . 'ENERGY SOURCE (MOLECULES) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,54)=
     . 'ENERGY SOURCE (TEST IONS) FROM PHOTON-PLASMA INTERACTION    '
      TXTTAL(1,55)=
     . 'ENERGY SOURCE (PHOTONS) FROM PHOTON-PLASMA INTERACTION      '
      TXTTAL(1,56)=
     . 'ENERGY SOURCE (BULK IONS) FROM PHOTON-PLASMA INTERACTION    '
C  TALLY NTALA=57 (SEE PARMMOD.F)
C        ADDITIONAL TRACKLENGTH-ESTIMATED TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 10A
      TXTTAL(1,NTALA)=
     . 'ADDITIONAL TALLIES, TRACKLENGTH ESTIMATOR, SUBR. UPTUSR.F   '
C  TALLY NTALC=58 (SEE PARMMOD.F)
C        ADDITIONAL COLLISION-ESTIMATED TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 10B
      TXTTAL(1,NTALC)=
     . 'ADDITIONAL TALLIES, COLLISION ESTIMATOR, SUBR. UPCUSR.F     '
C  TALLY NTALT=59 (SEE PARMMOD.F)
C        ADDITIONAL SNAPSHOT-ESTIMATED TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 13B
      TXTTAL(1,NTALT)=
     . 'ADDITIONAL TALLIES, SNAPSHOT ESTIMATOR, SUBR. UPNUSR.F      '
C  TALLY NTALM=60 (SEE PARMMOD.F)
C        ADDITIONAL TALLIES FOR INTERFACING TO OTHER CODES
C        TXTTAL MAY BE OVERWRITTEN IN SUBR. INFCOP
      TXTTAL(1,NTALM)=
     . 'ADDITIONAL TALLIES FOR INTERFACING, SUBR. INFCOP.F          '
C  TALLY NTALB=61 (SEE PARMMOD.F)
C        ADDITIONAL TALLIES FOR ITERATIVE MODE (BGK ITERATION)
      TXTTAL(1,NTALB)=
     . 'ADDITIONAL TALLIES FOR ITERATIVE MODE, SUBR. UPTBGK.F       '
C  TALLY NTALR=62 (SEE PARMMOD.F)
C        ADDITIONAL TALLIES, ALGEBRAIC EXPRESSION IN EXISTING TALLIES
C        TXTTAL IS OVERWRITTEN BY INPUT BLOCK 10C
      TXTTAL(1,NTALR)=
     . 'ADDITIONAL TALLIES, ALGEBRAIC EXPRESSIONS, INPUT BLOCK 10C  '
C  GENERATION LIMIT TALLIES
      TXTTAL(1,63)=
     . 'PARTICLE SINK (ATOMS) DUE TO GENERATION LIMIT               '
      TXTTAL(1,64)=
     . 'PARTICLE SINK (MOLECULES) DUE TO GENERATION LIMIT           '
      TXTTAL(1,65)=
     . 'PARTICLE SINK (TEST IONS) DUE TO GENERATION LIMIT           '
      TXTTAL(1,66)=
     . 'PARTICLE SINK (PHOTONS) DUE TO GENERATION LIMIT             '
      TXTTAL(1,67)=
     . 'ENERGY SINK (ATOMS) DUE TO GENERATION LIMIT                 '
      TXTTAL(1,68)=
     . 'ENERGY SINK (MOLECULES) DUE TO GENERATION LIMIT             '
      TXTTAL(1,69)=
     . 'ENERGY SINK (TEST IONS) DUE TO GENERATION LIMIT             '
      TXTTAL(1,70)=
     . 'ENERGY SINK (PHOTONS) DUE TO GENERATION LIMIT               '
      TXTTAL(1,71)=
     . 'MOMENTUM SINK (ATOMS) DUE TO GENERATION LIMIT               '
      TXTTAL(1,72)=
     . 'MOMENTUM SINK (MOLECULES) DUE TO GENERATION LIMIT           '
      TXTTAL(1,73)=
     . 'MOMENTUM SINK (TEST IONS) DUE TO GENERATION LIMIT           '
      TXTTAL(1,74)=
     . 'MOMENTUM SINK (PHOTONS) DUE TO GENERATION LIMIT             '
C  VOLUMETRIC PRIMARY SOURCE TALLIES  (E.G. RECOMBINATION)
      TXTTAL(1,75)=
     . 'PRIMARY PARTICLE SOURCE (ATOMS) FROM PLASMA INTERACTIONS    '
      TXTTAL(1,76)=
     . 'PRIMARY PARTICLE SOURCE (MOLECULES) FROM PLASMA INTERACTIONS'
      TXTTAL(1,77)=
     . 'PRIMARY PARTICLE SOURCE (TEST IONS) FROM PLASMA INTERACTIONS'
      TXTTAL(1,78)=
     . 'PRIMARY PARTICLE SOURCE (PHOTONS) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,79)=
     . 'PRIMARY PARTICLE SOURCE (BULK IONS) FROM PLASMA INTERACTIONS'
      TXTTAL(1,80)=
     . 'PRIMARY ENERGY SOURCE (ATOMS) FROM PLASMA INTERACTIONS      '
      TXTTAL(1,81)=
     . 'PRIMARY ENERGY SOURCE (MOLECULES) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,82)=
     . 'PRIMARY ENERGY SOURCE (TEST IONS) FROM PLASMA INTERACTIONS  '
      TXTTAL(1,83)=
     . 'PRIMARY ENERGY SOURCE (PHOTONS) FROM PLASMA INTERACTIONS    '
      TXTTAL(1,84)=
     . 'PRIMARY ENERGY SOURCE (BULK IONS) FROM PLASMA INTERACTIONS  '
C  MOMENTUM DENSITY TALLIES
      TXTTAL(1,85)=
     . 'MOMENTUM DENSITY, X-DIRECTION (ATOMS)                       '
      TXTTAL(1,86)=
     . 'MOMENTUM DENSITY, X-DIRECTION (MOLECULES)                   '
      TXTTAL(1,87)=
     . 'MOMENTUM DENSITY, X-DIRECTION (TEST IONS)                   '
      TXTTAL(1,88)=
     . 'MOMENTUM DENSITY, X-DIRECTION (PHOTONS)                     '

      TXTTAL(1,89)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (ATOMS)                       '
      TXTTAL(1,90)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (MOLECULES)                   '
      TXTTAL(1,91)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (TEST IONS)                   '
      TXTTAL(1,92)=
     . 'MOMENTUM DENSITY, Y-DIRECTION (PHOTONS)                     '

      TXTTAL(1,93)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (ATOMS)                       '
      TXTTAL(1,94)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (MOLECULES)                   '
      TXTTAL(1,95)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (TEST IONS)                   '
      TXTTAL(1,96)=
     . 'MOMENTUM DENSITY, Z-DIRECTION (PHOTONS)                     '
cdr  vectorial momentum source tallies, wrt. B field
      ICO=0
      KK=ICO*NPLS
      TXTTAL(KK+1,97)=
     . 'MOMENTUM SOURCE FROM ATOM-PLASMA INTERACTION, PARALLEL      '
      TXTTAL(KK+1,98)=
     . 'MOMENTUM SOURCE FROM MOLECULE-PLASMA INTERACTION, PARALLEL  '
      TXTTAL(KK+1,99)=
     . 'MOMENTUM SOURCE FROM TEST ION-PLASMA INTERACTION, PARALLEL  '
      TXTTAL(KK+1,100)=
     . 'MOMENTUM SOURCE FROM PHOTON-PLASMA INTERACTION, PARALLEL    '
      ICO=1
      KK=ICO*NPLS
      TXTTAL(KK+1,97)=
     . 'MOMENTUM SOURCE FROM ATOM-PLASMA INTERACTION, PERP.         '
      TXTTAL(KK+1,98)=
     . 'MOMENTUM SOURCE FROM MOLECULE-PLASMA INTERACTION, PERP.     '
      TXTTAL(KK+1,99)=
     . 'MOMENTUM SOURCE FROM TEST ION-PLASMA INTERACTION, PERP.     '
      TXTTAL(KK+1,100)=
     . 'MOMENTUM SOURCE FROM PHOTON-PLASMA INTERACTION, PERP.       '
      ICO=2
      KK=ICO*NPLS
      TXTTAL(KK+1,97)=
     . 'MOMENTUM SOURCE FROM ATOM-PLASMA INTERACTION, DIAMAG.       '
      TXTTAL(KK+1,98)=
     . 'MOMENTUM SOURCE FROM MOLECULE-PLASMA INTERACTION, DIAMAG.   '
      TXTTAL(KK+1,99)=
     . 'MOMENTUM SOURCE FROM TEST ION-PLASMA INTERACTION, DIAMAG.   '
      TXTTAL(KK+1,100)=
     . 'MOMENTUM SOURCE FROM PHOTON-PLASMA INTERACTION, DIAMAG.     '

C  RADIATION RATES
      TXTTAL(1,101)=
     . 'RADIATION RATE EMITTED FROM ATOMS                           '
      TXTTAL(1,102)=
     . 'RADIATION RATE EMITTED FROM MOLECULES                       '
      TXTTAL(1,103)=
     . 'RADIATION RATE EMITTED FROM TEST IONS                       '
C
      DO 1 J=1,NTALV
        DO I=2,N1MX  ! rather than n1mx we should use nfstpi(j)
          IF (J.NE.97 .AND. J.NE.98 .AND. J.NE.99 .AND. J.NE.100) THEN
cdr  scalar output tallies 1...96
            TEXT72=TXTTAL(1,J)
            TXTTAL(I,J)=TEXT72
          ELSE
cdr  vectorial momentum source tallies
            ICO=(I-1)/NPLS  ! = 0,1,2
            KK=ICO*NPLS
            TEXT72=TXTTAL(KK+1,J)
            TXTTAL(I,J)=TEXT72
          ENDIF
        END DO
    1 CONTINUE
C
      TXTUNT(1,1)='CM**-3                  '
      TXTUNT(1,2)='CM**-3                  '
      TXTUNT(1,3)='CM**-3                  '
      TXTUNT(1,4)='CM**-3                  '
      TXTUNT(1,5)='EV*CM**-3               '
      TXTUNT(1,6)='EV*CM**-3               '
      TXTUNT(1,7)='EV*CM**-3               '
      TXTUNT(1,8)='EV*CM**-3               '
      TXTUNT(1,9)='AMP*CM**-3              '
      TXTUNT(1,10)='AMP*CM**-3              '
      TXTUNT(1,11)='AMP*CM**-3              '
      TXTUNT(1,12)='AMP*CM**-3              '
      TXTUNT(1,13)='AMP*CM**-3              '
      TXTUNT(1,14)='AMP*CM**-3              '
      TXTUNT(1,15)='AMP*CM**-3              '
      TXTUNT(1,16)='AMP*CM**-3              '
      TXTUNT(1,17)='AMP*CM**-3              '
      TXTUNT(1,18)='AMP*CM**-3              '
      TXTUNT(1,19)='AMP*CM**-3              '
      TXTUNT(1,20)='AMP*CM**-3              '
      TXTUNT(1,21)='AMP*CM**-3              '
      TXTUNT(1,22)='AMP*CM**-3              '
      TXTUNT(1,23)='AMP*CM**-3              '
      TXTUNT(1,24)='AMP*CM**-3              '
      TXTUNT(1,25)='AMP*CM**-3              '
      TXTUNT(1,26)='AMP*CM**-3              '
      TXTUNT(1,27)='AMP*CM**-3              '
      TXTUNT(1,28)='AMP*CM**-3              '
      TXTUNT(1,29)='AMP*CM**-3              '
      TXTUNT(1,30)='AMP*CM**-3              '
      TXTUNT(1,31)='AMP*CM**-3              '
      TXTUNT(1,32)='AMP*CM**-3              '
      TXTUNT(1,33)='WATT*CM**-3             '
      TXTUNT(1,34)='WATT*CM**-3             '
      TXTUNT(1,35)='WATT*CM**-3             '
      TXTUNT(1,36)='WATT*CM**-3             '
      TXTUNT(1,37)='WATT*CM**-3             '
      TXTUNT(1,38)='WATT*CM**-3             '
      TXTUNT(1,39)='WATT*CM**-3             '
      TXTUNT(1,40)='WATT*CM**-3             '
      TXTUNT(1,41)='WATT*CM**-3             '
      TXTUNT(1,42)='WATT*CM**-3             '
      TXTUNT(1,43)='WATT*CM**-3             '
      TXTUNT(1,44)='WATT*CM**-3             '
      TXTUNT(1,45)='WATT*CM**-3             '
      TXTUNT(1,46)='WATT*CM**-3             '
      TXTUNT(1,47)='WATT*CM**-3             '
      TXTUNT(1,48)='WATT*CM**-3             '
      TXTUNT(1,49)='WATT*CM**-3             '
      TXTUNT(1,50)='WATT*CM**-3             '
      TXTUNT(1,51)='WATT*CM**-3             '
      TXTUNT(1,52)='WATT*CM**-3             '
      TXTUNT(1,53)='WATT*CM**-3             '
      TXTUNT(1,54)='WATT*CM**-3             '
      TXTUNT(1,55)='WATT*CM**-3             '
      TXTUNT(1,56)='WATT*CM**-3             '
C  ADDITIONAL TALLIES
      TXTUNT(1,NTALA)='TO BE READ              '
      TXTUNT(1,NTALC)='TO BE READ              '
      TXTUNT(1,NTALT)='TO BE READ              '
      TXTUNT(1,NTALM)='TO BE DEFINED IN INFCOP '
      TXTUNT(1,NTALB)='TO BE DEFINED IN BGK    '
      TXTUNT(1,NTALR)='TO BE READ              '
C  GENERATION LIMIT TALLIES
      TXTUNT(1,63)='AMP*CM**-3              '
      TXTUNT(1,64)='AMP*CM**-3              '
      TXTUNT(1,65)='AMP*CM**-3              '
      TXTUNT(1,66)='AMP*CM**-3              '
      TXTUNT(1,67)='WATT*CM**-3             '
      TXTUNT(1,68)='WATT*CM**-3             '
      TXTUNT(1,69)='WATT*CM**-3             '
      TXTUNT(1,70)='WATT*CM**-3             '
      TXTUNT(1,71)='G*CM/S*AMP*CM**-3       '
      TXTUNT(1,72)='G*CM/S*AMP*CM**-3       '
      TXTUNT(1,73)='G*CM/S*AMP*CM**-3       '
      TXTUNT(1,74)='G*CM/S*AMP*CM**-3       '
C  VOLUMETRIC PRIMARY SOURCE TALLIES  (E.G. RECOMBINATION)
      TXTUNT(1,75)='AMP*CM**-3              '
      TXTUNT(1,76)='AMP*CM**-3              '
      TXTUNT(1,77)='AMP*CM**-3              '
      TXTUNT(1,78)='AMP*CM**-3              '
      TXTUNT(1,79)='AMP*CM**-3              '
      TXTUNT(1,80)='WATT*CM**-3             '
      TXTUNT(1,81)='WATT*CM**-3             '
      TXTUNT(1,82)='WATT*CM**-3             '
      TXTUNT(1,83)='WATT*CM**-3             '
      TXTUNT(1,84)='WATT*CM**-3             '
C  MOMENTUM DENSITY TALLIES, X,Y,Z
      TXTUNT(1,85)='G*CM/SEC*CM**-3         '
      TXTUNT(1,86)='G*CM/SEC*CM**-3         '
      TXTUNT(1,87)='G*CM/SEC*CM**-3         '
      TXTUNT(1,88)='G*CM/SEC*CM**-3         '

      TXTUNT(1,89)='G*CM/SEC*CM**-3         '
      TXTUNT(1,90)='G*CM/SEC*CM**-3         '
      TXTUNT(1,91)='G*CM/SEC*CM**-3         '
      TXTUNT(1,92)='G*CM/SEC*CM**-3         '

      TXTUNT(1,93)='G*CM/SEC*CM**-3         '
      TXTUNT(1,94)='G*CM/SEC*CM**-3         '
      TXTUNT(1,95)='G*CM/SEC*CM**-3         '
      TXTUNT(1,96)='G*CM/SEC*CM**-3         '
C  MOMENTUM SOURCES (BULK IONS)
      TXTUNT(1,97)='G*CM/S*AMP*CM**-3       '
      TXTUNT(1,98)='G*CM/S*AMP*CM**-3       '
      TXTUNT(1,99)='G*CM/S*AMP*CM**-3       '
      TXTUNT(1,100)='G*CM/S*AMP*CM**-3       '
C  RADIATION RATES
      TXTUNT(1,101)='WATT*CM**-3             '
      TXTUNT(1,102)='WATT*CM**-3             '
      TXTUNT(1,103)='WATT*CM**-3             '
      DO 2 J=1,NTALV
        DO I=2,N1MX
          TEXT24=TXTUNT(1,J)
          TXTUNT(I,J)=TEXT24
        END DO
    2 CONTINUE

      IF (NADVI > 0) THEN
        TXTTAL(1:NADVI,NTALA) = TXTTLA(1:NADVI)
        TXTSPC(1:NADVI,NTALA) = TXTSCA(1:NADVI)
        TXTUNT(1:NADVI,NTALA) = TXTUTA(1:NADVI)
      END IF

      IF (NCLVI > 0) THEN
        TXTTAL(1:NCLVI,NTALC) = TXTTLC(1:NCLVI)
        TXTSPC(1:NCLVI,NTALC) = TXTSCC(1:NCLVI)
        TXTUNT(1:NCLVI,NTALC) = TXTUTC(1:NCLVI)
      END IF

      IF (NALVI > 0) THEN
        TXTTAL(1:NALVI,NTALR) = TXTTLR(1:NALVI)
        TXTSPC(1:NALVI,NTALR) = TXTSCR(1:NALVI)
        TXTUNT(1:NALVI,NTALR) = TXTUTR(1:NALVI)
      END IF

      IF (NSNVI > 0) THEN
        TXTTAL(1:NSNVI,NTALT) = TXTTLT(1:NSNVI)
        TXTSPC(1:NSNVI,NTALT) = TXTSCT(1:NSNVI)
        TXTUNT(1:NSNVI,NTALT) = TXTUTT(1:NSNVI)
      END IF

      CALL EIRENE_DEALLOC_CTEXT3

C  SURFACE-AVERAGED TALLIES

C  PARTICLE FLUXES, INCIDENT AND EMITTED
      TXTTLW(1,1)='PARTICLE FLUX, INCIDENT, ATOMS                   '
      TXTTLW(1,2)='PARTICLE FLUX, EMITTED, ATS. => ATOMS            '
      TXTTLW(1,3)='PARTICLE FLUX, EMITTED, MLS. => ATOMS            '
      TXTTLW(1,4)='PARTICLE FLUX, EMITTED, T.I. => ATOMS            '
      TXTTLW(1,5)='PARTICLE FLUX, EMITTED, PHS. => ATOMS            '
      TXTTLW(1,6)='PARTICLE FLUX, EMITTED, B.I. => ATOMS            '
      TXTTLW(1,7)='PARTICLE FLUX, INCIDENT, MOLECULES               '
      TXTTLW(1,8)='PARTICLE FLUX, EMITTED, ATS. => MOLECULES        '
      TXTTLW(1,9)='PARTICLE FLUX, EMITTED, MLS. => MOLECULES        '
      TXTTLW(1,10)='PARTICLE FLUX, EMITTED, T.I. => MOLECULES        '
      TXTTLW(1,11)='PARTICLE FLUX, EMITTED, PHS. => MOLECULES        '
      TXTTLW(1,12)='PARTICLE FLUX, EMITTED, B.I. => MOLECULES        '
      TXTTLW(1,13)='PARTICLE FLUX, INCIDENT, TEST IONS               '
      TXTTLW(1,14)='PARTICLE FLUX, EMITTED, ATS. => TEST IONS        '
      TXTTLW(1,15)='PARTICLE FLUX, EMITTED, MLS. => TEST IONS        '
      TXTTLW(1,16)='PARTICLE FLUX, EMITTED, T.I. => TEST IONS        '
      TXTTLW(1,17)='PARTICLE FLUX, EMITTED, PHS. => TEST IONS        '
      TXTTLW(1,18)='PARTICLE FLUX, EMITTED, B.I. => TEST IONS        '
      TXTTLW(1,19)='PARTICLE FLUX, INCIDENT, PHOTONS                 '
      TXTTLW(1,20)='PARTICLE FLUX, EMITTED, ATS. => PHOTONS          '
      TXTTLW(1,21)='PARTICLE FLUX, EMITTED, MLS. => PHOTONS          '
      TXTTLW(1,22)='PARTICLE FLUX, EMITTED, T.I. => PHOTONS          '
      TXTTLW(1,23)='PARTICLE FLUX, EMITTED, PHS. => PHOTONS          '
      TXTTLW(1,24)='PARTICLE FLUX, EMITTED, B.I. => PHOTONS          '
      TXTTLW(1,25)='PARTICLE FLUX, INCIDENT, BULK IONS               '
C  ENERGY FLUXES, INCIDENT AND EMITTED
      TXTTLW(1,26)='ENERGY FLUX, INCIDENT, ATOMS                     '
      TXTTLW(1,27)='ENERGY FLUX, EMITTED, ATS. => ATOMS              '
      TXTTLW(1,28)='ENERGY FLUX, EMITTED, MLS. => ATOMS              '
      TXTTLW(1,29)='ENERGY FLUX, EMITTED, T.I. => ATOMS              '
      TXTTLW(1,30)='ENERGY FLUX, EMITTED, PHS. => ATOMS              '
      TXTTLW(1,31)='ENERGY FLUX, EMITTED, B.I. => ATOMS              '
      TXTTLW(1,32)='ENERGY FLUX, INCIDENT, MOLECULES                 '
      TXTTLW(1,33)='ENERGY FLUX, EMITTED, ATS. => MOLECULES          '
      TXTTLW(1,34)='ENERGY FLUX, EMITTED, MLS. => MOLECULES          '
      TXTTLW(1,35)='ENERGY FLUX, EMITTED, T.I. => MOLECULES          '
      TXTTLW(1,36)='ENERGY FLUX, EMITTED, PHS. => MOLECULES          '
      TXTTLW(1,37)='ENERGY FLUX, EMITTED, B.I. => MOLECULES          '
      TXTTLW(1,38)='ENERGY FLUX, INCIDENT, TEST IONS                 '
      TXTTLW(1,39)='ENERGY FLUX, EMITTED, ATS. => TEST IONS          '
      TXTTLW(1,40)='ENERGY FLUX, EMITTED, MLS. => TEST IONS          '
      TXTTLW(1,41)='ENERGY FLUX, EMITTED, T.I. => TEST IONS          '
      TXTTLW(1,42)='ENERGY FLUX, EMITTED, PHS. => TEST IONS          '
      TXTTLW(1,43)='ENERGY FLUX, EMITTED, B.I. => TEST IONS          '
      TXTTLW(1,44)='ENERGY FLUX, INCIDENT, PHOTONS                   '
      TXTTLW(1,45)='ENERGY FLUX, EMITTED, ATS. => PHOTONS            '
      TXTTLW(1,46)='ENERGY FLUX, EMITTED, MLS. => PHOTONS            '
      TXTTLW(1,47)='ENERGY FLUX, EMITTED, T.I. => PHOTONS            '
      TXTTLW(1,48)='ENERGY FLUX, EMITTED, PHS. => PHOTONS            '
      TXTTLW(1,49)='ENERGY FLUX, EMITTED, B.I. => PHOTONS            '
      TXTTLW(1,50)='ENERGY FLUX, INCIDENT, BULK IONS                 '
C  SPUTTERED FLUXES,  EMITTED FROM SURFACE
      TXTTLW(1,51)='SPUTTERED FLUX BY INCIDENT ATS. => ATOMS         '
      TXTTLW(1,52)='SPUTTERED FLUX BY INCIDENT MLS. => ATOMS         '
      TXTTLW(1,53)='SPUTTERED FLUX BY INCIDENT T.I. => ATOMS         '
      TXTTLW(1,54)='SPUTTERED FLUX BY INCIDENT PHS. => ATOMS         '
      TXTTLW(1,55)='SPUTTERED FLUX BY INCIDENT B.I. => ATOMS         '
      TXTTLW(1,56)='SPUTTERED FLUX BY INCIDENT ATS. => MOLECULES     '
      TXTTLW(1,57)='SPUTTERED FLUX BY INCIDENT MLS. => MOLECULES     '
      TXTTLW(1,58)='SPUTTERED FLUX BY INCIDENT T.I. => MOLECULES     '
      TXTTLW(1,59)='SPUTTERED FLUX BY INCIDENT PHS. => MOLECULES     '
      TXTTLW(1,60)='SPUTTERED FLUX BY INCIDENT B.I. => MOLECULES     '
      TXTTLW(1,61)='SPUTTERED FLUX BY INCIDENT ATS. => TEST IONS     '
      TXTTLW(1,62)='SPUTTERED FLUX BY INCIDENT MLS. => TEST IONS     '
      TXTTLW(1,63)='SPUTTERED FLUX BY INCIDENT T.I. => TEST IONS     '
      TXTTLW(1,64)='SPUTTERED FLUX BY INCIDENT PHS. => TEST IONS     '
      TXTTLW(1,65)='SPUTTERED FLUX BY INCIDENT B.I. => TEST IONS     '
      TXTTLW(1,66)='SPUTTERED FLUX BY INCIDENT ATS. => PHOTONS       '
      TXTTLW(1,67)='SPUTTERED FLUX BY INCIDENT MLS. => PHOTONS       '
      TXTTLW(1,68)='SPUTTERED FLUX BY INCIDENT T.I. => PHOTONS       '
      TXTTLW(1,69)='SPUTTERED FLUX BY INCIDENT PHS. => PHOTONS       '
      TXTTLW(1,70)='SPUTTERED FLUX BY INCIDENT B.I. => PHOTONS       '
      TXTTLW(1,71)='SPUTTERED FLUX BY INCIDENT ATS. => BULK IONS     '
      TXTTLW(1,72)='SPUTTERED FLUX BY INCIDENT MLS. => BULK IONS     '
      TXTTLW(1,73)='SPUTTERED FLUX BY INCIDENT T.I. => BULK IONS     '
      TXTTLW(1,74)='SPUTTERED FLUX BY INCIDENT PHS. => BULK IONS     '
      TXTTLW(1,75)='SPUTTERED FLUX BY INCIDENT B.I. => BULK IONS     '
      TXTTLW(1,76)='SPUTTERED FLUX BY INCIDENT ATS., TOTAL           '
      TXTTLW(1,77)='SPUTTERED FLUX BY INCIDENT MLS., TOTAL           '
      TXTTLW(1,78)='SPUTTERED FLUX BY INCIDENT T.I., TOTAL           '
      TXTTLW(1,79)='SPUTTERED FLUX BY INCIDENT PHS., TOTAL           '
      TXTTLW(1,80)='SPUTTERED FLUX BY INCIDENT B.I., TOTAL           '
      TXTTLW(1,81)='SPUTTERED FLUX, TOTAL                            '
C  TALLY NTLSA=82 (SEE PARMMOD.F)
C   ADDIT. TALLIES, SUBR. UPSUSR.F
C   TXTTLW IS OVERWRITTEN BY INPUT BLOCK 10D
      TXTTLW(1,82)='ADDITIONAL SURFACE TALLY, SUBR. UPSUSR.F         '
C  TALLY NTLSR=83 (SEE PARMMOD.F)
C   ADDIT. TALLIES, ALGEBRAIC EXPRESSION IN EXISTING TALLIES
C   TXTTLW IS OVERWRITTEN BY INPUT BLOCK 10E
      TXTTLW(1,83)='ALGEBRAIC EXPRESSION IN SURFACE-AVERAGED TALLIES '
C  PUMPED FLUXES
      TXTTLW(1,84)='PUMPED FLUX BY SPECIES                           '
C
      DO J=1,NTALS
        TXTTLW(2:N2MX,J)=TXTTLW(1,J)
      END DO
C  particle fluxes
      TXTUNW(1,1)='AMP                     '
      TXTUNW(1,2)='AMP                     '
      TXTUNW(1,3)='AMP                     '
      TXTUNW(1,4)='AMP                     '
      TXTUNW(1,5)='AMP                     '
      TXTUNW(1,6)='AMP                     '
      TXTUNW(1,7)='AMP                     '
      TXTUNW(1,8)='AMP                     '
      TXTUNW(1,9)='AMP                     '
      TXTUNW(1,10)='AMP                     '
      TXTUNW(1,11)='AMP                     '
      TXTUNW(1,12)='AMP                     '
      TXTUNW(1,13)='AMP                     '
      TXTUNW(1,14)='AMP                     '
      TXTUNW(1,15)='AMP                     '
      TXTUNW(1,16)='AMP                     '
      TXTUNW(1,17)='AMP                     '
      TXTUNW(1,18)='AMP                     '
      TXTUNW(1,19)='AMP                     '
      TXTUNW(1,20)='AMP                     '
      TXTUNW(1,21)='AMP                     '
      TXTUNW(1,22)='AMP                     '
      TXTUNW(1,23)='AMP                     '
      TXTUNW(1,24)='AMP                     '
      TXTUNW(1,25)='AMP                     '
c  energy fluxes
      TXTUNW(1,26)='WATT                    '
      TXTUNW(1,27)='WATT                    '
      TXTUNW(1,28)='WATT                    '
      TXTUNW(1,29)='WATT                    '
      TXTUNW(1,30)='WATT                    '
      TXTUNW(1,31)='WATT                    '
      TXTUNW(1,32)='WATT                    '
      TXTUNW(1,33)='WATT                    '
      TXTUNW(1,34)='WATT                    '
      TXTUNW(1,35)='WATT                    '
      TXTUNW(1,36)='WATT                    '
      TXTUNW(1,37)='WATT                    '
      TXTUNW(1,38)='WATT                    '
      TXTUNW(1,39)='WATT                    '
      TXTUNW(1,40)='WATT                    '
      TXTUNW(1,41)='WATT                    '
      TXTUNW(1,42)='WATT                    '
      TXTUNW(1,43)='WATT                    '
      TXTUNW(1,44)='WATT                    '
      TXTUNW(1,45)='WATT                    '
      TXTUNW(1,46)='WATT                    '
      TXTUNW(1,47)='WATT                    '
      TXTUNW(1,48)='WATT                    '
      TXTUNW(1,49)='WATT                    '
      TXTUNW(1,50)='WATT                    '
c  sputter tallies
      TXTUNW(1,51)='AMP                     '
      TXTUNW(1,52)='AMP                     '
      TXTUNW(1,53)='AMP                     '
      TXTUNW(1,54)='AMP                     '
      TXTUNW(1,55)='AMP                     '
      TXTUNW(1,56)='AMP                     '
      TXTUNW(1,51)='AMP                     '
      TXTUNW(1,52)='AMP                     '
      TXTUNW(1,53)='AMP                     '
      TXTUNW(1,54)='AMP                     '
      TXTUNW(1,55)='AMP                     '
      TXTUNW(1,56)='AMP                     '
      TXTUNW(1,57)='AMP                     '
      TXTUNW(1,58)='AMP                     '
      TXTUNW(1,59)='AMP                     '
      TXTUNW(1,60)='AMP                     '
      TXTUNW(1,61)='AMP                     '
      TXTUNW(1,62)='AMP                     '
      TXTUNW(1,63)='AMP                     '
      TXTUNW(1,64)='AMP                     '
      TXTUNW(1,65)='AMP                     '
      TXTUNW(1,66)='AMP                     '
      TXTUNW(1,67)='AMP                     '
      TXTUNW(1,68)='AMP                     '
      TXTUNW(1,69)='AMP                     '
      TXTUNW(1,70)='AMP                     '
      TXTUNW(1,71)='AMP                     '
      TXTUNW(1,72)='AMP                     '
      TXTUNW(1,73)='AMP                     '
      TXTUNW(1,74)='AMP                     '
      TXTUNW(1,75)='AMP                     '
      TXTUNW(1,76)='AMP                     '
      TXTUNW(1,77)='AMP                     '
      TXTUNW(1,78)='AMP                     '
      TXTUNW(1,79)='AMP                     '
      TXTUNW(1,80)='AMP                     '
      TXTUNW(1,81)='AMP                     '
      TXTUNW(1,82)='TO BE READ              '    ! ADD. SURF. TALLY
      TXTUNW(1,83)='TO BE READ              '    ! ALG. SURF. TALLY
      TXTUNW(1,84)='AMP                     '    ! PUMPED FLUX
      DO J=1,NTALS
        TXTUNW(2:N2MX,J)=TXTUNW(1,J)
      END DO

      RETURN
      END SUBROUTINE EIRENE_SETTXT
C
      SUBROUTINE EIRENE_STTXT1
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
      USE EIRMOD_CLOGAU

      IMPLICIT NONE
      INTEGER :: IATM, IION, IPLS, IMOL, ISPZ, IPHOT,
     .           JATM, JION, JMOL, JPHOT, JPLS,
     .           N1,  N2,  N3,  N4,  N5,  N5V, N6,  N7,  N8,  N9,  N10,
     .           N11, N12, N13, N14, N15, N16, N17, N18, N19, N20, N21,
     .           N22, N23, N24, N25, N26, N27, N28, N29, N30, N31, N32,
     .           KK, ICO
      INTEGER :: IAD, EIRENE_INDIRECT_ADDRESS
      EXTERNAL :: EIRENE_TEXT_COMBINE, EIRENE_INDIRECT_ADDRESS
C
C
      NFSTVI(1)=NATMI
      NFSTVI(2)=NMOLI
      NFSTVI(3)=NIONI
      NFSTVI(4)=NPHOTI
      NFSTVI(5)=NATMI
      NFSTVI(6)=NMOLI
      NFSTVI(7)=NIONI
      NFSTVI(8)=NPHOTI
      NFSTVI(9)=1
      NFSTVI(10)=NATMI
      NFSTVI(11)=NMOLI
      NFSTVI(12)=NIONI
      NFSTVI(13)=NPHOTI
      NFSTVI(14)=NPLSI
      IF (NLSPCSCL_ATM) THEN
        NFSTVI(10)=NATM*NATMP
        NFSTVI(11)=NMOL*NATMP
        NFSTVI(12)=NION*NATMP
        NFSTVI(13)=NPHOT*NATMP
        NFSTVI(14)=NPLS*NATMP
        NEXTVI(10)=NATM
        NEXTVI(11)=NMOL
        NEXTVI(12)=NION
        NEXTVI(13)=NPHOT
        NEXTVI(14)=NPLS
      END IF
      NFSTVI(15)=1
      NFSTVI(16)=NATMI
      NFSTVI(17)=NMOLI
      NFSTVI(18)=NIONI
      NFSTVI(19)=NPHOTI
      NFSTVI(20)=NPLSI
      IF (NLSPCSCL_MOL) THEN
        NFSTVI(16)=NATM*NMOLP
        NFSTVI(17)=NMOL*NMOLP
        NFSTVI(18)=NION*NMOLP
        NFSTVI(19)=NPHOT*NMOLP
        NFSTVI(20)=NPLS*NMOLP
        NEXTVI(16)=NATM
        NEXTVI(17)=NMOL
        NEXTVI(18)=NION
        NEXTVI(19)=NPHOT
        NEXTVI(20)=NPLS
      END IF
      NFSTVI(21)=1
      NFSTVI(22)=NATMI
      NFSTVI(23)=NMOLI
      NFSTVI(24)=NIONI
      NFSTVI(25)=NPHOTI
      NFSTVI(26)=NPLSI
      IF (NLSPCSCL_ION) THEN
        NFSTVI(22)=NATM*NIONP
        NFSTVI(23)=NMOL*NIONP
        NFSTVI(24)=NION*NIONP
        NFSTVI(25)=NPHOT*NIONP
        NFSTVI(26)=NPLS*NIONP
        NEXTVI(22)=NATM
        NEXTVI(23)=NMOL
        NEXTVI(24)=NION
        NEXTVI(25)=NPHOT
        NEXTVI(26)=NPLS
      END IF
      NFSTVI(27)=1
      NFSTVI(28)=NATMI
      NFSTVI(29)=NMOLI
      NFSTVI(30)=NIONI
      NFSTVI(31)=NPHOTI
      NFSTVI(32)=NPLSI
      IF (NLSPCSCL_PHOT) THEN
        NFSTVI(28)=NATM*NPHOTP
        NFSTVI(29)=NMOL*NPHOTP
        NFSTVI(30)=NION*NPHOTP
        NFSTVI(31)=NPHOT*NPHOTP
        NFSTVI(32)=NPLS*NPHOTP
        NEXTVI(28)=NATM
        NEXTVI(29)=NMOL
        NEXTVI(30)=NION
        NEXTVI(31)=NPHOT
        NEXTVI(32)=NPLS
      END IF
      NFSTVI(33)=1
      NFSTVI(34)=1
      NFSTVI(35)=1
      NFSTVI(36)=1
      NFSTVI(37)=1
      NFSTVI(38)=NPLSI
      NFSTVI(39)=1
      NFSTVI(40)=1
      NFSTVI(41)=1
      NFSTVI(42)=1
      NFSTVI(43)=1
      NFSTVI(44)=NPLSI
      NFSTVI(45)=1
      NFSTVI(46)=1
      NFSTVI(47)=1
      NFSTVI(48)=1
      NFSTVI(49)=1
      NFSTVI(50)=NPLSI
      NFSTVI(51)=1
      NFSTVI(52)=1
      NFSTVI(53)=1
      NFSTVI(54)=1
      NFSTVI(55)=1
      NFSTVI(56)=NPLSI
C
      NFSTVI(NTALA)=NADVI   ! TALLY 57
      IF (NLEMIS) NFSTVI(NTALA)=NADVI+NADV_ADD
      NFSTVI(NTALC)=NCLVI   ! TALLY 58
      NFSTVI(NTALT)=NSNVI   ! TALLY 59
      NFSTVI(NTALM)=NCPVI   ! TALLY 60
C     NFSTVI(NTALB)=NBGVI   ! TALLY 61, BUT NBGVI IS DEFINED IN SUBR. XSECT...
      NFSTVI(NTALB)=0
      NFSTVI(NTALR)=NALVI   ! TALLY 62
      NFSTVI(63)=NATMI
      NFSTVI(64)=NMOLI
      NFSTVI(65)=NIONI
      NFSTVI(66)=NPHOTI
      NFSTVI(67)=NATMI
      NFSTVI(68)=NMOLI
      NFSTVI(69)=NIONI
      NFSTVI(70)=NPHOTI
      NFSTVI(71)=NATMI
      NFSTVI(72)=NMOLI
      NFSTVI(73)=NIONI
      NFSTVI(74)=NPHOTI
      NFSTVI(75)=NATMI
      NFSTVI(76)=NMOLI
      NFSTVI(77)=NIONI
      NFSTVI(78)=NPHOTI
      NFSTVI(79)=NPLSI
      NFSTVI(80)=1
      NFSTVI(81)=1
      NFSTVI(82)=1
      NFSTVI(83)=1
      NFSTVI(84)=NPLSI
C  MOMENTUM DENSITY TALLIES, X,Y,Z
      NFSTVI(85)=NATMI
      NFSTVI(86)=NMOLI
      NFSTVI(87)=NIONI
      NFSTVI(88)=NPHOTI
      NFSTVI(89)=NATMI
      NFSTVI(90)=NMOLI
      NFSTVI(91)=NIONI
      NFSTVI(92)=NPHOTI
      NFSTVI(93)=NATMI
      NFSTVI(94)=NMOLI
      NFSTVI(95)=NIONI
      NFSTVI(96)=NPHOTI
cdr vectorial momentum sources, for plasma species ipls, 1,...,nplsi
      NFSTVI(97)=3*NPLSI
      NFSTVI(98)=3*NPLSI
      NFSTVI(99)=3*NPLSI
      NFSTVI(100)=3*NPLSI

      NFSTVI(101)=NATMI
      NFSTVI(102)=NMOLI
      NFSTVI(103)=NIONI
C
C
      NFSTWI(1)=NATMI
      NFSTWI(2)=NATMI
      IF (NLSPCSCL_ATM) THEN
        NFSTWI(2)=NATM*NATMP
        NEXTWI(2)=NATM
      END IF
      NFSTWI(3)=NATMI
      IF (NLSPCSCL_MOL) THEN
        NFSTWI(3)=NATM*NMOLP
        NEXTWI(3)=NATM
      END IF
      NFSTWI(4)=NATMI
      IF (NLSPCSCL_ION) THEN
        NFSTWI(4)=NATM*NIONP
        NEXTWI(4)=NATM
      END IF
      NFSTWI(5)=NATMI
      IF (NLSPCSCL_PHOT) THEN
        NFSTWI(5)=NATM*NPHOTP
        NEXTWI(5)=NATM
      END IF
      NFSTWI(6)=NATMI
      NFSTWI(7)=NMOLI
      NFSTWI(8)=NMOLI
      IF (NLSPCSCL_ATM) THEN
        NFSTWI(8)=NMOL*NATMP
        NEXTWI(8)=NMOL
      END IF
      NFSTWI(9)=NMOLI
      IF (NLSPCSCL_MOL) THEN
        NFSTWI(9)=NMOL*NMOLP
        NEXTWI(9)=NMOL
      END IF
      NFSTWI(10)=NMOLI
      IF (NLSPCSCL_ION) THEN
        NFSTWI(10)=NMOL*NIONP
        NEXTWI(10)=NMOL
      END IF
      NFSTWI(11)=NMOLI
      IF (NLSPCSCL_PHOT) THEN
        NFSTWI(11)=NMOL*NPHOTP
        NEXTWI(11)=NMOL
      END IF
      NFSTWI(12)=NMOLI
      NFSTWI(13)=NIONI
      NFSTWI(14)=NIONI
      IF (NLSPCSCL_ATM) THEN
        NFSTWI(14)=NION*NATMP
        NEXTWI(14)=NION
      END IF
      NFSTWI(15)=NIONI
      IF (NLSPCSCL_MOL) THEN
        NFSTWI(15)=NION*NMOLP
        NEXTWI(15)=NION
      END IF
      NFSTWI(16)=NIONI
      IF (NLSPCSCL_ION) THEN
        NFSTWI(16)=NION*NIONP
        NEXTWI(16)=NION
      END IF
      NFSTWI(17)=NIONI
      IF (NLSPCSCL_PHOT) THEN
        NFSTWI(17)=NION*NPHOTP
        NEXTWI(17)=NION
      END IF
      NFSTWI(18)=NIONI
      NFSTWI(19)=NPHOTI
      NFSTWI(20)=NPHOTI
      IF (NLSPCSCL_ATM) THEN
        NFSTWI(20)=NPHOT*NATMP
        NEXTWI(20)=NPHOT
      END IF
      NFSTWI(21)=NPHOTI
      IF (NLSPCSCL_MOL) THEN
        NFSTWI(21)=NPHOT*NMOLP
        NEXTWI(21)=NPHOT
      END IF
      NFSTWI(22)=NPHOTI
      IF (NLSPCSCL_ION) THEN
        NFSTWI(22)=NPHOT*NIONP
        NEXTWI(22)=NPHOT
      END IF
      NFSTWI(23)=NPHOTI
      IF (NLSPCSCL_PHOT) THEN
        NFSTWI(23)=NPHOT*NPHOTP
        NEXTWI(23)=NPHOT
      END IF
      NFSTWI(24)=NPHOTI
      NFSTWI(25)=NPLSI

      NFSTWI(26)=NATMI
      NFSTWI(27)=NATMI
      NFSTWI(28)=NATMI
      NFSTWI(29)=NATMI
      NFSTWI(30)=NATMI
      NFSTWI(31)=NATMI
      NFSTWI(32)=NMOLI
      NFSTWI(33)=NMOLI
      NFSTWI(34)=NMOLI
      NFSTWI(35)=NMOLI
      NFSTWI(36)=NMOLI
      NFSTWI(37)=NMOLI
      NFSTWI(38)=NIONI
      NFSTWI(39)=NIONI
      NFSTWI(40)=NIONI
      NFSTWI(41)=NIONI
      NFSTWI(42)=NIONI
      NFSTWI(43)=NIONI
      NFSTWI(44)=NPHOTI
      NFSTWI(45)=NPHOTI
      NFSTWI(46)=NPHOTI
      NFSTWI(47)=NPHOTI
      NFSTWI(48)=NPHOTI
      NFSTWI(49)=NPHOTI
      NFSTWI(50)=NPLSI

      NFSTWI(51)=NATMI
      IF (NLSPCSPT_ATM) THEN
        NFSTWI(51)=NATM*NATMP
        NEXTWI(51)=NATM
      END IF
      NFSTWI(52)=NATMI
      IF (NLSPCSPT_MOL) THEN
        NFSTWI(52)=NATM*NMOLP
        NEXTWI(52)=NATM
      END IF
      NFSTWI(53)=NATMI
      IF (NLSPCSPT_ION) THEN
        NFSTWI(53)=NATM*NIONP
        NEXTWI(53)=NATM
      END IF
      NFSTWI(54)=NATMI
      IF (NLSPCSPT_PHOT) THEN
        NFSTWI(54)=NATM*NPHOTP
        NEXTWI(54)=NATM
      END IF
      NFSTWI(55)=NATMI
      IF (NLSPCSPT_PLS) THEN
        NFSTWI(55)=NATM*NPLSP
        NEXTWI(55)=NATM
      END IF
      NFSTWI(56)=NMOLI
      IF (NLSPCSPT_ATM) THEN
        NFSTWI(56)=NMOL*NATMP
        NEXTWI(56)=NMOL
      END IF
      NFSTWI(57)=NMOLI
      IF (NLSPCSPT_MOL) THEN
        NFSTWI(57)=NMOL*NMOLP
        NEXTWI(57)=NMOL
      END IF
      NFSTWI(58)=NMOLI
      IF (NLSPCSPT_ION) THEN
        NFSTWI(58)=NMOL*NIONP
        NEXTWI(58)=NMOL
      END IF
      NFSTWI(59)=NMOLI
      IF (NLSPCSPT_PHOT) THEN
        NFSTWI(59)=NMOL*NPHOTP
        NEXTWI(59)=NMOL
      END IF
      NFSTWI(60)=NMOLI
      IF (NLSPCSPT_PLS) THEN
        NFSTWI(60)=NMOL*NPLSP
        NEXTWI(60)=NMOL
      END IF
      NFSTWI(61)=NIONI
      IF (NLSPCSPT_ATM) THEN
        NFSTWI(61)=NION*NATMP
        NEXTWI(61)=NION
      END IF
      NFSTWI(62)=NIONI
      IF (NLSPCSPT_MOL) THEN
        NFSTWI(62)=NION*NMOLP
        NEXTWI(62)=NION
      END IF
      NFSTWI(63)=NIONI
      IF (NLSPCSPT_ION) THEN
        NFSTWI(63)=NION*NIONP
        NEXTWI(63)=NION
      END IF
      NFSTWI(64)=NIONI
      IF (NLSPCSPT_PHOT) THEN
        NFSTWI(64)=NION*NPHOTP
        NEXTWI(64)=NION
      END IF
      NFSTWI(65)=NIONI
      IF (NLSPCSPT_PLS) THEN
        NFSTWI(65)=NION*NPLSP
        NEXTWI(65)=NION
      END IF
      NFSTWI(66)=NPHOTI
      IF (NLSPCSPT_ATM) THEN
        NFSTWI(66)=NPHOT*NATMP
        NEXTWI(66)=NPHOT
      END IF
      NFSTWI(67)=NPHOTI
      IF (NLSPCSPT_MOL) THEN
        NFSTWI(67)=NPHOT*NMOLP
        NEXTWI(67)=NPHOT
      END IF
      NFSTWI(68)=NPHOTI
      IF (NLSPCSPT_ION) THEN
        NFSTWI(68)=NPHOT*NIONP
        NEXTWI(68)=NPHOT
      END IF
      NFSTWI(69)=NPHOTI
      IF (NLSPCSPT_PHOT) THEN
        NFSTWI(69)=NPHOT*NPHOTP
        NEXTWI(69)=NPHOT
      END IF
      NFSTWI(70)=NPHOTI
      IF (NLSPCSPT_PLS) THEN
        NFSTWI(70)=NPHOT*NPLSP
        NEXTWI(70)=NPHOT
      END IF
      NFSTWI(71)=NPLSI
      IF (NLSPCSPT_ATM) THEN
        NFSTWI(71)=NPLS*NATMP
        NEXTWI(71)=NPLS
      END IF
      NFSTWI(72)=NPLSI
      IF (NLSPCSPT_MOL) THEN
        NFSTWI(72)=NPLS*NMOLP
        NEXTWI(72)=NPLS
      END IF
      NFSTWI(73)=NPLSI
      IF (NLSPCSPT_ION) THEN
        NFSTWI(73)=NPLS*NIONP
        NEXTWI(73)=NPLS
      END IF
      NFSTWI(74)=NPLSI
      IF (NLSPCSPT_PHOT) THEN
        NFSTWI(74)=NPLS*NPHOTP
        NEXTWI(74)=NPLS
      END IF
      NFSTWI(75)=NPLSI
      IF (NLSPCSPT_PLS) THEN
        NFSTWI(75)=NPLS*NPLSP
        NEXTWI(75)=NPLS
      END IF
      NFSTWI(76)=1
      NFSTWI(77)=1
      NFSTWI(78)=1
      NFSTWI(79)=1
      NFSTWI(80)=1
      NFSTWI(81)=1
      NFSTWI(NTLSA)=NADSI    !  ADD SURF. TALLY
      NFSTWI(NTLSR)=NALSI    !  ALG. SURF. TALLY
      NFSTWI(NTALS)=NSPTOT   !  PUMPED FLUX
C
C  DEFINE INITIAL AND LAST "SPECIES INDEX" FOR VOLUME TALLIES, NSPAN, NSPEN
C  NEEDED FOR STANDARD DEVIATION MARKER LMETSP(..), which is reset for each history

      N1=NPHOTI
      N2=N1+NATMI
      N3=N2+NMOLI
      N4=N3+NIONI
      N5=N4+NPLSI
      N5V =N4+3*NPLSI
c  additional tallies
      N6=N5+NADVI
      IF (NLEMIS) N6 = N5+NADVI+NADV_ADD
      N7=N6+NALVI
      N8=N7+NCLVI
      N9=N8+NCPVI
      N10=N9+NBGVI
      N11=N10+NSNVI

c  for species specific rescaling
C  tallies resolved for atoms
      N12=N11+NATM*NATMP
      N13=N12+NMOL*NATMP
      N14=N13+NION*NATMP
      N15=N14+NPHOT*NATMP
      N16=N15+NPLS*NATMP
C  tallies resolved for molecules
      N17=N16+NATM*NMOLP
      N18=N17+NMOL*NMOLP
      N19=N18+NION*NMOLP
      N20=N19+NPHOT*NMOLP
      N21=N20+NPLS*NMOLP
C  tallies resolved for test ions
      N22=N21+NATM*NIONP
      N23=N22+NMOL*NIONP
      N24=N23+NION*NIONP
      N25=N24+NPHOT*NIONP
      N26=N25+NPLS*NIONP
C  tallies resolved for photons
      N27=N26+NATM*NPHOTP
      N28=N27+NMOL*NPHOTP
      N29=N28+NION*NPHOTP
      N30=N29+NPHOT*NPHOTP
      N31=N30+NPLS*NPHOTP

      NSPAN(1)=N1+1
      NSPAN(2)=N2+1
      NSPAN(3)=N3+1
      NSPAN(4)=1

      NSPAN(5)=N1+1
      NSPAN(6)=N2+1
      NSPAN(7)=N3+1
      NSPAN(8)=1

      NSPAN(9)=0
      IF (NLSPCSCL_ATM) THEN
        NSPAN(10)=N11+1
        NSPAN(11)=N12+1
        NSPAN(12)=N13+1
        NSPAN(13)=N14+1
        NSPAN(14)=N15+1
      ELSE
        NSPAN(10)=N1+1
        NSPAN(11)=N2+1
        NSPAN(12)=N3+1
        NSPAN(13)=1
        NSPAN(14)=N4+1
      END IF
      NSPAN(15)=0
      IF (NLSPCSCL_MOL) THEN
        NSPAN(16)=N16+1
        NSPAN(17)=N17+1
        NSPAN(18)=N18+1
        NSPAN(19)=N19+1
        NSPAN(20)=N20+1
      ELSE
        NSPAN(16)=N1+1
        NSPAN(17)=N2+1
        NSPAN(18)=N3+1
        NSPAN(19)=1
        NSPAN(20)=N4+1
      END IF
      NSPAN(21)=0
      IF (NLSPCSCL_ION) THEN
        NSPAN(22)=N21+1
        NSPAN(23)=N22+1
        NSPAN(24)=N23+1
        NSPAN(25)=N24+1
        NSPAN(26)=N25+1
      ELSE
        NSPAN(22)=N1+1
        NSPAN(23)=N2+1
        NSPAN(24)=N2+1
        NSPAN(25)=1
        NSPAN(26)=N4+1
      END IF
      NSPAN(27)=0
      IF (NLSPCSCL_PHOT) THEN
        NSPAN(28)=N26+1
        NSPAN(29)=N27+1
        NSPAN(30)=N28+1
        NSPAN(31)=N29+1
        NSPAN(32)=N30+1
      ELSE
        NSPAN(28)=N1+1
        NSPAN(29)=N2+1
        NSPAN(30)=N2+1
        NSPAN(31)=1
        NSPAN(32)=N4+1
      END IF
      NSPAN(33)=0
      NSPAN(34)=0
      NSPAN(35)=0
      NSPAN(36)=0
      NSPAN(37)=0
      NSPAN(38)=N4+1
      NSPAN(39)=0
      NSPAN(40)=0
      NSPAN(41)=0
      NSPAN(42)=0
      NSPAN(43)=0
      NSPAN(44)=N4+1
      NSPAN(45)=0
      NSPAN(46)=0
      NSPAN(47)=0
      NSPAN(48)=0
      NSPAN(49)=0
      NSPAN(50)=N4+1
      NSPAN(51)=0
      NSPAN(52)=0
      NSPAN(53)=0
      NSPAN(54)=0
      NSPAN(55)=0
      NSPAN(56)=N4+1
C  ADDITIONAL TALLIES

c  additional tracklength estimators
      NSPAN(NTALA)=N5+1
c  additional collision estimators
      NSPAN(NTALC)=N7+1
c  additional snapshot estimators
      NSPAN(NTALT)=N10+1
c  additional couple tallies
      NSPAN(NTALM)=N8+1
c  additional bgk tallies
      NSPAN(NTALB)=N9+1
c  additional algebraic tallies
      NSPAN(NTALR)=N6+1

C  GENERATION LIMIT TALLIES
      NSPAN(63)=N1+1
      NSPAN(64)=N2+1
      NSPAN(65)=N3+1
c
      NSPAN(66)=1
      NSPAN(67)=N1+1
      NSPAN(68)=N2+1
      NSPAN(69)=N3+1
c
      NSPAN(70)=1
      NSPAN(71)=N1+1
      NSPAN(72)=N2+1
      NSPAN(73)=N3+1
c
      NSPAN(74)=1
      NSPAN(75)=N1+1
      NSPAN(76)=N2+1
      NSPAN(77)=N3+1
c
      NSPAN(78)=1
      NSPAN(79)=N4+1
      NSPAN(80)=0
      NSPAN(81)=0
      NSPAN(82)=0
      NSPAN(83)=0
      NSPAN(84)=N4+1

      NSPAN(85)=N1+1
      NSPAN(86)=N2+1
      NSPAN(87)=N3+1
      NSPAN(88)=1
      NSPAN(89)=N1+1
      NSPAN(90)=N2+1
      NSPAN(91)=N3+1
      NSPAN(92)=1
      NSPAN(93)=N1+1
      NSPAN(94)=N2+1
      NSPAN(95)=N3+1
      NSPAN(96)=1
cdr VECTORIAL MOMENTUM SOURCE TALLIES
      NSPAN(97)=N4+1
      NSPAN(98)=N4+1
      NSPAN(99)=N4+1
      NSPAN(100)=N4+1

      NSPAN(101)=N1+1
      NSPAN(102)=N2+1
      NSPAN(103)=N3+1

      NSPEN(1)=N2
      NSPEN(2)=N3
      NSPEN(3)=N4
      NSPEN(4)=N1

      NSPEN(5)=N2
      NSPEN(6)=N3
      NSPEN(7)=N4
      NSPEN(8)=N1

      NSPEN(9)=0
      IF (NLSPCSCL_ATM) THEN
        NSPEN(10)=N12
        NSPEN(11)=N13
        NSPEN(12)=N14
        NSPEN(13)=N15
        NSPEN(14)=N16
      ELSE
        NSPEN(10)=N2
        NSPEN(11)=N3
        NSPEN(12)=N4
        NSPEN(13)=N1
        NSPEN(14)=N5
      ENDIF
      NSPEN(15)=0
      IF (NLSPCSCL_MOL) THEN
        NSPEN(16)=N17
        NSPEN(17)=N18
        NSPEN(18)=N19
        NSPEN(19)=N20
        NSPEN(20)=N21
      ELSE
        NSPEN(16)=N2
        NSPEN(17)=N3
        NSPEN(18)=N4
        NSPEN(19)=N1
        NSPEN(20)=N5
      ENDIF
      NSPEN(21)=0
      IF (NLSPCSCL_ION) THEN
        NSPEN(22)=N22
        NSPEN(23)=N23
        NSPEN(24)=N24
        NSPEN(25)=N25
        NSPEN(26)=N26
      ELSE
        NSPEN(22)=N2
        NSPEN(23)=N3
        NSPEN(24)=N4
        NSPEN(25)=N1
        NSPEN(26)=N5
      ENDIF
      NSPEN(27)=0
      IF (NLSPCSCL_PHOT) THEN
        NSPEN(28)=N27
        NSPEN(29)=N28
        NSPEN(30)=N29
        NSPEN(31)=N30
        NSPEN(32)=N31
      ELSE
        NSPEN(28)=N2
        NSPEN(29)=N3
        NSPEN(30)=N4
        NSPEN(31)=N1
        NSPEN(32)=N5
      ENDIF
      NSPEN(33)=0
      NSPEN(34)=0
      NSPEN(35)=0
      NSPEN(36)=0
      NSPEN(37)=0
      NSPEN(38)=N5
      NSPEN(39)=0
      NSPEN(40)=0
      NSPEN(41)=0
      NSPEN(42)=0
      NSPEN(43)=0
      NSPEN(44)=N5
      NSPEN(45)=0
      NSPEN(46)=0
      NSPEN(47)=0
      NSPEN(48)=0
      NSPEN(49)=0
      NSPEN(50)=N5
      NSPEN(51)=0
      NSPEN(52)=0
      NSPEN(53)=0
      NSPEN(54)=0
      NSPEN(55)=0
      NSPEN(56)=N5
C  ADDITIONAL TALLIES
      NSPEN(NTALA)=N6
      NSPEN(NTALC)=N8
      NSPEN(NTALT)=N11
      NSPEN(NTALM)=N9
      NSPEN(NTALB)=N10
      NSPEN(NTALR)=N7
C  GENERATION LIMIT TALLIES
      NSPEN(63)=N2
      NSPEN(64)=N3
      NSPEN(65)=N4
      NSPEN(66)=N1
      NSPEN(67)=N2
      NSPEN(68)=N3
      NSPEN(69)=N4
      NSPEN(70)=N1
      NSPEN(71)=N2
      NSPEN(72)=N3
      NSPEN(73)=N4
      NSPEN(74)=N1
      NSPEN(75)=N2
      NSPEN(76)=N3
      NSPEN(77)=N4
      NSPEN(78)=N1

      NSPEN(79)=N5
      NSPEN(80)=0
      NSPEN(81)=0
      NSPEN(82)=0
      NSPEN(83)=0

      NSPEN(84)=N5

      NSPEN(85)=N2
      NSPEN(86)=N3
      NSPEN(87)=N4
      NSPEN(88)=N1

      NSPEN(89)=N2
      NSPEN(90)=N3
      NSPEN(91)=N4
      NSPEN(92)=N1

      NSPEN(93)=N2
      NSPEN(94)=N3
      NSPEN(95)=N4
      NSPEN(96)=N1
cdr VECTORIAL MOMENTUM SOURCE TALLIES
      NSPEN(97)=N5V
      NSPEN(98)=N5V
      NSPEN(99)=N5V
      NSPEN(100)=N5V

      NSPEN(101)=N2
      NSPEN(102)=N3
      NSPEN(103)=N4

      DO IPHOT=1,NPHOTI
        ISPZ=IPHOT
        TXTSPC(IPHOT,4)=TEXTS(ISPZ)
        TXTSPC(IPHOT,8)=TEXTS(ISPZ)
        TXTSPC(IPHOT,13)=TEXTS(ISPZ)
        IF (NLSPCSCL_ATM) THEN
          DO IATM=1, NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IATM,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPC(IAD,13))
          END DO
        END IF
        TXTSPC(IPHOT,19)=TEXTS(ISPZ)
        IF (NLSPCSCL_MOL) THEN
          DO IMOL=1, NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IMOL,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPC(IAD,19))
          END DO
        END IF
        TXTSPC(IPHOT,25)=TEXTS(ISPZ)
        IF (NLSPCSCL_ION) THEN
          DO IION=1, NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IION,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPC(IAD,25))
          END DO
        END IF
        TXTSPC(IPHOT,31)=TEXTS(ISPZ)
        IF (NLSPCSCL_PHOT) THEN
          DO JPHOT=1, NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,JPHOT,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(JPHOT),
     .                               TXTSPC(IAD,31))
          END DO
        END IF
        TXTSPC(IPHOT,66)=TEXTS(ISPZ)
        TXTSPC(IPHOT,70)=TEXTS(ISPZ)
        TXTSPC(IPHOT,74)=TEXTS(ISPZ)
        TXTSPC(IPHOT,78)=TEXTS(ISPZ)
        TXTSPC(IPHOT,88)=TEXTS(ISPZ)
        TXTSPC(IPHOT,92)=TEXTS(ISPZ)
        TXTSPC(IPHOT,96)=TEXTS(ISPZ)
      END DO

      DO 10 IATM=1,NATMI
        ISPZ=NSPH+IATM
        TXTSPC(IATM,1)=TEXTS(ISPZ)
        TXTSPC(IATM,5)=TEXTS(ISPZ)
        TXTSPC(IATM,10)=TEXTS(ISPZ)
        IF (NLSPCSCL_ATM) THEN
          DO JATM=1, NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,JATM,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+JATM),
     .                               TXTSPC(IAD,10))
          END DO
        END IF
        TXTSPC(IATM,16)=TEXTS(ISPZ)
        IF (NLSPCSCL_MOL) THEN
          DO IMOL=1, NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IMOL,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPC(IAD,16))
          END DO
        END IF
        TXTSPC(IATM,22)=TEXTS(ISPZ)
        IF (NLSPCSCL_ION) THEN
          DO IION=1, NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IION,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPC(IAD,22))
          END DO
        END IF
        TXTSPC(IATM,28)=TEXTS(ISPZ)
        IF (NLSPCSCL_PHOT) THEN
          DO IPHOT=1, NPHOT
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IPHOT,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPC(IAD,28))
          END DO
        END IF
        TXTSPC(IATM,63)=TEXTS(ISPZ)
        TXTSPC(IATM,67)=TEXTS(ISPZ)
        TXTSPC(IATM,71)=TEXTS(ISPZ)
        TXTSPC(IATM,75)=TEXTS(ISPZ)
        TXTSPC(IATM,85)=TEXTS(ISPZ)
        TXTSPC(IATM,89)=TEXTS(ISPZ)
        TXTSPC(IATM,93)=TEXTS(ISPZ)
        TXTSPC(IATM,101)=TEXTS(ISPZ)
   10 CONTINUE
C
      DO 20 IMOL=1,NMOLI
        ISPZ=NSPA+IMOL
        TXTSPC(IMOL,2)=TEXTS(ISPZ)
        TXTSPC(IMOL,6)=TEXTS(ISPZ)
        TXTSPC(IMOL,11)=TEXTS(ISPZ)
        IF (NLSPCSCL_ATM) THEN
          DO IATM=1, NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IATM,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPC(IAD,11))
          END DO
        END IF
        TXTSPC(IMOL,17)=TEXTS(ISPZ)
        IF (NLSPCSCL_MOL) THEN
          DO JMOL=1, NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,JMOL,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+JMOL),
     .                               TXTSPC(IAD,17))
          END DO
        END IF
        TXTSPC(IMOL,23)=TEXTS(ISPZ)
        IF (NLSPCSCL_ION) THEN
          DO IION=1, NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IION,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPC(IAD,23))
          END DO
        END IF
        TXTSPC(IMOL,29)=TEXTS(ISPZ)
        IF (NLSPCSCL_PHOT) THEN
          DO IPHOT=1, NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IPHOT,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPC(IAD,29))
          END DO
        END IF
        TXTSPC(IMOL,64)=TEXTS(ISPZ)
        TXTSPC(IMOL,68)=TEXTS(ISPZ)
        TXTSPC(IMOL,72)=TEXTS(ISPZ)
        TXTSPC(IMOL,76)=TEXTS(ISPZ)
        TXTSPC(IMOL,86)=TEXTS(ISPZ)
        TXTSPC(IMOL,90)=TEXTS(ISPZ)
        TXTSPC(IMOL,94)=TEXTS(ISPZ)
        TXTSPC(IMOL,102)=TEXTS(ISPZ)
   20 CONTINUE
C
      DO 30 IION=1,NIONI
        ISPZ=NSPAM+IION
        TXTSPC(IION,3)=TEXTS(ISPZ)
        TXTSPC(IION,7)=TEXTS(ISPZ)
        TXTSPC(IION,12)=TEXTS(ISPZ)
        IF (NLSPCSCL_ATM) THEN
          DO IATM=1, NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IATM,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPC(IAD,12))
          END DO
        END IF
        TXTSPC(IION,18)=TEXTS(ISPZ)
        IF (NLSPCSCL_MOL) THEN
          DO IMOL=1, NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IMOL,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPC(IAD,18))
          END DO
        END IF
        TXTSPC(IION,24)=TEXTS(ISPZ)
        IF (NLSPCSCL_ION) THEN
          DO JION=1, NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,JION,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IMOL),
     .                               TXTSPC(IAD,24))
          END DO
        END IF
        TXTSPC(IION,30)=TEXTS(ISPZ)
        IF (NLSPCSCL_PHOT) THEN
          DO IPHOT=1, NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IPHOT,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPC(IAD,30))
          END DO
        END IF
        TXTSPC(IION,65)=TEXTS(ISPZ)
        TXTSPC(IION,69)=TEXTS(ISPZ)
        TXTSPC(IION,73)=TEXTS(ISPZ)
        TXTSPC(IION,77)=TEXTS(ISPZ)
        TXTSPC(IION,87)=TEXTS(ISPZ)
        TXTSPC(IION,91)=TEXTS(ISPZ)
        TXTSPC(IION,95)=TEXTS(ISPZ)
        TXTSPC(IION,103)=TEXTS(ISPZ)
   30 CONTINUE
C
      DO 40 IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTSPC(IPLS,14)=TEXTS(ISPZ)
        IF (NLSPCSCL_ATM) THEN
          DO IATM=1, NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IATM,NPLS)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPC(IAD,14))
          END DO
        END IF
        TXTSPC(IPLS,20)=TEXTS(ISPZ)
        IF (NLSPCSCL_MOL) THEN
          DO IMOL=1, NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IMOL,NPLS)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPC(IAD,20))
          END DO
        END IF
        TXTSPC(IPLS,26)=TEXTS(ISPZ)
        IF (NLSPCSCL_ION) THEN
          DO IION=1, NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IION,NPLS)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPC(IAD,26))
          END DO
        END IF
        TXTSPC(IPLS,32)=TEXTS(ISPZ)
        IF (NLSPCSCL_PHOT) THEN
          DO IPHOT=1, NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IPHOT,NPLS)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPC(IAD,32))
          END DO
        END IF
        TXTSPC(IPLS,38)=TEXTS(ISPZ)
        TXTSPC(IPLS,44)=TEXTS(ISPZ)
        TXTSPC(IPLS,50)=TEXTS(ISPZ)
        TXTSPC(IPLS,56)=TEXTS(ISPZ)
        TXTSPC(IPLS,79)=TEXTS(ISPZ)
        TXTSPC(IPLS,84)=TEXTS(ISPZ)
        DO ICO=0,2
          KK=ICO*NPLS
          TXTSPC(KK+IPLS,97)=TEXTS(ISPZ)
          TXTSPC(KK+IPLS,98)=TEXTS(ISPZ)
          TXTSPC(KK+IPLS,99)=TEXTS(ISPZ)
          TXTSPC(KK+IPLS,100)=TEXTS(ISPZ)
        ENDDO
   40 CONTINUE
C
      TXTSPC(1,9)='ELECTRONS               '
      TXTSPC(1,15)='ELECTRONS               '
      TXTSPC(1,21)='ELECTRONS               '
      TXTSPC(1,27)='ELECTRONS               '
      TXTSPC(1,33)='ELECTRONS               '
      TXTSPC(1,39)='ELECTRONS               '
      TXTSPC(1,45)='ELECTRONS               '
      TXTSPC(1,51)='ELECTRONS               '
C
      TXTSPC(1,34)='ATOMS                   '
      TXTSPC(1,40)='ATOMS                   '
      TXTSPC(1,46)='ATOMS                   '
      TXTSPC(1,52)='ATOMS                   '
      TXTSPC(1,80)='ATOMS                   '
C
      TXTSPC(1,35)='MOLECULES               '
      TXTSPC(1,41)='MOLECULES               '
      TXTSPC(1,47)='MOLECULES               '
      TXTSPC(1,53)='MOLECULES               '
      TXTSPC(1,81)='MOLECULES               '
C
      TXTSPC(1,36)='TEST IONS               '
      TXTSPC(1,42)='TEST IONS               '
      TXTSPC(1,48)='TEST IONS               '
      TXTSPC(1,54)='TEST IONS               '
      TXTSPC(1,82)='TEST IONS               '
C
      TXTSPC(1,37)='PHOTONS                 '
      TXTSPC(1,43)='PHOTONS                 '
      TXTSPC(1,49)='PHOTONS                 '
      TXTSPC(1,55)='PHOTONS                 '
      TXTSPC(1,83)='PHOTONS                 '
C
C
C
C  INITIALISE SPECIES ARRAYS FOR SURFACE TALLIES
c  n1 -- n5: same as for volume tallies
c     N1=NPHOTI
c     N2=N1+NATMI
c     N3=N2+NMOLI
c     N4=N3+NIONI
c     N5=N4+NPLSI   = NSPZTOT
cdr also pumped flux SPUMP: now 1:N5  (was: n7+1:n8)
c
      N6=N5+NADSI
      N7=N6+NALSI
C  for species specified scoring
C  for fluxes going into atoms
      N8=N7+NATM*NATMP
      N9=N8+NATM*NMOLP
      N10=N9+NATM*NIONP
      N11=N10+NATM*NPHOTP
C  for fluxes going into molecules
      N12=N11+NMOL*NATMP
      N13=N12+NMOL*NMOLP
      N14=N13+NMOL*NIONP
      N15=N14+NMOL*NPHOTP
C  for fluxes going into test ions
      N16=N15+NION*NATMP
      N17=N16+NION*NMOLP
      N18=N17+NION*NIONP
      N19=N18+NION*NPHOTP
C  for fluxes going into photons
      N20=N19+NPHOT*NATMP
      N21=N20+NPHOT*NMOLP
      N22=N21+NPHOT*NIONP
      N23=N22+NPHOT*NPHOTP

c  extra dimension for bulk ion resolved sputtering
c  for fluxes going into atoms
      N24=NATM*NPLSP
c  for fluxes going into molecules
      N25=N24+NMOL*NPLSP
c  for fluxes going into test ions
      N26=N25+NION*NPLSP
c  for fluxes going into photons
      N27=N23+N26+NPHOT*NPLSP
c  for fluxes going into bulk ions
      N28=N27+NPLS*NATMP
      N29=N28+NPLS*NMOLP
      N30=N29+NPLS*NIONP
      N31=N30+NPLS*NPHOTP
      N32=N31+NPLS*NPLSP

      NSPANW(1)=N1+1
      NSPANW(2)=N1+1
      IF (NLSPCSCL_ATM) NSPANW(2)=N7+1
      NSPANW(3)=N1+1
      IF (NLSPCSCL_MOL) NSPANW(3)=N8+1
      NSPANW(4)=N1+1
      IF (NLSPCSCL_ION) NSPANW(4)=N9+1
      NSPANW(5)=N1+1
      IF (NLSPCSCL_PHOT) NSPANW(5)=N10+1
      NSPANW(6)=N1+1
      NSPANW(7)=N2+1
      NSPANW(8)=N2+1
      IF (NLSPCSCL_ATM) NSPANW(8)=N11+1
      NSPANW(9)=N2+1
      IF (NLSPCSCL_MOL) NSPANW(9)=N12+1
      NSPANW(10)=N2+1
      IF (NLSPCSCL_ION) NSPANW(10)=N13+1
      NSPANW(11)=N2+1
      IF (NLSPCSCL_PHOT) NSPANW(11)=N14+1
      NSPANW(12)=N2+1
      NSPANW(13)=N3+1
      NSPANW(14)=N3+1
      IF (NLSPCSCL_ATM) NSPANW(14)=N15+1
      NSPANW(15)=N3+1
      IF (NLSPCSCL_MOL) NSPANW(15)=N16+1
      NSPANW(16)=N3+1
      IF (NLSPCSCL_ION) NSPANW(16)=N17+1
      NSPANW(17)=N3+1
      IF (NLSPCSCL_PHOT) NSPANW(17)=N18+1
      NSPANW(18)=N3+1
      NSPANW(19)=1
      NSPANW(20)=1
      IF (NLSPCSCL_ATM) NSPANW(20)=N19+1
      NSPANW(21)=1
      IF (NLSPCSCL_MOL) NSPANW(21)=N20+1
      NSPANW(22)=1
      IF (NLSPCSCL_ION) NSPANW(22)=N21+1
      NSPANW(23)=1
      IF (NLSPCSCL_PHOT) NSPANW(23)=N22+1
      NSPANW(24)=1
      NSPANW(25)=N4+1
      NSPANW(26)=N1+1
      NSPANW(27)=N1+1
      NSPANW(28)=N1+1
      NSPANW(29)=N1+1
      NSPANW(30)=N1+1
      NSPANW(31)=N1+1
      NSPANW(32)=N2+1
      NSPANW(33)=N2+1
      NSPANW(34)=N2+1
      NSPANW(35)=N2+1
      NSPANW(36)=N2+1
      NSPANW(37)=N2+1
      NSPANW(38)=N3+1
      NSPANW(39)=N3+1
      NSPANW(40)=N3+1
      NSPANW(41)=N3+1
      NSPANW(42)=N3+1
      NSPANW(43)=N3+1
      NSPANW(44)=1
      NSPANW(45)=1
      NSPANW(46)=1
      NSPANW(47)=1
      NSPANW(48)=1
      NSPANW(49)=1
      NSPANW(50)=N4+1
      NSPANW(51)=N1+1
      IF (NLSPCSPT_ATM) NSPANW(51)=N7+1
      NSPANW(52)=N1+1
      IF (NLSPCSPT_MOL) NSPANW(52)=N8+1
      NSPANW(53)=N1+1
      IF (NLSPCSPT_ION) NSPANW(53)=N9+1
      NSPANW(54)=N1+1
      IF (NLSPCSPT_PHOT) NSPANW(54)=N10+1
      NSPANW(55)=N1+1
      IF (NLSPCSPT_PLS) NSPANW(55)=N11+1
      NSPANW(56)=N2+1
      IF (NLSPCSPT_ATM) NSPANW(56)=N11+N24+1
      NSPANW(57)=N2+1
      IF (NLSPCSPT_MOL) NSPANW(57)=N12+N24+1
      NSPANW(58)=N2+1
      IF (NLSPCSPT_ION) NSPANW(58)=N13+N24+1
      NSPANW(59)=N2+1
      IF (NLSPCSPT_PHOT) NSPANW(59)=N14+N24+1
      NSPANW(60)=N2+1
      IF (NLSPCSPT_PLS) NSPANW(60)=N15+N24+1
      NSPANW(61)=N3+1
      IF (NLSPCSPT_ATM) NSPANW(61)=N15+N25+1
      NSPANW(62)=N3+1
      IF (NLSPCSPT_MOL) NSPANW(62)=N16+N25+1
      NSPANW(63)=N3+1
      IF (NLSPCSPT_ION) NSPANW(63)=N17+N25+1
      NSPANW(64)=N3+1
      IF (NLSPCSPT_PHOT) NSPANW(64)=N18+N25+1
      NSPANW(65)=N3+1
      IF (NLSPCSPT_PLS) NSPANW(65)=N19+N25+1
      NSPANW(66)=1
      IF (NLSPCSPT_ATM) NSPANW(66)=N19+N26+1
      NSPANW(67)=1
      IF (NLSPCSPT_MOL) NSPANW(67)=N20+N26+1
      NSPANW(68)=1
      IF (NLSPCSPT_ION) NSPANW(68)=N21+N26+1
      NSPANW(69)=1
      IF (NLSPCSPT_PHOT) NSPANW(69)=N22+N26+1
      NSPANW(70)=1
      IF (NLSPCSPT_PLS) NSPANW(70)=N23+N26+1
      NSPANW(71)=N4+1
      IF (NLSPCSPT_ATM) NSPANW(71)=N27+1
      NSPANW(72)=N4+1
      IF (NLSPCSPT_MOL) NSPANW(72)=N28+1
      NSPANW(73)=N4+1
      IF (NLSPCSPT_ION) NSPANW(73)=N29+1
      NSPANW(74)=N4+1
      IF (NLSPCSPT_PHOT) NSPANW(74)=N30+1
      NSPANW(75)=N4+1
      IF (NLSPCSPT_PLS) NSPANW(75)=N31+1
      NSPANW(76)=0
      NSPANW(77)=0
      NSPANW(78)=0
      NSPANW(79)=0
      NSPANW(80)=0
      NSPANW(81)=0
      NSPANW(82)=N5+1
      NSPANW(83)=N6+1
      NSPANW(84)=1   !PUMPED FLUX

      NSPENW(1)=N2
      NSPENW(2)=N2
      IF (NLSPCSCL_ATM) NSPENW(2)=N8
      NSPENW(3)=N2
      IF (NLSPCSCL_MOL) NSPENW(3)=N9
      NSPENW(4)=N2
      IF (NLSPCSCL_ION) NSPENW(4)=N10
      NSPENW(5)=N2
      IF (NLSPCSCL_PHOT) NSPENW(5)=N11
      NSPENW(6)=N2
      NSPENW(7)=N3
      NSPENW(8)=N3
      IF (NLSPCSCL_ATM) NSPENW(8)=N12
      NSPENW(9)=N3
      IF (NLSPCSCL_MOL) NSPENW(9)=N13
      NSPENW(10)=N3
      IF (NLSPCSCL_ION) NSPENW(10)=N14
      NSPENW(11)=N3
      IF (NLSPCSCL_PHOT) NSPENW(11)=N15
      NSPENW(12)=N3
      NSPENW(13)=N4
      NSPENW(14)=N4
      IF (NLSPCSCL_ATM) NSPENW(14)=N16
      NSPENW(15)=N4
      IF (NLSPCSCL_MOL) NSPENW(15)=N17
      NSPENW(16)=N4
      IF (NLSPCSCL_ION) NSPENW(16)=N18
      NSPENW(17)=N4
      IF (NLSPCSCL_PHOT) NSPENW(17)=N19
      NSPENW(18)=N4
      NSPENW(19)=N1
      NSPENW(20)=N1
      IF (NLSPCSCL_ATM) NSPENW(20)=N20
      NSPENW(21)=N1
      IF (NLSPCSCL_MOL) NSPENW(21)=N21
      NSPENW(22)=N1
      IF (NLSPCSCL_ION) NSPENW(22)=N22
      NSPENW(23)=N1
      IF (NLSPCSCL_PHOT) NSPENW(23)=N23
      NSPENW(24)=N1
      NSPENW(25)=N5
      NSPENW(26)=N2
      NSPENW(27)=N2
      NSPENW(28)=N2
      NSPENW(29)=N2
      NSPENW(30)=N2
      NSPENW(31)=N2
      NSPENW(32)=N3
      NSPENW(33)=N3
      NSPENW(34)=N3
      NSPENW(35)=N3
      NSPENW(36)=N3
      NSPENW(37)=N3
      NSPENW(38)=N4
      NSPENW(39)=N4
      NSPENW(40)=N4
      NSPENW(41)=N4
      NSPENW(42)=N4
      NSPENW(43)=N4
      NSPENW(44)=N1
      NSPENW(45)=N1
      NSPENW(46)=N1
      NSPENW(47)=N1
      NSPENW(48)=N1
      NSPENW(49)=N1
      NSPENW(50)=N5
      NSPENW(51)=N2
      IF (NLSPCSPT_ATM) NSPENW(51)=N8
      NSPENW(52)=N2
      IF (NLSPCSPT_MOL) NSPENW(52)=N9
      NSPENW(53)=N2
      IF (NLSPCSPT_ION) NSPENW(53)=N10
      NSPENW(54)=N2
      IF (NLSPCSPT_PHOT) NSPENW(54)=N11
      NSPENW(55)=N2
      IF (NLSPCSPT_PLS) NSPENW(55)=N11+N24
      NSPENW(56)=N3
      IF (NLSPCSPT_ATM) NSPENW(56)=N12+N24
      NSPENW(57)=N3
      IF (NLSPCSPT_MOL) NSPENW(57)=N13+N24
      NSPENW(58)=N3
      IF (NLSPCSPT_ION) NSPENW(58)=N14+N24
      NSPENW(59)=N3
      IF (NLSPCSPT_PHOT) NSPENW(59)=N15+N24
      NSPENW(60)=N3
      IF (NLSPCSPT_PLS) NSPENW(60)=N15+N25
      NSPENW(61)=N4
      IF (NLSPCSPT_ATM) NSPENW(61)=N16+N25
      NSPENW(62)=N4
      IF (NLSPCSPT_MOL) NSPENW(62)=N17+N25
      NSPENW(63)=N4
      IF (NLSPCSPT_ION) NSPENW(63)=N18+N25
      NSPENW(64)=N4
      IF (NLSPCSPT_PHOT) NSPENW(64)=N19+N25
      NSPENW(65)=N4
      IF (NLSPCSPT_PLS) NSPENW(65)=N19+N26
      NSPENW(66)=N1
      IF (NLSPCSPT_ATM) NSPENW(66)=N20+N26
      NSPENW(67)=N1
      IF (NLSPCSPT_MOL) NSPENW(67)=N21+N26
      NSPENW(68)=N1
      IF (NLSPCSPT_ION) NSPENW(68)=N22+N26
      NSPENW(69)=N1
      IF (NLSPCSPT_PHOT) NSPENW(69)=N23+N26
      NSPENW(70)=N1
      IF (NLSPCSPT_PLS) NSPENW(70)=N27
      NSPENW(71)=N5
      IF (NLSPCSPT_ATM) NSPENW(71)=N28
      NSPENW(72)=N5
      IF (NLSPCSPT_MOL) NSPENW(72)=N29
      NSPENW(73)=N5
      IF (NLSPCSPT_ION) NSPENW(73)=N30
      NSPENW(74)=N5
      IF (NLSPCSPT_PHOT) NSPENW(74)=N31
      NSPENW(75)=N5
      IF (NLSPCSPT_PLS) NSPENW(75)=N32
      NSPENW(76)=0
      NSPENW(77)=0
      NSPENW(78)=0
      NSPENW(79)=0
      NSPENW(80)=0
      NSPENW(81)=0
      NSPENW(82)=N6    ! ADD. SURF. TALLY: N5+1--N6
      NSPENW(83)=N7    ! ALG. SURF. TALLY: N6+1--N7
      NSPENW(84)=N5    ! PUMPED FLUX     : 1   --N5

      DO IPHOT=1,NPHOTI
        ISPZ=IPHOT
        TXTSPW(IPHOT,19)=TEXTS(ISPZ)
        TXTSPW(IPHOT,20)=TEXTS(ISPZ)
        IF (NLSPCSCL_ATM) THEN
          DO IATM=1,NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IATM,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPW(IAD,20))
          END DO
        END IF
        TXTSPW(IPHOT,21)=TEXTS(ISPZ)
        IF (NLSPCSCL_MOL) THEN
          DO IMOL=1,NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IMOL,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPW(IAD,21))
          END DO
        END IF
        TXTSPW(IPHOT,22)=TEXTS(ISPZ)
        IF (NLSPCSCL_ION) THEN
          DO IION=1,NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IION,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPW(IAD,22))
          END DO
        END IF
        TXTSPW(IPHOT,23)=TEXTS(ISPZ)
        IF (NLSPCSCL_PHOT) THEN
          DO JPHOT=1,NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,JPHOT,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(JPHOT),
     .                               TXTSPW(IAD,23))
          END DO
        END IF
        TXTSPW(IPHOT,24)=TEXTS(ISPZ)
        TXTSPW(IPHOT,44)=TEXTS(ISPZ)
        TXTSPW(IPHOT,45)=TEXTS(ISPZ)
        TXTSPW(IPHOT,46)=TEXTS(ISPZ)
        TXTSPW(IPHOT,47)=TEXTS(ISPZ)
        TXTSPW(IPHOT,48)=TEXTS(ISPZ)
        TXTSPW(IPHOT,49)=TEXTS(ISPZ)
        TXTSPW(IPHOT,66)=TEXTS(ISPZ)
        IF (NLSPCSPT_ATM) THEN
          DO IATM=1,NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IATM,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPW(IAD,66))
          END DO
        END IF
        TXTSPW(IPHOT,67)=TEXTS(ISPZ)
        IF (NLSPCSPT_MOL) THEN
          DO IMOL=1,NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IMOL,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPW(IAD,67))
          END DO
        END IF
        TXTSPW(IPHOT,68)=TEXTS(ISPZ)
        IF (NLSPCSPT_ION) THEN
          DO IION=1,NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IION,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPW(IAD,68))
          END DO
        END IF
        TXTSPW(IPHOT,69)=TEXTS(ISPZ)
        IF (NLSPCSPT_PHOT) THEN
          DO JPHOT=1,NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,JPHOT,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(JPHOT),
     .                               TXTSPW(IAD,69))
          END DO
        END IF
        TXTSPW(IPHOT,70)=TEXTS(ISPZ)
        IF (NLSPCSPT_PLS) THEN
          DO IPLS=1,NPLSI
            IAD = EIRENE_INDIRECT_ADDRESS(IPHOT,IPLS,NPHOT)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAMI+IPLS),
     .                               TXTSPW(IAD,70))
          END DO
        END IF
      END DO

      DO IATM=1,NATMI
        ISPZ=NSPH+IATM
        TXTSPW(IATM,1)=TEXTS(ISPZ)
        TXTSPW(IATM,2)=TEXTS(ISPZ)
        IF (NLSPCSCL_ATM) THEN
          DO JATM=1,NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,JATM,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+JATM),
     .                               TXTSPW(IAD,2))
          END DO
        END IF
        TXTSPW(IATM,3)=TEXTS(ISPZ)
        IF (NLSPCSCL_MOL) THEN
          DO IMOL=1,NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IMOL,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPW(IAD,3))
          END DO
        END IF
        TXTSPW(IATM,4)=TEXTS(ISPZ)
        IF (NLSPCSCL_ION) THEN
          DO IION=1,NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IION,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPW(IAD,4))
          END DO
        END IF
        TXTSPW(IATM,5)=TEXTS(ISPZ)
        IF (NLSPCSCL_PHOT) THEN
          DO IPHOT=1,NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IPHOT,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPW(IAD,5))
          END DO
        END IF
        TXTSPW(IATM,6)=TEXTS(ISPZ)
        TXTSPW(IATM,26)=TEXTS(ISPZ)
        TXTSPW(IATM,27)=TEXTS(ISPZ)
        TXTSPW(IATM,28)=TEXTS(ISPZ)
        TXTSPW(IATM,29)=TEXTS(ISPZ)
        TXTSPW(IATM,30)=TEXTS(ISPZ)
        TXTSPW(IATM,31)=TEXTS(ISPZ)
        TXTSPW(IATM,51)=TEXTS(ISPZ)
        IF (NLSPCSPT_ATM) THEN
          DO JATM=1,NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,JATM,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+JATM),
     .                               TXTSPW(IAD,51))
          END DO
        END IF
        TXTSPW(IATM,52)=TEXTS(ISPZ)
        IF (NLSPCSPT_MOL) THEN
          DO IMOL=1,NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IMOL,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPW(IAD,52))
          END DO
        END IF
        TXTSPW(IATM,53)=TEXTS(ISPZ)
        IF (NLSPCSPT_ION) THEN
          DO IION=1,NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IION,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPW(IAD,53))
          END DO
        END IF
        TXTSPW(IATM,54)=TEXTS(ISPZ)
        IF (NLSPCSPT_PHOT) THEN
          DO IPHOT=1,NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IPHOT,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPW(IAD,54))
          END DO
        END IF
        TXTSPW(IATM,55)=TEXTS(ISPZ)
        IF (NLSPCSPT_PLS) THEN
          DO IPLS=1,NPLSI
            IAD = EIRENE_INDIRECT_ADDRESS(IATM,IPLS,NATM)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAMI+IPLS),
     .                               TXTSPW(IAD,55))
          END DO
        END IF
      END DO
C
      DO IMOL=1,NMOLI
        ISPZ=NSPA+IMOL
        TXTSPW(IMOL,7)=TEXTS(ISPZ)
        TXTSPW(IMOL,8)=TEXTS(ISPZ)
        IF (NLSPCSCL_ATM) THEN
          DO IATM=1,NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IATM,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPW(IAD,8))
          END DO
        END IF
        TXTSPW(IMOL,9)=TEXTS(ISPZ)
        IF (NLSPCSCL_MOL) THEN
          DO JMOL=1,NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,JMOL,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+JMOL),
     .                               TXTSPW(IAD,9))
          END DO
        END IF
        TXTSPW(IMOL,10)=TEXTS(ISPZ)
        IF (NLSPCSCL_ION) THEN
          DO IION=1,NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IION,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPW(IAD,10))
          END DO
        END IF
        TXTSPW(IMOL,11)=TEXTS(ISPZ)
        IF (NLSPCSCL_PHOT) THEN
          DO IPHOT=1,NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IPHOT,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPW(IAD,11))
          END DO
        END IF
        TXTSPW(IMOL,12)=TEXTS(ISPZ)
        TXTSPW(IMOL,32)=TEXTS(ISPZ)
        TXTSPW(IMOL,33)=TEXTS(ISPZ)
        TXTSPW(IMOL,34)=TEXTS(ISPZ)
        TXTSPW(IMOL,35)=TEXTS(ISPZ)
        TXTSPW(IMOL,36)=TEXTS(ISPZ)
        TXTSPW(IMOL,37)=TEXTS(ISPZ)
        TXTSPW(IMOL,56)=TEXTS(ISPZ)
        IF (NLSPCSPT_ATM) THEN
          DO IATM=1,NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IATM,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPW(IAD,56))
          END DO
        END IF
        TXTSPW(IMOL,57)=TEXTS(ISPZ)
        IF (NLSPCSPT_MOL) THEN
          DO JMOL=1,NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,JMOL,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+JMOL),
     .                               TXTSPW(IAD,57))
          END DO
        END IF
        TXTSPW(IMOL,58)=TEXTS(ISPZ)
        IF (NLSPCSPT_ION) THEN
          DO IION=1,NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IION,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPW(IAD,58))
          END DO
        END IF
        TXTSPW(IMOL,59)=TEXTS(ISPZ)
        IF (NLSPCSPT_PHOT) THEN
          DO IPHOT=1,NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IPHOT,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPW(IAD,59))
          END DO
        END IF
        TXTSPW(IMOL,60)=TEXTS(ISPZ)
        IF (NLSPCSPT_PLS) THEN
          DO IPLS=1,NPLSI
            IAD = EIRENE_INDIRECT_ADDRESS(IMOL,IPLS,NMOL)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAMI+IPLS),
     .                               TXTSPW(IAD,60))
          END DO
        END IF
      END DO
C
      DO IION=1,NIONI
        ISPZ=NSPAM+IION
        TXTSPW(IION,13)=TEXTS(ISPZ)
        TXTSPW(IION,14)=TEXTS(ISPZ)
        IF (NLSPCSCL_ATM) THEN
          DO IATM=1,NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IATM,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPW(IAD,14))
          END DO
        END IF
        TXTSPW(IION,15)=TEXTS(ISPZ)
        IF (NLSPCSCL_MOL) THEN
          DO IMOL=1,NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IMOL,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPW(IAD,15))
          END DO
        END IF
        TXTSPW(IION,16)=TEXTS(ISPZ)
        IF (NLSPCSCL_ION) THEN
          DO JION=1,NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,JION,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+JION),
     .                               TXTSPW(IAD,16))
          END DO
        END IF
        TXTSPW(IION,17)=TEXTS(ISPZ)
        IF (NLSPCSCL_PHOT) THEN
          DO IPHOT=1,NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IPHOT,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPW(IAD,17))
          END DO
        END IF
        TXTSPW(IION,18)=TEXTS(ISPZ)
        TXTSPW(IION,38)=TEXTS(ISPZ)
        TXTSPW(IION,39)=TEXTS(ISPZ)
        TXTSPW(IION,40)=TEXTS(ISPZ)
        TXTSPW(IION,41)=TEXTS(ISPZ)
        TXTSPW(IION,42)=TEXTS(ISPZ)
        TXTSPW(IION,43)=TEXTS(ISPZ)
        TXTSPW(IION,61)=TEXTS(ISPZ)
        IF (NLSPCSPT_ATM) THEN
          DO IATM=1,NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IATM,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPW(IAD,61))
          END DO
        END IF
        TXTSPW(IION,62)=TEXTS(ISPZ)
        IF (NLSPCSPT_MOL) THEN
          DO IMOL=1,NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IMOL,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPW(IAD,62))
          END DO
        END IF
        TXTSPW(IION,63)=TEXTS(ISPZ)
        IF (NLSPCSPT_ION) THEN
          DO JION=1,NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,JION,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+JION),
     .                               TXTSPW(IAD,63))
          END DO
        END IF
        TXTSPW(IION,64)=TEXTS(ISPZ)
        IF (NLSPCSPT_PHOT) THEN
          DO IPHOT=1,NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IPHOT,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPW(IAD,64))
          END DO
        END IF
        TXTSPW(IION,65)=TEXTS(ISPZ)
        IF (NLSPCSPT_PLS) THEN
          DO IPLS=1,NPLSI
            IAD = EIRENE_INDIRECT_ADDRESS(IION,IPLS,NION)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAMI+IPLS),
     .                               TXTSPW(IAD,65))
          END DO
        END IF
      END DO
C
      DO IPLS=1,NPLSI
        ISPZ=NSPAMI+IPLS
        TXTSPW(IPLS,25)=TEXTS(ISPZ)
        TXTSPW(IPLS,50)=TEXTS(ISPZ)
        TXTSPW(IPLS,71)=TEXTS(ISPZ)
        IF (NLSPCSPT_ATM) THEN
          DO IATM=1,NATMI
            IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IATM,NPLS)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPH+IATM),
     .                               TXTSPW(IAD,71))
          END DO
        END IF
        TXTSPW(IPLS,72)=TEXTS(ISPZ)
        IF (NLSPCSPT_MOL) THEN
          DO IMOL=1,NMOLI
            IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IMOL,NPLS)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPA+IMOL),
     .                               TXTSPW(IAD,72))
          END DO
        END IF
        TXTSPW(IPLS,73)=TEXTS(ISPZ)
        IF (NLSPCSPT_ION) THEN
          DO IION=1,NIONI
            IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IION,NPLS)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAM+IION),
     .                               TXTSPW(IAD,73))
          END DO
        END IF
        TXTSPW(IPLS,74)=TEXTS(ISPZ)
        IF (NLSPCSPT_PHOT) THEN
          DO IPHOT=1,NPHOTI
            IAD = EIRENE_INDIRECT_ADDRESS(IPLS,IPHOT,NPLS)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(IPHOT),
     .                               TXTSPW(IAD,74))
          END DO
        END IF
        TXTSPW(IPLS,75)=TEXTS(ISPZ)
        IF (NLSPCSPT_PLS) THEN
          DO JPLS=1,NPLSI
            IAD = EIRENE_INDIRECT_ADDRESS(IPLS,JPLS,NPLS)
            CALL EIRENE_TEXT_COMBINE(TEXTS(ISPZ),TEXTS(NSPAMI+JPLS),
     .                               TXTSPW(IAD,75))
          END DO
        END IF
      END DO
C
      TXTSPW(1,42)='                        '
      TXTSPW(1,43)='                        '
      TXTSPW(1,44)='                        '
      TXTSPW(1,45)='                        '
C
      RETURN
      END SUBROUTINE EIRENE_STTXT1


      SUBROUTINE EIRENE_TEXT_COMBINE(TXT1,TXT2,TXTOUT)

      IMPLICIT NONE

      CHARACTER(8), INTENT(IN) :: TXT1, TXT2
      CHARACTER(24), INTENT(OUT) :: TXTOUT
      CHARACTER(8) :: TT1, TT2
      CHARACTER(50) :: TTOUT
      INTEGER :: LL1, LL2, LL

      LL1 = LEN_TRIM(ADJUSTL(TRIM(TXT1)))
      TT1(1:LL1) = ADJUSTL(TRIM(TXT1))

      LL2 = LEN_TRIM(ADJUSTL(TRIM(TXT2)))
      TT2(1:LL2) = ADJUSTL(TRIM(TXT2))

      TTOUT = TT1(1:LL1) // ' FROM ' // TT2(1:LL2)

      LL = MAX(24, LL1 + 6 + LL2)
      TXTOUT(1:LL) = TTOUT(1:LL)

      RETURN
      END SUBROUTINE EIRENE_TEXT_COMBINE
