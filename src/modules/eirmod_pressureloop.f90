MODULE EIRMOD_PRESSURELOOP
  USE EIRMOD_PRECISION
  USE EIRMOD_CGRID, ONLY: NSBOX_TAL

  IMPLICIT NONE

  TYPE pressureFeedback
    INTEGER:: id = 0 !ID of PFL
    LOGICAL:: active = .FALSE. !Indicates if the loop is active
    REAL(DP):: pressureRef = 0.0_DP !Reference pressure
    REAL(DP):: integralError = 0.0_DP !Accumulative error
    REAL(DP):: Kp = 0.0_DP, Ki = 0.0_DP !Proportional and Integral coefficients
    INTEGER:: cellRef = 0 !Id of reference cell
    REAL(DP):: delta = 0.0_DP !Increment of RECYCT obtained from the PI loop

  END TYPE pressureFeedback

  TYPE(pressureFeedback), ALLOCATABLE:: RPRESSFED(:)

  CONTAINS
    SUBROUTINE initPressureFeedback(PFL, id, cellRef, pressureRef)
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY : iunout

      IMPLICIT NONE
      INTEGER, INTENT(in):: id
      INTEGER, INTENT(in):: cellRef
      REAL(DP), INTENT(in):: pressureRef
      TYPE(pressureFeedback), INTENT(inout):: PFL
      EXTERNAL :: EIRENE_EXIT_OWN

      !Surface cell id
      PFL%id = id

      !Gets the cell from the reference pressure position 
      IF (cellRef > 0 .AND. cellRef <= NSBOX_TAL) THEN
        PFL%cellRef = cellRef
      ELSE
        WRITE(iunout, '(A)') 'Incorrect reference cell for PFL:'
        WRITE(iunout, '(A,I3)') 'id = ', id
        WRITE(iunout, '(A,I3)') 'cell = ', cellRef
        CALL EIRENE_EXIT_OWN(1)

      END IF

      PFL%active = .TRUE.

      PFL%pressureRef = pressureRef

      PFL%integralError = 0.0_DP
      
      PFL%delta = 0.0_DP

      !TODO: make these parameters user input (JSON)
      PFL%Kp = 1.D-2
      PFL%Ki = 1.D-5
      
    END SUBROUTINE initPressureFeedback

    SUBROUTINE calculatePressureRatio(PFL)
      USE EIRMOD_PRECISION
      USE EIRMOD_CESTIM
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA

      IMPLICIT NONE
      TYPE(pressureFeedback), INTENT(inout):: PFL
      REAL(DP):: pressure
      INTEGER:: i
      REAL(DP):: error
      REAL(DP):: ergVol2Pa = 1.D-1


      !Calculate pressure in reference cell (in Pa)
      !Atoms
      pressure = 0._DP
      DO i = 1, NATM
        IF (PDENA(i, PFL%cellRef) > EPS30) THEN
        pressure = pressure + EV_TO_J*1.D6*EDENA(i, PFL%cellRef)/1.5_DP -  &
                              ergVol2Pa*(VXDENA(i, PFL%cellRef)**2 +       &
                                         VYDENA(i, PFL%cellRef)**2 +       &
                                         VZDENA(i, PFL%cellRef)**2)/       &
                              (3._DP*AMUA*RMASSA(i)*PDENA(i, PFL%cellRef))

        END IF

      END DO

      !Molecules
      DO i =1, NMOL
        IF (PDENM(i, PFL%cellRef) > EPS30) THEN
        pressure = pressure + EV_TO_J*1.D6*EDENM(i, PFL%cellRef)/1.5_DP -  &
                              ergVol2Pa*(VXDENM(i, PFL%cellRef)**2 +       &
                                         VYDENM(i, PFL%cellRef)**2 +       &
                                         VZDENM(i, PFL%cellRef)**2)/       &
                              (3._DP*AMUA*RMASSM(i)*PDENM(i, PFL%cellRef))

        END IF


      END DO

      WRITE(iunout,*) PFL%id, pressure, PFL%cellRef

      !Calculate the pressure ratio required to update the RECYCT parameter
      error = PFL%pressureRef - pressure
      PFL%integralError = PFL%integralError + error
      PFL%delta = PFL%Kp*error + PFL%Ki*PFL%integralError
      WRITE(iunout,*) "PFL: surf, error, delta:", PFL%id, error, PFL%delta

    END SUBROUTINE calculatePressureRatio

    SUBROUTINE updatePressureFeedback(PFL)
      USE EIRMOD_CLGIN
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT

      IMPLICIT NONE
      TYPE(pressureFeedback), INTENT(inout):: PFL
      IF (PFL%active) THEN
        !Update pressure ratio
        CALL calculatePressureRatio(PFL)

        !Update RECYCT value
        RECYCT(:,PFL%id) = RECYCT(:,PFL%id) + PFL%delta

        !Limits RECYCT to [0 1]
        IF (ANY(RECYCT(:,PFL%id) > 1.D0)) THEN
          RECYCT(:,PFL%id) = 1.D0

        ELSEIF (ANY(RECYCT(:,PFL%id) < 0.0_DP)) THEN
          RECYCT(:,PFL%id) = 0.0_DP

        END IF

        WRITE(iunout,*) "PFL: surf, RECYCT: ", PFL%id, RECYCT(:,PFL%id)

      END IF

    END SUBROUTINE updatePressureFeedback

END MODULE EIRMOD_PRESSURELOOP
