MODULE solmod

    USE RANDGEN
    USE DIHEDRAL
    USE VECTOR_CLASS
    USE LIB

    IMPLICIT NONE

    CONTAINS

    FUNCTION SOLUTE_INIT(SOL, SYM, DIHEDRAL, DROTSOLV) RESULT(SOLROT)

        TYPE (VECTOR), DIMENSION(:), INTENT(IN)  :: SOL
        TYPE (VECTOR), DIMENSION(SIZE(SOL))      :: SOLROT
        CHARACTER*4, DIMENSION(:), INTENT(IN)    :: SYM
        INTEGER, DIMENSION(:,:), INTENT(IN)      :: DIHEDRAL
        DOUBLE PRECISION:: HOEK
        DOUBLE PRECISION:: HOEK_PRE
        DOUBLE PRECISION:: HOEK_POST
        DOUBLE PRECISION:: DROTSOLV
        INTEGER:: A1
        INTEGER:: A2
        INTEGER:: A3
        INTEGER:: A4
        INTEGER:: NATOM
        INTEGER:: I
        INTEGER:: J
        INTEGER:: M
        INTEGER:: NDIH
        INTEGER, DIMENSION(8) :: NEIGHBOURS

        NATOM = SIZE(SOL)
        NDIH = SIZE(DIHEDRAL) / 2

        I = INT(RAND() * NDIH) + 1 ! Willeukeurige binding selecteren
        A2 = DIHEDRAL(I,1)
        A3 = DIHEDRAL(I,2)

        A1 = FIND(A2, A3, SOL, NATOM)
        A4 = FIND(A3, A2, SOL, NATOM)

        ! DRAAIEN
        HOEK = RAND() * 2 * DROTSOLV - DROTSOLV
        HOEK_PRE = GETDIHEDRAL(SOL, A1, A2, A3, A4) * 180 / PI
        SOLROT = SETDIHEDRAL(SOL, SYM, A1, A2, A3, A4, HOEK * PI / 180)
        HOEK_POST = GETDIHEDRAL(SOLROT, A1, A2, A3, A4) * 180 / PI

        !WRITE (*,"(A,I2.2,A3,I2.2)") "Rotation around axis ", A2, " - ", A3
        WRITE (*,"(A,I2.2,A3,I2.2,A3,I2.2,A3,I2.2)") "Dihedral angle with ", A1, " - ", A2, " - ", A3, " - ", A4
        WRITE (*,"(A, F7.2, A, F7.2, A, F7.2)") "Dihedral change from ", HOEK_PRE, " with ", HOEK, " to ", HOEK_POST

    END FUNCTION SOLUTE_INIT

    FUNCTION SOLUTE_METROPOLIS(SOL, PRE_EN, POST_EN, TEMP) RESULT(isOK)
        TYPE (VECTOR), DIMENSION(:), INTENT(IN)  :: SOL
        DOUBLE PRECISION:: PRE_EN
        DOUBLE PRECISION:: POST_EN
        DOUBLE PRECISION:: TEMP
        DOUBLE PRECISION:: EXPONENT
        DOUBLE PRECISION:: RV
        DOUBLE PRECISION:: KANS
        DOUBLE PRECISION:: DELTA
        LOGICAL:: ISOK

        DELTA = POST_EN - PRE_EN
        EXPONENT = -1.D0 * DELTA  * 1000.D0 / (8.314D0 * TEMP)
        IF (EXPONENT .LT. -75.D0) THEN ! e^-75 < 3*10^-33: 0% kans anyway
        !write(500,*) "Large Exponent!", I
            KANS = 0.D0
            RV = 1.D0 ! Skip rand() voor cpu tijd besparing.
        ELSE IF(EXPONENT .GE. 0.D0) THEN ! Lager in energie, dus 100% kans
            KANS = 1.D0
            RV = 0.D0
        ELSE
            KANS = E ** EXPONENT
            RV = RAND()
            IF (KANS .GT. 1.D0) KANS = 1.D0
        END IF

        ! Bepaal if succesvol -> volgende config
        IF(RV .LE. KANS) THEN ! Succes!
            WRITE (*,*) "New solute accepted!", PRE_EN, POST_EN
            ISOK = .TRUE.
        ELSE ! Fail!
            WRITE (*,*) "New solute rejected!", PRE_EN, POST_EN
            ISOK = .FALSE.
        END IF

    END FUNCTION SOLUTE_METROPOLIS


    FUNCTION FIND(A1, A2, SOL, NATOM) RESULT(A)
    TYPE (VECTOR), DIMENSION(:), INTENT(IN)  :: SOL
        INTEGER:: A1
        INTEGER:: A2
        INTEGER:: A
        INTEGER:: NATOM
        INTEGER:: I

        A = 0
        DO I=1, NATOM
            IF (I .NE. A1 .AND. I .NE. A2) THEN
                IF (getDist(SOL(I), SOL(A1)) .LT. 1.90) THEN
                    A = I
                    RETURN
                END IF
            END IF
        END DO
        IF (A .EQ. 0) THEN
            WRITE (500,*) "No match find in dihedral-find with", A1, A2
            STOP
        END IF
    END FUNCTION FIND
END MODULE solmod
