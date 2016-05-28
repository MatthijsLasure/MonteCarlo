!====================================================================
! lib.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Library met nuttige functies
!====================================================================
MODULE lib
    USE vector_class
    USE randgen

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: PI = 4.D0 * DATAN(1.D0)
    DOUBLE PRECISION, PARAMETER :: TAU = 2.D0 * PI
    DOUBLE PRECISION, PARAMETER :: E = 2.71828182845904523536028747135266249775724709369995
    DOUBLE PRECISION, PARAMETER :: EPSILON0 = 8.854187817620 ! C^2 N^-1 m^-2
    DOUBLE PRECISION, PARAMETER :: HARTREE2KJMOL = 4.359744 * 6.02214085 * 100 ! Hartree -> KJ/mol

    CONTAINS

    ! RotMatrix
    !==========
    ! Roteer een set coördinaten rond de CoM met gegeven hoeken
    ! Hoek: in graden!
    FUNCTION RotMatrix(COM, RELPOS, HOEK) RESULT(AbsPos)
        TYPE (vector), INTENT(IN) :: COM
        TYPE (vector), INTENT(IN) :: HOEK
        TYPE (vector), DIMENSION(:), INTENT(IN) :: RELPOS
        TYPE (vector), DIMENSION(size(RELPOS)) :: ABSPOS

        ! Variabelen
        !===========
        DOUBLE PRECISION, DIMENSION(3,3) :: RM ! Rotatiematrix
        DOUBLE PRECISION, DIMENSION(3) :: TEMP ! Tijdelijke array
        TYPE (vector):: V ! Tijdelijke vector

        DOUBLE PRECISION:: COS1
        DOUBLE PRECISION:: COS2
        DOUBLE PRECISION:: COS3
        DOUBLE PRECISION:: SIN1
        DOUBLE PRECISION:: SIN2
        DOUBLE PRECISION:: SIN3
        INTEGER:: N ! Aantal atomen
        INTEGER:: I

        N = size(RELPOS)

        ! Bereken cos / sin
        !==================
        COS1 = cos(HOEK%x)
        COS2 = cos(HOEK%y)
        COS3 = cos(HOEK%z)
        SIN1 = sin(HOEK%x)
        SIN2 = sin(HOEK%y)
        SIN3 = sin(HOEK%z)

        ! Opstellen matrix
        !=================
        ! Rij 1
        RM(1,1) = COS1 * COS2 - SIN1 * SIN2 * COS3
        RM(1,2) = SIN1 * COS2 + COS1 * SIN2 * COS3
        RM(1,3) = SIN2 * SIN3
        ! Rij 2
        RM(2,1) = -1.D0 * COS1 * SIN2 - SIN1 * COS2 * COS3
        RM(2,2) = -1.D0 * SIN1 * SIN2 + COS1 * COS2 * COS3
        RM(2,3) = COS2 * SIN3
        ! Rij 3
        RM(3,1) = SIN1 * SIN3
        RM(3,2) = -1.D0 * COS1 * SIN3
        RM(3,3) = COS3

        ! Roteren + transleren atomen
        !============================
        DO I=1,N
            TEMP = getArray(RelPos(I))
            TEMP = matmul(RM, TEMP)
            V = fromArray(TEMP)
            ABSPOS(I) = COM + V
        END DO

    END FUNCTION RotMatrix

    ! findSym
    !========
    FUNCTION findSym(TYPE, SYM) RESULT(POS)
        CHARACTER*4:: TYPE
        INTEGER:: POS
        INTEGER:: I
        INTEGER:: N
        CHARACTER*4, DIMENSION(:) :: SYM

        POS = 0

        N = size(SYM)
        DO I=1,N
            IF( sym(I) .EQ. TYPE) THEN
                POS = I
                RETURN
            END IF
        END DO
        IF (POS .EQ. 0) THEN
            POS = 1
            write (500, *) "Error in findSym: did not find", TYPE
        END IF
    END FUNCTION findSym

    ! RandVec: random vector with max total length of MAX
    !====================================================
    FUNCTION randVec(MAX) RESULT(rv)
        DOUBLE PRECISION, INTENT(IN) :: MAX
        TYPE (vector):: RV

        RV%x = (RAND() - 0.5D0) * 2
        RV%y = (RAND() - 0.5D0) * 2
        RV%z = (RAND() - 0.5D0) * 2
        RV = RV / sqrt(3.D0)
        RV = RV * MAX
    END FUNCTION randVec

    ! RandVecHoek: maakt vector met willekeurige getallen tussen -max en +max
    !========================================================================
    FUNCTION randVecHoek(MAX) RESULT(rv)
        DOUBLE PRECISION, INTENT(IN) :: MAX
        TYPE (vector):: RV

        RV%x = (RAND() - 0.5D0) * 2 * MAX
        RV%y = (RAND() - 0.5D0) * 2 * MAX
        RV%z = (RAND() - 0.5D0) * 2 * MAX
    END FUNCTION randVecHoek

    END MODULE lib
