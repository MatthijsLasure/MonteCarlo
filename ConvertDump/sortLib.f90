module sortLib
    implicit none

    CONTAINS

    SUBROUTINE SORT(AR, IX)
        DOUBLE PRECISION, DIMENSION(:) :: AR
        INTEGER, DIMENSION(:) :: IX

        INTEGER :: I, J, K, N

        N = SIZE(AR)

        DO J = N, 2, -1
            DO I = 2,J - 1
                IF ( AR(I) .LT. AR(I-1)) CALL SWAP(AR, IX, I, I-1)
            END DO
        END DO

    END SUBROUTINE SORT

    SUBROUTINE SWAP(AR, IX, I1, I2)
        DOUBLE PRECISION, DIMENSION(:) :: AR
        INTEGER, DIMENSION(:) :: IX
        DOUBLE PRECISION :: ATEMP
        INTEGER :: ITEMP, I1, I2


        ATEMP = AR(I1)
        ITEMP = IX(I1)

        AR(I1) = AR(I2)
        IX(I1) = IX(I2)

        AR(I2) = ATEMP
        IX(I2) = ITEMP

    END SUBROUTINE SWAP
end module sortLib
