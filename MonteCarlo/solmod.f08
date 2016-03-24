module solmod

    USE RANDGEN
    USE DIHEDRAL
    USE VECTOR_CLASS
    USE LIB

    implicit none

    contains

    FUNCTION SOLUTE_INIT(SOL, DIHEDRAL) RESULT(SOLROT)

        TYPE (VECTOR), DIMENSION(:), INTENT(IN)  :: SOL
        TYPE (VECTOR), DIMENSION(SIZE(SOL))      :: SOLROT
        INTEGER, DIMENSION(:,:), INTENT(IN)      :: DIHEDRAL
        DOUBLE PRECISION :: HOEK, HOEK_PRE, HOEK_POST
        INTEGER :: A1, A2, A3, A4, NATOM, I, J, M, NDIH
        INTEGER, DIMENSION(8) :: NEIGHBOURS

        NATOM = SIZE(SOL)
        NDIH = SIZE(DIHEDRAL)

        I = INT(RAND() * NDIH) + 1 ! Willeukeurige binding selecteren
        A2 = DIHEDRAL(I,1)
        A3 = DIHEDRAL(I,2)

        A1 = FIND(A2, A3, SOL, NATOM)
        A4 = FIND(A3, A2, SOL, NATOM)

        ! DRAAIEN
        HOEK = RAND() * 360 - 180
        HOEK_PRE = GETDIHEDRAL(SOL, A1, A2, A3, A4) * 180 / PI
        SOLROT = SETDIHEDRAL(SOL, A1, A2, A3, A4, HOEK * PI / 180)
        HOEK_POST = GETDIHEDRAL(SOLROT, A1, A2, A3, A4) * 180 / PI

        WRITE (*,"(A, F5.2, A, F5.2, A, F5.2)") "Dihedral change from ", HOEK_PRE, " with ", HOEK, " to ", HOEK_POST

    END FUNCTION SOLUTE_INIT

    FUNCTION FIND(A1, A2, SOL, NATOM) RESULT(A)
    TYPE (VECTOR), DIMENSION(:), INTENT(IN)  :: SOL
        INTEGER :: A1, A2, A, NATOM
        INTEGER :: I

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
end module solmod
