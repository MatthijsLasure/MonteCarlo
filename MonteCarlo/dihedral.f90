MODULE dihedral

    !implicit none
    USE vector_class

    CONTAINS

    FUNCTION GETDIHEDRAL(MOL, I1, I2, I3, I4) RESULT(HOEK)
        TYPE (VECTOR), DIMENSION(:) :: MOL
        INTEGER:: I1
        INTEGER:: I2
        INTEGER:: I3
        INTEGER:: I4
        DOUBLE PRECISION:: HOEK
        DOUBLE PRECISION:: TEMP
        DOUBLE PRECISION:: SIGN
        TYPE (VECTOR):: V1
        TYPE (VECTOR):: V2
        TYPE (VECTOR):: V3
        TYPE (VECTOR):: V4
        TYPE (VECTOR):: B1
        TYPE (VECTOR):: B2
        TYPE (VECTOR):: B3
        TYPE (VECTOR):: N1
        TYPE (VECTOR):: N2

        V1 = MOL(I1)
        V2 = MOL(I2)
        V3 = MOL(I3)
        V4 = MOL(I4)

        B1 = V2 - V1
        B2 = V3 - V2
        B3 = V4 - V3

        N1 = B1 * B2
        N2 = B2 * B3

        TEMP = DOT(N1, N2) / (LENGTH(N1) * LENGTH(N2))
        HOEK = DACOS(TEMP)

        SIGN = DOT(B1, N2)
        IF (SIGN .LT. 0.D0) HOEK = -1.D0 * HOEK

    END FUNCTION GETDIHEDRAL

!====================================================================
!====================================================================

    FUNCTION SETDIHEDRAL(MOL, I1, I2, I3, I4, HOEK) RESULT(MOLOUT)
        TYPE (VECTOR), DIMENSION(:) :: MOL
        INTEGER, DIMENSION(SIZE(MOL)) :: FRAG ! Welk fragment
        TYPE (VECTOR), DIMENSION(SIZE(MOL)) :: MOLOUT
        INTEGER:: I
        INTEGER:: I1
        INTEGER:: I2
        INTEGER:: I3
        INTEGER:: I4
        INTEGER:: NATOM
        DOUBLE PRECISION:: HOEK
        DOUBLE PRECISION:: dist1
        DOUBLE PRECISION:: dist2
        DOUBLE PRECISION:: cosH
        DOUBLE PRECISION:: sinH
        TYPE (VECTOR):: V1
        TYPE (VECTOR):: V2
        TYPE (VECTOR):: V3
        TYPE (VECTOR):: V4
        TYPE (VECTOR):: V
        TYPE (VECTOR):: K ! Rotatie as
        TYPE (VECTOR):: ATOM
        TYPE (VECTOR):: R1
        TYPE (VECTOR):: R2
        LOGICAL:: mayGo = .TRUE.

        DOUBLE PRECISION, PARAMETER :: BONDL = 2.05
        DOUBLE PRECISION, PARAMETER :: BONDL2 = BONDL * BONDL

        NATOM = SIZE(MOL)

        V1 = MOL(I1)
        V2 = MOL(I2)
        V3 = MOL(I3)
        V4 = MOL(I4)

        ! Hoek
        cosH = DCOS(HOEK)
        sinH = DSIN(HOEK)

        ! Rotatie as
        K = V3 - V2
        K = norm(K)

        !Stel alle atomen in op fragment 99
        DO I = 1, NATOM
            FRAG(I) = 99
            MOLOUT(I) = MOL(I)
        END DO

        FRAG(I2) = 0 ! As waarrond gedraaid wordt
        FRAG(I3) = 0

        first: DO I = 1, NATOM
            IF (FRAG(I) .NE. 0) THEN
                R1 = MOL(I) - V2
                R2 = MOL(I) - V3
                dist1 = length_sq(R1)
                dist2 = length_sq(R2)

                ! Atoom hoort bij fragment 1
                IF ( (dist1 .LT. BONDL2) .AND. (dist2 .GE. 4.21D0) ) THEN
                    FRAG(I) = 1
                ! Atoom hoort bij fragment 2
                ELSEIF ( (dist1 .GE. 4.21D0) .AND. (dist2 .LT. BONDL2) ) THEN
                    FRAG(I) = 2
                END IF
            END IF
        END DO first

        WRITE (*,*) "findAtoms"

        ! Vind de aangrenzende atomen
        findAtoms: DO WHILE(mayGo)
            atomLoop: DO I = 1, NATOM
                IF (FRAG(I) .EQ. 99) THEN
                    DO J = 1, NATOM
                        DIST1 = HUGE(DIST1)
                        ! Geen vergelijking met zichzelf, en de andere is al toegewezen
                        IF (J .NE. I .AND. ( FRAG(J) .EQ. 1 .OR. FRAG(J) .EQ. 2 ) ) THEN
                            R1 = MOL(I) - MOL(J)
                            DIST1 = length_sq(R1)

                            ! I behoort tot hetzelfde fragment als J: opschrijven dus
                            IF (DIST1 .LT. BONDL2) THEN
                                FRAG(I) = FRAG(J)
                                CYCLE atomLoop
                            END IF
                        END IF
                        WRITE (500,*) I, J, FRAG(I), FRAG(J), DIST1
                    END DO
                END IF
            END DO atomLoop

            mayGo = .FALSE.
            ! Check of alles is toegewezen
            DO I = 1, NATOM
                IF (FRAG(I) .EQ. 99) mayGo = .TRUE.
            END DO
        END DO findAtoms
        DO I = 1,NATOM
            WRITE (*,*) I, FRAG(I)
        END DO

        WRITE (*,*) "rotlus"
        ! Alleen fragment 2 wordt geroteerd
        rotlus: DO I = 1, NATOM
            IF( FRAG(I) .EQ. 2) THEN
                ATOM = MOL(I)
                V = ATOM - V3
                ATOM = V * cosH + ( K * V) * sinH + K * DOT(K, V) * (1.D0 - cosH)
                !ATOM = ATOM + sinH * (K * ATOM) + (1.D0 - cosH) * (K * (K * ATOM))
                MOLOUT(I) = ATOM + V3
            END IF
        END DO rotlus

    END FUNCTION SETDIHEDRAL
END MODULE dihedral
