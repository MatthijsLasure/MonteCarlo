MODULE LennardJones

    USE vector_class
    USE lib

    IMPLICIT NONE

    CONTAINS

!====================================================================

    ! calcLJ: bereken alle interacties via LJ potentiaal
    ! Als j = 0: bereken solvent - solute
    ! V = E_iE_j Q_i*Q_j / 4pi*epsilon0*r_ij + E_iE_j 4 epsilon_ij * [ (sigma_ij/r_ij)^12 - (sigma_ij/r_ij)^6 ]
    !
    ! SYNTAX:
    ! mol1: eerste molecule met abs coord (gebruik RotMat)
    ! mol2: idem als mol1
    ! sym1: atoomtypes voor mol1
    ! sym2: atoomtypes voor mol2
    ! table_sym: atoomtypes voor parameters
    ! table_Q/e/s: parameters voor LJ

    SUBROUTINE calcLJ(MOL1, MOL2, SYM1, TABLE, EN, BOXL, BOXL2)
        ! INPUT
        TYPE (vector), DIMENSION(:), INTENT(IN)         :: MOL1, MOL2 ! absolute coords!
        CHARACTER*4, DIMENSION(:), INTENT(IN)           :: SYM1 ! Atoomtypes
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)    :: TABLE ! params
        DOUBLE PRECISION                                :: CONV, BOXL, BOXL2

        ! OUTPUT
        DOUBLE PRECISION    :: EN

        ! TEMP
        TYPE (vector)       :: R
        DOUBLE PRECISION    :: RV, E, S ! temp params
        DOUBLE PRECISION    :: SUML, SUMR ! sommen
        DOUBLE PRECISION, DIMENSION(2)       :: TEMP ! tijdelijke optelling
        INTEGER             :: I, J, N1, N2 ! loop vars, totale grootte arrays
        INTEGER             :: A, B ! welke atoomsoort

        ! Omzetting naar kJ/mol
        CONV = 6.023 * (1.60217646)**2
        CONV = CONV / 4
        CONV = CONV / PI
        CONV = CONV / 8.854187817620 ! J/mol * 10^-5
        CONV = CONV * 10**(-8) ! kJ/mol

        N1 = size(MOL1)
        N2 = size(MOL2)

        SUML = 0.D0
        SUMR = 0.D0

        DO I=1,N1
            IF(sym1(I) .NE. "H") THEN
                DO J=1,N2
                    IF(sym1(J) .NE. "H") THEN
                        R = mol2(J) - mol1(I)

                        TEMP = LJ_FORMULA(table(I,2), table(J,2), &
                        table(I,3), table(J,3), table(I,1), table(J,1), &
                        R, BOXL, BOXL2)
                        SUML = SUML + TEMP(1)
                        SUMR = SUMR + TEMP(2)
                    END IF
                END DO
            END IF
        END DO
        SUML = SUML * CONV

        EN = SUML + SUMR

    END SUBROUTINE calcLJ

!====================================================================
!====================================================================

SUBROUTINE CalcLJ_SV(MOL1, MOL2, SYM1, TABLE_DMSO, TABLE_SOL, EN, BOXL, BOXL2)
        ! INPUT
        TYPE (vector), DIMENSION(:), INTENT(IN)         :: MOL1, MOL2 ! absolute coords!
        CHARACTER*4, DIMENSION(:), INTENT(IN)           :: SYM1 ! Atoomtypes
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)    :: TABLE_DMSO, TABLE_SOL ! params
        DOUBLE PRECISION                                :: CONV, BOXL, BOXL2

        ! OUTPUT
        DOUBLE PRECISION    :: EN

        ! TEMP
        TYPE (vector)       :: R
        DOUBLE PRECISION    :: RV, E, S ! temp params
        DOUBLE PRECISION    :: SUML, SUMR ! sommen
        DOUBLE PRECISION, DIMENSION(2)       :: TEMP ! tijdelijke optelling
        INTEGER             :: I, J, N1, N2 ! loop vars, totale grootte arrays
        INTEGER             :: A, B ! welke atoomsoort

        ! Omzetting naar kJ/mol
        CONV = 6.023 * (1.60217646)**2
        CONV = CONV / 4
        CONV = CONV / PI
        CONV = CONV / 8.854187817620 ! J/mol * 10^-5
        CONV = CONV * 10**(-8) ! kJ/mol

        N1 = size(MOL1)
        N2 = size(MOL2)

        SUML = 0.D0
        SUMR = 0.D0

        DO I=1,N1
            IF(sym1(I) .NE. "H") THEN
                DO J=1,N2
                        R = mol2(J) - mol1(I)

                        TEMP = LJ_FORMULA(table_DMSO(I,2), table_SOL(J,2), &
                        table_DMSO(I,3), table_SOL(J,3), table_DMSO(I,1), table_SOL(J,1), &
                        R, BOXL, BOXL2)
                        SUML = SUML + TEMP(1)
                        SUMR = SUMR + TEMP(2)
                END DO
            END IF
        END DO
        SUML = SUML * CONV

        EN = SUML + SUMR

END SUBROUTINE CalcLJ_SV

!====================================================================
!====================================================================

FUNCTION LJ_FORMULA(E1, E2, S1, S2, Q1, Q2, RIN, BOXL, BOXL2)

    ! INPUT
    DOUBLE PRECISION, INTENT(IN)    :: E1, E2, S1, S2, Q1, Q2
    DOUBLE PRECISION, INTENT(IN)    :: BOXL, BOXL2
    TYPE (VECTOR), INTENT(IN)       :: RIN

    ! INTERN
    DOUBLE PRECISION, DIMENSION(2)  :: LJ_FORMULA
    DOUBLE PRECISION                :: E, S, RV
    TYPE (VECTOR)                   :: R

    E = sqrt(E1 * E2)
    S = (S1 + S2)/2
    R = RIN

    ! Minimal image convention
    R%X = R%X - BOXL2 * AINT(R%X / BOXL)
    R%Y = R%Y - BOXL2 * AINT(R%Y / BOXL)
    R%Z = R%Z - BOXL2 * AINT(R%Z / BOXL)

    RV = length(R)

    LJ_FORMULA(1) = Q1 * Q2 / RV
    LJ_FORMULA(2) = 4 * E * (S**12 / RV**12 - S**6 / RV**6)

END FUNCTION LJ_FORMULA

!====================================================================
!====================================================================

SUBROUTINE ASSOSIATE_DMSO(SYM, TABLE_SYM, TABLE_Q, TABLE_E, TABLE_S, TABLE)
    ! INPUT
    CHARACTER*4, DIMENSION(:), INTENT(IN)         :: SYM, TABLE_SYM ! Atoomtypes
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: TABLE_Q, TABLE_E, TABLE_S ! params

    ! OUTPUT
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: TABLE

    ! TEMP
    INTEGER :: N, I, A

    N = SIZE(SYM)

    DO I=1,N
        IF (SYM(I) .NE. "H") THEN
            A = findSym(SYM(I), TABLE_SYM)
            TABLE(I,1) = TABLE_Q(A)
            TABLE(I,2) = TABLE_E(A)
            TABLE(I,3) = TABLE_S(A)
        END IF
    END DO

END SUBROUTINE ASSOSIATE_DMSO

!====================================================================
!====================================================================

SUBROUTINE ASSOSIATE_SOLUTE(SYM, TABLE_SYM, TABLE_E, TABLE_S, TABLE)
    ! INPUT
    CHARACTER*4, DIMENSION(:), INTENT(IN)         :: SYM, TABLE_SYM ! Atoomtypes
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: TABLE_E, TABLE_S ! params

    ! OUTPUT
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: TABLE

    ! TEMP
    INTEGER :: N, I, A

    N = SIZE(SYM)

    DO I=1,N
        A = findSym(SYM(I), TABLE_SYM)
        TABLE(I,1) = 0.D0
        TABLE(I,2) = TABLE_E(A)
        TABLE(I,3) = TABLE_S(A)
    END DO

END SUBROUTINE ASSOSIATE_SOLUTE

END MODULE LennardJones
