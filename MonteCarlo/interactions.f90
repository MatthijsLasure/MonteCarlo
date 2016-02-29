!====================================================================
! interactions.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Bereken energie via LJ of QC
!====================================================================
MODULE interactions
    USE vector_class
    USE lib
    !implicit none

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

    ! Calls in subroutine CALCLJ: 
    ! => getConv (on line <49>)
    SUBROUTINE calcLJ(MOL1, MOL2, SYM1, SYM2, TABLE_SYM, TABLE_Q, TABLE_E, TABLE_S, EN, BOXL, BOXL2)
        ! INPUT
        TYPE (vector), DIMENSION(:), INTENT(IN) :: MOL1, MOL2 ! absolute coords!
        CHARACTER*4, DIMENSION(:), INTENT(IN) :: SYM1, SYM2, TABLE_SYM ! Atoomtypes
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: TABLE_Q, TABLE_E, TABLE_S ! params
        DOUBLE PRECISION :: CONV, BOXL, BOXL2

        ! OUTPUT
        DOUBLE PRECISION:: EN

        ! TEMP
        TYPE (vector):: R
        DOUBLE PRECISION:: RV, E, S ! temp params
        DOUBLE PRECISION:: SUML, SUMR ! sommen
        DOUBLE PRECISION:: TEMPL, TEMPR ! tijdelijke optelling
        INTEGER:: I, J, N1, N2 ! loop vars, totale grootte arrays
        INTEGER:: A, B ! welke atoomsoort

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
                A = findSym(sym1(I), TABLE_SYM)
                DO J=1,N2
                    IF(sym2(J) .NE. "H") THEN
                        B = findSym(sym2(J), TABLE_SYM)
                        E = sqrt(table_e(A) * table_e(B))
                        S = (table_s(A) + table_s(B))/2
                        !R = getDist(mol1(I), mol2(J))
                        R = mol2(J) - mol1(I)

                        ! Minimal image convention
                        R%X = R%X - BOXL2 * ANINT(R%X / BOXL)
                        R%Y = R%Y - BOXL2 * ANINT(R%Y / BOXL)
                        R%Z = R%Z - BOXL2 * ANINT(R%Z / BOXL)

                        RV = length(R)

                        TEMPL = table_Q(A) * table_Q(B) / RV
                        TEMPR = 4 * E * (S**12 / RV**12 - S**6 / RV**6)
                        SUML = SUML + TEMPL
                        SUMR = SUMR + TEMPR
                    END IF
                END DO
            END IF
        END DO
        SUML = SUML * CONV

        EN = SUML + SUMR

    END SUBROUTINE calcLJ

!====================================================================

    SUBROUTINE calcGa(i, j, mol1, mol2, sym1, sym2, en)
        ! INPUT
        TYPE (vector), DIMENSION(:), INTENT(IN) :: MOL1, MOL2 ! absolute coords!
        CHARACTER*4, DIMENSION(:), INTENT(IN) :: SYM1, SYM2 ! Atoomtypes
        INTEGER, INTENT(in) :: I, J

        ! OUTPUT
        DOUBLE PRECISION :: en

        ! INTERNAL VARS
        INTEGER :: K ! loop de loop
        CHARACTER*32 :: gauss_file, gauss_log
        CHARACTER*16 :: str_i, str_j
        INTEGER :: N1, N2

        905 FORMAT(A2, 3F16.8)

        write(str_i, *) I
        write(str_j, *) J

        gauss_file = "input-" // str_i // "-" // str_j // ".com"
        gauss_log = "output-" // str_i // "-" // str_j // ".log"
        N1 = size(mol1)
        N2 = size(mol2)

        ! OPEN FILE
        OPEN(unit=15, file=gauss_file)

        ! Print Mol1
        do K=21,N1
            write (15,905) sym1(K), mol1(K)%X, mol1(K)%Y, mol1(K)%Z
            write (15,905) sym2(K), mol2(K)%X, mol2(K)%Y, mol2(K)%Z
        end do

        CLOSE(15)

    END SUBROUTINE calcGa

END MODULE interactions
