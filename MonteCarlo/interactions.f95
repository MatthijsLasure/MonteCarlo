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
    SUBROUTINE calcLJ(MOL1, MOL2, SYM1, SYM2, TABLE_SYM, TABLE_Q, TABLE_E, TABLE_S, EN)
        ! INPUT
        TYPE (vector), DIMENSION(:), INTENT(IN) :: MOL1 ! absolute coords!
        TYPE (vector), DIMENSION(:), INTENT(IN) :: MOL2 ! absolute coords!
        CHARACTER*4, DIMENSION(:), INTENT(IN) :: SYM1 ! Atoomtypes
        CHARACTER*4, DIMENSION(:), INTENT(IN) :: SYM2 ! Atoomtypes
        CHARACTER*4, DIMENSION(:), INTENT(IN) :: TABLE_SYM ! Atoomtypes
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: TABLE_Q ! params
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: TABLE_E ! params
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: TABLE_S ! params
        DOUBLE PRECISION:: CONVERSION

        ! OUTPUT
        DOUBLE PRECISION:: EN

        ! TEMP
        TYPE (vector):: R
        DOUBLE PRECISION:: RV ! temp params
        DOUBLE PRECISION:: E ! temp params
        DOUBLE PRECISION:: S ! temp params
        DOUBLE PRECISION:: SUML ! sommen
        DOUBLE PRECISION:: SUMR ! sommen
        DOUBLE PRECISION:: TEMPL ! tijdelijke optelling
        DOUBLE PRECISION:: TEMPR ! tijdelijke optelling
        INTEGER:: I ! loop vars, totale grootte arrays
        INTEGER:: J ! loop vars, totale grootte arrays
        INTEGER:: N1 ! loop vars, totale grootte arrays
        INTEGER:: N2 ! loop vars, totale grootte arrays
        INTEGER:: A ! welke atoomsoort
        INTEGER:: B ! welke atoomsoort

        CALL getConv(CONVERSION)

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
                        !R = getDist(mol1(I), mol2(J))
                        R = mol1(I) - mol2(J)

                        if(.true.)then
                        ! Minimal image convention
                        R%X = R%X - BOXL * ANINT(R%X / BOXL)
                        R%Y = R%Y - BOXL * ANINT(R%Y / BOXL)
                        R%Z = R%Z - BOXL * ANINT(R%Z / BOXL)
                        end if
                        RV = length(R)

                        S = (table_s(A) + table_s(B))/2

                        TEMPL = table_Q(A) * table_Q(B) / RV
                        TEMPR = 4 * E * (S**12 / RV**12 - S**6 / RV**6)
                        SUML = SUML + TEMPL
                        SUMR = SUMR + TEMPR
                    END IF
                END DO
            END IF
        END DO
        SUML = SUML * CONVERSION

        EN = SUML + SUMR

    END SUBROUTINE calcLJ

!====================================================================

    SUBROUTINE calcGa(I)
    INTEGER:: I

    END SUBROUTINE calcGa

END MODULE interactions
