!====================================================================
! interactions.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Bereken energie via LJ of QC
!====================================================================
module interactions
    use vector_class
    use lib
    !implicit none
    contains

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

    subroutine calcLJ(mol1, mol2, sym1, sym2, table_sym, table_Q, table_e, table_s, en)
        ! INPUT
        TYPE (vector), dimension(:), INTENT(in) :: mol1, mol2 ! absolute coords!
        character*4, dimension(:), INTENT(in) :: sym1, sym2, table_sym ! Atoomtypes
        double precision, dimension(:), INTENT(in) :: table_Q, table_e, table_s ! params

        ! OUTPUT
        double precision :: en

        ! TEMP
        double precision :: r, e, s ! temp params
        double precision :: sumL, sumR ! sommen
        integer :: i, j, n1, n2 ! loop vars, totale grootte arrays
        integer :: a, b ! welke atoomsoort

        n1 = size(mol1)
        n2 = size(mol2)

        sumL = 0.D0
        sumR = 0.D0

        do i=1,n1
            a = findSym(sym1(i), table_sym)
            do j=i+1,n2
            b = findSym(sym2(j), table_sym)
                e = sqrt(table_e(a) * table_e(b))
                r = getDist(mol1(i), mol2(i))
                s = (table_s(a) + table_s(b))/2

                sumL = sumL + table_Q(a) * table_Q(b) / r
                sumR = sumR + 4 * e * (s**12 / r**12 - s**6 / r**6)
            end do
        end do
        sumL = sumL / (4 * pi * epsilon0)

        en = sumL + sumR

    end subroutine calcLJ

!====================================================================

    subroutine calcGa(i)
    integer :: i

    end subroutine calcGa

end module interactions
