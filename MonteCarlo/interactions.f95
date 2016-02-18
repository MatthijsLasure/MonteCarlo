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
    implicit none
    contains

!====================================================================

    ! calcLJ: bereken alle interacties via LJ potentiaal
    ! Als j = 0: bereken solvent - solute
    subroutine calcLJ(i, j, en)
        integer, intent(in) :: i,j
        double precision :: en


    end subroutine calcLJ

!====================================================================

    subroutine calcGa(i)
    integer :: i

    end subroutine calcGa

end module interactions
