!====================================================================
! lib.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Library met nuttige functies
!====================================================================
module lib
    use vector_class

    implicit none

    contains

    ! RotMatrix
    !==========
    ! Roteer een set coördinaten rond de CoM met gegeven hoeken
    ! Hoek: in graden!
    function RotMatrix(CoM, RelPos, hoek) result(AbsPos)
        TYPE (vector), INTENT(in) :: CoM, hoek
        TYPE (vector), DIMENSION(:), INTENT(in) :: RelPos
        TYPE (vector), DIMENSION(size(RelPos)) :: AbsPos

        ! Variabelen
        !===========
        double precision, dimension(3,3) :: RM ! Rotatiematrix
        double precision, dimension(3) :: temp ! Tijdelijke array
        TYPE (vector) :: v ! Tijdelijke vector

        double precision :: cos1, cos2, cos3
        double precision :: sin1, sin2, sin3
        integer :: n ! Aantal atomen
        integer :: i

        n = size(RelPos)

        ! Bereken cos / sin
        !==================
        cos1 = cos(hoek%x)
        cos2 = cos(hoek%y)
        cos3 = cos(hoek%z)
        sin1 = sin(hoek%x)
        sin2 = sin(hoek%y)
        sin3 = sin(hoek%z)

        ! Opstellen matrix
        !=================
        ! Rij 1
        RM(1,1) = cos1 * cos2 - sin1 * sin2 * cos3
        RM(1,2) = sin1 * cos2 + cos1 * sin2 * cos3
        RM(1,3) = sin2 * sin3
        ! Rij 2
        RM(2,1) = -1.D0 * cos1 * sin2 - sin1 * cos2 * cos3
        RM(2,2) = -1.D0 * sin1 * sin2 + cos1 * cos2 * cos3
        RM(2,3) = cos2 * sin3
        ! Rij 3
        RM(3,1) = sin1 * sin3
        RM(3,2) = -1.D0 * cos1 * sin3
        RM(3,3) = cos3

        ! Roteren + transleren atomen
        !============================
        do i=1,n
            temp = getArray(RelPos(i))
            temp = matmul(RM, temp)
            v = fromArray(temp)
            AbsPos(i) = CoM + v
        end do

    end function RotMatrix

    function findSym(type, sym) result(pos)
        character*4 :: type
        integer :: pos, i, n
        character*4, dimension(:) :: sym

        n = size(sym)
        do i=1,n
            if( sym(i) .EQ. type) then
                pos = i
                return
            end if
            pos = 0
        end do
    end function findSym
end module lib
