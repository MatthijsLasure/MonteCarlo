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
        TYPE (vector) :: AbsPos
    end function RotMatrix
end module lib
