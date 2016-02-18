!========================================================================================
! Vector_class.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Definiëer een type vector, met bijhorende functies
! vector: double precision x, y, z
! functies: +, -, * (met getal), / (met getal), . (scalair product), * (cross product)
! functies: length, length_sq, normalize (vec_norm), setlength
!
!========================================================================================

module vector_class
    implicit none
    public

    type vector
        double precision :: x, y, z
    end type vector

    ! Overload arimetric operators
    INTERFACE OPERATOR (+)
        MODULE PROCEDURE vecadd
    END INTERFACE

    INTERFACE OPERATOR (-)
        MODULE PROCEDURE vecsub
    END INTERFACE

    INTERFACE OPERATOR (*)
        MODULE PROCEDURE vecmult
    END INTERFACE

    INTERFACE OPERATOR (/)
        MODULE PROCEDURE vecdiv
    END INTERFACE

    INTERFACE OPERATOR (*)
        MODULE PROCEDURE cross
    END INTERFACE

    contains

    ! Arimetrische functies
    function vecadd(a, b) result(c)
        TYPE (vector), INTENT(in) :: a, b
        TYPE (vector) :: c
        c%x = a%x + b%x
        c%y = a%y + b%y
        c%z = a%z + b%z
    end function vecadd

    function vecsub(a, b) result(c)
        TYPE (vector), INTENT(in) :: a, b
        TYPE (vector) :: c
        c%x = a%x - b%x
        c%y = a%y - b%y
        c%z = a%z - b%z
    end function vecsub

    function vecmult(a, m) result(c)
        TYPE (vector), INTENT(in) :: a
        TYPE (vector) :: c
        double precision, INTENT(in) :: m
        c%x = m * a%x
        c%y = m * a%y
        c%z = m * a%z
    end function vecmult

    function vecdiv(a, m) result(c)
        TYPE (vector), INTENT(in) :: a
        double precision, INTENT(in) :: m
        TYPE (vector) :: c
        if (m /= 0) then
            c = vecmult(a, 1.D0 / m)
        else
            c = vecmult(a, 0.D0)
        end if
    end function vecdiv

    function cross(a, b) result(c)
        TYPE (vector), INTENT(in) :: a, b
        TYPE (vector) :: c
        c%x = a%y * b%z - a%z * b%y
        c%y = a%z * b%x - a%x * b%z
        c%z = a%x * b%y - a%y * b%x
    end function cross

    function dot(a, b) result(d)
        TYPE (vector) :: a, b
        double precision d
        d = 0.D0
        d = d + a%x + b%x
        d = d + a%y + b%y
        d = d + a%z + b%z
    end function dot

    ! Niet arimetrische functies
    function length_sq(v) result(r)
        TYPE (vector) :: v
        double precision :: r
        r = v%x ** 2 + v%y ** 2 + v%z ** 2
    end function

    function length(v) result(r)
        TYPE (vector) :: v
        double precision :: r
        r = sqrt(length_sq(v))
    end function length

    ! Normalize vector zodat r = 1
    function norm(v) result(w)
        TYPE (vector) :: v, w
        double precision :: r
        r = length(v)
        w = v / r
    end function norm

    function setlength(v, a) result(w)
        TYPE (vector) :: v, w
        double precision :: a
        w = norm(v)
        w = w * a
    end function setlength

    ! bereken afstand tussen 2 vectoren
    function getDist(a, b) result(r)
        TYPE (vector), INTENT(in) :: a, b
        double precision :: r
        r = length(b-a)
    end function getDist

    ! Geef een aray voor matrix multiplication
    function getArray(v) result(a)
        TYPE (vector), INTENT(in) :: v
        double precision, dimension(3) :: a
        a(1) = v%x
        a(2) = v%y
        a(3) = v%z
    end function getArray

    function fromArray(a) result(v)
        double precision, dimension(3), INTENT(in) :: a
        TYPE (vector) :: v
        v%x = a(1)
        v%y = a(2)
        v%z = a(3)
    end function fromArray

end module vector_class
