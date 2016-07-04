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

MODULE vector_class
    IMPLICIT NONE
    PUBLIC

    TYPE vector
        DOUBLE PRECISION :: x, y, z
    END TYPE vector

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

    CONTAINS

    ! Arimetrische functies
    FUNCTION vecadd(A, B) RESULT(c)
        TYPE (vector), INTENT(IN) :: A
        TYPE (vector), INTENT(IN) :: B
        TYPE (vector):: C
        C%x = A%x + B%x
        C%y = A%y + B%y
        C%z = A%z + B%z
    END FUNCTION vecadd

    FUNCTION vecsub(A, B) RESULT(c)
        TYPE (vector), INTENT(IN) :: A
        TYPE (vector), INTENT(IN) :: B
        TYPE (vector):: C
        C%x = A%x - B%x
        C%y = A%y - B%y
        C%z = A%z - B%z
    END FUNCTION vecsub

    FUNCTION vecmult(A, M) RESULT(c)
        TYPE (vector), INTENT(IN) :: A
        TYPE (vector):: C
        DOUBLE PRECISION, INTENT(IN) :: M
        C%x = M * A%x
        C%y = M * A%y
        C%z = M * A%z
    END FUNCTION vecmult

    FUNCTION vecdiv(A, M) RESULT(c)
        TYPE (vector), INTENT(IN) :: A
        DOUBLE PRECISION, INTENT(IN) :: M
        TYPE (vector):: C
        IF (M /= 0) THEN
            C = vecmult(A, 1.D0 / M)
        ELSE
            C = vecmult(A, 0.D0)
        END IF
    END FUNCTION vecdiv

    FUNCTION cross(A, B) RESULT(c)
        TYPE (vector), INTENT(IN) :: A
        TYPE (vector), INTENT(IN) :: B
        TYPE (vector):: C
        C%x = A%y * B%z - A%z * B%y
        C%y = A%z * B%x - A%x * B%z
        C%z = A%x * B%y - A%y * B%x
    END FUNCTION cross

    FUNCTION dot(A, B) RESULT(d)
        TYPE (vector):: A
        TYPE (vector):: B
        DOUBLE PRECISION:: D
        D = 0.D0
        D = D + A%x * B%x
        D = D + A%y * B%y
        D = D + A%z * B%z
    END FUNCTION dot

    ! Niet arimetrische functies
    FUNCTION length_sq(V) RESULT(r)
        TYPE (vector):: V
        DOUBLE PRECISION:: R
        R = V%x ** 2 + V%y ** 2 + V%z ** 2
    END FUNCTION

    FUNCTION length(V) RESULT(r)
        TYPE (vector):: V
        DOUBLE PRECISION:: R
        R = sqrt(length_sq(V))
    END FUNCTION length

    ! Normalize vector zodat r = 1
    FUNCTION norm(V) RESULT(w)
        TYPE (vector):: V
        TYPE (vector):: W
        DOUBLE PRECISION:: R
        R = length(V)
        W = V / R
    END FUNCTION norm

    FUNCTION setlength(V, A) RESULT(w)
        TYPE (vector):: V
        TYPE (vector):: W
        DOUBLE PRECISION:: A
        W = norm(V)
        W = W * A
    END FUNCTION setlength

    ! bereken afstand tussen 2 vectoren
    FUNCTION getDist(A, B) RESULT(r)
        TYPE (vector), INTENT(IN) :: A
        TYPE (vector), INTENT(IN) :: B
        DOUBLE PRECISION:: R
        R = length(B-A)
    END FUNCTION getDist

    ! bereken afstand tussen 2 vectoren
    FUNCTION getDistSq(A, B) RESULT(r)
        TYPE (vector), INTENT(IN) :: A
        TYPE (vector), INTENT(IN) :: B
        DOUBLE PRECISION:: R
        R = length_sq(B-A)
    END FUNCTION getDistSq

    ! Geef een aray voor matrix multiplication
    FUNCTION getArray(V) RESULT(a)
        TYPE (vector), INTENT(IN) :: V
        DOUBLE PRECISION, DIMENSION(3) :: A
        A(1) = V%x
        A(2) = V%y
        A(3) = V%z
    END FUNCTION getArray

    FUNCTION fromArray(A) RESULT(v)
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: A
        TYPE (vector):: V
        V%x = a(1)
        V%y = a(2)
        V%z = a(3)
    END FUNCTION fromArray

END MODULE vector_class
