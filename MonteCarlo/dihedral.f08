MODULE dihedral

    !implicit none
    use vector_class

    contains

    FUNCTION GETDIHEDRAL(MOL, I1, I2, I3, I4) RESULT(HOEK)
        TYPE (VECTOR), DIMENSION(:), ALLOCATABLE :: MOL
        INTEGER :: I1, I2, I3, I4
        DOUBLE PRECISION :: HOEK, TEMP, SIGN
        TYPE (VECTOR) :: V1, V2, V3, V4
        TYPE (VECTOR) :: B1, B2, B3, N1, N2

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
end MODULE dihedral
