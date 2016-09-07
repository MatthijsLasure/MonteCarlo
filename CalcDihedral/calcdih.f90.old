program calcdih

    USE vector_class
    USE DIHEDRAL
    USE LIB

    implicit none

    CHARACTER*1000                          :: what_char, sol_file, I1C, I2C, I3C, I4C
    INTEGER                                 :: I1, I2, I3, I4, NSOL, I
    TYPE(vector), DIMENSION(:), ALLOCATABLE :: solute
    CHARACTER*4, DIMENSION(:), ALLOCATABLE  :: sol_sym
    DOUBLE PRECISION                        :: HOEK

    CALL GET_COMMAND_ARGUMENT(1, what_char)
    CALL GET_COMMAND_ARGUMENT(2, sol_file)
    CALL GET_COMMAND_ARGUMENT(3, I1C)
    CALL GET_COMMAND_ARGUMENT(4, I2C)
    CALL GET_COMMAND_ARGUMENT(5, I3C)
    CALL GET_COMMAND_ARGUMENT(6, I4C)

    READ (I1C, "(I4.4)") I1
    READ (I2C, "(I4.4)") I2
    READ (I3C, "(I4.4)") I3
    READ (I4C, "(I4.4)") I4

    OPEN (UNIT=10, FILE=trim(sol_file))

    ! Wade through the box stuff
    READ (*,*) ! BOXL
    READ (*,*) NSOL
    READ (*,*) ! En
    READ (*,*) ! BOXNAME
    DO I=1,NSOL
        READ (*,*) ! CoM, HOEK
    END DO
    READ (*,*) ! SOLUTE

    ! solute.txt: conformatie opgeloste molecule (sol)

    READ (10, *) nSol ! Lees aantal atomen
    READ (10, *) ! Comment
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(sol_sym(NSOL))
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (10,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    CLOSE(10)

    HOEK = GETDIHEDRAL(solute, I1, I2, I3, I4) * 180 / PI
    WRITE (*,"(F15.10)") HOEK
end program calcdih
