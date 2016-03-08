program BoxScale
    implicit none

    CHARACTER*500 :: input, output, charfactor
    DOUBLE PRECISION :: factor

    DOUBLE PRECISION :: CoM1, CoM2, CoM3
    INTEGER :: nCoM, I
    DOUBLE PRECISION :: BOXL

    CALL GET_COMMAND_ARGUMENT(1, CHARFACTOR)
    CALL GET_COMMAND_ARGUMENT(2, INPUT)
    CALL GET_COMMAND_ARGUMENT(3, OUTPUT)

    read(CHARFACTOR, *) FACTOR

    open(10, file=input)
    open(11, file=output)
    read (10, *) BOXL
    read (10, *) nCoM

    write (11, *) BOXL * factor
    write (11, *) nCoM

    DO I=1,nCoM
        read (10, *) CoM1, CoM2, CoM3
        write (11, *) CoM1 * factor, CoM2 * factor, CoM3 * factor
    END Do

    close(10)
    close(11)
end program BoxScale
