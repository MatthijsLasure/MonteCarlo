program BoxScale
    implicit none

    CHARACTER*500 :: input, output, charfactor
    DOUBLE PRECISION :: factor

    DOUBLE PRECISION :: CoM1, CoM2, CoM3, h1, h2, h3
    INTEGER :: nCoM, I
    DOUBLE PRECISION :: BOXL

    WRITE (*,*) "************"
    WRITE (*,*) "* BoxScale *"
    WRITE (*,*) "************"

    CALL GET_COMMAND_ARGUMENT(1, CHARFACTOR)
    CALL GET_COMMAND_ARGUMENT(2, INPUT)
    CALL GET_COMMAND_ARGUMENT(3, OUTPUT)

    read(CHARFACTOR, *) FACTOR

    open(10, file=input)
    open(11, file=output)
    read (10, *) BOXL
    read (10, *) nCoM
    read (10, "(A,I10.10)") charfactor, I

    write (11, *) BOXL * factor
    write (11, *) nCoM
    write (11, "(A,I10.10)") trim(charfactor), I

    DO I=1,nCoM
        read (10, *) CoM1, CoM2, CoM3, h1, h2, h3
        write (11, *) CoM1 * factor, CoM2 * factor, CoM3 * factor, h1, h2, h3
    END Do

    close(10)
    close(11)

    WRITE (*,*) "BoxScale has finished!"
end program BoxScale
