program BoxScale
    implicit none

    CHARACTER*500 :: input, output, charfactor
    CHARACTER*4   :: sym
    DOUBLE PRECISION :: factor, en

    DOUBLE PRECISION :: CoM1, CoM2, CoM3, h1, h2, h3
    INTEGER :: nCoM, I, stat, I1, I2
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
    ! START BOX
    read (10, *) BOXL
    read (10, *) nCoM
    !read (10, *) en
    read (10, "(ES20.10)", IOSTAT=stat) en
    read (10, "(A,I10.10)") charfactor, I

    write (11, *) BOXL * factor
    write (11, *) nCoM
    !write (11, *) en
    if (stat .eq. 0) write (11, *) en
    write (11, "(A,I10.10)") trim(charfactor), I

    DO I=1,nCoM
        read (10, *) CoM1, CoM2, CoM3, h1, h2, h3
        write (11, *) CoM1 * factor, CoM2 * factor, CoM3 * factor, h1, h2, h3
    END Do

    ! SOLUTE
    !=======
    READ (10, *) ! Line with SOLUTE
    READ (10, *) nCoM ! Lees aantal atomen
    READ (10, "(A)") charfactor ! Comment line


    WRITE (11, *) "SOLUTE"
    WRITE (11, *) nCoM
    WRITE (11, "(A)") trim(charfactor)


    DO I=1, NCOM ! Lees de coördinaten uit
        READ (10,*) SYM, CoM1, CoM2, CoM3
        WRITE (11,*) SYM, CoM1, CoM2, CoM3
    END DO

    READ (10,*) en
    READ (10,*) nCoM

    WRITE (11,*) en
    WRITE (11,*) nCoM

    DO I=1, nCoM
        READ (10,"(I3, 1X, I3, 1X, F10.6)") I1, I2, CoM3
        IF (CoM3 .EQ. 0.D0) THEN
            WRITE (11,"(I3, 1X, I3)") I1, I2
        ELSE
            WRITE (11,"(I3, 1X, I3, 1X, F10.6)") I1, I2, CoM3
        END IF
    END DO

    close(10)
    close(11)

    WRITE (*,*) "BoxScale has finished!"
end program BoxScale
