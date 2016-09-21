program changeDih
    USE vector_class
    USE lib
    USE dihedral
    implicit none

    ! Variables
    CHARACTER*1000      :: BOX_FILE, OUT_FILE, SOL_FILE, DMSO_FILE, DUMP_FILE, TEMP, THISFLAG
    CHARACTER*500       :: ARG1, ARG2, ARG3, ARG4, BOXID, SOLID
    LOGICAL             :: SILENT = .FALSE., USE_IN = .TRUE., USE_OUT = .TRUE., USE_DMSO = .FALSE.
    LOGICAL             :: BOX = .FALSE.
    LOGICAL             :: RANKING = .FALSE.

    LOGICAL             :: DOROTATION = .FALSE.
    INTEGER             :: I, ISTATUS
    INTEGER             :: MODE, IOout = 11, IOin = 10
    INTEGER             :: A1, A2, A3, A4
    DOUBLE PRECISION    :: THISHOEK, HOEKROT, HOEK_POST

    ! DATA
    TYPE (VECTOR), DIMENSION(:), ALLOCATABLE :: COM, HOEK, SOLUTE, DMSO, MOL1, SOLROT
    CHARACTER*4, DIMENSION(:), ALLOCATABLE :: SOL_SYM, DMSO_SYM
    ! Array met bindingen die gebruikt mogen worden voor de dihedrale rotatie
    INTEGER, DIMENSION(:,:), ALLOCATABLE        :: DIHOEK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DROTSOLU_ARRAY
    INTEGER :: NDIHOEK
    DOUBLE PRECISION :: BOXL, EN, E_SOL
    INTEGER :: NCOM, NDMSO, NSOL
    CALL ARGLOAD()

    IF ( .NOT. SILENT) THEN
        WRITE (*,*) "#####################"
        WRITE (*,*) "# CalcDihedral v2.0 #"
        WRITE (*,*) "#-------------------#"
        WRITE (*,*) "#  Matthijs Lasure  #"
        WRITE (*,*) "#####################"
    END IF


    CALL read_norm()

    THISHOEK = GETDIHEDRAL(solute, A1, A2, A3, A4) * 180 / PI
    IF (SILENT) THEN
        WRITE (*,"(F10.5)") THISHOEK
    ELSE
        501 FORMAT("Dihedral angle of ", I3.3, " - ", I3.3, " - ", I3.3, &
        " - ", I3.3, " is ", F15.10, " degrees.")
        WRITE (*,501) A1, A2, A3, A4, THISHOEK
    END IF

    IF (DOROTATION) THEN
        SOLROT = SETDIHEDRAL(SOLUTE, A1, A2, A3, A4, HOEKROT * PI / 180)
        HOEK_POST = GETDIHEDRAL(SOLROT, A1, A2, A3, A4) * 180 / PI

        502 FORMAT("Rotated from ", F7.2, " degrees with ", F7.2, " degrees to ", F7.2, " degrees.")
        IF (.NOT. SILENT) WRITE (*,502) THISHOEK, HOEKROT, HOEK_POST
        IF ((HOEK_POST - THISHOEK - HOEKROT) .LT. 0.01) THEN
            IF (.NOT. SILENT) WRITE (*,*) "Check is OK, good rotation (error < 0.01)."
        ELSE
            IF (.NOT. SILENT) THEN
                WRITE (*,*) "Check is not OK! Something went wrong (error >= 0.01) :'(",HOEK_POST - THISHOEK + HOEKROT
            ELSE
                WRITE (0,*) "ERROR"
            END IF
        END IF
        CALL WRITE_NORM()
    END IF


CONTAINS
SUBROUTINE ARGLOAD()
    INTEGER :: SKIP = 0
    INTEGER :: FLAG, NUMFLAGS

    NUMFLAGS = COMMAND_ARGUMENT_COUNT()
    IF ( NUMFLAGS .LE. 1) THEN
        CALL HELP()
        STOP
    END IF

    flags: DO FLAG=1,NUMFLAGS
        CALL get_command_argument(FLAG, THISFLAG)
        IF (SKIP .GT. 0) THEN
            SKIP = SKIP - 1
            CYCLE
        END IF

        SELECT CASE (THISFLAG)
            CASE ("-r", "--rotate")
                DOROTATION = .TRUE.
                CALL get_command_argument(FLAG+1, ARG1)
                READ(ARG1, *) HOEKROT
                SKIP = 1

            CASE ("-b", "--box")
                CALL get_command_argument(FLAG+1, BOX_FILE)
                BOX = .TRUE.
                USE_IN = .FALSE.
                SKIP = 1
            CASE ("-o", "--output")
                CALL get_command_argument(FLAG+1, OUT_FILE)
                USE_OUT = .FALSE.
                SKIP = 1
            CASE ("-q", "--quiet")
                SILENT = .TRUE.
            CASE ("-n", "--atoms")
                CALL get_command_argument(FLAG+1, ARG1)
                CALL get_command_argument(FLAG+2, ARG2)
                CALL get_command_argument(FLAG+3, ARG3)
                CALL get_command_argument(FLAG+4, ARG4)

                READ (ARG1, "(I3.3)") A1
                READ (ARG2, "(I3.3)") A2
                READ (ARG3, "(I3.3)") A3
                READ (ARG4, "(I3.3)") A4

                SKIP = 4

            CASE DEFAULT
                WRITE (*,*) "Incorrect argument: ", THISFLAG
                CALL HELP()
                CALL EXIT(1)
        END SELECT
    END DO flags
END SUBROUTINE ARGLOAD


    SUBROUTINE read_norm

    INTEGER :: IOWORK = 10, ISTAT = 0

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    IF (USE_IN) THEN
        IOWORK = 5
    ELSE
        OPEN (UNIT=IOWORK, FILE=BOX_FILE, ACTION='READ', STATUS='OLD', IOSTAT=istat)
        IF (ISTAT .NE. 0) THEN
            WRITE (*,*) "Error opening box in ", TRIM(BOX_FILE)
            CALL EXIT(2)
        END IF
    END IF

    READ (IOWORK, *) BOXL ! Box grootte
    READ (IOWORK, *) nCoM ! Lees aantal moleculen
    READ (IOWORK,"(E20.10)", IOSTAT=istat) EN ! Energy from previous run

    READ (IOWORK, *) BOXID ! Box ID
    ALLOCATE(CoM(NCOM))
    ALLOCATE(hoek(NCOM))
    DO I= 1, NCOM
        READ(IOWORK,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO

    ! solute.txt: conformatie opgeloste molecule (sol)
    ! NOW INCORPERATED IN BOX.TXT
    READ (IOWORK, *) ! Line with SOLUTE
    READ (IOWORK, *) nSol ! Lees aantal atomen
    READ (IOWORK, "(A)") SOLID! Comment line
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(sol_sym(NSOL))
    DO I=1, NSOL ! Lees de coÃ¶rdinaten uit
        READ (IOWORK,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO

    CLOSE (IOWORK)


END SUBROUTINE read_norm

SUBROUTINE WRITE_NORM()
    INTEGER :: IOwork = 10, ISTAT = 0
    WRITE (*,*) "Writing data...", TRIM(OUT_FILE)

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=IOwork, FILE=OUT_FILE, IOSTAT=istat)
    IF ( ISTAT .NE. 0) THEN
        WRITE (0,"(A, I6.6)") "Error opening output file! #", ISTAT
        STOP -1
    END IF
    WRITE (IOwork, *) BOXL ! Box grootte
    WRITE (IOwork, *) NCOM ! Lees aantal moleculen
    WRITE (IOwork, *) EN
    WRITE (IOwork, *) BOXID
    DO I= 1, NCOM
        WRITE(IOwork,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO
    WRITE (IOwork, *) "SOLUTE"
    ! solute.txt: conformatie opgeloste molecule (sol)
    WRITE (IOwork, *) NSOL ! Schrijf aantal atomen
    WRITE (IOwork, *) trim(SOLID) ! Comment line
    DO I=1, NSOL ! Schrijf de coördinaten uit
        WRITE (IOwork,*) sol_sym(I), SOLROT(I)%X, SOLROT(I)%Y, SOLROT(I)%Z
    END DO
    WRITE (IOwork,*) E_SOL / HARTREE2KJMOL ! Schrijf energie in HARTREE
    WRITE (IOwork,*) NDIHOEK ! Schrijf aantal hoeken
    DO I=1, NDIHOEK
        IF ( DROTSOLU_ARRAY(I) .EQ. -1.D0) THEN
            WRITE (IOwork,"(I3, 1X, I3)") DIHOEK(I,1), DIHOEK(I,2)
        ELSE
            WRITE (IOwork,"(I3, 1X, I3, 1X, F10.6)") DIHOEK(I,1), DIHOEK(I,2), DROTSOLU_ARRAY(I)
        END IF
    END DO
    CLOSE(IOwork)
END SUBROUTINE WRITE_NORM

end program changeDih
