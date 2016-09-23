program boxbin
    USE VECTOR_CLASS
    USE LIB

    implicit none

    CHARACTER*255 :: DMSO_FILE, BOX_FILE, OUT_FILE, SOLNAME, BOX_ID
    TYPE (VECTOR), DIMENSION(:), ALLOCATABLE :: DMSO, COM, HOEK, SOLUTE
    CHARACTER*4, DIMENSION(:), ALLOCATABLE :: DMSO_SYM, SOL_SYM
    INTEGER :: NDMSO, NCOM, NSOL, ISTAT, IOWORK=10, I
    DOUBLE PRECISION :: E_DMSO, ENG, E_SOL, BOXL
    LOGICAL :: INBIN = .FALSE., DOPRINT = .TRUE.

    ! Array met bindingen die gebruikt mogen worden voor de dihedrale rotatie
    INTEGER, DIMENSION(:,:), ALLOCATABLE        :: DIHOEK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DROTSOLU_ARRAY
    INTEGER :: NDIHOEK

    DMSO_FILE = "na"
    OUT_FILE = "no"

    CALL GETARGS

    CALL READ_DMSO

    IF(INBIN) THEN
        CALL READ_BIN
    ELSE
        CALL READ_NORM
    END IF

    IF (DOPRINT) THEN
        CALL PRINTIT
    END IF

    IF (OUT_FILE .NE. "no") THEN
        CALL WRITE_BIN
    END IF

CONTAINS

SUBROUTINE GETARGS
    INTEGER :: NUMARGS, SKIP, I
    CHARACTER*255 :: ARG

    SKIP = 0
    NUMARGS = command_argument_count()

    IF (NUMARGS .EQ. 0) THEN
        CALL HELP
        STOP 1
    END IF

    DO I=1, NUMARGS
        IF (SKIP .GT. 0) THEN
            SKIP = SKIP - 1
            CYCLE
        END IF
        CALL get_command_argument(I, ARG)

        SELECT CASE (ARG)
            CASE ("-b", "--binary")
                INBIN = .TRUE.
                SKIP = 1
                CALL get_command_argument(I+1, BOX_FILE)
            CASE ("-t", "--text")
                INBIN = .FALSE.
                SKIP = 1
                CALL get_command_argument(I+1, BOX_FILE)
            CASE ("-o", "--output")
                SKIP = 1
                CALL get_command_argument(I+1, OUT_FILE)
            CASE ("-d", "--DMSO")
                SKIP = 1
                CALL get_command_argument(I+1, DMSO_FILE)
            CASE ("-q", "--silent")
                DOPRINT = .FALSE.
        END SELECT


    END DO



END SUBROUTINE GETARGS

SUBROUTINE READ_DMSO
! DMSO.txt: conformatie DMSO
    OPEN (UNIT=IOwork, FILE=DMSO_FILE, ACTION='READ', STATUS='OLD', IOSTAT=istat)
    IF (istat .NE. 0) THEN
        WRITE (*,*) "Error opening DMSO in ", TRIM(DMSO_FILE), " Code ", istat
        WRITE (*,*) "Using internal coordinates."

        ALLOCATE(DMSO(10))
        ALLOCATE(DMSO_SYM(10))
        CALL DMSO_init(nDMSO, DMSO_sym, DMSO, E_DMSO)
        !STOP
    ELSE
        READ (IOwork, *) nDMSO ! Lees aantal atomen
        READ (IOwork, *) ! Comment line
        ALLOCATE(DMSO(NDMSO)) ! Ken correcte groottes toe aan de arrays
        ALLOCATE(DMSO_sym(NDMSO))
        DO I=1, NDMSO
            READ (IOwork,*) DMSO_sym(I), DMSO(I)%X, DMSO(I)%Y, DMSO(I)%Z
        END DO
        READ (IOwork,*) E_DMSO
        CLOSE(IOwork)
    END IF
END SUBROUTINE READ_DMSO

! Read input files as normal files
SUBROUTINE read_norm

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=IOwork, FILE=BOX_FILE, ACTION='READ', STATUS='OLD', IOSTAT=istat)
    IF (istat .NE. 0) THEN
        WRITE (*,*) "READ_NORM: Error opening box in ", TRIM(BOX_FILE)
        STOP
    END IF
    READ (IOwork, *) BOXL ! Box grootte
    READ (IOwork, *) nCoM ! Lees aantal moleculen
    READ (IOwork,"(E20.10)", IOSTAT=istat) ENG ! Energy from previous run
    IF (istat .NE. 0) THEN
        WRITE (*,*) "No energy detected in input box. Will be calculated."
        WRITE (*,*) istat
        ENG = 123456.789
        BACKSPACE(IOwork)
    END IF

    READ (IOwork, *) BOX_ID ! Box ID
    ALLOCATE(CoM(NCOM))
    ALLOCATE(hoek(NCOM))
    DO I= 1, NCOM
        READ(IOwork,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO

    ! solute.txt: conformatie opgeloste molecule (sol)
    ! NOW INCORPERATED IN BOX.TXT
    READ (IOwork, *) ! Line with SOLUTE
    READ (IOwork, *) nSol ! Lees aantal atomen
    READ (IOwork, "(A)") SOLNAME ! Comment line
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(sol_sym(NSOL))
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (IOwork,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    READ (IOwork,*) E_sol
    READ (IOwork,*) NDIHOEK
    ALLOCATE(DIHOEK(NDIHOEK, 2))
    ALLOCATE(DROTSOLU_ARRAY(NDIHOEK))
    DO I=1, NDIHOEK
        READ (IOwork,"(I3, 1X, I3, 1X, F10.6)") DIHOEK(I,1), DIHOEK(I,2), DROTSOLU_ARRAY(I)
        IF (DROTSOLU_ARRAY(I) .EQ. 0.D0) DROTSOLU_ARRAY(I) = -1.D0
    END DO

    CLOSE (IOwork)


END SUBROUTINE read_norm

! Read input files as normal files
SUBROUTINE read_bin
    CHARACTER*6 :: SOLSOL = "SOLUTE"
    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=IOwork, FILE=BOX_FILE, ACTION='READ', STATUS='OLD', IOSTAT=istat, FORM='UNFORMATTED')
    IF (istat .NE. 0) THEN
        WRITE (*,*) "READ_BIN: Error opening box in ", TRIM(BOX_FILE)
        STOP
    END IF
    READ (IOwork) BOXL ! Box grootte
    READ (IOwork) nCoM ! Lees aantal moleculen
    READ (IOwork) ENG ! Energy from previous run

    READ (IOwork) BOX_ID     ! Box ID
    ALLOCATE(CoM(NCOM))
    ALLOCATE(hoek(NCOM))
    DO I= 1, NCOM
        READ(IOwork) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO

    ! solute.txt: conformatie opgeloste molecule (sol)
    ! NOW INCORPERATED IN BOX.TXT
    READ (IOwork) SOLSOL ! Line with SOLUTE
    READ (IOwork) nSol ! Lees aantal atomen
    READ (IOwork) SOLNAME ! Comment line
    ALLOCATE(SOLUTE(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(sol_sym(NSOL))
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (IOwork) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    READ (IOwork) E_sol
    READ (IOwork) NDIHOEK
    ALLOCATE(DIHOEK(NDIHOEK, 2))
    ALLOCATE(DROTSOLU_ARRAY(NDIHOEK))
    DO I=1, NDIHOEK
        READ (IOwork) DIHOEK(I,1), DIHOEK(I,2), DROTSOLU_ARRAY(I)
        IF (DROTSOLU_ARRAY(I) .EQ. 0.D0) DROTSOLU_ARRAY(I) = -1.D0
    END DO

    CLOSE (IOwork)


END SUBROUTINE read_bin

SUBROUTINE WRITE_BIN

    CHARACTER*6 :: SOLSOL = "SOLUTE"
    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=IOwork, FILE=OUT_FILE, ACTION='WRITE', IOSTAT=istat, FORM='UNFORMATTED')
    IF (istat .NE. 0) THEN
        WRITE (*,*) "WRITE_BIN: Error opening box in ", TRIM(OUT_FILE)
        STOP
    END IF
    WRITE (IOwork) BOXL ! Box grootte
    WRITE (IOwork) nCoM ! Lees aantal moleculen
    WRITE (IOwork) ENG ! Energy from previous run

    WRITE (IOwork) BOX_ID     ! Box ID
    DO I= 1, NCOM
        WRITE(IOwork) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO

    ! solute.txt: conformatie opgeloste molecule (sol)
    ! NOW INCORPERATED IN BOX.TXT
    WRITE (IOwork) SOLSOL ! Line with SOLUTE
    WRITE (IOwork) nSol ! Lees aantal atomen
    WRITE (IOwork) SOLNAME ! Comment line
    DO I=1, NSOL ! Lees de coördinaten uit
        WRITE (IOwork) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    WRITE (IOwork) E_sol
    WRITE (IOwork) NDIHOEK
    DO I=1, NDIHOEK
        WRITE (IOwork) DIHOEK(I,1), DIHOEK(I,2), DROTSOLU_ARRAY(I)
    END DO

    CLOSE (IOwork)

END SUBROUTINE WRITE_BIN

SUBROUTINE PRINTIT

    900 FORMAT (A, I6.6)
    901 FORMAT (A, F8.2)
    902 FORMAT (6F8.4)
    903 FORMAT (A, 3F8.4)
    904 FORMAT (2I3.3,F8.2)

    WRITE (*,901) "BoxL: ", BOXL
    WRITE (*,900) "nCoM: ", NCOM
    WRITE (*,901) "ENG:  ", ENG
    WRITE (*, "(A,A)") "BOX ID: ", BOX_ID
    WRITE (*,"(A)") "Coordinates + angles:"
    DO I=1,NCOM
        WRITE (*,902) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO
    WRITE (*,"(A)") "End coordinates"

    WRITE (*,900) "nSol: ", NSOL
    WRITE (*,"(A,A)") "SOL ID: ", SOLNAME
    DO I=1,NSOL
        WRITE (*,903) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    WRITE (*,901) "Esol: ", E_SOL
    WRITE (*,900) "Dihoek: ", NDIHOEK
    DO I=1,NDIHOEK
        WRITE (*,904) DIHOEK(I,1), DIHOEK(I,2), DROTSOLU_ARRAY(I)
    END DO

END SUBROUTINE PRINTIT

end program boxbin
