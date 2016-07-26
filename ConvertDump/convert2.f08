PROGRAM convert2
    USE lib
    USE sortLib
    IMPLICIT NONE

    ! Variables
    CHARACTER*1000      :: BOX_FILE, OUT_FILE, SOL_FILE, DMSO_FILE, DUMP_FILE, TEMP, THISFLAG
    CHARACTER*500       :: ARG1, ARG2, BOXID, SOLID
    LOGICAL             :: SILENT = .FALSE., USE_IN = .TRUE., USE_OUT = .TRUE., USE_DMSO = .FALSE.
    LOGICAL             :: BOX = .FALSE., SOL = .FALSE., DUMP = .FALSE., doH = .FALSE., CORNERS = .TRUE.
    LOGICAL             :: RANKING = .FALSE.
    INTEGER             :: I, ISTATUS
    INTEGER             :: MODE, IOout = 11, IOin = 10
    INTEGER             :: ATOM1, ATOM2

    ! DATA
    TYPE (VECTOR), DIMENSION(:), ALLOCATABLE :: COM, HOEK, SOLUTE, DMSO, MOL1
    CHARACTER*4, DIMENSION(:), ALLOCATABLE :: SOL_SYM, DMSO_SYM
    DOUBLE PRECISION :: BOXL, EN
    INTEGER :: NCOM, NDMSO, NSOL

    ! Process args
    CALL ARGLOAD()

    IF (USE_IN) IOin = 5
    IF (USE_OUT) THEN
        SILENT = .TRUE.
        IOout = 6
    END IF

    IF ( .NOT. SILENT) THEN
        WRITE (*,*) "+------------------+"
        WRITE (*,*) "| Converter for MC |"
        WRITE (*,*) "+------------------+"
    END IF

    ! Load data
    CALL read_norm()

    SELECT CASE (mode)
        ! DUMP
        CASE (0)
            CALL doDump()
        ! BOX
        CASE (1)
            CALL doBox(CoM, HOEK)
        ! Dump 2 dist
        CASE (2)
            CALL doDist()
        ! BOX SORT
        CASE (3)
            CALL boxDist()
    END SELECT

    CONTAINS
SUBROUTINE ARGLOAD()
    INTEGER :: SKIP = 0
    INTEGER :: FLAG, NUMFLAGS

    NUMFLAGS = COMMAND_ARGUMENT_COUNT()
    IF ( NUMFLAGS .EQ. 0) THEN
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
            CASE ("-m", "--mode")
                CALL get_command_argument(FLAG+1, ARG1)
                READ (ARG1, "(I2.2)") MODE
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
            CASE ("-d", "--dmso")
                CALL get_command_argument(FLAG+1, DMSO_FILE)
                USE_DMSO = .TRUE.
                SKIP = 1
            CASE ("-s", "--solute")
                CALL get_command_argument(FLAG+1, SOL_FILE)
                SOL = .TRUE.
                SKIP = 1
            CASE ("-D", "--dump")
                CALL get_command_argument(FLAG+1, DUMP_FILE)
                DUMP = .TRUE.
                USE_IN = .FALSE.
                SKIP = 1
            CASE ("-q", "--quiet")
                SILENT = .TRUE.
            CASE ("-n", "--atoms")
                CALL get_command_argument(FLAG+1, ARG1)
                CALL get_command_argument(FLAG+2, ARG2)

                READ (ARG1, "(I3.3)") ATOM1
                READ (ARG2, "(I3.3)") ATOM2
                SKIP = 2
            CASE ("-e", "--empty")
                CORNERS = .FALSE.
            CASE ("-h", "--hydrogen")
                doH = .TRUE.
            CASE ("-r", "--ranking")
                ranking = .TRUE.

            CASE DEFAULT
                WRITE (*,*) "Incorrect argument: ", THISFLAG
                CALL HELP()
                CALL EXIT(1)
        END SELECT
    END DO flags
END SUBROUTINE ARGLOAD

!========================================

SUBROUTINE read_norm

    INTEGER :: IOWORK = 10, ISTAT = 0
    DOUBLE PRECISION :: E_DMSO ! Delete afterwards

! DMSO.txt: conformatie DMSO
    IF (USE_DMSO) THEN
        OPEN (UNIT=IOWORK, FILE=DMSO_FILE, ACTION='READ', STATUS='OLD', IOSTAT=istat)
        IF (ISTAT .NE. 0) THEN
            IF (.NOT. SILENT) WRITE (*,*) "Error opening DMSO in ", TRIM(DMSO_FILE), ": code ", ISTAT
            IF (.NOT. SILENT) WRITE (*,*) "Using internal coordinates."

            ALLOCATE(DMSO(10))
            ALLOCATE(DMSO_SYM(10))
            ALLOCATE(MOL1(10))
            CALL DMSO_init(NDMSO, DMSO_SYM, DMSO, E_DMSO)
        ELSE
            READ (IOWORK, *) nDMSO ! Lees aantal atomen
            READ (IOWORK, *) ! Comment line
            ALLOCATE(DMSO(NDMSO)) ! Ken correcte groottes toe aan de arrays
            ALLOCATE(DMSO_sym(NDMSO))
            ALLOCATE(MOL1(NDMSO))
            DO I=1, NDMSO
                READ (IOWORK,*) DMSO_sym(I), DMSO(I)%X, DMSO(I)%Y, DMSO(I)%Z
            END DO
            READ (IOWORK,*) E_DMSO
            CLOSE(IOWORK)
        END IF
    ELSE
            IF (.NOT. SILENT) WRITE (*,*) "Using internal coordinates for DMSO."
            ALLOCATE(DMSO(10))
            ALLOCATE(DMSO_SYM(10))
            ALLOCATE(MOL1(10))
            CALL DMSO_init(NDMSO, DMSO_SYM, DMSO, E_DMSO)
    END IF

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    IF (BOX) THEN
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
        IOWORK = 10
    END IF

    IF (SOL) THEN
        ! solute.txt: conformatie opgeloste molecule (sol)
        OPEN (UNIT=10, FILE=TRIM(SOL_FILE), IOSTAT=istat)
        IF (ISTAT .NE. 0) THEN
            WRITE (*,*) "Error opening solute in ", TRIM(SOL_FILE)
            CALL EXIT(3)
        END IF
        READ (10, *) nSol ! Lees aantal atomen
        READ (10, *) SOLID! Comment
        ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
        ALLOCATE(sol_sym(NSOL))
        DO I=1, NSOL ! Lees de coördinaten uit
            READ (10,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
        END DO
    CLOSE(10)
    END IF

END SUBROUTINE read_norm

SUBROUTINE doDump()
    INTEGER :: NHYDROGEN, I, J, K, N, TIMESTEP
    TYPE (VECTOR) :: CoM, hoek

    IF (.NOT. SILENT) WRITE (*,*) "Processing dump. This may take a while..."

    NHYDROGEN = 0
    IF (doH) THEN
        DO I=1, NDMSO
            IF (DMSO_SYM(I) .EQ. "H") NHYDROGEN = NHYDROGEN + 1
        END DO
    END IF

    IF (.NOT. USE_IN)  OPEN (IOin, FILE=TRIM(DUMP_FILE))
    IF (.NOT. USE_OUT) OPEN (IOout, FILE=TRIM(OUT_FILE))

    READ (IOin, *) BOXL

    loop: DO
        READ (IOin, *, IOSTAT=ISTATUS) NCOM
        IF (ISTATUS < 0) exit
        READ (IOin, "(E20.10)", IOSTAT=ISTATUS) En
        IF (ISTATUS .NE. 0) THEN
            En = 0.D0
            BACKSPACE(10)
        END IF
        READ (10, "(A10,I20.20)") TEMP, TIMESTEP

        IF (CORNERS) THEN
            N = NCOM * (NDMSO - NHYDROGEN) + NSOL + 8
        ELSE
            N = NCOM * (NDMSO - NHYDROGEN) + NSOL
        END IF

        WRITE (IOout,*) N
        WRITE (IOout, "('ID ', I20.20, ' Energy: ', F20.16, ' kJ/mol')") TIMESTEP, En

        IF (CORNERS) THEN
            DO I=0,1
                DO J=0,1
                    DO K=0,1
                        WRITE (IOout,*) "Ne", (-1.D0)**I * BOXL, (-1.D0)**J * BOXL, (-1.D0)**K * BOXL
                    END DO
                END DO
            END DO
        END IF

        SOL_loop: DO I=1,NSOL
            WRITE (IOout,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
        END DO SOL_loop
        CoM_loop: DO I=1,NCOM
            READ (IOin, *) CoM%X, CoM%Y, CoM%Z, hoek%X, hoek%Y, hoek%Z
            MOL1 = RotMatrix(CoM, DMSO, hoek)
            DO J=1,NDMSO
                IF (DMSO_SYM(J) .NE. "H" .OR. doH) THEN
                    WRITE (IOout, *) DMSO_sym(J), MOL1(J)%X, MOL1(J)%Y, MOL1(J)%Z
                END IF
            END DO
        END DO CoM_loop
    END DO loop

    IF (.NOT. USE_IN)  CLOSE(IOin)
    IF (.NOT. USE_OUT) CLOSE(IOout)

END SUBROUTINE doDump


SUBROUTINE doDist()

    TYPE (VECTOR) :: vecA, vecB
    INTEGER :: I, IMIN, TIMESTEP
    DOUBLE PRECISION :: R, RMIN
    TYPE (VECTOR) :: CoM, hoek

    IF (.NOT. USE_IN)  OPEN (IOin, FILE=TRIM(DUMP_FILE))
    IF (.NOT. USE_OUT) OPEN (IOout, FILE=TRIM(OUT_FILE))

    READ (IOin, *) BOXL

    vecA = SOLUTE(ATOM1)
    DO
        READ (10, *, IOSTAT=ISTATUS) NCOM
        IF (ISTATUS < 0) exit
        READ (10, "(A10,I20.20)") TEMP, TIMESTEP
        RMIN = HUGE(R)
        IMIN = HUGE(I)
        DO I=1,NCOM
            !READ (10, "(6F24.16)") CoM%X, CoM%Y, CoM%Z, hoek%X, hoek%Y, hoek%Z
            READ (10,"(A)") TEMP
            !WRITE (*,*) trim(TEMP)
            READ (TEMP, *) CoM%X, CoM%Y, CoM%Z, hoek%X, hoek%Y, hoek%Z
            MOL1 = RotMatrix(CoM, DMSO, hoek)
            vecB = MOL1(ATOM2)
            R = getDist(vecA, vecB)
            IF (R .LT. RMIN) THEN
                RMIN = R
                IMIN = I
            END IF
        END DO
        WRITE (11, *) TIMESTEP, IMIN, RMIN
    END DO

    IF (.NOT. USE_IN)  CLOSE(IOin)
    IF (.NOT. USE_OUT) CLOSE(IOout)

END SUBROUTINE doDist

SUBROUTINE boxDist()

    DOUBLE PRECISION, DIMENSION(NCOM) :: DIST
    TYPE (VECTOR), DIMENSION(NCOM) :: COMS, HOEKS
    INTEGER, DIMENSION(NCOM) :: IX
    DOUBLE PRECISION :: R, RMIN
    INTEGER :: I, A, B, AMIN, BMIN
    TYPE (VECTOR) :: vecA, vecB

    ! Print header if necessary
    IF (.NOT. (SILENT .OR. RANKING)) WRITE (*,*) "Solv SoluNr SolvNr r"

    ! Init indiches
    DO I = 1,NCOM
        IX(I) = I
    END DO

    ! Get distances
    box: DO I=1, NCOM
        MOL1 = RotMatrix(CoM(I), DMSO, HOEK(I))
        RMIN = HUGE(RMIN)
        BMIN = HUGE(BMIN)
        AMIN = HUGE(AMIN)

        molecule: DO A=1,NSOL
            vecA = SOLUTE(A)
            DO B=1,nDMSO
                vecB = MOL1(B)
                R = getDistSq(vecA, vecB)
                IF (R .LT. RMIN) THEN
                    RMIN = R
                    BMIN = B
                    AMIN = A
                END IF
            END DO
        END DO molecule

        RMIN = SQRT(RMIN) ! Because it was squared before (save cpu)
        DIST(I) = RMIN

        IF (.NOT. (SILENT .OR. RANKING)) WRITE (*,*) I, AMIN, BMIN, RMIN
    END DO box

    ! Sort
    CALL SORT(DIST, IX)

    ! Sort the CoM & HOEK
    IF (RANKING .AND. .NOT. SILENT) WRITE (*,"(A)") "+--------+---------+"
    IF (RANKING .AND. .NOT. SILENT) WRITE (*,"(A)") "|  Place | Solvent |"
    IF (RANKING .AND. .NOT. SILENT) WRITE (*,"(A)") "+--------+---------+"
    DO I=1,NCOM
        COMS(I) = COM(IX(I))
        HOEKS(I) = HOEKS(IX(I))
        IF (RANKING .AND. .NOT. SILENT) WRITE (*,"('| ', I6, ' |  ', I6, ' |')") I, IX(I)
    END DO
    IF (RANKING .AND. .NOT. SILENT) WRITE (*,"(A)") "+--------+---------+"

    ! PRINT
    CALL doBox(CoMS, HOEKS)

END SUBROUTINE boxDist

SUBROUTINE doBox(CoM, HOEK)
    TYPE (VECTOR), DIMENSION(NCOM), INTENT(IN) :: CoM, HOEK
    INTEGER :: I, J, K, N, NHYDROGEN

    ! Printing time!
    NHYDROGEN = 0
    IF (doH) THEN
        DO I=1, NDMSO
            IF (DMSO_SYM(I) .EQ. "H") NHYDROGEN = NHYDROGEN + 1
        END DO
    END IF

    IF (.NOT. USE_OUT) OPEN (IOout, FILE=TRIM(OUT_FILE))

    IF (CORNERS) THEN
        N = NCOM * (NDMSO - NHYDROGEN) + NSOL + 8
    ELSE
        N = NCOM * (NDMSO - NHYDROGEN) + NSOL
    END IF

    WRITE (IOout,*) N
    WRITE (IOout, "(A, ' - ',A)") TRIM(BOXID), TRIM(SOLID)

    IF (CORNERS) THEN
        DO I=0,1
            DO J=0,1
                DO K=0,1
                    WRITE (IOout,*) "Ne", (-1.D0)**I * BOXL, (-1.D0)**J * BOXL, (-1.D0)**K * BOXL
                END DO
            END DO
        END DO
    END IF

    SOL_loop: DO I=1,NSOL
        WRITE (IOout,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO SOL_loop

    CoM_loop: DO I=1,NCOM
        MOL1 = RotMatrix(CoM(I), DMSO, HOEK(I))
        DO J=1,NDMSO
            IF (DMSO_SYM(J) .NE. "H" .OR. doH) THEN
                WRITE (IOout, *) DMSO_sym(J), MOL1(J)%X, MOL1(J)%Y, MOL1(J)%Z
            END IF
        END DO
    END DO CoM_loop

    IF (.NOT. USE_OUT) CLOSE(IOout)

END SUBROUTINE doBox

END PROGRAM convert2
