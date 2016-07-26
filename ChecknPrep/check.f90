PROGRAM CHECKNPREP
    USE LIB

    IMPLICIT NONE

    ! VARIABLES
    CHARACTER*1000  :: FILE, TEMP, THING, SOLNAME, BOXID, STR
    INTEGER         :: LINE = 0, I, IOwork=10, istat=0, NERR = 0, NWAR = 0

    ! Energiën van de moleculen, bekomen via extern programma
    DOUBLE PRECISION                :: E_DMSO, E_SOL, BOXL, PRE_ENG

    ! Vectoren voor moleculen & atoomtypes
    TYPE (vector), DIMENSION(:), ALLOCATABLE    :: DMSO, COM, SOLUTE, HOEK
    CHARACTER*4, DIMENSION(:), ALLOCATABLE      :: DMSO_SYM, SOL_SYM
    INTEGER                                     :: NDMSO, NCOM, NSOL, NPARAM, NPARSOL ! Aantal units

    ! Array met bindingen die gebruikt mogen worden voor de dihedrale rotatie
    INTEGER, DIMENSION(:,:), ALLOCATABLE        :: DIHOEK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DROTSOLU_ARRAY
    INTEGER :: NDIHOEK

    ! Arrays voor parameters van DMSO (Q, epsilon, sigma, mass)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: Q, EPSILON, SIGMA
    CHARACTER*4, DIMENSION(:), ALLOCATABLE          :: SYM
    ! Arrays voor paremeters van solute (epsilon, sigma)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: SOL_Q, SOL_EPSILON, SOL_SIGMA
    CHARACTER*4, DIMENSION(:), ALLOCATABLE          :: SOLPAR_SYM
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: TABLE_DMSO, TABLE_SOL

    ! Format for errors
    501 FORMAT('Error opening file, error ', I4.4)

    ! What file?
    IF (COMMAND_ARGUMENT_COUNT() .EQ. 0) THEN
        WRITE (*,*) "Please give a file"
        CALL EXIT(20)
    END IF
    CALL get_command_argument(1, FILE)
    CALL get_command_argument(2, THING)

    WRITE (*,"(A,A)") "Running on file: ", TRIM(FILE)

    ! OPEN FILE
    OPEN (UNIT=IOwork, FILE=FILE, ACTION='READ', STATUS='OLD', IOSTAT=istat)
    IF (istat .NE. 0) THEN
        WRITE (*,501) istat
        CALL EXIT(1)
    END IF


    ! READ
    SELECT CASE(THING)
        CASE ("dmso")
            CALL readDMSO()
        CASE ("box")
            CALL readBOX()
    END SELECT

    ! PRINT INFO
    IF (NERR .EQ. 0) THEN
        WRITE (*,*) "No read errors occured!"
    ELSEIF (NERR .EQ. 1) THEN
        WRITE (*,"(A)") "One error has occured. Please check the file."
    ELSE
        WRITE (*,"(I4.4, A)") NERR, " errors have occured. Please check the file."
        CALL EXIT(2)
    END IF

    CLOSE(IOwork)

    CONTAINS
SUBROUTINE readE(thing)

    CHARACTER(len=*) :: THING


    500 FORMAT('Error reading ', A, ' at line ', I4.4, ': ID ', I6.6)
    LINE = LINE + 1
    IF (istat .NE. 0) THEN
        WRITE (*,500) THING, LINE, ISTAT
        BACKSPACE(IOwork)
        READ (IOwork,"(A)", IOSTAT=istat) STR
        WRITE (*,"(A)") TRIM(STR)
        NERR = NERR + 1
    END IF

END SUBROUTINE readE

SUBROUTINE loopE(thing, N)
    CHARACTER(len=*) :: THING
    INTEGER :: N
    502 FORMAT(A, A, A, I4.4, A, I3.3, A, I3.3)
    IF (istat .EQ. 5001) THEN ! EOF
        WRITE(*,502) "End of file while reading ", thing, &
        ". At line ", LINE+1, ", was at element ", I, " of ", N
        CALL EXIT(3)
    ELSE
        CALL readE(thing)
    END IF
END SUBROUTINE loopE

SUBROUTINE readDMSO()
        READ (IOwork, *, IOSTAT=istat) nDMSO ! Lees aantal atomen
    CALL readE("nDMSO")

    READ (IOwork, *, IOSTAT=istat) TEMP! Comment line
    CALL readE("comment")

    ALLOCATE(DMSO(NDMSO)) ! Ken correcte groottes toe aan de arrays
    ALLOCATE(DMSO_sym(NDMSO))
    ALLOCATE(TABLE_DMSO(NDMSO, 3))
    DO I=1, NDMSO
        READ (IOwork,*, IOSTAT=istat) DMSO_sym(I), DMSO(I)%X, DMSO(I)%Y, DMSO(I)%Z
        CALL loopE("coordinates", NDMSO)
    END DO
    READ (IOwork,*, IOSTAT=istat) E_DMSO
    CALL readE("E_DMSO")

END SUBROUTINE readDMSO

SUBROUTINE readBOX()
    READ (IOwork, *, IOSTAT=istat) BOXL ! Box grootte
    CALL readE("BOXL")
    READ (IOwork, *, IOSTAT=istat) nCoM ! Lees aantal moleculen
    CALL readE("nCoM")
    READ (IOwork,"(E20.10)", IOSTAT=istat) PRE_ENG ! Energy from previous run
    CALL readE("box energy")
    IF (istat .NE. 0) THEN
        WRITE (*,*) "INFO: No energy detected in input box."
        WRITE (*,*) istat
        PRE_ENG = 123456.789
        BACKSPACE(IOwork)
    END IF

    IF (BOXL .EQ. 0.D0) THEN
        WRITE (*,*) "WARNING: BOXL is 0!"
        NWAR = NWAR + 1
    END IF
    IF (NCOM .EQ. 0) THEN
        WRITE (*,*) "WARNING: nCoM is 0!"
        NWAR = NWAR + 1
    END IF
    IF (PRE_ENG .EQ. 0.D0) THEN
        WRITE (*,*) "WARNING: ENG is 0!"
        NWAR = NWAR + 1
    END IF

    READ (IOwork, *, IOSTAT=istat)  BOXID    ! Box ID
    CALL readE("box ID")
    ALLOCATE(CoM(NCOM))
    ALLOCATE(hoek(NCOM))
    DO I= 1, NCOM
        READ(IOwork,*, IOSTAT=istat) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
        CALL loopE("Center of Mass + Hoek", NCOM)

        IF (CoM(I)%X .EQ. 0.D0 .OR. CoM(I)%Y .EQ. 0.D0 .OR. CoM(I)%Z .EQ. 0.D0 .OR. &
        HOEK(I)%X .EQ. 0.D0 .OR. HOEK(I)%Y .EQ. 0.D0 .OR. HOEK(I)%Z .EQ. 0.D0) THEN
        WRITE (*,*) "WARNING: coordinates are 0! Line: ", LINE
        NWAR = NWAR + 1
    END IF
    END DO

    ! solute.txt: conformatie opgeloste molecule (sol)
    ! NOW INCORPERATED IN BOX.TXT
    READ (IOwork, *, IOSTAT=istat) TEMP! Line with SOLUTE
    CALL readE("SOLUTE line")
    READ (IOwork, *, IOSTAT=istat) nSol ! Lees aantal atomen
    CALL readE("nSOL")
    READ (IOwork, "(A)", IOSTAT=istat) SOLNAME ! Comment line
    CALL readE("solute name")
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(sol_sym(NSOL))
    ALLOCATE(SOL_Q(NSOL))
    ALLOCATE(TABLE_SOL(NSOL, 3))

    IF (NSOL .EQ. 0) THEN
        WRITE (*,*) "WARNING: nSOL is 0!"
        NWAR = NWAR + 1
    END IF

    DO I=1, NSOL ! Lees de coördinaten uit
        READ (IOwork,*, IOSTAT=istat) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
        CALL loopE("solute", NSOL)
    END DO
    READ (IOwork,*, IOSTAT=istat) E_sol
    CALL readE("E_SOL")
    READ (IOwork,*, IOSTAT=istat) NDIHOEK
    CALL readE("N DIHOEK")
    ALLOCATE(DIHOEK(NDIHOEK, 2))
    ALLOCATE(DROTSOLU_ARRAY(NDIHOEK))

    IF (E_SOL .EQ. 0.D0) THEN
        WRITE (*,*) "WARNING: E_SOL is 0!"
        NWAR = NWAR + 1
    END IF
    IF (NCOM .EQ. 0) THEN
        WRITE (*,*) "WARNING: NDHIHOEK is 0!"
        NWAR = NWAR + 1
    END IF

    DO I=1, NDIHOEK
        READ (IOwork,"(I3, 1X, I3, 1X, F10.6)", IOSTAT=istat) DIHOEK(I,1), DIHOEK(I,2), DROTSOLU_ARRAY(I)
        CALL loopE("Dihedral angle", NDIHOEK)
        IF (DROTSOLU_ARRAY(I) .EQ. 0.D0) THEN
            WRITE (*,"('INFO: MaxRot is not defined at line ', I3.3)") LINE
        END IF

        IF (DIHOEK(I,1) .GT. NSOL .OR. DIHOEK(I,2) .GT. NSOL .OR. &
        DIHOEK(I,1) .EQ. 0 .OR. DIHOEK(I,2) .EQ. 0) THEN
            WRITE (*,"(A, I4.4)") "ERROR: Values are out of bounds! Line ", LINE
            BACKSPACE (IOwork)
            READ (IOwork, "(A)") STR
            WRITE (*,"('Read as: ',A)") TRIM(STR)
            WRITE (*,"('Values: ',I3.3, 1X, I3.3, 1X, F8.2)") DIHOEK(I,1), DIHOEK(I,2), DROTSOLU_ARRAY(I)
            NERR = NERR + 1
        END IF
    END DO


END SUBROUTINE readBOX
END PROGRAM CHECKNPREP
