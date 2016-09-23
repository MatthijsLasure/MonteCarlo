MODULE LIB_READ
    USE vector_class

    IMPLICIT NONE

!    CHARACTER*255 :: DMSO_FILE, BOX_FILE, OUT_FILE, SOLNAME, BOX_ID
!    TYPE (VECTOR), DIMENSION(:), ALLOCATABLE :: DMSO, COM, HOEK, SOLUTE
!    CHARACTER*4, DIMENSION(:), ALLOCATABLE :: DMSO_SYM, SOL_SYM
!    INTEGER :: NDMSO, NCOM, NSOL, ISTAT, IOWORK=10, I
!    DOUBLE PRECISION :: E_DMSO, ENG, E_SOL, BOXL
!    LOGICAL :: INBIN = .FALSE., DOPRINT = .TRUE.
!
!    ! Array met bindingen die gebruikt mogen worden voor de dihedrale rotatie
!    INTEGER, DIMENSION(:,:), ALLOCATABLE        :: DIHOEK
!    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DROTSOLU_ARRAY
!    INTEGER :: NDIHOEK
CONTAINS

SUBROUTINE READ_DMSO(DMSO_FILE, NDMSO, DMSO_SYM, DMSO, E_DMSO)

CHARACTER*255, INTENT(IN) :: DMSO_FILE
INTEGER, INTENT(OUT) :: NDMSO
CHARACTER*4, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: DMSO_SYM
TYPE(VECTOR), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: DMSO
DOUBLE PRECISION :: E_DMSO

INTEGER :: ISTAT, I, IOwork = 10


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
SUBROUTINE read_norm(BOX_FILE, BOXL, BOX_ID, ENG, NCOM, COM, HOEK, SOLNAME, &
NSOL, SOLUTE, SOL_SYM, E_SOL, NDIHOEK, DIHOEK, DROTSOLU_ARRAY)


    CHARACTER*255, INTENT(IN) :: BOX_FILE
    CHARACTER*255, INTENT(OUT) :: SOLNAME, BOX_ID
    TYPE (VECTOR), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: COM, HOEK, SOLUTE
    CHARACTER*4, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: SOL_SYM
    INTEGER, INTENT(OUT) :: NCOM, NSOL
    INTEGER :: ISTAT, IOWORK=10, I
    DOUBLE PRECISION, INTENT(OUT) :: ENG, E_SOL, BOXL

    ! Array met bindingen die gebruikt mogen worden voor de dihedrale rotatie
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)        :: DIHOEK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: DROTSOLU_ARRAY
    INTEGER, INTENT(OUT) :: NDIHOEK

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

    READ (IOwork, "(A255)") BOX_ID ! Box ID
    ALLOCATE(CoM(NCOM))
    ALLOCATE(hoek(NCOM))
    DO I= 1, NCOM
        READ(IOwork,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO

    ! solute.txt: conformatie opgeloste molecule (sol)
    ! NOW INCORPERATED IN BOX.TXT
    READ (IOwork, *) ! Line with SOLUTE
    READ (IOwork, *) nSol ! Lees aantal atomen
    READ (IOwork, "(A255)") SOLNAME ! Comment line
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
SUBROUTINE read_bin(BOX_FILE, BOXL, BOX_ID, ENG, NCOM, COM, HOEK, SOLNAME, &
NSOL, SOLUTE, SOL_SYM, E_SOL, NDIHOEK, DIHOEK, DROTSOLU_ARRAY)

    CHARACTER*255, INTENT(IN) :: BOX_FILE
    CHARACTER*255, INTENT(OUT) :: SOLNAME, BOX_ID
    TYPE (VECTOR), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: COM, HOEK, SOLUTE
    CHARACTER*4, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: SOL_SYM
    INTEGER, INTENT(OUT) :: NCOM, NSOL
    INTEGER :: ISTAT, IOWORK=10, I
    DOUBLE PRECISION, INTENT(OUT) :: ENG, E_SOL, BOXL

    ! Array met bindingen die gebruikt mogen worden voor de dihedrale rotatie
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)        :: DIHOEK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: DROTSOLU_ARRAY
    INTEGER, INTENT(OUT) :: NDIHOEK

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

SUBROUTINE WRITE_BIN(OUT_FILE, BOXL, BOX_ID, ENG, NCOM, COM, HOEK, SOLNAME, &
NSOL, SOLUTE, SOL_SYM, E_SOL, NDIHOEK, DIHOEK, DROTSOLU_ARRAY)

    CHARACTER*255, INTENT(IN) :: OUT_FILE
    CHARACTER*255, INTENT(IN) :: SOLNAME, BOX_ID
    TYPE (VECTOR), DIMENSION(:), ALLOCATABLE, INTENT(IN) :: COM, HOEK, SOLUTE
    CHARACTER*4, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: SOL_SYM
    INTEGER, INTENT(IN) :: NCOM, NSOL
    INTEGER :: ISTAT, IOWORK=10, I
    DOUBLE PRECISION, INTENT(IN) :: ENG, E_SOL, BOXL

    ! Array met bindingen die gebruikt mogen worden voor de dihedrale rotatie
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(IN)        :: DIHOEK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(IN) :: DROTSOLU_ARRAY
    INTEGER, INTENT(IN) :: NDIHOEK

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

    ! DMSO: return DMSO molecule
    !===========================
    SUBROUTINE DMSO_init(NDMSO, SYM, LOC, E)
        INTEGER             :: NDMSO
        DOUBLE PRECISION    :: E
        CHARACTER*4 , DIMENSION(:)      :: SYM
        TYPE (vector), DIMENSION(:)     :: LOC

        NDMSO = 10
        E =  -.0499674791107

        LOC(01) = VECTOR(-0.000001,  0.243032, -0.439459)
        LOC(02) = VECTOR(-1.363478, -0.822308,  0.180148)
        LOC(03) = VECTOR(-1.275543, -0.934699,  1.264845)
        LOC(04) = VECTOR(-2.300201, -0.313459, -0.059709)
        LOC(05) = VECTOR(-1.332307, -1.796078, -0.318322)
        LOC(06) = VECTOR( 1.363477, -0.822310,  0.180148)
        LOC(07) = VECTOR( 1.332384, -1.796034, -0.318417)
        LOC(08) = VECTOR( 2.300186, -0.313388, -0.059610)
        LOC(09) = VECTOR( 1.275493, -0.934809,  1.264829)
        LOC(10) = VECTOR( 0.000001,  1.508457,  0.386993)

        SYM(01) = "S"
        SYM(02) = "C"
        SYM(03) = "H"
        SYM(04) = "H"
        SYM(05) = "H"
        SYM(06) = "C"
        SYM(07) = "H"
        SYM(08) = "H"
        SYM(09) = "H"
        SYM(10) = "O"

    END SUBROUTINE DMSO_init

END MODULE LIB_READ
