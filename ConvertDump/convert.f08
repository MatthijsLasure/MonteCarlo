program convert
    use lib
    implicit none

    CHARACTER*1000   :: in_file, out_file, sol_file, DMSO_file, temp
    CHARACTER*20     :: what_char, A_char, B_char
    TYPE (vector)    :: CoM, hoek
    TYPE (vector)    :: vecA, vecB
    DOUBLE PRECISION :: R, RMIN
    INTEGER          :: IMIN
    TYPE(vector), DIMENSION(:), ALLOCATABLE :: solute, DMSO, pos, MOL1
    CHARACTER*4, DIMENSION(:), ALLOCATABLE      :: DMSO_SYM, SOL_SYM
    CHARACTER*4      :: symbol
    INTEGER          :: nDMSO, nSOL, NCOM = 30, N, TIMESTEP
    INTEGER          :: what, I, J, IOSTATUS, A, B, NHYDROGEN
    LOGICAL          :: doH = .FALSE.

    WRITE (*,*) "***************"
    WRITE (*,*) "* ConvertDump *"
    WRITE (*,*) "***************"

    IF( COMMAND_ARGUMENT_COUNT() .LT. 5) THEN
        WRITE (0,*) "Incorrect arguments!"
        WRITE (*,*) "ConvertDump what in out solute DMSO"
        STOP
    END IF

    CALL GET_COMMAND_ARGUMENT(1, what_char)
    CALL GET_COMMAND_ARGUMENT(2, in_file)
    CALL GET_COMMAND_ARGUMENT(3, out_file)
    CALL GET_COMMAND_ARGUMENT(4, sol_file)
    CALL GET_COMMAND_ARGUMENT(5, DMSO_file)

        ! Laden van configuraties
!====================================================================
!====================================================================

    WRITE (*,*) "Loading data..."

    ! DMSO.txt: conformatie DMSO
    OPEN (UNIT=10, FILE=trim(dmso_file))
    READ (10, *) nDMSO ! Lees aantal atomen
    ALLOCATE(DMSO(NDMSO)) ! Ken correcte groottes toe aan de arrays
    ALLOCATE(MOL1(NDMSO))
    ALLOCATE(DMSO_sym(NDMSO))
    DO I=1, NDMSO
        READ (10,*) DMSO_sym(I), DMSO(I)%X, DMSO(I)%Y, DMSO(I)%Z
    END DO
    CLOSE(10)

    ! solute.txt: conformatie opgeloste molecule (sol)
    OPEN (UNIT=10, FILE=trim(sol_file))
    READ (10, *) nSol ! Lees aantal atomen
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(sol_sym(NSOL))
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (10,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    CLOSE(10)

    WRITE (*,*) "Done loading data!"

!====================================================================
!====================================================================

    READ(what_char, *) what
    IF (what .EQ. 1) THEN
        doH = .TRUE.
        NHYDROGEN = 0
        DO I=1, NDMSO
            IF (DMSO_SYM(I) .EQ. "H") NHYDROGEN = NHYDROGEN + 1
        END DO
    END IF

    OPEN(10, FILE=trim(in_file))
    OPEN(11, FILE=trim(out_file))

    WRITE (*,*) "Processing dump. This may take a while..."

IF (what .LE. 1) THEN ! Just dump XYZ

    loop: DO
        READ (10, *, IOSTAT=IOSTATUS) N
        IF (IOSTATUS < 0) exit
        READ (10, "(A10,I20.20)") TEMP, TIMESTEP
        IF (doH) THEN
            N = NCOM * NDMSO + NSOL
        ELSE
            N = NCOM * (NDMSO - NHYDROGEN) + NSOL
        END IF
        WRITE (11,*) N
        WRITE (11, *) "Timestep:", TIMESTEP
        SOL_loop: DO I=1,NSOL
            WRITE (11,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
        END DO SOL_loop
        CoM_loop: DO I=1,NCOM
            READ (10, *) CoM%X, CoM%Y, CoM%Z, hoek%X, hoek%Y, hoek%Z
            MOL1 = RotMatrix(CoM, DMSO, hoek)
            DO J=1,NDMSO
                IF (DMSO_SYM(J) .NE. "H" .OR. doH) THEN
                    WRITE (11, *) DMSO_sym(J), MOL1(J)%X, MOL1(J)%Y, MOL1(J)%Z
                END IF
            END DO
        END DO CoM_loop
    END DO loop


ELSE IF (what .EQ. 2) THEN! Calculate distances

IF (COMMAND_ARGUMENT_COUNT() .LT. 7) THEN
    WRITE (0,*) "Not enough arguments!"
    WRITE (*,*) "ConvertDump 2 in out solute DMSO A B"
    STOP
END IF
CALL GET_COMMAND_ARGUMENT(6, A_char)
CALL GET_COMMAND_ARGUMENT(7, B_char)
READ(A_char, *) A
READ(B_char, *) B

vecA = SOLUTE(A)
    DO
        READ (10, *, IOSTAT=IOSTATUS) N
        IF (IOSTATUS < 0) exit
        READ (10, "(A10,I20.20)") TEMP, TIMESTEP
        RMIN = HUGE(R)
        IMIN = HUGE(I)
        DO I=1,NCOM
            READ (10, *) CoM%X, CoM%Y, CoM%Z, hoek%X, hoek%Y, hoek%Z
            MOL1 = RotMatrix(CoM, DMSO, hoek)
            vecB = MOL1(B)
            R = getDist(vecA, vecB)
            IF (R .LT. RMIN) THEN
                RMIN = R
                IMIN = I
            END IF
        END DO
        WRITE (11, *) TIMESTEP, IMIN, RMIN
    END DO
END IF

CLOSE(10)
CLOSE(11)

WRITE (*,*) "ConvertDump is ready!"

end program convert
