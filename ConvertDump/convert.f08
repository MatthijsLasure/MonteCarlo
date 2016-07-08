program convert
    use lib
    implicit none

    CHARACTER*1000   :: in_file, out_file, sol_file, DMSO_file, temp
    CHARACTER*20     :: what_char, A_char, B_char
    TYPE (vector)    :: CoM, hoek
    TYPE (vector)    :: vecA, vecB
    DOUBLE PRECISION :: R, RMIN, En
    INTEGER          :: IMIN, AMIN
    TYPE(vector), DIMENSION(:), ALLOCATABLE :: solute, DMSO, pos, MOL1, UCOM, UHOEK
    CHARACTER*4, DIMENSION(:), ALLOCATABLE      :: DMSO_SYM, SOL_SYM
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DIST
    CHARACTER*4      :: symbol
    INTEGER          :: nDMSO, nSOL, NCOM, N, TIMESTEP
    INTEGER          :: what, I, J, K, IOSTATUS, A, B, NHYDROGEN
    LOGICAL          :: doH = .FALSE.
    DOUBLE PRECISION :: BOXL

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

    READ(what_char, *) what

        ! Laden van configuraties
!====================================================================
!====================================================================

    WRITE (*,*) "Loading data..."

    ! DMSO.txt: conformatie DMSO
    OPEN (UNIT=10, FILE=trim(dmso_file))
    READ (10, *) nDMSO ! Lees aantal atomen
    READ (10, *) ! Comment
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
    READ (10, *) ! Comment
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(sol_sym(NSOL))
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (10,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    CLOSE(10)

    WRITE (*,*) "Done loading data!"

!====================================================================
!====================================================================

    NHYDROGEN = 0
    IF (what .EQ. 1) THEN
        doH = .TRUE.
    ELSE
        DO I=1, NDMSO
            IF (DMSO_SYM(I) .EQ. "H") NHYDROGEN = NHYDROGEN + 1
        END DO
    END IF


    OPEN(10, FILE=trim(in_file))
    OPEN(11, FILE=trim(out_file))

    WRITE (*,*) "Processing dump. This may take a while..."

    READ (10, *) BOXL

IF (what .LE. 1) THEN ! Just dump XYZ

    loop: DO
        READ (10, *, IOSTAT=IOSTATUS) NCOM
        IF (IOSTATUS < 0) exit
        READ (10, "(E20.10)", IOSTAT=IOSTATUS) En
        IF (IOSTATUS .NE. 0) En = 0.D0
        READ (10, "(A10,I20.20, F20.20)") TEMP, TIMESTEP

        N = NCOM * (NDMSO - NHYDROGEN) + NSOL + 8
        WRITE (11,*) N
        WRITE (11, "('ID ', I6.6, ' Energy: ', F20.16, ' kJ/mol')") TIMESTEP, En

        DO I=0,1
            DO J=0,1
                DO K=0,1
                    WRITE (11,*) "Ne", (-1.D0)**I * BOXL, (-1.D0)**J * BOXL, (-1.D0)**K * BOXL
                END DO
            END DO
        END DO

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
        READ (10, *, IOSTAT=IOSTATUS) NCOM
        IF (IOSTATUS < 0) exit
        READ (10, "(A10,I20.20)") TEMP, TIMESTEP
        RMIN = HUGE(R)
        IMIN = HUGE(I)
        DO I=1,NCOM
            !READ (10, "(6F24.16)") CoM%X, CoM%Y, CoM%Z, hoek%X, hoek%Y, hoek%Z
            READ (10,"(A)") TEMP
            !WRITE (*,*) trim(TEMP)
            READ (TEMP, *) CoM%X, CoM%Y, CoM%Z, hoek%X, hoek%Y, hoek%Z
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


ELSE IF (what .EQ. 3) THEN! Calculate all distances


    DO
        READ (10, *, IOSTAT=IOSTATUS) NCOM
        IF (IOSTATUS < 0) exit
        READ (10, "(A10,I20.20)") TEMP, TIMESTEP
            DO I=1,NCOM
                !READ (10, "(6F24.16)") CoM%X, CoM%Y, CoM%Z, hoek%X, hoek%Y, hoek%Z
                READ (10,"(A)") TEMP
                !WRITE (*,*) trim(TEMP)
                READ (TEMP, *) CoM%X, CoM%Y, CoM%Z, hoek%X, hoek%Y, hoek%Z
                MOL1 = RotMatrix(CoM, DMSO, hoek)

                RMIN = HUGE(R)
                IMIN = HUGE(I)
                AMIN = HUGE(A)
                DO A=1,nSOL
                    vecA = SOLUTE(A)
                    DO B=1,nDMSO
                        vecB = MOL1(B)
                        R = getDist(vecA, vecB)
                        IF (R .LT. RMIN) THEN
                            RMIN = R
                            IMIN = B
                            AMIN = A
                        END IF
                    END DO
                END DO
                WRITE (11, *) TIMESTEP, AMIN, I, IMIN, RMIN
            END DO
    END DO

END IF

CLOSE(10)
CLOSE(11)

WRITE (*,*) "ConvertDump is ready!"

end program convert
