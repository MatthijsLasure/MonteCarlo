program boxbin
    USE VECTOR_CLASS
    USE LIB_READ
    USE LIB

    implicit none

    CHARACTER*255 :: DMSO_FILE, BOX_FILE, OUT_FILE, SOLNAME, BOX_ID
    TYPE (VECTOR), DIMENSION(:), ALLOCATABLE :: DMSO, COM, HOEK, SOLUTE
    CHARACTER*4, DIMENSION(:), ALLOCATABLE :: DMSO_SYM, SOL_SYM
    INTEGER :: NDMSO, NCOM, NSOL, I
    DOUBLE PRECISION :: E_DMSO, ENG, E_SOL, BOXL
    LOGICAL :: INBIN = .FALSE., DOPRINT = .TRUE.

    ! Array met bindingen die gebruikt mogen worden voor de dihedrale rotatie
    INTEGER, DIMENSION(:,:), ALLOCATABLE        :: DIHOEK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DROTSOLU_ARRAY
    INTEGER :: NDIHOEK

    DMSO_FILE = "na"
    OUT_FILE = "no"

    CALL GETARGS

    CALL READ_DMSO(DMSO_FILE, NDMSO, DMSO_SYM, DMSO, E_DMSO)

    IF(INBIN) THEN
        CALL READ_BIN(BOX_FILE, BOXL, BOX_ID, ENG, NCOM, COM, HOEK, SOLNAME, &
NSOL, SOLUTE, SOL_SYM, E_SOL, NDIHOEK, DIHOEK, DROTSOLU_ARRAY)
    ELSE
        CALL READ_NORM(BOX_FILE, BOXL, BOX_ID, ENG, NCOM, COM, HOEK, SOLNAME, &
NSOL, SOLUTE, SOL_SYM, E_SOL, NDIHOEK, DIHOEK, DROTSOLU_ARRAY)
    END IF

    IF (DOPRINT) THEN
        CALL PRINTIT
    END IF

    IF (OUT_FILE .NE. "no") THEN
        CALL WRITE_BIN(OUT_FILE, BOXL, BOX_ID, ENG, NCOM, COM, HOEK, SOLNAME, &
NSOL, SOLUTE, SOL_SYM, E_SOL, NDIHOEK, DIHOEK, DROTSOLU_ARRAY)
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



SUBROUTINE PRINTIT

    900 FORMAT (A, I6.6)
    901 FORMAT (A, F8.2)
    902 FORMAT (6F8.4)
    903 FORMAT (A, 3F8.4)
    904 FORMAT (2I3.3,F8.2)

    WRITE (*,901) "BoxL: ", BOXL
    WRITE (*,900) "nCoM: ", NCOM
    WRITE (*,901) "ENG:  ", ENG
    WRITE (*, "(A,A)") "BOX ID: ", TRIM(BOX_ID)
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
