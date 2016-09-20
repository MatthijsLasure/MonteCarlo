program postProcess
    USE RANDGEN

    implicit none

    INTEGER :: NUMARG, I, J, K, ISTAT, IOwork = 10, RUNS, NUMLOG
    CHARACTER*255, DIMENSION(:), ALLOCATABLE :: LOGLIST
    CHARACTER*255, DIMENSION(:,:), ALLOCATABLE :: BOXES
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ENERGY
    INTEGER, DIMENSION(:), ALLOCATABLE :: OLD, NOW, S
    LOGICAL, DIMENSION(:), ALLOCATABLE :: ACC
    LOGICAL :: ISOK, DRYRUN = .FALSE.
    CHARACTER*8 :: PRESTRING
    CHARACTER*255 :: OUTDIR, OUTFILE
    CHARACTER*2000 :: FORM, ROW, STR, COMMAND
    DOUBLE PRECISION :: TEMPERATURE = 300.D0

#ifdef WINDOWS
    CHARACTER*5, PARAMETER :: COPYCOM = "COPY "
    CHARACTER*1, PARAMETER :: SEP = "\"
    CHARACTER*1, PARAMETER :: BAD = "/"
#else
    CHARACTER*3, PARAMETER :: COPYCOM = "cp "
    CHARACTER*1, PARAMETER :: SEP = "/"
    CHARACTER*1, PARAMETER :: BAD = "\"
#endif


    CALL random_seed(size = I)
    ALLOCATE (S(I))
    INQUIRE(FILE='seed',EXIST=ISOK)
    IF(ISOK) THEN
        OPEN (UNIT=IOwork, FILE='seed', ACTION='READ', STATUS='OLD', IOSTAT=istat)
        READ (IOwork, *) S
        CLOSE(IOwork)
        WRITE (*,"(A)") "File 'seed' was found, read the seed succesfully."
        WRITE (*,*) S
        CALL random_seed(put=S)
    ELSE
        CALL init_random_seed()
        CALL random_seed(get=S)
        OPEN (UNIT=IOwork, FILE='seed', ACTION='WRITE', STATUS='NEW', IOSTAT=istat)
        WRITE (IOwork, *) S
        CLOSE(IOwork)
        WRITE (*,"(A)") "File 'seed' was not found, a new seed was generated."
        WRITE (*,*) S
    END IF


    WRITE (*,"(5A)") "Will use '", COPYCOM, "' for copying and '", SEP, "' for directory separation."
    WRITE (*,"(A)") "To change, use '-D [LINUX/WINDOWS]' while compiling."

    NUMARG = COMMAND_ARGUMENT_COUNT()

    IF (NUMARG .EQ. 0) THEN
        CALL HELP()
        STOP
    END IF
    CALL GETARGS()

    ! Get number of runs in the first logfile
    CALL RECON()
    WRITE (*,'(A,I3.3,A,I3.3,A)') "Energy array will be ", RUNS, " runs from ", NUMLOG, " logs."

    ALLOCATE(BOXES(NUMLOG,RUNS))
    ALLOCATE(ENERGY(NUMLOG, RUNS))
    ALLOCATE(OLD(RUNS))
    ALLOCATE(NOW(RUNS))
    ALLOCATE(ACC(RUNS))

    DO I=1,NUMLOG
        CALL READLOG(LOGLIST(I), I)
    END DO


    ! PRINTDO
    FORM = "('| ', I4.4, ' | '"
    ROW = "+------+"
    DO I=1, NUMLOG
        FORM = TRIM(FORM) // ", F8.2, ' | '"
        ROW = TRIM(ROW) // "----------+"
    END DO
    FORM = TRIM(FORM) // ")"
    WRITE (*,"(A)") TRIM(ROW)
    DO I=1, RUNS
        WRITE (*,FORM) I,  ENERGY(:,I)
    END DO
    WRITE (*,"(A)") TRIM(ROW)

    ! END READING

    WRITE (*,"(A,F6.2,A)") "Temperature: ", TEMPERATURE, " K"
    K = 0
    CALL COMPARE()

    ! PRINTING
    901 FORMAT(A, I3.3, A, I3.3, A, F6.2, A)
    WRITE (*,901) "Accepted: ", K, " / ", RUNS, " (", FLOAT(K) / FLOAT(RUNS) * 100, " %)"
    WRITE (*,'(A)') "Legend: I old solute with new rejected. - old solute. + new solute accepted. o new solute rejected"
    DO I=1,RUNS
        STR = "('| ', I4.4,' |"
        DO J=1,NUMLOG
            IF (OLD(I) .EQ. J .AND. .NOT. ACC(I)) THEN
                STR = TRIM(STR) // "I"
            ELSEIF (NOW(I) .EQ. J .AND. .NOT. ACC(I)) THEN
                STR = TRIM(STR) // "o"
            ELSEIF ( OLD(I) .EQ. J) THEN
                STR = TRIM(STR) // "-"
            ELSEIF ( NOW(I) .EQ. J) THEN
                STR = TRIM(STR) // "+"
            ELSE
                STR = TRIM(STR) // "."
            END IF
        END DO
        STR = TRIM(STR) // "|')"
        WRITE (*,TRIM(STR)) I
    END DO

    ! COPY
    OUTDIR = "output"
#ifdef __GFORTRAN__
    INQUIRE(FILE=OUTDIR, EXIST=ISOK)
#endif
    IF(ISOK) THEN
        WRITE (*,"(A)") "Directory exists."
    ELSE
        WRITE (*,"(A)") "Directory does not exist, creating..."
        COMMAND="mkdir " // trim(OUTDIR)
        CALL SYSTEM(TRIM(COMMAND))
    END IF

    DO I=1,RUNS
        IF (ACC(I)) THEN ! If accepted
            J = NOW(I)
            WRITE(OUTFILE, "('box_', I3.3, '_', I3.3, '.txt')") I, J
            OUTFILE = TRIM(OUTDIR) // SEP // TRIM(OUTFILE)
            INQUIRE(FILE=TRIM(OUTFILE), EXIST=ISOK)
            COMMAND = COPYCOM // TRIM(BOXES(J, I))
            COMMAND = TRIM(COMMAND) // " " // TRIM(OUTFILE)
            IF (ISOK) THEN
                WRITE (*,"(A,A,A)") "File '", TRIM(OUTFILE), "' exists. Skipping this one."
            ELSE
                WRITE (*,"(A)") TRIM(COMMAND)
                IF(.NOT. DRYRUN) CALL SYSTEM(TRIM(COMMAND))
            END IF
        END IF
    END DO

CONTAINS

! =============================================================

SUBROUTINE GETARGS() ! TODO more complicated
    CHARACTER*255 :: MYARG, TEMP_C
    LOGICAL :: LIST
    INTEGER :: SKIP = 0, J = 0, START

    J = 0
    DO I=1,NUMARG
        IF (SKIP .GT. 0) THEN
            SKIP = SKIP - 1
            CYCLE
        END IF
        CALL GET_COMMAND_ARGUMENT(I, MYARG)
        SELECT CASE (MYARG)
            CASE ("-l")
                LIST = .TRUE.
                START = I + 1
            CASE ("-o")
                SKIP = 1
                CALL GET_COMMAND_ARGUMENT(I+1, OUTDIR)
                LIST = .FALSE.
            CASE ("-T")
                SKIP = 1
                CALL GET_COMMAND_ARGUMENT(I+1, TEMP_C)
                READ(TEMP_C, *) TEMPERATURE
                LIST = .FALSE.
            CASE ("-d")
                DRYRUN = .TRUE.
                LIST = .FALSE.
            CASE DEFAULT
                IF (LIST) THEN
                    !LOGLIST(J) = TRIM(MYARG)
                    J = J + 1
                ELSE
                    WRITE (0,"(A,A)") "Unknown argument: ", TRIM(MYARG)
                END IF
        END SELECT
    END DO
    NUMLOG = J

    ALLOCATE(LOGLIST(NUMLOG))
    J = 1
    DO I=START,NUMLOG+START-1
        CALL GET_COMMAND_ARGUMENT(I, MYARG)
        LOGLIST(J) = TRIM(MYARG)
        J = J + 1
    END DO
END SUBROUTINE GETARGS

! =============================================================

SUBROUTINE RECON()
    INTEGER :: I, J = 0, IOwork = 10
    INTEGER, DIMENSION(NUMLOG) :: LOGRUN

    DO I=1,NUMLOG
        J = 0
        OPEN (IOwork, FILE=LOGLIST(I))
        DO
            READ (IOwork,'(A8)', IOSTAT=ISTAT) PRESTRING
            IF (ISTAT .NE. 0) EXIT ! EOF
            IF (PRESTRING .EQ. ' Energy:') J = J + 1
        END DO
        CLOSE (IOwork)
        LOGRUN(I) = J
    END DO

    RUNS = MINVAL(LOGRUN)

    WRITE (*,'(A, I3.3, A)') "Found ", RUNS, " common runs."
    DO I=1, NUMLOG
        WRITE (*,'(A10, 1X, I3.3)') TRIM(LOGLIST(I)), LOGRUN(I)
    END DO

END SUBROUTINE RECON

! =============================================================

SUBROUTINE READLOG(LOGNAME, I)
    CHARACTER*255 :: LOGNAME, FULLSTRING, STUFF, BOX
    DOUBLE PRECISION :: E
    INTEGER :: I, J = 0, IOwork = 20

    OPEN (IOwork, FILE=LOGNAME)
    J = 0
    DO
        READ (IOwork,'(A255)', IOSTAT=ISTAT) FULLSTRING
        IF (ISTAT .NE. 0) EXIT ! EOF
        IF (FULLSTRING(1:7) .EQ. 'BOX    ') THEN
            J = J + 1
            IF (J .GT. RUNS) EXIT
            READ (FULLSTRING, "(7X, A100)") BOX
            BOX = Replace_Text(BOX, BAD, SEP)
            BOXES(I, J) = TRIM(BOX)
        END IF
        IF (FULLSTRING(1:8) .EQ. ' Energy:') THEN
            !J = J + 1
            READ(FULLSTRING, '(A8, F25.14)') STUFF, E
            ENERGY(I, J) = E
            !IF (J .GE. RUNS) EXIT
        END IF

    END DO

    CLOSE (IOwork)
END SUBROUTINE READLOG

! =============================================================

SUBROUTINE COMPARE()
    DOUBLE PRECISION :: RV, KANS, EXPONENT, DELTA
    INTEGER :: I, CURRENT, NEXT

    CURRENT = MINLOC(ENERGY(:,1),1)
    OLD(1) = CURRENT

    DO I=1,RUNS
        RV = RAND()
        IF (RV .LT. 0.5) THEN
            NEXT = CURRENT + 1
        ELSE
            NEXT = CURRENT - 1
        END IF
        IF (NEXT .GT. NUMLOG) NEXT = 1
        IF (NEXT .LE. 0) NEXT = NUMLOG

        DELTA = ENERGY(I, NEXT) - ENERGY(I, CURRENT)
        EXPONENT = -1.D0 * DELTA  * 1000.D0 / (8.314D0 * TEMPERATURE)
        IF (EXPONENT .LT. -75.D0) THEN ! e^-75 < 3*10^-33: 0% kans anyway
        !write(IOerr,*) "Large Exponent!", I
            KANS = 0.D0
            RV = 1.D0 ! Skip rand() voor cpu tijd besparing
        ELSE IF(EXPONENT .GE. 0.D0) THEN ! Lager in energie, dus 100% kans
            KANS = 1.D0
            RV = 0.D0
        ELSE
            KANS = EXP(EXPONENT)
            RV = RAND()
            IF (KANS .GT. 1.D0) KANS = 1.D0
        END IF

        IF (RV .LT. KANS) THEN ! Success!
            K = K + 1
            ACC(I) = .TRUE.
            OLD(I) = CURRENT
            NOW(I) = NEXT
            CURRENT = NEXT
        ELSE
            ACC(I) = .FALSE.
            OLD(I) = CURRENT
            NOW(I) = NEXT
        END IF
    END DO
END SUBROUTINE COMPARE

FUNCTION Replace_Text (s,text,rep)  RESULT(outs)
    CHARACTER(*)        :: s,text,rep
    CHARACTER(LEN(s)+100) :: outs     ! provide outs with extra 100 char len
    INTEGER             :: i, nt, nr

    outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
    DO
       i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
       outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
    END DO
END FUNCTION Replace_Text

! =============================================================

SUBROUTINE HELP()
WRITE (*,"(A)") "PostProcess: extract boxes from the simulations.                                  "
WRITE (*,"(A)") "                                                                                  "
WRITE (*,"(A)") "Must be used in the directory of the MonteCarlo program (uses relative box paths) "
WRITE (*,"(A)") "                                                                                  "
WRITE (*,"(A)") "Usage: PostProcess.exe -l <logs> [args]                                           "
WRITE (*,"(A)") "                                                                                  "
WRITE (*,"(A)") "-l <logs>                                                                         "
WRITE (*,"(A)") "List of log files to process.                                                     "
WRITE (*,"(A)") "                                                                                  "
WRITE (*,"(A)") "-T <temp>                                                                         "
WRITE (*,"(A)") "Temperature for the MC, in Kelvin.                                                "
WRITE (*,"(A)") "                                                                                  "
WRITE (*,"(A)") "-d                                                                                "
WRITE (*,"(A)") "Do not actually copy the files, just print information.                           "
WRITE (*,"(A)") "                                                                                  "
WRITE (*,"(A)") "NOTICE                                                                            "
WRITE (*,"(A)") "This program will save the seed upon the first run.                               "
WRITE (*,"(A)") "In this way, the results will be the same in the same folder.                     "
WRITE (*,"(A)") "Please do not delete the 'seed' file.                                             "
END SUBROUTINE HELP

end program postProcess
