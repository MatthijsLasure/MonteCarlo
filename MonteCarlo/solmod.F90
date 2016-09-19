MODULE solmod

    USE mcconstants
    USE RANDGEN
    USE DIHEDRAL
    USE VECTOR_CLASS
    USE LIB
    USE OMP_LIB

    IMPLICIT NONE

    CONTAINS

    FUNCTION SOLUTE_INIT(SOL, SYM, DIHEDRAL, DROTSOLV, DROTSOLU) RESULT(SOLROT)

        TYPE (VECTOR), DIMENSION(:), INTENT(IN)  :: SOL
        TYPE (VECTOR), DIMENSION(SIZE(SOL))      :: SOLROT
        CHARACTER*4, DIMENSION(:), INTENT(IN)    :: SYM
        INTEGER, DIMENSION(:,:), INTENT(IN)      :: DIHEDRAL
        DOUBLE PRECISION, INTENT(in) :: DROTSOLU
        DOUBLE PRECISION:: HOEK, DMAX
        DOUBLE PRECISION:: HOEK_PRE
        DOUBLE PRECISION:: HOEK_POST
        DOUBLE PRECISION, DIMENSION(:)           :: DROTSOLV
        INTEGER:: A1
        INTEGER:: A2
        INTEGER:: A3
        INTEGER:: A4
        INTEGER:: NATOM
        INTEGER:: I
        INTEGER:: J
        INTEGER:: M
        INTEGER:: NDIH
        INTEGER, DIMENSION(8) :: NEIGHBOURS

        NATOM = SIZE(SOL)
        NDIH = SIZE(DIHEDRAL) / 2

        I = INT(RAND() * NDIH) + 1 ! Willeukeurige binding selecteren
        A2 = DIHEDRAL(I,1)
        A3 = DIHEDRAL(I,2)

        A1 = FIND(A2, A3, SOL, NATOM)
        A4 = FIND(A3, A2, SOL, NATOM)

        ! DRAAIEN
        IF(DROTSOLV(I) .EQ. -1.D0) THEN
            DMAX = DROTSOLU
        ELSE
            DMAX = DROTSOLV(I)
        END IF
        HOEK = RAND() * 2 * DMAX - DMAX
        HOEK_PRE = GETDIHEDRAL(SOL, A1, A2, A3, A4) * 180 / PI
        SOLROT = SETDIHEDRAL(SOL, SYM, A1, A2, A3, A4, HOEK * PI / 180)
        HOEK_POST = GETDIHEDRAL(SOLROT, A1, A2, A3, A4) * 180 / PI

        !WRITE (*,"(A,I2.2,A3,I2.2)") "Rotation around axis ", A2, " - ", A3
        WRITE (*,"(A,I2.2,A3,I2.2,A3,I2.2,A3,I2.2)") "Dihedral angle with ", A1, " - ", A2, " - ", A3, " - ", A4
        WRITE (*,"(A, F7.2, A, F7.2, A, F7.2)") "Dihedral change from ", HOEK_PRE, " with ", HOEK, " to ", HOEK_POST

    END FUNCTION SOLUTE_INIT

    FUNCTION SOLUTE_METROPOLIS(PRE_EN, POST_EN, TEMP) RESULT(isOK)
        DOUBLE PRECISION:: PRE_EN
        DOUBLE PRECISION:: POST_EN
        DOUBLE PRECISION:: TEMP
        DOUBLE PRECISION:: EXPONENT
        DOUBLE PRECISION:: RV
        DOUBLE PRECISION:: KANS
        DOUBLE PRECISION:: DELTA
        LOGICAL:: ISOK

        DELTA = POST_EN - PRE_EN
        EXPONENT = -1.D0 * DELTA  * 1000.D0 / (8.314D0 * TEMP)
        IF (EXPONENT .LT. -75.D0) THEN ! e^-75 < 3*10^-33: 0% kans anyway
        !write(500,*) "Large Exponent!", I
            KANS = 0.D0
            RV = 1.D0 ! Skip rand() voor cpu tijd besparing.
        ELSE IF(EXPONENT .GE. 0.D0) THEN ! Lager in energie, dus 100% kans
            KANS = 1.D0
            RV = 0.D0
        ELSE
            KANS = E ** EXPONENT
            RV = RAND()
            IF (KANS .GT. 1.D0) KANS = 1.D0
        END IF

        ! Bepaal if succesvol -> volgende config
        906 FORMAT(A, E20.10, ' / ', ES20.10)
        IF(RV .LE. KANS) THEN ! Succes!
            WRITE (*,906) "New solute accepted! (Prev / Post)", PRE_EN, POST_EN
            ISOK = .TRUE.
        ELSE ! Fail!
            WRITE (*,906) "New solute rejected! (Prev / Post)", PRE_EN, POST_EN
            ISOK = .FALSE.
        END IF

    END FUNCTION SOLUTE_METROPOLIS


    FUNCTION FIND(A1, A2, SOL, NATOM) RESULT(A)
    TYPE (VECTOR), DIMENSION(:), INTENT(IN)  :: SOL
        INTEGER:: A1
        INTEGER:: A2
        INTEGER:: A
        INTEGER:: NATOM
        INTEGER:: I

        A = 0
        DO I=1, NATOM
            IF (I .NE. A1 .AND. I .NE. A2) THEN
                IF (getDist(SOL(I), SOL(A1)) .LT. 1.90) THEN
                    A = I
                    RETURN
                END IF
            END IF
        END DO
        IF (A .EQ. 0) THEN
            WRITE (IOerr,*) "No match find in dihedral-find with", A1, A2
            STOP
        END IF
    END FUNCTION FIND

!====================================================================
!====================================================================

FUNCTION calcSol(SYM, SOL, WORKDIR) RESULT(NEWE)

    IMPLICIT NONE

    ! Input parameters
    TYPE (vector), DIMENSION(:), INTENT(IN) :: SOL ! absolute coords!
    CHARACTER*4, DIMENSION(:), INTENT(IN) :: SYM ! Atoomtypes
    CHARACTER*(*), INTENT(IN) :: WORKDIR

    ! Output parameters
    !DOUBLE PRECISION, INTENT(OUT) :: NEWE
    DOUBLE PRECISION :: NEWE

    ! Internal variables
    CHARACTER*500  :: GAUSS_SCRATCH
    CHARACTER*500  :: GAUSS_IN
    CHARACTER*500  :: GAUSS_OUT
    CHARACTER*500  :: GAUSS_ERR
    CHARACTER*2000 :: GAUSS_COMM
    CHARACTER*100  :: DUMMYSTRING
    INTEGER        :: IOSTATUS ! Check for EOF
    INTEGER        :: ThreadNum
    INTEGER        :: FI, FO         ! File handle numbers for input- and output file/pipe

    ! Get our OpenMP thread number. This is used to avoid directory conflicts
    ThreadNum = OMP_get_thread_num()
#ifdef DEBUG
    WRITE (*,"(A,I3.3,A)") "DEBUG: calcSol thread", ThreadNum, ": Entered calcSol"
#endif
    FI = IOstartrange + 2*ThreadNum
    FO = FI + 1

    ! Prepare strings with file names and commands to execute.
    WRITE(GAUSS_SCRATCH, "(A,'/work-',I3.3)")  TRIM(WORKDIR), ThreadNum
#if defined(USE_PIPES)
    WRITE(GAUSS_IN, "(A,'/input.pipe')")  TRIM(GAUSS_SCRATCH)
    WRITE(GAUSS_OUT,"(A,'/output.pipe')") TRIM(GAUSS_SCRATCH)
#else
    WRITE(GAUSS_IN, "(A,'/input.mop')")   TRIM(GAUSS_SCRATCH)
    WRITE(GAUSS_OUT,"(A,'/output.txt')")  TRIM(GAUSS_SCRATCH)
#endif
    WRITE(GAUSS_ERR,"(A,'/error_calcGaEn.log')") TRIM(GAUSS_SCRATCH)

!#if defined(USE_PIPES)

    GAUSS_COMM = "cd " // TRIM(GAUSS_SCRATCH) // "; " // &
                 "mopac " // TRIM(GAUSS_IN) // " 2> /dev/null; " // &
                 "grep 'FINAL HEAT OF FORMATION =' input.pipe.out" // &
                 " > " // TRIM(GAUSS_OUT)
!#endif
#ifdef DEBUG
    WRITE (*,"(A,I3.3,A)") "DEBUG: calcSol thread", ThreadNum, ": Want to execute: " // TRIM(GAUSS_COMM)
#endif


#if defined(USE_PIPES)
    ! Call Gaussian. Since the command string starts Gaussian in the background,
    ! the SYSTEM call should return immediately. GAussian will be waiting in the background
    ! for us to open the pipes.
#ifdef DEBUG
    WRITE(*,'(A,I3.3,A)') "DEBUG: calcSol thread ", ThreadNum, ": Starting Gaussian in calcSol (to run in the background)"
#endif
!    CALL SYSTEM( TRIM(GAUSS_COMM) // " &" )
#ifdef DEBUG
    WRITE(*,'(A,I3.3,A)') "DEBUG: calcSol thread ", ThreadNum, &
      &                   ": Finishing system call for Gaussian; it should be running in the background"
#endif
#endif

    ! Open the input file/pipe
    OPEN(UNIT=FI, FILE=TRIM(GAUSS_IN), ACTION='WRITE')
    CALL calcSol_input(SOL, SYM, FI)
    CLOSE(FI)

    CALL SYSTEM( TRIM(GAUSS_COMM))! // " &" )

#if !defined(USE_PIPES)
    ! Using files instead of pipes: It is now the moment to launch Gaussian.
#ifdef DEBUG
    WRITE(*,'(A,I3.3,A)') "DEBUG: calcSol thread ", ThreadNum, ": Starting Gaussian in calcSol"
#endif
    CALL SYSTEM(GAUSS_COMM)
#ifdef DEBUG
    WRITE(*,'(A,I3.3,A)') "DEBUG: calcSol thread ", ThreadNum, &
      &                   ": Finishing system call for Gaussian, Gaussian should be completed"
#endif
#endif

    ! Open the output file/pipe.
    OPEN(UNIT=FO, FILE=TRIM(GAUSS_OUT), STATUS='OLD', ACTION='READ')

    ! Now get the result from the output pipe or output file.
    READ (FO, "(A63,F12.4)", IOSTAT=IOSTATUS) DUMMYSTRING, NEWE ! lees resultaat in (A22,F14.12)
    IF (IOSTATUS .NE. 0 ) THEN
        WRITE(*,"(A,I3.3,A,I3)") "! calcSol thread ", ThreadNum, &
          &                      "  : Unexpected error/end-of-file while reading " // TRIM(GAUSS_OUT) // &
          &                      ", READ returned code ", IOSTATUS
!$OMP CRITICAL (CS_IOERR)
        ! Safeguarding writing to the error file with a critical section so that another thread
        ! doesn't mix its I/O, assuming it also uses a critical section of course.
        WRITE(IOerr,"(A)") "========================================"
        WRITE(IOerr, *) "Error @ calculating solute"
        WRITE(IOerr,"(A,I3.3,A)") &
                           "calcSol thread ", ThreadNum, " failed reading output data from the execution of " // &
                           TRIM(GAUSS_COMM) // ", Input was:"
        WRITE(IOerr,"(A)") "BEGIN INPUT below this line"
        CALL calcSol_input(SOL, SYM, IOerr)
        WRITE(IOerr,"(A)") "END INPUT above this line"
        WRITE(IOerr,"(A)") "========================================"
!$OMP END CRITICAL (CS_IOERR)
        !CALL ABORT()
        NEWE = 10000.D0
        !RETURN
    END IF
#if defined(DEBUG) && defined(USE_PIPES)
    WRITE(*,'(A, I3.3, A, A21, F20.10)') "DEBUG: calcSol thread ", ThreadNum, ": Read from pipe: ", DUMMYSTRING, EN
#endif

    CLOSE(FO)

#ifdef DEBUG
    WRITE (*,"(A,I3.3,A)") "DEBUG: calcSol thread ", ThreadNum, ": Exiting calcSol"
#endif

CONTAINS

    SUBROUTINE calcSol_input(MOL, SYM, FI)

        IMPLICIT NONE

        ! Input arguments
        TYPE (vector), DIMENSION(:), INTENT(IN) :: MOL ! absolute coords!
        CHARACTER*4, DIMENSION(:), INTENT(IN)   :: SYM ! Atoomtypes
        INTEGER, INTENT(IN) :: FI

        ! Other variables
        INTEGER :: K

        ! Output formats
        905 FORMAT(A, 3F16.8)

        ! Do the actual IO
        WRITE (FI,"(A)") 'PM6 1SCF THREADS=1 '
        WRITE (FI,"(A)") 'Dual calculation   '
        WRITE (FI,"(A)") 'yes dual           '

        ! Print Mol1
        DO K=1,size(MOL)
            WRITE (FI,905) sym(K), mol(K)%X, mol(K)%Y, mol(K)%Z
        END DO

        WRITE (FI,"(A)") '                                        '

    END SUBROUTINE calcSol_input

END FUNCTION calcSol
END MODULE solmod
