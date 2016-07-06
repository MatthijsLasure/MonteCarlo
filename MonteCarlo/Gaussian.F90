!====================================================================
! interactions.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Bereken energie via LJ of QC
!====================================================================
MODULE gaussian

    USE MCconstants, ONLY : IOstartrange, IOerr
    USE vector_class
    USE lib
    USE OMP_LIB

    CONTAINS

!====================================================================

SUBROUTINE prepGaussian( WORKDIR )

    IMPLICIT NONE

    ! Input parameter
    CHARACTER*(*), INTENT(IN) :: WORKDIR

    ! Global work variables
    CHARACTER*1000 :: GLOB_COMM

    ! Per thread work variables
    INTEGER :: ThreadNum
    CHARACTER*500  :: THREAD_DIR
    CHARACTER*1000 :: THREAD_COMM

    ! First make sure the work directory exists
    GLOB_COMM = "mkdir -p " // TRIM(WORKDIR)
#ifdef DEBUG
    WRITE (*,'(A)') "DEBUG prepGaussian: Executing " // TRIM(GLOB_COMM)
#endif
    CALL SYSTEM(GLOB_COMM)

    ! Now the per thread part: Create thread work directory and
    ! the pipes. We do this in an OpenMP parallel section to
    ! make sure we got the right number of threads.
    !$OMP PARALLEL PRIVATE(ThreadNum,THREAD_DIR,THREAD_COMM)
    ThreadNum = OMP_get_thread_num()
    WRITE(THREAD_DIR,"(A,'/work-',I3.3)") TRIM(WORKDIR), ThreadNum
#if defined(USE_PIPES)
    THREAD_COMM = "/bin/mkdir -p " // TRIM(THREAD_DIR)  // " ; " // &
      &           "cd " // TRIM(THREAD_DIR)             // " ; " // &
!      &           "[ -e input.pipe ]  || mkfifo input.pipe ; " // &
      &           "touch input.pipe; " // &
      &           "touch output.pipe"
!      &           "[ -e output.pipe ] || mkfifo output.pipe"
#else
    THREAD_COMM = "mkdir -p " // TRIM(THREAD_DIR)
#endif
#ifdef DEBUG
    WRITE (*,'(A)') "DEBUG prepGaussian: Executing " // TRIM(THREAD_COMM)
#endif
    CALL SYSTEM(THREAD_COMM)

    ! Make the script
    WRITE(THREAD_COMM, "(A,'/do.sh')") THREAD_DIR
    OPEN(10, FILE=THREAD_COMM)
    WRITE (10, "(A)") "#!/bin/bash"
    WRITE (10, "(A)") "trap 'exit' SIGRTMAX-10"
    WRITE (10, "(A)") "echo $$ >> ../PIDS.txt"
    WRITE (10, "(A)") "while : do g09 < input.pipe > output.pipe; done"
    CLOSE(10)

    WRITE (THREAD_COMM, "(A)") "cd " // TRIM(THREAD_DIR) // " ; " // &
     &                         "chmod 777 do.sh; ./do.sh"

    CALL SYSTEM(THREAD_COMM)

    !$OMP END PARALLEL

#ifdef DEBUG
    ! DEEBUG: List the work directory.
    GLOB_COMM="ls -l " // TRIM(WORKDIR) // "/work-*"
    WRITE (*,'(A)') "DEBUG prepGaussian: Executing " // TRIM(GLOB_COMM)
    CALL SYSTEM(GLOB_COMM)
#endif

END SUBROUTINE prepGaussian

!====================================================================
!====================================================================

SUBROUTINE cleanGaussian(WORKDIR)

    IMPLICIT NONE

    ! Input parameter
    CHARACTER*(*), INTENT(IN) :: WORKDIR

    ! Global work variables
    CHARACTER*250 :: GLOB_COMM
    INTEGER       :: PID, STAT
    CHARACTER*500 :: FILE

    FILE = TRIM(WORKDIR) // "PIDS.txt"
    OPEN(10, FILE=FILE)
    DO
        READ(10,'(I5.5)', IOSTAT=STAT) PID
        IF (STAT .LT. 0) EXIT
        CALL KILL(PID, 54)
    END
    CLOSE(10)



    ! Remove all work-* subdirectories
    GLOB_COMM = "/bin/rm -rf " // TRIM(WORKDIR) // "/work-*"
#ifdef DEBUG
    WRITE (*,'(A)') "DEBUG cleanGaussian: Executing " // TRIM(GLOB_COMM)
#endif
    CALL SYSTEM(GLOB_COMM)

END SUBROUTINE cleanGaussian

!====================================================================
!====================================================================

SUBROUTINE calcGaEn(I, J, MOL1, MOL2, SYM1, SYM2, EN, WORKDIR)

    IMPLICIT NONE

    ! Input parameters
    TYPE (vector), DIMENSION(:), INTENT(IN) :: MOL1, MOL2 ! absolute coords!
    CHARACTER*4, DIMENSION(:), INTENT(IN) :: SYM1, SYM2 ! Atoomtypes
    INTEGER, INTENT(IN) :: I, J
    CHARACTER*(*), INTENT(IN) :: WORKDIR

    ! Output parameters
    DOUBLE PRECISION, INTENT(OUT) :: EN

    ! Internal variables
    CHARACTER*500  :: GAUSS_SCRATCH
    CHARACTER*500  :: GAUSS_IN
    CHARACTER*500  :: GAUSS_OUT
    CHARACTER*500  :: GAUSS_ERR
    CHARACTER*2000 :: GAUSS_COMM
    CHARACTER*100  :: DUMMYSTRING
    CHARACTER*16   :: STR_I, STR_J
    INTEGER        :: IOSTATUS ! Check for EOF
    INTEGER        :: ThreadNum
    INTEGER        :: FI, FO         ! File handle numbers for input- and output file/pipe

    ! Get our OpenMP thread number. This is used to avoid directory conflicts
    ThreadNum = OMP_get_thread_num()
#ifdef DEBUG
    WRITE (*,"(A,I3.3,A)") "DEBUG: calcGaEn thread", ThreadNum, ": Entered calcGaEn"
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
                 "mopac " // TRIM(GAUSS_IN) // " &> /dev/null; " // &
                 "grep 'FINAL HEAT OF FORMATION =' input.pipe.out" // &
                 " > " // TRIM(GAUSS_OUT)
!#endif
#ifdef DEBUG
    WRITE (*,"(A,I3.3,A)") "DEBUG: calcGaEn thread", ThreadNum, ": Want to execute: " // TRIM(GAUSS_COMM)
#endif


#if defined(USE_PIPES)
    ! Call Gaussian. Since the command string starts Gaussian in the background,
    ! the SYSTEM call should return immediately. GAussian will be waiting in the background
    ! for us to open the pipes.
#ifdef DEBUG
    WRITE(*,'(A,I3.3,A)') "DEBUG: calcGaEn thread ", ThreadNum, ": Starting Gaussian in calcGaEn (to run in the background)"
#endif
!    CALL SYSTEM( TRIM(GAUSS_COMM) // " &" )
#ifdef DEBUG
    WRITE(*,'(A,I3.3,A)') "DEBUG: calcGaEn thread ", ThreadNum, &
      &                   ": Finishing system call for Gaussian; it should be running in the background"
#endif
#endif

    ! Open the input file/pipe
    OPEN(UNIT=FI, FILE=TRIM(GAUSS_IN), ACTION='WRITE')
    CALL calcGaEn_input(MOL1, MOL2, SYM1, SYM2, FI)
    CLOSE(FI)

    CALL SYSTEM( TRIM(GAUSS_COMM))! // " &" )

#if !defined(USE_PIPES)
    ! Using files instead of pipes: It is now the moment to launch Gaussian.
#ifdef DEBUG
    WRITE(*,'(A,I3.3,A)') "DEBUG: calcGaEn thread ", ThreadNum, ": Starting Gaussian in calcGaEn"
#endif
    CALL SYSTEM(GAUSS_COMM)
#ifdef DEBUG
    WRITE(*,'(A,I3.3,A)') "DEBUG: calcGaEn thread ", ThreadNum, &
      &                   ": Finishing system call for Gaussian, Gaussian should be completed"
#endif
#endif

    ! Open the output file/pipe.
    OPEN(UNIT=FO, FILE=TRIM(GAUSS_OUT), STATUS='OLD', ACTION='READ')

    ! Now get the result from the output pipe or output file.
    READ (FO, "(A63,F12.4)", IOSTAT=IOSTATUS) DUMMYSTRING, EN ! lees resultaat in (A22,F14.12)
    REWIND(FO)
    READ(FO, *) DUMMYSTRING
    WRITE (*,*) "calcGaEn thread ", ThreadNum, " found interaction energy ", EN
    WRITE (*,"(A,I3.3,A,A)") "calcGaEn thread ", ThreadNum, "string: ", DUMMYSTRING
    IF (IOSTATUS .NE. 0 ) THEN
        WRITE(*,"(A,I3.3,A,I3)") "calcGaEn thread ", ThreadNum, &
          &                      "  : Unexpected error/end-of-file while reading " // TRIM(GAUSS_OUT) // &
          &                      ", READ returned code ", IOSTATUS
!$OMP CRITICAL (CS_IOERR)
        ! Safeguarding writing to the error file with a critical section so that another thread
        ! doesn't mix its I/O, assuming it also uses a critical section of course.
        WRITE(IOerr,"(A)") "========================================"
        WRITE(IOerr, *) "Error @ ", I, " | ", J
        WRITE(IOerr,"(A,I3.3,A)") &
                           "calcGaEn thread ", ThreadNum, " failed reading output data from the execution of " // &
                           TRIM(GAUSS_COMM) // ", Input was:"
        WRITE(IOerr,"(A)") "BEGIN INPUT below this line"
        CALL calcGaEn_input(MOL1, MOL2, SYM1, SYM2, IOerr)
        WRITE(IOerr,"(A)") "END INPUT above this line"
        WRITE(IOerr,"(A)") "========================================"
!$OMP END CRITICAL (CS_IOERR)
        !CALL ABORT()
        EN = 10000.D0
        !RETURN
    END IF
#if defined(DEBUG) && defined(USE_PIPES)
    WRITE(*,'(A, I3.3, A, A21, F20.10)') "DEBUG: calcGaEn thread ", ThreadNum, ": Read from pipe: ", DUMMYSTRING, EN
#endif

    CLOSE(FO)

#ifdef DEBUG
    WRITE (*,"(A,I3.3,A)") "DEBUG: calcGaEn thread ", ThreadNum, ": Exiting calcGaEn"
#endif

CONTAINS

    SUBROUTINE calcGaEn_input(MOL1, MOL2, SYM1, SYM2, FI)

        IMPLICIT NONE

        ! Input arguments
        TYPE (vector), DIMENSION(:), INTENT(IN) :: MOL1, MOL2 ! absolute coords!
        CHARACTER*4, DIMENSION(:), INTENT(IN)   :: SYM1, SYM2 ! Atoomtypes
        INTEGER, INTENT(IN) :: FI

        ! Other variables
        INTEGER :: K

        ! Output formats
        905 FORMAT(A, 3F16.8)

        ! Do the actual IO
        WRITE (FI,"(A)") 'PM6 1SCF         '
        WRITE (FI,"(A)") 'Dual calculation '
        WRITE (FI,"(A)") 'yes dual         '

        ! Print Mol1
        DO K=1,size(MOL1)
            WRITE (FI,905) sym1(K), mol1(K)%X, mol1(K)%Y, mol1(K)%Z
        END DO

        ! Print mol2
        DO K=1,size(MOL2)
            WRITE (FI,905) sym2(K), mol2(K)%X, mol2(K)%Y, mol2(K)%Z
        END DO


        WRITE (FI,"(A)") '                                        '

    END SUBROUTINE calcGaEn_input

END SUBROUTINE calcGaEn

!====================================================================
!====================================================================

SUBROUTINE DO_SOLUTE(SOL_SYM, SOL, SOL_Q, WORKDIR)

    IMPLICIT NONE

    CHARACTER*4, DIMENSION(:)       :: SOL_SYM
    TYPE (VECTOR), DIMENSION(:)     :: SOL
    DOUBLE PRECISION, DIMENSION(:)  :: SOL_Q
    CHARACTER*(*), INTENT(IN)       :: WORKDIR

    INTEGER        :: I
    CHARACTER*20   :: NCHAR, NCHAR2
    !CHARACTER*1000                  :: COMMAND1, COMMAND2, COMMAND3
    CHARACTER*1000 :: COMMAND
    CHARACTER*500  :: GAUSS_SCRATCH
    CHARACTER*500  :: GAUSS_IN
    CHARACTER*500  :: GAUSS_OUT
    CHARACTER*500  :: GAUSS_ERR
    INTEGER        :: ThreadNum
    INTEGER        :: FI, FO         ! File handle numbers for input- and output file/pipe
    INTEGER        :: IOSTATUS

    ! Get our OpenMP thread number. This is used to avoid directory conflicts
    ThreadNum = OMP_get_thread_num()
#ifdef DEBUG
    WRITE (*,"(A,I3.3,A)") "DEBUG: Thread", ThreadNum, ": Entered DO_SOLUTE"
#endif
    FI = IOstartrange + 2*ThreadNum
    FO = FI + 1

    ! Prepare strings with file names and commands to execute.
    WRITE(GAUSS_SCRATCH, "(A,'/work-',I3.3)") TRIM(WORKDIR), ThreadNum
#if defined(USE_PIPES)
    WRITE(GAUSS_IN, "(A,'/input.pipe')")  TRIM(GAUSS_SCRATCH)
    WRITE(GAUSS_OUT,"(A,'/output.pipe')") TRIM(GAUSS_SCRATCH)
#else
    WRITE(GAUSS_IN, "(A,'/solute_charge.com')")        TRIM(GAUSS_SCRATCH)
    WRITE(GAUSS_OUT,"(A,'/solute_charge_output.log')") TRIM(GAUSS_SCRATCH)
#endif
    WRITE(GAUSS_ERR,"(A,'/error_DO_SOLUTE.log')") TRIM(GAUSS_SCRATCH)

    WRITE(NCHAR, "(I0)")  SIZE(SOL)
    WRITE(NCHAR2, "(I0)") SIZE(SOL) + 2

    ! Prepare the commands to call Gaussian
    COMMAND = "export GAUSS_SCRDIR=" // TRIM(GAUSS_SCRATCH) // "; " // &
      &       "g09 <" // TRIM(GAUSS_IN) // " 2>" // TRIM(GAUSS_ERR) // " | " // &
      &       "grep -A" // TRIM(NCHAR2) // " 'Charges from ESP fit' | " // &
      &       "tail -" // TRIM(NCHAR) // " >" // TRIM(GAUSS_OUT)
#ifdef DEBUG
    WRITE(*,"(A)") "DEBUG: DO_SOLUTE: Will start Gaussian through: " // TRIM(COMMAND)
#endif

#if defined(USE_PIPES)
    ! Using pipes: Start Gaussian in the background so that it is waiting for the pipes
    ! to be opened by the Fortran program.
#ifdef DEBUG
    WRITE(*,"(A)") "DEBUG: DO_SOLUTE: Calling SYSTEM (Gaussian in the background)"
#endif
!    CALL SYSTEM( TRIM(COMMAND))! // " &" )
#ifdef DEBUG
    WRITE(*,"(A)") "DEBUG: DO_SOLUTE: Returned from SYSTEM (Gaussian in the background)"
#endif
#endif

    ! Open the input file/pipe
    OPEN(FI, FILE=GAUSS_IN, ACTION='WRITE')
    CALL solute_input(SOL_SYM, SOL, SOL_Q, FI)
    CLOSE(FI)
#if defined(DEBUG) && !defined(USE_PIPES)
    WRITE(*,"(A)") "DEBUG: DO_SOLUTE: Input file for Gaussian:"
    CALL SYSTEM( 'cat <' // TRIM(GAUSS_IN) )
    WRITE(*,"(A)") "DEBUG: DO_SOLUTE: END DEBUG INFO input file for Gaussian"
#endif

CALL SYSTEM( TRIM(COMMAND))! // " &" )

#if !defined(USE_PIPES)
    ! Not using pipes: We now run Gaussian.
#ifdef DEBUG
    WRITE(*,"(A)") "DEBUG: DO_SOLUTE: Calling SYSTEM (Gaussian, wait till completed)"
#endif
    CALL SYSTEM(COMMAND)
#ifdef DEBUG
    WRITE(*,"(A)") "DEBUG: DO_SOLUTE: Returned from SYSTEM (Gaussian completed)"
#endif
#endif

    ! Open the output file/pipe.
    OPEN(FO, FILE=GAUSS_OUT, STATUS='OLD', ACTION='READ')

    DO I = 1, SIZE(SOL)
        READ (FO, "(A12, F9.7)",IOSTAT=IOSTATUS) NCHAR, SOL_Q(I)
        IF (IOSTATUS .NE. 0 ) THEN
            WRITE(*,"(A,I3.3)") "DO_SOLUTE: Unexpected error/end-of-file while reading " // TRIM(GAUSS_OUT) // &
              &                 ", READ returned code ", IOSTATUS
!$OMP CRITICAL (CS_IOERR)
            ! Safeguarding writing to the error file with a critical section so that another thread
            ! doesn't mix its I/O, assuming it also uses a critical section of course.
            WRITE(IOerr,"(A)") "========================================"
            WRITE(IOerr,"(A,I3.3,A)") &
                           "DO_SOLUTE thread ", ThreadNum, " failed reading output data from the execution of " // &
                           TRIM(COMMAND) // ", Input was:"
            WRITE(IOerr,"(A)") "BEGIN INPUT below this line"
            CALL solute_input(SOL_SYM, SOL, SOL_Q, IOerr)
            WRITE(IOerr,"(A)") "END INPUT above this line"
            WRITE(IOerr,"(A)") "========================================"
!$OMP END CRITICAL (CS_IOERR)
            !CALL ABORT()
        END IF
#if defined(DEBUG) && defined(USE_PIPES)
        WRITE(*,'(A, A12, F10.7)') "DEBUG: DO_SOLUTE: Read from pipe: ", NCHAR, SOL_Q(I)
#endif
    END DO

    CLOSE(FO)

#ifdef DEBUG
    WRITE (*,"(A,I3.3,A)") "DEBUG: Thread", ThreadNum, ": Exiting DO_SOLUTE"
#endif

CONTAINS

    SUBROUTINE solute_input(SOL_SYM, SOL, SOL_Q, FI)

        IMPLICIT NONE

        ! Input arguments
        CHARACTER*4, DIMENSION(:)       :: SOL_SYM
        TYPE (VECTOR), DIMENSION(:)     :: SOL
        DOUBLE PRECISION, DIMENSION(:)  :: SOL_Q
        INTEGER, INTENT(IN)             :: FI

        ! Other variables
        CHARACTER*4 :: GAUSS_PROCS
        INTEGER     :: I

        ! Initialisations
        WRITE(GAUSS_PROCS,'(I4)') OMP_get_num_procs()

        ! Preamble
        WRITE (FI,"(A)") '%nproc=' // ADJUSTL(GAUSS_PROCS) // '           '
        WRITE (FI,"(A)") '%mem=12GB                                       '
        WRITE (FI,"(A)") '%CHK=solute_charge.chk                          '
        WRITE (FI,"(A)") '#P B3LYP/6-31G* POP=(CHELPG,DIPOLE, READRADII)  '
        WRITE (FI,"(A)") '                                                '
        WRITE (FI,"(A)") 'SOLUTE CHARGE CALCULATION                       '
        WRITE (FI,"(A)") '                                                '
        WRITE (FI,"(A)") '0 1                                             '

        ! Data
        DO I = 1, SIZE(SOL)
            WRITE (FI, *) SOL_SYM(I), SOL(I)%X, SOL(I)%Y, SOL(I)%Z
        END DO

        ! Termination line
        WRITE (FI,"(A)") '                                                '
        WRITE (FI,"(A)") 'Br 2.44                                         '
        WRITE (FI,"(A)") '                                                '

    END SUBROUTINE solute_input

END SUBROUTINE DO_SOLUTE

!====================================================================

END MODULE gaussian
