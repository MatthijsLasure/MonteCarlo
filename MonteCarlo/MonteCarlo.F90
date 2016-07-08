!====================================================================
! MonteCarlo.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Doe een Monte-Carlo simulatie met 1 solute in DMSO
! Eerst werken met LJ voor ruwe benadering (eg 10.000 stappen)
! Dan verfijnen met Gaussian (eg 250 stappen)
!====================================================================


PROGRAM MonteCarlo

    USE MCconstants
    USE vector_class
    USE LennardJones
    USE Gaussian
    USE lib
    USE randgen
    USE readconfig
    USE solmod
    USE iso_fortran_env

    IMPLICIT NONE


    ! Variabelen
    !===========
    INTEGER             :: RUN_ID = 0 ! Welke run
    DOUBLE PRECISION    :: BOXL, BOXL2, TEMPERATURE ! box grootte, halve box
    INTEGER             :: I, J, UNICORN, ISEED
    INTEGER             :: LJ_STEPS, GA_STEPS ! Aantal stappen per loop
    INTEGER             :: RSOLV ! Geselecteerde molecule voor MC
    INTEGER             :: NSUC = 0 ! Aantal succesvolle MC
    INTEGER             :: LJ_NADJ, LJ_NPRINT, GA_NADJ, GA_NPRINT ! Aanpassen dposMax / printen om de n cycli
    INTEGER             :: LJ_DUMP, GA_DUMP
    INTEGER             :: NACCEPT = 0
    DOUBLE PRECISION    :: RV, KANS,DELTA = 0, EXPONENT ! Random variabele en toebehoren voor Metropolis
    DOUBLE PRECISION    :: RATIO ! Percentage succes in NADJ trials
    DOUBLE PRECISION    :: PADJ ! Hoeveel de dposmax aangepast mag worden
    LOGICAL             :: REJECTED, DOROTSOLV
    INTEGER             :: PROC ! Aantal processoren voor gaussian
    DOUBLE PRECISION    :: BOXSCALE = 0.9D0 ! Schalen van de box
    DOUBLE PRECISION    :: DROTSOLV
    CHARACTER(LEN=30)   :: DATE
    INTEGER             :: STATUS

    ! FILES
    !======
    CHARACTER*500, DIMENSION(11)    :: FILES
    CHARACTER*500                   :: CONFILE, LJ_STEPS_TEMP, GA_STEPS_TEMP, ID_TEMP, SOLNAME ! Config
    CHARACTER*100                   :: DMSO_FILE, BOX_FILE, SOL_FILE, PARAM_FILE, PARSOL_FILE ! input files
    CHARACTER*100                   :: OUT_FILE, ERR_FILE, DUMP_FILE, SOLVSOLV_FILE, RESULT_FILE, SOLOUT_FILE ! output files
    CHARACTER*500                   :: WORKDIR

    ! Energiën van de moleculen, bekomen via extern programma
    DOUBLE PRECISION                :: E_DMSO, E_SOL

    ! Vectoren voor moleculen & atoomtypes
    TYPE (vector), DIMENSION(:), ALLOCATABLE    :: DMSO, COM, SOLUTE, HOEK
    ! Variabelen voor de vorige run van MC
    TYPE (vector), DIMENSION(:), ALLOCATABLE    :: COM_OLD, SOLUTE_OLD, HOEK_OLD
    CHARACTER*4, DIMENSION(:), ALLOCATABLE      :: DMSO_SYM, SOL_SYM
    INTEGER                                     :: NDMSO, NCOM, NSOL, NPARAM, NPARSOL ! Aantal units
    ! DMSO: relatieve coördinaten voor de atomen
    ! CoM: Centre of Mass: locaties van de DMSO moleculen
    ! solute: conformatie van de solute
    ! Hoek: oriëntatie van de DMSO moleculen -> RotMatrix

    ! Array met bindingen die gebruikt mogen worden voor de dihedrale rotatie
    INTEGER, DIMENSION(:,:), ALLOCATABLE        :: DIHOEK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DROTSOLV_ARRAY
    INTEGER :: NDIHOEK

    ! Tijdelijke vectoren voor mol
    TYPE (vector), DIMENSION(:), ALLOCATABLE :: MOL1, MOL2

    ! Debug
!====================================================================
    INTEGER:: K, L                                                     !
    TYPE (vector), DIMENSION(10) :: ABSPOS                          !
    INTEGER:: START                                                 !
    LOGICAL:: DODEBUG = .FALSE.                                          !
!====================================================================

    ! output calc
    DOUBLE PRECISION                                :: EN ! Energie van een run
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: SOLVENTSOLVENT, SSOLD
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: ENERGY, EOLD
    DOUBLE PRECISION                                :: TOTENG = 0.D0, PRE_ENG, POST_ENG
    DOUBLE PRECISION                                :: TOTENG_OLD = 0.D0

    ! Arrays voor parameters van DMSO (Q, epsilon, sigma, mass)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: Q, EPSILON, SIGMA
    CHARACTER*4, DIMENSION(:), ALLOCATABLE          :: SYM
    ! Arrays voor paremeters van solute (epsilon, sigma)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: SOL_Q, SOL_EPSILON, SOL_SIGMA
    CHARACTER*4, DIMENSION(:), ALLOCATABLE          :: SOLPAR_SYM
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: TABLE_DMSO, TABLE_SOL

    ! Maximale verarndering bij MC
    DOUBLE PRECISION    :: DPOSMAX, DHOEKMAX, DPOSMIN, DHOEKMIN

    CALL FDATE(DATE)

    ! Config
!====================================================================
!====================================================================
    WRITE (*,*) "****************************************"
    WRITE (*,*) "* Welcome to the MC simulation of DMSO *"
    WRITE (*,*) "* Author: Matthijs Lasure              *"
    WRITE (*,*) "****************************************"
    WRITE (*,*) "Program initiated @ ", DATE
    WRITE (*,*) "Variable init Done!"

    ! Prepare the work directoy for the calls to Gaussian.
    CALL get_environment_variable("MonteCarlo_WORKDIR", value=WORKDIR, status=STATUS)
    IF ( STATUS > 0 ) THEN
        WORKDIR = "/dev/shm/MonteCarlo-exp"
        WRITE (*,*)  "WARNING: We suggest to set the environment variable MonteCarlo_WORKDIR "// &
          &          "to point ot the work directory for Gaussian"
    END IF
    WRITE (*,*) "Using the work directory " // TRIM(WORKDIR)

    CALL prepGaussian( WORKDIR )
    WRITE (*,*) "Work directory initialized"

    ! Start reading the config file procedure.
    WRITE (*,*) "Fase 0 started!"
    WRITE (*,*) "Reading config..."

    ! Read command line
    CALL get_command_argument(1, CONFILE)
    WRITE (*,*) "Reading config from: ", TRIM(CONFILE)

    ! Read from config.ini
    CALL rConfig(CONFILE, LJ_STEPS, GA_STEPS, ISEED, LJ_NADJ, LJ_NPRINT, GA_NADJ, &
    GA_NPRINT, LJ_DUMP, GA_DUMP, DPOSMAX, DPOSMIN, DHOEKMAX, DHOEKMIN, PADJ, PROC, &
    FILES, TEMPERATURE, DOROTSOLV, DROTSOLV)

    IF (LJ_DUMP .EQ. 0) LJ_DUMP = huge(LJ_DUMP)
    IF (GA_DUMP .EQ. 0) GA_DUMP = huge(GA_DUMP)

    ! Override stuff with command line
    IF (COMMAND_ARGUMENT_COUNT() .GT. 1) THEN ! Seriële modus

        CALL GET_COMMAND_ARGUMENT(2, LJ_STEPS_TEMP)
        CALL GET_COMMAND_ARGUMENT(3, GA_STEPS_TEMP)
        CALL GET_COMMAND_ARGUMENT(4, ID_TEMP)
        READ (LJ_STEPS_TEMP, *) LJ_STEPS
        READ (GA_STEPS_TEMP, *) GA_STEPS
        READ (ID_TEMP, *) RUN_ID

        WRITE (*,*) "-----------------------"
        WRITE (*,*) "Serial Modus requested!"
        WRITE (*,*) "ID: ", trim(ID_TEMP)
        WRITE (*,*) "LJ: ", trim(LJ_STEPS_TEMP)
        WRITE (*,*) "GA: ", trim(GA_STEPS_TEMP)
        WRITE (*,*) "-----------------------"

        ! Files

        ! INPUT
        BOX_FILE= trim(files(1)) // "." // trim(ID_TEMP) // ".in"
        SOL_FILE = trim(files(3)) // "." // trim(ID_TEMP) // ".in"
        !write(box_file, 910) files(1), ".", RUN_ID, ".in"

        ! OUTPUT
        OUT_FILE = trim(files(5)) // "." // trim(ID_TEMP) // ".txt"
        ERR_FILE = trim(files(6)) // "." // trim(ID_TEMP) // ".txt"
        DUMP_FILE = trim(files(7)) // "." // trim(ID_TEMP) // ".txt"
        SOLVSOLV_FILE = trim(files(8)) // "." // trim(ID_TEMP) // ".txt"
        RESULT_FILE = trim(files(9)) // "." // trim(ID_TEMP) // ".out"
        SOLOUT_FILE = trim(files(11)) // "." // trim(ID_TEMP) // ".out"

    ELSE ! No command line given
        ! INPUT
        BOX_FILE = files(1)

        ! OUTPUT
        SOL_FILE = files(3)
        OUT_FILE = files(5)
        ERR_FILE = files(6)
        DUMP_FILE = files(7)
        SOLVSOLV_FILE = files(8)
        RESULT_FILE = files(9)
        SOLOUT_FILE = files(11)
    END IF

    DHOEKMAX = DHOEKMAX * PI
    DHOEKMIN = DHOEKMIN * PI

    ! Standaard shit, altijd hetzelfde
    DMSO_FILE = files(2)
    PARAM_FILE = files(4)
    PARSOL_FILE = files(10)

    ! Print info
    WRITE(*,"(A, A)") "BOX    ", trim(BOX_FILE)
    WRITE(*,"(A, A)") "DMSO   ", trim(DMSO_FILE)
    WRITE(*,"(A, A)") "SOL    ", trim(SOL_FILE)
    WRITE(*,"(A, A)") "PARAM  ", trim(PARAM_FILE)
    WRITE(*,"(A, A)") "OUT    ", trim(OUT_FILE)
    WRITE(*,"(A, A)") "ERR    ", trim(ERR_FILE)
    WRITE(*,"(A, A)") "DUMP   ", trim(DUMP_FILE)
    WRITE(*,"(A, A)") "SOLV   ", trim(SOLVSOLV_FILE)
    WRITE(*,"(A, A)") "RES    ", trim(RESULT_FILE)
    WRITE(*,"(A, A)") "SOLOUT ", trim(SOLOUT_FILE)


    ! START ERR/OUT
    OPEN(UNIT=IOout, FILE=OUT_FILE)
    WRITE (*,*) "Outputstream started in ", OUT_FILE
    OPEN(UNIT=IOerr, FILE=ERR_FILE)
    WRITE (*,*) "Errorstream started in ", ERR_FILE
    OPEN(UNIT=IOdump, FILE=DUMP_FILE)
    WRITE (*,*) "Dump file opened on ", DUMP_FILE

    WRITE (*,*) "Config loaded!"

    ! Get seed for randgen.
    CALL init_random_seed()
    WRITE (*,"(A5,I20)") "SEED ", ISEED

    ! Laden van configuraties
!====================================================================
!====================================================================

    WRITE (*,*) "Loading data..."

    ! DMSO.txt: conformatie DMSO
    OPEN (UNIT=IOwork, FILE=DMSO_FILE)
    READ (IOwork, *) nDMSO ! Lees aantal atomen
    READ (IOwork, *) ! Comment line
    ALLOCATE(DMSO(NDMSO)) ! Ken correcte groottes toe aan de arrays
    ALLOCATE(DMSO_sym(NDMSO))
    ALLOCATE(TABLE_DMSO(NDMSO, 3))
    DO I=1, NDMSO
        READ (IOwork,*) DMSO_sym(I), DMSO(I)%X, DMSO(I)%Y, DMSO(I)%Z
    END DO
    READ (IOwork,*) E_DMSO
    CLOSE(IOwork)

    ! solute.txt: conformatie opgeloste molecule (sol)
    OPEN (UNIT=IOwork, FILE=SOL_FILE)
    READ (IOwork, *) nSol ! Lees aantal atomen
    READ (IOwork, "(A)") SOLNAME ! Comment line
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(SOLUTE_OLD(NSOL))
    ALLOCATE(sol_sym(NSOL))
    ALLOCATE(SOL_Q(NSOL))
    ALLOCATE(TABLE_SOL(NSOL, 3))
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (IOwork,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    READ (IOwork,*) E_sol
    READ (IOwork,*) NDIHOEK
    ALLOCATE(DIHOEK(NDIHOEK, 2))
    ALLOCATE(DROTSOLV_ARRAY(NDIHOEK))
    DO I=1, NDIHOEK
        READ (IOwork,"(I3, 1X, I3, 1X, F10.6)") DIHOEK(I,1), DIHOEK(I,2), DROTSOLV_ARRAY(I)
        IF (DROTSOLV_ARRAY(I) .EQ. 0.D0) DROTSOLV_ARRAY(I) = DROTSOLV
    END DO
    CLOSE(IOwork)


    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=IOwork, FILE=BOX_FILE)
    READ (IOwork, *) BOXL ! Box grootte
    READ (IOwork, *) nCoM ! Lees aantal moleculen
    READ (IOwork, *)      ! Box ID
    ALLOCATE(CoM(NCOM))
    ALLOCATE(hoek(NCOM))
    ALLOCATE(CoM_old(NCOM))
    ALLOCATE(hoek_old(NCOM))
    DO I= 1, NCOM
        READ(IOwork,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO
    CLOSE(IOwork)

    ! param.txt: parameters voor LJ etc
    OPEN (UNIT=IOwork, FILE=PARAM_FILE)
    READ (IOwork, *) nParam ! Aantal beschikbare parameters
    READ (IOwork, *) ! skip comment line
    ALLOCATE(sym(NPARAM))
    ALLOCATE(Q(NPARAM))
    ALLOCATE(epsilon(NPARAM))
    ALLOCATE(sigma(NPARAM))
    DO I= 1,NPARAM
        READ (IOwork,*) sym(I), Q(I), epsilon(I), sigma(I)
    END DO
    CLOSE(IOwork)

    ! par_solute.txt: parameters voor LJ van het solute
    OPEN (UNIT=IOwork, FILE=PARSOL_FILE)
    READ (IOwork, *) NPARSOL
    READ (IOwork, *) ! Comment line
    ALLOCATE(SOLPAR_SYM(NPARSOL))
    ALLOCATE(SOL_EPSILON(NPARSOL))
    ALLOCATE(SOL_SIGMA(NPARSOL))
    DO I= 1,NPARSOL
        READ (IOwork,*) SOLPAR_SYM(I), SOL_EPSILON(I), SOL_SIGMA(I)
    END DO
    CLOSE(IOwork)

    WRITE (*,*) "Done loading data!"

    WRITE (*,*) "BOXL DPOSMAX DPOSMIN DHOEKMAX DHOEKMIN"
    WRITE (*,*) BOXL, DPOSMAX, DPOSMIN, DHOEKMAX, DHOEKMIN
    IF (DOROTSOLV) THEN
        WRITE (*,"(A,F13.10)") "Solvent rotation, max ", DROTSOLV
        DO I=1,NDIHOEK
            WRITE (*,"('Bond ', I3, ' - ', I3, ' MAX ', F10.6)") DIHOEK(I,1), DIHOEK(I,2), DROTSOLV_ARRAY(I)
        END DO
    END IF

    BOXL2 = BOXL * 2.D0
    E_SOL = E_SOL * HARTREE2KJMOL
    E_DMSO = E_DMSO * HARTREE2KJMOL

    WRITE (*,*) "E_SOL  (kJ/mol): ", E_SOL
    WRITE (*,*) "E_DMSO (kJ/mol): ", E_DMSO

!====================================================================
!====================================================================

    WRITE(*,*) "Working on assosiation tables..."
    CALL ASSOSIATE_DMSO(DMSO_SYM, SYM, Q, EPSILON, SIGMA, TABLE_DMSO)
    CALL ASSOSIATE_SOLUTE(SOL_SYM, SOLPAR_SYM, SOL_EPSILON, SOL_SIGMA, TABLE_SOL)

    ! Initiële berekening interacties
    !================================

    WRITE (*,*) "Initial energy calculations of the system"

    ! Arrays
    ALLOCATE(mol1(NDMSO))
    ALLOCATE(mol2(NDMSO))
    ALLOCATE(solventsolvent(NCOM, NCOM))
    ALLOCATE(ssold(NCOM, NCOM))
    ALLOCATE(energy(NCOM))
    ALLOCATE(eold(NCOM))

    ! Roteer de solute
    SOLUTE_OLD = SOLUTE
    IF (DOROTSOLV) SOLUTE = SOLUTE_INIT(SOLUTE, SOL_SYM, DIHOEK, DROTSOLV_ARRAY)

    ! Initiële berekening ladingen solute
    !====================================
    IF (LJ_STEPS .GT. 0) THEN
        WRITE (*,*) "Calculating partial charges on solute..."
        CALL DO_SOLUTE(SOL_SYM, SOLUTE, TABLE_SOL(:,1), WORKDIR)
        WRITE (*,*) "Done calculating"
    END IF

    IF (LJ_STEPS .GT. 0 .OR. LJ_STEPS .EQ. -5) THEN
        DO I=1,NCOM
            CALL calculateLJ(I) ! Bereken alle energiën!
            !CALL calculateGA(I, 0)
        END DO

        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)
        PRE_ENG = TOTENG


        ! Dump energiën
        !tot = 0.0D
        OPEN(UNIT=IOwork, FILE=SOLVSOLV_FILE)
        DO I=1,NCOM
            WRITE(IOwork,*) solventsolvent(I,:)
        END DO
        CLOSE(IOwork)
    END IF

!====================================================================
!====================================================================


    !WRITE (*,*) "Dump file opened on ", DUMP_FILE
    !OPEN(UNIT=IOdump, FILE=DUMP_FILE)
    WRITE (IOdump,*) BOXL
    CALL DUMP(0, IOdump)
    !CLOSE(IOdump)

    901 FORMAT(A12, 1X, A20, 1X, A20, 1X, A6, 1X, A6, 1X, A3, 1X, A6, 1X, A6, 1X, A6)
    902 FORMAT(I12.12, 1X, ES20.10, 1X, ES20.10, 1X, F6.4, 1X, F6.4, 1X, I3.3, 1X, F6.4, 1X, F6.4, 1X, F6.5)
    WRITE (IOout,901) "i", "TotEng", "TotEng_old", "kans", "rv", "rSolv", "pSuc", "ratio","dposmax"
    IF(LJ_STEPS .GT. 0) WRITE (IOout,902) 0, TOTENG, TOTENG, 0.D0, 0.D0, 0, REAL(0) / real(1), 0.D0, DPOSMAX

!====================================================================
!====================================================================

    WRITE (*,*) "Fase 0 done!"
    WRITE (*,*) "Fase 1 started!"
    WRITE (*,*) "Entering Lennard-Jones loop..."
    !OPEN(UNIT=IOdump, FILE=DUMP_FILE, ACCESS="APPEND")

    ! Loop 1: LJ
    !===========
    loop_LJ: DO UNICORN=1,LJ_STEPS

        CALL MCINIT(UNICORN) ! Verander 1 molecule + check voor OOB & PB

        ! Bereken veranderde interacties
        CALL calculateLJ(RSOLV)

        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)

        IF (MOD(UNICORN, LJ_DUMP) .EQ. 0) THEN
            CALL DUMP(UNICORN, IOdump)
        END IF

        ! Check for invalid shit
        IF(TOTENG > HUGE(TOTENG)) THEN
            WRITE (0,*) "TOTENG IS INFINITY @", UNICORN, TOTENG
            TOTENG = HUGE(TOTENG)-1
        END IF
        IF(TOTENG .NE. TOTENG) THEN
            WRITE (0,*) "TOTENG IS NaN @", UNICORN, TOTENG
            TOTENG = HUGE(TOTENG)-1
        END IF

        ! Doe Metropolis
        CALL METROPOLIS(UNICORN, LJ_NADJ, LJ_NPRINT, REJECTED)

    END DO loop_LJ

    !CLOSE(IOdump)
    WRITE (*,*) "Exiting Lennard-Jones loop..."
    WRITE (*,*) "Writing data..."

    ! Write box
    OPEN(UNIT=IOwork, FILE="backup.txt")
    WRITE(IOwork,*) BOXL
    WRITE(IOwork,*) NCOM
    DO I=1,NCOM
        WRITE(IOwork,*) CoM(I)%x, CoM(I)%y, CoM(I)%z, hoek(I)%x, hoek(I)%y, hoek(I)%z
    END DO
    CLOSE(IOwork)
    
    WRITE (*,*) "Fase 1 done!"
    
    !====================================================================
    !====================================================================

    WRITE (*,*) "Fase 2 started!"

    
    IF (GA_STEPS .GT. 0 .OR. LJ_STEPS .EQ. -1) THEN
        DO I=1,NCOM
            CALL calculateGA(I, 0)
        END DO
        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)
        WRITE (IOout,902) 0, TOTENG, TOTENG, 0.D0, 0.D0, 0, 0.D0, 0.D0, DPOSMAX
    END IF
    
    IF(LJ_STEPS .EQ. 0) PRE_ENG = TOTENG
    
    ! Dump energiën
    !tot = 0.0D
    OPEN(UNIT=IOwork, FILE="gaussian.txt")
    DO I=1,NCOM
        WRITE(IOwork,*) solventsolvent(I,:)
    END DO
    CLOSE(IOwork)

    WRITE (*,*) "Entering Gaussian loop..."

    ! Loop 2: Gaussian
    loop_Ga: DO UNICORN=1+LJ_STEPS,GA_STEPS+LJ_STEPS

        ! Doe MC
        CALL MCINIT(UNICORN)

        ! Bereken veranderde interacties
        CALL calculateGA(RSOLV, UNICORN)
        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)

        !OPEN(UNIT=IOdump, FILE=dump_file, ACCESS="APPEND")
        IF (MOD(UNICORN, GA_DUMP) .EQ. 0) THEN
            CALL DUMP(UNICORN, IOdump)
        END IF
        !CLOSE(IOdump)

        ! Doe Metropolis
        CALL METROPOLIS(UNICORN, GA_NADJ, GA_NPRINT, REJECTED)

    END DO loop_Ga

    !CLOSE(IOdump)
    WRITE (*,*) "Exiting Gaussian loop..."

    !TEMPDUMP
    DO I=1,NCOM
    MOL1 = RotMatrix(COM(I), DMSO, HOEK(I))
        DO K = 1, NDMSO
            DO L = 1, NSOL
                IF(getDist(MOL1(K), SOLUTE(L)) .LT. 1.8D0) THEN
                    CALL CalcLJ_SV(MOL1, SOLUTE, DMSO_SYM, TABLE_DMSO, TABLE_SOL, EN, BOXL, BOXL2)
                    WRITE (*,*) I, 0, K, L, getDist(MOL1(K), SOLUTE(L)), EN
                END IF
            END DO
        END DO
        DO J=1, NCOM
        IF (J .NE. I) THEN
            MOL2 = RotMatrix(COM(J), DMSO, HOEK(J))
            DO K = 1, NDMSO
                DO L = 1, NDMSO
                    IF(getDist(MOL1(K), MOL2(L)) .LT. 1.8D0) THEN
                        CALL calcLJ(MOL1, MOL2, DMSO_SYM, TABLE_DMSO, EN, BOXL, BOXL2)
                        WRITE (*,*) I, J, K, L, getDist(MOL1(K), MOL2(L)), EN
                    END IF
                END DO
            END DO
        END IF
        END DO
    END DO

    CALL system_clock(START)
    !write (IOerr,*) START

!====================================================================
!====================================================================
    ! Wegschrijven resultaten
    
    IF(DOROTSOLV) THEN
        POST_ENG = TOTENG
        IF(.NOT. SOLUTE_METROPOLIS(SOLUTE, PRE_ENG, POST_ENG, TEMPERATURE)) THEN
            SOLUTE = SOLUTE_OLD
        END IF
    END IF
    
    ! solute.txt: conformatie opgeloste molecule (sol)
    OPEN (UNIT=IOwork, FILE=SOLOUT_FILE)
    WRITE (IOwork, *) NSOL ! Lees aantal atomen
    WRITE (IOwork, *) trim(SOLNAME) ! Comment line
    DO I=1, NSOL ! Lees de co�rdinaten uit
        WRITE (IOwork,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    WRITE (IOwork,*) E_SOL / HARTREE2KJMOL
    WRITE (IOwork,*) NDIHOEK
    DO I=1, NDIHOEK
        WRITE (IOwork,"(I3, 1X, I3, 1X, F10.6)") DIHOEK(I,1), DIHOEK(I,2), DROTSOLV_ARRAY(I)
    END DO
    CLOSE(IOwork)
    
    CALL DUMP(UNICORN+1, IOdump)
    
    WRITE (*,*) "Fase 2 done!"
    WRITE (*,*) "Fase POST started!"
    
    WRITE (*,*) "Writing data..."
    
    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=IOwork, FILE=RESULT_FILE)
    WRITE (IOwork, *) BOXL ! Box grootte
    WRITE (IOwork, *) NCOM ! Lees aantal moleculen
    WRITE (IOwork, *) "BOX ", trim(ID_TEMP)
    DO I= 1, NCOM
        WRITE(IOwork,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO
    CLOSE(IOwork)
    
!====================================================================
!====================================================================
    
    WRITE (*,*) "Deallocating..."
    DEALLOCATE(DMSO)
    DEALLOCATE(DMSO_SYM)
    DEALLOCATE(TABLE_DMSO)
    DEALLOCATE(TABLE_SOL)
    DEALLOCATE (DIHOEK)
    DEALLOCATE(SOLUTE)
    DEALLOCATE (SOLUTE_OLD)
    DEALLOCATE(SOL_SYM)
    DEALLOCATE(CoM)
    DEALLOCATE(hoek)
    DEALLOCATE(COM_OLD)
    DEALLOCATE(HOEK_OLD)
    DEALLOCATE(SYM)
    DEALLOCATE(Q)
    DEALLOCATE(EPSILON)
    DEALLOCATE(SIGMA)
    DEALLOCATE(SOL_Q)
    DEALLOCATE(SOLPAR_SYM)
    DEALLOCATE(SOL_EPSILON)
    DEALLOCATE(SOL_SIGMA)
    DEALLOCATE(MOL1)
    DEALLOCATE(MOL2)
    DEALLOCATE(SOLVENTSOLVENT)
    DEALLOCATE(SSOLD)
    DEALLOCATE(ENERGY)
    DEALLOCATE(EOLD)
    
    
!====================================================================
!====================================================================
    
    CLOSE(IOout)
    CLOSE(IOerr)
    CLOSE(IOdump)
    
!====================================================================
!====================================================================
    CALL cleanGaussian(WORKDIR)

    CALL FDATE(DATE)
    
    WRITE (*,*) "Fase POST done!"
    WRITE (*,*) "We're done here. Signing off!"
    WRITE (*,*) "Program finished @ ", DATE

CONTAINS

!====================================================================
!====================================================================

SUBROUTINE MCINIT(I)

        INTEGER :: I

        ! Verhuis oude vars naar de _old vars
        COM_OLD = COM
        HOEK_OLD = HOEK
        SOLUTE_OLD = SOLUTE ! Wordt nog niet gevarieerd
        TOTENG_OLD = TOTENG
        EOLD = ENERGY
        SSOLD = SOLVENTSOLVENT

        ! Doe MC
        RSOLV = INT(RAND() * NCOM) + 1 ! Willekeurige DMSO molecule
        IF( RSOLV .EQ. 31) RSOLV = 30
        COM(RSOLV) = CoM(RSOLV) + randVec(DPOSMAX)
        HOEK(RSOLV) = hoek(RSOLV) + randVecHoek(DHOEKMAX)


        ! Check of hoeken nog binnen [-Pi, +Pi] liggen
        IF (.TRUE.) THEN
        HOEK(RSOLV)%X = HOEK(RSOLV)%X - TAU * DINT(HOEK(RSOLV)%X / PI)
        HOEK(RSOLV)%Y = HOEK(RSOLV)%Y - TAU * DINT(HOEK(RSOLV)%Y / PI)
        HOEK(RSOLV)%Z = HOEK(RSOLV)%Z - TAU * DINT(HOEK(RSOLV)%Z / PI)

        ! Periodic boundaries => Computer Simulation of Liquids, p30
        COM(RSOLV)%X = COM(RSOLV)%X - BOXL2 * DINT(COM(RSOLV)%X / BOXL)
        COM(RSOLV)%Y = COM(RSOLV)%Y - BOXL2 * DINT(COM(RSOLV)%Y / BOXL)
        COM(RSOLV)%Z = COM(RSOLV)%Z - BOXL2 * DINT(COM(RSOLV)%Z / BOXL)
        END IF

END SUBROUTINE MCINIT

!====================================================================
!====================================================================

SUBROUTINE METROPOLIS(I, NADJ, NPRINT, REJECTED)

        INTEGER :: I
        INTEGER :: NADJ, NPRINT
        LOGICAL, INTENT(OUT) :: REJECTED

        902 FORMAT(I12.12, 1X, ES20.10, 1X, ES20.10, 1X, F6.4, 1X, F6.4, 1X, I3.3, 1X, F6.4, 1X, F6.4, 1X, F6.5)

        DELTA = TOTENG - TOTENG_OLD
        EXPONENT = -1.D0 * DELTA  * 1000.D0 / (8.314D0 * TEMPERATURE)
        IF (EXPONENT .LT. -75.D0) THEN ! e^-75 < 3*10^-33: 0% kans anyway
        !write(IOerr,*) "Large Exponent!", I
            KANS = 0.D0
            RV = 1.D0 ! Skip rand() voor cpu tijd besparing
        ELSE IF(EXPONENT .GE. 0.D0) THEN ! Lager in energie, dus 100% kans
            KANS = 1.D0
            RV = 0.D0
        ELSE
            KANS = E ** EXPONENT
            RV = RAND()
            IF (KANS .GT. 1.D0) KANS = 1.D0
        END IF

        ! Bepaal if succesvol -> volgende config
        IF(RV .LE. KANS) THEN ! Succes!
            NSUC = NSUC + 1
            NACCEPT = NACCEPT + 1
            REJECTED = .FALSE.
        ELSE ! Fail!
            REJECTED = .TRUE.
        END IF

        RATIO = REAL(NACCEPT) / real(NADJ)

        ! Prints every NPRINT times
        IF(mod(I, NPRINT) .EQ. 0) THEN
            !CALL DUMP(I, IOdump)
            WRITE (IOout,902) I, TOTENG, TOTENG_OLD, KANS, RV, RSOLV, REAL(NSUC) / real(I), RATIO, DPOSMAX
        END IF

        ! Pas de dposMax aan indien ratio =/= 50%
        IF (mod(I, NADJ) .EQ. 0) THEN
            IF (RATIO .GT. 0.5) THEN
                DPOSMAX = DPOSMAX * (1.D0 + PADJ)
                DHOEKMAX = DHOEKMAX * (1.D0 + PADJ)
            ELSEIF (RATIO .LT. 0.5D0) THEN
                DPOSMAX = DPOSMAX * (1.D0 - PADJ)
                DHOEKMAX = DHOEKMAX * (1.D0 - PADJ)
            END IF
            DPOSMAX = DMAX1(DPOSMAX, DPOSMIN)
            DHOEKMAX = DMAX1(DHOEKMAX, DHOEKMIN)
            NACCEPT = 0
        END IF

        IF(REJECTED) THEN
            ! Verhuis oude vars terug naar de nieuwe
            COM = COM_OLD
            HOEK = HOEK_OLD
            TOTENG = TOTENG_OLD
            ENERGY = EOLD ! Resetten!!
            SOLVENTSOLVENT = SSOLD
        END IF

END SUBROUTINE METROPOLIS
!====================================================================
!====================================================================

SUBROUTINE calculateLJ(I)

    INTEGER:: I
    INTEGER:: J

    TYPE (vector) :: R

    ! Solvent solvent
    ! Notice: geen i loop: alleen molecule i is veranderd en moet opnieuw berekend worden
    SOLVENTSOLVENT(I,I) = 0.D0
    MOL1 = RotMatrix(CoM(I), DMSO, hoek(I))

    DO J=1,NCOM
        IF ( J .NE. I) THEN
            ! Check lengte, met MIC!
            R = CoM(J) - CoM(I)
            R%X = R%X - BOXL2 * DINT(R%X / BOXL)
            R%Y = R%Y - BOXL2 * DINT(R%Y / BOXL)
            R%Z = R%Z - BOXL2 * DINT(R%Z / BOXL)

            IF (length(R) .GT. 7.D0) THEN
                EN = 0.D0
            ELSE
                !MOL2 = RotMatrix(TEMPJ, DMSO, hoek(J))
                MOL2 = RotMatrix(CoM(J), DMSO, hoek(J))
                CALL calcLJ(MOL1, MOL2, DMSO_SYM, TABLE_DMSO, EN, BOXL, BOXL2)
            END IF
            SOLVENTSOLVENT(I,J) = EN
            SOLVENTSOLVENT(J,I) = EN
        END IF
    END DO

    ! Solvent-solute
    CALL CalcLJ_SV(MOL1, SOLUTE, DMSO_SYM, TABLE_DMSO, TABLE_SOL, EN, BOXL, BOXL2)
    ENERGY(I) = EN

END SUBROUTINE calculateLJ

!====================================================================
!====================================================================

SUBROUTINE calculateGA(I, LOOPNR)

    INTEGER:: I, J, LOOPNR ! gevraagde moleculen, welke loop
    INTEGER :: K, L, M, KMIN = 2, LMIN = 2, MMIN = 2 ! Minimal Image Convention
    DOUBLE PRECISION :: R, RMIN, EN = 0
    TYPE (vector) :: TEMPJ
    LOGICAL :: TOOFAR, CLIPPED
    LOGICAL :: MAYCONTINUE = .TRUE.

    ! Solvent solvent
    ! Notice: geen i loop: alleen molecule i is veranderd en moet opnieuw berekend worden
    SOLVENTSOLVENT(I,I) = 0.D0
    MOL1 = RotMatrix(CoM(I), DMSO, hoek(I))
    MAYCONTINUE = .TRUE.

    !J = 1
    ! EXECUTE
    ! Begin of parallel loop
    !================================================================
#ifdef DEBUG
    PRINT *, "Entering parallel section."
#endif

OPEN (5414, FILE="debug.xyz", access='append')
    !$OMP PARALLEL default(none) PRIVATE(EN,MOL2, RMIN, KMIN, LMIN, MMIN, R, K, L, M, TEMPJ, TOOFAR, CLIPPED) &
    !$OMP SHARED(NCOM, SOLVENTSOLVENT, I, COM, BOXL2, DMSO, HOEK, NDMSO, MOL1, LOOPNR, WORKDIR, PROC, DMSO_SYM, E_DMSO, E_SOL)
    !$OMP DO
    EXEC: DO J = 1, NCOM
        SOLVENTSOLVENT(I,J) = 0.D0
        SOLVENTSOLVENT(J,I) = 0.D0

        CLIPPED = .FALSE.
        TOOFAR = .FALSE.

        IF (J .NE. I) THEN ! Not with itself
            ! MINIMIZE DISTANCE
            RMIN = huge(EN)
            K_LOOP: DO K=-1,1
                DO L=-1,1
                    DO M=-1,1
                        TEMPJ%X = CoM(J)%X + FLOAT(K) * BOXL2
                        TEMPJ%Y = CoM(J)%Y + FLOAT(L) * BOXL2
                        TEMPJ%Z = CoM(J)%Z + FLOAT(M) * BOXL2

                        R = getDistSq(CoM(I), TEMPJ) ! Check
                        IF(R .LT. RMIN) THEN
                            RMIN = R
                            KMIN = K
                            LMIN = L
                            MMIN = M
                        END IF
                    END DO
                END DO
            END DO K_LOOP

            TEMPJ%X = CoM(J)%X + DBLE(KMIN) * BOXL2
            TEMPJ%Y = CoM(J)%Y + DBLE(LMIN) * BOXL2
            TEMPJ%Z = CoM(J)%Z + DBLE(MMIN) * BOXL2

            !write(IOerr,*) I, J, TEMPJ%X, TEMPJ%Y, TEMPJ%Z, CoM(J)%X, CoM(J)%Y, CoM(J)%Z, KMIN, LMIN, MMIN
            !write(IOerr,*) I, J, getDist(CoM(I), CoM(J)), sqrt(RMIN), KMIN, LMIN, MMIN

            IF (RMIN .GT. 49.D0) THEN
                TOOFAR = .TRUE.
                !WRITE (IOerr, *) "Molecules too far - ", I, J, " @ " , RMIN
            ELSE
                MOL2 = RotMatrix(TEMPJ, DMSO, hoek(J))
!                MOL2 = RotMatrix(CoM(J), DMSO, hoek(J))
!                DO K=1,NDMSO
!                    MOL2(K)%X = MOL2(K)%X + DBLE(KMIN) * BOXL2
!                    MOL2(K)%Y = MOL2(K)%Y + DBLE(LMIN) * BOXL2
!                    MOL2(K)%Z = MOL2(K)%Z + DBLE(MMIN) * BOXL2
!                END DO

                DO K=1,NDMSO
                    DO L=1,NDMSO
                        R = getDistSq(MOL1(K), MOL2(L))
                        IF (R .LT. 0.25) THEN ! Clipped
                         CLIPPED = .TRUE.
                         WRITE (IOerr, *) "Clipped @", LOOPNR, I, J, K, L, R
                        END IF
                    END DO
                END DO

                IF (CLIPPED) THEN
                    EN = 1000.D0
                ELSE
                    CALL calcGaEn(I, J, MOL1, MOL2, DMSO_SYM, DMSO_SYM, EN, WORKDIR)
                    if(en .EQ. 10000.D0) then
                        WRITE (IOerr,*) "+++ ", LOOPNR, I, J, " +++"
                        WRITE (IOerr,*) KMIN, LMIN, MMIN
                        WRITE(IOerr,*) COM(I)%X, COM(I)%Y, COM(I)%Z, HOEK(I)%X, HOEK(I)%Y, HOEK(I)%Z
                        WRITE(IOerr,*) TEMPJ%X, TEMPJ%Y, TEMPJ%Z, HOEK(J)%X, HOEK(J)%Y, HOEK(J)%Z
                        WRITE (IOerr,*) "+++ ", LOOPNR, " +++"

                        905 FORMAT(A, 3F16.8)

                        !$OMP CRITICAL (CS_IOERR)
                        WRITE (5414, *) 30
                        WRITE (5414, *) LOOPNR, I, J, KMIN, LMIN, MMIN
                        ! Print Mol1
                        DO K=1,size(MOL1)
                            WRITE (5414,905) DMSO_SYM(K), mol1(K)%X, mol1(K)%Y, mol1(K)%Z
                        END DO

                        ! Print mol2
                        DO K=1,size(MOL2)
                            WRITE (5414,905) DMSO_SYM(K), mol2(K)%X, mol2(K)%Y, mol2(K)%Z
                        END DO
                        MOL2 = RotMatrix(CoM(J), DMSO, hoek(J))
                        !MOL2 = RotMatrix(TEMPJ, DMSO, hoek(J))
                        DO K=1,size(MOL2)
                            WRITE (5414,905) DMSO_SYM(K), mol2(K)%X + DBLE(KMIN) * BOXL2, mol2(K)%Y + DBLE(LMIN) * BOXL2&
                            , mol2(K)%Z + DBLE(MMIN) * BOXL2
                        END DO
                        !$OMP END CRITICAL (CS_IOERR)
                       CALL ABORT()

                    END IF
                    !EN = EN - (E_DMSO + E_DMSO) * HARTREE2KJMOL
                    EN = EN - E_DMSO - E_DMSO
                    !EN = EN * HARTREE2KJMOL
                END IF

                SOLVENTSOLVENT(I,J) = EN
                SOLVENTSOLVENT(J,I) = EN
                !IF(CLIPPED(j)) MayContinue = .FALSE.
            END IF

        END IF  ! Not with itself
    END DO  EXEC
    !$OMP END DO
    !$OMP END PARALLEL
    close(5414)
    ! END OF PARALLEL LOOP
    !================================================================
#ifdef DEBUG
    PRINT *, "Exiting parallel section."
#endif

    ! Solvent-solute
    !CALL calcGa(I, 0, MOL1, SOLUTE, DMSO_SYM, SOL_SYM, Clipped(NCOM+1), PROC)
    !CALL execGa(I, 0, EN)
    !CALL grepit(I, 0, EN)
    CALL calcGaEn(I, 0, MOL1, SOLUTE, DMSO_SYM, SOL_SYM, EN, WORKDIR)
    !EN = EN - (E_SOL + E_DMSO) * HARTREE2KJMOL
    EN = EN - E_SOL - E_DMSO
    !EN = EN * HARTREE2KJMOL
    ENERGY(I) = EN

END SUBROUTINE CALCULATEGA

!====================================================================
!====================================================================

FUNCTION calcEnergy(ENERGY, SOLVENTSOLVENT) RESULT(OUT)
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ENERGY
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: SOLVENTSOLVENT
    DOUBLE PRECISION:: OUT
    INTEGER:: I, J

    ! Totale E
    OUT = 0.D0

    DO I=1, NCOM
        OUT = OUT + energy(I) ! solv - solu
        DO J=I+1, NCOM
            OUT = OUT + solventsolvent(I,J)
        END DO
    END DO

END FUNCTION calcEnergy

!====================================================================
!====================================================================

SUBROUTINE dump(I,IOunit)
        INTEGER :: I
        INTEGER :: IOunit

        WRITE (IOunit,*) NCOM
        WRITE (IOunit,*) "Timestep: ", I
        DO J=1,NCOM
            WRITE(IOunit, *) CoM(J)%X, CoM(J)%Y, CoM(J)%Z, hoek(J)%X, hoek(J)%Y, hoek(J)%Z
        END DO

END SUBROUTINE dump

!====================================================================
!====================================================================

! The end :)
END PROGRAM MonteCarlo
