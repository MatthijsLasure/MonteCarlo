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

    ! FILES
    !======
    CHARACTER*500, DIMENSION(11)    :: FILES
    CHARACTER*500                   :: CONFILE, LJ_STEPS_TEMP, GA_STEPS_TEMP, ID_TEMP, SOLNAME ! Config
    CHARACTER*100                   :: DMSO_FILE, BOX_FILE, SOL_FILE, PARAM_FILE, PARSOL_FILE ! input files
    CHARACTER*100                   :: OUT_FILE, ERR_FILE, DUMP_FILE, SOLVSOLV_FILE, RESULT_FILE, SOLOUT_FILE ! output files

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
    INTEGER, DIMENSION(:,:), ALLOCATABLE     :: DIHOEK
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
    WRITE (*,*) "Fase 0 started!"
    WRITE (*,*) "Reading config..."

    ! Read command line
    CALL get_command_argument(1, CONFILE)

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
    WRITE (*,*) "Outputstream started in ", OUT_FILE
    WRITE (*,*) "Errorstream started in ", ERR_FILE
    OPEN(UNIT=501, FILE=OUT_FILE)
    OPEN(UNIT=500, FILE=ERR_FILE)

    WRITE (*,*) "Config loaded! Writing time..."

    CALL system_clock (START)
    WRITE(500,*) START

    ! Get seed for randgen.
    CALL init_random_seed()
    WRITE (*,"(A5,I20)") "SEED ", ISEED

    ! Laden van configuraties
!====================================================================
!====================================================================

    WRITE (*,*) "Loading data..."

    ! DMSO.txt: conformatie DMSO
    OPEN (UNIT=10, FILE=DMSO_FILE)
    READ (10, *) nDMSO ! Lees aantal atomen
    READ (10, *) ! Comment line
    ALLOCATE(DMSO(NDMSO)) ! Ken correcte groottes toe aan de arrays
    ALLOCATE(DMSO_sym(NDMSO))
    ALLOCATE(TABLE_DMSO(NDMSO, 3))
    DO I=1, NDMSO
        READ (10,*) DMSO_sym(I), DMSO(I)%X, DMSO(I)%Y, DMSO(I)%Z
    END DO
    READ (10,*) E_DMSO
    CLOSE(10)

    ! solute.txt: conformatie opgeloste molecule (sol)
    OPEN (UNIT=10, FILE=SOL_FILE)
    READ (10, *) nSol ! Lees aantal atomen
    READ (10, "(A)") SOLNAME ! Comment line
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(SOLUTE_OLD(NSOL))
    ALLOCATE(sol_sym(NSOL))
    ALLOCATE(TABLE_SOL(NSOL, 3))
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (10,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    READ (10,*) E_sol
    READ (10,*) NDIHOEK
    ALLOCATE(DIHOEK(NDIHOEK, 2))
    DO I=1, NDIHOEK
        READ (10,*) DIHOEK(I,1), DIHOEK(I,2)
    END DO
    CLOSE(10)

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=10, FILE=BOX_FILE)
    READ (10, *) BOXL ! Box grootte
    READ (10, *) nCoM ! Lees aantal moleculen
    READ (10, *) ! Box ID
    ALLOCATE(CoM(NCOM))
    ALLOCATE(hoek(NCOM))
    ALLOCATE(CoM_old(NCOM))
    ALLOCATE(hoek_old(NCOM))
    DO I= 1, NCOM
        READ(10,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO
    CLOSE(10)

    ! param.txt: parameters voor LJ etc
    OPEN (UNIT=10, FILE=PARAM_FILE)
    READ (10, *) nParam ! Aantal beschikbare parameters
    READ (10, *) ! skip comment line
    ALLOCATE(sym(NPARAM))
    ALLOCATE(Q(NPARAM))
    ALLOCATE(epsilon(NPARAM))
    ALLOCATE(sigma(NPARAM))
    DO I= 1,NPARAM
        READ (10,*) sym(I), Q(I), epsilon(I), sigma(I)
    END DO
    CLOSE(10)

    ! par_solute.txt: parameters voor LJ van het solute
    OPEN (UNIT=10, FILE=PARSOL_FILE)
    READ (10, *) NPARSOL
    READ (10, *) ! Comment line
    ALLOCATE(SOL_Q(NPARSOL))
    ALLOCATE(SOLPAR_SYM(NPARSOL))
    ALLOCATE(SOL_EPSILON(NPARSOL))
    ALLOCATE(SOL_SIGMA(NPARSOL))
    DO I= 1,NPARSOL
        READ (10,*) SOLPAR_SYM(I), SOL_EPSILON(I), SOL_SIGMA(I)
    END DO
    CLOSE(10)

    WRITE (*,*) "Done loading data!"

    WRITE (*,*) "BOXL DPOSMAX DPOSMIN DHOEKMAX DHOEKMIN"
    WRITE (*,*) BOXL, DPOSMAX, DPOSMIN, DHOEKMAX, DHOEKMIN
    IF (DOROTSOLV) WRITE (*,"(A,F13.10)") "Solvent rotation, max ", DROTSOLV

    BOXL2 = BOXL * 2.D0

!====================================================================
!====================================================================

    WRITE(*,*) "Working on assosiation tables..."
    CALL ASSOSIATE_DMSO(DMSO_SYM, SYM, Q, EPSILON, SIGMA, TABLE_DMSO)
    CALL ASSOSIATE_SOLUTE(SOL_SYM, SOLPAR_SYM, SOL_Q, SOL_EPSILON, SOL_SIGMA, TABLE_SOL)

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
    IF (DOROTSOLV) SOLUTE = SOLUTE_INIT(SOLUTE, DIHOEK, DROTSOLV)

    ! Initiële berekening ladingen solute
    !====================================
    WRITE (*,*) "Calculating partial charges on solute..."
    CALL DO_SOLUTE(SOL_SYM, SOLUTE,SOL_Q)

    DO I=1,NCOM
        CALL calculateLJ(I) ! Bereken alle energiën!
        !CALL calculateGA(I, 0)
    END DO

    TOTENG = 0.D0
    TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)
    PRE_ENG = TOTENG


    ! Dump energiën
    !tot = 0.0D
    OPEN(UNIT=10, FILE=SOLVSOLV_FILE)
    DO I=1,NCOM
        WRITE(10,*) solventsolvent(I,:)
    END DO
    CLOSE(10)


!====================================================================
!====================================================================


WRITE (*,*) "Dump file opened on ", DUMP_FILE
OPEN(UNIT=20, FILE=DUMP_FILE)
WRITE (20,*) BOXL
CALL DUMP(0)
CLOSE(20)

901 FORMAT(A12, 1X, A20, 1X, A20, 1X, A6, 1X, A6, 1X, A3, 1X, A6, 1X, A6, 1X, A6)
902 FORMAT(I12.12, 1X, ES20.10, 1X, ES20.10, 1X, F6.4, 1X, F6.4, 1X, I3.3, 1X, F6.4, 1X, F6.4, 1X, F6.5)
WRITE (501,901) "i", "TotEng", "TotEng_old", "kans", "rv", "rSolv", "pSuc", "ratio","dposmax"
WRITE (501,902) 0, TOTENG, TOTENG, 0.D0, 0.D0, 0, REAL(0) / real(1), 0.D0, DPOSMAX

CALL system_clock(START)
WRITE(500, *) START

!====================================================================
!====================================================================

    WRITE (*,*) "Fase 0 done!"
    WRITE (*,*) "Fase 1 started!"
    WRITE (*,*) "Entering Lennard-Jones loop..."
    OPEN(UNIT=20, FILE=DUMP_FILE, ACCESS="APPEND")

    ! Loop 1: LJ
    !===========
    loop_LJ: DO UNICORN=1,LJ_STEPS

        CALL MCINIT(UNICORN) ! Verander 1 molecule + check voor OOB & PB

        ! Bereken veranderde interacties
        CALL calculateLJ(RSOLV)

        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)

        IF (MOD(UNICORN, LJ_DUMP) .EQ. 0) THEN
            CALL DUMP(UNICORN)
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

    !CLOSE(20)
    WRITE (*,*) "Exiting Lennard-Jones loop..."
    WRITE (*,*) "Writing data..."

! Write box
OPEN(UNIT=10, FILE="backup.txt")
WRITE(10,*) BOXL
WRITE(10,*) NCOM
DO I=1,NCOM
    WRITE(10,*) CoM(I)%x, CoM(I)%y, CoM(I)%z, hoek(I)%x, hoek(I)%y, hoek(I)%z
END DO
CLOSE(10)

WRITE (*,*) "Fase 1 done!"
WRITE (*,*) "Fase 2 started!"

!====================================================================
!====================================================================

CALL system_clock(START)
!WRITE(500,*) START

IF (GA_STEPS .GT. 0) THEN
    DO I=1,NCOM
        CALL calculateGA(I, 0)
    END DO
END IF

TOTENG = 0.D0
TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)
IF(LJ_STEPS .EQ. 0) PRE_ENG = TOTENG

! Dump energiën
!tot = 0.0D
!OPEN(UNIT=10, FILE="gaussian.txt")
!DO I=1,NCOM
!    WRITE(10,*) solventsolvent(I,:)
!END DO
!CLOSE(10)

    WRITE (*,*) "Entering Gaussian loop..."

    ! Loop 2: Gaussian
    loop_Ga: DO UNICORN=1+LJ_STEPS,GA_STEPS+LJ_STEPS

        ! Doe MC
        CALL MCINIT(UNICORN)

        ! Bereken veranderde interacties
        CALL calculateGA(RSOLV, UNICORN)
        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)

        !OPEN(UNIT=20, FILE=dump_file, ACCESS="APPEND")
        IF (MOD(UNICORN, GA_DUMP) .EQ. 0) THEN
            CALL DUMP(UNICORN)
        END IF
        !CLOSE(20)

        ! Doe Metropolis
        CALL METROPOLIS(UNICORN, GA_NADJ, GA_NPRINT, REJECTED)

    END DO loop_Ga

    CLOSE(20)
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
    !write (500,*) START

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
OPEN (UNIT=10, FILE=SOLOUT_FILE)
WRITE (10, *) NSOL ! Lees aantal atomen
WRITE (10, *) trim(SOLNAME) ! Comment line
DO I=1, NSOL ! Lees de coördinaten uit
    WRITE (10,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
END DO
WRITE (10,*) E_SOL
WRITE (10,*) NDIHOEK
DO I=1, NDIHOEK
    WRITE (10,*) DIHOEK(I,1), DIHOEK(I,2)
END DO
CLOSE(10)


OPEN(UNIT=20, FILE=DUMP_FILE, ACCESS="APPEND")
CALL DUMP(UNICORN+1)
CLOSE(20)

WRITE (*,*) "Fase 2 done!"
WRITE (*,*) "Fase POST started!"

WRITE (*,*) "Writing data..."

! box.txt: plaatsen van de moleculen (CoM, hoek)
OPEN (UNIT=10, FILE=RESULT_FILE)
WRITE (10, *) BOXL ! Box grootte
WRITE (10, *) NCOM ! Lees aantal moleculen
WRITE (10, *) "BOX ", trim(ID_TEMP)
DO I= 1, NCOM
    WRITE(10,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
END DO
CLOSE(10)

!====================================================================
!====================================================================

WRITE (*,*) "Deallocating..."
DEALLOCATE(DMSO)
DEALLOCATE(DMSO_SYM)
DEALLOCATE(TABLE_DMSO)
DEALLOCATE(TABLE_SOL)
DEALLOCATE(SOLUTE)
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

CLOSE(501)
CLOSE(500)

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
        SOLUTE_OLD = SOLUTE ! Wordt nog niet gevariëerd
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
        !write(500,*) "Large Exponent!", I
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
            !CALL DUMP(I)
            WRITE (501,902) I, TOTENG, TOTENG_OLD, KANS, RV, RSOLV, REAL(NSUC) / real(I), RATIO, DPOSMAX
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
    LOGICAL, DIMENSION(NCOM+1) :: TOOFAR, CLIPPED
    LOGICAL :: MAYCONTINUE = .TRUE.

    ! Solvent solvent
    ! Notice: geen i loop: alleen molecule i is veranderd en moet opnieuw berekend worden
    SOLVENTSOLVENT(I,I) = 0.D0
    MOL1 = RotMatrix(CoM(I), DMSO, hoek(I))
    MAYCONTINUE = .TRUE.

    !J = 1
    preploop: DO J = 1, NCOM
        CLIPPED(J) = .FALSE.
        TOOFAR(J) = .FALSE.

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

            TEMPJ%X = CoM(J)%X + FLOAT(KMIN) * BOXL2
            TEMPJ%Y = CoM(J)%Y + FLOAT(LMIN) * BOXL2
            TEMPJ%Z = CoM(J)%Z + FLOAT(MMIN) * BOXL2

            !write(500,*) I, J, TEMPJ%X, TEMPJ%Y, TEMPJ%Z, CoM(J)%X, CoM(J)%Y, CoM(J)%Z, KMIN, LMIN, MMIN
            !write(500,*) I, J, getDist(CoM(I), CoM(J)), sqrt(RMIN), KMIN, LMIN, MMIN

            IF (RMIN .GT. 49.D0) THEN
                TOOFAR(J) = .TRUE.
            ELSE
                MOL2 = RotMatrix(TEMPJ, DMSO, hoek(J))
                CALL calcGa(I, J, MOL1, MOL2, DMSO_SYM, DMSO_SYM, CLIPPED(J), PROC)
                !IF(CLIPPED(j)) MayContinue = .FALSE.
            END IF

        END IF  ! Not with itself
    END DO preploop


    ! EXECUTE
    ! Begin of parallel loop
    !================================================================

    !$OMP PARALLEL
    !$OMP DO SCHEDULE(DYNAMIC) PRIVATE(En)
    EXEC: DO J=1,NCOM
        IF (J .NE. I .AND. .NOT. TooFar(J) ) THEN
            EN = 0.D0
            CALL execGa(I, J, EN)
        END IF
    END DO EXEC
    !$OMP END DO
    !$OMP END PARALLEL

    ! END OF PARALLEL LOOP
    !================================================================

        ! Read in
        grep: DO J=1,NCOM
            EN = 0.D0
            IF (J .NE. I .AND. .NOT. TooFar(J)) THEN
                CALL GREPIT(I, J, EN)

                EN = EN - E_DMSO - E_DMSO
                EN = EN * HARTREE2KJMOL

                SOLVENTSOLVENT(I,J) = EN
                SOLVENTSOLVENT(J,I) = EN
            ELSE
                SOLVENTSOLVENT(I,J) = 0.D0
                SOLVENTSOLVENT(J,I) = 0.D0
            END IF
        END DO GREP

        ! Solvent-solute
        CALL calcGa(I, 0, MOL1, SOLUTE, DMSO_SYM, SOL_SYM, Clipped(NCOM+1), PROC)
        CALL execGa(I, 0, EN)
        CALL grepit(I, 0, EN)
        EN = EN - E_SOL - E_DMSO
        EN = EN * HARTREE2KJMOL
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

SUBROUTINE dump(I)
        INTEGER :: I

        WRITE (20,*) NCOM
        WRITE (20,*) "Timestep: ", I
        DO J=1,NCOM
            WRITE(20, *) CoM(J)%X, CoM(J)%Y, CoM(J)%Z, hoek(J)%X, hoek(J)%Y, hoek(J)%Z
        END DO

END SUBROUTINE dump

!====================================================================
!====================================================================

! The end :)
END PROGRAM MonteCarlo
