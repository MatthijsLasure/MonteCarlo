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
    USE OMP_LIB
    !use iso_fortran_env

    IMPLICIT NONE


    ! Variabelen
    !===========
    INTEGER             :: RUN_ID = 0 ! Welke run
    DOUBLE PRECISION    :: BOXL, BOXL2 ! box grootte, halve box
    INTEGER             :: I, J, UNICORN, ISEED
    INTEGER             :: LJ_STEPS, GA_STEPS ! Aantal stappen per loop
    INTEGER             :: RSOLV ! Geselecteerde molecule voor MC
    INTEGER             :: NSUC = 0 ! Aantal succesvolle MC
    INTEGER             :: LJ_NADJ, LJ_NPRINT, GA_NADJ, GA_NPRINT ! Aanpassen dposMax / printen om de n cycli
    INTEGER             :: NACCEPT = 0
    DOUBLE PRECISION    :: RV, KANS,DELTA = 0, EXPONENT ! Random variabele en toebehoren voor Metropolis
    DOUBLE PRECISION    :: RATIO ! Percentage succes in NADJ trials
    DOUBLE PRECISION    :: PADJ ! Hoeveel de dposmax aangepast mag worden
    DOUBLE PRECISION    :: BETA ! p = EXP(-BETA * DELTA / (RT)
    LOGICAL             :: REJECTED
    INTEGER             :: PROC ! Aantal processoren voor gaussian
    DOUBLE PRECISION    :: BOXSCALE = 0.9D0 ! Schalen van de box

    ! FILES
    !======
    CHARACTER*500, DIMENSION(10)    :: files
    CHARACTER*500                   :: confile, LJ_STEPS_TEMP, GA_STEPS_TEMP, ID_TEMP ! Config
    CHARACTER*100                   :: dmso_file, box_file, sol_file, param_file, parsol_file ! input files
    CHARACTER*100                   :: out_file, err_file, dump_file, solvsolv_file, result_file ! output files

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

    ! Tijdelijke vectoren voor mol
    TYPE (vector), DIMENSION(:), ALLOCATABLE :: MOL1, MOL2

    ! Debug
!====================================================================
    INTEGER:: K                                                     !
    TYPE (vector), DIMENSION(10) :: ABSPOS                          !
    INTEGER:: START                                                 !
LOGICAL:: DODEBUG = .FALSE.                                          !
!====================================================================

    ! output calc
    DOUBLE PRECISION                                :: EN ! Energie van een run
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: SOLVENTSOLVENT, SSOLD
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: ENERGY, EOLD
    DOUBLE PRECISION                                :: TOTENG = 0.D0
    DOUBLE PRECISION                                :: TOTENG_OLD = 0.D0

    ! Arrays voor parameters van DMSO (Q, epsilon, sigma, mass)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: Q, EPSILON, SIGMA, MASS
    CHARACTER*4, DIMENSION(:), ALLOCATABLE          :: SYM
    ! Arrays voor paremeters van solute (epsilon, sigma)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: SOL_Q, SOL_EPSILON, SOL_SIGMA
    CHARACTER*4, DIMENSION(:), ALLOCATABLE          :: SOLPAR_SYM
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: TABLE_DMSO, TABLE_SOL

    ! Maximale verarndering bij MC
    DOUBLE PRECISION    :: DPOSMAX, DHOEKMAX



    ! Config
!====================================================================
!====================================================================
    WRITE (*,*) "******************************************************************"
    WRITE (*,*) "* Welcome to the MC simulation of DMSO                           *"
    WRITE (*,*) "* Author: Matthijs Lasure, Matthijs.Lasure@student.uantwerpen.be *"
    WRITE (*,*) "******************************************************************"
    WRITE (*,*) "Variable init Done!"
    WRITE (*,*) "Fase 0 started!"
    WRITE (*,*) "Reading config..."

    ! Read command line
    CALL get_command_argument(1, confile)

    ! Read from config.ini
    CALL rConfig(confile, BOXL, LJ_STEPS, Ga_STEPS, iseed, DoDebug, LJ_nadj, LJ_nprint, GA_nadj, &
    GA_nprint, dposmax, dhoekmax, padj, beta, proc, files)

    ! BOXL = OBSOLETE: wordt ingelezen uit box.txt

    ! Override stuff with command line
    IF (COMMAND_ARGUMENT_COUNT() .GT. 1) THEN ! Seriële modus

        CALL GET_COMMAND_ARGUMENT(2, LJ_STEPS_TEMP)
        CALL GET_COMMAND_ARGUMENT(3, GA_STEPS_TEMP)
        CALL GET_COMMAND_ARGUMENT(4, ID_TEMP)
        read (LJ_STEPS_TEMP, *) LJ_STEPS
        read (GA_STEPS_TEMP, *) GA_STEPS
        read (ID_TEMP, *) RUN_ID

        WRITE (*,*) "------------------------"
        WRITE (*,*) "Serial Modus detected!"
        WRITE (*,*) "ID: ", trim(ID_TEMP)
        WRITE (*,*) "LJ: ", trim(LJ_STEPS_TEMP)
        WRITE (*,*) "GA: ", trim(GA_STEPS_TEMP)
        WRITE (*,*) "------------------------"

        ! Files

        ! INPUT
        box_file= trim(files(1)) // "." // trim(ID_TEMP) // ".in"
        !write(box_file, 910) files(1), ".", RUN_ID, ".in"

        ! OUTPUT
        out_file = trim(files(5)) // "." // trim(ID_TEMP) // ".txt"
        err_file = trim(files(6)) // "." // trim(ID_TEMP) // ".txt"
        dump_file = trim(files(7)) // "." // trim(ID_TEMP) // ".txt"
        solvsolv_file = trim(files(8)) // "." // trim(ID_TEMP) // ".txt"
        result_file = trim(files(9)) // "." // trim(ID_TEMP) // ".out"

    ELSE ! No command line given
        ! INPUT
        box_file = files(1)

        ! OUTPUT
        out_file = files(5)
        err_file = files(6)
        dump_file = files(7)
        solvsolv_file = files(8)
        result_file = files(9)
    END IF

    DHOEKMAX = DHOEKMAX * PI

    ! Standaard shit, altijd hetzelfde
    dmso_file = files(2)
    sol_file = files(3)
    param_file = files(4)
    parsol_file = files(10)

    ! START ERR/OUT
    WRITE (*,*) "Outputstream started in ", out_file
    WRITE (*,*) "Errorstream started in ", err_file
    OPEN(UNIT=501, FILE=out_file)
    OPEN(UNIT=500, FILE=err_file)

    WRITE (*,*) "Config loaded! Writing time..."

    CALL system_clock (START)
    WRITE(500,*) START
    CALL SRAND(REAL(iseed)) ! Prime randgen

    ! Laden van configuraties
!====================================================================
!====================================================================

    WRITE (*,*) "Loading data..."

    ! DMSO.txt: conformatie DMSO
    OPEN (UNIT=10, FILE=dmso_file)
    READ (10, *) nDMSO ! Lees aantal atomen
    ALLOCATE(DMSO(NDMSO)) ! Ken correcte groottes toe aan de arrays
    ALLOCATE(DMSO_sym(NDMSO))
    ALLOCATE(TABLE_DMSO(NDMSO, 3))
    DO I=1, NDMSO
        READ (10,*) DMSO_sym(I), DMSO(I)%X, DMSO(I)%Y, DMSO(I)%Z
    END DO
    READ (10,*) E_DMSO
    CLOSE(10)

    ! solute.txt: conformatie opgeloste molecule (sol)
    OPEN (UNIT=10, FILE=sol_file)
    READ (10, *) nSol ! Lees aantal atomen
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(sol_sym(NSOL))
    ALLOCATE(TABLE_SOL(NSOL, 3))
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (10,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    READ (10,*) E_sol
    CLOSE(10)

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=10, FILE=box_file)
    READ (10, *) BOXL ! Box grootte
    READ (10, *) nCoM ! Lees aantal moleculen
    ALLOCATE(CoM(NCOM))
    ALLOCATE(hoek(NCOM))
    ALLOCATE(CoM_old(NCOM))
    ALLOCATE(hoek_old(NCOM))
    DO I= 1, NCOM
        READ(10,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO
    CLOSE(10)

    BOXL2 = BOXL * 2.D0

    ! param.txt: parameters voor LJ etc
    OPEN (UNIT=10, FILE=param_file)
    READ (10, *) nParam ! Aantal beschikbare parameters
    READ (10, *) ! skip comment line
    ALLOCATE(sym(NPARAM))
    ALLOCATE(Q(NPARAM))
    ALLOCATE(epsilon(NPARAM))
    ALLOCATE(sigma(NPARAM))
    ALLOCATE(mass(NPARAM))
    DO I= 1,NPARAM
        READ (10,*) sym(I), Q(I), epsilon(I), sigma(I), mass(I)
    END DO
    CLOSE(10)

    ! par_solute.txt: parameters voor LJ van het solute
    OPEN (UNIT=10, FILE=parsol_file)
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

!====================================================================
!====================================================================

    ! Initiële berekening ladingen solute
    !====================================
    WRITE (*,*) "Calculating partial charges on solute..."
    CALL DO_SOLUTE(SOL_SYM, SOLUTE,SOL_Q)

    FLUSH(5)
    FLUSH(6)

    write(*,*) "Working on assosiation tables..."
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

    DO I=1,NCOM
        CALL calculateLJ(I) ! Bereken alle energiën!
    END DO

    TOTENG = 0.D0
    TOTENG =  calcEnergy(ENERGY, SOLVENTSOLVENT)


    ! Dump energiën
    !tot = 0.0D
    OPEN(UNIT=10, FILE=solvsolv_file)
    DO I=1,NCOM
        WRITE(10,*) solventsolvent(I,:)
    END DO
    CLOSE(10)


!====================================================================
!====================================================================


WRITE (*,*) "Dump file opened on ", dump_file
OPEN(UNIT=20, FILE=dump_file)
CALL DUMP(0)
CLOSE(20)

901 FORMAT(A12, 1X, A20, 1X, A20, 1X, A6, 1X, A6, 1X, A3, 1X, A6, 1X, A6, 1X, A6)
902 FORMAT(I12.12, 1X, ES20.10, 1X, ES20.10, 1X, F6.4, 1X, F6.4, 1X, I3.3, 1X, F6.4, 1X, F6.4, 1X, F6.5)
WRITE (501,901) "i", "TotEng", "TotEng_old", "kans", "rv", "rSolv", "pSuc", "ratio","dposmax"
WRITE (501,902) 0, TOTENG, 0.D0, 0.D0, 0.D0, 0, REAL(0) / real(1), 0.D0, dposmax

CALL system_clock(start)
!write(500, *) start

!====================================================================
!====================================================================

    WRITE (*,*) "Fase 0 done!"
    WRITE (*,*) "Fase 1 started!"
    WRITE (*,*) "Entering Lennard-Jones loop..."

    ! Loop 1: LJ
    !===========
    loop_LJ: DO UNICORN=1,LJ_STEPS

        CALL MCINIT(UNICORN) ! Verander 1 molecule + check voor OOB & PB

        ! Bereken veranderde interacties
        CALL calculateLJ(RSOLV)
        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)

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

    WRITE (*,*) "Exiting Lennard-Jones loop..."
    WRITE (*,*) "Writing data..."

! Write box
OPEN(UNIT=10, FILE="backup.txt")
WRITE(10,*) boxl
WRITE(10,*) nCoM
DO I=1,nCoM
    WRITE(10,*) CoM(I)%x, CoM(I)%y, CoM(I)%z, hoek(I)%x, hoek(I)%y, hoek(I)%z
END DO
CLOSE(10)

OPEN(UNIT=20, FILE=dump_file, ACCESS="APPEND")
CALL DUMP(UNICORN+1)
CLOSE(20)

    WRITE (*,*) "Fase 1 done!"
    WRITE (*,*) "Fase 2 started!"

!====================================================================
!====================================================================

CALL system_clock(start)
!WRITE(500,*) START

    WRITE (*,*) "Entering Gaussian loop..."

    ! Loop 2: Gaussian
    loop_Ga: DO UNICORN=1+LJ_STEPS,GA_STEPS+LJ_STEPS

        ! Doe MC
        CALL MCINIT(UNICORN)

        ! Bereken veranderde interacties
        CALL calculateGA(RSOLV, UNICORN)
        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)

        OPEN(UNIT=20, FILE=dump_file, ACCESS="APPEND")
        CALL dump(UNICORN)
        CLOSE(20)

        ! Doe Metropolis
        CALL METROPOLIS(UNICORN, GA_NADJ, GA_NPRINT, REJECTED)

    END DO loop_Ga

    CALL system_clock(START)
    !write (500,*) START

!====================================================================
!====================================================================
! Wegschrijven resultaten

    WRITE (*,*) "Exiting Gaussian loop..."

    WRITE (*,*) "Fase 2 done!"
    WRITE (*,*) "Fase POST started!"

    WRITE (*,*) "Writing data..."

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=10, FILE=result_file)
    WRITE (10, *) BOXL ! Box grootte
    WRITE (10, *) nCoM ! Lees aantal moleculen
    DO I= 1, NCOM
        WRITE(10,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO
    CLOSE(10)

!====================================================================
!====================================================================

CLOSE(501)
CLOSE(500)

WRITE (*,*) "Fase POST done!"
WRITE (*,*) "We're done here. Signing off!"

CONTAINS

!====================================================================
!====================================================================

SUBROUTINE MCINIT(I)

        INTEGER :: I

        ! Verhuis oude vars naar de _old vars
        COM_OLD = COM
        HOEK_OLD = HOEK
        !solute_old = solute ! Wordt nog niet gevariëerd
        TOTENG_OLD = TOTENG
        EOLD = ENERGY
        SSOLD = SOLVENTSOLVENT

        ! Doe MC
        RSOLV = INT(RAND() * NCOM) + 1 ! Willekeurige DMSO molecule
        COM(RSOLV) = CoM(RSOLV) + randVec(DPOSMAX)
        HOEK(RSOLV) = hoek(RSOLV) + randVec(DHOEKMAX)


        ! Check of hoeken nog binnen [-Pi, +Pi] liggen
        HOEK(RSOLV)%X = HOEK(RSOLV)%X - TAU * AINT(HOEK(RSOLV)%X / PI)
        HOEK(RSOLV)%Y = HOEK(RSOLV)%Y - TAU * AINT(HOEK(RSOLV)%Y / PI)
        HOEK(RSOLV)%Z = HOEK(RSOLV)%Z - TAU * AINT(HOEK(RSOLV)%Z / PI)

        ! Periodic boundaries => Computer Simulation of Liquids, p30
        COM(RSOLV)%X = COM(RSOLV)%X - BOXL2 * AINT(COM(RSOLV)%X / BOXL)
        COM(RSOLV)%Y = COM(RSOLV)%Y - BOXL2 * AINT(COM(RSOLV)%Y / BOXL)
        COM(RSOLV)%Z = COM(RSOLV)%Z - BOXL2 * AINT(COM(RSOLV)%Z / BOXL)

        ! Checks for OOB
        IF(COM(RSOLV)%X .GT. BOXL .OR. COM(RSOLV)%X .LT. -1.D0 * BOXL) THEN
            WRITE(500,*) "OOB on X with ", RSOLV, "@", I, ":", COM(RSOLV)%X
        END IF
        IF(COM(RSOLV)%Y .GT. BOXL .OR. COM(RSOLV)%Y .LT. -1.D0 * BOXL) THEN
            WRITE(500,*) "OOB on Y with ", RSOLV, "@", I, ":", COM(RSOLV)%Y
        END IF
        IF(COM(RSOLV)%Z .GT. BOXL .OR. COM(RSOLV)%Z .LT. -1.D0 * BOXL) THEN
            WRITE(500,*) "OOB on Z with ", RSOLV, "@", I, ":", COM(RSOLV)%Z
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
        EXPONENT = -1.D0 * BETA * DELTA  * 1000.D0 / (8.314D0 * 300.D0)
        IF (EXPONENT .LT. -75.D0) THEN ! e^-75 < 3*10^-33: 0% kans anyway
        !write(500,*) "Large Exponent!", I
            KANS = 0.D0
            RV = 1.D0 ! Skip rand() voor cpu tijd besparing
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
            WRITE (501,902) I, TOTENG, TOTENG_OLD, KANS, RV, RSOLV, REAL(NSUC) / real(I), ratio, dposmax
        END IF

        ! Pas de dposMax aan indien ratio =/= 50%
        IF (mod(I, NADJ) .EQ. 0) THEN
            IF (ratio .GT. 0.5) THEN
                dposMax = dposMax * (1.D0 + pAdj)
            ELSE
                dposMax = dposMax * (1.D0 - pAdj)
            END IF
            IF (dposMax .LT. 0.00001) dposMax = 0.00001
            NACCEPT = 0
        END IF

        IF(REJECTED) THEN
            ! Verhuis oude vars terug naar de nieuwe
            COM = COM_OLD
            HOEK = HOEK_OLD
            !solute = solute_old ! Wordt nog niet gevariëerd
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

    DO J=I+1,NCOM
        ! Check lengte, met MIC!
        R = CoM(J) - CoM(I)
        R%X = R%X - BOXL2 * ANINT(R%X / BOXL)
        R%Y = R%Y - BOXL2 * ANINT(R%Y / BOXL)
        R%Z = R%Z - BOXL2 * ANINT(R%Z / BOXL)

        IF (length(R) .GT. 7) THEN
            EN = 0.D0
        ELSE
            MOL2 = RotMatrix(CoM(J), DMSO, hoek(J))
            CALL calcLJ(MOL1, MOL2, DMSO_SYM, TABLE_DMSO, EN, BOXL, BOXL2)
        END IF
        SOLVENTSOLVENT(I,J) = EN
        SOLVENTSOLVENT(J,I) = EN
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
    DOUBLE PRECISION :: R, RMIN
    TYPE (vector) :: TEMPJ
    LOGICAL, DIMENSION(NCOM+1) :: TooFar, Clipped
    LOGICAL :: MayContinue = .TRUE.

    ! Solvent solvent
    ! Notice: geen i loop: alleen molecule i is veranderd en moet opnieuw berekend worden
    SOLVENTSOLVENT(I,I) = 0.D0
    MOL1 = RotMatrix(CoM(I), DMSO, hoek(I))

    J = I + 1
    DO WHILE (J .LE. NCOM .AND. MayContinue)
        ! MINIMIZE DISTANCE
        RMIN = huge(en)
        K_LOOP: DO K=-1,1
            DO L=-1,1
                DO M=-1,1
                    TEMPJ%X = CoM(J)%X + FLOAT(K) * BOXL2
                    TEMPJ%Y = CoM(J)%Y + FLOAT(L) * BOXL2
                    TEMPJ%Z = CoM(J)%Z + FLOAT(M) * BOXL2

                    R = getDist(CoM(I), TEMPJ) ! Check
                    IF(R .LT. RMIN) THEN
                        RMIN = R
                        KMIN = K
                        LMIN = L
                        MMIN = M
                    ELSE
                    END IF
                END DO
            END DO
        END DO K_LOOP

        TEMPJ%X = CoM(J)%X + FLOAT(KMIN) * BOXL2
        TEMPJ%Y = CoM(J)%Y + FLOAT(LMIN) * BOXL2
        TEMPJ%Z = CoM(J)%Z + FLOAT(MMIN) * BOXL2

        !write(500,*) I, J, TEMPJ%X, TEMPJ%Y, TEMPJ%Z, CoM(J)%X, CoM(J)%Y, CoM(J)%Z, KMIN, LMIN, MMIN

        CLIPPED(j) = .FALSE.
        TOOFAR(j) = .FALSE.

        IF (RMIN .GT. 7.D0) THEN
            EN = 0.D0
            TooFar(j) = .TRUE.
        ELSE
            MOL2 = RotMatrix(TEMPJ, DMSO, hoek(J))
            CALL calcGa(I, J, MOL1, MOL2, DMSO_SYM, DMSO_SYM, CLIPPED(j), proc)
            IF(CLIPPED(j)) MayContinue = .FALSE.
        END IF

        J = J + 1
    END DO


    ! EXECUTE

    ! Begin of parallel loop
    !================================================================

        IF (MayContinue) THEN

        !$OMP PARALLEL
        !$OMP DO SCHEDULE(GUIDED) PRIVATE(En)
        EXEC: DO J=I+1,NCOM
            IF (CLipped(J)) THEN ! it may not run
                En = HUGE(En) ! Set energy to infinity
                WRITE(500, *) "Clipped, infty", I, J, LOOPNR
            ELSEIF (TooFar(j)) THEN ! too far apart
                En = 0.D0 ! Set energy to 0
            ELSE
                CALL execGa(I, J, EN)
                EN = EN - E_DMSO - E_DMSO
                EN = EN * HARTREE2KJMOL
            END IF

            SOLVENTSOLVENT(I,J) = EN
            SOLVENTSOLVENT(J,I) = EN
        END DO EXEC
        !$OMP END DO
        !$OMP END PARALLEL

        ! Solvent-solute
        CALL calcGa(I, 0, MOL1, SOLUTE, DMSO_SYM, SOL_SYM, Clipped(J+1), proc)
        IF (Clipped(J+1)) THEN
            EN = huge(en)
            WRITE (500,*) "Clipped with solute, infty", I, 0, LOOPNR
        ELSE
            CALL execGa(I, 0, EN)
            EN = EN - E_SOL - E_DMSO
            EN = EN * HARTREE2KJMOL
        END IF

        ENERGY(I) = EN
    ELSE
        ENERGY(I) = HUGE(EN)
        SOLVENTSOLVENT(I,J) = HUGE(EN)
        SOLVENTSOLVENT(J,I) = HUGE(EN)
    END IF

END SUBROUTINE calculateGA

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

SUBROUTINE dump(i)
        INTEGER :: i

        !WRITE (20,*) NCOM*NDMSO+NSOL - 6*NCOM
        WRITE (20,*) NCOM+NSOL
        WRITE (20,*) "Timestep: ", I
        DO J=1, NSOL
            WRITE(20,*) sol_sym(J), solute(J)%x, solute(J)%y, solute(J)%z
        END DO
        DO J=1,NCOM
        !WRITE (20,*) "S", CoM(J)%x, CoM(J)%y, CoM(J)%z
        IF(.TRUE.) THEN
            ABSPOS = RotMatrix(CoM(J), DMSO, hoek(J))
            DO K=1, NDMSO
                IF (DMSO_sym(K) .NE. "H") THEN
                WRITE(20,*) DMSO_Sym(K), absPos(K)%x, absPos(K)%y, absPos(K)%z
                END IF
            END DO
            END IF
        END DO

END SUBROUTINE dump

!====================================================================
!====================================================================

! The end :)
END PROGRAM MonteCarlo
