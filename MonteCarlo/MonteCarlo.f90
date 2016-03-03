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

! Calls in program MONTECARLO: 
! => system_clock (on line <82>)
! ==> SRAND (on line <83>)
! ===> calculateLJ (on line <156>)
! ====> calculateLJ (on line <251>)
PROGRAM MonteCarlo

    USE vector_class
    USE interactions
    USE lib
    USE randgen
    USE readconfig
    !use iso_fortran_env

    IMPLICIT NONE


    ! Variabelen
    !===========

    DOUBLE PRECISION:: BOXL, BOXL2 ! box grootte, halve box
    INTEGER:: I, J, UNICORN, ISEED
    INTEGER:: LJ_STEPS, GA_STEPS ! Aantal stappen per loop
    INTEGER:: RSOLV ! Geselecteerde molecule voor MC
    INTEGER:: NSUC = 0 ! Aantal succesvolle MC
    INTEGER :: NADJ, NPRINT, NACCEPT = 0 ! Aanpassen dposMax / printen om de n cycli
    DOUBLE PRECISION:: RV, KANS,DELTA = 0, EXPONENT ! Random variabele en toebehoren voor Metropolis
    DOUBLE PRECISiON :: RATIO ! Percentage succes in NADJ trials
    DOUBLE PRECISION :: PADJ ! Hoeveel de dposmax aangepast mag worden
    DOUBLE PRECISION :: BETA ! p = EXP(-BETA * DELTA / (RT)

    ! FILES
    !======
    character*100 :: dmso_file, box_file, sol_file, param_file ! input files
    CHARACTER*100 :: out_file, err_file, solvsolv_file ! output files

    ! Energiën van de moleculen, bekomen via extern programma
    DOUBLE PRECISION:: E_DMSO, E_SOL

    ! Vectoren voor moleculen & atoomtypes
    TYPE (vector), DIMENSION(:), ALLOCATABLE :: DMSO, COM, SOLUTE, HOEK
    ! Variabelen voor de vorige run van MC
    TYPE (vector), DIMENSION(:), ALLOCATABLE :: COM_OLD, SOLUTE_OLD, HOEK_OLD
    CHARACTER*4, DIMENSION(:), ALLOCATABLE :: DMSO_SYM, SOL_SYM
    INTEGER:: NDMSO, NCOM, NSOL, NPARAM ! Aantal units
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
    DOUBLE PRECISION:: EN ! Energie van een run
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SOLVENTSOLVENT, SSOLD
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ENERGY, EOLD
    DOUBLE PRECISION:: TOTENG = 0.D0
    DOUBLE PRECISION:: TOTENG_OLD = 0.D0

    ! Arrays voor parameters van DMSO (Q, epsilon, sigma, mass)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Q, EPSILON, SIGMA, MASS
    CHARACTER*4, DIMENSION(:), ALLOCATABLE :: SYM

    ! Maximale verarndering bij MC
    DOUBLE PRECISION:: DPOSMAX, DHOEKMAX
    DPOSMAX = 0.35D0
    DHOEKMAX = 1 ! In aantal * Pi

    ! START ERR/OUT
    out_file = "out.txt"
    err_file = "err.txt"
    OPEN(UNIT=501, FILE=out_file)
    OPEN(UNIT=500, FILE=err_file)

    ! Config
!====================================================================
!====================================================================
    ! Read from config.ini
    CALL rConfig(BOXL, LJ_STEPS, Ga_STEPS, iseed, DoDebug, nadj, nprint, dposmax, dhoekmax, padj, beta)

    CALL system_clock (START)
    write(500,*) START
    CALL SRAND(REAL(iseed)) ! Prime randgen

    BOXL2 = BOXL * 2.D0
    DHOEKMAX = DHOEKMAX * PI

    ! INPUT
    dmso_file = "DMSO.txt"
    box_file = "box.txt"
    sol_file = "solute.txt"
    param_file = "param.txt"

    ! OUTPUT
    solvsolv_file = "solventsolvent.txt"


    ! Laden van configuraties
!====================================================================
!====================================================================

    ! DMSO.txt: conformatie DMSO
    OPEN (UNIT=10, FILE=dmso_file)
    READ (10, *) nDMSO ! Lees aantal atomen
    ALLOCATE(DMSO(NDMSO)) ! Ken correcte groottes toe aan de arrays
    ALLOCATE(DMSO_sym(NDMSO))
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
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (10,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    READ (10,*) E_sol
    CLOSE(10)

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=10, FILE=box_file)
    READ (10, *) !BOXL ! Box grootte
    READ (10, *) nCoM ! Lees aantal moleculen
    ALLOCATE(CoM(NCOM))
    ALLOCATE(hoek(NCOM))
    ALLOCATE(CoM_old(NCOM))
    ALLOCATE(hoek_old(NCOM))
    DO I= 1, NCOM
        READ(10,*) CoM(I)%X, CoM(I)%Y, CoM(I)%Z, hoek(I)%X, hoek(I)%Y, hoek(I)%Z
    END DO
    CLOSE(10)

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

!====================================================================
!====================================================================

    ! Initiële berekening interacties
    !================================

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
    OPEN(UNIT=10, FILE="solventsolvent.txt")
    DO I=1,NCOM
        WRITE(10,*) solventsolvent(I,:)
    END DO
    CLOSE(10)


!====================================================================
!====================================================================



OPEN(UNIT=20, FILE="DUMP.txt")
CALL DUMP(0)
CLOSE(20)

901 FORMAT(A12, 1X, A20, 1X, A20, 1X, A6, 1X, A6, 1X, A3, 1X, A6, 1X, A6, 1X, A6)
902 FORMAT(I12.12, 1X, ES20.10, 1X, ES20.10, 1X, F6.4, 1X, F6.4, 1X, I3.3, 1X, F6.4, 1X, F6.4, 1X, F6.5)
WRITE (501,901) "i", "TotEng", "TotEng_old", "kans", "rv", "rSolv", "pSuc", "ratio","dposmax"
WRITE (501,902) 0, TOTENG, 0.D0, 0.D0, 0.D0, 0, REAL(0) / real(1), 0.D0, dposmax

CALL system_clock(start)
write(500, *) start


    ! Loop 1: LJ
    !===========
    loop_LJ: DO UNICORN=1,LJ_STEPS

        CALL MCINIT(UNICORN) ! Verander 1 molecule + check voor OOB & PB

        ! Bereken veranderde interacties
        CALL calculateLJ(RSOLV)
        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)

        ! Check for invalid shit
        if(TOTENG > HUGE(TOTENG)) then
            write (0,*) "TOTENG IS INFINITY @", UNICORN, TOTENG
            TOTENG = HUGE(TOTENG)-1
        end if
        if(TOTENG .NE. TOTENG) then
            write (0,*) "TOTENG IS NaN @", UNICORN, TOTENG
            TOTENG = HUGE(TOTENG)-1
        end if

        ! Doe Metropolis
        CALL METROPOLIS(UNICORN)

    END DO loop_LJ

! Write box
open(unit=10, file="backup.txt")
write(10,*) boxl
write(10,*) nCoM
do I=1,nCoM
    write(10,*) CoM(I)%x, CoM(I)%y, CoM(I)%z, hoek(I)%x, hoek(I)%y, hoek(I)%z
end do
close(10)

OPEN(UNIT=20, FILE="DUMP.txt", ACCESS="APPEND")
CALL DUMP(UNICORN+1)
CLOSE(20)

!====================================================================
!====================================================================
CALL system_clock(start)
WRITE(500,*) START

    ! Loop 2: Gaussian
    loop_Ga: DO UNICORN=1,GA_STEPS
        ! Doe MC
        CALL MCINIT(UNICORN)

        ! Bereken veranderde interacties
        CALL calculateGA(RSOLV, UNICORN)
        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)

        call dump(UNICORN)

        ! Doe Metropolis
        CALL METROPOLIS(UNICORN)

    END DO loop_Ga

    CALL system_clock(START)
    write (500,*) START

!====================================================================
!====================================================================

!close(20)
close(501)
close(500)

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
        if(COM(RSOLV)%X .GT. BOXL .OR. COM(RSOLV)%X .LT. -1.D0 * BOXL) THEN
            write(500,*) "OOB on X with ", RSOLV, "@", I, ":", COM(RSOLV)%X
        END IF
        if(COM(RSOLV)%Y .GT. BOXL .OR. COM(RSOLV)%Y .LT. -1.D0 * BOXL) THEN
            write(500,*) "OOB on Y with ", RSOLV, "@", I, ":", COM(RSOLV)%Y
        END IF
        if(COM(RSOLV)%Z .GT. BOXL .OR. COM(RSOLV)%Z .LT. -1.D0 * BOXL) THEN
            write(500,*) "OOB on Z with ", RSOLV, "@", I, ":", COM(RSOLV)%Z
        END IF

END SUBROUTINE MCINIT

!====================================================================
!====================================================================

SUBROUTINE METROPOLIS(I)

        INTEGER :: I

        902 FORMAT(I12.12, 1X, ES20.10, 1X, ES20.10, 1X, F6.4, 1X, F6.4, 1X, I3.3, 1X, F6.4, 1X, F6.4, 1X, F6.5)

        DELTA = TOTENG - TOTENG_OLD
        EXPONENT = -1.D0 * BETA * DELTA  * 1000.D0 / (8.314D0 * 300.D0)
        if (EXPONENT .LT. -75.D0) then ! e^-75 < 3*10^-33: 0% kans anyway
        write(500,*) "Large Exponent!", I
            KANS = 0.D0
            RV = 1.D0 ! Skip rand() voor cpu tijd besparing
        else
            KANS = E ** EXPONENT
            RV = RAND()
            IF (KANS .GT. 1.D0) KANS = 1.D0
        end if

        ! Bepaal if succesvol -> volgende config
        IF(RV .LE. KANS) THEN ! Succes!
            NSUC = NSUC + 1
            NACCEPT = NACCEPT + 1
        ELSE ! Fail!
            ! Verhuis oude vars terug naar de nieuwe
            COM = COM_OLD
            HOEK = HOEK_OLD
            !solute = solute_old ! Wordt nog niet gevariëerd
            TOTENG = TOTENG_OLD
            ENERGY = EOLD ! Resetten!!
            SOLVENTSOLVENT = SSOLD
        END IF

        RATIO = REAL(NACCEPT) / real(NADJ)

        ! Prints every NPRINT times
        if(mod(I, NPRINT) .EQ. 0) then
            !CALL DUMP(I)
            WRITE (501,902) I, TOTENG, TOTENG_OLD, KANS, RV, RSOLV, REAL(NSUC) / real(I), ratio, dposmax
        end if

        ! Pas de dposMax aan indien ratio =/= 50%
        if (mod(I, NADJ) .EQ. 0) then
            if (ratio .GT. 0.5) then
                dposMax = dposMax * (1.D0 + pAdj)
            else
                dposMax = dposMax * (1.D0 - pAdj)
            end if
            if (dposMax .LT. 0.00001) dposMax = 0.00001
            NACCEPT = 0
        end if

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
            CALL calcLJ(MOL1, MOL2, DMSO_SYM, DMSO_SYM, SYM, Q, EPSILON, SIGMA, EN, BOXL, BOXL2)
        END IF
        SOLVENTSOLVENT(I,J) = EN
        SOLVENTSOLVENT(J,I) = EN
    END DO

    ! Solvent-solute
    CALL calcLJ(MOL1, SOLUTE, DMSO_SYM, SOL_SYM, SYM, Q, EPSILON, SIGMA, EN, BOXL, BOXL2)
    ENERGY(I) = EN

END SUBROUTINE calculateLJ

!====================================================================
!====================================================================

SUBROUTINE calculateGA(I, LOOPNR)

    INTEGER:: I, J, LOOPNR ! gevraagde moleculen, welke loop
    INTEGER :: K, L, M, KMIN = 2, LMIN = 2, MMIN = 2 ! Minimal Image Convention
    DOUBLE PRECISION :: R, RMIN
    TYPE (vector) :: TEMPJ

    ! Solvent solvent
    ! Notice: geen i loop: alleen molecule i is veranderd en moet opnieuw berekend worden
    SOLVENTSOLVENT(I,I) = 0.D0
    MOL1 = RotMatrix(CoM(I), DMSO, hoek(I))


    ! Begin of parallel loop
    !================================================================

    !$OMP PARALLEL
    !$OMP DO SCHEDULE(GUIDED)
    DO J=I+1,NCOM
        ! MINIMIZE DISTANCE
        RMIN = huge(en)
        K_LOOP: DO K=-1,1
            DO L=-1,1
                DO M=-1,1
                    TEMPJ%X = CoM(J)%X + FLOAT(K) * BOXL2
                    TEMPJ%Y = CoM(J)%Y + FLOAT(L) * BOXL2
                    TEMPJ%Z = CoM(J)%Z + FLOAT(M) * BOXL2

                    R = getDist(CoM(I), TEMPJ) ! Check
                    if(R .LT. RMIN) THEN
                        RMIN = R
                        KMIN = K
                        LMIN = L
                        MMIN = M
                    else
                    end if
                END DO
            END DO
        END DO K_LOOP

        TEMPJ%X = CoM(J)%X + FLOAT(KMIN) * BOXL2
        TEMPJ%Y = CoM(J)%Y + FLOAT(LMIN) * BOXL2
        TEMPJ%Z = CoM(J)%Z + FLOAT(MMIN) * BOXL2

        !write(500,*) I, J, TEMPJ%X, TEMPJ%Y, TEMPJ%Z, CoM(J)%X, CoM(J)%Y, CoM(J)%Z, KMIN, LMIN, MMIN


        IF (RMIN .GT. 7.D0) THEN
            EN = 0.D0
        ELSE
            MOL2 = RotMatrix(TEMPJ, DMSO, hoek(J))
            !CALL system("rm gauss/*")
            CALL calcGa(I, J, MOL1, MOL2, DMSO_SYM, DMSO_SYM, EN, LOOPNR)

            !WRITE(500, *) "Ga - ", I, "-", J

            EN = EN - E_DMSO - E_DMSO
            EN = EN * HARTREE2KJMOL
        END IF

        SOLVENTSOLVENT(I,J) = EN
        SOLVENTSOLVENT(J,I) = EN
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    ! Solvent-solute
    CALL calcGa(I, 0, MOL1, SOLUTE, DMSO_SYM, SOL_SYM, EN, LOOPNR)
    EN = EN - E_SOL - E_DMSO
    EN = EN * HARTREE2KJMOL
    ENERGY(I) = EN

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

subroutine dump(i)
        INTEGER :: i

        !WRITE (20,*) NCOM*NDMSO+NSOL - 6*NCOM
        WRITE (20,*) NCOM+NSOL
        WRITE (20,*) "Timestep: ", I
        DO J=1, NSOL
            WRITE(20,*) sol_sym(J), solute(J)%x, solute(J)%y, solute(J)%z
        END DO
        DO J=1,NCOM
        WRITE (20,*) "S", CoM(J)%x, CoM(J)%y, CoM(J)%z
        IF(.FALSE.) THEN
            ABSPOS = RotMatrix(CoM(J), DMSO, hoek(J))
            DO K=1, NDMSO
                IF (DMSO_sym(K) .NE. "H") THEN
                WRITE(20,*) DMSO_Sym(K), absPos(K)%x, absPos(K)%y, absPos(K)%z
                END IF
            END DO
            END IF
        END DO

end subroutine dump

!====================================================================
!====================================================================

! The end :)
END PROGRAM MonteCarlo
