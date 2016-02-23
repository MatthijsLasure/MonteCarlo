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

    IMPLICIT NONE


    ! Variabelen
    !===========

    DOUBLE PRECISION:: BOXL, BOXL2 ! box grootte, halve box
    INTEGER:: I, J, UNICORN, ISEED
    INTEGER:: LJ_STEPS, GA_STEPS ! Aantal stappen per loop
    INTEGER:: RSOLV ! Geselecteerde molecule voor MC
    INTEGER:: NSUC = 0 ! Aantal succesvolle MC
    DOUBLE PRECISION:: RV, KANS,DELTA = 0, EXPONENT ! Random variabele en toebehoren voor Metropolis
    CHARACTER(LEN=32), DIMENSION(5) :: ARG

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
    DOUBLE PRECISION:: EN
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SOLVENTSOLVENT, SSOLD
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ENERGY, EOLD
    DOUBLE PRECISION:: TOTENG = 0.D0
    DOUBLE PRECISION:: TOTENG_OLD = 0.D0
    DOUBLE PRECISION :: RADIUS

    ! Arrays voor parameters van DMSO (Q, epsilon, sigma, mass)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Q, EPSILON, SIGMA, MASS
    CHARACTER*4, DIMENSION(:), ALLOCATABLE :: SYM

    ! Maximale verarndering bij MC
    DOUBLE PRECISION:: DPOSMAX, DHOEKMAX
    DPOSMAX = 0.25D0
    DHOEKMAX = PI

    ! Config
!====================================================================
!====================================================================
    CALL rConfig(BOXL, LJ_STEPS, Ga_STEPS, iseed, DoDebug)

    CALL system_clock (START)
    write(*,*) START
    CALL SRAND(REAL(iseed)) ! Prime randgen

    BOXL2 = BOXL * 2

    ! TODO inlezen uit config.txt?


    !LJ_STEPS = 1000 ! Ruwe loop
    GA_STEPS = 1 ! Loop met Gaussian
    !BOXL = 15.7D0

    ! Laden van configuraties
!====================================================================
!====================================================================

    ! DMSO.txt: conformatie DMSO
    OPEN (UNIT=10, FILE="DMSO.txt")
    READ (10, *) nDMSO ! Lees aantal atomen
    ALLOCATE(DMSO(NDMSO)) ! Ken correcte groottes toe aan de arrays
    ALLOCATE(DMSO_sym(NDMSO))
    DO I=1, NDMSO
        READ (10,*) DMSO_sym(I), DMSO(I)%X, DMSO(I)%Y, DMSO(I)%Z
    END DO
    READ (10,*) E_DMSO
    CLOSE(10)

    ! solute.txt: conformatie opgeloste molecule (sol)
    OPEN (UNIT=10, FILE="solute.txt")
    READ (10, *) nSol ! Lees aantal atomen
    ALLOCATE(solute(NSOL)) ! Maak de arrays groot genoeg
    ALLOCATE(sol_sym(NSOL))
    DO I=1, NSOL ! Lees de coördinaten uit
        READ (10,*) sol_sym(I), solute(I)%X, solute(I)%Y, solute(I)%Z
    END DO
    READ (10,*) E_sol
    CLOSE(10)

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    OPEN (UNIT=10, FILE="box.txt")
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
    OPEN (UNIT=10, FILE="param.txt")
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

IF (DODEBUG) THEN
OPEN(UNIT=11, FILE="out/DUMP.txt")
! DUMP
!WRITE (11,*) NCOM*NDMSO+NSOL - 6*NCOM
WRITE (11,*) NCOM+NSOL
WRITE (11,*) "Timestep: ", 0
DO J=1, NSOL
    WRITE(11,*) sol_sym(J), solute(J)%x, solute(J)%y, solute(J)%z
END DO
DO J=1,NCOM
    WRITE (11,*) "S", CoM(J)%x, CoM(J)%y, CoM(J)%z
    IF(.FALSE.) THEN
        ABSPOS = RotMatrix(CoM(J), DMSO, hoek(J))
        DO K=1, NDMSO
            IF (DMSO_sym(K) .NE. "H") THEN
            WRITE(11,*) DMSO_Sym(K), absPos(K)%x, absPos(K)%y, absPos(K)%z
            END IF
        END DO
    END IF
END DO
END IF

!====================================================================
!====================================================================

WRITE (*,"(A, 1X, A, 1X, A, 1X, A, 1X, A, 1X, A, 1X, A)") "i", "TotEng", "TotEng_old", "kans", "rv", "rSolv", "pSuc"


    ! Loop 1: LJ
    !===========
    loop_LJ: DO UNICORN=1,LJ_STEPS
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
        HOEK(RSOLV)%X = HOEK(RSOLV)%X - TAU * ANINT(HOEK(RSOLV)%X / PI)
        HOEK(RSOLV)%Y = HOEK(RSOLV)%Y - TAU * ANINT(HOEK(RSOLV)%Y / PI)
        HOEK(RSOLV)%Z = HOEK(RSOLV)%Z - TAU * ANINT(HOEK(RSOLV)%Z / PI)

        IF(.TRUE.)THEN
        ! Periodic boundaries => Computer Simulation of Liquids, p30
        COM(RSOLV)%X = COM(RSOLV)%X - BOXL2 * ANINT(COM(RSOLV)%X / BOXL)
        COM(RSOLV)%Y = COM(RSOLV)%Y - BOXL2 * ANINT(COM(RSOLV)%Y / BOXL)
        COM(RSOLV)%Z = COM(RSOLV)%Z - BOXL2 * ANINT(COM(RSOLV)%Z / BOXL)
        END IF

        ! DUMP
        IF(DODEBUG) THEN
        !WRITE (11,*) NCOM*NDMSO+NSOL - 6*NCOM
        WRITE (11,*) NCOM+NSOL
        WRITE (11,*) "Timestep: ", UNICORN
        DO J=1, NSOL
            WRITE(11,*) sol_sym(J), solute(J)%x, solute(J)%y, solute(J)%z
        END DO
        DO J=1,NCOM
        WRITE (11,*) "S", CoM(J)%x, CoM(J)%y, CoM(J)%z
        IF(.FALSE.) THEN
            ABSPOS = RotMatrix(CoM(J), DMSO, hoek(J))
            DO K=1, NDMSO
                IF (DMSO_sym(K) .NE. "H") THEN
                WRITE(11,*) DMSO_Sym(K), absPos(K)%x, absPos(K)%y, absPos(K)%z
                END IF
            END DO
            END IF
        END DO
        END IF

        ! Bereken veranderde interacties
        CALL calculateLJ(RSOLV)
        TOTENG = 0.D0
        TOTENG = calcEnergy(ENERGY, SOLVENTSOLVENT)

        ! Doe Metropolis
        DELTA = TOTENG - TOTENG_OLD
        EXPONENT = -1.D0 * DELTA  * 1000.D0 / (8.315D0 * 300.D0)
        KANS = E ** EXPONENT
        IF (KANS .GT. 1.D0) KANS = 1.D0
        RV = RAND()

        WRITE (*,"(I12.12, 1X, ES20.10, 1X, ES20.10, 1X, F6.4, 1X, F6.4, 1X, I3.3, 1X, F6.4)") &
        UNICORN, TOTENG, TOTENG_OLD, KANS, RV, RSOLV, real(NSUC) / real(UNICORN)

        ! Bepaal if succesvol -> volgende config
        IF(RV .LE. KANS) THEN ! Succes!
            NSUC = NSUC + 1
        ELSE ! Fail!
            ! Verhuis oude vars terug naar de nieuwe
            COM = COM_OLD
            HOEK = HOEK_OLD
            !solute = solute_old ! Wordt nog niet gevariëerd
            TOTENG = TOTENG_OLD
            ENERGY = EOLD ! Resetten!!
            SOLVENTSOLVENT = SSOLD
        END IF

    END DO loop_LJ

    IF(DODEBUG) CLOSE(11)

    ! DUMP
OPEN(UNIT=11, FILE="DUMP.txt")
WRITE (11,*) NCOM*NDMSO+NSOL - 6*NCOM
!WRITE (11,*) NCOM+NSOL
WRITE (11,*) "Timestep: ", UNICORN
DO J=1, NSOL
    WRITE(11,*) sol_sym(J), solute(J)%x, solute(J)%y, solute(J)%z
END DO
DO J=1,NCOM
    !WRITE (11,*) "S", CoM(J)%x, CoM(J)%y, CoM(J)%z
    ABSPOS = RotMatrix(CoM(J), DMSO, hoek(J))
    DO K=1, NDMSO
        IF (DMSO_sym(K) .NE. "H") THEN
        WRITE(11,*) DMSO_Sym(K), absPos(K)%x, absPos(K)%y, absPos(K)%z
        END IF
    END DO
END DO
close(11)

! Write box
open(unit=10, file="backup.txt")
write(10,*) boxl
write(10,*) nCoM
do I=1,nCoM
    write(10,*) CoM(I)%x, CoM(I)%y, CoM(I)%z, hoek(I)%x, hoek(I)%y, hoek(I)%z
end do
close(10)

!====================================================================
!====================================================================

    ! Loop 2: Gaussian
    loop_Ga: DO I=1,GA_STEPS
        ! Doe MC

        ! Bereken veranderde interacties
        J = I
        ! Doe Metropolis

        ! Bepaal if succesvol -> volgende config
    END DO loop_Ga

    CALL system_clock(START)
    write (*,*) START

!====================================================================
!====================================================================
CONTAINS

!====================================================================
!====================================================================

! Calls in subroutine CALCULATELJ: 
! => calcLJ (on line <316>)
! ==> calcLJ (on line <324>)
SUBROUTINE calculateLJ(I)

    INTEGER:: I
    INTEGER:: J

    ! Solvent solvent
    ! Notice: geen i loop: alleen molecule i is veranderd en moet opnieuw berekend worden
    SOLVENTSOLVENT(I,I) = 0.D0
    MOL1 = RotMatrix(CoM(I), DMSO, hoek(I))
    DO J=I+1,NCOM
        MOL2 = RotMatrix(CoM(J), DMSO, hoek(J))
        CALL calcLJ(MOL1, MOL2, DMSO_SYM, DMSO_SYM, SYM, Q, EPSILON, SIGMA, EN, BOXL, BOXL2)

        SOLVENTSOLVENT(I,J) = EN
        SOLVENTSOLVENT(J,I) = EN
    END DO

    ! Solvent-solute
    !do i=1,nCoM ! Not needed!
        CALL calcLJ(MOL1, SOLUTE, DMSO_SYM, SOL_SYM, SYM, Q, EPSILON, SIGMA, EN, BOXL, BOXL2)
        ENERGY(I) = EN
    !end do


    !!write (*,*) totEng

END SUBROUTINE calculateLJ

FUNCTION calcEnergy(ENERGY, SOLVENTSOLVENT) RESULT(out)
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ENERGY
    DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: SOLVENTSOLVENT
    DOUBLE PRECISION:: OUT
    INTEGER:: I
    INTEGER:: J

    ! Totale E
    OUT = 0.D0

    DO I=1, NCOM
        OUT = OUT + energy(I) ! solv - solu
        DO J=I+1, NCOM
            OUT = OUT + solventsolvent(I,J)
        END DO
    END DO

END FUNCTION calcEnergy

END PROGRAM MonteCarlo
