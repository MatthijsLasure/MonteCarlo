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

program MonteCarlo

    use vector_class
    use interactions

    implicit none

    ! Pi
    double precision, parameter :: PI = 4.D0 * DATAN(1.D0)

    ! Variabelen
    !===========

    double precision :: box_length ! box grootte
    integer :: i, j
    integer :: LJ_steps, Ga_steps ! Aantal stappen per loop

    ! Energiën van de moleculen, bekomen via extern programma
    double precision :: E_DMSO, E_sol

    ! Vectoren voor moleculen & atoomtypes
    TYPE(vector), dimension(:), allocatable :: DMSO, CoM, solute, hoek
    character*4, dimension(:), allocatable :: DMSO_sym, sol_sym
    integer :: nDMSO, nCoM, nSol, nParam ! Aantal units
    ! DMSO: relatieve coördinaten voor de atomen
    ! CoM: Centre of Mass: locaties van de DMSO moleculen
    ! solute: conformatie van de solute
    ! Hoek: oriëntatie van de DMSO moleculen -> RotMatrix

    ! Arrays voor parameters van DMSO (Q, epsilon, sigma, mass)
    double precision, dimension(:), allocatable :: Q, epsilon, sigma, mass
    character*4, dimension(:), allocatable :: sym

    ! Config
    !=======

    ! TODO inlezen uit config.txt?
    LJ_steps = 10000 ! Ruwe loop
    Ga_steps = 1 ! Loop met Gaussian

    ! Laden van configuraties
    !========================

    ! DMSO.txt: conformatie DMSO
    open (unit=10, file="DMSO.txt")
    read (10, *) nDMSO ! Lees aantal atomen
    allocate(DMSO(nDMSO)) ! Ken correcte groottes toe aan de arrays
    allocate(DMSO_sym(nDMSO))
    do i=1, nDMSO
        read (10,*) DMSO_sym(i), DMSO(i)%x, DMSO(i)%y, DMSO(i)%z
    end do
    read (10,*) E_DMSO
    close(10)

    ! solute.txt: conformatie opgeloste molecule (sol)
    open (unit=10, file="solute.txt")
    read (10, *) nSol ! Lees aantal atomen
    allocate(solute(nSol)) ! Maak de arrays groot genoeg
    allocate(sol_sym(nSol))
    do i=1, nSol ! Lees de coördinaten uit
        read (10,*) sol_sym(i), solute(i)%x, solute(i)%y, solute(i)%z
    end do
    read (10,*) E_sol
    close(10)

    ! box.txt: plaatsen van de moleculen (CoM, hoek)
    open (unit=10, file="box.txt")
    read (10, *) box_length ! Box grootte
    read (10, *) nCoM ! Lees aantal moleculen
    allocate(CoM(nCoM))
    allocate(hoek(nCoM))
    do i= 1, nCoM
        read(10,*) CoM(i)%x, CoM(i)%y, CoM(i)%z, hoek(i)%x, hoek(i)%y, hoek(i)%z
    end do
    close(10)

    ! param.txt: parameters voor LJ etc
    open (unit=10, file="param.txt")
    read (10, *) nParam ! Aantal beschikbare parameters
    read (10, *) ! skip comment line
    allocate(sym(nParam))
    allocate(Q(nParam))
    allocate(epsilon(nParam))
    allocate(sigma(nParam))
    allocate(mass(nParam))
    do i= 1,nParam
        read (10,*) sym(i), Q(i), epsilon(i), sigma(i), mass(i)
    end do
    close(10)

    ! Initiële berekening interacties

    ! Loop 1: LJ
    loop_LJ: do i=1,LJ_steps
        ! Doe MC

        ! Bereken veranderde interacties

        ! Doe Metropolis

        ! Bepaal if succesvol -> volgende config
    end do loop_LJ

    ! Loop 2: Gaussian
    loop_Ga: do i=1,Ga_steps
        ! Doe MC

        ! Bereken veranderde interacties

        ! Doe Metropolis

        ! Bepaal if succesvol -> volgende config
    end do loop_Ga


end program MonteCarlo
