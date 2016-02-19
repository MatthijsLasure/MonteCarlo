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
    use lib

    implicit none

    ! Pi -> lib
    !double precision, parameter :: PI = 4.D0 * DATAN(1.D0)

    ! Variabelen
    !===========

    double precision :: box_length ! box grootte
    integer :: i, j
    integer :: LJ_steps, Ga_steps ! Aantal stappen per loop

    ! Energi�n van de moleculen, bekomen via extern programma
    double precision :: E_DMSO, E_sol

    ! Vectoren voor moleculen & atoomtypes
    TYPE (vector), dimension(:), allocatable :: DMSO, CoM, solute, hoek
    character*4, dimension(:), allocatable :: DMSO_sym, sol_sym
    integer :: nDMSO, nCoM, nSol, nParam ! Aantal units
    ! DMSO: relatieve co�rdinaten voor de atomen
    ! CoM: Centre of Mass: locaties van de DMSO moleculen
    ! solute: conformatie van de solute
    ! Hoek: ori�ntatie van de DMSO moleculen -> RotMatrix

    ! Tijdelijke vectoren voor mol
    TYPE (vector), dimension(:), allocatable :: mol1, mol2

    ! Debug
    !integer :: k

    ! output calc
    double precision :: en
    double precision, dimension(:,:), allocatable :: solventsolvent
    double precision, dimension(:), allocatable :: energy
    double precision :: TotEng

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
    do i=1, nSol ! Lees de co�rdinaten uit
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

    ! Initi�le berekening interacties
    !================================
    ! hoe kan dit in een functie?
    ! Vars naar lib???

    ! Arrays
    allocate(mol1(nDMSO))
    allocate(mol2(nDMSO))
    allocate(solventsolvent(nCoM, nCoM))
    allocate(energy(nCoM))

    call calculate

    open(unit=10, file="solventsolvent.txt")
    do i=1,nCoM
        write(10,*) solventsolvent(i,:)
    end do
    close(10)


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

contains

subroutine calculate

    ! Solvent-solvent
    do i=1,nCoM
        solventsolvent(i,i) = 0.D0
        mol1 = RotMatrix(CoM(i), DMSO, hoek(i))
        do j=i+1,nCoM
            mol2 = RotMatrix(CoM(j), DMSO, hoek(j))
            write (*,*) "calcLJ on ", i, " + ", j
            call calcLJ(mol1, mol2, DMSO_sym, DMSO_sym, sym, Q, epsilon, sigma, en)

            solventsolvent(i,j) = en
            solventsolvent(j,i) = en
        end do
    end do

    write (*,*) "Solvent-Solute"

    ! Solvent-solute
    do i=1,nCoM
        mol1 = RotMatrix(CoM(i), DMSO, hoek(i))
        write (*,*) "calcLJ on ", i, " + solute"
        call calcLJ(mol1, solute, DMSO_sym, sol_sym, sym, Q, epsilon, sigma, en)
        energy(i) = en
    end do

    ! Totale E
    totEng = 0.D0

    do i=1, nCoM
        totEng = totEng + energy(i) ! solv - solu
        do j=i+1, nCoM
            totEng = totEng + solventsolvent(i,j)
        end do
    end do
    write (*,*) totEng

end subroutine calculate
end program MonteCarlo
