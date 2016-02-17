!====================================================================
! MonteCarlo.f95
! Auteur: Matthijs Lasure
!
! Description:
! Doe een Monte-Carlo simulatie met 1 solute in DMSO
! Eerst werken met LJ voor ruwe benadering (eg 10.000 stappen)
! Dan verfijnen met Gaussian (eg 250 stappen)
!====================================================================

program MonteCarlo


    !use lib
    !use interactions
    !implicit none

    ! Variabelen
    double precision, parameter :: PI = 4.D0 * DATAN(1.D0)
    integer :: i, j
    integer :: LJ_steps, Ga_steps ! Aantal stappen per loop

    ! Config
    ! TODO inlezen uit config.txt?
    LJ_steps = 10000 ! Ruwe loop
    Ga_steps = 1 ! Loop met Gaussian

    ! Laden van configuraties
    ! DMSO.txt: conformatie DMSO
    ! solute.txt: conformatie opgeloste molecule
    ! box.txt: plaatsen van de moleculen (CoM)

    ! Initiële berekening interacties

    write (*,*) PI
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
