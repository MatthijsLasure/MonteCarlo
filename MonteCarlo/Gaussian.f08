!====================================================================
! interactions.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Bereken energie via LJ of QC
!====================================================================
MODULE gaussian

    USE vector_class
    USE lib

    CONTAINS

!====================================================================

SUBROUTINE calcGa(I, J, MOL1, MOL2, SYM1, SYM2, HASCLIPPED, PROC)
    ! INPUT
    TYPE (vector), DIMENSION(:), INTENT(IN) :: MOL1, MOL2 ! absolute coords!
    CHARACTER*4, DIMENSION(:), INTENT(IN) :: SYM1, SYM2 ! Atoomtypes
    INTEGER, INTENT(IN) :: I, J
    INTEGER, INTENT(IN) :: PROC ! aantal processoren voor gaussian

    ! OUTPUT
    LOGICAL, INTENT(OUT) :: HASCLIPPED ! Says if this set may run

    ! INTERNAL VARS
    INTEGER :: K, L ! loop de loop
    CHARACTER*100 :: GAUSS_FILE
    CHARACTER*16 :: STR_I, STR_J
    INTEGER :: N1, N2
    INTEGER :: IOSTATUS ! Check for EOF


    905 FORMAT(A, 3F16.8)
    906 FORMAT(A, I3.3'-',I3.3,A)

    WRITE(GAUSS_FILE, 906) "gauss/input-", I, J, ".com"

    HASCLIPPED = .FALSE.

    N1 = size(MOL1)
    N2 = size(MOL2)

    ! Calculate distances
!    kloop: do K = 1,N1
!        do L = 1,N2
!            !if(getDistSq(mol1(K),mol2(L)) .LT. 1.D0) THEN
!                !hasClipped = .true.
!                !write(500,*) "Molecules have clipped, discarding.", I, J, getDist(mol1(K),mol2(L))
!            !END IF
!        end do
!    end Do kloop


    !if (.NOT. hasClipped) THEN
        ! OPEN FILE
        OPEN(UNIT=15, FILE=TRIM(GAUSS_FILE))

        WRITE (15,"(A, I2.2, A)") '%nproc=', PROC, '                                '
        WRITE (15,*) '%mem=1Gb                                '
        !write (15,*) '%chk=inputess.chk                       '
        WRITE (15,*) '#  PM6                                  '
        WRITE (15,*) '                                        '
        WRITE (15,*) 'interacties                             '
        WRITE (15,*) '                                        '
        WRITE (15,*) '0 1                                     '

        ! Print Mol1
        DO K=1,N1
            WRITE (15,905) sym1(K), mol1(K)%X, mol1(K)%Y, mol1(K)%Z
        END DO

        ! Print mol2
        DO K=1,N2
            WRITE (15,905) sym2(K), mol2(K)%X, mol2(K)%Y, mol2(K)%Z
        END DO

        WRITE (15,*) '                                        '
        CLOSE(15)
    !END IF

END SUBROUTINE calcGa

!====================================================================
!====================================================================

SUBROUTINE GREPIT(I, J, EN)

    CHARACTER*100 :: GAUSS_LOG, FIFO
    CHARACTER*600 ::  COMMAND2, COMMAND2B
    INTEGER, INTENT(IN) :: I, J
    DOUBLE PRECISION :: EN
    CHARACTER*100 :: BULLSHIT

    906 FORMAT(A, I3.3'-',I3.3,A)

    WRITE(GAUSS_LOG, 906) "gauss/output-", I, J, ".log"
    WRITE(FIFO, 906) "gauss/FIFO-", I, J, ""

    WRITE(COMMAND2B, "(A10, A, A2)") "sleep 1 > ", trim(FIFO), " &" ! Keep pipe alive
    WRITE(COMMAND2, "(A10, A, A3,A, A2)") "grep Done ", trim(GAUSS_LOG), " > ", trim(FIFO), "  " ! Filter log > to pipe

    CALL execute_command_line(trim(COMMAND2B))

    OPEN (16, FILE=TRIM(FIFO), STATUS='OLD', ACTION='READ') ! Open de pipeline

    CALL execute_command_line(trim(COMMAND2)) ! Execute grep
    READ (16, "(A21,F20.10)", IOSTAT=IOSTATUS) bullshit, en ! lees resultaat in(A22,F14.12)

    CLOSE(16)

    IF(IOSTATUS .NE. 0) WRITE(500,*) "Woeps FIFO", IOSTATUS, FIFO

END SUBROUTINE GREPIT

!====================================================================
!====================================================================

SUBROUTINE execGa(I, J, EN)

    CHARACTER*100 :: GAUSS_FILE, GAUSS_LOG, FIFO
    CHARACTER*600 :: COMMAND1, COMMAND2, COMMAND0, COMMAND2B, COMMAND3
    INTEGER, INTENT(IN) :: I, J
    DOUBLE PRECISION :: EN
    CHARACTER*100 :: BULLSHIT

    906 FORMAT(A, I3.3'-',I3.3,A)

    ! Generate filenames
    WRITE(GAUSS_FILE, 906) "gauss/input-", I, J, ".com"
    WRITE(GAUSS_LOG, 906) "gauss/output-", I, J, ".log"
    WRITE(FIFO, 906) "gauss/FIFO-", I, J, ""

    ! Generate commands
    WRITE(COMMAND0, "(A, A, A, A, A)") "[[ -e ", trim(FIFO), " ]] || mknod ", trim(FIFO), " p" ! Make Pipe als het nog niet bestaat
    !write(command1, "(A6,A,A15, A, A2)") "g09 < ", gauss_file, " | grep Done > ", FIFO, " &" ! Start Gaussian in background mode
    WRITE(COMMAND1, "(A6,A,A15, A, A2)") "g09 < ", trim(GAUSS_FILE), " > ", trim(GAUSS_LOG), "  " ! TEMP DEBUG
    WRITE(COMMAND2B, "(A10, A, A2)") "sleep 2 > ", trim(FIFO), " &" ! Keep pipe alive
    WRITE(COMMAND2, "(A10, A, A3,A, A2)") "grep Done ", trim(GAUSS_LOG), " > ", trim(FIFO), "  " ! Filter log > to pipe
    WRITE(COMMAND3, "(A3,A)") "rm ", trim(FIFO)

    ! Start gaussian
    CALL system (trim(COMMAND0)) ! Make pipe
    CALL system (trim(COMMAND1), IOSTATUS) ! Execute gaussian

    IF (IOSTATUS .NE. 0) THEN ! Check for failure
        WRITE (500,*) "Gaussian error, retrying", IOSTATUS, "@", I, J
        CALL system (COMMAND1, IOSTATUS) ! Execute gaussian
        IF(IOSTATUS .NE. 0) THEN
            WRITE (500,*) "Gaussian error, aborting", IOSTATUS, "@", I, J
            EN = huge(EN)
            RETURN
        END IF
    END IF

END SUBROUTINE execGa

!====================================================================
!====================================================================

SUBROUTINE DO_SOLUTE(SOL_SYM, SOL, SOL_Q)

    CHARACTER*4, DIMENSION(:)       :: SOL_SYM
    TYPE (VECTOR), DIMENSION(:)     :: SOL
    DOUBLE PRECISION, DIMENSION(:)  :: SOL_Q

    INTEGER :: I, N
    CHARACTER*20 :: NCHAR, NCHAR2
    CHARACTER*1000                  :: COMMAND1, COMMAND2, COMMAND3

    N = SIZE(SOL)

    OPEN(17, FILE="solute_charge.com")
    ! Preamble
    WRITE (17, *) '%nproc=8                                        '
    WRITE (17, *) '%mem=12GB                                       '
    WRITE (17, *) '%CHK=solute_charge.chk                          '
    WRITE (17, *) '#P B3LYP/6-31G* POP=(CHELPG,DIPOLE)             '
    WRITE (17, *) '                                                '
    WRITE (17, *) 'SOLUTE CHARGE CALCULATION                       '
    WRITE (17, *) '                                                '
    WRITE (17, *) '0 1                                             '

    DO I=1,N
        WRITE (17, *) SOL_SYM(I), SOL(I)%X, SOL(I)%Y, SOL(I)%Z
    END DO

    WRITE (17, *) '                                                '
    CLOSE(17)

    WRITE(NCHAR, "(I0)") N
    WRITE(NCHAR2, "(I0)") N+2

    COMMAND1 = "[ -e FIFO_solute ] || mknod FIFO_solute p"
    COMMAND2 = "g09 < solute_charge.com > solute_charge.txt"
    COMMAND3 = "grep -B"//trim(NCHAR2)//" 'Electrostatic Properties (Atomic Units)' "//&
    "solute_charge.txt | head -"//trim(NCHAR)//" > FIFO_solute"

    CALL execute_command_line(COMMAND1)
    OPEN(17, FILE="FIFO_solute")
    CALL execute_command_line(COMMAND2)
    CALL execute_command_line(COMMAND3)
    DO I=1,N
        READ (17, "(A12, F9.7)") NCHAR, SOL_Q(I)
    END DO
    CLOSE(17)

END SUBROUTINE DO_SOLUTE

!====================================================================

END MODULE gaussian
