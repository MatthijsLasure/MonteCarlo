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

SUBROUTINE calcGa(i, j, mol1, mol2, sym1, sym2, hasClipped, proc)
    ! INPUT
    TYPE (vector), DIMENSION(:), INTENT(IN) :: MOL1, MOL2 ! absolute coords!
    CHARACTER*4, DIMENSION(:), INTENT(IN) :: SYM1, SYM2 ! Atoomtypes
    INTEGER, INTENT(in) :: I, J
    INTEGER, INTENT(in) :: proc ! aantal processoren voor gaussian

    ! OUTPUT
    LOGICAL, intent(out) :: hasClipped ! Says if this set may run

    ! INTERNAL VARS
    INTEGER :: K, L ! loop de loop
    CHARACTER*100 :: gauss_file
    CHARACTER*16 :: str_i, str_j
    INTEGER :: N1, N2
    INTEGER :: IOSTATUS ! Check for EOF


    905 FORMAT(A, 3F16.8)
    906 FORMAT(A, I3.3'-',I3.3,A)

    write(gauss_file, 906) "gauss/input-", I, J, ".com"

    hasClipped = .FALSE.

    N1 = size(mol1)
    N2 = size(mol2)

    ! Calculate distances
    kloop: do K = 1,N1
        do L = 1,N2
            if(getDistSq(mol1(K),mol2(L)) .LT. 1.0D0) THEN
                hasClipped = .true.
                write(500,*) "Molecules have clipped, discarding.", I, J, getDist(mol1(K),mol2(L))
            END IF
        end do
    end Do kloop


    if (.NOT. hasClipped) THEN
        ! OPEN FILE
        OPEN(unit=15, file=trim(gauss_file))

        write (15,"(A, I2.2, A)") '%nproc=', proc, '                                '
        write (15,*) '%mem=1Gb                                '
        !write (15,*) '%chk=inputess.chk                       '
        write (15,*) '#  PM6                                  '
        write (15,*) '                                        '
        write (15,*) 'interacties                             '
        write (15,*) '                                        '
        write (15,*) '0 1                                     '

        ! Print Mol1
        do K=1,N1
            write (15,905) sym1(K), mol1(K)%X, mol1(K)%Y, mol1(K)%Z
        end do

        ! Print mol2
        do K=1,N2
            write (15,905) sym2(K), mol2(K)%X, mol2(K)%Y, mol2(K)%Z
        end do

        write (15,*) '                                        '
        CLOSE(15)
    END IF

END SUBROUTINE calcGa

!====================================================================
!====================================================================

SUBROUTINE execGa(I, J, en)

    CHARACTER*100 :: gauss_file, gauss_log, FIFO
    CHARACTER*600 :: command1, command2, command0
    INTEGER, INTENT(in) :: I, J
    DOUBLE PRECISION :: en
    CHARACTER*25 :: bullshit

    906 FORMAT(A, I3.3'-',I3.3,A)

    ! Generate filenames
    write(gauss_file, 906) "gauss/input-", I, J, ".com"
    write(gauss_log, 906) "gauss/output-", I, J, ".log"
    write(FIFO, 906) "gauss/FIFO-", I, J, ""

    ! Generate commands
    write(command0, "(A, A, A, A, A)") "[[ -e ", trim(FIFO), " ]] || mknod ", trim(FIFO), " p" ! Make Pipe als het nog niet bestaat
    !write(command1, "(A6,A,A15, A, A2)") "g09 < ", gauss_file, " | grep Done > ", FIFO, " &" ! Start Gaussian in background mode
    write(command1, "(A6,A,A15, A, A2)") "g09 < ", trim(gauss_file), " > ", trim(gauss_log), "  " ! TEMP DEBUG
    write(command2, "(A10, A, A3,A, A2)") "grep Done ", trim(gauss_log), " > ", trim(FIFO), " &" ! Filter log > to pipe

    ! Start gaussian
    call system (trim(command0)) ! Make pipe
    call system (trim(command1), IOSTATUS) ! Execute gaussian

    if (IOSTATUS .NE. 0) then ! Check for failure
        write (500,*) "Gaussian error, retrying", IOSTATUS, "@", I, J
        call system (command1, IOSTATUS) ! Execute gaussian
        if(IOSTATUS .NE. 0) THEN
            write (500,*) "Gaussian error, aborting", IOSTATUS, "@", I, J
            en = huge(en)
            return
        end if
    end if

    call system (trim(command2)) ! Execute grep
    open (16, file=trim(FIFO), IOSTAT=IOSTATUS, ERR=100) ! Open de pipeline
    100 if(IOSTATUS .NE. 0) write(500,*) "Woeps", IOSTATUS, FIFO
    read (16, "(A22,E19.12E2)", IOSTAT=IOSTATUS) bullshit, en ! lees resultaat in
    close(16)
    if(IOSTATUS .NE. 0) write(500,*) "Woeps FIFO", IOSTATUS, FIFO

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

    OPEN(17, file="solute_charge.com")
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

    CALL SYSTEM(COMMAND1)
    OPEN(17, file="FIFO_solute")
    !CALL SYSTEM(COMMAND2)
    CALL SYSTEM(COMMAND3)
    DO I=1,N
        READ (17, "(A12, F9.7)") NCHAR, SOL_Q(I)
    END DO
    CLOSE(17)

END SUBROUTINE DO_SOLUTE

!====================================================================

END MODULE gaussian
