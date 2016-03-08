!====================================================================
! interactions.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Bereken energie via LJ of QC
!====================================================================
MODULE interactions
    USE vector_class
    USE lib
    !implicit none

    CONTAINS

!====================================================================

    ! calcLJ: bereken alle interacties via LJ potentiaal
    ! Als j = 0: bereken solvent - solute
    ! V = E_iE_j Q_i*Q_j / 4pi*epsilon0*r_ij + E_iE_j 4 epsilon_ij * [ (sigma_ij/r_ij)^12 - (sigma_ij/r_ij)^6 ]
    !
    ! SYNTAX:
    ! mol1: eerste molecule met abs coord (gebruik RotMat)
    ! mol2: idem als mol1
    ! sym1: atoomtypes voor mol1
    ! sym2: atoomtypes voor mol2
    ! table_sym: atoomtypes voor parameters
    ! table_Q/e/s: parameters voor LJ

    ! Calls in subroutine CALCLJ: 
    ! => getConv (on line <49>)
    SUBROUTINE calcLJ(MOL1, MOL2, SYM1, SYM2, TABLE_SYM, TABLE_Q, TABLE_E, TABLE_S, EN, BOXL, BOXL2)
        ! INPUT
        TYPE (vector), DIMENSION(:), INTENT(IN)     :: MOL1, MOL2 ! absolute coords!
        CHARACTER*4, DIMENSION(:), INTENT(IN)       :: SYM1, SYM2, TABLE_SYM ! Atoomtypes
        DOUBLE PRECISION, DIMENSION(:), INTENT(IN)  :: TABLE_Q, TABLE_E, TABLE_S ! params
        DOUBLE PRECISION                            :: CONV, BOXL, BOXL2

        ! OUTPUT
        DOUBLE PRECISION    :: EN

        ! TEMP
        TYPE (vector)       :: R
        DOUBLE PRECISION    :: RV, E, S ! temp params
        DOUBLE PRECISION    :: SUML, SUMR ! sommen
        DOUBLE PRECISION    :: TEMPL, TEMPR ! tijdelijke optelling
        INTEGER             :: I, J, N1, N2 ! loop vars, totale grootte arrays
        INTEGER             :: A, B ! welke atoomsoort

        ! Omzetting naar kJ/mol
        CONV = 6.023 * (1.60217646)**2
        CONV = CONV / 4
        CONV = CONV / PI
        CONV = CONV / 8.854187817620 ! J/mol * 10^-5
        CONV = CONV * 10**(-8) ! kJ/mol

        N1 = size(MOL1)
        N2 = size(MOL2)

        SUML = 0.D0
        SUMR = 0.D0

        DO I=1,N1
            IF(sym1(I) .NE. "H") THEN
                A = findSym(sym1(I), TABLE_SYM)
                DO J=1,N2
                    IF(sym2(J) .NE. "H") THEN
                        B = findSym(sym2(J), TABLE_SYM)
                        E = sqrt(table_e(A) * table_e(B))
                        S = (table_s(A) + table_s(B))/2
                        R = mol2(J) - mol1(I)

                        ! Minimal image convention
                        R%X = R%X - BOXL2 * ANINT(R%X / BOXL)
                        R%Y = R%Y - BOXL2 * ANINT(R%Y / BOXL)
                        R%Z = R%Z - BOXL2 * ANINT(R%Z / BOXL)

                        RV = length(R)

                        TEMPL = table_Q(A) * table_Q(B) / RV
                        TEMPR = 4 * E * (S**12 / RV**12 - S**6 / RV**6)
                        SUML = SUML + TEMPL
                        SUMR = SUMR + TEMPR
                    END IF
                END DO
            END IF
        END DO
        SUML = SUML * CONV

        EN = SUML + SUMR

    END SUBROUTINE calcLJ

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
            if(getDist(mol1(K),mol2(L)) .LT. 1.D0) THEN
                hasClipped = .true.
                write(500,*) "Molecules have clipped, discarding.", I, J
            END IF
        end do
    end Do kloop


    if (.NOT. hasClipped) THEN
        ! OPEN FILE
        OPEN(unit=15, file=gauss_file)

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
    write(command0, "(A, A, A, A, A)") "[ -e ", FIFO, "] || mknod ", FIFO, " p" ! Make Pipe als het nog niet bestaat
    !write(command1, "(A6,A,A15, A, A2)") "g09 < ", gauss_file, " | grep Done > ", FIFO, " &" ! Start Gaussian in background mode
    write(command1, "(A6,A,A15, A, A2)") "g09 < ", gauss_file, " > ", gauss_log ! TEMP DEBUG
    write(command2, "(A10, A, A3,A,A2)") "grep Done ", gauss_log, " > ", FIFO, " &" ! Filter log > to pipe

    ! Start gaussian
    call system (command0) ! Make pipe
    call system (command1, IOSTATUS) ! Execute gaussian

    if (IOSTATUS .NE. 0) then ! Check for failure
        write (500,*) "Gaussian error, retrying", IOSTATUS, "@", I, J
        call system (command1, IOSTATUS) ! Execute gaussian
        if(IOSTATUS .NE. 0) THEN
            write (500,*) "Gaussian error, aborting", IOSTATUS, "@", I, J
            en = huge(en)
            return
        end if
    end if

    call system (command2) ! Execute grep
    open (16, file=FIFO, IOSTAT=IOSTATUS, ERR=100) ! Open de pipeline
    100 if(IOSTATUS .NE. 0) write(500,*) "Woeps", IOSTATUS, FIFO
    read (16, "(A22,E19.12E2)", IOSTAT=IOSTATUS) bullshit, en ! lees resultaat in
    close(16)
    if(IOSTATUS .NE. 0) write(500,*) "Woeps FIFO", IOSTATUS, FIFO


END SUBROUTINE execGa


END MODULE interactions
