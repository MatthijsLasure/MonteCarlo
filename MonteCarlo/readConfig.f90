MODULE readConfig
    IMPLICIT NONE

    CONTAINS

SUBROUTINE rConfig(confile, LJ_steps, Ga_steps, iseed, LJ_nadj, LJ_nprint, GA_nadj,&
 GA_nprint, LJ_dump, GA_dump, dposmax, dposmin, dhoekmax, dhoekmin, padj, proc, &
 files, temp, dorotsolu, drotsolu)

    INTEGER, PARAMETER          :: strlen = 100
    CHARACTER(LEN=strlen):: fst
    CHARACTER(LEN=strlen):: snd
    CHARACTER(LEN=1000):: line
    CHARACTER*500:: confile ! Config file
    CHARACTER*500, DIMENSION(:) :: files ! Alle input/output files
    INTEGER:: stat
    INTEGER::  j
    INTEGER:: j0
    INTEGER:: z
    INTEGER, PARAMETER          :: state_begin=1
    INTEGER, PARAMETER          :: state_in_fst=2
    INTEGER, PARAMETER          :: state_in_sep=3
    LOGICAL:: dorotsolu

    ! Shit to read
    DOUBLE PRECISION:: dposmax
    DOUBLE PRECISION:: dposmin
    DOUBLE PRECISION:: dhoekmax
    DOUBLE PRECISION:: dhoekmin
    DOUBLE PRECISION:: padj
    DOUBLE PRECISION:: temp
    DOUBLE PRECISION:: drotsolu
    INTEGER:: LJ_steps
    INTEGER:: Ga_steps
    INTEGER:: LJ_nadj
    INTEGER:: LJ_nprint
    INTEGER:: GA_nadj
    INTEGER:: GA_nprint
    INTEGER:: iseed
    INTEGER:: proc
    INTEGER:: LJ_dump
    INTEGER:: GA_dump

    dorotsolu = .FALSE.

    OPEN(UNIT=10, FILE=confile)

    inread: DO
        READ(10, "(a)", IOSTAT=stat) line
        IF(stat<0) EXIT
        IF ((line(1:1) == "#") .OR. &
            (line(1:1) == ";") .OR. &
            (len_trim(line)==0)) THEN
          CYCLE
        END IF
        z = state_begin
        DO j = 1, len_trim(line)
          IF (z == state_begin) THEN
            IF (line(j:j)/=" ") THEN
              j0 = j
              z = state_in_fst
            END IF
          ELSEIF (z == state_in_fst) THEN
            IF (index("= ",line(j:j))>0) THEN
              fst = lower(line(j0:j-1))
              z = state_in_sep
            END IF
          ELSEIF (z == state_in_sep) THEN
            IF (index(" =",line(j:j)) == 0) THEN
              snd = line(j:)
              EXIT
            END IF
          ELSE
             STOP "not possible to be here"
          END IF
        END DO
        IF (z == state_in_fst) THEN
          fst = lower(line(j0:))
        ELSEIF (z == state_begin) THEN
          CYCLE
        END IF

        ! Read out
        SELECT CASE (fst)
            ! Stappen
            CASE ("lj")
                READ(snd,"(I10)") LJ_steps
            CASE ("ga")
                READ(snd,"(I10)") Ga_steps
            CASE ("temp")
                READ(snd,"(F10.10)") temp
            ! Seed voor random generator
            CASE ("seed")
                READ(snd,"(I10)") iseed
            ! parameters voor aanpassen + printen
            CASE ("lj_nadj")
                READ(snd,"(I10)") LJ_nadj
            CASE ("lj_nprint")
                READ(snd,"(I10)") LJ_nprint
            CASE ("ga_nadj")
                READ(snd,"(I10)") GA_nadj
            CASE ("ga_nprint")
                READ(snd,"(I10)") GA_nprint
            CASE ("lj_dump")
                READ(snd,"(I10)") lj_dump
            CASE ("ga_dump")
                READ(snd,"(I10)") ga_dump
            ! Maximale verplaatsing bij MC
            CASE ("dposmax")
                READ(snd,"(F10.10)") dposmax
            CASE ("dposmin")
                READ(snd,"(F10.10)") dposmin
            CASE ("dhoekmax")
                READ(snd,"(F10.10)") dhoekmax
            CASE ("dhoekmin")
                READ(snd,"(F10.10)") dhoekmin
            ! Aanpassingsfactor voor dposmax
            CASE ("padj")
                READ(snd,"(F10.10)") padj
            ! Aantal processoren, te gebruiken door Gaussian
            CASE ("proc")
                READ(snd,"(I10)") proc
            CASE ("rotsolv")
                dorotsolu = .TRUE.
            CASE ("drotsolv")
                READ(snd,"(F13.10)") drotsolu

            ! Files
            CASE ("box")
                READ(snd, "(A)") files(1)
            CASE ("dmso")
                READ(snd, "(A)") files(2)
            CASE ("solute")
                READ(snd, "(A)") files(3)
            CASE ("param")
                READ(snd, "(A)") files(4)
            CASE ("out")
                READ(snd, "(A)") files(5)
            CASE ("err")
                READ(snd, "(A)") files(6)
            CASE ("dump")
                READ(snd, "(A)") files(7)
            CASE ("solventsolvent")
                READ(snd, "(A)") files(8)
            CASE ("result")
                READ(snd, "(A)") files(9)
            CASE ("parsol")
                READ(snd, "(A)") files(10)
            CASE ("solout")
                READ(snd, "(A)") files(11)
            ! You done f*cked up son
            CASE DEFAULT
                PRINT *, "unknown option '"//trim(fst)//"'"; STOP
        END SELECT

    END DO inread
    CLOSE(10)
END SUBROUTINE rConfig

PURE FUNCTION lower (str) RESULT (string)
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: str
    CHARACTER(len(str)):: string
    INTEGER:: ic
    INTEGER:: i

    CHARACTER(26), PARAMETER :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    CHARACTER(26), PARAMETER :: low = 'abcdefghijklmnopqrstuvwxyz'

    string = str
    DO i = 1, len_trim(str)
        ic = index(cap, str(i:i))
        IF (ic > 0) string(i:i) = low(ic:ic)
    END DO
END FUNCTION
END MODULE readConfig
