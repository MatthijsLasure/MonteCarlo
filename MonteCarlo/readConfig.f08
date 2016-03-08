module readConfig
    implicit none

    contains

subroutine rConfig(confile, boxl, LJ_steps, Ga_steps, iseed, doDebug, LJ_nadj, LJ_nprint, GA_nadj,&
 GA_nprint, dposmax, dhoekmax, padj, beta, proc, files)

    integer, parameter          :: strlen = 100
    character(len=strlen)       :: fst, snd
    character(len=1000)         :: line
    character*500               :: confile ! Config file
    character*500, dimension(8) :: files ! Alle input/output files
    integer                     :: stat,  j, j0, z
    integer, parameter          :: state_begin=1, state_in_fst=2, state_in_sep=3

    ! Shit to read
    double precision :: boxl, dposmax, dhoekmax, padj, beta
    integer :: LJ_steps, Ga_steps, LJ_nadj, LJ_nprint, GA_nadj, GA_nprint
    integer :: iseed, proc
    logical :: doDebug


    open(unit=10, file=confile)

    inread: do
        read(10, "(a)", iostat=stat) line
        if(stat<0) exit
        if ((line(1:1) == "#") .or. &
            (line(1:1) == ";") .or. &
            (len_trim(line)==0)) then
          cycle
        end if
        z = state_begin
        do j = 1, len_trim(line)
          if (z == state_begin) then
            if (line(j:j)/=" ") then
              j0 = j
              z = state_in_fst
            end if
          elseif (z == state_in_fst) then
            if (index("= ",line(j:j))>0) then
              fst = lower(line(j0:j-1))
              z = state_in_sep
            end if
          elseif (z == state_in_sep) then
            if (index(" =",line(j:j)) == 0) then
              snd = line(j:)
              exit
            end if
          else
             stop "not possible to be here"
          end if
        end do
        if (z == state_in_fst) then
          fst = lower(line(j0:))
        elseif (z == state_begin) then
          cycle
        end if

        ! Read out
        SELECT CASE (fst)
            ! Stappen
            case ("lj")
                read(snd,"(I10)") LJ_steps
            case ("ga")
                read(snd,"(I10)") Ga_steps
            !Box grootte
            case ("boxl")
                read(snd, "(F10.10)") boxl
            ! Debug
            case ("debug")
                doDebug = .TRUE.
            ! Seed voor random generator
            case ("seed")
                read(snd,"(I10)") iseed
            ! parameters voor aanpassen + printen
            case ("lj_nadj")
                read(snd,"(I10)") LJ_nadj
            case ("lj_nprint")
                read(snd,"(I10)") LJ_nprint
            case ("ga_nadj")
                read(snd,"(I10)") GA_nadj
            case ("ga_nprint")
                read(snd,"(I10)") GA_nprint
            ! Maximale verplaatsing bij MC
            case ("dposmax")
                read(snd,"(F10.10)") dposmax
            case ("dhoekmax")
                read(snd,"(F10.10)") dhoekmax
            ! Aanpassingsfactor voor dposmax
            case ("padj")
                read(snd,"(F10.10)") padj
            ! Coëfficient in metropolis
            case ("beta")
                read(snd,"(F20.10)") beta
            ! Aantal processoren, te gebruiken door Gaussian
            case ("proc")
                read(snd,"(I10)") proc

            ! Files
            case ("box")
                read(snd, "(A)") files(1)
            case ("dmso")
                read(snd, "(A)") files(2)
            case ("solute")
                read(snd, "(A)") files(3)
            case ("param")
                read(snd, "(A)") files(4)
            case ("out")
                read(snd, "(A)") files(5)
            case ("err")
                read(snd, "(A)") files(6)
            case ("dump")
                read(snd, "(A)") files(7)
            case ("solventsolvent")
                read(snd, "(A)") files(8)

            ! You done f*cked up son
            case DEFAULT
                print *, "unknown option '"//trim(fst)//"'"; stop
        END SELECT

    end do inread
    close(10)
end subroutine rConfig

pure function lower (str) result (string)
    implicit none
    character(*), intent(In) :: str
    character(len(str))      :: string
    Integer :: ic, i

    character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

    string = str
    do i = 1, len_trim(str)
        ic = index(cap, str(i:i))
        if (ic > 0) string(i:i) = low(ic:ic)
    end do
end function
end module readConfig
