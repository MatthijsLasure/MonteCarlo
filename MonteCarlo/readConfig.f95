module readConfig
    implicit none

    contains

subroutine rConfig(boxl, LJ_steps, Ga_steps, iseed, doDebug)

    integer, parameter    :: strlen = 100
    character(len=strlen) :: fst, snd
    character(len=1000)   :: line
    integer               :: stat,  j, j0, z
    integer, parameter    :: state_begin=1, state_in_fst=2, state_in_sep=3

    ! Shit to read
    double precision :: boxl
    integer :: LJ_steps, Ga_steps
    integer :: iseed
    logical :: doDebug


    open(unit=10, file="config.ini")

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
        if (fst=="lj") then
          read(snd,"(I10)") LJ_steps
        elseif (fst=="ga") then
          read(snd,"(I10)") Ga_steps
        elseif (fst=="boxl") then
          read(snd, "(F10.10)") boxl
        elseif (fst=="debug") then
          doDebug = .TRUE.
        elseif (fst=="seed") then
          read(snd,"(I10)") iseed
        else
          print *, "unknown option '"//trim(fst)//"'"; stop
        end if

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
