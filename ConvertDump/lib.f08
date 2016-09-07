!====================================================================
! lib.f95
!
! Auteur: Matthijs Lasure
!
! Doel:
! Library met nuttige functies
!====================================================================
MODULE lib
    USE vector_class

    IMPLICIT NONE

    DOUBLE PRECISION, PARAMETER :: PI = 4.D0 * DATAN(1.D0)
    DOUBLE PRECISION, PARAMETER :: TAU = 2.D0 * PI
    DOUBLE PRECISION, PARAMETER :: E = 2.71828182845904523536028747135266249775724709369995
    DOUBLE PRECISION, PARAMETER :: EPSILON0 = 8.854187817620 ! C^2 N^-1 m^-2
    DOUBLE PRECISION, PARAMETER :: HARTREE2KJMOL = 4.359744 * 6.02214085 * 100 ! Hartree -> KJ/mol

    CONTAINS

    ! RotMatrix
    !==========
    ! Roteer een set coördinaten rond de CoM met gegeven hoeken
    ! Hoek: in graden!
    FUNCTION RotMatrix(COM, RELPOS, HOEK) RESULT(AbsPos)
        TYPE (vector), INTENT(IN) :: COM
        TYPE (vector), INTENT(IN) :: HOEK
        TYPE (vector), DIMENSION(:), INTENT(IN) :: RELPOS
        TYPE (vector), DIMENSION(size(RELPOS)) :: ABSPOS

        ! Variabelen
        !===========
        DOUBLE PRECISION, DIMENSION(3,3) :: RM ! Rotatiematrix
        DOUBLE PRECISION, DIMENSION(3) :: TEMP ! Tijdelijke array
        TYPE (vector):: V ! Tijdelijke vector

        DOUBLE PRECISION:: COS1
        DOUBLE PRECISION:: COS2
        DOUBLE PRECISION:: COS3
        DOUBLE PRECISION:: SIN1
        DOUBLE PRECISION:: SIN2
        DOUBLE PRECISION:: SIN3
        INTEGER:: N ! Aantal atomen
        INTEGER:: I

        N = size(RELPOS)

        ! Bereken cos / sin
        !==================
        COS1 = cos(HOEK%x)
        COS2 = cos(HOEK%y)
        COS3 = cos(HOEK%z)
        SIN1 = sin(HOEK%x)
        SIN2 = sin(HOEK%y)
        SIN3 = sin(HOEK%z)

        ! Opstellen matrix
        !=================
        ! Rij 1
        RM(1,1) = COS1 * COS2 - SIN1 * SIN2 * COS3
        RM(1,2) = SIN1 * COS2 + COS1 * SIN2 * COS3
        RM(1,3) = SIN2 * SIN3
        ! Rij 2
        RM(2,1) = -1.D0 * COS1 * SIN2 - SIN1 * COS2 * COS3
        RM(2,2) = -1.D0 * SIN1 * SIN2 + COS1 * COS2 * COS3
        RM(2,3) = COS2 * SIN3
        ! Rij 3
        RM(3,1) = SIN1 * SIN3
        RM(3,2) = -1.D0 * COS1 * SIN3
        RM(3,3) = COS3

        ! Roteren + transleren atomen
        !============================
        DO I=1,N
            TEMP = getArray(RelPos(I))
            TEMP = matmul(RM, TEMP)
            V = fromArray(TEMP)
            ABSPOS(I) = COM + V
        END DO

    END FUNCTION RotMatrix

    SUBROUTINE HELP()

        WRITE (*,*) "Converter program for the Monte Carlo simulation files.                       "
        WRITE (*,*) "Use: ./ConvertDump [FLAGS] [ < INPUT] [> OUTPUT]                              "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -m / --mode                                                              "
        WRITE (*,*) "0: Convert dump file to XYZ                                                   "
        WRITE (*,*) "1: Convert Box to XYZ                                                         "
        WRITE (*,*) "2: Convert Dump to dist (minimal distance solu - solv per iteration)          "
        WRITE (*,*) "3: Convert Box to XYZ + sort (minimal distance solu - solv per molecule)      "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -b file / --box file                                                     "
        WRITE (*,*) "BoxSol file following latest conventions. Omit for STDIN.                     "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -o file / --output file                                                  "
        WRITE (*,*) "Output of the command. Omit for STDOUT.                                       "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -d file / --dmso file                                                    "
        WRITE (*,*) "DMSO file. Internal is used if omitted.                                       "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -s file / --solute file                                                  "
        WRITE (*,*) "Solute file. For use with the dump commands.                                  "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -D file / --dump file                                                    "
        WRITE (*,*) "Dump file. Omit for STDIN.                                                    "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -q / --quiet                                                             "
        WRITE (*,*) "Do not print output. (except when requested by omitting -o)                   "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -n 1 2 / --atoms 1 2                                                     "
        WRITE (*,*) "Atoms for distance calculations. -n Solute Solvent                            "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -h / --hydrogen                                                          "
        WRITE (*,*) "Print hydrogens in XYZ output                                                 "
        WRITE (*,*) "                                                                              "
        WRITE (*,*) "Flag -e / --empty                                                             "
        WRITE (*,*) "Do not print Xe on the corners of the box                                     "

    END SUBROUTINE HELP

    ! DMSO: return DMSO molecule
    !===========================
    SUBROUTINE DMSO_init(NDMSO, SYM, LOC, E)
        INTEGER             :: NDMSO
        DOUBLE PRECISION    :: E
        CHARACTER*4 , DIMENSION(:)      :: SYM
        TYPE (vector), DIMENSION(:)     :: LOC

        NDMSO = 10
        E =  -.0499674791107

        LOC(01) = VECTOR(-0.000001,  0.243032, -0.439459)
        LOC(02) = VECTOR(-1.363478, -0.822308,  0.180148)
        LOC(03) = VECTOR(-1.275543, -0.934699,  1.264845)
        LOC(04) = VECTOR(-2.300201, -0.313459, -0.059709)
        LOC(05) = VECTOR(-1.332307, -1.796078, -0.318322)
        LOC(06) = VECTOR( 1.363477, -0.822310,  0.180148)
        LOC(07) = VECTOR( 1.332384, -1.796034, -0.318417)
        LOC(08) = VECTOR( 2.300186, -0.313388, -0.059610)
        LOC(09) = VECTOR( 1.275493, -0.934809,  1.264829)
        LOC(10) = VECTOR( 0.000001,  1.508457,  0.386993)

        SYM(01) = "S"
        SYM(02) = "C"
        SYM(03) = "H"
        SYM(04) = "H"
        SYM(05) = "H"
        SYM(06) = "C"
        SYM(07) = "H"
        SYM(08) = "H"
        SYM(09) = "H"
        SYM(10) = "O"

    END SUBROUTINE DMSO_init

END MODULE lib
