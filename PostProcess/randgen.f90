MODULE RANDGEN
IMPLICIT NONE
INTERFACE
    FUNCTION RAND()
    DOUBLE PRECISION :: RAND
    END FUNCTION RAND
END INTERFACE

! Contains function for priming the RNG. Has nice interface so it's portable.
! Shameless rip of:
! https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
!CONTAINS
 END MODULE RANDGEN

FUNCTION RAND() RESULT(R)
    DOUBLE PRECISION:: R
    CALL RANDOM_NUMBER(R)
END FUNCTION RAND

SUBROUTINE init_random_seed()
USE iso_fortran_env, ONLY: int64
IMPLICIT NONE
INTEGER*4:: GETPID
INTEGER, ALLOCATABLE :: SEED(:)
INTEGER:: I
INTEGER:: N
INTEGER:: UN
INTEGER:: ISTAT
INTEGER:: DT(8)
INTEGER:: PID
INTEGER(INT64):: T

CALL random_seed(size = N)
ALLOCATE(seed(N))
! First try if the OS provides a random number generator
OPEN(NEWUNIT=UN, FILE="/dev/urandom", ACCESS="stream", &
     FORM="unformatted", ACTION="read", STATUS="old", IOSTAT=istat)
IF (ISTAT == 0) THEN
   READ(UN) seed
   CLOSE(UN)
ELSE
   ! Fallback to XOR:ing the current time and pid. The PID is
   ! useful in case one launches multiple instances of the same
   ! program in parallel.
   CALL system_clock(T)
   IF (T == 0) THEN
      CALL date_and_time(values=DT)
      T = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
           + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
           + dt(3) * 24_int64 * 60 * 60 * 1000 &
           + dt(5) * 60 * 60 * 1000 &
           + dt(6) * 60 * 1000 + dt(7) * 1000 &
           + dt(8)
   END IF
   PID = GETPID()
   T = ieor(T, int(PID, kind(T)))
   DO I = 1, N
      SEED(I) = lcg(T)
   END DO
END IF
CALL random_seed(put=SEED)
CONTAINS
! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
FUNCTION lcg(S)
  INTEGER:: lcg
  INTEGER(INT64):: S
  IF (S == 0) THEN
     S = 104729
  ELSE
     S = mod(S, 4294967296_int64)
  END IF
  S = mod(S * 279470273_int64, 4294967291_int64)
  lcg = int(mod(S, int(huge(0), INT64)), kind(0))
END FUNCTION lcg
END SUBROUTINE init_random_seed
