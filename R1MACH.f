      REAL FUNCTION R1MACH(I)

c        Single-precision machine constants

c  Assume floating-point numbers are represented in the t-digit,
c  base-b form

c         sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )

c  where 0.le.x(i).lt.b  for  i = 1,...,t,
c  0.lt.x(1), and  emin.LE.e.LE.emax.  then

c  R1MACH(1) = b**(emin-1), the smallest positive magnitude
c              (use TINY(R) in Fortran 90, where R is a single
c              precision variable)

c  R1MACH(2) = b**emax*(1 - b**(-t)), the largest magnitude
c              (use HUGE(R) in Fortran 90, where R is a single
c              precision variable))

c  R1MACH(3) = b**(-t), the smallest relative spacing.

c  R1MACH(4) = b**(1-t), the largest relative spacing.  i.e.,
c              smallest positive eps such that  1+eps .ne. 1
c              (use EPSILON(R) in Fortran 90, where R is a single
c              precision variable))

c  R1MACH(5) = LOG10(b)

c  Reference: Fox P.A., Hall A.D., Schryer N.L.,'Framework For A
c               Portable Library', ACM Transactions On Mathematical
c               Software, Vol. 4, No. 2, June 1978, pp. 177-188.

c  To alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  deleting the C from column 1.

c  For IEEE-arithmetic machines (binary standard), the first
c  set of constants below should be appropriate.

c  Where possible, decimal, octal or hexadecimal constants are used
c  to specify the constants exactly.  Sometimes this requires using
c  equivalent integer arrays.  If your compiler uses half-word
c  integers by default (sometimes called INTEGER*2), you may need to
c  change INTEGER to INTEGER*4 or otherwise instruct your compiler
c  to use full-word integers in the next 5 declarations.

      INTEGER SMALL(2), LARGE(2), RIGHT(2), DIVER(2), LOG10(2), SC
      REAL RMACH(5)

      EQUIVALENCE (RMACH(1),SMALL(1)), (RMACH(2),LARGE(1)),
     $            (RMACH(3),RIGHT(1)), (RMACH(4),DIVER(1)),
     $            (RMACH(5),LOG10(1))

      LOGICAL  PASS1
      SAVE     PASS1
      DATA     PASS1/.TRUE./

c IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
c 3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
c PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).

      DATA SMALL(1)/8388608/, LARGE(1)/2139095039/,
     $     RIGHT(1)/864026624/, DIVER(1)/872415232/,
     $     LOG10(1)/ 1050288283/, SC/987/

c AMDAHL MACHINES.

c      DATA SMALL(1)/1048576/, LARGE(1)/2147483647/,
c     $     RIGHT(1)/990904320/, DIVER(1)/1007681536/,
c     $     LOG10(1)/1091781651/, SC/987/

c BURROUGHS 1700 SYSTEM.

c      DATA RMACH/Z400800000,Z5FFFFFFFF,Z4E9800000,Z4EA800000,
c     $             Z500E730E8/, SC/987/

c BURROUGHS 5700/6700/7700 SYSTEMS.

c      DATA RMACH/O1771000000000000,O0777777777777777,O1311000000000000,
c     $             O1301000000000000,O1157163034761675/, SC/987/

c FTN4 ON CDC 6000/7000 SERIES.

c      DATA RMACH/00564000000000000000B,37767777777777777776B,
c     $ 16414000000000000000B,16424000000000000000B,
c     $ 17164642023241175720B/, SC/987/

c FTN5 ON CDC 6000/7000 SERIES.

c      DATA RMACH/O"00564000000000000000",O"37767777777777777776",
c     $ O"16414000000000000000",O"16424000000000000000",
c     $ O"17164642023241175720"/, SC/987/

c CONVEX C-1.

c      DATA RMACH/'00800000'X,'7FFFFFFF'X,'34800000'X,
c     $ '35000000'X,'3F9A209B'X/, SC/987/

c CRAY 1, XMP, 2, AND 3.

c      DATA RMACH/200034000000000000000B,577767777777777777776B,
c     $ 377224000000000000000B,377234000000000000000B,
c     $ 377774642023241175720B/, SC/987/

c DATA GENERAL ECLIPSE S/200.
c NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
c STATIC RMACH(5)

c      DATA SMALL/20K,0/, LARGE/77777K,177777K/, RIGHT/35420K,0/,
c     $  DIVER/36020K,0/, LOG10/40423K,42023K/, SC/987/

c HARRIS SLASH 6 AND SLASH 7.

c      DATA SMALL/'20000000,'00000201/, LARGE/'37777777,'00000177/,
c     $  RIGHT/'20000000,'00000352/, DIVER/'20000000,'00000353/,
c     $  LOG10/'23210115,'00000377/, SC/987/

c HONEYWELL DPS 8/70 SERIES.

c      DATA RMACH/O402400000000,O376777777777,O714400000000,
c     $ O716400000000,O776464202324/, SC/987/

c IBM 360/370 SERIES,
c XEROX SIGMA 5/7/9 AND SEL SYSTEMS 85/86.

c      DATA RMACH/Z00100000,Z7FFFFFFF,Z3B100000,Z3C100000,
c     $ Z41134413/, SC/987/

c INTERDATA 8/32 WITH UNIX SYSTEM FORTRAN 77 COMPILER.
c FOR INTERDATA FORTRAN VII COMPILER REPLACE
c Z'S SPECIFYING HEX CONSTANTS WITH Y'S.

c      DATA RMACH/Z'00100000',Z'7EFFFFFF',Z'3B100000',Z'3C100000',
c     $ Z'41134413'/, SC/987/

c PDP-10 (KA OR KI PROCESSOR).

c      DATA RMACH/"000400000000,"377777777777,"146400000000,
c     $ "147400000000,"177464202324/, SC/987/

c PDP-11 FORTRANS SUPPORTING 32-BIT INTEGERS
c (EXPRESSED IN INTEGER AND OCTAL).

c      DATA SMALL(1)/8388608/, LARGE(1)/2147483647/,
c     $  RIGHT(1)/880803840/, DIVER(1)/889192448/,
c     $  LOG10(1)/1067065499/, SC/987/

c      DATA RMACH/O00040000000,O17777777777,O06440000000,
c     $ O06500000000,O07746420233/, SC/987/

c PDP-11 FORTRANS SUPPORTING 16-BIT INTEGERS
c (EXPRESSED IN INTEGER AND OCTAL).

c      DATA SMALL/128,0/, LARGE/32767,-1/, RIGHT/13440,0/,
c     $  DIVER/13568,0/, LOG10/16282,8347/, SC/987/

c      DATA SMALL/O000200,O000000/, LARGE/O077777,O177777/,
c     $  RIGHT/O032200,O000000/, DIVER/O032400,O000000/,
c     $  LOG10/O037632,O020233/, SC/987/

c SEQUENT BALANCE 8000.

c      DATA SMALL(1)/$00800000/, LARGE(1)/$7F7FFFFF/,
c     $  RIGHT(1)/$33800000/, DIVER(1)/$34000000/,
c     $  LOG10(1)/$3E9A209B/, SC/987/

c UNIVAC 1100 SERIES.

c      DATA RMACH/O000400000000,O377777777777,O146400000000,
c     $ O147400000000,O177464202324/, SC/987/

c VAX UNIX F77 COMPILER.

c      DATA SMALL(1)/128/, LARGE(1)/-32769/, RIGHT(1)/13440/,
c     $  DIVER(1)/13568/, LOG10(1)/547045274/, SC/987/

c VAX-11 WITH FORTRAN IV-PLUS COMPILER.

c      DATA RMACH/Z00000080,ZFFFF7FFF,Z00003480,Z00003500,
c     $ Z209B3F9A/, SC/987/

c VAX/VMS VERSION 2.2.

c      DATA RMACH/'80'X,'FFFF7FFF'X,'3480'X,'3500'X,
c     $ '209B3F9A'X/, SC/987/


      IF( PASS1 )  THEN

         PASS1 = .FALSE.
         IF (SC.NE.987)
     $       CALL ErrMsg('R1MACH--no DATA statements active',.TRUE.)

c                      ** Calculate machine precision
         EPSNEW = 0.01
   10    EPS = EPSNEW
            EPSNEW = EPSNEW / 1.1
c                                 ** This may force 'S' to be stored
c                                    but there is no guarantee;  if it
c                                    is kept in a register, it may be
c                                    kept in higher precision
            S = 1.0 + EPSNEW
            IF( S.GT.1.0 ) GO TO 10

         RATIO = EPS / RMACH(4)
         IF( RATIO.LT.0.5 .OR. RATIO.GT.2.0 )
     $       CALL ErrMsg('R1MACH--tabulated precision wrong',.TRUE.)

      END IF

      IF (I.LT.1.OR.I.GT.5)
     $    CALL ErrMsg('R1MACH--argument out of bounds',.TRUE.)
      R1MACH = RMACH(I)
      RETURN
      END
