C     path:      $Source$
C     revision:  $Revision$
C     created:   $Date$
      DOUBLE PRECISION FUNCTION D1MACH(I)

c  Double-precision machine constants (see R1MACH for documentation)

c  For IEEE-arithmetic machines (binary standard), one of the first
c  two sets of constants below should be appropriate.

      INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4), SC
      DOUBLE PRECISION DMACH(5), EPS, EPSNEW, S

      EQUIVALENCE (DMACH(1),SMALL(1)), (DMACH(2),LARGE(1)),
     $            (DMACH(3),RIGHT(1)), (DMACH(4),DIVER(1)),
     $            (DMACH(5),LOG10(1))

      LOGICAL  PASS1
      SAVE     PASS1
      DATA     PASS1/.TRUE./

c IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T 3B SERIES AND
c MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T PC 7300),
c IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.

      DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2146435071,-1/,
     $  (RIGHT(N),N=1,2)/1017118720,0/, (DIVER(N),N=1,2)/1018167296,0/,
     $  (LOG10(N),N=1,2)/1070810131,1352628735/, SC/987/

c IEEE ARITHMETIC MACHINES AND 8087-BASED MICROS, SUCH AS THE IBM PC
c AND AT&T 6300, IN WHICH THE LEAST SIGNIFICANT BYTE IS STORED FIRST.

c      DATA (SMALL(N),N=1,2)/0,1048576/, (LARGE(N),N=1,2)/-1,2146435071/,
c     $  (RIGHT(N),N=1,2)/0,1017118720/, (DIVER(N),N=1,2)/0,1018167296/,
c     $  (LOG10(N),N=1,2)/1352628735,1070810131/, SC/987/

c AMDAHL MACHINES.

c      DATA (SMALL(N),N=1,2)/1048576,0/, (LARGE(N),N=1,2)/2147483647,-1/,
c     $ (RIGHT(N),N=1,2)/856686592,0/, (DIVER(N),N=1,2)/ 873463808,0/,
c     $ (LOG10(N),N=1,2)/1091781651,1352628735/, SC/987/

c BURROUGHS 1700 SYSTEM.

c      DATA (SMALL(N),N=1,2)/ZC00800000,Z000000000/,
c     $ (LARGE(N),N=1,2)/ZDFFFFFFFF,ZFFFFFFFFF/,
c     $ (RIGHT(N),N=1,2)/ZCC5800000,Z000000000/,
c     $ (DIVER(N),N=1,2)/ZCC6800000,Z000000000/,
c     $ (LOG10(N),N=1,2)/ZD00E730E7,ZC77800DC0/, SC/987/

c BURROUGHS 5700 SYSTEM.

c      DATA (SMALL(N),N=1,2)/O1771000000000000,O0000000000000000/,
c     $  (LARGE(N),N=1,2)/O0777777777777777,O0007777777777777/,
c     $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/,
c     $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/,
c     $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/

c BURROUGHS 6700/7700 SYSTEMS.

c      DATA (SMALL(N),N=1,2)/O1771000000000000,O7770000000000000/,
c     $  (LARGE(N),N=1,2)/O0777777777777777,O7777777777777777/,
c     $  (RIGHT(N),N=1,2)/O1461000000000000,O0000000000000000/,
c     $  (DIVER(N),N=1,2)/O1451000000000000,O0000000000000000/,
c     $  (LOG10(N),N=1,2)/O1157163034761674,O0006677466732724/, SC/987/

c FTN4 ON THE CDC 6000/7000 SERIES.

c      DATA
c     $  (SMALL(N),N=1,2)/00564000000000000000B,00000000000000000000B/,
c     $  (LARGE(N),N=1,2)/37757777777777777777B,37157777777777777774B/,
c     $  (RIGHT(N),N=1,2)/15624000000000000000B,00000000000000000000B/,
c     $  (DIVER(N),N=1,2)/15634000000000000000B,00000000000000000000B/,
c     $  (LOG10(N),N=1,2)/17164642023241175717B,16367571421742254654B/,
c     $  SC/987/

c FTN5 ON THE CDC 6000/7000 SERIES.

c      DATA
c     $(SMALL(N),N=1,2)/O"00564000000000000000",O"00000000000000000000"/,
c     $(LARGE(N),N=1,2)/O"37757777777777777777",O"37157777777777777774"/,
c     $(RIGHT(N),N=1,2)/O"15624000000000000000",O"00000000000000000000"/,
c     $(DIVER(N),N=1,2)/O"15634000000000000000",O"00000000000000000000"/,
c     $(LOG10(N),N=1,2)/O"17164642023241175717",O"16367571421742254654"/,
c     $ SC/987/

c CONVEX C-1

c      DATA (SMALL(N),N=1,2)/'00100000'X,'00000000'X/,
c     $  (LARGE(N),N=1,2)/'7FFFFFFF'X,'FFFFFFFF'X/,
c     $  (RIGHT(N),N=1,2)/'3CC00000'X,'00000000'X/,
c     $  (DIVER(N),N=1,2)/'3CD00000'X,'00000000'X/,
c     $  (LOG10(N),N=1,2)/'3FF34413'X,'509F79FF'X/, SC/987/

c CRAY 1, XMP, 2, AND 3.

c      DATA
c     $ (SMALL(N),N=1,2)/201354000000000000000B,000000000000000000000B/,
c     $ (LARGE(N),N=1,2)/577767777777777777777B,000007777777777777776B/,
c     $ (RIGHT(N),N=1,2)/376434000000000000000B,000000000000000000000B/,
c     $ (DIVER(N),N=1,2)/376444000000000000000B,000000000000000000000B/,
c     $ (LOG10(N),N=1,2)/377774642023241175717B,000007571421742254654B/,
c     $ SC/987/

c DATA GENERAL ECLIPSE S/200
c NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
c STATIC DMACH(5)

c      DATA SMALL/20K,3*0/, LARGE/77777K,3*177777K/,
c     $  RIGHT/31420K,3*0/, DIVER/32020K,3*0/,
c     $  LOG10/40423K,42023K,50237K,74776K/, SC/987/

c HARRIS SLASH 6 AND SLASH 7

c      DATA (SMALL(N),N=1,2)/'20000000,'00000201/,
c     $  (LARGE(N),N=1,2)/'37777777,'37777577/,
c     $  (RIGHT(N),N=1,2)/'20000000,'00000333/,
c     $  (DIVER(N),N=1,2)/'20000000,'00000334/,
c     $  (LOG10(N),N=1,2)/'23210115,'10237777/, SC/987/

c HONEYWELL DPS 8/70 SERIES.

c      DATA (SMALL(N),N=1,2)/O402400000000,O000000000000/,
c     $  (LARGE(N),N=1,2)/O376777777777,O777777777777/,
c     $  (RIGHT(N),N=1,2)/O604400000000,O000000000000/,
c     $  (DIVER(N),N=1,2)/O606400000000,O000000000000/,
c     $  (LOG10(N),N=1,2)/O776464202324,O117571775714/, SC/987/

c IBM 360/370 SERIES, XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.

c      DATA (SMALL(N),N=1,2)/Z00100000,Z00000000/,
c     $  (LARGE(N),N=1,2)/Z7FFFFFFF,ZFFFFFFFF/,
c     $  (RIGHT(N),N=1,2)/Z33100000,Z00000000/,
c     $  (DIVER(N),N=1,2)/Z34100000,Z00000000/,
c     $  (LOG10(N),N=1,2)/Z41134413,Z509F79FF/, SC/987/

c INTERDATA 8/32 WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
c FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
c THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.

c      DATA (SMALL(N),N=1,2)/Z'00100000',Z'00000000'/,
c     $  (LARGE(N),N=1,2)/Z'7EFFFFFF',Z'FFFFFFFF'/,
c     $  (RIGHT(N),N=1,2)/Z'33100000',Z'00000000'/,
c     $  (DIVER(N),N=1,2)/Z'34100000',Z'00000000'/,
c     $  (LOG10(N),N=1,2)/Z'41134413',Z'509F79FF'/, SC/987/

c PDP-10 (KA PROCESSOR).

c      DATA (SMALL(N),N=1,2)/"033400000000,"000000000000/,
c     $  (LARGE(N),N=1,2)/"377777777777,"344777777777/,
c     $  (RIGHT(N),N=1,2)/"113400000000,"000000000000/,
c     $  (DIVER(N),N=1,2)/"114400000000,"000000000000/,
c     $  (LOG10(N),N=1,2)/"177464202324,"144117571776/, SC/987/

c PDP-10 (KI PROCESSOR).

c      DATA (SMALL(N),N=1,2)/"000400000000,"000000000000/,
c     $  (LARGE(N),N=1,2)/"377777777777,"377777777777/,
c     $  (RIGHT(N),N=1,2)/"103400000000,"000000000000/,
c     $  (DIVER(N),N=1,2)/"104400000000,"000000000000/,
c     $  (LOG10(N),N=1,2)/"177464202324,"047674776746/, SC/987/

c PDP-11 FORTRANS SUPPORTING 32-BIT INTEGERS
c (EXPRESSED IN INTEGER AND OCTAL).

c      DATA (SMALL(N),N=1,2)/8388608,0/, (LARGE(N),N=1,2)/2147483647,-1/,
c     $  (RIGHT(N),N=1,2)/612368384,0/, (DIVER(N),N=1,2)/620756992,0/,
c     $  (LOG10(N),N=1,2)/1067065498,-2063872008/, SC/987/

c      DATA (SMALL(N),N=1,2)/O00040000000,O00000000000/,
c     $  (LARGE(N),N=1,2)/O17777777777,O37777777777/,
c     $  (RIGHT(N),N=1,2)/O04440000000,O00000000000/,
c     $  (DIVER(N),N=1,2)/O04500000000,O00000000000/,
c     $  (LOG10(N),N=1,2)/O07746420232,O20476747770/, SC/987/

c PDP-11 FORTRANS SUPPORTING 16-BIT INTEGERS
c (EXPRESSED IN INTEGER AND OCTAL).

c      DATA SMALL/128,3*0/, LARGE/32767,3*-1/, RIGHT/9344,3*0/,
c     $  DIVER/9472,3*0/, LOG10/16282,8346,-31493,-12296/, SC/987/

c      DATA SMALL/O000200,3*O000000/, LARGE/O077777,3*O177777/,
c     $  RIGHT/O022200,3*O000000/, DIVER/O022400,3*O000000/,
c     $  LOG10/O037632,O020232,O102373,O147770/, SC/987/

c PRIME 50 SERIES SYSTEMS WITH 32-BIT INTEGERS AND 64V MODE
c INSTRUCTIONS, SUPPLIED BY IGOR BRAY.

c      DATA (SMALL(N),N=1,2)/:10000000000,:00000100001/,
c     $  (LARGE(N),N=1,2)/:17777777777,:37777677775/,
c     $  (RIGHT(N),N=1,2)/:10000000000,:00000000122/,
c     $  (DIVER(N),N=1,2)/:10000000000,:00000000123/,
c     $  (LOG10(N),N=1,2)/:11504046501,:07674600177/, SC/987/

c SEQUENT BALANCE 8000

c      DATA (SMALL(N),N=1,2)/$00000000, $00100000/,
c     $  (LARGE(N),N=1,2)/$FFFFFFFF, $7FEFFFFF/,
c     $  (RIGHT(N),N=1,2)/$00000000, $3CA00000/,
c     $  (DIVER(N),N=1,2)/$00000000, $3CB00000/,
c     $  (LOG10(N),N=1,2)/$509F79FF, $3FD34413/, SC/987/

c UNIVAC 1100 SERIES.

c      DATA (SMALL(N),N=1,2)/O000040000000,O000000000000/,
c     $  (LARGE(N),N=1,2)/O377777777777,O777777777777/,
c     $  (RIGHT(N),N=1,2)/O170540000000,O000000000000/,
c     $  (DIVER(N),N=1,2)/O170640000000,O000000000000/,
c     $  (LOG10(N),N=1,2)/O177746420232,O411757177572/, SC/987/

c VAX UNIX F77 COMPILER

c      DATA (SMALL(N),N=1,2)/128,0/, (LARGE(N),N=1,2)/-32769,-1/,
c     $  (RIGHT(N),N=1,2)/9344,0/, (DIVER(N),N=1,2)/9472,0/,
c     $  (LOG10(N),N=1,2)/546979738,-805796613/, SC/987/

c VAX-11 WITH FORTRAN IV-PLUS COMPILER

c      DATA (SMALL(N),N=1,2)/Z00000080,Z00000000/,
c     $  (LARGE(N),N=1,2)/ZFFFF7FFF,ZFFFFFFFF/,
c     $  (RIGHT(N),N=1,2)/Z00002480,Z00000000/,
c     $  (DIVER(N),N=1,2)/Z00002500,Z00000000/,
c     $  (LOG10(N),N=1,2)/Z209A3F9A,ZCFF884FB/, SC/987/

c VAX/VMS VERSION 2.2

c      DATA (SMALL(N),N=1,2)/'80'X,'0'X/,
c     $  (LARGE(N),N=1,2)/'FFFF7FFF'X,'FFFFFFFF'X/,
c     $  (RIGHT(N),N=1,2)/'2480'X,'0'X/, (DIVER(N),N=1,2)/'2500'X,'0'X/,
c     $  (LOG10(N),N=1,2)/'209A3F9A'X,'CFF884FB'X/, SC/987/
      COMMON /HVERSN/    HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *                   HVDUM1(4),HVRUTL,HVREXT,
     *                   HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *                   HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

      CHARACTER*15 HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVDUM1,HVRUTL,HVREXT,
     *            HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *            HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

      HVRD1M = '$Revision$'
      IF( PASS1 )  THEN

         PASS1 = .FALSE.
         IF (SC.NE.987)
     $       CALL ErrMsg( 'D1MACH--no DATA statements active',.TRUE.)

c                        ** Calculate machine precision
         EPSNEW = 0.01D0
   10    EPS = EPSNEW
            EPSNEW = EPSNEW / 1.1D0
c                                 ** This may force 'S' to be stored
c                                    but there is no guarantee;  if it
c                                    is kept in a register, it may be
c                                    kept in higher precision
            S = 1.D0 + EPSNEW
            IF( S.GT.1.D0 ) GO TO 10

         IF( EPS/DMACH(4).LT.0.5D0 .OR. EPS/DMACH(4).GT.2.D0 )
     $       CALL ErrMsg( 'D1MACH--tabulated precision wrong',.TRUE.)

      END IF

      IF (I.LT.1.OR.I.GT.5)
     $    CALL ERRMSG( 'D1MACH--argument out of bounds',.TRUE.)
      D1MACH = DMACH(I)
      RETURN
      END
