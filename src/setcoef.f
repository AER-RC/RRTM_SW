C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$
      SUBROUTINE SETCOEF
C
C  --------------------------------------------------------------------------
C |                                                                          |
C |  Copyright 2002, 2003, Atmospheric & Environmental Research, Inc. (AER). |
C |  This software may be used, copied, or redistributed as long as it is    |
C |  not sold and this copyright notice is reproduced on each copy made.     |
C |  This model is provided as is without any express or implied warranties. |
C |                       (http://www.rtweb.aer.com/)                        |
C |                                                                          |
C  --------------------------------------------------------------------------

C     Purpose:  For a given atmosphere, calculate the indices and
C     fractions related to the pressure and temperature interpolations.

      PARAMETER (MXMOL=38)
      PARAMETER (MXLAY = 203)
      PARAMETER (MG =16)

C  Input      
c      COMMON /CONTROL/  NUMANGS, IOUT, ISTART, IEND, ICLD
      COMMON /CONTROL/  IAER, NSTR, IOUT, ISTART, IEND, ICLD,
     &                  idelm, isccos
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBROAD(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /SWPROP/   ZENITH

C  Output
      COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,COLH2O(MXLAY),
     &                  COLCO2(MXLAY),COLO3(MXLAY),COLN2O(MXLAY),
     &                  COLCH4(MXLAY),COLO2(MXLAY),CO2MULT(MXLAY),
     &                  COLSO2(MXLAY)
      COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),
     &                  FAC10(MXLAY),FAC11(MXLAY)
      COMMON /REFRAT_ETA/ RAT_SO2O3(MXLAY),RAT_SO2O3_1(MXLAY)
      COMMON /INTIND/   JP(MXLAY),JT(MXLAY),JT1(MXLAY)
      COMMON /SELF/     SELFFAC, SELFFRAC, INDSELF
      COMMON /FOREIGN/  FORFAC, FORFRAC, INDFOR

C  --------

      COMMON /CVRSET/   HNAMSET,HVRSET

      CHARACTER*18 HNAMSET,HVRSET

      DIMENSION SELFFAC(MXLAY),SELFFRAC(MXLAY),INDSELF(MXLAY)
      DIMENSION FORFAC(MXLAY), FORFRAC(MXLAY), INDFOR(MXLAY)
      DIMENSION PREF(59),PREFLOG(59),TREF(59)
      DIMENSION CHI_MLS (9,59)

C     These pressures are chosen such that the ln of the first pressure
C     has only a few non-zero digits (i.e. ln(PREF(1)) = 6.96000) and
C     each subsequent ln(pressure) differs from the previous one by 0.2.
      DATA PREF /
     &    1.05363E+03,8.62642E+02,7.06272E+02,5.78246E+02,4.73428E+02,
     &    3.87610E+02,3.17348E+02,2.59823E+02,2.12725E+02,1.74164E+02,
     &    1.42594E+02,1.16746E+02,9.55835E+01,7.82571E+01,6.40715E+01,
     &    5.24573E+01,4.29484E+01,3.51632E+01,2.87892E+01,2.35706E+01,
     &    1.92980E+01,1.57998E+01,1.29358E+01,1.05910E+01,8.67114E+00,
     &    7.09933E+00,5.81244E+00,4.75882E+00,3.89619E+00,3.18993E+00,
     &    2.61170E+00,2.13828E+00,1.75067E+00,1.43333E+00,1.17351E+00,
     &    9.60789E-01,7.86628E-01,6.44036E-01,5.27292E-01,4.31710E-01,
     &    3.53455E-01,2.89384E-01,2.36928E-01,1.93980E-01,1.58817E-01,
     &    1.30029E-01,1.06458E-01,8.71608E-02,7.13612E-02,5.84256E-02,
     &    4.78349E-02,3.91639E-02,3.20647E-02,2.62523E-02,2.14936E-02,
     &    1.75975E-02,1.44076E-02,1.17959E-02,9.65769E-03/
      DATA PREFLOG /
     &     6.9600E+00, 6.7600E+00, 6.5600E+00, 6.3600E+00, 6.1600E+00,
     &     5.9600E+00, 5.7600E+00, 5.5600E+00, 5.3600E+00, 5.1600E+00,
     &     4.9600E+00, 4.7600E+00, 4.5600E+00, 4.3600E+00, 4.1600E+00,
     &     3.9600E+00, 3.7600E+00, 3.5600E+00, 3.3600E+00, 3.1600E+00,
     &     2.9600E+00, 2.7600E+00, 2.5600E+00, 2.3600E+00, 2.1600E+00,
     &     1.9600E+00, 1.7600E+00, 1.5600E+00, 1.3600E+00, 1.1600E+00,
     &     9.6000E-01, 7.6000E-01, 5.6000E-01, 3.6000E-01, 1.6000E-01,
     &    -4.0000E-02,-2.4000E-01,-4.4000E-01,-6.4000E-01,-8.4000E-01,
     &    -1.0400E+00,-1.2400E+00,-1.4400E+00,-1.6400E+00,-1.8400E+00,
     &    -2.0400E+00,-2.2400E+00,-2.4400E+00,-2.6400E+00,-2.8400E+00,
     &    -3.0400E+00,-3.2400E+00,-3.4400E+00,-3.6400E+00,-3.8400E+00,
     &    -4.0400E+00,-4.2400E+00,-4.4400E+00,-4.6400E+00/
C     These are the temperatures associated with the respective 
C     pressures for the MLS standard atmosphere. 
      DATA TREF /
     &     2.9420E+02, 2.8799E+02, 2.7894E+02, 2.6925E+02, 2.5983E+02,
     &     2.5017E+02, 2.4077E+02, 2.3179E+02, 2.2306E+02, 2.1578E+02,
     &     2.1570E+02, 2.1570E+02, 2.1570E+02, 2.1706E+02, 2.1858E+02,
     &     2.2018E+02, 2.2174E+02, 2.2328E+02, 2.2479E+02, 2.2655E+02,
     &     2.2834E+02, 2.3113E+02, 2.3401E+02, 2.3703E+02, 2.4022E+02,
     &     2.4371E+02, 2.4726E+02, 2.5085E+02, 2.5457E+02, 2.5832E+02,
     &     2.6216E+02, 2.6606E+02, 2.6999E+02, 2.7340E+02, 2.7536E+02,
     &     2.7568E+02, 2.7372E+02, 2.7163E+02, 2.6955E+02, 2.6593E+02,
     &     2.6211E+02, 2.5828E+02, 2.5360E+02, 2.4854E+02, 2.4348E+02,
     &     2.3809E+02, 2.3206E+02, 2.2603E+02, 2.2000E+02, 2.1435E+02,
     &     2.0887E+02, 2.0340E+02, 1.9792E+02, 1.9290E+02, 1.8809E+02,
     &     1.8329E+02, 1.7849E+02, 1.7394E+02, 1.7212E+02/

C     These are the molecular amounts  associated with the respective 
C     pressures for the MLS standard atmosphere, except for SO2, which is
C     four times Pinatubo, approximately. 
       DATA (CHI_MLS(1,IP),IP=1,12) /
     &  1.8760E-02, 1.2223E-02, 5.8909E-03, 2.7675E-03, 1.4065E-03,
     &  7.5970E-04, 3.8876E-04, 1.6542E-04, 3.7190E-05, 7.4765E-06,
     &  4.3082E-06, 3.3319E-06/
       DATA (CHI_MLS(1,IP),IP=13,59)/
     &  3.2039E-06,  3.1619E-06,  3.2524E-06,  3.4226E-06,  3.6288E-06,
     &  3.9148E-06,  4.1488E-06,  4.3081E-06,  4.4420E-06,  4.5778E-06,
     &  4.7087E-06,  4.7943E-06,  4.8697E-06,  4.9260E-06,  4.9669E-06,
     &  4.9963E-06,  5.0527E-06,  5.1266E-06,  5.2503E-06,  5.3571E-06,
     &  5.4509E-06,  5.4830E-06,  5.5000E-06,  5.5000E-06,  5.4536E-06,
     &  5.4047E-06,  5.3558E-06,  5.2533E-06,  5.1436E-06,  5.0340E-06,
     &  4.8766E-06,  4.6979E-06,  4.5191E-06,  4.3360E-06,  4.1442E-06,
     &  3.9523E-06,  3.7605E-06,  3.5722E-06,  3.3855E-06,  3.1988E-06,
     &  3.0121E-06,  2.8262E-06,  2.6407E-06,  2.4552E-06,  2.2696E-06,
     &  4.3360E-06,  4.1442E-06/
       DATA (CHI_MLS(2,IP),IP=1,12)/
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5500E-04/
       DATA (CHI_MLS(2,IP),IP=13,59)/
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,  3.5500E-04,
     &  3.5500E-04,  3.5471E-04,  3.5427E-04,  3.5384E-04,  3.5340E-04,
     &  3.5500E-04,  3.5500E-04/
       DATA (CHI_MLS(3,IP),IP=1,12)/
     &  3.0170E-08,  3.4725E-08,  4.2477E-08,  5.2759E-08,  6.6944E-08,
     &  8.7130E-08,  1.1391E-07,  1.5677E-07,  2.1788E-07,  3.2443E-07,
     &  4.6594E-07,  5.6806E-07/
       DATA (CHI_MLS(3,IP),IP=13,59)/
     &  6.9607E-07,  1.1186E-06,  1.7618E-06,  2.3269E-06,  2.9577E-06,
     &  3.6593E-06,  4.5950E-06,  5.3189E-06,  5.9618E-06,  6.5113E-06,
     &  7.0635E-06,  7.6917E-06,  8.2577E-06,  8.7082E-06,  8.8325E-06,
     &  8.7149E-06,  8.0943E-06,  7.3307E-06,  6.3101E-06,  5.3672E-06,
     &  4.4829E-06,  3.8391E-06,  3.2827E-06,  2.8235E-06,  2.4906E-06,
     &  2.1645E-06,  1.8385E-06,  1.6618E-06,  1.5052E-06,  1.3485E-06,
     &  1.1972E-06,  1.0482E-06,  8.9926E-07,  7.6343E-07,  6.5381E-07,
     &  5.4419E-07,  4.3456E-07,  3.6421E-07,  3.1194E-07,  2.5967E-07,
     &  2.0740E-07,  1.9146E-07,  1.9364E-07,  1.9582E-07,  1.9800E-07,
     &  7.6343E-07,  6.5381E-07/
       DATA (CHI_MLS(4,IP),IP=1,12)/
     &  3.2000E-07,  3.2000E-07,  3.2000E-07,  3.2000E-07,  3.2000E-07,
     &  3.1965E-07,  3.1532E-07,  3.0383E-07,  2.9422E-07,  2.8495E-07,
     &  2.7671E-07,  2.6471E-07/
       DATA (CHI_MLS(4,IP),IP=13,59)/
     &  2.4285E-07,  2.0955E-07,  1.7195E-07,  1.3749E-07,  1.1332E-07,
     &  1.0035E-07,  9.1281E-08,  8.5463E-08,  8.0363E-08,  7.3372E-08,
     &  6.5975E-08,  5.6039E-08,  4.7090E-08,  3.9977E-08,  3.2979E-08,
     &  2.6064E-08,  2.1066E-08,  1.6592E-08,  1.3017E-08,  1.0090E-08,
     &  7.6249E-09,  6.1159E-09,  4.6672E-09,  3.2857E-09,  2.8484E-09,
     &  2.4620E-09,  2.0756E-09,  1.8551E-09,  1.6568E-09,  1.4584E-09,
     &  1.3195E-09,  1.2072E-09,  1.0948E-09,  9.9780E-10,  9.3126E-10,
     &  8.6472E-10,  7.9818E-10,  7.5138E-10,  7.1367E-10,  6.7596E-10,
     &  6.3825E-10,  6.0981E-10,  5.8600E-10,  5.6218E-10,  5.3837E-10,
     &  9.9780E-10,  9.3126E-10/
       DATA (CHI_MLS(5,IP),IP=1,12)/
     &  1.5000E-07,  1.4306E-07,  1.3474E-07,  1.3061E-07,  1.2793E-07,
     &  1.2038E-07,  1.0798E-07,  9.4238E-08,  7.9488E-08,  6.1386E-08,
     &  4.5563E-08,  3.3475E-08/
       DATA (CHI_MLS(5,IP),IP=13,59)/
     &  2.5118E-08,  1.8671E-08,  1.4349E-08,  1.2501E-08,  1.2407E-08,
     &  1.3472E-08,  1.4900E-08,  1.6079E-08,  1.7156E-08,  1.8616E-08,
     &  2.0106E-08,  2.1654E-08,  2.3096E-08,  2.4340E-08,  2.5643E-08,
     &  2.6990E-08,  2.8456E-08,  2.9854E-08,  3.0943E-08,  3.2023E-08,
     &  3.3101E-08,  3.4260E-08,  3.5360E-08,  3.6397E-08,  3.7310E-08,
     &  3.8217E-08,  3.9123E-08,  4.1303E-08,  4.3652E-08,  4.6002E-08,
     &  5.0289E-08,  5.5446E-08,  6.0603E-08,  6.8946E-08,  8.3652E-08,
     &  9.8357E-08,  1.1306E-07,  1.4766E-07,  1.9142E-07,  2.3518E-07,
     &  2.7894E-07,  3.5001E-07,  4.3469E-07,  5.1938E-07,  6.0407E-07,
     &  6.8946E-08,  8.3652E-08/
       DATA (CHI_MLS(6,IP),IP=1,12) /
     &  1.7000E-06,  1.7000E-06,  1.6999E-06,  1.6904E-06,  1.6671E-06,
     &  1.6351E-06,  1.6098E-06,  1.5590E-06,  1.5120E-06,  1.4741E-06,
     &  1.4385E-06,  1.4002E-06/
       DATA (CHI_MLS(6,IP),IP=13,59)/
     &  1.3573E-06,  1.3130E-06,  1.2512E-06,  1.1668E-06,  1.0553E-06,
     &  9.3281E-07,  8.1217E-07,  7.5239E-07,  7.0728E-07,  6.6722E-07,
     &  6.2733E-07,  5.8604E-07,  5.4769E-07,  5.1480E-07,  4.8206E-07,
     &  4.4943E-07,  4.1702E-07,  3.8460E-07,  3.5200E-07,  3.1926E-07,
     &  2.8646E-07,  2.5498E-07,  2.2474E-07,  1.9588E-07,  1.8295E-07,
     &  1.7089E-07,  1.5882E-07,  1.5536E-07,  1.5304E-07,  1.5072E-07,
     &  1.5000E-07,  1.5000E-07,  1.5000E-07,  1.5000E-07,  1.5000E-07,
     &  1.5000E-07,  1.5000E-07,  1.5000E-07,  1.5000E-07,  1.5000E-07,
     &  1.5000E-07,  1.5000E-07,  1.5000E-07,  1.5000E-07,  1.5000E-07,
     &  1.5000E-07,  1.5000E-07/
       DATA (CHI_MLS(7,IP),IP=1,12)/
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    /
       DATA (CHI_MLS(7,IP),IP=13,59)/
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,  0.2090    ,
     &  0.2090    ,  0.2090    /
       DATA (CHI_MLS(8,IP),IP=1,12)/
     &  3.0000E-10,  3.0000E-10,  3.0000E-10,  3.0000E-10,  3.0000E-10,
     &  3.0000E-10,  3.0000E-10,  3.0000E-10,  3.0000E-10,  2.9902E-10,
     &  2.9312E-10,  2.7879E-10/
       DATA (CHI_MLS(8,IP),IP=13,59)/
     &  2.5954E-10,  2.4448E-10,  2.4633E-10,  2.6374E-10,  2.9746E-10,
     &  3.6442E-10,  5.2243E-10,  8.2054E-10,  1.1399E-09,  1.6843E-09,
     &  2.2826E-09,  3.2340E-09,  4.3246E-09,  5.6736E-09,  7.0639E-09,
     &  8.4849E-09,  9.6481E-09,  1.0674E-08,  1.1333E-08,  1.1673E-08,
     &  1.1787E-08,  1.1500E-08,  1.1137E-08,  1.0690E-08,  1.0506E-08,
     &  1.0337E-08,  1.0168E-08,  1.0135E-08,  1.0120E-08,  1.0105E-08,
     &  1.0130E-08,  1.0173E-08,  1.0216E-08,  1.0331E-08,  1.0589E-08,
     &  1.0847E-08,  1.1106E-08,  1.1857E-08,  1.2837E-08,  1.3816E-08,
     &  1.4795E-08,  1.6649E-08,  1.8938E-08,  2.1228E-08,  2.3518E-08,
     &  2.6383E-08,  3.5385E-08/
       DATA (CHI_MLS(9,IP),IP=1,12)/
     &  1.1998E-09,  1.0443E-09,  7.7068E-10,  5.3051E-10,  3.9263E-10,
     &  3.1387E-10,  2.6740E-10,  2.4259E-10,  2.2840E-10,  2.1866E-10,
     &  2.5645E-08,  1.2862E-07/
       DATA (CHI_MLS(9,IP),IP=13,59)/
     &  2.4639E-07,  2.9317E-07,  2.8936E-07,  2.8936E-07,  2.8936E-07,
     &  2.8936E-07,  2.8936E-07,  2.8936E-07,  2.8936E-07,  2.8936E-07,
     &  2.8936E-07,  2.8936E-07,  2.8936E-07,  2.8936E-07,  2.8936E-07,
     &  2.8936E-07,  2.8936E-07,  2.8939E-07,  2.8952E-07,  2.8983E-07,
     &  2.8988E-07,  2.7607E-07,  2.5116E-07,  2.1387E-07,  1.7593E-07,
     &  1.3795E-07,  9.9972E-08,  6.9612E-08,  4.0271E-08,  1.0930E-08,
     &  1.5803E-09,  1.2028E-09,  8.2528E-10,  5.3973E-10,  4.3789E-10,
     &  3.3606E-10,  2.3423E-10,  2.0212E-10,  2.0212E-10,  2.0212E-10,
     &  2.0212E-10,  2.0212E-10,  2.0212E-10,  2.0212E-10,  2.0212E-10,
     &  2.0212E-10,  2.0212E-10/

C ****************** START OF EXECUTABLE CODE ***************************
      HVRSET = '$Revision$'

      STPFAC = 296./1013.

      INDBOUND = TBOUND - 159.
      TBNDFRAC = TBOUND - INT(TBOUND)
      INDLEV0 = TZ(0) - 159.
      T0FRAC = TZ(0) - INT(TZ(0))

      LAYTROP = 0
      LAYSWTCH = 0
      LAYLOW = 0
      DO 7000 LAY = 1, NLAYERS
C        Find the two reference pressures on either side of the
C        layer pressure.  Store them in JP and JP1.  Store in FP the
C        fraction of the difference (in ln(pressure)) between these
C        two values that the layer pressure lies.
         PLOG = ALOG(PAVEL(LAY))
         JP(LAY) = INT(36. - 5*(PLOG+0.04))
         IF (JP(LAY) .LT. 1) THEN
            JP(LAY) = 1
         ELSEIF (JP(LAY) .GT. 58) THEN
            JP(LAY) = 58
         ENDIF
         JP1 = JP(LAY) + 1
         FP = 5. * (PREFLOG(JP(LAY)) - PLOG)

C        Determine, for each reference pressure (JP and JP1), which
C        reference temperature (these are different for each  
C        reference pressure) is nearest the layer temperature but does
C        not exceed it.  Store these indices in JT and JT1, resp.
C        Store in FT (resp. FT1) the fraction of the way between JT
C        (JT1) and the next highest reference temperature that the 
C        layer temperature falls.
         JT(LAY) = INT(3. + (TAVEL(LAY)-TREF(JP(LAY)))/15.)
         IF (JT(LAY) .LT. 1) THEN
            JT(LAY) = 1
         ELSEIF (JT(LAY) .GT. 4) THEN
            JT(LAY) = 4
         ENDIF
         FT = ((TAVEL(LAY)-TREF(JP(LAY)))/15.) - FLOAT(JT(LAY)-3)
         JT1(LAY) = INT(3. + (TAVEL(LAY)-TREF(JP1))/15.)
         IF (JT1(LAY) .LT. 1) THEN
            JT1(LAY) = 1
         ELSEIF (JT1(LAY) .GT. 4) THEN
            JT1(LAY) = 4
         ENDIF
         FT1 = ((TAVEL(LAY)-TREF(JP1))/15.) - FLOAT(JT1(LAY)-3)

         WATER = WKL(1,LAY)/COLDRY(LAY)
         SCALEFAC = PAVEL(LAY) * STPFAC / TAVEL(LAY)

C        If the pressure is less than ~100mb, perform a different
C        set of species interpolations.
         IF (PLOG .LE. 4.56) GO TO 5300
         LAYTROP =  LAYTROP + 1
         IF (PLOG .GE. 6.62) LAYLOW = LAYLOW + 1

C        Set up factors needed to separately include the water vapor
C        foreign-continuum in the calculation of absorption coefficient.
         FORFAC(LAY) = SCALEFAC / (1.+WATER)
         FACTOR = (332.0-TAVEL(LAY))/36.0
         INDFOR(LAY) = MIN(2, MAX(1, INT(FACTOR)))
         FORFRAC(LAY) = FACTOR - FLOAT(INDFOR(LAY))
C
C        Set up factors needed to separately include the water vapor
C        self-continuum in the calculation of absorption coefficient.
         SELFFAC(LAY) = WATER * FORFAC(LAY)
         FACTOR = (TAVEL(LAY)-188.0)/7.2
         INDSELF(LAY) = MIN(9, MAX(1, INT(FACTOR)-7))
         SELFFRAC(LAY) = FACTOR - FLOAT(INDSELF(LAY) + 7)

C        Calculate needed column amounts.
         COLH2O(LAY) = 1.E-20 * WKL(1,LAY)
         COLCO2(LAY) = 1.E-20 * WKL(2,LAY)
         COLO3(LAY) = 1.E-20 * WKL(3,LAY)
c         COLO3(LAY) = 0.
C         COLO3(LAY) = colo3(lay)/1.16
         COLN2O(LAY) = 1.E-20 * WKL(4,LAY)
         COLCH4(LAY) = 1.E-20 * WKL(6,LAY)
         COLO2(LAY) = 1.E-20 * WKL(7,LAY)
         COLSO2(LAY) = 1.E-20 * WKL(9,LAY)
         COLMOL(LAY) = 1.E-20 * COLDRY(LAY) + COLH2O(LAY)
c         colco2(lay) = 0.
c         colo3(lay) = 0.
c         coln2o(lay) = 0.
c         colch4(lay) = 0.
c         colo2(lay) = 0.
c         colmol(lay) = 0.
         IF (COLCO2(LAY) .EQ. 0.) COLCO2(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLN2O(LAY) .EQ. 0.) COLN2O(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLCH4(LAY) .EQ. 0.) COLCH4(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLO2(LAY) .EQ. 0.) COLO2(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLSO2(LAY) .EQ. 0.) COLSO2(LAY) = 1.E-32 * COLDRY(LAY)
C        Using E = 1334.2 cm-1.
         CO2REG = 3.55E-24 * COLDRY(LAY)
         CO2MULT(LAY)= (COLCO2(LAY) - CO2REG) *
     &        272.63*EXP(-1919.4/TAVEL(LAY))/(8.7604E-4*TAVEL(LAY))

         RAT_SO2O3(LAY)=CHI_MLS(9,JP(LAY))/CHI_MLS(3,JP(LAY))
         RAT_SO2O3_1(LAY)=CHI_MLS(9,JP(LAY)+1)/CHI_MLS(3,JP(LAY)+1)

         GO TO 5400

C        Above LAYTROP.
 5300    CONTINUE

C        Set up factors needed to separately include the water vapor
C        foreign-continuum in the calculation of absorption coefficient.
         FORFAC(LAY) = SCALEFAC / (1.+WATER)
         FACTOR = (TAVEL(LAY)-188.0)/36.0
         INDFOR(LAY) = 3
         FORFRAC(LAY) = FACTOR - 1.0
C
C        Calculate needed column amounts.
         COLH2O(LAY) = 1.E-20 * WKL(1,LAY)
         COLCO2(LAY) = 1.E-20 * WKL(2,LAY)
         COLO3(LAY) = 1.E-20 * WKL(3,LAY)
         COLN2O(LAY) = 1.E-20 * WKL(4,LAY)
         COLCH4(LAY) = 1.E-20 * WKL(6,LAY)
         COLO2(LAY) = 1.E-20 * WKL(7,LAY)
         COLSO2(LAY) = 1.E-20 * WKL(9,LAY)
         COLMOL(LAY) = 1.E-20 * COLDRY(LAY) + COLH2O(LAY)

         IF (COLCO2(LAY) .EQ. 0.) COLCO2(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLN2O(LAY) .EQ. 0.) COLN2O(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLCH4(LAY) .EQ. 0.) COLCH4(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLO2(LAY) .EQ. 0.) COLO2(LAY) = 1.E-32 * COLDRY(LAY)
         IF (COLSO2(LAY) .EQ. 0.) COLSO2(LAY) = 1.E-32 * COLDRY(LAY)
         CO2REG = 3.55E-24 * COLDRY(LAY)
         CO2MULT(LAY)= (COLCO2(LAY) - CO2REG) *
     &        272.63*EXP(-1919.4/TAVEL(LAY))/(8.7604E-4*TAVEL(LAY))

         RAT_SO2O3(LAY)=CHI_MLS(9,JP(LAY))/CHI_MLS(3,JP(LAY))
         RAT_SO2O3_1(LAY)=CHI_MLS(9,JP(LAY)+1)/CHI_MLS(3,JP(LAY)+1)

 5400    CONTINUE

C        We have now isolated the layer ln pressure and temperature,
C        between two reference pressures and two reference temperatures 
C        (for each reference pressure).  We multiply the pressure 
C        fraction FP with the appropriate temperature fractions to get 
C        the factors that will be needed for the interpolation that yields
C        the optical depths (performed in routines TAUGBn for band n).

         COMPFP = 1. - FP
         FAC10(LAY) = COMPFP * FT
         FAC00(LAY) = COMPFP * (1. - FT)
         FAC11(LAY) = FP * FT1
         FAC01(LAY) = FP * (1. - FT1)

 7000 CONTINUE

C    Adjust SO2 layer amounts
      
      sumso2 = 0.0
      do lay=nlayers,laytrop,-1
         !print *,lay,colso2(lay)
         sumso2 = sumso2+colso2(lay)/zenith
         colso2(lay) = colso2(lay)*(1.0-0.04*sumso2/1.0e-2)
         if (colso2(lay).le.0) colso2(lay) = 1.E-32 * COLDRY(LAY)
         !print *,colso2(lay)
      end do
         
        
      RETURN
      END

