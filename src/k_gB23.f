C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$
      PARAMETER (MG = 16)
      REAL KA(5,13,MG)
      DIMENSION SELFREF(10,MG),FORREF(3,MG)

      COMMON /HVRSN23/ HVRKG23
      COMMON /K23/ KA ,SELFREF, FORREF

      CHARACTER*15 HVRKG23

      DATA HVRKG23 /'$Revision$'/

C     The array KA contains absorption coefs at the 16 chosen g-values 
C     for a range of pressure levels > ~100mb and temperatures.  The first
C     index in the array, JT, which runs from 1 to 5, corresponds to 
C     different temperatures.  More specifically, JT = 3 means that the 
C     data are for the corresponding TREF for this  pressure level, 
C     JT = 2 refers to the temperatureTREF-15, JT = 1 is for TREF-30, 
C     JT = 4 is for TREF+15, and JT = 5 is for TREF+30.  The second 
C     index, JP, runs from 1 to 13 and refers to the corresponding 
C     pressure level in PREF (e.g. JP = 1 is for a pressure of 1053.63 mb).  
C     The third index, IG, goes from 1 to 16, and tells us which 
C     g-interval the absorption coefficients are for.
      DATA (KA(JT, 1, 1),JT=1,5) /
     &3.3078E-08,3.4034E-08,3.5124E-08,3.4187E-08,3.4744E-08/
      DATA (KA(JT, 2, 1),JT=1,5) /
     &2.5544E-08,2.5873E-08,2.6742E-08,2.7512E-08,2.7504E-08/
      DATA (KA(JT, 3, 1),JT=1,5) /
     &1.8549E-08,1.9611E-08,2.0840E-08,2.2548E-08,2.3069E-08/
      DATA (KA(JT, 4, 1),JT=1,5) /
     &2.8794E-08,3.0837E-08,3.2679E-08,3.4307E-08,3.6901E-08/
      DATA (KA(JT, 5, 1),JT=1,5) /
     &3.6776E-08,3.9144E-08,4.1300E-08,4.3264E-08,4.6626E-08/
      DATA (KA(JT, 6, 1),JT=1,5) /
     &5.9710E-08,6.2941E-08,6.5500E-08,6.7353E-08,6.8774E-08/
      DATA (KA(JT, 7, 1),JT=1,5) /
     &1.2143E-07,1.2932E-07,1.3250E-07,1.3526E-07,1.3849E-07/
      DATA (KA(JT, 8, 1),JT=1,5) /
     &1.2531E-07,1.3241E-07,1.3939E-07,1.4705E-07,1.5465E-07/
      DATA (KA(JT, 9, 1),JT=1,5) /
     &2.0209E-07,2.1134E-07,2.2163E-07,2.3098E-07,2.4004E-07/
      DATA (KA(JT,10, 1),JT=1,5) /
     &1.0750E-06,1.1204E-06,1.1575E-06,1.1923E-06,1.2227E-06/
      DATA (KA(JT,11, 1),JT=1,5) /
     &2.7782E-06,2.8204E-06,2.8406E-06,2.8380E-06,2.8440E-06/
      DATA (KA(JT,12, 1),JT=1,5) /
     &3.8510E-06,3.9934E-06,4.0697E-06,4.1102E-06,4.1571E-06/
      DATA (KA(JT,13, 1),JT=1,5) /
     &4.3157E-06,4.4488E-06,4.5799E-06,4.6585E-06,4.7223E-06/
      DATA (KA(JT, 1, 2),JT=1,5) /
     &8.4637E-07,8.6989E-07,9.0697E-07,9.0000E-07,9.1373E-07/
      DATA (KA(JT, 2, 2),JT=1,5) /
     &6.7062E-07,6.8649E-07,7.2334E-07,7.3645E-07,7.2978E-07/
      DATA (KA(JT, 3, 2),JT=1,5) /
     &5.2317E-07,5.3924E-07,5.5425E-07,5.8347E-07,5.7813E-07/
      DATA (KA(JT, 4, 2),JT=1,5) /
     &3.9868E-07,4.1431E-07,4.2761E-07,4.3892E-07,4.5982E-07/
      DATA (KA(JT, 5, 2),JT=1,5) /
     &3.2074E-07,3.3452E-07,3.4754E-07,3.5582E-07,3.7378E-07/
      DATA (KA(JT, 6, 2),JT=1,5) /
     &4.2465E-07,4.4058E-07,4.5605E-07,4.7192E-07,4.8493E-07/
      DATA (KA(JT, 7, 2),JT=1,5) /
     &4.7581E-07,5.0000E-07,5.2487E-07,5.4192E-07,5.5955E-07/
      DATA (KA(JT, 8, 2),JT=1,5) /
     &1.0592E-06,1.1093E-06,1.1483E-06,1.1923E-06,1.2169E-06/
      DATA (KA(JT, 9, 2),JT=1,5) /
     &5.0835E-06,5.1710E-06,5.2329E-06,5.2644E-06,5.2818E-06/
      DATA (KA(JT,10, 2),JT=1,5) /
     &7.7867E-06,8.2156E-06,8.6002E-06,8.9664E-06,9.2394E-06/
      DATA (KA(JT,11, 2),JT=1,5) /
     &8.9031E-06,9.3573E-06,9.7686E-06,1.0236E-05,1.0591E-05/
      DATA (KA(JT,12, 2),JT=1,5) /
     &9.8068E-06,1.0238E-05,1.0700E-05,1.1108E-05,1.1471E-05/
      DATA (KA(JT,13, 2),JT=1,5) /
     &1.1145E-05,1.1697E-05,1.2123E-05,1.2482E-05,1.2856E-05/
      DATA (KA(JT, 1, 3),JT=1,5) /
     &6.6049E-06,6.7547E-06,7.0104E-06,6.9745E-06,7.0687E-06/
      DATA (KA(JT, 2, 3),JT=1,5) /
     &5.4104E-06,5.5245E-06,5.7623E-06,5.8532E-06,5.8099E-06/
      DATA (KA(JT, 3, 3),JT=1,5) /
     &4.3608E-06,4.4784E-06,4.5778E-06,4.7665E-06,4.7295E-06/
      DATA (KA(JT, 4, 3),JT=1,5) /
     &3.5399E-06,3.6454E-06,3.7356E-06,3.8107E-06,3.9586E-06/
      DATA (KA(JT, 5, 3),JT=1,5) /
     &2.8576E-06,2.9552E-06,3.0325E-06,3.1028E-06,3.2276E-06/
      DATA (KA(JT, 6, 3),JT=1,5) /
     &2.1017E-06,2.1836E-06,2.2554E-06,2.3121E-06,2.3486E-06/
      DATA (KA(JT, 7, 3),JT=1,5) /
     &1.9384E-06,1.9914E-06,2.0428E-06,2.0879E-06,2.1191E-06/
      DATA (KA(JT, 8, 3),JT=1,5) /
     &3.2672E-06,3.3753E-06,3.4694E-06,3.5273E-06,3.5984E-06/
      DATA (KA(JT, 9, 3),JT=1,5) /
     &8.2257E-06,8.4975E-06,8.6969E-06,8.8995E-06,9.0123E-06/
      DATA (KA(JT,10, 3),JT=1,5) /
     &3.5363E-05,3.6204E-05,3.6916E-05,3.7708E-05,3.8294E-05/
      DATA (KA(JT,11, 3),JT=1,5) /
     &4.8837E-05,5.0952E-05,5.2710E-05,5.4347E-05,5.5976E-05/
      DATA (KA(JT,12, 3),JT=1,5) /
     &5.6059E-05,5.8301E-05,6.0100E-05,6.1776E-05,6.3179E-05/
      DATA (KA(JT,13, 3),JT=1,5) /
     &5.7871E-05,5.9505E-05,6.0930E-05,6.2289E-05,6.3278E-05/
      DATA (KA(JT, 1, 4),JT=1,5) /
     &2.7624E-05,2.8017E-05,2.9078E-05,2.8567E-05,2.8724E-05/
      DATA (KA(JT, 2, 4),JT=1,5) /
     &2.3107E-05,2.3464E-05,2.4357E-05,2.4519E-05,2.4024E-05/
      DATA (KA(JT, 3, 4),JT=1,5) /
     &1.9113E-05,1.9400E-05,1.9615E-05,2.0420E-05,1.9962E-05/
      DATA (KA(JT, 4, 4),JT=1,5) /
     &1.5873E-05,1.6138E-05,1.6349E-05,1.6528E-05,1.7176E-05/
      DATA (KA(JT, 5, 4),JT=1,5) /
     &1.3198E-05,1.3437E-05,1.3641E-05,1.3810E-05,1.4393E-05/
      DATA (KA(JT, 6, 4),JT=1,5) /
     &1.0951E-05,1.1172E-05,1.1352E-05,1.1506E-05,1.1631E-05/
      DATA (KA(JT, 7, 4),JT=1,5) /
     &8.6121E-06,8.8300E-06,9.0149E-06,9.1565E-06,9.2594E-06/
      DATA (KA(JT, 8, 4),JT=1,5) /
     &7.1478E-06,7.2918E-06,7.4035E-06,7.4959E-06,7.5566E-06/
      DATA (KA(JT, 9, 4),JT=1,5) /
     &1.6458E-05,1.7092E-05,1.7686E-05,1.7967E-05,1.8273E-05/
      DATA (KA(JT,10, 4),JT=1,5) /
     &4.7953E-05,4.9663E-05,5.1524E-05,5.2800E-05,5.4192E-05/
      DATA (KA(JT,11, 4),JT=1,5) /
     &9.4263E-05,9.5557E-05,9.6513E-05,9.7430E-05,9.7733E-05/
      DATA (KA(JT,12, 4),JT=1,5) /
     &1.2087E-04,1.2152E-04,1.2240E-04,1.2318E-04,1.2403E-04/
      DATA (KA(JT,13, 4),JT=1,5) /
     &1.2781E-04,1.2897E-04,1.3049E-04,1.3171E-04,1.3337E-04/
      DATA (KA(JT, 1, 5),JT=1,5) /
     &8.2859E-05,8.4817E-05,8.9056E-05,8.8057E-05,8.9410E-05/
      DATA (KA(JT, 2, 5),JT=1,5) /
     &7.0937E-05,7.2685E-05,7.6796E-05,7.7955E-05,7.6735E-05/
      DATA (KA(JT, 3, 5),JT=1,5) /
     &5.9876E-05,6.1448E-05,6.2837E-05,6.6293E-05,6.4996E-05/
      DATA (KA(JT, 4, 5),JT=1,5) /
     &5.0598E-05,5.2054E-05,5.3293E-05,5.4332E-05,5.7333E-05/
      DATA (KA(JT, 5, 5),JT=1,5) /
     &4.2742E-05,4.4035E-05,4.5164E-05,4.6134E-05,4.8935E-05/
      DATA (KA(JT, 6, 5),JT=1,5) /
     &3.5769E-05,3.6975E-05,3.8038E-05,3.8917E-05,3.9681E-05/
      DATA (KA(JT, 7, 5),JT=1,5) /
     &2.9747E-05,3.0824E-05,3.1756E-05,3.2589E-05,3.3314E-05/
      DATA (KA(JT, 8, 5),JT=1,5) /
     &2.1994E-05,2.2945E-05,2.3786E-05,2.4517E-05,2.5155E-05/
      DATA (KA(JT, 9, 5),JT=1,5) /
     &2.2298E-05,2.2688E-05,2.3127E-05,2.3803E-05,2.4335E-05/
      DATA (KA(JT,10, 5),JT=1,5) /
     &8.8898E-05,9.1280E-05,9.3333E-05,9.5224E-05,9.6298E-05/
      DATA (KA(JT,11, 5),JT=1,5) /
     &1.2299E-04,1.2407E-04,1.2536E-04,1.2642E-04,1.2813E-04/
      DATA (KA(JT,12, 5),JT=1,5) /
     &1.4539E-04,1.4851E-04,1.5022E-04,1.5157E-04,1.5204E-04/
      DATA (KA(JT,13, 5),JT=1,5) /
     &1.5949E-04,1.6239E-04,1.6467E-04,1.6667E-04,1.6801E-04/
      DATA (KA(JT, 1, 6),JT=1,5) /
     &2.5339E-04,2.5995E-04,2.7170E-04,2.6963E-04,2.7413E-04/
      DATA (KA(JT, 2, 6),JT=1,5) /
     &2.1908E-04,2.2404E-04,2.3522E-04,2.3902E-04,2.3717E-04/
      DATA (KA(JT, 3, 6),JT=1,5) /
     &1.8611E-04,1.9076E-04,1.9518E-04,2.0503E-04,2.0222E-04/
      DATA (KA(JT, 4, 6),JT=1,5) /
     &1.5769E-04,1.6210E-04,1.6597E-04,1.6948E-04,1.7817E-04/
      DATA (KA(JT, 5, 6),JT=1,5) /
     &1.3402E-04,1.3792E-04,1.4153E-04,1.4434E-04,1.5267E-04/
      DATA (KA(JT, 6, 6),JT=1,5) /
     &1.1390E-04,1.1743E-04,1.2065E-04,1.2316E-04,1.2543E-04/
      DATA (KA(JT, 7, 6),JT=1,5) /
     &9.6417E-05,9.9612E-05,1.0233E-04,1.0453E-04,1.0652E-04/
      DATA (KA(JT, 8, 6),JT=1,5) /
     &8.1395E-05,8.4205E-05,8.6346E-05,8.8406E-05,9.0122E-05/
      DATA (KA(JT, 9, 6),JT=1,5) /
     &4.7776E-05,4.8971E-05,4.9736E-05,4.9917E-05,5.0289E-05/
      DATA (KA(JT,10, 6),JT=1,5) /
     &1.0698E-04,1.0815E-04,1.0817E-04,1.0799E-04,1.0851E-04/
      DATA (KA(JT,11, 6),JT=1,5) /
     &2.0220E-04,2.0727E-04,2.1241E-04,2.1675E-04,2.1989E-04/
      DATA (KA(JT,12, 6),JT=1,5) /
     &2.3474E-04,2.3601E-04,2.3974E-04,2.4383E-04,2.4876E-04/
      DATA (KA(JT,13, 6),JT=1,5) /
     &2.3410E-04,2.3809E-04,2.4185E-04,2.4554E-04,2.4952E-04/
      DATA (KA(JT, 1, 7),JT=1,5) /
     &6.7024E-04,6.8026E-04,7.0419E-04,7.0159E-04,7.1089E-04/
      DATA (KA(JT, 2, 7),JT=1,5) /
     &5.8729E-04,5.9778E-04,6.2097E-04,6.2912E-04,6.2423E-04/
      DATA (KA(JT, 3, 7),JT=1,5) /
     &5.0967E-04,5.1900E-04,5.2765E-04,5.4794E-04,5.4266E-04/
      DATA (KA(JT, 4, 7),JT=1,5) /
     &4.4167E-04,4.5006E-04,4.5793E-04,4.6469E-04,4.8335E-04/
      DATA (KA(JT, 5, 7),JT=1,5) /
     &3.8096E-04,3.8881E-04,3.9576E-04,4.0259E-04,4.2086E-04/
      DATA (KA(JT, 6, 7),JT=1,5) /
     &3.2818E-04,3.3539E-04,3.4192E-04,3.4829E-04,3.5405E-04/
      DATA (KA(JT, 7, 7),JT=1,5) /
     &2.8259E-04,2.8946E-04,2.9584E-04,3.0203E-04,3.0742E-04/
      DATA (KA(JT, 8, 7),JT=1,5) /
     &2.4273E-04,2.4912E-04,2.5546E-04,2.6110E-04,2.6607E-04/
      DATA (KA(JT, 9, 7),JT=1,5) /
     &1.9937E-04,2.0653E-04,2.1314E-04,2.1968E-04,2.2520E-04/
      DATA (KA(JT,10, 7),JT=1,5) /
     &1.3306E-04,1.3331E-04,1.3393E-04,1.3538E-04,1.3616E-04/
      DATA (KA(JT,11, 7),JT=1,5) /
     &1.6236E-04,1.6154E-04,1.6187E-04,1.6113E-04,1.6209E-04/
      DATA (KA(JT,12, 7),JT=1,5) /
     &1.7872E-04,1.8355E-04,1.8612E-04,1.8792E-04,1.8745E-04/
      DATA (KA(JT,13, 7),JT=1,5) /
     &1.8970E-04,1.9384E-04,1.9773E-04,2.0261E-04,2.0377E-04/
      DATA (KA(JT, 1, 8),JT=1,5) /
     &1.8130E-03,1.8305E-03,1.8716E-03,1.8655E-03,1.8814E-03/
      DATA (KA(JT, 2, 8),JT=1,5) /
     &1.6420E-03,1.6600E-03,1.7006E-03,1.7156E-03,1.7108E-03/
      DATA (KA(JT, 3, 8),JT=1,5) /
     &1.4687E-03,1.4870E-03,1.5042E-03,1.5423E-03,1.5376E-03/
      DATA (KA(JT, 4, 8),JT=1,5) /
     &1.3068E-03,1.3248E-03,1.3421E-03,1.3592E-03,1.3960E-03/
      DATA (KA(JT, 5, 8),JT=1,5) /
     &1.1574E-03,1.1753E-03,1.1928E-03,1.2097E-03,1.2467E-03/
      DATA (KA(JT, 6, 8),JT=1,5) /
     &1.0167E-03,1.0342E-03,1.0515E-03,1.0681E-03,1.0840E-03/
      DATA (KA(JT, 7, 8),JT=1,5) /
     &8.8992E-04,9.0662E-04,9.2299E-04,9.3844E-04,9.5347E-04/
      DATA (KA(JT, 8, 8),JT=1,5) /
     &7.8445E-04,8.0031E-04,8.1639E-04,8.3131E-04,8.4631E-04/
      DATA (KA(JT, 9, 8),JT=1,5) /
     &6.9812E-04,7.1312E-04,7.2801E-04,7.4222E-04,7.5605E-04/
      DATA (KA(JT,10, 8),JT=1,5) /
     &3.2521E-04,3.3835E-04,3.5194E-04,3.6381E-04,3.7594E-04/
      DATA (KA(JT,11, 8),JT=1,5) /
     &3.1406E-04,3.2013E-04,3.2456E-04,3.3261E-04,3.3741E-04/
      DATA (KA(JT,12, 8),JT=1,5) /
     &2.8132E-04,2.8732E-04,2.9674E-04,3.0509E-04,3.1393E-04/
      DATA (KA(JT,13, 8),JT=1,5) /
     &2.5704E-04,2.6316E-04,2.7195E-04,2.7732E-04,2.8905E-04/
      DATA (KA(JT, 1, 9),JT=1,5) /
     &6.7370E-03,6.7873E-03,6.8896E-03,6.8819E-03,6.9268E-03/
      DATA (KA(JT, 2, 9),JT=1,5) /
     &6.3111E-03,6.3622E-03,6.4623E-03,6.5060E-03,6.5027E-03/
      DATA (KA(JT, 3, 9),JT=1,5) /
     &5.8834E-03,5.9361E-03,5.9858E-03,6.0811E-03,6.0806E-03/
      DATA (KA(JT, 4, 9),JT=1,5) /
     &5.4753E-03,5.5309E-03,5.5823E-03,5.6306E-03,5.7229E-03/
      DATA (KA(JT, 5, 9),JT=1,5) /
     &5.0781E-03,5.1373E-03,5.1895E-03,5.2391E-03,5.3328E-03/
      DATA (KA(JT, 6, 9),JT=1,5) /
     &4.6791E-03,4.7408E-03,4.7949E-03,4.8470E-03,4.8993E-03/
      DATA (KA(JT, 7, 9),JT=1,5) /
     &4.2724E-03,4.3381E-03,4.3951E-03,4.4503E-03,4.5061E-03/
      DATA (KA(JT, 8, 9),JT=1,5) /
     &3.8568E-03,3.9254E-03,3.9848E-03,4.0431E-03,4.1005E-03/
      DATA (KA(JT, 9, 9),JT=1,5) /
     &3.5657E-03,3.6362E-03,3.6994E-03,3.7608E-03,3.8215E-03/
      DATA (KA(JT,10, 9),JT=1,5) /
     &3.3774E-03,3.4501E-03,3.5171E-03,3.5811E-03,3.6450E-03/
      DATA (KA(JT,11, 9),JT=1,5) /
     &2.3923E-03,2.4623E-03,2.5263E-03,2.5890E-03,2.6517E-03/
      DATA (KA(JT,12, 9),JT=1,5) /
     &1.6959E-03,1.7542E-03,1.8077E-03,1.8614E-03,1.9187E-03/
      DATA (KA(JT,13, 9),JT=1,5) /
     &1.1732E-03,1.2232E-03,1.2712E-03,1.3221E-03,1.3677E-03/
      DATA (KA(JT, 1,10),JT=1,5) /
     &1.9604E-02,1.9698E-02,1.9938E-02,1.9854E-02,1.9950E-02/
      DATA (KA(JT, 2,10),JT=1,5) /
     &1.8714E-02,1.8803E-02,1.9035E-02,1.9131E-02,1.9102E-02/
      DATA (KA(JT, 3,10),JT=1,5) /
     &1.7676E-02,1.7785E-02,1.7904E-02,1.8189E-02,1.8156E-02/
      DATA (KA(JT, 4,10),JT=1,5) /
     &1.6662E-02,1.6773E-02,1.6908E-02,1.7056E-02,1.7303E-02/
      DATA (KA(JT, 5,10),JT=1,5) /
     &1.5655E-02,1.5775E-02,1.5942E-02,1.6103E-02,1.6359E-02/
      DATA (KA(JT, 6,10),JT=1,5) /
     &1.4694E-02,1.4838E-02,1.5029E-02,1.5200E-02,1.5320E-02/
      DATA (KA(JT, 7,10),JT=1,5) /
     &1.3797E-02,1.3959E-02,1.4174E-02,1.4350E-02,1.4471E-02/
      DATA (KA(JT, 8,10),JT=1,5) /
     &1.2902E-02,1.3089E-02,1.3313E-02,1.3489E-02,1.3626E-02/
      DATA (KA(JT, 9,10),JT=1,5) /
     &1.1897E-02,1.2105E-02,1.2333E-02,1.2512E-02,1.2656E-02/
      DATA (KA(JT,10,10),JT=1,5) /
     &1.1834E-02,1.2072E-02,1.2331E-02,1.2548E-02,1.2721E-02/
      DATA (KA(JT,11,10),JT=1,5) /
     &1.1416E-02,1.1699E-02,1.1968E-02,1.2141E-02,1.2311E-02/
      DATA (KA(JT,12,10),JT=1,5) /
     &1.0776E-02,1.1070E-02,1.1309E-02,1.1507E-02,1.1702E-02/
      DATA (KA(JT,13,10),JT=1,5) /
     &9.9577E-03,1.0263E-02,1.0492E-02,1.0719E-02,1.0951E-02/
      DATA (KA(JT, 1,11),JT=1,5) /
     &2.9783E-02,2.9883E-02,3.0210E-02,3.0141E-02,3.0248E-02/
      DATA (KA(JT, 2,11),JT=1,5) /
     &2.8562E-02,2.8743E-02,2.9186E-02,2.9339E-02,2.9258E-02/
      DATA (KA(JT, 3,11),JT=1,5) /
     &2.7212E-02,2.7429E-02,2.7654E-02,2.8083E-02,2.8023E-02/
      DATA (KA(JT, 4,11),JT=1,5) /
     &2.5949E-02,2.6197E-02,2.6424E-02,2.6623E-02,2.7069E-02/
      DATA (KA(JT, 5,11),JT=1,5) /
     &2.4686E-02,2.4942E-02,2.5176E-02,2.5403E-02,2.5908E-02/
      DATA (KA(JT, 6,11),JT=1,5) /
     &2.3430E-02,2.3686E-02,2.3923E-02,2.4158E-02,2.4433E-02/
      DATA (KA(JT, 7,11),JT=1,5) /
     &2.2171E-02,2.2424E-02,2.2653E-02,2.2909E-02,2.3203E-02/
      DATA (KA(JT, 8,11),JT=1,5) /
     &2.0928E-02,2.1171E-02,2.1407E-02,2.1699E-02,2.2018E-02/
      DATA (KA(JT, 9,11),JT=1,5) /
     &1.9076E-02,1.9320E-02,1.9548E-02,1.9858E-02,2.0150E-02/
      DATA (KA(JT,10,11),JT=1,5) /
     &1.9537E-02,1.9788E-02,2.0064E-02,2.0391E-02,2.0843E-02/
      DATA (KA(JT,11,11),JT=1,5) /
     &1.9137E-02,1.9444E-02,1.9793E-02,2.0268E-02,2.0695E-02/
      DATA (KA(JT,12,11),JT=1,5) /
     &1.8393E-02,1.8715E-02,1.9143E-02,1.9652E-02,2.0037E-02/
      DATA (KA(JT,13,11),JT=1,5) /
     &1.7255E-02,1.7680E-02,1.8170E-02,1.8641E-02,1.9010E-02/
      DATA (KA(JT, 1,12),JT=1,5) /
     &4.6641E-02,4.6796E-02,4.7107E-02,4.6977E-02,4.7093E-02/
      DATA (KA(JT, 2,12),JT=1,5) /
     &4.6819E-02,4.6956E-02,4.7337E-02,4.7426E-02,4.7339E-02/
      DATA (KA(JT, 3,12),JT=1,5) /
     &4.6276E-02,4.6366E-02,4.6462E-02,4.6935E-02,4.6862E-02/
      DATA (KA(JT, 4,12),JT=1,5) /
     &4.4986E-02,4.5103E-02,4.5287E-02,4.5535E-02,4.6152E-02/
      DATA (KA(JT, 5,12),JT=1,5) /
     &4.3367E-02,4.3552E-02,4.3761E-02,4.4035E-02,4.4674E-02/
      DATA (KA(JT, 6,12),JT=1,5) /
     &4.1584E-02,4.1793E-02,4.2020E-02,4.2344E-02,4.2701E-02/
      DATA (KA(JT, 7,12),JT=1,5) /
     &3.9785E-02,4.0007E-02,4.0269E-02,4.0635E-02,4.1003E-02/
      DATA (KA(JT, 8,12),JT=1,5) /
     &3.7918E-02,3.8155E-02,3.8468E-02,3.8839E-02,3.9210E-02/
      DATA (KA(JT, 9,12),JT=1,5) /
     &3.5060E-02,3.5328E-02,3.5698E-02,3.6105E-02,3.6551E-02/
      DATA (KA(JT,10,12),JT=1,5) /
     &3.5216E-02,3.5547E-02,3.5872E-02,3.6327E-02,3.6542E-02/
      DATA (KA(JT,11,12),JT=1,5) /
     &3.5158E-02,3.5374E-02,3.5869E-02,3.6242E-02,3.6691E-02/
      DATA (KA(JT,12,12),JT=1,5) /
     &3.4264E-02,3.4598E-02,3.5106E-02,3.5505E-02,3.6021E-02/
      DATA (KA(JT,13,12),JT=1,5) /
     &3.2716E-02,3.3195E-02,3.3830E-02,3.4228E-02,3.4925E-02/
      DATA (KA(JT, 1,13),JT=1,5) /
     &7.6084E-02,7.6052E-02,7.6051E-02,7.5851E-02,7.5753E-02/
      DATA (KA(JT, 2,13),JT=1,5) /
     &7.9580E-02,7.9564E-02,7.9664E-02,7.9619E-02,7.9414E-02/
      DATA (KA(JT, 3,13),JT=1,5) /
     &8.2218E-02,8.2302E-02,8.2367E-02,8.2543E-02,8.2246E-02/
      DATA (KA(JT, 4,13),JT=1,5) /
     &8.3613E-02,8.3740E-02,8.3824E-02,8.3866E-02,8.4103E-02/
      DATA (KA(JT, 5,13),JT=1,5) /
     &8.3913E-02,8.4123E-02,8.4289E-02,8.4403E-02,8.4840E-02/
      DATA (KA(JT, 6,13),JT=1,5) /
     &8.3159E-02,8.3442E-02,8.3718E-02,8.3891E-02,8.3993E-02/
      DATA (KA(JT, 7,13),JT=1,5) /
     &8.1401E-02,8.1826E-02,8.2202E-02,8.2451E-02,8.2656E-02/
      DATA (KA(JT, 8,13),JT=1,5) /
     &7.8949E-02,7.9505E-02,7.9978E-02,8.0351E-02,8.0690E-02/
      DATA (KA(JT, 9,13),JT=1,5) /
     &7.6002E-02,7.6671E-02,7.7279E-02,7.7752E-02,7.8257E-02/
      DATA (KA(JT,10,13),JT=1,5) /
     &6.9777E-02,7.0425E-02,7.1054E-02,7.1706E-02,7.2167E-02/
      DATA (KA(JT,11,13),JT=1,5) /
     &7.2929E-02,7.3732E-02,7.4323E-02,7.5246E-02,7.5786E-02/
      DATA (KA(JT,12,13),JT=1,5) /
     &7.3007E-02,7.4429E-02,7.4922E-02,7.5916E-02,7.6947E-02/
      DATA (KA(JT,13,13),JT=1,5) /
     &7.1376E-02,7.2507E-02,7.3710E-02,7.4716E-02,7.5702E-02/
      DATA (KA(JT, 1,14),JT=1,5) /
     &1.2585E-01,1.2569E-01,1.2576E-01,1.2571E-01,1.2570E-01/
      DATA (KA(JT, 2,14),JT=1,5) /
     &1.3868E-01,1.3853E-01,1.3849E-01,1.3824E-01,1.3803E-01/
      DATA (KA(JT, 3,14),JT=1,5) /
     &1.5142E-01,1.5135E-01,1.5112E-01,1.5093E-01,1.5074E-01/
      DATA (KA(JT, 4,14),JT=1,5) /
     &1.6359E-01,1.6359E-01,1.6339E-01,1.6311E-01,1.6285E-01/
      DATA (KA(JT, 5,14),JT=1,5) /
     &1.7462E-01,1.7472E-01,1.7456E-01,1.7421E-01,1.7393E-01/
      DATA (KA(JT, 6,14),JT=1,5) /
     &1.8403E-01,1.8447E-01,1.8450E-01,1.8424E-01,1.8395E-01/
      DATA (KA(JT, 7,14),JT=1,5) /
     &1.9179E-01,1.9239E-01,1.9256E-01,1.9253E-01,1.9252E-01/
      DATA (KA(JT, 8,14),JT=1,5) /
     &1.9772E-01,1.9863E-01,1.9901E-01,1.9923E-01,1.9947E-01/
      DATA (KA(JT, 9,14),JT=1,5) /
     &2.0154E-01,2.0279E-01,2.0355E-01,2.0417E-01,2.0473E-01/
      DATA (KA(JT,10,14),JT=1,5) /
     &1.8853E-01,1.9028E-01,1.9160E-01,1.9253E-01,1.9393E-01/
      DATA (KA(JT,11,14),JT=1,5) /
     &1.8013E-01,1.8167E-01,1.8320E-01,1.8375E-01,1.8507E-01/
      DATA (KA(JT,12,14),JT=1,5) /
     &1.9011E-01,1.9027E-01,1.9283E-01,1.9402E-01,1.9478E-01/
      DATA (KA(JT,13,14),JT=1,5) /
     &1.9594E-01,1.9738E-01,1.9911E-01,2.0124E-01,2.0282E-01/
      DATA (KA(JT, 1,15),JT=1,5) /
     &2.2369E-01,2.2259E-01,2.2155E-01,2.2059E-01,2.1997E-01/
      DATA (KA(JT, 2,15),JT=1,5) /
     &2.5602E-01,2.5478E-01,2.5377E-01,2.5306E-01,2.5237E-01/
      DATA (KA(JT, 3,15),JT=1,5) /
     &2.9258E-01,2.9107E-01,2.8998E-01,2.8920E-01,2.8830E-01/
      DATA (KA(JT, 4,15),JT=1,5) /
     &3.3067E-01,3.2888E-01,3.2753E-01,3.2646E-01,3.2566E-01/
      DATA (KA(JT, 5,15),JT=1,5) /
     &3.7114E-01,3.6880E-01,3.6713E-01,3.6598E-01,3.6499E-01/
      DATA (KA(JT, 6,15),JT=1,5) /
     &4.1494E-01,4.1167E-01,4.0935E-01,4.0779E-01,4.0636E-01/
      DATA (KA(JT, 7,15),JT=1,5) /
     &4.6115E-01,4.5729E-01,4.5455E-01,4.5230E-01,4.5004E-01/
      DATA (KA(JT, 8,15),JT=1,5) /
     &5.0906E-01,5.0463E-01,5.0137E-01,4.9843E-01,4.9550E-01/
      DATA (KA(JT, 9,15),JT=1,5) /
     &5.5829E-01,5.5330E-01,5.4936E-01,5.4557E-01,5.4168E-01/
      DATA (KA(JT,10,15),JT=1,5) /
     &6.0814E-01,6.0274E-01,5.9806E-01,5.9326E-01,5.8833E-01/
      DATA (KA(JT,11,15),JT=1,5) /
     &6.2954E-01,6.2588E-01,6.2076E-01,6.1665E-01,6.1184E-01/
      DATA (KA(JT,12,15),JT=1,5) /
     &6.2585E-01,6.2437E-01,6.1807E-01,6.1303E-01,6.0869E-01/
      DATA (KA(JT,13,15),JT=1,5) /
     &6.4856E-01,6.4505E-01,6.3861E-01,6.3277E-01,6.2702E-01/
      DATA (KA(JT, 1,16),JT=1,5) /
     &3.3327E-01,3.3385E-01,3.3538E-01,3.3638E-01,3.3736E-01/
      DATA (KA(JT, 2,16),JT=1,5) /
     &4.0916E-01,4.0842E-01,4.0848E-01,4.0865E-01,4.0854E-01/
      DATA (KA(JT, 3,16),JT=1,5) /
     &5.0099E-01,4.9888E-01,4.9727E-01,4.9588E-01,4.9422E-01/
      DATA (KA(JT, 4,16),JT=1,5) /
     &6.0389E-01,6.0029E-01,5.9704E-01,5.9367E-01,5.9071E-01/
      DATA (KA(JT, 5,16),JT=1,5) /
     &7.1868E-01,7.1337E-01,7.0835E-01,7.0318E-01,6.9852E-01/
      DATA (KA(JT, 6,16),JT=1,5) /
     &8.4815E-01,8.4138E-01,8.3446E-01,8.2728E-01,8.2023E-01/
      DATA (KA(JT, 7,16),JT=1,5) /
     &9.9512E-01,9.8644E-01,9.7698E-01,9.6695E-01,9.5712E-01/
      DATA (KA(JT, 8,16),JT=1,5) /
     &1.1606E+00,1.1485E+00,1.1354E+00,1.1218E+00,1.1077E+00/
      DATA (KA(JT, 9,16),JT=1,5) /
     &1.3444E+00,1.3282E+00,1.3102E+00,1.2917E+00,1.2735E+00/
      DATA (KA(JT,10,16),JT=1,5) /
     &1.5423E+00,1.5207E+00,1.4970E+00,1.4733E+00,1.4494E+00/
      DATA (KA(JT,11,16),JT=1,5) /
     &1.7462E+00,1.7138E+00,1.6827E+00,1.6518E+00,1.6199E+00/
      DATA (KA(JT,12,16),JT=1,5) /
     &1.9577E+00,1.9145E+00,1.8728E+00,1.8320E+00,1.7910E+00/
      DATA (KA(JT,13,16),JT=1,5) /
     &2.1716E+00,2.1171E+00,2.0639E+00,2.0113E+00,1.9587E+00/

      DATA FORREF/ 
     &     3.15770E-08, 6.71978E-08, 4.40649E-07,
     &     3.13674E-07, 2.85252E-07, 4.21024E-06,
     &     1.35818E-06, 1.45071E-06, 6.11285E-06,
     &     5.34065E-06, 5.86268E-06, 9.33970E-06,
     &     9.64007E-06, 1.07110E-05, 1.04486E-05,
     &     3.02775E-05, 3.57530E-05, 3.40724E-05,
     &     1.02437E-04, 1.08475E-04, 1.05245E-04,
     &     1.46054E-04, 1.41490E-04, 1.33071E-04,
     &     1.63978E-04, 1.50208E-04, 1.42864E-04,
     &     2.20412E-04, 1.82943E-04, 1.50941E-04,
     &     2.28877E-04, 1.97679E-04, 1.63220E-04,
     &     2.34177E-04, 2.17734E-04, 1.85038E-04,
     &     2.57187E-04, 2.41570E-04, 2.21178E-04,
     &     2.72455E-04, 2.70637E-04, 2.56269E-04,
     &     3.39445E-04, 3.00268E-04, 2.86574E-04,
     &     3.38841E-04, 3.55428E-04, 3.53794E-04/

C     The array SELFREF contains the coefficient of the water vapor
C     self-continuum (including the energy term).  The first index
C     refers to temperature in 7.2 degree increments.  For instance,
C     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
C     etc.  The second index runs over the g-channel (1 to 16).
      DATA (SELFREF(JT, 1),JT=1,10)  /
     & 1.00945E-05, 8.01113E-06, 6.35771E-06, 5.04554E-06, 4.00419E-06,
     & 3.17777E-06, 2.52191E-06, 2.00141E-06, 1.58834E-06, 1.26052E-06/
      DATA (SELFREF(JT, 2),JT=1,10)  /
     & 1.07573E-05, 9.99809E-06, 9.29245E-06, 8.63661E-06, 8.02706E-06,
     & 7.46053E-06, 6.93399E-06, 6.44460E-06, 5.98976E-06, 5.56702E-06/
      DATA (SELFREF(JT, 3),JT=1,10)  /
     & 3.50389E-05, 3.19234E-05, 2.90850E-05, 2.64989E-05, 2.41428E-05,
     & 2.19962E-05, 2.00404E-05, 1.82586E-05, 1.66351E-05, 1.51560E-05/
      DATA (SELFREF(JT, 4),JT=1,10)  /
     & 1.22993E-04, 1.10885E-04, 9.99691E-05, 9.01277E-05, 8.12551E-05,
     & 7.32559E-05, 6.60443E-05, 5.95426E-05, 5.36809E-05, 4.83963E-05/
      DATA (SELFREF(JT, 5),JT=1,10)  /
     & 2.06434E-04, 1.87435E-04, 1.70185E-04, 1.54522E-04, 1.40301E-04,
     & 1.27388E-04, 1.15664E-04, 1.05019E-04, 9.53540E-05, 8.65783E-05/
      DATA (SELFREF(JT, 6),JT=1,10)  /
     & 5.90645E-04, 5.33109E-04, 4.81177E-04, 4.34305E-04, 3.91998E-04,
     & 3.53812E-04, 3.19346E-04, 2.88238E-04, 2.60160E-04, 2.34817E-04/
      DATA (SELFREF(JT, 7),JT=1,10)  /
     & 1.63029E-03, 1.48773E-03, 1.35763E-03, 1.23891E-03, 1.13057E-03,
     & 1.03170E-03, 9.41483E-04, 8.59153E-04, 7.84023E-04, 7.15462E-04/
      DATA (SELFREF(JT, 8),JT=1,10)  /
     & 2.04528E-03, 1.89258E-03, 1.75128E-03, 1.62053E-03, 1.49954E-03,
     & 1.38758E-03, 1.28398E-03, 1.18812E-03, 1.09941E-03, 1.01733E-03/
      DATA (SELFREF(JT, 9),JT=1,10)  /
     & 2.10589E-03, 1.97078E-03, 1.84434E-03, 1.72601E-03, 1.61528E-03,
     & 1.51164E-03, 1.41466E-03, 1.32390E-03, 1.23896E-03, 1.15947E-03/
      DATA (SELFREF(JT,10),JT=1,10)  /
     & 2.45098E-03, 2.33745E-03, 2.22918E-03, 2.12592E-03, 2.02745E-03,
     & 1.93353E-03, 1.84397E-03, 1.75856E-03, 1.67710E-03, 1.59941E-03/
      DATA (SELFREF(JT,11),JT=1,10)  /
     & 2.67460E-03, 2.53325E-03, 2.39936E-03, 2.27255E-03, 2.15244E-03,
     & 2.03868E-03, 1.93093E-03, 1.82888E-03, 1.73222E-03, 1.64067E-03/
      DATA (SELFREF(JT,12),JT=1,10)  /
     & 3.04510E-03, 2.83919E-03, 2.64720E-03, 2.46820E-03, 2.30130E-03,
     & 2.14568E-03, 2.00059E-03, 1.86531E-03, 1.73918E-03, 1.62157E-03/
      DATA (SELFREF(JT,13),JT=1,10)  /
     & 3.38445E-03, 3.14719E-03, 2.92655E-03, 2.72139E-03, 2.53060E-03,
     & 2.35319E-03, 2.18822E-03, 2.03482E-03, 1.89217E-03, 1.75952E-03/
      DATA (SELFREF(JT,14),JT=1,10)  /
     & 3.88649E-03, 3.57018E-03, 3.27961E-03, 3.01269E-03, 2.76750E-03,
     & 2.54226E-03, 2.33535E-03, 2.14528E-03, 1.97068E-03, 1.81029E-03/
      DATA (SELFREF(JT,15),JT=1,10)  /
     & 4.12547E-03, 3.87413E-03, 3.63810E-03, 3.41646E-03, 3.20831E-03,
     & 3.01285E-03, 2.82930E-03, 2.65693E-03, 2.49506E-03, 2.34305E-03/
      DATA (SELFREF(JT,16),JT=1,10)  /
     & 5.34327E-03, 4.82967E-03, 4.36544E-03, 3.94583E-03, 3.56655E-03,
     & 3.22373E-03, 2.91387E-03, 2.63378E-03, 2.38062E-03, 2.15179E-03/
