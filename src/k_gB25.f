      PARAMETER (MG = 16)
      REAL KA(5,13,MG)

      COMMON /HVRSNB/ HVRKG(16:15+NBANDS)
      COMMON /K25/ KA 

      CHARACTER*8 HVRKG

      DATA HVRKG(25)  / '%I%' /

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
     &1.6461E-09,1.6782E-09,1.9339E-09,1.7100E-09,1.7045E-09/
      DATA (KA(JT, 2, 1),JT=1,5) /
     &2.8759E-09,2.9469E-09,3.3789E-09,3.4357E-09,2.8833E-09/
      DATA (KA(JT, 3, 1),JT=1,5) /
     &5.5148E-09,5.4808E-09,5.4190E-09,6.8260E-09,5.1972E-09/
      DATA (KA(JT, 4, 1),JT=1,5) /
     &9.5336E-09,9.4552E-09,9.3001E-09,9.0961E-09,1.4451E-08/
      DATA (KA(JT, 5, 1),JT=1,5) /
     &1.4930E-08,1.4736E-08,1.4432E-08,1.4074E-08,2.4102E-08/
      DATA (KA(JT, 6, 1),JT=1,5) /
     &2.2770E-08,2.2301E-08,2.1778E-08,2.1194E-08,2.0569E-08/
      DATA (KA(JT, 7, 1),JT=1,5) /
     &3.4699E-08,3.3951E-08,3.3124E-08,3.2144E-08,3.1220E-08/
      DATA (KA(JT, 8, 1),JT=1,5) /
     &6.2339E-08,6.0405E-08,5.9548E-08,5.8214E-08,5.6977E-08/
      DATA (KA(JT, 9, 1),JT=1,5) /
     &1.7411E-07,1.7654E-07,1.8315E-07,1.8100E-07,1.7839E-07/
      DATA (KA(JT,10, 1),JT=1,5) /
     &2.3526E-07,2.2729E-07,2.1947E-07,2.1188E-07,2.0454E-07/
      DATA (KA(JT,11, 1),JT=1,5) /
     &2.3535E-07,2.2737E-07,2.1956E-07,2.1196E-07,2.0461E-07/
      DATA (KA(JT,12, 1),JT=1,5) /
     &2.3539E-07,2.2740E-07,2.1959E-07,2.1199E-07,2.0465E-07/
      DATA (KA(JT,13, 1),JT=1,5) /
     &2.3543E-07,2.2744E-07,2.1962E-07,2.1202E-07,2.0467E-07/
      DATA (KA(JT, 1, 2),JT=1,5) /
     &6.2912E-09,6.1559E-09,8.4640E-09,5.9240E-09,5.8217E-09/
      DATA (KA(JT, 2, 2),JT=1,5) /
     &8.3749E-09,8.0756E-09,1.1623E-08,1.1272E-08,7.3636E-09/
      DATA (KA(JT, 3, 2),JT=1,5) /
     &1.3304E-08,1.2795E-08,1.2343E-08,2.1235E-08,1.1577E-08/
      DATA (KA(JT, 4, 2),JT=1,5) /
     &2.0704E-08,1.9736E-08,1.8900E-08,1.8228E-08,3.1601E-08/
      DATA (KA(JT, 5, 2),JT=1,5) /
     &3.1149E-08,2.9669E-08,2.8318E-08,2.7101E-08,4.9649E-08/
      DATA (KA(JT, 6, 2),JT=1,5) /
     &4.5713E-08,4.3519E-08,4.1488E-08,3.9918E-08,3.8291E-08/
      DATA (KA(JT, 7, 2),JT=1,5) /
     &7.7265E-08,7.3848E-08,7.0437E-08,6.7945E-08,6.6127E-08/
      DATA (KA(JT, 8, 2),JT=1,5) /
     &1.5754E-07,1.5664E-07,1.5378E-07,1.5027E-07,1.4633E-07/
      DATA (KA(JT, 9, 2),JT=1,5) /
     &1.6439E-07,1.4678E-07,1.2610E-07,1.1532E-07,1.0591E-07/
      DATA (KA(JT,10, 2),JT=1,5) /
     &1.4366E-07,1.3506E-07,1.2583E-07,1.1774E-07,1.1011E-07/
      DATA (KA(JT,11, 2),JT=1,5) /
     &1.4521E-07,1.3766E-07,1.3072E-07,1.2218E-07,1.1400E-07/
      DATA (KA(JT,12, 2),JT=1,5) /
     &1.4524E-07,1.3769E-07,1.3074E-07,1.2241E-07,1.1552E-07/
      DATA (KA(JT,13, 2),JT=1,5) /
     &1.4525E-07,1.3770E-07,1.3075E-07,1.2252E-07,1.1553E-07/
      DATA (KA(JT, 1, 3),JT=1,5) /
     &1.4060E-08,1.3587E-08,2.4644E-08,1.2716E-08,1.2367E-08/
      DATA (KA(JT, 2, 3),JT=1,5) /
     &1.7055E-08,1.6577E-08,3.2443E-08,3.1273E-08,1.5381E-08/
      DATA (KA(JT, 3, 3),JT=1,5) /
     &2.5414E-08,2.4672E-08,2.3874E-08,4.7281E-08,2.2346E-08/
      DATA (KA(JT, 4, 3),JT=1,5) /
     &3.9536E-08,3.8124E-08,3.6836E-08,3.5587E-08,7.2260E-08/
      DATA (KA(JT, 5, 3),JT=1,5) /
     &5.9488E-08,5.7630E-08,5.5623E-08,5.3878E-08,1.1230E-07/
      DATA (KA(JT, 6, 3),JT=1,5) /
     &9.9996E-08,9.6206E-08,9.3184E-08,9.0812E-08,8.9206E-08/
      DATA (KA(JT, 7, 3),JT=1,5) /
     &1.7678E-07,1.7554E-07,1.7358E-07,1.7091E-07,1.6830E-07/
      DATA (KA(JT, 8, 3),JT=1,5) /
     &1.8672E-07,1.7850E-07,1.6967E-07,1.6275E-07,1.5875E-07/
      DATA (KA(JT, 9, 3),JT=1,5) /
     &1.3558E-07,1.3493E-07,1.3633E-07,1.3799E-07,1.3932E-07/
      DATA (KA(JT,10, 3),JT=1,5) /
     &1.8883E-07,2.0452E-07,2.2206E-07,2.4347E-07,2.6091E-07/
      DATA (KA(JT,11, 3),JT=1,5) /
     &2.1296E-07,2.3580E-07,2.6439E-07,3.0148E-07,3.4942E-07/
      DATA (KA(JT,12, 3),JT=1,5) /
     &2.2072E-07,2.5535E-07,2.8661E-07,3.4814E-07,3.9337E-07/
      DATA (KA(JT,13, 3),JT=1,5) /
     &2.2515E-07,2.6161E-07,3.0833E-07,3.6527E-07,4.0123E-07/
      DATA (KA(JT, 1, 4),JT=1,5) /
     &3.2735E-08,3.1345E-08,5.8846E-08,2.8258E-08,2.7022E-08/
      DATA (KA(JT, 2, 4),JT=1,5) /
     &3.7754E-08,3.6873E-08,8.2776E-08,8.0947E-08,3.5190E-08/
      DATA (KA(JT, 3, 4),JT=1,5) /
     &7.6368E-08,7.5292E-08,7.4075E-08,1.0820E-07,7.3183E-08/
      DATA (KA(JT, 4, 4),JT=1,5) /
     &1.6392E-07,1.6130E-07,1.5926E-07,1.5700E-07,2.2041E-07/
      DATA (KA(JT, 5, 4),JT=1,5) /
     &2.9704E-07,2.8924E-07,2.8301E-07,2.7633E-07,4.2284E-07/
      DATA (KA(JT, 6, 4),JT=1,5) /
     &4.8466E-07,4.7240E-07,4.6143E-07,4.5012E-07,4.3867E-07/
      DATA (KA(JT, 7, 4),JT=1,5) /
     &7.1637E-07,6.9847E-07,6.7384E-07,6.5368E-07,6.3375E-07/
      DATA (KA(JT, 8, 4),JT=1,5) /
     &1.1904E-06,1.1714E-06,1.1524E-06,1.1354E-06,1.1172E-06/
      DATA (KA(JT, 9, 4),JT=1,5) /
     &2.1976E-06,2.1606E-06,2.1332E-06,2.0944E-06,2.0536E-06/
      DATA (KA(JT,10, 4),JT=1,5) /
     &2.1713E-06,2.1144E-06,2.0553E-06,1.9901E-06,1.9286E-06/
      DATA (KA(JT,11, 4),JT=1,5) /
     &2.1443E-06,2.0785E-06,2.0048E-06,1.9232E-06,1.8295E-06/
      DATA (KA(JT,12, 4),JT=1,5) /
     &2.1363E-06,2.0578E-06,1.9811E-06,1.8729E-06,1.7807E-06/
      DATA (KA(JT,13, 4),JT=1,5) /
     &2.1319E-06,2.0513E-06,1.9580E-06,1.8546E-06,1.7725E-06/
      DATA (KA(JT, 1, 5),JT=1,5) /
     &3.6050E-08,3.6125E-08,4.6253E-08,3.7280E-08,3.7359E-08/
      DATA (KA(JT, 2, 5),JT=1,5) /
     &6.5102E-08,6.4266E-08,6.8896E-08,6.5925E-08,6.1190E-08/
      DATA (KA(JT, 3, 5),JT=1,5) /
     &1.2173E-07,1.1889E-07,1.1625E-07,1.7574E-07,1.0921E-07/
      DATA (KA(JT, 4, 5),JT=1,5) /
     &2.0555E-07,1.9853E-07,1.9068E-07,1.8313E-07,3.0241E-07/
      DATA (KA(JT, 5, 5),JT=1,5) /
     &3.0900E-07,2.9996E-07,2.8857E-07,2.7772E-07,5.1631E-07/
      DATA (KA(JT, 6, 5),JT=1,5) /
     &4.3774E-07,4.2465E-07,4.0920E-07,3.9315E-07,3.7901E-07/
      DATA (KA(JT, 7, 5),JT=1,5) /
     &6.3869E-07,6.1654E-07,6.0324E-07,5.8966E-07,5.7948E-07/
      DATA (KA(JT, 8, 5),JT=1,5) /
     &9.8362E-07,9.6271E-07,9.4180E-07,9.2206E-07,9.1105E-07/
      DATA (KA(JT, 9, 5),JT=1,5) /
     &1.2061E-06,1.1895E-06,1.1564E-06,1.1296E-06,1.1110E-06/
      DATA (KA(JT,10, 5),JT=1,5) /
     &1.2958E-06,1.2694E-06,1.2425E-06,1.2153E-06,1.1880E-06/
      DATA (KA(JT,11, 5),JT=1,5) /
     &1.2962E-06,1.2698E-06,1.2429E-06,1.2156E-06,1.1883E-06/
      DATA (KA(JT,12, 5),JT=1,5) /
     &1.2964E-06,1.2701E-06,1.2431E-06,1.2158E-06,1.1885E-06/
      DATA (KA(JT,13, 5),JT=1,5) /
     &1.2966E-06,1.2702E-06,1.2433E-06,1.2160E-06,1.1886E-06/
      DATA (KA(JT, 1, 6),JT=1,5) /
     &7.3925E-08,7.0231E-08,2.1454E-07,6.3477E-08,6.0912E-08/
      DATA (KA(JT, 2, 6),JT=1,5) /
     &6.7794E-08,6.5807E-08,1.3854E-07,1.3061E-07,5.9361E-08/
      DATA (KA(JT, 3, 6),JT=1,5) /
     &9.8353E-08,9.5275E-08,9.2426E-08,1.5768E-07,8.7986E-08/
      DATA (KA(JT, 4, 6),JT=1,5) /
     &1.5855E-07,1.5394E-07,1.4948E-07,1.4655E-07,2.3172E-07/
      DATA (KA(JT, 5, 6),JT=1,5) /
     &2.7764E-07,2.6941E-07,2.6299E-07,2.5975E-07,4.0526E-07/
      DATA (KA(JT, 6, 6),JT=1,5) /
     &4.5469E-07,4.4417E-07,4.3276E-07,4.2440E-07,4.1489E-07/
      DATA (KA(JT, 7, 6),JT=1,5) /
     &7.1540E-07,7.1291E-07,7.0656E-07,6.9823E-07,6.8342E-07/
      DATA (KA(JT, 8, 6),JT=1,5) /
     &7.9651E-07,7.9807E-07,8.0621E-07,8.0941E-07,7.9835E-07/
      DATA (KA(JT, 9, 6),JT=1,5) /
     &1.8716E-07,1.6713E-07,1.4725E-07,1.3728E-07,1.1763E-07/
      DATA (KA(JT,10, 6),JT=1,5) /
     &9.2638E-08,8.6207E-08,8.0877E-08,7.0432E-08,6.4517E-08/
      DATA (KA(JT,11, 6),JT=1,5) /
     &1.3396E-07,1.2820E-07,1.2387E-07,1.0427E-07,9.4091E-08/
      DATA (KA(JT,12, 6),JT=1,5) /
     &1.4877E-07,1.4827E-07,1.4350E-07,1.2154E-07,1.0552E-07/
      DATA (KA(JT,13, 6),JT=1,5) /
     &1.5437E-07,1.5323E-07,1.4992E-07,1.2715E-07,1.0933E-07/
      DATA (KA(JT, 1, 7),JT=1,5) /
     &7.2717E-07,7.0656E-07,1.3933E-06,6.6449E-07,6.4269E-07/
      DATA (KA(JT, 2, 7),JT=1,5) /
     &5.2595E-07,5.0791E-07,1.1171E-06,1.0538E-06,4.5644E-07/
      DATA (KA(JT, 3, 7),JT=1,5) /
     &2.9919E-07,2.9227E-07,2.8284E-07,6.5215E-07,2.6347E-07/
      DATA (KA(JT, 4, 7),JT=1,5) /
     &2.7961E-07,2.7579E-07,2.7068E-07,2.6343E-07,4.1265E-07/
      DATA (KA(JT, 5, 7),JT=1,5) /
     &3.7031E-07,3.6318E-07,3.5475E-07,3.4488E-07,5.3740E-07/
      DATA (KA(JT, 6, 7),JT=1,5) /
     &5.3195E-07,5.2692E-07,5.2224E-07,5.1934E-07,5.1146E-07/
      DATA (KA(JT, 7, 7),JT=1,5) /
     &8.3043E-07,8.4552E-07,8.4833E-07,8.2800E-07,8.0930E-07/
      DATA (KA(JT, 8, 7),JT=1,5) /
     &1.4910E-06,1.5179E-06,1.5248E-06,1.5091E-06,1.4853E-06/
      DATA (KA(JT, 9, 7),JT=1,5) /
     &3.7340E-06,3.7823E-06,3.8311E-06,3.8453E-06,3.8567E-06/
      DATA (KA(JT,10, 7),JT=1,5) /
     &8.6791E-06,8.9697E-06,9.2118E-06,9.3991E-06,9.5564E-06/
      DATA (KA(JT,11, 7),JT=1,5) /
     &1.1878E-05,1.2201E-05,1.2588E-05,1.2897E-05,1.3151E-05/
      DATA (KA(JT,12, 7),JT=1,5) /
     &1.3192E-05,1.3732E-05,1.4137E-05,1.4465E-05,1.4643E-05/
      DATA (KA(JT,13, 7),JT=1,5) /
     &1.3716E-05,1.4229E-05,1.4617E-05,1.4944E-05,1.5182E-05/
      DATA (KA(JT, 1, 8),JT=1,5) /
     &3.9538E-06,3.8949E-06,5.6188E-06,3.7475E-06,3.6648E-06/
      DATA (KA(JT, 2, 8),JT=1,5) /
     &3.4231E-06,3.3633E-06,5.1877E-06,5.0048E-06,3.1425E-06/
      DATA (KA(JT, 3, 8),JT=1,5) /
     &2.8073E-06,2.7497E-06,2.6875E-06,4.4405E-06,2.5492E-06/
      DATA (KA(JT, 4, 8),JT=1,5) /
     &1.9229E-06,1.8818E-06,1.8382E-06,1.7896E-06,3.3073E-06/
      DATA (KA(JT, 5, 8),JT=1,5) /
     &1.1453E-06,1.1293E-06,1.1095E-06,1.0866E-06,1.9344E-06/
      DATA (KA(JT, 6, 8),JT=1,5) /
     &1.4565E-06,1.4517E-06,1.4369E-06,1.4141E-06,1.3944E-06/
      DATA (KA(JT, 7, 8),JT=1,5) /
     &2.3228E-06,2.2753E-06,2.2395E-06,2.2124E-06,2.1731E-06/
      DATA (KA(JT, 8, 8),JT=1,5) /
     &3.4877E-06,3.4362E-06,3.3796E-06,3.3389E-06,3.2924E-06/
      DATA (KA(JT, 9, 8),JT=1,5) /
     &6.3448E-06,6.3701E-06,6.3619E-06,6.2632E-06,6.1645E-06/
      DATA (KA(JT,10, 8),JT=1,5) /
     &1.2155E-05,1.1880E-05,1.1762E-05,1.1759E-05,1.1651E-05/
      DATA (KA(JT,11, 8),JT=1,5) /
     &1.4093E-05,1.3835E-05,1.3547E-05,1.3205E-05,1.2690E-05/
      DATA (KA(JT,12, 8),JT=1,5) /
     &1.4428E-05,1.4056E-05,1.3932E-05,1.3396E-05,1.2885E-05/
      DATA (KA(JT,13, 8),JT=1,5) /
     &1.5229E-05,1.4534E-05,1.3849E-05,1.3292E-05,1.2704E-05/
      DATA (KA(JT, 1, 9),JT=1,5) /
     &1.9250E-05,1.9148E-05,2.1702E-05,1.8906E-05,1.8761E-05/
      DATA (KA(JT, 2, 9),JT=1,5) /
     &1.8132E-05,1.8040E-05,2.0884E-05,2.0523E-05,1.7656E-05/
      DATA (KA(JT, 3, 9),JT=1,5) /
     &1.6928E-05,1.6843E-05,1.6742E-05,1.9715E-05,1.6470E-05/
      DATA (KA(JT, 4, 9),JT=1,5) /
     &1.5526E-05,1.5463E-05,1.5377E-05,1.5268E-05,1.8367E-05/
      DATA (KA(JT, 5, 9),JT=1,5) /
     &1.3545E-05,1.3511E-05,1.3455E-05,1.3362E-05,1.6722E-05/
      DATA (KA(JT, 6, 9),JT=1,5) /
     &9.7183E-06,9.7218E-06,9.7084E-06,9.6717E-06,9.6030E-06/
      DATA (KA(JT, 7, 9),JT=1,5) /
     &5.0307E-06,5.0984E-06,5.1628E-06,5.2093E-06,5.2354E-06/
      DATA (KA(JT, 8, 9),JT=1,5) /
     &4.5837E-06,4.5939E-06,4.5938E-06,4.5639E-06,4.5109E-06/
      DATA (KA(JT, 9, 9),JT=1,5) /
     &1.2254E-05,1.2319E-05,1.2397E-05,1.2584E-05,1.2620E-05/
      DATA (KA(JT,10, 9),JT=1,5) /
     &2.1545E-05,2.1836E-05,2.1718E-05,2.1511E-05,2.1211E-05/
      DATA (KA(JT,11, 9),JT=1,5) /
     &2.0079E-05,1.9539E-05,1.8859E-05,1.8393E-05,1.8181E-05/
      DATA (KA(JT,12, 9),JT=1,5) /
     &1.7115E-05,1.6357E-05,1.5410E-05,1.5220E-05,1.5207E-05/
      DATA (KA(JT,13, 9),JT=1,5) /
     &1.4935E-05,1.4679E-05,1.4593E-05,1.4448E-05,1.4436E-05/
      DATA (KA(JT, 1,10),JT=1,5) /
     &5.3569E-05,5.3042E-05,5.5454E-05,5.2098E-05,5.1678E-05/
      DATA (KA(JT, 2,10),JT=1,5) /
     &5.2196E-05,5.1739E-05,5.4777E-05,5.4075E-05,5.0624E-05/
      DATA (KA(JT, 3,10),JT=1,5) /
     &5.0339E-05,5.0046E-05,4.9769E-05,5.3168E-05,4.9370E-05/
      DATA (KA(JT, 4,10),JT=1,5) /
     &4.8505E-05,4.8316E-05,4.8143E-05,4.7993E-05,5.1621E-05/
      DATA (KA(JT, 5,10),JT=1,5) /
     &4.6313E-05,4.6267E-05,4.6119E-05,4.6064E-05,5.0279E-05/
      DATA (KA(JT, 6,10),JT=1,5) /
     &4.2662E-05,4.2818E-05,4.2935E-05,4.3007E-05,4.3099E-05/
      DATA (KA(JT, 7,10),JT=1,5) /
     &3.5762E-05,3.6149E-05,3.6450E-05,3.6639E-05,3.6887E-05/
      DATA (KA(JT, 8,10),JT=1,5) /
     &1.3516E-06,1.8607E-06,2.3061E-06,2.7339E-06,3.6516E-06/
      DATA (KA(JT, 9,10),JT=1,5) /
     &3.6432E-06,4.0739E-06,4.3830E-06,4.1136E-06,4.3128E-06/
      DATA (KA(JT,10,10),JT=1,5) /
     &6.2049E-06,6.9116E-06,7.3244E-06,6.5087E-06,7.8951E-06/
      DATA (KA(JT,11,10),JT=1,5) /
     &3.2156E-06,3.8834E-06,4.1231E-06,4.3386E-06,4.3405E-06/
      DATA (KA(JT,12,10),JT=1,5) /
     &2.2152E-06,2.6754E-06,3.1971E-06,3.4911E-06,3.7935E-06/
      DATA (KA(JT,13,10),JT=1,5) /
     &1.9792E-06,2.6543E-06,3.1511E-06,3.4597E-06,4.0624E-06/
      DATA (KA(JT, 1,11),JT=1,5) /
     &7.5384E-05,7.5103E-05,7.7406E-05,7.4222E-05,7.3734E-05/
      DATA (KA(JT, 2,11),JT=1,5) /
     &7.5458E-05,7.5244E-05,7.7778E-05,7.7018E-05,7.3942E-05/
      DATA (KA(JT, 3,11),JT=1,5) /
     &7.5023E-05,7.4844E-05,7.4477E-05,7.7271E-05,7.3633E-05/
      DATA (KA(JT, 4,11),JT=1,5) /
     &7.3633E-05,7.3539E-05,7.3257E-05,7.2934E-05,7.6232E-05/
      DATA (KA(JT, 5,11),JT=1,5) /
     &7.1348E-05,7.1322E-05,7.1227E-05,7.1069E-05,7.5258E-05/
      DATA (KA(JT, 6,11),JT=1,5) /
     &6.7784E-05,6.7873E-05,6.7974E-05,6.7924E-05,6.7903E-05/
      DATA (KA(JT, 7,11),JT=1,5) /
     &6.1855E-05,6.1922E-05,6.1973E-05,6.2206E-05,6.2496E-05/
      DATA (KA(JT, 8,11),JT=1,5) /
     &3.6622E-05,3.7413E-05,3.8740E-05,4.0550E-05,4.1833E-05/
      DATA (KA(JT, 9,11),JT=1,5) /
     &2.8544E-06,2.8831E-06,3.1445E-06,3.2900E-06,2.7967E-06/
      DATA (KA(JT,10,11),JT=1,5) /
     &5.3755E-06,4.2123E-06,5.1154E-06,6.3481E-06,5.4219E-06/
      DATA (KA(JT,11,11),JT=1,5) /
     &1.2605E-06,1.4078E-06,1.9167E-06,2.3729E-06,3.0161E-06/
      DATA (KA(JT,12,11),JT=1,5) /
     &1.1370E-06,9.1524E-07,1.1150E-06,1.4746E-06,2.0128E-06/
      DATA (KA(JT,13,11),JT=1,5) /
     &1.0511E-06,1.0014E-06,1.1405E-06,1.3852E-06,1.5576E-06/
      DATA (KA(JT, 1,12),JT=1,5) /
     &1.1184E-04,1.1117E-04,1.1327E-04,1.0989E-04,1.0910E-04/
      DATA (KA(JT, 2,12),JT=1,5) /
     &1.1379E-04,1.1322E-04,1.1555E-04,1.1462E-04,1.1135E-04/
      DATA (KA(JT, 3,12),JT=1,5) /
     &1.1508E-04,1.1459E-04,1.1421E-04,1.1671E-04,1.1339E-04/
      DATA (KA(JT, 4,12),JT=1,5) /
     &1.1596E-04,1.1563E-04,1.1538E-04,1.1511E-04,1.1770E-04/
      DATA (KA(JT, 5,12),JT=1,5) /
     &1.1597E-04,1.1581E-04,1.1569E-04,1.1553E-04,1.1890E-04/
      DATA (KA(JT, 6,12),JT=1,5) /
     &1.1443E-04,1.1445E-04,1.1443E-04,1.1443E-04,1.1438E-04/
      DATA (KA(JT, 7,12),JT=1,5) /
     &1.0852E-04,1.0888E-04,1.0912E-04,1.0934E-04,1.0942E-04/
      DATA (KA(JT, 8,12),JT=1,5) /
     &9.3194E-05,9.4766E-05,9.5355E-05,9.5090E-05,9.4926E-05/
      DATA (KA(JT, 9,12),JT=1,5) /
     &1.1836E-06,1.6115E-06,1.2883E-06,1.4202E-06,1.6541E-06/
      DATA (KA(JT,10,12),JT=1,5) /
     &1.8748E-06,3.4401E-06,3.9984E-06,4.4576E-06,3.3683E-06/
      DATA (KA(JT,11,12),JT=1,5) /
     &2.9890E-07,4.8741E-07,6.6276E-07,9.9698E-07,1.9230E-06/
      DATA (KA(JT,12,12),JT=1,5) /
     &1.5034E-07,3.9966E-07,5.6523E-07,7.0494E-07,1.0046E-06/
      DATA (KA(JT,13,12),JT=1,5) /
     &1.5016E-07,2.5751E-07,4.8928E-07,6.3534E-07,9.3575E-07/
      DATA (KA(JT, 1,13),JT=1,5) /
     &1.7305E-04,1.7234E-04,1.7389E-04,1.7055E-04,1.6974E-04/
      DATA (KA(JT, 2,13),JT=1,5) /
     &1.8170E-04,1.8075E-04,1.8265E-04,1.8138E-04,1.7772E-04/
      DATA (KA(JT, 3,13),JT=1,5) /
     &1.8990E-04,1.8892E-04,1.8776E-04,1.8950E-04,1.8494E-04/
      DATA (KA(JT, 4,13),JT=1,5) /
     &1.9649E-04,1.9552E-04,1.9424E-04,1.9281E-04,1.9464E-04/
      DATA (KA(JT, 5,13),JT=1,5) /
     &2.0197E-04,2.0109E-04,1.9993E-04,1.9856E-04,2.0092E-04/
      DATA (KA(JT, 6,13),JT=1,5) /
     &2.0595E-04,2.0549E-04,2.0452E-04,2.0331E-04,2.0199E-04/
      DATA (KA(JT, 7,13),JT=1,5) /
     &2.0703E-04,2.0710E-04,2.0649E-04,2.0552E-04,2.0428E-04/
      DATA (KA(JT, 8,13),JT=1,5) /
     &1.9874E-04,1.9767E-04,1.9696E-04,1.9655E-04,1.9591E-04/
      DATA (KA(JT, 9,13),JT=1,5) /
     &2.0434E-05,2.3398E-05,2.7400E-05,3.2409E-05,3.8451E-05/
      DATA (KA(JT,10,13),JT=1,5) /
     &1.8617E-06,9.9513E-07,1.0554E-06,1.6516E-06,3.7792E-06/
      DATA (KA(JT,11,13),JT=1,5) /
     &1.2517E-07,2.9518E-07,7.7058E-07,1.1660E-06,1.5349E-06/
      DATA (KA(JT,12,13),JT=1,5) /
     &1.2734E-07,3.6524E-07,6.6699E-07,1.0362E-06,1.4158E-06/
      DATA (KA(JT,13,13),JT=1,5) /
     &1.2431E-07,3.9389E-07,6.7331E-07,1.0292E-06,1.4448E-06/
      DATA (KA(JT, 1,14),JT=1,5) /
     &2.9365E-04,2.9046E-04,2.9008E-04,2.8509E-04,2.8286E-04/
      DATA (KA(JT, 2,14),JT=1,5) /
     &3.1990E-04,3.1668E-04,3.1617E-04,3.1332E-04,3.0885E-04/
      DATA (KA(JT, 3,14),JT=1,5) /
     &3.4787E-04,3.4432E-04,3.4112E-04,3.4052E-04,3.3589E-04/
      DATA (KA(JT, 4,14),JT=1,5) /
     &3.7401E-04,3.7027E-04,3.6696E-04,3.6394E-04,3.6355E-04/
      DATA (KA(JT, 5,14),JT=1,5) /
     &3.9840E-04,3.9446E-04,3.9082E-04,3.8763E-04,3.8760E-04/
      DATA (KA(JT, 6,14),JT=1,5) /
     &4.2165E-04,4.1729E-04,4.1335E-04,4.1006E-04,4.0721E-04/
      DATA (KA(JT, 7,14),JT=1,5) /
     &4.4257E-04,4.3782E-04,4.3364E-04,4.3014E-04,4.2736E-04/
      DATA (KA(JT, 8,14),JT=1,5) /
     &4.5299E-04,4.4953E-04,4.4586E-04,4.4260E-04,4.4006E-04/
      DATA (KA(JT, 9,14),JT=1,5) /
     &4.0190E-04,3.9751E-04,3.9238E-04,3.8812E-04,3.8612E-04/
      DATA (KA(JT,10,14),JT=1,5) /
     &6.4278E-06,1.8248E-06,1.6996E-06,3.1086E-07,1.6836E-07/
      DATA (KA(JT,11,14),JT=1,5) /
     &1.4350E-06,9.4778E-07,4.1349E-07,2.0817E-07,2.0238E-07/
      DATA (KA(JT,12,14),JT=1,5) /
     &1.6805E-06,1.5323E-06,6.2348E-07,9.9743E-08,1.2977E-07/
      DATA (KA(JT,13,14),JT=1,5) /
     &1.6858E-06,1.7103E-06,8.0574E-07,1.5825E-07,1.5032E-07/
      DATA (KA(JT, 1,15),JT=1,5) /
     &5.2181E-04,5.1578E-04,5.1251E-04,5.0356E-04,4.9731E-04/
      DATA (KA(JT, 2,15),JT=1,5) /
     &5.9491E-04,5.8822E-04,5.8413E-04,5.7646E-04,5.6692E-04/
      DATA (KA(JT, 3,15),JT=1,5) /
     &6.7653E-04,6.6881E-04,6.6126E-04,6.5540E-04,6.4461E-04/
      DATA (KA(JT, 4,15),JT=1,5) /
     &7.6388E-04,7.5456E-04,7.4556E-04,7.3649E-04,7.2840E-04/
      DATA (KA(JT, 5,15),JT=1,5) /
     &8.5507E-04,8.4417E-04,8.3378E-04,8.2338E-04,8.1349E-04/
      DATA (KA(JT, 6,15),JT=1,5) /
     &9.5034E-04,9.3798E-04,9.2553E-04,9.1287E-04,8.9957E-04/
      DATA (KA(JT, 7,15),JT=1,5) /
     &1.0496E-03,1.0352E-03,1.0206E-03,1.0054E-03,9.8958E-04/
      DATA (KA(JT, 8,15),JT=1,5) /
     &1.1507E-03,1.1337E-03,1.1169E-03,1.0991E-03,1.0806E-03/
      DATA (KA(JT, 9,15),JT=1,5) /
     &1.2408E-03,1.2207E-03,1.1996E-03,1.1773E-03,1.1531E-03/
      DATA (KA(JT,10,15),JT=1,5) /
     &1.2042E-04,1.1501E-04,1.1424E-04,1.1450E-04,1.3219E-04/
      DATA (KA(JT,11,15),JT=1,5) /
     &6.8914E-07,8.3960E-07,7.4591E-07,1.8660E-06,3.2503E-06/
      DATA (KA(JT,12,15),JT=1,5) /
     &3.5963E-08,4.6256E-07,5.6223E-07,9.8816E-07,9.2366E-07/
      DATA (KA(JT,13,15),JT=1,5) /
     &3.6605E-08,5.6591E-07,8.4008E-07,8.6042E-07,6.8452E-07/
      DATA (KA(JT, 1,16),JT=1,5) /
     &7.6517E-04,7.5944E-04,7.6010E-04,7.6100E-04,7.6498E-04/
      DATA (KA(JT, 2,16),JT=1,5) /
     &9.2375E-04,9.1357E-04,9.0997E-04,9.0997E-04,9.0993E-04/
      DATA (KA(JT, 3,16),JT=1,5) /
     &1.1142E-03,1.0974E-03,1.0835E-03,1.0789E-03,1.0748E-03/
      DATA (KA(JT, 4,16),JT=1,5) /
     &1.3278E-03,1.3025E-03,1.2802E-03,1.2631E-03,1.2539E-03/
      DATA (KA(JT, 5,16),JT=1,5) /
     &1.5712E-03,1.5343E-03,1.5017E-03,1.4713E-03,1.4497E-03/
      DATA (KA(JT, 6,16),JT=1,5) /
     &1.8525E-03,1.7982E-03,1.7525E-03,1.7101E-03,1.6714E-03/
      DATA (KA(JT, 7,16),JT=1,5) /
     &2.1731E-03,2.0986E-03,2.0340E-03,1.9757E-03,1.9210E-03/
      DATA (KA(JT, 8,16),JT=1,5) /
     &2.5325E-03,2.4346E-03,2.3473E-03,2.2687E-03,2.1950E-03/
      DATA (KA(JT, 9,16),JT=1,5) /
     &2.9269E-03,2.8006E-03,2.6863E-03,2.5805E-03,2.4878E-03/
      DATA (KA(JT,10,16),JT=1,5) /
     &2.9442E-03,2.7008E-03,2.3913E-03,2.1437E-03,1.8865E-03/
      DATA (KA(JT,11,16),JT=1,5) /
     &2.3220E-06,2.2310E-05,4.8349E-05,6.7183E-05,8.8908E-05/
      DATA (KA(JT,12,16),JT=1,5) /
     &2.2857E-06,1.1848E-05,4.2066E-05,6.7613E-05,8.6033E-05/
      DATA (KA(JT,13,16),JT=1,5) /
     &2.2823E-06,6.9105E-06,3.6212E-05,6.6247E-05,8.5488E-05/
