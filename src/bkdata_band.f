C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

      block data rrtm_band

      include 'param.f'

      COMMON /CONSTANTS/ PI,FLUXFAC,HEATFAC
      COMMON /FEATURES/  NG(IB1:IB2),NSPA(IB1:IB2),NSPB(IB1:IB2)
      COMMON /HVERSN/    HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *                   HVDUM1(4),HVRUTL,HVREXT,
     *                   HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *                   HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

      CHARACTER*15 HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVDUM1,HVRUTL,HVREXT,
     *            HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *            HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

      DATA HVRBKB /'$Revision$'/
      DATA WAVENUM1(16) /2600./,WAVENUM2(16) /3250./,DELWAVE(16) /650./
      DATA WAVENUM1(17) /3250./,WAVENUM2(17) /4000./,DELWAVE(17) /750./
      DATA WAVENUM1(18) /4000./,WAVENUM2(18) /4650./,DELWAVE(18) /650./
      DATA WAVENUM1(19) /4650./,WAVENUM2(19) /5150./,DELWAVE(19) /500./
      DATA WAVENUM1(20) /5150./,WAVENUM2(20) /6150./,DELWAVE(20) /1000./
      DATA WAVENUM1(21) /6150./,WAVENUM2(21) /7700./,DELWAVE(21) /1550./
      DATA WAVENUM1(22) /7700./,WAVENUM2(22) /8050./,DELWAVE(22) /350./
      DATA WAVENUM1(23) /8050./,WAVENUM2(23)/12850./,DELWAVE(23) /4800./
      DATA WAVENUM1(24)/12850./,WAVENUM2(24)/16000./,DELWAVE(24) /3150./
      DATA WAVENUM1(25)/16000./,WAVENUM2(25)/22650./,DELWAVE(25) /6650./
      DATA WAVENUM1(26)/22650./,WAVENUM2(26)/29000./,DELWAVE(26) /6350./
      DATA WAVENUM1(27)/29000./,WAVENUM2(27)/38000./,DELWAVE(27) /9000./
      DATA WAVENUM1(28)/38000./,WAVENUM2(28)/50000./,DELWAVE(28)/12000./
      DATA WAVENUM1(29)/820./,  WAVENUM2(29)/2600./, DELWAVE(29)/1780./

      DATA NG  /16,16,16,16,16,16,16,16,16,16,16,16,16,16/
      DATA NSPA /9, 9, 9, 9, 1, 9, 9, 1, 9, 1, 0, 1, 9, 1/
      DATA NSPB /1, 5, 1, 1, 1, 5, 1, 0, 1, 0, 0, 1, 5, 1/

C     HEATFAC is the factor by which one must multiply delta-flux/ 
C     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
C     the heating rate in units of degrees/day.  It is equal to 
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
      DATA HEATFAC /8.4391/

        end
