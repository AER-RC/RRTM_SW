C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$
      SUBROUTINE RTRDIS

C *** This program calculates the upward fluxes, downward fluxes,
C     and heating rates for an arbitrary atmosphere.  The input to
C     this program is the atmospheric profile and all Planck function
C     information.  First-order "numerical" quadrature is used for the 
C     angle integration, i.e. only one exponential is computed per layer
C     per g-value per band.

      IMPLICIT DOUBLE PRECISION (V)                                     
      PARAMETER (MXANG = 4)
      PARAMETER ( MCMU = 32, MUMU = 32,
     &          MPHI = 3)

      INCLUDE 'param.f'

      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /FEATURES/  NG(IB1:IB2),NSPA(IB1:IB2),NSPB(IB1:IB2)
      COMMON /CONTROL/   IAER, NSTR, IOUT, ISTART, IEND, ICLD,
     &                   idelm, isccos
      COMMON /SWPROP/    ZENITH, ALBEDO, ADJFLUX
      COMMON /SURFACE/   IREFLECT,SEMISS(NBANDS)
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /SOLARIN/   SFLUXZEN(MG)
      COMMON /PLANKG/    FRACS(MXLAY,MG)
      COMMON /TAUGCOM/   TAUG(MXLAY,MG)
      COMMON /SSAGCOM/   SSA(MXLAY,MG)
      COMMON /AERDAT/    ssaaer(mxlay,nbands), phase(mcmu,mxlay,nbands), 
     &                   tauaer(mxlay,nbands)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),TAUCLOUD(MXLAY,NBANDS),
     &                   SSACLOUD(MXLAY,NBANDS),
     &                   XMOM(0:16,MXLAY,NBANDS),
     &                   TAUCLDORIG(MXLAY,NBANDS)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY), TOTDFLUX(0:MXLAY),
     &                   DIFDOWN(0:MXLAY), DIRDOWN(0:MXLAY),
     &                   FNET(0:MXLAY), HTR(0:MXLAY)
      COMMON /HVERSN/    HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *                   HVDUM1(4),HVRUTL,HVREXT,
     *                   HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *                   HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

      CHARACTER*15 HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVDUM1,HVRUTL,HVREXT,
     *            HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *            HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR
                                       
      DIMENSION BBD1(MXLAY),BBD2(MXLAY),BBD3(MXLAY)
      DIMENSION WTNUM(MXANG)


      CHARACTER HEADER*127
      LOGICAL   DELTAM, LAMBER, ONLYFL, PLANK, USRANG, USRTAU
      INTEGER   IBCND, MAXCLY, MAXCMU, MAXPHI, MAXULV, MAXUMU, 
     &          NPHI, NSTR, NTAU, NUMU
      REAL      ACCUR, ALBEDO, BTEMP, FBEAM, FISOT, PHI0, TEMIS, TTEMP,
     &          UMU0, WVNMHI, WVNMLO
      LOGICAL   PRNT( 7 ), sccos
      REAL      ALBMED( MUMU ), DFDT( MXLAY ), TAUREV(MXLAY),
     &          FLUP( MXLAY ), HL( 0:MCMU ), PHI( MPHI ),
     &          PMOM( 0:MCMU, MXLAY ), RFLDIR( MXLAY ),
     &          RFLDN( MXLAY ), SSALB( MXLAY ), TEMPER( 0:MXLAY ),
     &          TRNMED( MUMU ), U0U( MUMU, MXLAY ), UAVG( MXLAY ),
     &          UMU( MUMU ), UTAU( MXLAY ),
     &          UU( MUMU, MXLAY, MPHI ),fldir(mxlay),fldn(mxlay)
      DIMENSION PHASERAY(0:MXSTR)


      DATA PRNT /.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,
     &     .FALSE.,.FALSE./

      HVRRTR = '$Revision$'  

      sccos = .false.
      if (isccos .gt. 0) sccos = .true.

      PHASERAY(0) = 1.0
      PHASERAY(1) = 0. 
      PHASERAY(2) = 0.1
      DO 120 ISTR = 3, NSTR
         PHASERAY (ISTR) = 0.
 120  CONTINUE

      SECZEN = 1. / ZENITH
      HEADER = ''
      USRTAU = .FALSE.
      USRANG = .FALSE.
      NPHI = 0
      IBCND = 0
      PHI0 = 0.
      IF (IREFLECT .EQ. 0) THEN
         LAMBER = .TRUE.
      ELSE
         LAMBER = .FALSE.
      ENDIF

      DELTAM = .FALSE.
      PLANK = .FALSE.
      ONLYFL = .TRUE.
      ACCUR = 0.0001
      MAXCLY = MXLAY
      MAXULV = MXLAY
      MAXUMU = MUMU
      MAXCMU = MCMU
      MAXPHI = MPHI
      DO 200 LAY = 0, NLAYERS
         IF (LAY .NE. 0) THEN
            SSALB(LAY) = 0.
            DO 180 IQ = 0, NSTR
               PMOM(IQ,LAY) = 0.
 180        CONTINUE
         ENDIF

         TOTUFLUX(LAY) = 0.0
         DIRDOWN(LAY) = 0.0
         DIFDOWN(LAY) = 0.0
 200  CONTINUE

C *** Loop over frequency bands.
      DO 6000 IBAND = ISTART, IEND
        
         IB = IBAND

         IF (IBAND .EQ. 16) THEN
            CALL TAUGB16
         ELSEIF (IBAND .EQ. 17) THEN
            CALL TAUGB17
         ELSEIF (IBAND .EQ. 18) THEN
            CALL TAUGB18
         ELSEIF (IBAND .EQ. 19) THEN
            CALL TAUGB19
         ELSEIF (IBAND .EQ. 20) THEN
            CALL TAUGB20
         ELSEIF (IBAND .EQ. 21) THEN
            CALL TAUGB21
         ELSEIF (IBAND .EQ. 22) THEN
            CALL TAUGB22
         ELSEIF (IBAND .EQ. 23) THEN
            CALL TAUGB23
         ELSEIF (IBAND .EQ. 24) THEN
            CALL TAUGB24
         ELSEIF (IBAND .EQ. 25) THEN
            CALL TAUGB25
         ELSEIF (IBAND .EQ. 26) THEN
            CALL TAUGB26
         ELSEIF (IBAND .EQ. 27) THEN
            CALL TAUGB27
         ELSEIF (IBAND .EQ. 28) THEN
            CALL TAUGB28
         ELSEIF (IBAND .EQ. 29) THEN
            CALL TAUGB29
         ENDIF

c  set albedo for this band
         ALBEDO = 1. - SEMISS(IBAND)

C ***    Loop over g-channels.
         IG = 1
 1000    CONTINUE
C ***    Downward radiative transfer.
         SOL = SFLUXZEN(IG)
         FBEAM = ADJFLUX * SOL
         UMU0 = abs(1./SECZEN)
         DO 3900 LAY = NLAYERS, 1, -1
            FORWAER = PHASE(1,LAY,IBAND)**NSTR
            PHMULT = 1. / (1. - FORWAER)
            SSADMAER = SSAAER(LAY,IBAND) * (1. - FORWAER) /
     &           (1. - FORWAER * SSAAER(LAY,IBAND))
            TAUDMAER = TAUAER(LAY,IBAND) * (1. - FORWAER *
     &           SSAAER(LAY,IBAND))
            TAUREV(NLAYERS-LAY+1) = TAUG(LAY,IG)  + TAUDMAER +
     &           TAUCLOUD(LAY,IB)
            SCATAER = SSADMAER * TAUDMAER
            SCATRAY = SSA(LAY,IG) * TAUG(LAY,IG)
            SCATCLD = SSACLOUD(LAY,IB) * TAUCLOUD(LAY,IB)
            SSALB(NLAYERS-LAY+1) = (SCATAER + SCATRAY + SCATCLD)/
     &           TAUREV(NLAYERS-LAY+1)
            IF (SSALB(NLAYERS-LAY+1) .GT. 1.0) THEN
               PRINT*,'WARNING SSALB > 1.0, LAYER ',LAY,
     &              SSALB(NLAYERS-LAY+1)
           ENDIF
            PMOM(0,NLAYERS-LAY+1) = 1.
            PHAER1 = PHMULT * (PHASE(1,LAY,IBAND) - FORWAER)
            PMOM(1,NLAYERS-LAY+1) = (SCATAER * PHAER1 +
     &           SCATCLD*XMOM(1,LAY,IB))/ 
     &           (SCATAER + SCATRAY + SCATCLD) 
            PHAER2 = PHMULT * (PHASE(2,LAY,IBAND) - FORWAER)
            PMOM(2,NLAYERS-LAY+1) = (SCATAER * PHAER2 + 
     &           SCATRAY*PHASERAY(2) + SCATCLD*XMOM(2,LAY,IB))
     &           / (SCATAER + SCATRAY + SCATCLD)
            DO 3850 K = 3, NSTR
               PHAERK = PHMULT * (PHASE(K,LAY,IBAND) - FORWAER)
               PMOM(K,NLAYERS-LAY+1) = (SCATAER * PHAERK +
     &              SCATCLD*XMOM(K,LAY,IB))/ 
     &              (SCATAER + SCATRAY + SCATCLD)
 3850       CONTINUE
 3900    CONTINUE
         CALL DISORT( NLAYERS, TAUREV, SSALB, PMOM, TEMPER, WVNMLO,
     &                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
     &                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
     &                   FISOT,LAMBER, ALBEDO, HL, BTEMP, TTEMP, TEMIS,
     &                   DELTAM, sccos, PLANK, ONLYFL, ACCUR, PRNT, 
     &                   HEADER, MAXCLY, MAXULV, MAXUMU, MAXCMU, MAXPHI, 
     &                   RFLDIR, RFLDN, FLUP, fldir, fldn, dirscale,
     &                   DFDT, UAVG, UU, U0U, ALBMED, TRNMED )
         DO 3950 LEV = NLAYERS, 0, -1
C     Sum up either actual downward fluxes (IDELM=0) or delta-M scaled
C     downward fluxes.  Since the delta-M scaling is being handled 
C     outside of DISORT, the RFL and FL output variables are equal. 
C     (FYI: Had DELTAM been set to 1 in the call to DISORT, the FL
C     outputs would be the deltaM scaled fluxes and the RFL the 
C     unscaled fluxes.
            IF (IDELM .EQ. 0) THEN
               TOTDFLX = FLDIR(NLAYERS-LEV+1) + FLDN(NLAYERS-LEV+1)
               IF (LEV .EQ. NLAYERS) THEN
                  DIRFLUX = FBEAM * UMU0
               ELSE
                  TAUORIG = TAUG(LEV+1,IG) + 
     &                 TAUAER(LEV+1,IBAND) + TAUCLDORIG(LEV+1,IB)
                  TRANS = EXP(-TAUORIG/UMU0)
                  DIRFLUX = DIRFLUX * TRANS
               ENDIF
               DIRDOWN(LEV) = DIRDOWN(LEV) + DIRFLUX
               DIFDOWN(LEV) = DIFDOWN(LEV) + (TOTDFLX - DIRFLUX)
            ELSE
               DIRDOWN(LEV) = DIRDOWN(LEV) + FLDIR(NLAYERS-LEV+1)
               DIFDOWN(LEV) = DIFDOWN(LEV) + FLDN(NLAYERS-LEV+1)
            ENDIF
            TOTUFLUX(LEV) = TOTUFLUX(LEV) + FLUP(NLAYERS-LEV+1)
 3950    CONTINUE

         IG = IG + 1
         IF (IG .LE. NG(IBAND)) GO TO 1000
 6000 CONTINUE

      DIFDOWN(NLAYERS) = 0.
      IF (ISCCOS .EQ. 2) DIRDOWN(NLAYERS) = DIRDOWN(NLAYERS)/DIRSCALE
      TOTDFLUX(NLAYERS) = DIRDOWN(NLAYERS)
      FNET(NLAYERS) = TOTDFLUX(NLAYERS) - TOTUFLUX(NLAYERS)
      HTR(NLAYERS) = 0.
      DO 3951 LEV = NLAYERS-1, 0, -1
         IF (ISCCOS .EQ. 2) DIRDOWN(LEV) =  DIRDOWN(LEV)/DIRSCALE
         TOTDFLUX(LEV) =  DIRDOWN(LEV) + DIFDOWN(LEV)
         FNET(LEV) = TOTDFLUX(LEV) - TOTUFLUX(LEV)
         HTR(LEV) = -HEATFAC * (FNET(LEV) -FNET(LEV+1)) /
     &        (PZ(LEV) - PZ(LEV+1))
 3951 CONTINUE

      RETURN
      END   
