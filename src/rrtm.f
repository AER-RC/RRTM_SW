C     path:      $Source: /storm/rc1/cvsroot/rc/rrtm_sw/src/rrtm.f,v $
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$
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

****************************************************************************
*                                                                          *
*                               RRTM_SW                                    *
*                                                                          *
*                                                                          *
*                                                                          *
*                   A RAPID RADIATIVE TRANSFER MODEL                       *
*                    FOR THE SOLAR SPECTRAL REGION                         *
*                                                                          *
*            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                  *
*                        131 HARTWELL AVE                                  *   
*                        LEXINGTON, MA 02421                               *
*                                                                          *
*                                                                          *
*                           ELI J. MLAWER                                  *   
*                         JENNIFER DELAMERE                                *
*                         STEVEN J. TAUBMAN~                               *
*                         SHEPARD A. CLOUGH                                *
*                                                                          *
*                                                                          *
*                         ~currently at GFDL                               *
*                                                                          *
*                                                                          *
*                                                                          *
*                       email:  mlawer@aer.com                             *
*                                                                          *
*        The authors wish to acknowledge the contributions of the          *
*        following people:  Patrick D. Brown, Michael J. Iacono,           *
*        Ronald E. Farren, Luke Chen, Robert Bergstrom.                    *
*                                                                          *
****************************************************************************

       PROGRAM RRTM_SW
                    
C *** This program is the driver for RRTM_SW, the AER rapid model.  
C     For each atmosphere the user wishes to analyze, this routine
C     a) calls READPROF to read in the atmospheric profile
C     b) calls SETCOEF to calculate various quantities needed for 
C        the radiative transfer algorithm
C     c) calls RTR or RTREG (depending on angular quadrature
C         method) to do the radiative transfer calculation
C     d) writes out the upward, downward, and net flux for each
C        level and the heating rate for each layer

      INCLUDE 	'param.f'

      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /FEATURES/  NG(IB1:IB2),NSPA(IB1:IB2),NSPB(IB1:IB2)
      COMMON /PRECISE/   ONEMINUS
      COMMON /CONTROL/   IAER, NSTR, IOUT, ISTART, IEND, ICLD,
     &                   idelm, isccos
      COMMON /SWPROP/    ZENITH, ALBEDO, ADJFLUX(NBANDS)
      COMMON /SURFACE/   IREFLECT,SEMISS(NBANDS)
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /IFIL/     IRD,IPR,IPU,IDUM(15)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),
     &     TAUCLOUD(MXLAY,NBANDS),SSACLOUD(MXLAY,NBANDS),
     &     xmom(0:16,MXLAY,NBANDS),taucldorig(mxlay,nbands)
      COMMON /OUTPUT/    TOTUFLUX(0:MXLAY,0:nbands), 
     &                   TOTDFLUX(0:MXLAY,0:nbands),
     &                   DIFDOWNFLUX(0:MXLAY,0:nbands), 
     &                   DIRDOWNFLUX(0:MXLAY,0:nbands),
     &                   FNET(0:MXLAY,0:nbands), 
     &                   HTR(0:MXLAY,0:nbands) 


      COMMON /CVRRTM/    HNAMRTM,HVRRTM
      COMMON /CVRSET/    HNAMSET,HVRSET
      COMMON /CVRATM/    HNAMATM,HVRATM
      COMMON /CVRUTL/    HNAMUTL,HVRUTL
      COMMON /CVRTAU/    HNAMTAU,HVRTAU
      COMMON /CVRCLD/    HNAMCLD,HVRCLD
      COMMON /CVREXT/    HNAMEXT,HVREXT
      COMMON /CVRRTR/    HNAMRTR,HVRRTR
      COMMON /CVRRDI/    HNAMRDI,HVRRDI
      COMMON /CVRERR/    HNAMERR,HVRERR
      COMMON /CVRLPK/    HNAMLPK,HVRLPK
      COMMON /CVRDIS/    HNAMDIS,HVRDIS

      COMMON /CVRSN16/   HNAMKG16,HVRKG16
      COMMON /CVRSN17/   HNAMKG17,HVRKG17
      COMMON /CVRSN18/   HNAMKG18,HVRKG18
      COMMON /CVRSN19/   HNAMKG19,HVRKG19
      COMMON /CVRSN20/   HNAMKG20,HVRKG20
      COMMON /CVRSN21/   HNAMKG21,HVRKG21
      COMMON /CVRSN22/   HNAMKG22,HVRKG22
      COMMON /CVRSN23/   HNAMKG23,HVRKG23
      COMMON /CVRSN24/   HNAMKG24,HVRKG24
      COMMON /CVRSN25/   HNAMKG25,HVRKG25
      COMMON /CVRSN26/   HNAMKG26,HVRKG26
      COMMON /CVRSN27/   HNAMKG27,HVRKG27
      COMMON /CVRSN28/   HNAMKG28,HVRKG28
      COMMON /CVRSN29/   HNAMKG29,HVRKG29

      CHARACTER*18 HVRRTM,HVRSET,HVRATM,HVRUTL,HVRTAU,
     *             HVRCLD,HVREXT,HVRRTR,HVRRDI,HVRERR,
     *             HVRLPK,HVRDIS

      CHARACTER*18 HNAMRTM,HNAMSET,HNAMATM,HNAMUTL,HNAMTAU,
     *             HNAMCLD,HNAMEXT,HNAMRTR,HNAMRDI,HNAMERR,
     *             HNAMLPK,HNAMDIS

      CHARACTER*18 HNAMKG16,HNAMKG17,HNAMKG18,HNAMKG19,
     *             HNAMKG20,HNAMKG21,HNAMKG22,HNAMKG23,
     *             HNAMKG24,HNAMKG25,HNAMKG26,HNAMKG27,
     *             HNAMKG28,HNAMKG29

      CHARACTER*18 HVRKG16,HVRKG17,HVRKG18,HVRKG19,
     *             HVRKG20,HVRKG21,HVRKG22,HVRKG23,
     *             HVRKG24,HVRKG25,HVRKG26,HVRKG27,
     *             HVRKG28,HVRKG29

      CHARACTER PAGE

      CHARACTER*50 OUTFORM(7)


c     Setup format statements for output

      DATA OUTFORM 
     1/'(1X,I3,3X,F7.6,4X,4(F10.4,4X),F11.6,4X,F10.5)',
     2 '(1X,I3,4X,F6.5,4X,4(F10.4,4X),F11.6,4X,F10.5)',
     3 '(1X,I3,4X,F6.4,4X,4(F10.4,4X),F11.6,4X,F10.5)',
     4 '(1X,I3,4X,F6.3,4X,4(F10.4,4X),F11.6,4X,F10.5)',
     5 '(1X,I3,4X,F6.2,4X,4(F10.4,4X),F11.6,4X,F10.5)',
     6 '(1X,I3,4X,F6.1,4X,4(F10.4,4X),F11.6,4X,F10.5)',
     7 '(1X,I3,4X,F6.1,4X,4(F10.4,4X),F11.6,4X,F10.5)'/

      PAGE = CHAR(12)

      HVRRTM = '$Revision$'      

      ONEMINUS = 1. - 1.E-6
c      PI = 2.*ASIN(1.)
      FLUXFAC = PI * 2.D4  

      IWR = 10

      
C     Multiple atmospheres not yet implemented. 
      NUMATMOS = 1
      DO 4000 IATMOS = 1, NUMATMOS

C ***    Input atmospheric profile from INPUT_RRTM.
         CALL READPROF
         
         IFLAG = IOUT

         IF (IFLAG .GT. 0 .AND. IFLAG .LE. IB2) THEN
            ISTART = IFLAG
            IEND = IFLAG
            IPBAND = IFLAG
         ELSE
            ISTART = IB1
            IEND = IB2
            IPBAND = 0
            IFLAG = IOUT
         ENDIF

C ***    Calculate information needed by the radiative transfer routine
C        that is specific to this atmosphere, especially some of the 
C        coefficients and indices needed to compute the optical depths
C        by interpolating data from stored reference atmospheres. 
         ICLDATM = 0
         IF (ICLD .EQ. 1) CALL CLDPROP(ICLDATM)

         CALL SETCOEF

C ***    Call the radiative transfer routine.
         CALL RTRDIS

         IF (IOUT .LT. 0) GO TO 4000

C ***    Process output for this atmosphere.
         OPEN (IWR,FILE='OUTPUT_RRTM',FORM='FORMATTED')

 1000    CONTINUE

         IF (IFLAG .GT. 0 .AND. IFLAG .LE. IB2) THEN
            ISTART = IFLAG
            IEND = IFLAG
            IPBAND = IFLAG
         ELSE
            IEND = IEND-1
            ISTART = IB2
            IPBAND = 0
         ENDIF

         if (isccos .eq. 1) then 
            write(iwr,9880) 
         elseif (isccos .eq. 2) then
            write(iwr,9881)
         else
            write(iwr,9879)
         endif

         if (idelm .eq. 0) then
            write(iwr,9883)
         else
            write(iwr,9882)
         endif

         WRITE(IWR,9899)WAVENUM1(ISTART),WAVENUM2(IEND)
         WRITE(IWR,9900)
         WRITE(IWR,9901)
C
         DO 3000 I = NLAYERS, 0, -1
            IF (PZ(I) .LT. 1.E-2) THEN
               INDFORM = 1
            ELSEIF (PZ(I) .LT. 1.E-1) THEN
               INDFORM = 2
            ELSEIF (PZ(I) .LT. 1.) THEN
               INDFORM = 3
            ELSEIF (PZ(I) .LT. 10.) THEN
               INDFORM = 4
            ELSEIF (PZ(I) .LT. 100.) THEN
               INDFORM = 5
            ELSEIF (PZ(I) .LT. 1000.) THEN
               INDFORM = 6
            ELSE
               INDFORM = 7
            ENDIF
            WRITE(IWR,OUTFORM(INDFORM)) I, PZ(I), 
     &           TOTUFLUX(I,IPBAND), DIFDOWNFLUX(I,IPBAND), 
     &           DIRDOWNFLUX(I,IPBAND), TOTDFLUX(I,IPBAND), 
     &           FNET(I,IPBAND), HTR(I,IPBAND)
 3000    CONTINUE
         WRITE(IWR,9903)PAGE
 3001    CONTINUE

         IF (IOUT .GE. 0 .AND. IOUT .LE. IB2) GO TO 3500
         IF (IFLAG .EQ. 98) THEN
            IFLAG = IB1
         ELSEIF (IFLAG .LT. IB2) THEN
            IFLAG = IFLAG + 1
         ELSE
            GO TO 3500
         ENDIF
         GO TO 1000
 3500    CONTINUE

C
C ***    Output module version numbers

         WRITE(IWR,9910) HNAMRTM, HVRRTM, HNAMATM, HVRATM,
     *                   HNAMSET, HVRSET, HNAMTAU, HVRTAU,
     *                   HNAMUTL, HVRUTL, HNAMCLD, HVRCLD,
     *                   HNAMRTR, HVRRTR, HNAMDIS, HVRDIS, 
     *                   HNAMRDI, HVRRDI, HNAMERR, HVRERR,
     *                   HNAMLPK, HVRLPK, HNAMEXT, HVREXT,
     *             HNAMKG16, HVRKG16, HNAMKG17, HVRKG17, 
     *             HNAMKG18, HVRKG18, HNAMKG19, HVRKG19, 
     *             HNAMKG20, HVRKG20, HNAMKG21, HVRKG21, 
     *             HNAMKG22, HVRKG22, HNAMKG23, HVRKG23, 
     *             HNAMKG24, HVRKG24, HNAMKG25, HVRKG25, 
     *             HNAMKG26, HVRKG26, HNAMKG27, HVRKG27, 
     *             HNAMKG28, HVRKG28, HNAMKG29, HVRKG29 

         CLOSE(IWR)

 4000 CONTINUE

 9879 format(1x)
 9880 format(1x,'All output fluxes adjusted to account for ins
     &trumental cosine response; heating rates invalid')
 9881 format(1x,'The output diffuse fluxes adjusted to account
     & for instrumental cosine response; heating rates invalid')

 9882 format(1x,'The downwelling direct and diffuse fluxes have been com
     &puted using the delta-M scaling approximation.') 
 9883 format(1x)

 9899 FORMAT(1X,'Wavenumbers: ',F6.0,' - ',F6.0,' cm-1')
 9900 FORMAT(1X,'LEVEL PRESSURE   UPWARD FLUX   DIFDOWN FLUX  DIRDOWN FL  
     &UX  DOWNWARD FLUX   NET FLUX    HEATING RATE')
 9901 FORMAT(1X,'         mb          W/m2          W/m2          W/m2
     &        W/m2          W/m2       degree/day')
 9902 FORMAT(1X,I3,3X,F11.6,4X,1P,2(G12.6,2X),G13.6,3X,G16.9,0P)
 9903 FORMAT(A)
 9910 FORMAT('  Modules and versions used in this calculation:',/,/,
     *         13(5X,a18,2X,A18,10X, a18,2X,A18,/))
      STOP
      END

C************************  SUBROUTINE READPROF  *****************************C

      SUBROUTINE READPROF                                                     
                                                                         
C     Read in atmospheric profile.

      IMPLICIT DOUBLE PRECISION (V)                                      
                                                                         
      INCLUDE 	'param.f'
      PARAMETER (MXMOL = 38)
      PARAMETER (MAXINPX=35)
      PARAMETER (MAXXSEC=4)
C      PARAMETER (MAXPROD = MXLAY*MAXXSEC)

      DIMENSION ALTZ(0:MXLAY),IXTRANS(14)
      DIMENSION SOLVAR(NBANDS)

      COMMON /CONTROL/  IAER, NSTR, IOUT, ISTART, IEND, ICLD,
     &                  idelm, isccos
      COMMON /CONSTANTS/FLUXFAC,HEATFAC
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
      COMMON /SWPROP/   ZENITH, ALBEDO, ADJFLUX(NBANDS)
      COMMON /SURFACE/  IREFLECT,SEMISS(NBANDS)
      COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /SPECIES/  COLDRY(MXLAY),WKL(MXMOL,MXLAY),WBRODL(MXLAY),
     &                  COLMOL(MXLAY),NMOL
      COMMON /IFIL/     IRD,IPR,IPU,IDUM(15)
      COMMON /XSECCTRL/ NXMOL,IXINDX(MAXINPX)
      COMMON /XSEC/     WX(MAXXSEC,MXLAY)
      COMMON /PATHX/    IXMAX,NXMOL0,IXINDX0(MAXINPX),WX0(MAXINPX,MXLAY)    
      COMMON /XRRTATM/  IXSECT

      CHARACTER*80 FORM1(0:1),FORM2(0:1),FORM3(0:1)
      CHARACTER*1 CTEST, CDOLLAR, CDUM

      DATA CDOLLAR /'$'/
      DATA IXTRANS /0,0,0,1,2,3,0,0,0,0,0,4,0,0/
C      DATA WX /MAXPROD*0.0/

      FORM1(0) = '(3F10.4,A3,I2,1X,2(F7.2,F8.3,F7.2))'
      FORM2(0) = '(3F10.4,A3,I2,23X,(F7.2,F8.3,F7.2))'
      FORM3(0) = '(8E10.3)'
      FORM1(1) = '(G15.7,G10.4,G10.4,A3,I2,1X,2(G7.2,G8.3,G7.2))'
      FORM2(1) = '(G15.7,G10.4,G10.4,A3,I2,23X,(G7.2,G8.3,G7.2))'
      FORM3(1) = '(8G15.7)'

C  Initialize molecular amount and cross section arrays to zero here.
      
      DO 1200 ILAY = 1,MXLAY
         DO 1100 ISP = 1,35
 1100       WKL(ISP,ILAY) = 0.0
         DO 1150 ISP = 1,MAXXSEC
 1150       WX(ISP,ILAY) = 0.0
 1200 CONTINUE

      IXMAX = MAXINPX
      IRD = 9
      OPEN (IRD,FILE='INPUT_RRTM',FORM='FORMATTED')

 1000 CONTINUE
      READ (IRD,9009,END=8800) CTEST
      IF (CTEST .NE. CDOLLAR) GO TO 1000
      READ (IRD,9011) IAER, IATM, ISCAT, ISTRM, IOUT, ICLD, IDELM, ICOS

      if (idelm.gt.1 .or. idelm.lt.0 .or. icos.gt.2 .or. icos.lt.0) then
         print *,'INVALID MEASUREMENT COMPARISON FLAG'
         stop
      endif
      isccos = icos

C     No cross-sections implemented in shortwave.
      IXSECT = 0

      IF (ISCAT .NE. 0) THEN
         PRINT *,' INVALID SCATTERING OPTION CHOSEN'
         STOP
      ENDIF

      IF (ISTRM .EQ. 0) THEN 
         NSTR = 4
      ELSE IF  (ISTRM .EQ. 1) THEN
         NSTR = 8
      ELSE IF  (ISTRM .EQ. 2) THEN
         NSTR = 16
      ELSE 
         PRINT *, 'INVALID VALUE FOR ISTRM'
         STOP
      ENDIF

      READ (IRD,9020) JULDAT, SZA, ISOLVAR, (SOLVAR(IB),IB=IB1,IB2)

      ZENITH = COS(SZA * PI / 180.)
      IF (JULDAT .EQ. 0) THEN
         ADJFLUX_JD = 1.
      ELSE
         ADJFLUX_JD = EARTH_SUN (JULDAT)
      ENDIF

C     If clouds are present, read in appropriate input file, IN_CLD_RRTM.
      IF (ICLD .EQ. 1) CALL READCLD

C     If aerosols are present, read in appropriate input file, IN_AER_RRTM.
      IF (IAER.EQ.10) CALL READAER 


      IF (ISOLVAR .EQ. 0) THEN
         DO 1400 IB = IB1,IB2
            ADJFLUX(IB) = ADJFLUX_JD
 1400    CONTINUE
      ELSEIF (ISOLVAR .EQ. 1) THEN
         DO 1450 IB=IB1,IB2
            ADJFLUX(IB) = ADJFLUX_JD * SOLVAR(IB1)
 1450    CONTINUE
      ELSEIF (ISOLVAR .EQ. 2) THEN
         DO 1475 IB=IB1,IB2
            ADJFLUX(IB) = ADJFLUX_JD * SOLVAR(IB)
 1475    CONTINUE
      ELSE
         PRINT *, 'ISOLVAR = ', ISOLVAR, ' NOT A VALID INPUT VALUE'
         STOP
      ENDIF

      READ (IRD,9012) IEMIS, IREFLECT, (SEMISS(IB),IB=IB1,IB2)
      IF (IEMIS .EQ. 0) THEN
         DO 1500 IB = IB1, IB2
            SEMISS(IB) = 1.
 1500    CONTINUE
      ELSEIF (IEMIS .EQ. 1) THEN
         DO 1600 IB = IB1, IB2
            SEMISS(IB) = SEMISS(IB1)
 1600    CONTINUE
      ELSEIF (IEMIS .EQ. 2) THEN
C          PRINT *, 'THESE ARE THE INPUT EMISSIVITY VALUES'
C          PRINT *, SEMISS(IB1:IB2)
      ELSE
          PRINT *, 'IEMIS = ', IEMIS, ' NOT A VALID INPUT VALUE'
          STOP
      ENDIF
     
      IF (IATM .EQ. 0) THEN
         READ (IRD,9013) IFORM,NLAYERS,NMOL
         IF (NMOL.EQ.0) NMOL = 7                                    
         READ (IRD,FORM1(IFORM)) PAVEL(1),TAVEL(1),SECNTK,CINP,
     &        IPTHAK,ALTZ(0),PZ(0),TZ(0),ALTZ(1),PZ(1),TZ(1)
         READ (IRD,FORM3(IFORM)) (WKL(M,1),M=1,7), WBRODL(1)
         IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,1),M=8,NMOL)

         DO 2000 L = 2, NLAYERS
            READ (IRD,FORM2(IFORM)) PAVEL(L),TAVEL(L),SECNTK,CINP,
     &           IPTHRK,ALTZ(L),PZ(L),TZ(L)
            READ (IRD,FORM3(IFORM)) (WKL(M,L),M=1,7), WBRODL(L)
            IF(NMOL .GT. 7) READ (IRD,FORM3(IFORM)) (WKL(M,L),M=8,NMOL)
 2000    CONTINUE                                                            
           
         IF (IXSECT .EQ. 1) THEN                                 
            READ (IRD,9300) NXMOL0
            NXMOL = NXMOL0
            CALL XSIDENT(IRD)
            READ (IRD,9301) IFORMX
C     
            DO 3000 L = 1, NLAYERS       
               READ (IRD,9010) CDUM
               READ (IRD, FORM3(IFORMX)) (WX0(M,L),M=1,7),WBRODX    
               IF (NXMOL0 .GT. 7) READ (IRD,FORM3(IFORMX)) 
     &              (WX0(M,L),M=8,NXMOL0)
 3000       CONTINUE
         ENDIF
      ELSE
         IPU = 7
         IPR = 66
         OPEN(UNIT=IPR,FILE='TAPE6',STATUS='UNKNOWN')
         CALL RRTATM
         IF (IXSECT .EQ. 1) THEN
            DO 3300 MX = 1, NXMOL0
               IXINDX(MX) = IXTRANS(IXINDX0(MX))
 3300       CONTINUE
         ENDIF
      ENDIF

C     Test for mixing ratio input.
      IMIX = 1
      DO 3500 M = 1, NMOL
         IF (WKL(M,1) .GT. 1.0) THEN
            IMIX = 0
            GO TO 3600
         ENDIF
 3500 CONTINUE
 3600 CONTINUE

      IF (IXSECT .EQ. 1) THEN
         IMIXX = 0
         IF (WX0(1,1) .LE. 1.0) IMIXX = 1
      ENDIF
      DO 5000 L = 1, NLAYERS
         SUMMOL = 0.0
         DO 4100 IMOL = 2, NMOL
            SUMMOL = SUMMOL + WKL(IMOL,L)
 4100    CONTINUE
         IF (IMIX .EQ. 1) THEN
            COLDRY(L) = WBRODL(L) / (1. - SUMMOL)
            DO 4200 IMOL = 1, NMOL
               WKL(IMOL,L) = COLDRY(L) * WKL(IMOL,L)
 4200       CONTINUE
         ELSE
            COLDRY(L) = WBRODL(L) + SUMMOL
         ENDIF
         IF (IXSECT .EQ. 1) THEN
            DO 4400 IX = 1, NXMOL0
               IF (IXINDX(IX) .NE. 0) THEN
                  IF (IMIXX .EQ. 1) THEN
                     WX(IXINDX(IX),L) = COLDRY(L) * WX0(IX,L) * 1.E-20
                  ELSE
                     WX(IXINDX(IX),L) = WX0(IX,L) * 1.E-20
                  ENDIF
               ENDIF
 4400       CONTINUE
         ENDIF
 5000 CONTINUE


      CLOSE(IRD)
      GO TO 9000

 8800 CONTINUE
      STOP ' INVALID INPUT_RRTM '

 9000 CONTINUE

 9009 FORMAT (A1,1X,I2,I2,I2)
 9010 FORMAT (A1)
 9011 FORMAT (18X,I2,29X,I1,32X,I1,1X,I1,2X,I3,4X,I1,3x,i1,i1)
 9012 FORMAT (11X,I1,2X,I1,14F5.3)
 9013 FORMAT (1X,I1,I3,I5)                                     
 9020 format (12X, I3, 3X, F7.4, 4X, I1,14F7.5)
 9300 FORMAT (I5)
 9301 FORMAT (1X,I1)

      RETURN
      END 

C************************  SUBROUTINE READCLD  *****************************C

      SUBROUTINE READCLD

C     Purpose:  To read in IN_CLD_RRTM_SW, the file that contains input 
C               cloud properties.

      INCLUDE 	'param.f'
      COMMON /CONTROL/   IAER, NSTR, IOUT, ISTART, IEND, ICLD,
     &                   idelm, isccos
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /CLOUDIN/   INFLAG,CLDDAT1(MXLAY),CLDDAT2(MXLAY),
     &                   ICEFLAG,LIQFLAG,CLDDAT3(MXLAY),CLDDAT4(MXLAY),
     &                   CLDDATMOM(0:16,MXLAY)
      COMMON /CLOUDDAT/  NCBANDS,CLDFRAC(MXLAY),
     &     TAUCLOUD(MXLAY,NBANDS),SSACLOUD(MXLAY,NBANDS),
     &     xmom(0:16,MXLAY,NBANDS)
      CHARACTER*1 CTEST, CPERCENT

      DATA CPERCENT /'%'/
      IRDCLD = 11

      OPEN(IRDCLD,FILE='IN_CLD_RRTM',FORM='FORMATTED')

C     Read in cloud input option.  
      READ(IRDCLD,9050) INFLAG, ICEFLAG, LIQFLAG

      DO 500 LAY = 1, NLAYERS
         CLDFRAC(LAY) = 0.
 500  CONTINUE

      IF (INFLAG .EQ. 0) THEN
 950     CONTINUE
c     For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
c     DAT1 are pertinent.  If CTEST = '%', then there are no more 
C     cloudy layers to process.
         READ (IRDCLD,9099,END=8950) CTEST,LAY,FRAC,
     &        DAT1,DAT2,(CLDDATMOM(ISTR,LAY),ISTR=0,NSTR)
         IF (CTEST .EQ. CPERCENT) GO TO 8950
         CLDFRAC(LAY) = FRAC
         CLDDAT1(LAY) = DAT1
         CLDDAT2(LAY) = DAT2
         GO TO 950
 8950    CONTINUE
      ELSE
 1000    CONTINUE
c     For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
c     DAT1 are pertinent.  If CTEST = '%', then there are no more 
C     cloudy layers to process.
         READ (IRDCLD,9100,END=9000) CTEST,LAY,FRAC,
     &        DAT1,DAT2,DAT3,DAT4
         IF (CTEST .EQ. CPERCENT) GO TO 9000
         CLDFRAC(LAY) = FRAC
         CLDDAT1(LAY) = DAT1
         CLDDAT2(LAY) = DAT2
         CLDDAT3(LAY) = DAT3
         CLDDAT4(LAY) = DAT4
         GO TO 1000
 9000    CONTINUE
      ENDIF

      CLOSE(IRDCLD)

 9050 FORMAT (3X,I2,4X,I1,4X,I1)
 9099 FORMAT (A1,1X,I3,19E10.5)
 9100 FORMAT (A1,1X,I3,5E10.5)
      RETURN
      END

C************************  SUBROUTINE READAER  *****************************C

      SUBROUTINE READAER 

C     Purpose:  To read in IN_AER_RRTM, the file that contains input
C               aerosol properties.

      INCLUDE 	'param.f'
      PARAMETER (MCMU = 32)
      real aerpar(3), ssa(nbands), asym(nbands), aod(mxlay),aod1(nbands)
      real rlambda(nbands), specfac(nbands)
      real rnu0(16:29),rnu1(23:26)
      real f1(23:26),od0(23:26),od1(23:26)
      integer lay(mxlay),ivec(mxlay)

      COMMON /CONTROL/   IAER, NSTR, IOUT, ISTART, IEND, ICLD,
     &                   idelm, isccos
      COMMON /PROFILE/   NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND
      COMMON /SWPROP/    ZENITH, ALBEDO, ADJFLUX(NBANDS)
      common /AERDAT/    ssaaer(mxlay,nbands), phase(mcmu,mxlay,nbands), 
     &                   tauaer(mxlay,nbands), phase_in(mcmu,nbands)
      CHARACTER*1 CTEST, CPERCENT

      DATA CPERCENT /'%'/

      data rnu0 /2903.,3601.,4310.,4892.,5623.,6872.,7872.,
     &               10590.,14420.,18970.,25015.,30390.,43507.,1412./

       data rnu1 /10530.7,14293.3,18678.0,24475.1/

       data f1/0.9929,0.9883,0.978,0.9696/
       data od0/0.1084,0.167,0.245,0.3611/
       data od1 /0.3144,0.4822,0.7013,1.0239/

      eps = 1.e-10
      IRDAER = 12
      OPEN(IRDAER,FILE='IN_AER_RRTM',FORM='FORMATTED')

      aod(:) = 0.0
      tauaer(:,ib1:ib2) = 0.0

C     Read in number of different aerosol models option.
      read (irdaer, 9010) naer
c      if (naer .gt. 4) then
c         print *, 'NAER (= ', naer, ') IS GREATER THAN 4'
c         stop
c      endif

c     For each aerosol read in optical properties and layer aerosol 
c     optical depths.
      do ia = 1, naer
	 read (irdaer, 9011) nlay, iaod, issa, iasym, (aerpar(i),i=1,3)

         if (iaod .eq. 0) then
c           Set defaults to get standard Angstrom relation.
            if (aerpar(2) .lt. eps) aerpar(2) = 1.

            do ib = ib1, ib2
	       rlambda(ib)=10000./rnu0(ib)
               specfac(ib) = (aerpar(2) + aerpar(3) * rlambda(ib)) /
     &              ((aerpar(2) + aerpar(3) - 1.) + 
     &              rlambda(ib)**aerpar(1))
            enddo
         endif

C        For this aerosol, read in layers and optical depth information.
C        Store a nonzero optical depth in aod to check for double
C        specification.
         do il = 1, nlay
            read(irdaer, 9012) lay(il), (aod1(ib), ib = ib1,ib2)
            if (aod(lay(il)) .lt. eps) then
               if (iaod .eq. 0) then
                  aod(lay(il)) = aod1(ib1)
                  do ib = ib1, ib2
                     tauaer(lay(il),ib) = aod(lay(il)) * specfac(ib)
                  enddo
               else
                  do ib = ib1, ib2
                     aod(lay(il)) = max(aod(lay(il)),aod1(ib))
                     tauaer(lay(il),ib) = aod1(ib)
                  enddo
               endif
            else
               print *,'LAYER ',lay(il),' HAS MORE THAN 
     &              ONE AEROSOL TYPE'
               stop
            endif
         enddo

c      Build vector of aerosol layer indices 

       do il=1,nlay
          ivec(il) = lay(il) 
       end do

c      Correct bands 23 through 26 for sza effect (negligible for others)
         do ib=23,26
            if (iaod.eq.0) then
                od = sum(tauaer(ivec(1:nlay),ib))/zenith
                rnu = rnu0(ib)+
     &           (rnu1(ib)-rnu0(ib))*(od-od0(ib))/(od1(ib)-od0(ib))
               rlambda_new=10000./rnu
               specfac_new = (aerpar(2)+aerpar(3)*rlambda_new) /
     &          ((aerpar(2)+aerpar(3)- 1.)+rlambda_new**aerpar(1))
               do il=1,nlay
                  tauaer(lay(il),ib) = tauaer(lay(il),ib)*
     &			specfac_new/specfac(ib)
               end do
            endif
         end do

c        For this aerosol, read and store optical properties
         read (irdaer, 9013) (ssa(ib), ib = ib1,ib2)

         DO 2000 IB = IB1, IB2
            do il = 1, nlay
               if (issa .eq. 0) then 
                  ssaaer(lay(il),IB) = ssa(ib1)
               else
                  ssaaer(lay(il),IB) = ssa(IB)
               endif
            enddo
 2000    CONTINUE

         if (iasym .lt. 2) then
            read (irdaer, 9013) (asym(ib), ib = ib1,ib2)

            DO 3000 IB = IB1, IB2
               do il = 1, nlay
                  do istr = 1,  nstr
                     if (iasym .eq. 0) then 
                        phase(istr,lay(il),IB) = asym(ib1)**istr
                     elseif (iasym .eq. 1) then
                        phase(istr,lay(il),IB) = asym(IB)**istr
                     endif
                  enddo
               enddo
 3000       CONTINUE
         else
	    do istr = 1, nstr
	       read (irdaer, 9013) (phase_in(istr,ib), 
     &                 ib = ib1,ib2)
	    enddo
            do il = 1, nlay
               do istr = 1, nstr
                  phase(istr,lay(il),ib) = phase_in(istr,ib)
               enddo
            enddo
         endif

      enddo    ! end of naer loop


 9000 CONTINUE
      CLOSE(IRDAER)

 9010 format (3x, i2)
 9011 format (2x, i3, 4x, i1, 4x, i1, 4x, i1, 3f8.3)
 9012 format (2x, i3, 14f7.4)
 9013 format (14f5.3)

      RETURN
      END


C************************  SUBROUTINE XSIDENT  *****************************C

      SUBROUTINE XSIDENT(IRD)
C                                                                         
C     This subroutine identifies which cross-sections are to be used.

      PARAMETER (MAXINPX=35)
      PARAMETER (MAXXSEC=4)

      IMPLICIT DOUBLE PRECISION (V)                                     ! 
C                                                                         
      COMMON /XSECCTRL/ NXMOL,IXINDX(MAXINPX)
C                                                                         
C     NXMOL     - number of cross-sections input by user
C     IXINDX(I) - index of cross-section molecule corresponding to Ith
C                 cross-section specified by user
C                 = 0 -- not allowed in RRTM
C                 = 1 -- CCL4
C                 = 2 -- CFC11
C                 = 3 -- CFC12
C                 = 4 -- CFC22
C                                                                         
C     XSNAME=NAMES, ALIAS=ALIASES OF THE CROSS-SECTION MOLECULES          
C                                                                         
      CHARACTER*10 XSNAME(MAXINPX),ALIAS(MAXXSEC,4),BLANK               
C                                                                         
      DATA (ALIAS(1,I),I=1,4)/                                           
     *    'CCL4      ', 'CCL3F     ', 'CCL2F2    ', 'CHCLF2    '/ 
      DATA (ALIAS(2,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'CFCL3     ', 'CF2CL2    ', 'CHF2CL    '/         
      DATA (ALIAS(3,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'CFC11     ', 'CFC12     ', 'CFC22     '/         
      DATA (ALIAS(4,I),I=1,4)/                                           
     *    ' ZZZZZZZZ ', 'F11       ', 'F12       ', 'F22       '/        

      DATA BLANK / '          '/                                          
C                                                                         
      DO 10 I = 1, NXMOL                                                 
         XSNAME(I) = BLANK                                                
   10 CONTINUE                                                            
C                                                                         
C     READ IN THE NAMES OF THE MOLECULES                                  
C                                                                         
      IF (NXMOL.GT.7) THEN                                               
         READ (IRD,'(7A10)') (XSNAME(I),I=1,7)                            
         READ (IRD,'(8A10)') (XSNAME(I),I=8,NXMOL)                       
      ELSE                                                                
         READ (IRD,'(7A10)') (XSNAME(I),I=1,NXMOL)                       
      ENDIF                                                               
C                                                                         
C     MATCH THE NAMES READ IN AGAINST THE NAMES STORED IN ALIAS           
C     AND DETERMINE THE INDEX VALUE.  
      IXMAX = 4                                                          
      DO 40 I = 1, NXMOL                                                 
C        Left-justify all inputed names.                                      
         CALL CLJUST (XSNAME(I),10)
         IXINDX(I) = 0
         DO 20 J = 1, IXMAX
            IF ((XSNAME(I).EQ.ALIAS(1,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(2,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(3,J)) .OR.                            
     *          (XSNAME(I).EQ.ALIAS(4,J))) THEN                           
               IXINDX(I) = J                                              
            ENDIF                                                         
   20    CONTINUE
   40 CONTINUE                                                            

      RETURN
      END

*****************************************************************
	real function earth_sun (idn)

C   function to calculate  correction factor of the Earth's orbit 

C     dn	 	: Julian day
C     earth_sun_ratio 	: square of the ratio of mean to actual Earth-Sun distance    


	PI   = 	3.141592654
	gamma = 2.*PI*(idn-1)/365. 

c   use Iqbal's equation 1.2.1 

	earth_sun= 1.000110 + .034221 * cos(gamma) + .001289*sin(gamma)+
     1                 .000719 *cos(2.*gamma) + .000077 * sin(2.*gamma)
     
	return 
	end

*****************************************************************
      BLOCK DATA

      include 'param.f'

      COMMON /CVRRTM/    HNAMRTM,HVRRTM
      COMMON /CVRSET/    HNAMSET,HVRSET
      COMMON /CVRATM/    HNAMATM,HVRATM
      COMMON /CVRUTL/    HNAMUTL,HVRUTL
      COMMON /CVRTAU/    HNAMTAU,HVRTAU
      COMMON /CVRCLD/    HNAMCLD,HVRCLD
      COMMON /CVREXT/    HNAMEXT,HVREXT
      COMMON /CVRRTR/    HNAMRTR,HVRRTR
      COMMON /CVRRDI/    HNAMRDI,HVRRDI
      COMMON /CVRERR/    HNAMERR,HVRERR
      COMMON /CVRLPK/    HNAMLPK,HVRLPK
      COMMON /CVRDIS/    HNAMDIS,HVRDIS

      COMMON /CVRSN16/   HNAMKG16,HVRKG16
      COMMON /CVRSN17/   HNAMKG17,HVRKG17
      COMMON /CVRSN18/   HNAMKG18,HVRKG18
      COMMON /CVRSN19/   HNAMKG19,HVRKG19
      COMMON /CVRSN20/   HNAMKG20,HVRKG20
      COMMON /CVRSN21/   HNAMKG21,HVRKG21
      COMMON /CVRSN22/   HNAMKG22,HVRKG22
      COMMON /CVRSN23/   HNAMKG23,HVRKG23
      COMMON /CVRSN24/   HNAMKG24,HVRKG24
      COMMON /CVRSN25/   HNAMKG25,HVRKG25
      COMMON /CVRSN26/   HNAMKG26,HVRKG26
      COMMON /CVRSN27/   HNAMKG27,HVRKG27
      COMMON /CVRSN28/   HNAMKG28,HVRKG28
      COMMON /CVRSN29/   HNAMKG29,HVRKG29

      COMMON /CONSTANTS/ FLUXFAC,HEATFAC
      COMMON /FEATURES/  NG(IB1:IB2),NSPA(IB1:IB2),NSPB(IB1:IB2)

      CHARACTER*18 HVRRTM,HVRSET,HVRATM,HVRUTL,HVRTAU,
     *             HVRCLD,HVREXT,HVRRTR,HVRRDI,HVRERR,
     *             HVRLPK,HVRDIS

      CHARACTER*18 HNAMRTM,HNAMSET,HNAMATM,HNAMUTL,HNAMTAU,
     *             HNAMCLD,HNAMEXT,HNAMRTR,HNAMRDI,HNAMERR,
     *             HNAMLPK,HNAMDIS

      CHARACTER*18 HNAMKG26

      CHARACTER*18 HVRKG26


      DATA HVRRTM / 'NOT USED' /,  
     *     HVRRTR / 'NOT USED' /,   HVRATM / 'NOT USED' /,
     *     HVRSET / 'NOT USED' /,   HVRTAU / 'NOT USED' /,
     *     HVRUTL / 'NOT USED' /,
     *     HVREXT / 'NOT USED' /, 
     *     HVRRDI / 'NOT USED' /,
     *     HVRERR / 'NOT USED' /,   HVRLPK / 'NOT USED' /,
     *     HVRCLD / 'NOT USED' /,
     *     HVRDIS / 'NOT USED' / 
 
      DATA HNAMRTM / '           rrtm.f:' /,
     *     HNAMSET / '        setcoef.f:' /,
     *     HNAMATM / '         rrtatm.f:' /,
     *     HNAMUTL / '       util_xxx.f:' /,
     *     HNAMTAU / '      taumoldis.f:' /,
     *     HNAMCLD / '        cldprop.f:' /,
     *     HNAMEXT / '          extra.f:' /,
     *     HNAMRTR / '         rtrdis.f:' /,
     *     HNAMRDI / '       RDI1MACH.f:' /,
     *     HNAMERR / '        ErrPack.f:' /,
     *     HNAMLPK / '         LINPAK.f:' /,
     *     HNAMDIS / '         disort.f:' /

      DATA HVRKG26 / ' ' /
      DATA HNAMKG26 / ' '/

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


      END
c**********************************************************************
      Block Data phys_consts
c
      COMMON /CONSTS/ PI,PLANCK,BOLTZ,CLIGHT,AVOGAD,ALOSMT,GASCON,
     *                RADCN1,RADCN2 
c
      DATA PI / 3.14159265 /
c
c    Constants from NIST 01/11/2002
c
      DATA PLANCK / 6.62606876E-27 /, BOLTZ  / 1.3806503E-16 /,
     *     CLIGHT / 2.99792458E+10 /, 
     *     AVOGAD / 6.02214199E+23 /, ALOSMT / 2.6867775E+19 /,
     *     GASCON / 8.314472  E+07 /
     *     RADCN1 / 1.191042722E-12 /, RADCN2 / 1.4387752    /
c
c     Pi was obtained from   PI = 2.*ASIN(1.)                             A03980
c
c     units are genrally cgs
c
c     The first and second radiation constants are taken from NIST.
c     They were previously obtained from the relations:
c                            RADCN1 = 2.*PLANCK*CLIGHT*CLIGHT*1.E-07      A03990
c                            RADCN2 = PLANCK*CLIGHT/BOLTZ                 A04000
      end
c

