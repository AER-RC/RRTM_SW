C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

	parameter (mxlay = 403, nbands = 29)
	parameter (ib1 = 16, ib2 = 29)
        parameter (mg = 16)
	parameter (mxstr = 16)

      COMMON /BANDS/     WAVENUM1(IB1:IB2),
     &                   WAVENUM2(IB1:IB2),
     &                   DELWAVE(IB1:IB2)
