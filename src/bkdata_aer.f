C     path:      $Source$
C     author:    $Author$
C     revision:  $Revision$
C     created:   $Date$

	block data aerosol
	
	include 'param.f'
	PARAMETER (MCMU = 32)
	PARAMETER (NPTS = MXLAY * NBANDS)
	PARAMETER (NP = NPTS * MCMU)
	common /AERDAT/    ssaaer(mxlay,nbands), phase(mcmu,mxlay,nbands), 
	1    tauaer(mxlay,nbands)

      COMMON /HVERSN/    HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *                   HVDUM1(4),HVRUTL,HVREXT,
     *                   HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *                   HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

      CHARACTER*15 HVRRTM,HVRRTR,HVRATM,HVRSET,HVRTAU,
     *            HVDUM1,HVRUTL,HVREXT,
     *            HVRD1M,HVRR1M,HVREPK,HVRLPK,HVRAER,HVRBKA,
     *            HVRBKB,HVRCLD,HVRDIS,HVRLAM,HVRPAR

	data HVRBKA /'$Revision$'/
	data ssaaer /npts*0./
	data phase /np*0./
	data tauaer /npts*0./

	end
