	block data aerosol
	
	include 'param.f'
	PARAMETER (MCMU = 32)
	PARAMETER (NPTS = MXLAY * NBANDS)
	PARAMETER (NP = NPTS * MCMU)
	common /AERDAT/    ssaaer(mxlay,nbands), phase(mcmu,mxlay,nbands), 
	1    tauaer(mxlay,nbands)

	
	data ssaaer /npts*0./
	data phase /np*0./
	data tauaer /npts*0./

	end
