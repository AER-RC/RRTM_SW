	subroutine get_lambda 
c
c  this routine calculates the average AOD assuming the Angstrom relation 
c  for aerosol optical depth; this average AOD is then used to calculate 
c  the corresponding wavelength, which will be used by aerprop to calculate 
c  the AOD for the band from the full Molineaux expression


	include 'param.f'
        real  rl1(nbands), rl2(nbands)

        common /AERIN/   aod(mxlay), aerpar(3,mxlay), 
     &                   rlambda(mxlay, nbands)
      	COMMON /PROFILE/ NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND


	do ib = ib1, ib2
	    rl1(ib) = 10000. / wavenum1(ib)	
	    rl2(ib) = 10000. / wavenum2(ib)	
	enddo

        do il = 1, nlayers
           if (aod(il) .gt. 0.) then
	      omaer  = 1. - aerpar(1,il)
	      factor = 1. / omaer
	      do ib = ib1, ib2
		 aodbar = factor *(rl2(ib)**omaer - rl1(ib)**omaer)/ 
     1                  (rl2(ib) - rl1(ib))
		 rlambda(il,ib) = (1. / aodbar) ** (1./ aerpar(1,il))
	      enddo
           endif
        enddo

	return
	end
