	subroutine aerprop


c   this routine calculates aerosol optical depth (AOD) from the full Molineaux formulation

	include 'param.f'

      	common /AERDAT/  ssaaer(mxlay,nbands), gaer(mxlay,nbands), 
     &                   tauaer(mxlay,nbands)
        common /AERIN/   aod(mxlay), aerpar(3,mxlay), 
     &                   rlambda(mxlay,nbands)
      	COMMON /PROFILE/ NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),
     &                   PZ(0:MXLAY),TZ(0:MXLAY),TBOUND

c  get best wavelength for each band for AOD calculations

        call get_lambda 

      	do il = 1,nlayers
	  if (aod(il) .gt. 0.) then
	     do ib = ib1, ib2
                tauaer(il,ib) = aod(il) *
     1               ( aerpar(2,il) + aerpar(3,il) * rlambda(il,ib)) /
     2               ((aerpar(2,il) + aerpar(3,il) - 1.) + 
     3               rlambda(il,ib)**aerpar(1,il))
	      enddo
	   endif
      	enddo

      return
      end
	    
