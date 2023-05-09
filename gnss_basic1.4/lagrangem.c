/*=================================================================
 *
 * MEX file for the lagrange interpolation.  Takes in an MX32 array of
 * GPS Satellite positions (ECEF), at an MX1 array of times.  Interplates the
 * position of the satellite in the correct column of gpsx, according
 * to the PRN given in the NX1 array prn, at the time given in the MX1
 * array datatime.
 *
 *
 *=================================================================*/

extern void lagrangec(
 		   double	gps_ecef_component[],
 		   double	gpstime[],
                   double       prnset[],
                   double       interptime[],
                   double       *interp_gps_ecef_component,
                   int          ninterp_pts,
                   int          ngps,
                   double       *test1
		   );

#include <math.h>
#include "mex.h"

/* Input Arguments - defn. of Right hand pointers */

#define	GPS_ECEF_IN	prhs[0]
#define	GPS_TIME_IN	prhs[1]
#define	PRNSET_IN	prhs[2]
#define GPS_INTERP_TIME_IN prhs[3]

/* Output Arguments - defn. of Left hand pointers */

#define	GPS_ECEF_INTERP_OUT	plhs[0]
#define	TEST1_OUT	plhs[1]

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *gpsx_inp, *gpsi_outp, *gpstime_inp, *gpsinterptime_inp, *prn_inp;
    double *test1;
    unsigned int ninput, noutput;
    
    ninput = mxGetM(GPS_TIME_IN); 
    noutput = mxGetM(GPS_INTERP_TIME_IN); 

    /* Create a matrix for the return argument */ 
 

    GPS_ECEF_INTERP_OUT = mxCreateDoubleMatrix(noutput, 1, mxREAL); 
    TEST1_OUT = mxCreateDoubleMatrix(noutput, 1, mxREAL); 
    
    
    /* Assign pointers to the various parameters */ 

    gpsx_inp = mxGetPr(GPS_ECEF_IN);
    gpsi_outp = mxGetPr(GPS_ECEF_INTERP_OUT);
    gpstime_inp = mxGetPr(GPS_TIME_IN);
    gpsinterptime_inp = mxGetPr(GPS_INTERP_TIME_IN);
    prn_inp = mxGetPr(PRNSET_IN);
    test1 = mxGetPr(TEST1_OUT);

    /* Do the actual computations in the subroutine lagrangec */

    lagrangec( gpsx_inp, gpstime_inp, prn_inp, gpsinterptime_inp, 
               gpsi_outp, noutput, ninput, test1);   

    return;
    
}


