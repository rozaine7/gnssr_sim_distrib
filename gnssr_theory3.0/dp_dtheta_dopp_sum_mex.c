/*=================================================================
 *
 * MEX file interface for dp/dtheta for the computation of ocean scattering
 * Derived from the MATLAB code to make the simulations run faster.
 * Kept general, however  the PDF description must be explicitly modified, 
 * 
 * Calls dp_dtheta_dopp_sumf.c - which actually does the computation
 * Requires function gram_charlierf.c - G/C PDF model.
 * 
 * MATLAB global variables, must be entered through the calling line:
 *   model_flag, Ti, VR, VG, alt, gamma, RE, RGPS
 *
 *=================================================================*/

extern void gram_charlierf(double tau_betau[],
                           double tau_betac[],
                           double PDF_params[],
                           double *PDF );


extern void dp_dtheta_dopp_sumf(
 		   double	theta_input[],
 		   double	delta_input[],
                   double       *ttilde_input,
                   int          mdelta,
                   int          ntheta,
                   double       *alpha,
                   double       PDF_params[],
                   double       *Ti,
                   double       *alt,
                   double       *gamma,
                   double       *RE,
                   double       *RGPS,
                   double       *fc,
                   double       VR[],
                   double       VG[],
  		   double	*dp_dopp_mat_real,
         	   double	*dp_dopp_mat_imag,
		   double       *tbx,
		   double       *tby,
                   double       *xell,
                   double       *yell,
                   double       *pdf_val,
                   double       *ftilde_step,
                   double       *ftilde_input,
                   double       correlation_spectrum_real[],
                   double       correlation_spectrum_imag[],
                   double       *wavelegth
		   );

#include <math.h>
#include "mex.h"

/* Input Arguments - defn. of Right hand pointers */

#define	THETA_IN	prhs[0]
#define	DELTA_IN	prhs[1]
#define	ALPHA_IN	prhs[2]
#define	PDF_IN  	prhs[3]
#define	TI_IN   	prhs[4]
#define	ALT_IN  	prhs[5]
#define	GAMMA_IN  	prhs[6]
#define	RE_IN   	prhs[7]
#define	RGPS_IN  	prhs[8]
#define FC_IN           prhs[9]
#define VR_IN           prhs[10]
#define VG_IN           prhs[11]
#define TTILDE_IN       prhs[12]
#define FTILDE_IN       prhs[13]
#define DFTILDE_IN       prhs[14]
#define LAMBDA_IN       prhs[15]

/* Output Arguments - defn. of Left hand pointers */

#define	DPREAL_OUT	plhs[0]
#define	DPIMAG_OUT	plhs[1]
#define	TX_OUT	plhs[2]
#define	TY_OUT	plhs[3] 
#define	X_OUT	plhs[4] 
#define	Y_OUT	plhs[5] 
#define	PDF_OUT	plhs[6] 
#define	SPREAL_OUT	plhs[7]
#define	SPIMAG_OUT	plhs[8]


#define PI 3.14159265

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
  double *theta_inp, *delta_inp, *ttilde_inp, *alpha_inp, *pdf_inp, *ti_inp;
    double *alt_inp, *gamma_inp, *re_inp, *rgps_inp;
    double *fc_inp, *vg_inp, *vr_inp;
    double *dpreal_outp, *dpimag_outp;
    double *tx_outp, *ty_outp, *x_outp, *y_outp, *pdf_outp;
    double *ftilde_inp, *dftilde_inp, *sreal_outp, *simag_outp;
    double *wavelength_inp;
    unsigned int ntheta, mdelta;
    int k, i; 
    
    /* Check for proper number of arguments */
    
    /*    if (nrhs != 2) { 

    } else if (nlhs > 1) {
	mexErrMsgTxt("Too many output arguments."); 
	} */

    
    /* Check the dimensions of Y.  Y can be 4 X 1 or 1 X 4. */ 
    
    ntheta = mxGetN(THETA_IN); 
    mdelta = mxGetM(DELTA_IN); 

    /*    printf("md = %5i, nd=%5i", mdelta, ntheta);  */

    /*    if (!mxIsDouble(Y_IN) || mxIsComplex(Y_IN) || 
	(MAX(m,n) != 4) || (MIN(m,n) != 1)) { 
	mexErrMsgTxt("YPRIME requires that Y be a 4 x 1 vector."); 
	} */
    
    /* Create a matrix for the return argument */ 

    DPREAL_OUT = mxCreateDoubleMatrix(mdelta, ntheta, mxREAL); 
    DPIMAG_OUT = mxCreateDoubleMatrix(mdelta, ntheta, mxREAL); 
    TX_OUT = mxCreateDoubleMatrix(mdelta, ntheta, mxREAL); 
    TY_OUT = mxCreateDoubleMatrix(mdelta, ntheta, mxREAL); 
    X_OUT = mxCreateDoubleMatrix(mdelta, ntheta, mxREAL); 
    Y_OUT = mxCreateDoubleMatrix(mdelta, ntheta, mxREAL); 
    PDF_OUT = mxCreateDoubleMatrix(mdelta, ntheta, mxREAL); 
    SPREAL_OUT = mxCreateDoubleMatrix(mdelta, ntheta, mxREAL); 
    SPIMAG_OUT = mxCreateDoubleMatrix(mdelta, ntheta, mxREAL); 

    
    /* Assign pointers to the various parameters */ 

    theta_inp = mxGetPr(THETA_IN);
    delta_inp = mxGetPr(DELTA_IN);
    ttilde_inp = mxGetPr(TTILDE_IN);
    ftilde_inp = mxGetPr(FTILDE_IN);
    dftilde_inp = mxGetPr(DFTILDE_IN);
    alpha_inp = mxGetPr(ALPHA_IN);
    pdf_inp = mxGetPr(PDF_IN);
    ti_inp  = mxGetPr(TI_IN);
    alt_inp = mxGetPr(ALT_IN);
    gamma_inp = mxGetPr(GAMMA_IN);
    re_inp = mxGetPr(RE_IN); 
    rgps_inp = mxGetPr(RGPS_IN);
    fc_inp = mxGetPr(FC_IN);
    vr_inp = mxGetPr(VR_IN);
    vg_inp = mxGetPr(VG_IN);
    dpreal_outp = mxGetPr(DPREAL_OUT);
    dpimag_outp = mxGetPr(DPIMAG_OUT);
    tx_outp = mxGetPr(TX_OUT);
    ty_outp = mxGetPr(TY_OUT); 
    x_outp = mxGetPr(X_OUT); 
    y_outp = mxGetPr(Y_OUT); 
    pdf_outp = mxGetPr(PDF_OUT); 
    sreal_outp = mxGetPr(SPREAL_OUT); 
    simag_outp = mxGetPr(SPIMAG_OUT); 
    wavelength_inp = mxGetPr(LAMBDA_IN); 

    /* Do the actual computations in the subroutine dp_dtheta_dopp_sumf */

    dp_dtheta_dopp_sumf( theta_inp,  delta_inp, ttilde_inp, mdelta, ntheta,  
                         alpha_inp, pdf_inp, ti_inp, alt_inp, 
                         gamma_inp, re_inp, rgps_inp, fc_inp,
                         vr_inp, vg_inp, dpreal_outp, dpimag_outp, 
                         tx_outp, ty_outp,
                         x_outp, y_outp, pdf_outp, dftilde_inp, ftilde_inp,
                         sreal_outp, simag_outp, wavelength_inp);

    return;
    
}


