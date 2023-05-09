/*
 * MATLAB Compiler: 2.0.1
 * Date: Fri Mar  1 11:27:28 2002
 * Arguments: "-x" "dp_dtheta_dopp_sum" 
 */

#ifndef MLF_V2
#define MLF_V2 1
#endif

#ifndef __dp_dtheta_dopp_sum_h
#define __dp_dtheta_dopp_sum_h 1

#include "matlab.h"

extern mxArray * mlfDp_dtheta_dopp_sum(mxArray * * tan_betax,
                                       mxArray * * tan_betay,
                                       mxArray * * sinc_lobes,
                                       mxArray * theta_input,
                                       mxArray * delta_input,
                                       mxArray * alpha,
                                       mxArray * PDF_params,
                                       mxArray * fc);
extern void mlxDp_dtheta_dopp_sum(int nlhs,
                                  mxArray * plhs[],
                                  int nrhs,
                                  mxArray * prhs[]);

#endif
