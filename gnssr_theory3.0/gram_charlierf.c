/*****************************************
 *
 * Gram-Charlier PDF 
 * PDF_params[7] will set how the anisotropy is accounted for.
 *              =0 - up wind/cross wind passed through.
 *              =1 - x,y slopes passed through, and complete 2X2 matrix *                  of variances used.
 * Note that the order of u and c have been reversed from older 
 * versions of the code !!! This is to make them agree with other
 * functions. 
 ******************************************/
#include <math.h>

extern void gram_charlierf(
                   double tan_betau, 
                   double tan_betac, 
                   double PDF_params[], 
                   double *pval)
{
  double sigmac, sigmau, c21, c03, c40, c22, c04, ksi, eta, pnorm;
  double normalpdf;
  double lambda1, lambda2, lambda3;
  double pi;
  pi = 3.141593E0;

  if (PDF_params[7] ==0)
    {
     sigmau = PDF_params[0];
     sigmac = PDF_params[1];
     c21 = PDF_params[2];
     c03 = PDF_params[3];
     c40 = PDF_params[4];
     c22 = PDF_params[5];
     c04 = PDF_params[6];
     ksi = tan_betau/sigmau;
     eta = tan_betac /sigmac;
     pnorm = exp( -(ksi*ksi + eta*eta)/2.0E0) *
     ( 1 - c21 *(ksi*ksi -1) * eta/2 - c03*(eta*eta*eta - 3*eta)/6 + 
      c40 * (ksi*ksi*ksi*ksi - 6 * ksi*ksi + 3)/24 + 
      c22 * (ksi*ksi - 1)*(eta*eta - 1)/4 + 
      c04 * (eta*eta*eta*eta - 6 * eta*eta + 3)/24 ); 
  
    *pval = pnorm / (2 * pi * sigmau * sigmac); 
    }
  else
    {
      lambda1 = PDF_params[0];
      lambda2 = PDF_params[1];
      lambda3 = PDF_params[2];
      pnorm = exp( -0.5E0 * ( lambda1 * tan_betau * tan_betau -
                  2.0E0 * lambda3 * tan_betau * tan_betac + 
			    lambda2 * tan_betac * tan_betac ));
      *pval = pnorm * sqrt( lambda1 * lambda2 - lambda3 * lambda3)/
	              ( 2.0E0 * pi); 

    }

}


