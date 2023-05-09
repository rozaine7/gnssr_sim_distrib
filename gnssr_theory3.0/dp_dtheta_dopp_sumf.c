/*********************************************************************
 *
 * MEX file version of GO code for differential reflected power at 
 * given delay (delta) and azimuth (theta)
 * Ver 2.2 - 1-17-2013 - included antenna gain pattern
 * Ver 2.1 - 09-20-07 - expanded to include the work from H. You's thesis - 
 *        to generate the second-order moments of the waveform and the 
 *        correlation spectrum.
 *
 * Ver 1.6 - 03-19-03 - corrected error in ellipse dimensions causing 
 *   problems in low grazing angle (earlier had center of ellipse corrected, 
 *   but was still using \delta \theta coordinates based upon assumption that 
 *   the ellipse center does not move from the specular point with increasing
 *   delay.  Also - first attempt to re-introduce Doppler into the C-MEX 
 *   version of the code.
 *  
 * Presently does not include option for Doppler, and is fixed to the 
 * G/C PDF
 *
 *********************************************************************/


#include <math.h>
#include "mex.h"


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
 		   double	dp_dopp_mat_real[],
                   double       dp_dopp_mat_imag[],
		   double       tbx[],
       	           double       tby[],
                   double       xell[],
                   double       yell[],
                   double       pdf_val[],
                   double       *ftilde_step,
                   double       *ftilde_input,
                   double       correlation_spectrum_real[],
                   double       correlation_spectrum_imag[],
                   double       *wavelength
		   )
{

  /*
     function for numerical integration around an ellipse of constant
     delay.  This function computes the differential power received 
     at the antenna from a surface area located on the ellipse with 
     delay delta, at  polar coordinate theta for a differential area 
     defined by dr, dtheta.
 
     Input parameters are the location of the reflecting facet (delta,
     theta),  reflection geometry: h=alt of receiver, alpha=direction 
     of wind vector,  gamma=elevation of GPS satellite
   
     Constants of the Gram-Charlier series approximation of the PDF

     j = index over delay (delta)
     i = index over azimuth (theta)

  */

  int i,j, nonzero_indicator;
  double xs_col, denominator, denominator2, J1;
  double h, singamma, cosgamma, C2, delta_prime, dS;
  double sintheta, costheta;
  double sinalpha, cosalpha;
  double r_col, x_col, y_col, R_col, ux_col, uy_col, uz_col;
  double vx_col, vy_col, vz_col, nx, ny, nz, nmag, cos_omega, cos_beta;
  double tan_betax, tan_betay, tan_betau, tan_betac, cross_section, dP;
  double PDF_value, R;
  double sinc, sinc_lobes, doppler_filter, pi, fd;
  double VGdotv, VRdotu;
  double sincthreshold, abssinc;
  double phase_rotation_real, phase_rotation_imag, indicator, indicator_step;
  double fdold, fnearest, dthetadf;
  double gain, bhatx, bhaty, bhatz, boresight, xyzmag, uxn_col, uyn_col, uzn_col;

  h = *alt;
  sincthreshold = 1.0E-14;
  sinalpha = sin(*alpha);
  cosalpha = cos(*alpha);
  singamma = sin(*gamma);
  cosgamma = cos(*gamma);

/*  indicator_step =  *ftilde_step*10.0; */
  indicator_step =  *ftilde_step; 
 /* Try making it as narrow as possible  
  indicator_step = 1.0; */
   nonzero_indicator = 0;
   /* *(correlation_spectrum_imag) = 0; */

   /* Changed order of the loops - to make calculation of the delta function in       
    the Doppler frequency easier */

  for (j=0; j<mdelta; j++)
 
    {
      fdold = 0.0;
      *(correlation_spectrum_real+j) = 0.0; 

      for( i=0; i<ntheta; i++)
	{
        costheta = cos( theta_input[i] );
        sintheta = sin( theta_input[i] );

        denominator2 =  costheta*costheta*singamma*singamma + sintheta*sintheta;
        denominator = sqrt(denominator2);

        J1 = cosgamma * costheta / ( singamma * denominator * denominator2);

 
         /* calc of ellipse center - allow for \delta dep */
        
	  xs_col = h * cosgamma/singamma + 
	    (cosgamma/(singamma * singamma)) * delta_input[j]; 

          delta_prime = delta_input[j] + h * singamma;
          C2 = delta_prime * delta_prime / ( singamma * singamma) - h*h;
          
          r_col = sqrt(C2)/ denominator; 

          x_col = xs_col + r_col *  costheta;  
          y_col = r_col * sintheta;
          R_col = sqrt( h*h + x_col*x_col + y_col*y_col);

	/*  unit vectors, u, v and n */

          ux_col =  -x_col/R_col; 
          uy_col = -y_col/R_col; 
          uz_col = h/R_col;

          vx_col = -cosgamma;
          vy_col = 0.0;
          vz_col = -singamma;
          
          /*  Calculation of antenna gain - bhat is unit vector along antenna
              boresight, right now assumes flat Earth */
           
          xyzmag = sqrt(ux_col*ux_col + uy_col*uy_col + uz_col*uz_col);
          uxn_col = ux_col / xyzmag;
          uyn_col = uy_col / xyzmag;
          uzn_col = uz_col / xyzmag;
          
          boresight = -bhatz*uxn_col - bhaty*uyn_col - bhatz*uzn_col; 
          /* mod by JLG 12/2010 to set gain to 1 - return to previous results */
          
          /* simple cardioid */
	     gain = 1.0+cos(boresight); 
         /*  gain = 1; */

        /* Doppler filter converted from MATLAB 3/2003 
         remove hard-coded 
         * nthg JLG 5/2016 */

    /*      wavelength = 0.19; */ 
          pi = 3.1415926535897;

          VGdotv = VG[0] * vx_col + VG[1]*vy_col + VG[2] * vz_col;
          VRdotu = VR[0] * ux_col + VR[1]*uy_col + VR[2] * uz_col;

          fdold = fd;

          fd = (VGdotv - VRdotu)/(*wavelength);

          sinc_lobes = (fd-(*fc)) * (*Ti);

        /* cut out only part of the range bin which is within a few 
         side-lobes of a Doppler bin - this code may be needed for 
         spacecraft simulation 
         if (model_flag(5) > 0)
          doppler_flag = logical( abs(sinc_lobes) < model_flag(5));
          doppler_filter =  sinc( sinc_lobes(doppler_flag) ).^2; 
          ux = ux_col(doppler_flag);
          uy = uy_col(doppler_flag);
          uz = uz_col(doppler_flag);
          vx = vx_col(doppler_flag);
          vy = vy_col(doppler_flag);
          vz = vz_col(doppler_flag);
          r = r_col(doppler_flag);
          x = x_col(doppler_flag);
          y = y_col(doppler_flag);
          R = R_col(doppler_flag);
          theta = theta_col(doppler_flag);
         else
         END of MATLAB code */
      
	  /* force sinc(0) = 1 to prevent problems when VR(1) = 0 */

	  abssinc = fabs(sinc_lobes);
      
	  /* printf("sinclobes %12.5f %18.12f \n", abssinc, sincthreshold);
	   */    

	 if(fabs(sinc_lobes) <= sincthreshold) 
	    { 
	                   sinc = 1.0; 
        }
	 else 
	    {
	         sinc =  sin( pi * sinc_lobes )/ (pi * sinc_lobes);
        }    
  
          doppler_filter = sinc * sinc; 
      
 /*      doppler_filter = 1; */

          nx = ( ux_col - vx_col );
          ny = ( uy_col - vy_col );
          nz = ( uz_col - vz_col );

          nmag = sqrt( nx*nx + ny*ny + nz*nz);

          nx = nx /nmag;
          ny = ny /nmag;
          nz = nz /nmag;

        /* Angles and slopes necessary for PDF and the 
           projected areas */

          cos_beta = nz;
          tan_betax = nx / nz;
          tan_betay = ny / nz;  /* These assume locally-flat Earth. */

          *(tbx+i*mdelta+j) = tan_betax;
          *(tby+i*mdelta+j) = tan_betay;

          *(xell+i*mdelta+j) = x_col;
          *(yell+i*mdelta+j) = y_col;

	   
        /* delete locally spherical Earth fixes ...
          *****************

          if (model_flag(2) == 1)
           surface_normalx1 = (r .* cos(theta));
           surface_normalx2 = zeros( ndelay*ntheta,1);
           surface_normalx3 = sqrt(RE^2 - surface_normalx1.^2);

           surface_normaly1 = zeros (ndelay*ntheta,1);
           surface_normaly2 = (r.*sin(theta));
           surface_normaly3 = sqrt(RE^2 - surface_normaly2.^2);

           tan_betax_correction = surface_normalx1./surface_normalx3;
           tan_betay_correction = surface_normaly2./surface_normaly3;

           tan_betax = tan_betax - tan_betax_correction;
           tan_betay = tan_betay - tan_betay_correction;
          end
          ***************
          END of SPHERICAL EARTH */

        /* rotate slopes to upwind/crosswind coordinates. 
           fixed 5/03 to use a RH defn of alpha */

          tan_betau = ( tan_betax * cosalpha + tan_betay * sinalpha );
          tan_betac = ( -tan_betax * sinalpha + tan_betay * cosalpha );

        /*  compute the scattering cross secion and projected areas */

          gram_charlierf(tan_betau,tan_betac,PDF_params,&PDF_value);

          *(pdf_val+i*mdelta+j) = PDF_value;

          cross_section = PDF_value / 
            ( 4 * cos_beta*cos_beta*cos_beta*cos_beta * R_col * R_col);
  
          dS = delta_prime / (singamma * singamma * denominator2) +
	    J1*sqrt( delta_input[j]*delta_input[j] + 
		     2*h*singamma*delta_input[j]);

          dP = cross_section * dS * gain;
     

	  /* Term added for the phase change over the time ttilde_input */

          phase_rotation_real = cos( 2.0*pi*(fd-(*fc)) * (*ttilde_input));
          phase_rotation_imag = sin( 2.0*pi*(fd-(*fc)) * (*ttilde_input));

          *(dp_dopp_mat_real+i*mdelta+j) = 
                     doppler_filter * dP *phase_rotation_real;
          *(dp_dopp_mat_imag+i*mdelta+j) = 
                     doppler_filter * dP *phase_rotation_imag; 

          /* Correlation spectrum includes the indicator function - this is 
             numerically approximated as a pulse with a width equal to the step
             of the correlation frequency, ftilde 
           
           See derivation notes fro 6/21/2010 */

          
          *(correlation_spectrum_imag+i*mdelta+j) = fd-(*fc);    
          
          if (( i>0) && (fd-(*fc) - *ftilde_input) * (fdold-(*fc) - *ftilde_input) < 0 )
            {        

             /* indicator = 1.0/fabs(fd-fdold); */
              dthetadf = fabs(( theta_input[i]-theta_input[i-1])/(fd - fdold));
             indicator = dthetadf;
              *(correlation_spectrum_real+j) = *(correlation_spectrum_real+j)+
                doppler_filter * dP * indicator; 

          }
          else
          { 
	          indicator = 0.0;
          } 

         
        
                   
     }
  }
}
