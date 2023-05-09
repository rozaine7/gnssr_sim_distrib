/*********************************************************************
 *
 * Lagrange interpolation -called by lagrangem.c
 *
 *********************************************************************/

#include <math.h>
#include "mex.h"


extern void lagrangec(
 		   double	gps_ecef_component[],
 		   double	gpstime[],
                   double         prnset[],
                   double       interptime[],
                   double       *interp_gps_ecef_component,
                   int          ninterp_pts,
                   int          ngps,
                   double       *test1
		   )
{

  int interp_pt, kgps, lgps,kstartgps, kstopgps, iprn;
  double lnumerator, ldenominator, sum_ecef_interp;
  double tbuffer, gpsecef;

  tbuffer = 3600.0;

  for ( interp_pt = 0; interp_pt < ninterp_pts; interp_pt++)
    {
      kstartgps = 0;
      kstopgps = 0;
      sum_ecef_interp = 0;
      if (prnset[ interp_pt] > 0)
	{
         iprn = ((int) prnset[ interp_pt]) - 1;

         while ( gpstime[kstartgps] <  (interptime[interp_pt] - tbuffer) )
         {
        	kstartgps++;
         } 
         kstopgps = kstartgps;

         while ( gpstime[kstopgps] <  (interptime[interp_pt] + tbuffer) )
         {
	     kstopgps++;
         }
         
         /* printf(" kstart = %4i  kstop = %4i \n", kstartgps, kstopgps);*/

         for (kgps=kstartgps; kgps < kstopgps; kgps++)
	 {
            lnumerator = 1;
            ldenominator = 1;
            gpsecef = gps_ecef_component[ngps*iprn + kgps];

        	 /*  printf("GPS ECEF: %12.1f\n", gpsecef);*/

            for (lgps=kstartgps; lgps <kstopgps; lgps++)
            {  
	      if (lgps != kgps)
		  {
        	lnumerator = lnumerator * 
                   ( interptime[interp_pt] - gpstime[lgps]);
            ldenominator = ldenominator * (gpstime[kgps]-gpstime[lgps]);
                  *test1 = prnset[ interp_pt];
                  
           /* printf(" kgps = %4i  lgps = %4i \n", kgps, lgps);*/
          }
        }

          sum_ecef_interp = sum_ecef_interp + 
        	   (lnumerator / ldenominator) * gpsecef;
          }

          *(interp_gps_ecef_component+interp_pt) = sum_ecef_interp;
    
	}
      else
	{
          *(interp_gps_ecef_component+interp_pt) = 0.0;
        }
    }
}
