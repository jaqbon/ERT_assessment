#include <math.h>
#include <stdio.h>


int main(void)
{
  const double deg_to_rad = M_PI/180.0;

#define NPTS (10)

  double points[NPTS][3] =
  {
    {121.39, 13.51, 1.494},
    {126.19, 12.02, 1.934},
    {130.27, 13.11, 2.148},
    {127.42, 10.09, 9.155},
    {126.14, 15.33, 2.221},
    {125.96, 14.00, 8.100},
    {123.15, 10.88, 2.039},
    {130.50, 11.18, 1.916},
    {129.08, 15.78, 3.729},
    {122.74, 15.82, 7.137}
  };

  int nlat = 50;
  int nlon = 70;

  double lat0 = 10.0;
  double lat1 = 16.0;

  double lon0 = 121.0;
  double lon1 = 131.0;


  double C0[NPTS], S0[NPTS];

  for (int ipts = 0; ipts < NPTS; ipts++)
  {
    C0[ipts] = cos(deg_to_rad*points[ipts][1]); // latitude
    S0[ipts] = sin(deg_to_rad*points[ipts][1]);
  }
    

  FILE *fp = fopen("interp.dat", "w");


  for (int ilat = 0; ilat < nlat; ilat++)
  {
    double lat = (lat0*(nlat-1-ilat) + lat1*ilat)/(double)(nlat-1);

    double C = cos(deg_to_rad*lat);
    double S = sin(deg_to_rad*lat);
    
    for (int ilon = 0; ilon < nlon; ilon++)
    {
      double lon = (lon0*(nlon-1-ilon) + lon1*ilon)/(double)(nlon-1);
      
      double sum = 0, norm = 0;

      for (int ipts = 0; ipts < NPTS; ipts++)
      {
        // distance on unit sphere
        double d = acos(S*S0[ipts] + C*C0[ipts]*cos(deg_to_rad*(lon - points[ipts][0])));

        if (d == 0)
        {
          sum = points[ipts][2];

          norm = 1;

          break;
        }

        d = d*d; // extend influence of points, smoother graph

        sum += points[ipts][2]/d;

        norm += 1/d;
      }

      double interp = sum/norm;

      fprintf(fp, "%.2f %.2f %.2f\n", lon, lat, interp);
    }

    fputc('\n', fp);
  }


  fclose(fp);


  return 0;
}

  
