#include <fftw3.h>
#include <math.h>
#include <stdio.h>


// dimensions of raw 0
#define NX0 (187)
#define NY0 (  9)

#define NXY0 (NX0*NY0)

// dimensions of raw 1
#define NX1 (409)
#define NY1 ( 11)

#define NXY1 (NX1*NY1)

// dimensions of padded array
// generous, but powers of 2 are easy to think about, numbers in this problem aren't too large
#define NX (1024) // need at least NX0 + NX1 - 1
#define NY (  32) // need at least NY0 + NY1 - 1

#define NXY (NX*NY)

// dimensions of fftw arrays
#define N1 NY
#define N2 (NX+2) // space for one extra complex number in each row

#define N (N1*N2) // total array size


double dx = 0.04;

double dy = 1.0/3000.0;


int ix(double x, int i0)
{
  return (i0 + (int) floor(x/dx + 0.5) + NX) % NX;
}

double xi(int i, int i0)
{
  return (i - i0)*dx;
}


int iy(double y, int i0)
{
  return (i0 + (int) floor(y/dy + 0.5) + NY) % NY;
}

double yi(int i, int i0)
{
  return (i - i0)*dy;
}


int main(void)
{
  double s0[NXY] = {0}, s1[NXY] = {0}, s0_for_plot[NXY];
  double S0[N], S1[N], S[N], s[NXY];
  FILE *fp;


//------------------------------------------------------------------------------
// read data

  fp = fopen("convolve_raw_0.dat", "r");

  for (int i = 0; i < NXY0; i++)
  {
    double x,y,z;

    fscanf(fp, "%lf %lf %lf\n", &x, &y, &z);

    s0[ix(x,0) + NX*iy(y,0)] = z; // center at corner for fft convolution

    s0_for_plot[ix(x,NX/2) + NX*iy(y,NY/3)] = z; // center at offset for plot
  }

  fclose(fp);


  fp = fopen("convolve_raw_1.dat", "r");

  for (int i = 0; i < NXY1; i++)
  {
    double x,y,z;

    fscanf(fp, "%lf %lf %lf\n", &x, &y, &z);

    s1[ix(x,NX/2) + NX*iy(y,NY/3)] = z; // center at offset for plot
  }

  fclose(fp);
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// write input for plots, both linear and log

  fp = fopen("raw0.dat", "w");

  for (int j = 0; j < NY; j++)
  {
    for (int i = 0; i < NX; i++)
    {
      double z = s0_for_plot[i + j*NX];

      fprintf(fp, "%f %f %f\n", xi(i,NX/2), yi(j,NY/3), z);
    }

    fputc('\n', fp);
  }

  fclose(fp);


  fp = fopen("raw1.dat", "w");

  for (int j = 0; j < NY; j++)
  {
    for (int i = 0; i < NX; i++)
    {
      double z = s1[i + j*NX];

      fprintf(fp, "%f %f %f\n", xi(i,NX/2), yi(j,NY/3), z);
    }

    fputc('\n', fp);
  }

  fclose(fp);

  fp = fopen("raw0_log.dat", "w");

  for (int j = 0; j < NY; j++)
  {
    for (int i = 0; i < NX; i++)
    {
      double z = s0_for_plot[i + j*NX];

      fprintf(fp, "%f %f %f\n", xi(i,NX/2), yi(j,NY/3), (z == 0.0 ? -10.0 : log10(z)));
    }

    fputc('\n', fp);
  }

  fclose(fp);


  fp = fopen("raw1_log.dat", "w");

  for (int j = 0; j < NY; j++)
  {
    for (int i = 0; i < NX; i++)
    {
      double z = s1[i + j*NX];

      fprintf(fp, "%f %f %f\n", xi(i,NX/2), yi(j,NY/3), (z == 0.0 ? -10.0 : log10(z)));
    }

    fputc('\n', fp);
  }

  fclose(fp);
//------------------------------------------------------------------------------


  // arrays are small, so don't bother with in-place
  fftw_plan pf0 = fftw_plan_dft_r2c_2d(NY, NX, s0, (fftw_complex*)S0, FFTW_ESTIMATE);
  fftw_plan pf1 = fftw_plan_dft_r2c_2d(NY, NX, s1, (fftw_complex*)S1, FFTW_ESTIMATE);

  fftw_plan pb0 = fftw_plan_dft_c2r_2d(NY, NX, (fftw_complex*)S0, s0, FFTW_ESTIMATE);
  fftw_plan pb1 = fftw_plan_dft_c2r_2d(NY, NX, (fftw_complex*)S1, s1, FFTW_ESTIMATE);

  fftw_plan pb = fftw_plan_dft_c2r_2d(NY, NX, (fftw_complex*)S, s, FFTW_ESTIMATE);
  

  fftw_execute(pf0);
  fftw_execute(pf1);


  int rl = 0, im = 1;

  for (; rl < N; rl += 2, im += 2)
  {
    S[rl] = S0[rl]*S1[rl] - S0[im]*S1[im];
    S[im] = S0[rl]*S1[im] + S0[im]*S1[rl];
  }


  fftw_execute(pb0);
  fftw_execute(pb1);

  fftw_execute(pb);
    

  double norm = 1.0/NXY;

  for (int i = 0; i < NXY; i++)
  {
    s0[i] *= norm;
    s1[i] *= norm;

    s[i] *= norm;
  }


//------------------------------------------------------------------------------
// write convolution result, both linear and log

  fp = fopen("conv.dat", "w");

  for (int j = 0; j < NY; j++)
  {
    for (int i = 0; i < NX; i++)
    {
      double z = s[i + j*NX];

      fprintf(fp, "%f %f %f\n", xi(i,NX/2), yi(j,NY/3), z); // centered at offset
    }

    fputc('\n', fp);
  }

  fclose(fp);


  fp = fopen("conv_log.dat", "w");

  for (int j = 0; j < NY; j++)
  {
    for (int i = 0; i < NX; i++)
    {
      double z = s[i + j*NX];

      fprintf(fp, "%f %f %f\n", xi(i,NX/2), yi(j,NY/3), (z <= 0.0 ? -10.0 : log10(z))); // centered at offset
    }

    fputc('\n', fp);
  }

  fclose(fp);
//------------------------------------------------------------------------------


  fftw_destroy_plan(pf0);
  fftw_destroy_plan(pf1);
  
  fftw_destroy_plan(pb0);
  fftw_destroy_plan(pb1);
  
  fftw_destroy_plan(pb);


  return 0;
}
