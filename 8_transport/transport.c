// based on (incorrectly, no doubt)
//
// LCPFCT - A Flux-Corrected Transport Algorithm for Solving Generalized Continuity Equations
//
// Boris, Landsberg, Oran, and Gardner
//
// https://apps.dtic.mil/sti/pdfs/ADA265011.pdf
//
//
// pages 20 and 21


#include "flux_corrected_transport.h"

#include <stdio.h>
#include <string.h>


#define N (1000)

void write_dat(double *r, FILE *fp)
{
  int i;

  for (i = 0; i < N && r[i] == 0.0; i++);

  fprintf(fp, "%d 0.0\n", i);

  for (; i < N && r[i] != 0.0; i++)
    fprintf(fp, "%d %g\n", i, r[i]);

  if (i != N)
    fprintf(fp, "%d 0.0\n", i-1);
}


int main(void)
{
  double r[N]={0}, rv[N]={0}, v[N];
  double r_half[N], rv_half[N], v_half[N];
  double r_full[N], rv_full[N];

  double Cmax = 0.25;
  double t = 0.0, dx = 1.0, dt, dt0 = 0.01;


  for (int i = 50; i < 150; i++)
  {
    r[i] = 25.0;

    rv[i] = 25.0*10.0;
  }


  FILE *fp = fopen("out.dat", "w");

  write_dat(r, fp);


  while (t < 10.0)
  {
    // calculate v
    for (int i = 0; i < N; i++)
      v[i] = (r[i] ? rv[i]/r[i] : 0.0);


    double vmax = 0;

    for (int i = 0; i < N; i++)
      if (v[i] > vmax)
        vmax = v[i];

    double dt_max = Cmax*dx/vmax;

    dt = (dt0 < dt_max ? dt0 : dt_max);


    // convect r half step
    flux_corr_method(r, v, N, 0.5*dt, dx, r_half);

    // convect rv half step
    flux_corr_method(rv, v, N, 0.5*dt, dx, rv_half);


    for (int i = 0; i < N ; i++)
      if (r_half[i] <= 1.0e-6)
        r_half[i] = 0.0;

    // calculate v_half
    for (int i = 0; i < N; i++)
      v_half[i] = (r_half[i] ? rv_half[i]/r_half[i] : 0.0);


    // convect r full step
    flux_corr_method(r, v_half, N, dt, dx, r_full);

    // convect rv full step
    flux_corr_method(rv, v_half, N, dt, dx, rv_full);


    memcpy(r, r_full, N*sizeof(double));

    memcpy(rv, rv_full, N*sizeof(double));


    for (int i = 0; i < N ; i++)
      if (r[i] <= 1.0e-6)
        r[i] = 0.0;


    t += dt;


    fputc('\n', fp);
    fputc('\n', fp);

    write_dat(r, fp);
  }


  fclose(fp);


  return 0;
}
