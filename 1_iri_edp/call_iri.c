#include "iri.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define OUTF(row, col) outf[((row)-1) + ((col)-1)*outf_stride]
#define OAR(row, col) oar[((row)-1) + ((col)-1)*oar_stride]


int main(void)
{
  FILE *fp;

  flogical jmag, jf[50];
  freal alati, along;
  finteger iyyyy, mmdd[2], iut, ivar;
  freal dhour[2];
  freal heibeg, heiend, heistp;
  freal hxx, h_tec_max;
  freal outf[20*1000], oar[100*1000], oarr[100];

  float edp_to_mhz = 1.602e-19/sqrt(9.11e-31*8.85e-12)/(2.0*M_PI)/1000.0;


  int outf_stride = 20;
  int oar_stride = 100;


  jmag = 0;

  for (int i = 0; i < 50; i++)
    jf[i] = 1;

  jf[ 3] =
  jf[ 4] =
  jf[ 5] =
  jf[20] =
  jf[22] =
  jf[27] =
  jf[28] =
  jf[29] =
  jf[32] =
  jf[34] =
  jf[38] =
  jf[39] =
  jf[46] = 0;

  for (int i = 0; i < 100; i++)
    oar[i] = oarr[i] = -1.0;

  alati =  37.8;
  along = -75.4; // iri_sub checks for < 0, adds 360

  iyyyy = 2021;

  mmdd[0] = 303;
  dhour[0] = 11.0;

  mmdd[1] = 304;
  dhour[1] = 23.0;

  iut = 1;

  int N = 51; // number of sample heights (<= 100)

  heibeg = 100.0;
  heiend = 600.0;
  heistp = (heiend - heibeg)/(N-1);

  hxx = heibeg;

  h_tec_max = 0.0;

  ivar = 1; // altitude


  int use_iri_sub = 1;


  read_ig_rz_();
  readapf107_();


  for (int i_date_time = 0; i_date_time < 2; i_date_time++)
  {
    if (use_iri_sub)
    {
     // add 25 for UTC (iri_web does this before calling iri_sub)
      freal hour = dhour[i_date_time] + 25.0;

      iri_sub_(jf, &jmag,
               &alati, &along,
               &iyyyy, &mmdd[i_date_time], &hour,
               &heibeg, &heiend, &heistp,
               outf, oarr);

      for (int i = 0; i < 100; i++)
        oar[i] = oarr[i];
    }
    else
    {
      iri_web_(&jmag, jf,
               &alati, &along,
               &iyyyy, &mmdd[i_date_time], &iut, &dhour[i_date_time],
               &hxx,
               &h_tec_max,
               &ivar,
               &heibeg, &heiend, &heistp,
               outf,
               oar);
    }


    char name[32];

    sprintf(name, "edp_%d_%.2f.dat", mmdd[i_date_time], dhour[i_date_time]);

    fp = fopen(name, "w");


    freal xcor = heibeg;

    for (int li = 1; li <= N; li++)
    {
      OAR( 1,li) = OAR( 1,1);
      OAR(37,li) = OAR(37,1);
      OAR(38,li) = OAR(38,1);

      finteger jne = (int)(OUTF(1,li)/1.e6 + 0.5);

      freal xner = OUTF(1,li)/OAR(1,li);

      finteger jtn = (int)(OUTF(2,li) + 0.5);
      finteger jti = (int)(OUTF(3,li) + 0.5);
      finteger jte = (int)(OUTF(4,li) + 0.5);

      freal scid = 1.0e-8;

      if (jf[21]) scid = 10.0;

      finteger jio  = (int)(OUTF( 5,li)*scid + 0.5);
      finteger jih  = (int)(OUTF( 6,li)*scid + 0.5);
      finteger jihe = (int)(OUTF( 7,li)*scid + 0.5);
      finteger jio2 = (int)(OUTF( 8,li)*scid + 0.5);
      finteger jino = (int)(OUTF( 9,li)*scid + 0.5);
      finteger jicl = (int)(OUTF(10,li)*scid + 0.5);
      finteger jin  = (int)(OUTF(11,li)*scid + 0.5);

      if (OUTF( 1,li) < 0) jne  = -1;
      if (OUTF( 1,li) < 0) xner = -1.0;
      if (OUTF( 2,li) < 0) jtn  = -1;
      if (OUTF( 3,li) < 0) jti  = -1;
      if (OUTF( 4,li) < 0) jte  = -1;
      if (OUTF( 5,li) < 0) jio  = -1;
      if (OUTF( 6,li) < 0) jih  = -1;
      if (OUTF( 7,li) < 0) jihe = -1;
      if (OUTF( 8,li) < 0) jio2 = -1;
      if (OUTF( 9,li) < 0) jino = -1;
      if (OUTF(10,li) < 0) jicl = -1;
      if (OUTF(11,li) < 0) jin  = -1;

      freal tec = OAR(37,li);
      finteger itopp;
      if (tec > 0.0)
      {
        tec /= 1.0e16;
        itopp = (int)(OAR(38,li) + 0.5);
      }
      else
      {
        tec = -1;
        itopp = -1;
      }

      // output full table, like iritest
      //
      //fprintf(fp, "%f %d %f %d %d %d %d %d %d %d %d %d %d %f %d\n",
      //      xcor, jne, xner, jtn, jti, jte, jio, jin, jih, jihe, jio2, jino, jicl, tec, itopp);

      fprintf(fp, "%.2f %.1f\n", (edp_to_mhz*sqrt((double)jne)), xcor);

      xcor += heistp;
    }

    fclose(fp);
  }


  return 0;
}
