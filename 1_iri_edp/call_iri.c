#include "iri.h"

#include <stdio.h>
#include <stdlib.h>

int main(void)
{
  flogical jmag, jf[50]={1}; // most defaults are true
  freal alati, along;
  finteger iyyyy, mmdd;
  freal dhour;
  freal heibeg, heiend, heistp;
  freal outf[20*1000], oarr[100];


  jmag = 0;

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
    oarr[i] = -1;


  alati =  37.8;
  along = -75.4;

  iyyyy = 2021;
  mmdd  = -62;// 303;
  dhour = 11.0 + 25.0; // add 25 for UTC

  int N = 51; // number of sample heights (<= 100)

  heibeg = 100.0;
  heiend = 600.0;
  heistp = (heiend - heibeg)/(N-1);


  read_ig_rz_();
  readapf107_();


  iri_sub_(jf, &jmag,
           &alati, &along,
           &iyyyy, &mmdd, &dhour,
           &heibeg, &heiend, &heistp,
           outf, oarr);


  FILE *fp = fopen("outf.dat", "w");

  for (int col = 0; col < N; col++)
  {
    for (int row = 0; row < 20; row++)
      fprintf(fp, "%.3g ", outf[row + col*20]);

    fputc('\n', fp);
  }

  fclose(fp);


  return 0;
}


/*
c
c Programs using subroutine IRI_SUB need to include (see IRITEST.FOR):
c
c		call read_ig_rz
c       call readapf107
c

       SUBROUTINE IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &    HEIBEG,HEIEND,HEISTP,OUTF,OARR)
C-----------------------------------------------------------------
C
C INPUT:  JF(1:50)      true/false switches for several options
C         JMAG          =0 geographic   = 1 geomagnetic coordinates
C         ALATI,ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
C         IYYYY         Year as YYYY, e.g. 1985
C         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
C         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL 
C                          HOURS
C         HEIBEG,       HEIGHT RANGE IN KM; maximal 100 heights, i.e.
C          HEIEND,HEISTP        int((heiend-heibeg)/heistp)+1.le.100

C  OUTPUT:  OUTF(1:20,1:1000)
C            OARR(1:100)   ADDITIONAL OUTPUT PARAMETERS         


*/
/*
C    JF switches to turn off/on (.true./.false.) several options
C
C    i       .true.                  .false.          standard version
C    -----------------------------------------------------------------
C    1    Ne computed            Ne not computed                     t
C    2    Te, Ti computed        Te, Ti not computed                 t
C    3    Ne & Ni computed       Ni not computed                     t
C    4    B0,B1 - Bil-2000       B0,B1 - other models jf(31)     false
C    5    foF2 - CCIR            foF2 - URSI                     false
C    6    Ni - DS-1995 & DY-1985 Ni - RBV-2010 & TBT-2015        false
C    7    Ne - Tops: f10.7<188   f10.7 unlimited                     t            
C    8    foF2 from model        foF2 or NmF2 - user input           t
C    9    hmF2 from model        hmF2 or M3000F2 - user input        t
C   10    Te - Standard          Te - Using Te/Ne correlation        t
C   11    Ne - Standard Profile  Ne - Lay-function formalism         t
C   12    Messages to unit 6     to messages.txt on unit 11          t
C   13    foF1 from model        foF1 or NmF1 - user input           t
C   14    hmF1 from model        hmF1 - user input (only Lay version)t
C   15    foE  from model        foE or NmE - user input             t
C   16    hmE  from model        hmE - user input                    t
C   17    Rz12 from file         Rz12 - user input                   t
C   18    IGRF dip, magbr, modip old FIELDG using POGO68/10 for 1973 t
C   19    F1 probability model   only if foF1>0 and not NIGHT        t
C   20    standard F1            standard F1 plus L condition        t
C (19,20) = (t,t) f1-prob, (t,f) f1-prob-L, (f,t) old F1, (f,f) no F1
C   21    ion drift computed     ion drift not computed          false
C   22    ion densities in %     ion densities in m-3                t
C   23    Te_tops (Bil-1985)     Te_topside (TBT-2012)           false
C   24    D-region: IRI-1990     FT-2001 and DRS-1995                t
C   25    F107D from APF107.DAT  F107D user input (oarr(41))         t
C   26    foF2 storm model       no storm updating                   t
C   27    IG12 from file         IG12 - user                         t
C   28    spread-F probability 	 not computed                    false
C   29    IRI01-topside          new options as def. by JF(30)   false
C   30    IRI01-topside corr.    NeQuick topside model   	     false 
C (29,30) = (t,t) IRIold, (f,t) IRIcor, (f,f) NeQuick, (t,f) IRIcor2
C   31    B0,B1 ABT-2009	     B0 Gulyaeva-1987 h0.5               t   
C (4,31) = (t,t) Bil-00, (f,t) ABT-09, (f,f) Gul-87, (t,f) not used
C   32    F10.7_81 from file     F10.7_81 - user input (oarr(46))    t
C   33    Auroral boundary model on/off  true/false	             false
C   34    Messages on            Messages off                        t
C   35    foE storm model        no foE storm updating           false
C   36    hmF2 w/out foF2_storm  with foF2-storm                     t
C   37    topside w/out foF2-storm  with foF2-storm                  t
C   38    turn WRITEs off in IRIFLIP   turn WRITEs on                t
C   39    hmF2 (M3000F2)         new models                      false
C   40    hmF2 AMTB-model        Shubin-COSMIC model             false
C (39,40) = (t,t) hmF2-old, (f,t) AMTB, (f,f) Shubin, (t,f) not used
C   41    Use COV=F10.7_365      COV=f(IG12) (IRI before Oct 2015)   t
C   42    Te with PF10.7 dep.	 w/o PF10.7 dependance               t
C   43    B0 from model          B0 user input in OARR(10)           t
C   44    B1 from model          B1 user input in OARR(35)           t
C   45    HNEA=65/80km dya/night HNEA user input in OARR(89)         t
C   46    HNEE=2000km 	         HNEE user input in OARR(90)         t
C   47    CGM computation on 	 CGM computation off             false
C      ....
C   50    
C   ------------------------------------------------------------------
*/
