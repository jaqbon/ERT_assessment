#ifndef IRI_H
#define IRI_H
//------------------------------------------------------------------------------

#include <stdint.h>

//------------------------------------------------------------------------------

typedef int32_t flogical;

typedef int32_t finteger;

typedef float freal;

//------------------------------------------------------------------------------

void read_ig_rz_(void);

void readapf107_(void);

void iri_sub_(flogical* jf,
              flogical* jmag,
              freal* alati,
              freal* along,
              finteger* iyyyy,
              finteger* mmdd,
              freal* dhour,
              freal* heibeg,
              freal* heiend,
              freal* heistp,
              freal outf[20*1000],
              freal oarr[100]);

//------------------------------------------------------------------------------
#endif
