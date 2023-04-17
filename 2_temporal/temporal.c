#include <ctype.h>
#include <stdlib.h>
#include <stdio.h> 
#include <string.h>


typedef struct
{
  int seconds; // sorted value

  float foF2;
  float hmF2;
}
record;


char* next_field(char *p)
{
  while (isspace(*p)) ++p;
  while (!isspace(*p)) ++p;

  return p;
}


int read_input(const char *filename, int *num, record **data)
{
  FILE *fp = fopen(filename, "r");
  
  if (!fp) return -1;


  char *p, buf[128]; // sufficient for this file
  int n; // local
  record *d; // local
  int rtn = 0;


  for (int loop = 0; loop < 2; loop++)
  {
    n = 0;

    fgets(buf, 128, fp); // skip header and blank line
    fgets(buf, 128, fp);

    while (p = fgets(buf, 128, fp))
    {
      if (loop)
      {
//        yyyy.MM.dd (DDD) HH:mm:ss C-score   foF2  foF1   foE  foEs   h`Es   hmF2   hmF1    hmE    B0   B1   D1
//        2020.01.31 (031) 22:30:00      100 4.913    --- 2.420 2.400 105.000 192.800    --- 103.900 63.400 1.800  

        p = buf + 12; // use known file structure, months and days written with leading 0s, so fixed size

        int day = strtol(p, &p, 10);

        int hour = strtol(++p, &p, 10);
        int min  = strtol(++p, &p, 10);
        int sec  = strtol(++p, &p, 10);

        d[n].seconds = sec + 60*(min + 60*(hour + 24*day)); // seconds since start of year

        p = buf + 35; // continue to use known file structure

        d[n].foF2 = (float)strtod(p, &p);

        p = next_field(p); // skip foF1
        p = next_field(p); // skip foE
        p = next_field(p); // skip foEs
        p = next_field(p); // skip h`Es

        d[n].hmF2 = (float)strtod(p, &p);
      }

      n++;
    }

    if (!loop)
    {
      if ((*num = n) == 0) { rtn = -2; break; }

      if ((*data = d = (record*) malloc(n*sizeof(record))) == NULL) { rtn = -3; break; }
    }

    rewind(fp);
  }


  fclose(fp);


  return rtn;
}


int comp(const void *rec1, const void *rec2)
{
  int s1 = ((record*)rec1)->seconds;
  int s2 = ((record*)rec2)->seconds;

  return (s1 > s2 ? 1 : s1 == s2 ? 0 : -1);
}


void median(const record *data, int i, record *m)
{
  m->seconds = data[i].seconds;

  if ((data[i-1].foF2 <= data[i].foF2) && (data[i-1].foF2 <= data[i+1].foF2))
    m->foF2 = (data[i].foF2 <= data[i+1].foF2 ? data[i].foF2 : data[i+1].foF2);
  else if ((data[i].foF2 <= data[i-1].foF2) && (data[i].foF2 <= data[i+1].foF2))
    m->foF2 = (data[i-1].foF2 <= data[i+1].foF2 ? data[i-1].foF2 : data[i+1].foF2);
  else
    m->foF2 = (data[i-1].foF2 <= data[i].foF2 ? data[i-1].foF2 : data[i].foF2);
  
  if ((data[i-1].hmF2 <= data[i].hmF2) && (data[i-1].hmF2 <= data[i+1].hmF2))
    m->hmF2 = (data[i].hmF2 <= data[i+1].hmF2 ? data[i].hmF2 : data[i+1].hmF2);
  else if ((data[i].hmF2 <= data[i-1].hmF2) && (data[i].hmF2 <= data[i+1].hmF2))
    m->hmF2 = (data[i-1].hmF2 <= data[i+1].hmF2 ? data[i-1].hmF2 : data[i+1].hmF2);
  else
    m->hmF2 = (data[i-1].hmF2 <= data[i].hmF2 ? data[i-1].hmF2 : data[i].hmF2);
}


int main(void)
{
  int rtn;
  int num;
  record *data = NULL;

  
  if (rtn = read_input("AU930_ROAM.TXT", &num, &data))
  {
    switch(rtn)
    {
      case -1:
        printf("can't open AU930_ROAM.TXT\n");
        break;

      case -2:
        printf("no data\n");
        break;

      case -3:
        printf("malloc error\n");
        break;
    }

    return rtn;
  }


  qsort(data, num, sizeof(record), comp);


  FILE *fp = fopen("median.txt", "w");

  if (fp)
  {
    fprintf(fp, "seconds foF2(in) foF2(out) hmF2(in) hmF2(out)\n\n");

    for (int i = 1; i < num-1; i++)
    {
      record m;

      median(data, i, &m);

      fprintf(fp, "%d %.3f %.3f %.3f %.3f\n", m.seconds, data[i].foF2, m.foF2, data[i].hmF2, m.hmF2);
    }

    fclose(fp);
  }
  else
    printf("can't open median.txt\n");

  if (data) free(data);


  return 0;
}
