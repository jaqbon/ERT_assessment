#include <math.h>
#include <stdio.h>


const double deg_to_rad = M_PI/180;
const double rad_to_deg = 180/M_PI;

const double earth_radius = 6378100.0; // meters


typedef struct
{
  double x;
  double y;
  double z;
}
vector;


double dot(const vector *v1, const vector *v2)
{
  return (v1->x*v2->x + v1->y*v2->y + v1->z*v2->z);
}


void gis_to_unit_vectors(double lon,
                         double lat,
                         vector *r, // radial
                         vector *n, // north tangent
                         vector *e) // east  tangent
{
  double l = deg_to_rad*lon;
  double L = deg_to_rad*lat;

  double c = cos(l);
  double s = sin(l);
  
  double C = cos(L);
  double S = sin(L);
  
  if (r)
  {
    r->x = C*c;
    r->y = C*s;
    r->z = S;
  }

  if (n)
  {
    n->x = -S*c;
    n->y = -S*s;
    n->z =  C;
  }

  if (e)
  {
    e->x = -s;
    e->y =  c;
    e->z =  0;
  }
}


void radial_vector_to_gis(const vector *r, double *lon, double *lat)
{
  *lon = rad_to_deg*atan2(r->y, r->x);

  *lat = rad_to_deg*atan2(r->z, sqrt(r->x*r->x + r->y*r->y));
}


int gis_to_radar(double *range,
                 double *bearing,
                 double lon0,
                 double lat0,
                 double lon1,
                 double lat1)
{
  vector r0, n0, e0, r1;

  gis_to_unit_vectors(lon0, lat0, &r0, &n0, &e0);

  gis_to_unit_vectors(lon1, lat1, &r1, NULL, NULL);

  *bearing = rad_to_deg*atan2(dot(&e0,&r1), dot(&n0,&r1));

  *range = earth_radius*acos(dot(&r0, &r1));

  return 0;
}


int radar_to_gis(double range,
                 double bearing,
                 double lon0,
                 double lat0,
                 double *lon1,
                 double *lat1)
{
  double r0_dot_r1 = cos(range/earth_radius);

  double mult = sqrt(1 - r0_dot_r1*r0_dot_r1);

  double n_dot_t = cos(deg_to_rad*bearing);
  double e_dot_t = sin(deg_to_rad*bearing);

  vector r0, n0, e0, r1;

  gis_to_unit_vectors(lon0, lat0, &r0, &n0, &e0);

  r1.x = r0.x*r0_dot_r1 + n0.x*mult*n_dot_t + e0.x*mult*e_dot_t;
  r1.y = r0.y*r0_dot_r1 + n0.y*mult*n_dot_t + e0.y*mult*e_dot_t;
  r1.z = r0.z*r0_dot_r1 + n0.z*mult*n_dot_t + e0.z*mult*e_dot_t;

  radial_vector_to_gis(&r1, lon1, lat1);
}



int main(void)
{
  double lat_wallops =  37;
  double lon_wallops = -75;

  double lat_pr =  18;
  double lon_pr = -66;

  double range, bearing;

  gis_to_radar(&range, &bearing,
               lon_wallops, lat_wallops,
               lon_pr, lat_pr);

  printf("range = %f, bearing = %f\n\n", range, bearing);

  // range = 2291214.151280, bearing = 154.963217

  // Google Maps says 2288670, so (range - GM)/GM = 0.0011115626106

  radar_to_gis(range, bearing,
               lon_wallops, lat_wallops,
               &lon_pr, &lat_pr);

  printf("lon_pr = %f, lat_pr = %f\n", lon_pr, lat_pr);


  return 0;
}

