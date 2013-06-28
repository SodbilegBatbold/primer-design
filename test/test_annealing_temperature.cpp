#include <cmath>
#include "test_annealing_temperature.h" 
#include "../pd5/marmur_temperature.h"
#include "../pd5/freier_temperature.h"
#include "../pd5/wallace_temperature.h"
#include "../pd5/primer3_Tm_temperature.h"

static char test_seq[] = "ATTGCGTATTATGC";

char * test_marmur() {
  marmur_temperature a;
  double t = a.calculate_temperature(test_seq);
  mu_assert((char *)"test_marmur error, temperature != 38", t == 38.0);
  return 0;
}


char * test_freier() {
  freier_temperature a;
  double t = a.calculate_temperature(test_seq);
  mu_assert((char *)"test_freier error, temperature != 26.3286", abs(t - 26.3286) < 0.1 );
  return 0;
}


char * test_wallace() {
  wallace_temperature a;
  double t = a.calculate_temperature(test_seq);
  mu_assert((char *)"test_Wallace error, temperature != 31.5143", abs(t - 31.5143) < 0.1);
  return 0;
}


char * test_primer3_Tm() {
  primer3_Tm_temperature a;
  double t = a.calculate_temperature(test_seq);
  mu_assert((char *)"test_primer3_Tm error, temperature != 34.5963", abs(t - 34.5963) < 0.1);
  return 0;
}
 

