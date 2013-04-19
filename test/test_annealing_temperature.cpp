#include <cmath>
#include "test_annealing_temperature.h" 
#include "../pd5/annealing_temperature.h"

static char test_seq[] = "ATTGCGTATTATGC";

char * test_marmur() {
  annealing_temperature a;
  double t = a.Marmur_method(test_seq);
  mu_assert((char *)"test_marmur error, temperature != 38", t == 38.0);
  return 0;
}


char * test_freier() {
  annealing_temperature a;
  double t = a.Freier_method(test_seq);
  mu_assert((char *)"test_freier error, temperature != 26.3286", abs(t - 26.3286) < 0.1 );
  return 0;
}


char * test_wallace() {
  annealing_temperature a;
  double t = a.Wallace_method(test_seq);
  mu_assert((char *)"test_Wallace error, temperature != 31.5143", abs(t - 31.5143) < 0.1);
  return 0;
}


char * test_primer3_Tm() {
  annealing_temperature a;
  double t = a.primer3_Tm(test_seq);
  mu_assert((char *)"test_primer3_Tm error, temperature != 34.5963", abs(t - 34.5963) < 0.1);
  return 0;
}
 

