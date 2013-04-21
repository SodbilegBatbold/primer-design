#include <cstdio>
#include "test_sequence_utils.h"
#include "test_annealing_temperature.h"
#include "test_DNAfind.h"
#include "test_dimerisation.h"
#include "test_primer.h"
#include "test_primer_pair.h"


int tests_run = 0;

static char * all_temp_tests() {
  mu_run_test(test_marmur);
  mu_run_test(test_freier);
  mu_run_test(test_wallace);
  mu_run_test(test_primer3_Tm);
  return 0;
}

static char * all_sequence_utils_tests() {
  mu_run_test(test_tail_check);
  mu_run_test(test_sticky_tail_check1);
  mu_run_test(test_sticky_tail_check2);
  mu_run_test(test_reverse_complement);
  mu_run_test(test_nucleotide_content); 
  mu_run_test(test_nucleotide_complement); 
  return 0;
}
 
static char * all_dna_find_tests() {
  mu_run_test(test_search_for_binding_sites);
  mu_run_test(test_search_for_pcr_products);
  mu_run_test(test_mismatch_search_for_pcr_products);
  return 0;
}

static char * all_dimerisation_tests() {
  mu_run_test(test_hairpin1);
  mu_run_test(test_hairpin2);
  mu_run_test(test_self_dimer);
  mu_run_test(test_pair_dimer);
  return 0;
}

static char * all_primer_tests() {
  mu_run_test(test_generate_100_candidate_primers);
  mu_run_test(test_generate_20mer_primers_at_start);
  return 0;
}

static char * all_primer_pair_tests() {
  mu_run_test(test_generate_38_candidate_pairs);
  mu_run_test(test_generate_sort_candidate_pairs);
  return 0;
}



static char * all_tests() {
  char * msg;
  msg = all_temp_tests();
  if(msg != 0) { return msg; }
  msg = all_sequence_utils_tests();
  if(msg != 0) { return msg; }
  msg = all_dna_find_tests();
  if(msg != 0) { return msg; }
  msg = all_dimerisation_tests();
  if(msg != 0) { return msg; }
  msg = all_primer_tests();
  if(msg != 0) { return msg; }
  msg = all_primer_pair_tests();
  if(msg != 0) { return msg; }
  return 0;
}


int main(int argc, char **argv) {
  char *result = all_tests();
  if (result != 0) {
    printf("%s\n", result);
  }
  else {
    printf("ALL TESTS PASSED\n");
  }
  printf("Tests run: %d\n", tests_run);
 
  return result != 0;
}
