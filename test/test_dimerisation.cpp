#include "test_dimerisation.h" 
#include "../pd5/dimerisation.h"

static char test_seq1[] = "ATAGGCGACGCTAA";
static char test_seq2[] = "TGGGGGGGCCCCCCCCCCC";
static char test_seq3[] = "TTAGCGTTTTAAAAAA";

char * test_hairpin1() {
  dimerisation d;
  d.hairpin(test_seq1);
  mu_assert((char *)"hairpin error, hairpin_score for seq1 != 0", d.hairpin_score == 0);
  return 0;
}


char * test_hairpin2() {
  dimerisation d;
  d.hairpin(test_seq2);
  mu_assert((char *)"hairpin error, hairpin_score for seq2 != 18", d.hairpin_score == 18);
  return 0;
}

char * test_self_dimer() {
  dimerisation d;
  d.self_dimer(test_seq2);
  mu_assert((char *)"self_dimer error, self_dimer_score for seq2 != 18", d.self_dimer_score == 18);
  return 0;
}

char * test_pair_dimer() {
  dimerisation d;
  d.pair_dimer(test_seq1,test_seq3);
  mu_assert((char *)"pair_dimer error, forward_dimer_score != 12", d.forward_dimer_score == 12);
  return 0;
}

