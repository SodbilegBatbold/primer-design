#include <string>

#include "test_sequence_utils.h" 
#include "../pd5/sequence_utils.h"

static char test_seq1[] = "ATTGCGTATTATGC";
static char test_seq2[] = "ATTGCGTATTATAT";

char * test_tail_check() {
  ofstream fout;
  fout.open("temp.txt");
  int hits = sequence_utils::tail_check(test_seq1,test_seq2,fout);
  fout.close();
  mu_assert((char *)"test_tail_check error, hits != 3", hits == 3);
  return 0;
}

char * test_sticky_tail_check1() {
  int ok = sequence_utils::sticky_tail_check(test_seq1);
  mu_assert((char *)"sticky_tail_check test_seq1 error, test_seq1 should have sticky tail", ok == TRUE);
  return 0;
}

char * test_sticky_tail_check2() {
  int ok = sequence_utils::sticky_tail_check(test_seq2);
  mu_assert((char *)"sticky_tail_check test_seq2 error, test_seq2 should not have sticky tail", ok == FALSE);
  return 0;
}


char * test_reverse_complement() {
  char rev[15];
  string correct("GCATAATACGCAAT");
  sequence_utils::reverse_complement(test_seq1,rev);
  mu_assert((char *)"reverse_complement error, test_seq1 was not revcomped correctly", correct.compare(rev) == 0);
  return 0;
}

char * test_nucleotide_content() {
  int nc = sequence_utils::nucleotide_content(GUANINE, test_seq1);
  mu_assert((char *)"nucleotide_content error, nc != 3", nc == 3);
  return 0;
}

char * test_nucleotide_complement() {
  int nc = sequence_utils::nucleotide_complement(GUANINE);
  mu_assert((char *)"nucleotide_complement error, nc != CYTOSINE", nc == CYTOSINE);
  return 0;
}


