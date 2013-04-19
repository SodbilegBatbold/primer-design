#include "test_DNAfind.h" 
#include "../pd5/dna_find.h"

static char test_seq1[] = "ATAGGCGACGCTAA";
static char test_seq2[] = "TTACATCCATATT";

char * test_search_for_binding_sites() {
  DNAfind d("test/example_genome.fa");
  int n = d.search_for_binding_sites(test_seq1);
  mu_assert((char *)"search_for_binding_sites error, n != 1", n == 1);
  return 0;
}

char * test_search_for_pcr_products() {
  DNAfind d("test/example_genome.fa");
  int n = d.search_for_pcr_products(test_seq1,test_seq2);
  mu_assert((char *)"search_for_pcr_products error, n != 1", n == 1);
  return 0;
}

char * test_mismatch_search_for_pcr_products() {
  DNAfind d("test/example_genome.fa");
  d.set_max_mismatches(3);
  int n = d.search_for_pcr_products(test_seq1,test_seq2);
  mu_assert((char *)"mismatched search_for_pcr_products error, n != 3", n == 3);
  return 0;
}

