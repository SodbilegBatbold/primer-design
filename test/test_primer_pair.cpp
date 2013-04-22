#include "test_primer_pair.h" 
#include "../pd5/primer_pair.h"

/**
 * Reads template sequence from a fasta-formatted file.
 */
static int readTemplate(const char* seqFile, char** seqP, char** seqName){ 
  ifstream fin(seqFile);
  char buffer[10000];
  *seqP = new char[10000];
  *seqName = new char[100];
  *seqP[0] = '\0';

  while(fin.getline(buffer, 9999)) 
  {
    if(buffer[0] == 0x3E)  // fasta header
    {
                strcpy(*seqName, strtok(buffer,"> \t\n\r"));
    }
    else 
        {
                strcat(*seqP, strtok(buffer,"> \t\n\r"));
    }
        
    fin.clear();
  }
  return 1;
}




char * test_generate_38_candidate_pairs() {
  char* seqP,*seqName;
  readTemplate("test/example_genome.fa", &seqP, &seqName);
  primer_pair pp;
  pp.forward_primer.start_location_range_begin = 1;
  pp.forward_primer.start_location_range_end   = 50;
  pp.reverse_primer.start_location_range_begin = 1000;
  pp.reverse_primer.start_location_range_end   = 950;
  pp.set_primer_length_range(18,25);
  pp.generate_candidates(seqP);
  pp.make_pair_candidates(pp.forward_primer.candidates_found, pp.reverse_primer.candidates_found);
  mu_assert((char *)"generate_100_candidate_pairs error, number_of_candidates != 38", pp.number_of_pair_candidates == 38);
  free(seqP);
  free(seqName);
  return 0;
}


char * test_generate_sort_candidate_pairs() {
  char* seqP,*seqName;
  readTemplate("test/example_genome.fa", &seqP, &seqName);

  DNAfind nsb("test/example_genome.fa");
  nsb.set_max_mismatches(0); 
  nsb.set_tail_length(12);
  nsb.set_max_viable_product_length(5000); 

  primer_pair pp;
  pp.forward_primer.start_location_range_begin = 1;
  pp.forward_primer.start_location_range_end   = 50;
  pp.reverse_primer.start_location_range_begin = 1000;
  pp.reverse_primer.start_location_range_end   = 950;
  pp.set_primer_length_range(18,25);
  pp.nsbP = &nsb;
  pp.generate_candidates(seqP);
  pp.candidate_analysis();
  pp.sort_individual_candidates("HAIRPIN,LENGTH,TEMPERATURE");
  pp.sort_pair_candidates("TM_DIFF,PRODUCTS");
  mu_assert((char *)"generate_sort_candidate_pairs error, good_pair_candidates != 36", pp.good_pair_candidates == 36);
  free(seqP);
  free(seqName);
  return 0;
}
