#include "test_primer.h" 
#include "../pd5/primer.h"

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




char * test_generate_100_candidate_primers() {
  primer p;
  p.start_location_range_begin = 1;
  p.start_location_range_end   = 50;
  p.max_number_candidates      = 100;
  char *seqName,*seqP;
  char filename[] = "test/example_genome.fa";
  readTemplate(filename,&seqP,&seqName);
  p.generate_candidates(seqP);
  mu_assert((char *)"generate_100_candidate_primers error, candidates_found != 100", p.candidates_found == 100);
  free(seqP);
  free(seqName);
  return 0;
}

char * test_generate_20mer_primers_at_start() {
  primer p;
  p.start_location_range_begin = 1;
  p.start_location_range_end   = 10;
  p.length_range_longest       = 20;
  p.length_range_shortest      = 20;
  p.max_number_candidates      = 100;
  char *seqName,*seqP;
  char filename[] = "test/example_genome.fa";
  readTemplate(filename,&seqP,&seqName);
  p.generate_candidates(seqP);
  mu_assert((char *)"generate_fewer_candidate_primers error, candidates_found == 5", p.candidates_found == 5);
  free(seqP);
  free(seqName);
  return 0;
}

