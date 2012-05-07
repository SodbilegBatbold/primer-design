#include <cmath>
#include "../display_utils.h"
#include "sputnik_ssr.h"
#include "../sequence_utils.h"
#include "../annealing_temperature.h"
#include "sputnik.h"

/**
 * lbound(a,b) is a macro that returns a, using b as a lower bound to
 * limit the return value
 */
#define lbound(a, b) (((a) < (b)) ? (b) : (a))
/**
 * ubound(a,b) is a macro that returns a, using b as an upper bound to
 * limit the return value
 */
#define ubound(a, b) (((a) > (b)) ? (b) : (a))

#define MAX_TEMP_DIFF 2

sputnik_ssr::sputnik_ssr(){
  fwdP = new primer;
  revP = new primer;
}

sputnik_ssr::~sputnik_ssr(){
    delete(fwdP);
    delete(revP);
    fwdP = NULL;
    revP = NULL;
}

int sputnik_ssr::has_repeat(char *seq) {
  sputnik sp;
  sp.Min_Unit_Length = 2; 
  sp.Max_Unit_Length = 5;
  sp.Exact_Match_Points = 1;
  sp.Error_Match_Points = -6;
  sp.Match_Min_Score = 12;
  sp.findingPrimers = False;
  sputnik::SeqStruct seq_struct;
  strcpy(seq_struct.seqStr, seq);
  strcpy(seq_struct.descStr, "dummy");
  seq_struct.seqLen = strlen(seq);
  Boolean has_repeats = sp.findRepeats(&seq_struct);
  if(has_repeats) {
    return 1;
  } else {
    return 0;
  }
}


int sputnik_ssr::find_primers(char* seq, int start, int end, int repeat_len) {
  //int good_candidates = 0;

  // This pattern is the repeating motif and shouldn't be found in the primer
  char repeat[repeat_len + 1];
  for (int i=start; i < start + repeat_len; i++) {
    repeat[i - start] = seq[i];
  }
  repeat[repeat_len] = 0;

  // Set locations for search
  int fwd5begin = start - 90;
  int fwd5end   = start - 30;
  if(fwd5begin < 1) { return 0; }
  int rev5begin = end + 90;
  int rev5end = end + 30;
  if(rev5begin > (int) strlen(seq)) { return 0; }


  revP->reverse_primer                = TRUE;
  revP->start_location_range_begin    = rev5begin;
  revP->start_location_range_end      = rev5end;
  revP->length_range_shortest         = 20;
  revP->length_range_longest          = 30;
  revP->optimum_primer_length         = revP->length_range_shortest;
  revP->optimum_Tm = 57;
  revP->seq_to_avoid                  = repeat;
  int okrev = revP->generate_candidates(seq);
  if(! okrev) {
    return 0;
  } 
  // First we calculate the basic cheap properties and sort the candidates
  for(int i = 0; i < revP->candidates_found; i++) {
    revP->hairpin(i);
    revP->self_dimer(i);
    revP->candidate[i].seqsim_matches = revP->blast_seq(i, seq); 
  }
  revP->priority[0] = SELF_DIMER;
  revP->priority[1] = SEQSIM_MATCH;
  revP->priority[2] = LENGTH;
  revP->priority[3] = SORT_END;

  
  if(! revP->rank_selection()) {
    return 0;
  }

  // Next we calculate temperature for the good ones, and resort.
  for(int i = 0; i < revP->good_candidates; i++) {
    revP->calculate_temperature(i); 
  }
  revP->candidates_found = revP->good_candidates;
  
  revP->priority[0] = TEMPERATURE;
  revP->priority[1] = SORT_END;
  
  if(! revP->rank_selection()) {
    return 0;
  }
  int bestRev = 0;

  int r_primer_end = revP->candidate[bestRev].location_5_prime_end;
  char rev_flank[r_primer_end - (end + 1) + 1];
  strncpy(rev_flank, seq + end, r_primer_end - end+1);
  rev_flank[r_primer_end - end + 1] = '\0';
  if(has_repeat(rev_flank)) {
    return 0;
  }

  
  char *revseq = revP->candidate[bestRev].sequence;
  int len_rev = strlen(revseq);
    

  // set the properties we need for a forward primer
  fwdP->reverse_primer             = FALSE;
  fwdP->start_location_range_begin = fwd5begin;
  fwdP->start_location_range_end   = fwd5end;
  fwdP->length_range_shortest = lbound(len_rev, revP->length_range_shortest);
  fwdP->length_range_longest  = ubound(len_rev, revP->length_range_longest);
  fwdP->optimum_primer_length = fwdP->length_range_shortest;
  fwdP->optimum_Tm = revP->candidate[bestRev].annealing_temperature;
  fwdP->required_GC_content   = sequence_utils::GC_content(revseq);
  fwdP->GC_tolerance          = 0;
  fwdP->seq_to_avoid          = repeat;

  int okfwd = fwdP->generate_candidates(seq);
  if(! okfwd) {
    return 0;
  } 

  // First we calculate the basic cheap properties and sort the candidates
  for(int i = 0; i < fwdP->candidates_found; i++) {
    fwdP->hairpin(i);
    fwdP->self_dimer(i);
    fwdP->primer_dimer_2(i, revseq);
    fwdP->candidate[i].seqsim_matches = fwdP->blast_seq(i, seq); 
  }
  fwdP->priority[0] = SELF_DIMER;
  fwdP->priority[1] = F_DIMER;
  fwdP->priority[2] = R_DIMER;
  revP->priority[3] = SEQSIM_MATCH;
  fwdP->priority[4] = LENGTH;
  fwdP->priority[5] = SORT_END;
    
  if(! fwdP->rank_selection()) {
    return 0;
  }


  for(int i = 0; i < fwdP->good_candidates; i++) {
    fwdP->calculate_temperature(i); 
  }
  fwdP->candidates_found = fwdP->good_candidates;
  fwdP->priority[0] = TEMPERATURE;
  fwdP->priority[1] = SORT_END;  
  if(! fwdP->rank_selection()) {
    return 0;
  }

  int bestFwd = 0;

  int f_primer_end = fwdP->candidate[bestFwd].location_5_prime_end;
  char fwd_flank[start - f_primer_end + 1];
  strncpy(fwd_flank, seq + f_primer_end, start - f_primer_end);
  fwd_flank[start - f_primer_end] = '\0';
  if(has_repeat(fwd_flank)) {
    return 0;
  }

  display_utils display;
  char * product = NULL;
  display.extract_product(seq, &fwdP->candidate[bestFwd], &revP->candidate[bestRev], product);    
  if (sequence_utils::nucleotide_content(ANYNUCLEOTIDE,product) > 0) {
    //std::cout << "Product would contains Ns" << std::endl;
    return 0;
  }
  if(fabs(revP->candidate[bestRev].annealing_temperature - fwdP->candidate[bestFwd].annealing_temperature) > MAX_TEMP_DIFF) {
    return 0;
  } 


  // Forward
  std::cout << "-----------------------------------------------" << std::endl;
  std::cout << "Forward primer" << std::endl;
  std::cout << "Sequence: " << fwdP->candidate[bestFwd].sequence << std::endl;
  std::cout << "Size: " << strlen(fwdP->candidate[bestFwd].sequence) << std::endl;
  std::cout << "Fwdflank: "<< fwd_flank << std::endl;
  std::cout << "Location 5',3': " << fwdP->candidate[bestFwd].location_5_prime_end + 1 << "," << fwdP->candidate[bestFwd].location_5_prime_end + strlen(fwdP->candidate[bestFwd].sequence) << std::endl;
  std::cout << "Hairpin score:            " << fwdP->candidate[bestFwd].hairpin << std::endl;
  std::cout << "Self dimer score:         " << fwdP->candidate[bestFwd].self_dimer << std::endl;
  std::cout << "Seq sim score:            " << fwdP->candidate[bestFwd].seqsim_matches << std::endl;
  std::cout << "Temperature:              " << fwdP->candidate[bestFwd].annealing_temperature << std::endl;
  // Reverse
  std::cout << "Reverse primer" << std::endl;
  std::cout << "Sequence: " << revseq  << std::endl;
  std::cout << "Size: " << strlen(revseq)  << std::endl;
  std::cout << "Revflank: "<< rev_flank << std::endl;
  std::cout << "Location 5',3': " << revP->candidate[bestRev].location_5_prime_end + 1 << "," << revP->candidate[bestRev].location_5_prime_end + 1 - strlen(revseq)<< std::endl;
  std::cout << "Hairpin score:            " << revP->candidate[bestRev].hairpin  << std::endl;
  std::cout << "Self dimer score:         " << revP->candidate[bestRev].self_dimer << std::endl;
  std::cout << "Seq sim score:            " << revP->candidate[bestRev].seqsim_matches << std::endl;
  std::cout << "Temperature:              " << revP->candidate[bestRev].annealing_temperature << std::endl;

  // Primer dimer
  std::cout << "Pair details" << std::endl;
  std::cout << "Primer dimer score (fwd): " << fwdP->candidate[bestFwd].forward_dimer  << std::endl;
  std::cout << "Primer dimer score (rev): " << fwdP->candidate[bestFwd].reverse_dimer << std::endl;

  // Product
  std::cout << "Product length: " << strlen(product) << std::endl;
  std::cout << "Product: " << product << std::endl;
  free(product);

  std::cout << std::endl;
    

  
  return 1;

} 
