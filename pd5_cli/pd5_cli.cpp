/***********************************************************************\
 
 primer_design_cli.cpp 
 
 Created by Michael C. Riley on 28/07/2010.
 Edited by Amanda Clare, 2011.
 *  Copyright 2011 Aberystwyth University. All rights reserved.
 
\***********************************************************************/

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <vector>
//#include "../pd5/primer.h"
#include "../pd5/primer_pair.h"
#include "../pd5/display_utils.h"
//#include "../pd5/sequence_utils.h"
//#include "../pd5/DNAfind.h"

using namespace std;

#define TRUE 1
#define FALSE 0


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

/**
 * OutputType allows the choice of plain text or HTML
 */
enum OutputType { TXT, HTML, CSV };



/** 
 * Describes the command line arguments available. Outputs to stdout.
 */
void usage() {
  std::cout << "primer_design_cli [OPTION]" << std::endl << std::endl;
  std::cout << "Find primers in a given target sequence. Options are as follows:" << std::endl;
  std::cout << "-h, --help" << std::endl;
  std::cout << "      Display this message." << std::endl;
  std::cout << "-f, --fwd fwd_args" << std::endl;
  std::cout << "      fwd_args should be a comma-spearated string of 4 integers, representing the target region and primer sizes: min_target_start,max_target_stop,min_primer_length,max_primer_length. For example \"1,500,20,40\"." << std::endl;
  std::cout << "-r, --rev rev_args" << std::endl;
  std::cout << "      rev_args should be a comma-spearated string of 4 integers, representing the target region and primer sizes: min_target_start,max_target_stop,min_primer_length,max_primer_length. For example \"1000,900,20,40\"." << std::endl;
  std::cout << "-s, --seq sequence" << std::endl;
  std::cout << "      Template sequence in which to find primers" << std::endl;
  std::cout << "-S, --Seqfile filename" << std::endl;
  std::cout << "      Fasta file containing template sequence in which to find primers" << std::endl;
  std::cout << "-n, --nsbfile filename" << std::endl;
  std::cout << "      Name of fasta file against which to search for non-specific binding." << std::endl;
  std::cout << "-o, --output t|h|c" << std::endl;
  std::cout << "      Output type: text, HTML or CSV" << std::endl;
  exit(0);
}



/**
 * Tokenises a string such as "1,20,30,60", splitting on commas, and
 * converting component parts to ints. Used to split command line
 * arguments supplied for forward and reverse primers.
 */
int splitargs(char *str, vector<int> &args) {
  char * tok = strtok (str,",");
  while (tok != NULL)
  {
    int num = atoi(tok); 
    args.push_back(num);
    tok = strtok (NULL, ",");
  }
  return TRUE;
}



/**
 * Reads template sequence from a fasta-formatted file.
 */
int readTemplate(char* seqFile,char** seqP, char** seqName){ 
  ifstream fin(seqFile);
  char buffer[256];
  *seqP = new char[7000];
  *seqName = new char[100];	  
  while(!fin.getline(buffer, 255).eof()) {
    if(buffer[0] == 0x3E)  // fasta header
      {
	strcpy(*seqName, strtok(buffer,"> \t\n"));
      }
    else {
      strcat(*seqP, buffer);
    }
    fin.clear();
  }
  return 1;
}




/**
 * Takes many arguments.
 * @param seq the template sequence
 * @param seqFile a filename of a file containing the template sequence
 * @param fwdargs_s the single string representing the parameters for the forward %primer
 * @param revargs_s the single string representing the parameters for the reverse %primer
 * @param output_type the type of output required (text or html, see OutputType for possible values)
 * @param nsbfile the filename of a fasta file in which to search for non-specific binding
 * \sa OutputType
 */
int new_process(char* template_sequence,
		char *seqFile, 
		char *fwdargs_s, 
		char *revargs_s, 
		OutputType output_type, 
		char* nsbfile)
{

  char *seqName = NULL;
  DNAfind *nsbP = NULL;
  
  cout << "New process\n";

  // Check that we have a genome file and, if so, set up DNAfind
  // and parameters for secondary binding detection

  if(nsbfile != NULL) 
  {
    nsbP = new DNAfind(nsbfile);
    //nsbP = &nsb;
    if(!nsbP->set_max_mismatches(0)) 
    {
      cout << "Could not set max_mismatches" << std::endl;
      return FALSE;
    }

    if(!nsbP->set_tail_length(12)) 
    {
      cout << "Could not set tail length" << std::endl;
      return FALSE;
    }

    if(!nsbP->set_max_viable_product_length(5000)) 
    {
      cout << "Could not set max amplicon length" << std::endl;
      return FALSE;
    }
  }

  /* Test and get the query sequence */

  if(seqFile != NULL) 
  {
    readTemplate(seqFile,&template_sequence, &seqName);
  }

  /* If all ok, make a pair of primers */

  if(template_sequence != NULL)
  {
    primer_pair pcr1;

    // Set parameters
    pcr1.set_target_location(200, 500);
    pcr1.set_primer_length_range(20, 20);

    // Get candidate primers
    pcr1.generate_candidates(template_sequence);

    // Analyse candidates and display
    pcr1.candidate_analysis();
    //pcr1.show_individual_candidates();

    // Sort individual candidates and display
    pcr1.sort_individual_candidates("HAIRPIN, SELF_DIMER, TEMPERATURE");
	
// Testing:
    //pcr1.show_individual_candidates();

    // Select and sort primer pairs
    pcr1.sort_pair_candidates("TM_DIFF, F_DIMER, R_DIMER, MOO_SORT");

    // For testing only: Display best 6 candidate pairs
    pcr1.show_best_pair_candidates(6);

	int best = 0; // 0 is the index for the best candidate

    display_utils display;
    //char * product = NULL;
	char product[] = "AGTCGTCGAGCTCGATGCTAGCTCGATCGAT";

    //display.extract_product(template_sequence, &fwd.candidate[bestFwd], &rev.candidate[bestRev], product);

	
    switch (output_type) 
    {

    case TXT: {
      // Forward
      if(seqName) {
	      cout << "Name: " << seqName << endl;
      }
      cout << "Forward primer" << endl;
      cout << "Sequence: " << pcr1.pair_candidate[best].forward_sequence << endl;
      cout << "Size: " << strlen(pcr1.pair_candidate[best].forward_sequence) << endl;
      cout << "Location 5',3': " << pcr1.pair_candidate[best].location_forward_5_prime_end + 1 << ","
           << pcr1.pair_candidate[best].location_forward_5_prime_end + strlen(pcr1.pair_candidate[best].forward_sequence) << endl;
      cout << "Hairpin score: " << pcr1.pair_candidate[best].forward_hairpin_score << endl;
      cout << "Self dimer score: " << pcr1.pair_candidate[best].forward_self_dimer_score << endl;
      cout << "Temperature: " << pcr1.pair_candidate[best].forward_annealing_temperature << endl;
      // Reverse
      cout << "Reverse primer" << endl;
      cout << "Sequence: " << pcr1.pair_candidate[best].reverse_sequence  << endl;
      cout << "Size: " << strlen(pcr1.pair_candidate[best].reverse_sequence)  << endl;
      cout << "Location 5',3': " << pcr1.pair_candidate[best].location_reverse_5_prime_end + 1 << ","
           << pcr1.pair_candidate[best].location_reverse_5_prime_end + 1 - strlen(pcr1.pair_candidate[best].reverse_sequence) << endl;
      cout << "Hairpin score: " << pcr1.pair_candidate[best].reverse_hairpin_score  << endl;
      cout << "Self dimer score: " << pcr1.pair_candidate[best].reverse_self_dimer_score << endl;
      cout << "Temperature: " << pcr1.pair_candidate[best].reverse_annealing_temperature << endl;

      // Primer dimer
      cout << "Pair details" << endl;
      cout << "Primer dimer score (fwd): " << pcr1.pair_candidate[best].forward_pair_dimer_score  << endl;
      cout << "Primer dimer score (rev): " << pcr1.pair_candidate[best].reverse_pair_dimer_score << endl;
      cout << "NSB score: " << pcr1.pair_candidate[best].number_of_pcr_products << endl;

      // Product
      cout << "Product length: " << strlen(product) << endl;
      cout << "Product: " << product << endl;
      break;
    }
    case HTML: {
      if(seqName) {
	std::cout << "<h3>Name</h3><p>" << seqName << "</p>" << std::endl;
      }
      // Forward
      /*std::cout << "<h3>Forward primer</h3>" << std::endl;
      std::cout << "<p>Sequence: " << pcr1.forward.candidate[bestFwd].sequence << "<br />" << std::endl;
      std::cout << "Size: " << strlen(pcr1.forward.candidate[bestFwd].sequence) << "<br />" << std::endl;
      std::cout << "Location 5': " << pcr1.forward.candidate[bestFwd].location_5_prime_end + 1 << "<br />" << std::endl;
      std::cout << "Hairpin score: " << pcr1.forward.candidate[bestFwd].hairpin << "<br />" << std::endl;
      std::cout << "Self dimer score: " << pcr1.forward.candidate[bestFwd].self_dimer << "<br />" << std::endl;
      std::cout << "Temperature: " << pcr1.forward.candidate[bestFwd].annealing_temperature << "</p>" << std::endl;

      // Reverse
      std::cout << "<h3>Reverse primer</h3>" << std::endl;
      std::cout << "<p>Sequence: " << pcr1.reverse.candidate[bestRev].sequence << "<br />" << std::endl;
      std::cout << "Size: " << strlen(pcr1.reverse.candidate[bestRev].sequence) << "<br />" << std::endl;
      std::cout << "Location 5': " << pcr1.reverse.candidate[bestRev].location_5_prime_end + 1 << "<br />" << std::endl;
      std::cout << "Hairpin score: " << pcr1.reverse.candidate[bestRev].hairpin << "<br />" << std::endl;
      std::cout << "Self dimer score: " << pcr1.reverse.candidate[bestRev].self_dimer << "<br />" << std::endl;
      std::cout << "Temperature: " << pcr1.reverse.candidate[bestRev].annealing_temperature << "</p>" << std::endl;

      // Primer dimer
      std::cout << "<h3>Pair details</h3>" << std::endl;
      std::cout << "<p>Primer dimer score (fwd): " << pcr1.forward.candidate[bestFwd].forward_dimer  << "<br />" << std::endl;
      std::cout << "Primer dimer score (rev): " << pcr1.forward.candidate[bestFwd].reverse_dimer << "<br />" << std::endl;
      std::cout << "NSB score: " << pcr1.forward.candidate[bestFwd].seqsim_matches << "</p>" << std::endl;

      // Product
      std::cout << "<p>Product length: " << strlen(product) << "<br />" << std::endl;
      std::cout << "Product: " << product << "</p>" << std::endl;

      // Location
      char * fancy_template = NULL;
      display.html_colour_sequence(seq, &forward.candidate[bestFwd], &rev.candidate[bestRev], fancy_template);
      //std::cout << "<p>Template: " << seq << "</p>" << std::endl;
      std::cout << "<p>Locations of primers within template: " << fancy_template << "</p>" << std::endl;
      free(fancy_template);*/
      break;
    }
    case CSV: {
      //std::cout << "Fwd sequence, Fwd size, Fwd loc 5', Fwd loc 3', Fwd hairpin, Fwd self dimer, Fwd temp, Rev sequence, Rev size, Rev loc 5', Rev loc 3', Rev hairpin, Rev self dimer, Rev temp, Pair dimer (fwd), Pair dimer (rev), Pair NSB, Product len, Product" << endl;
      if(seqName) {
	std::cout << seqName << ",";
      }
      // Forward
      /*std::cout << pcr1.forward.candidate[bestFwd].sequence << ",";
      std::cout << strlen(pcr1.forward.candidate[bestFwd].sequence) << ",";
      std::cout << pcr1.forward.candidate[bestFwd].location_5_prime_end + 1 << "," << fwd.candidate[bestFwd].location_5_prime_end + strlen(fwd.candidate[bestFwd].sequence) << ",";
      std::cout << pcr1.forward.candidate[bestFwd].hairpin << ",";
      std::cout << pcr1.forward.candidate[bestFwd].self_dimer << ",";
      std::cout << pcr1.forward.candidate[bestFwd].annealing_temperature << ",";
      // Reverse
      std::cout << revseq  << ",";
      std::cout << strlen(revseq)  << ",";
      std::cout << pcr1.reverse.candidate[bestRev].location_5_prime_end + 1 << "," << rev.candidate[bestRev].location_5_prime_end + 1 - strlen(revseq)<< ",";
      std::cout << pcr1.reverse.candidate[bestRev].hairpin  << ",";
      std::cout << pcr1.reverse.candidate[bestRev].self_dimer << ",";
      std::cout << pcr1.reverse.candidate[bestRev].annealing_temperature << ",";

      // Primer dimer
      std::cout << pcr1.forward.candidate[bestFwd].forward_dimer  << ",";
      std::cout << pcr1.forward.candidate[bestFwd].reverse_dimer << ",";
      std::cout << pcr1.forward.candidate[bestFwd].seqsim_matches << ",";

      // Product
      std::cout << strlen(product) << ",";
      std::cout << product << std::endl;*/

      break;
    }
    } 
    //free(product);
  }

  return(TRUE);

}




/**
 * Main method for primer_design_cli. 
 * Deals with basic command line argument processing.
 */
 
int main(int argc, char** argv)
{
  int c;
  char *seq                = NULL;
  char *seqFile            = NULL;
  char *revargs            = NULL;
  char *fwdargs            = NULL;
  char *output_type_string = NULL;
  char *nsbfile            = NULL;
  OutputType output_type   = TXT; 

  //cout << "Args: " << argc << endl;

  while(1)
    {
      static struct option long_options[] =
	{
	  {"help",     no_argument,   0, 'h'},
	  {"output",   no_argument,   0, 'o'},
	  {"fwd",  required_argument, 0, 'f'},
	  {"rev",  required_argument, 0, 'r'},
	  {"seq",  required_argument, 0, 's'},
	  {"Seqfile",  required_argument, 0, 'S'},
	  {"nsbfile",  required_argument, 0, 'n'},
	  {0, 0, 0, 0}
	};
      /* getopt_long stores the option index here. */
      int option_index = 0;
      
      c = getopt_long (argc, argv, "ho:f:r:s:S:n:",
		       long_options, &option_index);
     
      /* Detect the end of the options. */
      if (c == -1)
	break;
      
      switch (c) {
      case 0:
	/* If this option set a flag, do nothing else now. */
	if (long_options[option_index].flag != 0)
	  break;
	printf ("option %s", long_options[option_index].name);
	if (optarg)
	  printf (" with arg %s", optarg);
	printf ("\n");
	break;
	 
      case 'h':
	usage();
	break;
	
      case 'f':
	fwdargs = optarg;
	break;
	
      case 'r':
	revargs = optarg;
	break;
	
      case 's':
	seq = optarg;
	seqFile = NULL;
	break;

      case 'S':
	seq = NULL;
	seqFile = optarg;
	break;

      case 'n':
	nsbfile = optarg;
	break;

      case 'o':
	output_type_string = optarg;
	if( output_type_string[0] == 'h') {
	  output_type = HTML;
	} else if (output_type_string[0] == 'c') {
	  output_type = CSV;
	}
	break;
	
      case '?':
	/* getopt_long already printed an error message. */
	break;
	
      default:
	abort();
      }
    }
  new_process(seq, seqFile, fwdargs, revargs, output_type, nsbfile);

  exit(1);
}





