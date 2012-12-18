/***********************************************************************\
 
 primer_design_cli.cpp 
 
 Created by Michael C. Riley on 28/07/2010.
 Edited by Amanda Clare, 2011.
 Last edit by MCR 12/06/2012
 
 *  Copyright 2011, 2012 Aberystwyth University. All rights reserved.
 
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
  std::cout << "      fwd_args should be a comma-separated string of 5 integers, representing the target region and primer sizes: min_target_start,max_target_stop,min_primer_length,optimum_primer_length,max_primer_length. For example \"1,500,20,23,40\"." << std::endl;
  std::cout << "-r, --rev rev_args" << std::endl;
  std::cout << "      rev_args should be a comma-separated string of 5 integers, representing the target region and primer sizes: min_target_start,max_target_stop,min_primer_length,optimum_primer_length,max_primer_length. For example \"1000,900,20,23,40\"." << std::endl;
  std::cout << "-t, --tm tm_args" << std::endl;
  std::cout << "      tm_args should be a comma-separated string of 3 integers, representing the annealing temperature range: min_temperature,optimum_temperature,max_temperature. For example \"50,55,60\"." << std::endl;
  std::cout << "-s, --seq sequence" << std::endl;
  std::cout << "      Template sequence in which to find primers" << std::endl;
  std::cout << "-S, --Seqfile filename" << std::endl;
  std::cout << "      Fasta file containing template sequence in which to find primers" << std::endl;
  std::cout << "-n, --nsbfile filename" << std::endl;
  std::cout << "      Name of fasta file against which to search for non-specific binding." << std::endl;
  std::cout << "-o, --output filename" << std::endl;
  std::cout << "      Name of output file. If this option is not used, output will go to stdout." << std::endl;
  std::cout << "-O, --output-type t|h|c" << std::endl;
  std::cout << "      Output type: text, HTML or CSV" << std::endl;
  std::cout << "-V, --version" << std::endl;
  std::cout << "      Display version information" << std::endl;
  exit(0);
}

/** 
 * Describes the version information
 */
void version() {
  std::cout << "pd5_cli: Version 1.0.0" << std::endl;
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
int readTemplate(char* seqFile, char** seqP, char** seqName){ 
  ifstream fin(seqFile);
  char buffer[256];
  *seqP = new char[7000];
  *seqName = new char[100];
  
  while(fin.getline(buffer, 255)) 
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

int extract_product(const char* template_seq, primer_pair_data *myprimers, char * &product) 
{  
  int begin = myprimers->location_forward_5_prime_end;
  int end   = myprimers->location_reverse_5_prime_end;
  int length = 1 + end - begin;
  
  //cout << "extract product length " << length << endl;

  if (template_seq == 0 || strlen(template_seq) == 0 || (int)strlen(template_seq) < begin || (int)strlen(template_seq) < (begin + length - 1)) 
  {
    return 0;
  } 
  else 
  {
    product = (char *)(malloc (length + 1));
    for (int i = 0; i < length; i++)
	{
      product[i] = template_seq[begin + i];
    }
    product[length] = 0;
    return 1;
  }
}


/**
 * Takes many arguments.
 * @param seq the template sequence
 * @param seqFile a filename of a file containing the template sequence
 * @param fwdargs_s the single string representing the parameters for the forward %primer
 * @param revargs_s the single string representing the parameters for the reverse %primer
 * @param tm_args_s the single string representing the parameters for the annealing temperature
 * @param output_type the type of output required (text or html, see OutputType for possible values)
 * @param nsbfile the filename of a fasta file in which to search for non-specific binding
 * \sa OutputType
 */
int new_process(char* template_sequence,
		char *seqFile, 
		char *fwdargs_s, 
		char *revargs_s,
		char *tm_args_s,		
		OutputType output_type, 
		char* nsbFile,
		char* outFile)
{

  char *seqName = NULL;
  DNAfind *nsbP = NULL;
  
  //cout << "New process\n";

  // Check that we have a genome file and, if so, set up DNAfind
  // and parameters for secondary binding detection

  if(nsbFile != NULL) 
  {
    nsbP = new DNAfind(nsbFile);
    //nsbP = &nsb;
    if(!nsbP->set_max_mismatches(0)) 
    {
      cerr << "Could not set max_mismatches" << std::endl;
      return FALSE;
    }

    if(!nsbP->set_tail_length(12)) 
    {
      cerr << "Could not set tail length" << std::endl;
      return FALSE;
    }

    if(!nsbP->set_max_viable_product_length(5000)) 
    {
      cerr << "Could not set max amplicon length" << std::endl;
      return FALSE;
    }
  }

  /* Test and get the query sequence */

  if(seqFile != NULL) 
  {
    readTemplate(seqFile, &template_sequence, &seqName);
  }

  /* If all ok, make a pair of primers */

  if(template_sequence != NULL)
  {
    primer_pair pcr1;

	// Get and set parameters
	//pcr1.set_target_location(200, 500);
    //pcr1.set_primer_length_range(20, 20);
	
	if (revargs_s != NULL) 
	{
		std::vector<int> revargs;
		splitargs(revargs_s, revargs);
		
		pcr1.reverse_primer.set_primer_location_range(revargs[0]-1, revargs[1]-1);
		pcr1.reverse_primer.set_primer_length_range(revargs[2], revargs[4]); 
		pcr1.reverse_primer.optimum_primer_length = revargs[3];
    } 
	else 
	{
		pcr1.reverse_primer.set_primer_location_range(999, 999);
		pcr1.reverse_primer.set_primer_length_range(18, 30);
		pcr1.reverse_primer.optimum_primer_length = 20;
    }
	
	if (fwdargs_s != NULL) 
	{
		std::vector<int> fwdargs;
		splitargs(fwdargs_s, fwdargs);
	  
		pcr1.forward_primer.set_primer_location_range(fwdargs[0]-1, fwdargs[1]-1);
		pcr1.forward_primer.set_primer_length_range(fwdargs[2], fwdargs[4]);
		pcr1.forward_primer.optimum_primer_length = fwdargs[3];
    } 
	else 
	{
		pcr1.forward_primer.set_primer_location_range(0, 499);
		pcr1.forward_primer.set_primer_length_range(18, 30);
		pcr1.forward_primer.optimum_primer_length = 20;
    }
	
	if(tm_args_s != NULL)  // else uses default settings 50, 55, 60
	{
		vector<int> tm_args;
		splitargs(tm_args_s, tm_args);
		
		pcr1.set_Tm_range(tm_args[0], tm_args[1], tm_args[2]);	
	}
  
    //cout << "Get candidate primers\n";
    pcr1.generate_candidates(template_sequence);

    //cout << "Analyse candidates and display\n";
    pcr1.candidate_analysis();
    //pcr1.show_individual_candidates();

    //cout << "Sort individual candidates and display\n";
    pcr1.sort_individual_candidates("HAIRPIN, SELF_DIMER, TEMPERATURE");
	
    // Testing:
    //pcr1.show_individual_candidates();

    //cout << "Select and sort primer pairs\n";
    pcr1.sort_pair_candidates("TM_DIFF, F_DIMER, R_DIMER, MOO_SORT");

    // For testing only: Display best 6 candidate pairs
    //pcr1.show_best_pair_candidates(6);

    int best = 0; // 0 is the index for the best candidate primer pair

    display_utils display;
    char * product = NULL;
	//char product[] = "AGTCGTCGAGCTCGATGCTAGCTCGATCGAT";

    extract_product(template_sequence, &pcr1.pair_candidate[best], product);

    std::streambuf * buf;
    std::ofstream of;

    if(outFile != NULL) {
      of.open(outFile);
      buf = of.rdbuf();
    } else {
      buf = std::cout.rdbuf();
    }

    std::ostream out(buf);

    switch (output_type) 
    {

    case TXT: {
      // Forward
      if(seqName) {
	out << "Name: " << seqName << endl;
      }
      out << "Forward primer" << endl;
      out << "Sequence: " << pcr1.pair_candidate[best].forward_sequence << endl;
      out << "Size: " << strlen(pcr1.pair_candidate[best].forward_sequence) << endl;
      out << "Location 5',3': " << pcr1.pair_candidate[best].location_forward_5_prime_end + 1 << ","
           << pcr1.pair_candidate[best].location_forward_5_prime_end + strlen(pcr1.pair_candidate[best].forward_sequence) << endl;
      out << "Hairpin score: " << pcr1.pair_candidate[best].forward_hairpin_score << endl;
      out << "Self dimer score: " << pcr1.pair_candidate[best].forward_self_dimer_score << endl;
      out << "Temperature: " << pcr1.pair_candidate[best].forward_annealing_temperature << endl;
      // Reverse
      out << "Reverse primer" << endl;
      out << "Sequence: " << pcr1.pair_candidate[best].reverse_sequence  << endl;
      out << "Size: " << strlen(pcr1.pair_candidate[best].reverse_sequence)  << endl;
      out << "Location 5',3': " << pcr1.pair_candidate[best].location_reverse_5_prime_end + 1 << ","
           << pcr1.pair_candidate[best].location_reverse_5_prime_end + 1 - strlen(pcr1.pair_candidate[best].reverse_sequence) << endl;
      out << "Hairpin score: " << pcr1.pair_candidate[best].reverse_hairpin_score  << endl;
      out << "Self dimer score: " << pcr1.pair_candidate[best].reverse_self_dimer_score << endl;
      out << "Temperature: " << pcr1.pair_candidate[best].reverse_annealing_temperature << endl;

      // Primer dimer
      out << "Pair details" << endl;
      out << "Primer dimer score (fwd): " << pcr1.pair_candidate[best].forward_pair_dimer_score  << endl;
      out << "Primer dimer score (rev): " << pcr1.pair_candidate[best].reverse_pair_dimer_score << endl;
      out << "NSB score: " << pcr1.pair_candidate[best].number_of_pcr_products << endl;

      // Product
      out << "Product length: " << strlen(product) << endl;
      out << "Product: " << product << endl;
      break;
    }
    case HTML: {
      if(seqName) {
	out << "<h3>Name</h3><p>" << seqName << "</p>" << std::endl;
      }
      // Forward
      out << "<h3>Forward primer</h3>" << std::endl;
      out << "<p>Sequence: " << pcr1.pair_candidate[best].forward_sequence << "<br />" << std::endl;
      out << "Size: " << strlen(pcr1.pair_candidate[best].forward_sequence) << "<br />" << std::endl;
      out << "Location 5': " << pcr1.pair_candidate[best].location_forward_5_prime_end + 1 << "<br />" << std::endl;
      out << "Hairpin score: " << pcr1.pair_candidate[best].forward_hairpin_score << "<br />" << std::endl;
      out << "Self dimer score: " << pcr1.pair_candidate[best].forward_self_dimer_score << "<br />" << std::endl;
      out << "Temperature: " << pcr1.pair_candidate[best].forward_annealing_temperature << "</p>" << std::endl;

      // Reverse
      out << "<h3>Reverse primer</h3>" << std::endl;
      out << "<p>Sequence: " << pcr1.pair_candidate[best].reverse_sequence << "<br />" << std::endl;
      out << "Size: " << strlen(pcr1.pair_candidate[best].reverse_sequence) << "<br />" << std::endl;
      out << "Location 5': " << pcr1.pair_candidate[best].location_reverse_5_prime_end + 1 << "<br />" << std::endl;
      out << "Hairpin score: " << pcr1.pair_candidate[best].reverse_hairpin_score << "<br />" << std::endl;
      out << "Self dimer score: " << pcr1.pair_candidate[best].reverse_self_dimer_score << "<br />" << std::endl;
      out << "Temperature: " << pcr1.pair_candidate[best].reverse_annealing_temperature << "</p>" << std::endl;

      // Primer dimer
      out << "<h3>Pair details</h3>" << std::endl;
      out << "<p>Primer dimer score (fwd): " << pcr1.pair_candidate[best].forward_pair_dimer_score  << "<br />" << std::endl;
      out << "Primer dimer score (rev): " << pcr1.pair_candidate[best].reverse_pair_dimer_score << "<br />" << std::endl;
      //out << "NSB score: " << pcr1.pair_candidate[best].number_of_pcr_products << "</p>" << std::endl;

      // Product
      out << "<p>Product length: " << strlen(product) << "<br />" << std::endl;
      out << "Product: " << product << "</p>" << std::endl;

      // Location
      char * fancy_template = NULL;
      display.html_colour_sequence(template_sequence, &pcr1.pair_candidate[best], fancy_template);
      out << "<p>Template: " << template_sequence << "</p>" << std::endl;
      out << "<p>Locations of primers within template: " << fancy_template << "</p>" << std::endl;
      free(fancy_template);
      break;
    }
    case CSV: {
      out << "Fwd sequence, Fwd size, Fwd loc 5', Fwd loc 3', Fwd hairpin, Fwd self dimer, Fwd temp, Rev sequence, Rev size, Rev loc 5', Rev loc 3', Rev hairpin, Rev self dimer, Rev temp, Pair dimer (fwd), Pair dimer (rev), Pair NSB, Product len, Product" << endl;
			
      if(seqName) out << seqName << ",";
      
      // Forward
      out << pcr1.pair_candidate[best].forward_sequence << ",";
      out << strlen(pcr1.pair_candidate[best].forward_sequence) << ",";
      out << pcr1.pair_candidate[best].location_forward_5_prime_end + 1 << "," 
		<< pcr1.pair_candidate[best].location_forward_5_prime_end + strlen(pcr1.pair_candidate[best].forward_sequence) << ",";
      out << pcr1.pair_candidate[best].forward_hairpin_score << ",";
      out << pcr1.pair_candidate[best].forward_self_dimer_score << ",";
      out << pcr1.pair_candidate[best].forward_annealing_temperature << ",";
      // Reverse
      out << pcr1.pair_candidate[best].reverse_sequence  << ",";
      out << strlen(pcr1.pair_candidate[best].reverse_sequence)  << ",";
      out << pcr1.pair_candidate[best].location_reverse_5_prime_end + 1 << "," 
		<< pcr1.pair_candidate[best].location_reverse_5_prime_end + 1 - strlen(pcr1.pair_candidate[best].reverse_sequence)<< ",";
      out << pcr1.pair_candidate[best].reverse_hairpin_score  << ",";
      out << pcr1.pair_candidate[best].reverse_self_dimer_score << ",";
      out << pcr1.pair_candidate[best].reverse_annealing_temperature << ",";
      
      // Primer dimer
      out << pcr1.pair_candidate[best].forward_pair_dimer_score  << ",";
      out << pcr1.pair_candidate[best].reverse_pair_dimer_score << ",";
      out << pcr1.pair_candidate[best].number_of_pcr_products << ",";
      
      // Product
      out << strlen(product) << ",";
      out << product << std::endl;

      break;
     }
    }
	
    free(product);
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
  char *tm_args            = NULL;
  char *output_type_string = NULL;
  char *nsbFile            = NULL;
  char *outFile            = NULL;
  OutputType output_type   = TXT; 

  //cout << "Args: " << argc << endl;
	if(!argv[1])usage();

  while(1)
    {
      static struct option long_options[] =
	{
	  {"help",     no_argument,   0, 'h'},
	  {"output",   required_argument, 0, 'o'},
	  {"output-type", required_argument, 0, 'O'},
	  {"fwd",  required_argument, 0, 'f'},
	  {"rev",  required_argument, 0, 'r'},
	  {"tm",  required_argument, 0, 't'},
	  {"seq",  required_argument, 0, 's'},
	  {"Seqfile",  required_argument, 0, 'S'},
	  {"nsbfile",  required_argument, 0, 'n'},
	  {"version",  no_argument, 0, 'V'},
	  {0, 0, 0, 0}
	};
      /* getopt_long stores the option index here. */
      int option_index = 0;
      
      c = getopt_long (argc, argv, "ho:O:f:r:t:s:S:n:V",
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
	
      case 't':
	tm_args = optarg;
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
	nsbFile = optarg;
	break;

      case 'o':
	outFile = optarg;
	break;

      case 'O':
	output_type_string = optarg;
	if( output_type_string[0] == 'h') {
	  output_type = HTML;
	} else if (output_type_string[0] == 'c') {
	  output_type = CSV;
	} else if (output_type_string[0] == 't') {
	  output_type = TXT;
	} else {
	  cerr << "Error: Invalid output type (should be c, t or h), but was " << output_type_string << endl;
	  exit(1);
	}

	break;
	
      case 'V':
	version();
	break;

      case '?':
	/* getopt_long already printed an error message. */
	break;
	
      default:
	abort();
      }
    }
	
  int ok = new_process(seq, seqFile, fwdargs, revargs, tm_args, output_type, nsbFile, outFile);

  if (!ok) {
    exit(1);
  } else {
    exit(0);
  }
}





