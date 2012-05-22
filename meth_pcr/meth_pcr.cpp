 /*************************************************\
 
 meth_pcr.cpp - 8th April 2011
 
 \*************************************************/


/*
 A = 0x41 = 65
 T = 0x54 = 84
 G = 0x47 = 71
 C = 0x43 = 67
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include "../pd5/primer_pair.h"
#include "../pd5/DNAfind.h"
#include "../pd5/sequence_utils.h"

using namespace std;


// TESTING
#define TESTING

#define TRUE 1
#define FALSE 0

char* bisulphite_sequence_conversion(const char* sequence)
{
	int sequence_length = strlen(sequence);
	char* bisulphite_sequence = (char*) malloc(sequence_length + 1);
	
	for(int i = 0; i < sequence_length; i++)
	{
		if(sequence[i] == CYTOSINE && sequence[i + 1] != GUANINE)
			bisulphite_sequence[i] = THYMINE;
		else 
			bisulphite_sequence[i] = sequence[i];
	}
	
	return(bisulphite_sequence);
}

/****************************************************************************/

class FASTA
{
public:
	FASTA(){};
	~FASTA(){};
	int file_reader(const char* filename);
	void free_sequence(){free(sequence);};
	void free_header(){free(header);};
	char* header;
	char* sequence;
	
private:
	
};

/** This will get one sequence and one header from a FASTA format file. If there is more than one
 sequence then the contents of header and sequence will be those of the last sequence in the file.
 */

int FASTA::file_reader(const char* filename)
{
	ifstream sequence_fa(filename);
	char buffer[2056];
	bool FASTA_SEQUENCE = FALSE;
	
	if(!sequence_fa.is_open())
	{
		cout << "Unable to open " << filename << endl;
		return(0);
	}
	
	while(sequence_fa.getline(buffer, 2055))
	{
		if(buffer[0] == 0x3E)  // fasta header
		{
			header = (char*) malloc(strlen(buffer) + 1);
			strcpy(header, strtok(buffer,">\n"));
			FASTA_SEQUENCE = FALSE;
		}
		else if(FASTA_SEQUENCE)
		{
			sequence = (char*)realloc(sequence, strlen(sequence) + strlen(buffer) + 1);
			strcat(sequence, buffer);
		}
		else 
		{
			sequence = (char*) malloc(strlen(buffer) + 1);
			strcpy(sequence, buffer);
			FASTA_SEQUENCE = TRUE;
		}
	}
	return(1);
}


/****************************************************************************/


int main(int argc, char** argv)
{
	char* bisulphite_sequence;
	char filename[128];
	
   	//cout << "Methylation primer testing - 8/4/11\n";
	//cout << "Version 1.1\n";
	
	cout << argv[1] << endl << endl;
	
	if(!(argv[1]))
	{
		cout << "Use: exe <gene sequence file>\n";
		return(0);
	}
	else 
	{
		strcpy(filename, argv[1]);
	}
	
	/** Use FASTA object my_fa to read FASTA file given in first argument of command line */
	FASTA my_fa;	
	my_fa.file_reader(filename);

	/** Take the FASTA sequence in my_fa.sequence and convert it such that it had been bisulphite treated */	
	bisulphite_sequence = bisulphite_sequence_conversion(my_fa.sequence);
	
	// Free space
	my_fa.free_sequence();
	my_fa.free_header();
	
	/** Make an object of Helygen's primer_pair class and use it to find primers */
	primer_pair pcr1; 
	
	//get_pair(pcr1, 800, bisulphite_sequence);
	
	// Set parameters
	pcr1.set_target_location(800, 1100);
	pcr1.set_primer_length_range(20, 30);
	pcr1.set_flank_lengths(150, 150);
	pcr1.forward.optimum_primer_length = 25;
	pcr1.reverse.optimum_primer_length = 25;
	//set_primer3_method(tm_method_type method);
	
	if(!pcr1.set_Tm_range(50.0, 55.0, 60.0))
	{
		cout << "Unable to set Tm range\n";
		return(FALSE);
	}
	
	// Methylation specific params
	pcr1.forward.no_C_primer = TRUE;
	pcr1.reverse.no_G_primer = TRUE;
	
	// Get candidate primers (bisulphite_sequence is the dna template we wish to use)
	pcr1.generate_candidates(bisulphite_sequence);  
	
	// Analyse candidates and display
	
	// To compare with Methprimer we need to use Primer3's default settings, not the recommended settings.
	pcr1.forward.set_primer3_method((tm_method_type) 0);
	pcr1.forward.set_primer3_salt_correction_method((salt_correction_type) 0);
	pcr1.reverse.set_primer3_method((tm_method_type) 0);
	pcr1.reverse.set_primer3_salt_correction_method((salt_correction_type) 0);
	
	pcr1.candidate_analysis();
	//pcr1.show_individual_candidates();
	
	// Sort individual candidates and display
	pcr1.sort_individual_candidates("HAIRPIN, SELF_DIMER, LENGTH, TEMPERATURE");
	//pcr1.show_individual_candidates();
	
	// Select and sort primer pairs
	pcr1.sort_pair_candidates("TM_DIFF, F_DIMER, R_DIMER, MOO_SORT");
	
	// Display best 5 candidate pairs
	//pcr1.show_best_pair_candidates(5);
	
	// Output for LaTeX ********************
	
	cout << "\\scriptsize\n";
	cout << "\\begin{tabular}{@{} llrrrl @{}} % Column formatting, @{} suppresses leading/trailing space\n";
	cout << "No  &    Primer      &        Start  & Size  &  Tm  &        Sequence \\\\\n";
	
	std::cout.setf(ios::fixed);
	std::cout.setf(ios::showpoint);
	
	for(int x = 0; x < 5; x++)
	{
		std::cout.precision(1);
		std::cout << x + 1 << " & ";
		//std::cout << "No.\tPrimer\tHpin\tSdimer\tTm\tLocus\tLength\tSequence\n";
		

		std::cout << "Left primer & ";
		std::cout << pcr1.pair_candidate[x].location_forward_5_prime_end << " & ";
		std::cout << strlen(pcr1.pair_candidate[x].forward_sequence) << " & ";
		std::cout << pcr1.pair_candidate[x].forward_annealing_temperature << " & ";
		std::cout << pcr1.pair_candidate[x].forward_sequence << " \\\\";
		std::cout << std::endl;
		
		std::cout << " & ";
		std::cout << "Right primer & ";
		std::cout << pcr1.pair_candidate[x].location_reverse_5_prime_end << " & ";
		std::cout << strlen(pcr1.pair_candidate[x].reverse_sequence) << " & ";
		std::cout << pcr1.pair_candidate[x].reverse_annealing_temperature << " & ";
		std::cout << pcr1.pair_candidate[x].reverse_sequence << " \\\\";
		std::cout << std::endl;
	}
cout << "\\end{tabular}\n";
	cout << "\\normalsize\n";
	
	// End of Output for LaTeX ********************
	
	/*annealing_temperature myTm;
	cout << "Primer3 Query Tm = " << myTm.primer3_Tm("TAAATTATGTAATTGATAGAGGTAGAATG") << endl;
	cout << "Oligo Query Tm = " << myTm.oligo_Tm("TAAATTATGTAATTGATAGAGGTAGAATG") << endl;
	cout << "Wallace Query Tm = " << myTm.Wallace_method("TAAATTATGTAATTGATAGAGGTAGAATG") << endl;
	cout << "Freier Query Tm = " << myTm.Freier_method("TAAATTATGTAATTGATAGAGGTAGAATG") << endl;
	*/
	
    return(1);
}


