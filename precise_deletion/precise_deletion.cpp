
/*************************************************\
 
 precise_deletion.cpp (Ver 2.00 - 1/7/12)
 
 Updated from:
	precise_auto.cpp (Ver 1.00 - 10/2/11)
	s_primerII.cpp (Ver. 1.00 - 5/7/10) 
	sprimer.cpp (Ver. 1.00 - 7/6/10) to run on Beowulf
 
\*************************************************/

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <time.h>

using namespace std;

//#include "../pd5/primer_pair.h"
#include "../pd5/dna_find.h"
#include "../pd5/primer3_Tm_temperature.h"
//#include "../pd5/sequence_utils.h"

#include "../pd5/pd5.h"

// ORGANISM (POMBE, S_CERE, LACTIS)
#define POMBE

#define DB_OUT
#define DISPLAY_RESULTS

// OUTPUT FORMAT
#define XML 1

// Secondary binding search 
// BLAST, FASTA or NSB
#define NSB

// CLASSIFICATION OF RESULTS
#define PCR_RESULTS_CLASS "A"


// Error defines
#define OK 1
//#define ERROR 0
#define NO_PRIMER_A -1
#define NO_PRIMER_B -2
#define NO_PRIMERS_A_AND_B -9
#define NO_PRIMER_D -3
#define NO_PRIMER_E -4
#define NO_PCR2_FORWARD -5
#define NO_PCR2_REVERSE -6
#define NO_URA3 -7
#define CONF_PRIMER_ERROR -8
#define NO_PCR1_PRIMERS -10
#define NO_PLASMID_PRIMERS -11
#define NO_PCR2_PRIMERS -12

// TESTING
#define TESTING
#define Todo "To Do"

#define BEOWULF // To include/exclude changes required for distributed processing

#define TRUE 1
#define FALSE 0

//#define EXTERNAL_D_AND_E  // For D & E provided in a seperate file in FASTA format
#ifndef EXTERNAL_D_AND_E  // Otherwise, D & E are defined here

// Define fixed primers D and E here
#define PRIMER_D "GATCCCAATACAACAGAT" // This the original @ 97
//#define PRIMER_D "GTCTAGAGATCCCAATACA"
//#define PRIMER_D "CTAGAGATCCCAATACAAC" // This looks good for YGR125W (17/08/10)
//#define PRIMER_D "
//#define PRIMER_D "

#define PRIMER_E "CACACATTACTTGCCTC" // original
//#define PRIMER_E "AGCGACAAGAAGAGATAG" // (A)
//#define PRIMER_E "ACAATCATATGGGAGAAG" // (B)
//#define PRIMER_E "ACACTCCTCAGAAGCTC" // (C) @ 1443 Good for YGL202W (17/08/10)

#endif

// Cassette confirmation primers
/*
Lactis		Hpin	Sdim	BindA	BindB	Loc		Length	Tm		Sequence
FinCass		5.0		6.0		0		1		3107	29		54.6	TTCATCCTAAACCAAAAGTAAACAGTGTC
RinCass		7.0		9.0		0		1		3026	26		54.9	AAGGTACGCTTGTAGAATCCTTCTTC
 
 Pombe		Hpin	Sdim	BindA	BindB	Loc		Length	Tm		Sequence
 FinCass	8.0		8.0		1		1		5284	26		55.3	TACCGTCAAGCTACAATATGCATCTG
 RinCass	9.0		9.0		1		1		5273	23		55.3	CTGCGAATTTGCGATCCTCAAAG
 */

char m_complement[22] =  {"T G   C      N     A"};
int base_complement(int nucleotide)
{
	return(*(m_complement - 65 + nucleotide));	
}

int reverse_complement(const char* a_string, char* b_string)
{
	unsigned int i;
	for(i = 0; i < strlen(a_string); i++)
		b_string[i] = base_complement(a_string[strlen(a_string) - i - 1]);
	
	b_string[i] = 0;
	return TRUE;
}

/*
 *	Finds the location of a_string in b_string
 */

int findstring(const char* query_string, const char* main_string)
{
	unsigned int i, count;
	int location[30];
	bool found = FALSE;
	unsigned int main_string_length = strlen(main_string);
	unsigned int query_string_length = strlen(query_string);
	
	count = 0;
	
	for(i = 0; i < (main_string_length - query_string_length); i++)
	{
		found = TRUE;
		for(unsigned int j = 0; j < query_string_length; j++)
		{
			if(query_string[j] != main_string[i + j])
			{
				found = FALSE; 
				break;
			}
		}
		
		if(found)
		{
			location[count++] = i;
			cout << "Found " << query_string << " at " << i << endl;
			if(count >= 30)break;
		}
	}
	
	if(count == 0)
	{
		cout << "Could not find " << query_string << endl;
		return(-2);
	}
	else if(count == 1)
		return(location[0]);
	else 
	{
		cout << "Multiple hits for " << query_string << endl;
		return(-1);
	}
	
}


int start_codon_location(const char* orf)
{
	int start_location = 0;
	char start_codon[16] = {"ATG"};
	
	if(orf[1000] == start_codon[0] && orf[1001] == start_codon[1] && orf[1002] == start_codon[2])
	{
#ifdef TESTING
		cout << "Start location good\n";
#endif
		return(1000);
	}
	
	for(int i = 995; i < 1005; i++)
	{
	   	if(orf[i] == start_codon[0] && orf[i+1] == start_codon[1] && orf[i+2] == start_codon[2])
		{
			start_location = i;
			break;
		}
	}
	
	return(start_location);	
}

/**
	prokaryote start codons can be ATG, GTG, TTG
*/

int start_codon_location_for_prokaryotes(const char* orf)
{
	int start_location = 0;
	char start_codon[16] = {"ATG"};
	
	if(orf[1000] == start_codon[0] || orf[1000] == start_codon[1] || orf[1000] == start_codon[2])
	{
		if(orf[1001] == start_codon[1] && orf[1002] == start_codon[2])
		{
#ifdef TESTING
			cout << "Start location good\n";
#endif
			return(1000);
		}
	}
	
	for(int i = 995; i < 1005; i++)
	{
	   	if(orf[1000] == start_codon[0] || orf[1000] == start_codon[1] || orf[1000] == start_codon[2])
		{
			if(orf[1001] == start_codon[1] && orf[1002] == start_codon[2])
			{
				start_location = i;
				break;
			}
		}
	}
	
	return(start_location);	
}

int stop_codon_location(int start_codon_location, const char* orf)
{
	int stop_codon_location = 0;
	bool fnd = 0;
	
	// Find the stop codon
	
	for(unsigned int i = start_codon_location; i < strlen(orf); i += 3)
	{   if(orf[i] == 84)
		{
			if(orf[i+1] == 65 && orf[i+2] == 71){stop_codon_location = i; fnd = 1; break;}
			else if(orf[i+1] == 65 && orf[i+2] == 65){stop_codon_location = i; fnd = 1; break;}
			else if(orf[i+1] == 71 && orf[i+2] == 65){stop_codon_location = i; fnd = 1; break;}
		}
   	}
	
	if(fnd)
		return(stop_codon_location);
	else 
		return(0);
}

int check_sequence_contents(const char* sequence)
{
	long length = strlen(sequence);
	int letter[256];
	long letter_count[256];
	int tail = 0;
	int i, j;
	bool found = FALSE;
	
	for(i = 0; i < length ; i++)
	{
		found = FALSE;
		
		for(j = 0; j < tail; j++)
		{
			if((int)sequence[i] == letter[j])
			{
				letter_count[j]++;
				found = TRUE;
			}			
		}
		if(!found)
		{
			letter[tail] = (int)sequence[i];
			letter_count[tail] = 1;
			tail++;
		}
	}
	
	for(j = 0; j < tail; j++)
	{
		cout << (char)letter[j] << " = " << letter_count[j] << endl;
	}
		
	return(1);
}

/*
 The Whitehead Inst. uses the following parameters for their confirmation primers used with S. cerevisiae: -
 
 Primer length: 17 - 28 bases
 Primer Tm: 65 degrees C +/- 2 degrees C
 GC content: 30 - 70% (although Tm will have a large part in determining this)
 Primer A located 200 - 400 bases upstream of target ORF
 Primer B located downstream of primer A within the ORF to give a product of 300 - 1000 bases 
 Primer C located upstream of primer D within the ORF to give a product of 300 - 1000 bases
 Primer D located 200 - 400 bases downstream of target ORF
 
*/

#ifndef CONFIRM 

int precise_confirmation_primers(const char* orf_sequence,
								 const char* genome_file_name,
								 const char* plasmid_sequence,
								 const char* Primer_D,
								 const char* Primer_E,
								 const char* RC_Primer_C,
								 const char* Primer_R,
								 const char* Fwd_Cassette,
								 const char* Rev_Cassette,
								 char* cnf_results,
								 char* csv_results)
{
	int i = 0;
	char tempstr[15000]; // for message output cnf_results
	int orf_length = strlen(orf_sequence);
	char RC_orf_sequence[orf_length];
	int fwd_candidates, rev_candidates;
	char CPS_A[40];
	//char CPS_B[40];
	//char CPS_C[40];
	char CPS_D[40];
	char gene_present_confirmation_sequence[20000];
	char gene_deleted_confirmation_sequence[2000];
	
	cout << "Confirmation primers: \n";
	
	DNAfind genome_template(genome_file_name);
	
	if(!genome_template.set_max_mismatches(0)) 
		cout << "Could not set max_mismatches\n";
	
	if(!genome_template.set_tail_length(12)) 
		cout << "Could not set tail length\n";
	
	primer_pair confirmation; // This is for A and D, B and C are FinCass and RinCass
	
	confirmation.set_target_location(800, orf_length - 800);
	confirmation.set_flank_lengths(200, 200);
	confirmation.set_primer_length_range(17, 28);
	confirmation.set_Tm_range(50, 55, 60);
	
	confirmation.generate_candidates(orf_sequence);
	
	confirmation.candidate_analysis();
	confirmation.sort_individual_candidates("HAIRPIN, SELF_DIMER, TEMPERATURE");	
	
	if(confirmation.forward_primer.good_candidates > 5) 
		fwd_candidates = 6;
	else
		fwd_candidates = confirmation.forward_primer.good_candidates;
	
	if(confirmation.reverse_primer.good_candidates > 5)
		rev_candidates = 6;
	else
		rev_candidates = confirmation.reverse_primer.good_candidates;
	
	if(fwd_candidates < 1 || rev_candidates < 1)
		return(FALSE);
	
	cout << "Actual cands " << fwd_candidates << ", " << rev_candidates << endl;
	
	confirmation.make_pair_candidates(fwd_candidates, rev_candidates);
	
	for(int i = 0; i < confirmation.number_of_pair_candidates; i++)
	{
		//cout << pcr1.pair_candidate[i].forward_sequence << ", " << pcr1.pair_candidate[i].reverse_sequence << endl;
		confirmation.pair_candidate[i].number_of_pcr_products  =  
			genome_template.search_for_pcr_products(confirmation.pair_candidate[i].forward_sequence, 
													confirmation.pair_candidate[i].reverse_sequence);
		
		//cout << "Pair " << i << ": " << confirmation.pair_candidate[i].number_of_pcr_products << endl;
	}
	
	// sort pairs
	confirmation.sort_pair_candidates("TM_DIFF, F_DIMER, R_DIMER, PRODUCTS, MOO_SORT");
	
	strcpy(CPS_A, confirmation.pair_candidate[0].forward_sequence);
	strcpy(CPS_D, confirmation.pair_candidate[0].reverse_sequence);
	
	cout << "CPS_A: " << CPS_A << ", CPS_D: " << CPS_D << endl;
	
	//***************************************************
	
	strcpy(cnf_results,"\t\t\t<confirmationPrimerForward>");
	strcat(cnf_results, CPS_A);
	strcat(cnf_results, "</confirmationPrimerForward>\n");
	
	strcat(cnf_results, "\t\t\t<confirmationPrimerReverse>");
	strcat(cnf_results, CPS_D);
	strcat(cnf_results, "</confirmationPrimerReverse>\n");
	
	sprintf(tempstr, "%d, %d, ", (int)strlen(CPS_A), (int)strlen(CPS_D));
	strcpy(csv_results, tempstr);
	
	cout << "Length of ORF: " << strlen(orf_sequence) << endl;

	reverse_complement(orf_sequence, RC_orf_sequence);
	int forward_at = findstring(CPS_A, orf_sequence); //cout << "Forward strand\n";
	int reverse_at = findstring(CPS_D, RC_orf_sequence); //cout << "Reverse strand\n";
	
	cout << "Forward C primer at: " << forward_at << ", reverse C primer at: " << reverse_at << endl;
	
	if(forward_at > 1000){cout << "CPS_A location error\n"; return(0);}
	if(reverse_at > 1000){cout << "CPS_D location error\n"; return(0);}
	
	// FOR EUROFAN PRIMERs TEST
	// **************************************************
	
	//if(forward_at < -1 || reverse_at < -1) return(0);
	//else return(1);
	
	// **************************************************
	
	int A_to_start = 1000 - forward_at;
	int D_to_stop = 1000 - reverse_at;
	int gene_length = orf_length - 2000;
	
	// Make sequences
	int x = 0;
	int y = 0;

	
	for(i = forward_at; i < orf_length - reverse_at; i++)
	{
		gene_present_confirmation_sequence[x++] = orf_sequence[i];
		
		if(i < 1000 || i >= orf_length - 1000)
			gene_deleted_confirmation_sequence[y++] = orf_sequence[i];
	}
	gene_present_confirmation_sequence[x] = 0;
	gene_deleted_confirmation_sequence[y] = 0;
	
	/*	
	 fout << "CPS A to CPS D with gene present: " << A_to_start + D_to_stop + gene_length << " - " << gene_present_confirmation_sequence << endl;
	 fout << "CPS A to CPS D with gene deleted: " << A_to_start + D_to_stop << " - " << gene_deleted_confirmation_sequence << endl ;
	 */	
	
	
	sprintf(tempstr, "\t\t\t<genePresentLength>%d</genePresentLength>\n", A_to_start + D_to_stop + gene_length);
	strcat(cnf_results, tempstr);
	sprintf(tempstr, "\t\t\t<genePresentProduct>%s</genePresentProduct>\n", gene_present_confirmation_sequence);
	strcat(cnf_results, tempstr);
	
	sprintf(tempstr, "\t\t\t<geneDeletedLength>%d</geneDeletedLength>\n", A_to_start + D_to_stop);
	strcat(cnf_results, tempstr);
	sprintf(tempstr, "\t\t\t<geneDeletedProduct>%s</geneDeletedProduct>\n", gene_deleted_confirmation_sequence);
	strcat(cnf_results, tempstr);
	
	
	// Find the cassette between E and D
	char RC_plasmid_sequence[20000];
	reverse_complement(plasmid_sequence, RC_plasmid_sequence);
	
	int plasmid_length = strlen(plasmid_sequence);
		
	cout << "Looking for Primer D\n";
	int Primer_D_at = findstring(Primer_D, RC_plasmid_sequence); //cout << "Primer Forward strand\n";
	
	cout << "Looking for Primer E\n";
	int Primer_E_at = findstring(Primer_E, plasmid_sequence); //cout << "Primer Reverse strand\n";
	
	//cout << endl << ura3 << " Length is " << strlen(ura3) << endl;
	
	
	
	//cout << "Length of URA3: " << strlen(ura3) << endl;
	cout << "Looking for Primer RinKl\n";
	int Rinkl_at = findstring(Rev_Cassette, RC_plasmid_sequence); //cout << "Forward strand\n";
	cout << "Looking for Primer FinKl\n";
	int Finkl_at = findstring(Fwd_Cassette, plasmid_sequence); //cout << "Reverse strand\n";
	
	//cout << "Rev_Cassette: " << forward_at << endl ;
	//cout << "Fwd_Cassette: " << reverse_at << endl;
	
	int E_to_Rinkl = plasmid_length - Rinkl_at - Primer_E_at;
	int D_to_Finkl = plasmid_length - Finkl_at - Primer_D_at;
	
	//cout << "E to Rinkl: " << E_to_Rinkl << endl;
	//cout << "D to Finkl: " << D_to_Finkl << endl;
	
	// Get cassette sequence
	char cassette_sequence[20000];
	
	x = 0;
	
	for(i = Primer_E_at; i < (plasmid_length - Primer_D_at); i++) 
		cassette_sequence[x++] = RC_plasmid_sequence[i];
	cassette_sequence[x] = 0;
	
	// Get sequence from CPS A to the start of the gene
	char A_to_start_sequence[1000];
	
	x = 0;
	
	for(i = forward_at; i < 1000; i++) 
		A_to_start_sequence[x++] = orf_sequence[i];
	A_to_start_sequence[x] = 0;
	
	// Get the sequence from end of gene to CPS D
	char D_to_stop_sequence[1000];
	
	x = 0;
	
	for(i = orf_length - 1000; i < orf_length - reverse_at; i++) D_to_stop_sequence[x++] = orf_sequence[i];
	D_to_stop_sequence[x] = 0;
	
	/*	
	 fout << "CPS A to CPS D with URA3 present: " << A_to_start + D_to_stop + strlen(URA3_sequence) << " - ";
	 fout << A_to_start_sequence << URA3_sequence << D_to_stop_sequence << endl ;
	 
	 fout << "Finkl: " << Fwd_Cassette << endl;
	 fout << "Rinkl: " << Rev_Cassette << endl;
	 */
	sprintf(tempstr, "\t\t\t<CassettePresentLength>%d</CassettePresentLength>\n", A_to_start + D_to_stop + (int)strlen(cassette_sequence));
	strcat(cnf_results, tempstr);
	sprintf(tempstr, "\t\t\t<CassettePresentProduct>%s%s%s</CassettePresentProduct>\n", A_to_start_sequence, cassette_sequence, D_to_stop_sequence);
	strcat(cnf_results, tempstr);
	
	sprintf(tempstr, "\t\t\t<finkl>%s</finkl>\n\t\t\t<rinkl>%s</rinkl>\n", Fwd_Cassette, Rev_Cassette);
	strcat(cnf_results, tempstr);
	
	
	// Get Primer E to Rinkl sequence on RC_ura3
	char E_to_Rinkl_sequence[4000];
	
	x = 0;
	
	for(i = Primer_E_at; i < (plasmid_length - Rinkl_at); i++) E_to_Rinkl_sequence[x++] = RC_plasmid_sequence[i];
	E_to_Rinkl_sequence[x] = 0;
	
	sprintf(tempstr, "\t\t\t<cpsA_rinklLength>%d</cpsA_rinklLength>\n", A_to_start + 40 + E_to_Rinkl);
	strcat(cnf_results, tempstr);
	sprintf(tempstr, "\t\t\t<cpsA_rinklProduct>%s%s%s</cpsA_rinklProduct>\n", A_to_start_sequence, Primer_R, E_to_Rinkl_sequence);
	strcat(cnf_results, tempstr);
	
	// Get Primer D to Finkl sequence on ura3
	char D_to_Finkl_sequence[4000];
	
	x = 0;
	
	for(i = Primer_D_at; i < (plasmid_length - Finkl_at); i++) D_to_Finkl_sequence[x++] = plasmid_sequence[i];
	D_to_Finkl_sequence[x] = 0;
	/*	
	 fout << "CPS_D to Finkl: " << D_to_stop + 45 + D_to_Finkl << " - ";
	 */	
	char RC_D_to_stop_sequence[4000];
	reverse_complement(D_to_stop_sequence, RC_D_to_stop_sequence);
	/*	
	 fout << RC_D_to_stop_sequence << RC_Primer_C << D_to_Finkl_sequence << endl;
	 */	
	sprintf(tempstr, "\t\t\t<cpsD_finklLength>%d</cpsD_finklLength>\n", D_to_stop + 45 + D_to_Finkl);
	strcat(cnf_results, tempstr);
	sprintf(tempstr, "\t\t\t<cpsD_finklProduct>%s%s%s</cpsD_finklProduct>\n", RC_D_to_stop_sequence, RC_Primer_C, D_to_Finkl_sequence);
	strcat(cnf_results, tempstr);
	
	return(TRUE);	
}

#endif

#define BEST 0

int precise_process(const char* organism_name,
					const char* gene, 
				   const char* orf_sequence, 
				   const char* plasmid_sequence, 
				   const char* genome_file_name,
				   const char* plasmid_file_name,
					const char* Fwd_Cassette,
					const char* Rev_Cassette,
				   ofstream &fout)
{

	int i = 0;
	int j = 0;
	
#ifdef DISPLAY_RESULTS
	int examples_to_show;
#endif
	
	DNAfind genome_template(genome_file_name);
	
	if(!genome_template.set_max_mismatches(0)) 
		cout << "Could not set max_mismatches\n";
	
	if(!genome_template.set_tail_length(12)) 
		cout << "Could not set tail length\n";
	
	DNAfind plasmid_template(plasmid_file_name);
	
	if(!plasmid_template.set_max_mismatches(0)) 
		cout << "Could not set max_mismatches\n";
	
	if(!plasmid_template.set_tail_length(8)) 
		cout << "Could not set tail length\n";
	
// Check or find start codon
	
	int start = start_codon_location_for_prokaryotes(orf_sequence);
	
	if(start < 995) 
	{
		cout << "Start location not found \n";
		return(0);
	}
	
//  Make PCR1 primers
	
	primer_pair pcr1;
	
// Start by finding the best reverse primer
	pcr1.reverse_primer.set_primer_location_range(start - 1, start - 1); // Primer B starts immediately before start codon
	
#ifndef LACTIS  //default homopolymeric run check FALSE
	//pcr1.reverse_primer.homopolymeric_run_check = TRUE;
#endif

	pcr1.reverse_primer.set_primer_length_range(30, 60);
	pcr1.reverse_primer.optimum_primer_length = 40;
	
	pcr1.reverse_primer.generate_candidates(orf_sequence);
	
	pcr1.reverse_primer.analyse_all_candidates();	
	pcr1.reverse_primer.set_priorities("HAIRPIN, SELF_DIMER, LENGTH");
	pcr1.reverse_primer.sort_candidates();
	
	//pcr1.reverse_primer.show_all_single_candidates();
	

// Now set params for forward primer using reverse primer data	
	pcr1.forward_primer.set_primer_location_range(1, 500);
	pcr1.forward_primer.downstream_search = TRUE;	
	pcr1.forward_primer.set_primer_length_range(35, 45);
	pcr1.forward_primer.homopolymeric_run_check = TRUE; // We can afford to be choosy with the forward primer and reduce potential secondary products	
	
// To get a close Tm we constrain GC content in fwd primer to be close to the best rev primer
	//pcr1.forward_primer.required_GC_content = sequence_utils::GC_content(pcr1.reverse_primer.candidate[0].sequence);
	//pcr1.forward_primer.GC_tolerance = 2;
	
	pcr1.forward_primer.generate_candidates(orf_sequence);
	
	pcr1.forward_primer.analyse_all_candidates();
	
	//cout << "Unsorted forward primers\n";
	//pcr1.forward_primer.show_all_single_candidates();
	
	//cout << "Ideal temp = " << pcr1.reverse_primer.candidate[0].annealing_temperature << endl;
	
	pcr1.forward_primer.optimum_Tm = pcr1.reverse_primer.candidate[0].annealing_temperature;
	pcr1.forward_primer.max_Tm = pcr1.reverse_primer.candidate[0].annealing_temperature + 10;
	pcr1.forward_primer.min_Tm = pcr1.reverse_primer.candidate[0].annealing_temperature - 10;
	
	pcr1.forward_primer.set_priorities("HAIRPIN, SELF_DIMER, TEMPERATURE");
	pcr1.forward_primer.sort_candidates();

	//cout << "Sorted forward primers\n";
	//pcr1.forward_primer.show_all_single_candidates();
	
	//cout << "Fwd Candidates = " << pcr1.forward_primer.candidates_found << endl;
	//cout << "Rev Candidates = " << pcr1.reverse_primer.candidates_found << endl;

	cout << "Good pcr1 cands: " << pcr1.forward_primer.good_candidates << " and " << pcr1.reverse_primer.good_candidates << endl;

// Check that we have pcr1 primers
	if(pcr1.forward_primer.good_candidates < 1) return(NO_PRIMER_A);
	if(pcr1.reverse_primer.good_candidates < 1) return(NO_PRIMER_B);
	
	
// Make and test pairs	
	int fwd_candidates, rev_candidates;
	
	if(pcr1.forward_primer.good_candidates > 5) 
		fwd_candidates = 6;
	else
		fwd_candidates = pcr1.forward_primer.good_candidates;
	
	if(pcr1.reverse_primer.good_candidates > 5)
		rev_candidates = 6;
	else
		rev_candidates = pcr1.reverse_primer.good_candidates;
	
	//cout << "Actual cands " << fwd_candidates << ", " << rev_candidates << endl;
	
	pcr1.make_pair_candidates(fwd_candidates, rev_candidates);
	
	for(i = 0; i < pcr1.number_of_pair_candidates; i++)
	{
		//cout << pcr1.pair_candidate[i].forward_sequence << ", " << pcr1.pair_candidate[i].reverse_sequence << endl;
		pcr1.pair_candidate[i].number_of_pcr_products  =  genome_template.search_for_pcr_products(pcr1.pair_candidate[i].forward_sequence, pcr1.pair_candidate[i].reverse_sequence);
		pcr1.pair_candidate[i].number_of_pcr_products += plasmid_template.search_for_pcr_products(pcr1.pair_candidate[i].forward_sequence, pcr1.pair_candidate[i].reverse_sequence);
		
		//cout << "Pair " << i << ": " << pcr1.pair_candidate[i].number_of_pcr_products << endl;
	}	

	
	
	//temp_products = genome_template.search_for_pcr_products(pcr1.pair_candidate[0].forward_sequence, pcr1.pair_candidate[0].reverse_sequence);
	//cout << "Products = " << temp_products << endl;
	
	pcr1.sort_pair_candidates("TM_DIFF, F_DIMER, R_DIMER, PRODUCTS, MOO_SORT");
	
	if(pcr1.good_pair_candidates < 1)
	{
		return(NO_PCR1_PRIMERS);
	}
	
#ifdef DISPLAY_RESULTS
	 // Display best candidates
	 
	if(pcr1.good_pair_candidates > 5)
		examples_to_show = 6;
	else 
		examples_to_show = pcr1.good_pair_candidates;
	
	pcr1.show_best_pair_candidates(examples_to_show);
#endif
	
	//Make PCR1 product
	char pcr1_product[1001];
	j = 0;
	
	for(i = pcr1.pair_candidate[0].location_forward_5_prime_end; i < start; i++)
	{
		if(j >= 1000)
		{
			cout << "Error: pcr1 product over 1000 bases\n";
			return(0);
		}
		pcr1_product[j++] = orf_sequence[i];	
	}
	pcr1_product[j] = 0;

//  Make plasmid primers
		
	primer_pair plasmid;

#ifdef LACTIS
// For Lactis: pCS1966 0-4740, ermAM 2623-3360 & oroP 3786-4709
	plasmid.forward_primer.set_primer_location_range(2573, 2623);
	plasmid.reverse_primer.set_primer_location_range(4709, 4740);
#endif
	
#ifdef POMBE
// For Pombe: pFS118 sequence is 7282 bases long, promoter 4932-4981 URA4 4982-5776, terminator 5797-6211
	plasmid.forward_primer.set_primer_location_range(4832, 4932);
	plasmid.reverse_primer.set_primer_location_range(6211, 6311);
#endif
	
	plasmid.set_primer_length_range(18, 20);
	
	plasmid.generate_candidates(plasmid_sequence);
	plasmid.candidate_analysis();
	
	/*cout << "Plasmid results \n";
	plasmid.forward_primer.show_all_single_candidates();
	plasmid.reverse_primer.show_all_single_candidates();*/
	
	plasmid.sort_individual_candidates("HAIRPIN, SELF_DIMER"); // Temps will be low here - not worth sorting?
	
// Check that we have pcr1 primers
	if(plasmid.forward_primer.good_candidates < 1) return(NO_PRIMER_E);
	if(plasmid.reverse_primer.good_candidates < 1) return(NO_PRIMER_D);
	
// Make and test pairs	
	int plasmid_fwd_candidates, plasmid_rev_candidates;
	
	if(plasmid.forward_primer.good_candidates > 5) 
		plasmid_fwd_candidates = 6;
	else
		plasmid_fwd_candidates = plasmid.forward_primer.good_candidates;
	
	if(plasmid.reverse_primer.good_candidates > 5)
		plasmid_rev_candidates = 6;
	else
		plasmid_rev_candidates = plasmid.reverse_primer.good_candidates;
	
	plasmid.make_pair_candidates(plasmid_fwd_candidates, plasmid_rev_candidates);
	
	for(i = 0; i < plasmid.number_of_pair_candidates; i++)
	{
		//cout << pcr1.pair_candidate[i].forward_sequence << ", " << pcr1.pair_candidate[i].reverse_sequence << endl;
		plasmid.pair_candidate[i].number_of_pcr_products  =  genome_template.search_for_pcr_products(plasmid.pair_candidate[i].forward_sequence, plasmid.pair_candidate[i].reverse_sequence);
		plasmid.pair_candidate[i].number_of_pcr_products += plasmid_template.search_for_pcr_products(plasmid.pair_candidate[i].forward_sequence, plasmid.pair_candidate[i].reverse_sequence);
		
		//cout << "Pair " << i << ": " << pcr1.pair_candidate[i].number_of_pcr_products << endl;
	}
	
	plasmid.sort_pair_candidates("TM_DIFF, F_DIMER, R_DIMER, PRODUCTS, MOO_SORT");
	
	if(plasmid.good_pair_candidates < 1)
	{
		return(NO_PLASMID_PRIMERS);
	}
	
#ifdef DISPLAY_RESULTS
	// Display best candidate pairs
	if(plasmid.good_pair_candidates > 5)
		examples_to_show = 6;
	else 
		examples_to_show = plasmid.good_pair_candidates;
	
	plasmid.show_best_pair_candidates(examples_to_show);
#endif
	
//	Make C and R sequences
	
// 	C = 45 nt upstream of the stop codon and include stop
//  R = 40 nt downstream of the stop codon not including the stop codon
	char sequence_C[50];
	char sequence_R[50];
	char rc_sequence_C[50];
	int stop = stop_codon_location(start, orf_sequence);
	
	//cout << "Stop = " << stop << endl;
	
   	for(i = 0; i < 45; i++) sequence_C[i] = orf_sequence[i + stop - 42];
   	sequence_C[45] = 0;
	
   	for(i = 0; i < 40; i++) sequence_R[i] = orf_sequence[i + stop +3];
    sequence_R[40] = 0;
	
	//cout << "C = " << sequence_C << " and R = " << sequence_R << endl;
	

	
	// Get PCR1 product/PCR2 forward primer homology sequence
	char pcr1_homology[45];
#define PCR1_INIT_HOMOLOGY_LENGTH 25	
#define PCR1_HOMOLOGY_LENGTH 40
	
	for(i = 0; i < PCR1_INIT_HOMOLOGY_LENGTH; i++) 
		pcr1_homology[i] = pcr1_product[strlen(pcr1_product) - PCR1_INIT_HOMOLOGY_LENGTH + i];
	
	pcr1_homology[PCR1_INIT_HOMOLOGY_LENGTH] = 0; // terminate string
	
	// The homology sequence is only 25 bases unless the GC content < 40%, then make it 40 bases:
	
	int GC_number = 0;
	int pcr1_homology_length = strlen(pcr1_homology);
	
	for(i = 0; i < pcr1_homology_length; i++)
		if(pcr1_homology[i] == 67 || pcr1_homology[i] == 71) GC_number++;
	
	double GC_content = (double)GC_number/(double)pcr1_homology_length;
	
	if(GC_content < 0.4)
	{
		for(i = 0; i < PCR1_HOMOLOGY_LENGTH; i++) 
			pcr1_homology[i] = pcr1_product[strlen(pcr1_product) - PCR1_HOMOLOGY_LENGTH + i];
		pcr1_homology[PCR1_HOMOLOGY_LENGTH] = 0;
	}
	
	// Make Reverse complement Primer C
	reverse_complement(sequence_C, rc_sequence_C);
	
	//  Make PCR2 primers
	//cout << "PCR2 primers\n";
	
	//chimaera.forward_primer = B*RE and chimaera.reverse_primer = C'D. B* is actually 25/40 bases from the 3' end of the 
	// pcr1 product sense strand, which is the reverse comp of pcr1.reverse_primer (aka Primer B).
	
	primer_pair pcr2;

	for(i = 0; i < 6; i++)
	{
		// Make pcr2 forward   	
		strcpy(pcr2.forward_primer.candidate[i].sequence, pcr1_homology);
		strcat(pcr2.forward_primer.candidate[i].sequence, sequence_R);
		strcat(pcr2.forward_primer.candidate[i].sequence, plasmid.pair_candidate[i].forward_sequence);  
		
		// Make pcr2 reverse 	   	
		strcpy(pcr2.reverse_primer.candidate[i].sequence, rc_sequence_C);
		strcat(pcr2.reverse_primer.candidate[i].sequence, plasmid.pair_candidate[i].reverse_sequence);
		
	}
	
	// Need to tell primer class how many candidates have been manually set
	pcr2.forward_primer.candidates_found = 6;
	pcr2.reverse_primer.candidates_found = 6;
	
	pcr2.candidate_analysis();
	//pcr2.show_individual_candidates();
	pcr2.sort_individual_candidates("HAIRPIN, SELF_DIMER");

	cout << "Good pcr2 cands: " << pcr2.forward_primer.good_candidates << " and " << pcr2.reverse_primer.good_candidates << endl;

	// Check that we have pcr1 primers
	if(pcr2.forward_primer.good_candidates < 1) return(NO_PCR2_FORWARD);
	if(pcr2.reverse_primer.good_candidates < 1) return(NO_PCR2_REVERSE);
	
	// Make and test pairs	
	int pcr2_fwd_candidates, pcr2_rev_candidates;	
	
	if(pcr2.forward_primer.good_candidates > 5) 
		pcr2_fwd_candidates = 6;
	else
		pcr2_fwd_candidates = pcr2.forward_primer.good_candidates;
	
	if(pcr2.reverse_primer.good_candidates > 5)
		pcr2_rev_candidates = 6;
	else
		pcr2_rev_candidates = pcr2.reverse_primer.good_candidates;
	
	pcr2.make_pair_candidates(pcr2_fwd_candidates, pcr2_rev_candidates);
	
	for(i = 0; i < pcr2.number_of_pair_candidates; i++)
	{
		//cout << pcr2.pair_candidate[i].forward_sequence << ", " << pcr2.pair_candidate[i].reverse_sequence << endl;
		pcr2.pair_candidate[i].number_of_pcr_products  =  genome_template.search_for_pcr_products(pcr2.pair_candidate[i].forward_sequence, pcr2.pair_candidate[i].reverse_sequence);
		pcr2.pair_candidate[i].number_of_pcr_products += plasmid_template.search_for_pcr_products(pcr2.pair_candidate[i].forward_sequence, pcr2.pair_candidate[i].reverse_sequence);
		
		//cout << "Pair " << i << ": " << pcr2.pair_candidate[i].number_of_pcr_products << endl;
	}
	
	pcr2.sort_pair_candidates("TM_DIFF, F_DIMER, R_DIMER, PRODUCTS, MOO_SORT");
	
	if(pcr2.good_pair_candidates < 1)
	{
		return(NO_PCR2_PRIMERS);
	}
	
#ifdef DISPLAY_RESULTS
	 // Display best candidate pairs
	 
	if(pcr2.good_pair_candidates > 5)
		examples_to_show = 6;
	else 
		examples_to_show = pcr2.good_pair_candidates;
	
	pcr2.show_best_pair_candidates(examples_to_show);
#endif
	
	// Make cassette sequence
	// When we know which pcr2 candidates are to be used, only then can we determine the cassette
	// This requires specific chimaera handling code. For now we can use the plasmid components (D & E) of the pcr2
	// primers and find the plasmid sequence between them
	
	char primer_E[25];
	char primer_D[25];
	char rc_reverse_tail[25];
	int fwd_length = strlen(pcr2.pair_candidate[0].forward_sequence);
	int rev_length = strlen(pcr2.pair_candidate[0].reverse_sequence);
	int BR_length;
	
	// Length of the B*R component of the chimaera
	fwd_length > 90? BR_length = 80: BR_length = 65; 		
	
	for(i = 0; i < fwd_length - BR_length; i++) primer_E[i] = pcr2.pair_candidate[0].forward_sequence[BR_length + i];
		primer_E[fwd_length - BR_length] = 0;
	
	for(i = 0; i < rev_length - 45; i++) primer_D[i] = pcr2.pair_candidate[0].reverse_sequence[45 + i];
	primer_D[rev_length - 45] = 0;
	
	cout << "Primer E: " << primer_E << ", Primer D: " << primer_D << endl;
	
	int fwd_loc = findstring(primer_E, plasmid_sequence);
	
	reverse_complement(primer_D, rc_reverse_tail);
	
	int rev_loc = findstring(rc_reverse_tail, plasmid_sequence) + strlen(primer_D);
	
	cout << "RC primer D: " << rc_reverse_tail << endl;
	
	//cout << "Fwd loc = " << fwd_loc << " and " << rev_loc << endl;
		
	char cassette_sequence[5000];	
	
	j = 0;
	
	for(i = fwd_loc; i <= rev_loc; i++) cassette_sequence[j++] = plasmid_sequence[i];
	cassette_sequence[j] = 0;
	
	cout << cassette_sequence << endl;
	
	// Confirmation primers
	
	char cnf_results[30000] = "None";
	char csv_results[2000] = "None";
	
	/*if(XML)
	{
		if(!xml_confirmation_primers(gene, primer_d[idxP_D].sequence, primer_e[idxP_E].sequence, Primer_R, RC_Primer_C, cnf_results))
			return(CONF_PRIMER_ERROR);
	}
	else 
	{
		if(!confirmation_primers(gene, primer_d[idxP_D].sequence, primer_e[idxP_E].sequence, Primer_R, RC_Primer_C, cnf_results))
			return(CONF_PRIMER_ERROR);
	}*/
	
	if(!precise_confirmation_primers(orf_sequence, genome_file_name, plasmid_sequence, primer_D, primer_E, rc_sequence_C, sequence_R, Fwd_Cassette, Rev_Cassette, cnf_results, csv_results))
		return(CONF_PRIMER_ERROR);
	
	// RESULTS
	// **************************************
	time_t mtime;
	time(&mtime);
	char mytime[64];
	
	strcpy(mytime, ctime(&mtime));	
	mytime[strlen(mytime) - 1] = 0;	
	
#define PCR1_BEST 0
#define PCR2_BEST 0
	
	char output_file_name[32];
	sprintf(output_file_name, "%s.xml", gene);
	fout.open(output_file_name);

		
		//#ifndef BEOWULF		
		fout << "<organism>\n\t<name>" << organism_name << "</name>\n";
		//#endif
		fout << "\t<geneForDeletion>\n\t\t<name>" << gene << "</name>\n\t\t<date>" << mytime << "</date>\n";
		//fout << "\t\t<class>" << PCR_RESULTS_CLASS << "</class>\n";
		
		// PCR1	
		fout << "\t\t<primerPair>\n";
		fout << "\t\t\t<name>PCR1</name>\n";
		
		fout << "\t\t\t<forwardPrimer>\n";
		fout << "\t\t\t\t<sequence>"       << pcr1.pair_candidate[PCR1_BEST].forward_sequence << "</sequence>\n";
		fout << "\t\t\t\t<hairpinScore>"   << pcr1.pair_candidate[PCR1_BEST].forward_hairpin_score << "</hairpinScore>\n";
		fout << "\t\t\t\t<selfDimerScore>" << pcr1.pair_candidate[PCR1_BEST].forward_self_dimer_score << "</selfDimerScore>\n";
		/*fout << "\t\t\t\t<bindingSites>\n";
		fout << "\t\t\t\t\t<yeast>" << primer_a.candidate[idxP_A].binding_A << "</yeast>\n";
		fout << "\t\t\t\t\t<ura3>" << primer_a.candidate[idxP_A].binding_B << "</ura3>\n";
		fout << "\t\t\t\t</bindingSites>\n";*/
		fout << "\t\t\t\t<annealingTemperature>" << pcr1.pair_candidate[PCR1_BEST].forward_annealing_temperature << "</annealingTemperature>\n";
		fout << "\t\t\t</forwardPrimer>\n";
		
		fout << "\t\t\t<reversePrimer>\n";
		fout << "\t\t\t\t<sequence>"       << pcr1.pair_candidate[PCR1_BEST].reverse_sequence << "</sequence>\n";
		fout << "\t\t\t\t<hairpinScore>"   << pcr1.pair_candidate[PCR1_BEST].reverse_hairpin_score << "</hairpinScore>\n";
		fout << "\t\t\t\t<selfDimerScore>" << pcr1.pair_candidate[PCR1_BEST].reverse_self_dimer_score << "</selfDimerScore>\n";
		/*fout << "\t\t\t\t<bindingSites>\n";
		fout << "\t\t\t\t\t<yeast>" << primer_b.candidate[idxP_B].binding_A << "</yeast>\n";
		fout << "\t\t\t\t\t<ura3>" << primer_b.candidate[idxP_B].binding_B << "</ura3>\n";
		fout << "\t\t\t\t</bindingSites>\n";*/
		fout << "\t\t\t\t<annealingTemperature>" << pcr1.pair_candidate[PCR1_BEST].reverse_annealing_temperature << "</annealingTemperature>\n";
		fout << "\t\t\t</reversePrimer>\n";
		
		fout << "\t\t\t<primerDimerForward>" << pcr1.pair_candidate[PCR1_BEST].forward_pair_dimer_score << "</primerDimerForward>\n";
		fout << "\t\t\t<primerDimerReverse>" << pcr1.pair_candidate[PCR1_BEST].reverse_pair_dimer_score << "</primerDimerReverse>\n";
		
		fout << "\t\t\t<numberOfProducts>"   << pcr1.pair_candidate[PCR1_BEST].number_of_pcr_products << "</numberOfProducts>\n";
		
		fout << "\t\t\t<product>"       << pcr1_product << "</product>\n";
		fout << "\t\t\t<productLength>" << strlen(pcr1_product) << "</productLength>\n";
		
		fout << "\t\t</primerPair>\n";
		
		// PCR2
		fout << "\t\t<primerPair>\n";
		fout << "\t\t\t<name>PCR2</name>\n";
		
		fout << "\t\t\t<forwardPrimer>\n";
		fout << "\t\t\t\t<sequence>"       << pcr2.pair_candidate[PCR2_BEST].forward_sequence << "</sequence>\n";
		fout << "\t\t\t\t<hairpinScore>"   << pcr2.pair_candidate[PCR2_BEST].forward_hairpin_score << "</hairpinScore>\n";
		fout << "\t\t\t\t<selfDimerScore>" << pcr2.pair_candidate[PCR2_BEST].forward_self_dimer_score << "</selfDimerScore>\n";
		/*fout << "\t\t\t\t<bindingSites>\n";  
		fout << "\t\t\t\t\t<yeast>" << PCR2_forward_data_p->binding_A << "</yeast>\n";
		fout << "\t\t\t\t\t<ura3>" << PCR2_forward_data_p->binding_B << "</ura3>\n";		
		fout << "\t\t\t\t</bindingSites>\n";*/
		fout << "\t\t\t\t<annealingTemperature>" << pcr2.pair_candidate[PCR2_BEST].forward_annealing_temperature << "</annealingTemperature>\n";
		
		fout << "\t\t\t</forwardPrimer>\n";
		
		fout << "\t\t\t<reversePrimer>\n";
		fout << "\t\t\t\t<sequence>"       << pcr2.pair_candidate[PCR2_BEST].reverse_sequence << "</sequence>\n";
		fout << "\t\t\t\t<hairpinScore>"   << pcr2.pair_candidate[PCR2_BEST].reverse_hairpin_score << "</hairpinScore>\n";
		fout << "\t\t\t\t<selfDimerScore>" << pcr2.pair_candidate[PCR2_BEST].reverse_self_dimer_score << "</selfDimerScore>\n";
		/*fout << "\t\t\t\t<bindingSites>\n";  
		fout << "\t\t\t\t\t<yeast>" << PCR2_reverse_data_p->binding_A << "</yeast>\n";
		fout << "\t\t\t\t\t<ura3>" << PCR2_reverse_data_p->binding_B << "</ura3>\n";
		fout << "\t\t\t\t</bindingSites>\n";*/
		fout << "\t\t\t\t<annealingTemperature>" << pcr2.pair_candidate[PCR2_BEST].reverse_annealing_temperature << "</annealingTemperature>\n";
		
		fout << "\t\t\t</reversePrimer>\n";
		
		fout << "\t\t\t<primerDimerForward>" << pcr2.pair_candidate[PCR2_BEST].forward_pair_dimer_score << "</primerDimerForward>\n";
		fout << "\t\t\t<primerDimerReverse>" << pcr2.pair_candidate[PCR2_BEST].reverse_pair_dimer_score << "</primerDimerReverse>\n";
		
		fout << "\t\t\t<numberOfProducts>" << pcr2.pair_candidate[PCR2_BEST].number_of_pcr_products << "</numberOfProducts>\n";
		fout << "\t\t\t<primer_B>" << pcr1_homology << "</primer_B>\n";
		fout << "\t\t\t<primer_R>" << sequence_R  << "</primer_R>\n";
		fout << "\t\t\t<primer_E>" << primer_E << "</primer_E>\n";
		fout << "\t\t\t<primer_C>" << sequence_C  << "</primer_C>\n";
		fout << "\t\t\t<primer_D>" << primer_D << "</primer_D>\n";
		fout << "\t\t\t<product>" << pcr1_homology << sequence_R << cassette_sequence << sequence_C << "</product>\n";
		fout << "\t\t\t<productLength>" << strlen(pcr1_homology) + strlen(sequence_R) + strlen(cassette_sequence) + strlen(sequence_C) << "</productLength>\n";
		
		
		fout << "\t\t</primerPair>\n";
		// PCR3
		fout << "\t\t<pcr3Product>" << pcr1_product << sequence_R << cassette_sequence << sequence_C << "</pcr3Product>\n";
		fout << "\t\t<pcr3ProductLength>" << strlen(pcr1_product) + strlen(sequence_R) + strlen(cassette_sequence) + strlen(sequence_C) << "</pcr3ProductLength>\n"; 
		
		// CONFIRMATION RESULTS
		fout << "\t\t<confirmationData>\n" << cnf_results << "\t\t</confirmationData>\n";
		
		fout << "\t</geneForDeletion>\n";
		//#ifndef BEOWULF
		fout << "</organism>\n";
		//#endif

	fout.close();
	
	/*
	 To collate all .csv results use cat *.csv > all.csv
	 */

#ifdef DB_OUT
	char csv_file_name[32];
	sprintf(csv_file_name, "%s.csv", gene);	
	ofstream fcsv(csv_file_name);
	fcsv << gene << ", ";
	fcsv << pcr1.pair_candidate[PCR1_BEST].forward_hairpin_score << ", ";
	fcsv << pcr1.pair_candidate[PCR1_BEST].forward_self_dimer_score << ", ";
	fcsv << pcr1.pair_candidate[PCR1_BEST].reverse_hairpin_score << ", ";
	fcsv << pcr1.pair_candidate[PCR1_BEST].reverse_self_dimer_score << ", ";
	fcsv << pcr1.pair_candidate[PCR1_BEST].forward_pair_dimer_score << ", ";
	fcsv << pcr1.pair_candidate[PCR1_BEST].reverse_pair_dimer_score << ", ";
	fcsv << pcr1.pair_candidate[PCR1_BEST].number_of_pcr_products << ", ";
	fcsv << " - , ";
	fcsv << pcr2.pair_candidate[PCR2_BEST].forward_hairpin_score << ", ";
	fcsv << pcr2.pair_candidate[PCR2_BEST].forward_self_dimer_score << ", ";
	fcsv << pcr2.pair_candidate[PCR2_BEST].reverse_hairpin_score << ", ";
	fcsv << pcr2.pair_candidate[PCR2_BEST].reverse_self_dimer_score << ", ";
	fcsv << pcr2.pair_candidate[PCR2_BEST].forward_pair_dimer_score << ", ";
	fcsv << pcr2.pair_candidate[PCR2_BEST].reverse_pair_dimer_score << ", ";	
	fcsv << pcr2.pair_candidate[PCR2_BEST].number_of_pcr_products << ", ";
	fcsv << " - , ";
	fcsv << strlen(pcr1_homology) << ", ";
	fcsv << strlen(sequence_R) << ", ";	
	fcsv << strlen(sequence_C) << ", ";
	fcsv << strlen(primer_E) << ", ";
	fcsv << strlen(primer_D) << ", ";
	fcsv << " - , ";
	fcsv << csv_results;
	fcsv << endl;
	
	fcsv.close();
#endif	
	
	return(TRUE);
}

/**
 Note that the Fwd_Cassette_confirmation_primer and Rev_Cassette_confirmation_primer are obtained by building 
 (use ./build_cassette_confirmation) and running cassette_confirmation.cpp
 */

int main(int argc, char** argv)
{
	char genebuffer[32];
	char buffer[4001];
	char orf_sequence[50000];
	char plasmid_sequence[10000];
	char query_gene_id[32];
	char header_id[32];
	bool found = FALSE;
	int error = 1;
	char *token;
	
#ifdef POMBE
	// Pombe parameters
	char organism_name[128] = "Schizosaccharomyces pombe";
	char genes_1000_all[128] = "pombe_1000_all.fa";
	char genome_file_name[128] = "pombe_genome_09052011.fasta";
	char plasmid_file_name[128] = "pFS118.fa";
	char Fwd_Cassette_confirmation_primer[40] = "TACCGTCAAGCTACAATATGCATCTG";
	char Rev_Cassette_confirmation_primer[40] = "CTGCGAATTTGCGATCCTCAAAG";
#endif

#ifdef S_CERE
	// S. cere
	char organism_name[128] = "Saccharomyces cereviseae";
	char genes_1000_all[128] = "orf_genomic_1000_all.fasta";
	char genome_file_name[128] = "s_cere_genome.fa";
	char plasmid_file_name[128] = "rc_ura3.fa";
	char Fwd_Cassette_confirmation_primer[40] = "ACGTTCGTTCGACTGATGAG";
	char Rev_Cassette_confirmation_primer[40] = "CTCATCAGTCGAACGAACGT";
#endif
	
#ifdef LACTIS
	// Lactis
	char organism_name[128] = "Lactococcus lactis";
	char genes_1000_all[128] = "lactis_220612_1000_all.fa";
	char genome_file_name[128] = "l_lactis.fa";
	char plasmid_file_name[128] = "pCS1966.fa";
	char Fwd_Cassette_confirmation_primer[40] = "TTCATCCTAAACCAAAAGTAAACAGTGT";
	char Rev_Cassette_confirmation_primer[40] = "AAGGTACGCTTGTAGAATCCTTCTTC";
#endif

   	cout << "Precise_Deletion: Primer design\n";
	cout << "Version 1.3, 1st August 2012\n";
	
	if(!(argv[1] && argv[2]))
	{
		cout << "Use: exe <list of genes file> <missing genes output file>\n";
		return(0);
	}
	
	ofstream fout;
	
	// Get plasmid sequence
	
	ifstream plasmid(plasmid_file_name);
	if(!plasmid.is_open()) cout << "file opening error\n";
	
	plasmid_sequence[0] = 0; // init for strcat
	
	while(plasmid.getline(buffer, 4000))
	{
		if(buffer[0] != 0x3E)
		{
			token = strtok(buffer, " \t\n\r");
			if(token != NULL)strcat(plasmid_sequence, token);
		}
	}	
	plasmid.close();
	
	//cout << plasmid_sequence << endl << strlen(plasmid_sequence) << endl;
	
	// Open file with the list of genes to be primer'ed
	
	ifstream genelist(argv[1]);
	if(!genelist.is_open())
	{
		cout << "Unable to open " << argv[1] << endl;
		return(0);
	}
	
	// Open missing results file
	
	ofstream ferr(argv[2]);
	
	// For each gene get the sequence 1000 upstream, gene and 1000 downstream and name it orf
	
	ifstream fin;

	while(genelist.getline(genebuffer, 32))
	{
		//if(genebuffer){
		
		cout << genebuffer << endl;
		
		token = strtok(genebuffer,"> \t\n");
		
		if(token)
		{
			strcpy(query_gene_id, token);
		
			if(!(query_gene_id[0] == 37)) // 37 = % => comment
			{
				fin.open(genes_1000_all);
				
				if(!fin.is_open())
				{
					cout << "unable to open " << genes_1000_all << endl;
					break;
				}
				
				found = FALSE;
				orf_sequence[0] = 0; 
				
				while(fin.getline(buffer, 2055))
				{
					if(buffer[0] == 0x3E)  // fasta header
					{
						if(found)break;
						strcpy(header_id, strtok(buffer,"> \t\n"));
						if(!strcmp(header_id, query_gene_id))found = TRUE;			
					}
					else if(found) strcat(orf_sequence, buffer);
					
					//if(i++ > 1000)break;
				}  
				fin.close();
				
				if(!found)
				{
					cout << "Could not find " << query_gene_id << endl;
					//initialization_error = TRUE;
					//return(FALSE);
				}
				else
				{
					//cout << "Found " << header_id << endl;
					
					//error = auto_process(header_id, orf_sequence, genome_file_name, fout);
					error = precise_process(organism_name, header_id, orf_sequence, plasmid_sequence, genome_file_name, plasmid_file_name, Fwd_Cassette_confirmation_primer, Rev_Cassette_confirmation_primer, fout);
					
					switch(error)
					{
						case NO_PRIMER_B: ferr << header_id << ": No primer B\n"; break;
						case NO_PRIMER_A: ferr << header_id << ": No primer A\n"; break;
						case NO_PCR1_PRIMERS: ferr << header_id << "; No PCR1 primer pairs after sorting\n"; break;
						case NO_PCR2_FORWARD: ferr << header_id << ": No PCR2 forward primers \n"; break;
						case NO_PCR2_REVERSE: ferr << header_id << ": No PCR2 reverse primers \n"; break;
						case NO_PCR2_PRIMERS: ferr << header_id << ": No PCR2 primer pair after sorting\n"; break;
						case NO_PRIMER_D: ferr << header_id << ": No plasmid primer D \n"; break;
						case NO_PRIMER_E: ferr << header_id << ": No plasmid primer E \n"; break;
						case NO_PLASMID_PRIMERS: ferr << header_id << ": No plasmid primer pairs after sorting\n"; break;
						case NO_URA3: ferr << header_id << ": Unable to open ura3_rc.fa \n"; break;
						case CONF_PRIMER_ERROR : ferr << header_id << ": Confirmation primers error\n"; break;
					}
					
					
				}
			}
		} // end if(token)
	}
	
	fin.close();
	ferr.close();
	
    return(1);
}
