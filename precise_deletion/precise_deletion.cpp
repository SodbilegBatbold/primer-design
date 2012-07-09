/*************************************************\
 
 precise_auto.cpp (Ver 1.00 - 10/2/11)
 
 Updated from:
	s_primerII.cpp (Ver. 1.00 - 5/7/10) 
	sprimer.cpp (Ver. 1.00 - 7/6/10) to run on Beowulf
 
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

using namespace std;

//#include "primrose_utils.h"
#include "../pd5/primer_pair.h"
#include "../pd5/dna_find.h"
#include "../pd5/annealing_temperature.h"
#include "../pd5/sequence_utils.h"

// OUTPUT FORMAT
#define XML 1

// Secondary binding search 
// BLAST, FASTA or NSB
#define NSB

// CLASSIFICATION OF RESULTS
#define PCR_RESULTS_CLASS "A"


// Error defines

#define NO_PRIMER_B 0
#define NO_PRIMERS_A_AND_B -1
#define NO_PCR2_PRIMERS -2
#define NO_PRIMER_D -3
#define NO_PRIMER_E -4
#define NO_URA3 -5
#define CONF_PRIMER_ERROR -6

// TESTING
//#define TESTING
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

int findstring(const char* a_string, const char* b_string)
{
	unsigned int i, count;
	int location[30];
	bool found = FALSE;
	
	count = 0;
	
	for(i = 0; i < (strlen(b_string) - strlen(a_string)); i++)
	{
		found = TRUE;
		for(unsigned int j = 0; j < strlen(a_string); j++)
		{
			if(a_string[j] != b_string[i + j])
			{
				found = FALSE; 
				break;
			}
		}
		
		if(found)
		{
			location[count++] = i;
			//cout << "Found " << a_string << " at " << location << endl;
			if(count >= 30)break;
		}
	}
	
	if(count == 0)
	{
		cout << "Could not find " << a_string << endl;
		return(-2);
	}
	else if(count == 1)
		return(location[0]);
	else 
	{
		cout << "Multiple hits for " << a_string << endl;
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

int stop_codon_location(int start_codon_location, const char* orf)
{
	int stop_codon_location = 0;
	
	// Find the stop codon
	
	for(unsigned int i = start_codon_location; i < strlen(orf); i += 3)
	{   if(orf[i] == 84)
		{
			if(orf[i+1] == 65 && orf[i+2] == 71){stop_codon_location = i; break;}
			else if(orf[i+1] == 65 && orf[i+2] == 65){stop_codon_location = i; break;}
			else if(orf[i+1] == 71 && orf[i+2] == 65){stop_codon_location = i; break;}
		}
   	}	
	
	return(stop_codon_location);
	
}

/***********************************************************************\
 
 \***********************************************************************/

//int confirmation_primers(const char* query_gene, const char* Primer_D, const char* Primer_E, const char* Primer_R, const char* RC_Primer_C, ofstream &fout)
int confirmation_primers(const char* query_gene, const char* Primer_D, const char* Primer_E, const char* Primer_R, const char* RC_Primer_C, char* cnf_results)
{
	char buffer[2056];
	char tbuffer[2056];
   	int i = 0;
	char orf[20000];
	char ura3[20000];
	char RC_orf[20000];
	//char query_gene[32];
	char header_id[32];
	//char output[64];
	bool found;
	
	char CPS_A[40];
	char CPS_B[40];
	char CPS_C[40];
	char CPS_D[40];

	char temp_A[40];
	char temp_B[40];
	char temp_C[40];
	char temp_D[40];
	
	char gene_present_confirmation_sequence[20000];
	char gene_deleted_confirmation_sequence[2000];
	
	char tempstr[15000]; // for message output cnf_results
	
	cout << "Confirmation primers 26/07/10\n";
	
	ifstream fin("orf_genomic_1000_all.fasta");
	if(!fin.is_open())
	{
		cout << "unable to open orf_genomic_1000_all.fasta \n";
		return(0);
	}
	
	found = FALSE;
	orf[0] = 0; 
	
	while(fin.getline(buffer, 2055))
	{
		if(buffer[0] == 0x3E)  // fasta header
		{
			if(found)break;
			strcpy(header_id, strtok(buffer,"> \t\n"));
			if(!strcmp(header_id, query_gene))found = TRUE;			
		}
		else if(found) strcat(orf, buffer);
		
		//if(i++ > 1000)break;
	}  
	fin.close();
	
	if(!found)
	{
		cout <<  query_gene << "Not in orf_genomic_100_all.fasta" << endl;
		//initialization_error = TRUE;
		return(FALSE);
	}
	else
	{
		cout << "Found " << header_id << endl;
	}
	
	/*	
	 YDR093W  	              	         	
	 ATAGACAAACCAGATTCATCGCTTA     	
	 GTGTGCTAATATCTTTCATTGGCTT     	
	 AGATCGTTAGAGAAATGTGGTTACG     	
	 CCGTAAAGTAGTTACTTGCCACAAT
	 */	
	ifstream fin2("Eurofan_CPs.txt");
	if(!fin2.is_open())
	{
		//cout << "unable to open confirmation_primers_S_cere.txt \n";
		cout << "unable to open Eurofan_CPs.txt \n";
		return(0);
	}
	
	found = FALSE;
	
	while(fin2.getline(buffer, 2055))
	{
		strcpy(tbuffer, buffer);
		
		strtok(buffer,"\t\n\r");
		strcpy(temp_A, strtok(NULL,"\t\n\r"));
		strcpy(header_id, strtok(temp_A," "));
		
		//cout << "*" << header_id << "*" << endl;
		
		if(!strcmp(header_id, query_gene))
		{
			found = TRUE;
			//cout << "FOund: " << header_id << endl;
			strtok(tbuffer,"\t\n\r");
			strtok(NULL,"\t\n\r");
			//strtok(NULL,"\t\n\r");
			strcpy(temp_A, strtok(NULL,"\t\n\r"));
			strcpy(temp_B, strtok(NULL,"\t\n\r"));
			strcpy(temp_C, strtok(NULL,"\t\n\r"));
			strcpy(temp_D, strtok(NULL,"\t\n\r"));
			
			/*for(i = 0; i < 30; i++) temp[i] = tbuffer[i + 49];
			 strcpy(CPS_A, strtok(temp," \t\n\r"));
			 
			 for(i = 0; i < 30; i++) temp[i] = tbuffer[i + 145];
			 strcpy(CPS_D, strtok(temp," \t\n\r"));*/
			break;
		}
	} 
	
	if(!found) return(-1);
	
	strcpy(CPS_A, strtok(temp_A," "));
	strcpy(CPS_B, strtok(temp_B," "));
	strcpy(CPS_C, strtok(temp_C," "));
	strcpy(CPS_D, strtok(temp_D," "));
/*
	 fout << "CPS A (forward): " << CPS_A << endl;
	 //cout << CPS_B << ", ";
	 //cout << CPS_C << ", ";
	 fout << "CPS D (reverse): " << CPS_D << endl;
*/
	strcpy(cnf_results,"CPS A (forward): ");
	strcat(cnf_results, CPS_A);
	strcat(cnf_results, "\nCPS D (reverse): ");
	strcat(cnf_results, CPS_D);
	
	cout << "Length of ORF: " << strlen(orf) << endl;
	
	int forward_at = findstring(CPS_A, orf); //cout << "Forward strand\n";
	reverse_complement(orf, RC_orf);
	int reverse_at = findstring(CPS_D, RC_orf); //cout << "Reverse strand\n";
	
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
	int gene_length = strlen(orf) - 2000;
	
// Make sequences
	int x = 0;
	int y = 0;
	int orf_length = strlen(orf);
	
	for(i = forward_at; i < orf_length - reverse_at; i++)
	{
		gene_present_confirmation_sequence[x++] = orf[i];
		
		if(i < 1000 || i >= orf_length - 1000)
			gene_deleted_confirmation_sequence[y++] = orf[i];
	}
	gene_present_confirmation_sequence[x] = 0;
	gene_deleted_confirmation_sequence[y] = 0;
		
/*	
	fout << "CPS A to CPS D with gene present: " << A_to_start + D_to_stop + gene_length << " - " << gene_present_confirmation_sequence << endl;
	fout << "CPS A to CPS D with gene deleted: " << A_to_start + D_to_stop << " - " << gene_deleted_confirmation_sequence << endl ;
 */	
	
	
	sprintf(tempstr, "\nCPS A to CPS D with gene present: %d - %s\n", A_to_start + D_to_stop + gene_length, gene_present_confirmation_sequence);
	strcat(cnf_results, tempstr);
	sprintf(tempstr, "CPS A to CPS D with gene deleted: %d - %s\n", A_to_start + D_to_stop, gene_deleted_confirmation_sequence);
	strcat(cnf_results, tempstr);
	
	/***********************************************************************\
	 
	 Now need to find these:
	 
	 FinklURA3 5'-ACGTTCGTTCGACTGATGAG-3'
	 RinklURA3 5'-CTCATCAGTCGAACGAACGT-3'
	 
	 in the URA3 sequence between Primer D and Primer E
	 
	 PRIMER_D "GATCCCAATACAACAGAT" 
	 PRIMER_E "CACACATTACTTGCCTC"
	 
	 \***********************************************************************/
	
	// Get URA3
	
	fin.open("ura3_rc.fa");
	if(!fin.is_open())
	{
		cout << "unable to open ura3_rc.fa \n";
		return(0);
	}
	
	// Get Sequence
	
	x = 0; 
	
	while(fin.getline(buffer, 2055))
	{		
		if(buffer[0] != 0x3E)  // fasta header
		{
			int buffer_length = strlen(buffer);
			
			for(i = 0; i < buffer_length; i++)
			{
				if(buffer[i] > 31) // remove ctrl chars
					ura3[x++] = buffer[i];
			}
		}
	} ura3[x] = 0; // terminate string
	fin.close();
	
	// Find the cassette between E and D
	char RC_ura3[20000];
	
	int Primer_D_at = findstring(Primer_D, ura3); //cout << "Primer Forward strand\n";
	reverse_complement(ura3, RC_ura3);
	int Primer_E_at = findstring(Primer_E, RC_ura3); //cout << "Primer Reverse strand\n";
	
	//cout << endl << ura3 << " Length is " << strlen(ura3) << endl;
	
	char FinklURA3[40] = "ACGTTCGTTCGACTGATGAG";
	char RinklURA3[40] = "CTCATCAGTCGAACGAACGT";
	
	//cout << "Length of URA3: " << strlen(ura3) << endl;
	int Rinkl_at = findstring(RinklURA3, ura3); //cout << "Forward strand\n";
	reverse_complement(ura3, RC_ura3);
	int Finkl_at = findstring(FinklURA3, RC_ura3); //cout << "Reverse strand\n";
	
	//cout << "RinklURA3: " << forward_at << endl ;
	//cout << "FinklURA3: " << reverse_at << endl;
	
	int ura3_length = strlen(ura3);
	
	int E_to_Rinkl = ura3_length - Rinkl_at - Primer_E_at;
	int D_to_Finkl = ura3_length - Finkl_at - Primer_D_at;
	
	//cout << "E to Rinkl: " << E_to_Rinkl << endl;
	//cout << "D to Finkl: " << D_to_Finkl << endl;
	
// Get URA3 sequence
	char URA3_sequence[20000];
	
	x = 0;
	
	for(i = Primer_E_at; i < (ura3_length - Primer_D_at); i++) URA3_sequence[x++] = RC_ura3[i];
	URA3_sequence[x] = 0;

// Get sequence from CPS A to the start of the gene
	char A_to_start_sequence[1000];
	
	x = 0;
	
	for(i = forward_at; i < 1000; i++) A_to_start_sequence[x++] = orf[i];
	A_to_start_sequence[x] = 0;

// Get the sequence from end of gene to CPS D
	char D_to_stop_sequence[1000];
	
	x = 0;
	
	for(i = orf_length - 1000; i < orf_length - reverse_at; i++) D_to_stop_sequence[x++] = orf[i];
	D_to_stop_sequence[x] = 0;
	
/*	
	fout << "CPS A to CPS D with URA3 present: " << A_to_start + D_to_stop + strlen(URA3_sequence) << " - ";
	fout << A_to_start_sequence << URA3_sequence << D_to_stop_sequence << endl ;
	
	fout << "Finkl: " << FinklURA3 << endl;
	fout << "Rinkl: " << RinklURA3 << endl;
*/
	sprintf(tempstr, "CPS A to CPS D with URA3 present: %d - %s%s%s\n", A_to_start + D_to_stop + (int)strlen(URA3_sequence), A_to_start_sequence, URA3_sequence, D_to_stop_sequence);
	strcat(cnf_results, tempstr);
	sprintf(tempstr, "Finkl: %s\nRinkl: %s\n", FinklURA3, RinklURA3);
	strcat(cnf_results, tempstr);
	
// Get Primer E to Rinkl sequence on RC_ura3
	char E_to_Rinkl_sequence[4000];
	
	x = 0;
	
	for(i = Primer_E_at; i < (ura3_length - Rinkl_at); i++) E_to_Rinkl_sequence[x++] = RC_ura3[i];
	E_to_Rinkl_sequence[x] = 0;

/*
	fout << "CPS_A to Rinkl: " << A_to_start + 40 + E_to_Rinkl << " - ";
	fout << A_to_start_sequence << Primer_R << E_to_Rinkl_sequence << endl;
*/
	sprintf(tempstr, "CPS_A to Rinkl: %d - %s%s%s\n", A_to_start + 40 + E_to_Rinkl, A_to_start_sequence, Primer_R, E_to_Rinkl_sequence);
	strcat(cnf_results, tempstr);
	
	
// Get Primer D to Finkl sequence on ura3
	char D_to_Finkl_sequence[4000];
	
	x = 0;
	
	for(i = Primer_D_at; i < (ura3_length - Finkl_at); i++) D_to_Finkl_sequence[x++] = ura3[i];
	D_to_Finkl_sequence[x] = 0;
/*	
	fout << "CPS_D to Finkl: " << D_to_stop + 45 + D_to_Finkl << " - ";
*/	
	char RC_D_to_stop_sequence[4000];
	reverse_complement(D_to_stop_sequence, RC_D_to_stop_sequence);
/*	
	fout << RC_D_to_stop_sequence << RC_Primer_C << D_to_Finkl_sequence << endl;
*/	
	sprintf(tempstr, "CPS_D to Finkl: %d - %s%s%s\n", D_to_stop + 45 + D_to_Finkl, RC_D_to_stop_sequence, RC_Primer_C, D_to_Finkl_sequence);
	strcat(cnf_results, tempstr);
	
	return(1);
}

/***********************************************************************\
 
 \***********************************************************************/

//int confirmation_primers(const char* query_gene, const char* Primer_D, const char* Primer_E, const char* Primer_R, const char* RC_Primer_C, ofstream &fout)
int xml_confirmation_primers(const char* query_gene, const char* Primer_D, const char* Primer_E, const char* Primer_R, const char* RC_Primer_C, char* cnf_results)
{
	char buffer[2056];
	char tbuffer[2056];
   	int i = 0;
	char orf[20000];
	char ura3[20000];
	char RC_orf[20000];
	//char query_gene[32];
	char header_id[32];
	//char output[64];
	bool found;
	
	char CPS_A[40];
	char CPS_B[40];
	char CPS_C[40];
	char CPS_D[40];

	char temp_A[40];
	char temp_B[40];
	char temp_C[40];
	char temp_D[40];
	
	char gene_present_confirmation_sequence[20000];
	char gene_deleted_confirmation_sequence[2000];
	
	char tempstr[15000]; // for message output cnf_results
	
	cout << "Confirmation primers 26/07/10\n";
	
	ifstream fin("orf_genomic_1000_all.fasta");
	if(!fin.is_open())
	{
		cout << "unable to open 1000_all.fa \n";
		return(0);
	}
	
	found = FALSE;
	orf[0] = 0; 
	
	while(fin.getline(buffer, 2055))
	{
		if(buffer[0] == 0x3E)  // fasta header
		{
			if(found)break;
			strcpy(header_id, strtok(buffer,"> \t\n"));
			if(!strcmp(header_id, query_gene))found = TRUE;			
		}
		else if(found) strcat(orf, buffer);
		
		//if(i++ > 1000)break;
	}  
	fin.close();
	
	if(!found)
	{
		cout <<  query_gene << ": Not found in orf_genomic_1000_all.fasta" << endl;
		//initialization_error = TRUE;
		return(FALSE);
	}
	else
	{
		cout << "Found " << header_id << endl;
	}
	
	/*	
	 YDR093W  	              	         	
	 ATAGACAAACCAGATTCATCGCTTA     	
	 GTGTGCTAATATCTTTCATTGGCTT     	
	 AGATCGTTAGAGAAATGTGGTTACG     	
	 CCGTAAAGTAGTTACTTGCCACAAT
	 */	
	
	/// FOR S. CEREVISIAE *****************************
	
	ifstream fin2("Eurofan_CPs.txt");
	if(!fin2.is_open())
	{
		//cout << "unable to open confirmation_primers_S_cere.txt \n";
		cout << "unable to open Eurofan_CPs.txt \n";
		return(0);
	}
	
	found = FALSE;
	
	while(fin2.getline(buffer, 2055))
	{
		strcpy(tbuffer, buffer);
		
		strtok(buffer,"\t\n\r");
		strcpy(temp_A, strtok(NULL,"\t\n\r"));
		strcpy(header_id, strtok(temp_A," "));
		
		//cout << "*" << header_id << "*" << endl;
		
		if(!strcmp(header_id, query_gene))
		{
			found = TRUE;
			//cout << "FOund: " << header_id << endl;
			strtok(tbuffer,"\t\n\r");
			strtok(NULL,"\t\n\r");
			//strtok(NULL,"\t\n\r");
			strcpy(temp_A, strtok(NULL,"\t\n\r"));
			strcpy(temp_B, strtok(NULL,"\t\n\r"));
			strcpy(temp_C, strtok(NULL,"\t\n\r"));
			strcpy(temp_D, strtok(NULL,"\t\n\r"));
			
			//for(i = 0; i < 30; i++) temp[i] = tbuffer[i + 49];
			// strcpy(CPS_A, strtok(temp," \t\n\r"));
			 
			// for(i = 0; i < 30; i++) temp[i] = tbuffer[i + 145];
			// strcpy(CPS_D, strtok(temp," \t\n\r"));
			break;
		}
	} 
	 
	 if(!found) return(-1);
	 
	 strcpy(CPS_A, strtok(temp_A," "));
	 strcpy(CPS_B, strtok(temp_B," "));
	 strcpy(CPS_C, strtok(temp_C," "));
	 strcpy(CPS_D, strtok(temp_D," "));
	 
	
	/// FOR S. POMBE *****************************
	/*
	ifstream fin2("PombConf_v1.csv");
	if(!fin2.is_open())
	{
		//cout << "unable to open confirmation_primers_S_cere.txt \n";
		cout << "unable to open PombConf_v1.csv \n";
		return(0);
	}
	
	found = FALSE;
	
	while(fin2.getline(buffer, 2055))
	{
		strcpy(tbuffer, buffer);		
		strcpy(header_id, strtok(buffer,",\t\n\r"));
		
		//cout << "*" << header_id << "*" << endl;
		
		if(!strcmp(header_id, query_gene))
		{
			found = TRUE;
			//cout << "FOund: " << header_id << endl;
			strtok(tbuffer,",\t\n\r");
			strcpy(temp_A, strtok(NULL,",\t\n\r"));
			strtok(NULL,",\t\n\r");
			strtok(NULL,",\t\n\r");
			strtok(NULL,",\t\n\r");
			strtok(NULL,",\t\n\r");
			strtok(NULL,",\t\n\r");
			strtok(NULL,",\t\n\r");
			strcpy(temp_D, strtok(NULL,",\t\n\r"));
			
			break;
		}
	} 
	
	if(!found) return(-1);
	
	strcpy(CPS_A, strtok(temp_A," "));
	//strcpy(CPS_B, strtok(temp_B," "));
	//strcpy(CPS_C, strtok(temp_C," "));
	strcpy(CPS_D, strtok(temp_D," "));*/
	
	// END OF CONFIRMATION FILE READERs ********************
	
	
	
	strcpy(cnf_results,"\t\t\t<confirmationPrimerForward>");
	strcat(cnf_results, CPS_A);
	strcat(cnf_results, "</confirmationPrimerForward>\n");
	
	strcat(cnf_results, "\t\t\t<confirmationPrimerReverse>");
	strcat(cnf_results, CPS_D);
	strcat(cnf_results, "</confirmationPrimerReverse>\n");
	
	cout << "Length of ORF: " << strlen(orf) << endl;
	
	int forward_at = findstring(CPS_A, orf); //cout << "Forward strand\n";
	reverse_complement(orf, RC_orf);
	int reverse_at = findstring(CPS_D, RC_orf); //cout << "Reverse strand\n";
	
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
	int gene_length = strlen(orf) - 2000;
	
	// Make sequences
	int x = 0;
	int y = 0;
	int orf_length = strlen(orf);
	
	for(i = forward_at; i < orf_length - reverse_at; i++)
	{
		gene_present_confirmation_sequence[x++] = orf[i];
		
		if(i < 1000 || i >= orf_length - 1000)
			gene_deleted_confirmation_sequence[y++] = orf[i];
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
	
	/***********************************************************************\
	 
	 Now need to find these:
	 
	 FinklURA3 5'-ACGTTCGTTCGACTGATGAG-3'
	 RinklURA3 5'-CTCATCAGTCGAACGAACGT-3'
	 
	 in the URA3 sequence between Primer D and Primer E
	 
	 PRIMER_D "GATCCCAATACAACAGAT" 
	 PRIMER_E "CACACATTACTTGCCTC"
	 
	 \***********************************************************************/
	
	// Get URA3
	
	fin.open("ura3_rc.fa");
	if(!fin.is_open())
	{
		cout << "unable to open ura3_rc.fa \n";
		return(0);
	}
	
	// Get Sequence
	
	x = 0; 
	
	while(fin.getline(buffer, 2055))
	{		
		if(buffer[0] != 0x3E)  // fasta header
		{
			int buffer_length = strlen(buffer);
			
			for(i = 0; i < buffer_length; i++)
			{
				if(buffer[i] > 31) // remove ctrl chars
					ura3[x++] = buffer[i];
			}
		}
	} 
	ura3[x] = 0; // terminate string
	
	fin.close();
	
	// Find the cassette between E and D
	char RC_ura3[20000];
	int ura3_length = strlen(ura3);
	
	int Primer_D_at = findstring(Primer_D, ura3); //cout << "Primer Forward strand\n";
	reverse_complement(ura3, RC_ura3);
	int Primer_E_at = findstring(Primer_E, RC_ura3); //cout << "Primer Reverse strand\n";
	
	//cout << endl << ura3 << " Length is " << strlen(ura3) << endl;
	
	char FinklURA3[40] = "ACGTTCGTTCGACTGATGAG";
	char RinklURA3[40] = "CTCATCAGTCGAACGAACGT";
	
	//cout << "Length of URA3: " << strlen(ura3) << endl;
	int Rinkl_at = findstring(RinklURA3, ura3); //cout << "Forward strand\n";
	reverse_complement(ura3, RC_ura3);
	int Finkl_at = findstring(FinklURA3, RC_ura3); //cout << "Reverse strand\n";
	
	//cout << "RinklURA3: " << forward_at << endl ;
	//cout << "FinklURA3: " << reverse_at << endl;
	
	int E_to_Rinkl = ura3_length - Rinkl_at - Primer_E_at;
	int D_to_Finkl = ura3_length - Finkl_at - Primer_D_at;
	
	//cout << "E to Rinkl: " << E_to_Rinkl << endl;
	//cout << "D to Finkl: " << D_to_Finkl << endl;
	
	// Get URA3 sequence
	char URA3_sequence[20000];
	
	x = 0;
	
	for(i = Primer_E_at; i < (ura3_length - Primer_D_at); i++) 
		URA3_sequence[x++] = RC_ura3[i];
	URA3_sequence[x] = 0;
	
	// Get sequence from CPS A to the start of the gene
	char A_to_start_sequence[1000];
	
	x = 0;
	
	for(i = forward_at; i < 1000; i++) 
		A_to_start_sequence[x++] = orf[i];
	A_to_start_sequence[x] = 0;
	
	// Get the sequence from end of gene to CPS D
	char D_to_stop_sequence[1000];
	
	x = 0;
	
	for(i = orf_length - 1000; i < orf_length - reverse_at; i++) D_to_stop_sequence[x++] = orf[i];
	D_to_stop_sequence[x] = 0;
	
	/*	
	 fout << "CPS A to CPS D with URA3 present: " << A_to_start + D_to_stop + strlen(URA3_sequence) << " - ";
	 fout << A_to_start_sequence << URA3_sequence << D_to_stop_sequence << endl ;
	 
	 fout << "Finkl: " << FinklURA3 << endl;
	 fout << "Rinkl: " << RinklURA3 << endl;
	 */
	sprintf(tempstr, "\t\t\t<ura3PresentLength>%d</ura3PresentLength>\n", A_to_start + D_to_stop + (int)strlen(URA3_sequence));
	strcat(cnf_results, tempstr);
	sprintf(tempstr, "\t\t\t<ura3PresentProduct>%s%s%s</ura3PresentProduct>\n", A_to_start_sequence, URA3_sequence, D_to_stop_sequence);
	strcat(cnf_results, tempstr);
	
	sprintf(tempstr, "\t\t\t<finkl>%s</finkl>\n\t\t\t<rinkl>%s</rinkl>\n", FinklURA3, RinklURA3);
	strcat(cnf_results, tempstr);
	
	
	// Get Primer E to Rinkl sequence on RC_ura3
	char E_to_Rinkl_sequence[4000];
	
	x = 0;
	
	for(i = Primer_E_at; i < (ura3_length - Rinkl_at); i++) E_to_Rinkl_sequence[x++] = RC_ura3[i];
	E_to_Rinkl_sequence[x] = 0;
	
	sprintf(tempstr, "\t\t\t<cpsA_rinklLength>%d</cpsA_rinklLength>\n", A_to_start + 40 + E_to_Rinkl);
	strcat(cnf_results, tempstr);
	sprintf(tempstr, "\t\t\t<cpsA_rinklProduct>%s%s%s</cpsA_rinklProduct>\n", A_to_start_sequence, Primer_R, E_to_Rinkl_sequence);
	strcat(cnf_results, tempstr);
	
	// Get Primer D to Finkl sequence on ura3
	char D_to_Finkl_sequence[4000];
	
	x = 0;
	
	for(i = Primer_D_at; i < (ura3_length - Finkl_at); i++) D_to_Finkl_sequence[x++] = ura3[i];
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
	return(1);
}

int auto_process(const char* gene, const char* orf, const char* genome_file_name, ofstream &fout)
{
	char buffer[4000];
	char Primer_A[62];
	char Primer_C[47], Primer_R[42];
	char RC_Primer_B[62], RC_Primer_A[62];
	char PCR1_product[1000];
	
	int i, j;
	int idxP_A = 0;
	int idxP_B = 0;
	
	char ura3_rc[4000];
	char ura3[4000];
	char URA3_sequence[4000];
		
	ifstream fin("ura3_rc.fa");
	
#ifdef NSB	
	DNAfind s_cere_nsb(genome_file_name);
	
	if(!s_cere_nsb.set_max_mismatches(0)) 
		cout << "Could not set max_mismatches\n";
	
	if(!s_cere_nsb.set_tail_length(12)) 
		cout << "Could not set tail length\n";
	
	DNAfind ura_3_nsb("ura3_rc.fa");
	
	if(!ura_3_nsb.set_max_mismatches(0)) 
		cout << "Could not set max_mismatches\n";
	
	if(!ura_3_nsb.set_tail_length(8)) 
		cout << "Could not set tail length\n";
	
#endif
	
	if(!fin.is_open())
	{
		cout << "Unable to open ura3_rc.fa" << endl;
		return(NO_URA3);
	}
	
	// Get Sequence
	
	int x = 0; 
	
	while(fin.getline(buffer, 4000))
	{		
		if(buffer[0] != 0x3E)  // fasta header
		{
			int buffer_length = strlen(buffer);
			
			for(i = 0; i < buffer_length; i++)
			{
				if(buffer[i] > 31) // remove ctrl chars
					ura3_rc[x++] = buffer[i];
			}
		}
	} ura3_rc[x] = 0; // terminate string
	
	reverse_complement(ura3_rc, ura3);
	
	fin.close();
	
		
	//sequence_analysis primer_b_analysis;		
	//primer_candidate_generation primer_b(REVERSE);
	
	primer primer_b;
	
	int start = start_codon_location(orf);
	
	//primer_b.tail_complementarity_check = FALSE;
	primer_b.reverse_primer = TRUE;
	primer_b.start_location_range_begin = start - 1; // Primer B starts 
	primer_b.start_location_range_end = start - 1;   // immediately before start codon
	primer_b.length_range_shortest = 30;
	primer_b.length_range_longest = 60;
		
	if(!primer_b.generate_candidates(orf))return(NO_PRIMER_B);
	
	for(i = 0; i < primer_b.candidates_found; i++)
	{
		//cout << "Checking candidate " << i << endl;
		
		//primer_b.tail_complementarity_report(primer_b.candidate[i].sequence);
		
		primer_b.hairpin(i);
		primer_b.self_dimer(i);

		reverse_complement(primer_b.candidate[i].sequence, RC_Primer_B);
		
#ifdef BLAST
		primer_b.candidate[i].yeast_matches = sequence_utils::blast_yeast(RC_Primer_B, "RC_primer_B", primer_b.expectation);
		primer_b.candidate[i].plasmid_matches = sequence_utils::blast_plasmid(RC_Primer_B, "RC_primer_B", primer_b.expectation);
#endif		

#ifdef FASTA
		primer_b.candidate[i].yeast_matches = sequence_utils::fasta3(RC_Primer_B, "s_cere_genome.fa", primer_b.expectation);
		primer_b.candidate[i].plasmid_matches = sequence_utils::fasta3(RC_Primer_B, "ura3_rc.fa", primer_b.expectation);
#endif
		
#ifdef NSB
		primer_b.candidate[i].yeast_matches = s_cere_nsb.search_for_binding_sites(RC_Primer_B);
		primer_b.candidate[i].plasmid_matches = ura_3_nsb.search_for_binding_sites(RC_Primer_B);
#endif
	}	
	
	// Auto select Primer B
	// Set up sorting priorities - higher the array index, the higher the priority, so 0 is lowest.
	primer_b.priority[0] = LENGTH;
	primer_b.priority[1] = PLASMID_MATCH;
	primer_b.priority[2] = YEAST_MATCH;
	primer_b.priority[3] = HAIRPIN;
	primer_b.priority[4] = SELF_DIMER;
	primer_b.priority[5] = SORT_END;
	
	cout << "Rank selection - "  << primer_b.candidates_found << " candidates found\n";
	
	primer_b.rank_selection();
	primer_b.show_all_candidates();	
	
	bool A_B_OK = FALSE;
	idxP_B = 0;
	
	
// NOW DO PRIMER A	
	
	primer primer_a;
	
	//primer_a.tail_complementarity_check = FALSE;
	primer_a.start_location_range_begin = 0;
	primer_a.start_location_range_end = 500;
	primer_a.downstream_search = TRUE; 
	
	primer_a.max_number_candidates = 100;
	
	while(!A_B_OK)
	{	
		primer_a.length_range_shortest = strlen(primer_b.candidate[idxP_B].sequence);
		primer_a.length_range_longest = strlen(primer_b.candidate[idxP_B].sequence);
	
		// Set required GC content for primer A
		primer_a.required_GC_content = primer_b.candidate[idxP_B].get_GC_content();
		primer_a.GC_clamping = FALSE;
		//primer_a.GC_tolerance = 1;
	
		if(primer_a.generate_candidates(orf))
		{
#ifdef TESTING				
			cout << "Number of candidates = " << primer_a.candidates_found << endl;
#endif
			for(i = 0; i < primer_a.candidates_found; i++)
			{
				//cout << "Checking candidate " << i << endl;
				//primer_a.tail_complementarity_report(primer_a.candidate[i].sequence);
#ifdef TESTING				
				cout << "Sequence " << i << " = " << primer_a.candidate[i].sequence << endl;
#endif				
				primer_a.hairpin(i);	
				primer_a.self_dimer(i);
				primer_a.primer_dimer_2(i, primer_b.candidate[idxP_B].sequence);
			
				reverse_complement(primer_a.candidate[i].sequence, RC_Primer_A);
				reverse_complement(primer_b.candidate[idxP_B].sequence, RC_Primer_B);
				
#ifdef TESTING				
				cout << "DImer tests done \n";
#endif	
				
#ifdef BLAST
				primer_a.candidate[i].yeast_matches = sequence_utils::blast_yeast(RC_Primer_A, "RC_primer_A", primer_a.expectation);
				primer_a.candidate[i].plasmid_matches = sequence_utils::blast_plasmid(RC_Primer_A, "RC_primer_A", primer_a.expectation);
#endif
				
#ifdef FASTA
				primer_a.candidate[i].yeast_matches = sequence_utils::fasta3(RC_Primer_A, "s_cere_genome.fa", primer_a.expectation);
				primer_a.candidate[i].plasmid_matches = sequence_utils.fasta3(RC_Primer_A, "ura3_rc.fa", primer_a.expectation);
#endif
				
#ifdef NSB 
				primer_a.candidate[i].yeast_matches = s_cere_nsb.search_for_binding_sites(RC_Primer_A);
				primer_a.candidate[i].plasmid_matches = ura_3_nsb.search_for_binding_sites(RC_Primer_A);
#ifdef TESTING				
				cout << "Binding tests done \n";
#endif	
				
				// Now we have both primers we can check for non-specific products
				
				primer_a.candidate[i].products = s_cere_nsb.search_for_pcr_products(RC_Primer_A, RC_Primer_B);
				primer_a.candidate[i].products += ura_3_nsb.search_for_pcr_products(RC_Primer_A, RC_Primer_B);
				
#endif
				
#ifdef TESTING				
				cout << "Number of products for " << i << " = " << primer_a.candidate[i].products << endl;
#endif
			}
	
		// Auto select Primer A
		// Set up sorting priorities - higher the array index, the higher the priority, so 0 is lowest.
#ifdef TESTING			
			cout << "Auto select Primer A\n";
#endif
			
			if(primer_b.candidate[idxP_B].yeast_matches == 1)
			{
				primer_a.priority[0] = PLASMID_MATCH;
				primer_a.priority[1] = YEAST_MATCH;
				primer_a.priority[2] = HAIRPIN;
				primer_a.priority[3] = SELF_DIMER;
				primer_a.priority[4] = F_DIMER;
				primer_a.priority[5] = R_DIMER;
				primer_a.priority[6] = PRODUCTS;
				primer_a.priority[7] = MOO_SORT;
				primer_a.priority[8] = SORT_END;
		
			}
			else 
			{	
				primer_a.priority[0] = PLASMID_MATCH;
				primer_a.priority[1] = HAIRPIN;
				primer_a.priority[2] = SELF_DIMER;
				primer_a.priority[3] = F_DIMER;
				primer_a.priority[4] = R_DIMER;
				primer_a.priority[5] = YEAST_MATCH;
				primer_a.priority[6] = PRODUCTS;
				primer_a.priority[7] = MOO_SORT;
				primer_a.priority[8] = SORT_END;				
			}

			if(primer_a.rank_selection()) // We have at least one acceptable primer
			{
				primer_a.show_all_candidates();
		
				idxP_A = 0;
	
				strcpy(Primer_A, primer_a.candidate[idxP_A].sequence);
				cout << idxP_A << " is " << Primer_A << endl;
			
				A_B_OK = TRUE;
			}
			else 
				if(idxP_B++ >= primer_b.good_candidates)
				{
					cout << "Cannot find suitable primers A and B\n";
					return(NO_PRIMERS_A_AND_B);
				}			
		}
		else 
			if(idxP_B++ >= primer_b.good_candidates)
			{
				cout << "Cannot find suitable primers A and B\n";
				return(NO_PRIMERS_A_AND_B);
			}		
	}

	/*
	 * 	GET PRIMER C AND PRIMER R
	 */
	
	// 	C = 45 nt upstream of the stop codon and include stop
	//  R = 40 nt downstream of the stop codon not including the stop codon
	
	int stop = stop_codon_location(start, orf);
	// Copy primer C and primer R straight from the orf
	
   	for(i = 0; i < 45; i++) Primer_C[i] = orf[i + stop - 42];
   	Primer_C[45] = 0;
   	for(i = 0; i < 40; i++) Primer_R[i] = orf[i + stop +3];
   	Primer_R[40] = 0;
	
	cout << "Primer B: " << primer_b.candidate[idxP_B].sequence << endl;
	cout << "Primer A: " << primer_a.candidate[idxP_A].sequence << endl;
   	cout << "Primer C: " << Primer_C << endl;
	cout << "Primer R: " << Primer_R << endl;	
	
	//Make PCR1 product
	j = 0;
	
	for(i = primer_a.candidate[idxP_A].location_5_prime_end; i < start; i++)
		PCR1_product[j++] = orf[i];	
	PCR1_product[j] = 0;
	
// MAKE D AND E CANDIDATES
// *****************************************
	
	primer_data primer_d[6];
	primer_data primer_e[6];
		
	strcpy(primer_d[0].sequence, "GATCCCAATACAACAGAT"); primer_d[0].location_5_prime_end = 3892;// This the original @ 97
	strcpy(primer_d[1].sequence, "GTCTAGAGATCCCAATACA"); primer_d[1].location_5_prime_end = 3899;
	strcpy(primer_d[2].sequence, "CTAGAGATCCCAATACAAC"); primer_d[2].location_5_prime_end = 3897; // This looks good for YGR125W (17/08/10)
	strcpy(primer_d[3].sequence, ""); primer_d[3].location_5_prime_end = 0;
	strcpy(primer_d[4].sequence, ""); primer_d[4].location_5_prime_end = 0;
	strcpy(primer_d[5].sequence, ""); primer_d[5].location_5_prime_end = 0;
	
	strcpy(primer_e[0].sequence, "CACACATTACTTGCCTC"); primer_e[0].location_5_prime_end = 2736;// original
	strcpy(primer_e[1].sequence, "AGCGACAAGAAGAGATAG"); primer_e[1].location_5_prime_end = 2500;// (A)
	strcpy(primer_e[2].sequence, "ACAATCATATGGGAGAAG"); primer_e[2].location_5_prime_end = 2622;// (B)
	strcpy(primer_e[3].sequence, "ACACTCCTCAGAAGCTC"); primer_e[3].location_5_prime_end = 2545;// (C) @ 1443 Good for YGL202W (17/08/10)
	strcpy(primer_e[4].sequence, ""); primer_e[4].location_5_prime_end = 0;
	strcpy(primer_e[5].sequence, ""); primer_e[5].location_5_prime_end = 0;	

	
// PCR2 PRIMERS
// *****************************************
	
	cout << "PCR2 primers\n";
	
	char PCR2_forward[256];
	char PCR2_reverse[256];
	char RC_PCR2_forward[256];
	char RC_PCR2_reverse[256];
	char RC_Primer_C[47];
	char PCR1_homology[45];

	primer_data *PCR2_forward_data_p;
	primer_data *PCR2_reverse_data_p;

	primer pcr2f;
	primer pcr2r;
		
	// Get PCR1 product/PCR2 forward primer homology sequence
	
#define PCR1_INIT_HOMOLOGY_LENGTH 25	
#define PCR1_HOMOLOGY_LENGTH 40
	
	for(i = 0; i < PCR1_INIT_HOMOLOGY_LENGTH; i++) 
		PCR1_homology[i] = PCR1_product[strlen(PCR1_product) - PCR1_INIT_HOMOLOGY_LENGTH + i];
	PCR1_homology[PCR1_INIT_HOMOLOGY_LENGTH] = 0; // terminate string
	
	// So the homology sequence is only 25 bases unless the GC content < 40%, then make it 40 bases:
	
	int GC_number = 0;
	int PCR1_homology_length = strlen(PCR1_homology);
	
	for(i = 0; i < PCR1_homology_length; i++)
		if(PCR1_homology[i] == 67 || PCR1_homology[i] == 71) GC_number++;
	
	double GC_content = (double)GC_number/(double)PCR1_homology_length;
	
	if(GC_content < 0.4)
	{
		for(i = 0; i < PCR1_HOMOLOGY_LENGTH; i++) 
			PCR1_homology[i] = PCR1_product[strlen(PCR1_product) - PCR1_HOMOLOGY_LENGTH + i];
		PCR1_homology[PCR1_HOMOLOGY_LENGTH] = 0;
	}
	// Make Reverse complement Primer C
	reverse_complement(Primer_C, RC_Primer_C);
	
	bool D_E_OK = FALSE;
	bool FORWARD_HAIRPIN_OK;
	bool FORWARD_SELF_DIMER_OK;
	bool REVERSE_HAIRPIN_OK;
	bool REVERSE_SELF_DIMER_OK;
	bool PAIR_DIMER_OK;
	
	
	int idxP_D = 0; // idxP = index for primer...
	int idxP_E = 0;
	
	while(!D_E_OK)
	{
		// Make PRC2 forward   	
		strcpy(PCR2_forward, PCR1_homology);
		strcat(PCR2_forward, Primer_R);
		strcat(PCR2_forward, primer_e[idxP_E].sequence);  
   	
		// Make PCR2 reverse 	   	
		strcpy(PCR2_reverse, RC_Primer_C);
		strcat(PCR2_reverse, primer_d[idxP_D].sequence);
   	
		cout << "PCR2 forward: " << PCR2_forward << endl;
		cout << "PCR2 reverse: " << PCR2_reverse << endl;
	
		strcpy(pcr2f.candidate[0].sequence, PCR2_forward);
		strcpy(pcr2r.candidate[0].sequence, PCR2_reverse);
		pcr2f.hairpin(0);
		pcr2r.hairpin(0);

		pcr2f.hairpin(0);
		if(pcr2f.candidate[0].hairpin < 12) 
			FORWARD_HAIRPIN_OK = TRUE;
		else 
			FORWARD_HAIRPIN_OK = FALSE;
		
		pcr2f.self_dimer(0);
		if(pcr2f.candidate[0].self_dimer < 12) 
			FORWARD_SELF_DIMER_OK = TRUE;
		else 
			FORWARD_SELF_DIMER_OK = FALSE;
		
		pcr2r.hairpin(0);
		if(pcr2r.candidate[0].hairpin < 12) 
			REVERSE_HAIRPIN_OK = TRUE;
		else 
			REVERSE_HAIRPIN_OK = FALSE;
		
		pcr2r.self_dimer(0);
		if(pcr2r.candidate[0].self_dimer < 12) 
			REVERSE_SELF_DIMER_OK = TRUE;
		else 
			REVERSE_SELF_DIMER_OK = FALSE;

		pcr2r.primer_dimer_2(0, pcr2f.candidate[0].sequence);
		if(pcr2r.candidate[0].forward_dimer < 12 && pcr2r.candidate[0].reverse_dimer < 12) 
			PAIR_DIMER_OK = TRUE;
		else 
			PAIR_DIMER_OK = FALSE;
		
		reverse_complement(PCR2_forward, RC_PCR2_forward);
		reverse_complement(PCR2_reverse, RC_PCR2_reverse);
		strcpy(pcr2f.candidate[0].sequence, RC_PCR2_forward);
		strcpy(pcr2r.candidate[0].sequence, RC_PCR2_reverse);
		
#ifdef BLAST
		pcr2f.blast_plasmid(0, "PCR2_forward");	
		pcr2f.blast_yeast(0, "PCR2_forward");
	
		pcr2r.blast_plasmid(0, "PCR2_reverse");	
		pcr2r.blast_yeast(0, "PCR2_reverse");
#endif

#ifdef NSB
		pcr2f.candidate[0].yeast_matches = s_cere_nsb.search_for_binding_sites(RC_PCR2_forward);
		pcr2f.candidate[0].plasmid_matches = ura_3_nsb.search_for_binding_sites(RC_PCR2_forward);
		
		pcr2r.candidate[0].yeast_matches = s_cere_nsb.search_for_binding_sites(RC_PCR2_reverse);
		pcr2r.candidate[0].plasmid_matches = ura_3_nsb.search_for_binding_sites(RC_PCR2_reverse);
		
		pcr2r.candidate[0].products = s_cere_nsb.search_for_pcr_products(RC_PCR2_forward, RC_PCR2_reverse);
		pcr2f.candidate[0].products = ura_3_nsb.search_for_pcr_products(RC_PCR2_forward, RC_PCR2_reverse);

		int total_products = pcr2r.candidate[0].products + pcr2f.candidate[0].products;
#endif

		if(FORWARD_HAIRPIN_OK && 
		   FORWARD_SELF_DIMER_OK && 
		   REVERSE_HAIRPIN_OK && 
		   REVERSE_SELF_DIMER_OK && 
		   PAIR_DIMER_OK)
		{
			D_E_OK = TRUE;
		}
		else 			
		{
			if(!FORWARD_HAIRPIN_OK || !FORWARD_SELF_DIMER_OK) idxP_E++;
			else if(!REVERSE_HAIRPIN_OK || !REVERSE_SELF_DIMER_OK) idxP_D++;
			else if(!PAIR_DIMER_OK)
			{
				if(idxP_E < 3)
					idxP_E++;
				else 
					idxP_D++;
			}

			else 
			{
				if(pcr2r.candidate[0].forward_dimer > 11) idxP_D++;
				if(pcr2r.candidate[0].reverse_dimer > 11) idxP_E++;
#ifdef NSB
				if(total_products > 1) 
				{
					if(idxP_E < 3)
						idxP_E++;
					else 
						idxP_D++;
				}
#endif
			}
		}
		
		if(idxP_D > 2)return(NO_PRIMER_D);
		if(idxP_E > 3)return(NO_PRIMER_E);

	}
	
// Report D and E to console
	cout << "Primer D: " << primer_d[idxP_D].sequence << endl;
	cout << "Primer E: " << primer_e[idxP_E].sequence << endl;
	
// Make the URA3 cassette sequence
	x = 0;
	
	for(i = primer_e[idxP_E].location_5_prime_end; i <= primer_d[idxP_D].location_5_prime_end; i++) URA3_sequence[x++] = ura3[i];
	URA3_sequence[x] = 0; // terminate string

// Confirmation primers
	
	char cnf_results[20000] = "None";
	
	if(XML)
	{
		if(!xml_confirmation_primers(gene, primer_d[idxP_D].sequence, primer_e[idxP_E].sequence, Primer_R, RC_Primer_C, cnf_results))
			return(CONF_PRIMER_ERROR);
	}
	else 
	{
		if(!confirmation_primers(gene, primer_d[idxP_D].sequence, primer_e[idxP_E].sequence, Primer_R, RC_Primer_C, cnf_results))
			return(CONF_PRIMER_ERROR);
	}
	 
	
// RESULTS
// **************************************
	time_t mtime;
	time(&mtime);
	char output[64];
	char mytime[64];
	
	annealing_temperature Tm;
	
	strcpy(mytime, ctime(&mtime));
	
	mytime[strlen(mytime) - 1] = 0;
	
	PCR2_forward_data_p = &pcr2f.candidate[0];
	PCR2_reverse_data_p = &pcr2r.candidate[0];

	if(XML)
	{
		
//#ifndef BEOWULF		
		fout << "<organism>\n\t<name>Schizosaccharomyces pombe</name>\n";
//#endif
		fout << "\t<geneForDeletion>\n\t\t<name>" << gene << "</name>\n\t\t<date>" << mytime << "</date>\n";
		fout << "\t\t<class>" << PCR_RESULTS_CLASS << "</class>\n";
		
	// PCR1	
		fout << "\t\t<primerPair>\n";
		fout << "\t\t\t<name>PCR1</name>\n";
				
		fout << "\t\t\t<forwardPrimer>\n";
		fout << "\t\t\t\t<sequence>" << primer_a.candidate[idxP_A].sequence << "</sequence>\n";
		fout << "\t\t\t\t<hairpinScore>" << primer_a.candidate[idxP_A].hairpin << "</hairpinScore>\n";
		fout << "\t\t\t\t<selfDimerScore>" << primer_a.candidate[idxP_A].self_dimer << "</selfDimerScore>\n";
		fout << "\t\t\t\t<bindingSites>\n";
		fout << "\t\t\t\t\t<yeast>" << primer_a.candidate[idxP_A].yeast_matches << "</yeast>\n";
		fout << "\t\t\t\t\t<ura3>" << primer_a.candidate[idxP_A].plasmid_matches << "</ura3>\n";
		fout << "\t\t\t\t</bindingSites>\n";
		fout << "\t\t\t\t<annealingTemperature>" << Tm.primer3_Tm(primer_a.candidate[idxP_A].sequence) << "</annealingTemperature>\n";
		fout << "\t\t\t</forwardPrimer>\n";
		
		fout << "\t\t\t<reversePrimer>\n";
		fout << "\t\t\t\t<sequence>" << primer_b.candidate[idxP_B].sequence << "</sequence>\n";
		fout << "\t\t\t\t<hairpinScore>" << primer_b.candidate[idxP_B].hairpin << "</hairpinScore>\n";
		fout << "\t\t\t\t<selfDimerScore>" << primer_b.candidate[idxP_B].self_dimer << "</selfDimerScore>\n";
		fout << "\t\t\t\t<bindingSites>\n";
		fout << "\t\t\t\t\t<yeast>" << primer_b.candidate[idxP_B].yeast_matches << "</yeast>\n";
		fout << "\t\t\t\t\t<ura3>" << primer_b.candidate[idxP_B].plasmid_matches << "</ura3>\n";
		fout << "\t\t\t\t</bindingSites>\n";
		fout << "\t\t\t\t<annealingTemperature>" << Tm.primer3_Tm(primer_b.candidate[idxP_B].sequence) << "</annealingTemperature>\n";
		fout << "\t\t\t</reversePrimer>\n";
		
		fout << "\t\t\t<primerDimerForward>" << primer_a.candidate[idxP_A].forward_dimer << "</primerDimerForward>\n";
		fout << "\t\t\t<primerDimerReverse>" << primer_a.candidate[idxP_A].reverse_dimer << "</primerDimerReverse>\n";
		
		fout << "\t\t\t<numberOfProducts>" << primer_a.candidate[idxP_A].products << "</numberOfProducts>\n";
		
		fout << "\t\t\t<product>" << PCR1_product << "</product>\n";
		fout << "\t\t\t<productLength>" << start - primer_a.candidate[idxP_A].location_5_prime_end << "</productLength>\n";
		
		fout << "\t\t</primerPair>\n";
	
	// PCR2
		fout << "\t\t<primerPair>\n";
		fout << "\t\t\t<name>PCR2</name>\n";
		
		fout << "\t\t\t<forwardPrimer>\n";
		fout << "\t\t\t\t<sequence>" << PCR2_forward << "</sequence>\n";
		fout << "\t\t\t\t<hairpinScore>" << PCR2_forward_data_p->hairpin << "</hairpinScore>\n";
		fout << "\t\t\t\t<selfDimerScore>" << PCR2_forward_data_p->self_dimer << "</selfDimerScore>\n";
		fout << "\t\t\t\t<bindingSites>\n";  
		fout << "\t\t\t\t\t<yeast>" << PCR2_forward_data_p->yeast_matches << "</yeast>\n";
		fout << "\t\t\t\t\t<ura3>" << PCR2_forward_data_p->plasmid_matches << "</ura3>\n";		
		fout << "\t\t\t\t</bindingSites>\n";
		fout << "\t\t\t\t<annealingTemperature>" << Tm.primer3_Tm(PCR2_forward) << "</annealingTemperature>\n";
		
		fout << "\t\t\t</forwardPrimer>\n";
		
		fout << "\t\t\t<reversePrimer>\n";
		fout << "\t\t\t\t<sequence>" << PCR2_reverse << "</sequence>\n";
		fout << "\t\t\t\t<hairpinScore>" << PCR2_reverse_data_p->hairpin << "</hairpinScore>\n";
		fout << "\t\t\t\t<selfDimerScore>" << PCR2_reverse_data_p->self_dimer << "</selfDimerScore>\n";
		fout << "\t\t\t\t<bindingSites>\n";  
		fout << "\t\t\t\t\t<yeast>" << PCR2_reverse_data_p->yeast_matches << "</yeast>\n";
		fout << "\t\t\t\t\t<ura3>" << PCR2_reverse_data_p->plasmid_matches << "</ura3>\n";
		fout << "\t\t\t\t</bindingSites>\n";
		fout << "\t\t\t\t<annealingTemperature>" << Tm.primer3_Tm(PCR2_reverse) << "</annealingTemperature>\n";
		
		fout << "\t\t\t</reversePrimer>\n";
		
		fout << "\t\t\t<primerDimerForward>" << PCR2_reverse_data_p->forward_dimer << "</primerDimerForward>\n";
		fout << "\t\t\t<primerDimerReverse>" << PCR2_reverse_data_p->reverse_dimer << "</primerDimerReverse>\n";
		
		fout << "\t\t\t<numberOfProducts>" << PCR2_reverse_data_p->products + PCR2_forward_data_p->products << "</numberOfProducts>\n";
		fout << "\t\t\t<primer_B*>" << PCR1_homology << "</primer_B*>\n";
		fout << "\t\t\t<primer_R>" << Primer_R  << "</primer_R>\n";
		fout << "\t\t\t<primer_E>" << primer_e[idxP_E].sequence << "</primer_E>\n";
		fout << "\t\t\t<primer_C>" << Primer_C  << "</primer_C>\n";
		fout << "\t\t\t<primer_D>" << primer_d[idxP_D].sequence << "</primer_D>\n";
		fout << "\t\t\t<product>" << PCR1_homology << Primer_R << URA3_sequence << Primer_C << "</product>\n";
		fout << "\t\t\t<productLength>" << strlen(PCR1_homology) + strlen(Primer_R) + strlen(URA3_sequence) + strlen(Primer_C) << "</productLength>\n";
		
		
		fout << "\t\t</primerPair>\n";
	// PCR3
		fout << "\t\t<pcr3Product>" << PCR1_product << Primer_R << URA3_sequence << Primer_C << "</pcr3Product>\n";
		fout << "\t\t<pcr3ProductLength>" << strlen(PCR1_product) + strlen(Primer_R) + strlen(URA3_sequence) + strlen(Primer_C) << "</pcr3ProductLength>\n"; 
		
	// CONFIRMATION RESULTS
		fout << "\t\t<confirmationData>\n" << cnf_results << "\t\t</confirmationData>\n";
		
		fout << "\t</geneForDeletion>\n";
//#ifndef BEOWULF
		fout << "</organism>\n";
//#endif
	}
	else
	{
		sprintf(output,"%s.results", gene);
		ofstream fout(output);
		
		fout << gene << " RESULTS SUMMARY - " << ctime(&mtime) << "==================================================\n\n";
		fout << "Name - No.: Sequence, Hpin, Sdimer, F dimer, R dimer, Yeast NSB, URA3 NSB, Location\n";
		fout << "Primer B - ";
		primer_b.show_candidate(idxP_B, fout);
		fout << "Primer A - ";
		primer_a.show_candidate(idxP_A, fout);
		fout << endl;
		
		fout << "Primer C: " << Primer_C << endl;
		fout << "Primer R: " << Primer_R << endl;	
		fout << "Primer D: " << primer_d[idxP_D].sequence << endl;
		fout << "Primer E: " << primer_e[idxP_E].sequence << endl;
		
		fout << "PCR1 product: " << PCR1_product << endl;
		fout << "PCR1 product location: " << primer_a.candidate[idxP_A].location_5_prime_end << " to " << start - 1 << endl;
		fout << "PCR1 product length: " << start - primer_a.candidate[idxP_A].location_5_prime_end << endl;
		
		fout << "PCR2 forward primer: " << PCR2_forward << endl;
		fout << "PCR2 reverse primer: " << PCR2_reverse << endl;
		
		fout << "PCR2 product: " << PCR1_homology << Primer_R << URA3_sequence << Primer_C << endl;
		fout << "PCR2 product length: " << strlen(PCR1_homology) + strlen(Primer_R) + strlen(URA3_sequence) + strlen(Primer_C) << endl;
		
		fout << "PCR3 product: " << PCR1_product << Primer_R << URA3_sequence << Primer_C << endl;
		fout << "PCR3 product length: " << strlen(PCR1_product) + strlen(Primer_R) + strlen(URA3_sequence) + strlen(Primer_C) << endl;
		
		fout << endl << gene << " CONFIRMATION PRIMERS\n============================\n\n";
		//confirmation_primers(gene, primer_d[idxP_D].sequence, primer_e[idxP_E].sequence, Primer_R, RC_Primer_C, fout);
		fout << cnf_results << endl;

		fout << endl << gene << " FULL RESULTS\n====================\n\n"; 
		fout << "PRIMER B\n";
		fout << "Hairpin \n";
		primer_b.hairpin(idxP_B, fout);
		fout << "\nSelf dimer \n";
		primer_b.self_dimer(idxP_B, fout);
		
		reverse_complement(primer_b.candidate[idxP_B].sequence, RC_Primer_B);
		fout << "\nBLAST results - Yeast\n";
		//primer_b.blast_yeast(RC_Primer_B, "RC_primer_B", fout);
		fout << "==============\n";
		
		fout << "BLAST results - URA3\n";
		//primer_b.blast_plasmid(RC_Primer_B, "RC_primer_B", fout);
		fout << "==============\n";
		
		fout << "PRIMER A\n";
		fout << "Hairpin \n";
		primer_a.hairpin(idxP_A, fout);
		fout << "\nSelf dimer \n";
		primer_a.self_dimer(idxP_A, fout);
		fout << "\nPrimer dimer A/B\n";
		primer_a.primer_dimer_2(idxP_A, primer_b.candidate[idxP_B].sequence, fout);
		
		reverse_complement(primer_a.candidate[idxP_A].sequence, RC_Primer_A);
		fout << "\nBLAST results - Yeast\n";
		//primer_a.blast_yeast(RC_Primer_A, "RC_primer_A", fout);
		fout << "==============\n";
		
		fout << "BLAST results - URA3\n";
		//primer_a.blast_plasmid(RC_Primer_A, "RC_primer_A", fout);
		fout << "==============\n";
		
		fout << "PCR2 FORWARD\n";
		fout << "Hairpin\n";
		///PCR_analysis.hairpin(PCR2_forward, PCR2_forward_data, fout);				
		pcr2f.hairpin(0, fout);				
		fout << "\nSelf dimer\n";
		///PCR_analysis.self_dimer(PCR2_forward, PCR2_forward_data, fout);
		pcr2f.self_dimer(0, fout);
		
		reverse_complement(PCR2_forward, RC_PCR2_forward);
		
		fout << "\nBLAST results - Yeast\n";
		//PCR_analysis.blast_yeast(RC_Primer_temp, "PCR2_forward", fout);
		fout << "==============\n";
		
		fout << "BLAST results - URA3\n";
		//PCR_analysis.blast_plasmid(RC_Primer_temp, "PCR2_forward", fout);
		fout << "==============\n";
		
		fout << "PCR2 REVERSE\n";
		fout << "Hairpin\n";
		///PCR_analysis.hairpin(PCR2_reverse, PCR2_reverse_data, fout);	
		pcr2r.hairpin(0, fout);	
		fout << "\nSelf dimer\n";
		///PCR_analysis.self_dimer(PCR2_reverse, PCR2_reverse_data, fout);	
		pcr2f.self_dimer(0, fout);	
		fout << "\nPrimer dimer - reverse/forward\n";
		///PCR_analysis.primer_dimer_2(PCR2_reverse, PCR2_forward, PCR2_reverse_data, fout);
		pcr2r.primer_dimer_2(0, PCR2_forward, fout);
		
		reverse_complement(PCR2_reverse, RC_PCR2_reverse);
		
		fout << "\nBLAST results - Yeast\n";	
		//PCR_analysis.blast_yeast(RC_Primer_temp, "PCR2_reverse", fout);
		fout << "==============\n";
		
		fout << "BLAST results - URA3\n";
		//PCR_analysis.blast_plasmid(RC_Primer_temp, "PCR2_reverse", fout);
		fout << "==============\n";
		
		fout << endl << gene << " SUPPLEMENTARY\n=====================\n\n"; 
		
		fout << "Primer B candidates\n";
		primer_b.show_all_candidates(fout);
		
		fout << "\nPrimer A candidates\n";
		primer_a.show_all_candidates(fout);
		
		fout.close();
	}
	
	return(1);	
}

/*int new_process(const char* gene, const char* orf)
{
	primrose_pair my_primers;
	
	int start = start_codon_location(orf);
	
	// set both primers parameters
	
	my_primers.forward_start_location_range_begin;
	my_primers.forward_start_location_range_end;
	my_primers.forward_length_range_shortest;
	my_primers.forward_length_range_longest;
	//my_primers.forward_required_GC_content;
	//my_primers.forward_GC_tolerance;
	
	my_primers.reverse_start_location_range_begin = start - 1; // Primer B starts 
	my_primers.reverse_start_location_range_end = start - 1;   // immediately before start codon
	my_primers.reverse_length_range_shortest = 30;
	my_primers.reverse_length_range_longest = 60;
	//my_primers.reverse_required_GC_content;
	//my_primers.reverse_GC_tolerance;
	
	// Abort this! It is not possible to determine the forward primer parameters
	// until after the reverse primer has been selected (07/04/2011)
	
}*/

int main(int argc, char** argv)
{
	char genebuffer[32];
	char buffer[4000];
	char orf[50000];
	char query_gene[32];
	char header_id[32];
	bool found;
	int error;
	char *token;
	char output_file_name[32];
	//char genes_1000_all[128] = "pombe_1000_all.fa";
	//char genome_file_name[128] = "pombe_genome_09052011.fasta";
	char genes_1000_all[128] = "orf_genomic_1000_all.fasta";
	char genome_file_name[128] = "s_cere_genome.fa";
	
   	cout << "Precise_Deletion: Primer design\n";
	cout << "Version 1.2, 22nd November 2011\n";
	
	if(!(argv[1] && argv[2]))
	{
		cout << "Use: exe <list of genes file> <missing genes output file>\n";
		return(0);
	}
	
	ofstream fout;
	
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
			strcpy(query_gene, token);
		
			if(!(query_gene[0] == 37)) // 37 = % => comment
			{
				fin.open(genes_1000_all);
				
				if(!fin.is_open())
				{
					cout << "unable to open " << genes_1000_all << endl;
					break;
				}
				
				found = FALSE;
				orf[0] = 0; 
				
				while(fin.getline(buffer, 2055))
				{
					if(buffer[0] == 0x3E)  // fasta header
					{
						if(found)break;
						strcpy(header_id, strtok(buffer,"> \t\n"));
						if(!strcmp(header_id, query_gene))found = TRUE;			
					}
					else if(found) strcat(orf, buffer);
					
					//if(i++ > 1000)break;
				}  
				fin.close();
				
				if(!found)
				{
					cout << "Could not find " << query_gene << endl;
					//initialization_error = TRUE;
					//return(FALSE);
				}
				else
				{
					cout << "Found " << header_id << endl;
					
					sprintf(output_file_name, "%s.xml", header_id);
					fout.open(output_file_name);
					
					error = auto_process(header_id, orf, genome_file_name, fout);
					//error = new_process(header_id, orf);
					
					switch(error)
					{
						case NO_PRIMER_B: ferr << header_id << ": No primer B\n"; break;
						case NO_PRIMERS_A_AND_B: ferr << header_id << ": No primer A\n"; break;
						case NO_PCR2_PRIMERS: ferr << header_id << ": No PCR2 primers \n"; break;
						case NO_PRIMER_D: ferr << header_id << ": No PCR2 primer D \n"; break;
						case NO_PRIMER_E: ferr << header_id << ": No PCR2 primer E \n"; break;
						case NO_URA3: ferr << header_id << ": Unable to open ura3_rc.fa \n"; break;
						case CONF_PRIMER_ERROR : ferr << header_id << ": Confirmation primers error\n"; break;
					}
					
					fout.close();
				}
			}
		} // end if(token)
	}
	
	fin.close();
	ferr.close();
	
    return(1);
}
