/**
 *  second_site.cpp
 *  
 *
 *  Created by Mhr on 16/02/2012.
 *  Copyright 2012 Michael C. Riley. All rights reserved.
 *
 */

#include <cstdio>
#include <ctime>
#include "../pd5/DNAfind.h"

using namespace std;

int main(int argc, char** argv)
{
/** We are going to search the genome of S. cerevisiae for potential secondary products
 */
	int hits = 0;
	clock_t start, stop;
	
	DNAfind my_dnaf("s_cere_genome.fa");
	my_dnaf.set_max_mismatches(5);
	my_dnaf.set_tail_length(20);
	my_dnaf.set_max_viable_product_length(3500);
	
	my_dnaf.GC_array_optimisation = TRUE;
	my_dnaf.report_details = TRUE;
	
	// Testing YDR093W primers from /precise_deletion/First_11_primers.txt
	//
	// Testing YMR272C primers from /precise_deletion/First_11_primers.txt
	start = clock();
	//hits = my_dnaf.search_for_pcr_products("CTTCGTACTAGTTTTAGCTTACCACCGCTTTTAGTG", "GCTTACATTTGCAGACGTCAAAAAACCTGACAAACC");
	//hits = my_dnaf.search_for_pcr_products("ATATGTGTTGAGTTTATGCTATGATTCTTTTTAAAGAATTCTATTATTAAAAGGATTAGG", "GTGATCATTTCTGTATAAGTTTTCTCAACTGAATATATAAATATATATCATTACTATGCT");
	/*
	cout << "YLR079W, " << endl; 
	hits = my_dnaf.search_for_pcr_products("TTTCGTGTAATAGTCCCTGCTTCAAGGGTCAAAATTCCGTAGAGGTTTGGTTG", "TCATGATGGTCTCTGATGTACTAAAGGGGTTATCGGAGTCCATAACATGCGAT"); //(YLR079W - Primer B, Primer A)
	cout << "Hits = " << hits << endl;
	
	cout << "YMR272C, " << endl; 
	hits = my_dnaf.search_for_pcr_products("CTTCGTACTAGTTTTAGCTTACCACCGCTTTTAGTG", "GCTTACATTTGCAGACGTCAAAAAACCTGACAAACC"); //(Primer B, Primer A)
	cout << "Hits = " << hits << endl;
	
	cout << "YDR093W, " << endl; 
	hits = my_dnaf.search_for_pcr_products("ATATGTGTTGAGTTTATGCTATGATTCTTTTTAAAGAATTCTATTATTAAAAGGATTAGG", "GTGATCATTTCTGTATAAGTTTTCTCAACTGAATATATAAATATATATCATTACTATGCT");
	cout << "Hits = " << hits << endl;
	
	cout << "YGL202W, " << endl; 
	hits = my_dnaf.search_for_pcr_products("GATAGTAACGATCGGTTGTCCTATTCAATGTCTGTATC", "AAGACATTAACATCATTTTTCGTGCTGTTAGGGTGGTG");
	cout << "Hits = " << hits << endl;
	*/
	cout << "YHR163W, " << endl; 
	hits = my_dnaf.search_for_pcr_products("TTGTCGCGGACCAGGGGTAAAGCAATTGATGGTGCATTGCCTTC", "CAAGACGGTGCATTTCTGGGCTCCTACTTTGAAATGGGGTCTGG"); 
	cout << "Hits = " << hits << endl;
	
	cout << "YPR073C, " << endl; 
	hits = my_dnaf.search_for_pcr_products("TTTCGTGTTAACTTTCCCTTCTCAGTTTTCTATCGCTTATCAAAAATCACAGGGTTTC", "TTGTTTCTCTACCAGATTTCTTCACTTTTGACAAACAACCCCTAGAATACAGACTATC");
	cout << "Hits = " << hits << endl;
	
	cout << "YGR125W, " << endl; 
	hits = my_dnaf.search_for_pcr_products("ATTTCTAAATTAAATTTTTGAGATTTTAAACTGTCTTGCTTATAGACC", "AATCTATATAGATATATGTATATAATAAAAGTTTTCAATCTGCTGACC");
	cout << "Hits = " << hits << endl;
	
	cout << "YNL169C, " << endl; 
	hits = my_dnaf.search_for_pcr_products("GCTGGCTTTGCTTTTCCTTCTTCTTCAAAAAATAACGACCAAG", "ACCACTAGACCTGCAAACTCGACTTTCAATGGTTTAAATGAAG"); 
	cout << "Hits = " << hits << endl;
	
	cout << "YJL134W, " << endl; 
	hits = my_dnaf.search_for_pcr_products("TTTTTATCTATACAAATGGCTTGTGTCCGCTGGTCTGGTAGG", "CACATTGTCGTCAAAAGGTAATCTCCCCAGAAGTTTTCACTC");
	cout << "Hits = " << hits << endl;
	
	cout << "YPL262W, " << endl; 
	hits = my_dnaf.search_for_pcr_products("GATTGTATCTCTTATCCGTTTAATAGTTAGACTTTATGG", "TTATCCTTCTTATTCTTCTTTGCCATTAATCACTTCTAC"); 
	cout << "Hits = " << hits << endl;
	/*
	cout << "YDR074W, " << endl; 
	hits = my_dnaf.search_for_pcr_products("TTCGGCACAGAAATAGTGACAGGCAGTGTTATTTTGG", "CTATCTTCTCTTGGCAAGAGTGTTCGAGATATGGTAG");
	cout << "Hits = " << hits << endl;
	
	cout << "YMR205C, " << endl; 
	hits = my_dnaf.search_for_pcr_products("TGCGTATGGTTAGTTCTTGGCCAATGCTAAACTAGTCTCTAAATC", "CCATCTCCACGCTAAAAGGAAAACGAAGATTAAGATAAAGTTGGG");
	cout << "Hits = " << hits << endl;
	*/
	
	stop = clock();
	
	cout << "Processing time: " << stop - start << " microseconds" << endl;
	//cout << "Hits = " << hits << endl;
	
	/*										//	Primer B - sense						Primer A - antisense 
	my_dnaf.GC_array_optimisation = FALSE;
	
	start = clock();
	hits = my_dnaf.search_for_pcr_products("CTTCGTACTAGTTTTAGCTTACCACCGCTTTTAGTG", "GCTTACATTTGCAGACGTCAAAAAACCTGACAAACC");
	//hits = my_dnaf.search_for_pcr_products("ATATGTGTTGAGTTTATGCTATGATTCTTTTTAAAGAATTCTATTATTAAAAGGATTAGG", "GTGATCATTTCTGTATAAGTTTTCTCAACTGAATATATAAATATATATCATTACTATGCT");
	stop = clock();
	
	cout << "Time: " << stop - start << endl;
	cout << "Hits = " << hits << endl;
	*/
    return(1);
}

