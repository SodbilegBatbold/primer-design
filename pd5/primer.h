/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 primer.h
 
 Created by:	Michael C. Riley and Amanda Clare
 
 Date:			11/02/2011
 
 Created by Michael Riley on 11/02/2011 - mhr"at"aber.ac.uk
 
 Copyright (c) 2010, 2011 Aberystwyth University. 
 All rights reserved.
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 
 *******************************************************************/

/** \file primer.h
 \brief Individual primer design
 
 This file contains the methods and attributes required for individual primer design
 */

#ifndef PRIMER_H
#define PRIMER_H

#define nucleotide_A 0x41
#define nucleotide_C 0x43
#define nucleotide_G 0x47
#define nucleotide_T 0x54
#define nucleotide_N 0x4E

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "annealing_temperature.h"
#include "dimerisation.h"
#include "constraints.h"

#define TRUE 1
#define FALSE 0

#define ERROR -1
/**
 * Priorities to be used when sorting candidate primers.
 */
enum Priority { SORT_END, 
				YEAST_MATCH, 
				SELF_DIMER, 
				HAIRPIN, 
				LENGTH, 
				F_DIMER, 
				R_DIMER, 
				MOO_SORT,
				PLASMID_MATCH, 
				PRODUCTS, 
				SEQSIM_MATCH, 
				TEMPERATURE,
				TM_DIFF};

#define ADENINE 65
#define THYMINE 84
#define GUANINE 71
#define CYTOSINE 67
#define URACIL 85
#define ANYNUCLEOTIDE 78

using namespace std;

//! Individual primer data class
/**
 * Class holding the data for a single primer. 
 */

class primer_data
{
public:
	primer_data()
	{
	// initialise vars
		sequence[0] = 0;
		location_5_prime_end = 0;
		primer_length = 0;
		sticky_tail = FALSE;
		tail_check = 0;
		hairpin = 0;
		self_dimer = 0;
		primer_dimer = 0;
		forward_dimer = 0;
		reverse_dimer = 0;
		yeast_matches = 0;
		plasmid_matches = 0;
		seqsim_matches = 0;
		moo = 0.0;
		products = 0;
		annealing_temperature = 0;
	}

	char sequence[128];       ///< Primer sequence. Limited to 128 chars. We do not want primers over 100 nt in length anyway 
	int location_5_prime_end; ///< Location of the 5' end. Locations begin at 0, not 1. 
	double hairpin;           ///< The score for hairpin potential. 
	double self_dimer;        ///< The score for self-dimer potential.  
	int yeast_matches;        ///< The number of hits to the yeast genome 
	int plasmid_matches;      ///< The number of hits to the plasmid 
	int seqsim_matches;       ///< The number of sequence similarity hits 
	double moo;				///< multi objective optimisation score
	double annealing_temperature; ///< The calculated annealing temperature 
	int primer_length;		///< Sequence length of primer
	
	// Pair specific attributes (Maintained here for backward compatibility, but now located in primer_pair_data)
	double primer_dimer;      ///< Highest binding score obtained by primer pair convolution (Most common method).
	double forward_dimer;     ///< The highest score for forward primer tail binding to the reverse primer. 
	double reverse_dimer;     ///< The highest score for reverse primer tail binding to the forward primer.  
	int products;             ///< Number of potential products made by the primers. If > 1, then we have secondary products. 
	
	// Deprecated
	bool sticky_tail;         ///< Whether we allow sticky tails or not. 
	int tail_check;           ///< ?? 

	
	// Methods
	int self_dimerisation(void);
	int get_GC_content(void)
	{
		int GC_number = 0;
		int length = strlen(sequence);
		
		if(sequence[0] == 0) 
			return(ERROR);
		else
		{	
			for(int i = 0; i < length; i++) 
				if(sequence[i] == CYTOSINE || sequence[i] == GUANINE) GC_number++;
			
			return(GC_number);
		}
	}
};

//! Individual primer design class
/**
 The main class for individual primer design
 */

class primer: public annealing_temperature, constraints
{
public:
	
	primer();
	~primer()
	{	
				
	}

// Candidate generation

	/**
	 * generate_candidates: Use this function for most purposes.
	 * @param template_seq the template sequence in which you want
	 * to seach for primers.
	 */
	int generate_candidates(const char* template_seq)
	{
		return(generate_candidate_primers(template_seq));
	}
					
	int generate_candidate_primers(const char* template_sequence);
	
	int homopolymeric_run_detection(const char* sequence);

	/**
	 *\defgroup CandidateDisplay Candidate display
	 * @{
	 */
	int show_candidate(int candidate_index);
	int show_candidate(int candidate_index, ofstream &fout);
	int show_all_candidates();
	int show_all_candidates(ofstream &fout);
	int show_all_single_candidates();
	/** @} */

	/**
	 *\defgroup CandidateSelection Candidate selection
	 * @{
	 */

	/** 
	 * Returns the index of the first candidate primer which is
	 * good enough. Returns -1 if no such primer exists. 
	 * \sa rank_selection */
	int auto_selection(void);

	/** 
	 * Reorders the candidate list into ranked order. Returns 1 if 
	 * some successful candidates, 0 if no successful candidates. 
	 * \sa auto_selection */
	int sort_candidates(void);
	int rank_selection(void)
	{
		return(sort_candidates());
	}
	
	

	// int manual_selection(int candidate_number);
	/** @} */

	int set_primer_location_range(int begin, int end);
	int set_primer_length_range(int begin, int end);

	/**
	 *\defgroup PrimerParams User defined primer parameters
	 * @{
	 */
	
	bool verbose;
	
	int start_location_range_begin; ///< Beginning of range for 5' end 
	int start_location_range_end;   ///< End of range for 5' end 
	int length_range_shortest;      ///< Shortest allowable length 
	int length_range_longest;       ///< Longest allowable length 
	int required_GC_content;        ///< Absolute number of GCs required (useful when searching for the second primer in a pair). GC content not checked if 0 
	int GC_tolerance;               ///< Number that GC content can vary by. For example, 1 => GC content can vary by +/- 1. Percentage is not reliable here, use integer. 
	
	int optimum_primer_length;      ///< Preferred length for a primer if possible 
	
	double optimum_Tm;  ///< Preferred temp for a primer if possible 
	double max_Tm;					///< Maximum temp for a primer
	double min_Tm;					///< Minimum temp for a primer
	
	int candidates_found;           ///< Number of candidates meeting basic criteria. 
	int good_candidates;            ///< Number of candidates left after scoring and eliminating 
	
	bool reverse_primer;            ///< Whether this is  the reverse primer */
	bool downstream_search;			///< Set TRUE to search downstream towards target
	bool GC_clamping;               ///< ?? Default TRUE
	bool no_sticky_tails;           ///< ?? Default FALSE
	bool no_G_primer;				///< Candidate primers should contain no Gs
	bool no_C_primer;				///< Candidate primers should contain no Cs 
	bool tail_complementarity_check; ///< Used in the test_candidates method to opt for a tail complentarity check
    
	int set_sequence_to_avoid(const char* subsequence)
	{
		strcpy(seq_to_avoid, subsequence);
		avoid_subsequence_check = TRUE;
		return(TRUE);
	}
	
	bool homopolymeric_run_check; ///< Returns true if it find homopolymeric runs >= homopolymeric_run_length_limit
	int set_homopolymeric_run_length_limit(int limit) 
	{
		homopolymeric_run_length_limit = limit; ///< Default is 5
		return(TRUE);
	}
	
	// Auto selection parameters
	int yeast_nsb_limit;
	/**
	 * Default sorting priorities (can have up to 32 including SORT_END). Fill this
	 * with Priority elements in the order required. Finish with
	 * SORT_END ??.
	 */
#define PRIORITY_LIST_LENGTH_MAX 32
	
	Priority priority[PRIORITY_LIST_LENGTH_MAX];
	int set_priorities(const char* priority_list);
	int max_number_candidates;

	/** @} */

	/**
	 * The array of candidate primers. Maximum 102 candidates allowed.
	 */
	primer_data candidate[102];

	/**
	 *\defgroup SequenceAnalysis Sequence Analysis
	 * @{
	 */
	
// Sequence analysis
	double hairpin(int candidate_number);
	double hairpin(int candidate_number, ofstream &fout);
	int tail_check(int candidate_number, const char* b_sequence, ofstream &fout);
	int sticky_tail_check(int candidate_number);
	
	int tail_complementarity_report(const char* sequence);
	int tail_complementarity(const char* sequence);
	int self_dimer(int candidate_number);
	int self_dimer(int candidate_number, ofstream &fout);
	/** Deprecated. Do not use (only works for short primers) */
	int primer_dimer(int candidate_number_a, const char* b_sequence, ofstream &fout);
	int primer_dimer_2(int candidate_number, const char* b_sequence);
	int primer_dimer_2(int candidate_number, const char* b_sequence, ofstream &fout);
	int blast_yeast(int candidate_number, const char* name);
	int blast_plasmid(int candidate_number, const char* name);
	int blast_yeast(int candidate_number, const char* name, ofstream &fout);
	int blast_plasmid(int candidate_number, const char* name, ofstream &fout);
	int blast_seq(int candidate_number, const char* sequence);
	int blast_db(int candidate_number, const char* db_name);

	int calculate_temperature(int candidate_number);

	double expectation;
	
// FASTA implementation 17/2/11 & 24/2/11
	int Smith_Waterman(int candidate_number, const char* library);
	int fasta3(int candidate_number, const char* library);
	
// Multi Objective Optimisation weightings
	double hairpin_weighting;
	double self_dimer_weighting;
	double forward_pair_dimer_weighting;
	double reverse_pair_dimer_weighting;
	double annealing_temperature_weighting;
	double primer_length_weighting;


	/** @} */

private:
// Candidate generation
	//char complement[22];
	
	int test_candidate(const char* sequence, int location_3_prime_end, const char* template_sequence);
	int sort_yeast_matches(int data_size);
	int sort_plasmid_matches(int data_size);
	int sort_self_dimer(int data_size);
	int sort_hairpin(int data_size);
	int sort_f_dimer(int data_size);
	int sort_r_dimer(int data_size);
	int sort_moo(int data_size);
	int sort_products(int data_size);
	int sort_length(int data_size);	
	int sort_seqsim_matches(int data_size);
	int sort_Tm_optimum(int data_size);
	int sort_Tm_ascending(int data_size);
	int sort_Tm_descending(int data_size);
	
	bool avoid_subsequence_check;		///< Set to TRUE to avoid including seq_to_avoid in primer sequence
	char* seq_to_avoid;             ///< A sequence that must not be present in the primer 

	int homopolymeric_run_length_limit; //< Default = 5
	

};
#endif
