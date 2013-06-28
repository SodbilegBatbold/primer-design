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

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "dimerisation.h"
#include "primer_data.h"
#include "global_defs.h"

/**
 * Priorities to be used when sorting candidate primers.
 */
enum Priority { SORT_END, 
				SELF_DIMER, 
				HAIRPIN, 
				LENGTH, 
				F_DIMER, 
				R_DIMER, 
				MOO_SORT,
				BINDING_A, 
				BINDING_B, 
				PRODUCTS, 
				SEQSIM_MATCH, 
				TEMPERATURE,
				TM_DIFF};



//! Individual primer design class
/**
 The main class for individual primer design
 */

class primer
{
public:
	
	primer();
	~primer()
	{	
				
	}



	/** \deprecated 
	 * Do not use this function, use primer::generate_candidate_primers 
	 */
	int generate_candidates(const char* template_seq)
	{
		return(generate_candidate_primers(template_seq));
	}
			
	/** Use this function for most primer generation purposes.
	 * @param template_seq the template sequence in which you want
	 * to seach for primers. 
	 */
	int generate_candidate_primers(const char* template_sequence);
	

	/** @name CandidateDisplay 
	 * Primer candidate display: Various methods for displaying the primers
	 */
	///@{
	int show_candidate(int candidate_index);
	int show_candidate(int candidate_index, ofstream &fout);
	int show_all_candidates();
	int show_all_candidates(ofstream &fout);
	int show_all_single_candidates();
	///@}

	/** @name CandidateSelection 
	 * Methods for candidate selection
	 */
	///@{

	/** 
	 * Returns the index of the first candidate primer which is
	 * good enough. Returns ERROR if no such primer exists. 
	 * \sa rank_selection */
	int auto_selection(void);

	/** 
	 * Reorders the candidate list into ranked order. Returns TRUE if 
	 * some successful candidates, FALSE if no successful candidates. 
	 * \sa auto_selection */
	int sort_candidates(void);

	/** Reorders the candidate list into ranked order. Returns the index of 
	 * the first candidate primer which is good enough (this will always be 0).
	 * Number of good candidates is stored in \ref good_candidates.
	 * Returns ERROR if no such primer exists. 
	 * \sa auto_selection, sort_candidates  
	 */
	int rank_selection(void)
	{
	  if (!sort_candidates()) {return ERROR;} 
	  else { return 0; } // if sort is successful then first good candidate is in index 0; 
	}
	
       
	///@}

	/** Set the range of locations for the primer's 5' end */
	int set_primer_location_range(int begin, int end);
	/** Set the range of lengths for the primer */
	int set_primer_length_range(int begin, int end);

	/** Defines a sequence that must not be present in the primer */
	int set_sequence_to_avoid(const char* subsequence)
	{
		strcpy(seq_to_avoid, subsequence);
		avoid_subsequence_check = TRUE;
		return(TRUE);
	}

	/** Default limit is 5 */
	int set_homopolymeric_run_length_limit(int limit) 
	{
		homopolymeric_run_length_limit = limit; 
		return(TRUE);
	}

	/** Sets the primer::priority */
	int set_priorities(const char* priority_list); 



	/** @name SequenceAnalysis 
	 * Sequence Analysis
	 */
	///@{

	/** Calculates hairpin, self dimer and Tm for each candidate */
	int analyse_all_candidates(void); 

	/** Calculates hairpin score */
	double hairpin(int candidate_number);
	// double hairpin(int candidate_number, ofstream &fout);

	/** Calls sequence_utils::tail_check */
	int tail_check(int candidate_number, const char* b_sequence, ofstream &fout);
	/** Calls sequence_utils::sticky_tail_check */
	int sticky_tail_check(int candidate_number);

	/** As for primer::tail_complementarity but outputs detail to stdout */
	int tail_complementarity_report(const char* sequence);
	/** Checks if the tail of the sequence will bind to its complement, returns TRUE or FALSE */
	int tail_complementarity(const char* sequence);

	/** Calculates primer_data::self_dimer score and returns TRUE/FALSE */
	int self_dimer(int candidate_number);
	/** Calculates primer_data::self_dimer score, outputs to fout and returns TRUE/FALSE */
	int self_dimer(int candidate_number, ofstream &fout);

	/** \deprecated Do not use (only works for short primers) */
	int primer_dimer(int candidate_number_a, const char* b_sequence, ofstream &fout);
	/** Calculates and saves primer dimer forward and reverse scores. Returns TRUE/FALSE */
	int primer_dimer_2(int candidate_number, const char* b_sequence);
	/** Calculates and saves primer dimer forward and reverse scores, and outputs to fout. Returns TRUE/FALSE */
	int primer_dimer_2(int candidate_number, const char* b_sequence, ofstream &fout);

	/** \Deprecated */
	int blast_yeast(int candidate_number, const char* name);
	/** \Deprecated */
	int blast_plasmid(int candidate_number, const char* name);
	/** \Deprecated */
	int blast_yeast(int candidate_number, const char* name, ofstream &fout);
	/** \Deprecated */
	int blast_plasmid(int candidate_number, const char* name, ofstream &fout);

	/** Use BLAST to find the number of hits of the candidate against the given sequence, with e-value at primer::expectation or better. Calls sequence_utils::blast_seq */
	int blast_seq(int candidate_number, const char* sequence);
	/** Use BLAST to find the number of hits of the candidate against the given database, with e-value at primer::expectation or better. Calls sequence_utils::blast_db */
	int blast_db(int candidate_number, const char* db_name);
	
// FASTA implementation 17/2/11 & 24/2/11
	/** See sequence_utils::Smith_Waterman */
	int Smith_Waterman(int candidate_number, const char* library);
	/** See  sequence_utils::fasta3 */
	int fasta3(int candidate_number, const char* library);

	double expectation;     ///< E-value for use in BLAST searches.

	/** Uses a subclass of annealing_temperature to calculate and save the temperature. */
	int calculate_temperature(int candidate_number);

	/** Checks for homopolymeric run in sequence. Returns TRUE if found one that is longer 
	    than allowed length, FALSE otherwise. */
	int homopolymeric_run_detection(const char* sequence);


	///@} 

	/**
	 * The array of candidate primers. Maximum 102 candidates allowed.
	 */
	primer_data candidate[102];

	int candidates_found;           ///< Number of candidates meeting basic criteria. 
	int good_candidates;            ///< Number of candidates left after scoring and eliminating 


	/** @name PrimerParams 
	 * User defined primer parameters
	 */
	///@{

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
	
	
	bool reverse_primer;            ///< Whether this is  the reverse primer 
	bool downstream_search;			///< Set TRUE to search downstream towards target
	bool GC_clamping;               ///< Whether to require primer to end in G or C (but not allow more than 2 GC in final 4 bases. Default TRUE
	//bool no_sticky_tails;           ///< Default FALSE. AFC: NOT USED/CHECKED, DOES NOT CHANGE ANYTHING.
	bool no_G_primer;				///< Candidate primers should contain no Gs
	bool no_C_primer;				///< Candidate primers should contain no Cs 
	bool tail_complementarity_check; ///< Used in the test_candidates method to opt for a tail complentarity check
    
	bool homopolymeric_run_check; ///< Set to true to find homopolymeric runs >= homopolymeric_run_length_limit
	bool homopolymeric_run_nr_tail_check; ///< True will check for poly runs on the template at the 3' tail location of the primer

	int yeast_nsb_limit; ///< AFC is this used? Can we generalise it from yeast? 
	int max_number_candidates;

#define PRIORITY_LIST_LENGTH_MAX 32
	Priority priority[PRIORITY_LIST_LENGTH_MAX]; ///< Default sorting priorities (can have up to 32 including SORT_END). Fill this with Priority elements in the order required. Finish with SORT_END.


	///@}

	/** @name MOOWeightings 
	 * Multi Objective Optimisation weightings 
	 */

	///@{

	double hairpin_weighting;
	double self_dimer_weighting;
	double forward_pair_dimer_weighting;
	double reverse_pair_dimer_weighting;
	double annealing_temperature_weighting;
	double primer_length_weighting;

	///@}







private:
// Candidate generation
	//char complement[22];
	
	int test_candidate(const char* sequence, int location_3_prime_end, const char* template_sequence);
	int sort_binding_A(int data_size);
	int sort_binding_B(int data_size);
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

	int homopolymeric_run_length_limit; ///< Default = 5
	

};
#endif
