
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 primer_pair.h
 
 Created by:	Michael C. Riley
				Amanda Clare
 
 Date:			05/04/2011
 
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
 
 /** \mainpage PD5 software library for primer design applications
 
 \section intro_sec Introduction
 
 The PD5 software library is a versatile collection of software modules suitable for use in primer 
 design applications. It was developed in response to a growing need for increasingly complex primer
 design applications.
 
 The main class for designing primer pairs is primer_pair. The \ref primer class can be used
 where an individual primer needs to be assessed or designed. However, the analysis classes (\ref dimerisation, 
 annealing_temperature, DNAfind) are designed so that they can be used in isolation
 
 */

/** \file primer_pair.h
 \brief The main class for designing pairs of primers

 */


#ifndef PRIMER_PAIR_H
#define PRIMER_PAIR_H

#include "primer.h"
#include "primer_pair_data.h"
#include "dna_find.h"
#include "global_defs.h"




//! Primer pair design class
/**
 The main class for designing pairs of primers

 Having selected the location of the target and/or regions in which the primers may be located, 
 there are three steps :- 
 (1) Generate a number of candidates subject to hard constraints.
 (2) Analyse and score candidate primers
 (3) Sort according to score for most optimum pairs. 
 
 Hard constraints: primer location, primer length

 */

class primer_pair
{
public:
	primer_pair();
	~primer_pair(){};
	
	// Two instances of the primer class
	primer forward_primer;
	primer reverse_primer;
	
	// Array of optimum pair results
	primer_pair_data pair_candidate[38];
	
	
	// int get_primers(const char* dna_template);

	/** AFC: to check */
	int pair_dimerisation(void);

	/** Sets the begin and end locations of the region for amplification */
	int set_target_location(int begin, int end);
	/** Sets the allowed size of flanking regions that can be amplified either side of the target region */
	int set_flank_lengths(int upstream_flank_length, int downstream_flank_length);
	/** Sets the allowed length of the primers themselves */
	int set_primer_length_range(int shortest_length, int longest_length);
	
	/** Sets the maximum, minimum and optimum annealing temperatures for primers.
		Default 50, 55, 60.
	 */
	int set_Tm_range(double minimum, double optimum, double maximum);
	
	/** Generates up to 100 candidates passing hard constraints
	 */
	int generate_candidates(const char* template_sequence);
	
	/** For each candidate primer the scores for potential hairpin formation and 
	 dimerisation are calculated along with the annealing temperature. This must 
	 be called before sorting and selection methods
	 */
	int candidate_analysis(void);
	
	/** Sort the pair candidates by increasing Tm difference */
	int sort_Tm_difference(int data_size);
	/** Sort the pair candidates by increasing number of PCR products */
	int sort_products(int data_size);
	/** Sort the pair candidates by increasing pair dimer score of the forward primer */
	int sort_f_pair_dimer(int data_size);
	/** Sort the pair candidates by increasing pair dimer score of the reverse primer */
	int sort_r_pair_dimer(int data_size);
	
	/** Sort the individual candidates by increasing pair dimer score of the reverse primer 
	    (that is, sort the forward primer candidates, then sort the reverse primer candidates) 
	*/
	int sort_individual_candidates(const char* priority_list);
	/** Show all the forward primer candidates, then show all the reverse primer candidates */
	int show_individual_candidates(void);
	
	/** Fill a pair_candidates array with the pairs made by the cross product of the forward 
	    candidates and the reverse candidates. That is, each forward candidate is paired with 
	    every reverse candidate 
	*/
	int make_pair_candidates(int number_of_forward_candidates, int number_of_reverse_candidates);
	/** The total number of pairs of candidates that were generated (including ones that are not good) */
	int number_of_pair_candidates;
	/** The number of pairs of candidates being considered that are still "good" (that meet the 
	    criteria given) 
	*/
	int good_pair_candidates;
	/** Sort the primer pair candidates by the chosen priorities. If the pair 
	    candidates array was not already populated with pairs, get the best 6 
	    (or less if 6 not available) of each forward and reverse primers to 
	    make 36 primer pairs before sorting.
	*/
	int sort_pair_candidates(const char* priority_list);
	/** Write the details of the var best candidates to stdout */
	int show_best_pair_candidates(int var);
	
	//Make parameters public for now 6/4/11
	
	int forward_start_location_range_begin;
	int forward_start_location_range_end;
	int forward_length_range_shortest;
	int forward_length_range_longest;
	int forward_required_GC_content;
	int forward_GC_tolerance;
	
	int reverse_start_location_range_begin;
	int reverse_start_location_range_end;
	int reverse_length_range_shortest;
	int reverse_length_range_longest;
	int reverse_required_GC_content;
	int reverse_GC_tolerance;
	
	/** @name moopair 
	 * Primer pair multi objective optimisation 
	 */
	///@{
	/** Uses multi objective optimisation to sort for optimum primer pairs */
	int sort_moo(int data_size);
	
	//! Multi Objective Optimisation weighting
	double forward_pair_dimer_weighting;
	
	//! Multi Objective Optimisation weighting
	double reverse_pair_dimer_weighting;
	
	//! Multi Objective Optimisation weighting
	double Tm_diff_weighting;
	 
	///@}

	
private:
	int upstream_flank_length;
	int downstream_flank_length;
	int target_5_prime;
	int target_3_prime;
	
	bool pair_candidates_exist;
	
};

#endif 
