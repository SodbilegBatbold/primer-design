
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

/** \file primer_pair.h
 \brief Main classes primer pair design
 Primer pair
 
 Having selected the location of the target and/or regions in which the primers maybe located:-
 
 Three steps:-
 Generating a number of candidates subject to hard constraints.
 Analyse and score candidate primers
 Sort according to score for most optimum pairs
 
 Hard constraints
 primer location
 length
 
 */

#ifndef PRIMER_PAIR_H
#define PRIMER_PAIR_H

#include "primer.h"
#include "DNAfind.h"

#define NO_TEMPLATE_ERROR -32

//! Primer pair data class
/**
 The main class containing primer pair data
 */
class primer_pair_data
{
public:
	char forward_sequence[128];
	char reverse_sequence[128]; 
	int forward_index;
	int reverse_index;
	int location_forward_5_prime_end;
	int location_reverse_5_prime_end;
	double forward_hairpin_score;
	double reverse_hairpin_score;
	double forward_self_dimer_score;
	double reverse_self_dimer_score;
	double forward_annealing_temperature;
	double reverse_annealing_temperature;
	double forward_pair_dimer_score; ///< forward 3' tail on reverse primer
	double reverse_pair_dimer_score; ///< reverse 3' tail on forward primer
	int number_of_pcr_products;
	double annealing_temperature_difference;
	double moo_score;
	
	// Methods
	int pair_dimerisation(void);
	int pcr_products(DNAfind &mynsb);
};

//! Primer pair design class
/**
 The main class for designing pairs of primers
 */

class primer_pair
{
public:
	primer_pair();
	~primer_pair(){};
	
	// Make two instances of the primer class
	primer forward;
	primer reverse;
	
	// Array of optimum pair results
	primer_pair_data pair_candidate[38];
	
	int get_primers(const char* dna_template);
	int pair_dimerisation(void);
	int set_target_location(int begin, int end);
	int set_flank_lengths(int upstream_flank_length, int downstream_flank_length);
	int set_primer_length_range(int shortest_length, int longest_length);
	
	/** Sets the maximum, minimum and optimum annealing temperatures for primers
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
	
	int sort_Tm_difference(int data_size);
	int sort_f_pair_dimer(int data_size);
	int sort_r_pair_dimer(int data_size);
	
	int sort_individual_candidates(const char* priority_list);
	int sort_pair_candidates(const char* priority_list);
	int show_individual_candidates(void);
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
	
	/** \defgroup moopair Primer pair multi objective optimisation 
	 @{ */
	/** Uses multi objective optimisation to sort for optimum primer pairs */
	int sort_moo(int data_size);
	
	//! Multi Objective Optimisation weighting
	double forward_pair_dimer_weighting;
	
	//! Multi Objective Optimisation weighting
	double reverse_pair_dimer_weighting;
	
	//! Multi Objective Optimisation weighting
	double Tm_diff_weighting;
	 
	 /** @} */
// Accessor methods
	
private:
	int upstream_flank_length;
	int downstream_flank_length;
	int target_5_prime;
	int target_3_prime;
	
};


#endif

/** \mainpage PD5 software library for primer design applications
 
 \section intro_sec Introduction
 
 The PD5 software library is a versatile collection of software modules suitable for use in primer 
 design applications. It was developed in response to a growing need for increasingly complex primer
 design applications.
 
 The main class for designing primer pairs is primer_pair. The primer class can be used
 where an individual primer needs to be assessed or designed. However, the analysis classes (dimerisation, 
 annealing_temperature, DNAfind) are designed so that they can be used in isolation
 
 */