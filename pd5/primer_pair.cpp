
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 primer_pair.cpp
 
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

/** \file primer_pair.cpp
 \brief Main classes primer pair design
 
 
 */

#include <cmath> // for fabs

#include "primer_pair.h"


primer_pair::primer_pair(void)
{
	reverse_primer.reverse_primer = TRUE;
	
	upstream_flank_length = 100;
	downstream_flank_length = 100;
	
	// Init target begin/end locations
	target_5_prime = 0;
	target_3_prime = 0;
	
	//moo_score = 0.0;
	
	// Default weightings

	forward_pair_dimer_weighting = 1.0;
	reverse_pair_dimer_weighting = 1.0;
	Tm_diff_weighting = 1.0;
	
	pair_candidates_exist = FALSE;
}

int primer_pair::set_target_location(int begin, int end)
{
	if(begin < end)
	{
		target_5_prime = begin;
		target_3_prime = end;
	}
	else 
	{
		target_5_prime = end;
		target_3_prime = begin;
	}
	
	// We actually set the regions where the 5' end of the primers can be located
	
	forward_primer.set_primer_location_range(target_5_prime - upstream_flank_length, target_5_prime);
	reverse_primer.set_primer_location_range(target_3_prime, target_3_prime + downstream_flank_length);
	
	return(TRUE);
}

int primer_pair::set_flank_lengths(int upstream_flank_length, int downstream_flank_length)
{
	this->upstream_flank_length = upstream_flank_length;
	this->downstream_flank_length = downstream_flank_length;
	
	// Now reset primer location parameters to suit if target location set
	
	if(target_5_prime && target_3_prime)
	{
		forward_primer.set_primer_location_range(target_5_prime - upstream_flank_length, target_5_prime);
		reverse_primer.set_primer_location_range(target_3_prime, target_3_prime + downstream_flank_length);
	}
	
	return(TRUE);
}

/* 
int primer_pair::get_primers(const char* dna_template)
{
	bool quiet = FALSE;
#ifdef TESTING	
	std::cout << "Getting primers\n";
#endif
	
	// Check for template
	if(!dna_template)
	{
		// No template problem to deal with or:-
		return(NO_TEMPLATE_ERROR);
	}
	
	// Check parameters - assume ok for now 6/4/11
	return(FALSE)
}
*/

int primer_pair::pair_dimerisation(void)
{
	dimerisation pair_dimer;
	
	// Uses primer_data - DOES IT? AFC: LOOKS LIKE IT SETS PARAMS IN THE LOCAL pair_dimer
	pair_dimer.pair_dimer(forward_primer.candidate[0].sequence, reverse_primer.candidate[0].sequence);
	
	// Uses primer_pair_data
	//pair_dimer.pair_dimer(pair_candidate[0].forward_sequence, pair_candidate[0].reverse_sequence);
	
	return(TRUE);
}

int primer_pair::set_primer_length_range(int shortest_length, int longest_length)
{
	// This sets both primers to have the same length range
	// Use set_primer_length_range() individually for differing ranges
	
	forward_primer.set_primer_length_range(shortest_length, longest_length);
	reverse_primer.set_primer_length_range(shortest_length, longest_length);
		
	return(TRUE);
}

int primer_pair::generate_candidates(const char* template_sequence)
{
	forward_primer.generate_candidate_primers(template_sequence);
	reverse_primer.generate_candidate_primers(template_sequence);
	
	return(TRUE);
}

int primer_pair::candidate_analysis(void)
{
	if(forward_primer.candidates_found)
	{
		for(int idx = 0; idx < forward_primer.candidates_found; idx++)
		{
			forward_primer.hairpin(idx);
			forward_primer.self_dimer(idx);
			//forward_primer.candidate[idx].annealing_temperature = forward_primer.primer3_Tm(forward_primer.candidate[idx].sequence);
			forward_primer.calculate_temperature(idx);
		}
	}
	else {
		return(FALSE);
	}
	
	if(reverse_primer.candidates_found)
	{
		for(int idx = 0; idx < reverse_primer.candidates_found; idx++)
		{
			reverse_primer.hairpin(idx);
			reverse_primer.self_dimer(idx);
			//reverse_primer.candidate[idx].annealing_temperature = reverse_primer.primer3_Tm(reverse_primer.candidate[idx].sequence);
			reverse_primer.calculate_temperature(idx);
		}
	}
	else {
		return(FALSE);
	}

	
	return(TRUE);
}


/*******************************************
 *
 *			SORTING
 *
 *******************************************/

double s_sigmoid(double var)
{
	// double gain = 1.0;
	//return(1.0/(1.0 + (exp(-(gain * var)))));
	return(1.0/(1.0 + (exp(-var))));
}

inline void swap(int *a, int *b)
{
	int temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

int primer_pair::sort_Tm_difference(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(pair_candidate[i].annealing_temperature_difference > pair_candidate[i + 1].annealing_temperature_difference)
			{
				swap(pair_candidate[i], pair_candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}		
	return(1);
}

int primer_pair::sort_products(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(pair_candidate[i].number_of_pcr_products > pair_candidate[i + 1].number_of_pcr_products)
			{
				swap(pair_candidate[i], pair_candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}		
	return(1);
}
int primer_pair::sort_f_pair_dimer(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(pair_candidate[i].forward_pair_dimer_score > pair_candidate[i + 1].forward_pair_dimer_score)
			{
				swap(pair_candidate[i], pair_candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}

int primer_pair::sort_r_pair_dimer(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(pair_candidate[i].reverse_pair_dimer_score > pair_candidate[i + 1].reverse_pair_dimer_score)
			{
				swap(pair_candidate[i], pair_candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}

int primer_pair::sort_moo(int data_size)
{
	double max_score = 20.0; // This is true for all dimerisation scores at present (8/3/11)
	double max_Tm_diff = 10.0; // May need some work
	
	for(int i = 0; i < data_size; i++)
	{
		// Objective function
		pair_candidate[i].moo_score =	\
		forward_pair_dimer_weighting * s_sigmoid((double)pair_candidate[i].forward_pair_dimer_score - max_score) + \
		reverse_pair_dimer_weighting * s_sigmoid((double)pair_candidate[i].reverse_pair_dimer_score - max_score) + \
		Tm_diff_weighting * s_sigmoid((double)pair_candidate[i].annealing_temperature_difference - max_Tm_diff);
	}
	
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(pair_candidate[i].moo_score > pair_candidate[i + 1].moo_score)
			{
				swap(pair_candidate[i], pair_candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(TRUE);
}



int primer_pair::sort_individual_candidates(const char* priority_list)
{	
	forward_primer.set_priorities(priority_list);
	forward_primer.sort_candidates();
	
	reverse_primer.set_priorities(priority_list);
	reverse_primer.sort_candidates();
	
	return(TRUE);
}

int primer_pair::make_pair_candidates(int number_of_forward_candidates, int number_of_reverse_candidates)
{
	number_of_pair_candidates = (number_of_forward_candidates * number_of_reverse_candidates);
	
	int pair_idx = 0;
	
	for(int fwd_index = 0; fwd_index < number_of_forward_candidates; fwd_index++)
	{
		for(int rev_index = 0; rev_index < number_of_reverse_candidates; rev_index++)
		{
			strcpy(pair_candidate[pair_idx].forward_sequence, forward_primer.candidate[fwd_index].sequence);
			strcpy(pair_candidate[pair_idx].reverse_sequence, reverse_primer.candidate[rev_index].sequence);
			pair_candidate[pair_idx].forward_index = fwd_index;
			pair_candidate[pair_idx].reverse_index = rev_index;
			pair_candidate[pair_idx].location_forward_5_prime_end = forward_primer.candidate[fwd_index].location_5_prime_end;
			pair_candidate[pair_idx].location_reverse_5_prime_end = reverse_primer.candidate[rev_index].location_5_prime_end;
			pair_candidate[pair_idx].forward_hairpin_score = forward_primer.candidate[fwd_index].hairpin;
			pair_candidate[pair_idx].reverse_hairpin_score = reverse_primer.candidate[rev_index].hairpin;
			pair_candidate[pair_idx].forward_self_dimer_score = forward_primer.candidate[fwd_index].self_dimer;
			pair_candidate[pair_idx].reverse_self_dimer_score = reverse_primer.candidate[rev_index].self_dimer;
			pair_candidate[pair_idx].forward_annealing_temperature = forward_primer.candidate[fwd_index].annealing_temperature;
			pair_candidate[pair_idx].reverse_annealing_temperature = reverse_primer.candidate[rev_index].annealing_temperature;
			
			pair_candidate[pair_idx].annealing_temperature_difference = 
			fabs(forward_primer.candidate[fwd_index].annealing_temperature - reverse_primer.candidate[rev_index].annealing_temperature);
			
			pair_idx++;
		}
	}
	
	pair_candidates_exist = TRUE;
	return(TRUE);
}

int primer_pair::sort_pair_candidates(const char* priority_list)
{
	// If not already done, get the best 6 (or less if 6 not available) 
	// of each forward and reverse primers and make 36 primer pairs for sorting
	
	int i;
	int number_of_forward_candidates;
	int number_of_reverse_candidates;
	
	if(!pair_candidates_exist)
	{
		forward_primer.candidates_found > 5? 
			number_of_forward_candidates = 6: 
			number_of_forward_candidates = forward_primer.candidates_found;
		
		reverse_primer.candidates_found > 5? 
			number_of_reverse_candidates = 6: 
			number_of_reverse_candidates = reverse_primer.candidates_found;
		
		number_of_pair_candidates = (number_of_forward_candidates * number_of_reverse_candidates);
		
		int pair_idx = 0;
		
		for(int fwd_index = 0; fwd_index < number_of_forward_candidates; fwd_index++)
		{
			for(int rev_index = 0; rev_index < number_of_reverse_candidates; rev_index++)
			{
				strcpy(pair_candidate[pair_idx].forward_sequence, forward_primer.candidate[fwd_index].sequence);
				strcpy(pair_candidate[pair_idx].reverse_sequence, reverse_primer.candidate[rev_index].sequence);
				pair_candidate[pair_idx].forward_index = fwd_index;
				pair_candidate[pair_idx].reverse_index = rev_index;
				pair_candidate[pair_idx].location_forward_5_prime_end = forward_primer.candidate[fwd_index].location_5_prime_end;
				pair_candidate[pair_idx].location_reverse_5_prime_end = reverse_primer.candidate[rev_index].location_5_prime_end;
				pair_candidate[pair_idx].forward_hairpin_score = forward_primer.candidate[fwd_index].hairpin;
				pair_candidate[pair_idx].reverse_hairpin_score = reverse_primer.candidate[rev_index].hairpin;
				pair_candidate[pair_idx].forward_self_dimer_score = forward_primer.candidate[fwd_index].self_dimer;
				pair_candidate[pair_idx].reverse_self_dimer_score = reverse_primer.candidate[rev_index].self_dimer;
				pair_candidate[pair_idx].forward_annealing_temperature = forward_primer.candidate[fwd_index].annealing_temperature;
				pair_candidate[pair_idx].reverse_annealing_temperature = reverse_primer.candidate[rev_index].annealing_temperature;
				
				pair_candidate[pair_idx].annealing_temperature_difference = 
				fabs(forward_primer.candidate[fwd_index].annealing_temperature - reverse_primer.candidate[rev_index].annealing_temperature);
				
				pair_idx++;
			}
		}
	}
	
	// Analyse pairs
	
	for(i = 0; i < number_of_pair_candidates; i++)
	{
		pair_candidate[i].pair_dimerisation();
		
		// Initialising 
		//pair_candidate[i].number_of_pcr_products = 0;
		//pair_candidate[i].number_of_pcr_products = search_for_pcr_products(pcr1.pair_candidate[i].forward_sequence, pcr1.pair_candidate[i].reverse_sequence);
		
	}
	
	// Organise sort priorities
	
	//int priority_index = 0;
	int data_size = number_of_pair_candidates;
	
	char buffer[1024]; 
	char *token;
	
	strcpy(buffer, priority_list);
	token = strtok(buffer, ", \t\n\r");
	
	while(token)
	{
		if(!strcmp(token, "TM_DIFF"))
		{
			sort_Tm_difference(data_size);
		}
		else if(!strcmp(token, "F_DIMER"))
		{
			sort_f_pair_dimer(data_size);
			
			for(i = 0; i < data_size; i++)
				if(pair_candidate[i].forward_pair_dimer_score >= 12)data_size = i;
		}
		else if(!strcmp(token, "R_DIMER"))
		{
			sort_r_pair_dimer(data_size);
			
			for(i = 0; i < data_size; i++)
				if(pair_candidate[i].reverse_pair_dimer_score >= 12)data_size = i;
		}
		else if(!strcmp(token, "PRODUCTS"))
		{
			sort_products(data_size);
			for(i = 0; i < data_size; i++)
				if(pair_candidate[i].number_of_pcr_products > 1)data_size = i;
		}
		else if(!strcmp(token, "MOO_SORT"))
		{
			sort_moo(data_size);
		}
		else 
		{
			cout << "Priority list error\n";
		}
		
		token = strtok(NULL, ", \t\n\r");
		
	}
	good_pair_candidates = data_size;
	
	if(data_size < 1) return(FAIL);
	else return(PASS);
}

int primer_pair::show_individual_candidates(void)
{
	forward_primer.show_all_single_candidates();
	reverse_primer.show_all_single_candidates();
	
	return(TRUE);
}


int primer_pair::show_best_pair_candidates(int var)
{
	std::cout.setf(ios::fixed);
	std::cout.setf(ios::showpoint);
	
	for(int x = 0; x < var; x++)
	{
		std::cout.precision(1);
		std::cout << "Candidate pair " << x + 1 << std::endl;
		std::cout << "No.\tPrimer\tHpin\tSdimer\tTm\tLocus\tLength\tSequence\n";
		
		std::cout << pair_candidate[x].forward_index << "\t";
		std::cout << "Forward\t";
		std::cout << pair_candidate[x].forward_hairpin_score << "\t";
		std::cout << pair_candidate[x].forward_self_dimer_score << "\t";
		std::cout << pair_candidate[x].forward_annealing_temperature << "\t";
		std::cout << pair_candidate[x].location_forward_5_prime_end << "\t";
		std::cout << strlen(pair_candidate[x].forward_sequence) << "\t";
		std::cout << pair_candidate[x].forward_sequence << "\t";
		std::cout << std::endl;
		
		std::cout << pair_candidate[x].reverse_index << "\t";
		std::cout << "Reverse\t";
		std::cout << pair_candidate[x].reverse_hairpin_score << "\t";
		std::cout << pair_candidate[x].reverse_self_dimer_score << "\t";
		std::cout << pair_candidate[x].reverse_annealing_temperature << "\t";
		std::cout << pair_candidate[x].location_reverse_5_prime_end << "\t";
		std::cout << strlen(pair_candidate[x].reverse_sequence) << "\t";
		std::cout << pair_candidate[x].reverse_sequence << "\t";
		std::cout << std::endl << std::endl;
		
		std::cout << "Products ";
		if(pair_candidate[x].number_of_pcr_products)std::cout << pair_candidate[x].number_of_pcr_products;
		else std::cout << "n/a";
		std::cout << "\t";
		std::cout.precision(2);
		std::cout << "Tm diff. " << pair_candidate[x].annealing_temperature_difference << "\t";
		std::cout.precision(1);
		std::cout << "Fwd tail dimer " << pair_candidate[x].forward_pair_dimer_score << "\t";
		std::cout << "Rev tail dimer " << pair_candidate[x].reverse_pair_dimer_score << "\t";
		std::cout << std::endl << std::endl << "-----------------------------------" << std::endl << std::endl;
	}
	
	return(TRUE);
}



// Sets the required minimum, optimum and maximum annealing temperatures for both forward and reverse primers.

int primer_pair::set_Tm_range(double minimum, double optimum, double maximum)
{
	if(minimum < 0.0)
		return(FALSE);
	else
	{
		forward_primer.min_Tm = minimum;
		reverse_primer.min_Tm = minimum;
	}
	
	if(optimum > maximum || optimum < minimum)
		return(FALSE);
	else
	{
		forward_primer.optimum_Tm = optimum;
		reverse_primer.optimum_Tm = optimum;
	}
	
	if(maximum > 100.0)
		return(FALSE);
	else
	{
		forward_primer.max_Tm = maximum;
		reverse_primer.max_Tm = maximum;
	}
	
	return(TRUE);
}






