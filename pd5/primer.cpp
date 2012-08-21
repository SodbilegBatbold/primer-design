
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 primer.cpp
 
 Created by:	Michael C. Riley and Amanda Clare
 
 Date:			11/02/2011
 
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

/** \file primer.cpp
 \brief Individual primer design
 
 This file contains the methods and attributes required for individual primer design
 */

#include <cmath> // for fabs

#include "primer.h"

#define TRUE 1
#define FALSE 0

#define LOCATION_ERROR 0
#define SET_LOCATION_ERROR 0

//#define DEBUG

primer::primer()
{
	reverse_primer = FALSE;
	verbose = FALSE;
	// Initialise this look-up string for use by base_complement
	//strcpy(complement, "T G   C      N     A");
		
	expectation = 0.5;
	
	// DEFAULT VALUES
	max_number_candidates = 100;
	downstream_search = FALSE;
	GC_clamping = TRUE;
	no_sticky_tails = FALSE;
	no_G_primer = FALSE;
	no_C_primer = FALSE;
	tail_complementarity_check = TRUE;
	required_GC_content = 0; // GC content not checked if 0
	GC_tolerance = 0; // ie 1 => GC content can vary by +/- 1. Percentage is not reliable here, use integer.
	
	length_range_shortest = 18;      
	length_range_longest = 25; 
	
	yeast_nsb_limit = 1; //
	
	optimum_primer_length = 21;
	optimum_Tm = 55.0;
	max_Tm = 60.0;
	min_Tm = 50.0;
	
	avoid_subsequence_check = FALSE;
	seq_to_avoid = NULL;
	
	homopolymeric_run_check = FALSE;
	homopolymeric_run_nr_tail_check = TRUE;
	homopolymeric_run_length_limit = 5; //< popularly < 4-5 nt (Buck et al, BioTechniques 27:528-536 (September 1999))

	
	// Default sorting priorities (can have up to 10)
	priority[0] = BINDING_B;
	priority[1] = BINDING_A;
	priority[2] = SELF_DIMER;
	priority[3] = SORT_END;
	priority[4] = SORT_END;
	
	// Default weightings
	hairpin_weighting = 1.0;
	self_dimer_weighting = 1.0;
	forward_pair_dimer_weighting = 1.0;
	reverse_pair_dimer_weighting = 1.0;
	annealing_temperature_weighting = 1.0;
	primer_length_weighting = 2.0;
	
	good_candidates = 0;
	
}


int primer::set_primer_location_range(int begin, int end)
{
	if(begin == 0) begin = 1;  // Assume user enters 0 means the beginning of the sequence
	
	if(begin < 1 || end < 1)return(SET_LOCATION_ERROR);

	if(reverse_primer) // begin > end
	{
		if(end < begin)
		{
			start_location_range_begin = begin;
			start_location_range_end = end;
		}
		else 
		{
			start_location_range_begin = end;
			start_location_range_end = begin;
		}
	}
	else // begin < end		
	{
		if(end >= begin)
		{
			start_location_range_begin = begin;
			start_location_range_end = end;
		}
		else 
		{
			start_location_range_begin = end;
			start_location_range_end = begin;
		}
	}
	return(TRUE);
}

int primer::set_primer_length_range(int begin, int end)
{
	if(end >= begin)
	{
		length_range_shortest = begin;
		length_range_longest = end;
	}
	else 
	{
		length_range_shortest = end;
		length_range_longest = begin;
	}
	
	return(TRUE);
}


/** CANDIDATE GENERATION */

int primer::generate_candidate_primers(const char* template_sequence)
{
	int i, j, k;
	int n = 0;
	char possible_candidate[length_range_longest + 1];
	
	// Check for invalid primer location region
	if((unsigned int)start_location_range_begin > strlen(template_sequence))return(LOCATION_ERROR);
	if((unsigned int)(start_location_range_end + length_range_longest) > strlen(template_sequence))return(LOCATION_ERROR);
	
	// Start locations are physical locations starting at 1, not 0, so adjust for array indexing
	start_location_range_begin--;
	start_location_range_end--;
	   
	// Find possible primers and test
	if(downstream_search) // default search is in the upstream direction
	{
		if(reverse_primer)
		{
	#ifdef DEBUG
			std::cout << "Reverse downstream\n";
	#endif
			
			for(i = start_location_range_begin; i >= start_location_range_end; i--)
			{
				for(j = length_range_shortest; j <= length_range_longest; j++)
				{
					for(k = 0; k < j; k++)
					  possible_candidate[k] = sequence_utils::nucleotide_complement(template_sequence[i - k]);
					
					possible_candidate[k] = 0;  // end the string
					
					// CANDIDATE TESTING
					if(test_candidate(possible_candidate, i + strlen(possible_candidate), template_sequence))
					{
						candidate[n].location_5_prime_end = i + 1; // + 1 for physical location  
						candidate[n].primer_length = strlen(possible_candidate);
						strcpy(candidate[n++].sequence, possible_candidate);
					}
					if(n >= max_number_candidates)break;
				}
				if(n >= max_number_candidates)break;
				
			}
		}
		else
		{
	#ifdef DEBUG
			std::cout << " Forward downstream\n";
	#endif
			
			for(i = start_location_range_begin; i <= start_location_range_end; i++)
			{
				for(j = length_range_shortest; j <= length_range_longest; j++)
				{
					for(k = 0; k < j; k++)
						possible_candidate[k] = template_sequence[i + k]; 
					
					possible_candidate[k] = 0;  // end the string					
#ifdef DEBUG
					std::cout << possible_candidate << endl;
#endif
					
					// CANDIDATE TESTING
					if(test_candidate(possible_candidate, i + strlen(possible_candidate), template_sequence))
					{
						candidate[n].location_5_prime_end = i + 1; // + 1 for physical location  
						candidate[n].primer_length = strlen(possible_candidate);
						strcpy(candidate[n++].sequence, possible_candidate);
					}
					if(n >= max_number_candidates)break;					
				}
				if(n >= max_number_candidates)break;
			}
		}
	}
	else // upstream search
	{
		if(reverse_primer)
		{
#ifdef DEBUG
			std::cout << "Reverse upstream\n" << start_location_range_end << " and " << start_location_range_begin << endl;
#endif
			
			for(i = start_location_range_end; i <= start_location_range_begin; i++)
			{
				for(j = length_range_shortest; j <= length_range_longest; j++)
				{
					for(k = 0; k < j; k++)
						possible_candidate[k] = sequence_utils::nucleotide_complement(template_sequence[i - k]);
					
					possible_candidate[k] = 0;  // end the string
#ifdef DEBUG
					std::cout << possible_candidate << endl;
#endif
					
					// CANDIDATE TESTING
					if(test_candidate(possible_candidate, i + strlen(possible_candidate), template_sequence))
					{
						candidate[n].location_5_prime_end = i + 1; // + 1 for physical location 
						candidate[n].primer_length = strlen(possible_candidate);
						strcpy(candidate[n++].sequence, possible_candidate);
					}
					if(n >= max_number_candidates)break;
				}
				if(n >= max_number_candidates)break;
				
			}
		}
		else
		{
#ifdef DEBUG
			std::cout << " Forward upstream\n";
#endif
			
			for(i = start_location_range_end; i >= start_location_range_begin; i--)
			{
				for(j = length_range_shortest; j <= length_range_longest; j++)
				{
					for(k = 0; k < j; k++)
						possible_candidate[k] = template_sequence[i + k];
					
					possible_candidate[k] = 0;  // end the string
#ifdef DEBUG
					std::cout << possible_candidate << endl;
#endif
					
					// CANDIDATE TESTING
					if(test_candidate(possible_candidate, i + strlen(possible_candidate), template_sequence))
					{
						candidate[n].location_5_prime_end = i + 1; // + 1 for physical location
						candidate[n].primer_length = strlen(possible_candidate);
						strcpy(candidate[n++].sequence, possible_candidate);
					}
					if(n >= max_number_candidates)break;
					
				}
				if(n >= max_number_candidates)break;
			}
		}
	}
	
	
	candidates_found = n;
	
	if(n == 0)
		return(0); // i.e. no candidates found
	else 
		return(1);
}


int primer::test_candidate(const char* sequence, int location_3_prime_end, const char* template_sequence)
{
	bool Good_primer = TRUE;
	int tail_end = strlen(sequence) - 1;
	int i = 0;
	char poly_sequence[32]; // sequence on template about the 3' end of each primer. Used for homopolymeric run detection
	int poly_sequence_req_length = 11;
	
	//std::cout << "End is " << sequence[strlen(sequence) - 1] << std::endl;
	if(GC_clamping)
	{		
		if(!(sequence[tail_end] == nucleotide_G || sequence[tail_end] == nucleotide_C))
			Good_primer = FALSE;
		
		int tail_GC_content = 0;
		
		for(i = tail_end - 4; i < tail_end; i++)
		{
			if(sequence[i] == nucleotide_G || sequence[i] == nucleotide_C)
				tail_GC_content++;
		}
		
		if(tail_GC_content > 2)
			Good_primer = FALSE;
	}
	
	if(tail_complementarity_check)
	{
		if(tail_complementarity(sequence))
			Good_primer = FALSE;
	}
	
	if(required_GC_content)
	{
	  if(sequence_utils::GC_content(sequence) < required_GC_content - GC_tolerance || sequence_utils::GC_content(sequence) > required_GC_content + GC_tolerance) 
			Good_primer = FALSE;
	}
	
	if(no_G_primer)
	{
	  if(sequence_utils::nucleotide_content(GUANINE, sequence) > 0)
			Good_primer = FALSE;
	}

	if(no_C_primer)
	{
		if(sequence_utils::nucleotide_content(CYTOSINE, sequence) > 0)
			Good_primer = FALSE;
	}

	if(sequence_utils::nucleotide_content(ANYNUCLEOTIDE, sequence) > 0)
	  Good_primer = FALSE;
	
	if(avoid_subsequence_check) 
	{
	  if(strstr(sequence, seq_to_avoid)) 
	    Good_primer = FALSE;
	}
	
/** Checks for homopolymeric runs within the primer sequence. If they exist it increases the likelihood
 of secondary binding/products.
 */	
	if(homopolymeric_run_check)
	{
		if(homopolymeric_run_detection(sequence))
			Good_primer = FALSE;
	}

/** Checks for homopolymeric runs on the template at the primer sequence 3' end. If they exist it increases the likelihood
	 of mispriming. This should not be optional.
	 */	
	if(homopolymeric_run_nr_tail_check)
	{		
		for(i = 0; i < poly_sequence_req_length; i++) poly_sequence[i] = template_sequence[location_3_prime_end + i - 6];
		poly_sequence[poly_sequence_req_length] = 0;
		
		if(homopolymeric_run_detection(poly_sequence))
			Good_primer = FALSE;
	}
	

	return(Good_primer);
	
}

int primer::homopolymeric_run_detection(const char* sequence)
{
	if(sequence == NULL)return(FALSE);
	
	int sequence_length = strlen(sequence);
	int i, homopolymeric_length = 0;
	
	for(i = 1; i < sequence_length; i++)
	{
		if(sequence[i - 1] == sequence[i])
			homopolymeric_length++;
		else
			homopolymeric_length = 0;
		
		if(homopolymeric_length >= homopolymeric_run_length_limit)
			return(TRUE);
	}
	return(FALSE);
}

int primer::analyse_all_candidates(void)
{
	
	for(int i = 0; i < candidates_found; i++)
	{			
		hairpin(i);
		self_dimer(i);
		calculate_temperature(i);
	}
	
	return(0);
}

int primer::show_candidate(int i)
{
	std::cout << "Sequence, Hpin, Sdimer, F dimer, R dimer, Yeast NSB, URA3 NSB, Location\n";
	
	std::cout << i << ": " << candidate[i].sequence << ", ";
	//std::cout << candidate[i].sticky_tail << ", ";
	//std::cout << candidate[i].tail_check << ", ";
	std::cout << candidate[i].hairpin << ", ";
	std::cout << candidate[i].self_dimer << ", ";
	std::cout << candidate[i].forward_dimer << ", ";
	std::cout << candidate[i].reverse_dimer << ", ";
	std::cout << candidate[i].binding_A << ", ";
	std::cout << candidate[i].binding_B << ", ";
	std::cout << candidate[i].location_5_prime_end << ", ";
	std::cout << std::endl;
	
	
	return(0);
}

int primer::show_candidate(int i, ofstream &fout)
{
	//fout << "Sequence, Hpin, Sdimer, F dimer, R dimer, Yeast NSB, URA3 NSB\n";
	
	fout << i << ": " << candidate[i].sequence << ", ";
	fout << candidate[i].hairpin << ", ";
	fout << candidate[i].self_dimer << ", ";
	fout << candidate[i].forward_dimer << ", ";
	fout << candidate[i].reverse_dimer << ", ";
	fout << candidate[i].binding_A << ", ";
	fout << candidate[i].binding_B << ", ";
	fout << candidate[i].location_5_prime_end << ", ";
	fout << std::endl;
	
	return(0);
}

int primer::show_all_candidates(void)
{
	std::cout << "Sequence, Hpin, Sdimer, F dimer, R dimer, Yeast NSB, URA3 NSB, Location, Length, Products, Tm\n";
	
	for(int i = 0; i < candidates_found; i++)
	{			
		std::cout << i << ": " << candidate[i].sequence << ", ";
		//std::cout << candidate[i].sticky_tail << ", ";
		//std::cout << candidate[i].tail_check << ", ";
		std::cout << candidate[i].hairpin << ", ";
		std::cout << candidate[i].self_dimer << ", ";
		std::cout << candidate[i].forward_dimer << ", ";
		std::cout << candidate[i].reverse_dimer << ", ";
		std::cout << candidate[i].binding_A << ", ";
		std::cout << candidate[i].binding_B << ", ";
		std::cout << candidate[i].location_5_prime_end << ", ";
		std::cout << strlen(candidate[i].sequence) << ", ";
		std::cout << candidate[i].products << ", ";
		std::cout << candidate[i].annealing_temperature << ", ";
		std::cout << std::endl;
	}
	
	return(0);
}

int primer::show_all_candidates(ofstream &fout)
{
	fout << "Sequence, Hpin, Sdimer, F dimer, R dimer, Yeast NSB, URA3 NSB, Location\n";
	
	for(int i = 0; i < candidates_found; i++)
	{			
		fout << i << ": " << candidate[i].sequence << ", ";
		fout << candidate[i].hairpin << ", ";
		fout << candidate[i].self_dimer << ", ";
		fout << candidate[i].forward_dimer << ", ";
		fout << candidate[i].reverse_dimer << ", ";
		fout << candidate[i].binding_A << ", ";
		fout << candidate[i].binding_B << ", ";
		fout << candidate[i].location_5_prime_end << ", ";
		fout << std::endl;
		
	}
	
	return(0);
}

int primer::show_all_single_candidates(void)
{
	std::cout << "No.\tHpin\tSdimer\tNSB1\tNSB2\tLoc 5'\tLen\tTm\tSequence\n";
	
	for(int i = 0; i < candidates_found; i++)
	{	
		//std::cout.setf(ios::fixed);
		//std::cout.setf(ios::showpoint);
		std::cout << fixed << showpoint;
		std::cout.precision(1);
		std::cout << i << "\t";
		std::cout << candidate[i].hairpin << "\t";
		std::cout << candidate[i].self_dimer << "\t";
		std::cout << candidate[i].binding_A << "\t";
		std::cout << candidate[i].binding_B << "\t";
		std::cout << candidate[i].location_5_prime_end << "\t";
		std::cout << strlen(candidate[i].sequence) << "\t";
		std::cout << candidate[i].annealing_temperature << "\t";
		std::cout << candidate[i].sequence << "\t";
		std::cout << std::endl;
		
		// std::cout.precision(6); //default
	}
	
	return(0);
}

int primer::auto_selection(void)
{
	int selection;
	bool FOUND = FALSE;
	
	for(int i = 0; i < candidates_found; i++)
	{
		
		if((candidate[i].hairpin       < 11) &&
		   (candidate[i].self_dimer    < 11) &&
		   (candidate[i].forward_dimer    < 11) &&
		   (candidate[i].reverse_dimer    < 11) &&
		   (candidate[i].binding_A <= yeast_nsb_limit))
		{
#ifdef DEBUG
			std::cout << i << ": " << candidate[i].sequence << ", ";
			//std::cout << candidate[i].sticky_tail << ", ";
			//std::cout << candidate[i].tail_check << ", ";
			std::cout << candidate[i].hairpin << ", ";
			std::cout << candidate[i].self_dimer << ", ";
			std::cout << candidate[i].forward_dimer << ", ";
			std::cout << candidate[i].reverse_dimer << ", ";
			std::cout << candidate[i].binding_A << ", ";
			std::cout << candidate[i].binding_B << ", ";
			std::cout << std::endl;
#endif
			
			selection = i;
			
			FOUND = TRUE;
		}
		
		if(FOUND)break;
	}
	if(FOUND)
		return(selection);
	else
		return(-1);
	
}

int primer::set_priorities(const char* priority_list)
{
	char buffer[1024]; 
	char *token;
	 int priority_index = 0;
	 
	strcpy(buffer, priority_list);
	 token = strtok(buffer, ", \t\n\r");
	 
	 while(token)
	 {
		 if(!strcmp(token, "HAIRPIN"))
		 {
			 priority[priority_index] = HAIRPIN;
		 }
		 else if(!strcmp(token, "SELF_DIMER"))
		 {
			 priority[priority_index] = SELF_DIMER;
		 }
		 else if(!strcmp(token, "F_DIMER"))
		 {
			 priority[priority_index] = F_DIMER;
		 }
		 else if(!strcmp(token, "R_DIMER"))
		 {
			 priority[priority_index] = R_DIMER;
		 }
		 else if(!strcmp(token, "LENGTH"))
		 {
			 priority[priority_index] = LENGTH;
		 }
		 else if(!strcmp(token, "PRODUCTS"))
		 {
			 priority[priority_index] = PRODUCTS;
		 }
		 else if(!strcmp(token, "MOO_SORT"))
		 {
			 priority[priority_index] = MOO_SORT;
		 }
		 else if(!strcmp(token, "TEMPERATURE"))
		 {
			 priority[priority_index] = TEMPERATURE;
		 }
		 else if(!strcmp(token, "BINDING_A"))
		 {
			 priority[priority_index] = BINDING_A;
		 }
		 else if(!strcmp(token, "BINDING_B"))
		 {
			 priority[priority_index] = BINDING_B;
		 }
		 else 
		 {
			 cout << "Priority list error\n";
		 }
		 
		 token = strtok(NULL, ", \t\n\r");
	 
		 priority_index++;
	 }
	 
	 priority[priority_index] = SORT_END;
	 
	
	return(TRUE);
}

double sigmoid(double var)
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


int primer::sort_binding_A(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].binding_A > candidate[i + 1].binding_A)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}		
	return(1);
}

int primer::sort_binding_B(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].binding_B > candidate[i + 1].binding_B)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}		
	return(1);
}

int primer::sort_seqsim_matches(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].seqsim_matches > candidate[i + 1].seqsim_matches)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}		
	return(1);
}


int primer::sort_self_dimer(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].self_dimer > candidate[i + 1].self_dimer)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}


int primer::sort_hairpin(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].hairpin > candidate[i + 1].hairpin)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}


int primer::sort_f_dimer(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].forward_dimer > candidate[i + 1].forward_dimer)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}

int primer::sort_r_dimer(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].reverse_dimer > candidate[i + 1].reverse_dimer)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}

int primer::sort_moo(int data_size)
{
	double max_score = 20.0; // This is true for all dimerisation scores at present (8/3/11)
	
	for(int i = 0; i < data_size; i++)
	{
		// Objective function
		candidate[i].moo =	hairpin_weighting * sigmoid((double)candidate[i].hairpin - max_score) + \
		self_dimer_weighting * sigmoid((double)candidate[i].self_dimer - max_score) + \
		/*forward_pair_dimer_weighting * sigmoid((double)candidate[i].forward_dimer - max_score) + \
		reverse_pair_dimer_weighting * sigmoid((double)candidate[i].reverse_dimer - max_score) + \*/
		annealing_temperature_weighting * sigmoid(candidate[i].annealing_temperature - 40) + \
		primer_length_weighting * sigmoid(candidate[i].primer_length - 13);
	}
	
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].moo > candidate[i + 1].moo)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(TRUE);
}


int primer::sort_products(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].products > candidate[i + 1].products)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}

int primer::sort_length(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(abs((int)strlen(candidate[i].sequence) - optimum_primer_length) > abs((int)strlen(candidate[i + 1].sequence) - optimum_primer_length))
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}

int primer::sort_Tm_ascending(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].annealing_temperature > candidate[i + 1].annealing_temperature)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}

int primer::sort_Tm_descending(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(candidate[i].annealing_temperature < candidate[i + 1].annealing_temperature)
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}


int primer::sort_Tm_optimum(int data_size)
{
	bool swapped = TRUE;
	
	while(swapped)
	{
		swapped = FALSE;
		
		for(int i = 0; i < data_size - 1; i++)
		{
			if(fabs(candidate[i].annealing_temperature - optimum_Tm) > fabs(candidate[i + 1].annealing_temperature - optimum_Tm))
			{
				swap(candidate[i], candidate[i + 1]);
				swapped = TRUE;
			}
		}
	}
	return(1);
}


int primer::sort_candidates(void)
{
	int i, j;
	//bool swapped = TRUE;
	int data_size = candidates_found;

	for(j = 0; j < PRIORITY_LIST_LENGTH_MAX; j++)
	{
		if(priority[j])
		{
			switch(priority[j])
			{
				case BINDING_A:	
					sort_binding_A(data_size);									
					for(i = 0; i < data_size; i++)
						if(candidate[i].binding_A >= 4)data_size = i;
					break;
					
				case BINDING_B:	
					sort_binding_B(data_size);									
					for(i = 0; i < data_size; i++)
						if(candidate[i].binding_B >= 1)data_size = i;
					break;
					
				case SELF_DIMER:	
					sort_self_dimer(data_size);
					for(i = 0; i < data_size; i++)
						if(candidate[i].self_dimer >= 12)data_size = i;
#ifdef DEBUG
					cout << "Good candidates after Self dimer = " << data_size << endl;
#endif
					break;
					
				case HAIRPIN:	
					sort_hairpin(data_size);
					for(i = 0; i < data_size; i++)
						if(candidate[i].hairpin >= 12)data_size = i;
#ifdef DEBUG
					cout << "Good candidates after Hairpin = " << data_size << endl;
#endif
					break;
					
				case LENGTH:	// Primer length	
					sort_length(data_size);
					break;

				case TEMPERATURE:
					sort_Tm_ascending(data_size);
					for(i = 0; i < data_size; i++)
						if(candidate[i].annealing_temperature >= max_Tm)data_size = i;
					sort_Tm_descending(data_size);
					for(i = 0; i < data_size; i++)
						if(candidate[i].annealing_temperature <= min_Tm)data_size = i;
					sort_Tm_optimum(data_size);
#ifdef DEBUG
					cout << "Good candidates after Temperature = " << data_size << endl;
#endif
					break;
					
				case F_DIMER:	
					sort_f_dimer(data_size);
					for(i = 0; i < data_size; i++)
						if(candidate[i].forward_dimer >= 12)data_size = i;
					break;
					
				case R_DIMER:	
					sort_r_dimer(data_size);
					for(i = 0; i < data_size; i++)
						if(candidate[i].reverse_dimer >= 12)data_size = i;
					break;
					
				case PRODUCTS:	
					sort_products(data_size);
					for(i = 0; i < data_size; i++)
						if(candidate[i].products >= 2)data_size = i;
					break;
					
				case SEQSIM_MATCH:	
					sort_seqsim_matches(data_size);									
					for(i = 0; i < data_size; i++)
						if(candidate[i].seqsim_matches > 1)data_size = i; 
					break;	
					
				case MOO_SORT:
					sort_moo(data_size);
					break;
					
				case TM_DIFF:
					// Not valid for single primer
					break;
					
					
				case SORT_END:
				        j = PRIORITY_LIST_LENGTH_MAX;
					break;			
			}
		}
		else break;
	}
	
	good_candidates = data_size;
	
	if(data_size < 1) return(0);
	else return(1);
}





// CANDIDATE ANALYSIS

double primer::hairpin(int candidate_number) 
{	
	dimerisation hairpin;
	
	hairpin.hairpin(candidate[candidate_number].sequence);
	candidate[candidate_number].hairpin = hairpin.hairpin_score;
	return hairpin.hairpin_score;
}

/*
 double primer::hairpin(int candidate_number, ofstream &fout)
{	
	dimerisation hairpin;
	
	hairpin.hairpin(candidate[candidate_number].sequence, fout);
	candidate[candidate_number].hairpin = hairpin.hairpin_score;
	return hairpin.hairpin_score;
}
 */


int primer::tail_check(int candidate_number, const char* b_sequence, ofstream &fout) 
{
  return sequence_utils::tail_check(candidate[candidate_number].sequence, b_sequence, fout);
}


int primer::sticky_tail_check(int candidate_number)
{
  return  sequence_utils::sticky_tail_check(candidate[candidate_number].sequence);
}

int primer::tail_complementarity_report(const char* sequence)
{
	verbose = TRUE;
	tail_complementarity(sequence);
	verbose = FALSE;
	
	return(TRUE);
}


int primer::tail_complementarity(const char* sequence)
{
	int seqlen = strlen(sequence);
	
	//std::cout << "Tail complementarity check\n";
	
	// Length 3
	if(sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 2]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 3]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 4]))
	{
		if(verbose)
		{
			std::cout << "Length 3 complementarity at " << seqlen - 4;
			std::cout << " in " << sequence << endl;
		}
		else 
		{
			return(TRUE);
		}

	}
	
	/*if(sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 3]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 4]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 5]))
	{
		std::cout << "Length 3 complementarity at " << seqlen - 5;
		std::cout << " in " << sequence << endl;
	}*/
	
	if(sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 4]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 5]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 6]))
	{
		if(verbose)
		{
			std::cout << "Length 3 complementarity at " << seqlen - 6;
			std::cout << " in " << sequence << endl;
		}
		else 
		{
			return(TRUE);
		}
	}
	
	if(sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 5]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 6]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 7]))
	{
		if(verbose)
		{
			std::cout << "Length 3 complementarity at " << seqlen - 7;
			std::cout << " in " << sequence << endl;
		}
		else 
		{
			return(TRUE);
		}
	}
	
	// Length 4
	if(sequence[seqlen - 4] == sequence_utils::nucleotide_complement(sequence[seqlen - 1]) && \
	   sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 2]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 3]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 4]))
	{
		if(verbose)
		{
			std::cout << "Length 4 complementarity at " << seqlen - 4;
			std::cout << " in " << sequence << endl;
		}
		else 
		{
			return(TRUE);
		}
	}
	
	/*if(sequence[seqlen - 4] == sequence_utils::nucleotide_complement(sequence[seqlen - 2]) && \
	   sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 3]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 4]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 5]))
	{
		std::cout << "Length 4 complementarity at " << seqlen - 5;
		std::cout << " in " << sequence << endl;
	}*/
	
	if(sequence[seqlen - 4] == sequence_utils::nucleotide_complement(sequence[seqlen - 3]) && \
	   sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 4]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 5]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 6]))
	{
		if(verbose)
		{
			std::cout << "Length 4 complementarity at " << seqlen - 6;
			std::cout << " in " << sequence << endl;
		}
		else 
		{
			return(TRUE);
		}
	}
	
	/*if(sequence[seqlen - 4] == sequence_utils::nucleotide_complement(sequence[seqlen - 4]) && \
	   sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 5]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 6]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 7]))
	{
		std::cout << "Length 4 complementarity at " << seqlen - 7;
		std::cout << " in " << sequence << endl;
	}*/
	
	// Length 5
	if(sequence[seqlen - 5] == sequence_utils::nucleotide_complement(sequence[seqlen - 2]) && \
	   sequence[seqlen - 4] == sequence_utils::nucleotide_complement(sequence[seqlen - 3]) && \
	   sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 4]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 5]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 6]))
	{
		if(verbose)
		{
			std::cout << "Length 4 complementarity at " << seqlen - 4;
			std::cout << " in " << sequence << endl;
		}
		else 
		{
			return(TRUE);
		}
	}
	
	/*if(sequence[seqlen - 5] == sequence_utils::nucleotide_complement(sequence[seqlen - 3]) && \
	   sequence[seqlen - 4] == sequence_utils::nucleotide_complement(sequence[seqlen - 4]) && \
	   sequence[seqlen - 3] == sequence_utils::nucleotide_complement(sequence[seqlen - 5]) && \
	   sequence[seqlen - 2] == sequence_utils::nucleotide_complement(sequence[seqlen - 6]) && \
	   sequence[seqlen - 1] == sequence_utils::nucleotide_complement(sequence[seqlen - 7]))
	{
		std::cout << "Length 4 complementarity at " << seqlen - 5;
		std::cout << " in " << sequence << endl;
	}*/
	

	return(FALSE); // ie no tail complementarity detected
}

// SELF DIMERISATION METHODS

int primer::self_dimer(int candidate_number)
{
	dimerisation dimer;
	
	dimer.self_dimer(candidate[candidate_number].sequence);
	candidate[candidate_number].self_dimer = dimer.self_dimer_score;
	
	return(TRUE);	
}

int primer::self_dimer(int candidate_number, ofstream &fout)
{
	dimerisation dimer;
	
	dimer.self_dimer(candidate[candidate_number].sequence, fout);
	candidate[candidate_number].self_dimer = dimer.self_dimer_score;
	
	return(TRUE);	
}



int primer::primer_dimer(int candidate_number_a, const char* b_sequence, ofstream &fout)
{
  return sequence_utils::primer_dimer(candidate[candidate_number_a].sequence, b_sequence, fout);
}

int primer::primer_dimer_2(int candidate_number_a, const char* b_sequence)
{
	dimerisation dimer;
	
	if(dimer.primer_dimer_V2(candidate[candidate_number_a].sequence, b_sequence))
	{
		candidate[candidate_number_a].forward_dimer = dimer.forward_dimer_score;
		candidate[candidate_number_a].reverse_dimer = dimer.reverse_dimer_score;
	}
	
	return(0);
		
}



int primer::primer_dimer_2(int candidate_number_a, const char* b_sequence, ofstream &fout)
{
	dimerisation dimer;
	
	if(dimer.primer_dimer_V2(candidate[candidate_number_a].sequence, b_sequence, fout))
	{
		candidate[candidate_number_a].forward_dimer = dimer.forward_dimer_score;
		candidate[candidate_number_a].reverse_dimer = dimer.reverse_dimer_score;
	}
	
	return(0);
	
}



int primer::blast_yeast(int candidate_number, const char* name)
{
  char *s = candidate[candidate_number].sequence;
  char *r = NULL;
  int result = 0;

  r = new char[1+strlen(s)];
  sequence_utils::reverse_complement(s,r);
  result = sequence_utils::blast_yeast(r, name, expectation);
  delete[] r;
  
  return result;
}

int primer::blast_plasmid(int candidate_number, const char* name)
{
  char *s = candidate[candidate_number].sequence;
  char *r = NULL;
  int result = 0;

  r = new char[1+strlen(s)];
  sequence_utils::reverse_complement(s,r);
  result = sequence_utils::blast_plasmid(r, name, expectation);
  delete[] r;
  
  return result;
}


// BLAST - FILE REPORTING VERSIONS

int primer::blast_yeast(int candidate_number, const char* name, ofstream &fout)
{
  char *s = candidate[candidate_number].sequence;
  char *r = NULL;
  int result = 0;

  r = new char[1+strlen(s)];
  sequence_utils::reverse_complement(s,r);
  result = sequence_utils::blast_yeast(r, name, expectation, fout);
  delete[] r;
  
  return result;
}

int primer::blast_plasmid(int candidate_number, const char* name, ofstream &fout)
{
  char *s = candidate[candidate_number].sequence;
  char *r = NULL;
  int result = 0;

  r = new char[1+strlen(s)];
  sequence_utils::reverse_complement(s,r);
  result = sequence_utils::blast_plasmid(r, name, expectation, fout);
  delete[] r;
  
  return result;
}




int primer::blast_seq(int candidate_number, const char* sequence)
{
  char *s = candidate[candidate_number].sequence;
  char *r = NULL;
  int result = 0;

  r = new char[1+strlen(s)];
  sequence_utils::reverse_complement(s,r);
  result = sequence_utils::blast_seq(r, sequence, expectation);
  delete[] r;
  
  return result;
}



int primer::blast_db(int candidate_number, const char* db_name)
{
  char *s = candidate[candidate_number].sequence;
  char *r = NULL;
  int result = 0;

  r = new char[1+strlen(s)];
  sequence_utils::reverse_complement(s,r);
  result = sequence_utils::blast_db(r, db_name, expectation);
  delete[] r;
  return result;

}


int primer::Smith_Waterman(int candidate_number, const char* library_filename)
{
  return sequence_utils::Smith_Waterman(candidate[candidate_number].sequence, library_filename, expectation);
}



int primer::fasta3(int candidate_number, const char* library_filename)
{
  return sequence_utils::fasta3(candidate[candidate_number].sequence, library_filename, expectation);
}



int primer::calculate_temperature(int candidate_number)
{
  annealing_temperature temp;
  candidate[candidate_number].annealing_temperature = temp.primer3_Tm(candidate[candidate_number].sequence);
  return 1;
}
