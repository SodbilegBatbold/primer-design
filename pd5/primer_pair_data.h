
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 primer_pair_data.h
 
 Created by:	Michael C. Riley
				Amanda Clare
 
 Date:			05/04/2011
 
 Copyright (c) 2010, 2011, 2012 Aberystwyth University. 
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

/** \file primer_pair_data.h
 \brief Data pertaining to a pair of primers

 */

#ifndef PRIMER_PAIR_DATA_H
#define PRIMER_PAIR_DATA_H

#include "dimerisation.h"
#include "dna_find.h"

#define TRUE 1
#define FALSE 0

//! Primer pair data class
/**
 The main class containing primer pair data
 */

class primer_pair_data
{	
public:
	//primer_pair_data();
	//~primer_pair_data(){};
	
	char forward_sequence[128]; ///< Forward primer sequence. Limited to 128 chars. We do not want primers over 100 nt in length anyway 
	char reverse_sequence[128]; ///< Reverse primer sequence. Limited to 128 chars. We do not want primers over 100 nt in length anyway 
	int forward_index;   ///< Which forward primer this was, out of the array of candidates
	int reverse_index;   ///< Which reverse primer this was, out of the array of candidates
	int location_forward_5_prime_end;  ///< Location of the 5' end. Locations begin at 0, not 1.
	int location_reverse_5_prime_end;  ///< Location of the 5' end. Locations begin at 0, not 1.
	double forward_hairpin_score;  ///< The score for hairpin potential.
	double reverse_hairpin_score;  ///< The score for hairpin potential.
	double forward_self_dimer_score;  ///< The score for self-dimer potential.  
	double reverse_self_dimer_score;  ///< The score for self-dimer potential.  
	double forward_annealing_temperature; ///< The calculated annealing temperature 
	double reverse_annealing_temperature;  ///< The calculated annealing temperature 
	double forward_pair_dimer_score; ///< Dimer score of forward 3' tail on reverse primer
	double reverse_pair_dimer_score; ///< Dimer score of reverse 3' tail on forward primer
	int number_of_pcr_products;  ///< The number of potential products that this pair could make
	double annealing_temperature_difference; ///< The calculated difference between annealing temperatures 
	double moo_score;   ///< Multi objective optimisation score
	
	// Methods
	/** Calls pair_dimer on forward_sequence and reverse_sequence and sets the 
	    forward_pair_dimer_score and reverse_pair_dimer_score. Returns TRUE on success 
	    and FALSE on failure. */
	int pair_dimerisation(void);
	/** Sets number_of_pcr_products by calling DNAfind::search_for_pcr_products. Returns 
	    TRUE on success and FALSE on failure. */
	int pcr_products(DNAfind &mynsb);
};

#endif
