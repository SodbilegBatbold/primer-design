
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 primer_data.h
 
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

/** \file primer_data.h
 \brief Class holding the data for a single primer


 */



#ifndef PRIMER_DATA_H
#define PRIMER_DATA_H

#include <cstring>
#include "dimerisation.h"

#define TRUE 1
#define FALSE 0

#define ADENINE 65
#define THYMINE 84
#define GUANINE 71
#define CYTOSINE 67
#define URACIL 85
#define ANYNUCLEOTIDE 78

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
		binding_A = 0;
		binding_B = 0;
		seqsim_matches = 0;
		moo = 0.0;
		products = 0;
		annealing_temperature = 0;
	}
	
	char sequence[128];       ///< Primer sequence. Limited to 128 chars. We do not want primers over 100 nt in length anyway 
	int location_5_prime_end; ///< Location of the 5' end. Locations begin at 0, not 1. 
	double hairpin;           ///< The score for hairpin potential. 
	double self_dimer;        ///< The score for self-dimer potential.   
	int binding_A;	///< The number of binding hits to the A genome (1 > indicates secondary binding)
	int binding_B;	///< The number of binding hits to the B genome (1 > indicates secondary binding)
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
	int get_GC_content(void);
};

#endif
