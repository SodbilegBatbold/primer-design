/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 dimerisation.h
 
 Created by Michael Riley
			Amanda Clare
 
 Date:	10/02/2011
 
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
 THEORY OF LIABILITY, WHETHER IN CONTRACT, ST
 RICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 
 *******************************************************************/

/** \file dimerisation.h
 \brief Dimerisation analysis
 
 This file contains the methods and attributes required for dimerisation analysis
 */

#ifndef DIMERISATION_H
#define DIMERISATION_H

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "sequence_utils.h"

#define ERROR -1
/*
#define TRUE 1
#define FALSE 0


#define ADENINE 65
#define THYMINE 84
#define GUANINE 71
#define CYTOSINE 67
#define URACIL 85
#define ANYNUCLEOTIDE 78
*/

using namespace std;

//! Methods for dimerisation analysis
/**
 Class contains various methods for determining potential dimerisation
 */

class dimerisation: public sequence_utils
{
public:
	dimerisation();
	~dimerisation(){};
	
	/** \brief This is the current recommended method for hairpin prediction. It uses 3' tail 
	    binding affinity
	 */
	int hairpin(const char* sequence);
	
	/**
	 \brief This is the current recommended method for self dimerisation prediction. It uses 3' 
	 tail binding affinity

	 Sets self_dimer_score and self_dimer_location. Returns TRUE
	 */
	int self_dimer(const char* a_sequence);
	/**
	 \brief As for dimerisation::self_dimer(const char* a_sequence) but also outputs to fout.

	 Sets self_dimer_score and self_dimer_location. Returns TRUE
	 */
	int self_dimer(const char* a_sequence, ofstream &fout);
	
	/**
	 \brief This is the current recommended method for pair dimerisation prediction. 
	 It uses 3' tail binding affinity

	 Sets forward_dimer_score, reverse_dimer_score, forward_dimer_location and 
	 reverse_dimer_location. Returns TRUE
	 */
	int pair_dimer(const char* a_sequence, const char* b_sequence);
	
	/**
	 \brief Deprecated
	 
	 Deprecated method for determining potential dimerisation. 
	 */	
	int primer_dimer_V2(const char* a_sequence, const char* b_sequence);

	/**
	 \brief Deprecated
	 
	 Deprecated method for determining potential dimerisation. 
	 */	
	int primer_dimer_V2(const char* a_sequence, const char* b_sequence, ofstream &fout);
	
	/**
	 \brief Percentage complementarity method
	 
	 Old method for determining potential dimerisation. Displays in file all complementary matches over
	 a predefined percentage (default 80%).
	 
	 \sa pair_dimer, set_percentage_complementarity
	 */
	int percentage_complementarity(const char* a_sequence, const char* b_sequence, ofstream &fout);

	/**
	 \brief Percentage complementarity accessor. var must be between 0 and 100.
	*/
	int set_percentage_complementarity(int var);
	
	/** 
	 \brief Thermodynamic analysis for dimerisation prediction
	 
	 This is a container for a third party thermodynamic method for dimerisation prediction. Tests (Dec 2011) 
	 carried out using several currently available thermodynamic methods were disappointing. Currently, the 
	 thermodynamic method has not been made available pending further tests.
	 */
	int thermodynamic_dimerisation(const char* a_sequence, const char* b_sequence);
	
	//char* match;
	
	int dimer_tail_length;
	int high_scoring_limit;
	double GC_score;
	double AT_score;
	double high_score_factor;
	
	double hairpin_score;
	int hairpin_location;
	double self_dimer_score;
	int self_dimer_location;
	double forward_dimer_score; ///< Score for the A sequence (forward primer) tail matched to the entire B sequence
	int forward_dimer_location;
	double reverse_dimer_score; ///< Score for the B sequence (reverse primer) tail matched to the entire A sequence
	int reverse_dimer_location;
	
private:
	bool compare_sequences(const char* a_seq, const char* b_seq, double complementarity, char* match);
	//int nucleotide_complement(int nucleotide);
	int percentage_complementarity_threshold;
};
	
#endif	
