/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 sequence_utils.h
 
 Created by:	Michael C. Riley and Amanda Clare
 
 Date:			21/04/2011
 
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

/** \file sequence_utils.h
 \brief Utilities for sequence analysis
 
 
 */


#ifndef SEQUENCE_UTILS_H
#define SEQUENCE_UTILS_H

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "global_defs.h"



using namespace std;

//! Methods for sequence analysis and manipulation
/**
 * Class for general sequence utilities
 */

class sequence_utils
{
 public:

	sequence_utils()
	{
	};
	
	~sequence_utils(){};
	// int hairpin(const char* sequence, primer_data_structure &score);
	// int hairpin(const char* sequence, primer_data_structure &score, ofstream &fout);
	/** Compares last 15 bases of each sequence for tail complementarity. Returns number of hits */ 
	static int tail_check(const char* a_sequence, const char* b_sequence, ofstream &fout);
	/** Returns TRUE if sequence ends in GC or CG, else FALSE */
	static int sticky_tail_check(const char* query_sequence);
	
	//int self_dimer(const char* a_sequence, primer_data_structure &score);
	//int self_dimer(const char* a_sequence, primer_data_structure &score, ofstream &fout);
	/** Deprecated. Do not use (only works for short primers) */
	static int primer_dimer(const char* a_sequence, const char* b_sequence, ofstream &fout);
	/** Returns the maximum dimer score from forward and reverse primers. 
	    Better to use a dimerisation or primer object instead of this method, 
	    as those would provide more detail.
	    \sa dimerisation::primer_dimer_V2, primer::primer_dimer_2.*/
	static int primer_dimer_2(const char* a_sequence, const char* b_sequence);
	/** Returns the maximum dimer score from forward and reverse primers. Outputs to fout.
	    Better to use a dimerisation or primer object instead of this method, 
	    as those would provide more detail.
	    \sa dimerisation::primer_dimer_V2, primer::primer_dimer_2. */
	static int primer_dimer_2(const char* a_sequence, const char* b_sequence, ofstream &fout);

	/// deprecated
	static int blast_yeast(const char* sequence, const char* name, double expectation);
	/// deprecated
	static int blast_plasmid(const char* sequence, const char* name, double expectation);
	/// deprecated
	static int blast_yeast(const char* sequence, const char* name, double expectation, ofstream &fout);
	/// deprecated
	static int blast_plasmid(const char* sequence, const char* name, double expectation, ofstream &fout);

	/** Use blastn to find the number of hits of the candidate against the given 
	    sequence, with e-value at primer::expectation or better. blastn must be 
	    installed and in your PATH. Returns the number of hits. */
	static int blast_seq(const char* query, const char* sequence, double expectation);

	/** Use blastn to find the number of hits of the candidate against the given 
	    database, with e-value at primer::expectation or better. blastn must be 
	    installed and in your PATH. The database must exist and be pre-formatted. 
	    Returns the number of hits. */
	static int blast_db(const char* query, const char* db_name, double expectation);


	//static int GC_clustering(const char* sequence);

	
// FASTA implementation 17/2/11 & 24/2/11
	/** Uses ssearch35 to search for hits of query against library. Returns the number of hits. */
	static int Smith_Waterman(const char* query, const char* library, double expectation);
	/** Uses fasta35 to search for hits of query against library. Returns the number of hits. */
	static int fasta3(const char* query, const char* library, double expectation);

	/** Returns the number of hits parsed from the file-based output of fasta. */ 
	static int fasta_results_parser(const char* fasta_output_filename);

	/** Returns the number of Gs and Cs in a sequence */
	static int GC_content(const char* sequence);
	/** Returns the count of a particular nucleotide in the sequence. Use #def'd constants, 
	    ADENINE etc for the nucleotide */
	static int nucleotide_content(int nucleotide, const char* sequence);
	/** Returns the complement of a particular nucleotide. Use #def'd constants, 
	    ADENINE etc for the nucleotide */
	static int nucleotide_complement(int nucleotide);
	/** Calculate the reverse complement */ 
	static int reverse_complement(const char* orig_string, char* out_string);
	
private:

	static char complement[22]; 
	//char *rc_sequence;
	
	
// Sequence analysis	
	static int compare_sequences(const char* a_seq, const char* b_seq, double complementarity, ofstream &fout);
	static int blast_results(const char* filename, int length, double expectation);
	static int blast_results(const char* filename, int length, double expectation, ofstream &fout);
	static int blast_db_results(const char* filename, int length, double expectation);
	
};




#endif
