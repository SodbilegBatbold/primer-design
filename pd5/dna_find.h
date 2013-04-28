/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 Created by:	Michael C. Riley and Amanda Clare
 
 Date:			27/03/2011.
 
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

/** \file dna_find.h
 \brief Secondary binding analysis methods
 
 
 */

#ifndef DNAFIND_H
#define DNAFIND_H
 
#include <fstream>
#include <iostream> 
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "global_defs.h"
#include "nsb.h"


//#define DNAFIND_DEBUG

using namespace std;



//! Secondary binding class
/**
 Methods for detecting potential secondary primer binding on the template and potential secondary products.
 
 Example use for genome of S. cerevisiae:
 
     DNAfind my_dnafind("s_cere_genome.fa");
     my_dnafind.set_max_mismatches(5);
     my_dnafind.set_tail_length(20);
     my_dnafind.set_max_viable_product_length(3500);
  
     int products = my_dnafind.search_for_pcr_products(primer_A, primer_B); 
 
     cout << "Number of products = " << products << endl;
 */

class DNAfind : public nsb
{
public:
	//! Constructor
	/** Note that the instantiation of a DNAfind object requires the file name of the genome sequence for the 
	 organism under investigation. Must be in FASTA format with each chromosome separately delineated by >. */
        DNAfind(const char* filename);
	
	//!Destructor
	~DNAfind()
	{
		//if(rc_sequence) delete[] rc_sequence;
		fin.close();
	};
	
	//location_data forward_primer_match_locations[1000]; ///< Calculated results stored here
	//location_data reverse_primer_match_locations[1000]; ///< Calculated results stored here
	
	/** \brief Secondary binding site prediction */
	/** Finds all potential binding site locations. Returns number of binding sites found, locations 
	 of all sites found are in forward_primer_match_locations and reverse_primer_match_locations
	 \sa set_max_mismatches(), set_tail_length() */
	int search_for_binding_sites(const char* sequence);
	
	/** \brief Secondary product prediction */
	/** Finds all potential products. Returns number of products found. 
	 \sa set_max_mismatches(), set_tail_length(), set_max_viable_amplicon_length(); */
	int search_for_pcr_products(const char *forward_sequence, const char *reverse_sequence);
	
	//char* reverse_complement(const char* sequence);
	
	/** \brief Accessor method */
	/** Sets the number of mismatches allowed in the tail of the primer sequence. 
	 Recommend < 30% of tail_length. 0 => exact match. */
	int set_max_mismatches(int var);
	
	/** \brief Accessor method */
	/** Set length of the 3' tail of the primer used when matching up potential for binding. 
	 Max length = 20, default = 12 and recommended min = 8. */
	int set_tail_length(int var);
	
	/** \brief Accessor method */
	/** Sets the maximum length of the secondary PCR product to detect by search_for_pcr_products(). 
	 Max = 5000, default = 3500 */
	int set_max_viable_product_length(int var);
	
	//int set_optimisation(int var);
	bool GC_array_optimisation; ///< Optimisation parameter
	bool report_details;   ///< Outputs to stdout
	
private:
	//ifstream fin;
	
	string reverse_complement(string sequence);
    int reverse_complement(const char* orig_string, char* out_string);
	
	int first_GC_opt(const char* sequence, 
					 const char* dna_template,
					 location_data sequence_match[],
					 int match_count);
	
	int GC_array_opt(const char* sequence,
                     bool is_rc_sequence,
					 const char* dna_template,
					 location_data sequence_match[],
					 int match_count);
	
	int find_sequence(const char* sequence, 
					  const char* dna_template);
	
	int find_sequence(const char* sequence, 
					  const char* dna_template, 
					  location_data sequence_match[]);
	
	int process_chromosome(const char *chromosome, 
						   const char *forward_sequence, 
						   const char *reverse_sequence);
		
	int matches;
	int max_mismatches;
	int max_viable_product_length;
	int tail_length;

};

#endif
