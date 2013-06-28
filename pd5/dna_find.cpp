/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 DNAfind.cpp
 
 Created by:	Michael C. Riley and Amanda Clare
 
 Date:			27/03/2011
 
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

/** \file dna_find.cpp
 \brief Secondary binding analysis methods
 
 
 */

#include "dna_find.h"

#define REVERSE 0;
#define FORWARD 1;

using namespace std;

DNAfind::DNAfind(const char* filename) 
try : nsb(filename)
{
  //Set defaults		
  max_mismatches = 0;
  tail_length = 12; // Max 20
  max_viable_product_length = 3500;
  GC_array_optimisation = TRUE;
  report_details = FALSE;
  //cout << "tail "<< tail_length << endl;
  //cout << "mismatches "<< max_mismatches << endl;
} 
catch (...) {
  throw;
 }



int DNAfind::find_sequence(const char* sequence, 
						   const char* dna_template)
{
	long i;
	int index = 0;
	int errors, j;
	int First_GC;
	long location[1001];
	int count = 1000;
	bool found = FALSE;
	long seq_len = strlen(sequence);
	long dna_len = strlen(dna_template);
	
// Find a location for G or C in the query sequence - we can use this to speed up the search
	
	First_GC = 0; // Use 0 default if no Gs or Cs found in sequence
	
	for(i = 0; i < seq_len; i++)
	{	
		if(sequence[i] == GUANINE || sequence[i] == CYTOSINE)
		{
			First_GC = i; 
			break;
		}
	}
		
	for(i = 0; i < (dna_len - seq_len); i++)
	{
		if(sequence[First_GC] == dna_template[i + First_GC])
		{
			found = TRUE;
			errors = 0;
			for(j = 0; j < seq_len; j++)
			{
				if(sequence[j] != dna_template[i + j])
				{
					if(errors++ >= max_mismatches){found = FALSE; break;}
				}
			}
			
			if(found)
			{
				location[index] = i;
#ifdef DNAFIND_DEBUG				
				cout << "Found " << sequence << " at w " << location[index] << endl;
#endif				
				index++;
				if(index > count)
				{
#ifdef DNAFIND_DEBUG					
					cout << "Error: more than " << count << " matches\n";
#endif
					return(-1);
				}
				
				found = FALSE;
				//if(count >= 30)break;
			}
		}
	}
	
	char rc_sequence[seq_len + 1];	
	sequence_utils::reverse_complement(sequence, rc_sequence);
	
	//find a G or C
	First_GC = 0; // Use 0 default if no Gs or Cs found in sequence
	
	for(i = 0; i < seq_len; i++)if(rc_sequence[i] == GUANINE || rc_sequence[i] == CYTOSINE){First_GC = i; break;}
	
	for(i = 0; i < (dna_len - seq_len); i++)
	{
		if(rc_sequence[First_GC] == dna_template[i + First_GC])
		{
			found = TRUE;
			errors = 0;
			for(j = 0; j < seq_len; j++)
			{
				if(rc_sequence[j] != dna_template[i + j])
				{
					if(errors++ >= max_mismatches){found = FALSE; break;};
				}
			}
			
			if(found)
			{
				location[index] = i + seq_len;
#ifdef DNAFIND_DEBUG
				cout << "Found " << rc_sequence << " at c " << location[index] << endl;
#endif
				index++;
				if(index > count)
				{
#ifdef DNAFIND_DEBUG					
					cout << "Error: more than " << count << " matches\n";
#endif
					return(-1);
				}
				
				found = FALSE;
				//if(count >= 30)break;
			}
		}
	}
	
	matches += index;
	
	if(index == 0)
	{
#ifdef DNAFIND_DEBUG
		cout << "Could not find " << sequence << endl;
#endif		
		return(-2);
	}
	else if(index == 1)
		return(location[0]);
	else 
	{
#ifdef DNAFIND_DEBUG
		cout << "Multiple hits for " << sequence << endl;
#endif
		return(-1);
	}
	
}

char complement[22] = {"T G   C      N     A"};

int nucleotide_complement(int nucleotide)
{
	return(*(complement - 65 + nucleotide));
}

string DNAfind::reverse_complement(string sequence)
{
	int i;
	int sequence_length = sequence.length();
	string rc_sequence; // = new string[sequence_length + 1];
	
	cout << sequence << endl;
	
	for(i = 0; i < sequence_length; i++)
		rc_sequence[i] = nucleotide_complement(sequence[sequence_length - i - 1]);
	
	rc_sequence[i] = 0;
	
	cout << rc_sequence << endl;
	
	return(rc_sequence);	
}

//***************************************************************

int match_sequence(long location, string sequence, string dna_template)
{
	int i, seq_len = sequence.length();
	
	
	for(i = 0; i < seq_len; i++)
	{
		if(sequence[i] != dna_template[location + i])
			return(FALSE);
	}
	return(TRUE);
}

long find_subset(string sequence, string dna_template)
{
	int x = 0;
	int opt_array_limit = 5;
	int subset[opt_array_limit];
	long i, location;
	bool found = FALSE;
	
	long seq_len = sequence.length();
	long dna_len = dna_template.length();
	
	/** Find min freq subset in query sequence
	 */
	
	for(i = 0; i < seq_len; i++)
	{	
		if(sequence[i] == GUANINE || sequence[i] == CYTOSINE)
		{ 
			subset[x] = i;
			if(++x > opt_array_limit)break;
		}
	}
	
	/** Set size of subset */
	int offset_array_size = x;
	
	/** Search dna template for subset */
	for(location = 0; location < (dna_len - seq_len); location++)
	{
		found = TRUE;
		for(x = 0; x < offset_array_size; x++)
		{
			if(sequence[subset[x]] != dna_template[location + subset[x]])
			{	
				found = FALSE;
				break;
			}
		}
		if(found)
		{
			//cout << "*"; // match and evaluate BAV
			
			if(match_sequence(location, sequence, dna_template))
			{
				cout << "\nMatch at " << location << ": " << sequence << endl;
			}
		}
	}
	
	return(location);
}


/** Speed options:
 first_GC_opt() looks for the first S (S = G or C) starting from the 5' end of the sequence
 and notes the location on the sequence as an offset. Any sequence window on the DNA template
 that has an S at the offset is tagged as "subset_found". 
 
 GC_array_opt() looks for up to the first (opt_array_limit + 1) Ss starting from the 3' end and 
 therefore GC_array_opt() is much stricter on GC content at the 3' tail. So much so that 
 BAV testing is largely superfluous and therefore GC_array_opt() is recommended and is default.
 
 Test in second_site.cpp (without BAV evaluation) using:
	hits = my_dnaf.search_for_pcr_products("CTTCGTACTAGTTTTAGCTTACCACCGCTTTTAGTG", "GCTTACATTTGCAGACGTCAAAAAACCTGACAAACC");
 Parameters:
	my_dnaf.set_max_mismatches(4);
	my_dnaf.set_tail_length(16);
	my_dnaf.set_max_viable_product_length(3500);
 
 GC_array_opt() result:
	Time: 800212
	Hits = 2
 first_GC_opt()result:
	Time: 1191999
	Hits = 179
 
 
 */

int DNAfind::GC_array_opt(const char* sequence, 
						  const char* dna_template,
						  location_data sequence_match[],
						  int match_count)
{
	bool found = FALSE;
	bool subset_found = FALSE;
	bool high_BAV = FALSE; // Binding Affinity Value (BAV) exceeds threshold if TRUE
	int i, j;
	long location;
	
	long seq_len = strlen(sequence);
	long dna_len = strlen(dna_template);
	
	/** Find min freq subset in query sequence
	 */
	int x = 0;
	int opt_array_limit = 3; ///< Pattern of (opt_array_limit + 1) Ss will be used
	int subset[opt_array_limit + 1];
	
	for(i = seq_len; i >= 0; i--)
	{	
		if(sequence[i] == GUANINE || sequence[i] == CYTOSINE)
		{
			subset[x] = i;
			if(x++ > opt_array_limit)break;
		}
	}
	
	/** Set size of subset */
	int offset_array_size = x;
	
	/** Search dna template for subset */
	for(location = 0; location < (dna_len - seq_len); location++)
	{
		subset_found = TRUE;
		
		for(x = 0; x < offset_array_size; x++)
		{
			if(sequence[subset[x]] != dna_template[location + subset[x]])
			{	
				subset_found = FALSE;
				break;
			}
		}
		
		if(subset_found)
		{
			found = TRUE;
			int errors = 0;
			for(j = 0; j < seq_len; j++)
			{
				if(sequence[j] != dna_template[location + j])
				{
					if(errors++ >= max_mismatches){found = FALSE; break;}
				}
			}
		
		
		/** If a match is found, evaluate the binding affinity and if the BA value exceeds
		 the threshold, set high_BA to TRUE*/
			high_BAV = TRUE; // evaluate_BAV_P3(); evaluate_BAV_Hgen(); evaluate_BAV_Wallace();
		
			if(found && high_BAV)
			{
				sequence_match[match_count].location = location;
				sequence_match[match_count].sense_strand = TRUE;
			
				//cout << "Found " << sequence << " at " << location << endl;
				match_count++;
				if(match_count > 999)return(match_count);				
				found = FALSE;
			}
		}
	}	
	return(match_count);	
}


int DNAfind::first_GC_opt(const char* sequence, 
						  const char* dna_template,
						  location_data sequence_match[],
						  int match_count)
{
	bool found = FALSE;
	bool subset_found = FALSE;
	bool high_BAV = FALSE; // Binding Affinity Value (BAV) exceeds threshold if TRUE
	int i, j;
	long location;
	
	long seq_len = strlen(sequence);
	long dna_len = strlen(dna_template);
	
/** Optimise search for genomes with low GC content */
	int First_GC = 0; // Use 0 default if no Gs or Cs found in sequence	
	for(i = 0; i < seq_len; i++)
		if(sequence[i] == GUANINE || sequence[i] == CYTOSINE){First_GC = i; break;}
	
/** Search the DNA template (most likely a chromosome) for a match/close match to the sequence*/
	for(location = 0; location < (dna_len - seq_len); location++)
	{
		subset_found = TRUE;
		if(sequence[First_GC] != dna_template[location + First_GC])
			subset_found = FALSE;
		
		if(subset_found)
		{
			found = TRUE;
			int errors = 0;
			for(j = 0; j < seq_len; j++)
			{
				if(sequence[j] != dna_template[location + j])
				{
					if(errors++ >= max_mismatches){found = FALSE; break;}
				}
			}
			
/** If a match is found, evaluate the binding affinity and if the BA value exceeds
	the threshold, set high_BA to TRUE*/
			high_BAV = TRUE;// evaluate_BAV_P3(); evaluate_BAV_Hgen(); evaluate_BAV_Wallace();
			
			if(found && high_BAV)
			{
				sequence_match[match_count].location = location;
				sequence_match[match_count].sense_strand = TRUE;
				
				//cout << "Found " << sequence << " at " << location << endl;
				match_count++;
				if(match_count > 999)return(match_count);				
				found = FALSE;
			}
		}
	}
	return(match_count);
}

/** 
 find_sequence() looks for a match on the DNA template for a match to the sequence.
 Primer binding sites are complement
 */

int DNAfind::find_sequence(const char* sequence, 
						   const char* dna_template, 
						   location_data sequence_match[])
{
	int seq_len = strlen(sequence);
	int w_hits, all_hits;
	
/** Reverse complement sequence to search antisense strand */
	char rc_sequence[seq_len + 1];	
	sequence_utils::reverse_complement(sequence, rc_sequence);
	
	/** Search strands W and then C */
	if(!GC_array_optimisation)
	{
		w_hits = first_GC_opt(sequence, dna_template, sequence_match, 0);
		all_hits = first_GC_opt(rc_sequence, dna_template, sequence_match, w_hits);
	}
	else
	{
		w_hits = GC_array_opt(sequence, dna_template, sequence_match, 0);
		all_hits = GC_array_opt(rc_sequence, dna_template, sequence_match, w_hits);
	}
	
/** Adjust location results for antisense strand */
	for(int i = w_hits; i < all_hits; i++)
	{
		sequence_match[i].sense_strand = FALSE;
		sequence_match[i].location += seq_len; // 5' location is the other end for Crick results
	}

#ifdef DNAFIND_DEBUG
	for(int i = 0; i < all_hits; i++)
	{
		cout << "No " << i << ": " << sequence_match[i].location << ", " << sequence_match[i].sense_strand << endl;
	}	
#endif		
	return(all_hits);
}


int DNAfind::search_for_binding_sites(const char *sequence)
{	
	//fasta reader
	char header_id[32];
	char *s_buffer, *buffer, *token;
	int i = 0;
	long location;
	bool chromosome_for_processing = FALSE;	
	
	char tail_sequence[20];

	// For normal sequence
	int sequence_length = strlen(sequence);
	for(i = 0; i < tail_length; i++) 
		tail_sequence[i] = sequence[i + sequence_length - tail_length];

	/* For reverse complement sequence
	for(i = 0; i < tail_length; i++)
		tail_sequence[i] = sequence[i];
     */
	// Terminate sequence string
	tail_sequence[tail_length] = 0;


#ifdef DNAFIND_DEBUG
	cout << "Sequence " << sequence << endl;
	cout << "Tail is " << tail_sequence << endl;
#endif
	
	string t_buffer;
	string chromosome;
	string temp;
	
	matches = 0;
	
	fin.clear();
	fin.seekg(0, ios::beg);

	while(getline(fin, t_buffer))
	{	
		buffer = new char [t_buffer.size() + 1];
		strcpy(buffer, t_buffer.c_str());
		
		if(buffer[0] == 0x3E)  // fasta header
		{
			if(chromosome_for_processing)
			{
#ifdef DNAFIND_DEBUG				
				cout << "chromo for processing\n";
#endif				
				// Convert c_string to char*
				s_buffer = new char [chromosome.size() + 1];
				strcpy (s_buffer, chromosome.c_str());
				
				location = find_sequence(tail_sequence, s_buffer);
				
				delete[] s_buffer;
				
				//cout << "Found " << location << endl;
				chromosome.clear();
				chromosome_for_processing = FALSE;
			}
			strcpy(header_id, strtok(buffer, " >\n\t\r"));
			chromosome_for_processing = TRUE;
			//cout << header_id << endl;
		}
		else 
		{
			token = strtok(buffer, "\n\t\r");			
			chromosome += temp.assign(token); // Convert char* to c_string
		}
		delete[] buffer;
		//if(i++ > 2) break;
	} 
	
	if(chromosome_for_processing)
	{
#ifdef DNAFIND_DEBUG		
		cout << "Last one\n";
#endif		
		// Convert c_string to char*
		s_buffer = new char [chromosome.size() + 1];
		strcpy (s_buffer, chromosome.c_str());
		
		location = find_sequence(tail_sequence, s_buffer);
		
		delete[] s_buffer;
		
		//cout << "Found " << location << endl;
		chromosome.clear();
		chromosome_for_processing = FALSE;
	}
	
	// matches is a member var updated by find_sequence()
	return(matches);
}

/** Search the chromosome for every match for both forward and reverse sequences
 */

int DNAfind::process_chromosome(const char *chromosome, 
									   const char *forward_sequence, 
									   const char *reverse_sequence)
{	
	int products = 0;
	int i, j;
		 
	 int f_count = find_sequence(forward_sequence, chromosome, forward_primer_match_locations);
	 int r_count = find_sequence(reverse_sequence, chromosome, reverse_primer_match_locations);
	
#ifdef DNAFIND_DEBUG
	cout << "Number of forward primer binding sites = " << f_count << " and for the reverse primer = " << r_count << endl;
	
	for(i = 0; i < f_count; i++) 
	{	
		cout << "Forward no. " << i << " is located at: " << forward_primer_match_locations[i].location << " on "; 
		
		if(forward_primer_match_locations[i].sense_strand)
			cout << "sense strand\n";
		else 
			cout << "antisense strand\n";
	}
	
	for(i = 0; i < r_count; i++)
	{
		cout << "Reverse no. " << i << " is located at: " << reverse_primer_match_locations[i].location << " on ";
		
		if(reverse_primer_match_locations[i].sense_strand)
			cout << "sense strand\n";
		else 
			cout << "antisense strand\n";
	}
#endif
	
	/** Sort data into sense strand and antisense strand locations */
	
	int *sense_loci =  new int[f_count + r_count];
	int *sense_primer =  new int[f_count + r_count];
	int *antisense_loci = new int[f_count + r_count];
	int *antisense_primer = new int[f_count + r_count];
	
	int number_of_sense_loci = 0;
	int number_of_antisense_loci = 0;
	
	for(i = 0; i < f_count; i++)
	{
		if(forward_primer_match_locations[i].sense_strand)
		{
			sense_loci[number_of_sense_loci] = forward_primer_match_locations[i].location;
			sense_primer[number_of_sense_loci++] = FORWARD;
		}
		else
		{
			antisense_loci[number_of_antisense_loci] = forward_primer_match_locations[i].location;
			antisense_primer[number_of_antisense_loci++] = FORWARD;
		}
	}
	
	for(i = 0; i < r_count; i++)
	{
		if(reverse_primer_match_locations[i].sense_strand)
		{
			sense_loci[number_of_sense_loci] = reverse_primer_match_locations[i].location;
			sense_primer[number_of_sense_loci++] = REVERSE;
		}
		else
		{
			antisense_loci[number_of_antisense_loci] = reverse_primer_match_locations[i].location;
			antisense_primer[number_of_antisense_loci++] = REVERSE;
		}
	}	
	
	 
	/** Now search downstream (5' to 3') of each antisense strand location for a sense 
	 location within max_viable_product_length (default 3500).
	 */
	
	char rc_forward_sequence[21];
	char rc_reverse_sequence[21];

	for(i = 0; i < number_of_sense_loci ; i++)
	{
		for(j = 0; j < number_of_antisense_loci; j++)
		{
			if(antisense_loci[j] > sense_loci[i] &&
				antisense_loci[j] - sense_loci[i] < max_viable_product_length)
			{
				products++;
				
				if(report_details)
				{
					cout << sense_loci[i] << ">>>                                    <<<" << antisense_loci[j] << endl;
					if(sense_primer[i]) cout << forward_sequence << ">>>           ";
					else cout << reverse_sequence << ">>>           ";
					
					for(int k = 0 - tail_length; k < 0; k++)
						cout << chromosome[antisense_loci[j] + k];
					cout << endl;
					
					for(int k = 0; k < tail_length; k++)
						cout << chromosome[sense_loci[i] + k];
					//cout << endl;
					
					sequence_utils::reverse_complement(forward_sequence, rc_forward_sequence);
					sequence_utils::reverse_complement(reverse_sequence, rc_reverse_sequence);
					
					if(antisense_primer[j]) cout << "           <<<" << rc_forward_sequence;
					else cout << "           <<<" << rc_reverse_sequence;
					
					cout << endl << "Product length " << antisense_loci[j] - sense_loci[i];
					
				}
				 
			}
		}		
	}	
	delete[] sense_loci;
	delete[] antisense_loci;
	delete[] sense_primer;
	delete[] antisense_primer;
	return(products);
}


/** DNAfind::search_for_pcr_products:
 Takes the forward and reverse primer sequences and searches for primer binding sites that
 would result in a PCR reaction generating a product no larger that max_viable_product_length.
 */

int DNAfind::search_for_pcr_products(const char *forward_primer, 
									 const char *reverse_primer)
{	
	//fasta reader using strings
	
	int i;
	char header_id[32];
	char *buffer, *token;
	bool chromosome_for_processing = FALSE;
	int number_of_products_found = 0;	
	char *str_chromosome;	
	string t_buffer;
	string chromosome;
	string temp;
	
	// Make Tails: Make sure that the sequence is not shorter than DNAfind::tail_length 
	// otherwise make the tails the same as primer length.
	
	int fwd_seq_len = strlen(forward_primer);
	int rev_seq_len = strlen(reverse_primer);
	int fwd_tail_length;
	int rev_tail_length;
	
	if(DNAfind::tail_length > fwd_seq_len)
		fwd_tail_length = fwd_seq_len;
	else
		fwd_tail_length = tail_length;
	
	if(DNAfind::tail_length > rev_seq_len)
		rev_tail_length = rev_seq_len;
	else
		rev_tail_length = tail_length;
	
	char forward_sequence_tail[21]; // tails no longer than 20 bp
	char reverse_sequence_tail[21];
	
	// Make 3' tail sequences
	for(i = 0; i < fwd_tail_length; i++)
		forward_sequence_tail[i] = forward_primer[fwd_seq_len - fwd_tail_length + i];
	forward_sequence_tail[fwd_tail_length] = 0;
	
	for(i = 0; i < rev_tail_length; i++)
		reverse_sequence_tail[i] = reverse_primer[rev_seq_len - rev_tail_length + i];
	reverse_sequence_tail[rev_tail_length] = 0;
	
	// Read genome.fa and process each chromosome
	 
	fin.clear();
	fin.seekg(ios::beg);
	
	while(getline(fin, t_buffer))
	{	
		buffer = new char [t_buffer.size() + 1];
		strcpy(buffer, t_buffer.c_str());
		
		if(buffer[0] == 0x3E)  // fasta header
		{
			if(chromosome_for_processing)
			{				
				// Convert c_string to char*
				str_chromosome = new char [chromosome.size() + 1];
				strcpy (str_chromosome, chromosome.c_str());
				
				if(report_details)
				{
					int i = process_chromosome(str_chromosome, forward_sequence_tail, reverse_sequence_tail);
					if(i) cout << " on sequence " << header_id << endl << endl;
						
					number_of_products_found += i;
				}
				else number_of_products_found += process_chromosome(str_chromosome, forward_sequence_tail, reverse_sequence_tail);
				
				chromosome.clear();
				delete[] str_chromosome;
				chromosome_for_processing = FALSE;
			}
			
			strcpy(header_id, strtok(buffer, " >\n\t\r"));
			chromosome_for_processing = TRUE;					
		}
		else 
		{
			token = strtok(buffer, "\n\t\r");			
			chromosome += temp.assign(token); // Convert char* to c_string
		}
		delete[] buffer;
	} 
	
	if(chromosome_for_processing)
	{	
		// Convert c_string to char*
		str_chromosome = new char [chromosome.size() + 1];
		strcpy (str_chromosome, chromosome.c_str());
		
		if(report_details)
		{
			int i = process_chromosome(str_chromosome, forward_sequence_tail, reverse_sequence_tail);
			if(i) cout << " on sequence " << header_id << endl << endl;
			
			number_of_products_found += i;
		}
		else number_of_products_found += process_chromosome(str_chromosome, forward_sequence_tail, reverse_sequence_tail);
		
		chromosome.clear();
		delete[] str_chromosome;
		chromosome_for_processing = FALSE;
	}
	
	return(number_of_products_found);
}

// Accessor Methods

int DNAfind::set_max_mismatches(int var)
{
	if(var < 0 || var > tail_length) // Actually best if < 30% of tail_length. 0 => exact match.
		return(FALSE);
	else 
	{
		max_mismatches = var;
		return(TRUE);
	}
}

int DNAfind::set_tail_length(int var)
{
	if(var < 0 || var > 20) 
		return(FALSE);
	else 
	{
		tail_length = var;
		return(TRUE);
	}		
}

int DNAfind::set_max_viable_product_length(int var)
{
	if(var < 0 || var > 5000) 
		return(FALSE);
	else 
	{
		max_viable_product_length = var;
		return(TRUE);
	}	
}



