/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 dimerisation.cpp
 
 Created by Michael Riley on 10/08/2010 - mhr"at"aber.ac.uk
 
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

/** \file dimerisation.cpp
 \brief Dimerisation analysis
 
 This file contains the methods and attributes required for dimerisation analysis
 */

#include "dimerisation.h"


//#define TESTING

dimerisation::dimerisation(void)
{
// DEFAULT SETTINGS
	dimer_tail_length = 8;
	high_scoring_limit = 2;
	
	GC_score = 2;
	AT_score = 1;
	high_score_factor = 2;
	
	hairpin_score = 0;
	self_dimer_score = 0;
	forward_dimer_score = 0;
	reverse_dimer_score = 0;
	
}

bool dimerisation::compare_sequences(const char* a_seq, const char* b_seq, double complementarity, char* match)
{
	int sequence_length = strlen(a_seq);
	int i = 0;
	bool HIT = FALSE;
	
	if(sequence_length != (int)strlen(b_seq))
	{
#ifdef TESTING		
		std::cout << "Sequences are different lengths\n"; 
#endif		
		return(ERROR);
	}
	
	int number_of_matches = 0;
	for(i = 0; i < sequence_length; i++)
	{ 
		if(a_seq[i] == nucleotide_complement(b_seq[i]))
		{
			match[i] = 124; // In DOS | = 179, in File text | = 124
			number_of_matches++;
		}
		else match[i] = 32;  // a space
	}	
	match[i] = 0; // terminate string
	
	if(number_of_matches > sequence_length * complementarity) HIT = TRUE;
	
	return(HIT);
}

int dimerisation::hairpin(const char* sequence)
{	
	int len = strlen(sequence);
	char tail[dimer_tail_length + 1];
	double score;
	int i, j;
	
	// flip last (dimer_tail_length) bases at 3' end of a_sequence
	
	for(i = 0; i < dimer_tail_length; i++) tail[i] = sequence[len - 1 - i]; tail[i] = 0;
	
#ifdef TESTING
	std::cout << sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < len - (2 * dimer_tail_length); j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > hairpin_score)
		{
			hairpin_score = score;
			hairpin_location = j;
		}
	}
	
#ifdef TESTING	
	std::cout << "A: " << hairpin_score << " at " << hairpin_location << endl;
#endif	
	
	return(TRUE);
	
}


int dimerisation::self_dimer(const char* a_sequence, ofstream &fout)
{		
	int a_len = strlen(a_sequence);
	char tail[dimer_tail_length + 1];	
	int i, j;
	double score;
	
// Initialise
	
	self_dimer_score = 0;
	
// Get and flip last 8 bases at 3' end of a_sequence
	
	for(i = 0; i < dimer_tail_length; i++) tail[i] = a_sequence[a_len - 1 - i]; tail[i] = 0;
	
#ifdef TESTING
	std::cout << a_sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < a_len - dimer_tail_length; j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(a_sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > self_dimer_score)
		{
			self_dimer_score = score;
			self_dimer_location = j;
		}
	}

#ifdef TESTING	
	std::cout << "A: " << self_dimer_score << " at " << self_dimer_location << endl;
#endif	
	
	fout << "tail binds: " << self_dimer_score << " at " << self_dimer_location << endl;
	
	return(TRUE);
}

int dimerisation::self_dimer(const char* a_sequence)
{	
	int a_len = strlen(a_sequence);
	char tail[dimer_tail_length + 1];	
	int i, j;
	double score;
	
// Initialise
	
	self_dimer_score = 0;
	
	// Get and flip last 8 bases at 3' end of a_sequence
	for(i = 0; i < dimer_tail_length; i++) tail[i] = a_sequence[a_len - 1 - i]; tail[i] = 0;

#ifdef TESTING
	std::cout << a_sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < a_len - dimer_tail_length; j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(a_sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > self_dimer_score)
		{
			self_dimer_score = score;
			self_dimer_location = j;
		}
	}
	
#ifdef TESTING	
	std::cout << "A: " << self_dimer_score << " at " << self_dimer_location << endl;
#endif
	
	// Re-Initialise
	
	self_dimer_score = 0;
	
	// Get last 8 bases at 3' end of a_sequence
	for(i = 0; i < dimer_tail_length; i++) tail[i] = a_sequence[a_len - dimer_tail_length + i]; tail[i] = 0;
	
#ifdef TESTING
	std::cout << a_sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < a_len - dimer_tail_length; j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(a_sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > self_dimer_score)
		{
			self_dimer_score = score;
			self_dimer_location = j;
		}
	}
	
#ifdef TESTING	
	std::cout << "A: " << self_dimer_score << " at " << self_dimer_location << endl;
#endif
	
	return(TRUE);
}

int dimerisation::pair_dimer(const char* a_sequence, const char* b_sequence)
{
	int a_len = strlen(a_sequence);
	int b_len = strlen(b_sequence);
	char tail[dimer_tail_length + 1];	
	int i, j;
	double score;
	
	// Initialise
	
	forward_dimer_score = 0;
	reverse_dimer_score = 0;
	
	// Get and flip last 8 bases at 3' end of a_sequence
	for(i = 0; i < dimer_tail_length; i++) tail[i] = a_sequence[a_len - 1 - i]; tail[i] = 0;
	
#ifdef TESTING	
	std::cout << a_sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < b_len - dimer_tail_length; j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(b_sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > forward_dimer_score)
		{
			forward_dimer_score = score;
			forward_dimer_location = j;
		}
	}
	
	// Now do sequence B
	
	for(i = 0; i < dimer_tail_length; i++) tail[i] = b_sequence[b_len - 1 - i]; tail[i] = 0;
#ifdef TESTING	
	std::cout << b_sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < a_len - dimer_tail_length; j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(a_sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > reverse_dimer_score)
		{
			reverse_dimer_score = score;
			reverse_dimer_location = j;
		}
	}
	
#ifdef TESTING	
	std::cout << "A: " << forward_dimer_score << " at " << forward_dimer_location;
	std::cout << ", B: " << reverse_dimer_score<< " at " << reverse_dimer_location << endl;
#endif
	
	return(TRUE);
}


int dimerisation::primer_dimer_V2(const char* a_sequence, const char* b_sequence)
{
	int a_len = strlen(a_sequence);
	int b_len = strlen(b_sequence);
	char tail[dimer_tail_length + 1];	
	int i, j;
	double score;
	
	// Initialise
	
	forward_dimer_score = 0;
	reverse_dimer_score = 0;

	// Get and flip last 8 bases at 3' end of a_sequence
	for(i = 0; i < dimer_tail_length; i++) tail[i] = a_sequence[a_len - 1 - i]; tail[i] = 0;

#ifdef TESTING	
	std::cout << a_sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < b_len - dimer_tail_length; j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(b_sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > forward_dimer_score)
		{
			forward_dimer_score = score;
			forward_dimer_location = j;
		}
	}
	
	// Now do sequence B
	
	for(i = 0; i < dimer_tail_length; i++) tail[i] = b_sequence[b_len - 1 - i]; tail[i] = 0;
#ifdef TESTING	
	std::cout << b_sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < a_len - dimer_tail_length; j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(a_sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > reverse_dimer_score)
		{
			reverse_dimer_score = score;
			reverse_dimer_location = j;
		}
	}
	
#ifdef TESTING	
	std::cout << "A: " << forward_dimer_score << " at " << forward_dimer_location;
	std::cout << ", B: " << reverse_dimer_score<< " at " << reverse_dimer_location << endl;
#endif
	
	return(TRUE);
}

int dimerisation::primer_dimer_V2(const char* a_sequence, const char* b_sequence, ofstream &fout)
{
	int a_len = strlen(a_sequence);
	int b_len = strlen(b_sequence);
	char tail[dimer_tail_length + 1];	
	int i, j;
	double score;
	
	// Initialise
	
	forward_dimer_score = 0;
	reverse_dimer_score = 0;	

	// Get and flip last 8 bases at 3' end of a_sequence
	for(i = 0; i < dimer_tail_length; i++) tail[i] = a_sequence[a_len - 1 - i]; tail[i] = 0;

#ifdef TESTING	
	fout << a_sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < b_len - dimer_tail_length; j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(b_sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > forward_dimer_score)
		{
			forward_dimer_score = score;
			forward_dimer_location = j;
		}
	}
	
	// Repeat for sequence B
	
	for(i = 0; i < dimer_tail_length; i++) tail[i] = b_sequence[b_len - 1 - i]; tail[i] = 0;

#ifdef TESTING	
	fout << b_sequence << endl << tail << endl;
#endif
	
	for(j = 0; j < a_len - dimer_tail_length; j++)
	{
		score = 0;
		
		for(i = 0; i < dimer_tail_length; i++)
		{
			if(tail[i] == nucleotide_complement(a_sequence[j + i]))
			{
				if(tail[i] == CYTOSINE || tail[i] == GUANINE)
				{
					if(i < high_scoring_limit)score += GC_score * high_score_factor;
					else score += GC_score;
				}
				else
				{
					if(i < high_scoring_limit)score += AT_score * high_score_factor;
					else score += AT_score;
				}
			}
		}
		if(score > reverse_dimer_score)
		{
			reverse_dimer_score = score;
			reverse_dimer_location = j;
		}
	}
	
	
	fout << "tail A on B: " << forward_dimer_score << " at " << forward_dimer_location;
	fout << ", tail B on A: " << reverse_dimer_score<< " at " << reverse_dimer_location << endl;
	
	
	return(TRUE);
}


int dimerisation::percentage_complementarity(const char* a_sequence, const char* b_sequence, ofstream &fout)
{	
	int a_len = strlen(a_sequence);
	int b_len = strlen(b_sequence);
	char seq_a[a_len];
	char seq_b[b_len];
	char flip_seq[b_len];
	int i, j, k, pos, number_of_hits = 0; 
	int min_window_length = 3;
	
	int maxlen = 1000;
	
	if(a_len > b_len) maxlen = a_len;
	else maxlen = b_len;
	
	char match[maxlen];
	
	char msg[250];
	char spaces[250];
	
	double percentage = percentage_complementarity_threshold/100;
	/*
	for(pos = min_window_length - b_len; pos <= a_len - min_window_length; pos++)
	{
		j = 0;
		for(i = 0; i < b_len; i++)
		{
			if((pos + i) >= 0 && (pos + i) < a_len)
			{
				seq_a[j] = a_sequence[pos + i];
				seq_b[j] = b_sequence[i];
				j++;
			}			
		}
		seq_a[j] = 0x00;
		seq_b[j] = 0x00;
		if(compare_sequences(seq_a, seq_b, 0.5)) number_of_hits++;
	}*/
	
	// Flip B sequence and fully convolve
	
	for(j = 0; j < b_len; j++) flip_seq[j] = b_sequence[b_len - 1 - j];
	flip_seq[j] = 0;
	
	
	for(pos = min_window_length - b_len; pos <= a_len - min_window_length; pos++)
	{
		j = 0;
		for(i = 0; i < b_len; i++)
		{
			if((pos + i) >= 0 && (pos + i) < a_len)
			{
				seq_a[j] = a_sequence[pos + i];
				seq_b[j] = flip_seq[i];
				j++;
			}			
		}
		seq_a[j] = 0x00;
		seq_b[j] = 0x00;
		if(compare_sequences(seq_a, seq_b, percentage, match))
		{
			sprintf(msg, "Pos = %d\n", pos);
			cout << msg << endl;
			if(pos > 0)
			{
				for(i = 0; i < pos; i++)spaces[i] = 32; spaces[i] = 0;
				
				cout << a_sequence << endl;
				sprintf(msg, "%s%s", spaces, match);
				cout << msg << endl;
				sprintf(msg, "%s%s", spaces, flip_seq);
				cout << msg << endl;								
			}
			else 
			{
				k = 0;
				for(i = 0; i > pos; i--)spaces[k++] = 32; spaces[k] = 0;
				
				sprintf(msg, "%s%s", spaces, a_sequence);
				cout << msg << endl;	
				sprintf(msg, "%s%s", spaces, match);
				cout << msg << endl;
				cout << flip_seq << endl;
			}

			number_of_hits++;
		}
	}
	
	return(number_of_hits);
}

int dimerisation::set_percentage_complementarity(int var)
{
	if(var >= 0 && var <= 100)
	{
		percentage_complementarity_threshold = var;
		return(TRUE);
	}
	else 
	{
		return(FALSE);
	}	
}

