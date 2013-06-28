
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 sequence_utils.cpp
 
 Created by Michael Riley, Amanda Clare on 11/02/2011
 
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

/** \file sequence_utils.cpp
 \brief Utilities for sequence analysis
 
 
 */

#include "sequence_utils.h"
#include "dimerisation.h"

#define TRUE 1
#define FALSE 0

// initialise statics
char sequence_utils::complement[22] = {"T G   C      N     A"};

int sequence_utils::nucleotide_complement(int nucleotide)
{
	return(*(complement - 65 + nucleotide));
}

/*char* sequence_utils::reverse_complement(const char* sequence)
{
	int i;
	int sequence_length = strlen(sequence);
	rc_sequence = new char[sequence_length + 1];
	
	for(i = 0; i < sequence_length; i++)
		rc_sequence[i] = nucleotide_complement(sequence[sequence_length - i - 1]);
	
	rc_sequence[i] = 0;
	
	return(rc_sequence);	
}*/

int sequence_utils::reverse_complement(const char* orig_string, char* out_string)
{
	int i;
	int orig_string_length = strlen(orig_string);
	
	for(i = 0; i < orig_string_length; i++)
	{	
		out_string[i] = nucleotide_complement(orig_string[orig_string_length - i - 1]);
	}	
	out_string[i] = 0;
	
	//cout << out_string << endl;
	return(1);
}

int sequence_utils::GC_content(const char* sequence)
{
	unsigned int GC_number = 0;	
	
	for(unsigned int i = 0; i < strlen(sequence); i++) 
		if(sequence[i] == CYTOSINE || sequence[i] == GUANINE) GC_number++;
	
	return(GC_number);
}

int sequence_utils::nucleotide_content(int nucleotide, const char* sequence)
{
	unsigned int nucleotide_number = 0;	
	
	for(unsigned int i = 0; i < strlen(sequence); i++) 
		if(sequence[i] == nucleotide) nucleotide_number++;
	
	return(nucleotide_number);
}







int sequence_utils::compare_sequences(const char* a_seq, const char* b_seq, double complementarity, ofstream &fout)
{
        unsigned int sequence_length = strlen(a_seq), i;
	char match[sequence_length];
	bool HIT = FALSE;
	
	if(sequence_length != strlen(b_seq))
	{
#ifdef TESTING
		std::cout << "Sequences are different lengths\n";
#endif
		return(-1);
	}
	
	// Display in text: matches
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
	
	match[i] = 0;
	
	if(number_of_matches > sequence_length * complementarity) 
	{		
		fout << a_seq << endl;
		fout << match << endl;
		fout << b_seq << endl;
		
		HIT = TRUE;
	}
	return(HIT);
}

int sequence_utils::tail_check(const char* a_sequence, const char* b_sequence,  ofstream &fout)
{
	int a_length = strlen(a_sequence);
	int b_length = strlen(b_sequence);
	char seq_a[16];
	char seq_b[16];
	int i, j, number_of_hits = 0;	
	
	for(j = 15; j > 2; j--)
	{
		for(i = 0; i < j; i++)
		{
			seq_a[i] = a_sequence[i + a_length - j];
			seq_b[i] = b_sequence[b_length - i - 1];
		}
		seq_a[i] = 0x00;
		seq_b[i] = 0x00;
		
		//cout << seq_a << endl << seq_b << endl;
		if(compare_sequences(seq_a, seq_b, 0.5, fout)) number_of_hits++; 
	}
 	return(number_of_hits);
	
}

int sequence_utils::sticky_tail_check(const char* query_sequence)
{
	int length = strlen(query_sequence);
	bool HIT = FALSE;
	
	if(query_sequence[length - 2] == nucleotide_C && query_sequence[length - 1] == nucleotide_G) HIT = TRUE;
	if(query_sequence[length - 2] == nucleotide_G && query_sequence[length - 1] == nucleotide_C) HIT = TRUE;
	
	return(HIT);
}


int sequence_utils::primer_dimer(const char* a_sequence, const char* b_sequence, ofstream &fout)
{	
	int a_len = strlen(a_sequence);
	int b_len = strlen(b_sequence);
	char seq_a[a_len];
	char seq_b[b_len];
	char flip_seq[b_len];
	int i, j, pos, number_of_hits = 0; 
	int min_window_length = 3;
	
	char msg[250];
	
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
		if(compare_sequences(seq_a, seq_b, 0.5, fout)) number_of_hits++;
	}
	
	// Now flip sequence and fully convolve
	
	for(j = 0; j < b_len; j++) flip_seq[j] = b_sequence[b_len - 1 - j];
	
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
		if(compare_sequences(seq_a, seq_b, 0.5, fout))
		{
			sprintf(msg, "Pos = %d\n", pos);
			fout << msg << endl;
			number_of_hits++;
		}
	}
	
	return(number_of_hits);
}


int sequence_utils::primer_dimer_2(const char* a_sequence, const char* b_sequence)
{
	dimerisation dimer;

	if(dimer.primer_dimer_V2(a_sequence, b_sequence)) {
	  int max_score = dimer.forward_dimer_score;
	  if (dimer.reverse_dimer_score > max_score) {
	    max_score = dimer.reverse_dimer_score;
	  }
	  return(max_score);
	} else {
	  return(ERROR); 
	}
}


int sequence_utils::primer_dimer_2(const char* a_sequence, const char* b_sequence, ofstream &fout)
{
	dimerisation dimer;
	
	if(dimer.primer_dimer_V2(a_sequence, b_sequence, fout)) {
	  int max_score = dimer.forward_dimer_score;
	  if (dimer.reverse_dimer_score > max_score) {
	    max_score = dimer.reverse_dimer_score;
	  }
	  return(max_score);
	} else {
	  return(ERROR); 
	}
}


int sequence_utils::blast_yeast(const char* sequence, const char* name, double expectation)
{
  char filename[64], command[256]; //, temp[256];		
   	ofstream fseq; // for sequence to be BLASTed
	int number_of_hits = 0;
   	

	
	sprintf(filename, "%s.fa", name);
	fseq.open(filename);
	fseq << ">" << name << endl << sequence << endl;
	fseq.close();
	
//	for(int i = 1; i <= 17; i++)
//	{
		//sprintf(command, "./blastn -subject chr%d_s_cerevisiae.fa -query %s.fa -evalue 0.5 -word_size 7 -dust no -ungapped -out test%d.txt",i,name,i);	
	sprintf(command, "./blastn -subject s_cere_genome.fa -query %s.fa -evalue 0.5 -word_size 7 -dust no -ungapped -out %s.out",name,name);	
		system(command);
		
		//cout << "Chromosome " << i << endl << "=============\n\n";
		//sprintf(temp, "\nCHROMOSOME %d", i);
		//analysis_text.add_content(temp);
		
		sprintf(filename, "%s.out", name);
		number_of_hits = blast_results(filename, strlen(sequence), expectation);
//	}
	
	//system("c:/blast/bin/bl2seq -p blastn -g F -j test_primer.fa -i test_primer.fa -o p_dimer.res");
	//fclose(stdout);
	return(number_of_hits);	
}

int sequence_utils::blast_plasmid(const char* sequence, const char* name, double expectation)
{
	char filename[64], command[256];		
   	ofstream fout;
	int number_of_hits = 0;
	//char rc_sequence[strlen(sequence)];
	
	sprintf(filename, "%s.fa", name);
	fout.open(filename);
	fout << ">" << name << endl << sequence << endl;
	fout.close();
	
	sprintf(command, "./blastn -subject ura3_rc.fa -query %s.fa -evalue 0.5 -word_size 7 -dust no -ungapped  -out %s.res",name,name);
	//sprintf(command, "c:/blast/bin/bl2seq -p blastn -g F -j ura3.fa -i %s.fa -o %s.res",name, name);	
	system(command);
	
	//cout << "Plasmid \n" << "=========\n\n";
	//analysis_text.add_content("PLASMID");
	
	sprintf(filename, "%s.res", name);
	number_of_hits = blast_results(filename, strlen(sequence), expectation);
	
	return(number_of_hits);	
}

int sequence_utils::blast_results(const char* filename, int length, double expectation)
{	
	char buffer[256];
	char temp[256];
	//char e_value[256];
	char* token;
	int number_of_hits = 0;
	//char* token2;
   	ifstream fin(filename);
	bool HIT, EVA;
   	
   	while(fin.getline(buffer, 255))
   	{
	 	if(!strncmp(buffer," Score", 6))
	 	{ 
		 	token = strstr(buffer, "Expect") + 9;
			//cout << "token: " << token << endl;
		 	if(atof(token) < expectation) 
		 	{
			 	
		 		// Get lines 4-6 after this E value from BLAST results
				HIT = FALSE; EVA = TRUE;
			 	for(int i = 0; i < 7; i++) 
			 	{
				
				 	fin.getline(buffer, 255);
					
					if(i == 3)
					{
						strcpy(temp, buffer);
						strtok(temp, " ");
						token = strtok(NULL, " ");
						
						// Set OK if matched sequence begins within 3 bp of query sequence
						// noting that by BLASTing the reverse complement, the 3' end is 
						// at the beginning of the sequence
						if(atol(token) < 3)HIT = TRUE;
						
						strtok(NULL, " ");
						token = strtok(NULL, " ");
						
						// Set OK if matched sequence ends within 3 bp of query sequence
						//if(atol(token) > length - 3)OK = 1;
					}
					if(HIT && EVA) // We have a hit
					{
						number_of_hits++;
						EVA = FALSE;
					}
			 	
			 	}
		 	} 
	 	}	 		  	
   	} 
   	return(number_of_hits);	
}

// BLAST - FILE REPORTING VERSIONS

int sequence_utils::blast_yeast(const char* sequence, const char* name, double expectation, ofstream &fout)
{
        char filename[64], command[256]; //, temp[256];		
   	ofstream fseq; // for sequence to be BLASTed
	int number_of_hits = 0;
	int hits;
   	
   	//freopen("primer_result.txt","a", stdout);
	
	sprintf(filename, "%s.fa", name);
	fseq.open(filename);
	fseq << ">" << name << endl << sequence << endl;
	fseq.close();
	
	for(int i = 1; i <= 17; i++)
	{
		sprintf(command, "./blastn -subject chr%d_s_cerevisiae.fa -query %s.fa -evalue 0.5 -word_size 7 -dust no -ungapped -out test%d.txt",i,name,i);	
		system(command);
		
		fout << "CHROMOSOME " << i << endl;
		
		sprintf(filename, "test%d.txt", i);
		hits = blast_results(filename, strlen(sequence), expectation, fout);
		number_of_hits += hits;
		
		if(!hits) fout << endl; // Just a bit of file formatting
	}
	
	//system("c:/blast/bin/bl2seq -p blastn -g F -j test_primer.fa -i test_primer.fa -o p_dimer.res");
	//fclose(stdout);
	return(number_of_hits);	
}

int sequence_utils::blast_plasmid(const char* sequence, const char* name, double expectation, ofstream &fout)
{
	char filename[64], command[256];		
   	ofstream fseq;
	int number_of_hits = 0;
	//char rc_sequence[strlen(sequence)];
	
	sprintf(filename, "%s.fa", name);
	fseq.open(filename);
	fseq << ">" << name << endl << sequence << endl;
	fseq.close();
	
	sprintf(command, "./blastn -subject ura3_rc.fa -query %s.fa -evalue 0.5 -word_size 7 -dust no -ungapped  -out %s.res",name,name);
	//sprintf(command, "c:/blast/bin/bl2seq -p blastn -g F -j ura3.fa -i %s.fa -o %s.res",name, name);	
	system(command);
	
	fout << "PLASMID \n";
	
	sprintf(filename, "%s.res", name);
	number_of_hits = blast_results(filename, strlen(sequence), expectation, fout);
	
	if(!number_of_hits) fout << endl; // Just a bit of file formatting
	
	return(number_of_hits);	
}

int sequence_utils::blast_results(const char* filename, int length, double expectation, ofstream &fout)
{	
	char buffer[256];
	char temp[256];
	char e_value[256];
	char* token;
	int number_of_hits = 0;
	//char* token2;
   	ifstream fin(filename);
	bool HIT, EVA;
   	
   	while(fin.getline(buffer, 255))
   	{
	 	if(!strncmp(buffer," Score", 6))
	 	{ 
		        // cout << buffer << endl;
		 	token = strstr(buffer, "Expect") + 9;
		 	if(atof(token) < expectation) 
		 	{
			 	sprintf(e_value, "E = %s\tCandidate primer length = %d ", token, length);
			 	
			 	//analysis_text.add_content("\n");
			 	
		 		// Get lines 4-6 after this E value from BLAST results
				HIT = FALSE; EVA = TRUE;
			 	for(int i = 0; i < 7; i++) 
			 	{
				 	fin.getline(buffer, 255);
					
					if(i == 3)
					{
						strcpy(temp, buffer);
						strtok(temp, " ");
						token = strtok(NULL, " ");
						
						// Set OK if matched sequence begins within 3 bp of query sequence
						// noting that by BLASTing the reverse complement, the 3' end is 
						// at the beginning of the sequence
						if(atol(token) < 3)HIT = TRUE;
						
						strtok(NULL, " ");
						token = strtok(NULL, " ");
						//cout << "This should be a number: " << token << endl;
						//cout << "This is the length of the sequence: " << length << endl;
						
						// Set OK if matched sequence ends within 3 bp of query sequence
						//if(atol(token) > length - 3)OK = 1;
					}
					if(HIT && EVA) // We have a hit
					{
						fout << e_value << endl;
						number_of_hits++;
						EVA = FALSE;
					}
				 	if(i > 2 && HIT) fout << buffer << endl;
			 	}
		 	}
	 	}	 		  	
   	} 
  	//cout << number_of_hits << endl;
   	return(number_of_hits);	
}



int sequence_utils::blast_seq(const char* query, const char* sequence, double expectation)
{
	char filename[64], command[256];		
   	ofstream fseq; // for sequence to be BLASTed
	int number_of_hits = 0;
   	

       	sprintf(filename, "query.fa");
	fseq.open(filename);
	fseq << ">query" << endl << query << endl;
	fseq.close();

	sprintf(filename, "seq.fa");
	fseq.open(filename);
	fseq << ">seq" << endl << sequence << endl;
	fseq.close();
	
	sprintf(command, "blastn -query query.fa -subject seq.fa -evalue 0.5 -out blastseqout.txt -dust no -word_size 7 -ungapped");	
	system(command);
		
	number_of_hits += blast_results("blastseqout.txt", strlen(query), expectation);
	//std::cout << "Num hits seq:" << number_of_hits << std::endl;
	//std::cout << query << std::endl;
	//std::cin.get();
	return(number_of_hits);	
}


int sequence_utils::blast_db(const char* query, const char* db_name, double expectation)
{
	char filename[64], command[256];		
   	ofstream fseq; // for sequence to be BLASTed
	int number_of_hits = 0;
   	
       	sprintf(filename, "query.fa");
	fseq.open(filename);
	fseq << ">query" << endl << query << endl;
	fseq.close();

	sprintf(command, "blastn -db %s -query query.fa -evalue 0.5 -word_size 7 -dust no -ungapped  -out blastseqout.txt",db_name);
	
	system(command);
		
	number_of_hits += blast_db_results("blastseqout.txt", strlen(query), expectation);
	//std::cout << "Num hits db:" << number_of_hits << std::endl;
	//std::cout << query << std::endl;
	//std::cin.get();
	return(number_of_hits);	
}






int sequence_utils::blast_db_results(const char* filename, int length, double expectation) 
{
  char buffer[256];
  char temp[256];
  //char e_value[256];
  char* token;
  int number_of_hits = 0;
  //char* token2;
  ifstream fin(filename);
  bool HIT, EVA;
  
  // TODO
  while(fin.getline(buffer, 255))
    {
      if(!strncmp(buffer," Score", 6))
	{ 
	  char* eloc = strstr(buffer, "Expect");
	  eloc = strstr(eloc, " = ");
	  token = eloc + 3;
	  //cout << "token: " << token << endl;
	  if(atof(token) < expectation) 
	    {
	      
	      // Get lines 4-6 after this E value from BLAST results
	      HIT = FALSE; EVA = TRUE;
	      for(int i = 0; i < 7; i++) 
		{
		  
		  fin.getline(buffer, 255);
		  
		  if(i == 3)
		    {
		      strcpy(temp, buffer);
		      strtok(temp, " ");
		      token = strtok(NULL, " ");
		      
		      // Set OK if matched sequence begins within 3 bp of query sequence
		      // noting that by BLASTing the reverse complement, the 3' end is 
						// at the beginning of the sequence
		      if(atol(token) < 3)HIT = TRUE;
		      
		      strtok(NULL, " ");
		      token = strtok(NULL, " ");
		      
		      // Set OK if matched sequence ends within 3 bp of query sequence
		      //if(atol(token) > length - 3)OK = 1;
		    }
		  if(HIT && EVA) // We have a hit
		    {
		      number_of_hits++;
		      EVA = FALSE;
		    }
		  
		}
	    } 	 		  	
   	} 
    }
  return(number_of_hits);
}




/*
int sequence_utils::GC_clustering(const char* sequence)
{
	
	return(1);
}
*/



// FASTA: Smith Waterman algorithm implementation 17/2/11
int sequence_utils::Smith_Waterman(const char* query, const char* library_filename, double expectation)
{
	char command[512];
	int number_of_hits;
	char options[100];
	
// Make fasta format file for the query sequence	
	ofstream fout("fasta_input.fa");
	fout << ">fasta query sequence\n" << query << endl;
	fout.close();

// Set options
	sprintf(options,"-q -H");

	sprintf(command, "ssearch35 %s -E %E -O fasta_output.txt fasta_input.fa %s", options, expectation, library_filename);	
	system(command);
	
	number_of_hits = fasta_results_parser("fasta_output.txt");
	
	return(number_of_hits);
}


// FASTA: fasta3 (Pearson & Lipman) 24/2/11

int sequence_utils::fasta3(const char* query, const char* library_filename, double expectation)
{
	char command[512];
	int number_of_hits;
	char options[100];
	
	// Make fasta format file for the query sequence	
	ofstream fout("fasta_input.fa");
	fout << ">fasta query sequence\n" << query << endl;
	fout.close();
	
	// Set options
	sprintf(options,"-q -H");
	
	sprintf(command, "fasta35 %s -E %E -O fasta_output.txt fasta_input.fa %s", options, expectation, library_filename);	
	system(command);
	
	number_of_hits = fasta_results_parser("fasta_output.txt");
	
	return(number_of_hits);
}

int sequence_utils::fasta_results_parser(const char* fasta_output_filename)
{
	char buffer[256];
	char identity[16], location[64];
	int number_of_hits = 0;
	// char* token;
	
	ifstream fin(fasta_output_filename);
	
	while(fin.getline(buffer, 255))
   	{
	 	//cout << buffer << endl;
		
		if(buffer[0] == 0x3E && buffer[1] == 0x3E)
		{
			//cout << buffer << endl;
#ifdef TESTING
			cout << strtok(buffer, "> \n\r\t");
#endif			
			fin.getline(buffer, 255);
			//cout << buffer << endl;
			
#ifdef TESTING
			if(!strncmp (buffer,"rev-comp",8)) cout << ", anti-sense, ";
			   else cout << ", sense, ";
#endif			   
			fin.getline(buffer, 255);
			//cout << buffer << endl;
			strtok(buffer,";%\n\r");
			strcpy(identity, strtok(NULL,";% \n\r"));
			strtok(NULL,"(\n\r");
			strtok(NULL,"(\n\r");
			strcpy(location, strtok(NULL,"()\n\r"));
#ifdef TESTING
			cout << identity << " at " << location << endl;
#endif			
			number_of_hits++;
		}
   	} 
	
#ifdef TESTING
  	cout << number_of_hits << endl;
#endif
	
   	return(number_of_hits);	
}





	
