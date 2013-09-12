
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 genblast.cpp
 
 Created by:	Michael C. Riley and Amanda Clare
 
 Date:			04/05/2012
 
 Copyright (c) 2012, 2013 Aberystwyth University. 
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

/** \file genblast.cpp 
 \brief Wrapper class for BLAST 
 
 */


#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "genblast.h"

genblast::genblast(void)
{
// Defaults
    evalue = 1.0e-6;
	wordsize = 7;
	report = 0;
	
	min_match_length = 9;
	identity = 89;
	tail_loc = 1;	
}

int genblast::blastn_results_parser(const char* filename)
{       
    char buffer[256], a_buffer[256], b_buffer[256], c_buffer[256];
    char* token;	
	int hit_count = 0;
	char query[64];
    long query_length = 0;
	char subject[64];
	
    ifstream fin(filename);
        
    while(fin.getline(buffer, 255))
    {
        if(!strncmp(buffer, "Query=", 6))
        { 
            token = buffer + 7;
			strcpy(query, token);
			
			fin.getline(buffer, 255);			
            
			if(!strncmp(buffer, "Length", 6))
			{
					token = buffer + 7;
					query_length = atol(token);
			}
			else
			{
				fin.getline(buffer, 255);
				if(!strncmp(buffer, "Length", 6))
				{
					token = buffer + 7;
					query_length = atol(token);
				}
				else
					cout << "Error - Query length missed\n";
			}

			while(fin.getline(buffer, 255))
			{
				if(!strncmp(buffer, "Subject", 7))
				{
					token = strtok(buffer, " =(),\n\t\r"); // "Subject"
					token = strtok(NULL, " =(),\n\t\r"); 
					strcpy(subject, token);

					fin.getline(buffer, 255);
										
					if(!strncmp(buffer, "Length", 6))
					{
						token = buffer + 7;
                        tmp.subject_length = atol(token);
					}
					else
					{
						fin.getline(buffer, 255);
						if(!strncmp(buffer, "Length", 6))
						{
							token = buffer + 7;
                            tmp.subject_length = atol(token);
						}
					}					
					break;
				}
			}
			
// Now get hits

			//hit_count = 0;
			while(fin.getline(buffer, 255))
			{
				if(!strncmp(buffer, "*****", 5))
				{
					// No hits
					//cout << "No hits\n";
				}
				if(!strncmp(buffer, " Score", 6))
				{
					fin.getline(buffer, 255);
					token = strtok(buffer, " =(),/\n\t\r"); // "Identities"
					
					token = strtok(NULL, " =(),/\n\t\r"); // identity
					tmp.identity = atoi(token);

					token = strtok(NULL, " =(),/\n\t\r"); // match length
                    tmp.match_length = atoi(token);

					token = strtok(NULL, " =(),/%\n\t\r"); // score value
                    tmp.percent_identity = atof(token);

					// Next line example:- Strand=Plus/Plus, Strand=Plus/Minus
					
                    fin.getline(buffer, 255);
					
					// Blank line
					fin.getline(buffer, 255);
					
					// Next line - Eg: Query  43        GGCCCCTCGAGGCGGTTCGACGA  65
					fin.getline(buffer, 255);
					strcpy(a_buffer, buffer);
					
					token = strtok(buffer, " =(),/\n\t\r"); // "Query"
					token = strtok(NULL, " =(),/\n\t\r"); // start location
                    tmp.query_start_location = atoi(token);
					
					token = strtok(NULL, " =(),/\n\t\r"); // query_hit_sequence
					strcpy(tmp.query_hit_sequence, token);
		
					token = strtok(NULL, " =(),/%\n\t\r"); // end location
                    tmp.query_end_location = atoi(token);
					
					// Eg:                  ||||| |||||| ||||||||||
					fin.getline(buffer, 255);
					strcpy(b_buffer, buffer);
					
					token = strtok(buffer, " =(),/\n\t\r");
					strcpy(tmp.identity_sequence, token);
					
					// Eg: Sbjct  20307025  GGCCCGTCGAGGTGGTTCGACGA  20307003
					fin.getline(buffer, 255);
					strcpy(c_buffer, buffer);
                     
					token = strtok(buffer, " =(),/\n\t\r"); // "Sbjct"
					token = strtok(NULL, " =(),/\n\t\r"); // start location
                    tmp.subject_start_location = atoi(token);
					
					token = strtok(NULL, " =(),/\n\t\r"); // query_hit_sequence
					strcpy(tmp.subject_hit_sequence, token);
					
					token = strtok(NULL, " =(),/%\n\t\r"); // end location
                    tmp.subject_end_location = atoi(token);
                    
                    if(tmp.match_length > min_match_length && tmp.percent_identity > identity && tmp.query_end_location >= query_length - tail_loc)
					{
						hit_count++;
						if(report)
						{
							cout << a_buffer << endl;
							cout << b_buffer << endl;
							cout << c_buffer << endl << endl;
						}
					}
					
					
					if(tmp.match_length == min_match_length && tmp.percent_identity == 100 && tmp.query_end_location >= query_length - tail_loc)
					{
						hit_count++;
						if(report)
						{
							cout << a_buffer << endl;
							cout << b_buffer << endl;
							cout << c_buffer << endl << endl;
						}
					}
					 
					/*if(tmp.identity > 8 && tmp.query_end_location >= query_length - 1)
					{
						hit_count++;
						if(report)
						{
							cout << a_buffer << endl;
							cout << b_buffer << endl;
							cout << c_buffer << endl << endl;
						}
					}
					*/
				}
				
				if(!strncmp(buffer, "Lambda", 6)) // End of record for sequence/chromosome
				{
					fin.getline(buffer, 255);
					// tokenise(" \t") buffer to get Lambda, K and H values
					break;
				}

			}
			//cout << "We have " << hit_count << " for " << subject << endl;
			hits_total = hit_count;
        }                               
    } 
    return(1); 
}


int genblast::blastn_sequence(const char* sequence, const char* template_filename, const char* name, const char* blast_filename)
{
    char query_filename[64];                
    ofstream fseq; // for sequence to be BLASTed
     
// Make sequence into a FASTA file for input to blastn using 'name' for the filename
// and for the FASTA header
    sprintf(query_filename, "%s.fa", name);
    fseq.open(query_filename);
    fseq << ">" << name << endl << sequence << endl;
    fseq.close();
  
	blastn_file(query_filename, template_filename, name, blast_filename);
	
	return(1);
}

int genblast::blastn_file(const char* query_filename, const char* template_filename, const char* name, const char* blast_filename)
{
	char command[1024];
	
    sprintf(command,
	  "%s -subject %s -query %s -evalue %f -word_size %i -out %s.out -ungapped",
		blast_filename,	
	  template_filename, 
	  query_filename,
	  evalue,
	  wordsize,
	  name);
	  
    system(command);

// Parse results 
	char out_filename[32];
	
	sprintf(out_filename, "%s.out", name);
	blastn_results_parser(out_filename);

	return(1); 
}

    

