/********************************************************************
 
 PD5: a general purpose library for primer design.
 
 cassette_confirmation.cpp
 
 Created by:	Michael C. Riley
 		Amanda Clare
 
 Date:			8/8/2012
 
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

/** In this app we require a forward and reverse primer from within the cassette sequence
 used in precise deletion. These primers are not to be designed to work as a pair, so the analyses 
 required are just for hairpin formation, self dimerisation, annealing temperature and secondary binding
 to the host plasmid sequence and the target organism genome sequence.
 
 (This code is specifically for use with precise deletion and is in progress 8/8/12)
 
 Use: cassette_confirmation.exe <plasmid file name> <template file name> 
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include "../pd5/primer.h"
#include "../pd5/dna_find.h"

int main(int argc, char** argv)
{	
	char buffer[4001];
	char *token;
	int i;
	primer forward_cassette_primer;
	primer reverse_cassette_primer;
	
	reverse_cassette_primer.reverse_primer = TRUE;

	if(argc != 3) {
	  cout << argv[0] << ": Use - exe <plasmid file name> <template file name>\n";
	  exit(0);
	}
	
	// Get plasmid sequence
	char plasmid_sequence[10000];
	
	ifstream plasmid(argv[1]);
	if(!plasmid.is_open())
	{
		cout << "file opening error\n";
		return(0);
	}
	
	plasmid_sequence[0] = 0; // init for strcat
	
	while(plasmid.getline(buffer, 4000))
	{
		if(buffer[0] != 0x3E)
		{
			token = strtok(buffer, " \t\n\r");
			if(token != NULL)strcat(plasmid_sequence, token);
		}
	}	
	plasmid.close();
	
	//cout << plasmid_sequence << endl;
	
	DNAfind genome_template(argv[2]);
	
	if(!genome_template.set_max_mismatches(0)) 
		cout << "Could not set max_mismatches\n";
	
	if(!genome_template.set_tail_length(12)) 
		cout << "Could not set tail length\n";
	
	DNAfind plasmid_template(argv[1]);
	
	if(!plasmid_template.set_max_mismatches(0)) 
		cout << "Could not set max_mismatches\n";
	
	if(!plasmid_template.set_tail_length(8)) 
		cout << "Could not set tail length\n";

	// Set parameters
	// pFS118 cassette 4932-6211, so conf primer region 5332, 5432
	// pCS1966 cassette 2623-4709, so conf primer region 3023, 3123
	
	forward_cassette_primer.set_primer_location_range(3023, 3123);
	forward_cassette_primer.set_primer_length_range(20, 30);
	forward_cassette_primer.homopolymeric_run_check = TRUE;
	//forward_cassette_primer.downstream_search = TRUE;
	forward_cassette_primer.min_Tm = 50.0;
	forward_cassette_primer.optimum_Tm = 55.0;
	forward_cassette_primer.max_Tm = 60.0;
	
	reverse_cassette_primer.set_primer_location_range(3023, 3123);
	reverse_cassette_primer.set_primer_length_range(20, 30);
	reverse_cassette_primer.homopolymeric_run_check = TRUE;
	reverse_cassette_primer.min_Tm = 50.0;
	reverse_cassette_primer.optimum_Tm = 55.0;
	reverse_cassette_primer.max_Tm = 60.0;
	
	
	// Get candidate primers 
	if(!forward_cassette_primer.generate_candidate_primers(plasmid_sequence))
	{
		cout << "No forward primers found\n";
		return(0);
	}
	
	if(!reverse_cassette_primer.generate_candidate_primers(plasmid_sequence))
	{
		cout << "No reverse primers found\n";
		return(0);
	}
	
	// Analyse candidates and display
	forward_cassette_primer.analyse_all_candidates();
	reverse_cassette_primer.analyse_all_candidates();
													  
	for(i = 0; i < forward_cassette_primer.candidates_found; i++)
	{
		//cout << forward_cassette_primer.candidate[i].sequence << endl;
		forward_cassette_primer.candidate[i].binding_A = genome_template.search_for_binding_sites(forward_cassette_primer.candidate[i].sequence);
		forward_cassette_primer.candidate[i].binding_B = plasmid_template.search_for_binding_sites(forward_cassette_primer.candidate[i].sequence);
	}
	
	for(i = 0; i < reverse_cassette_primer.candidates_found; i++)
	{
		//cout << forward_cassette_primer.candidate[i].sequence << endl;
		reverse_cassette_primer.candidate[i].binding_A = genome_template.search_for_binding_sites(reverse_cassette_primer.candidate[i].sequence);
		reverse_cassette_primer.candidate[i].binding_B = plasmid_template.search_for_binding_sites(reverse_cassette_primer.candidate[i].sequence);
	}
	
	forward_cassette_primer.show_all_single_candidates();
	reverse_cassette_primer.show_all_single_candidates();

	// Sort individual candidates and display
	/*forward_cassette_primer.set_priorities("BINDING_B, BINDING_A, HAIRPIN, SELF_DIMER, TEMPERATURE");
	reverse_cassette_primer.set_priorities("BINDING_B, BINDING_A, HAIRPIN, SELF_DIMER, TEMPERATURE");

	forward_cassette_primer.sort_candidates();
	reverse_cassette_primer.sort_candidates();
	
	forward_cassette_primer.show_all_single_candidates();
	reverse_cassette_primer.show_all_single_candidates();*/

	return(1);
}	
