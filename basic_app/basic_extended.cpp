/********************************************************************
 
 HELYGEN: a general purpose library for primer design.
 
 helygen_basic.cpp
 
 Created by:	Michael C. Riley
				Amanda Clare
 
 Date:			05/11/2011
 
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


#include "../pd5/primer_pair.h"

int main(int argc, char** argv)
{	
	primer_pair pcr1;

	if(argc != 2) {
	  cout << argv[0] << ": Please provide a DNA string as argument" << endl;
	  exit(1);
	}
	
	// Set parameters for both primers
	pcr1.set_target_location(200, 500); // No default
	pcr1.set_flank_lengths(200, 200); // upstream, downstream. Default 100, 100
	pcr1.set_primer_length_range(20, 20); // shortest acceptable primer length, longest. Default 18, 25
	
	// Set parameters for primers individually
	pcr1.forward.
	
	pcr1.reverse.
	
	// Get candidate primers (argv[1] = dna template)
	pcr1.generate_candidates(argv[1]);  
	
	// Analyse candidates and display
	pcr1.candidate_analysis();
	pcr1.show_individual_candidates();
	
	// Sort individual candidates and display
	pcr1.sort_individual_candidates("HAIRPIN, SELF_DIMER, TEMPERATURE");
	pcr1.show_individual_candidates();
	
	// Select and sort primer pairs
	pcr1.sort_pair_candidates("TM_DIFF, F_DIMER, R_DIMER, MOO_SORT");
	
	// Display best 6 candidate pairs
	pcr1.show_best_pair_candidates(6); 
	
	return(1);
}	
