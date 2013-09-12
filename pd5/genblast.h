
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 genblast.h
 
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

/** \file genblast.h 
 \brief Wrapper class for BLAST 
 
 */

using namespace std;
class BLAST_results_data
{
public:
	char query[64];
	long query_length;
	char subject[64];
	long subject_length;
    //char chromosome[16];
	double e_result;
	double bit_score;
	int score;
	int identity;
	int match_length;
	double percent_identity;
	long query_start_location;
	char query_hit_sequence[1024];
	long query_end_location;
	char identity_sequence[1024];
	long subject_start_location;
	char subject_hit_sequence[1024];
	long subject_end_location;
};

//! Wrapper for BLAST use and results analysis
/**
 Class contains methods for running BLAST on a sequence string or a sequence from file and a further method for parsing the BLAST results file.
 */

class genblast
{
public:
	genblast();
	~genblast(){};
	
	int blastn_results_parser(const char* filename);
	int blastn_sequence(const char* sequence, const char* template_filename, const char* name, const char* blast_filename);
	int blastn_file(const char* query_filename, const char* template_filename, const char* name, const char* blast_filename);
	
	// Blast parameters	
	double evalue; ///< E value used by BLAST
    int wordsize;  ///< Word size used by BLAST
	
	// Blast results
	
	BLAST_results_data tmp;
	
	// BLAST sequence filters
	int min_match_length; ///< BLAST results filter: minimum length of BLAST hit sequence.
	int identity; ///< BLAST results filter: Percentage of matching bases.
	int tail_loc; ///< BLAST results filter: how far from 3' tail can sequence be.
	int hits_total; ///< Number of BLAST hits after filtering
	
	// Options
	bool report; ///< Set to true to send BLAST filtered sequence matches to stdout. (Default false)
	
private:
	
};

