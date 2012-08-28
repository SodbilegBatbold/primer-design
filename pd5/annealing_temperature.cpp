
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 annealing_temperature.cpp
 
 Created by:	Michael C. Riley and Amanda Clare
 
 Date:			19/02/2011
 
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

/** \file annealing_temperature.cpp 
 \brief Annealing temperature calculations/analysis 
 
 
 */


#include "annealing_temperature.h"


annealing_temperature::annealing_temperature(void)
{
	// Primer3 defaults
	dna_conc = 50;
	salt_conc = 50;
	divalent_conc = 0;
	dntp_conc = 0;
	nn_max_len = 35;
	tm_method = santalucia_auto; // breslauer_auto OR santalucia_auto
	salt_corrections = santalucia; //Santa Lucia
	
}


// Frier Eqn (or poss Baldino_89, but paper unavailable)	

double annealing_temperature::Freier_method(const char* sequence)
{
	double Tm = 0; // Tm is in Celsius
	int GC_content = 0;
	int length = strlen(sequence);
	
	// GC content
	for(int i = 0; i < length; i++) if(sequence[i] == 71 || sequence[i] == 67) GC_content++;
	
	// Equation
	Tm = 59.9 + (41.0 * GC_content/length) - (675.0/length);
	
	return(Tm);
}

// Marmur (also referred to as the Wallace eqn from Wallace_79)
// Restricted to short sequences (< 13 bp)

double annealing_temperature::Marmur_method(const char* sequence)
{
	double Tm = 0; // Tm is in Celsius
	int GC_content = 0;
	int length = strlen(sequence);
	
	// GC content
	for(int i = 0; i < length; i++) if(sequence[i] == 71 || sequence[i] == 67) GC_content++;
	
	// Equation
	Tm = (4 * GC_content) + (2 * (length - GC_content));
	
	return(Tm);
}

// Wallace formula

double annealing_temperature::Wallace_method(const char* sequence)
{
	double Tm = 0; // Tm is in Celsius
	int GC_content = 0;
	int length = strlen(sequence);
	
	// GC content
	for(int i = 0; i < length; i++) if(sequence[i] == 71 || sequence[i] == 67) GC_content++;
	
	// Equation
	Tm = 64.9 + (41.0 * (GC_content - 16.4)/length);
	
	return(Tm);
}

// Santa Lucia - Just a placeholder for now, use primer3_Tm instead. 13/7/11

double annealing_temperature::Santa_Lucia_method(const char* primer_sequence, const char* template_sequence)
{
	double Tm = 0; // Tm is in Celsius

	
	return(Tm);
}

int annealing_temperature::set_primer3_parameter(primer3_parameter parameter, double value)
{
	switch(parameter)
	{
		case DNA_CONC:
			dna_conc = value;
			break;
		case SALT_CONC:
			salt_conc = value;
			break;
		case DIVALENT_CONC:
			divalent_conc = value;
			break;
		case DNTP_CONC:
			dntp_conc = value;
			break;
		case NN_MAX_LEN:
			nn_max_len = (int)value;
			break;
		default:
			return(FALSE);
	}
	return(TRUE);
}
						   
						   
int annealing_temperature::set_primer3_method(tm_method_type method)
{
	tm_method = method;
	return(TRUE);
}

int annealing_temperature::set_primer3_salt_correction_method(salt_correction_type method)
{
	salt_corrections = method;
	return(TRUE);
}

double annealing_temperature::primer3_Tm(const char* primer_sequence)
{
	double Tm;
	
// Function call to Primer3's seqtm() in oligotm.cpp
	Tm = seqtm(primer_sequence,
			   dna_conc,
			   salt_conc,
			   divalent_conc,
			   dntp_conc,
			   nn_max_len,
			   tm_method,
			   salt_corrections);

	return(Tm);

}

double annealing_temperature::oligo_Tm(const char* primer_sequence)
{
	double Tm;
	
	// Function call to Primer3's seqtm() in oligotm.cpp
	Tm = oligotm(primer_sequence,
			   dna_conc,
			   salt_conc,
			   divalent_conc,
			   dntp_conc,
			   tm_method,
			   salt_corrections);
	
	return(Tm);
	
}
