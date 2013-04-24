
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 primer3_Tm_temperature.h
 
 Created by:	Michael C. Riley and Amanda Clare
 
 Date:			19/04/2013
 
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

/** \file primer3_Tm_temperature.h 
 \brief Annealing temperature calculations/analysis 
 
 A very basic method suitable for short sequences no more than 13 bp. Not recommended,
 but included for research/historical purposes. 
 
 */


#ifndef PRIMER3_TM_TEMPERATURE_H
#define PRIMER3_TM_TEMPERATURE_H

#include <cstring>
#include <cstdlib>
#include "sequence_utils.h" 
#include "annealing_temperature.h"
#include "oligotm.h"

//! Primer3_Tm method for primer annealing temperature calculations

class primer3_Tm_temperature : public annealing_temperature
{
public:
	
  primer3_Tm_temperature();
  ~primer3_Tm_temperature(){};
		

  /** Calls Primer3's annealing temperature calculation, which is based on thermodynamic methods.
      Note that the thermodynamic algorithm is not used on sequences > 32 bp. */ 	
  double calculate_temperature(const char* sequence);

	
  /**
     \defgroup Primer3 Primer3's annealing temperature methods
     @{
  */
  

  double oligo_Tm(const char* sequence);
  
  /// primer3 parameters
  typedef enum primer3_parameter 
  {
    DNA_CONC		= 0,
    SALT_CONC		= 1,
    DIVALENT_CONC   = 2,
    DNTP_CONC		= 3,
    NN_MAX_LEN		= 4,
  } primer3_parameter;
  
  /** For setting any of the 5 numerical values in Primer3's annealing temperature algorithm.
      Eg. set_primer3_parameter(DNA_CONC, 70); */
  int set_primer3_parameter(primer3_parameter parameter, double value);
  
  /** For setting either Santa Lucia or Breslauer method. \sa oligoTm.h (for typedef) */
  int set_primer3_method(tm_method_type method);
  
  /** For setting Schildkraut, Santa Lucia or Owczarzy salt correction method. \sa oligoTm.h (for typedef) */
  int set_primer3_salt_correction_method(salt_correction_type method);	
  /** @} */
	

private:

  /// primer3_Tm attributes
	double dna_conc;
	double salt_conc;
	double divalent_conc;
	double dntp_conc;
	int    nn_max_len;
	tm_method_type tm_method;
	salt_correction_type salt_corrections;

};

#endif

