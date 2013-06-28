
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 wallace_temperature.h
 
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

/** \file wallace_temperature.h 
 \brief Annealing temperature calculations/analysis 
 
 A very basic method suitable for short sequences no more than 13 bp. Not recommended,
 but included for research/historical purposes. 
 
 */


#ifndef WALLACE_TEMPERATURE_H
#define WALLACE_TEMPERATURE_H

#include <cstring>
#include <cstdlib>

#include "annealing_temperature.h"

//! Wallace method for primer annealing temperature calculations

class wallace_temperature : public annealing_temperature
{
public:
	
  wallace_temperature() {};
  ~wallace_temperature(){};
		
  /** Non linear equation method based on sequence length and GC content. 
   * Does not compensate for GC clustering. Maps primer3's Santa Lucia algorithm well 
   * for non clustered sequences and handles sequences > 32 bp. Returns the temperature.
   */
  double calculate_temperature(const char* sequence);
	


};

#endif

