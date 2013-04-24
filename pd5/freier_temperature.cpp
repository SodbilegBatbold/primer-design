
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 freier_temperature.cpp
 
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

/** \file freier_temperature.cpp 
 \brief Freier annealing temperature calculations/analysis 

 */


#include "freier_temperature.h"


// Frier Eqn (or poss Baldino_89, but paper unavailable)        
double freier_temperature::calculate_temperature(const char* sequence)
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


