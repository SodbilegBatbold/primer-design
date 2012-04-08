/**
   constraints.h
   
 
   Created by Mhr on 16/02/2012.
   Copyright 2012 University of Wales, Aberystwyth. All rights reserved.
 
 */

/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 Michael C. Riley and Amanda Clare
 
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

/** \file sequence_utils.h
 \brief Utilities for sequence analysis
 
 
 */


#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>


#define nucleotide_A 0x41
#define nucleotide_C 0x43
#define nucleotide_G 0x47
#define nucleotide_T 0x54
#define nucleotide_N 0x4E


#define TRUE 1
#define FALSE 0


#define ADENINE 65
#define THYMINE 84
#define GUANINE 71
#define CYTOSINE 67
#define URACIL 85
#define ANYNUCLEOTIDE 78

using namespace std;

//! Primer constraint testing
/**
 * 
 */

class constraints
{
public:
	constraints()
	{
		// Initialise defaults
		homopolymeric_run_length_limit = 5; //< popularly < 4-5 nt (Buck et al, BioTechniques 27:528-536 (September 1999))
	};
	~constraints(){};
	
	/** Detects runs of the same nucleotide longer than homopolymeric_run_length_limit in sequence.
	 Returns true if a long run is detected
	 */
	int homopolymeric_run_detection(const char* sequence);
	
	
private:
	int homopolymeric_run_length_limit; //< Default = 5
	
};

#endif
