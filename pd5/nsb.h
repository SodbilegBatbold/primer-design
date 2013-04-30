/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 Created by:	Michael C. Riley and Amanda Clare
 
 Date:			27/03/2011.
 
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

/** \file nsb.h
 \brief Abstract class for secondary binding analysis methods
 
 
 */

#ifndef NSB_H
#define NSB_H


#include <fstream>

using namespace std;


//! Secondary binding data
/**
 Data on the location and nature of secondary binding occurrences
 */


class location_data
{
public:
	long location;
	bool sense_strand;
	char matching_sequence[21]; 
	// int matching_sequence_score;
};

//! Abstract secondary binding class
/**
 * Abstract class for detecting potential secondary primer binding on the template and 
 * potential secondary products.
 */
 class nsb {

 public:

  //! Constructor
  /** Note that the instantiation of an nsb object requires a file name. Opens the file. 
   */
  nsb(const char* filename);

  //! Destructor. Closes the file.
  virtual ~nsb() = 0;
  

  location_data forward_primer_match_locations[1012]; ///< Calculated results stored here
  location_data reverse_primer_match_locations[1012]; ///< Calculated results stored here

  /** \brief Secondary binding site prediction */
  /** Finds all potential binding site locations. Returns number of binding sites found, locations 
   *  of all sites found are in forward_primer_match_locations and reverse_primer_match_locations
   */
  virtual int search_for_binding_sites(const char* sequence) = 0;

  /** \brief Secondary product prediction */
  /** Finds all potential products. Returns number of products found. 
   */
  virtual int search_for_pcr_products(const char *forward_sequence, const char *reverse_sequence) = 0;

 protected:
	ifstream fin;	

 private:
	nsb(); // we're not allowing derived classes to use this constructor
};


#endif
