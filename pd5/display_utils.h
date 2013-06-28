/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 display_utils.h
 
 Created by Michael Riley and Amanda Clare on 01/06/2011 
 
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

/** \file display_utils.h
 \brief Utilities for displaying results
 
 
 */

#ifndef DISPLAY_UTILS_H
#define DISPLAY_UTILS_H

//#include <fstream>
//#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "primer.h"
#include "primer_pair.h"


#define TRUE 1
#define FALSE 0


using namespace std;

//! Utilities for displaying results and data
/**
 Includes a method for displaying results in html
 */

class display_utils
{
 public:
  display_utils();
  ~display_utils() {};
  
  /**
     \brief Populates \a product with the product that would be created from \a template_seq using the \a fwd and \a rev primers. 

     Returns TRUE if successful and FALSE otherwise.  
   */
  int extract_product(const char* template_seq, primer_data *fwd, primer_data *rev, char * &product);
  
  /**
     \brief Populates \a output with the \a template_seq where the \a primers are highlighted in red. 

     Returns the length of the output.
   */
  int html_colour_sequence(const char* template_seq, primer_pair_data *primers, char * &output);
  /**
     \brief Populates \a output with the \a template_seq where the \a fwd and \a rev primers are highlighted in red. 

     Returns the length of the output.
   */
  int html_colour_sequence(const char* template_seq, primer_data *fwd, primer_data *rev, char * &output);

};



#endif
