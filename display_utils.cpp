/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 display_utils.cpp
 
 Created by Michael Riley and Amanda Clare on 01/06/2011 - mhr/afc"at"aber.ac.uk
 
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

/** \file display_utils.cpp
 \brief Utilities for displaying results
 
 
 */

#include <string.h>
#include "display_utils.h"
#include "primer.h"

#define TRUE 1
#define FALSE 0


display_utils::display_utils() {
}

int display_utils::extract_product(const char* template_seq, primer_data *fwd, primer_data *rev, char * &product) {
  unsigned int begin = fwd->location_5_prime_end;
  unsigned int end   = rev->location_5_prime_end;
  unsigned int len   = 1 + end - begin;

  if (template_seq == 0 || strlen(template_seq) == 0 || strlen(template_seq) < begin || strlen(template_seq) < (begin+len-1)) {
    return 0; 
  } else {
    product = (char *)(malloc (len + 1));
    for (int i=0; i< len; i++){ 
      product[i] = template_seq[begin+i];
    }
    product[len] = '\0';
    return 1;
  }
}




int display_utils::html_colour_sequence(const char* template_seq, primer_data *fwd, primer_data *rev, char * &output) {
    int fwd_5_prime = fwd->location_5_prime_end;
    int fwd_3_prime = fwd_5_prime + strlen(fwd->sequence);
    int rev_5_prime = rev->location_5_prime_end;
    int rev_3_prime = rev_5_prime - strlen(rev->sequence) + 1;
    char bold_start_tag[] = "<font color=\"FF0000\">";
    char bold_end_tag[]   = "</font>";
    int size_of_output = strlen(template_seq) + (2*strlen(bold_start_tag)) + (2*strlen(bold_end_tag)) +1;
    output = (char*) malloc(size_of_output);
    
    char *output_loc = output;
    strncpy(output_loc, template_seq, fwd_5_prime);
    output_loc += fwd_5_prime;

    strncpy(output_loc, bold_start_tag, strlen(bold_start_tag));
    output_loc += strlen(bold_start_tag);

    strncpy(output_loc, template_seq + fwd_5_prime, strlen(fwd->sequence));
    output_loc += strlen(fwd->sequence);

    strncpy(output_loc, bold_end_tag, strlen(bold_end_tag));
    output_loc += strlen(bold_end_tag);

    strncpy(output_loc, template_seq + fwd_3_prime, rev_3_prime - fwd_3_prime);
    output_loc += rev_3_prime - fwd_3_prime;

    strncpy(output_loc, bold_start_tag, strlen(bold_start_tag));
    output_loc += strlen(bold_start_tag);

    strncpy(output_loc, template_seq + rev_3_prime, strlen(rev->sequence));
    output_loc += strlen(rev->sequence);

    strncpy(output_loc, bold_end_tag, strlen(bold_end_tag));
    output_loc += strlen(bold_end_tag);

    strncpy(output_loc, template_seq + rev_5_prime + 1, strlen(template_seq));
    output_loc = output + size_of_output - 1; 

    *output_loc = '\0';
    return strlen(output);
}

