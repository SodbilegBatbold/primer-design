
/********************************************************************
 
 PD5: a general purpose library for primer design app development.
 
 primer_pair_data.cpp
 
 Created by:	Michael C. Riley
				Amanda Clare
 
 Date:			05/04/2011
 
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

/** \file primer_pair_data.cpp
 \brief Data pertaining to a pair of primers

 */


#include "primer_pair_data.h"

int primer_pair_data::pair_dimerisation(void)
{
	dimerisation pair_dimer;
	
	if(forward_sequence && reverse_sequence)
	{
		pair_dimer.pair_dimer(forward_sequence, reverse_sequence);
		forward_pair_dimer_score = pair_dimer.forward_dimer_score;
		reverse_pair_dimer_score = pair_dimer.reverse_dimer_score;
	}
	else
	{
		return(FALSE);
	}	
	
	return(TRUE);
}

int primer_pair_data::pcr_products(nsb &mynsb)
{
	if(forward_sequence && reverse_sequence)
	{
		//mynsb.search_nsb_product(forward_sequence, reverse_sequence);
		number_of_pcr_products = mynsb.search_for_pcr_products(forward_sequence, reverse_sequence);
	}
	else
	{
		return(FALSE);
	}	
	
	return(TRUE);
	
}
