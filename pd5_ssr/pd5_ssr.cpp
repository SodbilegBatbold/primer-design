/* #define DEBUG_SPUTNIK 1 */


/*
  find repeats in fasta format seq file 
  allows for indels, returns score.  

  beta version.  caveat emptor.  

  chrisa  29-Jul-94

  chris abajian
  University of Washington
  Dept. of Molecular Biotechnology  FJ-20
  Fluke Hall, Mason Road
  Seattle WA 98195
*/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include "sputnik_ssr.h"
#include "sputnik.h"



int main(int argc, char* argv[])
{
   sputnik::SeqStructPtr seqP;
   int count;
   sputnik sp;

   if (argc != 2)
     {
        fprintf(stderr,"Usage: %s <fasta format sequence file name>\n", argv[0]);
        exit(1);
     }

   sp.openFile(argv[1]);

   sp.initBuffer();

   count = 0;
   while (! sp.endOfFile) 
	   if ((seqP = sp.getSeq(argv[1])))
       {
#ifdef DEBUG_SPUTNIK
	 fprintf(stderr,"processing sequence %d\n", count++);
	 dumpSeq(seqP);
#endif
	 sp.findRepeats(seqP);
          free((void *)seqP);
       }

   return 0;
}
