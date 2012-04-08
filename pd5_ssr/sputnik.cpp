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
#include "sputnik_helygen.h"
#include "sputnik.h"

/* #define DEBUG_SPUTNIK */


char *sputnik::repeatName[6] = {
  "***ERROR***",    /* bad programmer!  no latte! */
  "single nucleotide",
  "dinucleotide",
  "trinucleotide",
  "tetranucleotide",
  "pentanucleotide"
};


/*
 ************************************************************
 * these routines are used to read and parse the fasta format
 * sequence file
 ************************************************************
 */

sputnik::sputnik() {
  /* search params and definitions */
  Min_Unit_Length = 2; /* start search with dinucleotide repeats */
  /* will search for di, tri, tetra ... <n>nucleotide repeats up to 
     this value for n */
  Max_Unit_Length = 4;  /* up to and including pentanucleotides */
  /* this is the point score for each exact match */
  Exact_Match_Points = 1;
  /* this is the point score for a mismatch, insertion or deletion */
  Error_Match_Points = -10;
  /* this is the minimum score required to be considered a match */
  Match_Min_Score = 18;
  /* this is the low score at which we stop trying */
  Match_Fail_Score = -1;
  /* this is the max recursion depth we try to recover errors */
  Max_Recursion = 5;

  findingPrimers = True;
}

void sputnik::fillBuf()
{
   size_t result;

   result = read(fd, (void *)readBuf, BUF_SIZE);
   if (result == -1)
     {
        fprintf(stderr,"error reading file! errno = %d\n",errno);
        exit(1);
     }
   else if (result == 0)
     endOfFile = True;
   else
     {
        curBufLen = result;
        curBufPos = 0;
     }
}  /* readBuf */


/* returns True on success */
Boolean sputnik::getChar(char *achar)
{
   if (havePutBack)
     {
        *achar = putBack;
        havePutBack = False;
        return(True);
     }

   if (curBufPos == curBufLen)
     fillBuf();

   if (endOfFile)
     return (False);

   *achar = readBuf[curBufPos++];
   return (True);
}


void sputnik::putCharBack(char c)
{
   havePutBack = True;
   putBack = c;
}


void sputnik::openFile(char *fn)
{
   /* open the specified file */
   fd = open(fn, O_RDONLY);
   if (fd == -1)
     {
        fprintf(stderr,"unable to open file %s\n", fn);
        exit(1);
     }
}

/* should call this once for each file read */
void sputnik::initBuffer()
{
   /* initialize length and pointer */
   curBufPos = 0;
   curBufLen = 0;
   havePutBack = False;
   endOfFile = False;
}

void sputnik::addCharToLine(char c, char *line, int *lineLen)
{
   if (*lineLen < MAX_DESCRIPTION_LEN)
     line[(*lineLen)++] = c;
   else {
     fprintf(stderr,"warning: description line truncated\n");
     /* fprintf(stderr,"%s\n",line); */
   }
}


/*
 *********************************************************************
 * these routines are (more) specific to reading the fasta file format
 *********************************************************************
 */


/* 
 * pick up a non-blank line from the file, presumably description.
 * truncates all leading blanks and/or blank lines 
 */
Boolean sputnik::getNonBlankLine(char *line)
{
   Boolean stop, nonBlank;
   char c;
   int lineLen;

   lineLen = 0;
   stop = False;
   nonBlank = False;  /* will be set by any non whitespace char */
   while ((! endOfFile) && (! stop))
     if (getChar(&c))
       if (c == '\n')
         stop = nonBlank; /* stop if have anything. don't save eol char. */
       else {
         if (nonBlank)
           /* add it to line no matter what */
           addCharToLine(c,line,&lineLen);
         else if ((c != ' ') && (c != '\t'))
           {
              /* only non whitespace will start the line */
              nonBlank = True;
              addCharToLine(c,line,&lineLen);
           }
       }
   /*printf("line: %s\n",line);*/
   return True;
}


/* load the sequence struct with comment line and bases */
sputnik::SeqStructPtr sputnik::getSeq(char *fname)
{
   SeqStructPtr newSeqP;
   Boolean endOfSeq;
   char c;

   if (endOfFile) return ((SeqStructPtr )0);   /* bombproofing */
   /* malloc a new seq */
   if (! (newSeqP = (SeqStructPtr )malloc(sizeof(SeqStruct)) ) )
     {
        fprintf(stderr,"unable to malloc() memory for sequence.\n");
        exit(1);
     }
   /* clear mem */
   memset( (void *)newSeqP, '\0', sizeof(SeqStruct));

   /* pick up description line */
   if (! getNonBlankLine(newSeqP->descStr) )
     {
        free(newSeqP);
        return ((SeqStructPtr )0);   
     }

   /* did it start correctly ? */
   if (newSeqP->descStr[0] != '>')
     {
        fprintf(stderr,"format error in input file:  missing '>'\n");
	/* fprintf(stderr,"%s\n",newSeqP->descStr); */
        exit(1);
     }

   endOfSeq = False;
   while ((!endOfFile) && (!endOfSeq))
     {
        if (getChar(&c))
          {
             if (c == '>')
               {
                  /* hit new sequence */
                  endOfSeq = True;
                  putCharBack(c);
               }
             else if (((c >= 'A') && (c <= 'Z')) ||
                      ((c >= 'a') && (c <= 'z')))/* bogus test, chris */
               /* have nucleotide */
               newSeqP->seqStr[newSeqP->seqLen++] = toupper(c);
             else if ((c != '\n') && (c != ' ') && (c != '\t'))
               {
                  /* wierd shit in file.  bail. */
                  fprintf(stderr,"bad char in sequence, file %s: %c\n",fname,c);
                  exit(1);
               } 
          }
     }     

   if (! newSeqP->seqLen)
     {
        fprintf(stderr,"? Null sequence encountered in file %s (ignored)\n",fname);
        fprintf(stderr,"  %s\n", newSeqP->descStr);
        free(newSeqP);
        return ((SeqStructPtr )0);   
     }

   return(newSeqP);
}  /* getSeq */


/* for debugging.  dump entire seq to stdout. */
#ifdef DEBUG_SPUTNIK
void dumpSeq(SeqStructPtr seqP)
{
   int i, charsOnLine;

   fprintf(stdout,"%s\n", seqP->descStr);
   fprintf(stdout,"Sequence (length = %d):\n", seqP->seqLen);
   i = 0;
   charsOnLine = 0;
   while (i < seqP->seqLen)
     {
        if (charsOnLine == MAX_OUT_LINE_CHARS)
          {
             fprintf(stdout,"\n");
             charsOnLine = 1;
          }
        else
          charsOnLine++;
        fprintf(stdout,"%c", seqP->seqStr[i++]);
     } 
   fprintf(stdout,"\n");
} /* dumpSeq */
#endif /* DEBUG_SPUTNIK */






/* dump the matched seq & stats to stdout */
void sputnik::dumpMatch(SeqStructPtr seqP, 
               MatchStructPtr matchP,
               Boolean anyMatchThisSeq)
{
   int i, charsOnLine;

   //if (! anyMatchThisSeq)
   fprintf(stdout,"%s\n", seqP->descStr);

   fprintf(stdout,"%s %d : %d -- length %d score %d\n",
           repeatName[matchP->testLen], 
           matchP->curPos+1,
           matchP->testPos,
           matchP->testPos - matchP->curPos,
           matchP->curScore);

#ifdef DEBUG_SPUTNIK
   fprintf(stdout,"mis = %d, del = %d, ins = %d\n", 
           matchP->missense,
           matchP->deletions,
           matchP->insertions);
#endif

   i = matchP->curPos;
   charsOnLine = 0;
   while (i < matchP->testPos)
     {
        if (charsOnLine == MAX_OUT_LINE_CHARS)
          {
             fprintf(stdout,"\n");
             charsOnLine = 1;
          }
        else
          charsOnLine++;
        fprintf(stdout,"%c", seqP->seqStr[i++]);
     } 
   fprintf(stdout,"\n");

#ifdef DEBUG_SPUTNIK
   i = 0;
   charsOnLine = 0;
   while (i < (matchP->testPos - matchP->curPos))
     {
        if (charsOnLine == MAX_OUT_LINE_CHARS)
          {
             fprintf(stdout,"\n");
             charsOnLine = 1;
          }
        else
          charsOnLine++;
        if (matchP->errCodes[i] == '\0')
          fprintf(stdout," ");
        else
          fprintf(stdout,"%c", matchP->errCodes[i]);
        i++;
     } 
   fprintf(stdout,"\n");
#endif


}  /* dumpMatch */


Boolean sputnik::testForNRepeat(SeqStructPtr seqP,
                       MatchStructPtr matchP)
{
   MatchStruct curMatch, recover, bestSoFar, bestOfABadLot;

   /* save matchP in case we fail altogether. */
   copyMSPtr(&curMatch, matchP);
   /* keep track of the best score and return that if over thresh. */
   copyMSPtr(&bestSoFar, matchP);

   while ( (curMatch.testPos < seqP->seqLen)           /* anything to test */
          && (curMatch.curScore > Match_Fail_Score) )  /* above fail threshold */
     {
        /* test a base */
        if (seqP->seqStr[curMatch.curPos+curMatch.testCtr] 
            == seqP->seqStr[curMatch.testPos])
          {
             /* we matched.  this is easy. */
             curMatch.curScore += Exact_Match_Points;  /* score your points */
             curMatch.testPos++; /* advance the downstream test position */
             bumpTestCtr(&curMatch); /* advance pos in the (presumed) repeating seq */
          }
        else if (seqP->seqStr[curMatch.testPos] == 'N')
          {
             /* don't call it wrong, but no credit either */
             curMatch.testPos++; /* advance the downstream test position */
             bumpTestCtr(&curMatch); /* advance pos in the (presumed) repeating seq */
          }
        else
          {
             /* no match.  take the score penalty, but keep going (maybe). */
             curMatch.curScore += Error_Match_Points;
             curMatch.testPos++; /* advance the downstream test position */
             bumpTestCtr(&curMatch);  /* advance pos in seq */
             /* is the score too bad to continue, or are we
                already too deep? */
             if ( (curMatch.curScore > Match_Fail_Score)
                  && (curMatch.depth < Max_Recursion) )
               {
                  /* try simple missense */
                  copyMSPtr(&recover,&curMatch);
                  if ((recover.testPos - recover.curPos) < MAX_ERRCODES)
                    recover.errCodes[recover.testPos - recover.curPos -1] = 'M';
                  recover.missense++;
                  recover.depth++;
                  (void )testForNRepeat(seqP,&recover);
                  copyMSPtr(&bestOfABadLot,&recover);

                  /* try deletion */
                  copyMSPtr(&recover,&curMatch);
                  if ((recover.testPos - recover.curPos) < MAX_ERRCODES)
                    recover.errCodes[recover.testPos - recover.curPos -1] = 'D';
                  recover.testPos--; /* DON'T advance downstream */
                  recover.deletions++;
                  recover.depth++;
                  (void )testForNRepeat(seqP,&recover);
                  if (recover.curScore > bestOfABadLot.curScore)
                    copyMSPtr(&bestOfABadLot,&recover);

                  /* try insertion */
                  copyMSPtr(&recover,&curMatch);
                  if ((recover.testPos - recover.curPos) < MAX_ERRCODES)
                    recover.errCodes[recover.testPos - recover.curPos -1] = 'I';
                  /* RETEST for this base in the repeating seq */
                  if (recover.testCtr == 0)
                    recover.testCtr = recover.testLen - 1;
                  else
                    recover.testCtr--;
                  recover.insertions++;
                  recover.depth++;
                  (void )testForNRepeat(seqP,&recover);
                  if (recover.curScore > bestOfABadLot.curScore)
                    copyMSPtr(&bestOfABadLot,&recover);

                  /* take the best of a bad lot */
                  bestOfABadLot.depth--;  /* dec recursion counter */ 
                  copyMSPtr(&curMatch, &bestOfABadLot);
               }  /* it was worth carrying on */
          }  /* no match, found best of bad lot */

        /* whatever happened, the best we could do is now in matchP */
        if (curMatch.curScore > bestSoFar.curScore)
          copyMSPtr(&bestSoFar, &curMatch);

     }  /* while loop to test a single base */

   /* for whatever reason, we've stopped searching for more of this
      putative repeat.  if there were any matches that passed
      the global threshold, return the best of them. note that this
      has the effect of NOT advancing the pointer(s) if nothing
      rang the bell.  remember that we will test the same position
      for ntide repeats of several different lengths. */
   if (bestSoFar.curScore > Match_Min_Score)
     {
        copyMSPtr(matchP, &bestSoFar);
        return(True);
     }
   return(False);   /* the whole thing was a waste of time */
}  /* testForNRepeat */


/* 
 * returns True if the sequence we want to look for repeats of is
 *
 * a) all the same base (i.e. 'AAA' or 'GG').  This filters out
 *    single nucleotide repeats
 *
 * b) conains 'N'.  we search against these, but don't use them
 *    as wildcards.
 */
Boolean sputnik::ignoreSeq(SeqStructPtr seqP,
                  MatchStructPtr matchP)
{
   int i;

   /* firstly, never search for any pattern that contains N */
   for (i = 0; i < matchP->testLen; i++)
     if (seqP->seqStr[matchP->curPos+i] == 'N')
       return(True);

   /* now test for mononucleotide repeat.  other tests may get 
      added, in which case this one will beed to be changed. */
   for (i = 1; i < matchP->testLen; i++)
      if (seqP->seqStr[matchP->curPos] != seqP->seqStr[matchP->curPos+i])
        return(False);  /* they're not all the same */
   return (True);  /* they ARE all same */
}


Boolean sputnik::findRepeats(SeqStructPtr seqP)
{
   int curPos;
   Boolean anyMatchThisSeq, matchAtThisPos;
   MatchStruct match;

   memset( (char *)&match, 0, sizeof(MatchStruct) );  /* clear match struct */

   anyMatchThisSeq = False; /* avoid dumping description more than once. */
   /* loop on all positions in the sequence.  note that a match
      will advance curPos past all matching chars to the first
      unmatched char. */
   while ( match.curPos <= seqP->seqLen)
     {
        /* now loop on all the different lengths of repeats we're
           looking for (i.e. di, tri, tetra nucleotides.  if we
           find a match at a shorter repeat length, forego testing
           for longer lengths. */
        match.testLen = Min_Unit_Length;
        matchAtThisPos = False;
        while ((match.testLen <= Max_Unit_Length) && (!matchAtThisPos))
          {
             /* initialize the state of the match */
             match.curScore = 0;  /* no points yet */
             match.testCtr = 0; /* no chars tested yet */
             match.testPos = match.curPos + match.testLen;
             match.insertions = 0;
             match.deletions = 0;
             match.missense = 0;
             /* there are some things we don't want to test for */
             if (! ignoreSeq(seqP,&match))
               matchAtThisPos = testForNRepeat(seqP, &match);
             else
               matchAtThisPos = False;
             if (! matchAtThisPos) match.testLen++;
          }

        if (matchAtThisPos)
          {
	    if(findingPrimers) {
	      sputnik_helygen primer_factory;
	      int good_primers = primer_factory.find_primers(seqP->seqStr, match.curPos+1, match.testPos, match.testLen);
	      if(good_primers) {
		
		dumpMatch(seqP,&match,anyMatchThisSeq);
	      }
	    } else {
	      //dumpMatch(seqP,&match,anyMatchThisSeq);
	      return True;
	    }
	    anyMatchThisSeq |= matchAtThisPos;
	    match.curPos = match.testPos;
          }
        else
          match.curPos++;  /* no, so advance to next base. */
     }
   return anyMatchThisSeq;
}


