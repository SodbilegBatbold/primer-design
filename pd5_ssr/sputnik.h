#ifndef SPUTNIK_H
#define SPUTNIK_H

/* trivial defs */
#ifndef True
#define True 1
#endif
#ifndef False
#define False 0
#endif



typedef int Boolean;

/* size of buffer for reads. */ 
#define BUF_SIZE 1024*10   /* 10K */
/* max size of description line (begins with ">") */
#define MAX_DESCRIPTION_LEN 1024
/* max sequence length */
#define MAX_SEQUENCE_LEN 1024*800 /* 800K */
/* max number of sequence chars dumped to line */
#define MAX_OUT_LINE_CHARS 60

/* for debugging only */
#define MAX_ERRCODES 1024



using namespace std;

class sputnik
{
public:	
        sputnik();



	//char *repeatName[MAX_UNIT_LENGTH+1];
	static char *repeatName[6];
	char readBuf[BUF_SIZE];
	Boolean endOfFile;
	int curBufLen;
	int curBufPos;
	int fd;
	Boolean havePutBack;
	char putBack;
	
	 
	/* search params and definitions */
	int Min_Unit_Length; /* start search with dinucleotide repeats */
	/* will search for di, tri, tetra ... <n>nucleotide repeats up to 
	   this value for n */
	int  Max_Unit_Length;  /* up to and including pentanucleotides */
	/* this is the point score for each exact match */
	int  Exact_Match_Points;
	/* this is the point score for a mismatch, insertion or deletion */
	int  Error_Match_Points;
	/* this is the minimum score required to be considered a match */
	int Match_Min_Score;
	/* this is the low score at which we stop trying */
	int  Match_Fail_Score;
	/* this is the max recursion depth we try to recover errors */
	int  Max_Recursion;


	Boolean findingPrimers;



	/* struct for indiv sequence in a file */
	typedef struct ss
	{
	  char descStr[MAX_DESCRIPTION_LEN];
	  char seqStr[MAX_SEQUENCE_LEN];
	  unsigned int seqLen;
	} SeqStruct, *SeqStructPtr;
	
	
	/*
	 * this structure describes the current state of a comparison.
	 * it gets passed down to recursive calls of the find repeat
	 * call so it can know when to bail out of an unsuccessful
	 * search, or return the size/state of a successful hit, etc.
	 */
	typedef struct ms
	{
	  int curPos;         /* putative pattern starts here */
	  int testPos;        /* start testing here */
	  int testLen;        /* di, tri, tetra, etc. */
	  int testCtr;        /* # chars in testLen already tested. mod counter */
	  int curScore;       /* current score */
	  int missense;       /* keep track of ins, del, err */
	  int insertions;     
	  int deletions;
	  int depth;          /* how deep is recursion for this match */
	  char errCodes[MAX_ERRCODES];
	} MatchStruct, *MatchStructPtr;
	/* a utility macro to copy one testStruct to another */
	
#define copyMSPtr(dest,source) memcpy((char *)dest,(char *)source,sizeof(MatchStruct))
/* a utility macro to increment the modular testCtr */
#define bumpTestCtr(msp) (msp)->testCtr++; if ((msp)->testCtr==(msp)->testLen) (msp)->testCtr=0;


	void fillBuf();
	Boolean getChar(char *achar);
	void putCharBack(char c);
	void openFile(char *fn);
	void initBuffer();
	void addCharToLine(char c, char *line, int *lineLen);
	Boolean getNonBlankLine(char *line);
	SeqStructPtr getSeq(char *fname);
	void dumpMatch(SeqStructPtr seqP, 
		       MatchStructPtr matchP,
		       Boolean anyMatchThisSeq);
	Boolean testForNRepeat(SeqStructPtr seqP,
			       MatchStructPtr matchP);
	Boolean ignoreSeq(SeqStructPtr seqP,
			  MatchStructPtr matchP);
	Boolean findRepeats(SeqStructPtr seqP);
};

#endif
