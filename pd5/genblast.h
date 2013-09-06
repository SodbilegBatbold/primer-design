/*
 *  genblast.h
 *  
 *
 *  Created by Mhr on 04/05/2012.
 *  Copyright 2012 University of Wales, Aberystwyth. All rights reserved.
 *
 */

using namespace std;
class results_data
{
public:
	char query[64];
	long query_length;
	char subject[64];
	long subject_length;
    //char chromosome[16];
	double e_result;
	double bit_score;
	int score;
	int identity;
	int match_length;
	double percent_identity;
	long query_start_location;
	char query_hit_sequence[1024];
	long query_end_location;
	char identity_sequence[1024];
	long subject_start_location;
	char subject_hit_sequence[1024];
	long subject_end_location;
};

class genblast
{
public:
	genblast();
	~genblast(){};
	
	int blastn_results_parser(const char* filename);
	int blastn_sequence(const char* sequence, const char* template_filename, const char* name, const char* blast_filename);
	int blastn_file(const char* query_filename, const char* template_filename, const char* name, const char* blast_filename);
	
	// Blast parameters	
	double evalue; 
    int wordsize;
	
	// Blast results
	int hits_total;
	results_data tmp;
	
	// BLAST sequence filters
	int min_match_length;
	int identity;
	int tail_loc;
	
	// Options
	bool report;
	
private:
	
};

