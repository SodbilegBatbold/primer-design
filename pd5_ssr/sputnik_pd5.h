#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "../primer.h"

#ifndef SPUTNIK_HELYGEN_H
#define SPUTNIK_HELYGEN_H

using namespace std;

class sputnik_helygen
{
 public:
  
  primer *fwdP;
  primer *revP;
  
  sputnik_helygen();
  ~sputnik_helygen();
  
  int find_primers(char*, int start, int end, int repeat_len);

 private:
  int has_repeat(char *seq);
};

#endif
