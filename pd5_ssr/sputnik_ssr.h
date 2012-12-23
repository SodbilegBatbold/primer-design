#include <cstring>
#include <cstdio>
#include <cstdlib>
#include "../pd5/primer.h"

#ifndef SPUTNIK_SSR_H
#define SPUTNIK_SSR_H

using namespace std;

class sputnik_ssr
{
 public:
  
  primer *fwdP;
  primer *revP;
  
  sputnik_ssr();
  ~sputnik_ssr();
  
  int find_primers(char*, int start, int end, int repeat_len);

 private:
  int has_repeat(char *seq);
};

#endif
