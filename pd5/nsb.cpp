#include "nsb.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>


nsb::~nsb() {
  fin.close();
};

// This constructor is private
nsb::nsb() {
};

// Public constructor
nsb::nsb(const char* filename) {
  fin.open(filename);
  if(!fin.is_open()) cout << "Cannot open file: " << filename << endl;
}
