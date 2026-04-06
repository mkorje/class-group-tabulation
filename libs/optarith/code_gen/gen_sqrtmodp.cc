/**
 * Generates sqrtmodp_list.c
 */
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <stdlib.h>

#include "primes.h"

using namespace std;

int line_elem = 0;
int para_elem = 0;
int total_elem = 0;
const int elem_per_line = 10; 
unsigned int max_p = 1000;

template<typename T>
string to_string(T t) {
  stringstream ss;
  ss << t;
  return ss.str();
}

void new_para() {
  line_elem = 0;
  para_elem = 0;
}

template<class T>
void print_elem(ofstream& f, T i) {
  if (para_elem > 0) {
    f << ',';
    if (line_elem % elem_per_line == 0) {
      f << "\n  ";
    } else {
      f << ' ';
    }
  }
  f << i;
  line_elem ++;
  para_elem ++;
  total_elem ++;
}

int main(int argc, char** argv) {
  if (argc != 2 && argc != 4) {
    cout << "Usage: " << argv[0] << " max_prime [output_c output_h]" << endl;
    cout << endl;
    return 0;
  }
  max_p = atoi(argv[1]);

  string output_c = "sqrtmodp_list.c";
  string output_h = "sqrtmodp_list.h";
  if (argc == 4) {
    output_c = argv[2];
    output_h = argv[3];
  }

  cout << "Generating " << output_c << " and " << output_h << endl;
	
  ofstream fc, fh;
  fc.open(output_c.c_str());
  fh.open(output_h.c_str());
  if (!fc.is_open() || !fh.is_open()) {
    cerr << "Unable to open output files." << endl;
    return 1;
  }
  fh << "#pragma once" << endl;
  fh << "#ifndef SQRTMODP_LIST__INCLUDED" << endl;
  fh << "#define SQRTMODP_LIST__INCLUDED" << endl;
  fh << endl;
  fh << "extern const int* sqrtmodp[];" << endl;
  fh << "extern const int sqrtmodp_maxp;" << endl;
  fh << endl;
	
  for (unsigned int prime_index = 0; prime_index < prime_list_count && prime_list[prime_index] <= max_p;  prime_index ++) {
    unsigned int p = prime_list[prime_index];
    fh << "extern const int sqrtmod_" << p << "[];" << endl;
    fc << "const int sqrtmod_" << p << "[] = {" << endl;
    fc << "  ";
    new_para();

    int* sqrt_table = new int[p];
    sqrt_table[0] = 0;
    for (unsigned int i = 1; i < p; i++) {
      sqrt_table[i] = -1;
    }

    for (unsigned int r = 1; r < p; r++) {
      unsigned int s = ((unsigned long long) r * r) %p;
      if (sqrt_table[s] == -1) {
        sqrt_table[s] = r;
      }
    }

    for (unsigned int s = 0; s < p; s++) {
      print_elem(fc, sqrt_table[s]);
    }
    delete[] sqrt_table;

    fc << "\n};\n";
    fc << endl;
  }

  // reference table
  fc << "const int sqrtmodp_maxp = " << max_p << ";" << endl;
  fc << "const int* sqrtmodp[] = {\n  ";
  unsigned int i = 0;
  new_para();
  for (unsigned int prime_index = 0; prime_index < prime_list_count && prime_list[prime_index] <= max_p;  prime_index ++) {
    unsigned int p = prime_list[prime_index];
    while (i < p) {
      print_elem(fc, 0);
      i ++;
    }
    string s = "sqrtmod_";
    s += to_string(p);
    print_elem(fc, s);
    i ++;
  }
  while (i <= max_p) {
    print_elem(fc, 0);
    i ++;
  }
  fc << "\n};\n";
  fc << endl;
  fh << endl;
  fh << "#endif  // SQRTMODP_LIST__INCLUDED" << endl;
  fh << endl;
  fc.close();
  fh.close();
  return 0;
}
