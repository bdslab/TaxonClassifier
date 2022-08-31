%module rna_stem_align
%include "std_string.i"
%include "std_vector.i"

namespace std {
  %template(vector_i) vector<int>;
  %template(vector_i2) vector< vector<int> >;
  %template(vector_s) vector<string>;
  %template(vector_s2) vector < vector<string> >;
}

%{
#include "rna_stem_align.h"
%}

%include "rna_stem_align.h";
