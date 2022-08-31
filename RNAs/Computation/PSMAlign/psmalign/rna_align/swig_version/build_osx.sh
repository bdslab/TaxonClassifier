g++ -O -c RNA_BASE.cpp
g++ -O -c RNA_DP_METHOD1.cpp

swig -perl5 -c++ rna_stem_align.i 
g++ -c `perl -MConfig -e 'print join(" ", @Config{qw(ccflags optimize cccdlflags)}, "-I$Config{archlib}/CORE")'` rna_stem_align.cpp rna_stem_align_WRAP.cxx
g++ `perl -MConfig -e 'print $Config{lddlflags}'` rna_stem_align.o rna_stem_align_WRAP.o RNA_BASE.o RNA_DP_METHOD1.o -o rna_stem_align.dylib
