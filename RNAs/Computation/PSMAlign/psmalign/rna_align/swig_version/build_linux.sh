g++ -fPIC -O -c RNA_BASE.cpp
g++ -fPIC -O -c RNA_DP_METHOD1.cpp

swig -perl5 -c++ rna_stem_align.i
g++ -fPIC -c rna_stem_align.cpp
g++ -fPIC -c rna_stem_align_wrap.cxx `perl -MConfig -e 'print "-I$Config{archlib}/CORE"'` -Dbool=char
g++ -shared rna_stem_align.o rna_stem_align_wrap.o RNA_BASE.o RNA_DP_METHOD1.o -o rna_stem_align.so
