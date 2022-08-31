rna_align : RNA_BASE.o RNA_DP.o RNA_DP_METHOD1.o makefile
	g++ -O -o dp  RNA_BASE.o RNA_DP.o RNA_DP_METHOD1.o

RNA_BASE.o : RNA_BASE.cpp RNA_BASE.h
	g++ -O -c RNA_BASE.cpp

RNA_DP_METHOD1.o : RNA_DP_METHOD1.cpp RNA_DP_METHOD1.h
	g++ -O -c RNA_DP_METHOD1.cpp

RNA_DP.o : RNA_DP.cpp RNA_BASE.h RNA_DP.h RNA_DP_METHOD1.h
	g++ -O -c RNA_DP.cpp

