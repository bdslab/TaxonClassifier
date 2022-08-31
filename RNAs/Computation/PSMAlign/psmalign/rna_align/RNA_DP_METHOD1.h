// RNA_DP_METHOD1.h: interface for the RNA_DP_METHOD1 class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RNA_DP_METHOD1_H__D514358D_80BE_11D4_9347_00E09806FDAE__INCLUDED_)
#define AFX_RNA_DP_METHOD1_H__D514358D_80BE_11D4_9347_00E09806FDAE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <vector>
#include "RNA_BASE.h"

using namespace std;


class RNA_DP_MATRIX1 : public RNA_DP_MATRIX
{
public:
	RNA_DP_MATRIX1();
	virtual ~RNA_DP_MATRIX1();
public:
	bool legal_segment(int i1,int i2,RNA_SEQ &rna1);
	bool legal_segment(int i1,int i2);
	int allocate_memory(RNA_SEQ &rna1,RNA_SEQ &rna2);
};

typedef vector<int> DEL_QUE;

class RNA_DP_METHOD1 : public RNA_DP_PROPERTY  
{
public:
	RNA_DP_MATRIX1 *p_matrix;
public:
	RNA_DP_METHOD1();
	virtual ~RNA_DP_METHOD1();
	bool recur(RNA_DP_SIZE &dpsize,RNA_DP_RECURLIST& candlist);
	float DP(RNA_DP_SIZE& dpsize);
	float basic_cost(int l1,int r1,int l2,int r2);
	bool basic_cost(RNA_DP_SIZE &dpsize, float *cost, DEL_QUE *delque1=NULL,DEL_QUE *delque2=NULL);
public:
	void write_align(ostream &out,int len_per_line);
	int back_track(RNA_DP_SIZE dpsize,DEL_QUE *delque1, DEL_QUE *delque2);
//	float min_cost(RNA_DP_SIZE& dpsize);//useless
};

#endif // !defined(AFX_RNA_DP_METHOD1_H__D514358D_80BE_11D4_9347_00E09806FDAE__INCLUDED_)
