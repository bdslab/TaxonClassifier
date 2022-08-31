// RNA_BASE.h: interface for the RNA_BASE class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RNA_BASE_H__D514358A_80BE_11D4_9347_00E09806FDAE__INCLUDED_)
#define AFX_RNA_BASE_H__D514358A_80BE_11D4_9347_00E09806FDAE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define MAXCOST 1e30
#define BEYONDMAXCOST 1e31
#define UNTOUCHEDCOST 1e32

class RNA_BASE  
{
	unsigned char m_base; //0,1,2,3 w.r.t. a,c,g,u.
	int m_pairoff; //pairoff==0: this is a single base;
				//pairoff!=0: this is paired up with a base at pairoff.
	int m_label;  //use this to determine if they are under the cover of
				//a same arc.
public:
	RNA_BASE();
	virtual ~RNA_BASE();

public:  //set or get the attributes.
	unsigned char getbase(){
		return m_base;
	}	//return the base (which is the number)
	unsigned char setbase(char ch);
	void setpairoff(int offset){
		m_pairoff=offset;
	}
	void setlabel(int label){
		m_label=label;
	}
	int getlabel(){ return m_label;};
	int getpairoff(){
		return m_pairoff;
	}
public:  //service
	char getchar();  //return the actual character of the code
	bool is_left(){
		return (m_pairoff>0);
	}
	bool is_right(){
		return m_pairoff<0;
	}
	bool is_single(){
		return m_pairoff==0;
	}
};

//////////////////////////////
//  A sequence of RNA bases
//////////////////////////////
class RNA_SEQ:public vector<RNA_BASE>
{
public:
	RNA_SEQ();
	virtual ~RNA_SEQ();

private:
	int m_maxlabel;
	int newlabel();
	int label(int startpos,int endpos, int startlabel);  
	//label the sequence recursively, return the last pos that is labelled.
public:
	bool is_pair(int l,int r);
	int output(ofstream& outfile);
	//int create(ifstream &infile);
	int create(string stem_base_seq_str, vector< vector < int > > stem_base_pairs);
	int label();//label the sequence with 1,2,3,.... 
				//such that all bases that are under a same pair is labelled with a
				//unique number.
};

//////////////////////////
// the property of the Dynamic Programming.
struct RNA_DP_SIZE
{
	bool is_interval1;
	bool is_interval2;  //l2,r2 indicate an interval
	int l1,l2,r1,r2;    //l1,l2,r1,r2 = -1 if they are space.
};

struct RNA_DP_RECURINFO
{
//	bool end_of_recur; //if(end_of_recur) then component is invalid.
//	float cost_of_partition;
	RNA_DP_SIZE component[2];
};

typedef	vector<RNA_DP_RECURINFO> RNA_DP_RECURLIST;

struct RNA_DP_CELL
{
	float cost;
//	int min_recur_index;
};

class RNA_DP_MATRIX
{
public:
	int m_dim1,m_dim2,m_dim3,m_dim4;
	RNA_DP_CELL ****m_buf;
public:
	RNA_DP_MATRIX();
	~RNA_DP_MATRIX();
public:
	RNA_DP_CELL &getcell(RNA_DP_SIZE &size){
		if(size.l1==(size.r1+1))
			return m_buf[1][0][size.l2+1][size.r2+1];
		else
		return m_buf[size.l1+1][size.r1+1][size.l2+1][size.r2+1];
	}
	RNA_DP_CELL &getcell(int i1,int i2,int i3,int i4)
	{
		return m_buf[i1][i2][i3][i4];
	}
	virtual int allocate_memory(RNA_SEQ &rna1,RNA_SEQ &rna2);
};

class RNA_DP_PROPERTY  
{
public:
	RNA_SEQ* p_rna1;
	RNA_SEQ* p_rna2;
	RNA_DP_MATRIX* p_matrix;
public:
	RNA_DP_PROPERTY();
	virtual ~RNA_DP_PROPERTY();
public:
	virtual bool recur(RNA_DP_SIZE &dpsize,RNA_DP_RECURLIST& candlist);
	//set the vector candlist to contain all possible recursions.
	//if dpsize is a basic component, then return a vector of
	// length 1 with end_of_recur==true.
//	virtual float min_cost(RNA_DP_SIZE& dpsize);
};

#endif // !defined(AFX_RNA_BASE_H__D514358A_80BE_11D4_9347_00E09806FDAE__INCLUDED_)
