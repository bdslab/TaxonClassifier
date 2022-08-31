// RNA_BASE.cpp: implementation of the RNA_BASE class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"


#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
using namespace std;


//#include "RNA_DP.h"
#include "RNA_BASE.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RNA_BASE::RNA_BASE()
{

}

RNA_BASE::~RNA_BASE()
{

}

unsigned char RNA_BASE::setbase(char ch){
	switch(ch){
	case 'a':;
	case 'A': m_base=0;break;
	case 'c':;
	case 'C': m_base=1;break;
	case 'g':;
	case 'G': m_base=2;break;
	case 'u':;
	case 'U': m_base=3;break;
	default: m_base=4;break;
	}
	return m_base;
}

char RNA_BASE::getchar(){
	switch(m_base){
	case 0: return 'A';
	case 1: return 'C';
	case 2: return 'G';
	case 3: return 'U';
	default: return '?';
	}
}

//////////////////////////////////////////////////////////////////////
// RNA_SEQ Class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RNA_SEQ::RNA_SEQ()
{

}

RNA_SEQ::~RNA_SEQ()
{

}

int RNA_SEQ::newlabel()
{
	return (++m_maxlabel);
}

int RNA_SEQ::label(int startpos,int endpos,int startlabel)
{
	int k=startpos;
	RNA_BASE *p=&(*this)[startpos];
	if(p->is_left()){
		int kl=k;
		int kr=kl+(*this)[kl].getpairoff();
		(*this)[kl].setlabel(startlabel);
		(*this)[kr].setlabel(startlabel);
		if(kr-kl>1)	label(kl+1,kr-1,newlabel());
		if(endpos>kr) label(kr+1,endpos,newlabel());
		return 1;
	}
	else for(k=startpos;k<=endpos;k++)
	{
		RNA_BASE *p=&(*this)[k];
		if(p->is_single()){
			p->setlabel(startlabel);
		}
		else if(p->is_left()){
			int kl=k;
			int kr=kl+(*this)[kl].getpairoff();
			label(kl,kr,newlabel());
			k=kr;
		}
	}
	return 1;
}

int RNA_SEQ::label()
{
	m_maxlabel=1;
	label(0,size()-1,m_maxlabel);
	return m_maxlabel;
}


int RNA_SEQ::create(ifstream &infile)
{
	//1. read the sequence.
	string str;
	RNA_BASE base;
	while(infile)
	{
		getline(infile,str);
		if(str[0]=='>') break;
		int i;
		for(i=0;i<str.size();i++)
		{
			base.setbase(str[i]);
			if(base.getbase()>=4) continue; //not a proper character.
			push_back(base);
		}
	}

	//2. read the pairs.
	int i;
	for(i=0;i<size();i++)
	{
		(*this)[i].setpairoff(0);
	}
	while(infile)
	{
		getline(infile,str);
		str+='\0';
		int l,r,n,k;
		k=str.find_first_of("0123456789");
		char * tmpp;
		l=strtol(str.data()+k,&tmpp,10); //this is the useless number.

		k=str.find_first_of("0123456789",tmpp-str.data());
		l=strtol(str.data()+k,&tmpp,10);
		k=str.find_first_of("0123456789",tmpp-str.data());
		r=strtol(str.data()+k,&tmpp,10);
		k=str.find_first_of("0123456789",tmpp-str.data());
		n=strtol(str.data()+k,&tmpp,10);

		int i;
		for(i=0;i<n;i++){
			(*this)[l-1+i].setpairoff(r-l-i-i); //suppose the input is 1 based.
			(*this)[r-1-i].setpairoff(l-r+i+i);
		}
	}

	return size();
}

int RNA_SEQ::output(ofstream &outfile)
{
	int i;
	for(i=0;i<size();i++)
	{
		outfile<<(*this)[i].getchar()<<",";
	}
	outfile<<endl;
	for(i=0;i<size();i++)
	{
		outfile<<(*this)[i].getpairoff()<<",";
	}
	outfile<<endl;
	for(i=0;i<size();i++)
	{
		outfile<<(*this)[i].getlabel()<<",";
	}
	return 1;
}

//////////////////////////////////////////////////////////////////////
// RNA_DP_PROPERTY Class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RNA_DP_PROPERTY::RNA_DP_PROPERTY()
{

}

RNA_DP_PROPERTY::~RNA_DP_PROPERTY()
{

}

bool RNA_DP_PROPERTY::recur(RNA_DP_SIZE &dpsize,RNA_DP_RECURLIST& candlist)
{
	RNA_DP_RECURINFO onecand;
	return false;
	//set candlist here.
}

/*
float RNA_DP_PROPERTY::min_cost(RNA_DP_SIZE& dpsize)
{
	//if already computed, return the storage:
	RNA_DP_CELL *thecell=&p_matrix->getcell(dpsize.l1,dpsize.l2,dpsize.r1,dpsize.r2);
	if(thecell->cost<BEYONDMAXCOST)
		return thecell->cost;

	//generate the candlist
	RNA_DP_RECURLIST candlist;
	recur(dpsize,candlist);

	//calculate the minimum cost of the dp;
	float mincost=MAXCOST;
	int min_index=0;
	int i;
	for(i=0;i<candlist.size();i++)
	{
		RNA_DP_RECURLIST::iterator p=&candlist[i];
		float cost=0;
		if(p->end_of_recur){
			cost=p->cost_of_partition;
		}
		else{
			cost=p->cost_of_partition+min_cost(p->component[0])
				+min_cost(p->component[1]);
		}
		if(cost<mincost){
			mincost=cost;
			min_index=i;
		}
	}

	//now minp is the minimum partition and mincost is the mincost
	//record it to p_matrix;
	thecell->cost=mincost;
	thecell->min_recur_index=min_index;  //use this to do backtracking.

	return mincost;
}
*/

///////////////////////////////////
//RNA_DP_MATRIX CLASS
///////////////////////////////////

RNA_DP_MATRIX::RNA_DP_MATRIX()
{
		m_dim1=m_dim2=m_dim3=m_dim4=0;
		m_buf=NULL;
};

RNA_DP_MATRIX::~RNA_DP_MATRIX()
{
	int i1,i2,i3;
	if(m_buf==NULL) return;
	for(i1=0;i1<m_dim1;i1++)
	{
		if(m_buf[i1]==NULL) continue;
		for(i2=0;i2<m_dim2;i2++)
		{
			if(m_buf[i1][i2]==NULL) continue;
			for(i3=0;i3<m_dim3;i3++)
			{
				if(m_buf[i1][i2][i3]==NULL) continue;
				delete m_buf[i1][i2][i3];
				m_buf[i1][i2][i3]=NULL;
			}
			delete m_buf[i1][i2];
			m_buf[i1][i2]=NULL;
		}
		delete m_buf[i1];
		m_buf[i1]=NULL;
	}
	delete m_buf;
	m_buf=NULL;
}


int RNA_DP_MATRIX::allocate_memory(RNA_SEQ &rna1,RNA_SEQ &rna2)
{
	m_dim1=m_dim2=rna1.size()+2;
	m_dim3=m_dim4=rna2.size()+2;
	int i1,i2,i3;
	if((m_buf=new RNA_DP_CELL***[m_dim1])==NULL) return -1;
	for(i1=0;i1<m_dim1;i1++)
	{
		if((m_buf[i1]=new RNA_DP_CELL**[m_dim2])==NULL) return -1;
		for(i2=0;i2<m_dim2;i2++)
		{
			if((m_buf[i1][i2]=new RNA_DP_CELL*[m_dim3])==NULL) return -1;
			for(i3=0;i3<m_dim3;i3++)
			{
				if((m_buf[i1][i2][i3]=new RNA_DP_CELL[m_dim4])==NULL) return -1;
			}
		}
	}
	return 1;
}


bool RNA_SEQ::is_pair(int l, int r)
{
	return ((l!=r)&&(((*this)[l].getpairoff()+l)==r));
}
