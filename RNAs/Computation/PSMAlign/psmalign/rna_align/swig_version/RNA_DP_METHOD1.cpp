// RNA_DP_METHOD1.cpp: implementation of the RNA_DP_METHOD1 class.
//
//////////////////////////////////////////////////////////////////////

//#include "stdafx.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
using namespace std;

#include "RNA_DP.h"
#include "RNA_DP_METHOD1.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RNA_DP_METHOD1::RNA_DP_METHOD1()
{

}

RNA_DP_METHOD1::~RNA_DP_METHOD1()
{

}

#define c_a0 1.5	//(a,u) v.s. [a,u]  break an arc
#define c_a1 1.75	//(a,u) v.s. [a,_]  break an arc and remove a base
#define c_a2 2	//(a,u) v.s. [_,_]  remove a pair
//#define c_a3 2  //(a,u) v.s. [g,_]  break an arc, remove and replace a base.
#define c_rep2 1   //(a,u) v.s. (g,c)  replace a pair
#define c_rep1 1   //a v.s. g     replace a base
#define c_del  1   //a v.s. _     delete/insert a base
#define x_del  ((c_a2-c_a0)/2.0)   //delete a arc related base.
#define c_break_penalty (c_a0/2.0)

#define c_a3 (c_del + c_rep1 + c_break_penalty)  //(a,u) v.s. [g,_]  break an arc, remove and replace a base.
#define c_a4 (c_a0+c_rep1+c_rep1)  //(a,u) v.s. [g,g]  break an arc, and replace two bases.
#define c_a5 (c_a0+c_rep1)  //(a,u) v.s. [a,g]

//suppose no pair is replaced at only one of the two bases.


float RNA_DP_METHOD1::basic_cost(int l1,int r1,int l2,int r2)
{
	return 1;
}

bool RNA_DP_METHOD1::basic_cost(RNA_DP_SIZE &dpsize, float *cost, int stem1_upstream_end, int stem1_downstream_start, int stem2_upstream_end,
				int stem2_downstream_start, DEL_QUE *delque1, DEL_QUE *delque2)
{
	int len1,len2;
//	bool adjacent_pair1,adjacent_pair2;
//	adjacent_pair1=( (dpsize.r1-dpsize.l1==1)&&p_rna1->is_pair(dpsize.l1,dpsize.r1));
//	adjacent_pair2=( (dpsize.r2-dpsize.l2==1)&&p_rna2->is_pair(dpsize.l2,dpsize.r2));
	if( (!dpsize.is_interval1)&&(!dpsize.is_interval2) )
	{//there are four possibilities:
		//1. pair vs (a,u), pair vs (a,_), pair vs (_,a), pair vs (_,_);
		if((dpsize.l2==-1)&&(dpsize.r2==-1)){
			//pair vs (_,_);
			*cost=c_a2;
			//record the deletion
			if(delque1!=NULL){
				delque1->push_back(dpsize.l1);
				delque1->push_back(dpsize.r1);
			}
		}
		else if((dpsize.l2==-1)&&(dpsize.r2!=-1)){
			//pair vs (_,a)
		  if (dpsize.r2 < stem2_downstream_start) {
		    *cost += MAXCOST;
		  }
		  else {
			if((*p_rna1)[dpsize.r1].getbase()==(*p_rna2)[dpsize.r2].getbase())
				*cost=c_a1;
			else
				*cost=c_a3;

			//if(!(*p_rna2)[dpsize.r2].is_single()){
			if((*p_rna2)[dpsize.r2].is_right())
				*cost+=c_break_penalty;
			else if((*p_rna2)[dpsize.r2].is_left())
				*cost+=MAXCOST;

			//record the deletion
			if(delque1!=NULL){
				delque1->push_back(dpsize.l1);
			}
		  }
		}
		else if(dpsize.r2==-1){
			//pair vs (a,_)
		  if (dpsize.l2 > stem2_upstream_end) {
		    *cost += MAXCOST;
		  }
		  else {
			if((*p_rna1)[dpsize.l1].getbase()==(*p_rna2)[dpsize.l2].getbase())
				*cost=c_a1;
			else
				*cost=c_a3;

                        //if(!(*p_rna2)[dpsize.l2].is_single()){
                        if((*p_rna2)[dpsize.l2].is_left())
                                *cost+=c_break_penalty;
                        else if((*p_rna2)[dpsize.l2].is_right())
                                *cost+=MAXCOST;

			//record the deletion
			if(delque1!=NULL){
				delque1->push_back(dpsize.r1);
			}
		  }
		}
		else{//pair vs (a,u)
			//if(l2,r2) is a pair
			if((*p_rna2)[dpsize.l2].getpairoff()+dpsize.l2==dpsize.r2)
			{
			  /*
				if((*p_rna2)[dpsize.l2].getbase()==(*p_rna1)[dpsize.l1].getbase())
					*cost=0;
				else
					*cost=c_rep2;
			  */
			  *cost = 0;
			  if((*p_rna2)[dpsize.l2].getbase()!=(*p_rna1)[dpsize.l1].getbase())
			    *cost += c_rep2 / 2.0;

			  if((*p_rna2)[dpsize.r2].getbase()!=(*p_rna1)[dpsize.r1].getbase())
			    *cost += c_rep2 / 2.0;
			}
			//if (l2,r2) are not pair
			else{
			  if (dpsize.l2 > stem2_upstream_end || dpsize.r2 < stem2_downstream_start) {
			    *cost += MAXCOST;
			  }
			  else {
				if((*p_rna2)[dpsize.l2].getbase()!=(*p_rna1)[dpsize.l1].getbase())
				{
					if((*p_rna2)[dpsize.r2].getbase()!=(*p_rna1)[dpsize.r1].getbase())
						*cost=c_a4;
					else
						*cost=c_a5;
				}
				else{
//					int a=(*p_rna2)[dpsize.r2].getbase();
//					int b=(*p_rna1)[dpsize.r1].getbase();
					if((*p_rna2)[dpsize.r2].getbase()!=(*p_rna1)[dpsize.r1].getbase())
						*cost=c_a5;
					else
						*cost=c_a0;
				}

	                        //if(!(*p_rna2)[dpsize.r2].is_single()){
        	                if((*p_rna2)[dpsize.r2].is_right())
                	                *cost+=c_break_penalty;
                        	else if((*p_rna2)[dpsize.r2].is_left())
                                	*cost+=MAXCOST;

	                        //if(!(*p_rna2)[dpsize.l2].is_single()){
        	                if((*p_rna2)[dpsize.l2].is_left())
                	                *cost+=c_break_penalty;
                        	else if((*p_rna2)[dpsize.l2].is_right())
                                	*cost+=MAXCOST;
			  }
			}
		}
		return true;
	}
	else if(((len1=dpsize.r1-dpsize.l1)<=0)&&((len2=dpsize.r2-dpsize.l2)<=0))
	//2. seg1 has at most 1 base and seg2 has at most 1 base.
	{
		if(len1==0){
			if(len2<0) 
			{ // base v.s. _
				*cost=c_del;
				//record the deletion
				if(delque1!=NULL){
					delque1->push_back(dpsize.l1);
				}
			}
			else{ //base vs. base
			  if ((dpsize.l1 > stem1_upstream_end && dpsize.l2 < stem2_downstream_start) || (dpsize.l1 < stem1_downstream_start && dpsize.l2 > stem2_upstream_end)) {
			    *cost += MAXCOST;
			  }
			  else {
				if((*p_rna1)[dpsize.l1].getbase()==(*p_rna2)[dpsize.l2].getbase())
					*cost=0;
				else
					*cost=c_rep1;

				if(!(*p_rna2)[dpsize.l2].is_single()){
					*cost+=c_break_penalty;
					if(!(*p_rna1)[dpsize.l1].is_single())
					{
						bool left1=(*p_rna1)[dpsize.l1].is_left();
						bool left2=(*p_rna1)[dpsize.l2].is_left();
						if( (left1&&(!left2)) || (left2 &&(!left1)) )
							*cost+=MAXCOST;
					}
				}
			  }
			}
		}
		else{//len1<0  
			if(len2<0) *cost=0;  //- v.s. -
			else{   // - v.s. base.
				if((*p_rna2)[dpsize.l2].is_single())
					*cost=c_del;
				else
					*cost=x_del+c_break_penalty;

				//record the deletion
				if(delque2!=NULL){
					delque2->push_back(dpsize.l2);
				}
				
			}
		}

		return true;
	}
//	else if(adjacent_pair1&&(dpsize.r2-dpsize.l2<0))
//	{
//		*cost=c_a2;
//		return true;
//	}
	else
		return false;
}

bool RNA_DP_METHOD1::recur(RNA_DP_SIZE &dpsize,RNA_DP_RECURLIST& candlist)
{
	RNA_DP_RECURINFO onecand;
	RNA_DP_SIZE size1,size2;
	//set candlist here.
	//////////////////////
	// if dpsize.l1
	////////////////////////////////////
	//1. decrease the size of the second rna by 1
	if(dpsize.r2>=dpsize.l2){
		//reduce at the right side.
		size1=dpsize;
		size1.r2--;

		size2.is_interval1=true;
		size2.is_interval2=true;
		size2.l1=dpsize.r1+1;
		size2.r1=dpsize.r1;
		size2.l2=dpsize.r2;
		size2.r2=dpsize.r2;

//		onecand.end_of_recur=false;

		onecand.component[0]=size1;
		onecand.component[1]=size2;

		candlist.push_back(onecand);
		//reduce at the left side.
		size1=dpsize;
		size1.l2++;

		size2.is_interval1=true;
		size2.is_interval2=true;
		size2.l1=dpsize.l1;
		size2.r1=dpsize.l1-1;
		size2.l2=dpsize.l2;
		size2.r2=dpsize.l2;

//		onecand.end_of_recur=false;

		onecand.component[0]=size1;
		onecand.component[1]=size2;

		candlist.push_back(onecand);

	};
	//2. decrease the size of the first rna.
	//2.0 if dpsize.r1 is single
	if(dpsize.l1<dpsize.r1)
	{
	if((*p_rna1)[dpsize.r1].is_single())
	{
		//base match _
		size1.is_interval1=true;
		size1.is_interval2=true;
		size2.is_interval1=true;
		size2.is_interval2=true;

		size1.l1=dpsize.l1;
		size1.r1=dpsize.r1-1;
		size1.l2=dpsize.l2;
		size1.r2=dpsize.r2;
		size2.l1=dpsize.r1;
		size2.r1=dpsize.r1;
		size2.l2=dpsize.r2+1;
		size2.r2=dpsize.r2;

//		onecand.end_of_recur=false;
//		onecand.cost_of_partition=0;
		onecand.component[0]=size1;
		onecand.component[1]=size2;

		candlist.push_back(onecand);

		//base match base
		if(dpsize.r2>=dpsize.l2)
		{
			size1.r2--;
			size2.l2--;
			onecand.component[0]=size1;
			onecand.component[1]=size2;

			candlist.push_back(onecand);
		}
	}
	//2.1 if dpsize.l1 and dpsize.r1 are a pair
	else if(((*p_rna1)[dpsize.l1].getpairoff()+dpsize.l1)==dpsize.r1)
	{
		size1.is_interval1=true;
		size1.is_interval2=true;
		size1.l1=dpsize.l1+1;
		size1.r1=dpsize.r1-1;

		size2.is_interval1=false;
		size2.is_interval2=false;
		size2.l1=dpsize.l1;
		size2.r1=dpsize.r1;

		//2.1.1 pair match pair.
		if(dpsize.r2-dpsize.l2>=1){
			size1.l2=dpsize.l2+1;
			size1.r2=dpsize.r2-1;
			size2.l2=dpsize.l2;
			size2.r2=dpsize.r2;

//			onecand.end_of_recur=false;
			onecand.component[0]=size1;
			onecand.component[1]=size2;

			candlist.push_back(onecand);
		}

		//2.1.2 pair match empty.
		if(true){
			size1.l2=dpsize.l2;
			size1.r2=dpsize.r2;
			size2.l2=-1;
			size2.r2=-1;

//			onecand.end_of_recur=false;
			onecand.component[0]=size1;
			onecand.component[1]=size2;
	
			candlist.push_back(onecand);
		}
		//2.1.3 pair match (a,_) and (_,a).
		if(dpsize.r2-dpsize.l2>=0){
			//2.1.3.1 pair match (a,_).
			size1.l2=dpsize.l2+1;
			size1.r2=dpsize.r2;
			size2.l2=dpsize.l2;
			size2.r2=-1;

//			onecand.end_of_recur=false;
			onecand.component[0]=size1;
			onecand.component[1]=size2;

			candlist.push_back(onecand);

			//2.1.3.2 pair match (_,a).
			size1.l2=dpsize.l2;
			size1.r2=dpsize.r2-1;
			size2.l2=-1;
			size2.r2=dpsize.r2;

//			onecand.end_of_recur=false;
			onecand.component[0]=size1;
			onecand.component[1]=size2;

			candlist.push_back(onecand);
		}
	}
	//2.2 if dpsize.l1 and dpsize.r1 are not pair.
	//    and dpsize.r1 is not single.
	else
	{
		//find a partition of rna1
		int i=(*p_rna1)[dpsize.r1].getpairoff()+dpsize.r1-1;
		if(i>=dpsize.l1){
			size1.is_interval1=true;
			size1.is_interval2=true;
			size1.l1=dpsize.l1;
			size1.r1=i;
			size1.l2=dpsize.l2;
			//size1.r2 to be defined later
			size2.is_interval1=true;
			size2.is_interval2=true;
			size2.l1=i+1;
			size2.r1=dpsize.r1;
			//size2.l2 to be defined later
			size2.r2=dpsize.r2;

			int k;
			for(k=dpsize.l2-1;k<=dpsize.r2;k++)
			{
				size1.r2=k;
				size2.l2=k+1;
				onecand.component[0]=size1;
				onecand.component[1]=size2;
				candlist.push_back(onecand);
			}
		}
	}
	}

	if(candlist.empty()) return false;
	else return true;
}



float RNA_DP_METHOD1::DP(RNA_DP_SIZE& dpsize, int stem1_upstream_end, int stem1_downstream_start, int stem2_upstream_end, int stem2_downstream_start)
{
	//if already computed, return the storage:
	RNA_DP_CELL *thecell=&p_matrix->getcell(dpsize);
//	if(thecell->cost<BEYONDMAXCOST)
//		return thecell->cost;

	float cost;
	if(basic_cost(dpsize,&cost, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start)){
		thecell->cost=cost;
//		thecell->min_recur_index=-1;
		return cost;
	}

	//generate the candlist
	RNA_DP_RECURLIST candlist;
	recur(dpsize,candlist);

	//calculate the minimum cost of the dp;
	float mincost=MAXCOST;
	int min_index=0;
	int i;
	for(i=0;i<candlist.size();i++)
	{
	  //	  	RNA_DP_RECURLIST::iterator p=&candlist[i];
	  RNA_DP_RECURINFO* p = &candlist[i];
		float cost,cost0,cost1;

		if(!basic_cost(p->component[0],&cost0, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start))
			cost0=p_matrix->getcell(p->component[0]).cost;

		if(!basic_cost(p->component[1],&cost1, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start))
			cost1=p_matrix->getcell(p->component[1]).cost;

		cost=cost0+cost1;

		if(cost<mincost){
			mincost=cost;
			min_index=i;
		}
	}

	//now minp is the minimum partition and mincost is the mincost
	//record it to p_matrix;
	thecell->cost=mincost;
//	thecell->min_recur_index=min_index;  //use this to do backtracking.

	return mincost;
}

//////////////////////////////////////////////////////////////////////
// RNA_DP_MATRIX1 Class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

RNA_DP_MATRIX1::RNA_DP_MATRIX1()
{

}

RNA_DP_MATRIX1::~RNA_DP_MATRIX1()
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

bool RNA_DP_MATRIX1::legal_segment(int i1,int i2,RNA_SEQ &rna1)
{
	if(i1==i2+1) return (i1==0);
	else if(i1>(i2+1)) return false;
	else if((i1 >= 0) && (i2 >= 0) && rna1.is_pair(i1,i2)) return true;
	
	bool legal_i1=false;
	if(i1==0) legal_i1=true;
	else if((i1 > 0) &&rna1[i1-1].is_left()) legal_i1=true;

	if(!legal_i1) return false;
	int j;
	//if there is an arc start between i1 and i2 but end
	//beyond i1 and i2, then this is illegal.
	for(j=i1;j<=i2;j++)
	{
		int k=j+rna1[j].getpairoff();
		if((k>i2)||(k<i1)) return false;
	}
	return true;
}

bool RNA_DP_MATRIX1::legal_segment(int i1,int i2)
{
	if((i1<-1)||(i2<-1)||(i2>=m_dim2)||(i1>=m_dim2))
		return false;
	return (m_buf[i1+1][i2+1]!=NULL);
}

int RNA_DP_MATRIX1::allocate_memory(RNA_SEQ &rna1,RNA_SEQ &rna2)
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
			if(!legal_segment(i1-1,i2-1,rna1))
			{
				m_buf[i1][i2]=NULL;
				continue;
			}
			//else
			if((m_buf[i1][i2]=new RNA_DP_CELL*[m_dim3])==NULL) return -1;
			else for(i3=0;i3<m_dim3;i3++)
			{
				if((m_buf[i1][i2][i3]=new RNA_DP_CELL[m_dim4])==NULL) return -1;
			}
		}
		//cout<<"allocated memory "<<i1<<"of "<<m_dim1<<endl;
	}
	return 1;
}



///////////////////
/*
float RNA_DP_METHOD1::min_cost(RNA_DP_SIZE& dpsize)
{
	//if already computed, return the storage:
	RNA_DP_CELL *thecell=&p_matrix->getcell(dpsize.l1,dpsize.l2,dpsize.r1,dpsize.r2);
//	if(thecell->cost<BEYONDMAXCOST)
//		return thecell->cost;

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

int RNA_DP_METHOD1::back_track(RNA_DP_SIZE dpsize, int stem1_upstream_end, int stem1_downstream_start, int stem2_upstream_end,
			       int stem2_downstream_start, DEL_QUE *delque1, DEL_QUE *delque2)
{
	RNA_DP_CELL *thecell=&p_matrix->getcell(dpsize);

	float cost;
	if(basic_cost(dpsize,&cost, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start, delque1, delque2)){
		return 1;
	}

	//generate the candlist
	RNA_DP_RECURLIST candlist;
	recur(dpsize,candlist);

	//calculate the minimum cost of the dp;
///////////
//	int min_index=thecell->min_recur_index;
        float mincost=MAXCOST;
        int min_index=0;
        int i;
	int len1_1, len1_2, len2_1, len2_2;
        for(i=0;i<candlist.size();i++)
        {
	  //                RNA_DP_RECURLIST::iterator p=&candlist[i];
	  RNA_DP_RECURINFO* p = &candlist[i];
                float cost,cost0,cost1;
                
                if(!basic_cost(p->component[0],&cost0, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start))
                        cost0=p_matrix->getcell(p->component[0]).cost;
        
                if(!basic_cost(p->component[1],&cost1, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start))
                        cost1=p_matrix->getcell(p->component[1]).cost;

                cost=cost0+cost1;
        
                if(cost<mincost){
                        mincost=cost;
                        min_index=i;
                }
		else if (cost == mincost) {
		  RNA_DP_SIZE size1 = p->component[0];
		  RNA_DP_SIZE size2 = p->component[1];
		  if (size2.l2 == -1 || size2.r2 == -1)
		    min_index = i;
		  else if (size1.l2 == -1 || size1.r2 == -1)
		    min_index = i;
		  /*
		  else if (((len1_1 = size1.r1 - size1.l1) < 0) || ((len1_2 = size1.r2 - size1.l2) < 0) ||
			   ((len2_1 = size2.r1 - size2.l1) < 0) || ((len2_2 = size2.r2 - size2.l2) < 0))
		    min_index = i;
		  */
		}
        }
//////////////
//	RNA_DP_RECURLIST::iterator p=&candlist[min_index];
	  RNA_DP_RECURINFO* p = &candlist[min_index];

	float cost0;
	if(!basic_cost(p->component[0],&cost0, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start, delque1, delque2))
	  back_track(p->component[0], stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start, delque1, delque2);

	if(!basic_cost(p->component[1],&cost0, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start, delque1, delque2))
	  back_track(p->component[1], stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start, delque1, delque2);

	return 1;
}

/*void RNA_DP_METHOD1::write_align(ostream &out,int len_per_line)
{
	DEL_QUE delque1,delque2;
	delque1.clear();
	delque2.clear();

	RNA_DP_SIZE size;
	size.is_interval1=true;
	size.is_interval2=true;
	size.l1=0;
	size.r1=p_rna1->size()-1;
	size.l2=0;
	size.r2=p_rna2->size()-1;
	back_track(size,&delque1,&delque2);
	/////////////////////
	int length=delque1.size()+delque2.size()+p_rna1->size()+p_rna2->size();
	cout<<endl<<"length="<<length;
	length/=2;
	int i,i1,i2;
	i1=i2=0;
	string out1;
	string out2;
	string outb1;
	string outb2;

	int j1,j2;
	j1=j2=0;

	for(i=0;i<length;i++)
	{
		if(find(delque1.begin(),delque1.end(),i1)!=delque1.end()){
			out1+=(*p_rna1)[i1].getchar();
			out2+='-';
			i1++;
		}
		else if(find(delque2.begin(),delque2.end(),i2)!=delque2.end()){
			out2+=(*p_rna2)[i2].getchar();
			out1+='-';
			i2++;
		}
		else{
			out1+=(*p_rna1)[i1].getchar();
			out2+=(*p_rna2)[i2].getchar();
			i1++;
			i2++;
		}

		//write bracks into outb1,outb2:
		if(out1[out1.length()-1]!='-'){
			if((*p_rna1)[i1-1].is_left()) 
				outb1+='(';
			else if((*p_rna1)[i1-1].is_right())
				outb1+=')';
			else outb1+=' ';
		}
		else outb1+=' ';

		if(out2[out2.length()-1]!='-'){
			if((*p_rna2)[i2-1].is_left()) 
				outb2+='(';
			else if((*p_rna2)[i2-1].is_right())
				outb2+=')';
			else outb2+=' ';
		}
		else outb2+=' ';
	}

//	cout<<out1.length()<<","<<out2.length()<<","<<outb1.length()<<endl;
	for(i=0;i<length;i+=len_per_line)
	{
		out<<outb1.substr(i,len_per_line)<<endl;
		out<<outb2.substr(i,len_per_line)<<endl;
		out<<out1.substr(i,len_per_line)<<endl;
		out<<out2.substr(i,len_per_line)<<endl;
		out<<endl<<endl;
	}
}
*/

void mark_occur_pos(int &first_occur_pos, int &last_occur_pos, int occur_pos) {
  if (first_occur_pos == -1)
    first_occur_pos = occur_pos;

  last_occur_pos = occur_pos;
}

string convert_int_to_str(int value){
  stringstream ss;
  ss << value;

  return ss.str();
}

string convert_float_to_str(float value) {
  stringstream ss;
  ss << value;

  return ss.str();
}

vector < vector<string> > RNA_DP_METHOD1::write_align(float alignment_cost, int stem1_upstream_end, int stem1_downstream_start,
						      int stem2_upstream_end, int stem2_downstream_start)
{
	DEL_QUE delque1,delque2;
	delque1.clear();
	delque2.clear();

	RNA_DP_SIZE size;
	size.is_interval1=true;
	size.is_interval2=true;
	size.l1=0;
	size.r1=p_rna1->size()-1;
	size.l2=0;
	size.r2=p_rna2->size()-1;
	back_track(size, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start, &delque1,&delque2);
	/////////////////////
	int length=delque1.size()+delque2.size()+p_rna1->size()+p_rna2->size();
	length/=2;
	int i,i1,i2;
	i1=i2=0;

	string out1u, out1d, out2u, out2d; 
	string outb1u, outb1d, outb2u, outb2d;

	int first_open_brac_pos1, last_open_brac_pos1, first_close_brac_pos1, last_close_brac_pos1;
	int first_open_brac_pos2, last_open_brac_pos2, first_close_brac_pos2, last_close_brac_pos2;
	first_open_brac_pos1 = -1;
	first_close_brac_pos1 = -1;
	first_open_brac_pos2 = -1;
	first_close_brac_pos2 = -1;

	int first_base_pos1, last_base_pos1;
	int first_base_pos2, last_base_pos2;
	first_base_pos1 = -1;
	first_base_pos2 = -1;

	int written_rna1_upstream_len = 0;
	int written_rna2_upstream_len = 0;

	for(i=0;i<length;i++)
	{
		if(find(delque1.begin(),delque1.end(),i1)!=delque1.end()){		      
		        if ((*p_rna1)[i1].is_left()) {
			  outb1u += '(';
			  mark_occur_pos(first_open_brac_pos1, last_open_brac_pos1, i);
			}
			else if ((*p_rna1)[i1].is_right()) {
			  outb1d += ')';
			  mark_occur_pos(first_close_brac_pos1, last_close_brac_pos1, i);
			}
			//else if (outb1d != "" || written_rna1_upstream_len > stem1_upstream_end)
			else if (written_rna1_upstream_len > stem1_upstream_end)
			  outb1d += '.';
			else
			  outb1u += '.';

			if (written_rna1_upstream_len > stem1_upstream_end) {
			  out1d += (*p_rna1)[i1].getchar();
			  out2d += '-';
			  outb2d += '.';
			}
			else {
			  out1u += (*p_rna1)[i1].getchar();
			  out2u += '-';
			  outb2u += '.';
			  written_rna1_upstream_len++;
			}

			mark_occur_pos(first_base_pos1, last_base_pos1, i);
			i1++;
		}
		else if(find(delque2.begin(),delque2.end(),i2)!=delque2.end()){
		        if ((*p_rna2)[i2].is_left()) {
			  outb2u += '(';
			  mark_occur_pos(first_open_brac_pos2, last_open_brac_pos2, i);
			}
			else if ((*p_rna2)[i2].is_right()) {
			  outb2d += ')';
			  mark_occur_pos(first_close_brac_pos2, last_close_brac_pos2, i);
			}
			//else if (outb2d != "" || written_rna2_upstream_len > stem2_upstream_end)
			else if (written_rna2_upstream_len > stem2_upstream_end)
			  outb2d += '.';
			else
			  outb2u += '.';

			if (written_rna2_upstream_len > stem2_upstream_end) {
			  out2d += (*p_rna2)[i2].getchar();
			  out1d += '-';
			  outb1d += '.';
			}
			else {
			  out2u += (*p_rna2)[i2].getchar();
			  out1u += '-';
			  outb1u += '.';
			  written_rna2_upstream_len++;
			}

			mark_occur_pos(first_base_pos2, last_base_pos2, i);
			i2++;
		}
		else{
			if ((*p_rna1)[i1].is_left()) {
			  outb1u += '(';
			  mark_occur_pos(first_open_brac_pos1, last_open_brac_pos1, i);
			}
			else if ((*p_rna1)[i1].is_right()) {
			  outb1d += ')';
			  mark_occur_pos(first_close_brac_pos1, last_close_brac_pos1, i);
			}
			//else if (outb1d != "" || written_rna1_upstream_len > stem1_upstream_end || written_rna2_upstream_len > stem2_upstream_end)
			else if (written_rna1_upstream_len > stem1_upstream_end || written_rna2_upstream_len > stem2_upstream_end)
			  outb1d += '.';
			else
			  outb1u += '.';

			if ((*p_rna2)[i2].is_left()) {
			  outb2u += '(';
			  mark_occur_pos(first_open_brac_pos2, last_open_brac_pos2, i);
			}
			else if ((*p_rna2)[i2].is_right()) {
			  outb2d += ')';
			  mark_occur_pos(first_close_brac_pos2, last_close_brac_pos2, i);
			}
			//else if (outb2d != "" || written_rna1_upstream_len > stem1_upstream_end || written_rna2_upstream_len > stem2_upstream_end)
			else if (written_rna1_upstream_len > stem1_upstream_end || written_rna2_upstream_len > stem2_upstream_end)
			  outb2d += '.';
			else
			  outb2u += '.';

			if (written_rna1_upstream_len > stem1_upstream_end || written_rna2_upstream_len > stem2_upstream_end) {
			  out1d += (*p_rna1)[i1].getchar();
			  out2d += (*p_rna2)[i2].getchar();
			}
			else {
			  out1u += (*p_rna1)[i1].getchar();
			  out2u += (*p_rna2)[i2].getchar();
			  written_rna1_upstream_len++;
			  written_rna2_upstream_len++;
			}

			mark_occur_pos(first_base_pos1, last_base_pos1, i);
			mark_occur_pos(first_base_pos2, last_base_pos2, i);
			i1++;
			i2++;
		}
	}

	vector < vector<string> > align_results(5);

	vector<string> align_cost_result(1);
	align_cost_result[0] = convert_float_to_str(alignment_cost);
	align_results[0] = align_cost_result;

	int outb1u_len = outb1u.length();
	int outb2u_len = outb2u.length();
	vector<string> align_upstream1(4);
	align_upstream1[0] = outb1u;
	align_upstream1[1] = out1u;
	if (first_open_brac_pos2 > -1) {
	  align_upstream1[2] = convert_int_to_str(first_open_brac_pos2);
	  align_upstream1[3] = convert_int_to_str(outb2u_len - last_open_brac_pos2 - 1);
	}
	else {
	  align_upstream1[2] = convert_int_to_str(first_base_pos2);
	  align_upstream1[3] = convert_int_to_str(outb2u_len - last_base_pos2 - 1);
	}

	align_results[1] = align_upstream1;

	vector<string> align_downstream1(4);
	align_downstream1[0] = outb1d;
	align_downstream1[1] = out1d;
	if (first_close_brac_pos2 > -1) {
	  align_downstream1[2] = convert_int_to_str(first_close_brac_pos2 - outb2u_len);
	  align_downstream1[3] = convert_int_to_str(outb2d.length() + outb2u_len - last_close_brac_pos2 - 1);
	}

	align_results[2] = align_downstream1;

	vector<string> align_upstream2(4);
	align_upstream2[0] = outb2u;
	align_upstream2[1] = out2u;
	if (first_open_brac_pos1 > -1) {
	  align_upstream2[2] = convert_int_to_str(first_open_brac_pos1);
	  align_upstream2[3] = convert_int_to_str(outb1u_len - last_open_brac_pos1 - 1);
	}
	else {
	  align_upstream2[2] = convert_int_to_str(first_base_pos1);
	  align_upstream2[3] = convert_int_to_str(outb1u_len - last_base_pos1 - 1);
	}

	align_results[3] = align_upstream2;

	vector<string> align_downstream2(4);
	align_downstream2[0] = outb2d;
	align_downstream2[1] = out2d;
	if (first_close_brac_pos1 > -1) {
	  align_downstream2[2] = convert_int_to_str(first_close_brac_pos1 - outb1u_len);
	  align_downstream2[3] = convert_int_to_str(outb1d.length() + outb1u_len - last_close_brac_pos1 - 1);
	}

	align_results[4] = align_downstream2;

	return align_results;
}
