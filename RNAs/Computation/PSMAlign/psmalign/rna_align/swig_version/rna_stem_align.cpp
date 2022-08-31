#include <iostream>
#include <vector>

using namespace std;

#include "RNA_DP_METHOD1.h"

vector < vector<string> > align_stem(string stem1_base_seq_str, string stem2_base_seq_str, vector < vector < int > > stem1_base_pairs, vector < vector < int > > stem2_base_pairs,
				     int stem1_upstream_end, int stem1_downstream_start, int stem2_upstream_end, int stem2_downstream_start)
{
  RNA_SEQ rna1, rna2;

  rna1.create(stem1_base_seq_str, stem1_base_pairs);
  rna2.create(stem2_base_seq_str, stem2_base_pairs);

  RNA_DP_MATRIX1 dpmatrix;
  RNA_DP_METHOD1 dpmethod;
	
  dpmethod.p_matrix=&dpmatrix;
  dpmethod.p_rna1=&rna1;
  dpmethod.p_rna2=&rna2;

  if (dpmatrix.allocate_memory(rna1, rna2) < 0)
  {
    cout << "not enough memory!" << endl;
    vector < vector<string> > empty_results(5);
    return empty_results;
  }

  int i,j,k,i1,k1;
  k = rna1.size();
  k1 = rna2.size();

  for (k =- 1; k < (int)rna1.size(); k++)
  {
      for (i = 0; i < (int)rna1.size() - k; i++)
      {
	  if ((k == -1) && (i > 0))
	    continue;
	  if (dpmatrix.legal_segment(i, i + k))
	  {
	      for (k1 = -1; k1 < (int)rna2.size(); k1++)
	      {
		  for (i1 = 0; i1 < (int)rna2.size() - k1; i1++)
		    {
		      RNA_DP_SIZE size;
		      size.is_interval1 = true;
		      size.is_interval2 = true;
		      size.l1 = i;
		      size.r1 = i+k;
		      size.l2 = i1;
		      size.r2 = i1+k1;
		      dpmethod.DP(size, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start);
		    }
		}
	  }
      }
  }

  RNA_DP_SIZE size;
  size.l1=0;
  size.r1=rna1.size()-1;
  size.l2=0;
  size.r2=rna2.size()-1;

  float alignment_cost = dpmatrix.getcell(size).cost;
  
  return dpmethod.write_align(alignment_cost, stem1_upstream_end, stem1_downstream_start, stem2_upstream_end, stem2_downstream_start);
}
