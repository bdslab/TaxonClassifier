// RNA_DP.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
using namespace std;

#include "RNA_DP_METHOD1.h"
//#include "RNA_DP.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// The one and only application object

/*
CWinApp theApp;

int mymain(int argc, char * argv[]);

int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRetCode = 0;

	// initialize MFC and print and error on failure
	if (!AfxWinInit(::GetModuleHandle(NULL), NULL, ::GetCommandLine(), 0))
	{
		// TODO: change error code to suit your needs
		cerr << _T("Fatal Error: MFC initialization failed") << endl;
		nRetCode = 1;
	}
	else
	{
		// TODO: code your application's behavior here.
		CString strHello;
		strHello.LoadString(IDS_HELLO);
		cout << (LPCTSTR)strHello << endl;

		mymain(argc,argv);
	}

	return nRetCode;
}
*/
//int mymain(int argc, char * argv[])
int main(int argc, char * argv[])
{
	if(argc!=4){
		cout<<"Usage: rna_align rna1 rna2 output!"<<endl;
		return 0;
	}
	char *frna1=argv[1];
	char *frna2=argv[2];
	char *outputfile=argv[3];

	RNA_SEQ rna1,rna2;

	ifstream inf(frna1);
	rna1.create(inf);
	inf.close();
	ifstream inf2(frna2);
	rna2.create(inf2);
	inf2.close();

	ofstream outf(outputfile);

//	rna1.label();
//	rna1.output(outf);
//	outf.close();

	//
	RNA_DP_MATRIX1 dpmatrix;
	RNA_DP_METHOD1 dpmethod;
	
	dpmethod.p_matrix=&dpmatrix;
	dpmethod.p_rna1=&rna1;
	dpmethod.p_rna2=&rna2;

	if(dpmatrix.allocate_memory(rna1,rna2)<0)
	{
		cout<<"not enough memory!"<<endl;
		return 0;
	}
	else{
		cout<<"memory allocated succesfully!"<<endl;
	}
	int i,j,k,i1,k1;
	k=rna1.size();
	k1=rna2.size();
	for(k=-1;k<(int)rna1.size();k++)
	{
		for(i=0;i<(int)rna1.size()-k;i++)
		{
			if((k==-1)&&(i>0)) continue;
			if(dpmatrix.legal_segment(i,i+k))
			{
				cout<<"["<<i+1<<","<<i+1+k<<"]";
				cout<<";  ";
				for(k1=-1;k1<(int)rna2.size();k1++)
				{
					for(i1=0;i1<(int)rna2.size()-k1;i1++)
					{
						RNA_DP_SIZE size;
						size.is_interval1=true;
						size.is_interval2=true;
						size.l1=i;
						size.r1=i+k;
						size.l2=i1;
						size.r2=i1+k1;
						dpmethod.DP(size);
					}
				}
			}
		}
		cout<<endl;
	}

	RNA_DP_SIZE size;
	size.l1=0;
	size.r1=rna1.size()-1;
	size.l2=0;
	size.r2=rna2.size()-1;

	cout<<dpmatrix.getcell(size).cost;

	dpmethod.write_align(outf,100);
	outf.close();

	return 0;
}


