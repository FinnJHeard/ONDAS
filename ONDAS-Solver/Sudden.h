// Sudden.h: interface for the CSudden class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SUDDEN_H__FBBF483A_FE81_404E_87D9_C87F7E2D43B7__INCLUDED_)
#define AFX_SUDDEN_H__FBBF483A_FE81_404E_87D9_C87F7E2D43B7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Boundary.h"
//#include "Node.h"
#include "PathLine.h"
#include "Pipe.h"
#include "Properties.h"

class CSudden : public CBoundary  
{
public:
	CSudden();
	virtual ~CSudden();
	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rSPIPES, int** &rSPIPES_ENDS, double* &rENDCORR, 
					int id, bool EX, int npipes, int assyid, string parent_assy_res_dir, string calling_object_str);
	void RunBoundary(CProperties* pPpt, double time, int timestep);
	void General(CProperties* pPpt, double time, int timestep);
	void Enlargement(CProperties* pPpt, double &rlambda_in_c2, double &rAA_c2, double &rlambda_out1, double &rlambda_out2, int* &rpipe_flow, double time, int timestep);
	void Contraction(CProperties* pPpt, double &rlambda_in_c2, double &rAA_c2, double &rlambda_out1, double &rlambda_out2, int* &rpipe_flow);

private:
	int UPSTREAM;
	int DOWNSTREAM;
	bool LEFT_TO_RIGHT;
	int TYPE; // ENLARGEMENT, CONTRACTION, or EQUAL_AREAS
	double area_ratio;
	double del1, del2, V1;// Variables for validation of sudden enlargement model
};

#endif // !defined(AFX_SUDDEN_H__FBBF483A_FE81_404E_87D9_C87F7E2D43B7__INCLUDED_)
