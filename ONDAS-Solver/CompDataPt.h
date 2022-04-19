// CompDataPt.h: interface for the CCompDataPt class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_COMPDATAPT_H__D82F43AB_72C3_4EA3_9CC9_1DFE8AFDD0BF__INCLUDED_)
#define AFX_COMPDATAPT_H__D82F43AB_72C3_4EA3_9CC9_1DFE8AFDD0BF__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Properties.h"

class CCompDataPt  
{
public:
	CCompDataPt();
	virtual ~CCompDataPt();

	// Copy constructor
	CCompDataPt(const CCompDataPt& inCompDataPt);

	CCompDataPt& operator=(const CCompDataPt& inCompDataPt);	// Overload =
	CCompDataPt operator+ (const CCompDataPt& rhsCompDataPt);	// Overload +
	CCompDataPt operator- (const CCompDataPt& rhsCompDataPt);	// Overload -
	CCompDataPt operator/ (const CCompDataPt& denCompDataPt);	// Overload /
	CCompDataPt operator* (const CCompDataPt& rhsCompDataPt);	// Overload *
	CCompDataPt operator* (const double value);

	void ProcessCompDataPtVar(CProperties* pPpt, double* f, double* C, double T);
	
	// Basic data set - common
	// =======================	
	double sp;
	double pr;
	double mfp;
	double eta;
	
	// Additional data set - for constant inlet pressure
	// =================================================
	double T2_over_T1;
	double M2;
	double lambda_in_star;	// Outlet side (2), i.e. lambda_in2_star
	double lambda_out_star;	// Outlet side (2), i.e. lambda_out2_star
	double AA2_over_AA1;

	// The following are the same for all points on a curve
	double lambda_in_star_max_comp;
//	double M2_max;
	double KS_comp;
	double lambda_in_star_min_comp; 

	// Additional data set - for variable inlet pressure - same as for turbine
	// =======================================================================
	double lambda_in_star_ratio;	// = lambda_in_star1_over_lambda_in_star2;
	double G1;						// Mass flow parameter G1
	double p1_over_p2;
//	double T2_over_T1; // Defined above
};

#endif // !defined(AFX_COMPDATAPT_H__D82F43AB_72C3_4EA3_9CC9_1DFE8AFDD0BF__INCLUDED_)
