// DataPoint.h: interface for the CDataPoint class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DATAPOINT_H__E2ED6DF7_A5FA_4AEF_A3C9_9ADAE661082C__INCLUDED_)
#define AFX_DATAPOINT_H__E2ED6DF7_A5FA_4AEF_A3C9_9ADAE661082C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Properties.h"

class CDataPoint  
{
public:
	CDataPoint();
	virtual ~CDataPoint();

	// Copy constructor
	CDataPoint(const CDataPoint& inDataPt);

	CDataPoint& operator=(const CDataPoint& inDataPt);	// Overload =
	CDataPoint operator+ (const CDataPoint& rhsDataPt);	// Overload +
	CDataPoint operator- (const CDataPoint& rhsDataPt);	// Overload -
	CDataPoint operator/ (const CDataPoint& denDataPt);	// Overload /
	CDataPoint operator* (const CDataPoint& rhsDataPt);	// Overload *
	CDataPoint operator* (const double value);

	void ProcessDataPtVar(CProperties* pPpt, double* f, double* C, bool TURBINE, double T);
	void PrintToScreen(CProperties* pPpt, bool VAR);
	
	// Basic data set - common
	// =======================	
	double sp;
	double pr;
	double mfp;
	double eta;
	
	// Additional data set - for constant turbine outlet/compressor inlet pressure
	// ===========================================================================
	double lambda_in_star;	// Outlet side (2), i.e. lambda_in2_star
	double lambda_out_star;	// Outlet side (2), i.e. lambda_out2_star
	
	// The following are the same for all points on a curve
	double lambda_in_star_min;
	double lambda_in_star_max;
	double Ks;

	// Compressor:
	double T2_over_T1;
	double M2;
	double AA2_over_AA1;

	// The following are the same for all points on a curve
	double lambda_in_star_max_comp;
	double lambda_in_star_min_comp; 
//	double M2_max;
	double KS_comp;
	
	// Additional data set - for variable turbine outlet/compressor inlet pressure
	// ===========================================================================
	double lambda_in_star_ratio;	// = lambda_in_star1_over_lambda_in_star2;
	double G1;						// Mass flow parameter G1
	double p1_over_p2;
//	double T2_over_T1; // Defined above

	// Variable compressor power parameters
	// ====================================
	double power_param;	// W_CI/(p1.sqrt(T1))
	double torque_param;// L_CI/p1

	// Common
	// ======
	double query_param;	// Identical to lambda_in_star_ratio (variable) or lambda_in_star (constant)
};

#endif // !defined(AFX_DATAPOINT_H__E2ED6DF7_A5FA_4AEF_A3C9_9ADAE661082C__INCLUDED_)
