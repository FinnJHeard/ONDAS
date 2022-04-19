// TurbDataPt.h: interface for the CTurbDataPt class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TURBDATAPT_H__2E3EEA87_61AD_4D89_BD71_349864A6CD43__INCLUDED_)
#define AFX_TURBDATAPT_H__2E3EEA87_61AD_4D89_BD71_349864A6CD43__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Properties.h"

class CTurbDataPt  
{
public:
	CTurbDataPt();
	virtual ~CTurbDataPt();

	// Copy constructor
	CTurbDataPt(const CTurbDataPt& inTurbDataPt);

	CTurbDataPt& operator=(const CTurbDataPt& inTurbDataPt);	// Overload =
	CTurbDataPt operator+ (const CTurbDataPt& rhsTurbDataPt);	// Overload +
	CTurbDataPt operator- (const CTurbDataPt& rhsTurbDataPt);	// Overload -
	CTurbDataPt operator/ (const CTurbDataPt& denTurbDataPt);	// Overload /
	CTurbDataPt operator* (const CTurbDataPt& rhsTurbDataPt);	// Overload *
	CTurbDataPt operator* (const double value);

	void ProcessTurbDataPtVar(CProperties* pPpt, double* f, double* C, double T);

	
	// Basic data set - common
	// =======================	
	double sp;
	double mfp;
	double pr;
	double eta;
	
	// Additional data set - for constant outlet pressure
	// ==================================================
	double lambda_in_star;	// Inlet side (1), i.e. lambda_in1_star
	double lambda_out_star; // Inlet side (1), i.e. lambda_out1_star

	// Additional data set - for variable outlet pressure
	// ==================================================
	double G1;
	double p1_over_p2;
	double T2_over_T1;
	double lambda_in_star_ratio;
};

#endif // !defined(AFX_TURBDATAPT_H__2E3EEA87_61AD_4D89_BD71_349864A6CD43__INCLUDED_)
