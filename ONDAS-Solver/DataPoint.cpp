// DataPoint.cpp: implementation of the CDataPoint class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "DataPoint.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CDataPoint::CDataPoint()
{

}

CDataPoint::~CDataPoint()
{

}

// Copy constructor
CDataPoint::CDataPoint(const CDataPoint& inDataPt)
{
	sp = inDataPt.sp;
	mfp = inDataPt.mfp;
	pr = inDataPt.pr;
	eta = inDataPt.eta;

	power_param = inDataPt.power_param;
	torque_param = inDataPt.torque_param;

	T2_over_T1 = inDataPt.T2_over_T1;

	M2 = inDataPt.M2;
	lambda_in_star = inDataPt.lambda_in_star;
	lambda_out_star = inDataPt.lambda_out_star;
	AA2_over_AA1 = inDataPt.AA2_over_AA1;
	
	G1 = inDataPt.G1;
	p1_over_p2 = inDataPt.p1_over_p2;
	lambda_in_star_ratio = inDataPt.lambda_in_star_ratio;

	lambda_in_star_min = inDataPt.lambda_in_star_min;
	lambda_in_star_max = inDataPt.lambda_in_star_max;
	Ks = inDataPt.Ks;

	query_param = inDataPt.query_param;
}

CDataPoint& CDataPoint::operator=(const CDataPoint& inDataPt)
{
	if(this != &inDataPt)
	{
		sp = inDataPt.sp;
		mfp = inDataPt.mfp;
		pr = inDataPt.pr;
		eta = inDataPt.eta;

		power_param = inDataPt.power_param;
		torque_param = inDataPt.torque_param;

		T2_over_T1 = inDataPt.T2_over_T1;

		M2 = inDataPt.M2;
		lambda_in_star = inDataPt.lambda_in_star;
		lambda_out_star = inDataPt.lambda_out_star;
		AA2_over_AA1 = inDataPt.AA2_over_AA1;

		G1 = inDataPt.G1;
		p1_over_p2 = inDataPt.p1_over_p2;
		lambda_in_star_ratio = inDataPt.lambda_in_star_ratio;

		lambda_in_star_min = inDataPt.lambda_in_star_min;
		lambda_in_star_max = inDataPt.lambda_in_star_max;
		Ks = inDataPt.Ks;

		query_param = inDataPt.query_param;
	}
	return *this;
}

CDataPoint CDataPoint::operator+ (const CDataPoint& rhsDataPt)
{
	// this + rhsDataPt
		
	CDataPoint temp;
	temp = *this; // Set temp equal to this, then add rhs values

	temp.sp += rhsDataPt.sp;
	temp.mfp += rhsDataPt.mfp;
	temp.pr += rhsDataPt.pr;
	temp.eta += rhsDataPt.eta;

	temp.power_param += rhsDataPt.power_param;
	temp.torque_param += rhsDataPt.torque_param;

	temp.T2_over_T1 += rhsDataPt.T2_over_T1;

	temp.M2 += rhsDataPt.M2;
	temp.lambda_in_star += rhsDataPt.lambda_in_star;
	temp.lambda_out_star += rhsDataPt.lambda_out_star;
	temp.AA2_over_AA1 += rhsDataPt.AA2_over_AA1;

	temp.G1 += rhsDataPt.G1;
	temp.p1_over_p2 += rhsDataPt.p1_over_p2;
	temp.lambda_in_star_ratio += rhsDataPt.lambda_in_star_ratio;

	temp.lambda_in_star_min += rhsDataPt.lambda_in_star_min;
	temp.lambda_in_star_max += rhsDataPt.lambda_in_star_max;
	temp.Ks += rhsDataPt.Ks;
	
	temp.query_param += rhsDataPt.query_param;

	return temp;
}

CDataPoint CDataPoint::operator- (const CDataPoint& rhsDataPt)
{
	// this - rhsDataPt
		
	CDataPoint temp;
	temp = *this; // Set temp equal to this, then subtract rhs values

	temp.sp -= rhsDataPt.sp;
	temp.mfp -= rhsDataPt.mfp;
	temp.pr -= rhsDataPt.pr;
	temp.eta -= rhsDataPt.eta;

	temp.power_param -= rhsDataPt.power_param;
	temp.torque_param -= rhsDataPt.torque_param;

	temp.T2_over_T1 -= rhsDataPt.T2_over_T1;

	temp.M2 -= rhsDataPt.M2;
	temp.lambda_in_star -= rhsDataPt.lambda_in_star;
	temp.lambda_out_star -= rhsDataPt.lambda_out_star;
	temp.AA2_over_AA1 -= rhsDataPt.AA2_over_AA1;

	temp.G1 -= rhsDataPt.G1;
	temp.p1_over_p2 -= rhsDataPt.p1_over_p2;
	temp.lambda_in_star_ratio -= rhsDataPt.lambda_in_star_ratio;

	temp.lambda_in_star_min -= rhsDataPt.lambda_in_star_min;
	temp.lambda_in_star_max -= rhsDataPt.lambda_in_star_max;
	temp.Ks -= rhsDataPt.Ks;

	temp.query_param -= rhsDataPt.query_param;

	return temp;
}

CDataPoint CDataPoint::operator/ (const CDataPoint& denDataPt)
{
	// this / denDataPt
		
	CDataPoint temp;
	temp = *this; // Set temp equal to this, then divide by denominator values

	temp.sp /= denDataPt.sp;
	temp.mfp /= denDataPt.mfp;
	temp.pr /= denDataPt.pr;
	temp.eta /= denDataPt.eta;

	temp.power_param /= denDataPt.power_param;
	temp.torque_param /= denDataPt.torque_param;

	temp.T2_over_T1 /= denDataPt.T2_over_T1;

	temp.M2 /= denDataPt.M2;
	temp.lambda_in_star /= denDataPt.lambda_in_star;
	temp.lambda_out_star /= denDataPt.lambda_out_star;
	temp.AA2_over_AA1 /= denDataPt.AA2_over_AA1;

	temp.G1 /= denDataPt.G1;
	temp.p1_over_p2 /= denDataPt.p1_over_p2;
	temp.lambda_in_star_ratio /= denDataPt.lambda_in_star_ratio;

	temp.lambda_in_star_min /= denDataPt.lambda_in_star_min;
	temp.lambda_in_star_max /= denDataPt.lambda_in_star_max;
	temp.Ks /= denDataPt.Ks;

	temp.query_param /= denDataPt.query_param;

	return temp;
}

CDataPoint CDataPoint::operator* (const CDataPoint& rhsDataPt)
{
	// this * rhsDataPt
		
	CDataPoint temp;
	temp = *this; // Set temp equal to this, then multiply by rhs values

	temp.sp *= rhsDataPt.sp;
	temp.mfp *= rhsDataPt.mfp;
	temp.pr *= rhsDataPt.pr;
	temp.eta *= rhsDataPt.eta;

	temp.power_param *= rhsDataPt.power_param;
	temp.torque_param *= rhsDataPt.torque_param;

	temp.T2_over_T1 *= rhsDataPt.T2_over_T1;

	temp.M2 *= rhsDataPt.M2;
	temp.lambda_in_star *= rhsDataPt.lambda_in_star;
	temp.lambda_out_star *= rhsDataPt.lambda_out_star;
	temp.AA2_over_AA1 *= rhsDataPt.AA2_over_AA1;

	temp.G1 *= rhsDataPt.G1;
	temp.p1_over_p2 *= rhsDataPt.p1_over_p2;
	temp.lambda_in_star_ratio *= rhsDataPt.lambda_in_star_ratio;

	temp.lambda_in_star_min *= rhsDataPt.lambda_in_star_min;
	temp.lambda_in_star_max *= rhsDataPt.lambda_in_star_max;
	temp.Ks *= rhsDataPt.Ks;
	
	temp.query_param *= rhsDataPt.query_param;

	return temp;
}

CDataPoint CDataPoint::operator* (const double value)
{
	// this * some double
		
	CDataPoint temp;
	temp = *this; // Set temp equal to this, then multiply by rhs values

	temp.sp *= value;
	temp.mfp *= value;
	temp.pr *= value;
	temp.eta *= value;

	temp.power_param *= value;
	temp.torque_param *= value;

	temp.T2_over_T1 *= value;
	
	temp.M2 *= value;
	temp.lambda_in_star *= value;
	temp.lambda_out_star *= value;
	temp.AA2_over_AA1 *= value;

	temp.G1 *= value;
	temp.p1_over_p2 *= value;
	temp.lambda_in_star_ratio *= value;

	temp.lambda_in_star_min *= value;
	temp.lambda_in_star_max *= value;
	temp.Ks *= value;

	temp.query_param *= value;

	return temp;
}

void CDataPoint::ProcessDataPtVar(CProperties* pPpt, double* f, double* C, bool TURBINE, double T)
// ============================================================ //
// Processes pr and mfp data into operating parameters			//
// for variable outlet pressure turbine/inlet pressure			//
// compressor maps.												//
// This function is called from									//
// CTurbocharger::ProcessMapVarComp/ProcessMapVarTurb			//
// ============================================================ //
{
	double c, M1, fM1, del_M1;
	double p01_over_p2, p01_over_p1, T01_over_T1, T2_over_T01;

	int UPSTREAM = 0;
	int DOWNSTREAM = 1;

	if(TURBINE)	p01_over_p2 = pr;
	else p01_over_p2 = 1/pr;

	if(mfp<=0.0) // Zero or reverse flow
	{
		G1 = 0;
		T2_over_T1 = 1.0;
		if(fabs(mfp) < 1e-6)
		{
			p1_over_p2 = p01_over_p2;
			lambda_in_star_ratio = pow(p1_over_p2, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)))
											* ( (1 + C[UPSTREAM]*G1)/(1 - C[DOWNSTREAM]*G1*p1_over_p2*sqrt(T2_over_T1)) );
		}
		/*
		if(fabs(G1)<1e-6)
		{
			cout << "G1 = " << G1 << endl;
			cout << "lambda_in_star_ratio = " << lambda_in_star_ratio << endl;
		}
		*/
	}
	else if (mfp>0.0)
	{
		c = pow( ( (mfp/f[UPSTREAM])*sqrt(pPpt->R_air/pPpt->gammaAir(T)) ), (2*(pPpt->gammaAir(T)-1)/(pPpt->gammaAir(T)+1)) );
		M1 = 0.001;

		do
		{
			fM1 = pow( c + ((pPpt->gammaAir(T)-1)/2)*c*pow(M1,2), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1)) );
			del_M1 = /*fabs*/(M1 - fM1);
			M1 = fM1;
		}
		//while(del_M1>0.0001);
		while(fabs(del_M1)>0.0001);

		p01_over_p1 = pow( 1 + ((pPpt->gammaAir(T)-1)/2)*pow(M1,2), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1) );
		T01_over_T1 = 1 + ((pPpt->gammaAir(T)-1)/2)*pow(M1,2);
		T2_over_T01 = 1 - eta*(1 - pow(1/p01_over_p2, (pPpt->gammaAir(T)-1)/pPpt->gammaAir(T)));

		//if(INLET_STAG) G1 = mfp*p01_over_p1*sqrt(1/T01_over_T1);
		//else G1 = mfp;
		G1 = mfp*p01_over_p1*sqrt(1/T01_over_T1);
//cout << "mfp = " << mfp << endl;
//cout << "p01_over_p1 = " << p01_over_p1 << endl;
//cout << "T01_over_T1 = " << T01_over_T1 << endl;
//cout << "G1 = " << G1 << endl;
		p1_over_p2 = (1/p01_over_p1)*p01_over_p2;
		T2_over_T1 = T2_over_T01*T01_over_T1;

		lambda_in_star_ratio = pow(p1_over_p2, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)))
										* ( (1 + C[UPSTREAM]*G1)/(1 - C[DOWNSTREAM]*G1*p1_over_p2*sqrt(T2_over_T1)) );
	}

	// In all cases:
	query_param = this->lambda_in_star_ratio;
}

void CDataPoint::PrintToScreen(CProperties* pPpt, bool VAR)
// ============================================================ //
// ============================================================ //
{
	cout << "Data point:\n";
	cout << "sp\t\t\t\t=\t" << sp << "\n";
	cout << "pr\t\t\t\t=\t" << pr << "\n";
	cout << "mfp\t\t\t\t=\t" << mfp << "\n";
	cout << "eta\t\t\t\t=\t" << eta << "\n";
	cout << "power_param\t\t\t=\t" << eta << "\n";
	cout << "torque_param\t\t\t=\t" << eta << "\n";
	cout << "\n";
	if(VAR)
	{
		cout << "lambda_in_star_ratio\t\t=\t" << lambda_in_star_ratio << "\n";
		cout << "G1\t\t\t\t=\t" << G1 << "\n";
		cout << "p1_over_p2\t\t\t=\t" << p1_over_p2 << "\n";
		cout << "T2_over_T1\t\t\t=\t" << T2_over_T1 << "\n";
	}
	else
	{
//		cout << "T2_over_T1\t\t\t=\t" << T2_over_T1 << "\n";
//		cout << "M2\t\t\t\t=\t" << M2 << "\n";
		cout << "lambda_in_star\t\t\t=\t" << lambda_in_star << "\n";
		cout << "lambda_out_star\t\t\t=\t" << lambda_out_star << "\n";
		cout << endl;
		cout << "lambda_in_star_min\t\t=\t" << lambda_in_star_min << "\n";
		cout << "lambda_in_star_max\t\t=\t" << lambda_in_star_max << "\n";
		cout << "Ks\t\t\t\t=\t" << Ks << "\n";
//		cout << "AA2_over_AA1\t\t\t=\t" << AA2_over_AA1 << "\n";
//		cout << "\n";
	}
	cout << "query_param\t\t\t=\t" << query_param << "\n";
}
