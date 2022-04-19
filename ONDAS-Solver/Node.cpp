// Node.cpp: implementation of the CNode class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Node.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CNode::CNode()
{
	CL1 = new double [2];
	CL2 = new double [2];
	AA = new double [2];
	W = new double [3];
	F = new double [3];
	C = new double [3];
	W_pred = new double [3];
	F_pred = new double [3];
	C_pred = new double [3];
	S = new double [3];
	S_pred = new double [3];
	W_prev = new double [3];
}

CNode::~CNode()
{
	delete [] CL1;
	delete [] CL2;
	delete [] AA;
	delete [] W;
	delete [] F;
	delete [] C;
	delete [] W_pred;
	delete [] F_pred;
	delete [] C_pred;
	delete [] S;
	delete [] S_pred;
	delete [] W_prev;
}

// Copy constructor
CNode::CNode(const CNode& inNode)
{
	CL1 = new double [2];
	CL2 = new double [2];
	AA = new double [2];
	W = new double [3];
	F = new double [3];
	C = new double [3];
	W_pred = new double [3];
	F_pred = new double [3];
	C_pred = new double [3];
	S = new double [3];
	S_pred = new double [3];
	W_prev = new double [3];

	CL1[0] = inNode.CL1[0];	
	CL1[1] = inNode.CL1[1];
	CL2[0] = inNode.CL2[0];	
	CL2[1] = inNode.CL2[1];
	AA[0] = inNode.AA[0]; 
	AA[1] = inNode.AA[1];
	
	for(int k=0; k<3; ++k)
	{
		W[k] = inNode.W[k];
		F[k] = inNode.F[k];
		C[k] = inNode.C[k];
		W_pred[k] = inNode.W_pred[k];
		F_pred[k] = inNode.F_pred[k];
		C_pred[k] = inNode.C_pred[k];
		S[k] = inNode.S[k];
		S_pred[k] = inNode.S_pred[k];
		W_prev[k] = inNode.W_prev[k];
	}

	cfa = inNode.cfa;
	cfa_pred = inNode.cfa_pred;
	cfa_delx = inNode.cfa_delx;
	cfa_delx_pred = inNode.cfa_delx_pred;
	d = inNode.d;
	dddx = inNode.dddx;
	d2ddx2 = inNode.d2ddx2;
	d_pred = inNode.d_pred;
	dddx_pred = inNode.dddx_pred;
	d2ddx2_pred = inNode.d2ddx2_pred;
	DELX_R = inNode.DELX_R;
	DELX_L = inNode.DELX_L;
	dfdx = inNode.dfdx;
	dfdx_pred = inNode.dfdx_pred;
	f = inNode.f;
	f_dash = inNode.f_dash;
	f_pred = inNode.f_pred;
	f_prev = inNode.f_prev;
	vol = inNode.vol;
	vol_prev = inNode.vol_prev;
	x = inNode.x;
	X = inNode.X;

	// Primitives
	U = inNode.U;
	A = inNode.A;
	M = inNode.M;
	p_dash = inNode.p_dash;
	T = inNode.T;
	rho = inNode.rho;
	rho_prev = inNode.rho_prev;
	mdot = inNode.mdot;
	Re = inNode.Re;

	// Other
	bc = inNode.bc;
	side = inNode.side;
	rflow = inNode.rflow;

	STEADY = inNode.STEADY;

	// Haemodynamics
	B = inNode.B;
	B_pred = inNode.B_pred;
	b = inNode.b;
	b_pred = inNode.b_pred;
	dbdr0 = inNode.dbdr0;
	dbdr0_pred = inNode.dbdr0_pred;
	d0 = inNode.d0;
	d0_pred = inNode.d0_pred;
	dd0dx = inNode.dd0dx;
	dd0dx_pred = inNode.dd0dx_pred;
	f0 = inNode.f0;
	f0_pred = inNode.f0_pred;
	p0_dash = inNode.p0_dash;
}

CNode& CNode::operator=(const CNode& inNode)
{
	if(this != &inNode)
	{
		CL1[0] = inNode.CL1[0];	
		CL1[1] = inNode.CL1[1];
		CL2[0] = inNode.CL2[0];	
		CL2[1] = inNode.CL2[1];
		AA[0] = inNode.AA[0]; 
		AA[1] = inNode.AA[1];
		
		for(int k=0; k<3; ++k)
		{
			W[k] = inNode.W[k];
			F[k] = inNode.F[k];
			C[k] = inNode.C[k];
			W_pred[k] = inNode.W_pred[k];
			F_pred[k] = inNode.F_pred[k];
			C_pred[k] = inNode.C_pred[k];
			S[k] = inNode.S[k];
			S_pred[k] = inNode.S_pred[k];
			W_prev[k] = inNode.W_prev[k];
		}

		x = inNode.x;
		X = inNode.X;
		DELX_R = inNode.DELX_R;
		DELX_L = inNode.DELX_L;
		d = inNode.d;
		dddx = inNode.dddx;
		d2ddx2 = inNode.d2ddx2;
		d_pred = inNode.d_pred;
		dddx_pred = inNode.dddx_pred;
		d2ddx2_pred = inNode.d2ddx2_pred;
		f = inNode.f;
		f_pred = inNode.f_pred;
		cfa = inNode.cfa;
		cfa_pred = inNode.cfa_pred;
		cfa_delx = inNode.cfa_delx;
		cfa_delx_pred = inNode.cfa_delx_pred;
		f_prev = inNode.f_prev;
		f_dash = inNode.f_dash;
		dfdx = inNode.dfdx;
		dfdx_pred = inNode.dfdx_pred;
		vol = inNode.vol;
		vol_prev = inNode.vol_prev;

		// Primitives
		U = inNode.U;
		A = inNode.A;
		M = inNode.M;
		p_dash = inNode.p_dash;
		T = inNode.T;
		rho = inNode.rho;
		rho_prev = inNode.rho_prev;
		mdot = inNode.mdot;
		Re = inNode.Re;

		// Other
		bc = inNode.bc;
		side = inNode.side;
		rflow = inNode.rflow;

		STEADY = inNode.STEADY;

		// Haemodynamics
		B = inNode.B;
		B_pred = inNode.B_pred;
		b = inNode.b;
		b_pred = inNode.b_pred;
		dbdr0 = inNode.dbdr0;
		dbdr0_pred = inNode.dbdr0_pred;
		d0 = inNode.d0;
		d0_pred = inNode.d0_pred;
		dd0dx = inNode.dd0dx;
		dd0dx_pred = inNode.dd0dx_pred;
		f0 = inNode.f0;
		f0_pred = inNode.f0_pred;
		p0_dash = inNode.p0_dash;
	}
	return *this;
}

CNode CNode::operator+ (const CNode& inNode)
{
	// this + inNode
		
	CNode temp;
	temp = *this; // Set temp equal to this, then add rhs values

	temp.CL1[0] += inNode.CL1[0];
	temp.CL1[1] += inNode.CL1[1];
	temp.CL2[0] += inNode.CL2[0];
	temp.CL2[1] += inNode.CL2[1];
	temp.AA[0] += inNode.AA[0];
	temp.AA[1] += inNode.AA[1];

	for(int k=0; k<3; ++k)
	{
		temp.W[k] += inNode.W[k];
		temp.F[k] += inNode.F[k];
		temp.C[k] += inNode.C[k];
		temp.W_pred[k] += inNode.W_pred[k];
		temp.F_pred[k] += inNode.F_pred[k];
		temp.C_pred[k] += inNode.C_pred[k];
		temp.S[k] += inNode.S[k];
		temp.S_pred[k] += inNode.S_pred[k];
		temp.W_prev[k] += inNode.W_prev[k];
	}

	temp.x += inNode.x;
	temp.X += inNode.X;
	temp.DELX_R += inNode.DELX_R;
	temp.DELX_L += inNode.DELX_L;
	temp.d += inNode.d;
	temp.dddx += inNode.dddx;
	temp.d2ddx2 += inNode.d2ddx2;
	temp.d_pred += inNode.d_pred;
	temp.dddx_pred += inNode.dddx_pred;
	temp.d2ddx2_pred += inNode.d2ddx2_pred;
	temp.f += inNode.f;
	temp.f_pred += inNode.f_pred;
	temp.cfa += inNode.cfa;
	temp.cfa_pred += inNode.cfa_pred;
	temp.cfa_delx += inNode.cfa_delx;
	temp.cfa_delx_pred += inNode.cfa_delx_pred;
	temp.f_prev += inNode.f_prev;
	temp.f_dash += inNode.f_dash;
	temp.dfdx += inNode.dfdx;
	temp.dfdx_pred += inNode.dfdx_pred;
	temp.vol += inNode.vol;
	temp.vol_prev += inNode.vol_prev;

	// Primitives
	temp.U += inNode.U;
	temp.A += inNode.A;
	temp.M += inNode.M;
	temp.p_dash += inNode.p_dash;
	temp.T += inNode.T;
	temp.rho += inNode.rho;
	temp.rho_prev += inNode.rho_prev;
	temp.mdot += inNode.mdot;
//	temp.Re += inNode.Re; // Calculate Re from interpolated properties

	// Other
	// N/A: bc, side, rflow
	
	// Haemodynamics
	temp.B += inNode.B;
	temp.B_pred += inNode.B_pred;
	temp.b += inNode.b;
	temp.b_pred += inNode.b_pred;
	temp.dbdr0 += inNode.dbdr0;
	temp.dbdr0_pred += inNode.dbdr0_pred;
	temp.d0 += inNode.d0;
	temp.d0_pred += inNode.d0_pred;
	temp.dd0dx += inNode.dd0dx;
	temp.dd0dx_pred += inNode.dd0dx_pred;
	temp.f0 += inNode.f0;
	temp.f0_pred += inNode.f0_pred;
	temp.p0_dash += inNode.p0_dash;

	return temp;
}

CNode CNode::operator- (const CNode& inNode)
{
	// this - inNode
		
	CNode temp;
	temp = *this; // Set temp equal to this, then subtract rhs values

	temp.CL1[0] -= inNode.CL1[0];
	temp.CL1[1] -= inNode.CL1[1];
	temp.CL2[0] -= inNode.CL2[0];
	temp.CL2[1] -= inNode.CL2[1];
	temp.AA[0] -= inNode.AA[0];
	temp.AA[1] -= inNode.AA[1];

	for(int k=0; k<3; ++k)
	{
		temp.W[k] -= inNode.W[k];
		temp.F[k] -= inNode.F[k];
		temp.C[k] -= inNode.C[k];
		temp.W_pred[k] -= inNode.W_pred[k];
		temp.F_pred[k] -= inNode.F_pred[k];
		temp.C_pred[k] -= inNode.C_pred[k];
		temp.S[k] -= inNode.S[k];
		temp.S_pred[k] -= inNode.S_pred[k];
		temp.W_prev[k] -= inNode.W_prev[k];
	}

	temp.x -= inNode.x;
	temp.X -= inNode.X;
	temp.DELX_R -= inNode.DELX_R;
	temp.DELX_L -= inNode.DELX_L;
	temp.d -= inNode.d;
	temp.dddx -= inNode.dddx;
	temp.d2ddx2 -= inNode.d2ddx2;
	temp.d_pred -= inNode.d_pred;
	temp.dddx_pred -= inNode.dddx_pred;
	temp.d2ddx2_pred -= inNode.d2ddx2_pred;
	temp.f -= inNode.f;
	temp.f_pred -= inNode.f_pred;
	temp.cfa -= inNode.cfa;
	temp.cfa_pred -= inNode.cfa_pred;
	temp.cfa_delx -= inNode.cfa_delx;
	temp.cfa_delx_pred -= inNode.cfa_delx_pred;
	temp.f_prev -= inNode.f_prev;
	temp.f_dash -= inNode.f_dash;
	temp.dfdx -= inNode.dfdx;
	temp.dfdx_pred -= inNode.dfdx_pred;
	temp.vol -= inNode.vol;
	temp.vol_prev -= inNode.vol_prev;

	// Primitives
	temp.U -= inNode.U;
	temp.A -= inNode.A;
	temp.M -= inNode.M;
	temp.p_dash -= inNode.p_dash;
	temp.T -= inNode.T;
	temp.rho -= inNode.rho;
	temp.rho_prev -= inNode.rho_prev;
	temp.mdot -= inNode.mdot;
//	temp.Re -= inNode.Re;

	// Other
	// N/A: bc, side, rflow

	// Haemodynamics
	temp.B -= inNode.B;
	temp.B_pred -= inNode.B_pred;
	temp.b -= inNode.b;
	temp.b_pred -= inNode.b_pred;
	temp.dbdr0 -= inNode.dbdr0;
	temp.dbdr0_pred -= inNode.dbdr0_pred;
	temp.d0 -= inNode.d0;
	temp.d0_pred -= inNode.d0_pred;
	temp.dd0dx -= inNode.dd0dx;
	temp.dd0dx_pred -= inNode.dd0dx_pred;
	temp.f0 -= inNode.f0;
	temp.f0_pred -= inNode.f0_pred;
	temp.p0_dash -= inNode.p0_dash;
	
	return temp;
}

CNode CNode::operator/ (const CNode& inNode)
{
	// this / inNode
		
	CNode temp;
	temp = *this; // Set temp equal to this, then divide by denominator values

	temp.CL1[0] /= inNode.CL1[0];
	temp.CL1[1] /= inNode.CL1[1];
	temp.CL2[0] /= inNode.CL2[0];
	temp.CL2[1] /= inNode.CL2[1];
	temp.AA[0] /= inNode.AA[0];
	temp.AA[1] /= inNode.AA[1];

	for(int k=0; k<3; ++k)
	{
		temp.W[k] /= inNode.W[k];
		temp.F[k] /= inNode.F[k];
		temp.C[k] /= inNode.C[k];
		temp.W_pred[k] /= inNode.W_pred[k];
		temp.F_pred[k] /= inNode.F_pred[k];
		temp.C_pred[k] /= inNode.C_pred[k];
		temp.S[k] /= inNode.S[k];
		temp.S_pred[k] /= inNode.S_pred[k];
		temp.W_prev[k] /= inNode.W_prev[k];
	}

	temp.x /= inNode.x;
	temp.X /= inNode.X;
	temp.DELX_R /= inNode.DELX_R;
	temp.DELX_L /= inNode.DELX_L;
	temp.d /= inNode.d;
	temp.dddx /= inNode.dddx;
	temp.d2ddx2 /= inNode.d2ddx2;
	temp.d_pred /= inNode.d_pred;
	temp.dddx_pred /= inNode.dddx_pred;
	temp.d2ddx2_pred /= inNode.d2ddx2_pred;
	temp.f /= inNode.f;
	temp.f_pred /= inNode.f_pred;
	temp.cfa /= inNode.cfa;
	temp.cfa_pred /= inNode.cfa_pred;
	temp.cfa_delx /= inNode.cfa_delx;
	temp.cfa_delx_pred /= inNode.cfa_delx_pred;
	temp.f_prev /= inNode.f_prev;
	temp.f_dash /= inNode.f_dash;
	temp.dfdx /= inNode.dfdx;
	temp.dfdx_pred /= inNode.dfdx_pred;
	temp.vol /= inNode.vol;
	temp.vol_prev /= inNode.vol_prev;

	// Primitives
	temp.U /= inNode.U;
	temp.A /= inNode.A;
	temp.M /= inNode.M;
	temp.p_dash /= inNode.p_dash;
	temp.T /= inNode.T;
	temp.rho /= inNode.rho;
	temp.rho_prev /= inNode.rho_prev;
	temp.mdot /= inNode.mdot;
//	temp.Re /= inNode.Re;
	
	// Other
	// N/A: bc, side, rflow
	
	// Haemodynamics
	temp.B /= inNode.B;
	temp.B_pred /= inNode.B_pred;
	temp.b /= inNode.b;
	temp.b_pred /= inNode.b_pred;
	temp.dbdr0 /= inNode.dbdr0;
	temp.dbdr0_pred /= inNode.dbdr0_pred;
	temp.d0 /= inNode.d0;
	temp.d0_pred /= inNode.d0_pred;
	temp.dd0dx /= inNode.dd0dx;
	temp.dd0dx_pred /= inNode.dd0dx_pred;
	temp.f0 /= inNode.f0;
	temp.f0_pred /= inNode.f0_pred;
	temp.p0_dash /= inNode.p0_dash;


	return temp;
}

CNode CNode::operator* (const CNode& inNode)
{
	// this * inNode
		
	CNode temp;
	temp = *this; // Set temp equal to this, then multiply by rhs values

	temp.CL1[0] *= inNode.CL1[0];
	temp.CL1[1] *= inNode.CL1[1];
	temp.CL2[0] *= inNode.CL2[0];
	temp.CL2[1] *= inNode.CL2[1];
	temp.AA[0] *= inNode.AA[0];
	temp.AA[1] *= inNode.AA[1];

	for(int k=0; k<3; ++k)
	{
		temp.W[k] *= inNode.W[k];
		temp.F[k] *= inNode.F[k];
		temp.C[k] *= inNode.C[k];
		temp.W_pred[k] *= inNode.W_pred[k];
		temp.F_pred[k] *= inNode.F_pred[k];
		temp.C_pred[k] *= inNode.C_pred[k];
		temp.S[k] *= inNode.S[k];
		temp.S_pred[k] *= inNode.S_pred[k];
		temp.W_prev[k] *= inNode.W_prev[k];
	}

	temp.x *= inNode.x;
	temp.X *= inNode.X;
	temp.DELX_R *= inNode.DELX_R;
	temp.DELX_L *= inNode.DELX_L;
	temp.d *= inNode.d;
	temp.dddx *= inNode.dddx;
	temp.d2ddx2 *= inNode.d2ddx2;
	temp.d_pred *= inNode.d_pred;
	temp.dddx_pred *= inNode.dddx_pred;
	temp.d2ddx2_pred *= inNode.d2ddx2_pred;
	temp.f *= inNode.f;
	temp.f_pred *= inNode.f_pred;
	temp.cfa *= inNode.cfa;
	temp.cfa_pred *= inNode.cfa_pred;
	temp.cfa_delx *= inNode.cfa_delx;
	temp.cfa_delx_pred *= inNode.cfa_delx_pred;
	temp.f_prev *= inNode.f_prev;
	temp.f_dash *= inNode.f_dash;
	temp.dfdx *= inNode.dfdx;
	temp.dfdx_pred *= inNode.dfdx_pred;
	temp.vol *= inNode.vol;
	temp.vol_prev *= inNode.vol_prev;

	// Primitives
	temp.U *= inNode.U;
	temp.A *= inNode.A;
	temp.M *= inNode.M;
	temp.p_dash *= inNode.p_dash;
	temp.T *= inNode.T;
	temp.rho *= inNode.rho;
	temp.rho_prev *= inNode.rho_prev;
	temp.mdot *= inNode.mdot;
//	temp.Re *= inNode.Re;

	// Other
	// N/A: bc, side, rflow
	
	// Haemodynamics
	temp.B *= inNode.B;
	temp.B_pred *= inNode.B_pred;
	temp.b *= inNode.b;
	temp.b_pred *= inNode.b_pred;
	temp.dbdr0 *= inNode.dbdr0;
	temp.dbdr0_pred *= inNode.dbdr0_pred;
	temp.d0 *= inNode.d0;
	temp.d0_pred *= inNode.d0_pred;
	temp.dd0dx *= inNode.dd0dx;
	temp.dd0dx_pred *= inNode.dd0dx_pred;
	temp.f0 *= inNode.f0;
	temp.f0_pred *= inNode.f0_pred;
	temp.p0_dash *= inNode.p0_dash;

	return temp;
}

CNode CNode::operator* (const double value)
{
	// this * some double
		
	CNode temp;
	temp = *this; // Set temp equal to this, then multiply by rhs values

	temp.CL1[0] *= value;
	temp.CL1[1] *= value;
	temp.CL2[0] *= value;
	temp.CL2[1] *= value;
	temp.AA[0] *= value;
	temp.AA[1] *= value;

	for(int k=0; k<3; ++k)
	{
		temp.W[k] *= value;
		temp.F[k] *= value;
		temp.C[k] *= value;
		temp.W_pred[k] *= value;
		temp.F_pred[k] *= value;
		temp.C_pred[k] *= value;
		temp.S[k] *= value;
		temp.S_pred[k] *= value;
		temp.W_prev[k] *= value;
	}

	temp.x *= value;
	temp.X *= value;
	temp.DELX_R *= value;
	temp.DELX_L *= value;
	temp.d *= value;
	temp.dddx *= value;
	temp.d2ddx2 *= value;
	temp.d_pred *= value;
	temp.dddx_pred *= value;
	temp.d2ddx2_pred *= value;
	temp.f *= value;
	temp.f_pred *= value;
	temp.cfa *= value;
	temp.cfa_pred *= value;
	temp.cfa_delx *= value;
	temp.cfa_delx_pred *= value;
	temp.f_prev *= value;
	temp.f_dash *= value;
	temp.dfdx *= value;
	temp.dfdx_pred *= value;
	temp.vol *= value;
	temp.vol_prev *= value;

	// Primitives
	temp.U *= value;
	temp.A *= value;
	temp.M *= value;
	temp.p_dash *= value;
	temp.T *= value;
	temp.rho *= value;
	temp.rho_prev *= value;
	temp.mdot *= value;
//	temp.Re *= value;

	// Other
	// N/A: bc, side, rflow

	// Haemodynamics
	temp.B *= value;
	temp.B_pred *= value;
	temp.b *= value;
	temp.b_pred *= value;
	temp.dbdr0 *= value;
	temp.dbdr0_pred *= value;
	temp.d0 *= value;
	temp.d0_pred *= value;
	temp.dd0dx *= value;
	temp.dd0dx_pred *= value;
	temp.f0 *= value;
	temp.f0_pred *= value;
	temp.p0_dash *= value;
	
	return temp;
}

void CNode::Print()
{
	int k;
	cout << endl;

	cout << "CNode::Print()" << endl;
	cout << "CL1[0] = " << CL1[0] << endl;
	cout << "CL2[0] = " << CL2[0] << endl;
	cout << "CL1[0] = " << CL1[0] << endl;
	cout << "CL2[1] = " << CL2[1] << endl;
	cout << "AA[0] = " << AA[0] << endl;
	cout << "AA[1] = " << AA[1] << endl;

	for(k=0; k<3; ++k) cout << "W[" << k << "] = " << W[k] << endl;
	for(k=0; k<3; ++k) cout << "F[" << k << "] = " << F[k] << endl;
	for(k=0; k<3; ++k) cout << "C[" << k << "] = " << C[k] << endl;
	for(k=0; k<3; ++k) cout << "W_pred[" << k << "] = " << W_pred[k] << endl;
	for(k=0; k<3; ++k) cout << "F_pred[" << k << "] = " << F_pred[k] << endl;
	for(k=0; k<3; ++k) cout << "C_pred[" << k << "] = " << C_pred[k] << endl;
	for(k=0; k<3; ++k) cout << "S[" << k << "] = " << S[k] << endl;
	for(k=0; k<3; ++k) cout << "S_pred[" << k << "] = " << S_pred[k] << endl;
	for(k=0; k<3; ++k) cout << "W_prev[" << k << "] = " << W_prev[k] << endl;

	cout << "x = " << x << endl;
	cout << "X = " << X << endl;
	cout << "DELX_R = " << DELX_R << endl;
	cout << "DELX_L = " << DELX_L << endl;
	cout << "d = " << d << endl;
	cout << "dddx = " << dddx << endl;
	cout << "d2ddx2 = " << d2ddx2 << endl;
	cout << "d_pred = " << d_pred << endl;
	cout << "ddx_pred = " << dddx_pred << endl;
	cout << "d2ddx2_pred = " << d2ddx2_pred << endl;
	cout << "f = " << f << endl;
	cout << "cfa = " << cfa << endl;
	cout << "cfa_pred = " << cfa_pred << endl;
	cout << "cfa_delx = " << cfa_delx << endl;
	cout << "cfa_delx_pred = " << cfa_delx_pred << endl;
	cout << "f_pred = " << f_pred << endl;
	cout << "f_prev = " << f_prev << endl;
	cout << "f_dash = " << f_dash << endl;
	cout << "dfdx = " << dfdx << endl;
	cout << "dfdx_pred = " << dfdx_pred << endl;
	cout << "vol = " << vol << endl;
	cout << "vol_prev = " << vol_prev << endl;

	// Primitives
	cout << "U = " << U << endl;
	cout << "A = " << A << endl;
	cout << "M = " << M << endl;
	cout << "p_dash = " << p_dash << endl;
	cout << "T = " << T << endl;
	cout << "rho = " << rho << endl;
	cout << "rho_prev = " << rho_prev << endl;
	cout << "mdot = " << mdot << endl;
//	cout << "Re = " << Re << endl;
	cout << endl;

	// Haemodynamics
	cout << "B = " << B << endl;
	cout << "B_pred = " << B_pred << endl;
	cout << "d0 = " << d0 << endl;
	cout << "d0_pred = " << d0_pred << endl;
	cout << "dd0dx = " << dd0dx << endl;
	cout << "dd0dx_pred = " << dd0dx_pred << endl;
	cout << "f0 = " << f0 << endl;
	cout << "f0_pred = " << f0_pred << endl;
	cout << "p0_dash = " << p0_dash << endl;
	cout << "b = " << b << endl;
	cout << "b_pred = " << b_pred << endl;
	cout << "dbdr0 = " << dbdr0 << endl;
	cout << "dbdr0_pred = " << dbdr0_pred << endl;
}