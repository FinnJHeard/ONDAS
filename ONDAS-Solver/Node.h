// Node.h: interface for the CNode class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_NODE_H__95E9F4DB_1C9D_4626_AE7E_0D043583F5CD__INCLUDED_)
#define AFX_NODE_H__95E9F4DB_1C9D_4626_AE7E_0D043583F5CD__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CNode  
{
public:
	CNode();
	virtual ~CNode();

	// Copy constructor
	CNode(const CNode& inNode);

	CNode& operator=(const CNode& inNode);	// Overload =
	CNode operator+ (const CNode& rhsNode);	// Overload +
	CNode operator- (const CNode& rhsNode);	// Overload -
	CNode operator/ (const CNode& denNode);	// Overload /
	CNode operator* (const CNode& rhsNode);	// Overload *
	CNode operator* (const double value);

	void Print();

private:

public:

	// Configuration variables (common)
	// ----------------------------------------------------------------------------------------------------
	int bc;					// Record boundary attachement type (INTERIOR for non-boundary nodes)
	int side;				// LEFT or RIGHT if a boundary node, INSIDE otherwise
	int rflow;				// If a junction node, an alias of its pipe odd_end_flow or even_end_flow
	
	// Geometry variables (common)
	// ----------------------------------------------------------------------------------------------------
	double cfa;				// Cell frontal area for this node (m^2)
	double cfa_pred;		// Predicted cell frontal area for this node (m^2)
	double cfa_delx;		// Cell length for this node (m^2)
	double cfa_delx_pred;	// Predicted cell length for this node (m^2)
	double d;				// Pipe diameter at node (m)
	double dddx;			// First derivative of diameter variation at node (m/m)
	double d2ddx2;			// Second derivative of diameter variation at node
	double d_pred;			// Predicted pipe diameter at node (m)
	double dddx_pred;		// Predicted first derivative of diameter variation at node (m/m)
	double d2ddx2_pred;		// Predicted second derivative of diameter variation at node
	double DELX_R, DELX_L;	// Non-dimensional node spacing, to the left and right of the node, respectively
	double dfdx;			// Area gradient at node (m^2/m)
	double dfdx_pred;		// Predicted area gradient at node (m^2/m)
	double f;				// Cross-sectional area at node (m^2)
	double f_dash;			// Non-dimensional cross-sectional area at node
	double f_pred;			// Predicted cross-sectional area at node (m^2)
	double f_prev;			// Cross-sectional area at node at previous time step (m^2)
	double vol;				// Volume (m^3)
	double vol_prev;		// Volume at previous time step (m^3)
	double x;				// Effective distance from odd end of pipe (m)
	double X;				// Effective non-dimensional distance from odd end

	// Mesh Method of Characteristics variables
	// ----------------------------------------------------------------------------------------------------
	double *CL1, *CL2;		// 2 item list for new and old Riemann invariants
	double *AA;				// Entropy level (used in non-homentropic MOC)
	double U, A;			// Non-dimensional velocity, speed of sound
	double M;				// Mach number;
	double p_dash;			// Non-dimensional static pressure;
	double T;				// Temperature (K)
	double mdot;			// Mass flow rate (kg.s^-1)
	
	// W_alpha_beta schemes variables
	// ----------------------------------------------------------------------------------------------------
	double *W, *F, *C;					// Solution vector W, flux vetor F, souce vector C
	double *W_pred, *F_pred, *C_pred;	// Predictor solution vector W_pred, predictor flux vector F_pred, predictor source vector C_pred
	double *S, *S_pred;					// Artificial viscosity and predictor artificial viscosity vector
	double *W_prev;						// Stores solution at previous time step

	// Primitive variables
	// ----------------------------------------------------------------------------------------------------
	double rho;
	double rho_prev;

	// Derived variables
	// ----------------------------------------------------------------------------------------------------
	double Re;				// Reynold's number

	// Screen output
	// ----------------------------------------------------------------------------------------------------
	bool CHOKED;			// Set true by the attached boundary method when choking occurs
	bool STEADY;			// Nodal velocity is steady between time steps

	// Haemodynamics variables
	// ----------------------------------------------------------------------------------------------------
	double B;				// Function of fluid density, unpressurized area and current area
	double B_pred;			// Function of fluid density, unpressurized area and current area evaluated at the predictor step
	double d0;				// Undeformed vessel diameter at node (m)
	double d0_pred;			// Undeformed vessel diameter at the predictor(m)
	double dd0dx;			// First derivative of undeformed diameter variation at node (m/m)
	double dd0dx_pred;		// First derivative of undeformed diameter variation at the predictor (m/m)
	double f0;				// Undeformed vessel cross-sectional area at the node (m^2)
	double f0_pred;			// Undeformed vessel cross-sectional area at the predictor (m^2)
	double p0_dash;			// Undeformed non-dimensional pressure at the node
	double b;				// Approximation of 4/3 * Eh/r0 based on empirical constants
	double b_pred;			// Approximation of 4/3 * Eh/r0 based on empirical constants evaluated at the predictor step
	double dbdr0;			// First derivative of b with respect to the inital radius of the vessel
	double dbdr0_pred;		// First derivative of b with respect to the inital radius the vessel evaluated at the predictor

	// ADDING MEMBER VARIABLES? REMEMBER TO UPDATE CLASS OPERATORS!
};

#endif // !defined(AFX_NODE_H__95E9F4DB_1C9D_4626_AE7E_0D043583F5CD__INCLUDED_)
