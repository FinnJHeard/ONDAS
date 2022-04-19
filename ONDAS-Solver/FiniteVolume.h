// FiniteVolume.h: interface for the CFiniteVolume class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FINITEVOLUME_H__BE19E240_E38E_41FA_9BE4_A6FE0D611AD0__INCLUDED_)
#define AFX_FINITEVOLUME_H__BE19E240_E38E_41FA_9BE4_A6FE0D611AD0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CFiniteVolume  
{
public:
	CFiniteVolume();
	virtual ~CFiniteVolume();



public:
/*
	// Configuration variables (common)
	// ================================
	int bc;			// Record boundary attachement type (INTERIOR for non-boundary nodes)
	int side;		// LEFT or RIGHT if a boundary node, INSIDE otherwise
	int rflow;		// If a junction node, an alias of its pipe odd_end_flow or even_end_flow
	
	// Geometry variables (common)
	// ===========================
	double x;		// Effective distance from odd end of pipe (m)
	double X;		// Effective non-dimensional distance from odd end
	double DELX_R, DELX_L; // Non-dimensional node spacing, to the left and right of the node, respectively
	double d;		// Pipe diameter at node (m)
	double f;		// Cross-sectional area at node (m^2)
	double f_dash;	// Non-dimensional cross-sectional area at node

	// Mesh Method of Characteristics variables
	// ========================================
	double *CL1, *CL2;		// 2 item list for new and old Riemann values
	double *AA;				// Entropy level (used in non-homentropic MOC)
	double U, A;			// Non-dimensional velocity, speed of sound
	double p_dash;			// Non-dimesional static pressure;
	double T;				// Temperature (K)
	double mdot;			// Mass flow rate (kg.s^-1)
*/

	// Filling and emptying variables
	// ==============================
	double l;		// Length of finute volume (m)
	double V;		// Volume of finite volume (m^3)
	double m;		// Mass contained within finite volume (kg)
	
	double T0;		// Stagnation (==static) temperature of finite volume (K)
	double p0;		// Stagnation (==static) pressure of finite volume (Pa)
	double p0_old;
	double h0;
};

#endif // !defined(AFX_FINITEVOLUME_H__BE19E240_E38E_41FA_9BE4_A6FE0D611AD0__INCLUDED_)
