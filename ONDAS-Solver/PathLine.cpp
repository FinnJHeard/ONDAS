// PathLine.cpp: implementation of the CPathLine class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "PathLine.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPathLine::CPathLine()
{

}

CPathLine::~CPathLine()
{

}

// Copy constructor
CPathLine::CPathLine(const CPathLine& inPathLine)
{
	ID = inPathLine.ID;
	AA = inPathLine.AA;
	XK = inPathLine.XK;
	XK_old = inPathLine.XK_old;
	lambda_K = inPathLine.lambda_K;
	beta_K = inPathLine.beta_K;
	new_pathline = inPathLine.new_pathline;
}


CPathLine& CPathLine::operator=(const CPathLine& inPathLine)
{
	if(this != &inPathLine)
	{
		ID = inPathLine.ID;
		AA = inPathLine.AA;
		XK = inPathLine.XK;
		XK_old = inPathLine.XK_old;
		lambda_K = inPathLine.lambda_K;
		beta_K = inPathLine.beta_K;

		new_pathline = inPathLine.new_pathline;
	}
	return *this;
}

void CPathLine::Print()
{
	cout << "Path line " << this->ID << ":" << endl;
	cout << "XK_old = " << this->XK_old << endl;
	cout << "XK = " << this->XK << endl;
	cout << "AA = " << this->AA << endl;
	cout << "lambda_K = " << this->lambda_K << endl;
	cout << "beta_K = " << this->beta_K << endl;
}
