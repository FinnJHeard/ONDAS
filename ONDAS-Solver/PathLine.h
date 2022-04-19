// PathLine.h: interface for the CPathLine class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PATHLINE_H__9393FEE7_2DCC_4048_B75C_AE5113B96AC8__INCLUDED_)
#define AFX_PATHLINE_H__9393FEE7_2DCC_4048_B75C_AE5113B96AC8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CPathLine  
{
public:
	CPathLine();
	virtual ~CPathLine();

	CPathLine(const CPathLine& inPathLine); // Copy constructor

	CPathLine& CPathLine::operator=(const CPathLine& inPathLine); // Overload =

	void Print(void);

private:

public:
	int ID;					// ID
	double AA;				// Entropy level
	double XK;				// Distance of the path line from the odd end
	double XK_old;			// Old distance of the path line from the odd end
	double lambda_K;		// Riemann value
	double beta_K;			// Riemann value

	bool new_pathline;
};

#endif // !defined(AFX_PATHLINE_H__9393FEE7_2DCC_4048_B75C_AE5113B96AC8__INCLUDED_)
