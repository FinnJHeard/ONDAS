// EndCap.h: interface for the CEndCap class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ENDCAP_H__62FC5EB2_097B_467B_B8A3_20DA8CC4CD94__INCLUDED_)
#define AFX_ENDCAP_H__62FC5EB2_097B_467B_B8A3_20DA8CC4CD94__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Boundary.h"
//#include "PathLine.h"
//#include "Pipe.h"
//#include "Properties.h"

class CEndCap : public CBoundary  
{
public:
	CEndCap();
	virtual ~CEndCap();
	CEndCap(const CEndCap& inEndCap);
	CEndCap& operator=(const CEndCap& inEndCap);
	void Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rEndCapPIPES, int** &rEndCapPIPES_ENDS, double* &rENDCORR, int id, bool ex, int npipes, int assyid, string parent_assy_res_dir, string calling_object_str);
	void RunBoundary();
	void PrintToFile(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca);
private:
};

#endif // !defined(AFX_ENDCAP_H__62FC5EB2_097B_467B_B8A3_20DA8CC4CD94__INCLUDED_)
