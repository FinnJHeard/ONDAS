// Joiner.h: interface for the CJoiner class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_JOINER_H__CE2333F6_292A_441B_8EDE_23347A51D101__INCLUDED_)
#define AFX_JOINER_H__CE2333F6_292A_441B_8EDE_23347A51D101__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "Properties.h"
//#include "FiniteVolume.h"
#include "Node.h"
//#include "PathLine.h"


class CJoiner  
{
public:
	CJoiner();
	virtual ~CJoiner();


public:
	CNode Node;					// The single node

};

#endif // !defined(AFX_JOINER_H__CE2333F6_292A_441B_8EDE_23347A51D101__INCLUDED_)
