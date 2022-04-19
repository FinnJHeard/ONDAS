// Time.h: interface for the CTime class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TIME_H
#define TIME_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <time.h>
#include <sys/timeb.h>

#include "Properties.h"

class CTime  
{
public:
	CTime();
	virtual ~CTime();
	double m_time;
	double del_t;
	void Increment(CProperties* pPpt, double IncreaseInTime);
	void StartTimer(CProperties* pPpt);
	void StopTimer(CProperties* pPpt);
	void OutStartTime(CProperties* pPpt);
	void OutStopTime(CProperties* pPpt);
	void OutTimeLapsed(CProperties* pPpt);

	void RestartCPUTime(CProperties* pPpt);
	void StopEvaluateCPUTime(CProperties* pPpt);
	void Pause(CProperties* pPpt, double howLong);
	void CheckIfXSecondsElapsed(CProperties* pPpt, double XSeconds);

private:
	struct timeb tpointer;
    time_t tstart, tstop, current_stop_time, last_true_time;
	int msstart, msstop;
	double stotal, mstotal;
	
	time_t cpu_restart_time, cpu_stop_time;
    int cpu_restart_time_ms, cpu_stop_time_ms;
    double cpu_time, cpu_time_ms;

public:
	double elapsed_time;
	bool XSecondsHaveElapsed;
};

#endif 
