// Time.cpp: implementation of the CTime class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Time.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTime::CTime()
{
	del_t = 100;

	msstart = 0;
	msstop = 0;
	stotal = 0;
	mstotal = 0;
	tstart = 0;
	tstop = 0;

	m_time = 0;
	cpu_time = 0;
	cpu_time_ms = 0;
	cpu_restart_time = 0;
	cpu_restart_time_ms = 0;
	cpu_stop_time = 0;
	cpu_stop_time_ms = 0;

	current_stop_time = 0;
	last_true_time = 0;
	elapsed_time = 0;
	XSecondsHaveElapsed = false;
}

CTime::~CTime()
{
}

void CTime::Increment(CProperties* pPpt, double IncreaseInTime)
{
	if(pPpt->SHOW_calls){pPpt->Out("MyTime.Increment\n");}
	m_time += IncreaseInTime;
}

void CTime::StartTimer(CProperties* pPpt)
// Only called at the very start of a simulation
{
	if(pPpt->SHOW_calls){pPpt->Out("MyTime.StartTimer\n");}
	ftime(&tpointer);
	tstart = tpointer.time;
	msstart = tpointer.millitm;
}

void CTime::StopTimer(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out("MyTime.StopTimer\n");}
	ftime(&tpointer);
	tstop = tpointer.time;
	msstop = tpointer.millitm;
}

void CTime::OutStartTime(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out("MyTime.OutStartTime\n");}
	Line(pPpt, 80, "%");
	pPpt->Out(" Start time\t\t\t\t\t\t"); pPpt->Out(ctime(&tstart));
	Line(pPpt, 80, "%");
}

void CTime::OutStopTime(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out("MyTime.OutStopTime\n");}
	pPpt->Out("\n");
	Line(pPpt, 80, "%");
	pPpt->Out("Job complete\n============\n");
	pPpt->Out(" End time\t\t\t"); pPpt->Out(ctime(&tstop));
	OutTimeLapsed(pPpt);
	Line(pPpt, 80, "%");
}

void CTime::OutTimeLapsed(CProperties* pPpt)
{
	stotal = difftime(tstop,tstart);
	mstotal = stotal*1000 + (msstop - msstart);
	pPpt->Out(" CPU time accumulated\t\t"); pPpt->Out(cpu_time); pPpt->Out(" s ("); pPpt->Out(cpu_time_ms); pPpt->Out(" ms)"); pPpt->Out("\n");
	pPpt->Out(" Real time elapsed\t\t"); pPpt->Out(stotal); pPpt->Out(" s ("); pPpt->Out(mstotal); pPpt->Out(" ms)"); pPpt->Out("\n");
}

void CTime::RestartCPUTime(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out("MyTime.RestartCPUTime\n");}
	ftime(&tpointer);
	cpu_restart_time = tpointer.time;
	cpu_restart_time_ms = tpointer.millitm;
}

void CTime::StopEvaluateCPUTime(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out("MyTime.StopEvaluateCPUTime\n");}
	ftime(&tpointer);
    cpu_stop_time = tpointer.time;
	cpu_stop_time_ms = tpointer.millitm;
	cpu_time += difftime(cpu_stop_time, cpu_restart_time);
	cpu_time_ms = cpu_time*1000 + cpu_stop_time_ms - cpu_restart_time_ms;
}

void CTime::Pause(CProperties* pPpt, double howLong)
{
	if(pPpt->SHOW_calls){pPpt->Out("MyTime.Pause\n");}
	RestartCPUTime(pPpt);
    StopEvaluateCPUTime(pPpt);
	ftime(&tpointer);
	time_t t_check;
    int t_check_ms;
    do
	{
        for(int i=0; i<1e7; ++i); // Creates a very small delay (~30ms)
		ftime(&tpointer); // Continue to check current time
        t_check = tpointer.time;
        t_check_ms = tpointer.millitm;
	}
    while(difftime(t_check,cpu_stop_time)*1000 + (t_check_ms - cpu_stop_time_ms)  < howLong*1000);
	RestartCPUTime(pPpt);
}

void CTime::CheckIfXSecondsElapsed(CProperties* pPpt, double XSeconds)
// Used to control timing of simulation output to screen
{
	if(pPpt->SHOW_calls){pPpt->Out("MyTime.CheckIfXSecondsElapsed\n");}
	ftime(&tpointer);
	current_stop_time = tpointer.time;
	elapsed_time = difftime(current_stop_time,tstart);
	if(difftime(current_stop_time,last_true_time)>=XSeconds)
	{
		XSecondsHaveElapsed = true;
		last_true_time = current_stop_time;
	}
	else XSecondsHaveElapsed = false;
}