// ====================================================================================================
// Main.cpp : Defines the entry point for the console application.
// ====================================================================================================
#include "stdafx.h"
#include "Globals.h"
#include "Tools.h"
#include "Properties.h"
#include "Assembly.h"
#include <windows.h>
#include <iostream>
#include <exception> // Standard library exceptions
using namespace std;

// Declare functions
// ----------------------------------------------------------------------------------------------------
void PrintToScreen(double test_var);
void Backup();
void Restore();
void CheckConfiguration(CProperties* pPpt, vector<string> &rConfStrs);
void RecordTappings(double test_var, double time_step_length_exh, double time_step_length_int);
void PrintToFiles();
void CloseFiles();
void FillingAndEmptying(double DELZe, double DELZi);
void ResetInitialConditions();

// Declare time object, property set and assembly array
// ----------------------------------------------------------------------------------------------------
CTime MyTime;
CProperties Ppt;
CAssembly *Assembly;

// Declare timing variables
// ----------------------------------------------------------------------------------------------------
int timestep;			// Timestep number
double *Ze, *Zi;		// Dimensionless time
double *TIMEe;			// Actual exhaust time (s)
double *TIMEi;			// Actual intake time (s)
double time_step_length_exh;
double time_step_length_int;
double factor;			// Timestep adjustment factor
bool RESTORE;			// Flag
int restores;			// Number of times simulation has been restored to previous values
bool ITERATION_COMPLETE;// Flag.

// Simulation case .out file
// ----------------------------------------------------------------------------------------------------
//ofstream outfile;
	
// Start of main function
// ----------------------------------------------------------------------------------------------------
int main(){
	try{			
		int p, i;
		// Start timer		 
		// ----------------------------------------------------------------------------------------------------
		MyTime.StartTimer(&Ppt);

		// Set code version (ymmdd) - must match FILE_VERSION in MAIN.txt		 
		// ----------------------------------------------------------------------------------------------------
		int CODE_VERSION = 210428;
		string dat_dir = "..\\dat." + IntToString(CODE_VERSION) + "\\";

		// Set location of data files and read case name from CASE_NAME.txt
		// ----------------------------------------------------------------------------------------------------
		Ppt.assemblies_dir = dat_dir + "assy\\";		// Set path to assemblies directory
		Ppt.cases_dir = dat_dir + "case\\";				// Set path to cases directory
		Ppt.param_dir = dat_dir + "param\\";			// Set path to parameters directory
		int max_num = 200;								// Maximum number of parameters permitted
		string case_name_dir_file = "CASE_NAME.txt";	// States directory case name that contains directory files
		Ppt.ReadCaseName(ConstructString(&Ppt, Ppt.cases_dir, case_name_dir_file), max_num); // Sets case_name
		Ppt.case_dir = Ppt.cases_dir + Ppt.case_name_slash;	// Set path to directory containing results directories for each assembly	

		// Open new simulation case .out file without overwriting existing ones
		// ----------------------------------------------------------------------------------------------------
		LPSECURITY_ATTRIBUTES attr; attr = NULL;
		int fk=1; bool EXISTS; FILE *fp;
		string str_outfile, str_resdir;
		do{
			EXISTS = false;
			str_outfile = Ppt.case_dir + Ppt.case_name;
			if(fk<10) str_outfile += ".00"; else{if(fk<100) str_outfile += ".0"; else if(fk<1000) str_outfile += ".";}
			str_resdir = str_outfile;
			str_outfile += IntToString(fk) + ".out"; str_resdir += IntToString(fk) + ".res";
			fp = fopen(StringToChar(str_outfile),"r");
			if(fp==NULL);
			else{EXISTS = true; fclose(fp); ++fk;} 
			if(fk>999){cout << "Maximum number of .out files allowed has been reached. Delete some and try again. Exiting.\n"; exit(1);} 
		}while(EXISTS);
		Ppt.OUTFILE = fopen(StringToChar(str_outfile), "w");
		CreateDirectory(StringToChar(str_resdir), attr);	// Create a numbered results directory within case directory
		Ppt.case_res_dir = StringToChar(str_resdir + "\\"); // Set path to numbered results directory within case directory

		// Read simulation properties from the main_input file in the case directory
		// ----------------------------------------------------------------------------------------------------
		string main_input = "MAIN_INPUT.txt";			// Set main input file name
		Ppt.ReadInput(ConstructString(&Ppt, Ppt.case_dir, main_input), max_num);
		
		// Print start of program header to screen
		// ----------------------------------------------------------------------------------------------------
		Message(&Ppt, CODE_VERSION); 
		
		// Output simulation start time
		// ----------------------------------------------------------------------------------------------------
		Ppt.Out("\n");
		MyTime.OutStartTime(&Ppt); 
		Ppt.Out("\n");
		MyTime.Pause(&Ppt, 0.5);
		
		// List case name and simulation properties
		// ----------------------------------------------------------------------------------------------------
		if(Ppt.SHOW_globals) Ppt.ListProperties(CODE_VERSION); 
		
		// Dimension assemblies
		// ----------------------------------------------------------------------------------------------------
		Assembly = new CAssembly[Ppt.NASSEMBLY];

		// Initialise all assemblies 
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Ppt.NASSEMBLY; ++i) Assembly[i].Initialise(&Ppt, i);

		// Read configuration files for required boundary types
		// ----------------------------------------------------------------------------------------------------
		vector<string> ConfStrs;
		for(i=0; i<Ppt.NASSEMBLY; ++i) Assembly[i].ConfigureAll(&Ppt, ConfStrs);
		
		// Check no duplicate boundary conditions exist
		// ----------------------------------------------------------------------------------------------------
		CheckConfiguration(&Ppt, ConfStrs);
		
		// Initialise boundary conditions (and their backups) followed by initialisation of the pipes
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Ppt.NASSEMBLY; ++i) Assembly[i].InitialiseBCs(&MyTime, &Ppt, Assembly);
		
		// Configure boundary nodes (can now attach boundaries to pipe boundary nodes)
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Ppt.NASSEMBLY; ++i) Assembly[i].ConfigureBCs(&Ppt);

		// Print individual object properties to screen
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Ppt.NASSEMBLY; ++i) Assembly[i].ListPropertiesAll(&Ppt);
				
		// Check and print boundary connections
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Ppt.NASSEMBLY; ++i) Assembly[i].PrintBoundaryConnections(&Ppt);
		
		// If simulation type is set to periodic, but there are no cylinders, reset it to continuous
		// ----------------------------------------------------------------------------------------------------
		Ppt.CONTINUOUS = true; i=0;
		while(Ppt.CONTINUOUS && i<Ppt.NASSEMBLY){if(Assembly[i].NCYLS!=0) Ppt.CONTINUOUS = false; ++i;}

		// Dimension timing variables
		// ----------------------------------------------------------------------------------------------------
		Ze = new double [2];
		Zi = new double [2];
		TIMEe = new double [2];
		TIMEi = new double [2];

		// Zero timing variables and print start time to screen
		// ----------------------------------------------------------------------------------------------------
		Ze[R] = 0;
		Zi[R] = 0;
		Ze[R+1] = 0;
		Zi[R+1] = 0;
		TIMEe[R] = 0;
		TIMEi[R] = 0;
		TIMEe[R+1] = 0;
		TIMEi[R+1] = 0;
		factor = 1;
		time_step_length_exh = 0;
		time_step_length_int = 0;

		RESTORE = false;
		restores = 0;
		double DELZtemp;
		double DELZe=1000, DELZe1;
		double DELZi=1000, DELZi1;
/*
		// Update domain initial conditions
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Ppt.NASSEMBLY; ++i)
		{		
			// Pipes
			// ----------------------------------------------------------------------------------------------------
			for(p=0; p<Assembly[i].NEXPIPES; ++p) Assembly[i].Exhaust[p].Update(&Ppt, DELZe); 
			for(p=0; p<Assembly[i].NINPIPES; ++p) Assembly[i].Intake[p].Update(&Ppt, DELZi);

			// Anechoic buffers and dampers
			// ----------------------------------------------------------------------------------------------------
			for(p=0; p<Assembly[i].NEXANEC; ++p)
			{
				Assembly[i].EX_ANEC[p].Buffer.Update(&Ppt, DELZe);
				Assembly[i].EX_ANEC[p].Damper.Update(&Ppt, DELZe);
			}
			for(p=0; p<Assembly[i].NINANEC; ++p)
			{
				Assembly[i].IN_ANEC[p].Buffer.Update(&Ppt, DELZi);
				Assembly[i].IN_ANEC[p].Damper.Update(&Ppt, DELZi);
			}
		}
//*/
		// Record initial domain data
		// ----------------------------------------------------------------------------------------------------
		double test_var;
		if(Ppt.CONTINUOUS) test_var = TIMEe[R+1]; else test_var = Assembly[0].Eng[0].ca_elapsed;

		RecordTappings(test_var, time_step_length_exh, time_step_length_int);	// Record initial tappings
		PrintToFiles();				// Print to all files
		
		// Ask user to begin simulation
		// ----------------------------------------------------------------------------------------------------
		PrintToScreen(test_var);
		Ppt.Out("Continue? (y/n):");
		char cont;
		bool stop = false;
		double until = 0;
		cin >> cont;
		Ppt.Out("\n");
		if(cont=='n' || cont=='N')
		{
			Ppt.Out("User answers NO\n");
			stop=true;
		}
		else 
		{
			Ppt.Out("User answers YES\n");
			if(Ppt.CONTINUOUS) Ppt.Out("Go to time (s): ");
			else
			{
				Ppt.Out("Go to crank angle ("); Ppt.Out(Deg()); Ppt.Out("CA): ");
			}
			cin >> until;
			MyTime.RestartCPUTime(&Ppt);
		}

		int a;
		bool CONTINUING = false; // Flags whether the simulation has been paused on the previous step
		
		// Start of main loop
		// ----------------------------------------------------------------------------------------------------
		timestep = 0;
		Ppt.STOP = false;
		while(TIMEe[R+1]<Ppt.ZMAX && !stop) {
			MyTime.CheckIfXSecondsElapsed(&Ppt, Ppt.outputInterval);
			
			if (MyTime.XSecondsHaveElapsed) { // Another interval has elapsed so print something to screen

				if (Ppt.CONTINUOUS) SimulationOutput(&Ppt, timestep, TIMEe[R + 1], MyTime.elapsed_time, CONTINUING);
				else SimulationOutput(&Ppt, timestep, TIMEe[R + 1], MyTime.elapsed_time, CONTINUING, Assembly[0].Eng[0].ca_elapsed);
			}
			
			timestep += 1;
			ITERATION_COMPLETE = false;
			RESTORE = false;
			bool done = false;		
			do {
				// Calculate timestep
				// ----------------------------------------------------------------------------------------------------
				DELZe = 1000;
				DELZi = 1000;
	
				if(timestep==1) { // Use minimum first time step
					DELZi = (Ppt.delt_min/Ppt.xref)*Ppt.AREFi;
					DELZe = (Ppt.delt_min/Ppt.xref)*Ppt.AREFe;
				}
				else {
					if(!Ppt.CONSTANT_DELT) { // If not constant time step
						for(a=0; a<Ppt.NASSEMBLY; ++a) {
							// Pipes
							for(p=0; p<Assembly[a].NEXPIPES; ++p) {
								DELZe1 = Assembly[a].Exhaust[p].TimeStepMOC(&Ppt);
								if(DELZe1<DELZe) DELZe = DELZe1;
							}
							for(p=0; p<Assembly[a].NINPIPES; ++p) {
								DELZi1 = Assembly[a].Intake[p].TimeStepMOC(&Ppt);
								if(DELZi1<DELZi) DELZi = DELZi1;
							}

							// Anechoic buffers and dampers
							for(p=0; p<Assembly[a].NEXANEC; ++p) {
								DELZe1 = Assembly[a].EX_ANEC[p].Buffer.TimeStepMOC(&Ppt);
								if(DELZe1<DELZe) DELZe = DELZe1;
								DELZe1 = Assembly[a].EX_ANEC[p].Damper.TimeStepMOC(&Ppt);
								if(DELZe1<DELZe) DELZe = DELZe1;
							}
							for(p=0; p<Assembly[a].NINANEC; ++p) {
								DELZi1 = Assembly[a].IN_ANEC[p].Buffer.TimeStepMOC(&Ppt);
								if(DELZi1<DELZi) DELZi = DELZi1;
								DELZi1 = Assembly[a].IN_ANEC[p].Damper.TimeStepMOC(&Ppt);
								if(DELZi1<DELZi) DELZi = DELZi1;
							}
						}
		
						DELZtemp = (Ppt.AREFi/Ppt.AREFe)*DELZe;
						if(DELZi>DELZtemp) DELZi = DELZtemp;
						else DELZe = (Ppt.AREFe/Ppt.AREFi)*DELZi;
		
						if((DELZe/Ppt.AREFe)*Ppt.xref > Ppt.delt_max) {
							DELZe = (Ppt.delt_max/Ppt.xref)*Ppt.AREFe;
							DELZi = (Ppt.delt_max/Ppt.xref)*Ppt.AREFi;
						}
						else {
							if((DELZe/Ppt.AREFe)*Ppt.xref < Ppt.delt_min) {
								DELZe = (Ppt.delt_min/Ppt.xref)*Ppt.AREFe;
								DELZi = (Ppt.delt_min/Ppt.xref)*Ppt.AREFi;
							}
						}

						if((DELZi/Ppt.AREFi)*Ppt.xref > Ppt.delt_max) {
							DELZi = (Ppt.delt_max/Ppt.xref)*Ppt.AREFi;
							DELZe = (Ppt.delt_max/Ppt.xref)*Ppt.AREFe;
						}
						else {
							if((DELZi/Ppt.AREFi)*Ppt.xref < Ppt.delt_min) {
								DELZi = (Ppt.delt_min/Ppt.xref)*Ppt.AREFi;
								DELZe = (Ppt.delt_min/Ppt.xref)*Ppt.AREFe;
							}
						}
					}
					else {
						DELZi = (Ppt.delt_max/Ppt.xref)*Ppt.AREFi;
						DELZe = (Ppt.delt_max/Ppt.xref)*Ppt.AREFe;
					}
				}
DELZe *= factor;
DELZi *= factor;
	
				// Calculate time elapsed at next step
				Ze[R+1] = Ze[R] + DELZe;
				Zi[R+1] = Zi[R] + DELZi;
				TIMEe[R+1] = (Ze[R+1]/Ppt.AREFe)*Ppt.xref;
				TIMEi[R+1] = (Zi[R+1]/Ppt.AREFi)*Ppt.xref;

				// Copy/backup objects in case of rollback
				// ----------------------------------------------------------------------------------------------------
				Backup();

				// Propagate charcateristics through each pipe in every assembly
				// ----------------------------------------------------------------------------------------------------
				for(a=0; a<Ppt.NASSEMBLY; ++a)
				{					
					// Pipes
					// ----------------------------------------------------------------------------------------------------
					for (p = 0; p < Assembly[a].NEXPIPES; ++p) {
						// Only propagate if any requested pipe delay has elapsed
						if(TIMEe[R] >= Assembly[a].Exhaust[p].delay) Assembly[a].Exhaust[p].RunPropagation(&Ppt, DELZe, timestep, RESTORE);
					}
					for (p = 0; p < Assembly[a].NINPIPES; ++p) {
						// Only propagate if any requested pipe delay has elapsed
						if(TIMEi[R] >= Assembly[a].Intake[p].delay) Assembly[a].Intake[p].RunPropagation(&Ppt, DELZi, timestep, RESTORE);
					}
	
					// Anechoic buffers and dampers
					// ----------------------------------------------------------------------------------------------------
					for(p=0; p<Assembly[a].NEXANEC; ++p)
					{
						Assembly[a].EX_ANEC[p].Buffer.RunPropagation(&Ppt, DELZe, timestep, RESTORE);
						Assembly[a].EX_ANEC[p].Damper.RunPropagation(&Ppt, DELZe, timestep, RESTORE);
					}
					for(p=0; p<Assembly[a].NINANEC; ++p)
					{
						Assembly[a].IN_ANEC[p].Buffer.RunPropagation(&Ppt, DELZi, timestep, RESTORE);
						Assembly[a].IN_ANEC[p].Damper.RunPropagation(&Ppt, DELZi, timestep, RESTORE);
					}
				}

				if(Ppt.ALTERNATE_MAC && (Ppt.beta==0 || Ppt.beta==1)) Ppt.beta = !Ppt.beta; // Alternate forward pred./backward corr., backward pred./forward corr. - put here so only done once for all pipes

				// Boundary calculations for all assemblies
				// ----------------------------------------------------------------------------------------------------
				for(a=0; a<Ppt.NASSEMBLY; ++a) Assembly[a].RunBoundary(&Ppt, MyTime, DELZe, DELZi, TIMEe, TIMEi, timestep, RESTORE);

				// Results for each assembly
				// ----------------------------------------------------------------------------------------------------
				// Calculate secondary values (pressures etc.) at all nodes
				for(a=0; a<Ppt.NASSEMBLY; ++a)
				{				
					// Pipes
					// ----------------------------------------------------------------------------------------------------
					for(p=0; p<Assembly[a].NEXPIPES; ++p) Assembly[a].Exhaust[p].Update(&Ppt, DELZe);
					for(p=0; p<Assembly[a].NINPIPES; ++p) Assembly[a].Intake[p].Update(&Ppt, DELZi);
	
					// Anechoic buffers and dampers
					// ----------------------------------------------------------------------------------------------------
					for(p=0; p<Assembly[a].NEXANEC; ++p) {
						Assembly[a].EX_ANEC[p].Buffer.Update(&Ppt, DELZe);
						Assembly[a].EX_ANEC[p].Damper.Update(&Ppt, DELZe);
					}
					for(p=0; p<Assembly[a].NINANEC; ++p) {
						Assembly[a].IN_ANEC[p].Buffer.Update(&Ppt, DELZi);
						Assembly[a].IN_ANEC[p].Damper.Update(&Ppt, DELZi);
					}
				}	
///*
				// Restore if required
				// ----------------------------------------------------------------------------------------------------
				if(RESTORE)
				{
					++restores;
					Restore();	// Restore object backup values
				
					// Reduce timestep for next time via factor (reset to unity if no further restore required
				//	factor *= 0.5;
					done = true;
				}
				else
				{
					RESTORE = false;
					factor = 1;
					done = true;
				}
//*/
			}while(!done);

			// Update timing variables on successful completion of the iteration
			// ----------------------------------------------------------------------------------------------------
			time_step_length_exh = TIMEe[R+1] - TIMEe[R];
			time_step_length_int = TIMEi[R + 1] - TIMEi[R];
			Ze[R] = Ze[R+1];
			Zi[R] = Zi[R+1];
			TIMEe[R] = TIMEe[R+1];
			TIMEi[R] = TIMEi[R+1];
			MyTime.Increment(&Ppt, TIMEe[R+1]);	// Add length of timestep to accumulative time // BUT TIMEe[R+1] is the current time, not the time step length

//cout << "time_step_length_exh = " << time_step_length_exh << endl;
//cout << "TIMEe[R+1] = " << TIMEe[R + 1] << endl;
 
			ITERATION_COMPLETE = true;
			
			// Test if simulation has reached desired point; print data to screen for all objects present
			// ----------------------------------------------------------------------------------------------------
			if(Ppt.CONTINUOUS) test_var = TIMEe[R+1]; else test_var = Assembly[0].Eng[0].ca_elapsed;
			if(test_var >= until 
			|| Ppt.STOP
		//	|| (RESTORE
			)
			{
				Ppt.STOP = false;
				if(Ppt.BEEP) printf("%c",toascii(7));
				MyTime.StopEvaluateCPUTime(&Ppt);
				PrintToScreen(test_var);	
				Ppt.Out("Continue? (y/n):");
				cin >> cont;
				Ppt.Out("\n");
				if(cont=='n' || cont=='N')
				{
					Ppt.Out("User answers NO\n");
					stop=true;
				}
				else 
				{
					Ppt.Out("User answers YES\n");
					if(Ppt.CONTINUOUS) Ppt.Out("Go to time (s): ");
					else 
					{
						Ppt.Out("Go to crank angle ("); Ppt.Out(Deg()); Ppt.Out("CA): ");
					}
					cin >> until;
					MyTime.RestartCPUTime(&Ppt);
					CONTINUING = true;
				}
			}

			// Print to all files
			// ----------------------------------------------------------------------------------------------------
			RecordTappings(test_var, time_step_length_exh, time_step_length_int);
			PrintToFiles();

			// Check if APLDev requires to reset initial gas flow conditions
			// ----------------------------------------------------------------------------------------------------
			bool RESET_INITIAL_CONDITIONS_THIS_TIME = false;
			for(a=0; a<Ppt.NASSEMBLY; ++a)
			{
				int ap;
				for(ap=0; ap<Assembly[a].NEXAPLDev; ++ap)
				{
					if(Assembly[a].APLDevE[ap].Get_RESET_INITIAL_CONDITIONS())
					{
						RESET_INITIAL_CONDITIONS_THIS_TIME = true;
						Assembly[a].APLDevE[ap].Set_RESET_INITIAL_CONDITIONS(false); // Reset the APLDev flag
					}
				}
				for(ap=0; ap<Assembly[a].NINAPLDev; ++ap)
				{
					if(Assembly[a].APLDevI[ap].Get_RESET_INITIAL_CONDITIONS())
					{
						RESET_INITIAL_CONDITIONS_THIS_TIME = true;
						Assembly[a].APLDevI[ap].Set_RESET_INITIAL_CONDITIONS(false); // Reset the APLDev flag
					}
				}
			}
			if(RESET_INITIAL_CONDITIONS_THIS_TIME)
			{
				ResetInitialConditions(); // Re-apply the initial gas conditions everywhere
				//PrintToScreen(test_var); // Print simulation data screen to check that reset occurred
			}
		}
	
		// End of main loop
		// ----------------------------------------------------------------------------------------------------	
		
		// Stop timing, sound bell
		// ----------------------------------------------------------------------------------------------------
		MyTime.StopTimer(&Ppt);
		MyTime.OutStopTime(&Ppt);
	
		// Close files
		// ----------------------------------------------------------------------------------------------------
		CloseFiles();
	}
	catch(exception& e)
	{
		Ppt.Out("standard exception: ");
		Ppt.Out(e.what());
		Ppt.Out("\n");
		exit(1);
	}
	catch(int e)
	{
		cout << "an exception occurred. exception no. " << e << endl;
		exit(1);
	}
	catch(char* e)
	{
		cout << "an exception occurred. exception description: " << e << endl; 
		exit(1);
	}
	catch(...)
	{
		cout << "an unknown exception occurred." << endl;
		exit(1);
	}
	return 0;
}

void PrintToScreen(double test_var)
// ============================================================ //
// Prints current values of object variables to screen			//
// ============================================================ //
{
	if(Ppt.SHOW_calls) Ppt.Out("PrintToScreen\n");
	Ppt.Out("\n");;
	std::string temp;
	temp = ITERATION_COMPLETE ? "Simulation paused; object print out at end of iteration no." 
							  : "Simulation paused; object print out BEFORE end of iteration no.";
	temp += IntToString(timestep);
	temp += ":";
	Ppt.Out(Underline(StringToChar(temp), "=", true));
	Ppt.Out("\n");

	int p, i, a;

	for(a=0; a<Ppt.NASSEMBLY; ++a)
	{
		// Pipes
		// ----------------------------------------------------------------------------------------------------
		for(p=0; p<Assembly[a].NEXPIPES; ++p)
		{
			if(Assembly[a].Exhaust[p].SHOW_DATA) 
			{
				Assembly[a].Exhaust[p].PrintToScreen(&Ppt); Ppt.Out("\n");
			}
		}
		for(p=0; p<Assembly[a].NINPIPES; ++p)
		{
			if(Assembly[a].Intake[p].SHOW_DATA) 
			{
				Assembly[a].Intake[p].PrintToScreen(&Ppt); Ppt.Out("\n");
			}
		}

		for(p=0; p<Assembly[a].NEXVOLUMES; ++p) Assembly[a].ExVolume[p].PrintToScreen(&Ppt);
		for(p=0; p<Assembly[a].NINVOLUMES; ++p) Assembly[a].InVolume[p].PrintToScreen(&Ppt);
		
		for(i=0; i<Assembly[a].NENGINES; ++i){Assembly[a].Eng[i].PrintToScreen(); Ppt.Out("\n");}
		for(i=0; i<Assembly[a].NCYLS; ++i){Assembly[a].Cyl[i].PrintToScreen(&Ppt); Ppt.Out("\n");}
		for(i=0; i<Assembly[a].NEXJUNCS; ++i){Assembly[a].JE[i].PrintToScreen(&Ppt); Ppt.Out("\n");}
		
		// Anechoic ends
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Assembly[a].NEXANEC; ++i) Assembly[a].EX_ANEC[i].PrintToScreen(&Ppt);
		for(i=0; i<Assembly[a].NINANEC; ++i) Assembly[a].IN_ANEC[i].PrintToScreen(&Ppt);

		// End environments
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Assembly[a].NEXEND; ++i) Assembly[a].EX_END[i].PrintToScreen(&Ppt);
		for(i=0; i<Assembly[a].NINEND; ++i) Assembly[a].IN_END[i].PrintToScreen(&Ppt);

		// Adiabatic pressure loss devices
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Assembly[a].NEXAPLDev; ++i) Assembly[a].APLDevE[i].PrintToScreen(&Ppt);
		for(i=0; i<Assembly[a].NINAPLDev; ++i) Assembly[a].APLDevI[i].PrintToScreen(&Ppt);
		
		for(i=0; i<Assembly[a].NTURBINE; ++i){Assembly[a].Turbine[i].PrintToScreen(&Ppt); Ppt.Out("\n");}
		
		// Transmissive boundaries
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<Assembly[a].NEXTRANSM; ++i)
		{
			if(Assembly[a].TransmE[i].PRINT_TO_SCREEN)
				Assembly[a].TransmE[i].PrintToScreen(&Ppt);
		}
		for(i=0; i<Assembly[a].NINTRANSM; ++i)
		{
			if(Assembly[a].TransmI[i].PRINT_TO_SCREEN)
				Assembly[a].TransmI[i].PrintToScreen(&Ppt);
		}
	}
	
	Ppt.Out(Underline("Simulation status" ,"=", true));
	if(ITERATION_COMPLETE)
	{
		if(RESTORE) 
		{
			Ppt.Out(" Iteration no."); Ppt.Out(timestep); Ppt.Out(" RESTORED\n");
		}
		else
		{
			Ppt.Out(" Iteration no."); Ppt.Out(timestep); Ppt.Out(" COMPLETE\n");
		}
//		cout << " Timestep adjustment factor for next iteration, factor = " << factor << " \n";
//		cout << " No. of times simulation has been restored, restores = " << restores << " \n";
		Ppt.Out(" Total time\t\t\t=\t"); Ppt.Out(TIMEe[R+1]); Ppt.Out(" s");
		if(!Ppt.CONTINUOUS)
		{ 
			Ppt.Out("\n"); Ppt.Out(" Elapsed crank angle\t\t=\t"); Ppt.Out(test_var); Ppt.Out(Deg());
		}
		Ppt.Out("\n");
	}
	else
	{
		Ppt.Out(" INSIDE iteration no."); Ppt.Out(timestep); Ppt.Out(" \n");
//		cout << " Timestep adjustment factor, factor = " << factor;
		Ppt.Out("\n");
//		cout << " No. of times simulation has been restored, restores = " << restores << " \n";
//		cout << " Total time = " << TIMEe[R+1] << " s" << endl;
	}
	
	temp = ITERATION_COMPLETE ? "Simulation paused; object print out at end of iteration no." 
							  : "Simulation paused; object print out BEFORE end of iteration no.";
	Ppt.Out(Line(StringToChar(temp),"="));
	Ppt.Out("\n");
}

void Backup()
// ============================================================ //
// Makes copies of objects prior to propagation and				//
// boundary condition updates									//
// ============================================================ //
{
	if(Ppt.SHOW_calls) cout << "Backup" << endl;

/*
	// Backup pipe values before propagation for possible rollback
	for(int p=0; p<Ppt.NEXPIPES; ++p) Exhaust[p].Backup();
	for(p=0; p<Ppt.NINPIPES; ++p) Intake[p].Backup();

	// Backup engines before propagation for possible rollback
	for(int i=0; i<Ppt.NENGINES; ++i) Eng_backup[i] = Eng[i];

	// Backup cylinders before propagation for possible rollback
	for(i=0; i<Ppt.NCYLS; ++i) Cyl_backup[i] = Cyl[i];

	// Backup junctions before propagation for possible rollback
	for(i=0; i<Ppt.NEXJUNCS; ++i) JE_backup[i] = JE[i];
	for(i=0; i<Ppt.NINJUNCS; ++i) JI_backup[i] = JI[i];
*/

	int a, p;
	for(a=0; a<Ppt.NASSEMBLY; ++a)
	{
		///*
		// Backup pipe values before propagation for possible rollback
		for(p=0; p<Assembly[a].NEXPIPES; ++p) Assembly[a].Exhaust[p].Backup(&Ppt);
		for(p=0; p<Assembly[a].NINPIPES; ++p) Assembly[a].Intake[p].Backup(&Ppt);
		//*/
		///*
		// Backup engines before propagation for possible rollback
//		for(i=0; i<Assembly[a].NENGINES; ++i) Assembly[a].Eng_backup[i] = Assembly[a].Eng[i];
		//*/
		///*
		// Backup cylinders before propagation for possible rollback
//		for(i=0; i<Assembly[a].NCYLS; ++i) Assembly[a].Cyl_backup[i] = Assembly[a].Cyl[i];
		//*/
		//*/
		// Backup junctions before propagation for possible rollback
//		for(i=0; i<Assembly[a].NEXJUNCS; ++i) Assembly[a].JE_backup[i] = Assembly[a].JE[i];
//		for(i=0; i<Assembly[a].NINJUNCS; ++i) Assembly[a].JI_backup[i] = Assembly[a].JI[i];
		//*/
	}
}

void Restore()
// ============================================================ //
// Restores the copied values to the working objects, i.e.		//
// moves the simulation back one step							//
// ============================================================ //
{
	if(Ppt.SHOW_calls) cout << "Restore" << endl;

/*
	// Restore pipes
	for(int p=0; p<Ppt.NEXPIPES; ++p) Exhaust[p].Restore();
	for(p=0; p<Ppt.NINPIPES; ++p) Intake[p].Restore();

///*
	// Restore engines
	for(int i=0; i<Ppt.NENGINES; ++i)
	{
		Eng[i] = Eng_backup[i];
//		Eng[i].rolledback = true;
	}

	// Restore cylinders
	for(i=0; i<Ppt.NCYLS; ++i)
	{
		Cyl[i] = Cyl_backup[i];
//		Cyl[i].rolledback = true;
	}

	// Restore junctions
	for(i=0; i<Ppt.NEXJUNCS; ++i)
	{
		JE[i] = JE_backup[i];
		JE[i].rolledback = true;
	}
	for(i=0; i<Ppt.NINJUNCS; ++i)
	{ 
		JI[i] = JI_backup[i];
		JI[i].rolledback = true;
	}


	// Restore timestep and times
	timestep-=1;
	Ze[R+1] = Ze[R];
	Zi[R+1] = Zi[R];
	TIMEe[R+1] = (Ze[R]/Ppt.AREFe)*Ppt.xref;
	TIMEi[R+1] = (Zi[R]/Ppt.AREFi)*Ppt.xref;0
//	TIMEe[R+1] = (Ze[R]/Ppt.AREF)*Ppt.xref;
//	TIMEi[R+1] = (Zi[R]/Ppt.AREF)*Ppt.xref;
*/
}

void CheckConfiguration(CProperties* pPpt, vector<string> &rConfStrs)
{
	if(Ppt.SHOW_calls){pPpt->Out("Main.CheckConfiguration\n");}
	//using std::cout;
	for(int v=0; v<int(rConfStrs.size()); ++v)
	{
		int w = v+1; // Need only compare after index v
		while(w<int(rConfStrs.size()))
		{
			/*
			std::cout << "Comparing ConfStrs[v=" << v << "] = ";
			std::cout << rConfStrs[v];
			std::cout << " against ConfStrs[w=" << w << "] = ";
			std::cout << rConfStrs[w];
			std::cout << std::endl;
			//*/
			if(rConfStrs[v]==rConfStrs[w])
			{
				std::cout << "Boundary condition configuration duplicated at: ";
				std::cout << rConfStrs[v];
				std::cout << std::endl;
				std::cout << "Remove duplicate boundary and try again. Exiting..." << std::endl << std::endl;
				exit(1);
			}
			++w;
		}
	}
}

void RecordTappings(double test_var, double time_step_length_exh, double time_step_length_int)
{
	if(Ppt.SHOW_calls) cout << "RecordTappings" << endl;

	int a, p;
	// Interpolate and record data for each tapping location of each domain
	for(a=0; a<Ppt.NASSEMBLY; ++a)
	{
		// Pipes
		// ----------------------------------------------------------------------------------------------------
		for(p=0; p<Assembly[a].NEXPIPES; ++p) Assembly[a].Exhaust[p].RecordTappings(&Ppt, timestep, test_var);
		for(p=0; p<Assembly[a].NINPIPES; ++p) Assembly[a].Intake[p].RecordTappings(&Ppt, timestep, test_var);

		// Anechoic buffers and dampers
		// ----------------------------------------------------------------------------------------------------
		for(p=0; p<Assembly[a].NEXANEC; ++p)
		{
			Assembly[a].EX_ANEC[p].Buffer.RecordTappings(&Ppt, timestep, test_var);
			Assembly[a].EX_ANEC[p].Damper.RecordTappings(&Ppt, timestep, test_var);
		}
		for(p=0; p<Assembly[a].NINANEC; ++p)
		{
			Assembly[a].IN_ANEC[p].Buffer.RecordTappings(&Ppt, timestep, test_var);
			Assembly[a].IN_ANEC[p].Damper.RecordTappings(&Ppt, timestep, test_var);
		}
	}

	// Interpolate and record data for each tapping location for each transmissive boundary
	for(a=0; a<Ppt.NASSEMBLY; ++a)
	{
		for(p=0; p<Assembly[a].NEXTRANSM; ++p) Assembly[a].TransmE[p].RecordTappings(&Ppt, timestep, test_var, Assembly[a].Exhaust);
		for(p=0; p<Assembly[a].NINTRANSM; ++p) Assembly[a].TransmI[p].RecordTappings(&Ppt, timestep, test_var, Assembly[a].Intake);
//		for(p=0; p<Assembly[a].NEXTRANSM; ++p) Assembly[a].TransmE[p].RecordTappingsOld(&Ppt, timestep, test_var, Assembly[a].Exhaust);
//		for(p=0; p<Assembly[a].NINTRANSM; ++p) Assembly[a].TransmI[p].RecordTappingsOld(&Ppt, timestep, test_var, Assembly[a].Intake);
	}

	// AVT (exhaust only at moment)
	// Interpolate and record data for the tapping at the AVT target location
	for (a = 0; a < Ppt.NASSEMBLY; ++a) {
		for (p = 0; p < Assembly[a].NCYLS; ++p) {

			for (int ev = 0; ev < Assembly[a].Cyl[p].EnginePtr->NEXVALVES; ++ev) {

				Assembly[a].Cyl[p].ExhaustValve[ev].RecordAVT(&Ppt, timestep, test_var, time_step_length_exh);// , Assembly[a].Exhaust); // test_var can be time or crank angle
			}

			// NEED TO ADD FOR INTAKE VALVES
		}
	}
}

void PrintToFiles()
{
	if(Ppt.SHOW_calls) cout << "PrintToFiles" << endl;
	
	int a; // assembly
	int b; // boundary
	int p; // pipe

	for(a=0; a<Ppt.NASSEMBLY; ++a)
	{
		// Pipes - tappings and movie data
		// ----------------------------------------------------------------------------------------------------
		for(p=0; p<Assembly[a].NEXPIPES; ++p)
		{
			Assembly[a].Exhaust[p].PrintToFileLocations(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
			Assembly[a].Exhaust[p].PrintToFileMovie(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
			if(Assembly[a].Exhaust[p].METHOD==3) Assembly[a].Exhaust[p].PrintToFileFandE(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca); 
		}
		for(p=0; p<Assembly[a].NINPIPES; ++p) 
		{
			Assembly[a].Intake[p].PrintToFileLocations(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
			Assembly[a].Intake[p].PrintToFileMovie(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
			if(Assembly[a].Intake[p].METHOD==3) Assembly[a].Intake[p].PrintToFileFandE(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca); 
		}

		// Closed ends
		// ----------------------------------------------------------------------------------------------------
		for (b = 0; b < Assembly[a].NEXENDCAP; ++b) {
			Assembly[a].CE[b].PrintToFile(&Ppt, timestep, TIMEe[R + 1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
		}
		for (b = 0; b < Assembly[a].NINENDCAP; ++b) {
			Assembly[a].CI[b].PrintToFile(&Ppt, timestep, TIMEe[R + 1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
		}

		// End environments
		// ----------------------------------------------------------------------------------------------------
		for(b=0; b<Assembly[a].NEXEND; ++b) {
			Assembly[a].EX_END[b].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca, Assembly[a].Exhaust);		
		}
		for(b=0; b<Assembly[a].NINEND; ++b) {
			Assembly[a].IN_END[b].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca, Assembly[a].Intake);		
		}

		// Anechoic buffers and dampers - tappings and movie data
		// ----------------------------------------------------------------------------------------------------
		for(p=0; p<Assembly[a].NEXANEC; ++p) {
			// Buffers
			Assembly[a].EX_ANEC[p].Buffer.PrintToFileLocations(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
			Assembly[a].EX_ANEC[p].Buffer.PrintToFileMovie(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
			if(Assembly[a].EX_ANEC[p].Buffer.METHOD==3) Assembly[a].EX_ANEC[p].Buffer.PrintToFileFandE(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca); 
		
			// Dampers
			Assembly[a].EX_ANEC[p].Damper.PrintToFileLocations(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
			Assembly[a].EX_ANEC[p].Damper.PrintToFileMovie(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
			if(Assembly[a].EX_ANEC[p].Damper.METHOD==3) Assembly[a].EX_ANEC[p].Damper.PrintToFileFandE(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca); 	
		}
		for(p=0; p<Assembly[a].NINANEC; ++p) {
			// Buffers
			Assembly[a].IN_ANEC[p].Buffer.PrintToFileLocations(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
			Assembly[a].IN_ANEC[p].Buffer.PrintToFileMovie(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
			if(Assembly[a].IN_ANEC[p].Buffer.METHOD==3) Assembly[a].IN_ANEC[p].Buffer.PrintToFileFandE(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca); 
		
			// Dampers
			Assembly[a].IN_ANEC[p].Damper.PrintToFileLocations(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
			Assembly[a].IN_ANEC[p].Damper.PrintToFileMovie(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
			if(Assembly[a].IN_ANEC[p].Damper.METHOD==3) Assembly[a].IN_ANEC[p].Damper.PrintToFileFandE(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca); 
		}
		
		// Volumes and valves
		for(p=0; p<Assembly[a].NEXVOLUMES; ++p) {
			Assembly[a].ExVolume[p].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
		}
		for(p=0; p<Assembly[a].NINVOLUMES; ++p) {
			Assembly[a].InVolume[p].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
		}
		
		// Cylinders and valves
		for(p=0; p<Assembly[a].NCYLS; ++p) {
			// Cylinders
			Assembly[a].Cyl[p].PrintToFileMovie(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);

			// Valves		
			int v;
			for(v=0; v<Assembly[a].Cyl[p].EnginePtr->NEXVALVES; ++v)
				Assembly[a].Cyl[p].ExhaustValve[v].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
		
			for(v=0; v<Assembly[a].Cyl[p].EnginePtr->NINVALVES; ++v)
				Assembly[a].Cyl[p].IntakeValve[v].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed, Assembly[0].Eng[0].ca);
		}
		
		// Junctions
		for(p=0; p<Assembly[a].NEXJUNCS; ++p) Assembly[a].JE[p].PrintToFile(TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
		for(p=0; p<Assembly[a].NINJUNCS; ++p) Assembly[a].JI[p].PrintToFile(TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);	
		/*
		// Print flow type to file
		for(f=0; f<Ppt.NEXJUNCS; ++f)
		{
			//fprintf(Ppt.EX_J_FILE[f],"%f\t%i\n", TIMEe[R+1], JE[f].flow_type_winterbone);
			if(fabs(JE[f].K_loss[0][0] - 9) < 1e-6) // If we are getting K9
				fprintf(Ppt.EX_J_FILE[f],"%f\t%f\n", JE[f].W[0], JE[f].K_loss[0][1]);
		}
		*/

		// APLDev boundaries
		for(p=0; p<Assembly[a].NEXAPLDev; ++p) Assembly[a].APLDevE[p].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
		for(p=0; p<Assembly[a].NINAPLDev; ++p) Assembly[a].APLDevI[p].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
		
		// Turbines
		for(p=0; p<Assembly[a].NTURBINE; ++p) {
			Assembly[a].Turbine[p].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
			Assembly[a].Turbine[p].PrintToMovieFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
		}

		// Transmissive boundaries
		for(p=0; p<Assembly[a].NEXTRANSM; ++p) Assembly[a].TransmE[p].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
		for(p=0; p<Assembly[a].NINTRANSM; ++p) Assembly[a].TransmI[p].PrintToFile(&Ppt, timestep, TIMEe[R+1], Assembly[0].Eng[0].ca_elapsed);
	}
}

void CloseFiles()
{
	if(Ppt.SHOW_calls) cout << "CloseFiles" << endl;

	int a, p;
	for(a=0; a<Ppt.NASSEMBLY; ++a)
	{
		// Pipes
		for(p=0; p<Assembly[a].NEXPIPES; ++p) Assembly[a].Exhaust[p].CloseFiles(&Ppt);
		for(p=0; p<Assembly[a].NINPIPES; ++p) Assembly[a].Intake[p].CloseFiles(&Ppt);
		
		// Volumes and valves
		for(p=0; p<Assembly[a].NEXVOLUMES; ++p) Assembly[a].ExVolume[p].CloseFiles();
		for(p=0; p<Assembly[a].NINVOLUMES; ++p) Assembly[a].InVolume[p].CloseFiles();

		// Cylinders and valves
		for(p=0; p<Assembly[a].NCYLS; ++p)
		{
			// Cylinders
			Assembly[a].Cyl[p].CloseFiles();

			// Valves		
			int v;
			for(v=0; v<Assembly[a].Cyl[p].EnginePtr->NEXVALVES; ++v)
				Assembly[a].Cyl[p].ExhaustValve[v].CloseFiles(&Ppt);

			for(v=0; v<Assembly[a].Cyl[p].EnginePtr->NINVALVES; ++v)
				Assembly[a].Cyl[p].IntakeValve[v].CloseFiles(&Ppt);
		}
		
		// Junctions
		for(p=0; p<Assembly[a].NEXJUNCS; ++p) Assembly[a].JE[p].CloseFiles(&Ppt);
		for(p=0; p<Assembly[a].NINJUNCS; ++p) Assembly[a].JI[p].CloseFiles(&Ppt);
		
		// APLDev boundaries
		for(p=0; p<Assembly[a].NEXAPLDev; ++p) Assembly[a].APLDevE[p].CloseFiles(&Ppt);
		for(p=0; p<Assembly[a].NINAPLDev; ++p) Assembly[a].APLDevI[p].CloseFiles(&Ppt);

		// End environments
		for(p=0; p<Assembly[a].NEXEND; ++p) Assembly[a].EX_END[p].CloseFiles(&Ppt);
		for(p=0; p<Assembly[a].NINEND; ++p) Assembly[a].IN_END[p].CloseFiles(&Ppt);

		// Turbines
		for(p=0; p<Assembly[a].NTURBINE; ++p) Assembly[a].Turbine[p].CloseFiles(&Ppt);
		
		// Transmissive boundaries
		for(p=0; p<Assembly[a].NEXTRANSM; ++p) Assembly[a].TransmE[p].WriteFiles();
		for(p=0; p<Assembly[a].NINTRANSM; ++p) Assembly[a].TransmI[p].WriteFiles();
		
		for(p=0; p<Assembly[a].NEXTRANSM; ++p) Assembly[a].TransmE[p].CloseFiles();
		for(p=0; p<Assembly[a].NINTRANSM; ++p) Assembly[a].TransmI[p].CloseFiles();
	}
	// Close simulation case .out file
	fclose(Ppt.OUTFILE);
}

///*
void FillingAndEmptying(double DELZe, double DELZi)
{
	int a, p, end;

	double start_beta = 0.9;
	double start_tol = 1;//1e-6;//1e-0;//1e-6;// %

	// Reset each pipe end CONVERGED flag
	// ==================================
	for(a=0; a<Ppt.NASSEMBLY; ++a)
	{
		for(p=0; p<Assembly[a].NEXPIPES; ++p)
		{
			for(end=0; end<2; ++end)
			{
				Assembly[a].Exhaust[p].A_throat[end] = 1e6; // To avoid convergence on first loop
				Assembly[a].Exhaust[p].U_throat[end] = 1e6; // To avoid convergence on first loop
				Assembly[a].Exhaust[p].CONVERGED[end] = false;
				Assembly[a].Exhaust[p].beta[end] = start_beta;
				Assembly[a].Exhaust[p].tol[end] = start_tol;
				Assembly[a].Exhaust[p].direction_str[end] = "NOFLOW";
			}
		}
		
		for(p=0; p<Assembly[a].NINPIPES; ++p)
		{
			for(end=0; end<2; ++end)
			{
				Assembly[a].Intake[p].A_throat[end] = 1e6; // To avoid convergence on first loop
				Assembly[a].Intake[p].U_throat[end] = 1e6; // To avoid convergence on first loop
				Assembly[a].Intake[p].CONVERGED[end] = false;
				Assembly[a].Intake[p].beta[end] = start_beta;
				Assembly[a].Intake[p].tol[end] = start_tol;
				Assembly[a].Intake[p].direction_str[end] = "NOFLOW";
			}
		}
	}

	int COUNTER = 0;
	bool CONVERGED;
	do
	{
		CONVERGED = true;

		// Boundary calculations for all assemblies
		// ==========================================================================================
		for(a=0; a<Ppt.NASSEMBLY; ++a) Assembly[a].RunBoundary(&Ppt, MyTime, DELZe, DELZi, TIMEe, TIMEi, timestep, RESTORE);

		// Refine incoming characteristic at boundaries
		// ==========================================================================================
		for(a=0; a<Ppt.NASSEMBLY; ++a)
		{
			for(p=0; p<Assembly[a].NEXPIPES; ++p)
			{
				//if(!(Assembly[a].Exhaust[p].ODD_CONVERGED && Assembly[a].Exhaust[p].EVEN_CONVERGED))
				if(!(Assembly[a].Exhaust[p].CONVERGED[0] && Assembly[a].Exhaust[p].CONVERGED[1]))
				{
					Assembly[a].Exhaust[p].FillingAndEmptyingAdjust(&Ppt, DELZe, COUNTER, TIMEe[1]);
					//Assembly[a].Exhaust[p].FillingAndEmptyingAdjust2(&Ppt, DELZe, COUNTER, TIMEe[1]);
					//Assembly[a].Exhaust[p].FillingAndEmptyingAdjust3(&Ppt, DELZe, COUNTER, TIMEe[1]);
					CONVERGED = false;
				}
			}

			for(p=0; p<Assembly[a].NINPIPES; ++p)
			{
				//if(!(Assembly[a].Intake[p].ODD_CONVERGED && Assembly[a].Intake[p].EVEN_CONVERGED))
				if(!(Assembly[a].Intake[p].CONVERGED[0] && Assembly[a].Intake[p].CONVERGED[1]))
				{
					Assembly[a].Intake[p].FillingAndEmptyingAdjust(&Ppt, DELZi, COUNTER, TIMEi[1]);
					//Assembly[a].Intake[p].FillingAndEmptyingAdjust2(&Ppt, DELZi, COUNTER, TIMEi[1]);
					//Assembly[a].Intake[p].FillingAndEmptyingAdjust3(&Ppt, DELZi, COUNTER, TIMEi[1]);
					CONVERGED = false;
				}
			}				
		}
		++COUNTER;

		//cout << "COUNTER = " << COUNTER << endl;
	}while(!CONVERGED);
}
//*/
void ResetInitialConditions()
// Re-applies initial gas flow conditions from input file to all pipes in all assemblies
{
	for(int a=0; a<Ppt.NASSEMBLY; ++a)
	{
		int p;
		for(p=0; p<Assembly[a].NEXPIPES; ++p) Assembly[a].Exhaust[p].InitialConditions(&Ppt);
		for(p=0; p<Assembly[a].NINPIPES; ++p) Assembly[a].Intake[p].InitialConditions(&Ppt);
	}
	cout << "Resetting initial conditions" << endl << endl;
}