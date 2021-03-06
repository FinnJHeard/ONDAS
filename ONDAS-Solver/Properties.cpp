// Properties.cpp: implementation of the CProperties class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Properties.h"
#include "Tools.h"

#include <string.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CProperties::CProperties()
{
	// Initialise labels for friction factor methods
	CONSTANT_FF = 0;
	SWAMEEJAIN_FF = 1;
	HAALAND_FF = 2;
	COLEBROOK_FF = 3;
	RICARDO_FF = 4;

	// Initialise labels for viscosity methods
	CONSTANT_VISCOSITY = 0;
	SUTHERLAND_VISCOSITY = 1;
	BLAIR_VISCOSITY = 2;
	RICARDO_VISCOSITY = 3;

	// Initialise labels for heat transfer methods
	REYNOLDS_ANALOGY = 0;
	NUSSELT_RELATION = 1;
	RICARDO_HT = 2;

	K_entropy = 1;
}

CProperties::~CProperties()
{
	delete [] labels;
	delete [] values;
}

void CProperties::Out(std::string text)
{
	cout << text;
	fprintf(OUTFILE,"%s", text.c_str());
}

void CProperties::Out(char* text)
{
	cout << text;
	fprintf(OUTFILE,"%s", text);
}

void CProperties::Out(const char* text)
{
	cout << text;
	fprintf(OUTFILE,"%s", text);
}

void CProperties::Out(char text)
{
	cout << text;
	fprintf(OUTFILE,"%c", text);
}

void CProperties::Out(int text)
{
	cout << text;
	fprintf(OUTFILE,"%d", text);
}

void CProperties::Out(double text)
{
	cout << text;
	fprintf(OUTFILE,"%f", text);
}

void CProperties::Out(double text, int precision)
{
	int temp_precision = cout.precision(); 
	cout << setprecision(precision);
	cout << text;
	//fprintf(OUTFILE,"%f", text);
	fprintf(OUTFILE, ConstructString(this, ConstructString(this, "%.", IntToString(precision)), "f"), text);
	cout << setprecision(temp_precision);
}

void CProperties::ReadCaseName(char* InputFile, int num_parameters)
{
	int last_entry;
	CommonReadInput(InputFile, num_parameters, this->labels, this->values, this->strings, last_entry);
	for(int r=0; r<last_entry+1; ++r) if(strcmp(labels[r], "case_name") == 0) case_name = strings[r];

	std::string s;
	s = case_name;
	
	std::string sz;
	for(int i=0; i<int(s.length()); ++i) sz += s[i]; 
	sz += "\\"; // Add directory slash

	case_name_slash = new char[sz.length() + 1];
	strcpy(case_name_slash, sz.c_str());
}

void CProperties::ReadInput(char* InputFile, int num_parameters)
{
    int last_entry;
	CommonReadInput(InputFile, num_parameters, this->labels, this->values, this->strings, last_entry);

	// Default settings in case not specified
	strFileExt = ".orf";
	rho_blood = 1060;
	ViscosityBlood = 3e-3;
	k1 = 2e6;
	k2 = -22.53e2;
	k3 = 8.65e4;
	BL_RELATIVE = 1;
	BL_rel = 0.2;
	BL_abs = 1e-3;
	
	for(int r=0; r<last_entry+1; ++r)
	{
		// ====================================================================================================
		// Main input parameter file
		// ====================================================================================================

		// Case description
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "strDesc") == 0) strDesc = strings[r];

		// Code version
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "FILE_VERSION") == 0) FILE_VERSION = int(values[r]);
		
		// Parameter file control
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NUM_PARAMS") == 0) NUM_PARAMS = int(values[r]);
		if(strcmp(labels[r], "SHOW_globals") == 0) SHOW_globals = DoubleToBool(values[r]);

		// Simulation control
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "ZMAX") == 0) ZMAX = values[r];
		if(strcmp(labels[r], "courant") == 0) courant = values[r];
		if(strcmp(labels[r], "delt_min") == 0) delt_min = values[r];
		if(strcmp(labels[r], "delt_max") == 0) delt_max = values[r];
		if(strcmp(labels[r], "CONSTANT_DELT") == 0) CONSTANT_DELT = DoubleToBool(values[r]);
		if(strcmp(labels[r], "CONTINUOUS") == 0) CONTINUOUS = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SUPERSONIC") == 0) SUPERSONIC = DoubleToBool(values[r]);
		if(strcmp(labels[r], "BEEP") == 0) BEEP = DoubleToBool(values[r]);
		if(strcmp(labels[r], "HAEMODYNAMICS") == 0) HAEMODYNAMICS= DoubleToBool(values[r]);
		
		// Propagation method
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "DEF_METHOD") == 0) DEF_METHOD = int(values[r]);
		
		// If mesh method of characteristics (MMOC), i.e., DEF_METHOD = 1 
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "HOMENTROPIC") == 0) HOMENTROPIC = DoubleToBool(values[r]);
		
		// If non-homentropic MMOC, i.e., HOMENTROPIC = 0
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NUM_PATH_MULT") == 0) NUM_PATH_MULT = int(values[r]);
		
		// If W(alpha,beta) scheme, i.e., DEF_METHOD = 2: MacCormack (1,0), Two-step L-W (0.5,0.5), Lerat & Peyret optimum (1+sqrt(5)/2,0.5)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "alpha") == 0) alpha = values[r];
		if(strcmp(labels[r], "beta") == 0) beta = values[r];
		if(strcmp(labels[r], "COMBINED_WAB_MOC") == 0) COMBINED_WAB_MOC = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SOURCES") == 0) SOURCES = DoubleToBool(values[r]);
		if(strcmp(labels[r], "WORK") == 0) WORK = DoubleToBool(values[r]);
		if(strcmp(labels[r], "ALTERNATE_MAC") == 0) ALTERNATE_MAC = DoubleToBool(values[r]);
		if(strcmp(labels[r], "VISC") == 0) VISC = DoubleToBool(values[r]);
		if(strcmp(labels[r], "Cx") == 0) Cx = values[r];
		if(strcmp(labels[r], "Cx_alpha") == 0) Cx_alpha = values[r];
		if(strcmp(labels[r], "TVD") == 0) TVD = DoubleToBool(values[r]);

		// (3) Filling and emptying
		// ----------------------------------------------------------------------------------------------------
		
		// Gas properties
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "COMPRESSIBLE") == 0) COMPRESSIBLE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "R_universal") == 0) R_universal = values[r];
		if(strcmp(labels[r], "M_air") == 0) M_air = values[r];
		if(strcmp(labels[r], "constProps") == 0) constProps = DoubleToBool(values[r]);
		if(strcmp(labels[r], "Cp_air") == 0) Cp_air = values[r];

		// Source term models
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "FRICTION_MODEL") == 0) FRICTION_MODEL = int(values[r]);
		if(strcmp(labels[r], "ff_const") == 0) ff_const = values[r];
		if(strcmp(labels[r], "VISCOSITY_AIR") == 0) VISCOSITY_AIR = int(values[r]);
		if(strcmp(labels[r], "mu_air_const") == 0) mu_air_const = values[r];
		if(strcmp(labels[r], "HEAT_TRANSFER") == 0) HEAT_TRANSFER = int(values[r]);
		if(strcmp(labels[r], "tol_steady") == 0) tol_steady = values[r];

		// Haemodynamic properties (k values from Olufsen (1999), "Structured tree outflow condition for blood flow in larger systemic arteries")
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "rho_blood") == 0) rho_blood = values[r];
		if (strcmp(labels[r], "ViscosityBlood") == 0) ViscosityBlood = values[r];
		if (strcmp(labels[r], "k1") == 0) k1 = values[r];
		if (strcmp(labels[r], "k2") == 0) k2 = values[r];
		if (strcmp(labels[r], "k3") == 0) k3 = values[r];
		if (strcmp(labels[r], "BL_RELATIVE") == 0) BL_RELATIVE = DoubleToBool(values[r]);

		// Boundary layer thickness 
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "BL_rel") == 0) BL_rel = values[r];
		if (strcmp(labels[r], "BL_abs") == 0) BL_abs = values[r];

		// Default values
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "discret") == 0) discret = values[r]/1000;
		if(strcmp(labels[r], "min_meshes") == 0) min_meshes = int(values[r]);
		if(strcmp(labels[r], "epsilon") == 0) epsilon = values[r]/1000;
		if(strcmp(labels[r], "CFTRANS") == 0) CFTRANS = values[r];
		if(strcmp(labels[r], "HGTRANS") == 0) HGTRANS = values[r];
		if(strcmp(labels[r], "CPTRANS") == 0) CPTRANS = values[r];
		if(strcmp(labels[r], "freq") == 0) freq = int(values[r]);

		// Reference values
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "PREF") == 0) PREF = values[r];
		if(strcmp(labels[r], "xref") == 0) xref = values[r];
		if(strcmp(labels[r], "fref") == 0) fref = values[r];

		// Boundary conditions
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "PARENT") == 0) PARENT = int(values[r]);
		if(strcmp(labels[r], "STRLEN") == 0) STRLEN = int(values[r]);
		if(strcmp(labels[r], "ZERO_TOL") == 0) ZERO_TOL = values[r];

		// Screen output
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "outputInterval") == 0) outputInterval = values[r];
		if(strcmp(labels[r], "SHOW_conn") == 0) SHOW_conn = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_config") == 0) SHOW_config = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_params") == 0) SHOW_params = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_calls") == 0) SHOW_calls = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_x") == 0) SHOW_x = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_X") == 0) SHOW_X = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_DELX") == 0) SHOW_DELX = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_d") == 0) SHOW_d = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_dddx") == 0) SHOW_dddx = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_d2ddx2") == 0) SHOW_d2ddx2 = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_f") == 0) SHOW_f = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_cfa") == 0) SHOW_cfa = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_cfa_delx") == 0) SHOW_cfa_delx = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_dfdx") == 0) SHOW_dfdx = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_W") == 0) SHOW_W = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_F") == 0) SHOW_F = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_C") == 0) SHOW_C = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_W_pred") == 0) SHOW_W_pred = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_F_pred") == 0) SHOW_F_pred = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_C_pred") == 0) SHOW_C_pred = DoubleToBool(values[r]);
		if(strcmp(labels[r], "SHOW_pathlines") == 0) SHOW_pathlines = DoubleToBool(values[r]);

		// File output
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "strFileExt") == 0) {
			strFileExt = ".";
			strFileExt += strings[r];
		}

		// Engines, cylinders & valves
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "CYL_DIR") == 0) CYL_DIR = strings[r];
		if(strcmp(labels[r], "POPPET_H_TOL1") == 0) POPPET_H_TOL1 = values[r];
		if(strcmp(labels[r], "POPPET_H_TOL2") == 0) POPPET_H_TOL2 = values[r];
		if(strcmp(labels[r], "VT_DIR") == 0) VT_DIR = strings[r];

		// End environments
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "ENDENV_DIR") == 0) ENDENV_DIR = strings[r];
		if(strcmp(labels[r], "USE_PHI") == 0) USE_PHI = DoubleToBool(values[r]);
		if(strcmp(labels[r], "NHI_TOL") == 0) NHI_TOL = values[r];
		if(strcmp(labels[r], "NOZZ_TOL") == 0) NOZZ_TOL = values[r];

		// Junctions
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "JUNC_DIR") == 0) JUNC_DIR = strings[r];
		if(strcmp(labels[r], "max_branches") == 0) max_branches = int(values[r]);

		// Adiabatic pressure loss devices
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "APLD_DIR") == 0) APLD_DIR = strings[r];
		if(strcmp(labels[r], "LOSS_TOL") == 0) LOSS_TOL = values[r];

		// Sudden area changes
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "tolSuddenEnlrgN") == 0) tolSuddenEnlrgN = values[r];
		if(strcmp(labels[r], "tolSuddenEnlrgLin") == 0) tolSuddenEnlrgLin = values[r];

		// Turbines
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "TURB_DIR") == 0) TURB_DIR = strings[r];

		// Transmissive boundaries
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "TRANSM_DIR") == 0) TRANSM_DIR = strings[r];

		// Assemblies
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "NASSEMBLY") == 0) NASSEMBLY = int(values[r]);
	}

	// Checks
	// ----------------------------------------------------------------------------------------------------
	if(delt_max < delt_min) {
		Out("delt_max < delt_min, exiting...\n"); exit(1);
	}
	if(DEF_METHOD!=1) HOMENTROPIC = false; // If not using MMOC must also be non-homentropic
	if(DEF_METHOD==1 && HOMENTROPIC && SHOW_pathlines) SHOW_pathlines = false; // Pathlines not present in homentropic MMOC
	if(DEF_METHOD==2 && TVD && courant>0.7){
		SimulationWarning(this, "A maximum Courant number of 0.7 will be enforced when using TVD criteria");
		courant=0.7;
	}
	SetParams(InputFile);
	tol_steady_multiplier = 1.0;
}

void CProperties::SetParams(char* InputFile)
{	
    if(SHOW_calls){Out("Ppt.SetParams("); Out(InputFile); Out(")\n");}
	// Gas properties
	// ----------------------------------------------------------------------------------------------------
	R_air = R_universal/M_air;
	if(constProps){ // Use value of Cp_air provided by user
		Cv_air = Cp_air - R_air;	
		gamma_air = Cp_air/Cv_air;
	}
	else
	{
		Cp_air = cpAir(TREFe);
		Cv_air = cvAir(TREFe);
		gamma_air = gammaAir(TREFe);
	}
      
    // Initial values in case of no pipes
    // ----------------------------------------------------------------------------------------------------
	TREFe = 300;
	TREFi = 300;
	AREFe = sqrt(gammaAir(TREFe)*R_air*TREFe);
	AREFi = sqrt(gammaAir(TREFi)*R_air*TREFi);

	// Constants
    // ----------------------------------------------------------------------------------------------------
	AA			= (3.0-gammaAir(TREFe))/(2.0*(gammaAir(TREFe)-1.0));
	BB			= (gammaAir(TREFe)+1.0)/(2.0*(gammaAir(TREFe)-1.0));
	Q			= (gammaAir(TREFe)-1.0)/(2.0*gammaAir(TREFe));
	QI			= 1.0/Q;
//	Eta			= (gammaAir()-1.0)/(2.0*gammaAir());
	sigma = 5.67e-8;		// Stefan-Boltzman constant (for cylinder radiative heat transfer)
	
	// Now run some checks on the compatibility of the above numbers

	// Multiple the junction configuration string length per pipe, by the number of pipes
//	JUNSTRLEN = JUNSTRLEN * max_branches;
/*	
	// Make sure enough pipes exist to initialise all the junctions
	if(NEXPIPES<2*NEXJUNCS + 1 && NEXJUNCS!=0)
	{
		cout << NEXPIPES << " exhaust pipes is insufficient to set up " << NEXJUNCS << " exhaust junctions" << endl;
		cout << 2*NEXJUNCS + 1 << " is the minimum required number of exhaust pipes for this configuration" << endl;
		cout << "Adjust " << InputFile << " to satisfy these requirements, and re-run" << endl;
		exit(1);
	}
	if(NINPIPES<2*NINJUNCS + 1 && NINJUNCS!=0)
	{
		cout << NINPIPES << " intake pipes is insufficient to set up " << NINJUNCS << " intake junctions" << endl;
		cout << 2*NINJUNCS + 1 << " is the minimum required number of intake pipes for this configuration" << endl;
		cout << "Adjust " << InputFile << " to satisfy these requirements, and re-run" << endl;
		exit(1);
	}
*/
/*
	// If simulation type is set to periodic, but there are no cylinders, reset it to continuous
	if(NCYLS==0) {
		CONTINUOUS = true;
	}
*/
	if(HOMENTROPIC) {
		// The reference temperature defines the temperature in the pipe
//		TREFe = TRIe*pow(PREF/PRIe, 2*Q);
//		TREFi = TRIi*pow(PREF/PRIi, 2*Q);
		// This will then ignore any TREF set in the input file!
	}

	// Read junction loss coefficients (should move this so that it is called only when junctions are invoked)
	ReadLossCoefficients(ConstructString(this, ConstructString(this, this->param_dir, this->JUNC_DIR), "losscoeffs.txt"));
}

void CProperties::ListProperties(int CODE_VERSION)
{
	if(SHOW_calls){Out("Ppt.ListProperties\n");}

	// ====================================================================================================
	// Main input parameter file
	// ====================================================================================================

	// Case description
	// ----------------------------------------------------------------------------------------------------
	std::string case_str;
	case_str = "Running Case: ";
	case_str += case_name;
	char* border_char; border_char = new char [2]; 
	border_char[0] = '=';//176;//254;
	border_char[1] = '\0';
	Out("\n");
	//Out(Underline(StringToChar(case_str), border_char, true));
	Out(Underline(StringToChar(case_str), border_char, true, 80, 1));
	Out("\n");
	Out("\t");
	int c=0;
	bool SPACE = false;
	while(strDesc[c]!=NULL)
	{
		if(strDesc[c]=='_') {
			if(SPACE) {
				Out("\n\t"); 
				SPACE = false;
			}
			else Out(' '); 
		}
		else Out(strDesc[c]);
		++c;
		if(c%60==0) SPACE = true; // Insert newline roughly every 60 characters
	}
	Out("\n\n\n");

	Out(Underline("Global parameters", "=", true, 80, 1));
	Out("\n");

	// Code version
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Code version", "-", "\t"));
	if(CODE_VERSION == FILE_VERSION) {
		Out("\tCode version ("); Out(CODE_VERSION); Out(") matches file version ("); Out(FILE_VERSION); Out(")\n");
	}
	else {
		if(CODE_VERSION > FILE_VERSION) {
			//cout << "\tCode version (" << CODE_VERSION << ") newer than file version (" << FILE_VERSION << "). Exiting." << endl;
			Out("\tCode version ("); Out(CODE_VERSION); Out(") newer than file version ("); Out(FILE_VERSION); Out("). Exiting.\n");
			exit(1);
		}
		else {
			//cout << "\tCode version (" << CODE_VERSION << ") older than file version (" << FILE_VERSION << "). Exiting." << endl;
			Out("\tCode version ("); Out(CODE_VERSION); Out(") older than file version ("); Out(FILE_VERSION); Out("). Exiting.\n");
			exit(1);
		}
	}
	Out("\n");

	// Parameter file control
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Parameter file control", "-", "\t"));
	Out("\tMaximum number of parameters, NUM_PARAMS\t=\t"); Out(NUM_PARAMS); Out("\n");
	if(SHOW_globals) {
		Out("\tShowing global parameters, SHOW_globals\t\t=\t"); Out(TrueOrFalse(SHOW_globals)); Out("\n");
	}
	else {
		Out("\tNot showing global parameters, SHOW_globals\t=\t"); Out(TrueOrFalse(SHOW_globals)); Out("\n");
	}
	Out("\n");

	// Simulation control
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Simulation control", "-", "\t"));
	Out("\tMaximum simulation time, ZMAX\t\t\t=\t"); Out(ZMAX); Out(" s\n");

	if(!CONSTANT_DELT) {
		Out("\tVariable timestep; Courant number, courant\t=\t"); Out(courant); Out("\n");
		Out("\tMinimum and first time step, delt_min\t\t=\t"); Out(delt_min); Out(" s\n");
		Out("\tMaximum timestep, delt_max\t\t\t=\t"); Out(delt_max); Out(" s\n");
	}
	else {
		Out("\tMinimum and first time step, delt_min\t\t=\t"); Out(delt_min); Out(" s\n");
		Out("\tOtherwise constant timestep, delt_max\t\t=\t"); Out(delt_max); Out(" s\n");
	}
	
	Out("\tSimulation time axis type, CONTINUOUS\t\t=\t"); Out(TrueOrFalse(CONTINUOUS)); Out("\n");
	if(CONTINUOUS) {
		Out("\t- continuous simulation; progress measured in s\n");
	}
	else {
		Out("\t- periodic simulation; progress measured in "); Out(Deg()); Out("CA\n");
	}
	
	if(SUPERSONIC) {
		Out("\tSupersonic flow permitted, SUPERSONIC\t\t=\t"); Out(TrueOrFalse(SUPERSONIC)); Out("\n");
	}
	else {
		Out("\tSupersonic flow not permitted, SUPERSONIC\t=\t"); Out(TrueOrFalse(SUPERSONIC)); Out("\n");
	}

	if(BEEP){Out("\tBeep for user input, BEEP\t\t\t=\t"); Out(TrueOrFalse(BEEP)); Out("\n");}
	else{Out("\tNo beep for user input, BEEP\t\t\t=\t"); Out(TrueOrFalse(BEEP)); Out("\n");}

	if (HAEMODYNAMICS) { Out("\tUsing haemodynamic models, HAEMODYNAMICS\t=\t"); Out(TrueOrFalse(HAEMODYNAMICS)); Out("\n"); }
	else { Out("\tNot using haemodynamic models, HAEMODYNAMICS\t=\t"); Out(TrueOrFalse(HAEMODYNAMICS)); Out("\n"); }

	if(SHOW_calls) {Out("\tEchoing function calls, SHOW_calls\t\t=\t"); Out(TrueOrFalse(SHOW_calls)); Out("\n");}
	else {Out("\tNot echoing function calls, SHOW_calls\t\t=\t"); Out(TrueOrFalse(SHOW_calls)); Out("\n");}

	Out("\n");

	// Propagation method
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Propagation method (default - individual domains may use different scheme)", "-", "\t"));

	// If mesh method of characteristics (MMOC), i.e., DEF_METHOD = 1 
	// ----------------------------------------------------------------------------------------------------
	if(DEF_METHOD==1) {
		Out("\tMesh Method of Characteristics (");
		if(HOMENTROPIC) Out("homentropic)\n");
		else {	
			// If non-homentropic MMOC, i.e., HOMENTROPIC = 0
			// ----------------------------------------------------------------------------------------------------
			Out("non-homentropic)\n");
			Out("\tPathlines multiplier, NUM_PATH_MULT\t\t=\t"); Out(NUM_PATH_MULT); Out("\n");
		}
	}
	else {
		// If W(alpha,beta) scheme, i.e., DEF_METHOD = 2: MacCormack (1,0), Two-step L-W (0.5,0.5), Lerat & Peyret optimum (1+sqrt(5)/2,0.5)
		// ----------------------------------------------------------------------------------------------------
		if(DEF_METHOD==2)
		{	
			// Character codes: 224 = alpha, 225 = beta
			// Lucida console: αβ 
			Out("\tW"); Out('a'); Out(char(225)); Out(" scheme ("); Out('a'); Out("="); Out(alpha); Out(", "); Out(char(225)); Out("="); Out(beta); Out(")"); 
			if(alpha==0.5 && beta==0.5) Out("\t\t\t=\tTwo-step Lax-Wendroff");
			else {
				if(alpha==1 && beta==0) {
					Out("\t\t\t\t=\tMacCormack method");
					if(ALTERNATE_MAC) Out(" (alternate)");
				}
				else {
					if(fabs(alpha-(1 + sqrt(5.)/2))<1e-6 && beta==0.5) Out("\t\t\t=\tLerat & Peyret \"optimum\"");
					else {Out("\t\t\t\t=\tUncategorised W"); Out(char(224)); Out(char(225)); Out(" scheme");}
				}
			}
			Out("\n\t\t\t\t\t\t\t\t");

			if(SOURCES)
			{
				if(WORK) Out(" with source terms (inc. work)");
				else Out(" with source terms (but no work)");
			}			
			else Out(" no source terms (except area variation)");
			Out("\n\t\t\t\t\t\t\t\t");

			if(TVD) Out(" + TVD");
			else Out(" no TVD");
			Out("\n");	

			if(COMBINED_WAB_MOC) {
				Out("\tContinuing to use pathlines, COMBINED_WAB_MOC\t=\t"); Out(TrueOrFalse(COMBINED_WAB_MOC)); Out("\n");
				Out("\tPathlines multiplier, NUM_PATH_MULT\t\t=\t"); Out(NUM_PATH_MULT); Out("\n");
			}
			else {Out("\tExtrapolating solution vector, COMBINED_WAB_MOC\t=\t"); Out(TrueOrFalse(COMBINED_WAB_MOC)); Out("\n");}
			
			if(VISC) {
				Out("\tApplying artificial viscosity, VISC\t\t=\t"); Out(TrueOrFalse(VISC)); Out("\n");
				Out("\tArtificial viscosity coefficient, Cx\t\t=\t"); Out(Cx); Out("\n");
				Out("\tArtificial viscosity coefficient, Cx_alpha\t=\t"); Out(Cx_alpha); Out("\n");
			}
			else {Out("\tNo artificial viscosity, VISC\t\t\t=\t"); Out(TrueOrFalse(VISC)); Out("\n");}
		}
		else {
			// (3) Filling and emptying
			// ----------------------------------------------------------------------------------------------------
			if(DEF_METHOD==3) Out("\tFilling and emptying\n");
			else Out("unknown method\n");
		}
	}	
	Out("\n");

	if (HAEMODYNAMICS) {
		// Haemodynamic properties (k values from Olufsen (1999), "Structured tree outflow condition for blood flow in larger systemic arteries")
		// ----------------------------------------------------------------------------------------------------
		Out(Underline("Haemodynamic properties (k values from Olufsen (1999), 'Structured tree outflow condition for blood flow in larger systemic arteries')", "-", "\t"));
		Out("\tDensity of blood, rho_blood\t\t\t=\t"); Out(rho_blood); Out(" kg.m^-3\n");
		Out("\tDynamic viscosity of blood, ViscosityBlood\t=\t"); Out(ViscosityBlood); Out(" Pa.s\n");
		Out("\tFitting coefficient k1 for large vessels, k1\t=\t"); Out(k1); Out(" kg.s^-1.m^-1\n");
		Out("\tFitting coefficient k2 for large vessels, k2\t=\t"); Out(k2); Out(" m^-1\n");
		Out("\tFitting coefficient k3 for large vessels, k3\t=\t"); Out(k3); Out(" kg.s^-1.m^-1\n");

		if (BL_RELATIVE) {
			Out("\tUsing relative BL thickness, BL_RELATIVE\t=\t"); Out(TrueOrFalse(BL_RELATIVE)); Out("\n");
			Out("\tBL thickness fraction of diameter, BL_rel\t=\t"); Out(BL_rel); Out(" \n");
		}
		else {
			Out("\tUsing fixed BL thickness, BL_RELATIVE\t\t=\t"); Out(TrueOrFalse(BL_RELATIVE)); Out("\n");
			Out("\tFixed absolute BL thickness, BL_abs\t\t=\t"); Out(BL_abs); Out(" \n");
		}
		Out("\n");
	}
	else {
		// Gas properties 
		// ----------------------------------------------------------------------------------------------------
		Out(Underline("Gas properties", "-", "\t"));
		if (COMPRESSIBLE) { Out("\tUsing compressible relation, COMPRESSIBLE\t=\t"); Out(TrueOrFalse(COMPRESSIBLE)); Out("\n"); }
		else { Out("\tUsing incompressible relation, COMPRESSIBLE\t=\t"); Out(TrueOrFalse(COMPRESSIBLE)); Out("\n"); }
		Out("\tUniversal gas constant, R_universal\t\t=\t"); Out(R_universal); Out(" J.mol^-1.K^-1\n");
		Out("\tMolecular mass of air, M_air\t\t\t=\t"); Out(M_air); Out(" kg.mol^-1\n");
		Out("\tGas constant for air, R_air\t\t\t=\t"); Out(R_air); Out(" J.kg^-1.K^-1\n");
		if (constProps) {
			Out("\tUsing constant gas properties below, constProps\t=\t"); Out(TrueOrFalse(constProps)); Out("\n");
			Out("\tSpec. heat of air, constant pressure, Cp_air\t=\t"); Out(Cp_air); Out(" J.kg^-1.K^-1\n");
			Out("\tSpec. heat of air, constant volume, Cv_air\t=\t"); Out(Cv_air); Out(" J.kg^-1.K^-1\n");
			Out("\tRatio of specific heats, gamma_air\t\t=\t"); Out(gamma_air); Out("\n");
		}
		else { Out("\tUsing Zucrow & Hoffman function, constProps\t=\t"); Out(TrueOrFalse(constProps)); Out("\n"); }
		Out("\n");
	}

	Out(Underline("Source term models", "-", "\t"));
	Out("\tFriction factor model, FRICTION_MODEL\t\t=\t");
	if(FRICTION_MODEL==CONSTANT_FF) {
		Out("CONSTANT_FF, constant friction factor\n");
		Out("\t- const. friction factor, ff_const\t\t=\t"); Out(ff_const); Out("\n");
	}
	else {
		if(FRICTION_MODEL==SWAMEEJAIN_FF) Out("SWAMEEJAIN_FF, Swamee-Jain\n");
		else {
			if(FRICTION_MODEL==HAALAND_FF) Out("HAALAND_FF, Haaland\n");
			else {
				if(FRICTION_MODEL==COLEBROOK_FF) Out("COLEBROOK_FF, Colebrook-White implicit\n");
				else {
					if(FRICTION_MODEL==RICARDO_FF) Out("RICARDO_FF, Ricardo WAVE\n");
					else {
						Out("unknown, use default (SWAMEEJAIN_FF)\n");
						FRICTION_MODEL=SWAMEEJAIN_FF;
					}
				}
			}
		}
	}

	Out("\tAir viscosity function, VISCOSITY_AIR\t\t=\t");
	if(VISCOSITY_AIR==CONSTANT_VISCOSITY) {
		Out("CONSTANT_VISCOSITY, using constant viscosity\n"); 
		Out("\tConstant air viscosity, mu_air_const\t\t=\t"); Out(mu_air_const); Out(" kg.m^-1.s^-1\n");
	}
	else {
		if(VISCOSITY_AIR==SUTHERLAND_VISCOSITY) Out("SUTHERLAND_VISCOSITY, Sutherland's equation\n");
		else {
			if(VISCOSITY_AIR==BLAIR_VISCOSITY) Out("BLAIR_VISCOSITY, Blair's equation\n");
			else {
				if(VISCOSITY_AIR==RICARDO_VISCOSITY) Out("RICARDO_VISCOSITY, Ricardo WAVE equation\n");
				else {
					Out("unknown, use default (SUTHERLAND_VISCOSITY)\n");
					VISCOSITY_AIR=SUTHERLAND_VISCOSITY;
				}
			}
		}
	}

	Out("\tHeat transfer model, HEAT_TRANSFER\t\t=\t");
	if(HEAT_TRANSFER==REYNOLDS_ANALOGY) Out("REYNOLDS_ANALOGY, using Reynolds' analogy for convective heat transfer from pipes\n");
	else {
		if(HEAT_TRANSFER==NUSSELT_RELATION)
			Out("NUSSELT_RELATION, using Nusselt relation for convective heat transfer from pipes\n");
		else {
			if(HEAT_TRANSFER==RICARDO_HT) Out("RICARDO_HT, Ricardo WAVE\n");
			else {
				Out("unknown, use default (REYNOLDS_ANALOGY)\n");
				HEAT_TRANSFER=REYNOLDS_ANALOGY;
			}
		}
	}
	Out("\n");
	Out("\tSteady flow tolerance, tol_steady\t\t=\t"); Out(tol_steady); Out("%\n");
	Out("\n");

	// Default values
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Default values", "-", "\t"));
	Out("\tDefault meshing parameters:\n");
	Out("\t- Target discretization length, discret\t\t=\t"); Out(discret*1000); Out(" mm\n");
	Out("\t- Minimum number of meshes, min_meshes\t\t=\t"); Out(min_meshes); Out("\n");
	Out("\tDefault roughness height and enhancement factors:\n");
	Out("\t- Roughness height, epsilon\t\t\t=\t"); Out(epsilon*1000); Out(" mm\n");
	Out("\t- Friction factor, CFTRANS\t\t\t=\t"); Out(CFTRANS); Out("\n");
	Out("\t- Heat transfer, HGTRANS\t\t\t=\t"); Out(HGTRANS); Out("\n");
	Out("\t- Bend pressure loss, CPTRANS\t\t\t=\t"); Out(CPTRANS); Out("\n");
	if(freq==1) Out("\tDefault sampling rate, freq\t\t\t=\tonce per timestep\n");
	else {Out("\tDefault sampling rate, freq\t\t\t=\tonce per "); Out(freq); Out(" timesteps\n");}
	Out("\n");

	// Reference values
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Reference values", "-", "\t"));
	Out("\tPressure, pref\t\t\t\t\t=\t"); Out(PREF); Out(" bar\n");
	Out("\tLength, xref\t\t\t\t\t=\t"); Out(xref); Out(" m\n");	
	Out("\tArea, fref\t\t\t\t\t=\t"); Out(fref); Out(" m^2\n");
	Out("\n");

	// Constants
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Constants", "-", "\t"));
	Out("\t(3 - k)/(2(k - 1)), AA\t\t\t\t=\t"); Out(AA); Out("\n");
	Out("\t(k + 1)/(2(k - 1)), BB\t\t\t\t=\t"); Out(BB); Out("\n");
	Out("\t(k - 1)/(2k), Q\t\t\t\t\t=\t"); Out(Q); Out("\n");
	Out("\t(2k)/(k - 1), QI\t\t\t\t=\t"); Out(QI); Out("\n");
	Out("\n");

	// Boundary conditions
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Boundary condition parameters", "-", "\t"));
	Out("\tSpecial integer identifier, PARENT\t\t=\t"); Out(PARENT); Out("\n");
	Out("\tConfiguration string length, STRLEN\t\t=\t"); Out(STRLEN); Out("\n");
	if(SHOW_conn) {
		Out("\tShowing pipe b.c. connections, SHOW_conn\t=\t"); Out(TrueOrFalse(SHOW_conn)); Out("\n");
	}
	else {
		Out("\tNot showing pipe b.c. connections, SHOW_conn\t=\t"); Out(TrueOrFalse(SHOW_conn)); Out("\n");
	}
	if(SHOW_config) {
		Out("\tShowing b.c. config. interp., SHOW_config\t=\t"); Out(TrueOrFalse(SHOW_config)); Out("\n");
	}
	else {
		Out("\tNot showing b.c. config. interp., SHOW_config\t=\t"); Out(TrueOrFalse(SHOW_config)); Out("\n");
	}
	if(SHOW_params) {
		Out("\tShowing object parameters, SHOW_params\t\t=\t"); Out(TrueOrFalse(SHOW_params)); Out("\n");
	}
	else {
		Out("\tNot showing object parameters, SHOW_params\t=\t"); Out(TrueOrFalse(SHOW_params)); Out("\n");	
	}
	Out("\tTolerance for effective zero, ZERO_TOL\t\t=\t"); Out(ZERO_TOL); Out("\n");
	Out("\n");

	// Screen output
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Screen output", "-", "\t"));
	Out("\tInterval between screen outputs, outputInterval\t=\t"); Out(outputInterval); Out(" s\n");
	if(SHOW_x) {
		Out("\tShowing pipe x (m) values, SHOW_x\t\t=\t"); Out(TrueOrFalse(SHOW_x)); Out("\n");
	}
	else {
		Out("\tNot showing pipe x (m) values, SHOW_x\t\t=\t"); Out(TrueOrFalse(SHOW_x)); Out("\n");
	}
	if(SHOW_X) {
		Out("\tShowing pipe X values, SHOW_X\t\t\t=\t"); Out(TrueOrFalse(SHOW_X)); Out("\n");
	}
	else {
		Out("\tNot showing pipe X values, SHOW_X\t\t=\t"); Out(TrueOrFalse(SHOW_X)); Out("\n");
	}
	if(SHOW_DELX) {
		Out("\tShowing pipe DELX values, SHOW_DELX\t\t=\t"); Out(TrueOrFalse(SHOW_DELX)); Out("\n");
	}
	else {
		Out("\tNot showing pipe DELX values, SHOW_DELX\t\t=\t"); Out(TrueOrFalse(SHOW_DELX)); Out("\n");
	}
	if(SHOW_d) {
		Out("\tShowing pipe d (m) values, SHOW_d\t\t=\t"); Out(TrueOrFalse(SHOW_d)); Out("\n");
	}
	else {
		Out("\tNot showing pipe d (m) values, SHOW_d\t\t=\t"); Out(TrueOrFalse(SHOW_d)); Out("\n");
	}
	if(SHOW_dddx) {
		Out("\tShowing pipe dddx (m) values, SHOW_dddx\t\t=\t"); Out(TrueOrFalse(SHOW_dddx)); Out("\n");
	}
	else {
		Out("\tNot showing pipe dddx (m) values, SHOW_dddx\t=\t"); Out(TrueOrFalse(SHOW_dddx)); Out("\n");
	}
	if(SHOW_d2ddx2)	{
		Out("\tShowing pipe d2ddx2 (m) values, SHOW_d2ddx2\t=\t"); Out(TrueOrFalse(SHOW_d2ddx2)); Out("\n");
	}
	else {
		Out("\tNot showing pipe d2ddx2 (m) values, SHOW_d2ddx2\t=\t"); Out(TrueOrFalse(SHOW_d2ddx2)); Out("\n");
	}
	if(SHOW_f) {
		Out("\tShowing pipe f (m^2) values, SHOW_f\t\t=\t"); Out(TrueOrFalse(SHOW_f)); Out("\n");
	}
	else {
		Out("\tNot showing pipe f (m^2) values, SHOW_f\t\t=\t"); Out(TrueOrFalse(SHOW_f)); Out("\n");
	}
	if(SHOW_cfa) {
		Out("\tShowing pipe cfa (m^2) values, SHOW_cfa\t\t=\t"); Out(TrueOrFalse(SHOW_cfa)); Out("\n");
	}
	else {
		Out("\tNot showing pipe cfa (m^2) values, SHOW_cfa\t=\t"); Out(TrueOrFalse(SHOW_cfa)); Out("\n");
	}
	if(SHOW_cfa_delx) {
		Out("\tShowing pipe cfa_delx (m) values, SHOW_cfa_delx\t=\t"); Out(TrueOrFalse(SHOW_cfa_delx)); Out("\n");
	}
	else {
		Out("\tNot showing cfa_delx (m) values, SHOW_cfa_delx\t=\t"); Out(TrueOrFalse(SHOW_cfa_delx)); Out("\n");
	}
	if(SHOW_dfdx) {
		Out("\tShowing pipe dfdx (m) values, SHOW_dfdx\t\t=\t"); Out(TrueOrFalse(SHOW_dfdx)); Out("\n");
	}
	else {
		Out("\tNot showing pipe dfdx (m) values, SHOW_dfdx\t=\t"); Out(TrueOrFalse(SHOW_dfdx)); Out("\n");
	}
	if(SHOW_vol) {
		Out("\tShowing node volume (m^3) values, SHOW_vol\t=\t"); Out(TrueOrFalse(SHOW_vol)); Out("\n");
	}
	else {
		Out("\tNot showing node volume (m^3) values, SHOW_vol\t=\t"); Out(TrueOrFalse(SHOW_vol)); Out("\n");
	}
	if(SHOW_W) {
		Out("\tShowing pipe solution vector W[], SHOW_W\t=\t"); Out(TrueOrFalse(SHOW_W)); Out("\n");
	}
	else {
		Out("\tNot showing pipe solution vector W[], SHOW_W\t=\t"); Out(TrueOrFalse(SHOW_W)); Out("\n");
	}
	if(SHOW_F) {
		Out("\tShowing pipe flux vector F[], SHOW_F\t\t=\t"); Out(TrueOrFalse(SHOW_F)); Out("\n");
	}
	else {
		Out("\tNot showing pipe flux vector F[], SHOW_F\t=\t"); Out(TrueOrFalse(SHOW_F)); Out("\n");
	}
	if(SHOW_C) {
		Out("\tShowing pipe source vector C[], SHOW_C\t\t=\t"); Out(TrueOrFalse(SHOW_C)); Out("\n");
	}
	if(SHOW_W_pred)	{
		Out("\tShowing pred sol vector W_pred[], SHOW_W_pred\t=\t"); Out(TrueOrFalse(SHOW_W_pred)); Out("\n");
	}
	else {
		Out("\tNot showing W_pred[] vector, SHOW_W_pred\t=\t"); Out(TrueOrFalse(SHOW_W_pred)); Out("\n");
	}
	if(SHOW_F_pred)	{
		Out("\tShowing pred flux vector F_pred[], SHOW_F_pred\t=\t"); Out(TrueOrFalse(SHOW_F_pred)); Out("\n");
	}
	else {
		Out("\tNot showing F_pred[] vector, SHOW_F_pred\t=\t"); Out(TrueOrFalse(SHOW_F_pred)); Out("\n");
	}
	if(SHOW_C_pred) {
		Out("\tShowing pred source vec. C_pred[], SHOW_C_pred\t=\t"); Out(TrueOrFalse(SHOW_C_pred)); Out("\n");
	}
	else {
		Out("\tNot showing C_pred[] vector, SHOW_C_pred\t=\t"); Out(TrueOrFalse(SHOW_C_pred)); Out("\n");
	}
	if(SHOW_pathlines) {
		Out("\tShowing pathlines, SHOW_pathlines\t\t=\t"); Out(TrueOrFalse(SHOW_pathlines)); Out("\n");
	}
	else {
		Out("\tNot showing pathlines, SHOW_pathlines\t\t=\t"); Out(TrueOrFalse(SHOW_pathlines)); Out("\n");
	}
	Out("\n");

	// File output
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("File output", "-", "\t"));
	Out("\tFile extension for results files, strFileExt\t=\t"); Out(strFileExt); Out("\n");
	Out("\n");

	// Engines, cylinders & valves
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Engine parameters", "-", "\t"));
	Out("\tCylinder data file directory, CYL_DIR\t\t=\t"); Out(CYL_DIR); Out("\n");
	
	if(!HOMENTROPIC) {
		Out("\n");
		Out(Underline("Valve parameters", "-", "\t"));
		Out("\tValve timing files directory, VT_DIR\t\t=\t"); Out(VT_DIR); Out("\n");
	}
	else {
		Out("\n");
		Out(Underline("Valve parameters", "-", "\t"));
		Out("\tHom. tol. 1, POPPET_H_TOL1\t\t\t=\t"); Out(POPPET_H_TOL1); Out("\n");
		Out("\tHom. tol. 2, POPPET_H_TOL2\t\t\t=\t"); Out(POPPET_H_TOL2); Out("\n");
		Out("\tValve timing files directory, VT_DIR\t\t=\t"); Out(VT_DIR); Out("\n");
	}
	Out("\n");

	// End environments
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("End environment parameters", "-", "\t"));
	Out("\tEnd environment file directory, ENDENV_DIR\t=\t"); Out(ENDENV_DIR); Out("\n");
	if(!HOMENTROPIC) {
		Out("\tLoop tolerance, NHI_TOL\t\t\t\t=\t"); Out(NHI_TOL); Out("\n");
	}
	Out("\tLoop tolerance, NOZZ_TOL\t\t\t=\t"); Out(NOZZ_TOL); Out("\n");
	if(USE_PHI) Out("\tArea ratio does matter under inflow, USE_PHI\t=\t");
	else Out("\tArea ratio doesn't matter under inflow, USE_PHI\t=\t");
	Out(TrueOrFalse(USE_PHI)); Out("\n");
	Out("\n");

	// Junctions
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Junctions", "-", "\t"));
	Out("\tJunction file directory, JUNC_DIR\t\t=\t"); Out(JUNC_DIR); Out("\n");
	Out("\tMax. branches per junction, max_branches\t=\t"); Out(max_branches); Out("\n");
	Out("\n");

	// Adiabatic pressure loss devices
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Adiabatic pressure loss devices", "-", "\t"));
	Out("\tAPLD loss file directory, APLD_DIR\t\t=\t"); Out(APLD_DIR); Out("\n");
	Out("\tLoss loop tolerance, LOSS_TOL\t\t\t=\t"); Out(LOSS_TOL); Out("\n");
	Out("\n");

	// Sudden area changes
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Sudden area changes", "-", "\t"));
	Out("\tTolerance on N convergence, tolSuddenEnlrgN\t=\t"); Out(tolSuddenEnlrgN); Out("\n");
	Out("\tTolerance on lambda_in, tolSuddenEnlrgLin\t=\t"); Out(tolSuddenEnlrgLin); Out("\n");
	Out("\n");

	// Turbines
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Turbines", "-", "\t"));
	Out("\tTurbine map directory, TURB_DIR\t\t\t=\t"); Out(TURB_DIR); Out("\n");
	Out("\n");
	
	// Transmissive boundaries
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Transmissive boundaries", "-", "\t"));
	Out("\tTransmissive files directory, TRANSM_DIR\t=\t"); Out(TRANSM_DIR); Out("\n");
	Out("\n");

	// Assemblies
	// ----------------------------------------------------------------------------------------------------
	Out(Underline("Assemblies", "-", "\t"));
	Out("\tNo. of assemblies, NASSEMBLY\t\t\t=\t"); Out(NASSEMBLY); Out("\n");
	Out("\n");
}

void CProperties::Configure(char *InputFile, char** &rCONFP, double* &rENDCORR, int** &rASSYS, int** &rPIPES, 
								 int** &rPIPES_ENDS, int* &rNUM_PIPES, int number_of, int strlen, char* bc_type, int npipes, 
								 vector<string> &rConfStrs, int assyid)
{
	if(SHOW_calls){Out("Ppt.Configure ("); Out(InputFile); Out(")\n");}
	int no_identified = 0;
	ReadConfigFile(InputFile, rCONFP, rENDCORR, number_of, strlen, bc_type, no_identified, npipes);
	ConvertStrings(InputFile, rCONFP, rASSYS, rPIPES, rPIPES_ENDS, rNUM_PIPES, number_of, bc_type, no_identified, npipes, rConfStrs, assyid);
}

// Reads configuration files for closed, open and nozzle boundaries
void CProperties::ReadConfigFile(char *InputFile, char** &rCONFP, double* &rENDCORR, int number_of, int strlen_per_pipe, char* bc_type,	int &rno_identified, int npipes)
{
	FILE *stream;
	stream = fopen(InputFile, "r");

	rCONFP = new char*[number_of];
	rENDCORR = new double[number_of];
	int q;
	int strlen = strlen_per_pipe*npipes; // String length per pipe * number of pipes joined at this boundary

//cout << InputFile << ": " << "number_of = " << number_of << endl;
//cout << "strlen_per_pipe = " << strlen_per_pipe << endl;
//cout << "npipes = " << npipes << endl;
//cout << "strlen = " << strlen << endl;
//cout << endl;

	if(SHOW_config)
	{
		Out("\n");
		Out(Underline(InputFile, "-", "\t"));
	}

	if(stream == NULL)
	{
		Out("Properties:ReadConfigFile: Error opening input file ");
		Out(InputFile);
		Out(". Exiting program...\n");
		exit(1);
	}
	else
	{
		// Set pointer to beginning of file
		fseek(stream, 0L, SEEK_SET);
		q=0;
		do
		{
			// Must create a new string holder each loop, otherwise rJCONFP[q]
			// will always point to the last string read
			char* s;
			s = new char [strlen];
			float f;

			if(q>=number_of)
			{
				if(SHOW_config)
				{
					Out("\tWarning: more "); Out(bc_type); Out(" configurations exist in "); Out(InputFile); Out("\n");
					Out("\tthan have been specified in the input file.\n");
				}
				break;
			}	

			fscanf(stream, "%s", s);	// This sees the label e.g. "CE0"
			if(SHOW_config){Out("\t");Out(s);Out("\t");}
			fscanf(stream, "%s", s);	// This sees the configuration e.g. "0e1o"
			if(SHOW_config){Out("\t");Out(s);Out("\t");}
			fscanf(stream, "%f", &f);	// This sees the end correction parameter e.g. "0.355"
			if(SHOW_config){Out("\t");Out(f);Out("\t");Out("\n");}
			rCONFP[q] = s;
			rENDCORR[q] = f;
			++q;
		}while(fscanf(stream, "%l")!=EOF);
		rno_identified = q;
		fclose(stream);
	}

	if(SHOW_config)
	{
		for(q=0; q<number_of; ++q) 
		{
			for(int r=0; r<strlen; ++r)
			{
				Out("\t"); Out(bc_type); Out(" rCONFP["); Out(q); Out("]["); Out(r); Out("] = "); Out(rCONFP[q][r]); Out("\n");
			}
		}
	}
}

void CProperties::ConvertStrings(char *InputFile, char** &rCONFP, int** &rASSYS, int** &rPIPES, 
								 int** &rPIPES_ENDS, int* &rNUM_PIPES, int number_of, char* bc_type, int no_identified, int npipes, 
								 vector<string> &rConfStrs, int assyid)
{
	// Convert configuration strings into usuable integers
	int q, r;

	// number_of is total number of this type of boundary
	// npipes is total number of connections at boundary[q out of number_of]
	
	// Dimension the real arrays
	rASSYS = new int* [number_of];
	rPIPES = new int* [number_of];		
	rPIPES_ENDS = new int* [number_of];
	rNUM_PIPES = new int [number_of];	// Lists the number of pipes connected at each boundary

	// The assembly numbers for each connection (temporary):
	int** rASSYStemp;
	rASSYStemp = new int* [number_of];
	for(q=0; q<number_of; ++q) rASSYStemp[q] = new int [npipes];

	// The pipe numbers for each connection (temporary):
	int** rPIPEStemp;
	rPIPEStemp = new int* [number_of];
	for(q=0; q<number_of; ++q) rPIPEStemp[q] = new int [npipes]; 

	// The end of each pipe to use for each connection (temporary):
	int** rPIPES_ENDStemp;
	rPIPES_ENDStemp = new int* [number_of];
	for(q=0; q<number_of; ++q) rPIPES_ENDStemp[q] = new int [npipes];

	int pipe_counter = 0;
	int assy_counter = 0;
	bool units, tens, hundreds;
	units=false; tens=false; hundreds=false;

	int digit, temp_number;
	bool found_int;
	bool stop;
	bool PIPE;

	for(q=0; q<number_of; ++q)
	{
		r=0;
		pipe_counter=0;
		assy_counter=0;
		stop=false;
		std::string temp;
		PIPE = true; // Is the b.c. for a pipe (default) or a loop?
		do
		{
			digit = rCONFP[q][r];
//			cout << "digit = " << digit << endl;
			if(digit>=48 && digit<58) found_int=true;
			else found_int=false;

			if(found_int)
			{
				if(tens) hundreds=true;
				else
				{
					if(units) tens=true;
					else units=true;
				}
			}
			else
			{
//				cout << "rCONFP[q][r] = '" << rCONFP[q][r] << "'" << endl;
				if(digit == 0/*rCONFP[q][r] == ' '*/)
				{ 
					// If we come to a space, the string length is
					// longer than absolutely necessary for this particular configuration
					stop = true;
//					cout << "stop = true" << endl;
				}
				else
				{
					// Then we have found a number, last digit at rCONFP[q][r-1]
					// NB: Integers start at ASCII character no. 48 
					if(hundreds)
						temp_number = 100*(rCONFP[q][r-3] - 48) 
												+ 10*(rCONFP[q][r-2] - 48) 
												+ (rCONFP[q][r-1] - 48);
					/*
						rPIPEStemp[q][pipe_counter] = 100*(rCONFP[q][r-3] - 48) 
												+ 10*(rCONFP[q][r-2] - 48) 
												+ (rCONFP[q][r-1] - 48);
					*/
					else
					{
						if(tens)
							temp_number = 10*(rCONFP[q][r-2] - 48) 
													+ (rCONFP[q][r-1] - 48);
						/*
							rPIPEStemp[q][pipe_counter] = 10*(rCONFP[q][r-2] - 48) 
													+ (rCONFP[q][r-1] - 48);
						*/
						else
							if(units)
							{
							//	cout << "rCONFP[q][r-1] = " << rCONFP[q][r-1] << endl;
								temp_number = (rCONFP[q][r-1] - 48);
								//rPIPEStemp[q][pipe_counter] = (rCONFP[q][r-1] - 48);
							}
					}

					// If temp_number is special identifier, then boundary is located on parent assembly
					if(temp_number==PARENT) temp_number=assyid;
	
					// Record which end of this pipe should be used
					if(rCONFP[q][r]=='A' || rCONFP[q][r]=='a')
					{
						rASSYStemp[q][assy_counter] = temp_number;
						++assy_counter; // Move on to the next assembly
					}
					else
					{
						if(rCONFP[q][r]=='E' || rCONFP[q][r]=='e')
						{
							rPIPEStemp[q][pipe_counter] = temp_number;
							rPIPES_ENDStemp[q][pipe_counter] = EVEN;
							++pipe_counter; // Move on to the next pipe making up the b.c.
						}
						else
						{	
							if(rCONFP[q][r]=='O' || rCONFP[q][r]=='o')
							{
								rPIPEStemp[q][pipe_counter] = temp_number;
								rPIPES_ENDStemp[q][pipe_counter] = ODD;
								++pipe_counter; // Move on to the next pipe making up the b.c.
							}
							else
							{				
								cout << "\t" << rCONFP[q][r] 
									<< " is not a valid pipe end specifier for " << bc_type << " ["
									<< q << "], pipe [" << pipe_counter << "]" << endl;
								cout << "\t" << "Please check the configuration file " 
									<< InputFile << endl;
								exit(1);
							}
						}
					}

					// Reset for next number
					units=false; tens=false; hundreds=false;

					// Add part string to config string
					temp += rCONFP[q][r];
					temp += IntToString(temp_number);	
				}
			}
			++r;
		}while(pipe_counter<npipes && !stop);

		// Add configuration string to list of configuration strings
		rConfStrs.push_back(temp);
		if(SHOW_config) std::cout << std::endl << "\tFinal configuration string: " << temp << std::endl << std::endl;

		// Now we know how many pipe connections per bc there are here
//		cout << "pipe_counter = " << pipe_counter << endl;
//		char pause;
//		cin >> pause;
		rASSYS[q] = new int [assy_counter];
		rPIPES[q] = new int [pipe_counter];
		rNUM_PIPES[q] = pipe_counter; // Stores the number of pipe connections here
		rPIPES_ENDS[q] = new int [pipe_counter]; // The end of each pipe to use for each junction:
	
		// Copy the first pipe_counter values from the temp arrays to the real arrays
		for(r=0; r<pipe_counter; ++r)
		{
			rASSYS[q][r] = rASSYStemp[q][r];
			rPIPES[q][r] = rPIPEStemp[q][r];	// Remember the [0] stores the npipes
			rPIPES_ENDS[q][r] = rPIPES_ENDStemp[q][r];
		}

		if(SHOW_config)
		{
			cout << endl;
			cout << Underline(bc_type, "-", q, "\t");
			if(rNUM_PIPES[q]>1) cout << "\t" << rNUM_PIPES[q] << " pipes join at this boundary" << endl;
			else
			{
				if(rNUM_PIPES[q]==1) cout << "\t" << rNUM_PIPES[q] << " pipe joins at this boundary" << endl;
				else cout << "\tNo pipes join at this boundary" << endl;
			}

			for(r=0; r<pipe_counter; ++r)
			{
				// Header
				if(bc_type=="cylinder" || bc_type=="CYLINDER" || bc_type=="Cylinder")
				{
					cout << endl;
					cout << Underline("Pipe connection", "-", r, "\t\t");
				}
				else 
				{
					if(bc_type=="EXHAUST SUDDEN AREA CHANGE" || bc_type=="INTAKE SUDDEN AREA CHANGE"
						|| bc_type=="EXHAUST ADIABATIC PRESSURE LOSS DEVICE" || bc_type=="INTAKE ADIABATIC PRESSURE LOSS DEVICE")
					{
						cout << endl;
						if(r==0) cout << Underline("Left side pipe connection", "-", "\t\t");
						else cout << Underline("Right side pipe connection", "-", "\t\t");
					}
					else // Junction
					{
						cout << endl; 
						cout << Underline("Pipe connection", "-", r, "\t\t");
					}
				}
			// Data
/*
			if((bc_type=="cylinder" || bc_type=="CYLINDER" || bc_type=="Cylinder")
				&& r==1 && !IPIPE) // Then no intake pipe intake on cylinders
			{
				cout << "  No intake pipe specified.\n";
				if(!IAIR) cout << "  No intake valve specified.\n";
				else cout << "  Using air receiver as inlet to intake side of cylinder.\n";
			}
			else
*/
				{
//					cout << "\t\trASSYS[" << q << "][" << r << "] = " << rASSYStemp[q][r] << endl;
//					cout << "\t\trPIPES[" << q << "][" << r << "] = " << rPIPES[q][r] << endl;
					if(rPIPES_ENDS[q][r]==ODD)
						cout << "\t\trPIPES_ENDS[" << q << "][" << r << "] = ODD\n";
					else
					{
						if(rPIPES_ENDS[q][r]==EVEN)
						{
							cout << "\t\trPIPES_ENDS[" << q << "][" << r << "] = EVEN\n";
						}
						else
						{
							cout << "\t\trPIPES_ENDS[" << q << "][" << r << "] = UNKNOWN\n";
						}
					}
				}
			}
		}
	}

	if(SHOW_config)
	{
		cout << endl;
		if(no_identified==1)
			cout << "\t" << no_identified << " " << bc_type << " boundary has been included in this assembly.\n";
		else
		{
			if(no_identified==0)
				cout << "\t" << "No " << bc_type << " boundaries have been included in this assembly.\n";
			else
				cout << "\t" << no_identified << " " << bc_type << " boundaries have been included in this assembly.\n";
		}
	}
}

void CProperties::SetupFiles()
{

/*
	EX_P_FILE_M = new FILE** [NEXPIPES];	// One set of files for each exhaust pipe
	for(f=0; f<NEXPIPES; ++f)
		EX_P_FILE_M[f] = new FILE* [NMEASUREMENTS];

	char*** EX_P_FileNames_M;
	CreateFileNames2(EX_P_FileNames_M, NEXPIPES, NMEASUREMENTS, "ex_");
	for(f=0; f<NEXPIPES; ++f)
	{
		for(n=0; n<NMEASUREMENTS; ++n)
		{
			EX_P_FILE_M[f][n] = fopen(EX_P_FileNames_M[f][n], "w");
		}
	}

	IN_P_FILE_M = new FILE** [NINPIPES];	// One set of files for each intake pipe
	for(f=0; f<NINPIPES; ++f)
		IN_P_FILE_M[f] = new FILE* [NMEASUREMENTS];

	char*** IN_P_FileNames_M;
	CreateFileNames2(IN_P_FileNames_M, NINPIPES, NMEASUREMENTS, "in_");
	for(f=0; f<NINPIPES; ++f)
	{
		for(n=0; n<NMEASUREMENTS; ++n)
		{
			IN_P_FILE_M[f][n] = fopen(IN_P_FileNames_M[f][n], "w");
		}
	}
*/

/*
	EX_J_FILE = new FILE* [NEXJUNCS];		// One for each exhaust junction
	char** EX_J_FileNames;
	CreateFileNames(EX_J_FileNames, NEXJUNCS, "junction");
	for(f=0; f<NEXJUNCS; ++f) EX_J_FILE[f] = fopen(EX_J_FileNames[f], "w");
*/
//	EX_TRANSM_FILE = new FILE* [NEXTRANSM];		// One for each transmissive boundary
//	char** EX_TRANS_FileNames;
//	CreateFileNames(EX_TRANS_FileNames, NEXTRANSM, "transmE"/* at least 7 letters */);
//	for(f=0; f<NEXTRANSM; ++f) EX_TRANSM_FILE[f] = fopen(EX_TRANS_FileNames[f], "w");

//	valvefileptr = fopen("valve_loop.txt", "w");
//	sudfileptr = fopen("sudden.txt", "w");
}

void CProperties::CreateFileNames(char** &rFileName, int num_objects, char* object_type)
{
	rFileName = new char* [num_objects];
	for(int f=0; f<num_objects; ++f)
	{
		rFileName[f] = new char [13];
		rFileName[f][0] = object_type[0]; rFileName[f][1] = object_type[1]; rFileName[f][2] = object_type[2];
		rFileName[f][3] = object_type[3]; rFileName[f][4] = object_type[4]; rFileName[f][5] = object_type[5]; 
		rFileName[f][6] = object_type[6];
		rFileName[f][7] = f + 48;
		rFileName[f][8] = '.'; rFileName[f][9] = 't'; rFileName[f][10] = 'x'; rFileName[f][11] = 't';
		rFileName[f][12] = '\0';
	}
}

void CProperties::CreateFileNames2(char*** &rFileName, int num_objects, int num_locs, char* object_type)
{
	int f, n;
	rFileName = new char** [num_objects];
	for(f=0; f<num_objects; ++f)
	{
		rFileName[f] = new char* [num_locs];
		for(n=0; n<num_locs; ++n)
		{
			rFileName[f][n] = new char [13];
			rFileName[f][n][0] = object_type[0]; rFileName[f][n][1] = object_type[1]; rFileName[f][n][2] = object_type[2];
			rFileName[f][n][3] = 'p';
			rFileName[f][n][4] = f + 48;
			rFileName[f][n][5] = '_';
			rFileName[f][n][6] = 'l';
			rFileName[f][n][7] = n + 48;
			rFileName[f][n][8] = '.'; rFileName[f][n][9] = 't'; rFileName[f][n][10] = 'x'; rFileName[f][n][11] = 't';
			rFileName[f][n][12] = '\0';
		}
	}
}

void CProperties::ReadLossCoefficients(char *InputFile)
{
    if(SHOW_calls){Out("Ppt.ReadLossCoefficients("); Out(InputFile); Out(")\n");}
	
	float C, C_old, W, W_old, M, L;
	int graph, num_graphs, curve, point;
	int *num_curves;
	int **num_points;

	int max_graphs = 12;
	int max_curves_per_graph = 12;
	num_curves = new int [max_graphs];
	num_points = new int* [max_graphs];

	for(graph=0; graph<max_graphs; ++graph)
		num_points[graph] = new int [max_curves_per_graph];

	FILE *stream;
	stream = fopen(InputFile, "r");
		
	if(stream == NULL)
	{
		printf("Properties::ReadLossCoefficients: Error opening input file...\n");
		exit(1);
	}
	else
	{
		// Set pointer to beginning of file
		fseek(stream, 0L, SEEK_SET);
		// Go through the file to count points per curve, etc.
		graph=0;
		curve=0;
		point=0;
		C=1;
		W=0;
		do
		{
			C_old = C;
			W_old = W;
			fscanf(stream, "%f", &C );
			fscanf(stream, "%f", &W );
			fscanf(stream, "%f", &M );
			fscanf(stream, "%f", &L );
		
			if((fabs(C_old - C)>1e-4))
			{
				// Changing graphs
				num_points[graph][curve] = point;
				//cout << num_points[graph][curve] << " points on curve " << curve << endl;
				++curve;
				num_curves[graph] = curve;
				//cout << num_curves[graph] << " curves on graph " << graph << endl;
				point=0;
				curve=0;
				++graph;
			}
			else
			{
				// Changing curves
				if((fabs(W_old - W)>1e-4))
				{ 
					num_points[graph][curve] = point;
					//num_points = point;
					//cout << num_points[graph][curve] << " points on curve " << curve << endl;
					point=0;
					++curve;					
				}
			}
			//cout << "graph = " << graph << endl;
			//cout << "curve = " << curve << endl;
			//cout << "L" << C << ": Curve W = " << W << ", point " << point << ": M = " << M << ", L = " << L << endl;
			++point;			
		}while(fscanf(stream, "%l")!=EOF);
		num_points[graph][curve] = point;
		//cout << num_points[graph][curve] << " points on curve " << curve << endl;
		++curve;
		num_curves[graph] = curve;
		//cout << num_curves[graph] << " curves on graph " << graph << endl;
		++graph;
		num_graphs = graph;
		//cout << "num_graphs = " << num_graphs << endl;	
		/*
		for(graph=0; graph<num_graphs; ++graph)
		{
			cout << "num_curves[" << graph << "] = " << num_curves[graph] << endl;
			for(curve=0; curve<num_curves[graph]; ++curve)
				cout << "num_points[" << graph << "][" << curve << "] = " << num_points[graph][curve] << endl;
		}
		*/

		// Can now dimension the matrix
		Loss = new double*** [num_graphs];
		for(graph=0; graph<num_graphs; ++graph) // For each graph
		{
			Loss[graph] = new double** [num_curves[graph]];
			for(curve=0; curve<num_curves[graph]; ++curve) // For each curve
			{
				Loss[graph][curve] = new double* [num_points[graph][curve]];
				for(point=0; point<num_points[graph][curve]; ++point)
					Loss[graph][curve][point] = new double [4]; // 4 pieces of data per point (C, W, Ma, L)
			}
		}

		// Now go through and copy the data to the relevant place
		// Set pointer to beginning of file
		fseek(stream, 0L, SEEK_SET);

		graph=0;
		curve=0;
		point=0;
		C=1;
		W=0;
		do
		{
			C_old = C;
			W_old = W;
			fscanf(stream, "%f", &C );
			fscanf(stream, "%f", &W );
			fscanf(stream, "%f", &M );
			fscanf(stream, "%f", &L );
		
			if((fabs(C_old - C)>1e-4))
			{ 
				point=0;
				curve=0;
				++graph;
			}
			else
			{
				if((fabs(W_old - W)>1e-4))
				{ 
					point=0;
					++curve;
				}
			}
			//cout << "curve = " << curve << endl;
			//cout << "L" << C << ": Curve W = " << W << ", point " << point << ": M = " << M << ", L = " << L << endl;
			Loss[graph][curve][point][0] = C;
			Loss[graph][curve][point][1] = W;
			Loss[graph][curve][point][2] = M;
			Loss[graph][curve][point][3] = L;
			++point;
		}while(fscanf(stream, "%l")!=EOF);
			
		fclose(stream);
	}

/*
	cout << setprecision(6);
	for(graph=0; graph<num_graphs; ++graph) // For each graph
	{
		for(curve=0; curve<num_curves[graph]; ++curve) // For each curve
		{
			for(point=0; point<num_points[graph][curve]; ++point) // For each point
			{
				//cout << "Graph " << graph << ", curve " << curve << ", point " << point << ":\n";
				cout << Loss[graph][curve][point][0] << "\t" 
						<< Loss[graph][curve][point][1] << "\t" 
							<< Loss[graph][curve][point][2] << "\t" 
								<< Loss[graph][curve][point][3] << "\t" << endl; 
			}
			cout << endl;
		}
		cout << endl;
	}
*/
}


double CProperties::ViscosityAir(double T)
// ====================================================	//
// Gas viscosity mu as a function of temperature		//
// using Sutherland's or Blair's equation				//
// ====================================================	//
{
    if(SHOW_calls){Out("Ppt.ViscosityAir("); Out(T); Out(")\n");}

	if(VISCOSITY_AIR==SUTHERLAND_VISCOSITY)
	{
		// Constants for air
		double C1 = 1.458e-6;	// (kg.m^-1.s^-1.K^-0.5)
		double C2 = 110.4;		// (K)
		return (C1*pow(T, 1.5))/(T + C2);
	}
	else
	{
		if(VISCOSITY_AIR==BLAIR_VISCOSITY)	
		// Gas viscosity mu as a function of temperature using Blair's curve fitted equation	
			return 7.457e-6 + 4.1547e-8*T - 7.4793e-12*pow(T,2);
		else
		{
			if(VISCOSITY_AIR==RICARDO_VISCOSITY) return 0.33e-6*pow(T,0.7);
			else
			{
				if(VISCOSITY_AIR==CONSTANT_VISCOSITY) return mu_air_const;
				else
				{
					cout << "VISCOSITY_AIR not specified correctly!\n";
					exit(1);
				}
			}
		}
	}
}
double CProperties::cpAir(double T)
// ====================================================	//
// Specific heat of air at constant pressure as a		//
// function of temperature in Kelvin using				//
// Zucrow & Hoffman equation.							//
// ====================================================	//
{
    if(SHOW_calls){Out("Ppt.cpAir("); Out(T); Out(")\n");}

    if(constProps) return Cp_air;
	else
	{
		if(T<1000) return (3.65359 - 1.33736e-3*T + 3.29421e-6*pow(T,2) - 1.91142e-9*pow(T,3) + 0.275462e-12*pow(T,4))*(R_universal/M_air);
		else return (3.04473 + 1.33805e-3*T - 0.488256e-6*pow(T,2) + 0.0855475e-9*pow(T,3) - 0.00570132e-12*pow(T,4))*(R_universal/M_air);
	}
}
double CProperties::cvAir(double T)
// ====================================================	//
// Specific heat of air at constant volume as a 
// function of cpAir()									//
// ====================================================	//
{
    if(SHOW_calls){Out("Ppt.cvAir("); Out(T); Out(")\n");}

	if(constProps) return Cv_air;
	else return cpAir(T) - R_air;
}

double CProperties::gammaAir(double T)
// ====================================================	//
// Ratio of specific heats as a function of cpAir()		//
// ====================================================	//
{
    if(SHOW_calls){Out("Ppt.gammaAir("); Out(T); Out(")\n");}

	if(constProps) return gamma_air;
	else return cpAir(T)/(cpAir(T) - R_air);
}

double CProperties::ViscosityBloodFunc(double T)
// ====================================================	//
// Blood viscosity 										//
//														//
// ====================================================	//
{
	if (SHOW_calls) { Out("Ppt.ViscosityBlood("); Out(T); Out(")\n"); }

	return ViscosityBlood;
}

double CProperties::BLThickness(double diameter)
// ====================================================	//
// Boundary layer thickness								//
//														//
// ====================================================	//
{
	if (SHOW_calls) { Out("Ppt.BLThickness("); Out(diameter); Out(")\n"); }

	if (BL_RELATIVE) return BL_rel * diameter; // If relative boundary thickness return value based on pipe diameter
	else return BL_abs; // Otherwise returned the fixed absolute value [m]
}
