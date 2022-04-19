// Turbine.cpp: implementation of the CTurbine class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Tools.h"
#include "Turbine.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CTurbine::CTurbine()
{

}

CTurbine::~CTurbine()
{

}

// Configuration
// ====================================================================================================
void CTurbine::Initialise(CProperties* pPpt, CPipe** pPipes, CPipe* &rPipe, int** &rTURBPIPES, int** &rTURBPIPES_ENDS, int* &rTURBPIPES_NUMS, double* &rENDCORR, int id, bool ex, std::string param_dir, int assyid, string parent_assy_res_dir, string calling_object_str)
// ============================================================ //
// Reads the object parameter, opens the results file,			//
// initialises variables and joins the pipe pointers to the		//
// relevant pipes.												//
// This function is called from the main function.				//
// ============================================================ //
{
	if (pPpt->SHOW_calls) { pPpt->Out("CTurbine.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }
	
	std::string bcname_str = "TURBINE";
	ReadInput(pPpt, ConstructString(pPpt, param_dir, bcname_str, id));

	// Enumerate labels
	INLET = 0;							// As opposed to ONE_SIDE; always at [0], always required
	OUTLET = 1;							// Used for cases with VARIABLE outlet pressure
	SP = 1;	MFP = 2; PR = 3; ETA = 4;	// SAE map column labels

	// Parameter definitions and units
	sp_str = "N/T01^0.5";
	mfp_str = "m T01^0.5/p01";	
	pr_str = "p01/p2";				
	eta_str = "eta_TS";
	//sp_units_str = "(r/min)/K^0.5";
	mfp_units_str = "kg/s K^0.5/kPa";
	
	// Join the pipe pointer(s) to the relevant pipe(s)
	InitialiseGen(pPpt, pPipes, rPipe, rTURBPIPES, rTURBPIPES_ENDS, rENDCORR, id, ex, rTURBPIPES_NUMS[id]/*VARIABLE ? 2 : 1*/, assyid, "CTurbine", parent_assy_res_dir);

	// Boundary name
	NAME = new int [rTURBPIPES_NUMS[ID]];
	NAME[0] = TURBINLET;				// Alway required and used

	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Initialise\n");}

	// Set turbine type as constant or variable depending on the number of pipes joined here
	if(rTURBPIPES_NUMS[ID]==2) {
		VARIABLE = true;
		NAME[1] = TURBOUTLET;

		//if(!RADIAL && EQUIVALENT_AREA) {
		if (EQUIVALENT_AREA) {
			EQUIVALENT_AREA = false;
			pPpt->Out("CTurbine::Initialise Turbine ["); pPpt->Out(ID); pPpt->Out("]; the pipe configuration shows a variable outlet pressure turbine - hence must use the single unique curve method. Equivalent areas can only be used for a constant outlet pressure configuration.");
			pPpt->Out("\n");
			pPpt->Out("(paused)");
			pPpt->Out("\n");
			char pause; cin >> pause;
		}
	}
	else {
		if(rTURBPIPES_NUMS[ID]==1) VARIABLE = false;
		else {
			cout << "CTurbine::Initialise Turbine [" << ID << "]; number of joining pipes specified by the configuration string (" 
				<< rTURBPIPES_NUMS[ID] << 
				") inconsistent with available configuration options - only 1 or 2 pipes can join here."  << endl;
			exit(1);
		}
	}

/*
	// Check number of pipes against desired configuration
	if(VARIABLE)
	{
		NAME[1] = TURBOUTLET;
		if(rTURBPIPES_NUMS[ID] != 2)
		{
			cout << "CTurbine::Initialise Turbine [" << ID << "]; number of joining pipes specified by the configuration string (" 
				<< rTURBPIPES_NUMS[ID] << 
				") inconsistent with the desired configuration in the parameter file" 
				<< ConstructString(pPpt, pPpt->param_dir, bcname_str, ID) << " which specifies variable outlet pressure - this requires 2 pipes.\n";
			exit(1);
		}
	}
	else
	{
		// Constant outlet pressure selected in parameter file
		if(rTURBPIPES_NUMS[ID] != 1)
		{
			cout << "CTurbine::Initialise Turbine [" << ID << "]; number of joining pipes specified by the configuration string (" 
				<< rTURBPIPES_NUMS[ID] << 
				") inconsistent with the desired configuration in the parameter file" 
				<< ConstructString(pPpt, pPpt->param_dir, bcname_str, ID) << " which specifies constant outlet pressure - this requires 1 pipe.\n";
			exit(1);
		}
	}
*/

//	pipe_flow = new int [NPIPES];
//	pipe_flow_old = new int [NPIPES];

	// Common initialisations
	// ----------------------
	ON_MAP = true;
	ON_CURVE = true;
	REVERSE_FLOW = false;
	CHOKED_FLOW = false;
	query = 1;

	if(VARIABLE) {
		// Variable pressure initialisations
		// ---------------------------------
		num_extra_reverse = 5;
		num_extra_choked = 10;
		C = new double [2];
	}
	else {
		// Constant pressure initialisations
		// ---------------------------------
		num_extra_reverse = 0;
		num_extra_choked = 0;
	}
	
	// Work and mass flow initialisations
	// ----------------------------------
	if(FIXED_SPEED) NT = fixedSpeed/60; // Convert to s^-1
	else NT = 0;
	dNTdt = 0;

	// Zero running and cycle/total values
	// -----------------------------------
	// Instantaneous
	eT = 0;
	eA = 0;
	eta_TS = 1;
	m_dot = 0;

	eA_over_m_dot = 0;
	W_TI = 0;
	L_TI = 0;

	W_B = 0;
	L_B = 0;
	eta_M = 1;
	
	// Most recent completed cycle
	W_Th = 0;
	W_TA = 0;
	eta_T = 1;
	m_dot_T = 0;

	// Running totals
	work_per_cycle_Th = 0;
	work_per_cycle_A = 0;
	mass_flow_per_cycle = 0;

	// Set up files
	std::string res_str = "res_turbine";
	SetupFiles(pPpt, res_str);
}

void CTurbine::InitialiseMap(CProperties* pPpt)
// ============================================================ //
// Loads turbine map in SAE format								//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".InitialiseMap\n");}

	// Load turbine map - placed here since requires use of pBN->f, for example
	bool SINGLE_CURVE;
	//if(!RADIAL) SINGLE_CURVE = true; else SINGLE_CURVE = false; // If RADIAL or AXIAL
	SINGLE_CURVE = false;

	T01 = pBN[INLET]->T + 0.5*pow(pBN[INLET]->U*pPipe[INLET]->AREF,2)/pPpt->Cp_air;
	T1 = pBN[INLET]->T;
	
	//LoadMap(pPpt, ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TURB_DIR), MAP_FILE), num_curves, num_points, Raw, SINGLE_CURVE);
	LoadMapSAE(pPpt, ConstructString(pPpt, ConstructString(pPpt, pPpt->param_dir, pPpt->TURB_DIR), MAP_FILE)); 

	std::string res_str;
	if (EX) res_str = "processed_map_ex_turbine"; else res_str = "processed_map_in_turbine";
	char* proc_map_file_str = StringToChar(res_str + IntToString(ID) + ".txt"); //"proc_map_turb.txt";
	if(VARIABLE) ProcessMapVar(pPpt, num_extra_reverse, num_extra_choked, proc_map_file_str, true, T1);
	else ProcessMapConst(pPpt, T1);

	// Initialise values in the interpolated data point
	double speed_param;
	//if(INLET_STAG) speed_param = this->NT/sqrt(this->T01);
	//else speed_param = this->NT/sqrt(this->T1);
	speed_param = this->NT/sqrt(this->T01);

	if(VARIABLE) {
		query = ( (*(pCLIN[INLET]))[1]/(pBN[INLET])->AA[1] )
				/ ( (*(pCLIN[OUTLET]))[1]/(pBN[OUTLET])->AA[1] ); // i.e. x = lambda_in_star[INLET]/lambda_in_star[OUTLET];
	}
	else {
		query = (*(pCLIN[ONE_SIDE]))[1]/(pBN[ONE_SIDE])->AA[1]; // i.e. lambda_in_star1
	}

	Interpolate(pPpt, speed_param, query);
}

void CTurbine::LoadMap(CProperties* pPpt, char *InputFile, int* &rnum_curves, int** &rnum_points, double**** &rRaw, bool SINGLE_CURVE)
// ============================================================ //
// Reads turbine map											//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".LoadMap\n");}
	
	float C, C_old, S, S_old, MFP, PR, eta; // C, the type of turbine? i.e. different graphs
	int graph, num_graphs, curve, point;

	int max_graphs = 1;//3;
	int max_curves_per_graph = 10;

	rnum_curves = new int [max_graphs];
	rnum_points = new int* [max_graphs];

	for(graph=0; graph<max_graphs; ++graph)	rnum_points[graph] = new int [max_curves_per_graph];

	FILE *stream;
	stream = fopen(InputFile, "r");
		
	if(stream == NULL)
	{
		pPpt->Out("CTurbine::LoadMap: Turbine ["); pPpt->Out(ID); pPpt->Out("]: Error opening input file...\n");
		exit(1);
	}
	else
	{
		//pPpt->Out("CTurbine::LoadMap: InputFile="); pPpt->Out(InputFile); pPpt->Out("\n");
		// Set pointer to beginning of file
		fseek(stream, 0L, SEEK_SET);
		// Go through the file to count points per curve, etc.
		graph=0;
		curve=0;
		point=0;

		fscanf(stream, "%f", &C );
		//cout << "First C = " << C << endl; // The first value of C
		fscanf(stream, "%f", &S );
		//cout << "First S = " << S << endl; // The first value of S
		fseek(stream, 0L, SEEK_SET);// Reset pointer to beginning of file

		do
		{
			C_old = C;
			S_old = S;
			fscanf(stream, "%f", &C );
			fscanf(stream, "%f", &S );
			fscanf(stream, "%f", &MFP );
			fscanf(stream, "%f", &PR );
			fscanf(stream, "%f", &eta );
		
			if((fabs(C_old - C)>1e-4))
			{
				// Changing graphs
				rnum_points[graph][curve] = point;//cout << rnum_points[graph][curve] << " points on curve " << curve << endl;
				++curve;
				rnum_curves[graph] = curve;//cout << rnum_curves[graph] << " curves on graph " << graph << endl;
				point=0;
				curve=0;
				++graph;
			}
			else
			{
				// Changing curves
//				if(TURBINE_TYPE != AXIAL && (fabs(S_old - S)>1e-4)) // Only one curve on an axial map
				if(!SINGLE_CURVE && (fabs(S_old - S)>1e-4)) // Only one curve on an axial map
				{ 
					rnum_points[graph][curve] = point;//cout << rnum_points[graph][curve] << " points on curve " << curve << endl;
					//num_points = point;
					point=0;
					++curve;					
				}
			}

//			cout << "graph = " << graph << endl;
//			cout << "curve = " << curve << endl;
//			cout << "Type" << C << ": Curve S = " << S << ", point " << point << ": MFP = " << MFP << ", PR = " << PR << endl;
			++point;			
		}while(fscanf(stream, "%l")!=EOF);

		rnum_points[graph][curve] = point;//cout << num_points[graph][curve] << " points on curve " << curve << endl;
		++curve;
		rnum_curves[graph] = curve;//cout << num_curves[graph] << " curves on graph " << graph << endl;
		++graph;
		num_graphs = graph;//cout << "num_graphs = " << num_graphs << endl;	
		
//		for(graph=0; graph<num_graphs; ++graph)
//		{
//			cout << "rnum_curves[" << graph << "] = " << rnum_curves[graph] << endl;
//			for(curve=0; curve<rnum_curves[graph]; ++curve)
//				cout << "rnum_points[" << graph << "][" << curve << "] = " << rnum_points[graph][curve] << endl;
//		}

		// Can now dimension the matrix
		rRaw = new double*** [num_graphs];
		for(graph=0; graph<num_graphs; ++graph) // For each graph
		{
			rRaw[graph] = new double** [rnum_curves[graph]];
			for(curve=0; curve<rnum_curves[graph]; ++curve) // For each curve
			{
				rRaw[graph][curve] = new double* [rnum_points[graph][curve]];
				for(point=0; point<rnum_points[graph][curve]; ++point)
					rRaw[graph][curve][point] = new double [5]; // 5 pieces of data per point (Type, S, MFP, PR, eta)
			}
		}

		// Now go through and copy the data to the relevant place
		// Set pointer to beginning of file
		fseek(stream, 0L, SEEK_SET);

		graph=0;
		curve=0;
		point=0;
		fscanf(stream, "%f", &C );//cout << "First C = " << C << endl;
		fscanf(stream, "%f", &S );//cout << "First S = " << S << endl;

		fseek(stream, 0L, SEEK_SET);// Reset pointer to beginning of file
		do
		{
			C_old = C;
			S_old = S;
			fscanf(stream, "%f", &C );
			fscanf(stream, "%f", &S );
			fscanf(stream, "%f", &MFP );
			fscanf(stream, "%f", &PR );
			fscanf(stream, "%f", &eta );
		
			if((fabs(C_old - C)>1e-4))
			{ 
				point=0;
				curve=0;
				++graph;
			}
			else
			{
//				if(TURBINE_TYPE != AXIAL && (fabs(S_old - S)>1e-4)) // Only one curve on an axial map
				if(!SINGLE_CURVE && (fabs(S_old - S)>1e-4)) // Only one curve on an axial map
				{ 
					point=0;
					++curve;
				}
			}

			//cout << "curve = " << curve << endl;
			//cout << "Type" << C << ": Curve S = " << S << ", point " << point << ": Re = " << Re << ", K = " << K << endl;
			rRaw[graph][curve][point][0] = C;
			rRaw[graph][curve][point][1] = S;
			rRaw[graph][curve][point][2] = MFP;
			rRaw[graph][curve][point][3] = PR;
			rRaw[graph][curve][point][4] = eta;
			++point;
		}while(fscanf(stream, "%l")!=EOF);
		fclose(stream);
	}
/*
	cout << setprecision(6);
	for(graph=0; graph<num_graphs; ++graph) // For each graph
	{
		for(curve=0; curve<rnum_curves[graph]; ++curve) // For each curve
		{
			for(point=0; point<rnum_points[graph][curve]; ++point) // For each point
			{
				//cout << "Graph " << graph << ", curve " << curve << ", point " << point << ":\n";
				cout << rRaw[graph][curve][point][0] << "\t" 
						<< rRaw[graph][curve][point][1] << "\t" 
							<< rRaw[graph][curve][point][2] << "\t" 
								<< rRaw[graph][curve][point][3] << "\t" << endl; 
			}
			cout << endl;
		}
		cout << endl;
	}
//*/
}

void CTurbine::LoadMapSAE(CProperties* pPpt, char *InputFile)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".LoadMapSAE\n");}

	FILE *stream;
	float fp;
	int max = 500;
	char *s;
	s = new char [max];
	
	int row, col;
	int rows, cols;
	int numHeaderRows = 7;

	stream = fopen(InputFile, "r");

	if(stream == NULL)
	{
		pPpt->Out(Identify()); pPpt->Out(".LoadMapSAE:"); pPpt->Out("Error opening input file: "); pPpt->Out(InputFile); pPpt->Out("\n");
		exit(1);
	}
	else
	{
		fseek(stream, 0L, SEEK_SET);		// Go to start of file
		
		fpos_t start_of_string;
		fpos_t temp_pos;
			
		row=0;
		do
		{
			fgetpos(stream, &start_of_string);
			fgets(s, max, stream);			// Read a line of the file at a time
			++row;
		}while(row<numHeaderRows);			// Data in SAE format starts on line 8

		rows=0;
		do
		{
			fgetpos(stream, &start_of_string);
			fgets(s, max, stream);
			++rows;
		}while(fscanf(stream, "%l")!=EOF);

		cols = 4;							// Number of data columns (speed, MFP, PR, eta)
		rawSAE = new double* [rows];
		for(row=0; row<rows; ++row) rawSAE[row] = new double [cols];
		//cout << "rows = " << rows << ". cols = " << cols << ".\n";
		fseek(stream, 0L, SEEK_SET);		// Go back to start of file

		row=0;
		do
		{
			fgetpos(stream, &start_of_string);
			fgets(s, max, stream);			// Read a line of the file at a time
			if(row==0) strNomenclature = s;	// Save turbine nomenclature from first line of SAE file
			++row;
		}while(row<numHeaderRows);			// Data in SAE format starts on line 8

		
		int strptr;
		int number_start;
		int number_pos;	
		bool found_digit;
		bool found_number;
		bool not_a_number; 
		bool NEWLINE = false;
		row=0;
		col=0;
		do
		{
			strptr = 0;
			number_start = 0;
			number_pos = 0;
			found_digit = false;
			found_number = false;
			not_a_number = false;

			fgetpos(stream, &start_of_string);
			fgets(s, max, stream);				
						
			while(s[strptr] != NULL)
			{
				if(( (int(s[strptr]) >=48 && int(s[strptr])<58) 
					|| s[strptr]=='-' || s[strptr]=='.' 
					|| (found_digit && (s[strptr]=='E' || s[strptr]=='e')) )
					&& !not_a_number)
					// A number 0-9 or '-' or decimal point or 'E' or 'e'
				{
					//pPpt->Out("Number\n");
					if(!found_digit)
					{
						if(s[strptr] == 'E' || s[strptr] == 'e') // Cannot start with E or e
						{
							//pPpt->Out("E/e\n");
							not_a_number = true;
							--strptr;	// Roll back
						}
						else
						{
							found_digit = true;
							number_start = strptr;
							number_pos = 0;
						}
					}
					else ++number_pos;	// Continuation of a number

					NEWLINE = false;
				}
				else // Not a number
				{
					//cout << "Not a number.\n";
					if(found_digit)		// Come to end of number
					{
						//cout << "End of a number\n";
						found_digit = false;
						// number runs from s[number_start] to s[number_start+number_pos]
						
						// Now move pointer to before this number
						temp_pos = start_of_string + number_start;
						fsetpos(stream, &temp_pos);

						// Read number into float
						fscanf(stream, "%f", &fp);
						
						rawSAE[row][col] = fp;
						//pPpt->Out("rawSAE["); pPpt->Out(row); pPpt->Out("]["); pPpt->Out(col); pPpt->Out("] = "); pPpt->Out(rawSAE[row][col]); pPpt->Out("\n");
						++col; // Tracks columns (speed, mfp, pr, eta)
						NEWLINE = false;
					}
					else
					{
						if(s[strptr]!=' ')
						{
							if (s[strptr] == '\n') {
								//pPpt->Out(Identify()); pPpt->Out(".LoadMapSAE: "); pPpt->Out("New line character (\n) found in "); pPpt->Out(InputFile); pPpt->Out("\n");
								//pPpt->Out("New line character\n");
								NEWLINE = true;
							}
							else {
								pPpt->Out(Identify()); pPpt->Out(".LoadMapSAE: "); pPpt->Out("Unrecognized character ("); pPpt->Out(s[strptr]);
								pPpt->Out(") found in "); pPpt->Out(InputFile); pPpt->Out("\n");
								exit(1);
							}
						}
						else {
							//cout << "Found a space\n";
							NEWLINE = false;
						}
					}
				}
				//cout << "Before ++strptr: ";
				//if (s[strptr] == NULL) cout << "s[strptr=" << strptr << "] == NULL\n\n";
				//else cout << "s[strptr=" << strptr << "] == " << s[strptr] << "\n\n";
				++strptr;
				//cout << "After ++strptr: ";
				//if (s[strptr] == NULL) cout << "s[strptr=" << strptr << "] == NULL\n\n";
				//else cout << "s[strptr=" << strptr << "] == " << s[strptr] << "\n\n";
			};

			// If there was a number right before the NULL, must do this:
			if(found_digit)		// Come to end of number
			{
				found_digit = false;
				// number runs from s[number_start] to s[number_start+number_pos]

				// Now move pointer to before this number
				temp_pos = start_of_string + number_start;
				fsetpos(stream, &temp_pos);

				// Read number into float
				fscanf(stream, "%f", &fp);
						
				rawSAE[row][col] = fp;
				//pPpt->Out("Right before NULL:\n");
				//pPpt->Out("Found number: "); pPpt->Out(fp); pPpt->Out("\n");
				//pPpt->Out("rawSAE["); pPpt->Out(row); pPpt->Out("]["); pPpt->Out(col); pPpt->Out("] = "); pPpt->Out(rawSAE[row][col]); pPpt->Out("\n");
			}
			//cout << "Now doing ++row. row was = " << row;
			if (!NEWLINE) ++row;
			//cout << ". row now = " << row << "\n";
			col=0;
		}while(fscanf(stream, "%l")!=EOF);
		fclose(stream);
	}
	delete [] s;

	// Now interpret rawSAE into separate speed lines
	speeds = 0;
	col = 0;
	for(row=0; row<rows; ++row)
	{
		if(row==0) speeds = 1;
		else
		{
			if(rawSAE[row][col] > ((1.0 + (tolSpeed/100))*rawSAE[row-1][col])) ++speeds;
		}
	}
	
	SAE = new double** [speeds];
	numPointsSAE = new int [speeds];
	int speed=0;
	int data_point_counter = 0;
	int tempPoint;
	for(row=0; row<rows; ++row)
	{
		if(row!=0)
		{
			if(rawSAE[row][col] > ((1.0 + (tolSpeed/100))*rawSAE[row-1][col]))
			{
				numPointsSAE[speed] = data_point_counter;
				SAE[speed] = new double* [numPointsSAE[speed]];
				for(tempPoint=0; tempPoint<numPointsSAE[speed]; ++tempPoint) SAE[speed][tempPoint] = new double [4];				
				data_point_counter = 0;
				++speed;
			}
		}
		++data_point_counter;
		
		if(row==rows-1)
		{
			numPointsSAE[speed] = data_point_counter;
			SAE[speed] = new double* [numPointsSAE[speed]];
			for(tempPoint=0; tempPoint<numPointsSAE[speed]; ++tempPoint) SAE[speed][tempPoint] = new double [4];
		}
	}

	speed=0;
	tempPoint=0;
	int totalPoints = 0;
	for(row=0; row<rows; ++row)
	{
		if( (row - totalPoints) > (numPointsSAE[speed] - 1) )
		{
			totalPoints += numPointsSAE[speed];
			++speed;
			tempPoint=0;
		}
		SAE[speed][tempPoint][0] = rawSAE[row][0]/60;	// Speed: convertSAE rpm to rps
		SAE[speed][tempPoint][1] = rawSAE[row][1]/1000;	// MFP: convert SAE kPa format to pure Pa
		SAE[speed][tempPoint][2] = rawSAE[row][2];		// PR
		SAE[speed][tempPoint][3] = rawSAE[row][3];		// eta
		++tempPoint;
	}

	// Convert SAE array into existing array format
	int num_graphs = 1;
	num_curves = new int [num_graphs];
	num_points = new int* [num_graphs];
	int graph;
	Raw = new double*** [num_graphs];
	for(graph=0; graph<num_graphs; ++graph)
	{
		num_curves[graph] = speeds;
		num_points[graph] = new int [speeds];
		Raw[graph] = new double** [speeds];
		for(speed=0; speed<speeds; ++speed)
		{
			num_points[graph][speed] = numPointsSAE[speed];
			Raw[graph][speed] = new double* [numPointsSAE[speed]];
			for(tempPoint=0; tempPoint<numPointsSAE[speed]; ++tempPoint)
			{
				Raw[graph][speed][tempPoint] = new double [5]; // 5 pieces of data per point (Type, S, MFP, PR, eta)
				Raw[graph][speed][tempPoint][0] = 1.0;
				Raw[graph][speed][tempPoint][1] = SAE[speed][tempPoint][0];
				Raw[graph][speed][tempPoint][2] = SAE[speed][tempPoint][1];
				Raw[graph][speed][tempPoint][3] = SAE[speed][tempPoint][2];
				Raw[graph][speed][tempPoint][4] = SAE[speed][tempPoint][3];
			}
		}
	}
}

void CTurbine::ProcessMapVar(CProperties* pPpt, int num_extra_reverse, int num_extra_choked, char* proc_map_file_str, bool TURBINE, double T)
// ============================================================ //
// For variable outlet pressure turbines.						//
// Processes input PR and MFP data into Riemann variables.		//
// The first processed point is set as lambda_in_star_min and	//
// the last as lambda_in_star_max.								//
// This function is called from CTurbine::Configure				//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ProcessMapVar\n");}

	int UPSTREAM = INLET;
	int DOWNSTREAM = OUTLET;

	int POINT_RAW;
	int ANY_POINT = 0;
	double f [2]; // Inlet and exit duct area

	f[INLET] = pBN[INLET]->f;
	f[OUTLET] = pBN[OUTLET]->f;
	C[INLET] = ((pPpt->gammaAir(T)-1)/2)*sqrt(pPpt->R_air/pPpt->gammaAir(T))*(1/f[INLET]);
	C[OUTLET] = ((pPpt->gammaAir(T)-1)/2)*sqrt(pPpt->R_air/pPpt->gammaAir(T))*(1/f[OUTLET]);
	
	int AXIAL = 1;
	int TYPE = AXIAL;
	int GRAPH = TYPE - 1;
	Data = new CDataPoint* [num_curves[GRAPH]];

	int CURVE, POINT;
 	for (CURVE=0; CURVE<num_curves[GRAPH]; ++CURVE) Data[CURVE] = new CDataPoint [num_extra_reverse + num_points[GRAPH][CURVE] + num_extra_choked];
 	
	x_min		= new double [num_curves[GRAPH]];
	Pt			= new double [num_curves[GRAPH]];
	theta_t		= new double [num_curves[GRAPH]];
	G1_max		= new double [num_curves[GRAPH]];
	K1			= new double [num_curves[GRAPH]];
	K2			= new double [num_curves[GRAPH]];
	K3			= new double [num_curves[GRAPH]];
	K4			= new double [num_curves[GRAPH]];
	x_max		= new double [num_curves[GRAPH]];

	double tempT = 300;

	for (CURVE=0; CURVE<num_curves[GRAPH]; ++CURVE) {
		// For the 'real' points only
		for (POINT = (0 + num_extra_reverse); POINT < (num_extra_reverse + num_points[GRAPH][CURVE]) ; ++POINT) {
			POINT_RAW = POINT - num_extra_reverse;
			Data[CURVE][POINT].sp = Raw[GRAPH][CURVE][POINT_RAW][SP];
			Data[CURVE][POINT].mfp = Raw[GRAPH][CURVE][POINT_RAW][MFP];
			Data[CURVE][POINT].pr = Raw[GRAPH][CURVE][POINT_RAW][PR];
			Data[CURVE][POINT].eta = Raw[GRAPH][CURVE][POINT_RAW][ETA];

			Data[CURVE][POINT].power_param = (1/Data[CURVE][POINT].eta)*(Data[CURVE][POINT].mfp)*pPpt->cpAir(tempT)
											 *(pow(Data[CURVE][POINT].pr, (pPpt->gammaAir(T)-1)/pPpt->gammaAir(T)) - 1);

			if(Data[CURVE][POINT].sp > 1e-12) {
				Data[CURVE][POINT].torque_param = (1/Data[CURVE][POINT].eta)*(Data[CURVE][POINT].mfp)*(1/Data[CURVE][POINT].sp)*(pPpt->cpAir(tempT)/(2*PI))
											 *(pow(Data[CURVE][POINT].pr, (pPpt->gammaAir(T)-1)/pPpt->gammaAir(T)) - 1);
			}
			else Data[CURVE][POINT].torque_param = 0;
		}
	}

	// Create/open the file to which to print the processed map and print some header information
	FILE* mapfileptr = fopen(ConstructString(pPpt, RES_DIR, proc_map_file_str), "w");
	string str_preamble = "Processed turbine map file for " + CharToString(Identify());
	fprintf(mapfileptr, Underline(str_preamble, "*", true)); fprintf(mapfileptr, "\n");

	for (CURVE=0; CURVE<num_curves[GRAPH]; ++CURVE) {
		fprintf(mapfileptr, Underline("SPEED LINE " + IntToString(CURVE+1) + " OF " + IntToString(num_curves[GRAPH]), "=")); fprintf(mapfileptr, "\n");
		fprintf(mapfileptr, Underline("REAL/MEASURED POINTS", "-"));

		fprintf(mapfileptr,"\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
			"POINT", "SPEED PARAM", "MFP     ", "PR      ", "ETA_TS   ", "G1      ", "p1_over_p2", "T2_over_T1", "lambda_in_star_ratio");

		// Process each real data point on this curve
		for (POINT = (0 + num_extra_reverse); POINT < (num_extra_reverse + num_points[GRAPH][CURVE]); ++POINT) {
			Data[CURVE][POINT].ProcessDataPtVar(pPpt, f, C, TURBINE, T);
/*
			if(Data[CURVE][POINT].G1 == 0)
			{
				cout << "Data[CURVE][POINT].G1 = " << Data[CURVE][POINT].G1 << endl;
				cout << "Data[CURVE][POINT].p1_over_p2 = " << Data[CURVE][POINT].p1_over_p2 << endl;
				cout << "Data[CURVE][POINT].lambda_in_star_ratio = " << Data[CURVE][POINT].lambda_in_star_ratio << endl;
			}
//*/			
			// Write the processed real data points to file
			fprintf(mapfileptr,"%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", 
				POINT, Data[CURVE][POINT].sp, Data[CURVE][POINT].mfp*1e3, Data[CURVE][POINT].pr, Data[CURVE][POINT].eta, 
				Data[CURVE][POINT].G1, Data[CURVE][POINT].p1_over_p2, Data[CURVE][POINT].T2_over_T1, Data[CURVE][POINT].lambda_in_star_ratio);
/*
cout << "Data[CURVE][POINT].sp = " << Data[CURVE][POINT].sp << endl;
cout << "Data[CURVE][POINT].mfp = " << Data[CURVE][POINT].mfp << endl;
cout << "Data[CURVE][POINT].pr = " << Data[CURVE][POINT].pr << endl;
cout << "Data[CURVE][POINT].eta = " << Data[CURVE][POINT].eta << endl;
cout << "Data[CURVE][POINT].power_param = " << Data[CURVE][POINT].power_param << endl;
cout << "Data[CURVE][POINT].torque_param = " << Data[CURVE][POINT].torque_param << endl;
cout << "Data[CURVE][POINT].G1 = " << Data[CURVE][POINT].G1 << endl;
cout << "Data[CURVE][POINT].p1_over_p2 = " << Data[CURVE][POINT].p1_over_p2 << endl;
cout << "Data[CURVE][POINT].T2_over_T1 = " << Data[CURVE][POINT].T2_over_T1 << endl;
cout << "Data[CURVE][POINT].lambda_in_star_ratio = " << Data[CURVE][POINT].lambda_in_star_ratio << endl;
cout << endl;
//*/
		}
		fprintf(mapfileptr, "\n\n");

		// Reverse flow data points - for each curve
		// =========================================
		double G1_test = 0;
		double G1_test_min = -0.2e-3;
		double del_G1_test = fabs(0 - G1_test_min)/(num_extra_reverse - 1);
		double T2_over_T1_temp = 1.0;
		double lambda_in_star_ratio_temp;
	
		// Assume the first 'real' data point is that for zero flow
		POINT = 0  + num_extra_reverse;
		double pr_min = Data[CURVE][POINT].p1_over_p2;
		pr_min = 1;
		x_min[CURVE] = pow(pr_min, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)));

		fprintf(mapfileptr, Underline("REVERSE/ZERO FLOW POINTS", "-"));
		fprintf(mapfileptr,"\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
			"POINT", "           ", "        ", "        ", "        ", " G1      ", "p1_over_p2", "T2_over_T1", "lambda_in_star_ratio");
		
		POINT = 0  + num_extra_reverse - 1; // Move to the last (reverse/zero flow) point, i.e. the one just before the first 'real' point
		do {
			lambda_in_star_ratio_temp = (x_min[CURVE]*(1 + C[UPSTREAM]*G1_test))/(1 - C[DOWNSTREAM]*G1_test*pr_min);
				
			// Store this new data point
			Data[CURVE][POINT].sp = Raw[GRAPH][CURVE][ANY_POINT][SP];
			Data[CURVE][POINT].mfp = 0;
			
			//double tempc, tempM1, tempfM1, tempdel_M1, tempp01_over_p1;
			//tempc = pow( ( (Data[CURVE][POINT].mfp/f[UPSTREAM])*sqrt(pPpt->R_air/pPpt->gammaAir(T)) ), (2*(pPpt->gammaAir(T)-1)/(pPpt->gammaAir(T)+1)) );
			//tempM1 = 0.001;
			//do
			//{
			//	tempfM1 = pow( tempc + ((pPpt->gammaAir(T)-1)/2)*tempc*pow(tempM1,2), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1)) );
			//	tempdel_M1 = /*fabs*/(tempM1 - tempfM1);
			//	tempM1 = tempfM1;
			//}
			//while(fabs(tempdel_M1)>0.0001);
			//tempp01_over_p1 = pow( 1 + ((pPpt->gammaAir(T)-1)/2)*pow(tempM1,2), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1) );
			//Data[CURVE][POINT].pr = p1_over_p2_test*tempp01_over_p1;

			Data[CURVE][POINT].pr = pr_min; // MFP is not necessarily 0 here
					
			Data[CURVE][POINT].eta = 1; // Correct?

			Data[CURVE][POINT].power_param = (1/Data[CURVE][POINT].eta)*(Data[CURVE][POINT].mfp)*pPpt->cpAir(tempT)
											 *(pow(Data[CURVE][POINT].pr, (pPpt->gammaAir(T)-1)/pPpt->gammaAir(T)) - 1);

			if(Data[CURVE][POINT].sp > 1e-12) {
				Data[CURVE][POINT].torque_param = (1/Data[CURVE][POINT].eta)*(Data[CURVE][POINT].mfp)*(1/Data[CURVE][POINT].sp)*(pPpt->cpAir(tempT)/(2*PI))
											 *(pow(Data[CURVE][POINT].pr, (pPpt->gammaAir(T)-1)/pPpt->gammaAir(T)) - 1);
			}
			else Data[CURVE][POINT].torque_param = 0;

			Data[CURVE][POINT].G1 = G1_test;
			Data[CURVE][POINT].T2_over_T1 = T2_over_T1_temp;
			Data[CURVE][POINT].p1_over_p2 = pr_min;
			//Data[CURVE][POINT].p1_over_p2 = 1.0;
			Data[CURVE][POINT].lambda_in_star_ratio = lambda_in_star_ratio_temp;
			Data[CURVE][POINT].query_param = Data[CURVE][POINT].lambda_in_star_ratio;

			// Write the processed reverse/zero flow data points to file
			fprintf(mapfileptr,"%d\t%s\t%s\t%s\t%s\t%.10f\t%.10f\t%.10f\t%.10f\n", 
				POINT, "          ", "          ", "          ", "          ", 
				Data[CURVE][POINT].G1, Data[CURVE][POINT].p1_over_p2, Data[CURVE][POINT].T2_over_T1, Data[CURVE][POINT].lambda_in_star_ratio);
		
			G1_test -= del_G1_test;
			--POINT;
/*
cout << "Data[CURVE][POINT].sp = " << Data[CURVE][POINT].sp << endl;
cout << "Data[CURVE][POINT].mfp = " << Data[CURVE][POINT].mfp << endl;
cout << "Data[CURVE][POINT].pr = " << Data[CURVE][POINT].pr << endl;
cout << "Data[CURVE][POINT].eta = " << Data[CURVE][POINT].eta << endl;
cout << "Data[CURVE][POINT].power_param = " << Data[CURVE][POINT].power_param << endl;
cout << "Data[CURVE][POINT].torque_param = " << Data[CURVE][POINT].torque_param << endl;
cout << "Data[CURVE][POINT].G1 = " << Data[CURVE][POINT].G1 << endl;
cout << "Data[CURVE][POINT].p1_over_p2 = " << Data[CURVE][POINT].p1_over_p2 << endl;
cout << "Data[CURVE][POINT].T2_over_T1 = " << Data[CURVE][POINT].T2_over_T1 << endl;
cout << "Data[CURVE][POINT].lambda_in_star_ratio = " << Data[CURVE][POINT].lambda_in_star_ratio << endl;
cout << endl;
//*/
		}while(POINT>=0);
		fprintf(mapfileptr, "\n\n"); // New line

		// Last and extra data points - choked flow data - for each curve
		// ==============================================================
		POINT = num_extra_reverse + num_points[GRAPH][CURVE] - 1; // Move to the last real point on the curve
		
		double maxMFP = Data[CURVE][POINT].mfp;	// Note maximum MFP from real data
		Pt[CURVE] = Data[CURVE][POINT].p1_over_p2; // The transition pressure ratio
		theta_t[CURVE] = Data[CURVE][POINT].T2_over_T1;
		G1_max[CURVE] = Data[CURVE][POINT].G1;
		K1[CURVE] = 1 + C[UPSTREAM]*G1_max[CURVE];
		K2[CURVE] = 1 - C[UPSTREAM]*G1_max[CURVE];
		K3[CURVE] = C[UPSTREAM]*G1_max[CURVE];
		K4[CURVE] = C[DOWNSTREAM]*G1_max[CURVE];		
		x_max[CURVE] = K1[CURVE]*pow(Pt[CURVE], (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T))) / (1 - K4[CURVE]*Pt[CURVE]*sqrt(theta_t[CURVE]));

		// Now process choked flow extra data, starting with p1_over_p2 at Pt + 1 increment, increasing by a suitable increment, del_p1_over_p2_test, each loop
		double p1_over_p2_test = Pt[CURVE];
		double p1_over_p2_max = Pt[CURVE]*1.25;

		// Given the extra number of data points available for choked data, work out the appropriate increment
		double del_p1_over_p2_test = (p1_over_p2_max - Pt[CURVE])/num_extra_choked;

		fprintf(mapfileptr, Underline("EXTRA/CHOKED POINTS", "-"));
		fprintf(mapfileptr,"\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
			"POINT", "           ", "        ", "        ", "        ", " G1      ", "p1_over_p2", "T2_over_T1", "lambda_in_star_ratio");
		
		++POINT;
		do {
			p1_over_p2_test += del_p1_over_p2_test;

			lambda_in_star_ratio_temp = (K1[CURVE]*pow(p1_over_p2_test, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T))))
														/(1 - K4[CURVE]*p1_over_p2_test*sqrt(theta_t[CURVE]));
			
			Data[CURVE][POINT].sp = Raw[GRAPH][CURVE][ANY_POINT][SP];
			Data[CURVE][POINT].mfp = maxMFP; // Set choked data point MFP to that of the last real data point
			

			double tempc, tempM1, tempfM1, tempdel_M1, tempp01_over_p1;
			tempc = pow( ( (Data[CURVE][POINT].mfp/f[UPSTREAM])*sqrt(pPpt->R_air/pPpt->gammaAir(T)) ), (2*(pPpt->gammaAir(T)-1)/(pPpt->gammaAir(T)+1)) );
			tempM1 = 0.001;
			do {
				tempfM1 = pow( tempc + ((pPpt->gammaAir(T)-1)/2)*tempc*pow(tempM1,2), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1)) );
				tempdel_M1 = /*fabs*/(tempM1 - tempfM1);
				tempM1 = tempfM1;
			}
			while(fabs(tempdel_M1)>0.0001);
			tempp01_over_p1 = pow( 1 + ((pPpt->gammaAir(T)-1)/2)*pow(tempM1,2), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1) );
			Data[CURVE][POINT].pr = p1_over_p2_test*tempp01_over_p1;
					
			Data[CURVE][POINT].eta = 1; // Correct??

			Data[CURVE][POINT].power_param = (1/Data[CURVE][POINT].eta)*(Data[CURVE][POINT].mfp)*pPpt->cpAir(tempT)
											 *(pow(Data[CURVE][POINT].pr, (pPpt->gammaAir(T)-1)/pPpt->gammaAir(T)) - 1);

			if(Data[CURVE][POINT].sp > 1e-12) {
				Data[CURVE][POINT].torque_param = (1/Data[CURVE][POINT].eta)*(Data[CURVE][POINT].mfp)*(1/Data[CURVE][POINT].sp)*(pPpt->cpAir(tempT)/(2*PI))
											 *(pow(Data[CURVE][POINT].pr, (pPpt->gammaAir(T)-1)/pPpt->gammaAir(T)) - 1);
			}
			else Data[CURVE][POINT].torque_param = 0;

			//Data[CURVE][POINT].G1 = rG1_max;
			//Data[CURVE][POINT].T2_over_T1 = rtheta_t;
			Data[CURVE][POINT].G1 = G1_max[CURVE];
			Data[CURVE][POINT].T2_over_T1 = theta_t[CURVE];
			Data[CURVE][POINT].p1_over_p2 = p1_over_p2_test;
			Data[CURVE][POINT].lambda_in_star_ratio = lambda_in_star_ratio_temp;
			Data[CURVE][POINT].query_param = Data[CURVE][POINT].lambda_in_star_ratio;

			// Write the processed extra/choked data points to file
			fprintf(mapfileptr,"%d\t\t\t\t\t\t\t%.10f\t%.10f\t%.10f\t%.10f\n", 
				POINT,
				Data[CURVE][POINT].G1, Data[CURVE][POINT].p1_over_p2, Data[CURVE][POINT].T2_over_T1, Data[CURVE][POINT].lambda_in_star_ratio);

			++POINT;
		}while(POINT < (num_extra_reverse + num_points[GRAPH][CURVE] + num_extra_choked));
		fprintf(mapfileptr, "\n\n");
	
		// Print finalised form of data to file, for each curve
		// ====================================================
		//fprintf(mapfileptr,"\nFINAL LIST - CURVE ");
		//fprintf(mapfileptr,"%d\n", CURVE);
		//fprintf(mapfileptr,"====================\n");
		fprintf(mapfileptr, Underline("FINAL LIST - SPEED LINE " + IntToString(CURVE+1) + " OF " + IntToString(num_curves[GRAPH]), "-"));

		fprintf(mapfileptr,	"\nPOINT\tG1\tp1_over_p2\tT2_over_T1\tlambda_in_star_ratio\n");
		for(POINT=0; POINT < (num_extra_reverse + num_points[GRAPH][CURVE] + num_extra_choked); ++POINT)
			fprintf(mapfileptr,"%d\t%.10f\t%.10f\t%.10f\t%.10f\n", 
				POINT,
				Data[CURVE][POINT].G1, Data[CURVE][POINT].p1_over_p2, Data[CURVE][POINT].T2_over_T1, Data[CURVE][POINT].lambda_in_star_ratio);
		fprintf(mapfileptr, "\n\n\n");
	}
	fclose(mapfileptr);
}

void CTurbine::ProcessMapConst(CProperties* pPpt, double T)
// ============================================================ //
// For constant outlet pressure turbines only:					//
// Processes input PR and MFP data into Riemann variables.		//
// The first processed point is set as lambda_in_star_min and	//
// the last as lambda_in_star_max.								//
// This function is called from CTurbine::Configure				//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ProcessMapConst\n");}

	int AXIAL = 1;
	int TYPE = AXIAL;
	int GRAPH = TYPE - 1;
	Data = new CDataPoint* [num_curves[GRAPH]];

	int CURVE, POINT;
 	for(CURVE=0; CURVE<num_curves[GRAPH]; ++CURVE) Data[CURVE] = new CDataPoint [num_points[GRAPH][CURVE]];
 	
	for(CURVE=0; CURVE<num_curves[GRAPH]; ++CURVE) {
		for (POINT=0; POINT<num_points[GRAPH][CURVE]; ++POINT) {
			Data[CURVE][POINT].sp = Raw[GRAPH][CURVE][POINT][SP];
			Data[CURVE][POINT].mfp =  Raw[GRAPH][CURVE][POINT][MFP];
			Data[CURVE][POINT].pr = Raw[GRAPH][CURVE][POINT][PR];
			Data[CURVE][POINT].eta = Raw[GRAPH][CURVE][POINT][ETA];
		}
	}
/*
	// Minimum lambda_in_star must be at zero mfp, though the data doesn't necessarily have to have a zero mfp point
	CURVE=0; POINT=0;
	lambda_in_star_min = pow(Data[CURVE][POINT].pr, (pPpt->gammaAir()-1)/(2*pPpt->gammaAir()));
	// This will extrapolate the first data point to get a pr value for zero mfp (if mfp here is not zero), in order to calculate lambda_in_star_min
*/
	//std::cout << ConstructString(pPpt, RES_DIR, "proc_map_turb.txt") << std::endl;
	FILE* mapfileptr = fopen(ConstructString(pPpt, RES_DIR, "proc_map_turb.txt"), "w");

	for(CURVE=0; CURVE<num_curves[GRAPH]; ++CURVE) {
		fprintf(mapfileptr,"%s\t%s\t%s\t%s\t%s\t%s\n", "SPEED PARAM", "MFP", "PR", "ETA_TS", "LAMBDA_IN_STAR", "LAMBDA_OUT_STAR");
//		fprintf(mapfileptr,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", 0, 0, 0, lambda_in_star_min, 0);
		
		for(POINT=0; POINT<num_points[GRAPH][CURVE]; ++POINT)	// Process each data point
		{
			ProcessDataPtConst(pPpt, Data[CURVE][POINT].pr, 
											Data[CURVE][POINT].mfp, 
												Data[CURVE][POINT].lambda_in_star, 
													Data[CURVE][POINT].lambda_out_star, T);
			// In all cases:
			Data[CURVE][POINT].query_param = Data[CURVE][POINT].lambda_in_star;

			// Record this point on the map being read in
			fprintf(mapfileptr,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", Data[CURVE][POINT].sp, Data[CURVE][POINT].mfp, Data[CURVE][POINT].pr, Data[CURVE][POINT].eta, Data[CURVE][POINT].lambda_in_star, Data[CURVE][POINT].lambda_out_star);
			//printf("%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n", Data[CURVE][POINT].sp, Data[CURVE][POINT].mfp, Data[CURVE][POINT].pr, Data[CURVE][POINT].eta, Data[CURVE][POINT].lambda_in_star, Data[CURVE][POINT].lambda_out_star);
		}	

		// Minimum lambda_in_star must be at zero mfp, though the data doesn't necessarily have to have a zero mfp point
		POINT=0;
		double lambda_in_star_min_temp = pow(Data[CURVE][POINT].pr, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)));
		// This will extrapolate the first data point to get a pr value for zero mfp (if mfp here is not zero), in order to calculate lambda_in_star_min

// Or assume minimum PR is 1, where MFP = 0
lambda_in_star_min_temp = 1.0;

		// Last data point on curve is set as max
		POINT = num_points[GRAPH][CURVE] - 1; // Move to the last point on the curve
		double lambda_in_star_max_temp = Data[CURVE][POINT].lambda_in_star;	
		// Constant Ks, for choked flow, calculated from lata data point
		double Ks_temp = (Data[CURVE][POINT].lambda_out_star)
								/(Data[CURVE][POINT].lambda_in_star);

		for(POINT=0; POINT<num_points[GRAPH][CURVE]; ++POINT)	// Process each data point
		{
			Data[CURVE][POINT].lambda_in_star_min = lambda_in_star_min_temp;
			Data[CURVE][POINT].lambda_in_star_max = lambda_in_star_max_temp;
			Data[CURVE][POINT].Ks = Ks_temp;
		}
	}
	fclose(mapfileptr);
/*
	// Last data point is set as max
	CURVE = num_curves[GRAPH] - 1; // Move to the last curve
	POINT = num_points[GRAPH][CURVE] - 1; // Move to the last point on the last curve
	lambda_in_star_max = Data[CURVE][POINT].lambda_in_star;	

	// Constant Ks, for choked flow, calculated from lata data point
	Ks = (Data[num_curves[GRAPH]-1][num_points[GRAPH][CURVE]-1].lambda_out_star)
			/(Data[num_curves[GRAPH]-1][num_points[GRAPH][CURVE]-1].lambda_in_star);
*/
}

void CTurbine::ProcessDataPtConst(CProperties* pPpt, double pr, double mfp, double &rlambda_in_star1, double &rlambda_out_star1, double T)
// ============================================================ //
// Processes input PR and MFP data into Riemann variables		//
// for a single data point.										//
// This function is called from CTurbine::ProcessMapConst		//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ProcessDataPtConst\n");}

	double C, M, fM, del_M;
	double inlet_ratio, inlet_to_exit_ratio;

	if(mfp<=0.0) // Reverse or no flow
	{
		double lambda_in_star_min = pow(pr, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T)));
		rlambda_in_star1 = lambda_in_star_min;
		rlambda_out_star1 = rlambda_in_star1; // Treat as closed end
	}
	else
	{
		C = pow( ( (mfp/pBN[INLET]->f)*sqrt(pPpt->R_air/pPpt->gammaAir(T)) ), (2*(pPpt->gammaAir(T)-1)/(pPpt->gammaAir(T)+1)) );
		M = 0.001;

		do
		{
			fM = pow( C + ((pPpt->gammaAir(T)-1)/2)*C*pow(M,2), (pPpt->gammaAir(T)+1)/(2*(pPpt->gammaAir(T)-1)) );
			del_M = /*fabs*/(M - fM);
			M = fM;
		}
		//while(del_M>0.0001);
		while(fabs(del_M)>1e-12);

		inlet_ratio = pow( 1 + ((pPpt->gammaAir(T)-1)/2)*pow(M,2), pPpt->gammaAir(T)/(pPpt->gammaAir(T)-1) ); //p01/p1
		inlet_to_exit_ratio = (1/inlet_ratio)*pr; // (p1/p2) = (p1/p01)*(p01/p2)
		
		rlambda_in_star1 = pow(inlet_to_exit_ratio, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T))) * (1 + ((pPpt->gammaAir(T)-1)/2)*M);
		rlambda_out_star1 = pow(inlet_to_exit_ratio, (pPpt->gammaAir(T)-1)/(2*pPpt->gammaAir(T))) * (1 - ((pPpt->gammaAir(T)-1)/2)*M);
	}
}

void CTurbine::RunBoundary(CProperties* pPpt, int timestep, double time, double DELZ)
// ============================================================ //
// The entry point into this boundary method. Called once per	//
// iteration per object.										//
// This function is called from the main function.				//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".RunBoundary\n");}

	double del_t;
	del_t = (DELZ/(EX ? pPpt->AREFe : pPpt->AREFi))*pPpt->xref;
	
	// Update inlet thermodynamic properties
	T01 = pBN[INLET]->T + 0.5*pow(pBN[INLET]->U*pPipe[INLET]->AREF,2)/pPpt->Cp_air;
	T1 = pBN[INLET]->T;

	if(VARIABLE) VariablePressure(pPpt/*, INLET, OUTLET,*//* true, C,*//* time*/); // Variable outlet pressure	
	else {	// Constant outlet pressure methods
		if (EQUIVALENT_AREA) EquivArea(pPpt, timestep, time); // Equivalent area method (intended for axial turbines) and imposes constant outlet pressure
		else ConstantPressure(pPpt);
	}

	InstantaneousWorkAndMassFlow(pPpt, del_t);
	InstantaneousTransmissionLoss(pPpt, del_t);		
	Matching(pPpt, del_t); // Resolve torque balance to update turbocharger speed
}

void CTurbine::VariablePressure(CProperties* pPpt/*, int& rINLET, int& rOUTLET,*//* bool TURBINE,*//* double* &rC,*//* double time*/)
// ============================================================ //
// Turbine model permitting variable							//
// downstream conditions e.g., variable back/outlet pressure.	//
// ============================================================ //																	//	
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".VariablePressure\n");}

	char pause;

	double *lambda_in_n, *AA_n, *lambda_in_star, *lambda_out_star, *lambda_out, *AA_c, *lambda_in_c, *lambda_in_c_old;
	int* pipe_flow;

	lambda_in_n = new double [NPIPES];
	AA_n = new double [NPIPES];
	lambda_in_star = new double [NPIPES];
	lambda_out_star = new double [NPIPES];
	lambda_out = new double [NPIPES];
	AA_c = new double [NPIPES];
	lambda_in_c = new double [NPIPES];
	lambda_in_c_old = new double [NPIPES];
	pipe_flow = new int [NPIPES];
	
	double x, error, limit;
	limit = 1e-6;
	
	int counter;

	int UPSTREAM, DOWNSTREAM;
	int INLET = 0;
	int OUTLET = 1;

	// Enter inital characteristic values
	// ----------------------------------------------------------------------------------------------------
	lambda_in_n[INLET] = (*(pCLIN[INLET]))[1];
	lambda_in_n[OUTLET] = (*(pCLIN[OUTLET]))[1];
	AA_n[INLET] = (pBN[INLET])->AA[1];
	AA_n[OUTLET] = (pBN[OUTLET])->AA[1];
	lambda_in_star[INLET] = lambda_in_n[INLET]/AA_n[INLET];
	lambda_in_star[OUTLET] = lambda_in_n[OUTLET]/AA_n[OUTLET];


	// Direction test
	// ----------------------------------------------------------------------------------------------------
	if(lambda_in_star[INLET] - lambda_in_star[OUTLET] > pPpt->ZERO_TOL) {
		// Conventional flow
		UPSTREAM = INLET;
		DOWNSTREAM = OUTLET;

		if(REVERSE_FLOW) {
//			cout << "Flow in variable pressure device is changing direction, to conventional flow\n";
//			cin >> pause;
		}
		REVERSE_FLOW = false;
	}
	else {
		if(lambda_in_star[OUTLET] - lambda_in_star[INLET] > pPpt->ZERO_TOL) {
			// Reverse flow
			UPSTREAM = OUTLET;
			DOWNSTREAM = INLET;

			if(!REVERSE_FLOW) {
//				cout << "Flow in variable pressure device is changing direction, to reverse flow\n";
//				cin >> pause;
			}
			REVERSE_FLOW = true;
		}
		else { // No flow
			// No flow but use the conventional direction
			UPSTREAM = INLET;
			DOWNSTREAM = OUTLET;
			REVERSE_FLOW = false;

			// No flow - treat both as closed ends
			pipe_flow[INLET] = NOFLOW;
			pipe_flow[OUTLET] = NOFLOW;
//			(*(pCLOUT[INLET]))[1] = (*(pCLIN[INLET]))[1];
//			(*(pCLOUT[OUTLET]))[1] = (*(pCLIN[OUTLET]))[1];
//			// No need to update lambda_in or AA
//			return;
		}
	}

	// Set flow directions (this is now a non-noflow case)
	// ----------------------------------------------------------------------------------------------------
	pipe_flow[UPSTREAM] = OUTFLOW;
	pipe_flow[DOWNSTREAM] = INFLOW;

	// Main algorithm loop
	// ----------------------------------------------------------------------------------------------------
	// Use uncorrected values as initial guess
	lambda_in_c[INLET] = lambda_in_n[INLET];
	lambda_in_c[OUTLET] = lambda_in_n[OUTLET];
	counter = 0;
	do {
		++counter;
		// Determine speed parameter and ratio x
		// ----------------------------------------------------------------------------------------------------
		x = lambda_in_star[INLET]/lambda_in_star[OUTLET];

		// Interpolate from the processed data
		// ----------------------------------------------------------------------------------------------------
		double speed_param;	
		//if(INLET_STAG) speed_param = NT/sqrt(T01);
		//else speed_param = NT/sqrt(T1);
		speed_param = NT/sqrt(T01);
		Interpolate(pPpt, speed_param, x);

		double lambda_in_star_ratio, G1, p1_over_p2, T2_over_T1;
		lambda_in_star_ratio = interp_data_pt.lambda_in_star_ratio;
		G1 = interp_data_pt.G1;
		p1_over_p2 = interp_data_pt.p1_over_p2;
		T2_over_T1 = interp_data_pt.T2_over_T1;
		
		if(!REVERSE_FLOW) // INLET lambda_in_star won't change
			lambda_in_star[OUTLET] = lambda_in_star[INLET]/lambda_in_star_ratio;
		else // OUTLET lambda_in_star won't change
			lambda_in_star[INLET] = lambda_in_star_ratio*lambda_in_star[OUTLET];

		// Evaluate equations
		// ----------------------------------------------------------------------------------------------------
		lambda_out_star[INLET] = lambda_in_star[INLET]*( (1 - C[INLET]*G1)/(1 + C[INLET]*G1) );
		lambda_out_star[OUTLET] = lambda_out_star[INLET]
									*pow( (1/p1_over_p2), (pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/(2*pPpt->gammaAir(pBN[UPSTREAM]->T)))
										*( (1 + C[OUTLET]*p1_over_p2*sqrt(T2_over_T1)*G1)
									      /(1 - C[INLET]*G1) );

		if(!REVERSE_FLOW) {
			AA_c[OUTLET] = AA_n[INLET]*sqrt(T2_over_T1)*pow(p1_over_p2, (pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/(2*pPpt->gammaAir(pBN[UPSTREAM]->T)));
			AA_c[INLET] = AA_n[INLET]; // No correction necessary for the inlet side
		}
		else { // Assume a throttling process across the turbine for reversed flow
			AA_c[OUTLET] = AA_n[OUTLET] ; // No correction necessary for the exit side
			AA_c[INLET] = AA_c[OUTLET]*pow(p1_over_p2, (pPpt->gammaAir(pBN[UPSTREAM]->T)-1)/(2*pPpt->gammaAir(pBN[UPSTREAM]->T)));
		}

		// Calculate corrected characteristics
		// ----------------------------------------------------------------------------------------------------
		lambda_out[INLET] = lambda_out_star[INLET]*AA_c[INLET];
		lambda_out[OUTLET] = lambda_out_star[OUTLET]*AA_c[OUTLET];

/*
		if(!REVERSE_FLOW)
		{	// INLET = UPSTREAM, OUTLET = DOWNSTREAM 
			lambda_in_c[INLET] = lambda_in_n[INLET];

			// OUTLET side need lambda_in correction as it's an INFLOW into the pipe
			lambda_in_c_old[OUTLET] = lambda_in_c[OUTLET];
			lambda_in_c[OUTLET] = lambda_in_n[OUTLET] + ((lambda_in_c[OUTLET] + lambda_out[OUTLET])/2)*((AA_c[OUTLET] - AA_n[OUTLET])/AA_c[OUTLET]);

			error = lambda_in_c[OUTLET] - lambda_in_c_old[OUTLET];
		}
		else // REVERSE FLOW
		{	// OUTLET = UPSTREAM, INLET = DOWNSTREAM 
			lambda_in_c[OUTLET] = lambda_in_n[OUTLET];

			// INLET side need lambda_in correction as it's an INFLOW into the pipe
			lambda_in_c_old[INLET] = lambda_in_c[INLET];
			lambda_in_c[INLET] = 	lambda_in_n[INLET] + ((lambda_in_c[INLET] + lambda_out[INLET])/2)*((AA_c[INLET] - AA_n[INLET])/AA_c[INLET]);

			error = lambda_in_c[INLET] - lambda_in_c_old[INLET];
		}
*/
		lambda_in_c[UPSTREAM] = lambda_in_n[UPSTREAM]; // UPSTREAM side needs no lambda_in correction as it's an OUTFLOW from the pipe
		lambda_in_c_old[DOWNSTREAM] = lambda_in_c[DOWNSTREAM]; // DOWNSTREAM side need lambda_in correction as it's an INFLOW into the pipe
		lambda_in_c[DOWNSTREAM] = lambda_in_n[DOWNSTREAM] + ((lambda_in_c[DOWNSTREAM] + lambda_out[DOWNSTREAM])/2)*((AA_c[DOWNSTREAM] - AA_n[DOWNSTREAM])/AA_c[DOWNSTREAM]);
		error = lambda_in_c[DOWNSTREAM] - lambda_in_c_old[DOWNSTREAM];

		if(counter>=100) {
			cout << "CTurbine::VariablePressure: Turbine [" << this->ID << "]: counter = " << counter << endl;
			cin >> pause;
		}
	}while(fabs(error) > limit);

	bool* CHOKED; CHOKED = new bool [NPIPES];
	for(int p=0; p<NPIPES; ++p) CHOKED[p] = false;

	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out, AA_c, pipe_flow, CHOKED); // Update boundary
}

void CTurbine::ConstantPressure(CProperties* pPpt)
// ============================================================ //
// A main calculation function for this boundary condition.		//
// This function is called from CTurbine::RunBoundary			//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ConstantPressure\n");}

	double AA_1, lambda_in_star1;
	double *lambda_in1, *AA_2, *lambda_out;
	int *pipe_flow;
	lambda_in1 = new double [NPIPES]; lambda_out = new double [NPIPES]; AA_2 = new double [NPIPES];
	pipe_flow = new int [NPIPES];

	// At each time step, lambda_in_1 and AA_2 are known (from the incoming pipe
	// characteristic and the constant back pressure value, respectively)

	// Set default flow directions
	// ===========================
	pipe_flow[ONE_SIDE] = OUTFLOW;
	REVERSE_FLOW = false;
	CHOKED_FLOW = false;

	// Enter inital characteristic values
	// ==================================
	lambda_in1[ONE_SIDE] = (*(pCLIN[ONE_SIDE]))[1];
	AA_1 = (pBN[ONE_SIDE])->AA[1];
	//lambda_in_star1 = lambda_in1[ONE_SIDE]/AA_1;
	AA_2[ONE_SIDE] = pow( (Pb/pPpt->PREF), (pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)/(2*pPpt->gammaAir(pBN[ONE_SIDE]->T)) );
	lambda_in_star1 = lambda_in1[ONE_SIDE]/AA_2[ONE_SIDE]; // Must use constant outlet entropy
/*
	// Test for choked flow
	// ====================
	if(lambda_in_star1 - interp_data_pt.lambda_in_star_max > 1e-6) 
	// Use old interpolated point value here (valid since speed parameter is same?)
	{
		CHOKED_FLOW = true;
		//cout << "CTurbine::ConstantOutletPressure: Turbine [" << this->ID << "]: choked flow\n";
		pipe_flow[ONE_SIDE] = OUTFLOW;				
		lambda_out = interp_data_pt.Ks*lambda_in1;
	}
	else
	{
		// Test for reverse flow (treat as open end)
		// =========================================
		if(lambda_in_star1 < interp_data_pt.lambda_in_star_min) // Use old interpolated point value here
		{
			REVERSE_FLOW = true;
			pipe_flow[ONE_SIDE] = INFLOW;
			lambda_out = 2*interp_data_pt.lambda_in_star_min*AA_2 - lambda_in1;

//			lambda_in1 = interp_data_pt.lambda_in_star_min*AA_2;
		}
		else
*/
		// Conventional flow (interpolate from the processed data)
		// =======================================================
		{
			pipe_flow[ONE_SIDE] = OUTFLOW;
			double speed_param;
			//if(INLET_STAG) speed_param = this->NT/sqrt(this->T01);
			//else speed_param = this->NT/sqrt(this->T1);
			speed_param = NT/sqrt(T01);

//cout << "speed_param = " << speed_param << endl;	
//cout << "lambda_in_star1 = " << lambda_in_star1 << endl;			
			Interpolate(pPpt, speed_param, lambda_in_star1);

//cout << "lambda_in_star1 = " << lambda_in_star1 << endl;
//interp_data_pt.PrintToScreen(pPpt, VARIABLE);

///*
			// Test for choked flow
			// ====================
			if(lambda_in_star1 - interp_data_pt.lambda_in_star_max > 1e-6)
			//if(lambda_in_star1 > interp_data_pt.lambda_in_star_max)
			{
				CHOKED_FLOW = true;
//cout << "CTurbine::ConstantOutletPressure: Turbine [" << this->ID << "]: choked flow\n";
				pipe_flow[ONE_SIDE] = OUTFLOW;				
				lambda_out[ONE_SIDE] = interp_data_pt.Ks*lambda_in1[ONE_SIDE];
				//lambda_out[ONE_SIDE] = interp_data_pt.Ks*lambda_in_star1*AA_2[ONE_SIDE];
			}
			else
			{
///*
				// Test for reverse flow (treat as open end)
				// =========================================
				if(lambda_in_star1<interp_data_pt.lambda_in_star_min)
				{
					REVERSE_FLOW = true;
//cout << "CTurbine::ConstantOutletPressure: Turbine [" << this->ID << "]: reverse flow\n";
					pipe_flow[ONE_SIDE] = INFLOW;
					lambda_out[ONE_SIDE] = 2*interp_data_pt.lambda_in_star_min*AA_2[ONE_SIDE] - lambda_in1[ONE_SIDE];

					//pipe_flow[ONE_SIDE] = NOFLOW;
					//lambda_out[ONE_SIDE] = lambda_in1[ONE_SIDE];
				}
				else
//*/
//*/
					// Conventional flow (used the interpolated point)
					// ===============================================
					lambda_out[ONE_SIDE] = interp_data_pt.lambda_out_star*AA_2[ONE_SIDE];
			}
		}
//	}
	bool* CHOKED; CHOKED = new bool [NPIPES]; CHOKED[ONE_SIDE] = false;
	common_UPDATE_NH(pPpt, lambda_in1, lambda_out, AA_2, pipe_flow, CHOKED);
}

void CTurbine::EquivArea(CProperties* pPpt, int timestep, double time)
// ====================================================================================================
// Routine treating turbine as an equivalent area nozzle
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".EquivArea\n");}

	double lambda_in_n, lambda_out_n, AA_n;
	double *lambda_in_c, *lambda_out_c, *AA_c;
	int *pipe_flow;
	lambda_in_c = new double [NPIPES]; lambda_out_c = new double [NPIPES]; AA_c = new double [NPIPES];	pipe_flow = new int [NPIPES];
	bool* CHOKED; CHOKED = new bool [NPIPES]; 
	CHOKED[ONE_SIDE] = false; // Reset here for all cases

	lambda_in_n = (*(pCLIN[ONE_SIDE]))[R+1];
	lambda_out_n = (*(pCLOUT[ONE_SIDE]))[R+1];
	AA_n = pBN[ONE_SIDE]->AA[R+1];

	lambda_in_c[ONE_SIDE] = lambda_in_n;
	lambda_out_c[ONE_SIDE] = lambda_out_n;
	AA_c[ONE_SIDE] = AA_n;
	
	double lambda_in_star1, AA_2n;

	// At each time step, lambda_in_n and AA_2n are known (from the incoming pipe
	// characteristic and the constant back pressure value, respectively)

	// Set default flow directions
	// ----------------------------------------------------------------------------------------------------
	pipe_flow[ONE_SIDE] = OUTFLOW;
	REVERSE_FLOW = false;
	CHOKED[ONE_SIDE] = false;

	// Enter inital characteristic values
	// ----------------------------------------------------------------------------------------------------
	lambda_in_star1 = lambda_in_n/AA_n;
	AA_2n = pow( (Pb/pPpt->PREF), (pPpt->gammaAir(pBN[ONE_SIDE]->T)-1)/(2*pPpt->gammaAir(pBN[ONE_SIDE]->T)) );

	// Conventional flow (interpolate from the processed data then calculate equivalent area)
	// ----------------------------------------------------------------------------------------------------
	{
		pipe_flow[ONE_SIDE] = OUTFLOW;
		
		double speed_param;
		//if(INLET_STAG) speed_param = NT/sqrt(T01);
		//else speed_param = NT/sqrt(T1);
		speed_param = NT/sqrt(T01);			
		Interpolate(pPpt, speed_param, lambda_in_star1);
//interp_data_pt.PrintToScreen(pPpt, VARIABLE);
		
		// Calculation of equivalent area
		double p01_over_p2;
		double T1_temp, k;
		double M, f_of_M, F_N;
		
//		if(INLET_STAG) // True - map uses p01 (and T01)
//		{
			p01_over_p2 = interp_data_pt.pr; // p2==p02==constant pressure for this type of turbine
//		}
//		else // False - map uses p1 (and T1) so interp_data_pt.pr = p1/p2 = p1/Pb here
//		{
			//double p1_temp = interp_data_pt.pr*Pb;
			//p1_temp = pow((lambda_in1 + lambda_out1)/(2*AA_1), (2*pPpt->gammaAir())/(pPpt->gammaAir()-1))*pPpt->PREF;
			//double lambda_out1 = lambda_out_n;	// Previous result
			T1_temp = pow((lambda_in_n + lambda_out_n)/2, 2)*(EX ? pPpt->TREFe : pPpt->TREFi);
			k = pPpt->gammaAir(T1_temp);
			//T01_temp= (pow((lambda_in_n + lambda_out_n)/2, 2) + ((pPpt->gammaAir(T1_temp)-1)/2)*pow((lambda_in_n - lambda_out_n)/(pPpt->gammaAir(T1_temp)-1), 2))*((EX ? pPpt->TREFe : pPpt->TREFi));
			//double T01_temp, rho1_temp, u1_temp, p01_temp, m_dot_temp;
			//rho1_temp = p1_temp/(pPpt->R_air*T1_temp);
			//u1_temp = (lambda_in_n - lambda_out_n)/(pPpt->gammaAir(T1_temp) - 1)*(EX ? pPpt->AREFe : pPpt->AREFi);
			//p01_temp = p1_temp + 0.5*rho1_temp*pow(u1_temp,2);
			//p01_over_p2 = p01_temp/this->Pb;
//		}
//		p01_over_p2 = interp_data_pt.pr;

//cout << "p01_over_p2 = " << p01_over_p2 << endl;
//cout << "T1_temp = " << T1_temp << endl;

		M = pow( (2/(pPpt->gammaAir(T1_temp)-1))*( pow( p01_over_p2, (pPpt->gammaAir(T1_temp)-1)/pPpt->gammaAir(T1_temp) ) - 1 ) , 0.5);

		if(M < 1.0 && M >= 0.0) 
			f_of_M = pow(pPpt->gammaAir(T1_temp)/pPpt->R_air, 0.5) 
			* (M/( pow(1 + ((pPpt->gammaAir(T1_temp)-1)/2)*pow(M,2), (pPpt->gammaAir(T1_temp)+1)/(2*(pPpt->gammaAir(T1_temp)-1))) ));
		else
		{
			if(M > 1.0)
				f_of_M = pow(pPpt->gammaAir(T1_temp)/pPpt->R_air, 0.5)*pow(2/(pPpt->gammaAir(T1_temp)+1), (pPpt->gammaAir(T1_temp)+1)/(2*(pPpt->gammaAir(T1_temp)-1)));
			else
			{
				cout << "Turbine[" << this->ID << "]:CTurbine::EquivArea Mach number is negative, M = " << M << "\n";
				char pause; cin >> pause;
			}
		}
		
/*
		if(INLET_STAG) // Can use MFP directly
			F_N = interp_data_pt.mfp*(1/f_of_M); // Equivalent area
		else
		{
			// MFP is then m_dot*sqrt(T1)/p1
			m_dot_temp = interp_data_pt.mfp*p1_temp/sqrt(T1_temp);
			F_N = (m_dot_temp*sqrt(T01_temp)/p01_temp) * (1/f_of_M); // Equivalent area
		}
*/
		F_N = interp_data_pt.mfp*(1/f_of_M); // Equivalent area
		double eq_phi = F_N/pBN[ONE_SIDE]->f_dash*pPpt->fref;
		
		// Test for choked flow
		// ----------------------------------------------------------------------------------------------------
		if(lambda_in_star1 - interp_data_pt.lambda_in_star_max > 1e-6)
		{
			CHOKED[ONE_SIDE] = true;
			REVERSE_FLOW = false;
			pipe_flow[ONE_SIDE] = OUTFLOW;				
			lambda_out_c[ONE_SIDE] = interp_data_pt.Ks*lambda_in_n;
		}
		else
		{
			// Test for reverse flow (treat as open inflow end)
			// ----------------------------------------------------------------------------------------------------
			if(lambda_in_star1<interp_data_pt.lambda_in_star_min)
			{
				// INFLOW
				CHOKED[ONE_SIDE] = false;
				REVERSE_FLOW = true;
				
				pipe_flow[ONE_SIDE] = INFLOW;
				
				bool SONIC;
				common_NHI_code(pPpt, lambda_in_n, lambda_out_n, AA_n, 
					lambda_in_c[ONE_SIDE], lambda_out_c[ONE_SIDE], AA_c[ONE_SIDE], 
					(pPpt->USE_PHI ? eq_phi : 1.0), Pb, T0, 
					CHOKED[ONE_SIDE], SONIC, T1_temp, timestep, time, true/*true=constant pressure (valve) model or false=pressure loss (port) model*/, 1e-6/*inflow main loop tolerance*/);
			}
			else
			{
				// Equivalent area method - call common nozzle fuction
				// ----------------------------------------------------------------------------------------------------
				lambda_out_c[ONE_SIDE] = common_NHN_code(pPpt, lambda_in_n, AA_n, eq_phi, Pb, CHOKED[ONE_SIDE], T1_temp, time);
			}
		}
	}
	common_UPDATE_NH(pPpt, lambda_in_c, lambda_out_c, AA_c, pipe_flow, CHOKED);
}

void CTurbine::Interpolate(CProperties* pPpt, double speed, double query_temp)
// ============================================================ //
// Interpolates the turbine map to set a new CDataPt			//
// for the current iteration from which the required			//
// properties (e.g. lambda_out_star) can be extracted.			//
// This function is called from									//
// CTurbine::VariablePressure									//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Interpolate\n");}

	query = query_temp;

	int curve, any_point, point_curve_above, point_curve_below;
	int GRAPH = 0;
	
	// Find the right curve
	curve = 0;
	any_point = 0;

	// Initial assumptions
	ON_MAP = true;
	ON_CURVE = true;
	CHOKED_FLOW = false;

	bool STOP;

	if (SINGLE_SPEED) {
		// Single speed line only variables
		int point_curve;

		// Find the two points on the single speed line that lambda_in_star1 lies between and interpolate
		curve = 0; curve_below = curve; curve_above = curve;

		point_curve = 0;
		ON_CURVE = true; // Assume this to be true at the start

		STOP = false;
		while (Data[curve][point_curve].query_param < query && !STOP) {
			if ((point_curve + 1) >= num_extra_reverse + num_points[GRAPH][curve] + num_extra_choked)
			{
				ON_CURVE = false;
				STOP = true;
			}
			else ++point_curve;
		}
		
		if (point_curve == 0) {
			ON_CURVE = false;
			point_curve_below = 0;
			point_curve_above = 1;
		}
		else { // lambda_in_star1 lies between point_curve and point_curve-1
			point_curve_above = point_curve;
			point_curve_below = point_curve - 1;
		}

		// Interpolate these points
		if (fabs(Data[curve][point_curve_above].query_param - Data[curve][point_curve_below].query_param) < 1e-6) {
			interp_data_pt = Data[curve][point_curve_below];
		}
		else {
			interp_data_pt = Data[curve][point_curve_below] +
				((Data[curve][point_curve_above] - Data[curve][point_curve_below])
					*
					((query - Data[curve][point_curve_below].query_param)
						/ (Data[curve][point_curve_above].query_param - Data[curve][point_curve_below].query_param)));
		}
	}
	else {
		// Multiple speed line interpolation variables
		int point_curve_above_above, point_curve_above_below, point_curve_below_above, point_curve_below_below;
		CDataPoint data_point_curve_above, data_point_curve_below;
		
		ON_MAP = true;  // Assume this to be true at the start
		STOP = false;
		
		while(Data[curve][any_point].sp < speed && !STOP) {
			if( (curve + 1) >= num_curves[GRAPH]) {
				STOP = true;
				ON_MAP = false;	// Since current speed is greater than that of the highest speed curve
			}
			else ++curve;	
		}

		if (curve == 0) {
			ON_MAP = false; // Current speed is lower than that of the lowest speed line
			curve_below = 0;
			curve_above = 1;
		}
		else {
			// The speed required is between speed lines [curve_above] and [curve_below]
			curve_above = curve;
			curve_below = curve - 1;
		}
		
		// Next find the two points on each of the curves that lambda_in_star1 lies between and interpolate each
		
		// 'curve above'
		point_curve_above = 0;
		ON_CURVE = true; // Assume this to be true at the start
		STOP = false;

		while(Data[curve_above][point_curve_above].query_param < query && !STOP) {
			if ((point_curve_above + 1) >= num_extra_reverse + num_points[GRAPH][curve_above] + num_extra_choked)
			{
				ON_CURVE = false;
				STOP = true;
			}
			else ++point_curve_above;
		}
		
		if (point_curve_above==0) {
			ON_CURVE = false;
			point_curve_above_below = 0;
			point_curve_above_above = 1;
		}
		else { // On this speed line lambda_in_star1 lies between point_curve_above and point_curve_above-1
			point_curve_above_above = point_curve_above;
			point_curve_above_below = point_curve_above - 1;
		}

		// Interpolate between the points on the speed line [curve above]

		// If operation lies below the first point on the speed line, just use the first point to prevent extrapolation
		if (query <= Data[curve_above][point_curve_above_below].query_param) data_point_curve_above = Data[curve_above][point_curve_above_below];
		else {
			if (fabs(Data[curve_above][point_curve_above_above].query_param - Data[curve_above][point_curve_above_below].query_param) < 1e-6)
			{
				data_point_curve_above = Data[curve_above][point_curve_above_below];
			}
			else {
				data_point_curve_above = Data[curve_above][point_curve_above_below] +
					((Data[curve_above][point_curve_above_above] - Data[curve_above][point_curve_above_below])
						*
						((query - Data[curve_above][point_curve_above_below].query_param)
							/ (Data[curve_above][point_curve_above_above].query_param - Data[curve_above][point_curve_above_below].query_param)));
			}
		}

		// 'curve below'
		point_curve_below = 0;
		ON_CURVE = true; // Assume this to be true at the start
		STOP = false;

		while(Data[curve_below][point_curve_below].query_param < query && !STOP) {
			if ((point_curve_below + 1) >= num_extra_reverse + num_points[GRAPH][curve_below] + num_extra_choked)
			{
				ON_CURVE = false;
				STOP = true;
			}
			else ++point_curve_below;
		}

		if(point_curve_below==0) {
			ON_CURVE = false;
			point_curve_below_below = 0;
			point_curve_below_above = 1;
		}
		else { // On this speed line lambda_in_star1 lies between point_curve_below and point_curve_below-1
			point_curve_below_above = point_curve_below;
			point_curve_below_below = point_curve_below - 1;
		}

		// Interpolate between the points on the speed line [curve below]

		// If operation lies below the first point on the speed line, just use the first point to prevent extrapolation
		if (query <= Data[curve_below][point_curve_below_below].query_param) data_point_curve_below = Data[curve_below][point_curve_below_below];
		else {
			if (fabs(Data[curve_below][point_curve_below_above].query_param - Data[curve_below][point_curve_below_below].query_param) < 1e-6)
			{
				data_point_curve_below = Data[curve_below][point_curve_below_below];
			}
			else
			{
				data_point_curve_below = Data[curve_below][point_curve_below_below] +
					((Data[curve_below][point_curve_below_above] - Data[curve_below][point_curve_below_below])
						*
						((query - Data[curve_below][point_curve_below_below].query_param)
							/ (Data[curve_below][point_curve_below_above].query_param - Data[curve_below][point_curve_below_below].query_param)));
			}
		}

		// Now interpolate between speed lines

		any_point = 0; // All points on one curve have the same speed
		if(fabs(Data[curve_above][any_point].sp - Data[curve_below][any_point].sp) < 1e-6) interp_data_pt = data_point_curve_below;
		else {
			interp_data_pt = data_point_curve_below +
				((data_point_curve_above - data_point_curve_below)
				*
				((speed - Data[curve_below][any_point].sp)
					/(Data[curve_above][any_point].sp - Data[curve_below][any_point].sp)) );
		}	

		if (VARIABLE) { // Calculate interpolated x_max based on speed
			double interp_x_max
				= x_max[curve_below]
				+
				((speed - Data[curve_below][any_point].sp)
					*
					((x_max[curve_above] - x_max[curve_below])
						/ (Data[curve_above][any_point].sp - Data[curve_below][any_point].sp)));

			// Check for special conditions
			if (query > interp_x_max) CHOKED_FLOW = true;
		}
	}
	
/*
	cout << "Interpolated data point:" << endl;
	interp_data_pt.PrintToScreen(pPpt, VARIABLE);
	cout << endl << endl;
//exit(1);
//*/
}

void CTurbine::InstantaneousWorkAndMassFlow(CProperties* pPpt, double del_t)
// ============================================================ //
// Calculates instantaneous work/power and mass flow rate		//
// for a variable outlet pressure turbine.						//
// This function is called from CTurbine::RunBoundary			//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".InstantaneousWorkAndMassFlow\n");}

	// Dimension characteristics and entropy level
	double lambda_in [2];
	double lambda_out [2];
	double AA [2];
	
	int INLET = 0;
	int OUTLET = 1;

	// Enter characteristic values (common)
	lambda_in[INLET] = (*(pCLIN[INLET]))[1];
	//lambda_out[INLET] = (*(pCLOUT[INLET]))[1];
	lambda_out[INLET] = (*(pCLOUT[INLET]))[0];
	AA[INLET] = (pBN[INLET])->AA[1];

	if(VARIABLE) // Varible outlet pressure; need to take account of variable outlet conditions
	{
		double p2_over_pref, p0_over_pref;
		double A0, A2;
	
		lambda_in[OUTLET] = (*(pCLIN[OUTLET]))[1];
		lambda_out[OUTLET] = (*(pCLOUT[OUTLET]))[1];
		AA[OUTLET] = (pBN[OUTLET])->AA[1];

		// Calculations
		A0 = sqrt( pow((lambda_in[INLET] + lambda_out[INLET])/2, 2) 
				+ ((pPpt->gammaAir(this->T1)-1)/2)*pow((lambda_in[INLET] - lambda_out[INLET])/(pPpt->gammaAir(this->T1)-1), 2) );
		p2_over_pref = pow( (lambda_in[OUTLET] + lambda_out[OUTLET])/(2*AA[OUTLET]), (2*pPpt->gammaAir(this->T1))/(pPpt->gammaAir(this->T1)-1) );
		p0_over_pref = pow(A0/AA[INLET], 2/(pPpt->gammaAir(this->T1)-1));
		A2 = A0*pow(p2_over_pref/p0_over_pref, (pPpt->gammaAir(this->T1)-1)/(2*pPpt->gammaAir(this->T1)));

		// Instantaneous work_per_unit_mass
		eA_over_m_dot = (pow(A0,2)/(pPpt->gammaAir(this->T1)-1))*(1 - pow(A2/A0, 2))*pow((EX ? pPpt->AREFe : pPpt->AREFi),2);
	}
	else // Constant outlet pressure
	{
		double A1, A2;	// Inlet, constant outlet
		double U1;		// Inlet n.d. velocity
		double AA_common;

		A1 = (lambda_in[INLET] + lambda_out[INLET])/2;
		U1 = (lambda_in[INLET] - lambda_out[INLET])/(pPpt->gammaAir(this->T1)-1);

		// A2 (turbine exit) is related to the constant exhaust pressure, P2 = Pb, and the exit entropy level
		AA_common = AA[INLET];	// Assuming isentropic flow AA_2 = AA_1 = AA
		A2 = AA_common*pow(Pb/pPpt->PREF, (pPpt->gammaAir(this->T1)-1)/(2*pPpt->gammaAir(this->T1)));

/*
cout << "A1= " << A1 << "\n";
cout << "U1= " << U1 << "\n";
cout << "AA_common= " << AA_common << "\n";
cout << "A2= " << A2 << "\n";
cout << "(pow(A1,2) - pow(A2,2))= " << (pow(A1,2) - pow(A2,2)) << "\n";
cout << "eA_over_m_dot= " << eA_over_m_dot << "\n";
*/

		// Instantaneous work_per_unit_mass
		if(REVERSE_FLOW)
		{
			if(fabs(A1 - A2) < 1e-6) eA_over_m_dot = 0;
			else eA_over_m_dot = pow((EX ? pPpt->AREFe : pPpt->AREFi),2)*( (1/(pPpt->gammaAir(this->T1)-1))*(pow(A1,2) - pow(A2,2)) + pow(U1,2)/2);
			
			// This case is treated as an open end (which is not correct since the turbine will
			// act as a restriction), so the m_dot will be too high, and slow the speed down too much.
			// Thus ignore the effect of reverse flow on this power, and rely on the bearing power
			// to slow the shaft down.

			// Better to treat as a nozzle using the equivalent area at the current speed?

			// Or better to set up reverse flow charcteristics like in the variable case?
		}
		else
		{
			if(fabs(A1 - A2) < 1e-6) eA_over_m_dot = pow((EX ? pPpt->AREFe : pPpt->AREFi),2)*( (1/(pPpt->gammaAir(this->T1)-1))*(0) + pow(U1,2)/2);
			else eA_over_m_dot = pow((EX ? pPpt->AREFe : pPpt->AREFi),2)*( (1/(pPpt->gammaAir(this->T1)-1))*(pow(A1,2) - pow(A2,2)) + pow(U1,2)/2);
		}
	}

	// Common (between variable and constant outlet pressure)
	// ======================================================

	// Calculate mass flow rate based on turbine inlet conditions
	m_dot = ((4*pPpt->gammaAir(this->T1))/(pPpt->gammaAir(this->T1)-1))
					*( (pPpt->PREF*1e5)*pBN[INLET]->f/(EX ? pPpt->AREFe : pPpt->AREFi))  // NB *1e5 since PREF in bar
					*((lambda_in[INLET] - lambda_out[INLET])/pow(lambda_in[INLET] + lambda_out[INLET], 2))
					*pow((lambda_in[INLET] + lambda_out[INLET])/(2*AA[INLET]), (2*pPpt->gammaAir(this->T1))/(pPpt->gammaAir(this->T1)-1));

	// Calculate T01/T1 after convergence (i.e. here) or during convergence???? 
	// Instantaneous stagnation temperature T01 ahead of the turbine
	T01 = (pow((lambda_in[INLET] + lambda_out[INLET])/2, 2) + ((pPpt->gammaAir(this->T1)-1)/2)*pow((lambda_in[INLET] - lambda_out[INLET])/(pPpt->gammaAir(this->T1)-1), 2))*((EX ? pPpt->TREFe : pPpt->TREFi));
	
	// Instantaneous static temperature T ahead of the turbine
	T1 = pow((lambda_in[INLET] + lambda_out[INLET])/2, 2)*((EX ? pPpt->TREFe : pPpt->TREFi));

	// Instantaneous power
	eA = m_dot*eA_over_m_dot;
	//eA = fabs(m_dot)*eA_over_m_dot;
	eta_TS = interp_data_pt.eta;
	eT = eA/eta_TS;
	
	// Update running totals for the current cycle
	work_per_cycle_Th += eT*del_t;
	work_per_cycle_A += eA*del_t;
	mass_flow_per_cycle += m_dot*del_t;

	// Instantaneous matching
	W_TI = eA;
//	if(fabs(NT)>1e-16)
	if(NT!=0) L_TI = W_TI/(2*PI*NT); //L_TI = W_TI/(2*PI*fabs(NT)); 
	else {
		if(fabs(W_TI)>1e-16) L_TI = 1e-6; //-1e-6; // "Stick" torque
		else L_TI = 0;
	}
}

void CTurbine::CycleWorkAndMassFlow(CProperties* pPpt, double Nc)
// ============================================================ //
// Calculates the turbine work and mass flow					//
// once per engine cycle.										//
// This function is called from the main function when the		//
// relevant engine signals a new cycle.							//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".CycleWorkAndMassFlow\n");}

	// Note Nc is number of cycles per second

	// Cycle/total values (calculated once per cycle)
	W_Th = Nc*work_per_cycle_Th;
	W_TA = Nc*work_per_cycle_A;
	eta_T = W_TA/W_Th;
	m_dot_T = Nc*mass_flow_per_cycle;

	// Reset running totals for the next cycle
	work_per_cycle_Th = 0;
	work_per_cycle_A = 0;
	mass_flow_per_cycle = 0;
}

void CTurbine::InstantaneousTransmissionLoss(CProperties* pPpt, double del_t)
// ============================================================ //
// Calculates instantaneous transmission power loss through		//
// the bearings, W_B.											//
// This function is called from CTurbine::RunBoundary			//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".InstantaneousTransmissionLoss\n");}

	double mu = 0.08; // A function of speed; zero at zero speed
	mu = NT*0.08; // ??
	double P = 100;
	double dm = 0.01;

	L_B = 0.5e-3 * mu * P * dm; // Should be positive
	W_B = L_B*(2*PI*fabs(NT));	// Same sign as torque
}

void CTurbine::Matching(CProperties* pPpt, double del_t)
// ============================================================ //
// Calculates the rotational acceleration of the turbocharger	//
// rotor assembly and the new shaft speed every per step,		//
// following the turbine and compressor calculations.			//
// This function is called from CTurbocharger::Runboundary		//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Matching\n");}

	// Instantaneous matching	
	if(W_TI!=0) eta_M = (W_TI - W_B)/W_TI; // eq. 9.175
//	if(W_TI!=0) eta_M = (W_TI - W_B)/fabs(W_TI); // eq. 9.175
	else eta_M = 1; // ??
	// This could be negative if the instantaneous tranmission power loss is greater than the instantaneous turbine power

	double L_CI = 0;
	if (FIXED_SPEED) L_CI = L_TI - L_B;			// If a fixed speed simulation, absorb the turbine torque so there is no acceleration
	dNTdt = (1/(2*PI*I))*(L_TI - L_B - L_CI);	// Rotational acceleration eq. 9.180
	// L_B always opposes rotation, so the value of L_B should always be positive, or could use fabs here
	
	if(!FIXED_SPEED) NT += dNTdt*del_t;			// Update rotational speed
}

// Read/write functions
// ====================================================================================================
void CTurbine::ReadInput(CProperties* pPpt, char *InputFile)
// ============================================================ //
// Reads turbine parameters from the relevant file				//
// and assigns the values to the appropriate member variables.	//
// This function is called from CTurbine::Initialise			//
// ============================================================ //
{
	//if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ReadInput\n");} // Can't call this - Identify() not valid yet...

	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, this->labels, this->values, this->strings, last_entry);
	for(int r=0; r<last_entry+1; ++r)
	{
		// ====================================================================================================
		// Parameter file for Turbine [0]
		// ====================================================================================================
		
		// Optional object description (max. 500 characters - use underscores for spaces)
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "strDesc") == 0) strDesc = strings[r];
/*
		// Type
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "RADIAL") == 0) RADIAL = DoubleToBool(values[r]);
		
			// If this is an axial flow turbine, i.e., RADIAL == 0 == false
			// ----------------------------------------------------------------------------------------------------	
			if(strcmp(labels[r], "EQUIVALENT_AREA") == 0) EQUIVALENT_AREA = DoubleToBool(values[r]);

			// Else this is a radial turbine, i.e., RADIAL == 1 == true
			// ----------------------------------------------------------------------------------------------------
			if(strcmp(labels[r], "SINGLE_SPEED") == 0) SINGLE_SPEED = DoubleToBool(values[r]);
		
		if(strcmp(labels[r], "FIXED_SPEED") == 0) FIXED_SPEED = DoubleToBool(values[r]);

			// If speed is fixed, i.e., FIXED_SPEED == 1 == true
			// ----------------------------------------------------------------------------------------------------	
			if(strcmp(labels[r], "fixedSpeed") == 0) fixedSpeed = values[r];
*/
		// Type
		// ----------------------------------------------------------------------------------------------------
		if (strcmp(labels[r], "EQUIVALENT_AREA") == 0) EQUIVALENT_AREA = DoubleToBool(values[r]);
		if (strcmp(labels[r], "FIXED_SPEED") == 0) FIXED_SPEED = DoubleToBool(values[r]);

			// If speed is fixed, i.e., FIXED_SPEED == 1 == true
			// ----------------------------------------------------------------------------------------------------	
			if (strcmp(labels[r], "fixedSpeed") == 0) fixedSpeed = values[r];

		if(strcmp(labels[r], "Pb") == 0) Pb = values[r];
		if(strcmp(labels[r], "T0") == 0) T0 = values[r];
		if(strcmp(labels[r], "I") == 0) I = values[r];
			
		// Map
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "MAP_FILE") == 0) MAP_FILE = strings[r];
		if(strcmp(labels[r], "SINGLE_SPEED") == 0) SINGLE_SPEED = DoubleToBool(values[r]);
		if(strcmp(labels[r], "tolSpeed") == 0) tolSpeed = values[r];
		if(strcmp(labels[r], "ECHO_MAP") == 0) ECHO_MAP = DoubleToBool(values[r]);
		if(strcmp(labels[r], "PRINT_PROC_MAP") == 0) PRINT_PROC_MAP = DoubleToBool(values[r]);

		// Measurements
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "USE_DEF_FREQ") == 0) USE_DEF_FREQ = DoubleToBool(values[r]);
		if(strcmp(labels[r], "freq") == 0) freq = int(values[r]);
		if(strcmp(labels[r], "PRINT_DEBUG_FILE") == 0) PRINT_DEBUG_FILE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "PRINT_MOVIE_FILE") == 0) PRINT_MOVIE_FILE = DoubleToBool(values[r]);

		// Deprecated
		// ----------------------------------------------------------------------------------------------------
		//if(strcmp(labels[r], "INLET_STAG") == 0) INLET_STAG = DoubleToBool(values[r]);
		//if(strcmp(labels[r], "OUTLET_STAG") == 0) OUTLET_STAG = DoubleToBool(values[r]);
		//if(strcmp(labels[r], "MFP_P") == 0) MFP_P = int(values[r]);

		// ====================================================================================================
		// End of file
		// ====================================================================================================
	}
	// Set some derived parameters
	if(USE_DEF_FREQ) freq = pPpt->freq;		// Use the default sampling rate
}

void CTurbine::ListProperties(CProperties* pPpt)
// ============================================================ //
// Lists turbine properties on screen.							//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ListProperties\n");}

	// ====================================================================================================
	// Parameter file for Turbine [0]
	// ====================================================================================================

	pPpt->Out(Underline(Identify(), "=", "\t", strDesc));
	pPpt->Out("\n");
/*
	// Type
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out(Underline("Type", "-", "\t"));
	if(RADIAL) {
		pPpt->Out("\tRadial flow turbine, RADIAL\t\t\t=\t"); pPpt->Out(TrueOrFalse(RADIAL)); pPpt->Out("\n");
		if(SINGLE_SPEED){pPpt->Out("\tConstant speed running, SINGLE_SPEED\t\t=\t"); pPpt->Out(TrueOrFalse(SINGLE_SPEED)); pPpt->Out("\n");}
	}
	else {
		pPpt->Out("\tAxial flow turbine, RADIAL\t\t\t=\t"); pPpt->Out(TrueOrFalse(RADIAL)); pPpt->Out("\n");
		if(EQUIVALENT_AREA){pPpt->Out("\tEquivalent area method, EQUIVALENT_AREA\t\t=\t"); pPpt->Out(TrueOrFalse(EQUIVALENT_AREA)); pPpt->Out("\n");}
		else{pPpt->Out("\tSimple unique curve, EQUIVALENT_AREA\t\t=\t"); pPpt->Out(TrueOrFalse(EQUIVALENT_AREA)); pPpt->Out("\n");}
	}
*/

	// Type
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out(Underline("Type", "-", "\t"));
	if (EQUIVALENT_AREA) { pPpt->Out("\tEquivalent area method, EQUIVALENT_AREA\t\t=\t"); pPpt->Out(TrueOrFalse(EQUIVALENT_AREA)); pPpt->Out("\n"); }
	if(FIXED_SPEED) {
		pPpt->Out("\tFix rotor speed, FIXED_SPEED\t\t\t=\t"); pPpt->Out(TrueOrFalse(FIXED_SPEED)); pPpt->Out("\n");
		
			// If speed is fixed, i.e., FIXED_SPEED == 1 == true
			// ----------------------------------------------------------------------------------------------------	
			pPpt->Out("\tFixed rotor speed, fixedSpeed\t\t\t=\t"); pPpt->Out(fixedSpeed); pPpt->Out(" r/min\n");
	}
	else {
		pPpt->Out("\tFree rotor speed, FIXED_SPEED\t\t\t=\t"); pPpt->Out(TrueOrFalse(FIXED_SPEED)); pPpt->Out("\n");
	}

	if(VARIABLE){pPpt->Out("\tVariable outlet pressure, VARIABLE\t\t=\t"); pPpt->Out(TrueOrFalse(VARIABLE)); pPpt->Out("\n");}
	else {
		pPpt->Out("\tConstant outlet pressure, VARIABLE\t\t=\t"); pPpt->Out(TrueOrFalse(VARIABLE)); pPpt->Out("\n");
		pPpt->Out("\tOutlet pressure, Pb\t\t\t\t=\t"); pPpt->Out(Pb); pPpt->Out(" bar\n");
		pPpt->Out("\tOutlet temperature, T0\t\t\t\t=\t"); pPpt->Out(T0); pPpt->Out(" K\n");
	}	
	pPpt->Out("\tTurbine wheel inertia, I\t\t\t=\t"); pPpt->Out(I); pPpt->Out(" kg.m^2\n");
	pPpt->Out("\n");

	// Map
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out(Underline("Map", "-", "\t"));
	pPpt->Out("\tTurbine map filename, MAP_FILE\t\t\t=\t"); pPpt->Out(MAP_FILE); pPpt->Out("\n");
	pPpt->Out("\tNomenclature, strNomenclature\t\t\t=\t"); pPpt->Out(strNomenclature); //pPpt->Out("\n");
	if (SINGLE_SPEED) { pPpt->Out("\Single speed line, SINGLE_SPEED\t\t=\t"); pPpt->Out(TrueOrFalse(SINGLE_SPEED)); pPpt->Out("\n"); }
	pPpt->Out("\tVariation in speed on same line, tolSpeed\t=\t"); pPpt->Out(tolSpeed); pPpt->Out("%\n");
	
	if(ECHO_MAP)
	{
		pPpt->Out("\n");
		pPpt->Out("\tPrinting SAE map to screen, ECHO_MAP\t\t=\t"); pPpt->Out(TrueOrFalse(ECHO_MAP)); pPpt->Out("\n");
		pPpt->Out("\n");
		pPpt->Out("\tN/T0^0.5 [rpm/K^0.5]\tMFP [kg/s K^0.5/kPa]\tPR (t-s)\teta (t-s)\n");
		pPpt->Out("\t--------------------\t--------------------\t--------\t---------\n");
		for(int speed=0; speed<speeds; ++speed)
		{
			for(int point=0; point<numPointsSAE[speed]; ++point)
			{
				pPpt->Out("\t");
				pPpt->Out(SAE[speed][point][0]*60);			// Convert back to SAE format (r/s to r/min)
				pPpt->Out("\t\t\t"); 
					pPpt->Out(SAE[speed][point][1]*1000);	// Convert back to SAE format (1/Pa to 1/kPa) 
					pPpt->Out("\t\t\t");
						pPpt->Out(SAE[speed][point][2]); 
						pPpt->Out("\t\t");
							pPpt->Out(SAE[speed][point][3]); 
							pPpt->Out("\n");
			}
		}
		pPpt->Out("\n");
	}
	else{pPpt->Out("\tNot printing SAE map to screen, ECHO_MAP\t=\t"); pPpt->Out(TrueOrFalse(ECHO_MAP)); pPpt->Out("\n");}
	
	if(PRINT_PROC_MAP)
	{
		pPpt->Out("\n");
		pPpt->Out("\tPrint processed map to screen, PRINT_PROC_MAP\t=\t"); pPpt->Out(TrueOrFalse(PRINT_PROC_MAP)); pPpt->Out("\n");

		int POINT;
		int temp_precision = cout.precision(); cout << setprecision(6);
		int GRAPH = 0;//TYPE - 1;
		//int CURVE = 0;
		for(int CURVE=0; CURVE<num_curves[GRAPH]; ++CURVE)
 		{
			string str_curve = "Curve [" + IntToString(CURVE) + "], " + IntToString(CURVE+1) + " of " + IntToString(num_curves[GRAPH]);
			pPpt->Out("\n");
			pPpt->Out(Underline(StringToChar(str_curve), "-", "\t"));	
			pPpt->Out("\n");	
			pPpt->Out("\t"); pPpt->Out(sp_str); pPpt->Out("\t"); pPpt->Out(mfp_str); pPpt->Out("\t"); pPpt->Out(pr_str); pPpt->Out("\t"); pPpt->Out(eta_str); pPpt->Out("\t"); pPpt->Out("power_param"); pPpt->Out("\t"); pPpt->Out("torque_param"); pPpt->Out("\tG1\t\tT2_over_T1\tp1_over_p2\tlambda_in_star_ratio\tquery_param\n");
			pPpt->Out("\t----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			for(POINT=0; POINT<(num_extra_reverse+0); ++POINT)
			{
				pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].sp); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].mfp); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].pr); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].eta); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].power_param); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].torque_param); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].G1); pPpt->Out("       \t");
				pPpt->Out(Data[CURVE][POINT].T2_over_T1); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].p1_over_p2); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].lambda_in_star_ratio); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].query_param); pPpt->Out("\n");
			}
			pPpt->Out("\t----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			for(POINT=(num_extra_reverse+0); POINT<(num_extra_reverse + num_points[GRAPH][CURVE]); ++POINT)
			{
				pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].sp); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].mfp); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].pr); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].eta); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].power_param); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].torque_param); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].G1); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].T2_over_T1); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].p1_over_p2); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].lambda_in_star_ratio); pPpt->Out("\t\t\t");
				pPpt->Out(Data[CURVE][POINT].query_param); pPpt->Out("\n");
			}
			pPpt->Out("\t----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			for(POINT=(num_extra_reverse + num_points[GRAPH][CURVE]); POINT<(num_extra_reverse + num_points[GRAPH][CURVE] + num_extra_choked); ++POINT)
			{
				pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].sp); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].mfp); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].pr); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].eta); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].power_param); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].torque_param); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].G1); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].T2_over_T1); pPpt->Out("\t");
				pPpt->Out(Data[CURVE][POINT].p1_over_p2); pPpt->Out("\t\t");
				pPpt->Out(Data[CURVE][POINT].lambda_in_star_ratio); pPpt->Out("\t\t\t");
				pPpt->Out(Data[CURVE][POINT].query_param); pPpt->Out("\n");
			}
			pPpt->Out("\t----------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
			if(VARIABLE)
			{
				pPpt->Out("\n");
	//			pPpt->Out("\tnum_extra_reverse = "); pPpt->Out(num_extra_reverse); pPpt->Out("\n");
	//			pPpt->Out("\tnum_extra_choked = "); pPpt->Out(num_extra_choked); pPpt->Out("\n");
				pPpt->Out("\tx_min\t=\t"); pPpt->Out(x_min[CURVE]); pPpt->Out("\n");
				pPpt->Out("\tPt\t=\t"); pPpt->Out(Pt[CURVE]); pPpt->Out("\n");
				pPpt->Out("\ttheta_t\t=\t"); pPpt->Out(theta_t[CURVE]); pPpt->Out("\n");
				pPpt->Out("\tG1_max\t=\t"); pPpt->Out(G1_max[CURVE]); pPpt->Out("\n");
				pPpt->Out("\tK1\t=\t"); pPpt->Out(K1[CURVE]); pPpt->Out("\n");
				pPpt->Out("\tK2\t=\t"); pPpt->Out(K2[CURVE]); pPpt->Out("\n");
				pPpt->Out("\tK3\t=\t"); pPpt->Out(K3[CURVE]); pPpt->Out("\n");
				pPpt->Out("\tK4\t=\t"); pPpt->Out(K4[CURVE]); pPpt->Out("\n");
				pPpt->Out("\tx_max\t=\t"); pPpt->Out(x_max[CURVE]); pPpt->Out("\n");
	//			pPpt->Out("\n");
			}
			pPpt->Out("\n");
		}
		if(VARIABLE)
		{
			//pPpt->Out("\n");
			pPpt->Out(Underline("All curves", "-", "\t"));
			pPpt->Out("\tnum_extra_reverse = "); pPpt->Out(num_extra_reverse); pPpt->Out("\n");
			pPpt->Out("\tnum_extra_choked = "); pPpt->Out(num_extra_choked); pPpt->Out("\n");
	//		pPpt->Out("\tx_min = "); pPpt->Out(x_min); pPpt->Out("\n");
	//		pPpt->Out("\tPt = "); pPpt->Out(Pt); pPpt->Out("\n");
	//		pPpt->Out("\ttheta_t = "); pPpt->Out(theta_t); pPpt->Out("\n");
	//		pPpt->Out("\tG1_max = "); pPpt->Out(G1_max); pPpt->Out("\n");
	//		pPpt->Out("\tK1 = "); pPpt->Out(K1); pPpt->Out("\n");
	//		pPpt->Out("\tK2 = "); pPpt->Out(K2); pPpt->Out("\n");
	//		pPpt->Out("\tK3 = "); pPpt->Out(K3); pPpt->Out("\n");
	//		pPpt->Out("\tK4 = "); pPpt->Out(K4); pPpt->Out("\n");
	//		pPpt->Out("\tx_max = "); pPpt->Out(x_max); pPpt->Out("\n");
			pPpt->Out("\n");
		}
		cout << setprecision(temp_precision);
	}
	else{pPpt->Out("\tPrint processed map to screen, PRINT_PROC_MAP\t=\t"); pPpt->Out(TrueOrFalse(PRINT_PROC_MAP)); pPpt->Out("\n\n");}
	 
	// Measurements
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out(Underline("Measurements", "-", "\t"));
	if(USE_DEF_FREQ)
	{
		if(freq==1) pPpt->Out("\tUsing default sampling rate, pPpt->freq\t=\tonce per timestep\n");
		else{pPpt->Out("\tUsing default sampling rate, pPpt->freq)\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");}
	}
	else
	{
		if(freq==1) pPpt->Out("\tUsing local sampling rate, freq\t\t\t=\tonce per timestep\n");
		else{pPpt->Out("\tUsing local sampling rate, freq\t\t\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");}
	}
	if(PRINT_DEBUG_FILE) pPpt->Out("\tPrinting debug file\n"); else pPpt->Out("\tNot printing debug file\n");
	if(PRINT_MOVIE_FILE) pPpt->Out("\tPrinting movie file\n"); else pPpt->Out("\tNot printing movie file\n");
	pPpt->Out("\n");
	pPpt->Out("\n");

	// ====================================================================================================
	// End of file
	// ====================================================================================================
}

void CTurbine::PrintToScreen(CProperties* pPpt)
// ============================================================ //
// Prints instantaneous and cycle turbine data to screen.		//
// This function is called from the main function.				//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintToScreen\n");}
	//cout << setprecision(8);
	pPpt->Out(Underline("Turbine", "=", ID));
//	interp_data_pt.PrintToScreen(pPpt, true);
	pPpt->Out("\n");

	//pPpt->Out(Underline("Dynamics", "-"));
	pPpt->Out("Rotor speed, NT\t\t\t\t\t\t=\t"); pPpt->Out(NT*60); pPpt->Out(" rpm ("); pPpt->Out(NT); pPpt->Out(" rps)\n");
	pPpt->Out("Acceleration, dNTdt\t\t\t\t\t=\t"); pPpt->Out(dNTdt); pPpt->Out(" rev/s^2"); if(FIXED_SPEED) pPpt->Out(" (fixed speed simulation)"); pPpt->Out("\n");
	pPpt->Out("Inlet total temperature, T01\t\t\t\t=\t"); pPpt->Out(T01); pPpt->Out(" K\n");
	pPpt->Out("Implied speed parameter, NT/sqrt(T01)\t\t\t=\t"); pPpt->Out((NT*60)/sqrt(T01)); pPpt->Out(" rpm/K^0.5 ("); pPpt->Out(NT/sqrt(T01)); pPpt->Out(" rps/K^0.5)\n");
	pPpt->Out("\n");

	pPpt->Out(Underline("Operating point (interpolated)", "-"));
	if(ON_MAP) pPpt->Out("On map"); else pPpt->Out("Off map"); pPpt->Out(", ");
	if(ON_CURVE) pPpt->Out("on curve"); else pPpt->Out("off curves"); pPpt->Out(", ");
	if(REVERSE_FLOW) pPpt->Out("reverse flow"); else pPpt->Out("forward flow");
	if(CHOKED_FLOW) pPpt->Out(", choked flow");
	pPpt->Out("\n");
	pPpt->Out("Speed parameter, "); pPpt->Out(sp_str); pPpt->Out("\t\t\t\t=\t"); pPpt->Out(interp_data_pt.sp*60); pPpt->Out(" rpm/K^0.5 ("); pPpt->Out(interp_data_pt.sp); pPpt->Out(" rps/K^0.5)\n");
	//cout << "Turbine-EquivArea-2" << endl;
	//cout << "curve_below = " << curve_below << endl;
	//cout << "curve_above = " << curve_above << endl;
	//exit(1);
	pPpt->Out("- hence interpolating between curves ["); pPpt->Out(curve_below); pPpt->Out("] ("); pPpt->Out(Data[curve_below][0].sp*60); pPpt->Out(") and ["); pPpt->Out(curve_above); pPpt->Out("] ("); pPpt->Out(Data[curve_above][0].sp*60); pPpt->Out(")\n");
	pPpt->Out("Mass flow parameter, "); pPpt->Out(mfp_str); pPpt->Out("\t\t\t=\t"); pPpt->Out(interp_data_pt.mfp); pPpt->Out(" "); pPpt->Out(mfp_units_str); pPpt->Out("\n");
	pPpt->Out("Pressure ratio, "); pPpt->Out(pr_str); pPpt->Out("\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.pr); pPpt->Out(" "); pPpt->Out("\n");
	pPpt->Out("Efficiency, "); pPpt->Out(eta_str); pPpt->Out("\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.eta*100); pPpt->Out("%\n");
	pPpt->Out("\n");
	
	if(VARIABLE) {
		pPpt->Out("Variable outlet pressure turbine parameters (interpolated point):\n");
		pPpt->Out("G1\t\t\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.G1); pPpt->Out("\n"); //pPpt->Out(" ?\n");
		pPpt->Out("p1_over_p2\t\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.p1_over_p2); pPpt->Out("\n");
		pPpt->Out("T2_over_T1\t\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.T2_over_T1); pPpt->Out("\n");
		pPpt->Out("lambda_in_star_ratio\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.lambda_in_star_ratio); pPpt->Out("\n");
		pPpt->Out("query_param\t\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.query_param); pPpt->Out("\n");
	}
	else {
		pPpt->Out("Constant outlet pressure turbine parameters (interpolated point):\n");
		pPpt->Out("lambda_in_star\t\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.lambda_in_star); pPpt->Out("\n");
		pPpt->Out("query_param\t\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.query_param); pPpt->Out("\n");
		pPpt->Out("lambda_out_star\t\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.lambda_out_star); pPpt->Out("\n");
		pPpt->Out("lambda_in_star_min\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.lambda_in_star_min); pPpt->Out("\n");
		pPpt->Out("lambda_in_star_max\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.lambda_in_star_max); pPpt->Out("\n");
		pPpt->Out("Ks\t\t\t\t\t\t\t=\t"); pPpt->Out(interp_data_pt.Ks); pPpt->Out("\n");
	}
	pPpt->Out("\n");
	pPpt->Out("Most recent completed cycle data:\n");
	pPpt->Out("Cycle theoretical turbine power, W_Th\t\t\t=\t"); pPpt->Out(W_Th/1000); pPpt->Out(" kW\n");
	pPpt->Out("Cycle actual turbine power, W_TA\t\t\t=\t"); pPpt->Out(W_TA/1000); pPpt->Out(" kW\n");
	pPpt->Out("Cycle turbine average efficiciency, eta_T\t\t=\t"); pPpt->Out(eta_T*100); pPpt->Out("%\n");
	pPpt->Out("Cycle turbine mass flow rate, m_dot_T\t\t\t=\t"); pPpt->Out(m_dot_T); pPpt->Out(" kg.s^-1\n");
	pPpt->Out("\n");
	pPpt->Out("Current instantaneous data:\n");
	pPpt->Out("Instantaneous theoretical turbine power, eT\t\t=\t"); pPpt->Out(eT); pPpt->Out(" W\n");
	pPpt->Out("Instantaneous actual turbine power, eA\t\t\t=\t"); pPpt->Out(eA); pPpt->Out(" W\n");
	pPpt->Out("Instantaneous turbine efficiency, eta_TS\t\t=\t"); pPpt->Out(eta_TS*100); pPpt->Out("%\n");
	pPpt->Out("Instantaneous turbine mass flow rate, m_dot\t\t=\t"); pPpt->Out(m_dot); pPpt->Out(" kg.s^-1\n");
	pPpt->Out("\n");
	pPpt->Out("Instantaneous work per unit mass, eA_over_m_dot\t\t=\t"); pPpt->Out(eA_over_m_dot); pPpt->Out(" J.kg^-1\n");
	pPpt->Out("Instantaneous turbine power, W_TI\t\t\t=\t"); pPpt->Out(W_TI); pPpt->Out(" W\n");
	pPpt->Out("Instantaneous turbine torque, L_TI\t\t\t=\t"); pPpt->Out(L_TI); pPpt->Out(" N.m\n");
	pPpt->Out("Instantaneous turbine inlet stag. temp., T01\t\t=\t"); pPpt->Out(T01); pPpt->Out(" K\n");
	pPpt->Out("Instantaneous turbine speed parameter, NT/T01^0.5\t=\t"); pPpt->Out(NT*60/sqrt(T01)); pPpt->Out(" rpm/K^0.5 ("); pPpt->Out(NT/sqrt(T01)); pPpt->Out(" rps/K^0.5)\n");
	pPpt->Out("\n");
	pPpt->Out("Instantaneous bearing power, W_B\t\t\t=\t"); pPpt->Out(W_B); pPpt->Out(" W\n");
	pPpt->Out("Instantaneous bearing torque, L_B\t\t\t=\t"); pPpt->Out(L_B); pPpt->Out(" N.m\n");
	pPpt->Out("Instantaneous mechanical efficiency, eta_M\t\t=\t"); pPpt->Out(eta_M*100); pPpt->Out("%\n");
	pPpt->Out("\n");
	pPpt->Out("Turbine current running totals:\n");
	pPpt->Out("Theoretical work this cycle, work_per_cycle_Th\t\t=\t"); pPpt->Out(work_per_cycle_Th); pPpt->Out(" J\n");
	pPpt->Out("Actual work this cycle, work_per_cycle_A\t\t=\t"); pPpt->Out(work_per_cycle_A); pPpt->Out(" J\n");
	pPpt->Out("Mass flow this cycle, mass_flow_per_cycle\t\t=\t"); pPpt->Out(mass_flow_per_cycle); pPpt->Out(" kg.s^-1\n");
	pPpt->Out("\n");
}

void CTurbine::PrintToFile(CProperties* pPpt, int timestep, double time, double ca)
// ============================================================ //
// Prints instantaneous turbine data to file.					//
// This function is called from the main function.				//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintToFile\n");}

	if (timestep == 0) {

		std::string temp_str = "Object results file for ";
		temp_str += Identify();
		fprintf(OUTPUT_FILE, "%s\n", Underline(StringToChar(temp_str), "-"));

		fprintf(OUTPUT_FILE, "%s\t%s%c%s", "Time (s)", "ca (", Deg(), "CA)");

		// Print inlet boundary conditions to file
		fprintf(OUTPUT_FILE, "\t%s\t%s\t%s\t%s", "Inlet static pressure (bar)", "Inlet static temperature (K)", "Inlet density (kg/m^3)", "Inlet velocity (m/s)");

		// Print outlet boundary conditions to file, or stagnation values if constant pressure outlet
		if (VARIABLE) fprintf(OUTPUT_FILE, "\t%s\t%s\t%s\t%s", "Outlet static pressure (bar)", "Outlet static temperature (K)", "Outlet density (kg/m^3)", "Outlet velocity (m/s)");
		else fprintf(OUTPUT_FILE, "\t%s\t%s", "Constant outlet pressure, Pb (bar)", "Constant outlet temperature, T0 (K)");

		// Print turbine operational parameters to file
		fprintf(OUTPUT_FILE, "\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
			"eT (W)", "eA (W)", "m_dot (kg.s^-1)", "NT (s^-1)", "L_TI (N.m)", "PR", "MFP [kg.s^-1.K^0.5.kPa^-1]",
			"eta_TS", "query", "G1 x 1e3", "ON_MAP", "ON_CURVE", "REVERSE_FLOW", "CHOKED_FLOW");

		// Print lambda_in and lambda_out to file (inlet boundary)
		fprintf(OUTPUT_FILE, "\t%s\t%s", "lambda_in (turbine entry)", "lambda_out (turbine entry)");

		// Print lambda_in and lambda_out to file (outlet boundary, if present)
		if(VARIABLE) fprintf(OUTPUT_FILE, "\t%s\t%s", "lambda_in (turbine exit)", "lambda_out (turbine exit)");

		// Print lambda_in_star and lambda_out_star to file (inlet boundary)
		fprintf(OUTPUT_FILE, "\t%s\t%s", "lambda_in_star (turbine entry)", "lambda_out_star (turbine entry)");

		// Print lambda_in_star and lambda_out_star to file (outlet boundary, if present)
		if (VARIABLE) fprintf(OUTPUT_FILE, "\t%s\t%s", "lambda_in_star (turbine exit)", "lambda_out_star (turbine exit)");

		fprintf(OUTPUT_FILE, "\n");
	}

	if (timestep % freq == 0) { // Print data at the specified sampling frequency
		
		fprintf(OUTPUT_FILE, "%f\t%f", time, ca);

		// Print inlet boundary conditions to file
		fprintf(OUTPUT_FILE, "\t%f\t%f\t%f\t%f", pBN[0]->p_dash * pPpt->PREF, pBN[0]->T, pBN[0]->rho, pBN[0]->U * pPipe[0]->AREF);

		// Print outlet boundary conditions to file, or stagnation values if constant pressure outlet
		if (VARIABLE) fprintf(OUTPUT_FILE, "\t%f\t%f\t%f\t%f", pBN[0]->p_dash * pPpt->PREF, pBN[0]->T, pBN[0]->rho, pBN[0]->U * pPipe[0]->AREF);
		else fprintf(OUTPUT_FILE, "\t%f\t%f", this->Pb, this->T0);

		// Print turbine operational parameters to file
		fprintf(OUTPUT_FILE, "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%s\t%s\t%s",
			eT, eA, m_dot, NT, L_TI, interp_data_pt.pr, interp_data_pt.mfp * 1e3, // MFP: convert to per kPa from per Pa 
			eta_TS, query, interp_data_pt.G1 * 1e3, TrueOrFalse(ON_MAP), TrueOrFalse(ON_CURVE), TrueOrFalse(REVERSE_FLOW), TrueOrFalse(CHOKED_FLOW));

		// Print lambda_in and lambda_out to file (inlet boundary)
		fprintf(OUTPUT_FILE, "\t%f\t%f", (*(pCLIN[0]))[R + 1], (*(pCLOUT[0]))[R + 1]);

		// Print lambda_in and lambda_out to file (outlet boundary, if present)
		if (VARIABLE) fprintf(OUTPUT_FILE, "\t%f\t%f", (*(pCLIN[1]))[R + 1], (*(pCLOUT[1]))[R + 1]);

		// Print lambda_in_star and lambda_out_star to file (inlet boundary)
		fprintf(OUTPUT_FILE, "\t%f\t%f", CLIN_STAR[0][R + 1], CLOUT_STAR[0][R + 1]);

		// Print lambda_in_star and lambda_out_star to file (outlet boundary, if present)
		if (VARIABLE) fprintf(OUTPUT_FILE, "\t%f\t%f", CLIN_STAR[1][R + 1], CLOUT_STAR[1][R + 1]);

		fprintf(OUTPUT_FILE, "\n");
	}
}

void CTurbine::PrintToMovieFile(CProperties* pPpt, int timestep, double time, double ca)
// ============================================================ //
// Prints instantaneous turbine data to the movie file.			//
// This function is called from the main function.				//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintToMovieFile\n");}

	if(PRINT_MOVIE_FILE)
	{
		fprintf(MOVIE_FILE,"%.6f\t", time);
		fprintf(MOVIE_FILE,"%.6f\t", ca);
		fprintf(MOVIE_FILE,"%.6f\t", interp_data_pt.pr);
		fprintf(MOVIE_FILE,"%.6f\n", interp_data_pt.mfp*1e3); // MFP: convert to per kPa from per Pa
	}
}
