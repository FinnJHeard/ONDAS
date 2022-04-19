// Pipe.cpp: implementation of the CPipe class.
//
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Globals.h"
#include "Pipe.h"
#include "Tools.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CPipe::CPipe()
{

}

CPipe::~CPipe()
{
	delete [] Node;
	delete [] PathLine;
	delete [] Node_Backup;
	delete [] PathLine_Backup;
	delete [] Measurements;
	delete [] MeasureNode;
	delete [] Results;
}

void CPipe::ReadInput(CProperties* pPpt, char *InputFile)
{
	int last_entry;
	CommonReadInput(InputFile, pPpt->NUM_PARAMS, this->labels, this->values, this->strings, last_entry);

	// Default settings in case not specified
	delay = 0;
	DIAMETER = 0;
	AREA = 0;
	STATIC_PRESSURE = 0;
	TEMPERATURE = 0;
	DENSITY = 0;
	VELOCITY = 0;
	MACH_NUMBER = 0;
	MASS_FLOW_RATE = 0;
	REYNOLDS_NO = 0;

	for(int r=0; r<last_entry+1; ++r)
	{
		// Optional object description (max. 500 characters - use underscores for spaces)
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "strDesc") == 0) strDesc = strings[r];

		// Propagation method
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "METHOD") == 0) METHOD = int(values[r]);
		if(strcmp(labels[r], "delay") == 0) delay = values[r];
		
		// Geometry
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "length") == 0) length = values[r]/1000; // Conversion from mm to m
		if(strcmp(labels[r], "d_odd") == 0) d_odd = values[r]/1000;
		if(strcmp(labels[r], "d_even") == 0) d_even = values[r]/1000;
		if(strcmp(labels[r], "bend_angle") == 0) bend_angle = values[r];
		if(strcmp(labels[r], "n_int_points") == 0)
		{
			n_int_points = int(values[r]);
			xi = new double [n_int_points+2]; // Internal points plus both ends
			di = new double [n_int_points+2]; // Internal points plus both ends
			xi[0] = 0; di[0] = d_odd;
			xi[n_int_points+2-1] = length; di[n_int_points+2-1] = d_even;
		}
		if(strcmp(labels[r], "xi1") == 0 && n_int_points>0) xi[1] = values[r]/1000; if(strcmp(labels[r], "di1") == 0 && n_int_points>0) di[1] = values[r]/1000;
		if(strcmp(labels[r], "xi2") == 0 && n_int_points>1) xi[2] = values[r]/1000; if(strcmp(labels[r], "di2") == 0 && n_int_points>1) di[2] = values[r]/1000;
		if(strcmp(labels[r], "xi3") == 0 && n_int_points>2) xi[3] = values[r]/1000; if(strcmp(labels[r], "di3") == 0 && n_int_points>2) di[3] = values[r]/1000;
		if(strcmp(labels[r], "xi4") == 0 && n_int_points>3) xi[4] = values[r]/1000; if(strcmp(labels[r], "di4") == 0 && n_int_points>3) di[4] = values[r]/1000;
		if(strcmp(labels[r], "LINEAR") == 0) LINEAR = DoubleToBool(values[r]);
		if(strcmp(labels[r], "LINEAR_F") == 0) LINEAR_F = DoubleToBool(values[r]);
		if(strcmp(labels[r], "x_int") == 0) x_int = values[r]/1000;
		if(strcmp(labels[r], "d_int") == 0) d_int = values[r]/1000;

		// Meshing
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "GLOBAL_MP") == 0) GLOBAL_MP = DoubleToBool(values[r]);
		if(strcmp(labels[r], "discret") == 0) discret = values[r]/1000;
		if(strcmp(labels[r], "min_meshes") == 0) min_meshes = int(values[r]);
		
		// Coefficients
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "GLOBAL_EF") == 0) GLOBAL_EF = DoubleToBool(values[r]);
		if(strcmp(labels[r], "epsilon") == 0) epsilon = values[r]/1000;
		if(strcmp(labels[r], "CFTRANS") == 0) CFTRANS = values[r];
		if(strcmp(labels[r], "HGTRANS") == 0) HGTRANS = values[r];
		if(strcmp(labels[r], "CPTRANS") == 0) CPTRANS = values[r];
	
		// Initialization
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "Tw") == 0) Tw = values[r];
		if(strcmp(labels[r], "nsections") == 0)
		{
			nsections = int(values[r]);
			
			// Dimension intialisation parameters using nsections
			xri = new double [nsections];
			pri = new double [nsections];
			Tri = new double [nsections];
			vri = new double [nsections];
			
			xri[nsections-1] = 1.0; // Always add an imaginary division at the end of the pipe
		}

		// Following allows for up to 5 sections == 4 divisions
		if(strcmp(labels[r], "xri0") == 0 && nsections>1) xri[0] = values[r];
		if(strcmp(labels[r], "xri1") == 0 && nsections>2) xri[1] = values[r];
		if(strcmp(labels[r], "xri2") == 0 && nsections>3) xri[2] = values[r];
		if(strcmp(labels[r], "xri3") == 0 && nsections>4) xri[3] = values[r];

		if(strcmp(labels[r], "pri0") == 0) pri[0] = values[r];
		if(strcmp(labels[r], "Tri0") == 0) Tri[0] = values[r];
		if(strcmp(labels[r], "vri0") == 0) vri[0] = values[r];

		if(strcmp(labels[r], "pri1") == 0 && nsections>1) pri[1] = values[r];
		if(strcmp(labels[r], "Tri1") == 0 && nsections>1) Tri[1] = values[r];
		if(strcmp(labels[r], "vri1") == 0 && nsections>1) vri[1] = values[r];

		if(strcmp(labels[r], "pri2") == 0 && nsections>2) pri[2] = values[r];
		if(strcmp(labels[r], "Tri2") == 0 && nsections>2) Tri[2] = values[r];
		if(strcmp(labels[r], "vri2") == 0 && nsections>2) vri[2] = values[r];

		if(strcmp(labels[r], "pri3") == 0 && nsections>3) pri[3] = values[r];
		if(strcmp(labels[r], "Tri3") == 0 && nsections>3) Tri[3] = values[r];
		if(strcmp(labels[r], "vri3") == 0 && nsections>3) vri[3] = values[r];

		if(strcmp(labels[r], "pri4") == 0 && nsections>4) pri[4] = values[r];
		if(strcmp(labels[r], "Tri4") == 0 && nsections>4) Tri[4] = values[r];
		if(strcmp(labels[r], "vri4") == 0 && nsections>4) vri[4] = values[r];
		
		// Pipe measurements
		// --------------------------------------------------------------------------------
		if(strcmp(labels[r], "USE_DEF_FREQ") == 0) USE_DEF_FREQ = DoubleToBool(values[r]);
		if(strcmp(labels[r], "freq") == 0) freq = int(values[r]);	
		if(strcmp(labels[r], "print_from_time") == 0) print_from_time = values[r];
		
		if(strcmp(labels[r], "ntappings") == 0)
		{
			ntappings = int(values[r]);
			loc_measure = new double [ntappings];
		}
		
		if(strcmp(labels[r], "loc_measure0") == 0 && ntappings>0) loc_measure[0] = values[r];
		if(strcmp(labels[r], "loc_measure1") == 0 && ntappings>1) loc_measure[1] = values[r];
		if(strcmp(labels[r], "loc_measure2") == 0 && ntappings>2) loc_measure[2] = values[r];
		if(strcmp(labels[r], "loc_measure3") == 0 && ntappings>3) loc_measure[3] = values[r];
		if(strcmp(labels[r], "loc_measure4") == 0 && ntappings>4) loc_measure[4] = values[r];
		if(strcmp(labels[r], "loc_measure5") == 0 && ntappings>5) loc_measure[5] = values[r];
		if(strcmp(labels[r], "loc_measure6") == 0 && ntappings>6) loc_measure[6] = values[r];
		if(strcmp(labels[r], "loc_measure7") == 0 && ntappings>7) loc_measure[7] = values[r];
		if(strcmp(labels[r], "loc_measure8") == 0 && ntappings>8) loc_measure[8] = values[r];
		if(strcmp(labels[r], "loc_measure9") == 0 && ntappings>9) loc_measure[9] = values[r];
		if(strcmp(labels[r], "loc_measure10") == 0 && ntappings>10) loc_measure[10] = values[r];

		if(strcmp(labels[r], "max_pts") == 0) max_pts = int(values[r]);

		// Measurements to record
//		if(strcmp(labels[r], "num_props_measured") == 0) num_props_measured = int(values[r]);
		if(strcmp(labels[r], "DIAMETER") == 0) DIAMETER = DoubleToBool(values[r]);
		if(strcmp(labels[r], "AREA") == 0) AREA = DoubleToBool(values[r]);
		if(strcmp(labels[r], "STATIC_PRESSURE") == 0) STATIC_PRESSURE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "TEMPERATURE") == 0) TEMPERATURE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "DENSITY") == 0) DENSITY = DoubleToBool(values[r]);
		if(strcmp(labels[r], "VELOCITY") == 0) VELOCITY = DoubleToBool(values[r]);
		if(strcmp(labels[r], "MACH_NUMBER") == 0) MACH_NUMBER = DoubleToBool(values[r]);
		if(strcmp(labels[r], "MASS_FLOW_RATE") == 0) MASS_FLOW_RATE = DoubleToBool(values[r]);
		if(strcmp(labels[r], "REYNOLDS_NO") == 0) REYNOLDS_NO = DoubleToBool(values[r]);

		// Screen output
		// ----------------------------------------------------------------------------------------------------
		if(strcmp(labels[r], "SHOW_DATA") == 0) SHOW_DATA = DoubleToBool(values[r]);
	}

	// Checks
	// ----------------------------------------------------------------------------------------------------
	if(METHOD==0) METHOD=pPpt->DEF_METHOD;	// Use the default method for this pipe
	else{
		if(METHOD!=1) pPpt->HOMENTROPIC = false; // If not using MMOC must also be non-homentropic
		if(METHOD==1 && pPpt->HOMENTROPIC && pPpt->SHOW_pathlines) pPpt->SHOW_pathlines = false; // Pathlines not present in homentropic MMOC
		if(METHOD==2 && pPpt->TVD && pPpt->courant>0.7){
			SimulationWarning(pPpt, "A maximum Courant number of 0.7 will be enforced when using TVD criteria");
			pPpt->courant=0.7;
		}
	}

	// Meshing
	// ----------------------------------------------------------------------------------------------------
	if(GLOBAL_MP){ // Apply global meshing parameters
		discret = pPpt->discret;
		min_meshes = pPpt->min_meshes;
	}
	
	// Coefficients
	// ----------------------------------------------------------------------------------------------------
	if(GLOBAL_EF){ // Apply global enhancement factors
		epsilon = pPpt->epsilon;
		CFTRANS = pPpt->CFTRANS;
		HGTRANS = pPpt->HGTRANS;
		CPTRANS = pPpt->CPTRANS;
	}

	// Measurements
	// ----------------------------------------------------------------------------------------------------
	if(USE_DEF_FREQ) freq = pPpt->freq;		// Use the default sampling rate for this pipe
	num_props_measured = int(DIAMETER) + int(AREA) + int(STATIC_PRESSURE) + int(TEMPERATURE) + int(DENSITY) + int(VELOCITY) + int(MACH_NUMBER) + int(MASS_FLOW_RATE) + int(REYNOLDS_NO);
}

void CPipe::Initialise(CProperties* pPpt, int assyid, int id, bool ex, CEngine* EngPtr, std::string param_dir, string parent_assy_res_dir, string calling_object_str)
{	
	if (pPpt->SHOW_calls) { pPpt->Out("CPipe.Initialise ("); pPpt->Out(calling_object_str); pPpt->Out(")\n"); }

	int m, S;

	BUFFER = false;
	DAMPER = false;

	// Identification
	// ----------------------------------------------------------------------------------------------------
	AssyID = assyid;
	ID = id;
	if(ex) EX = true; else EX = false;	
	pEng = EngPtr;
	std::string obname_str = "PIPE";
	//cout << ConstructString(pPpt, param_dir, obname_str, EX, ID) << endl;	//	char pause; cin >> pause;
	
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Initialise\n");}

	ReadInput(pPpt, ConstructString(pPpt, param_dir, obname_str, EX, ID));

	// Tappings
	// ----------------------------------------------------------------------------------------------------
	Measurements = new double [ntappings];
	MeasureNode = new CNode [ntappings];
	for(m=0; m<ntappings; ++m) Measurements[m] = loc_measure[m];

	Results = new double ** [ntappings];
	for(m=0; m<ntappings; ++m)
	{ 
		Results[m] = new double * [max_pts];
		for(int d=0; d<max_pts; ++d) Results[m][d] = new double [num_props_measured+1+1];
	}
	rowcounter=0;
	sample_factor=1;

	Pressure		= new double [ntappings];
	Temperature		= new double [ntappings];
	Velocity		= new double [ntappings];
	MassFlowRate	= new double [ntappings];
	Re				= new double [ntappings];
	PR				= new double [ntappings];
	MFP				= new double [ntappings];

	for(m=0; m<ntappings; ++m)
	{
		Pressure[m] = 1;
		Temperature[m] = 300;
		Velocity[m] = 0;
		MassFlowRate[m] = 0;
		Re[m] = 0;
		PR[m] = 1;
		MFP[m] = 0;
	}

	SetupFiles(pPpt, parent_assy_res_dir);// Setup and open results file
	
	// Mesh generation
	// ----------------------------------------------------------------------------------------------------					
	odd_end_flow = NOFLOW; even_end_flow = NOFLOW;
	MMOC=1; W_ALPHA_BETA=2; FandE = 3; JOINER=4;	// Numerate METHOD labels

	// Set end corrections
	end_corr_odd	= end_corr_odd_p*(d_odd);
	end_corr_even	= end_corr_even_p*(d_even);
	eff_length		= length + (end_corr_odd + end_corr_even);

	// Determine correct number of meshes based on the desired discretisation length
	if((eff_length/discret) - int(eff_length/discret)>1e-6) meshes = int(eff_length/discret) + 1;
	else													meshes = int(eff_length/discret);

	// Enforce a minimum number of meshes
	if(meshes<min_meshes)									meshes = min_meshes; 

	if(METHOD==3) meshes = 1;	// For filling and emptying, only makes sense to have one mesh (i.e. one volume)
	else
	{	
		if(METHOD==4)
		{
			end_corr_odd	= 0;
			end_corr_even	= 0;
			length			= 0;
			eff_length		= 0;
			meshes			= 0;
		}
	}

	XPIPE				= length/pPpt->xref;	// N.D. pipe physical length ()
	N					= meshes + 1;			// For a pipe with seperate ends (even if pipe is controlled by a loop)
	Node				= new CNode[N];			// Dimension Node vector
	Node_Backup			= new CNode[N];			// Dimension backup Node vector

	// Pipe bend loss
	// --------------
	if(bend_angle > 0)
	{
		double bend_radius, d_bar, C_pA, C_pB, C_pB1, C_pB2;
		
		bend_radius = eff_length/(bend_angle*PI/180); // Pipe radius based on pipe length and bend angle (converted to radians)
		d_bar = (d_odd + d_even)/2; // Average of pipe end diameters

		if(bend_radius/d_bar <= 1) C_pA = 0.21*pow(bend_radius/d_bar, -2.5);
		else C_pA = 0.21*pow(bend_radius/d_bar, -0.25);

		C_pB1 = 0.9*sin(bend_angle*PI/180);
		C_pB2 = 0.7 + 0.35*sin(bend_angle*PI/180)/90;

		if(bend_angle < 70) C_pB = C_pB1;
		else
		{
			if(bend_angle <= 100) C_pB = C_pB1*((100 - bend_angle)/30) + C_pB2*((bend_angle - 70)/30);
			else C_pB = C_pB2;
		}
		C_p = C_pA*C_pB;
	}
	else C_p = 0;


	// Applicable to all pipes
	// -----------------------
	if(meshes>0) xmesh	= eff_length/meshes;	// Mesh effective length (m)
	else		 xmesh	= eff_length;			// Special case for a single node
	XMESH				= xmesh/pPpt->xref;		// N.D. mesh effective length ()
	XPIPE_EFF			= eff_length/pPpt->xref;// N.D. pipe effective length ()
	C					= (d_even - d_odd)/XPIPE_EFF; // C = dD/dX required for ONEoD_dDdX in n-h MOC 
	for(S=0; S<N; ++S)
	{
		Node[S].bc = INTERIOR;
		Node[S].side = INSIDE;
	}
	Node[0].side = ODD;
	Node[N-1].side = EVEN;

	if(pPpt->NUM_PATH_MULT*(2*N - 1) > 2*N - 1) num_pathlines = pPpt->NUM_PATH_MULT*(2*N - 1);
	else num_pathlines = 2*N - 1;

	if(!pPpt->HOMENTROPIC) // Pathlines exist only for non-homentropic cases
	{
		PathLine = new CPathLine[num_pathlines];
		PathLine_Backup = new CPathLine[num_pathlines];
	}

	if(METHOD!=4)
	{
		// Need to define x (m), X () for all nodes
		for(S=0; S<N; ++S)
		{
			Node[S].x = -end_corr_odd + S*xmesh;
			Node[S].X = Node[S].x/pPpt->xref;
		}
/*
		// Pipe gradient calculation - don't need to consider end corrections as we want the physical gradient of the pipe
		if(LINEAR && LINEAR_F)
		{
			double f_odd = PI*pow(d_odd,2)/4;
			double f_even = PI*pow(d_even,2)/4;
			double dfdx_temp = (f_even - f_odd)/this->length;
			double x = this->length;
			this->vol = (dfdx_temp/2)*pow(x,2) + f_odd*x;
		}
		else
		{		
			// Calculate quadratic variation of pipe diameter given 3 points: odd, even and interior
			if(LINEAR) // If LINEAR is selected just set interior point accordingly
			{
				x_int = length/2;
				d_int = (d_odd + d_even)/2;
			}
			else
			{
				if(x_int>length || x_int<0)
				{
					cout << Identify() << ".Initialise: ";
					cout << "Interior point must lie between odd and even ends of pipe. Exiting." << endl;
					exit(1);
				}
			}

			double y1, y2, y3;
			y1 = d_odd;
			y2 = d_int;
			y3 = d_even;

			// dia = ax^2 + bx + c
			double x1, x2, x3;
			x1 = 0; 
			x2 = x_int; 
			x3 = this->length;
			this->a = ((y3-y2)/((x3-x2)*(x3-x1))) - ((y1-y2)/((x1-x2)*(x3-x1)));
			//this->b = ((y3-y2) + a*(pow(x2,2)-pow(x3,2)))/(x2-x3);
			this->b = ((y1-y2) + a*(pow(x2,2)-pow(x1,2)))/(x1-x2);
			//this->c = y1 - a*pow(x1,2) - b*x1;
			this->c = y2 - a*pow(x2,2) - b*x2;
			//this->c = y3 - a*pow(x3,2) - b*x3;

			double x = this->length;
			this->vol = (PI/4)*( (pow(a,2)/5)*pow(x,5) + (a*b/2)*pow(x,4) + ((2*a*c + pow(b,2))/3)*pow(x,3) + b*c*pow(x,2) + pow(c,2)*x );

			// Still required for MOC:
			// dia = AX^2 + BX + C
			{
			double X1, X2, X3;
			X1 = 0/pPpt->xref; 
			X2 = x_int/pPpt->xref; 
			X3 = this->length/pPpt->xref;
			this->A = ((y3-y2)/((X3-X2)*(X3-X1))) - ((y1-y2)/((X1-X2)*(X3-X1)));
			this->B = ((y1-y2) + A*(pow(X2,2)-pow(X1,2)))/(X1-X2);//this->B = ((y3-y2) + A*(pow(X2,2)-pow(X3,2)))/(X2-X3);
			this->C = y2 - A*pow(X2,2) - B*X2;//this->C = y1 - A*pow(X1,2) - B*X1;//this->C = y3 - A*pow(X3,2) - B*X3;
			}
//cout << "dia = " << a << "x^2 + " << b << "x + " << c << endl << endl;
		}
//*/
		for(S=0; S<N; ++S)
		{
			/*
			double x = Node[S].x;
			Node[S].d = d(x);
			Node[S].f = f(x);
			Node[S].f_dash = Node[S].f/pPpt->fref;
			
			Node[S].dfdx = dfdx(x);					// Use only in linear pipes		

			//if(S!=0 && S!=N-1) Node[S].dfdx = (dfdx(x) + dfdx(Node[S+1].x))/2;	
			//if(S!=0 && S!=N-1) Node[S].dfdx = (dfdx(Node[S-1].x) + dfdx(x) + dfdx(Node[S+1].x))/3;
			//if(S!=0 && S!=N-1) Node[S].dfdx = dfdx(x + 0.00005*xmesh);

//			// More fundamental geometry variation (not necessary for linear pipes)
//			// For circular cross-section:
//			Node[S].dddx = dddx(pPpt, x);					// First derivative of diameter variation
//			Node[S].d2ddx2 = d2ddx2(pPpt, x);				// Second derivative of diameter variation
			*/



/*
			// SET LINEARLY VARYING F, i.e., SO dfdx = constant for the pipe
			Node[S].dfdx = ((PI/4)*(pow(d_even,2) - pow(d_odd,2))/length);
			Node[S].f = (PI/4)*pow(d_odd,2) + Node[S].dfdx*Node[S].x;
			Node[S].f_dash = Node[S].f/pPpt->fref;
			Node[S].d = sqrt((4/PI)*Node[S].f);
//*/
///*
			// For each nodal location interpolate from the list of internal diameters
			int i=1;
			while(Node[S].x > xi[i] && i<(n_int_points+2)-1) ++i;
			// Node[S] lies between xi[i] and xi[i-1]
			Node[S].d = ((di[i] - di[i-1])/(xi[i]- xi[i-1])) // Gradient
						*(Node[S].x - xi[i-1])				 // Difference	
						+ di[i-1];

//cout << "Node[S=" << S << "].d = " << Node[S].d << endl;

			Node[S].f = PI*pow(Node[S].d,2)/4;
			Node[S].f_dash = Node[S].f/pPpt->fref;
			Node[S].dddx = ((di[i] - di[i-1])/(xi[i]- xi[i-1])); // Diameter gradient - actual, correct dddx at this node (as long as node doesn't fall precisely at the same location as an internal diameter specification)	
///*
			// f = (pi/4)*d^2 so df/dd = (pi/2)*d
			// Then df/dx = dd/dx * df/dd = dd/dx * (pi/2)*d
			Node[S].dfdx = Node[S].dddx * ((PI / 2) * Node[S].d);
			//Node[S].dfdx = (((PI * pow(di[i], 2) / 4) - (PI * pow(di[i - 1], 2) / 4)) / (xi[i] - xi[i - 1])); // Area gradient - actual, correct dfdx at this node (as long as node doesn't fall precisely at the same location as an internal diameter specification)
//*/
/*			
			cout << "Node[S=" << S << "].x=" << Node[S].x*1000 << "mm, lies between xi[i=" << i << "]=" << xi[i]*1000 << "mm and xi[i-1=" << i-1 << "]=" << xi[i-1]*1000 << "mm" << endl;
			cout << "Node[S=" << S << "].d = (" << (di[i] - di[i-1])*1000 << ")/(" << (xi[i]- xi[i-1])*1000 << ")*" << (Node[S].x - xi[i-1])*1000 << " + " << di[i-1]*1000 << " = " << Node[S].d*1000 << endl;
			cout << endl;
//*/		
			// Haemodynamics: set undeformed diameter, first derivative of undeformed diameter, and undeformed area
			Node[S].d0 = Node[S].d;
			Node[S].dd0dx = Node[S].dddx;
			Node[S].f0 = Node[S].f;

			//Haemodynamics: calculate b and dbdr0
			Node[S].b = 4.0 / 3.0 * (pPpt->k1 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2.0) + pPpt->k3);
			Node[S].dbdr0 = 4.0 / 3.0 * (pPpt->k1 * pPpt->k2 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2.0));
		}

		for(S=0; S<N; ++S) // After setting the diameter and area at each node, must set dfdx in a separate loop 
		{
			//if(S==0 || pPpt->beta==0) Node[S].dfdx = (Node[S+1].f - Node[S].f)/(Node[S+1].x - Node[S].x);
			//else Node[S].dfdx = (1/pPpt->beta)*((Node[S].f - Node[S-1].f)/(Node[S].x - Node[S-1].x)) - ((1 - pPpt->beta)/pPpt->beta)*Node[S-1].dfdx;
			// This is the expression that must be used with general Wab schemes to generate dfdx in order to get zero flow in a stagnant pipe with area variation

			//Node[S].dfdx = (PI/2)*Node[S].d*((d_even - d_odd)/length);

/*
			if(S==0) Node[S].dfdx = (Node[S+1].f - Node[S].f)/(Node[S+1].x - Node[S].x);
			else {
				Node[S].dfdx = 2 * ((Node[S].f - Node[S - 1].f) / (Node[S].x - Node[S - 1].x)) - Node[S - 1].dfdx;
			}
//*/

/*
			if(pPpt->alpha==0.5 && pPpt->beta==0.5)
			{
				if(S==0) Node[S].dfdx = (Node[S+1].f - Node[S].f)/(Node[S+1].x - Node[S].x);
				else Node[S].dfdx = 2*((Node[S].f - Node[S-1].f)/(Node[S].x - Node[S-1].x)) - Node[S-1].dfdx;
			}
			else
			{
				if(pPpt->alpha==1 && pPpt->beta==0)
				{
					if(S==N-1) Node[S].dfdx = 0.5*(Node[S].f/(Node[S].x - Node[S-1].x));
					else Node[S].dfdx = 0.5*(Node[S+1].f/(Node[S+1].x - Node[S].x));
				}
				else
				{
					if(pPpt->alpha==1 && pPpt->beta==1)
					{
						if(S==0) Node[S].dfdx = (Node[S+1].f - Node[S].f)/(Node[S+1].x - Node[S].x);
						else
						{
							if(S==1) Node[S].dfdx = (Node[S+1].f - Node[S].f)/(Node[S+1].x - Node[S].x);
							if(S==N-1); // Need not do anything
							else Node[S+1].dfdx = (2*Node[S+1].f - Node[S-1].f)/(Node[S+1].x - Node[S].x) - Node[S-1].dfdx;
						}
					}
				}
			}
//*/
		}

/*
		for(S=0; S<N; ++S)
		{
			if(S==0)
			{
				//Node[S].dfdx = 0;
				Node[S].dfdx = (Node[S+1].f - Node[S].f)/(Node[S+1].x - Node[S].x);
			}
			else 
			{
				if(S==N-1)
				{
					//Node[S].dfdx = 0;
					Node[S].dfdx = (Node[S].f - Node[S-1].f)/(Node[S].x - Node[S-1].x);
				}
				else 
					//Node[S].dfdx = (Node[S+1].f - Node[S-1].f)/(Node[S+1].x - Node[S-1].x); // about same as function
					//Node[S].dfdx = (Node[S+1].f - Node[S].f)/(Node[S+1].x - Node[S].x); // Bad
					//Node[S].dfdx = (Node[S].f - Node[S-1].f)/(Node[S].x - Node[S-1].x); // Bad
					Node[S].dfdx = ((Node[S+1].f - Node[S].f)/(Node[S+1].x - Node[S].x) +
									(Node[S].f - Node[S-1].f)/(Node[S].x - Node[S-1].x))/2;
			}


			// Integrate diameter between limits
			//double dFdx = (PI/2)*(pow(b,2)*Node[S].x + b*c);
		}
//*/

		for(S=0; S<N; ++S)
		{
			// Cell frontal area for circular cross-section
//			double x_left_boundary, x_right_boundary;
/*		
			if(S==0) 
			{
				x_left_boundary = Node[S].x;
				x_right_boundary =  (Node[S+1].x + Node[S].x)/2;
			}
			else
			{
				if(S==N-1)
				{
					x_left_boundary = (Node[S].x + Node[S-1].x)/2;
					x_right_boundary = Node[S].x;
				}
				else
				{
					x_left_boundary = (Node[S].x + Node[S-1].x)/2;
					x_right_boundary =  (Node[S+1].x + Node[S].x)/2;
				}
			}
//*/
/*
			if(S==0) 
			{
				x_right_boundary =  (Node[S+1].x + Node[S].x)/2;
				x_left_boundary = Node[S].x - (x_right_boundary - Node[S].x);
			}
			else
			{
				if(S==N-1)
				{
					x_left_boundary = (Node[S].x + Node[S-1].x)/2;
					x_right_boundary = Node[S].x + (Node[S].x - x_left_boundary);
				}
				else
				{
					x_left_boundary = (Node[S].x + Node[S-1].x)/2;
					x_right_boundary = (Node[S+1].x + Node[S].x)/2;
				}
			}

			Node[S].cfa = (PI/4)*( pow(d(x_right_boundary),2) - pow(d(x_left_boundary),2) );
			//if(S<N-1) Node[S].cfa = (PI/4)*( pow(Node[S+1].d,2) - pow(Node[S].d,2) ); else Node[S].cfa = 0;
			Node[S].cfa_delx = x_right_boundary - x_left_boundary;
*/
/*
cout << "Node[S].d = " << Node[S].d << endl;
cout << "Node[S].dddx = " << Node[S].dddx << endl;
cout << "Node[S].d2ddx2 = " << Node[S].d2ddx2 << endl;
cout << "Node[S].f = " << Node[S].f << endl;
cout << "Node[S].cfa = " << Node[S].cfa << endl;
cout << "Node[S].cfa_delx = " << Node[S].cfa_delx << endl;
cout << "Node[S].dfdx = " << Node[S].dfdx << endl;
cout << endl;
//*/
		}
	}
	else
	{
		// METHOD = 4 (JOINER)

		// Need to define x (m), X () for all nodes
		for(S=0; S<N; ++S)
		{
			Node[S].x = -end_corr_odd + S*xmesh;
			Node[S].X = Node[S].x/pPpt->xref;
		}

//cout << "d_odd = " << d_odd << endl;
//cout << "d_even = " << d_even << endl;
//exit(1);

		double d_mean = 0.5*(d_odd + d_even);
		d_int = d_mean;

		for(S=0; S<N; ++S)
		{
			Node[S].d = d_mean;
			Node[S].f = (PI/4)*pow(d_mean,2);
			Node[S].f_dash = Node[S].f/pPpt->fref;
		}
	}

	// Check that x value of last node equals (physical length of pipe plus even end correction)
	if(fabs(Node[N-1].x - (length + end_corr_even)) > 1e-9)
	{
		cout << Identify() << ": geometry mismatch:" << endl;
		cout << "Node[N-1].x = " << Node[N-1].x << " m" << endl;
		cout << "(length + end_corr_even) = " << (length + end_corr_even) << " m" << endl;
		cout << "Press any key to exit.\n";
		char cont;
		cin >> cont;
		exit(1);
	}
	
	// Pipe volume - sum individual mesh volumes (which are truncated cones)
	// NB assumes circular cross section
	int mesh;
	vol = 0;
	for(mesh=0; mesh<N-1; ++mesh) vol += (PI/12)*(Node[mesh+1].x - Node[mesh].x)*(pow(Node[mesh].d + Node[mesh+1].d,2) - (Node[mesh].d*Node[mesh+1].d));

	// Prepare to set initial Riemann values
	int section;
	double k;
	for(S=0; S<N; ++S)
	{
		// Node spacing
		if(S<N-1) Node[S].DELX_R = Node[S+1].X - Node[S].X;	else Node[S].DELX_R = 0;		
		if(S>0) Node[S].DELX_L = Node[S].X - Node[S-1].X; else Node[S].DELX_L = 0;
		if(Node[S].DELX_R < 0) 
		{
			cout << Identify() << ": DELX_R negative!" << endl;
			cout << "end_corr_odd_p = " << end_corr_odd_p << endl;
			cout << "end_corr_even_p = " << end_corr_even_p << endl;
			cout<<"Press any key to exit.\n";
			char cont;
			cin >> cont;
			cout << endl;
			exit(1);
		}
		if(Node[S].DELX_L < 0) 
		{
			cout << Identify() << ": DELX_L negative!" << endl; 
			cout << "end_corr_odd_p = " << end_corr_odd_p << endl;
			cout << "end_corr_even_p = " << end_corr_even_p << endl;
			cout<<"Press any key to exit.\n";
			char cont;
			cin >> cont;
			cout << endl;
			exit(1);
		}

		// Determine which section this node is in
		section = 0;
		while(eff_length*xri[section] < Node[S].x && section < nsections-1)
		{
			//cout << "eff_length*xri[section=" << section << "] = " << eff_length*xri[section] << endl;
			//cout << "Node[S=" << S << "].x = " << Node[S].x << endl;
			++section;
		}
			
		if(pPpt->HOMENTROPIC)
		{		
			// For homentropic, need to set TREF in order to show correct TRI; see p271 Simple Homentropic Program
			// Cannot have multiple section temps under homentropic flow though due to there only being one TREF
			double TREF_temp;
			if(pPpt->CONTINUOUS)
			{
				k = pPpt->gammaAir(Tri[section]);
				TREF_temp = pow(pPpt->PREF/pri[section], (k-1)/k)*Tri[section];			
			}
			else
			{
				//	TREF_temp = pEng->TCR*pow(pPpt->PREF/pEng->PCR, (pPpt->gammaAir()-1)/pPpt->gammaAir());
				//	TREF_temp = pEng->TCR;
				TREF_temp = 300;
			}
			
			// Need to have one AREF for entire exhaust system, and another separate one for entire intake system
			if(EX)
			{
				pPpt->TREFe = TREF_temp;
				pPpt->AREFe = sqrt(pPpt->gammaAir(pPpt->TREFe)*pPpt->R_air*pPpt->TREFe);
				this->AREF = pPpt->AREFe;
			}
			else
			{
				pPpt->TREFi = TREF_temp;
				pPpt->AREFi = sqrt(pPpt->gammaAir(pPpt->TREFi)*pPpt->R_air*pPpt->TREFi);
				this->AREF = pPpt->AREFi;
			}
			/*
			cout << "TREF_temp = " << TREF_temp << endl;
			cout << "pPpt->TREFe = " << pPpt->TREFe << endl;
			cout << "pPpt->TREFi = " << pPpt->TREFi << endl;
			cout << "pPpt->AREFe = " << pPpt->AREFe << endl;
			cout << "pPpt->AREFi = " << pPpt->AREFi << endl;
			*/
		}
		else
		{	
			// Choose TREF as the initial temperature in the first pipe (first section) of that system and calculate AREF		
			if(EX)
			{
				if(ID==0 && AssyID==0) // If this is the first exhaust pipe in the first assembly
				{
					pPpt->TREFe = Tri[/*section=*/0];
					pPpt->AREFe = sqrt(pPpt->gammaAir(pPpt->TREFe)*pPpt->R_air*pPpt->TREFe);
				}
				this->AREF = pPpt->AREFe;
			}
			else
			{
				if(ID==0 && AssyID==0) // If this is the first intake pipe in the first assembly
				{
					pPpt->TREFi = Tri[/*section=*/0];
					pPpt->AREFi = sqrt(pPpt->gammaAir(pPpt->TREFi)*pPpt->R_air*pPpt->TREFi);
				}
				this->AREF = pPpt->AREFi;
			}
		}
	}

	if(!pPpt->HOMENTROPIC) // Requirements for non-homentropic flow
	{
		// Pathlines
		for(int k=0; k<num_pathlines; ++k)
		{
			// Just need to position the pathlines initially and provide an entropy level value
			// PathLines() will calculate the remaining member variables
			PathLine[k].ID = k;
			PathLine[k].XK = k*(XPIPE_EFF/double(num_pathlines)) + 0.5*(XPIPE_EFF/double(num_pathlines)); // if XK has dimension (m)
			PathLine[k].XK_old = PathLine[k].XK;
		}
	}

	// W_alpha_beta schemes
	// ====================
	for(S=0; S<N; ++S)
	{
		// Allocate space for solution, flux, source vectors
		Node[S].W = new double [3];
		Node[S].F = new double [3];
		Node[S].C = new double [3];

		Node[S].W_pred = new double [3];
		Node[S].F_pred = new double [3];
		Node[S].C_pred = new double [3];

		Node[S].S = new double [3];
		Node[S].S_pred = new double [3];

		Node[S].W_prev = new double [3];
	}

	// Filling and emptying
	// ====================
	FV = new CFiniteVolume[N-1]; // i.e. same as the number of meshes
	for(S=0; S<N-1; ++S)
	{
		FV[S].l = Node[S+1].x - Node[S].x;

		// Volume assuming linear variation of diameter between nodes
		if(Node[S+1].d==Node[S].d) FV[S].V = (PI/4)*pow(Node[S].d,2)*FV[S].l;
		else FV[S].V = (PI/12)*(FV[S].l/(Node[S+1].d - Node[S].d))*(pow(Node[S+1].d,3) - pow(Node[S].d,3));
	}

	this->InitialConditions(pPpt); // Set the initial gas flow conditions for this pipe

	A_throat = new double [2];
	U_throat = new double [2];
	A_throat_old = new double [2];
	U_throat_old = new double [2];
	A_error = new double [2];
	U_error = new double [2];
	A_error_old = new double [2];
	U_error_old = new double [2];
	CONVERGED = new bool [2];
	beta = new double [2];
	tol = new double [2];
	direction_str = new char* [2];
	ps_throat = new double [2];
	ps_throat_old = new double [2];
}

void CPipe::InitialConditions(CProperties* pPpt)
// Sets or resets initial gas flow conditions based on input file
{	
	double CLI, k;
	int section;
	for(int S=0; S<N; ++S)
	{		
		// Determine which section this node is in
		section = 0;
		while(eff_length*xri[section] < Node[S].x && section < nsections-1) ++section;

		if(pPpt->HOMENTROPIC) {
			CLI = pow(pri[section]/pPpt->PREF, pPpt->Q);
			// Sets pressure == PRI; see p271 Simple Homentropic Program

			Node[S].CL1[0] = CLI;
			Node[S].CL2[0] = CLI;
			Node[S].CL1[1] = Node[S].CL1[0];
			Node[S].CL2[1] = Node[S].CL2[0];
		}
		else {	
			k = pPpt->gammaAir(Tri[section]);
			// For variable initial velocity
			Node[S].CL1[0] = sqrt(Tri[section]/(EX ? pPpt->TREFe : pPpt->TREFi)) + ((k-1)/2)*(vri[section]/(EX ? pPpt->AREFe : pPpt->AREFi));
			Node[S].CL2[0] = sqrt(Tri[section]/(EX ? pPpt->TREFe : pPpt->TREFi)) - ((k-1)/2)*(vri[section]/(EX ? pPpt->AREFe : pPpt->AREFi));
			Node[S].CL1[1] = Node[S].CL1[0];
			Node[S].CL2[1] = Node[S].CL2[0];
		}

		// Generate initial A and U
//cout << "Node[S=" << S << "].CL1[0] = " << Node[S].CL1[0] << endl;
//cout << "Node[S=" << S << "].CL2[0] = " << Node[S].CL2[0] << endl;
		Node[S].A = (Node[S].CL1[0] + Node[S].CL2[0])/2;
		Node[S].T = pow(Node[S].A,2)*(EX ? pPpt->TREFe : pPpt->TREFi);
		k = pPpt->gammaAir(Node[S].T);
		Node[S].U = (Node[S].CL1[0] - Node[S].CL2[0])/(k-1);
		if(fabs(Node[S].U) < pPpt->ZERO_TOL) Node[S].U = 0;
        Node[S].M = fabs(Node[S].U/Node[S].A);
		Node[S].p_dash = pri[section]/pPpt->PREF;
//if(ID==1) cout << "pri[section] = " << pri[section] << endl;
//cout << "Node[S=" << S << "].p_dash = " << Node[S].p_dash << endl;
		Node[S].p0_dash = Node[S].p_dash;

		if (pPpt->HAEMODYNAMICS) {
			Node[S].rho = pPpt->rho_blood;
		}
		else Node[S].rho = ((Node[S].p_dash*pPpt->PREF)*1e5)/(pPpt->R_air*Node[S].T);
		
		Node[S].rho_prev = Node[S].rho;
		//Node[S].mdot = ( ( (Node[S].p_dash*pPpt->PREF*1e5)/(pPpt->R_air*Node[S].T) ) * ( fabs(Node[S].U)*AREF ) * Node[S].f );
		Node[S].mdot = Node[S].rho*(Node[S].U*AREF)*Node[S].f;
		Node[S].CHOKED = false;
	}

	if(!pPpt->HOMENTROPIC) // Requirements for non-homentropic flow
	{
		for(int S=0; S<N; ++S) // Initial entropy levels
		{
			k = pPpt->gammaAir(Node[S].T);
			Node[S].AA[0] = Node[S].A/pow(Node[S].p_dash, (k-1)/(2*k));
			Node[S].AA[1] = Node[S].AA[0];
		}

		// Pathlines
		for(int k=0; k<num_pathlines; ++k)
		{
			// Must select pathline AA from the surrounding nodes
			int S=0;
			while(Node[S].X<PathLine[k].XK && S+1!=N)++S;

			if(N>1) {
				if(S==0) S=1; // Else S-1 undefined
				if(S==N) S=N-1;	// Else S undefined
				
				PathLine[k].AA = ((PathLine[k].XK - Node[S-1].X)/(Node[S].X - Node[S-1].X))
					*(Node[S].AA[0] - Node[S-1].AA[0])
					+ Node[S-1].AA[0];
			}
			else { // For single node, join pipes
				S=0;
				PathLine[k].AA = Node[S].AA[0];
			}
		}
	}

	// W_alpha_beta schemes
	// ====================
	int S=0;

//	// Generate remaining primitives
//	// -----------------------------
//	for(S=0; S<N; ++S)
//	{		
//		Node[S].rho = ((Node[S].p_dash*pPpt->PREF)*1e5)/(pPpt->R_air*Node[S].T);
//		Node[S].rho_prev = Node[S].rho;
//	}

	// Generate initial solution vector
	// --------------------------------
	for(S=0; S<N; ++S) {	
		k = pPpt->gammaAir(Node[S].T);
		Node[S].W[0] = Node[S].rho*Node[S].f;
		Node[S].W[1] = Node[S].rho*Node[S].f*(Node[S].U*AREF);

		// Need haemodynamic version for:
		Node[S].W[2] = Node[S].rho * (((Node[S].p_dash*pPpt->PREF)*1e5)/(Node[S].rho*(k-1)) + 0.5*pow((Node[S].U*AREF), 2)) * Node[S].f; 
	}

	// Filling and emptying
	// ====================
	for(int S=0; S<N-1; ++S) {
		//FV[S].p0 = Node[S+1].p_dash*pPpt->PREF*1e5 + 0.5*Node[S+1].rho*pow(Node[S+1].U*AREF,2);
		FV[S].p0/*Pa*/ = TotalPressureBar(pPpt, Node[S+1].p_dash*pPpt->PREF, Node[S+1].T, Node[S+1].U*AREF)*1e5;
		FV[S].p0_old = FV[S].p0; // Save previous value (for evaluating STEADY conditions)
		FV[S].T0 = Node[S+1].T + pow(Node[S+1].U*(EX ? pPpt->AREFe : pPpt->AREFi),2)/(2*pPpt->cpAir(Node[S+1].T));
		
		// Use perfect gas equation to obtain initial mass and density
		FV[S].m = (FV[S].p0*FV[S].V)/(pPpt->R_air*FV[S].T0);
		// Initial enthalpy
		FV[S].h0 = pPpt->cpAir(FV[S].T0)*FV[S].T0;
	}
}

void CPipe::PrintToFileLocations(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca)
// ============================================================ //
// Prints instantaneous pipe object data to file.				//
// This function is called from the main function.				//
// ============================================================ //
{
	if(ntappings>0)
	{
		int f, m;
		f = 0;
	
		if(timestep==0)
		{
			if(EX) fprintf(FILE_LOC,"%s", Underline("Object results for exhaust pipe", "-", ID));
			else fprintf(FILE_LOC,"%s", Underline("Object results for intake pipe", "-", ID));

			fprintf(FILE_LOC,"\t\t");
			if(!pPpt->CONTINUOUS) fprintf(FILE_LOC,"\t\t"); // Periodic
				
			for(m=0; m<ntappings; ++m) 
				//fprintf(FILE_LOC,"%s%i%s%.1f%s\t\t\t\t\t\t\t", "Loc. [", m, "] at ", loc_measure[m]*100, "% length of pipe");
				fprintf(FILE_LOC,"%s%i%s%.1f%s\t\t\t\t\t\t", "Loc. [", m, "] at ", loc_measure[m]*100, "% length of pipe");
			fprintf(FILE_LOC,"\n");
		
			fprintf(FILE_LOC,"%s", "Time(s)");
			if(!pPpt->CONTINUOUS) // Periodic
			{
				fprintf(FILE_LOC,"\t%s%c%s\t%s%c%s", "Elapsed (", Deg(), "CA)", "Crank angle (", Deg(), "CA)");
			}

			for(m=0; m<ntappings; ++m)
				//fprintf(FILE_LOC,"\t\t%s\t%s\t%s\t%s\t%s\t%s", 
				//"P (bar)", "T (K)", "u (m/s)", "MFR (kg/s)", "PR", "MFP");
				fprintf(FILE_LOC,"\t\t%s\t%s\t%s\t%s\t%s",
				"P (bar)", "T (K)", "u (m/s)", "MFR (kg/s)", "Re");
			fprintf(FILE_LOC,"\n");
		}

		if(time>=print_from_time)// Only start recording data as specified
		{
			if(timestep%freq==0) // Print data at the specified sampling frequency
			{
				fprintf(FILE_LOC,"%f", time);

				if(!pPpt->CONTINUOUS) // Periodic
				{
					fprintf(FILE_LOC,"\t%f\t%f", ca_elapsed, ca);
				}

				for(m=0; m<ntappings; ++m)
					//fprintf(FILE_LOC,"\t\t%f\t%f\t%f\t%f\t%f\t%f", 
					//Pressure[m], Temperature[m], Velocity[m], MassFlowRate[m], PR[m], MFP[m]);
					fprintf(FILE_LOC,"\t\t%f\t%f\t%f\t%f\t%f", 
					Pressure[m], Temperature[m], Velocity[m], MassFlowRate[m], Re[m]);
				fprintf(FILE_LOC,"\n");
			}
		}
	}
}

void CPipe::PrintToFileFandE(CProperties* pPpt, int timestep, double time, double ca_elapsed, double ca)
// ============================================================ //
// Prints instantaneous pipe object data to file.				//
// This function is called from the main function.				//
// ============================================================ //
{
	if(ntappings>0)
	{
		//int f, n, m;
		//f = 0;
		int S;
		
		if(timestep==0)
		{
			if(EX) fprintf(FILE_FandE,"%s", Underline("FILLING AND EMPTYING RESULTS FILE FOR EXHAUST PIPE", "-", ID));
			else fprintf(FILE_FandE,"%s", Underline("FILLING AND EMPTYING RESULTS FILE FOR INTAKE PIPE", "-", ID));

			fprintf(FILE_FandE,"\t\t");
			if(!pPpt->CONTINUOUS) fprintf(FILE_FandE,"\t\t"); // Periodic
			
			for(S=0; S<N-1; ++S) // For each mesh (should only be 1 for FandE)
				fprintf(FILE_FandE,"%s%i%s%i\t\t\t\t\t\t\t", "Filling and emptying volume properties for this pipe, mesh no. ", S+1, " of ", N-1);
			fprintf(FILE_FandE,"\n");
		
			fprintf(FILE_FandE,"%s", "Time(s)");
			if(!pPpt->CONTINUOUS) // Periodic
			{
				fprintf(FILE_FandE,"\t%s%c%s\t%s%c%s", "Elapsed (", Deg(), "CA)", "Crank angle (", Deg(), "CA)");
			}

//			for(S=0; S<N-1; ++S) // For each mesh (should only be 1 for FandE)
//				fprintf(FILE_FandE,"\t\t%s\t%s\t%s\t%s", 
//								   "Mass (kg)", "p (bar)", "T (K)", "rho (kg/m^3)");

			for(S=0; S<N-1; ++S) // For each mesh (should only be 1 for FandE)
				fprintf(FILE_FandE,"\t\t%s\t%s\t%s\t%s", 
								   "Mass (kg)", "p0 (bar)", "T0 (K)", "rho0 (kg/m^3)");

			fprintf(FILE_FandE,"\n");
		}

		if(time>=print_from_time)// Only start recording data as specified
		{
			if(timestep%freq==0) // Print data at the specified sampling frequency
			{
				fprintf(FILE_FandE,"%f", time);

				if(!pPpt->CONTINUOUS) // Periodic
				{
					fprintf(FILE_FandE,"\t%f\t%f", ca_elapsed, ca);
				}

				for(S=0; S<N-1; ++S) // For each mesh (should only be 1 for FandE)
					fprintf(FILE_FandE,"\t\t%f\t%f\t%f\t%f", 
									   FV[S].m, FV[S].p0/1e5, FV[S].T0, FV[S].p0/(pPpt->R_air*FV[S].T0));

				fprintf(FILE_FandE,"\n");
			}
		}
	}
}

void CPipe::CloseFiles(CProperties* pPpt)
{
	if(ntappings>0)
	{
		fclose(FILE_LOC);
		if(METHOD==FandE) fclose(FILE_FandE);
	}
	for(int f=0; f<num_props_measured; ++f)	fclose(FILE_OVERALL_MOV[f]);
}


void CPipe::Backup(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Backup\n");}

///*
	for(int n=0; n<N; ++n)
	{
		Node_Backup[n] = Node[n];
	}
//*/
///*
	for(int p=0; p<num_pathlines; ++p)
	{
		PathLine_Backup[p] = PathLine[p];
	}
//*/

}

void CPipe::Restore(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Restore\n");}

	for(int n=0; n<N; ++n)
	{
		Node[n] = Node_Backup[n];
	}

	for(int p=0; p<num_pathlines; ++p)
	{
		PathLine[p] = PathLine_Backup[p];
	}
}

void CPipe::RunPropagation(CProperties* pPpt, double DELZ, int timestep, bool &rRESTORE)
// ============================================================ //
// Controls the execution of the selected propagation method	//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".RunPropagation");}

	if(METHOD==MMOC) {
		if(pPpt->HOMENTROPIC) HomentropicMOC(pPpt, DELZ);
		else NonHomentropicMOC(pPpt, DELZ, timestep, rRESTORE);
	}
	else {
		if(METHOD==W_ALPHA_BETA) {
			if(pPpt->COMBINED_WAB_MOC) {
				PathLines(pPpt, DELZ);								// Determine path lines at new time step
				PathLinesAtDuctEnds(pPpt, timestep);				// Test to see if path lines enter or leave pipe, adjust accordingly
				InterpolateEntropyAtMeshPoints(pPpt);				// Calculate the new entropy level at mesh points
				W_alpha_beta(pPpt, DELZ);							// Run the propagation scheme
				W_alpha_beta_derive_lambdas(pPpt);					// Derive the interior characteristics
				MOC_Non_Homentropic_Boundary(pPpt, DELZ, rRESTORE); // Propagate characteristics to the boundaries
			}
			else {
				W_alpha_beta(pPpt, DELZ);							// Run the propagation scheme
				W_alpha_beta_prep_bcs(pPpt);						// Infer conditions at boundaries
			}
		}
		else {
			if(METHOD==FandE) {
///*
				PathLines(pPpt, DELZ);					// Determine path lines at new time step
				PathLinesAtDuctEnds(pPpt, timestep);	// Test to see if path lines enter or leave pipe, adjust accordingly
				InterpolateEntropyAtMeshPoints(pPpt);	// Calculate the new entropy level at mesh points
//*/							
				//FillingAndEmptying(pPpt, DELZ);		// Run the filling and emptying scheme
				FillingAndEmptyingSimple(pPpt);
			}
		}
	}
}

double CPipe::TimeStepMOC(CProperties* pPpt)
// For both homentropic and non-homentropic
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".TimeStepMOC\n");}

		double DELZ = 1000; // Reset to large value
		double DELZ1;		// Temporary holder
		double DELZ1_R, DELZ1_L;		// Temporary holders

		for(int S=0; S<this->N; ++S)
		{
			// The max time step when going right:
			if(S<N-1)
			{
				DELZ1_R = (Node[S].DELX_R)/( (Node[S].CL1[R] + Node[S].CL2[R])/2.0 
								+ fabs(Node[S].CL1[R] - Node[S].CL2[R])/(pPpt->gammaAir(Node[S].T)-1) );
			}
			else DELZ1_R = 1000;

			// The max time step when going left:
			if(S>0)
			{
			DELZ1_L = (Node[S].DELX_L)/( (Node[S].CL1[R] + Node[S].CL2[R])/2.0 
							+ fabs(Node[S].CL1[R] - Node[S].CL2[R])/(pPpt->gammaAir(Node[S].T)-1) );
			}
			else DELZ1_L = 1000;

			// Take the smallest:
			if(DELZ1_R<=DELZ1_L) DELZ1 = DELZ1_R;
			else DELZ1 = DELZ1_L;
			
			if(DELZ1<=DELZ) DELZ = DELZ1;
		}
		return pPpt->courant*DELZ;	// Apply Courant number
}

void CPipe::HomentropicMOC(CProperties* pPpt, double DELZ)
// Homentropic Method of Characteristics
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".HomentropicMOC\n");}
	int S;
	// Calculate Riemann variables at mesh points
	for(S=1; S<N; ++S)
		// Calculating new value using values from the left of node S, because CL1s move to the right, over distance DELX_L of node S
		Node[S].CL1[R+1] = Node[S].CL1[R] + (DELZ/Node[S].DELX_L)*
								(pPpt->BB*Node[S-1].CL1[R] - pPpt->AA*Node[S-1].CL2[R])*
								(Node[S-1].CL1[R] - Node[S].CL1[R]);
	for(S=0; S<N-1; ++S)
		// Calculating new value using values from the right of node S, because CL2s move to the left, over distance DELX_R of node S
		Node[S].CL2[R+1] = Node[S].CL2[R] + (DELZ/Node[S].DELX_R)*
								(pPpt->BB*Node[S+1].CL2[R] - pPpt->AA*Node[S+1].CL1[R])*
								(Node[S+1].CL2[R] - Node[S].CL2[R]);
}

void CPipe::NonHomentropicMOC(CProperties* pPpt, double DELZ, int timestep, bool &rRESTORE) 
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".NonHomentropicMOC\n");}
	PathLines(pPpt, DELZ);						// Determine path lines at new time step
	PathLinesAtDuctEnds(pPpt, timestep);		// Test to see if path lines enter or leave pipe, adjust accordingly
	InterpolateEntropyAtMeshPoints(pPpt);		// Calculate the new entropy level at mesh points
	MOC_Non_Homentropic(pPpt, DELZ, rRESTORE);	// Calculate lambda_I', lambda_II' for at Z + DELZ for all mesh points
}

void CPipe::PathLines(CProperties* pPpt, double DELZ)
// Update the XK distance and state value AA for the next time step
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PathLines\n");}

	double d, U, A, Tg, f, del_AA, del_AA_heat_transfer;
	int k, S;
	// NB: CL index; doesn't matter whether 0 or 1 - they are the same at this point
	for(k=0; k<this->num_pathlines; ++k)	// For number of pathlines K
	{
		S=0;
		// Loop to find the surrounding two nodes
		while(this->Node[S+1].X<this->PathLine[k].XK && S<this->N-2) ++S;
		// Then the pathline is between nodes S and S+1

		// lambda_K at this time step, eq. 7.72 of Benson
		// Everything is at time Z
//		this->PathLine[k].lambda_K = 
//			this->Node[S+1].CL1[R] - 
//				((this->Node[S+1].X - this->PathLine[k].XK)/this->DELX)
//				*(this->Node[S+1].CL1[R] - this->Node[S].CL1[R]);
		
		// The DELX value here is the X of S+1 minus the X of S, surely
		this->PathLine[k].lambda_K = 
			this->Node[S+1].CL1[R] - 
				((this->Node[S+1].X - this->PathLine[k].XK)/(Node[S+1].X-Node[S].X))
				*(this->Node[S+1].CL1[R] - this->Node[S].CL1[R]);

		// beta_K at this time step, eq. 7.73 of Benson
		// Everything is at time Z
//		this->PathLine[k].beta_K = 
//			this->Node[S+1].CL2[R] - 
//				((this->Node[S+1].X - this->PathLine[k].XK)/this->DELX)
//				*(this->Node[S+1].CL2[R] - this->Node[S].CL2[R]);

		// The DELX value here is the X of S+1 minus the X of S, surely
		this->PathLine[k].beta_K = 
			this->Node[S+1].CL2[R] - 
				((this->Node[S+1].X - this->PathLine[k].XK)/(Node[S+1].X-Node[S].X))
				*(this->Node[S+1].CL2[R] - this->Node[S].CL2[R]);

///*
		if(PathLine[k].lambda_K<0)
		{
			cout << "Pipe[" << ID << "]:" << endl;
			cout << "PathLine[k=" << k << "].lambda_K = " << PathLine[k].lambda_K << endl;
			cout << "Node[S+1=" << S+1 << "].CL1[R] = " << Node[S+1].CL1[R] << endl;
			cout << "Node[S=" << S << "].CL1[R] = " << Node[S].CL1[R] << endl;
			cout << "Node[S+1=" << S+1 << "].X = " << Node[S+1].X << endl;
			cout << "Node[S=" << S << "].X = " << Node[S].X << endl;
			cout << "PathLine[k=" << k << "].XK = " << PathLine[k].XK << endl;
			cout << endl;
		}
		if(PathLine[k].beta_K<0)
		{
			cout << "Pipe[" << ID << "]:" << endl;
			cout << "PathLine[k=" << k << "].beta_K = " << PathLine[k].beta_K << endl;
			cout << "Node[S+1=" << S+1 << "].CL2[R] = " << Node[S+1].CL2[R] << endl;
			cout << "Node[S=" << S << "].CL2[R] = " << Node[S].CL2[R] << endl;
			cout << "Node[S+1=" << S+1 << "].X = " << Node[S+1].X << endl;
			cout << "Node[S=" << S << "].X = " << Node[S].X << endl;
			cout << "PathLine[k=" << k << "].XK = " << PathLine[k].XK << endl;
			cout << endl;
		}
//*/

		// Use the interpolated lambda_K and beta_K to get A fot Tg, and U, for q
		A = (this->PathLine[k].lambda_K + this->PathLine[k].beta_K)/2;
		Tg = pow(A, 2)*(EX ? pPpt->TREFe : pPpt->TREFi);
		U = (this->PathLine[k].lambda_K - this->PathLine[k].beta_K)/(pPpt->gammaAir(Tg)-1);
		
		// Interpolate between S and S+1 for d at K
//if(this->PathLine[k].XK < 0) this->PathLine[k].XK = this->Node[S].X;
//		d = this->Node[S+1].d - ((this->Node[S+1].X - this->PathLine[k].XK)/this->DELX)
//								*(this->Node[S+1].d - this->Node[S].d);
		d = this->Node[S+1].d - ((this->Node[S+1].X - this->PathLine[k].XK)/(Node[S+1].X-Node[S].X))
								*(this->Node[S+1].d - this->Node[S].d);
/*
if(d<0)
{
//	d = (this->Node[S+1].d + this->Node[S].d)/2; // Just take average

	cout << "PathLines: d = " << d << endl;
	cout << "this->Node[S+1=" << S+1 << "].d = " << this->Node[S+1].d << endl;
	cout << "this->Node[S=" << S << "].d = " << this->Node[S].d << endl;
	cout << "this->Node[S+1=" << S+1 << "].X = " << this->Node[S+1].X << endl;
	cout << "this->Node[S=" << S << "].X = " << this->Node[S].X << endl;
	cout << "this->PathLine[k=" << k << "].XK = " << this->PathLine[k].XK << endl;
	cout << endl;
}
*/
		double Re_temp = Node[S].Re;
		f = FrictionFactor(pPpt, d, Re_temp, "PathLines", this->CFTRANS);

		// Calculate (del_AA)heat transfer - eq. 7.95 of Benson
		double u_temp = U*AREF;
		double unit_rho = 1; // Use 1 such that hc returned is hc/rho
		double hc_per_unit_rho = ConvectiveHTCoefficient(pPpt, f, unit_rho, u_temp, Tg, Re_temp, (Node[S+1].X-Node[S].X)*pPpt->xref);
		double q = ((4*hc_per_unit_rho)/(unit_rho*d))*(this->Tw - Tg);
		del_AA_heat_transfer = ((pPpt->gammaAir(Tg)-1)/2)*(this->PathLine[k].AA/pow(A,2))
								*((q*pPpt->xref)/pow(AREF,3))
								*DELZ;

		// del_AA, to obtain the state AA_dash at next time step, derived from eq. 7.76 of Benson
		del_AA = ((pPpt->gammaAir(Tg)-1)/2)*( (4*this->PathLine[k].AA)/pow(this->PathLine[k].lambda_K + this->PathLine[k].beta_K, 2) )
				 *(2*f*pPpt->xref/d)*pow(fabs((this->PathLine[k].lambda_K - this->PathLine[k].beta_K)/(pPpt->gammaAir(Tg)-1)), 3)
				 *DELZ
				+ del_AA_heat_transfer;

		// Results:

		// XK at next time step, eq. 7.75 of Benson
		this->PathLine[k].XK_old = this->PathLine[k].XK;	// Record current value
		this->PathLine[k].XK = 
			this->PathLine[k].XK + 
			((this->PathLine[k].lambda_K - this->PathLine[k].beta_K)/(pPpt->gammaAir(Tg)-1))*DELZ;
/*
		if(this->PathLine[k].XK<0)
		{
			cout << "PathLine[k=" << k << "].XK = " << PathLine[k].XK << endl;
			cout << "PathLine[k=" << k << "].lambda_K = " << PathLine[k].lambda_K << endl;
			cout << "PathLine[k=" << k << "].beta_K = " << PathLine[k].beta_K << endl;
			cout << endl;
		}
*/
		// State at next time step
		this->PathLine[k].AA = this->PathLine[k].AA + del_AA;
	}
}

void CPipe::PathLinesAtDuctEnds(CProperties* pPpt, int timestep)
// Test for inflow/outflow/noflow at both ends and take appropriate action
// to remove a path line and insert a new path line in the case of inflow,
// and to replace path lines leaving the pipe for outflow
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PathLinesAtDuctEnds\n");}

	bool odd_end_adjusted = false;
	bool even_end_adjusted = false;

	// Go through and make all pathlines "old"
	for(int k=0; k<this->num_pathlines; ++k) this->PathLine[k].new_pathline = false;

	// Odd end ----------------------------------------
	// ------------------------------------------------
	if(odd_end_flow == INFLOW)
	{
		// Inflow into pipe from boundary at odd end
		// Create new pathline and remove the least effective one
		CPathLine NewPathLine;
		NewPathLine = this->PathLine[0]; // At time Z, current/old timestep
		// Assign the values of the first path line to the new one as a first guess
		// The new AA value will be ajusted later by the boundary
		
//NewPathLine.XK = 0;			// It must have position at odd end of pipe
								// NewPathLine should already have ID==0
NewPathLine.XK = Node[0].X;

		NewPathLine.XK_old = NewPathLine.XK; // N/A really, but otherwise would contain rubbish

		this->RemovePathLine(pPpt, ODD, timestep); 
		// Removes the least useful path line and shift path lines accordingly

		NewPathLine.new_pathline = true;
		this->PathLine[0] = NewPathLine;
		// Insert new path line at the vacated position

		odd_end_adjusted = true;
	}
	else
	{
		if(odd_end_flow == OUTFLOW)
		{
			//cout << "Pipe[" << this->ID << "] - odd end OUTFLOW" << endl;
			// Outflow from pipe at odd end
			// Replace leaving pathlines with dummies
			CPathLine NewPathLine;
			// Some values for all new dummies created:
//NewPathLine.XK = 0;	// At new time step
NewPathLine.XK = Node[0].X;
			NewPathLine.XK_old = NewPathLine.XK;
			// Interpolate between the last path line to leave at Z + DELZ
			// and the last path lines remaining in the pipe
			
			// Identify which path lines have left at time Z + DELZ
			// by testing their XK value

			int k=0; int count=0;
			// Loop to find the last path line to leave the pipe
			// XK values must now be at time Z + DELZ
			
			//while(this->PathLine[k].XK<0)
			while(this->PathLine[k].XK<Node[0].X)
			{
				++k;
				++count;	// Number of path lines leaving the pipe
			};
			// Then the interpolation should be between 
			// path lines k-1 (has left pipe) and k (last inside pipe)
			if(k>0 && k<this->num_pathlines)	// For extra safety
			{
				NewPathLine.AA = this->PathLine[k-1].AA + 
								((NewPathLine.XK - this->PathLine[k-1].XK)
								/(this->PathLine[k].XK - this->PathLine[k-1].XK))
								*(this->PathLine[k].AA - this->PathLine[k-1].AA);
			}
/*
			else	// No path lines outside pipe
			{
				if(this->m_EX) 
					cout << "Outflow but no path lines leaving odd end of exhaust pipe " << this->ID << endl;
				else
					cout << "Outflow but no path lines leaving odd end of intake pipe " << this->ID << endl;
			}
*/
			// The first "count" path lines in the list should be replaced with NewPathLine
			// More than one path line leaving leads to dummy copies of NewPathLine
			// These will be removed first by the removal routine
			for(k=0; k<count; ++k)
			{
				NewPathLine.ID = k;	// Must assign an ID
				NewPathLine.new_pathline = true;
				this->PathLine[k] = NewPathLine;
			}
			odd_end_adjusted = true;
		}
		else this->odd_end_flow = NOFLOW; // PathLine[0].XK == PathLine[0].XK_old
	}

	// Even end ----------------------------------------
	// -------------------------------------------------
	if(even_end_flow == INFLOW)
	{
		// Inflow into pipe from boundary at even end
		// Create new pathline and remove the least effective one
		CPathLine NewPathLine;
		NewPathLine = this->PathLine[num_pathlines-1]; // At time Z, current/old timestep
		// Assign the values of the last path line to the new one as a first guess
		// The new AKK value will be adjusted later by the boundary
		
//NewPathLine.XK = this->XPIPE;		// It must have position at even end of pipe
											// NewPathLine should already have ID==num_pathlines-1
NewPathLine.XK = Node[N-1].X;


		NewPathLine.XK_old = NewPathLine.XK;	// N/A really, but otherwise would contain rubbish

		RemovePathLine(pPpt, EVEN, timestep); 
		// Removes the least useful path line and shift path lines accordingly

		NewPathLine.new_pathline = true;
		this->PathLine[num_pathlines-1] = NewPathLine;
		// Insert new path line at the vacated position

		even_end_adjusted = true;
	}
	else
	{
		if(even_end_flow == OUTFLOW)
		{
			//cout << "Pipe[" << this->ID << "] - even end OUTFLOW" << endl;
			// Outflow from pipe at even end
			// Replace leaving pathlines with dummies
			CPathLine NewPathLine;
			// Some values for all new dummies created:
//NewPathLine.XK = this->XPIPE;	// At new time step
NewPathLine.XK = Node[N-1].X;

			NewPathLine.XK_old = NewPathLine.XK;

			// Interpolate between the last path line to leave at Z + DELZ
			// and the last path lines remaining in the pipe
				
			// Identify which path lines have left at time Z + DELZ
			// by testing their XK value

			int k=this->num_pathlines-1; int count=0;
			// Loop to find the last path line to leave the pipe
			// XK values must now be at time Z + DELZ
			
			//while(this->PathLine[k].XK>this->XPIPE)
			while(this->PathLine[k].XK>Node[N-1].X)
			{
				--k;
				++count;	// Number of path lines leaving the pipe
			};
			// Then the interpolation should be between 
			// path lines k (last inside pipe) and k+1 (has left pipe)
			if(k<this->num_pathlines-1 && k>=0)	// For extra safety
			{
				NewPathLine.AA = this->PathLine[k].AA + 
								((NewPathLine.XK - this->PathLine[k].XK)
								/(this->PathLine[k+1].XK - this->PathLine[k].XK))
								*(this->PathLine[k+1].AA - this->PathLine[k].AA);
			}
/*
			else	// No path lines outside pipe
			{
				if(this->m_EX) 
					cout << "Outflow but no path lines leaving even end of exhaust pipe " << this->ID << endl;
				else
					cout << "Outflow but no path lines leaving even end of intake pipe " << this->ID << endl;
			}
*/
			// The last "count" path lines in the list should be replaced with NewPathLine
			// More than one path line leaving leads to dummy copies of NewPathLine
			// These will be removed first by the removal routine
			for(k=this->num_pathlines-1; k>(this->num_pathlines-1-count); --k)
			{
				NewPathLine.ID = k;	// Must assign an ID
				NewPathLine.new_pathline = true;
				this->PathLine[k] = NewPathLine;
			}
			even_end_adjusted = true;
		}
		else this->even_end_flow = NOFLOW;
	}

/*
	if(!odd_end_adjusted)
	{
		if(this->m_EX)
			cout << "No adjustment at odd end of exhaust pipe " << this->ID << endl;
		else
			cout << "No adjustment at odd end of intake pipe " << this->ID << endl;
	}

	if(!even_end_adjusted)
	{
		if(this->m_EX)
			cout << "No adjustment at even end of exhaust pipe " << this->ID << endl;
		else
			cout << "No adjustment at even end of intake pipe " << this->ID << endl;
	}
*/
}

void CPipe::InterpolateEntropyAtMeshPoints(CProperties* pPpt)
// Linearly interpolate path line characteristics either side of each mesh point
// Call this after path lines have been updated for the next time step
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".InterpolateEntropyAtMeshPoints\n");}

	int S, k;
	bool located, located2;

	for(S=0; S<this->N; ++S)	// For each node
	{
		k=0;
		located = false;
		// Loop to find the surrounding two path lines
		while(!located)
		{
			if( (this->PathLine[k+1].XK>this->Node[S].X ||
					fabs(this->PathLine[k+1].XK - this->Node[S].X) < 1e-6)
				 && this->PathLine[k+1].XK>this->PathLine[k].XK )
				located = true;
			else
			{
				if(k==this->num_pathlines-2)
					// Then the Node[S].X value is out of the range
					// of pathlines. We must find the two closest
					// pathlines of non-equal XK value from which
					// to extrapolate
				{
					located2 = false;
					while(!located2)
					{
						if(this->PathLine[k+1].XK>this->PathLine[k].XK) located2 = true;
						else --k;
					}
					located = true;
				}
				else ++k;
			}
		};
		// Then the node is between path lines k and k+1

		// Where the AAs are those for the new time step, the new entropy level
		// at the mesh point is:

		// Until here AA[1] and AA[0] are the same
		this->Node[S].AA[R+1] = (this->Node[S].X - this->PathLine[k].XK)
							  /(this->PathLine[k+1].XK - this->PathLine[k].XK)
							  *(this->PathLine[k+1].AA - this->PathLine[k].AA)
							  + this->PathLine[k].AA;
/*
		if((this->PathLine[k+1].XK - this->PathLine[k].XK) == 0)
		{
			cout << "Node[S].X = " << Node[S].X << endl;
			cout << "PathLine[k+1].XK = " << PathLine[k+1].XK << endl;
			cout << "PathLine[k].XK = " << PathLine[k].XK << endl;
		}
*/
		// Hence calculate del_AA if you need to?
		// Actually no, but need to remember new (above, R=1) and old (R=0) 
		// values for the del_lambda_entropy equation.
		// Hence only call MOC_Non_Homentropic after this is done, else you won't no the values
	}
}

void CPipe::MOC_Non_Homentropic(CProperties* pPpt, double DELZ, bool &rRESTORE)
// Non-Homentropic Method of Characteristics
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".MOC_Non_Homentropic\n");}
	int S;
	for(S=1; S<N; ++S)		// For lambda_I characteristics
		MOC_Non_Homentropic_lambdaI(pPpt, DELZ, S, rRESTORE);	// Uses values at [S] and [S-1]
	for(S=0; S<N-1; ++S)	// For lambda_II characteristics
		MOC_Non_Homentropic_lambdaII(pPpt, DELZ, S, rRESTORE);	// Uses values at [S] and [S+1]
}

void CPipe::MOC_Non_Homentropic_Interior(CProperties* pPpt, double DELZ, bool &rRESTORE)
// Non-Homentropic Method of Characteristics
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".MOC_Non_Homentropic_Interior\n");}
	int S;
//	for(S=1; S<N; ++S)	// For lambda_I characteristics
	for(S=1; S<N-1; ++S)	// For lambda_I characteristics
	{
		MOC_Non_Homentropic_lambdaI(pPpt, DELZ, S, rRESTORE);
		// Uses values at [S] and [S-1] to give Node[S].CL1[R+1]
	}

//	for(S=0; S<N-1; ++S)	// For lambda_II characteristics
	for(S=1; S<N-1; ++S)	// For lambda_II characteristics
	{
		MOC_Non_Homentropic_lambdaII(pPpt, DELZ, S, rRESTORE);
		// Uses values at [S] and [S+1]	to give Node[S].CL2[R+1]
	}
}

void CPipe::MOC_Non_Homentropic_Boundary(CProperties* pPpt, double DELZ, bool &rRESTORE)
// Non-homentropic Method of Characteristics
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".MOC_Non_Homentropic_Boundary\n");}

	// For lambda_I characteristic at the even end
	int S=N-1;
	MOC_Non_Homentropic_lambdaI(pPpt, DELZ, S, rRESTORE); 
	// Uses values at [S] and [S-1] to give Node[S].CL1[R+1]
	
	// For lambda_II characteristic at the odd end
	S=0;
	MOC_Non_Homentropic_lambdaII(pPpt, DELZ, S, rRESTORE); 
	// Uses values at [S] and [S+1]	to give Node[S].CL2[R+1]
}

void CPipe::MOC_Non_Homentropic_lambdaI(CProperties* pPpt, double DELZ, int S, bool &rRESTORE)
// Non-homentropic Method of Characteristics
// Propagates lambdaI for a single node
// For lambda_I characteristics
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".MOC_Non_Homentropic_lambdaI\n");}

	char pause;

	double A, B, lambda_A, lambda_B, beta_A, beta_B;
	double XW_sign, dDdX_sign;
	double AAR, AAQ, AAR_dash;
	double DELX, XR, delXdX;
	double f, Tg;

	bool BC_even_lambda1;
	BC_even_lambda1 = false;

	// Apply subsonic entry parameters for lambda_I
	// If actually supersonic but U is +ve, these parameters do not change
	A = pPpt->AA;
	B = pPpt->BB;
	lambda_A = Node[S].CL1[R];
	lambda_B = Node[S-1].CL1[R];
	beta_A = Node[S].CL2[R];
	beta_B = Node[S-1].CL2[R];
	XW_sign = 1; // XP_sign
	dDdX_sign = 1;
	AAR = Node[S].AA[R];
	AAQ = Node[S-1].AA[R];

	DELX = Node[S].DELX_L; // See just above, we use values at S-1 and S, so we want DELX_L
	if(Node[S].DELX_L==0)
	{
		cout << "Node[" << S << "].DELX_L = " << Node[S].DELX_L << endl;
		exit(1);
	}
	XR = Node[S].X;

	if((B*lambda_A - A*beta_A)<0) 
	// Then supersonic and lambda_I switches to an (l, l+1) domain
	// This implies U is -ve
	{
		A = -pPpt->AA;
		B = -pPpt->BB;
		XW_sign = -1;

		if(S==N-1)	// [S+1] out of range - must be supplied by the even end boundary
		{
			BC_even_lambda1 = true; // Flag the even end b.c. to supply required values
			// Extrapolate from interior values
			lambda_B = Node[S].CL1[R] +
					   (Node[S].X - Node[S-1].X)
					 *((Node[S].CL1[R] - Node[S-1].CL1[R])/(Node[S].X - Node[S-1].X));
			beta_B = Node[S].CL2[R] +
					   (Node[S].X - Node[S-1].X)
					 *((Node[S].CL2[R] - Node[S-1].CL2[R])/(Node[S].X - Node[S-1].X));
			AAQ = Node[S].AA[R] +
					   (Node[S].X - Node[S-1].X)
					 *((Node[S].AA[R] - Node[S-1].AA[R])/(Node[S].X - Node[S-1].X));
		}
		else
		{
			BC_even_lambda1 = false;
			lambda_B = Node[S+1].CL1[R];
			beta_B = Node[S+1].CL2[R];
			AAQ = Node[S+1].AA[R];
		}
	}
		
//	if(!BC_even_lambda1)
	{
		// Calculate delX/dX using eq. 7.71 of Benson - this must be positive
		delXdX = (B*lambda_A - A*beta_A)
				/((DELX/DELZ) + B*(lambda_A - lambda_B) - A*(beta_A - beta_B));

double Re_temp = Node[S].Re;
double d_temp = Node[S].d;
if(d_temp<0) cout << "MOC_Non_Homentropic_lambdaII: d_temp = " << d_temp << endl;


		f = FrictionFactor(pPpt, d_temp, Re_temp, "MOC_Non_Homentropic_lambdaI", this->CFTRANS);

		// Heat transfer - gas temperature (K)
		Tg = pow(Node[S].A,2)*(EX ? pPpt->TREFe : pPpt->TREFi);

		AAR_dash = Node[S].AA[R+1]; // Must run the entropy interpolation before here

		// Set appropriate characteristic to newly calculated lambda_R_dash
		Node[S].CL1[R+1] = 
			MOC_Non_Homentropic_Calcs(pPpt, S, lambda_A, lambda_B, beta_A, beta_B, delXdX, 
											AAR, AAQ, AAR_dash, DELX, XR, XW_sign, dDdX_sign, 
											Node[S].d, Tg, DELZ, f, Re_temp);
//if(S>=49 && S<=52)
if(false)
{
	cout << "S=" << S << endl;
	cout << "lambda_A = " << lambda_A << endl;
	cout << "lambda_B = " << lambda_B << endl;
	cout << "beta_A = " << beta_A << endl;
	cout << "beta_B = " << beta_B << endl;
	cout << "delXdX = " << delXdX << endl;
	cout << "AAR = " << AAR << endl;
	cout << "AAQ = " << AAQ << endl;
	cout << "AAR_dash = " << AAR_dash << endl;
	cout << "DELX = " << DELX << endl;
	cout << "XR = " << XR << endl;
	cout << "XW_sign = " << XW_sign << endl;
	cout << "dDdX_sign = " << dDdX_sign << endl;
	cout << "Node[S].d = " << Node[S].d << endl;
	cout << "Tg = " << Tg << endl;
	cout << "DELZ = " << DELZ << endl;
	cout << "f = " << f << endl;
	cout << "Re_temp = " << Re_temp << endl;
	cout << "Node[S=" << S << "].CL1[R+1] = " << Node[S].CL1[R+1] << endl << endl;
}

		// Set new entropy level
//		this->Node[S].AA[R+1] = AAR_dash;	// REMOVE when ready
	}
	
	if(BC_even_lambda1 && !pPpt->SUPERSONIC)
	{
		if(EX) cout << "Supersonic at odd end, Exhaust[" << ID << "]\n"; 
		else cout << "Supersonic at odd end, Intake[" << ID << "]\n";
		cin >> pause;
		rRESTORE = true;	
	}
}


void CPipe::MOC_Non_Homentropic_lambdaII(CProperties* pPpt, double DELZ, int S, bool &rRESTORE)
// Non-homentropic Method of Characteristics
// Propagates lambdaII for a singles node
// For lambda_II characteristics
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".MOC_Non_Homentropic_lambdaII\n");}

	char pause;

	double A, B, lambda_A, lambda_B, beta_A, beta_B;
	double XW_sign, dDdX_sign;
	double AAR, AAQ, AAR_dash;
	double DELX, XR, delXdX;
	double f, Tg;

	bool BC_odd_lambda11;

	BC_odd_lambda11 = false;
	// Apply subsonic entry parameters for lambda_II
	// If actually supersonic but U is -ve, these parameters do not change
	A = pPpt->AA;
	B = pPpt->BB;
	lambda_A = Node[S].CL2[R];
	lambda_B = Node[S+1].CL2[R];
	beta_A = Node[S].CL1[R];
	beta_B = Node[S+1].CL1[R];
	XW_sign = -1;
	dDdX_sign = -1;
	AAR = Node[S].AA[R]; 
	AAQ = Node[S+1].AA[R]; 

	DELX = Node[S].DELX_R; // See just above, we use values at S and S+1, so we want DELX_R
	if(Node[S].DELX_R==0)
	{
		cout << "Node[" << S << "].DELX_R = " << Node[S].DELX_R << endl;
		exit(1);
	}
	XR = Node[S].X;

	if((B*lambda_A - A*beta_A)<0) 
	// Then supersonic and lambda_II switches to an (l-1, l) domain
	// This implies U is +ve
	{
		A = -pPpt->AA;
		B = -pPpt->BB;
		XW_sign = 1;
		
		if(S==0)	// [S-1] out of range - must be supplied by the odd end boundary
		{
			BC_odd_lambda11 = true; // Flag the odd end b.c. to supply required values
			// Extrapolate from interior values
			lambda_B = Node[S].CL2[R] +
					   (Node[S].X - Node[S+1].X)
					 *((Node[S].CL2[R] - Node[S+1].CL2[R])/(Node[S].X - Node[S+1].X));
			beta_B = Node[S].CL1[R] +
					   (Node[S].X - Node[S+1].X)
					 *((Node[S].CL1[R] - Node[S+1].CL1[R])/(Node[S].X - Node[S+1].X));
			AAQ = Node[S].AA[R] +
					   (Node[S].X - Node[S+1].X)
					 *((Node[S].AA[R] - Node[S+1].AA[R])/(Node[S].X - Node[S+1].X));
		}
		else
		{
			BC_odd_lambda11 = false; // Flag the even end b.c. to supply required values
			lambda_B = Node[S-1].CL2[R];
			beta_B = Node[S-1].CL1[R];
			AAQ = Node[S-1].AA[R];
		}
	}

//	if(!BC_odd_lambda11)
	{
		// Calculate delX/dX using eq. 7.71 of Benson - this must be positive
		delXdX = (B*lambda_A - A*beta_A)
				/((DELX/DELZ) + B*(lambda_A - lambda_B) - A*(beta_A - beta_B));

double Re_temp = Node[S].Re;
double d_temp = Node[S].d;
if(d_temp<0) cout << "MOC_Non_Homentropic_lambdaII: d_temp = " << d_temp << endl;

		f = FrictionFactor(pPpt, d_temp, Re_temp, "MOC_Non_Homentropic_lambdaII", this->CFTRANS);

		// Heat transfer - gas temperature (K)
		Tg = pow(Node[S].A,2)*(EX ? pPpt->TREFe : pPpt->TREFi);

		AAR_dash = Node[S].AA[1]; // Must run the entropy interpolation before here

		// Set appropriate characteristic to newly calculated lambda_R_dash
		Node[S].CL2[R+1] = 
			MOC_Non_Homentropic_Calcs(pPpt, S, lambda_A, lambda_B, beta_A, beta_B, delXdX, 
											AAR, AAQ, AAR_dash, DELX, XR, XW_sign, dDdX_sign, 
											Node[S].d, Tg, DELZ, f, Re_temp);
		// Set new entropy level
//		this->Node[S].AA[R+1] = AAR_dash;	// NOT NECESSARY remove when ready
	}

	if(BC_odd_lambda11 && !pPpt->SUPERSONIC)
	{
		if(EX) cout << "Supersonic at odd end, Exhaust[" << ID << "]\n"; 
		else cout << "Supersonic at odd end, Intake[" << ID << "]\n";
		cin >> pause;
		rRESTORE = true;
	}
}


double CPipe::MOC_Non_Homentropic_Calcs(CProperties* pPpt, int S, double lambda_A, double lambda_B, 
								double beta_A, double beta_B, double delXdX, 
								double AAR, double AAQ, double AAR_dash, double DELX,
								double XR, double XW_sign, double dDdX_sign, 
								double d, double Tg, double DELZ, double f, double Re)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".MOC_Non_Homentropic_Calcs\n");}

	double lambda_P, beta_P, AAP;
	double XP, ONEoD_dDdX;
	double d_lambda_area, d_lambda_entropy, d_lambda_friction, d_lambda_heat_transfer;
	double d_lambda;

	// Calculate lambda_P - eq. 7.84 of Benson
	lambda_P = lambda_A - delXdX*(lambda_A - lambda_B);
	// Calculate beta_P - eq. 7.85 of Benson
	beta_P = beta_A - delXdX*(beta_A - beta_B);
	// Calculate AAP - eq. 7.86 of Benson
	AAP = AAR - delXdX*(AAR - AAQ);
	// Calculate XP - eq. 7.90 of Benson
	XP = XW_sign*(S*DELX - delXdX*DELX);			// Here S is equivalent to l-1 in eq. 7.90

	// Calculate (1/D)*dD/dX - eq. 7.89 of Benson (for linear case)
/*
	if(LINEAR) ONEoD_dDdX = (2*C)/(2*d_odd + C*(XR + XP));
	else
	{
		double X_AV = (XR + XP)/2;
		ONEoD_dDdX = (2*a + b)/(a*pow(X_AV,2) + b*X_AV + c);
	}
*/
//	double X_AV = (XR + XP)/2;
//	ONEoD_dDdX = (2*A + B)/(A*pow(X_AV,2) + B*X_AV + C);
	

	ONEoD_dDdX = (2*C)/(2*d_odd + C*(XR + XP));

/*
if(S==51)
{
//cout << "lambda_P = " << lambda_P << endl;
//cout << "beta_P = " << beta_P << endl;
//cout << "dDdX_sign = " << dDdX_sign << endl;
cout << "XR = " << XR << endl;
cout << "XP = " << XP << endl;
cout << "C = " << C << endl;
cout << "ONEoD_dDdX = " << ONEoD_dDdX << endl;
//cout << "DELZ = " << DELZ << endl;
}
//*/

	// Calculate component d_lambdas

	// Calculate d_lambda_area - eq. 7.87 of Benson
	d_lambda_area = -0.5*(lambda_P + beta_P)*(lambda_P - beta_P)*(dDdX_sign*ONEoD_dDdX)*DELZ;
	// NB: (dDdX_sign*ONEoD_dDdX) is ultimately +ve if diameter increases in direction 
	// of travel of the characteristic. e.g. lambda_I, subsonic, travels rightwards, and 
	// pipe increases in diameter going rightwards. Hence () should be +ve hence 
	// dDdX_sign is +ve.

	// Calculate d_lambda_entropy - eq. 7.81 of Benson
	d_lambda_entropy = 0.5*(lambda_P + beta_P)*(AAR_dash - AAP)/AAP;

		
	// Calculate d_lambda_friction - eq. 7.82 of Benson
	if(fabs(lambda_P - beta_P) < 1e-9) d_lambda_friction = 0;
	else
		d_lambda_friction = -(pPpt->gammaAir(Tg)-1)*(f*DELX/d)*(pPpt->xref/DELX)
						*pow((lambda_P - beta_P)/(pPpt->gammaAir(Tg)-1), 2)
						*((lambda_P - beta_P)/fabs(lambda_P - beta_P))
						*((3*beta_P - lambda_P)/(lambda_P + beta_P))
						*DELZ;

	// Calculate d_lambda_heat_transfer - eq. 7.83 of Benson
	double u_temp = (lambda_P - beta_P)/(pPpt->gammaAir(Tg)-1) * AREF;
	double unit_rho = 1; // Use 1 such that hc returned is hc/rho
	double hc_per_unit_rho = ConvectiveHTCoefficient(pPpt, f, unit_rho, u_temp, Tg, Re, DELX*pPpt->xref);
	double q = ((4*hc_per_unit_rho)/(unit_rho*d))*(Tw - Tg);
	d_lambda_heat_transfer = (pow(pPpt->gammaAir(Tg)-1, 2)/2)
							*q*pPpt->xref*(1/pow(AREF, 3))
							*(2/(lambda_P + beta_P))
							*DELZ;

	// Calculate combined d_lambda
	d_lambda = d_lambda_area + d_lambda_entropy + d_lambda_friction + d_lambda_heat_transfer;

	// Calculate new lambda_R at next time step and return it
	return (lambda_P + d_lambda);
}


void CPipe::RemovePathLine(CProperties* pPpt, int SIDE, int timestep)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".RemovePathLine\n");}
	int k;
	double AA_temp, del_AA, del_AA_min;
	int remove;
	// Reset temp values
	del_AA_min = 1.0e6;

	// Identify the most useless path line
	remove = 1; // Initialize
	for(k=1; k<this->num_pathlines-1; ++k)
	{
		// Interpolate AA between k-1 and k+1 
		// to calculate imaginary AA at path line k as if
		// pathline k were not there
		AA_temp = this->PathLine[k-1].AA + 
					((this->PathLine[k].XK - this->PathLine[k-1].XK)
					/(this->PathLine[k+1].XK - this->PathLine[k-1].XK))
					*(this->PathLine[k+1].AA - this->PathLine[k-1].AA);
		// Calculate the difference
		del_AA = fabs(AA_temp - this->PathLine[k].AA);

		if(del_AA<del_AA_min)
		{ 
			del_AA_min = del_AA;
			// Mark pathline with the smallest difference to be removed 
			remove = k;
		}
	}

	// Need to shift path lines in the list appropriately
	// If new path line is from the odd/LH end, shift all path lines
	// above "remove" down a place, allowing the new path line to be inserted at [0]
	// If new path line is from the even/RH end, shift all path lines
	// below "remove" up a place, allowing the new path line to be 
	// inserted at [num_pathlines-1]

	if(SIDE==ODD)
	{
		// Must start with the lowest in the list above "remove" to preserve data
		for(k=remove; k>0; --k)
		{
			this->PathLine[k] = this->PathLine[k-1];
			// Copy all data, but must retain ID accordingly to position in list
			this->PathLine[k].ID = k;
		}
	}
	else // (SIDE==EVEN)
	{
		// Must start with the path line just after "remove" to preserve data
		for(k=remove; k<this->num_pathlines-1; ++k)
		{
			this->PathLine[k] = this->PathLine[k+1];
			// Copy all data, but must retain ID accordingly to position in list
			this->PathLine[k].ID = k;
		}
	}
}

void CPipe::W_alpha_beta(CProperties* pPpt, double DELZ)
// ====================================================================================================
// General W_alpha_beta classification scheme					
// ----------------------------------------------------------------------------------------------------					
// Lerat and Peyret general class W_alpha_beta schemes			
// - Explicit conservative										
// - Two-step predictor-corrector								
// - Three-point space stencil, two-level point stencil			
// - Second order accuracy in time and space					
//																
// Scheme								alpha			beta	
// ------								-----			----	
// Two-step Lax-Wendroff (LW2)			0.5				0.5		
// MacCormack method (MAC):										
// - forward pred., rearward corr.		1				0		
// - rearward pred., forward corr.		1				1		
// Lerat + Peyret optimum (LP)			1 + sqrt(5)/2	0.5		
//																
// If in doubt, specify LW2!									
// ====================================================================================================
{
	if (pPpt->SHOW_calls) { pPpt->Out(Identify()); pPpt->Out(".W_alpha_beta\n"); }

	int S, K;
	double k;
	double f_of_R0; // Haemodynamics

	// Calculate time step in s
	// ----------------------------------------------------------------------------------------------------
	double delt = (DELZ / AREF) * pPpt->xref;

	// Predictor step
	// ----------------------------------------------------------------------------------------------------

	// Calculate flux vector F from known solution vector W
	// ----------------------------------------------------------------------------------------------------
	
	for (S = 0; S < N; ++S) { // For all nodes - requires separate loop since S+1 is referenced in predictor equation
		
		Node[S].F[0] = Node[S].W[1];

		k = pPpt->gammaAir(Node[S].T);

		if (pPpt->HAEMODYNAMICS) {
/*		
			f_of_R0 = pPpt->k1 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2) + pPpt->k3;

			//Node[S].B = (pPpt->k1 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2) + pPpt->k3) / Node[S].rho * sqrt(Node[S].f0 * Node[S].f);
			//Node[S].B = f_of_R0 / Node[S].rho * sqrt(Node[S].f0 * Node[S].f);
			//Node[S].B = (Node[S].f0 / Node[S].rho) * ( (2/9)*f_of_R0 * (1 - pow(Node[S].f/Node[S].f0,3)) + (Node[S].p0_dash*pPpt->PREF*1e5) );
			
			Node[S].B = (Node[S].f0 / Node[S].rho) * ( (2/9)*f_of_R0 * (pow(Node[S].f/Node[S].f0,3) - 1) + (Node[S].p0_dash*pPpt->PREF*1e5) );

			//Node[S].B = (Node[S].f0 / Node[S].rho) * (Node[S].p0_dash * pPpt->PREF * 1e5);

			Node[S].F[1] = pow(Node[S].W[1], 2) / Node[S].W[0] 
				+ Node[S].B * Node[S].rho;

//*/
			Node[S].B = Node[S].b / Node[S].rho * sqrt(Node[S].f0 * Node[S].f) - Node[S].b / Node[S].rho * Node[S].f0;//Calculate B at this timestep

			Node[S].F[1] = pow(Node[S].W[1], 2) / Node[S].W[0] + Node[S].B * Node[S].rho; // Calculate the flux

		}
		else {
			Node[S].F[1] = pow(Node[S].W[1], 2) / Node[S].W[0]
				+ (k - 1) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0]));
		}

/*
cout << "pow(Node[S].W[1], 2) / Node[S].W[0] = " << pow(Node[S].W[1], 2) / Node[S].W[0] << endl;
cout << "Node[S=" << S << "].B * rho = " << Node[S].B * Node[S].rho << endl;
cout << "(k - 1) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0])) = " << (k - 1) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0])) << endl;
cout << "Node[S].F[1] = pow(Node[S].W[1], 2) / Node[S].W[0] + Node[S].B * Node[S].rho = " << pow(Node[S].W[1], 2) / Node[S].W[0] + Node[S].B * Node[S].rho << endl;
cout << "Node[S].F[1] = pow(Node[S].W[1], 2) / Node[S].W[0] + (k - 1) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0])) = " << pow(Node[S].W[1], 2) / Node[S].W[0] + (k - 1) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0])) << endl;
cout << endl;
//*/

//Node[S].F[1] = pow(Node[S].W[1], 2) / Node[S].W[0];
//Node[S].F[1] = pow(Node[S].W[1], 2) / Node[S].W[0] + (k - 1) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0]));
//Node[S].F[1] = pow(Node[S].W[1], 2) / Node[S].W[0] + Node[S].B * Node[S].rho;
//Node[S].F[1] = pow(Node[S].W[1], 2) / Node[S].W[0] + (k - 1) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0])) + Node[S].B * Node[S].rho;

		Node[S].F[2] = (Node[S].W[1] * Node[S].W[2] / Node[S].W[0])
			+ (Node[S].W[1] / Node[S].W[0]) * (k - 1) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0]));
	}

	// Calculate source vector C via temporary primitives derived from known solution vector W
	// ----------------------------------------------------------------------------------------------------
	if (pPpt->SOURCES) {

		double rho_temp, p_temp, u_temp, T_temp, mu_air_temp, G;
		double Re_temp, ff_temp, ff_temp_C_p;
		double hc, q, w;

		//for(S=1; S<N; ++S) // For interior nodes, AND S==N-1 (C at S==N-1 required for predictor)
		for (S = 0; S < N; ++S) { // For all nodes
			k = pPpt->gammaAir(Node[S].T);

			// Continuity source term is always zero
			// ----------------------------------------------------------------------------------------------------
			Node[S].C[0] = 0;

			//// Momentum source due to change in area and friction (-p*dF/dx + rho*G*F)
			//// ----------------------------------------------------------------------------------------------------
			rho_temp = Node[S].W[0] / Node[S].f;
			p_temp = ((k - 1) / Node[S].f) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0]));
			u_temp = Node[S].W[1] / Node[S].W[0];
			T_temp = p_temp / (rho_temp * pPpt->R_air);
			mu_air_temp = pPpt->ViscosityAir(T_temp);

			//// Circular section pipe
			//// ----------------------------------------------------------------------------------------------------
			Re_temp = rho_temp * fabs(u_temp) * Node[S].d / mu_air_temp;
			ff_temp = FrictionFactor(pPpt, Node[S].d, Re_temp, "W_alpha_beta", CFTRANS);
			ff_temp_C_p = (Node[S].d / (((Node[S].DELX_L + Node[S].DELX_R) / 2) * pPpt->xref)) * C_p * CPTRANS;
			G = 0.5 * u_temp * fabs(u_temp) * (ff_temp + ff_temp_C_p) * (4 / Node[S].d);

			// Haemodynamics momentum sources (cu + dB/dr0 * dr0/dx)
			//----------------------------------------------------------------------------------------------------
			if (pPpt->HAEMODYNAMICS) {
				double delta = pPpt->BLThickness(Node[S].d);//double delta = 0.2 * d_pred_temp; //Boundary layer thickness, value entirely arbatory at the moment need to find a proper model for this!!
				double viscosity = pPpt->ViscosityBloodFunc(Node[S].T);//double viscosity = 0.003; //Fairly typical value in Pa.s.  Need to find a better model as non newtonian.
/*
				rho_temp = Node[S].W[0] / Node[S].f;

				Node[S].C[1] = -PI * Node[S].d * u_temp * viscosity / delta +
					(2 * sqrt(Node[S].f) * (sqrt(PI) * 4 / 3 * (pPpt->k1 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2) + pPpt->k3) + sqrt(Node[S].f0) * 4 / 3 * (pPpt->k1 * pPpt->k2 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2)))
						- Node[S].f * 4 / 3 * (pPpt->k1 * pPpt->k2 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2))) * Node[S].dd0dx / (2);
//*/			
				Node[S].C[1] = -PI * Node[S].d * u_temp * viscosity / delta +
					(2 * sqrt(Node[S].f * Node[S].f0) - Node[S].f - Node[S].f0) * Node[S].dbdr0 + 2 * sqrt(PI) * Node[S].b * (sqrt(Node[S].f) - sqrt(Node[S].f0)) * Node[S].dd0dx / 2;

			}
			else Node[S].C[1] = (-1) * p_temp * Node[S].dfdx + Node[S].W[0] * G; // Otherwise use (-p*dF/dx + rho*G*F)

//Node[S].C[1] = (-1) * p_temp * Node[S].dfdx + Node[S].W[0] * G;

			// Source due to work input and heat transfer (rho*w*F - rho*q*F)
			// ----------------------------------------------------------------------------------------------------

			// Source due to heat transfer (... - rho*q*F)
			// ----------------------------------------------------------------------------------------------------
			hc = ConvectiveHTCoefficient(pPpt, ff_temp, rho_temp, u_temp, T_temp, Re_temp, Node[S].DELX_R * pPpt->xref);
			q = ((4 * hc) / (rho_temp * Node[S].d)) * (Tw - T_temp); // General equation for heat transfer rate per unit mass

			// Source due to work input (rho*w*F - ...)
			// ----------------------------------------------------------------------------------------------------
			w = 0;

			Node[S].C[2] = Node[S].W[0] * (w - q);
		}
	}
	else { // No source terms except those due to area change  
		for (S = 0; S < N; ++S) { // For all nodes
			k = pPpt->gammaAir(Node[S].T);
			Node[S].C[0] = 0;
			double p_temp = ((k - 1) / Node[S].f) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0]));
			Node[S].C[1] = (-1) * p_temp * Node[S].dfdx;
			Node[S].C[2] = 0;
		}
	}

	// Artificial viscosity based on current pressure values
	// ----------------------------------------------------------------------------------------------------
	if (pPpt->VISC) {
		double* p_visc;
		p_visc = new double[N];
		for (S = 0; S < N; ++S) {
			k = pPpt->gammaAir(Node[S].T);
			p_visc[S] = ((k - 1) / Node[S].f) * (Node[S].W[2] - 0.5 * (pow(Node[S].W[1], 2) / Node[S].W[0]));
		}

		for (S = 1; S < N - 1; ++S) { // Calculate artificial viscosity S - note uses both forward and rearward values
			for (K = 0; K < 3; ++K)
				Node[S].S[K] = pPpt->Cx * (fabs((p_visc[S + 1]) - 2 * (p_visc[S]) + (p_visc[S - 1]))
					/ ((p_visc[S + 1]) + 2 * (p_visc[S]) + (p_visc[S - 1])))
				* (Node[S + 1].W[K] - 2 * Node[S].W[K] + Node[S - 1].W[K]);
		}
		S = 0; for (K = 0; K < 3; ++K) Node[S].S[K] = 0;
		S = N - 1; for (K = 0; K < 3; ++K) Node[S].S[K] = 0;
		delete[] p_visc;
	}
	else for (S = 0; S < N; ++S) for (K = 0; K < 3; ++K) Node[S].S[K] = 0;

	// Predictor equation - calculate predicted solution vector W_pred(i + beta, n + alpha)
	// ----------------------------------------------------------------------------------------------------
	for (S = 0; S < N - 1; ++S) {
		for (K = 0; K < 3; ++K) {

			Node[S].W_pred[K] = (1 - pPpt->beta) * Node[S].W[K]

				+ pPpt->beta * Node[S + 1].W[K]

				- pPpt->alpha * (delt / (Node[S].DELX_R * pPpt->xref)) * (Node[S + 1].F[K] - Node[S].F[K])		// Flux terms

				- pPpt->alpha * delt * (pPpt->beta * Node[S + 1].C[K] + (1 - pPpt->beta) * Node[S].C[K])		// Source terms

				+ Node[S].S[K];																					// Viscosity terms		
		}
	}

	// Next estimate boundary W_pred values
	// ----------------------------------------------------------------------------------------------------
	S = N - 1;
	
	// Set boundary predicted solution equal to existing solution
	// ----------------------------------------------------------------------------------------------------
	Node[S].W_pred[0] = Node[S].W[0];
	Node[S].W_pred[1] = Node[S].W[1];
	Node[S].W_pred[2] = Node[S].W[2];
	//*/
	/*
	// Estimate boundary W_pred values by extrapolating W_pred to the boundary (precarious if gradients are high)
	// ----------------------------------------------------------------------------------------------------
		Node[S].W_pred[0] = (((Node[S-1].W_pred[0] - Node[S-2].W_pred[0])/(Node[S-1].x - Node[S-2].x))
							*(Node[S].x - Node[S-1].x))
							+ Node[S-1].W_pred[0];
		Node[S].W_pred[1] = (((Node[S-1].W_pred[1] - Node[S-2].W_pred[1])/(Node[S-1].x - Node[S-2].x))
							*(Node[S].x - Node[S-1].x))
							+ Node[S-1].W_pred[1];
		Node[S].W_pred[2] = (((Node[S-1].W_pred[2] - Node[S-2].W_pred[2])/(Node[S-1].x - Node[S-2].x))
							*(Node[S].x - Node[S-1].x))
							+ Node[S-1].W_pred[2];
	*/

	// Predictor equation for node geometry values
	// ----------------------------------------------------------------------------------------------------
	
	// Haemo: geometry can vary in time so need to calculated f_pred, d_pred & dddx_pred from predicted soloution vector
	if (pPpt->HAEMODYNAMICS) 
	{
/*
		// Need to calculate f_pred and d_pred first for all nodes in this separate loop
		for (S = 0; S < N; ++S) {
			Node[S].f_pred = Node[S].W_pred[0] / Node[S].rho;
			Node[S].d_pred = sqrt(Node[S].f_pred * 4 / PI);
		}

		// Estimate dddx_pred and dfdx_pred for internal nodes using a central differencing scheme...
		for (S = 1; S < N - 1; ++S) { 
			Node[S].dddx_pred = (Node[S + 1].d_pred - Node[S - 1].d_pred) / (Node[S + 1].x - Node[S - 1].x);
			Node[S].dfdx_pred = (Node[S + 1].f_pred - Node[S - 1].f_pred) / (Node[S + 1].x - Node[S - 1].x);
		}
		// ...and hence must set end nodes separately
		S = 0;
		Node[S].dddx_pred = (Node[S + 1].d_pred - Node[S].d_pred) / (Node[S + 1].x - Node[S].x);
		Node[S].dfdx_pred = (Node[S + 1].f_pred - Node[S].f_pred) / (Node[S + 1].x - Node[S].x);
		S = N - 1;
		Node[S].dddx_pred = (Node[S].d_pred - Node[S - 1].d_pred) / (Node[S].x - Node[S - 1].x);
		Node[S].dfdx_pred = (Node[S].f_pred - Node[S - 1].f_pred) / (Node[S].x - Node[S - 1].x);

		S = 0; Node[S].dddx_pred = (Node[S + 1].d_pred - Node[S].d_pred) / (Node[S + 1].x - Node[S].x);
		
		for (S = 1; S < N - 1; ++S) Node[S].dddx_pred = (Node[S].d_pred - Node[S - 1].d_pred) / (Node[S].x - Node[S - 1].x);

		for (S = 0; S < N; ++S) Node[S].dfdx_pred = Node[S].dddx_pred * ((PI / 2) * Node[S].d_pred); // f = (pi/4)*d^2 so df/dd = (pi/2)*d then df/dx = dd/dx * df/dd = dd/dx * (pi/2)*d
	
for (S = 0; S < N - 1; ++S) {
Node[S].f_pred = pPpt->beta * Node[S + 1].f + (1 - pPpt->beta) * Node[S].f;
Node[S].dfdx_pred = pPpt->beta * Node[S + 1].dfdx + (1 - pPpt->beta) * Node[S].dfdx;
}
S = N - 1;
Node[S].f_pred = Node[S].f;
Node[S].dfdx_pred = Node[S].dfdx;

//*/

// Need to calculate f_pred and d_pred first for all nodes in this separate loop
		for (S = 0; S < N; ++S) {
			Node[S].f_pred = Node[S].W_pred[0] / Node[S].rho;
			Node[S].d_pred = sqrt(Node[S].f_pred * 4 / PI);

			//For f0 d0 which are invariant in time only beta needs to be considered
			Node[S].f0_pred = pPpt->beta * Node[S + 1].f0 + (1 - pPpt->beta) * Node[S].f0;
			Node[S].d0_pred = sqrt(Node[S].f0_pred * 4 / PI);
		}
		//Must set the value seporately for the end node the predicted node lies outside the domain for beta > 0
		S = N - 1;
		Node[S].d_pred = Node[S].d;
		Node[S].f_pred = Node[S].f;
		Node[S].d0_pred = Node[S].d0;
		Node[S].f0_pred = Node[S].f0;


		// Estimate dddx_pred, dd0dx_pred and dfdx_pred for internal nodes using a central differencing scheme...
		for (S = 1; S < N - 1; ++S) {
			Node[S].dddx_pred = (Node[S + 1].d_pred - Node[S - 1].d_pred) / (Node[S + 1].x - Node[S - 1].x);
			Node[S].dfdx_pred = (Node[S + 1].f_pred - Node[S - 1].f_pred) / (Node[S + 1].x - Node[S - 1].x);
			Node[S].dd0dx_pred = (Node[S + 1].d0_pred - Node[S - 1].d0_pred) / (Node[S + 1].x - Node[S - 1].x);

		}
		// ...and hence must set end nodes separately
		S = 0;
		Node[S].dddx_pred = (Node[S + 1].d_pred - Node[S].d_pred) / (Node[S + 1].x - Node[S].x);
		Node[S].dd0dx_pred = (Node[S + 1].d0_pred - Node[S].d0_pred) / (Node[S + 1].x - Node[S].x);
		Node[S].dfdx_pred = (Node[S + 1].f_pred - Node[S].f_pred) / (Node[S + 1].x - Node[S].x);
		S = N - 1;
		Node[S].dddx_pred = (Node[S].d_pred - Node[S - 1].d_pred) / (Node[S].x - Node[S - 1].x);
		Node[S].dd0dx_pred = (Node[S].d0_pred - Node[S - 1].d0_pred) / (Node[S].x - Node[S - 1].x);
		Node[S].dfdx_pred = (Node[S].f_pred - Node[S - 1].f_pred) / (Node[S].x - Node[S - 1].x);
	}
	else // Assume geometry is invariant in time, hence need only consider beta (i.e., in space)
	{ 	
		for (S = 0; S < N - 1; ++S) { // For appropriate nodes

			Node[S].f_pred = pPpt->beta * Node[S + 1].f + (1 - pPpt->beta) * Node[S].f;
			//Node[S].f_pred = Node[S+1].f;
			//Node[S].f_pred = Node[S].f;
			//Node[S].f_pred = (PI/4)*pow(d_odd,2) + Node[S].dfdx*(pPpt->beta*Node[S+1].x + (1-pPpt->beta)*Node[S].x); // For constant dfdx

			Node[S].dfdx_pred = pPpt->beta * Node[S + 1].dfdx + (1 - pPpt->beta) * Node[S].dfdx;
		}

		S = N - 1;
		//Node[S].d_pred = Node[S].d;
		Node[S].f_pred = Node[S].f;
		Node[S].dfdx_pred = Node[S].dfdx;
	}

	// Corrector step
	// ----------------------------------------------------------------------------------------------------
	
	double rho_pred_temp;
	
	// Derive corresponding predicted flux vector F_pred from predicted solution vector W_pred
	// ----------------------------------------------------------------------------------------------------
	for (S = 0; S < N; ++S) { // Really for interior nodes only, as W_pred is strictly only known for the interior
				
		Node[S].F_pred[0] = Node[S].W_pred[1];	
		
		k = pPpt->gammaAir(Node[S].T);
/*
		if (pPpt->HAEMODYNAMICS) {

			rho_pred_temp = Node[S].W_pred[0] / Node[S].f_pred;
			//rho_pred_temp = Node[S].W_pred[0] / Node[S].f0;

			f_of_R0 = pPpt->k1 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2) + pPpt->k3;

			//Node[S].B_pred = (pPpt->k1 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2) + pPpt->k3) / rho_pred_temp * sqrt(Node[S].f0 * Node[S].f_pred);
			//Node[S].B_pred = f_of_R0 / rho_pred_temp * sqrt(Node[S].f0 * Node[S].f_pred);
			//Node[S].B_pred = (Node[S].f0 / rho_pred_temp) * ((2 / 9) * f_of_R0 * (1 - pow(Node[S].f_pred / Node[S].f0, 3)) + (Node[S].p0_dash * pPpt->PREF * 1e5));
			
			Node[S].B_pred = (Node[S].f0 / rho_pred_temp) * ((2 / 9) * f_of_R0 * (pow(Node[S].f_pred / Node[S].f0, 3) - 1) + (Node[S].p0_dash * pPpt->PREF * 1e5));

			//Node[S].B_pred = (Node[S].f0 / rho_pred_temp) * (Node[S].p0_dash * pPpt->PREF * 1e5);

			//Node[S].B_pred = (Node[S].f0 / rho_pred_temp) * (Node[S].p0_dash * pPpt->PREF * 1e5);
// f_pred is in here

			Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] 
				+ Node[S].B_pred * rho_pred_temp;

			//Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0]
			//	+ (Node[S].f0) * (Node[S].p0_dash * pPpt->PREF * 1e5);
		}
		else {
			Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0]
				+ (k - 1) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0]));
		}

cout << "pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] = " << pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] << endl;
cout << "Node[S=" << S << "].B_pred * rho_pred_temp = " << Node[S].B_pred * rho_pred_temp << endl;
cout << "(k - 1) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0])) = " << (k - 1) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0])) << endl;
cout << "Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] + Node[S].B_pred * rho_pred_temp = " << pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] + Node[S].B_pred * rho_pred_temp << endl;
cout << "Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] + (k - 1) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0])) = " << pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] + (k - 1) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0])) << endl;
cout << endl;

//Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0];
//Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] + (k - 1) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0]));
//Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] + Node[S].B_pred * rho_pred_temp;
//Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] + (k - 1) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0])) + Node[S].B_pred * rho_pred_temp;

*/		

		if (pPpt->HAEMODYNAMICS) {

			rho_pred_temp = Node[S].W_pred[0] / Node[S].f_pred; //Im not sure this is nessassary as rho should be constant... tjf

			Node[S].b_pred = 4.0 / 3.0 * (pPpt->k1 * pow(constant_e, pPpt->k2 * Node[S].d0_pred / 2.0) + pPpt->k3);
			Node[S].dbdr0_pred = 4.0 / 3.0 * (pPpt->k1 * pPpt->k2 * pow(constant_e, pPpt->k2 * Node[S].d0_pred / 2.0));

			Node[S].B_pred = Node[S].b_pred / Node[S].rho * sqrt(Node[S].f0_pred * Node[S].f_pred) - Node[S].b_pred / Node[S].rho * Node[S].f0_pred;		//Calculate B_pred at this timestep

			//cout << Node[S].B_pred << "\t" << Node[S].f_pred << "\t" << Node[S].f0_pred << "\n";

			Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0] + Node[S].B_pred * Node[S].rho; //Calculate the flux
		}
		else {
			Node[S].F_pred[1] = pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0]
				+ (k - 1) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0]));
		}

		Node[S].F_pred[2] = (Node[S].W_pred[1] * Node[S].W_pred[2] / Node[S].W_pred[0])
			+ (Node[S].W_pred[1] / Node[S].W_pred[0]) * (k - 1) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0]));
	}
	// Calculate predicted source vector C_pred via temporary primitives derived from predicted solution vector W_pred
	// ----------------------------------------------------------------------------------------------------
	if(pPpt->SOURCES){
		double rho_pred_temp, p_pred_temp, u_pred_temp, T_pred_temp, mu_air_pred_temp, d_pred_temp, G_pred;
		double Re_pred_temp, ff_pred_temp, ff_pred_temp_C_p; 
		double hc_pred, q_pred, w_pred;

		for (S = 0; S < N - 1; ++S) { // For interior nodes, AND S==0 (C at S==0 required for corrector)
			k = pPpt->gammaAir(Node[S].T);

			// Continuity source term is always zero
			// ----------------------------------------------------------------------------------------------------
			Node[S].C_pred[0] = 0;

			// Momentum source due to change in area and friction (-p*dF/dx + rho*G*F)
			// ----------------------------------------------------------------------------------------------------
			rho_pred_temp = Node[S].W_pred[0] / Node[S].f_pred;
			p_pred_temp = ((k - 1) / Node[S].f_pred) * (Node[S].W_pred[2] - 0.5 * (pow(Node[S].W_pred[1], 2) / Node[S].W_pred[0]));
			u_pred_temp = Node[S].W_pred[1] / Node[S].W_pred[0];
			T_pred_temp = p_pred_temp / (rho_pred_temp * pPpt->R_air);
			mu_air_pred_temp = pPpt->ViscosityAir(T_pred_temp);

			// Circular section pipe
			// ----------------------------------------------------------------------------------------------------
			d_pred_temp = sqrt(Node[S].f_pred * 4 / PI);
			Re_pred_temp = rho_pred_temp * fabs(u_pred_temp) * d_pred_temp / mu_air_pred_temp;
			ff_pred_temp = FrictionFactor(pPpt, d_pred_temp, Re_pred_temp, "W_alpha_beta", CFTRANS);
			ff_pred_temp_C_p = (d_pred_temp / (((Node[S].DELX_L + Node[S].DELX_R) / 2) * pPpt->xref)) * C_p * CPTRANS;
			G_pred = 0.5 * u_pred_temp * fabs(u_pred_temp) * (ff_pred_temp + ff_pred_temp_C_p) * (4 / d_pred_temp);

			if (pPpt->HAEMODYNAMICS) {
/*
				//Hemodynamics momentum sources (cu + dB/dr0 * dr0/dx)
				//----------------------------------------------------------------------------------------------------
				double delta = pPpt->BLThickness(d_pred_temp);//double delta = 0.2 * d_pred_temp; //Boundary layer thickness, value entirely arbatory at the moment need to find a proper model for this!!
				double viscosity = pPpt->ViscosityBloodFunc(Node[S].T);//double viscosity = 0.003; //Fairly typical value in Pa.s.  Need to find a better model as non newtonian.
				rho_pred_temp = Node[S].W_pred[0] / Node[S].f_pred;

				Node[S].C_pred[1] = -PI * d_pred_temp * u_pred_temp * viscosity / delta
					+ (2 * sqrt(Node[S].f_pred) * (sqrt(PI) * 4 / 3 * (pPpt->k1 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2) + pPpt->k3) + sqrt(Node[S].f0) * 4 / 3 * (pPpt->k1 * pPpt->k2 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2)))
						- Node[S].f_pred * 4 / 3 * (pPpt->k1 * pPpt->k2 * pow(constant_e, pPpt->k2 * Node[S].d0 / 2))) * Node[S].dd0dx / 2;
//*/
				// Hemodynamics momentum sources (cu + dB/dr0 * dr0/dx)
				//----------------------------------------------------------------------------------------------------
				double delta = pPpt->BLThickness(d_pred_temp);//double delta = 0.2 * d_pred_temp; //Boundary layer thickness, value entirely arbatory at the moment need to find a proper model for this!!
				double viscosity = pPpt->ViscosityBloodFunc(Node[S].T);//double viscosity = 0.003; //Fairly typical value in Pa.s.  Need to find a better model as non newtonian.
				rho_pred_temp = Node[S].W_pred[0] / Node[S].f_pred;

				Node[S].C_pred[1] = -PI * d_pred_temp * u_pred_temp * viscosity / delta +
					((2 * sqrt(Node[S].f_pred * Node[S].f0_pred) - Node[S].f_pred - Node[S].f0_pred) * Node[S].dbdr0_pred + 2 * sqrt(PI) * Node[S].b_pred * (sqrt(Node[S].f_pred) - sqrt(Node[S].f0_pred))) * Node[S].dd0dx_pred / 2;

			}
			else { // Otherwise use  (-p*dF/dx + rho*G*F)
				Node[S].C_pred[1] = (-1) * p_pred_temp * Node[S].dfdx_pred + Node[S].W_pred[0] * G_pred;
			}

//Node[S].C_pred[1] = (-1) * p_pred_temp * Node[S].dfdx_pred + Node[S].W_pred[0] * G_pred;

			// Source due to work input and heat transfer (rho*w*F - rho*q*F)
			// ----------------------------------------------------------------------------------------------------
			
			// Source due to heat transfer (... - rho*q*F)
			// ----------------------------------------------------------------------------------------------------
			hc_pred = ConvectiveHTCoefficient(pPpt, ff_pred_temp, rho_pred_temp, u_pred_temp, T_pred_temp, Re_pred_temp, Node[S].DELX_R*pPpt->xref);
			q_pred = ((4*hc_pred)/(rho_pred_temp*d_pred_temp))*(Tw - T_pred_temp); // General equation for heat transfer rate per unit mass
				
			// Source due to work input (rho*w*F - ...)
			// ----------------------------------------------------------------------------------------------------
			w_pred = 0;
			
			Node[S].C_pred[2] = Node[S].W_pred[0]*(w_pred - q_pred);
		}
	}
	else{ // No source terms except those due to area change
		for(S=0; S<N; ++S){ // For all nodes
			k = pPpt->gammaAir(Node[S].T);
			Node[S].C_pred[0] = 0;
			double p_pred_temp = ((k-1)/Node[S].f_pred)*( Node[S].W_pred[2] - 0.5*(pow(Node[S].W_pred[1], 2)/Node[S].W_pred[0]) );		
			Node[S].C_pred[1] = (-1)*p_pred_temp*Node[S].dfdx_pred;
			Node[S].C_pred[2] = 0;
		}
	}

	// Predicted artificial viscosity S_pred based on predicted pressure values, hence requires a W_pred at boundaries
	// ----------------------------------------------------------------------------------------------------
	if(pPpt->VISC){
		double* p_visc_pred;
		p_visc_pred = new double [N];
		for(S=0; S<N; ++S){
			k = pPpt->gammaAir(Node[S].T);
			p_visc_pred[S] = ((k-1)/Node[S].f)*( Node[S].W_pred[2] - 0.5*(pow(Node[S].W_pred[1], 2)/Node[S].W_pred[0]) );
		}

		for(S=1; S<N-1; ++S){
			// Calculate artificial viscosity S - note uses both forward and rearward values
			for(K=0; K<3; ++K)
				Node[S].S_pred[K] = pPpt->Cx*(fabs((p_visc_pred[S+1]) - 2*(p_visc_pred[S]) + (p_visc_pred[S-1]))
										/((p_visc_pred[S+1]) + 2*(p_visc_pred[S]) + (p_visc_pred[S-1])))
									*(Node[S+1].W_pred[K] - 2*Node[S].W_pred[K] + Node[S-1].W_pred[K]);
		}
		S=0; for(K=0; K<3; ++K) Node[S].S_pred[K] = 0;
		S=N-1; for(K=0; K<3; ++K) Node[S].S_pred[K] = 0;
		delete [] p_visc_pred;
	}
	else for(S=0; S<N; ++S) for(K=0; K<3; ++K) Node[S].S_pred[K] = 0;

	// Save previous solution
	// ----------------------------------------------------------------------------------------------------
	for(S=0; S<N; ++S) for(K=0; K<3; ++K) Node[S].W_prev[K] = Node[S].W[K];

	// Corrector equation - calculate new solution vector W(i, n + 1)
	// ----------------------------------------------------------------------------------------------------
	if(pPpt->TVD) TVD(pPpt, delt);
	else{
		for(S=1; S<N-1; ++S){
			for(K=0; K<3; ++K){
				Node[S].W[K] = Node[S].W[K]

							- (delt/(2*pPpt->alpha*(Node[S].DELX_R*pPpt->xref)))	// NB: DELX_R
								* (pPpt->alpha - pPpt->beta) 
								* (Node[S+1].F[K] - Node[S].F[K])

							+ (delt/(2*pPpt->alpha*(Node[S].DELX_L*pPpt->xref)))	// NB: DELX_L
								* (1 - pPpt->alpha - pPpt->beta) 
								* (Node[S].F[K] - Node[S-1].F[K])

							- (delt/(2*pPpt->alpha*(Node[S].DELX_L*pPpt->xref)))	// NB: DELX_L, but...
								* (Node[S].F_pred[K] - Node[S-1].F_pred[K])			// ...DELX_R used to obtain F_pred

							+ (delt/(2*pPpt->alpha))*(
								 (
									(1 - pPpt->alpha - pPpt->beta)*Node[S+1].C[K]
								  +	(2*pPpt->beta - 1)*Node[S].C[K]
								  + (1 - pPpt->alpha - pPpt->beta)*Node[S-1].C[K]
								 )
								-
								 ((1 - pPpt->beta)*Node[S].C_pred[K] + pPpt->beta*Node[S-1].C_pred[K])
													)
								
								+ Node[S].S_pred[K];
			}
		}
	}

}

void CPipe::TVD(CProperties* pPpt, double del_t)
// ==================================================================================================== 
// Additional code to determine TVD flux limiting terms for a W_alpha_beta method 
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".TVD\n");}

	int i, K, S;

	//	Copy solution W in order to calculate differences
	// ----------------------------------------------------------------------------------------------------
	double** W_old;
	W_old = new double* [N];
	for(S=0; S<N; ++S){
		W_old[S] = new double [3];
		for(K=0; K<3; ++K) W_old[S][K] = Node[S].W[K];	
	}
	double** delW; delW = new double* [4]; for(i=0; i<4; ++i) delW[i] = new double [3];
	double* r; r = new double [4];
	double k, u_temp, rho_temp, p_temp, T_temp, a_temp;
	double lambda_1, lambda_2, lambda_3, lambda_max;
	double del_x = xmesh;

	for(S=1; S<N-1; ++S){ // r values are only required for internal mesh points (the solution can be propagated to the boundaries using MMOC)
		// Calculate the eigenvalues at each node; lambda_1 = u + a, lambda_2 = u - a, lambda_3 = u;
		// ----------------------------------------------------------------------------------------------------
		k = pPpt->gammaAir(Node[S].T);
		u_temp = (Node[S].W[1]/Node[S].W[0]);
		rho_temp = Node[S].W[0]/Node[S].f;
		p_temp = (((k-1)/Node[S].f)*( Node[S].W[2] - 0.5*(pow(Node[S].W[1], 2)/Node[S].W[0]) ));
		T_temp = (p_temp)/(rho_temp*pPpt->R_air);
		k = pPpt->gammaAir(T_temp);
		a_temp = sqrt(k*pPpt->R_air*T_temp);
		lambda_1 = u_temp + a_temp; lambda_2 = u_temp - a_temp; lambda_3 = u_temp;

		if(fabs(lambda_1)>fabs(lambda_2)){
			if(fabs(lambda_1)>fabs(lambda_3)) lambda_max = lambda_1;
			else lambda_max = lambda_3;
		}
		else{
			if(fabs(lambda_2)>fabs(lambda_3)) lambda_max = lambda_2;
			else lambda_max = lambda_3;
		}

		// Local artficial viscosity
		// ----------------------------------------------------------------------------------------------------
		for(i=0; i<4; ++i){
			for(K=0; K<3; ++K){
				if(S==1 && i==0){ // At the first internal mesh point assume required virtual gradient is identical to that between first real nodes
					delW[i][K] = W_old[S - 1 + i + 1][K] - W_old[S - 2 + i + 1][K];	// Just use W_old[S-1.5] = W_old[S-0.5]
				}
				else{
					if(S==N-2 && i==3) // At the last internal mesh point assume required virtual gradient is identical to that between last real nodes
						delW[i][K] = W_old[S - 1 + i - 1][K] - W_old[S - 2 + i - 1][K];	// Just use W_old[S+1.5] = W_old[S+0.5]
					else 
						delW[i][K] = W_old[S - 1 + i][K] - W_old[S - 2 + i][K];
				}
			}
		}
		// delW[0][k] = W_old[S-1] - W_old[S-2] = W_old[S-1.5]
		// delW[1][k] = W_old[S] - W_old[S-1] = W_old[S-0.5]
		// delW[2][k] = W_old[S+1] - W_old[S] = W_old[S+0.5]
		// delW[3][k] = W_old[S+2] - W_old[S+1] = W_old[S+1.5]

		r[0] = (delW[0][0]*delW[1][0] + delW[0][1]*delW[1][1] + delW[0][2]*delW[1][2])
				/(pow(delW[1][0],2) + pow(delW[1][1],2) + pow(delW[1][2],2)); // r+,i-1

		r[1] = (delW[1][0]*delW[2][0] + delW[1][1]*delW[2][1] + delW[1][2]*delW[2][2])
				/(pow(delW[1][0],2) + pow(delW[1][1],2) + pow(delW[1][2],2)); // r-,i

		r[2] = (delW[1][0]*delW[2][0] + delW[1][1]*delW[2][1] + delW[1][2]*delW[2][2])
				/(pow(delW[2][0],2) + pow(delW[2][1],2) + pow(delW[2][2],2)); // r+,i

		r[3] = (delW[2][0]*delW[3][0] + delW[2][1]*delW[3][1] + delW[2][2]*delW[3][2])
				/(pow(delW[2][0],2) + pow(delW[2][1],2) + pow(delW[2][2],2)); // r-,i+1
		
		for(K=0; K<3; ++K){
			Node[S].W[K] = Node[S].W[K]

						- (del_t/(2*pPpt->alpha*(Node[S].DELX_R*pPpt->xref)))	// NB: DELX_R
							* (pPpt->alpha - pPpt->beta) 
							* (Node[S+1].F[K] - Node[S].F[K])

						+ (del_t/(2*pPpt->alpha*(Node[S].DELX_L*pPpt->xref)))	// NB: DELX_L
							* (1 - pPpt->alpha - pPpt->beta) 
							* (Node[S].F[K] - Node[S-1].F[K])

						- (del_t/(2*pPpt->alpha*(Node[S].DELX_L*pPpt->xref)))	// NB: DELX_L, but...
							* (Node[S].F_pred[K] - Node[S-1].F_pred[K])			// ...DELX_R used to obtain F_pred

						- (del_t/(2*pPpt->alpha))*(
							-(
								(1 - pPpt->alpha - pPpt->beta)*Node[S+1].C[K]
							  +	(2*pPpt->beta - 1)*Node[S].C[K]
							  + (1 - pPpt->alpha - pPpt->beta)*Node[S-1].C[K]
							 )
							+
							 ((1 - pPpt->beta)*Node[S].C_pred[K] + pPpt->beta*Node[S-1].C_pred[K])
												)
							
						+ Node[S].S_pred[K]

						// Append the dissipative term
						// ----------------------------------------------------------------------------------------------------
						+ ( G_flux_limiter(pPpt, lambda_max, del_t, del_x, r[2]) 
							+ G_flux_limiter(pPpt, lambda_max, del_t, del_x, r[3]) )*delW[2][K]
						- ( G_flux_limiter(pPpt, lambda_max, del_t, del_x, r[0]) 
							+ G_flux_limiter(pPpt, lambda_max, del_t, del_x, r[1]) )*delW[1][K];
		}
	}
	for(S=0; S<N; ++S) delete [] W_old[S]; delete [] W_old;
	for(i=0; i<4; ++i) delete [] delW[i]; delete [] delW;
	delete [] r;
}

double CPipe::G_flux_limiter(CProperties* pPpt, double lambda_max, double del_t, double del_x, double r)
// ====================================================================================================
// The Davis flux limiter; places upper limit on Courant number of 0.7
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".G_flux_limiter\n");}

	double v, phi_of_r, C_of_v;
	v = fabs(lambda_max)*del_t/del_x;

	// The Davis defined flux limiter:
	if(r>0){
		// phi_of_r is min(2r, 1)
		if(2*r < 1)	phi_of_r = 2*r;
		else phi_of_r = 1;
	}
	else phi_of_r = 0;

	// C_of_v defined as:
	if(v<=0.5) C_of_v = v*(1 - v);
	else C_of_v = 0.25;

	return 0.5*C_of_v/*0--0.25*/*(1 - phi_of_r)/*0--1*/; // 0--0.125
}

void CPipe::W_alpha_beta_derive_lambdas(CProperties* pPpt)
// ==================================================================================================== 
// Generate characteristics following execution of W_alpha_beta scheme; interior nodes only
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".W_alpha_beta_derive_lambdas\n");}

	double k, rho_temp, U_temp, p_dash_temp, A_temp, T_temp;
	for(int S=1; S<N-1; ++S)
	{
		k = pPpt->gammaAir(Node[S].T);
		rho_temp = Node[S].W[0]/Node[S].f;
		U_temp = (Node[S].W[1]/Node[S].W[0])/AREF;	
		p_dash_temp = (	(((k-1)/Node[S].f)*( Node[S].W[2] - 0.5*(pow(Node[S].W[1], 2)/Node[S].W[0]) ))
						/1e5 // To get bar
						)/pPpt->PREF;	// Both in bar

		// Make use of entropy level derived from pathlines
		A_temp = pow(p_dash_temp, (k-1)/(2*k))*Node[S].AA[R+1];
		T_temp = pow(A_temp, 2)*(EX ? pPpt->TREFe : pPpt->TREFi);
		k = pPpt->gammaAir(T_temp);

		// Contruct characteristics at new time level
		Node[S].CL1[R+1] = A_temp + ((k-1)/2)*U_temp;
		Node[S].CL2[R+1] = A_temp - ((k-1)/2)*U_temp;
	}
}

void CPipe::W_alpha_beta_prep_bcs(CProperties* pPpt)
// ============================================================ //
// Prepare boundary nodes following interior propagation		//
// ============================================================ //
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".W_alpha_beta_prep_bcs\n");}

	int S, K;

	// Extrapolate solution vector to the boundary nodes
	// -------------------------------------------------
	S=0;
	for(K=0; K<3; ++K) Node[S].W[K] = (Node[S+1].W[K]/Node[S+1].f)*Node[S].f; // Set equal to adjacent node, remembering to correct for a possible difference in node area
	
	S=N-1;
	for(K=0; K<3; ++K) Node[S].W[K] = (Node[S-1].W[K]/Node[S-1].f)*Node[S].f; // Set equal to adjacent node, remembering to correct for a possible difference in node area

	// Derive characteristics and entropy values at ALL nodes at the new time level
	// ----------------------------------------------------------------------------
	double k, rho_temp, p_temp, u_temp, T_temp, p_dash_temp, U_temp, A_temp;
	for(S=0; S<N; ++S)
	{
		rho_temp = Node[S].W[0]/Node[S].f;
		k = pPpt->gammaAir(Node[S].T);
		p_temp = ((k-1)/Node[S].f)*( Node[S].W[2] - 0.5*(pow(Node[S].W[1], 2)/Node[S].W[0]) );
		u_temp = Node[S].W[1]/Node[S].W[0];
		T_temp = p_temp/(rho_temp*pPpt->R_air);
		k = pPpt->gammaAir(T_temp);
		p_dash_temp = (p_temp/1e5)/pPpt->PREF;
		U_temp = u_temp/AREF;
		A_temp = sqrt(T_temp/(EX ? pPpt->TREFe : pPpt->TREFi));
	
		// Construct characteristics and entropy level
		Node[S].CL1[R+1] = A_temp + ((k-1)/2)*U_temp;
		Node[S].CL2[R+1] = A_temp - ((k-1)/2)*U_temp;
		Node[S].AA[R+1] = A_temp/pow(p_dash_temp, (k-1)/(2*k));
	}
}

double CPipe::phi(CProperties* pPpt, double x, double epsilon)
// ============================================================ //
// Used in Glimm's method 										//
// ----------------------										//
// Implemented from Sod, JCP 27, 1-31 (1978)					//
// ============================================================ //
{
	double gamma = pPpt->gammaAir(300);
	if(fabs(1 - x) >= epsilon)
	{
		if(x >= 1) return sqrt( ((gamma + 1)/2)*x + (gamma - 1)/2 );
		else return ( ((gamma - 1)/2)*(1 - x) ) / ( sqrt(gamma)*(1 - pow(x, (gamma - 1)/(2*gamma))) );
	}
	else return sqrt(gamma);
}

void CPipe::RecordTappings(CProperties* pPpt, int timestep, double time)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".RecordTappings\n");}

	int n, r, d, S, columncounter;
	bool located;

//	if((timestep % sample_factor) == 0)
//	{
		for(n=0; n<ntappings; ++n)
		{
			// Interpolate at each tapping location
			if(Measurements[n]<0 || Measurements[n]>1) {cout << "Measurement location outside of pipe\n"; break;}

			// Locate the two nodes either side of the location of this tapping
			S=0;
			if(N>1)
			{
				located = false;
				while(!located)
				{
					//if(Node[S+1].X_REAL > (Measurements[n]*XPIPE_REAL) ||
					//		fabs(Node[S+1].X_REAL - (Measurements[n]*XPIPE_REAL)) < 1e-6)
					if(Node[S+1].x > (Measurements[n]*length) ||
						fabs(Node[S+1].x - (Measurements[n]*length)) < 1e-6)
						located = true;
					else ++S;
				};
		
				// Then the location is between nodes S and S+1
				//cout << "Measurements[n] = " << Measurements[n] << endl;
				//cout << "XPIPE_REAL = " << XPIPE_REAL << endl;
				//cout << "Node[S].X_REAL = " << Node[S].X_REAL << endl;
				
				//			MeasureNode[n] =  Node[S] +
				//								(Node[S+1] - Node[S])*
				//								(((Measurements[n]*XPIPE_REAL) - Node[S].X_REAL)
				//								/(Node[S+1].X_REAL - Node[S].X_REAL));

				// Then the location is between nodes S and S+1
				MeasureNode[n] =  Node[S] +
					(Node[S+1] - Node[S])*
					(((Measurements[n]*length) - Node[S].x)
					/(Node[S+1].x - Node[S].x));
				
			}
			else // For single node or joiner pipes
			{
				MeasureNode[n] =  Node[S];			
			}
			
			// Now update derived properties for this node
			MeasureNode[n].Re = (MeasureNode[n].rho*(fabs(MeasureNode[n].U)*AREF)*MeasureNode[n].d)/pPpt->ViscosityAir(MeasureNode[n].T);
		}
//	}
	
	if((timestep % sample_factor) == 0)
	{
		// Check whether we are above the limit of number of points to record
		if(rowcounter>=this->max_pts)
		{
			// Go through results array and remove every second row, shift up
			for(n=0; n<ntappings; ++n)
				for(r=1; r <= int((this->max_pts-1)/2); ++r)
					for(d=0; d<= this->num_props_measured+1; ++d) Results[n][r][d] = Results[n][2*r][d];

			rowcounter = int((this->max_pts-1)/2) + 1;
			sample_factor *= 2; // Reduce sampling rate by half
		}

		for(n=0; n<ntappings; ++n)
		{
			// Then extract the desired properties from the interpolated node
			columncounter=0;
			if(this->num_props_measured>0)
			{
//				Results[n][rowcounter][columncounter] = MeasureNode[n].X_REAL*pPpt->xref;
				Results[n][rowcounter][columncounter] = MeasureNode[n].x;
				++columncounter;
			}
			if(this->num_props_measured>0)
			{
				Results[n][rowcounter][columncounter] = time;
				++columncounter;
			}
			if(DIAMETER && columncounter <= this->num_props_measured + 1) {
				Results[n][rowcounter][columncounter] = MeasureNode[n].d;
				++columncounter;
			}
			if(AREA && columncounter <= this->num_props_measured + 1) {
				Results[n][rowcounter][columncounter] = MeasureNode[n].f;
				++columncounter;
			}
			if(STATIC_PRESSURE && columncounter<=this->num_props_measured + 1) {
				Results[n][rowcounter][columncounter] = MeasureNode[n].p_dash*pPpt->PREF;
				++columncounter;
			}
			if(TEMPERATURE && columncounter<=this->num_props_measured + 1) {
				Results[n][rowcounter][columncounter] = MeasureNode[n].T; 
				++columncounter;
			}
			if(REYNOLDS_NO && columncounter<=this->num_props_measured+1) {
				Results[n][rowcounter][columncounter] = MeasureNode[n].Re;
				++columncounter;
			}
		}
		++ rowcounter;
	}

	for(int m=0; m<ntappings; ++m)
	{
		Pressure[m] = MeasureNode[m].p_dash*pPpt->PREF;
		Temperature[m] = MeasureNode[m].T;
		Velocity[m] = MeasureNode[m].U*AREF;
		MassFlowRate[m] = MeasureNode[m].mdot;
		Re[m] = MeasureNode[m].Re;
//		MassFlowRate[m] = fabs(MeasureNode[m].mdot);
/*
		double mdot_temp, T1_temp, C1_temp, p1_temp, p2_temp;
		double T01_temp = 0;
		double p01_temp = 0;
*/
/*
		if(this->end==ODD)
		{
			mdot_temp = this->pPipe->Node[0].mdot;
			T1_temp = this->pPipe->Node[0].T;
			C1_temp = this->pPipe->Node[0].U*AREF;
			p1_temp = this->pPipe->Node[0].p_dash*pPpt->PREF;
//			p2_temp = this->pPipe->Node[this->pPipe->N-1].p_dash*pPpt->PREF;
		}
		else
		{
			mdot_temp = this->pPipe->Node[this->pPipe->N-1].mdot;
			T1_temp = this->pPipe->Node[this->pPipe->N-1].T;
			C1_temp = this->pPipe->Node[this->pPipe->N-1].U*AREF;
			p1_temp = this->pPipe->Node[this->pPipe->N-1].p_dash*pPpt->PREF;
//			p2_temp = this->pPipe->Node[0].p_dash*pPpt->PREF;
		}
*/
/*
		T01_temp = T1_temp*(1 + ((pPpt->gammaAir()-1)/2)*(pow(C1_temp,2)/(pPpt->gammaAir()*pPpt->R_air*T1_temp)));
		p01_temp = p1_temp*( pow(T01_temp/T1_temp, pPpt->gammaAir()/(pPpt->gammaAir()-1)) );
*/
		PR[m] = 1;//p01_temp/p2_temp;
		MFP[m] = 0;//mdot_temp*sqrt(T01_temp)/p01_temp;	
	}
}

void CPipe::ListProperties(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ListProperties\n");}
	
	int i, s;


	if(BUFFER) pPpt->Out(Underline("Buffer pipe", "=", "\t", strDesc));
	else
	{
		if(DAMPER) pPpt->Out(Underline("Damper pipe", "=", "\t", strDesc));
		else pPpt->Out(Underline(Identify(), "=", "\t", strDesc));
	}
	pPpt->Out("\n");

/*
	if(BUFFER) pPpt->Out(Underline("Buffer pipe", "=", "\t"));
	else
	{
		if(DAMPER) pPpt->Out(Underline("Damper pipe", "=", "\t"));
		else pPpt->Out(Underline(Identify(), "=", "\t"));
	}
	pPpt->Out("\n");
*/

	// Propagation method
	// ====================================================================================================
	pPpt->Out(Underline("Propagation method", "-", "\t"));

	// (1) Mesh Method of Characteristics (MMOC)
	// ----------------------------------------------------------------------------------------------------
	if(METHOD==1)
	{
		pPpt->Out("\tMesh Method of Characteristics (");
		if(pPpt->HOMENTROPIC) pPpt->Out("homentropic)\n");

		// Non-homentropic MMOC
		// ----------------------------------------------------------------------------------------------------
		else
		{
			pPpt->Out("non-homentropic)\n");
			pPpt->Out("\tPathlines multiplier, NUM_PATH_MULT\t\t=\t"); pPpt->Out(pPpt->NUM_PATH_MULT); pPpt->Out("\n");
		}
	}
	else
	{
		// (2) W_alpha_beta schemes
		// ----------------------------------------------------------------------------------------------------
		if(METHOD==2)
		{	
			// Character codes: 224 = alpha, 225 = beta
			pPpt->Out("\tW"); pPpt->Out('a'); pPpt->Out(char(225)); pPpt->Out(" scheme ("); pPpt->Out('a'); pPpt->Out("="); pPpt->Out(pPpt->alpha);pPpt->Out(", ");
			pPpt->Out(char(225)); pPpt->Out("="); pPpt->Out(pPpt->beta); pPpt->Out(")"); 

			if(pPpt->alpha==0.5 && pPpt->beta==0.5) pPpt->Out("\t\t\t=\tTwo-step Lax-Wendroff");
			else
			{
				if(pPpt->alpha==1 && pPpt->beta==0)
				{
					pPpt->Out("\t\t\t\t=\tMacCormack method");
					if(pPpt->ALTERNATE_MAC) pPpt->Out(" (alternate)");
				}
				else
				{
					if(fabs(pPpt->alpha-(1 + sqrt(5.)/2))<1e-6 && pPpt->beta==0.5) pPpt->Out("\t\t\t=\tLerat & Peyret \"optimum\"");
					else {pPpt->Out("\t\t\t\t=\tUncategorised W"); pPpt->Out(char(224)); pPpt->Out(char(225)); pPpt->Out(" scheme");}
				}
			}
			pPpt->Out("\n\t\t\t\t\t\t\t\t");
			if(pPpt->SOURCES)
			{
				if(pPpt->WORK) pPpt->Out(" with source terms (inc. work)");
				else pPpt->Out(" with source terms (but no work)");
			}			
			else pPpt->Out(" no source terms (except area variation)");
			pPpt->Out("\n\t\t\t\t\t\t\t\t");
			if(pPpt->TVD) pPpt->Out(" + TVD");
			else pPpt->Out(" no TVD");
			pPpt->Out("\n");
			if(pPpt->COMBINED_WAB_MOC)
			{
				pPpt->Out("\tContinuing to use pathlines, COMBINED_WAB_MOC\t=\t"); pPpt->Out(TrueOrFalse(pPpt->COMBINED_WAB_MOC)); pPpt->Out("\n");
				pPpt->Out("\tPathlines multiplier, NUM_PATH_MULT\t\t=\t"); pPpt->Out(pPpt->NUM_PATH_MULT); pPpt->Out("\n");
			}
			else
			{
				pPpt->Out("\tExtrapolating solution vector, COMBINED_WAB_MOC\t=\t"); pPpt->Out(TrueOrFalse(pPpt->COMBINED_WAB_MOC)); pPpt->Out("\n");
			}

			if(pPpt->VISC)
			{
				pPpt->Out("\tApplying artificial viscosity, VISC\t\t=\t"); pPpt->Out(TrueOrFalse(pPpt->VISC)); pPpt->Out("\n");
				pPpt->Out("\tArtificial viscosity coefficient, Cx\t\t=\t"); pPpt->Out(pPpt->Cx); pPpt->Out("\n");
				pPpt->Out("\tArtificial viscosity coefficient, Cx_alpha\t=\t"); pPpt->Out(pPpt->Cx_alpha); pPpt->Out("\n");
			}
			else {pPpt->Out("\tNo artificial viscosity, VISC\t\t\t=\t"); pPpt->Out(TrueOrFalse(pPpt->VISC)); pPpt->Out("\n");}
		}
		else
		{
			// (3) Filling and emptying
			// ----------------------------------------------------------------------------------------------------
			if(METHOD==3) pPpt->Out("\tFilling and emptying\n");
			else pPpt->Out("unknown method\n");
		}
	}
	pPpt->Out("\tDelay start of propagation by, delay\t\t=\t"); pPpt->Out(delay); pPpt->Out(" s\n");
	pPpt->Out("\n");

	// Geometry
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out(Underline("Geometry", "-", "\t"));
	pPpt->Out("\tPhysical length, length\t\t\t\t=\t"); pPpt->Out(length*1000); pPpt->Out(" mm\n");
	if(LINEAR)
  {
		if(LINEAR_F)
		{
			pPpt->Out("\tLinear pipe area variation:\n");
			pPpt->Out("\t- diameter odd end, d_odd\t\t\t=\t"); pPpt->Out(d_odd*1000); pPpt->Out(" mm\n");
			pPpt->Out("\t- diameter even end, d_even\t\t\t=\t"); pPpt->Out(d_even*1000); pPpt->Out(" mm\n");
		}
		else
		{
			pPpt->Out("\tLinear pipe diameter variation:\n");
			pPpt->Out("\t- diameter odd end, d_odd\t\t\t=\t"); pPpt->Out(d_odd*1000); pPpt->Out(" mm\n");
			pPpt->Out("\t- diameter even end, d_even\t\t\t=\t"); pPpt->Out(d_even*1000); pPpt->Out(" mm\n");
			//pPpt->Out("\t- D(X) = "); pPpt->Out(d_odd); pPpt->Out(" + "); pPpt->Out(C); pPpt->Out("X\n");
		}
		pPpt->Out("\tPipe turning angle, bend_angle\t\t\t=\t"); pPpt->Out(bend_angle); pPpt->Out(" degrees");
		if(bend_angle==0) pPpt->Out(" - straight pipe");
		pPpt->Out("\n");
		pPpt->Out("\tBend pressure loss coefficient due, C_p\t\t=\t"); pPpt->Out(C_p); pPpt->Out("\n");
		//pPpt->Out("\n");	
		if(n_int_points>0)
		{
			//pPpt->Out(Underline("Pipe diameters", "-", "\t"));
			pPpt->Out("\tNo. of internal diameters, n_int_points\t\t=\t"); pPpt->Out(n_int_points); pPpt->Out("\n");
		}
		for(i=0; i<n_int_points+2; ++i)
		{
			pPpt->Out("\t- distance = "); pPpt->Out(xi[i]*1000); pPpt->Out(" mm, diameter = "); pPpt->Out(di[i]*1000); pPpt->Out(" mm\n");
		}
		//pPpt->Out("\n");
	}
	else
	{
		pPpt->Out("\tQuadratic pipe diameter variation:\n");
		pPpt->Out("\t- currently disabled. Exiting.\n");
		exit(1);
		pPpt->Out("\t- dia. odd end, d_odd\t\t\t\t=\t"); pPpt->Out(d_odd*1000); pPpt->Out(" mm\n");
		pPpt->Out("\t- dia. interior point, d_int\t\t\t=\t"); pPpt->Out(d_int*1000); pPpt->Out(" mm\n");
		pPpt->Out("\t- distance of int. point from odd end, x_int\t=\t"); pPpt->Out(x_int*1000); pPpt->Out(" mm\n");
		pPpt->Out("\t- dia. even end, d_even\t\t\t=\t"); pPpt->Out(d_even*1000); pPpt->Out(" mm\n");
		//pPpt->Out("\t- D(X) = "); pPpt->Out(A); pPpt->Out("X^2 + "); pPpt->Out(B); pPpt->Out("X + "); pPpt->Out(C); pPpt->Out("\n");
	}
	pPpt->Out("\tVolume, vol\t\t\t\t\t=\t"); pPpt->Out(vol*1000); pPpt->Out(" litres\n");
	pPpt->Out("\n");

	// Meshing
	// ----------------------------------------------------------------------------------------------------
	if(GLOBAL_MP) pPpt->Out("\tUsing GLOBAL meshing parameters:");
	else pPpt->Out("\tUsing LOCAL meshing parameters factors:");
	pPpt->Out("\n");
	pPpt->Out("\t- Target discretization length, discret\t\t=\t"); pPpt->Out(discret*1000); pPpt->Out(" mm\n");
	pPpt->Out("\t- Minimum number of meshes, min_meshes\t\t=\t"); pPpt->Out(min_meshes); pPpt->Out("\n");
	pPpt->Out("\tActual discretization length, xmesh\t\t=\t"); pPpt->Out(xmesh*1000); pPpt->Out(" mm\n");
	pPpt->Out("\tNumber of meshes, meshes\t\t\t=\t"); pPpt->Out(meshes); pPpt->Out("\n");
	pPpt->Out("\n");

	// Coefficients
	// ----------------------------------------------------------------------------------------------------
	if(GLOBAL_EF) pPpt->Out("\tUsing GLOBAL roughness height and enhancement factors:");
	else pPpt->Out("\tUsing LOCAL roughness height and enhancements factors:");
	pPpt->Out("\n");
	pPpt->Out("\t- Roughness height, epsilon\t\t\t=\t"); pPpt->Out(epsilon*1000); pPpt->Out(" mm\n");
	pPpt->Out("\t- Friction enhancement factor, CFTRANS\t\t=\t"); pPpt->Out(CFTRANS); pPpt->Out("\n");
	pPpt->Out("\t- Heat transfer enhancement factor, HGTRANS\t=\t"); pPpt->Out(HGTRANS); pPpt->Out("\n");
	pPpt->Out("\t- Pressure loss enhancement factor, CPTRANS\t=\t"); pPpt->Out(CPTRANS); pPpt->Out("\n");
	pPpt->Out("\n");	

	// Initialization
	// ----------------------------------------------------------------------------------------------------
	pPpt->Out("\tWall temperature, Tw\t\t\t\t=\t"); pPpt->Out(Tw); pPpt->Out(" K\n\n");
	if(nsections>0)
	{
		pPpt->Out(Underline("Divisions", "-", "\t"));
		pPpt->Out("\tNo. of sections (nsections)\t\t\t=\t"); pPpt->Out(nsections); pPpt->Out("\n");
		if(nsections>1) pPpt->Out("\tAs a fraction of physical length from odd end:\n");
	}
	for(s=0; s<nsections-1; ++s)
	{
		pPpt->Out("\t- division "); pPpt->Out(s+1); pPpt->Out(" of "); pPpt->Out(nsections-1); pPpt->Out(" @ "); pPpt->Out(xri[s]); pPpt->Out("\n");
	}
	pPpt->Out("\n");
	for(s=0; s<nsections; ++s)
	{
		pPpt->Out("\tSection "); pPpt->Out(s+1); pPpt->Out(" of "); pPpt->Out(nsections); pPpt->Out(":\n");
		pPpt->Out("\tInitial pressure, pri["); pPpt->Out(s); pPpt->Out("]\t\t\t=\t"); pPpt->Out(pri[s]); pPpt->Out(" bar\n");
		pPpt->Out("\tInitial temperature, Tri["); pPpt->Out(s); pPpt->Out("]\t\t\t=\t"); pPpt->Out(Tri[s]); pPpt->Out(" K\n");
		pPpt->Out("\tInitial velocity, vri["); pPpt->Out(s); pPpt->Out("]\t\t\t=\t"); pPpt->Out(vri[s]); pPpt->Out(" m.s^-1\n");
		pPpt->Out("\n");
	}
	pPpt->Out(Underline("Measurements", "-", "\t"));
	if(USE_DEF_FREQ)
	{
		if(freq==1) 
		{
			pPpt->Out("\tUsing default sampling rate (pPpt->freq)\t=\tonce per timestep\n");
		}
		else
		{
			pPpt->Out("\tUsing default sampling rate (pPpt->freq)\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");
		}
	}
	else
	{
		if(freq==1) pPpt->Out("\tUsing local sampling rate (freq)\t\t=\tonce per timestep\n");
		pPpt->Out("\tUsing local sampling rate (freq)\t\t=\tonce per "); pPpt->Out(freq); pPpt->Out(" timesteps\n");
	}
	pPpt->Out("\tRecording data onwards from "); pPpt->Out(print_from_time); pPpt->Out("s\n");
	pPpt->Out("\n");
	
	if(ntappings>0)
	{
		pPpt->Out(Underline("Measuring locations", "-", "\t"));
		pPpt->Out("\tNo. measuring locations (ntappings)\t\t=\t"); pPpt->Out(ntappings); pPpt->Out("\n");
		pPpt->Out("\tAs a fraction of pipe physical length, measured from the odd end:\n");
	
		for(int m=0; m<ntappings; ++m)
		{
			pPpt->Out("\t- location "); pPpt->Out(m+1); pPpt->Out(" of "); pPpt->Out(ntappings); pPpt->Out(" @ "); pPpt->Out(loc_measure[m]); pPpt->Out("\n");
		}
		pPpt->Out("\n");
	}

	if(num_props_measured>0)
	{
		pPpt->Out("\n\tMax. data points per location (max_pts)\t\t=\t"); pPpt->Out(max_pts); pPpt->Out("\n");
		pPpt->Out("\tRecording:\n");
		if(DIAMETER) pPpt->Out("\t- diameter\n");
		if(AREA) pPpt->Out("\t- area\n");
		if(STATIC_PRESSURE) pPpt->Out("\t- static pressure\n");
		if(TEMPERATURE) pPpt->Out("\t- temperature\n");
		if(DENSITY) pPpt->Out("\t- density\n");
		if(VELOCITY) pPpt->Out("\t- velocity\n");
		if(MACH_NUMBER) pPpt->Out("\t- Mach number\n");
		if(MASS_FLOW_RATE) pPpt->Out("\t- mass flow rate\n");
		if(REYNOLDS_NO) pPpt->Out("\t- Reynolds number\n");
		pPpt->Out("\n");
	}
	pPpt->Out(Underline("Screen output", "-", "\t"));
	if(SHOW_DATA)
	{
		pPpt->Out("\tPrinting pipe data to screen (SHOW_DATA)\t=\t"); pPpt->Out(TrueOrFalse(SHOW_DATA)); pPpt->Out("\n");
	}
	else
	{
		pPpt->Out("\tNot printing pipe data to screen (SHOW_DATA)\t=\t"); pPpt->Out(TrueOrFalse(SHOW_DATA)); pPpt->Out("\n");
	}
	pPpt->Out("\n");
	pPpt->Out("\n");
}

/*
void CPipe::PrintToFile(CProperties* pPpt, int timestep, double time, double ca)
{
	int f, n;
	if(timestep==0)
	{
		for(f=0; f<num_props_measured; ++f)
		{
			fprintf(FILE_OVERALL[f],"%s\t", "TIMEe[R+1]");
			if(!pPpt->CONTINUOUS) fprintf(FILE_OVERALL[f],"%s\t", "Eng[0].ca_elapsed");
			else fprintf(FILE_OVERALL[f],"%f\t", 0);
			for(n=0; n<N-1; ++n)
			{
				fprintf(FILE_OVERALL[f],"%f\t", Node[n].x);
			}
			fprintf(FILE_OVERALL[f],"%f\n", Node[N-1].x);
		}
	}

	for(f=0; f<num_props_measured; ++f)
	{
		fprintf(FILE_OVERALL[f],"%.6f\t", time);
		if(!pPpt->CONTINUOUS) fprintf(FILE_OVERALL[f],"%.6f\t", ca);
		else fprintf(FILE_OVERALL[f],"%f\t", 0);
	}
	f=0;
	if(STATIC_PRESSURE)
	{
		for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL[f],"%f\t", Node[n].p_dash*pPpt->PREF);
		fprintf(FILE_OVERALL[f],"%f\n", Node[N-1].p_dash*pPpt->PREF);
		++f;
	}
	if(TEMPERATURE)
	{
		for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL[f],"%f\t", Node[n].T);
		fprintf(FILE_OVERALL[f],"%f\n", Node[N-1].T);
		++f;
	}
	if(VELOCITY)
	{
		for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL[f],"%f\t", Node[n].U*AREF);
		fprintf(FILE_OVERALL[f],"%f\n", Node[N-1].U*AREF);
		++f;
	}
	if(MACH_NUMBER)
	{
		for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL[f],"%f\t", Node[n].M);
		fprintf(FILE_OVERALL[f],"%f\n", Node[N-1].M);
		++f;
	}
	if(MASS_FLOW_RATE)
	{
		for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL[f],"%f\t", Node[n].mdot);
		fprintf(FILE_OVERALL[f],"%f\n", Node[N-1].mdot);
		++f;
	}
	if(REYNOLDS_NO)
	{
		for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL[f],"%f\t", Node[n].Re);
		fprintf(FILE_OVERALL[f],"%f\n", Node[N-1].Re);
		++f;
	}
}
*/

void CPipe::PrintToFileMovie(CProperties* pPpt, int timestep, double time, double ca)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintToFileMovie\n");}
	int f, n;
	if(num_props_measured>0)
	{
		if(timestep==0)
		{
			for(f=0; f<num_props_measured; ++f)
			{
				fprintf(FILE_OVERALL_MOV[f],"%f\t", 0); // To act as a space above the time column
				fprintf(FILE_OVERALL_MOV[f],"%f\t", 0); // To act as a space above the ca column
				//fprintf(FILE_OVERALL_MOV[f],"%s\t", "Time (s)"); // To act as a space above the time column
				//fprintf(FILE_OVERALL_MOV[f],"%c%s\t", Deg(), "CA"); // To act as a space above the ca column
				for(n=0; n<N-1; ++n)
				{
					fprintf(FILE_OVERALL_MOV[f],"%f\t", Node[n].x);
				}
				fprintf(FILE_OVERALL_MOV[f],"%f\n", Node[N-1].x);
			}
		}

		if(time>=print_from_time)// Only start recording data as specified
		{
			if(timestep%freq==0) // Print data at the specified sampling frequency
			{
				for(f=0; f<num_props_measured; ++f)
				{
					fprintf(FILE_OVERALL_MOV[f],"%.6f\t", time);
					if(!pPpt->CONTINUOUS) fprintf(FILE_OVERALL_MOV[f],"%.6f\t", ca);
					else fprintf(FILE_OVERALL_MOV[f],"%f\t", 0);
				}
				f=0;
				if(DIAMETER) {
					for (n = 0; n < N - 1; ++n) fprintf(FILE_OVERALL_MOV[f], "%f\t", Node[n].d);
					fprintf(FILE_OVERALL_MOV[f], "%f\n", Node[N - 1].d);
					++f;
				}
				if(AREA) {
					for (n = 0; n < N - 1; ++n) fprintf(FILE_OVERALL_MOV[f], "%f\t", Node[n].f);
					fprintf(FILE_OVERALL_MOV[f], "%f\n", Node[N - 1].f);
					++f;
				}
				if(STATIC_PRESSURE) {
					for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL_MOV[f],"%f\t", Node[n].p_dash*pPpt->PREF);
					fprintf(FILE_OVERALL_MOV[f],"%f\n", Node[N-1].p_dash*pPpt->PREF);
					++f;
				}
				if(TEMPERATURE) {
					for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL_MOV[f],"%f\t", Node[n].T);
					fprintf(FILE_OVERALL_MOV[f],"%f\n", Node[N-1].T);
					++f;
				}
				if(DENSITY) {
					for (n = 0; n < N - 1; ++n) fprintf(FILE_OVERALL_MOV[f], "%f\t", Node[n].rho);
					fprintf(FILE_OVERALL_MOV[f], "%f\n", Node[N - 1].rho);
					++f;
				}
				if(VELOCITY) {
					for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL_MOV[f],"%f\t", Node[n].U*AREF);
					fprintf(FILE_OVERALL_MOV[f],"%f\n", Node[N-1].U*AREF);
					++f;
				}
				if(MACH_NUMBER) {
					for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL_MOV[f],"%f\t", Node[n].M);
					fprintf(FILE_OVERALL_MOV[f],"%f\n", Node[N-1].M);
					++f;
				}
				if(MASS_FLOW_RATE) {
					for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL_MOV[f],"%f\t", Node[n].mdot);
					fprintf(FILE_OVERALL_MOV[f],"%f\n", Node[N-1].mdot);
					++f;
				}
				if(REYNOLDS_NO) {
					for(n=0; n<N-1; ++n) fprintf(FILE_OVERALL_MOV[f],"%f\t", Node[n].Re);
					fprintf(FILE_OVERALL_MOV[f],"%f\n", Node[N-1].Re);
					++f;
				}
			}
		}
	}
}

void CPipe::PrintToScreen(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintToScreen\n");}

	int k, S;
	if(BUFFER) pPpt->Out(Underline("Buffer pipe", "=", "", strDesc));
	else
	{
		if(DAMPER) pPpt->Out(Underline("Damper pipe", "=", "", strDesc));
		else pPpt->Out(Underline(Identify(), "=", "", strDesc));
	}
	pPpt->Out("\n");

	if(pPpt->SHOW_x || pPpt->SHOW_X || pPpt->SHOW_d || pPpt->SHOW_dddx || pPpt->SHOW_d2ddx2 || pPpt->SHOW_f || pPpt->SHOW_cfa || pPpt->SHOW_dfdx)
	{
		pPpt->Out(Underline("Geometry", "-"));	
		if(pPpt->SHOW_x)
		{ 
			pPpt->Out("x (mm):\t\t");
			for(S=0; S<N; ++S) { pPpt->Out(Node[S].x*1000); pPpt->Out(" "); } 
			pPpt->Out("\n\n"); 
		}
		if(pPpt->SHOW_X)
		{ 
			pPpt->Out("X (-):\t\t");
			for(S=0; S<N; ++S) { pPpt->Out(Node[S].X); pPpt->Out(" "); }
			pPpt->Out("\n\n");
		}
		if (pPpt->SHOW_d)
		{
			pPpt->Out("d0 (mm):\t");
			for (S = 0; S < N; ++S) { pPpt->Out(Node[S].d0*1000); pPpt->Out(" "); }
			pPpt->Out("\n\n");
		}
		if (pPpt->SHOW_dddx)
		{
			pPpt->Out("dd0dx (mm/mm):\t");
			for (S = 0; S < N; ++S) { pPpt->Out(Node[S].dd0dx); pPpt->Out(" "); }
			pPpt->Out("\n\n");
		}
		if(pPpt->SHOW_d)
		{ 
			pPpt->Out("d (mm):\t\t");
			for(S=0; S<N; ++S) { pPpt->Out(Node[S].d*1000); pPpt->Out(" ");	}
			pPpt->Out("\n\n");
		}
		if(pPpt->SHOW_dddx)
		{ 
			pPpt->Out("dddx (mm/mm):\t"); 
			for(S=0; S<N; ++S) { pPpt->Out(Node[S].dddx); pPpt->Out(" "); }
			pPpt->Out("\n\n");
		}
		if(pPpt->SHOW_d2ddx2)
		{ 
			pPpt->Out("d2ddx2:\t"); 
			for (S = 0; S < N; ++S) { pPpt->Out(Node[S].d2ddx2); pPpt->Out(" "); }
			pPpt->Out("\n\n");
		}
		if (pPpt->SHOW_f)
		{
			pPpt->Out("f (mm^2):\t");
			for (S = 0; S < N; ++S) { pPpt->Out(Node[S].f * 1e6); pPpt->Out(" "); }
			pPpt->Out("\n\n");
		}
		if (pPpt->SHOW_dfdx)
		{
			pPpt->Out("dfdx (mm^2/mm):\t");
			for (S = 0; S < N; ++S) { pPpt->Out(Node[S].dfdx * 1000); pPpt->Out(" "); }
			pPpt->Out("\n\n");
		}
		if(pPpt->SHOW_cfa)
		{ 
			pPpt->Out("cfa (m^2):\t");
			for(S=0; S<N; ++S) { pPpt->Out(Node[S].cfa ); pPpt->Out(" "); }
			pPpt->Out("\n\n");
		}
		if(pPpt->SHOW_cfa_delx)
		{ 
			pPpt->Out("cfa_delx (m):\t");
			for(S=0; S<N; ++S) { pPpt->Out(Node[S].cfa_delx); pPpt->Out(" "); }
			pPpt->Out("\n\n");
		}
	}
	if(pPpt->SHOW_F)
	{
		pPpt->Out(Underline("Flux vector", "-"));
		pPpt->Out("F[0]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].F[0]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("F[1]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].F[1]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("F[2]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].F[2]); pPpt->Out(" ");} pPpt->Out("\n\n");
	}
	if(pPpt->SHOW_C)
	{
		pPpt->Out(Underline("Source vector", "-"));
		pPpt->Out("C[0]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].C[0]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("C[1]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].C[1]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("C[2]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].C[2]); pPpt->Out(" ");} pPpt->Out("\n\n");
	}
	if(pPpt->SHOW_W_pred)
	{
		pPpt->Out(Underline("Predicted solution vector", "-"));
		pPpt->Out("W_pred[0]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].W_pred[0]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("W_pred[1]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].W_pred[1]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("W_pred[2]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].W_pred[2]); pPpt->Out(" ");} pPpt->Out("\n\n");
	}
	if(pPpt->SHOW_F_pred)
	{
		pPpt->Out(Underline("Predicted flux vector", "-"));
		pPpt->Out("F_pred[0]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].F_pred[0]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("F_pred[1]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].F_pred[1]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("F_pred[2]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].F_pred[2]); pPpt->Out(" ");} pPpt->Out("\n\n");
	}
	if(pPpt->SHOW_C_pred)
	{
		pPpt->Out(Underline("Predicted source vector", "-"));
		pPpt->Out("C_pred[0]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].C_pred[0]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("C_pred[1]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].C_pred[1]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("C_pred[2]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].C_pred[2]); pPpt->Out(" ");} pPpt->Out("\n\n");
	}
	if(pPpt->SHOW_W)
	{
		pPpt->Out(Underline("Solution vector", "-"));
		pPpt->Out("W[0]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].W[0]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("W[1]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].W[1]); pPpt->Out(" ");} pPpt->Out("\n\n");
		pPpt->Out("W[2]:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].W[2]); pPpt->Out(" ");} pPpt->Out("\n\n");
	}
	pPpt->Out(Underline("Flow", "-"));	
	pPpt->Out("Direction:\t");
	for(S=0; S<N; ++S)
	{
		if(Node[S].CHOKED) pPpt->Out("c");
		if(fabs(Node[S].U*AREF)==0) pPpt->Out("|");
		else 
		{
			if(fabs(Node[S].U)<1) // If subsonic
			{
				if(Node[S].U*AREF<0) pPpt->Out("<");
				else pPpt->Out(">");
			}
			else // Supersonic
			{
				if(Node[S].U<0) pPpt->Out("<s");
				else pPpt->Out("s>");
			}
		}
		pPpt->Out(" ");
	}
	pPpt->Out("\n");
	pPpt->Out("u (m/s):\t");
	for(S=0; S<N; ++S){pPpt->Out(Node[S].U*AREF); pPpt->Out(" ");} pPpt->Out("\n\n");
	pPpt->Out("m (kg/s):\t");
	for(S=0; S<N; ++S){pPpt->Out(Node[S].mdot); pPpt->Out(" ");} pPpt->Out("\n\n");
	pPpt->Out("Mach ():\t"); 
	for(S=0; S<N; ++S){pPpt->Out(Node[S].M); pPpt->Out(" ");} pPpt->Out("\n\n");
	if (pPpt->HAEMODYNAMICS) {
		pPpt->Out("Undeformed\n");
		pPpt->Out("pressure,\n");
		pPpt->Out("p0 (bar):\t");
		for (S = 0; S < N; ++S) { pPpt->Out(Node[S].p0_dash * pPpt->PREF); pPpt->Out(" "); } pPpt->Out("\n\n");
	}
	pPpt->Out("p (bar):\t"); 
	for(S=0; S<N; ++S){pPpt->Out(Node[S].p_dash*pPpt->PREF); pPpt->Out(" ");} pPpt->Out("\n\n");
	pPpt->Out("Temp. (K):\t");	
	for(S=0; S<N; ++S){pPpt->Out(Node[S].T); pPpt->Out(" ");} pPpt->Out("\n\n");
	pPpt->Out("rho (kg.m^-3):  ");	
	for(S=0; S<N; ++S){pPpt->Out(Node[S].rho); pPpt->Out(" ");} pPpt->Out("\n\n");
//	pPpt->Out("Re:\t"); 
//	for(S=0; S<N; ++S){pPpt->Out(Node[S].Re); pPpt->Out(" ");} pPpt->Out("\n\n");
	pPpt->Out("lambda_I:\t"); 
	for(S=0; S<N; ++S){pPpt->Out(Node[S].CL1[1]); pPpt->Out(" ");} pPpt->Out("\n\n");
	pPpt->Out("lambda_II:\t"); 
	for(S=0; S<N; ++S){pPpt->Out(Node[S].CL2[1]); pPpt->Out(" ");} pPpt->Out("\n\n");
	if(!pPpt->HOMENTROPIC)
	{pPpt->Out("AA at node:\t"); for(S=0; S<N; ++S){pPpt->Out(Node[S].AA[1]); pPpt->Out(" ");} pPpt->Out("\n\n");}
	
	if(METHOD==FandE)
	{
		//cout << Underline(StringToChar(temp), "-");
		pPpt->Out(Underline("Filling and emptying", "-"));
		pPpt->Out("length (m):\t"); 
		for(S=0; S<N-1; ++S) { pPpt->Out(FV[S].l); pPpt->Out(" "); } pPpt->Out("\n\n");
		pPpt->Out("vol. (m^3):\t"); 
		for(S=0; S<N-1; ++S) { pPpt->Out(FV[S].V); pPpt->Out(" "); } pPpt->Out("\n\n");
		pPpt->Out("mass (kg):\t"); 
		for(S=0; S<N-1; ++S) { pPpt->Out(FV[S].m); pPpt->Out(" "); } pPpt->Out("\n\n");
//		cout << "p (bar):\t"; for(S=0; S<N-1; ++S) cout << FV[S].p << " "; cout << endl; cout << endl;
//		cout << "Temp. (K):\t"; for(S=0; S<N-1; ++S) cout << FV[S].T << " "; cout << endl; cout << endl;
//		cout << "rho (kg/m^3):\t"; for(S=0; S<N-1; ++S) cout << FV[S].rho << " "; cout << endl; cout << endl;

		pPpt->Out("h0 (kJ/kg):\t"); 
		for(S=0; S<N-1; ++S) { pPpt->Out(FV[S].h0/1000); pPpt->Out(" "); } pPpt->Out("\n\n");
		pPpt->Out("p0 (bar):\t"); 
		for(S=0; S<N-1; ++S) { pPpt->Out(FV[S].p0/1e5); pPpt->Out(" "); } pPpt->Out("\n\n");
		pPpt->Out("T0. (K):\t"); 
		for(S=0; S<N-1; ++S) { pPpt->Out(FV[S].T0); pPpt->Out(" "); } pPpt->Out("\n\n");
//		cout << "rho0 (kg/m^3):\t"; for(S=0; S<N-1; ++S) cout << FV[S].p0/(pPpt->R_air*FV[S].T0) << " "; cout << endl; cout << endl;
	}

	if(pPpt->SHOW_pathlines)
	{
		pPpt->Out(Underline("Pathlines", "-"));	
		pPpt->Out("Pathline no.\t");
		for(k=0; k<num_pathlines; ++k) {
			if (PathLine[k].new_pathline) { pPpt->Out(PathLine[k].ID); pPpt->Out("n "); }
			else { pPpt->Out(PathLine[k].ID); pPpt->Out(" "); }
		}
		pPpt->Out("\n\n");
		pPpt->Out("X:\t\t");	
		for (k = 0; k < N; ++k) { pPpt->Out(Node[k].X); pPpt->Out(" "); }
		pPpt->Out("\n\n");
		pPpt->Out("XK:\t\t"); 
		for (k = 0; k < num_pathlines; ++k) { pPpt->Out(PathLine[k].XK); pPpt->Out(" "); }
		pPpt->Out("\n\n");
		pPpt->Out("AAK:\t\t"); 
		for (k = 0; k < num_pathlines; ++k) { pPpt->Out(PathLine[k].AA); pPpt->Out(" "); }
		pPpt->Out("\n\n");
	}
}


void CPipe::Update(CProperties* pPpt, double DELZ)
// ====================================================================================================
// Boundary methods have just been run, providing new time level values for CL1, CL2, AA at [R+1]
// Derive physical variables from characteristics
// ====================================================================================================
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Update\n");}

	double k;
	int S;
	for(S=0; S<N; ++S) {
		k = pPpt->gammaAir(Node[S].T);
		Node[S].A = (Node[S].CL1[R+1] + Node[S].CL2[R+1])/2;
		Node[S].T = pow(Node[S].A,2)*(EX ? pPpt->TREFe : pPpt->TREFi);
		Node[S].U = (Node[S].CL1[R+1] - Node[S].CL2[R+1])/(k-1);
		Node[S].M = fabs(Node[S].U/Node[S].A);
/*
		if(pPpt->HOMENTROPIC) {
			Node[S].p_dash = pow((Node[S].CL1[R+1] + Node[S].CL2[R+1])/2, pPpt->QI);
			Node[S].AA[R+1] = Node[S].A/pow(Node[S].p_dash, (k-1)/(2*k));
		}
		else Node[S].p_dash = pow(Node[S].A/Node[S].AA[R+1], (2*k)/(k-1));		
		
		if (pPpt->HAEMODYNAMICS) {
			// Fluid incompressible so density does not update
			// But walls are elastic so update new area and diameter
			Node[S].f = Node[S].W[0] / Node[S].rho;
			Node[S].d = sqrt(Node[S].f * 4 / PI);
			// dfdx and dddx are updated further below since all node f and d values are required to be set first

			// Update mass flow rate and Re number
			Node[S].mdot = Node[S].rho * (Node[S].U * AREF) * Node[S].f;
			Node[S].Re = (Node[S].rho * (fabs(Node[S].U) * AREF) * Node[S].d) / pPpt->ViscosityBlood;
		}
		else { 
			// Update density, mass flow rate, and Re number
			Node[S].rho = ((Node[S].p_dash*pPpt->PREF)*1e5)/(pPpt->R_air*Node[S].T);
			Node[S].mdot = Node[S].rho * (Node[S].U * AREF) * Node[S].f;
			Node[S].Re = (Node[S].rho * (fabs(Node[S].U) * AREF) * Node[S].d) / pPpt->ViscosityAir(Node[S].T);
		}
//*/
///*
		if (pPpt->HAEMODYNAMICS) {
			// Fluid incompressible so density does not update
			// But walls are elastic so update new area and diameter
			Node[S].f = Node[S].W[0] / Node[S].rho;
			Node[S].d = sqrt(Node[S].f * 4 / PI);
			// dfdx and dddx are updated further below since all node f and d values are required to be set first

			// Update mass flow rate and Re number
			Node[S].mdot = Node[S].rho * (Node[S].U * AREF) * Node[S].f;
			Node[S].Re = (Node[S].rho * (fabs(Node[S].U) * AREF) * Node[S].d) / pPpt->ViscosityBlood;

			// Calculate pressure from haemodynamic state equation based on updated area
			Node[S].p_dash = Node[S].p0_dash + Node[S].b * (1 - sqrt(Node[S].f0 / Node[S].f)) / (pPpt->PREF * 1e5);
		}
		else {
			// Calcuate pressure
			if (pPpt->HOMENTROPIC) {
				Node[S].p_dash = pow((Node[S].CL1[R + 1] + Node[S].CL2[R + 1]) / 2, pPpt->QI);
				Node[S].AA[R + 1] = Node[S].A / pow(Node[S].p_dash, (k - 1) / (2 * k));
			}
			else Node[S].p_dash = pow(Node[S].A / Node[S].AA[R + 1], (2 * k) / (k - 1));

			// Update density, mass flow rate, and Re number
			Node[S].rho = ((Node[S].p_dash * pPpt->PREF) * 1e5) / (pPpt->R_air * Node[S].T);
			Node[S].mdot = Node[S].rho * (Node[S].U * AREF) * Node[S].f;
			Node[S].Re = (Node[S].rho * (fabs(Node[S].U) * AREF) * Node[S].d) / pPpt->ViscosityAir(Node[S].T);
		}
//*/
		if(METHOD==W_ALPHA_BETA && (S==0 || S==N-1))// && pPpt->COMBINED_WAB_MOC)
		{
			// Construct solution vector for boundary nodes
			// --------------------------------------------
			Node[S].W[0] = Node[S].rho*Node[S].f;
			Node[S].W[1] = Node[S].rho*Node[S].f*(Node[S].U*AREF);
			// Need haemodynamic version for:
			Node[S].W[2] = Node[S].rho * (((Node[S].p_dash * pPpt->PREF) * 1e5) / (Node[S].rho * (k - 1)) + 0.5 * pow((Node[S].U * AREF), 2)) * Node[S].f;
		}
		
		// Record previous Riemann values and entropy levels
		// -------------------------------------------------
		Node[S].CL1[R] = Node[S].CL1[R+1];
		Node[S].CL2[R] = Node[S].CL2[R+1];
		Node[S].AA[R] = Node[S].AA[R+1];
	}	
	
	if (pPpt->HAEMODYNAMICS) {
///*
		// Calculate dddx and dfdx for internal nodes using a central differencing scheme
		for (S = 1; S < N - 1; ++S) {
			Node[S].dddx = (Node[S + 1].d - Node[S - 1].d) / (Node[S + 1].x - Node[S - 1].x);
			Node[S].dfdx = (Node[S + 1].f - Node[S - 1].f) / (Node[S + 1].x - Node[S - 1].x);
		}

		// Must set end nodes separately
		S = 0;
		Node[S].dddx = (Node[S + 1].d - Node[S].d) / (Node[S + 1].x - Node[S].x);
		Node[S].dfdx = (Node[S + 1].f - Node[S].f) / (Node[S + 1].x - Node[S].x);

		S = N - 1;
		Node[S].dddx = (Node[S].d - Node[S - 1].d) / (Node[S].x - Node[S - 1].x);
		Node[S].dfdx = (Node[S].f - Node[S - 1].f) / (Node[S].x - Node[S - 1].x);
//*/
	}

	// If F&E, calculate finite volume properties for the end of this step, now knowing the mass flow rates at the throats
	if(METHOD==FandE) {
///*
		// Calculate time increment
		// ========================
		double del_t = (DELZ/AREF)*pPpt->xref;
		
		for(int S=0; S<N-1; ++S)	// For each pipe finite volume
		{
			// Conservation of mass
			// ====================
			double mdot_odd = Node[S].mdot;
			double mdot_even = Node[S+1].mdot;
			double mdot = mdot_odd - mdot_even;	

			// Conservation of energy
			// ======================
			double R_gas = pPpt->R_air;
			
			// Initialise volume variables
			double m_temp = FV[S].m;
			double V_temp = FV[S].V;
			double p_temp = FV[S].p0;
			double T_temp = FV[S].T0;
			double rho_temp = p_temp/(R_gas*T_temp);
			double u_temp = pPpt->cvAir(T_temp)*T_temp;
			double h_temp = pPpt->cpAir(T_temp)*T_temp;
			
			// Enthalpy fluxes - Perfect gas relationship
			double h_odd, h_even;
			double u_odd, u_even;

			// Odd end enthalpy
			if(Node[S].U > 0)
			{
				//h_odd = pPpt->cpAir()*Node[S].T;	// Inflow - use node stagnation/static temperature
				h_odd = pPpt->cpAir(Node[S].T)*Node[S].T + pow(Node[S].U*AREF,2)/2; 				
				
				u_odd = pPpt->cvAir(Node[S].T)*Node[S].T; 
				//u_odd = pPpt->cvAir()*Node[S].T + pow(Node[S].U*AREF,2)/2; 
			}
			else
			{
				//h_odd = pPpt->cpAir()*T_temp;	// Outflow - use volume stagnation temperature
				h_odd = pPpt->cpAir(Node[S].T)*Node[S].T;
				//h_odd = pPpt->cpAir()*Node[S].T + pow(Node[S].U*AREF,2)/2; 

				u_odd = pPpt->cvAir(T_temp)*T_temp;
			}

			// Even end enthalpy
			if(Node[S+1].U < 0)
			{
				// Inflow - use node stagnation/static temperature
				
				h_even = pPpt->cpAir(Node[S+1].T)*Node[S+1].T;
				//h_even = pPpt->cpAir()*Node[S+1].T + pow(Node[S+1].U*AREF,2)/2; 
				
				u_even = pPpt->cvAir(Node[S+1].T)*Node[S+1].T; 
				//u_even = pPpt->cvAir()*Node[S+1].T + pow(Node[S+1].U*AREF,2)/2; 
			}
			else
			{
				// Outflow - use volume stagnation temperature
				//h_even = pPpt->cpAir()*T_temp;
				h_even = pPpt->cpAir(Node[S+1].T)*Node[S+1].T;
				//h_even = pPpt->cpAir()*Node[S+1].T + pow(Node[S+1].U*AREF,2)/2; 

				u_even = pPpt->cvAir(T_temp)*T_temp;
			}

			//double ke_odd = 1.5*(pow(Node[S].U*AREF,2)/2)*mdot_odd;

			// Enthalpy flux
			double sum_mdot_h = mdot_odd*h_odd - mdot_even*h_even;
		
			// Rate of change of internal energy
			double mdotu = mdot*u_temp;
			//double mdotu = mdot_odd*u_odd - mdot_even*u_even;
			
			// Partial derivatives with p				
			double pdrhopdp = 1/(R_gas*T_temp);
			double pdupdp = 0;
			double pdhpdp = 0;
			double pdRpdp = 0;		// Assume constant R

			// Partial derivatives with T
			double pdrhopdT = -1*(p_temp/(R_gas*pow(T_temp,2)));
			double pdupdT = pPpt->cvAir(T_temp);
			double pdhpdT = pPpt->cpAir(T_temp);
			double pdRpdT = 0;		// Assume constant R

			// Partial derivatives with rho
			double pdRpdrho = 0;	// Assume constant R
			
			// Partial derivatives with phi (all zero)
			double pdrhopdphi = 0;
			double pdupdphi = 0;
			double pdhpdphi = 0;
			double pdRpdphi = 0;	// Assume constant R
			
			// Further simplifications
			double Vdot = 0;	// Constant volume (since a manifold, not a cylinder)
			double phidot = 0;	// Ignore mass fraction effects
			double Qwdot = 0;	// No wall heat transfer			

			// Terms B, C & D; use instead of A', B' & C'
			//double B = -1*R_gas*T_temp*(Vdot/V_temp) + (1/m_temp)*(Qwdot + sum_mdot_h - mdotu - ke_odd);
			double B = -1*R_gas*T_temp*(Vdot/V_temp) + (1/m_temp)*(Qwdot + sum_mdot_h - mdotu);
			double C = 1 + (T_temp/R_gas)*pdRpdT;
			double D = 1 - (p_temp/R_gas)*pdRpdp;

			// Terms A', B' & C'; use instead of B, C & D
			double A_dash = pdhpdT + (pdrhopdT/pdrhopdp)*((1/rho_temp) - pdhpdp);
			double B_dash = (1 - rho_temp*pdhpdp)/pdrhopdp;
			double C_dash = pdhpdphi + (pdrhopdphi/pdrhopdp)*((1/rho_temp) - pdhpdp);
			
			// Rate of change of temperature with time, Tdot
			double Tdot = (B - (p_temp/D)*pdupdp*( mdot/m_temp - Vdot/V_temp + pdRpdphi*(phidot/R_gas) ) - pdupdphi*phidot)/(pdupdT + (C/D)*(p_temp/T_temp)*pdupdp);
//			double Tdot = (B_dash/A_dash)*( (mdot/m_temp)*(1 - h_temp/B_dash) - Vdot/V_temp - (C_dash/B_dash)*phidot + (1/(B_dash*m_temp))*(sum_mdot_h - Qwdot) );
			
			//Tdot *= 0.99;
			// Rate of change of pressure with time, pdot
			double pdot = (rho_temp/pdrhopdp)*(-1*Vdot/V_temp - (1/rho_temp)*pdrhopdT*Tdot - (1/rho_temp)*pdrhopdphi*phidot + mdot/m_temp);
			
			// Calculate the new volume properties at the end of the current step
			FV[S].m += mdot*del_t;
			FV[S].T0 += Tdot*del_t;
			FV[S].p0_old = FV[S].p0; // Save previous value (for evaluating STEADY conditions)
			FV[S].p0 += pdot*del_t;
			
			//double rho_temp2 = (pdot*del_t)/(pPpt->R_air*(Tdot*del_t));
			//FV[S].p0 += 0.5*rho_temp*pow(Node[S].U*AREF,2);// + 0.5*Node[S+1].rho*pow(Node[S+1].U*AREF,2);
		}
//*/
	}
}


void CPipe::PrintBoundaries(CProperties* pPpt)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".PrintBoundaries\n");}
	int temp_precision = cout.precision(); cout << setprecision(6/*10*/);
	if(EX) pPpt->Out(Underline("Exhaust Pipe", "=", ID, "\t"));
	else pPpt->Out(Underline("Intake Pipe", "=", ID, "\t"));
	pPpt->Out("\n");
	pPpt->Out("\tODD END BOUNDARY\t\t\t\t=\t");
	pPpt->Out(GetBoundaryName(Node[0].bc));
	pPpt->Out("\n\tend correction parameter (Lc/d)\t\t\t=\t"); pPpt->Out(end_corr_odd_p); pPpt->Out("\n");
	pPpt->Out("\n");
	pPpt->Out("\tEVEN END BOUNDARY\t\t\t\t=\t");
	pPpt->Out(GetBoundaryName(Node[N-1].bc));
	pPpt->Out("\n\tend correction parameter (Lc/d)\t\t\t=\t"); pPpt->Out(end_corr_even_p); pPpt->Out("\n");
	pPpt->Out("\n");
/*
	int S;
	if(pPpt->SHOW_x)
	{
		cout << Underline("Node x values (m)", "-", "\t");
		for(S=0; S<N; ++S) cout << "\tNode[" << S << "].x = " << Node[S].x << endl;
		cout << endl;
	}
	if(pPpt->SHOW_X)
	{
		cout << Underline("Node X values ()", "-", "\t");
		for(S=0; S<N; ++S) cout << "\tNode[" << S << "].X = " << Node[S].X << endl;
		cout << endl;
	}
	if(pPpt->SHOW_DELX)
	{
		cout << Underline("DELX_L<-[S]->DELX_R (*xref (m))", "-", "\t");
		for(S=0; S<N; ++S) cout << "\t" << Node[S].DELX_L*pPpt->xref << " <-[" << S << "]-> " << Node[S].DELX_R*pPpt->xref << endl;
		cout << endl;
	}
	if(pPpt->SHOW_d)
	{
		cout << Underline("Node diameter (d) values (m)", "-", "\t");
		for(S=0; S<N; ++S) cout << "\tNode[" << S << "].d = " << Node[S].d << endl;
		cout << endl;
	}
	if(pPpt->SHOW_dddx)
	{
		cout << Underline("Node dddx values ()", "-", "\t");
		for(S=0; S<N; ++S) cout << "\tNode[" << S << "].dddx = " << Node[S].dddx << endl;
		cout << endl;
	}
	if(pPpt->SHOW_d2ddx2)
	{
		cout << Underline("Node d2ddx2 values ()", "-", "\t");
		for(S=0; S<N; ++S) cout << "\tNode[" << S << "].d2ddx2 = " << Node[S].d2ddx2 << endl;
		cout << endl;
	}
	if(pPpt->SHOW_f)
	{
		cout << Underline("Node area (f) values (m^2)", "-", "\t");
		for(S=0; S<N; ++S) cout << "\tNode[" << S << "].f = " << Node[S].f << endl;
		cout << endl;
	}
	if(pPpt->SHOW_dfdx)
	{
		cout << Underline("Node dfdx (dfdx) values (m)", "-", "\t");
		for(S=0; S<N; ++S) cout << "\tNode[" << S << "].dfdx = " << Node[S].dfdx << endl;
		cout << endl;
	}
	if(PARALLEL)
//	if(pPpt->SHOW_gap)
	{
		cout << Underline("Node gap values (m)", "-", "\t");
		for(S=0; S<N; ++S) cout << "\tNode[" << S << "].gap = " << Node[S].gap << endl;
		cout << endl;
	}
	if(PARALLEL)
//	if(pPpt->SHOW_vtan)
	{
		cout << Underline("Node vtan values (m/s)", "-", "\t");
		for(S=0; S<N; ++S) cout << "\tNode[" << S << "].vtan = " << Node[S].vtan << endl;
		cout << endl;
	}
//*/
	cout << setprecision(temp_precision);
}

void CPipe::SetupFiles(CProperties* pPpt, string parent_assy_res_dir)
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".SetupFiles\n");}

	// Construct the correct assembly results directory
	// ================================================
	string res_str, fande_str;
/*
	std::string dir_str, res_str, fande_str;
	dir_str = "res_assembly";

	if(AssyID>=10) dir_str += int(AssyID/10) + 48;
	dir_str += (AssyID - int(AssyID/10)*10) + 48;
	dir_str += "\\";
	//RES_DIR = StringToChar(pPpt->case_dir + dir_str); // Set the results directory for this object
	RES_DIR = StringToChar(pPpt->case_res_dir + dir_str); // Set the results directory for this object
	//cout << dir_str << std::endl;
	//cout << RES_DIR << endl;
*/
	RES_DIR = StringToChar(parent_assy_res_dir);
	//cout << RES_DIR << endl;
	
	// Construct names and open files for the individual pipe tappings results
	// =======================================================================
	if(ntappings>0) {
		if(EX) res_str = "res_ex"; else res_str = "res_in";
		res_str += "_p";						// 'p' for pipe
		if(ID>=10) res_str += int(ID/10) + 48;
		res_str += (ID - int(ID/10)*10) + 48;	// Add the pipe ID
		fande_str = res_str;
		res_str += "_loc";						// Add the remainder of the file name
		res_str += pPpt->strFileExt;			// Add the file extension
		fande_str += "_fande";					// Add the remainder of the file name
		fande_str += pPpt->strFileExt;			// Add the file extension
		FILE_LOC = fopen(ConstructString(pPpt, RES_DIR, res_str), "w");
		if(METHOD==FandE) FILE_FandE = fopen(ConstructString(pPpt, RES_DIR, fande_str), "w");
	}

	// Construct names and open files for results across the whole pipe
	// ================================================================
	if(num_props_measured>0) {
		FILE_OVERALL_MOV = new FILE* [num_props_measured];
		if(EX) res_str = "res_ex"; else res_str = "res_in";
		res_str += "_p";
		if(ID>=10) res_str += int(ID/10) + 48;
		res_str += (ID - int(ID/10)*10) + 48; 
		std::string base_str = res_str;
		std::string mov_str;
		
		int n=0;
		if(DIAMETER) {
			res_str = base_str + "_dia";
			mov_str = res_str + pPpt->strFileExt;
			res_str += pPpt->strFileExt;
			FILE_OVERALL_MOV[n] = fopen(ConstructString(pPpt, RES_DIR, mov_str), "w");
			++n;
		}
		if(AREA) {
			res_str = base_str + "_area";
			mov_str = res_str + pPpt->strFileExt;
			res_str += pPpt->strFileExt;
			FILE_OVERALL_MOV[n] = fopen(ConstructString(pPpt, RES_DIR, mov_str), "w");
			++n;
		}
		if(STATIC_PRESSURE) {
			res_str = base_str + "_press";
			mov_str = res_str + pPpt->strFileExt;
			res_str += pPpt->strFileExt;
			FILE_OVERALL_MOV[n] = fopen(ConstructString(pPpt, RES_DIR, mov_str), "w");
			++n;
		}
		if(TEMPERATURE) {
			res_str = base_str + "_temp";
			mov_str = res_str + pPpt->strFileExt;
			res_str += pPpt->strFileExt;
			FILE_OVERALL_MOV[n] = fopen(ConstructString(pPpt, RES_DIR, mov_str), "w");
			++n;
		}
		if(DENSITY){
			res_str = base_str + "_dens";
			mov_str = res_str + pPpt->strFileExt;
			res_str += pPpt->strFileExt;
			FILE_OVERALL_MOV[n] = fopen(ConstructString(pPpt, RES_DIR, mov_str), "w");
			++n;
		}
		if(VELOCITY) {
			res_str = base_str + "_vel";
			mov_str = res_str + pPpt->strFileExt;
			res_str += pPpt->strFileExt;
			FILE_OVERALL_MOV[n] = fopen(ConstructString(pPpt, RES_DIR, mov_str), "w");
			++n;
		}
		if(MACH_NUMBER) {
			res_str = base_str + "_mach";
			mov_str = res_str + pPpt->strFileExt;
			res_str += pPpt->strFileExt;
			FILE_OVERALL_MOV[n] = fopen(ConstructString(pPpt, RES_DIR, mov_str), "w");
			++n;
		}
		if(MASS_FLOW_RATE) {
			res_str = base_str + "_mdot";
			mov_str = res_str + pPpt->strFileExt;
			res_str += pPpt->strFileExt;
			FILE_OVERALL_MOV[n] = fopen(ConstructString(pPpt, RES_DIR, mov_str), "w");
			++n;
		}
		if(REYNOLDS_NO) {
			res_str = base_str + "_reyn";
			mov_str = res_str + pPpt->strFileExt;
			res_str += pPpt->strFileExt;
			FILE_OVERALL_MOV[n] = fopen(ConstructString(pPpt, RES_DIR, mov_str), "w");
			++n;
		}
	}
}

/*
double CPipe::FrictionFactor(CProperties* pPpt, double d, double Re, char* callingfunc, double cftrans)
// ====================================================	//
// Swamee and Jain relationship for friction factor		//
// For Reynold's No. 5x10^3 <= Re <= 10^8				//
// For relative roughness 10^-6 <= k/D <= 10^-2			//
// ====================================================	//
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".FrictionFactor\n");}
//cout << callingfunc << ": d = " << d << endl;
	double ff;

	if(cftrans==0) ff = 0;
	else
	{
		if(d<1e-9) ff = 0.298932973;
		else
		{
			//if(epsilon/d>1e-2 || epsilon/d<1e-6) cout << "CLoop::PipeWallFrictionFactor (" << callingfunc << "): relative roughness (epsilon/d = " << epsilon << "/" << d << " = " << epsilon/d << ") for friction factor calculation out of normal range\n";

			// Determine regime based on Re no.
			// --------------------------------
			if(Re<3e3 || Re>1e8)
			{
				if(Re<3e3)
				{
//cout << "Re = " << Re << " => too low. Moody diagram" << endl;
					if(Re>=1) ff = -0.079644324*log10(Re) + 0.298932973;	// Derived from Moody diagram
					else ff = 0.298932973;
				}
				else // Re>1e8
				{
					if(Re>1e8)
					{
//cout << "Re = " << Re << " => too high. Set Re = 1e8 then use S&J equation" << endl;
						Re = 1e8; // OK since Moody diagram shows flat line for Re greater than this
						
						if(pPpt->FRICTION_MODEL==pPpt->CONSTANT_FF) ff = pPpt->ff_const;
						else
						{
							if(pPpt->FRICTION_MODEL==pPpt->SWAMEEJAIN_FF) ff = SwameeJain(d, Re);
							else
							{
								if(pPpt->FRICTION_MODEL==pPpt->HAALAND_FF) ff = Haaland(d, Re);
								else
								{
									if(pPpt->FRICTION_MODEL==pPpt->COLEBROOK_FF) ff = Colebrook(d, Re);
									else ff = SwameeJain(d, Re);
								}
							}
						}
					}
				}
				//cout << "Warning: Re (" << Re << ") for friction factor calculation out of normal range\n";
			}
			else
			{
//cout << "Re = " << Re << " => in range for S&J equation." << endl;

				if(pPpt->FRICTION_MODEL==pPpt->CONSTANT_FF) ff = pPpt->ff_const;
				else
				{
					if(pPpt->FRICTION_MODEL==pPpt->SWAMEEJAIN_FF) ff = SwameeJain(d, Re);
					else
					{
						if(pPpt->FRICTION_MODEL==pPpt->HAALAND_FF) ff = Haaland(d, Re);
						else
						{
							if(pPpt->FRICTION_MODEL==pPpt->COLEBROOK_FF) ff = Colebrook(d, Re);
							else ff = SwameeJain(d, Re);
						}
					}
				}
			}
		}

		// Apply friction enhancement factor
		// ---------------------------------		
		ff *= cftrans;
//cout << "cftrans = " << cftrans << endl;
//cout << "ff = " << ff << endl;
	}

	// Return result
	// -------------
	return ff;
}
*/

double CPipe::FrictionFactor(CProperties* pPpt, double d, double Re, char* callingfunc, double cftrans)
// ====================================================	//
// Turbulent pipe flow relations for friction factor	//
// - for Reynold's No. 5x10^3 <= Re <= 10^8				//
// - for relative roughness 10^-6 <= epsilon/d <= 10^-2	//
// ====================================================	//
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".FrictionFactor\n");}
//cout << callingfunc << ": d = " << d << endl;
	double ff;

	if(cftrans==0 || Re==0) ff = 0;
	else
	{
		if(pPpt->FRICTION_MODEL==pPpt->RICARDO_FF)
		{
			if(Re < 3000) ff = 64/Re; // Poiseuille law for laminar flow
			else
			{
				if(Re <= 4000) ff = 0.316/pow(Re,0.25); // Blasius law for smooth pipe
				else ff = Ricardo(Re);
			}
		}
		else
		{
			if(Re < 2000) ff = 64/Re; // Poiseuille law for laminar flow
			else
			{
				if(Re < 5000) ff = 0.316/pow(Re,0.25); // Blasius law for smooth pipe
				else
				{
					if(pPpt->FRICTION_MODEL==pPpt->CONSTANT_FF) ff = pPpt->ff_const;
					else
						{
						if(pPpt->FRICTION_MODEL==pPpt->SWAMEEJAIN_FF) ff = SwameeJain(d, Re);
						else
						{
							if(pPpt->FRICTION_MODEL==pPpt->HAALAND_FF) ff = Haaland(d, Re);
							else
							{
								if(pPpt->FRICTION_MODEL==pPpt->COLEBROOK_FF) ff = Colebrook(d, Re);
								else ff = SwameeJain(d, Re);
							}
						}
					}
				}
			}
		}
	}
	// Apply friction enhancement factor
	ff *= cftrans;
	return ff;
}

double CPipe::Colebrook(double d, double Re)
// Colebrook-White (1939) implicit relationship for friction factor
{
	double f = 1; // Initial guess
	double f_prev;
	int counter = 0;
	do{
		++counter;
		f_prev = f;
		f = pow(1 / (-2.0*log10( (this->epsilon/d)/3.7 + 2.51/(Re*sqrt(f)) )), 2);
/*
cout << "counter = " << counter << endl;
cout << "f_prev = " << f_prev << endl;
cout << "f = " << f << endl;
cout << endl;
//*/
	}
	//while(fabs(f - f_prev) > 1e-6);
	while((fabs(f - f_prev)/f)*100 > 0.01); // In %
	return f;
}

double CPipe::ConvectiveHTCoefficient(CProperties* pPpt, double f, double rho, double u, double T, double Re, double del_x)
// ========================================================================================================	//
// Returns convective heat transfer coefficient hc using various methods									//
// ========================================================================================================	//
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".ConvectiveHTCoefficient\n");}

	if(pPpt->HEAT_TRANSFER==pPpt->REYNOLDS_ANALOGY)	// Reynolds' analogy convective heat transfer coefficient
	{
		//double cp = (pPpt->gammaAir()*pPpt->R_air)/(pPpt->gammaAir()-1); // For an ideal gas
		//return (f/2)*rho*fabs(u)*cp*HGTRANS;
		return (f/2)*rho*fabs(u)*pPpt->cpAir(T)*HGTRANS;
	}
	else
	{
		if(pPpt->HEAT_TRANSFER==pPpt->NUSSELT_RELATION)
		{
			if(Re<1e3) Re=5e4;//1e3;

			double Nu = 0.0483*pow(Re, 0.783); // For straight sections of exhaust pipe (Heywood)
		//	double Nu = 0.022*pow(Re, 0.8);
			double k_air = 0.025;
			double L = del_x;
			return (Nu*k_air/L)*HGTRANS;	// Since Nu = hc*L/k;
		}
		else
		{
			if(pPpt->HEAT_TRANSFER==pPpt->RICARDO_HT)
			{
				//double cp = (pPpt->gammaAir()*pPpt->R_air)/(pPpt->gammaAir()-1); // For an ideal gas
				double k_air = 0.025;
				//double Pr = cp*pPpt->ViscosityAir(T)/k_air;
				double Pr = pPpt->cpAir(T)*pPpt->ViscosityAir(T)/k_air;
				//return (f/4/*==Cf/2*/)*rho*fabs(u)*cp*pow(Pr,-2/3)*HGTRANS;
				return (f/4/*==Cf/2*/)*rho*fabs(u)*pPpt->cpAir(T)*pow(Pr,-2/3)*HGTRANS;
			}
			else
			{
				cout << "unknown heat transfer model\n";
				exit(1);
			}
		}
	}
}

bool CPipe::Steady(CProperties* pPpt)
// ========================================================================================================	//
// Returns true if all pipe nodes show steady conditions													//
// ========================================================================================================	//
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".Steady\n");}
	
	bool STEADY = true; // Assume pipe is steady to start with
	
	if(METHOD==this->FandE) // Filling and emptying - check upon consecutive pressure reading in the volume
	{
		for(int S=0; S<N-1; ++S)
		{
/*
			if(fabs((FV[S].p - FV[S].p_old)/FV[S].p)*100 > pPpt->tol_steady*pPpt->tol_steady_multiplier)
			{
				STEADY = false;
				S = N-1; // End loop
			}
*/
			if(fabs((FV[S].p0 - FV[S].p0_old)/FV[S].p0)*100 > pPpt->tol_steady*pPpt->tol_steady_multiplier)
			{
				STEADY = false;
				S = N-1; // End loop
			}
		}
	}
	else
	{
		double u_temp_new, u_temp_old;
		for(int S=0; S<N; ++S)
		{
			Node[S].STEADY = true; // Assume node is steady to start with

			u_temp_new = ((Node[S].CL1[1] - Node[S].CL2[1])/(pPpt->gammaAir(Node[S].T) - 1))*AREF;
			u_temp_old = ((Node[S].CL1[0] - Node[S].CL2[0])/(pPpt->gammaAir(Node[S].T) - 1))*AREF;

			//if(fabs((Node[S].CL1[1] - Node[S].CL1[0])/Node[S].CL1[1])*100 > pPpt->tol_steady*pPpt->tol_steady_multiplier || fabs((Node[S].CL2[1] - Node[S].CL2[0])/Node[S].CL2[1])*100 > pPpt->tol_steady*pPpt->tol_steady_multiplier)
			if(fabs((u_temp_new - u_temp_old)/u_temp_old)*100 > pPpt->tol_steady*pPpt->tol_steady_multiplier)
			{
				Node[S].STEADY = false;
				STEADY = false;
				S = N-1; // End loop
			}
			/*
			cout << "Node[" << S << "].CL1[1] = " << Node[S].CL1[1] << endl;
			cout << "Node[" << S << "].CL1[0] = " << Node[S].CL1[0] << endl;
			cout << "Node[" << S << "].CL2[1] = " << Node[S].CL2[1] << endl;
			cout << "Node[" << S << "].CL2[0] = " << Node[S].CL2[0] << endl;
			*/
		}
	}
	return STEADY;
}


void CPipe::FillingAndEmptying(CProperties* pPpt, double DELZ)
// Filling and emptying
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".FillingAndEmptying\n");}

	// Calculate time increment
	// ========================
//	double del_t = (DELZ/AREF)*pPpt->xref;
	
	double T0_throat;	// Throat stagnation temperature (K)
	double p0_throat;	// Throat stagnation pressure (Pa)
	double ps_throat;	// Throat static pressure (Pa)

	int direction;

	for(int S=0; S<N; ++S)		// For each node
	{		
		if(S>0 && S<N-1)	// If not a boundary node
		{
			// Establish flow direction based on volume stagnation pressures
			if(FV[S-1].p0 > FV[S].p0)
			{
				direction = 1;
				T0_throat = FV[S-1].T0;	
				p0_throat = FV[S-1].p0;					
				ps_throat = FV[S].p0;
			}
			else
			{
				direction = -1;
				T0_throat = FV[S].T0;						
				p0_throat = FV[S].p0;					
				ps_throat = FV[S-1].p0;
			}
			
			if(ps_throat > p0_throat) // Indicates wrong direction
			{
				if(direction==1)
				{
					direction = -1;
					T0_throat = FV[S].T0;						
					p0_throat = FV[S].p0;					
					ps_throat = FV[S-1].p0;
				}
				else
				{
					direction = 1;
					T0_throat = FV[S-1].T0;	
					p0_throat = FV[S-1].p0;					
					ps_throat = FV[S].p0;
				}
			}
			
			bool REVERSED = false;
			double* AU;
			AU = FillingAndEmptyingEquations(pPpt, T0_throat, p0_throat, ps_throat, REVERSED);
			
			// Derive the characteristics
			Node[S].CL1[R+1] = AU[0] + ((pPpt->gammaAir(T0_throat)-1)/2)*direction*AU[1];
			Node[S].CL2[R+1] = AU[0] - ((pPpt->gammaAir(T0_throat)-1)/2)*direction*AU[1];	
			
			delete [] AU;
		}
		else // Boundary node
		{
/*
			// Provide initial guess for lambda_in
			// ===================================
			int V;
			if(S==0) V=S; else V=S-1;

			// Assume pressure at boundary is stagnation pressure p0, i.e. ps==p0 and U=0, CL1==CL2==A
			// This is true in the main body of the volume but not at the boundaries - hence only initial guess
//			if(S!=0) 
			Node[S].CL1[R+1] = pow((FV[V].p0/1e5)/pPpt->PREF, (pPpt->gammaAir()-1)/(2*pPpt->gammaAir()))*Node[S].AA[R+1];
//			if(S==0) 
			Node[S].CL2[R+1] = pow((FV[V].p0/1e5)/pPpt->PREF, (pPpt->gammaAir()-1)/(2*pPpt->gammaAir()))*Node[S].AA[R+1];	
			
//			Node[S].CL1[R+1] = Node[S].CL1[R];
//			Node[S].CL2[R+1] = Node[S].CL2[R];

*/
/*
			if(S==0)// Odd end boundary node
			{
				T0_throat = FV[S].T0;
				p0_throat = FV[S].p0;
			}
			else	// Even end boundary node
			{
				T0_throat = FV[S-1].T0;
				p0_throat = FV[S-1].p0;
			}

			double TREF = EX ? pPpt->TREFe : pPpt->TREFi;
			double A_throat = sqrt(T0_throat/TREF);
			double U_throat = 0;

			// Derive the characteristics
			Node[S].CL1[R+1] = A_throat + ((pPpt->gammaAir()-1)/2)*U_throat;
			Node[S].CL2[R+1] = A_throat - ((pPpt->gammaAir()-1)/2)*U_throat;
*/
/*
			double A_throat, A_throat_old, U_throat;
			double AA_throat = Node[S].AA[R+1];
			double TREF = EX ? pPpt->TREFe : pPpt->TREFi;

			if(S==0)// Odd end boundary node
			{
				T0_throat = FV[S].T0;
				p0_throat = FV[S].p0;
			}
			else	// Even end boundary node
			{
				T0_throat = FV[S-1].T0;
				p0_throat = FV[S-1].p0;
			}

			// Guess an initial value for A_throat
			A_throat = Node[S].A;

			// Recurrence relation for A

//			cout << "T0_throat = " << T0_throat << endl;
//			cout << "p0_throat = " << p0_throat << endl;
//			cout << "AA_throat = " << AA_throat << endl;
//			cout << "A_throat initial = " << A_throat << endl;
//			cout << endl;

			do
			{
				A_throat_old = A_throat;

				A_throat = AA_throat * pow( p0_throat/
											((1 + (pPpt->cpAir()/pPpt->R_air)*(T0_throat/(pow(A_throat,2)*TREF) - 1)) 
											  * (pPpt->PREF*1e5)) ,(pPpt->gammaAir()-1)/(2*pPpt->gammaAir()));


				A_throat = sqrt( 
					(1/TREF)*
					(T0_throat/
					(1 +
					 ((pPpt->R_air/pPpt->cpAir())*
					  (p0_throat/(pow(A_throat/AA_throat, (2*pPpt->gammaAir())/(pPpt->gammaAir()-1))*pPpt->PREF*1e5) - 1))
					)
					)
					);


//				cout << "A_throat_old = " << A_throat_old << endl;
//				cout << "A_throat = " << A_throat << endl;
//				cout << endl;

			}while((fabs(A_throat - A_throat_old)/A_throat_old)*100 > 0.001); // % 

			// Know A, derive U
			U_throat = sqrt(2*pPpt->cpAir()*(T0_throat - pow(A_throat,2)*TREF))/AREF;

			// Assume flow is toward the boundary

			if(S==0)// Odd end boundary node
			{
				direction = -1;
				Node[S].CL2[R+1] = A_throat - ((pPpt->gammaAir()-1)/2)*direction*U_throat;
			}
			else	// Even end boundary node
			{
				direction = 1;
				Node[S].CL1[R+1] = A_throat + ((pPpt->gammaAir()-1)/2)*direction*U_throat;
			}
//*/
		}
	}
}

double* CPipe::FillingAndEmptyingEquations(CProperties* pPpt, double T0, double p0, double ps, bool& rREVERSED)
// Filling and emptying eqiations
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".FillingAndEmptyingEquations\n");}

	double M, Ts, A, U;

	/*
	u = sqrt(
	((2*pPpt->gammaAir())
	/(pPpt->gammaAir()-1))
	*(ps/rhos)*( pow(p0/ps, (pPpt->gammaAir()-1)/pPpt->gammaAir()) - 1 )
	);
	// Subsonic isentropic compressible flow equation for velocity.
	*/
	
	//if(fabs((p0 - ps)/1e5) < 1e-5/*1e-6*/) M = 0; // To avoid sqrt(-ve);
	if(p0/ps < 1)// || fabs((p0 - ps)/1e5) < 1e-6)
	{
		M = 0; // To avoid sqrt(-ve);
		//					cout << "fabs((p0 - ps)/1e5) = " << fabs((p0 - ps)/1e5) << endl;
		//					cout << endl;
		rREVERSED = true;
	}
	else M = sqrt(
		(2/(pPpt->gammaAir(T0)-1))
		*( pow(p0/ps, (pPpt->gammaAir(T0)-1)/pPpt->gammaAir(T0)) - 1 )
		);
	// Subsonic isentropic compressible flow equation for Mach no. (no static density required)
	//		if(M > 1) M = 1; // Choked condition
	
	// Use stagnation temperature and Mach no. to obtain static temperature
	Ts = T0/(1 + ((pPpt->gammaAir(T0)-1)/2)*pow(M,2));
	
	// Now calculate A
	A = sqrt(pPpt->gammaAir(Ts)*pPpt->R_air*Ts)/AREF;
	
	// Hence U
	U = M*A;

	double* AU;
	AU = new double [2];
	AU[0] = A;
	AU[1] = U;

	return AU;
}

void CPipe::FillingAndEmptyingAdjust3(CProperties* pPpt, double DELZ, int& rCOUNTER, double TIME)
// Filling and emptying
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".FillingAndEmptyingAdjust3\n");}

	double p0 [2];
	double T0 [2];
	p0[0] = 1.25e5;
	T0[0] = 500;
	p0[1] = 1e5;
	T0[1] = 300;

	double A_lambda, U_lambda;
	double ps_lambda, p0_lambda, Ts_lambda, T0_lambda;
	double CL1_old, CL2_old;

	//double ps_throat;	// Throat static pressure (Pa);
	double ps_throat_old;
	double T0_throat;	// Throat stagnation temperature (K)
	double p0_throat;	// Throat stagnation pressure (Pa)
	double A_throat, U_throat;
	int direction;
	char* direction_str;
	int V;

	int end;
	for(int S=0; S<N; S+=N-1)		// For each boundary node
	{
		if(S==0) end = 0;	// Enumerate labels
		else end = 1;

		//cout << "p0[S] = " << p0[S] << endl;
		//cout << "T0[S] = " << T0[S] << endl;

		if(!CONVERGED[end]) 
		{
			if(S==0)
			{
				V=S;		// Volume adjacent to odd end boundary	
				//if(rCOUNTER>100)
				//	tol[S] *= (double(rCOUNTER)/100);
			}
			else
			{
				V=S-1;		// Volume adjacent to even end boundary
				//if(rCOUNTER>100)
				//	tol[S] *= (double(rCOUNTER)/100);
			}

			// Having run the boundary method, get updated CL1, CL2, AA
			A_lambda = (Node[S].CL1[R+1] + Node[S].CL2[R+1])/2;
			Ts_lambda = pow(A_lambda,2)*(EX ? pPpt->TREFe : pPpt->TREFi);
			U_lambda = (Node[S].CL1[R+1] - Node[S].CL2[R+1])/(pPpt->gammaAir(Ts_lambda)-1);
			ps_lambda = pow(A_lambda/Node[S].AA[R+1], (2*pPpt->gammaAir(Ts_lambda))/(pPpt->gammaAir(Ts_lambda)-1))*pPpt->PREF*1e5;
			//p0_lambda = ps_lambda + 0.5*(ps_lambda/(pPpt->R_air*Ts_lambda))*pow(U_lambda*AREF,2);
			p0_lambda = TotalPressurePa(pPpt, ps_lambda, Ts_lambda, U_lambda*AREF);
			//T0_lambda = Ts_lambda + pow(U_lambda*AREF,2)/(2*pPpt->cpAir(Ts_lambda));
			T0_lambda = TotalTemperature(pPpt, Ts_lambda, U_lambda*AREF);

		
			if(p0[S] > FV[V].p0) // INFLOW
			{
				direction_str = "INFLOW";
				if(end==0) direction = 1; else direction = -1;
				p0_throat = p0[S];
				T0_throat = T0[S];
				if(rCOUNTER==0)
				{
					ps_throat[S] = FV[V].p0; // This is an intial guess - the lower bound of ps_throat
					ps_throat_old = ps_throat[S];	
				}
				else
				{
					ps_throat_old = ps_throat[S];
					
					//ps_throat[S] += 0.00001*(ps_lambda - ps_throat[S]);
					
					ps_throat[S] = ps_lambda;

				}
			}
			else
			{
				if(p0[S] < FV[V].p0) // OUTFLOW
				{
					direction_str = "OUTFLOW";
					if(end==0) direction = -1; else direction = 1;
					p0_throat = FV[V].p0;
					T0_throat = FV[V].T0;
					if(rCOUNTER==0) ps_throat[S] = p0[S];
				}
				else // NOFLOW
				{
					direction_str = "NOFLOW";
					// Either direction will do
					if(end==0) direction = 1; else direction = -1;
					p0_throat = p0[S];
					T0_throat = T0[S];
					if(rCOUNTER==0) ps_throat[S] = FV[V].p0;
				}
			}
		
			bool REVERSED = false;
			double* AU;
			AU = FillingAndEmptyingEquations(pPpt, T0_throat, p0_throat, ps_throat[S], REVERSED);
			A_throat = AU[0];
			U_throat = direction*AU[1];

			CL1_old = Node[S].CL1[R+1];
			CL2_old = Node[S].CL2[R+1];
			Node[S].CL1[R+1] = A_throat + ((pPpt->gammaAir(T0_throat)-1)/2)*U_throat;
			Node[S].CL2[R+1] = A_throat - ((pPpt->gammaAir(T0_throat)-1)/2)*U_throat;

			// Need to achieve ps_throat at Node[S], so calculate AA to do this
			//Node[S].AA[R+1] = A_throat/pow((p0_throat/1e5)/pPpt->PREF, (pPpt->gammaAir()-1)/(2*pPpt->gammaAir()));
			//Node[S].p_dash = pow(Node[S].A/Node[S].AA[R+1], (2*pPpt->gammaAir())/(pPpt->gammaAir()-1));

			//CONVERGED[end] = true;


///*
			if(S==0)
			{
			int coutprecision = cout.precision();
			cout.precision(16);
			cout << "Time = " << TIME << endl;
			cout << "Node S = " << S << endl;
			cout << "End = " << end << endl;
			cout << "rCOUNTER = " << rCOUNTER << endl;
			cout << endl;
			
			if(direction_str=="INFLOW") cout << "p0[S] (" << p0[S] << ") > FV[V].p0 (" << FV[V].p0 << ") --> " << direction_str << endl;
			else 
			{
				if(direction_str=="OUTFLOW") cout << "FV[V].p0 (" << FV[V].p0 << ") > p0[S] (" << p0[S] << ") --> " << direction_str << endl;
				else cout << "p0[S] (" << p0[S] << ") == FV[V].p0 (" << FV[V].p0 << ") --> " << direction_str << endl;
			}
			cout << endl;

			cout << "b.c. updated Node[S].CL1[R+1] = " << CL1_old << endl;
			cout << "b.c. updated Node[S].CL2[R+1] = " << CL2_old << endl;
			cout << endl;

			cout << "p0_throat = " << p0_throat << endl;
			cout << "T0_throat = " << T0_throat << endl;
			cout << "old ps_throat = " << ps_throat_old << endl;
			cout << "ps_lambda = " << ps_lambda << endl;
			cout << "new ps_throat = " << ps_throat[S] << endl;
			cout << endl;

			cout << "A_lambda = " << A_lambda << endl;
			cout << "A_throat = " << A_throat << endl;
			cout << endl;
			cout << "U_lambda = " << U_lambda << endl;
			cout << "U_throat = " << U_throat << endl;
			cout << endl;
			cout << "New Node[S].CL1[R+1] = " << Node[S].CL1[R+1] << endl;
			cout << "New Node[S].CL2[R+1] = " << Node[S].CL2[R+1] << endl;
			cout << endl;

			Node[S].CL1[R+1] = A_throat + ((pPpt->gammaAir(T0_throat)-1)/2)*U_throat;
			Node[S].CL2[R+1] = A_throat - ((pPpt->gammaAir(T0_throat)-1)/2)*U_throat;
			cout << endl;

			cout.precision(coutprecision);
			}
//*/
		}
	}
}

void CPipe::FillingAndEmptyingAdjust2(CProperties* pPpt, double DELZ, int& rCOUNTER, double TIME)
// Filling and emptying
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".FillingAndEmptyingAdjust2\n");}

	int V;
//	int direction;
	
	double A_lambda, U_lambda;
//	double ps_lambda, p0_lambda, Ts_lambda, T0_lambda;
	
//	double ps_throat;	// Throat static pressure (Pa);
//	double T0_throat;	// Throat stagnation temperature (K)
//	double p0_throat;	// Throat stagnation pressure (Pa)

	int end;
	for(int S=0; S<N; S+=N-1)		// For each boundary node
	{
		// Enumerate labels
		if(S==0) end = 0;
		else end = 1;

		if(!CONVERGED[end]) 
		{
			if(S==0)
			{
				V=S;		// Volume adjacent to odd end boundary
				
				if(rCOUNTER>100)
					tol[end] *= (double(rCOUNTER)/100);				
			}
			else
			{
				V=S-1;		// Volume adjacent to even end boundary
				
				if(rCOUNTER>100)
					tol[end] *= (double(rCOUNTER)/100);
			}
			
			// Having run the boundary method, use updated CL1, CL2, AA
			A_lambda = (Node[S].CL1[R+1] + Node[S].CL2[R+1])/2;
			U_lambda = (Node[S].CL1[R+1] - Node[S].CL2[R+1])/(pPpt->gammaAir(Node[S].T)-1);
			
			// Compare
			A_error_old[end] = A_error[end];
			U_error_old[end] = U_error[end];
			A_error[end] = A_lambda - A_throat[end];
			U_error[end] = U_lambda - U_throat[end];
/*			
			if(S==0)
			{
				int coutprecision = cout.precision();
				cout.precision(16);
				cout << "A_throat_old = " << A_throat[end] << endl;
				cout << "A_error = " << A_error[end] << endl;
				cout << endl;
				cout.precision(coutprecision);
			}
//*/
			if( (fabs(A_error[end]) < tol[end] && fabs(U_error[end]) < tol[end]) && rCOUNTER!=0 )
			{
				CONVERGED[end] = true;
			}
			else
			{			
				if(rCOUNTER==0)		// Start by assuming zero velocity, i.e. U==0
				{	
					A_throat[end] = pow((FV[V].p0/1e5)/pPpt->PREF, (pPpt->gammaAir(Node[S].T)-1)/(2*pPpt->gammaAir(Node[S].T)))*Node[S].AA[R+1];
					U_throat[end] = 0;
					
					A_throat_old[end] = A_throat[end];
					U_throat_old[end] = U_throat[end];
				}
				else
				{
					A_throat_old[end] = A_throat[end];
					U_throat_old[end] = U_throat[end];
				/*	
					//A_throat[end] = A_lambda;
					if(A_error[end]/A_error_old[end] < 0) // There has been a sign change
						beta[end] *= 0.5;
					A_throat[end] += beta[end]*(A_lambda - A_throat[end]);
				*/	
					if(U_error[end]/U_error_old[end] < 0) // There has been a sign change
						beta[end] *= 0.5;
					U_throat[end] += beta[end]*(U_lambda - U_throat[end]);
				}

				double Ts_throat;
				if((U_throat[end] > 0 && end == 0) || (U_throat[end] < 0 && end == 1)) // Inflow
				{
					double ps_throat = FV[V].p0;
					Ts_throat = FV[V].T0;

					A_throat[end] = sqrt(Ts_throat/(EX ? pPpt->TREFe : pPpt->TREFi));

					//	Node[S].T = pow(Node[S].A,2)*(EX ? pPpt->TREFe : pPpt->TREFi);

				}
				else // Outflow
				{
					double p0_throat = FV[V].p0;
					double T0_throat = FV[V].T0;

					Ts_throat = T0_throat - pow(U_throat[end]*AREF,2)/(2*pPpt->cpAir(T0_throat));

					A_throat[end] = sqrt(Ts_throat/(EX ? pPpt->TREFe : pPpt->TREFi));
				}
				
				// Construct new characteristics
				Node[S].CL1[R+1] = A_throat[end] + ((pPpt->gammaAir(Ts_throat)-1)/2)*U_throat[end];
				Node[S].CL2[R+1] = A_throat[end] - ((pPpt->gammaAir(Ts_throat)-1)/2)*U_throat[end];
			}
/*			
			//if(S==0)
			{
				int coutprecision = cout.precision();
				cout.precision(16);
				cout << "COUNTER = " << rCOUNTER << endl;
				//if(rCOUNTER==0)	
					cout << "Time = " << TIME << endl;
				cout << "Node S = " << S << endl;
				if(end==0) cout << "ODD end (" << end << ")" << endl;
				else cout << "EVEN end (" << end << ")" << endl;
				cout << "beta = " << beta[end] << endl;
				cout << "tol = " << tol << endl;
				cout << endl;
				if(!CONVERGED[end])
				{
					cout << "A_throat_old = " << A_throat_old[end] << endl;
					cout << "A_lambda = " << A_lambda << endl;
					//cout << "A_lambda - A_throat_old  = " << A_lambda - A_throat_old[end] << endl;
					cout << endl;
					cout << "A_error = " << A_error[end] << endl;
					//cout << "A_error_old = " << A_error_old[end] << endl;
					cout << endl;
					cout << "U_throat_old = " << U_throat_old[end] << endl;
					cout << "U_lambda = " << U_lambda << endl;
					//cout << "U_lambda - U_throat_old  = " << U_lambda - U_throat_old[end] << endl;
					cout << endl;
					cout << "U_error = " << U_error[end] << endl;
					//cout << "U_error_old = " << U_error_old[end] << endl;
					cout << endl;
					cout << "new A_throat = " << A_throat[end] << endl;
					cout << "new U_throat = " << U_throat[end] << endl;
					cout << endl;
				}
				else
				{
					cout << "Final A_error = " << A_error[end] << endl;
					//cout << endl;
					cout << "Final U_error = " << U_error[end] << endl;
					cout << endl;
					cout << "Final A_throat = " << A_throat[end] << endl;
					cout << "Final U_throat = " << U_throat[end] << endl;
					cout << endl;
				}

				if(CONVERGED[end]) cout << "CONVERGED!" << endl;
				else cout << "NOT converged" << endl;
				cout << endl;
				cout << endl;
				cout.precision(coutprecision);
			}
//*/
		}

/*
		//if(S==0 && TIME>=0.5121)
		//if(S==0)
		{
			if(CONVERGED[end])
			{
				int coutprecision = cout.precision();
				cout.precision(16);
				cout << "COUNTER = " << rCOUNTER << endl;
				//if(rCOUNTER==0)	
					cout << "Time = " << TIME << endl;
				cout << "Node S = " << S << endl;
				if(end==0) cout << "ODD end (" << end << ")" << endl;
				else cout << "EVEN end (" << end << ")" << endl;
				cout << "A_throat = " << A_throat[end] << endl;
				cout << "U_throat = " << U_throat[end] << endl;
				cout << endl;
				cout << "CONVERGED!" << endl;
				cout << endl;
				cout << endl;
				cout.precision(coutprecision);
			}
		}
//*/
	}
}

void CPipe::FillingAndEmptyingSimple(CProperties* pPpt)//, double DELZ, int& rCOUNTER, double TIME)
// Filling and emptying
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".FillingAndEmptyingSimple\n");}

	int S, end, V;

	// Use lambda_in, lambda_out at boundary nodes to obtain values for A and U
	for(S=0; S<N; S+=N-1)		// For each boundary node
	{
		end = (S==0?0:1);
		if(end==0) V=S;	// Volume adjacent to odd end boundary
		else V=S-1;		// Volume adjacent to even end boundary
		
		double del_p_ideal = 0.5*Node[S].rho*pow(Node[S].U*(EX ? pPpt->AREFe : pPpt->AREFi),2)/1e5; // In bar

		double pX_temp;
		
		pX_temp = FV[V].p0/1e5; // pX_temp in bar

/*
		if( (end==0 && Node[S].U > 0) || (end==1 && Node[S].U < 0) ) // If flow is into this pipe
		{
			pX_temp = FV[V].p0/1e5 - del_p_ideal;

			// Expand to the corrected pressure
		}
		else pX_temp = FV[V].p0/1e5;
*/
		
		// Treat inflow as flow into a reservoir at stagnation pressure p0, i.e. ps==p0 and U=0
		Node[S].CL1[R+1] = pow((pX_temp)/pPpt->PREF, (pPpt->gammaAir(Node[S].T)-1)/(2*pPpt->gammaAir(Node[S].T)))*Node[S].AA[R+1];
		Node[S].CL2[R+1] = pow((pX_temp)/pPpt->PREF, (pPpt->gammaAir(Node[S].T)-1)/(2*pPpt->gammaAir(Node[S].T)))*Node[S].AA[R+1];	
	}
}

void CPipe::FillingAndEmptyingAdjust(CProperties* pPpt, double DELZ, int& rCOUNTER, double TIME)
// Filling and emptying
{
	if(pPpt->SHOW_calls){pPpt->Out(Identify()); pPpt->Out(".FillingAndEmptyingAdjust\n");}

	int end, V;
		
	double A_lambda, U_lambda;
	double ps_lambda, p0_lambda, Ts_lambda, T0_lambda;

	double CL1_old, CL2_old;

	// Use lambda_in, lambda_out at boundary nodes to obtain values for A and U
	for(int S=0; S<N; S+=N-1)		// For each boundary node
	{
		end = (S==0?0:1);
		
		if(!CONVERGED[end])
		{
			if(end==0)
			{
				V=S;	// Volume adjacent to odd end boundary
			}
			else
			{
				V=S-1;		// Volume adjacent to even end boundary
			}

			// Having run the boundary method, use updated CL1, CL2, AA
			A_lambda = (Node[S].CL1[R+1] + Node[S].CL2[R+1])/2;
			Ts_lambda = pow(A_lambda,2)*(EX ? pPpt->TREFe : pPpt->TREFi);
			U_lambda = (Node[S].CL1[R+1] - Node[S].CL2[R+1])/(pPpt->gammaAir(Ts_lambda)-1);
			ps_lambda = pow(A_lambda/Node[S].AA[R+1], (2*pPpt->gammaAir(Ts_lambda))/(pPpt->gammaAir(Ts_lambda)-1))*pPpt->PREF*1e5;
			//p0_lambda = ps_lambda + 0.5*(ps_lambda/(pPpt->R_air*Ts_lambda))*pow(U_lambda*AREF,2);
			p0_lambda = TotalPressurePa(pPpt, ps_lambda, Ts_lambda, U_lambda*AREF);
			//T0_lambda = Ts_lambda + pow(U_lambda*AREF,2)/(2*pPpt->cpAir(Ts_lambda));
			T0_lambda = TotalTemperature(pPpt, Ts_lambda, U_lambda*AREF);

			CL1_old = Node[S].CL1[R+1];// Save old values
			CL2_old = Node[S].CL2[R+1];

			
			if((end==0 && U_lambda > 0) || (end!=0 && U_lambda < 0))
			{
				// Flow is into the FV
				// -------------------
				if(direction_str[end] != "INFLOW") rCOUNTER = 0; // If switched directions
				direction_str[end] = "INFLOW";
			
///*
				// Treat inflow as flow into a reservoir at stagnation pressure p0, i.e. ps==p0 and U=0
				Node[S].CL1[R+1] = pow((FV[V].p0/1e5)/pPpt->PREF, (pPpt->gammaAir(Ts_lambda)-1)/(2*pPpt->gammaAir(Ts_lambda)))*Node[S].AA[R+1];
				
				
				Node[S].CL2[R+1] = pow((FV[V].p0/1e5)/pPpt->PREF, (pPpt->gammaAir(Ts_lambda)-1)/(2*pPpt->gammaAir(Ts_lambda)))*Node[S].AA[R+1];	
				//Node[S].CL2[R+1] = pow(( (FV[V].p0 - 0.5*Node[S].rho*pow(U_lambda*AREF,2)) /1e5)/pPpt->PREF, (pPpt->gammaAir()-1)/(2*pPpt->gammaAir()))*Node[S].AA[R+1];	
				


				CONVERGED[end] = true;
//*/	
/*
			
				// Set p0_throat = p0_lambda, T0_throat = T0_lambda, converge on ps_throat
				p0_throat = p0_lambda;
				T0_throat = T0_lambda;						
				
				if(rCOUNTER==0) ps_throat[end] = p0_throat;//FV[V].p0;//p0_throat*0.99;

				if(end==0)// Odd end boundary
				{	
					direction = 1;
				}
				else	// Even end boundary
				{
					direction = -1;
				}

				bool REVERSED = false;
				double* AU;
				AU = FillingAndEmptyingEquations(pPpt, T0_throat, p0_throat, ps_throat[end], REVERSED);
				A_throat[end] = AU[0];
				U_throat[end] = direction*AU[1];

				// Convergence test
				// ----------------
				double error = fabs(ps_throat[end] - ps_lambda);
								
				if(rCOUNTER>100)
				{
					tol[end] *= (double(rCOUNTER)/100);
				}
					
				if(error < tol[end])
				{
					CONVERGED[end] = true;								
				}
				else
				{				
					CL1_old = Node[S].CL1[R+1];
					Node[S].CL1[R+1] = A_throat[end] + ((pPpt->gammaAir()-1)/2)*U_throat[end];
					CL2_old = Node[S].CL2[R+1];
					Node[S].CL2[R+1] = A_throat[end] - ((pPpt->gammaAir()-1)/2)*U_throat[end];
					
					ps_throat_old[end] = ps_throat[end];
					ps_throat[end] -= (1-beta[end])*(ps_throat[end] - ps_lambda);
					
					while(ps_throat[end] > p0_throat)
					{
						//cout << "ps_throat[end] = " << ps_throat[end] << endl;
						//cout << "p0_throat = " << p0_throat << endl;
						//cout << "beta[end] = " << beta[end] << endl;
						//cout << endl;
						
						beta[end] = sqrt(beta[end]);
						ps_throat[end] = ps_throat_old[end];
						ps_throat[end] -= (1-beta[end])*(ps_throat[end] - ps_lambda);
					}
				}
//*/
			}
			else													
			{
				if((end==0 && U_lambda < 0) || (end!=0 && U_lambda > 0))
				{
					// Flow is toward the boundary
					// ---------------------------
					if(direction_str[end] != "OUTFLOW") rCOUNTER = 0; // If switched directions
					direction_str[end] = "OUTFLOW";
			
///*
					// Treat inflow as flow into a reservoir at stagnation pressure p0, i.e. ps==p0 and U=0
					Node[S].CL1[R+1] = pow((FV[V].p0/1e5)/pPpt->PREF, (pPpt->gammaAir(Ts_lambda)-1)/(2*pPpt->gammaAir(Ts_lambda)))*Node[S].AA[R+1];
					Node[S].CL2[R+1] = pow((FV[V].p0/1e5)/pPpt->PREF, (pPpt->gammaAir(Ts_lambda)-1)/(2*pPpt->gammaAir(Ts_lambda)))*Node[S].AA[R+1];	
				
					CONVERGED[end] = true;
//*/	
/*
					p0_throat = FV[V].p0; // As per Connor & Swain method
					T0_throat = FV[V].T0; // As per Connor & Swain method						
					
					if(rCOUNTER==0) ps_throat[end] = p0_throat*0.99;

					if(end==0)// Odd end boundary
					{	
						direction = -1;
					}
					else	// Even end boundary
					{
						direction = 1;
					}
					
					bool REVERSED = false;
					double* AU;
					AU = FillingAndEmptyingEquations(pPpt, T0_throat, p0_throat, ps_throat[end], REVERSED);
					A_throat[end] = AU[0];
					U_throat[end] = direction*AU[1];	
					
					// Convergence test
					// ----------------
					double error = fabs(ps_throat[end] - ps_lambda);
					
					if(rCOUNTER>100)
					{
						tol[end] *= (double(rCOUNTER)/100);
					}
					
					if(error < tol[end])
					{
						CONVERGED[end] = true;								
					}
					else
					{				
						CL1_old = Node[S].CL1[R+1];
						Node[S].CL1[R+1] = A_throat[end] + ((pPpt->gammaAir()-1)/2)*U_throat[end];
						CL2_old = Node[S].CL2[R+1];
						Node[S].CL2[R+1] = A_throat[end] - ((pPpt->gammaAir()-1)/2)*U_throat[end];
						
						ps_throat_old[end] = ps_throat[end];
						ps_throat[end] -= (1-beta[end])*(ps_throat[end] - ps_lambda);
						
						while(ps_throat[end] > p0_throat)
						{
							
							//cout << "ps_throat[end] = " << ps_throat[end] << endl;
							//cout << "p0_throat = " << p0_throat << endl;
							//cout << "beta[end] = " << beta[end] << endl;
							//cout << endl;
							
							beta[end] = sqrt(beta[end]);
							ps_throat[end] = ps_throat_old[end];
							ps_throat[end] -= (1-beta[end])*(ps_throat[end] - ps_lambda);
						}
					}
//*/					
				}
				else	// No flow
				{
					direction_str[end] = "NOFLOW";
					Node[S].CL1[R+1] = pow((FV[V].p0/1e5)/pPpt->PREF, (pPpt->gammaAir(Ts_lambda)-1)/(2*pPpt->gammaAir(Ts_lambda)))*Node[S].AA[R+1];
					Node[S].CL2[R+1] = pow((FV[V].p0/1e5)/pPpt->PREF, (pPpt->gammaAir(Ts_lambda)-1)/(2*pPpt->gammaAir(Ts_lambda)))*Node[S].AA[R+1];					
					CONVERGED[end] = true;
				}
			}
		}	
/*
//if(!CONVERGED[end])
if(direction_str[end]=="INFLOW")
{
	int coutprecision = cout.precision();
	cout.precision(16);
	cout << "S = " << S << endl;
	cout << "end = " << end << endl;
	cout << "TIME = " << TIME << endl;
	cout << "rCOUNTER = " << rCOUNTER << endl;
	cout << endl;
	cout << "beta = " << beta[end] << endl;
	cout << "tol = " << tol[end] << endl;
	cout << endl;

	if(direction_str[end]=="INFLOW")
	{
		cout << "ps_throat = " << ps_throat_old[end] << endl;
		cout << "ps_lambda = " << ps_lambda << endl;
		cout << "new ps_throat = " << ps_throat[end] << endl;
		cout << "FV[V].p0 = " << FV[V].p0 << endl;
		cout << "fabs(ps_throat - ps_lambda) = " << fabs(ps_throat[end] - ps_lambda) << endl;

		//cout << "p0_lambda = " << p0_lambda << endl;
		//cout << "p0_throat = " << p0_throat << endl;
		//cout << "fabs(p0_throat - p0_lambda) = " << fabs(p0_throat - p0_lambda) << endl;
		//cout << "ps_throat (constant) = " << ps_throat[end] << endl;
		//cout << endl;
		//cout << "T0_lambda = " << T0_lambda << endl;
		//cout << "T0_throat = " << T0_throat << endl;
		//cout << "fabs(T0_throat - T0_lambda) = " << fabs(T0_throat - T0_lambda) << endl;
	}

	if(direction_str[end]=="OUTFLOW")
	{
		cout << "ps_lambda = " << ps_lambda << endl;
		cout << "ps_throat = " << ps_throat[end] << endl;
		cout << "FV[V].p0 = " << FV[V].p0 << endl;
		cout << "fabs(ps_throat - ps_lambda) = " << fabs(ps_throat[end] - ps_lambda) << endl;
	}
	cout << endl;
	
	cout << "A_lambda\t\t=\t" << A_lambda << endl;
	cout << "A_throat\t\t=\t" << A_throat[end] << endl;
	cout << endl;
	cout << "U_lambda\t\t=\t" << U_lambda << endl;
	cout << "U_throat\t\t=\t" << U_throat[end] << endl;
	cout << endl;
	
	if(CONVERGED[end]) cout << "CONVERGED!" << endl;
	else cout << "NOT CONVERGED" << endl;
	cout << endl;
	cout << "Flow direction = " << direction_str[end] << endl;
	cout << endl;
	cout << endl;
	cout.precision(coutprecision);
}
//*/

	}
}

char* CPipe::Identify()
// ============================================================ //
// Returns identification of the current object					//
// ============================================================ //
{
	std::string sz;
	sz = "Assembly [";
	sz += IntToString(AssyID);
	sz += "], ";
	if(EX) sz += "Exhaust "; else sz += "Intake ";
	sz += "Pipe [";
	sz += IntToString(ID);
	sz += "]";

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str());
	return szz;
	delete [] szz;
}

double CPipe::d(double x)
// ============================================================ //
// Returns diameter d for the required location x				//
// ============================================================ //
{
	if(LINEAR && LINEAR_F) // If constant area variation (dfdx = constant)
	{
		double f_odd = PI*pow(d_odd,2)/4;
		double f_even = PI*pow(d_even,2)/4;
		double dfdx = (f_even - f_odd)/this->length;
		return sqrt((dfdx*x + f_odd)*4/PI);
	}
	else // For linear and quadratic pipe diameter variation (ax^2 + bx + c), a, b, and c are member variables calculated in the initialisation
	{
		cout << "Quadratic variation disabled. Exiting." << endl;
		exit(1);
		//return a*pow(x,2) + b*x + c;		
	}
}

double CPipe::dddx(CProperties* pPpt, double x)
// ============================================================ //
// Returns diameter gradient dddx for the required location x	//
// ============================================================ //
{
	// Both cases below assume circular cross-section:
	double dddx;
	if(LINEAR && LINEAR_F) // If constant area variation (dfdx = constant)
	{
		double f_odd = PI*pow(d_odd,2)/4;
		double f_even = PI*pow(d_even,2)/4;
		double dfdx = (f_even - f_odd)/this->length;

		dddx = 0.5*pow((dfdx*x + f_odd)*(4/PI),-0.5)*((4/PI)*dfdx); // Expression for dddx
	}
	else // For linear and quadratic pipe diameter variation (ax^2 + bx + c), a, b, and c are member variables calculated in the initialisation
	{
		cout << "Quadratic variation disabled. Exiting." << endl;
		exit(1);
		//dddx = 2*a*x + b; // Expression for dddx		
	}
	
	if(dddx < pPpt->ZERO_TOL) return 0; else return dddx;
}

double CPipe::d2ddx2(CProperties* pPpt, double x)
// ============================================================ //
// Returns second derivative of diameter variation				//
// for the required location x									//
// ============================================================ //
{
	// Both cases below assume circular cross-section:
	double d2ddx2;
	if(LINEAR && LINEAR_F) // If constant area variation (dfdx = constant)
	{
		double f_odd = PI*pow(d_odd,2)/4;
		double f_even = PI*pow(d_even,2)/4;
		double dfdx = (f_even - f_odd)/this->length;

		d2ddx2 = 0.5*(4/PI)*dfdx * ( -0.5*pow((dfdx*x + f_odd)*(4/PI),-1.5)*((4/PI)*dfdx) ); // Expression for d2ddx2
		cout << "d2ddx2 = " << d2ddx2 << "\n";
	}
	else // For linear and quadratic pipe diameter variation (ax^2 + bx + c), a, b, and c are member variables calculated in the initialisation
	{
		cout << "Quadratic variation disabled. Exiting." << endl;
		exit(1);
		//d2ddx2 = 2*a; // Expression for d2ddx2		
	}

	if(d2ddx2 < pPpt->ZERO_TOL) return 0; else return d2ddx2;
}

double CPipe::f(double x)
// ============================================================ //
// Returns area f for the required location x					//
// ============================================================ //
{
	if(LINEAR && LINEAR_F) // If constant area variation (dfdx = constant)
	{
		double f_odd = PI*pow(d_odd,2)/4;
		double f_even = PI*pow(d_even,2)/4;
		double dfdx = (f_even - f_odd)/this->length;
		return dfdx*x + f_odd;
	}
	else // For linear and quadratic pipe diameter variation (ax^2 + bx + c), a, b, and c are member variables calculated in the initialisation
	{
		// For pipe diameter = ax^2 + bx + c (for linear, a=0), this integrates to give a circular cross-sectional area:
		cout << "Quadratic variation disabled. Exiting." << endl;
		exit(1);
		//return (PI/4)*( pow(a,2)*pow(x,4) + 2*a*b*pow(x,3) + (2*a*c + pow(b,2))*pow(x,2) + 2*b*c*x + pow(c,2) );
	}
}

double CPipe::dfdx(double x)
// ============================================================ //
// Returns area gradient dfdx for the required location x		//
// ============================================================ //
{
	if(LINEAR && LINEAR_F) // If constant area variation (dfdx = constant)
	{
		double f_odd = PI*pow(d_odd,2)/4;
		double f_even = PI*pow(d_even,2)/4;
		return (f_even - f_odd)/this->length;
	}
	else // For linear and quadratic pipe diameter variation (ax^2 + bx + c), a, b, and c are member variables calculated in the initialisation
	{
		// The expression for pipe circular cross-sectional area differentiates to give the area gradient:
		//return (PI/2)*( 2*pow(a,2)*pow(x,3) + 3*a*b*pow(x,2) + (2*a*c + pow(b,2))*x + b*c );
		cout << "Quadratic variation disabled. Exiting." << endl;
		exit(1);
		//return (PI/2)*d(x)*(2*a*x + b);
	}
}