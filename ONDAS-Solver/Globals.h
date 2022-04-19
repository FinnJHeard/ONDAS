//#include "Time.h"

// Enumerate BOUNDARY CONDITIONS
#define UNKNOWN 0
#define INTERIOR 1
#define RESERVOIR 2
#define JUNCTION 3
#define OPEN 4
#define CLOSED 5
#define NOZZLE 6
#define INFLOWEND 7
#define CYLINDER 8
#define VALVE 9 //#define EXHVALVE 9
//#define INTVALVE 10
#define SUDDEN_LEFT 11
#define SUDDEN_RIGHT 12
#define APLDEV_LEFT 13
#define APLDEV_RIGHT 14
#define TURB_INLET_OR_OUTLET 15
#define TURBINLET 16
#define TURBOUTLET 17
#define COMPINLET 18
#define COMPOUTLET 19
#define TRANSMISSIVE 20
#define END_ENVIRONMENT 21
#define ASSEMBLY 22
#define VOLUMEVALVE 23
#define ANECHOIC 24

#define EXHAUST 25
#define INTAKE 26
#define ODD 33
#define INSIDE 38
#define EVEN 44

#define INFLOW 61			// Inflow INTO pipe from boundary
#define OUTFLOW 62			// Outflow from pipe towards boundary
#define NOFLOW 63

#define HCPJID 71
#define NHCPJID 72
#define NHPLTJID 73
#define NHPLPCID 74
#define NHPLJID 75

#define ENLARGEMENT 81		// Flow pertains to sudden enlargement
#define CONTRACTION 82		// Flow pertains to sudden contraction
#define EQUAL_AREAS 83		// Neither enalargement nor contraction

// Cylinder model
// ==============
#define CONST_P 90			// Constant conditions cylinder
#define CONST_VOL 91			// Constant volume, variable conditions cylinder
#define READ_P_T 92			// P & T read from file
#define RESET_P_T 93		// P & T reset at EVO or TDC
#define WATSON 94			// Watson single-zone

// Cylinder heat transfer model
// ============================
#define NONE 100	
#define ANNAND 101		
#define WOSCHNI 102	

// LABELS
const int R = 0;			// As in TIME[R+1] - or Timestep[R+1] (new or old)
const int ONE_SIDE = 0;		// Used in one-sided boundary conditions, e.g. nozzle	
const int LEFT_SIDE = 0;	// Used in double-sided boundary conditions, e.g. apldev 
const int RIGHT_SIDE = 1;
const int EXHAUST_SIDE = 0;			// To be sure of identifying the correct item in the two item configuration vectors
const int INTAKE_SIDE = 1;

// APL Device labels
const int GAUZE = 1;
const int THROTTLE = 2;
const int EGR = 3;

// DEFINITION OF PI
//const double PI = 3.14159265358979323;
const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

// DEFINITION OF e
const double constant_e = 2.718281828;