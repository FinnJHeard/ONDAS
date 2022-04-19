
#include "Properties.h"
#include "Time.h"
#include <limits>

void PrintAllCharacterCodes(CProperties* pPpt);
void Message(CProperties* pPpt, int build);
void SimulationOutput(CProperties* pPpt, int timestep, double time, double elapsedTime, bool& rCONTINUING);
void SimulationOutput(CProperties* pPpt, int timestep, double time, double elapsedTime, bool& rCONTINUING, double ca_elapsed);
void SimulationWarning(CProperties* pPpt,char* strWarning);
void Line(CProperties* pPpt, int length, int code);
void Line(CProperties* pPpt, int length, char character);
void Line(CProperties* pPpt, int length, char* character);
//void Line(CProperties &rPpt, int length, int code);
//void ReadInput1(char *InputFile, double* &rParams);
void GetHeader(char* &rheader, char* label, int n);

void CommonReadInput(char *InputFile, int num_parameters, char** &rlabels, double* &rvalues, char** &rstrings, int &rlast_entry);

char* ConstructString(CProperties* pPpt, std::string dir_str, std::string bcname_str, bool ex, int id);
char* ConstructString(CProperties* pPpt, std::string dir_str, std::string bcname_str, int id);
char* ConstructString(CProperties* pPpt, std::string dir_str, std::string name_str, int id, std::string ext_str);
char* ConstructString(CProperties* pPpt, std::string dir_str, std::string name_str, bool DOT_TEXT);
char* ConstructString(CProperties* pPpt, std::string dir_str, std::string name_str);

char* GetBoundaryName(int NAME);

char* StringToChar(std::string sz);
std::string IntToString(int num);
std::string CharToString(char* charstar);

char* Line(char* title_str, char* line_char);
char* Underline(char* title_str, char* line_char);
char* Underline(char* title_str, char* line_char, bool OVERLINE);
char* Underline(char* title_str, char* line_char, bool OVERLINE, int line_length, int indent_tabs);
char* Underline(char* title_str, char* line_char, bool OVERLINE, char* prefix);
char* Underline(char* title_str, char* line_char, int id);
char* Underline(char* title_str, char* line_char, char* prefix);
char* Underline(char* title_str, char* line_char, char* prefix, char* postfix);
char* Underline(char* title_str, char* line_char, int id, char* prefix);
char* Underline(string title_str, char* line_char);
char* Underline(string title_str, char* line_char, bool OVERLINE);
//char* Underline(string title_str, char* line_char, bool OVERLINE, char* prefix);
//char* Underline(string title_str, char* line_char, int id);
//char* Underline(string title_str, char* line_char, char* prefix);
//char* Underline(string title_str, char* line_char, char* prefix, char* postfix);
//char* Underline(string title_str, char* line_char, int id, char* prefix);

double AminusB(CProperties* pPpt, double A, double B);
void Interpolate(double xLower, double yLower, double xUpper, double yUpper, double xInterp, double& rYInterp);

double TotalPressurePa(CProperties* pPpt, double ps_Pa, double Ts, double u);
double TotalPressureBar(CProperties* pPpt, double ps_bar, double Ts, double u);
double TotalTemperature(CProperties* pPpt, double Ts, double u);
double StaticTemperature(CProperties* pPpt, double T0, double u);

char* InflowOrOutflow(int flowDir);
int InflowOrOutflowNumber(int flowDir);

inline bool DoubleToBool(double inDouble){if(fabs(inDouble - 1.0) < 1e-6) return true; else return false;}
inline char* TrueOrFalse(bool True){if(True) return "true"; else return "false";}
inline char* Choked(bool CHOKED){if(CHOKED) return "CHOKED"; else return "NOT CHOKED";}
inline int ChokedNo(bool CHOKED){return CHOKED;}
inline bool IsEven(int integer){if(integer%2==0) return true; else return false;}
inline char Deg(){return char(248);}

void PrecisionInformation();
