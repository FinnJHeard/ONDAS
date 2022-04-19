// ====================================================================================================100
// ===============================================================================================95
// ==========================================================================================90
// ============================================================60
// ----------------------------------------------------------------------------------------------------100
// --------------------------------------------------------------------------------80                    
// --------------------------------------------------50

#include "Time.h"
#include "StdAfx.h"

#include "Globals.h"
#include "Tools.h"

#include <sstream>  //this header automatically includes iostream and is needed for ostringstream 
#include <string> 

//using namespace std;

//#include <math.h>

void PrintAllCharacterCodes(CProperties* pPpt)
{
	for(int c=0; c<255; ++c){cout << c << "\t=\t"; Line(pPpt, 1, c);}
	return;
}

void Message(CProperties* pPpt, int build)
{
	// Codes
	// 1 - smiley face
	// 2 - smiley face
	// 3 - heart
	// 4 - diamond
	// 5 - club
	// 6 - spade
	// 22 - thick
	// 45 - dash
	// 176 - shading diagonally down
	// 177 - shading diagonally down
	// 178 - shading diagonally up
	// 224 - alpha
	// 225 - beta
	// 248 - degrees
	// 254 - blocks
/*
	for(int c=0; c<255; ++c){cout << c << "\t=\t"; Line(pPpt, 1, c);}
	exit(1);
//*/
	int length = 80;//50
	//Line(pPpt, length, 254);
	Line(pPpt, length, "%");
/*
	cout << endl;
	cout << "    ...    ...    ...     OOO   N   N  DDDD   AAAA   SSSSS  ...    ...    ...   " << endl;
	cout << "   ....   ....   ....   O   O  NN  N  D   D  A   A  S      ....   ....   ....   " << endl;
	cout << ".......................O   O  N N N  D   D  AAAAA  SSSSS........................" << endl;
	cout << "  ....   ....   ....  O   O  N  NN  D   D  A   A      S   ....   ....   ....    " << endl;
	cout << "  ...    ...    ...   OOO   N   N  DDDD   A   A  SSSSS    ...    ...    ...     " << endl;
	cout << endl;
	cout << "              ...:::ONe-Dimensional wave Action Simulator:::...                     " << endl;
	cout << endl;
*/
	pPpt->Out("\n");
	pPpt->Out("    ...    ...    ...     OOO   N   N  DDDD   AAAA   SSSSS  ...    ...    ...   "); pPpt->Out("\n");
	pPpt->Out("   ....   ....   ....   O   O  NN  N  D   D  A   A  S      ....   ....   ....   "); pPpt->Out("\n");
	pPpt->Out(".......................O   O  N N N  D   D  AAAAA  SSSSS........................"); pPpt->Out("\n");
	pPpt->Out("  ....   ....   ....  O   O  N  NN  D   D  A   A      S   ....   ....   ....    "); pPpt->Out("\n");
	pPpt->Out("  ...    ...    ...   OOO   N   N  DDDD   A   A  SSSSS    ...    ...    ...     "); pPpt->Out("\n");
	pPpt->Out("\n");
	pPpt->Out("              ...:::ONe-Dimensional wave Action Simulator:::...                     "); pPpt->Out("\n");
	pPpt->Out("\n");
	
	//Line(pPpt, length, 254);
	Line(pPpt, length, "%");
	//cout << endl;
	pPpt->Out("\n");

	int year, month, day;
	year = 2000 + (build - build%10000)/10000;
	month = (build%10000 - (build%10000)%100)/100;
	day = (build%10000)%100;

	char* str_month;
	switch(month)
	{
		case 1:
			str_month = "January";
			break;
		case 2:
			str_month = "February";
			break;
		case 3:
			str_month = "March";
			break;
		case 4:
			str_month = "April";
			break;
		case 5:
			str_month = "May";
			break;
		case 6:
			str_month = "June";
			break;
		case 7:
			str_month = "July";
			break;
		case 8:
			str_month = "August";
			break;
		case 9:
			str_month = "September";
			break;
		case 10:
			str_month = "October";
			break;
		case 11:
			str_month = "November";
			break;
		case 12:
			str_month = "December";
			break;
		default:
			str_month = "Invalid";
			break;
	}
	
	pPpt->Out("  This is ONDAS v1.0.0, build date "); pPpt->Out(day); pPpt->Out(" "); pPpt->Out(str_month); pPpt->Out(", "); pPpt->Out(year); pPpt->Out("\n");

	// ==========================================
	// Revision record
	// ==========================================
	// Code version		Version		Build Date
	// ------------		-------		----------
	// 90527			v0.5.1		May  27, 2009
	// 100913			v0.5.1		Sept 13, 2010
	// 140227			v0.6.0		Feb  27, 2014
	// 200721			v1.0.0		21 Jul, 2020
	//
	// ==========================================

	pPpt->Out("\n");
	pPpt->Out("  By:"); pPpt->Out("\n\n");
	pPpt->Out("  Dr Aaron Costall, Imperial College London, Mechanical Engineering Dept."); pPpt->Out("\n");
	pPpt->Out("  a.costall@imperial.ac.uk"); pPpt->Out("\n");
	pPpt->Out("\n");

	//Line(pPpt, length, 254);
	//Line(pPpt, length, "%");
//	for(int i=0; i<3e8; ++i){;}// Pause for a second
}
void SimulationOutput(CProperties* pPpt, int timestep, double time, double elapsedTime, bool& rCONTINUING)
{
	int temp_precision = cout.precision(); cout << setprecision(6);
	if(timestep==0 || rCONTINUING)
	{
		rCONTINUING = false; // Reset
		pPpt->Out("\n\n");
		pPpt->Out("......................................................ONDAS......................................................\n");
	}
	pPpt->Out(".                                                                                                               .\n");
	pPpt->Out(".                                                                                                               .\n");
	pPpt->Out("."); 
		pPpt->Out("\tTIMESTEP: "); 
		if(timestep<1000000) pPpt->Out(" ");
		if(timestep<100000) pPpt->Out(" ");
		if(timestep<10000) pPpt->Out(" ");
		if(timestep<1000) pPpt->Out(" ");
		if(timestep<100) pPpt->Out(" ");
		if(timestep<10) pPpt->Out(" ");
		//if(timestep<1) pPpt->Out(" ");
		pPpt->Out(timestep);
		pPpt->Out("\t\t");

		pPpt->Out("SIMULATION TIME: "); 
		if(time==0) pPpt->Out("0.000000");
		else pPpt->Out(time);
		pPpt->Out("s\t\t");

		pPpt->Out("ELAPSED TIME: "); 
		if(elapsedTime<1000) pPpt->Out(" ");
		if(elapsedTime<100) pPpt->Out(" ");
		if(elapsedTime<10) pPpt->Out(" ");
		//if(elapsedTime<1) pPpt->Out(" ");
		pPpt->Out(elapsedTime); 
		pPpt->Out("s\t\t");
		pPpt->Out(".\n");
	pPpt->Out(".                                                                                                               .\n");
	pPpt->Out(".                                                                                                               .\n");
	pPpt->Out("......................................................ONDAS......................................................\n");
	cout << setprecision(temp_precision);

	/*pPpt->Out("Timestep: "); pPpt->Out(timestep); pPpt->Out(";\t\tsimulation time: "); pPpt->Out(time);
	pPpt->Out("s;\t\telapsed time: "); pPpt->Out(elapsedTime); pPpt->Out("s\n");*/
}

void SimulationOutput(CProperties* pPpt, int timestep, double time, double elapsedTime, bool& rCONTINUING, double ca_elapsed)
{
	int temp_precision = cout.precision(); cout << setprecision(6);
	if (timestep == 0 || rCONTINUING)
	{
		rCONTINUING = false; // Reset
		pPpt->Out("\n\n");
		pPpt->Out("......................................................ONDAS......................................................\n");
	}
	pPpt->Out(".                                                                                                               .\n");
	pPpt->Out(".                                                                                                               .\n");
	pPpt->Out(".");
	pPpt->Out("  TIMESTEP: "); //pPpt->Out("\tTIMESTEP: ");
	if (timestep < 1000000) pPpt->Out(" ");
	if (timestep < 100000) pPpt->Out(" ");
	if (timestep < 10000) pPpt->Out(" ");
	if (timestep < 1000) pPpt->Out(" ");
	if (timestep < 100) pPpt->Out(" ");
	if (timestep < 10) pPpt->Out(" ");
	//if(timestep<1) pPpt->Out(" ");
	pPpt->Out(timestep);
	pPpt->Out("\t"); //pPpt->Out("\t\t");

	pPpt->Out("SIM TIME: ");
	if (time == 0) pPpt->Out("0.000000");
	else pPpt->Out(time);
	pPpt->Out("s\t"); //pPpt->Out("s\t\t");

	pPpt->Out("SIM ELAPSED CA: ");
	if (time == 0) pPpt->Out("0.000000");
	else pPpt->Out(ca_elapsed);
	pPpt->Out(" deg\t");

	pPpt->Out("REAL TIME ELAPSED: "); //pPpt->Out("ELAPSED TIME: ");
	if (elapsedTime < 1000) pPpt->Out(" ");
	if (elapsedTime < 100) pPpt->Out(" ");
	if (elapsedTime < 10) pPpt->Out(" ");
	//if(elapsedTime<1) pPpt->Out(" ");
	pPpt->Out(elapsedTime);
	pPpt->Out("s\t"); //pPpt->Out("s\t\t");
	pPpt->Out(".\n");
	pPpt->Out(".                                                                                                               .\n");
	pPpt->Out(".                                                                                                               .\n");
	pPpt->Out("......................................................ONDAS......................................................\n");
	cout << setprecision(temp_precision);

	/*pPpt->Out("Timestep: "); pPpt->Out(timestep); pPpt->Out(";\t\tsimulation time: "); pPpt->Out(time);
	pPpt->Out("s;\t\telapsed time: "); pPpt->Out(elapsedTime); pPpt->Out("s\n");*/
}

void SimulationWarning(CProperties* pPpt, char* strWarning)
{
	int temp_precision = cout.precision(); cout << setprecision(6);
	//if(timestep==0 || rCONTINUING)
	{
		//rCONTINUING = false; // Reset
		//pPpt->Out("\n\n");
		//pPpt->Out(".....................................................WARNING.....................................................\n");
		pPpt->Out("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	pPpt->Out("!\n");
	pPpt->Out("!\n");
/*	
	pPpt->Out("."); 
		pPpt->Out("\tTIMESTEP: "); 
		if(timestep<1000000) pPpt->Out(" ");
		if(timestep<100000) pPpt->Out(" ");
		if(timestep<10000) pPpt->Out(" ");
		if(timestep<1000) pPpt->Out(" ");
		if(timestep<100) pPpt->Out(" ");
		if(timestep<10) pPpt->Out(" ");
		//if(timestep<1) pPpt->Out(" ");
		pPpt->Out(timestep);
		pPpt->Out("\t\t");

		pPpt->Out("SIMULATION TIME: "); 
		if(time==0) pPpt->Out("0.000000");
		else pPpt->Out(time);
		pPpt->Out("s\t\t");

		pPpt->Out("ELAPSED TIME: "); 
		if(elapsedTime<1000) pPpt->Out(" ");
		if(elapsedTime<100) pPpt->Out(" ");
		if(elapsedTime<10) pPpt->Out(" ");
		//if(elapsedTime<1) pPpt->Out(" ");
		pPpt->Out(elapsedTime); 
		pPpt->Out("s\t\t");
		pPpt->Out(".\n");
		pPpt->Out(".");
*/
	pPpt->Out("!");
	pPpt->Out("\tWARNING: "); 
		pPpt->Out(strWarning); 
		//pPpt->Out("\t\t");
		//pPpt->Out(".");
		pPpt->Out("\n");
	pPpt->Out("!\n");
	pPpt->Out("!\n");
	pPpt->Out("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	cout << setprecision(temp_precision);
}

void Line(CProperties* pPpt, int length, int code)
{
	for(int dia=0; dia<length; ++dia) pPpt->Out(char(code));
	pPpt->Out("\n");
}
void Line(CProperties* pPpt, int length, char character)
{
	for(int dia=0; dia<length; ++dia) pPpt->Out(character);
	pPpt->Out("\n");
}

void Line(CProperties* pPpt, int length, char* character)
{
	for(int dia=0; dia<length; ++dia) pPpt->Out(character);
	pPpt->Out("\n");
}

/*
void ReadInput1(char *InputFile, double* &rParams)
{
	FILE *stream;
	float fp;
	char s[81];
	int q;

	stream = fopen(InputFile, "r");
	
	if(stream == NULL)
	{
		printf("Tools: ReadInput1: Error opening input file...\n");
		exit(1);
	}
	else
	{
//		fprintf(stream, "%s %f", "a-string", 3.593);

		// Set pointer to beginning of file:
		fseek(stream, 0L, SEEK_SET);

		q=0;
		do
		{
			// Read data back from file:
			fscanf(stream, "%s", s );
			fscanf(stream, "%f", &fp );
	
			cout << s << "\t = " << fp << endl;
			rParams[q] = fp;
			++q;
		}while(fscanf(stream, "%l")!=EOF);

		fclose(stream);
	}
}
*/

void GetHeader(char* &rheader, char* label, int n)
// Used to print column headers in output files
{
	//char header[10];
//	strcpy(rheader, "Node_");
	strcpy(rheader, label);
	if(n>=0 && n<=9)
	{ 
		rheader[5] = n + 48;
		rheader[6] = '\0';
	}
	else
	{
		if(n>=10 && n<=99)
		{
			rheader[5] = int(n/10) + 48;
			rheader[6] = n - 10*(int(n/10)) + 48;
			rheader[7] = '\0';
		}
		else
		{
			if(n>=100 && n<=999)
			{
				rheader[5] = int(n/100) + 48;
				rheader[6] = int((n - 100*(int(n/100)))/10) + 48;
				rheader[7] = n - 10*int((n - 100*(int(n/100)))/10)
								- 100*(int(n/100)) + 48;
				rheader[8] = '\0';
			}
			else cout << "More than 999 nodes per pipe!\n"; 
		}
	}
//	return header;
}

void CommonReadInput(char *InputFile, int num_parameters, char** &rlabels, double* &rvalues, char** &rstrings, int &rlast_entry)
{
	FILE *stream;
	float fp;
	int max = 500;
	char *s;
	s = new char [max];
	int q;
	int r, last_entry;

	stream = fopen(InputFile, "r");
	//cout << "Tools: CommonReadInput: input file: " << InputFile << "\n";	
	if(stream == NULL)
	{
		cout << "Tools: CommonReadInput: Error opening input file: " << InputFile << "\n";
		exit(1);
	}
	else
	{
		fseek(stream, 0L, SEEK_SET);

		q=0;

		rlabels = new char* [num_parameters];
		rvalues = new double [num_parameters];
		rstrings = new char* [num_parameters];
		for(r=0; r<num_parameters; ++r)
		{
			rlabels[r] = "empty";
			rvalues[r] = 0;
			rstrings[r] = "empty";
		}

	//	for(r=0; r<10; ++r) unord[r] = new char*[2];
	
		r=0;
		do
		{
			fpos_t start_of_string;
			fpos_t temp_pos;
			fpos_t end_of_string;

			fgetpos(stream, &start_of_string);
			//printf( "Before reading string, pos =  %ld \n", start_of_string);
			fgets(s, max, stream);	// Read a line of the file at a time

			int strptr = 0;
			int number_start = 0;
			int number_pos = 0;
			int word_start = 0;
			int word_pos = 0;

			bool found_digit = false;
			bool found_number = false;
			char* the_number;

			bool found_letter = false;
			bool found_word = false;
			char* the_word;

			bool not_a_number = false;
			bool comment = false;
			bool found_comment = false;

			bool space_or_tab_or_newline = true;
			bool exitt=false;

			while(s[strptr] != NULL && !found_comment)
			{
				if(( (int(s[strptr]) >=48 && int(s[strptr])<58) 
					|| s[strptr]=='-' || s[strptr]=='.' 
					|| (found_digit && (s[strptr]=='E' || s[strptr]=='e')) )
					&& !not_a_number && space_or_tab_or_newline)
					// A number 0-9 or '-' or decimal point or 'E' or 'e'
				{
					//cout << "a number" << endl;

					if(found_letter)		// Come to end of a word
					{
						found_letter = false;
						// Word runs from s[word_start] to s[word_start+word_pos]
						char *temp;
						temp = new char[word_pos+1];
						int ch;
						for(ch=0; ch<(word_pos+1); ++ch)
						{
							temp[ch] = s[word_start+ch];
						}
						temp[ch] = '\0';	// Terminates the string
//cout << "Found word: " << temp << endl;

						// Now if we have already found a word in the string prior to this
						// word, assume: word(label) = word(string), else discard
						if(found_word)
						{
							char* the_word_string;
							the_word_string = temp;

//cout << "Found parameter: " << the_word << " = " << the_word_string << endl;

							// Record in array
							rlabels[r] = the_word;
							rstrings[r] = the_word_string;
							last_entry = r;
							++r;
						
							found_word = false;		// Reset
						}
						else
						{
							found_word = true;
							the_word = temp;
						}
					}


					if(!found_digit)
					{
						if(s[strptr] == 'E' || s[strptr] == 'e') // Cannot start with E or e
						{
			//				cout << "E/e" << endl;
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
				}
				else // Not a number
				{
					if(found_digit)		// Come to end of number
					{
						found_digit = false;
						// number runs from s[number_start] to s[number_start+number_pos]
						char *temp;
						temp = new char[number_pos+1];
						int ch;
						for(ch=0; ch<(number_pos+1); ++ch)
						{
							temp[ch] = s[number_start+ch];
						}
						temp[ch] = '\0';	// Terminates the string
//cout << "Found number: " << temp << endl;

					//	found_number = true;
						the_number = temp;

						// Now if we have found a word in the string prior to this
						// number, assume: word = number, else discard
						if(found_word)
						{
//cout << "Found parameter: " << the_word << " = " << the_number << endl;

							fgetpos(stream, &end_of_string);
						//	printf("After reading string, end_of_string =  %ld \n", end_of_string);

							fsetpos(stream, &start_of_string);	// Rewind to before the current string
						//	fgetpos(stream, &temp_pos);
						//	printf("After rewinding, pos = %ld \n", temp_pos);

							// Now move pointer to before the relevant number
							temp_pos = start_of_string + number_start;
							fsetpos(stream, &temp_pos);
						//	cout << "positions moved = " << number_start << endl;
						//	fgetpos(stream, &temp_pos);
						//	printf("After moving, pos = %ld \n", temp_pos);

							// Read number into float
							fscanf(stream, "%f", &fp);
			//				cout << "fp = " << fp << endl;

						//	fgetpos(stream, &temp_pos);
						//	printf("After fscanf ing, pos = %ld \n", temp_pos);

							// Must reset the pointer to the end of this string
							fsetpos(stream, &end_of_string);
						//	fgetpos(stream, &temp_pos);
						//	printf("After resetting to end of string, pos = %ld \n", temp_pos);
						
							// Record in array
							rlabels[r] = the_word;
							rvalues[r] = fp;
							last_entry = r;
							++r;
							
							found_word = false;		// Reset
						}

					}

					if(s[strptr] == '/' && comment)
					{
						// Start of comment; ignore and go to next string line
						comment = false; // Reset
						found_comment = true; // This will exit the loop
					}
					else
					{
						if(s[strptr] == '/') comment=true;
						else
						{
							//if(s[strptr] != ' ' && s[strptr] != '\t')
							if(s[strptr] != ' ' && s[strptr] != '\t' && s[strptr] != '\n')
							{
								space_or_tab_or_newline = false;
								//cout << "a letter" << endl;
								if(!found_letter)
								{
									found_letter = true;
									word_start = strptr;
									word_pos = 0;
								}
								else ++word_pos;	// Continuation of a word
							}
							else
							{
							//	if(s[strptr] == ' ') cout << "a space" << endl;
							//	else
							//		if(s[strptr] == '\t') cout << "a tab" << endl;

								space_or_tab_or_newline = true;

								if(found_letter)	
								// Come to end of a word i.e. a space or a tab or a new line
								{
									found_letter = false;
									if(s[word_start] != '=') // Ignore =
									{
										// Word runs from s[word_start] to s[word_start+word_pos]
										char *temp;
										temp = new char[word_pos+1];
										int ch;
										for(ch=0; ch<(word_pos+1); ++ch)
										{
											temp[ch] = s[word_start+ch];
										}
										temp[ch] = '\0';	// Terminates the string
//cout << "Found word: " << temp << endl;

										// Now if we have already found a word in the string prior to this
										// word, assume: word(label) = word(string), else discard
										if(found_word)
										{
											char* the_word_string;
											the_word_string = temp;

//cout << "Found parameter: " << the_word << " = " << the_word_string << endl;

											// Record in array
											rlabels[r] = the_word;
											rstrings[r] = the_word_string;
											last_entry = r;
											++r;
						
											found_word = false;		// Reset
										}
										else
										{
											found_word = true;
											the_word = temp;
										}
									}
								}
							}
						}
					}
				}
//if(found_comment) cout << "Remainder of line was a comment" << endl;
				++strptr;
			};


			// If there was a number right before the NULL, must do this:
			if(found_digit)		// Come to end of number
			{
				found_digit = false;
				// number runs from s[number_start] to s[number_start+number_pos]
				char *temp;
				temp = new char [number_pos+1];
				int ch;
				for(ch=0; ch<(number_pos+1); ++ch)
				{
					temp[ch] = s[number_start+ch];
				}
				temp[ch] = '\0';	// Terminates the string

				//cout << endl;
				//cout << "number_pos+1 = " << number_pos+1 << endl;
	//			cout << "Right before NULL:" << endl;
	//			cout << "Found number: " << temp << endl;

			//	found_number = true;
				the_number = temp;

				// Now if we have found a word in the string prior to this
				// number, assume: word = number, else discard
				if(found_word)
				{
//cout << "Found parameter: " << the_word << " = " << the_number << endl;

					fgetpos(stream, &end_of_string);
				//	printf("After reading string, end_of_string =  %ld \n", end_of_string);

					fsetpos(stream, &start_of_string);	// Rewind to before the current string
				//	fgetpos(stream, &temp_pos);
				//	printf("After rewinding, pos = %ld \n", temp_pos);

					// Now move pointer to before the relevant number
					temp_pos = start_of_string + number_start;
					fsetpos(stream, &temp_pos);
				//	cout << "positions moved = " << number_start << endl;
				//	fgetpos(stream, &temp_pos);
				//	printf("After moving, pos = %ld \n", temp_pos);

					// Read number into float
					fscanf(stream, "%f", &fp);
	//				cout << "fp = " << fp << endl;

				//	fgetpos(stream, &temp_pos);
				//	printf("After fscanf ing, pos = %ld \n", temp_pos);

					// Must reset the pointer to the end of this string
					fsetpos(stream, &end_of_string);
				//	fgetpos(stream, &temp_pos);
				//	printf("After resetting to end of string, pos = %ld \n", temp_pos);
						
					// Record in array
					rlabels[r] = the_word;
					rvalues[r] = fp;
					last_entry = r;
					++r;
						
					found_word = false;		// Reset
			//		found_number = false;
				}
				//	else
				//		found_number = false;	// Discard
			}

			// If there was a word right before the NULL, must do this:
			if(found_letter)		// Come to end of number
			{
				found_letter = false;
				if(s[word_start] != '=') // Ignore =
				{
					// Word runs from s[word_start] to s[word_start+word_pos]
					char *temp;
					temp = new char[word_pos+1];
					int ch;
					for(ch=0; ch<(word_pos+1); ++ch)
					{
						temp[ch] = s[word_start+ch];
					}
					temp[ch] = '\0';	// Terminates the string
					//cout << "word_pos+1 = " << word_pos+1 << endl;
//cout << "Right before NULL:" << endl;
//cout << "Found word: " << temp << "END" << endl;

					// Now if we have already found a word in the string prior to this
					// word, assume: word(label) = word(string), else discard
					if(found_word)
					{
						char* the_word_string;
						the_word_string = temp;

//cout << "Found parameter: " << the_word << " = " << the_word_string << "END" << endl;

						// Record in array
						rlabels[r] = the_word;
						rstrings[r] = the_word_string;
						last_entry = r;
						++r;
						
						found_word = false;		// Reset
					}
					else
					{
						found_word = true;
						the_word = temp;
					}
				}
			}
			++q;
		}while(fscanf(stream, "%l")!=EOF);

		fclose(stream);
/*
		for(r=0; r<last_entry+1; ++r)
		{ 
			cout << rlabels[r] << "\t" << rvalues[r] << "\t" << rstrings[r] << endl;
		}
//*/
	}
	rlast_entry = last_entry;

	// Free dynamic memory
	delete [] s;
}

char* ConstructString(CProperties* pPpt, std::string dir_str, std::string name_str, bool ex, int id)
// For input parameter files of boundaries with normal EX or IN, e.g. nozzle
{
	std::string s;
	s = dir_str;
	if(ex) s += "EX_"; else s += "IN_";
	s += name_str; 

	if(id>=10) s += int(id/10) + 48;
	s += (id - int(id/10)*10) + 48; 
	s += ".txt";

	char *sz;
	sz = new char[s.length() + 1];
	strcpy(sz, s.c_str());

	return sz;
}
char* ConstructString(CProperties* pPpt, std::string dir_str, std::string name_str, int id)
// For input parameter files of boundaries with no EX or IN, e.g. cylinder, or results files
{
	std::string s;
	s = dir_str;
	s += name_str; 
	
	if(id>=10) s += int(id/10) + 48;
	s += (id - int(id/10)*10) + 48; 
	s += ".txt";
	
	char *sz;
	sz = new char[s.length() + 1];
	strcpy(sz, s.c_str());

	return sz;
}
char* ConstructString(CProperties* pPpt, std::string dir_str, std::string name_str, int id, std::string ext_str)
// For input parameter files of boundaries with no EX or IN, e.g. cylinder, or results files, specifying file extension
{
	std::string s;
	s = dir_str;
	s += name_str; 
	
	if(id>=10) s += int(id/10) + 48;
	s += (id - int(id/10)*10) + 48; 
	s += ext_str;
	
	char *sz;
	sz = new char[s.length() + 1];
	strcpy(sz, s.c_str());

	return sz;
}
char* ConstructString(CProperties* pPpt, std::string dir_str, std::string name_str, bool DOT_TEXT)
// For input parameter files of boundaries with no EX or IN, and ID included in name_str, e.g. engine
{
	std::string s;
	s = dir_str;
	s += name_str;
	if(DOT_TEXT) s += ".txt"; // Add .txt if desired
	char *sz;

	sz = new char[s.length() + 1];
	strcpy(sz, s.c_str());

	return sz;
}
char* ConstructString(CProperties* pPpt, std::string dir_str, std::string name_str)
// For input parameter files of boundaries with no EX or IN, and ID included in name_str, e.g. engine
{
	std::string s;
	s = dir_str;
	s += name_str;
	char *sz;

	sz = new char[s.length() + 1];
	strcpy(sz, s.c_str());

	return sz;
}
char* GetBoundaryName(int NAME)
{
	char* temp;
	switch(NAME)
	{
	case UNKNOWN:
		temp = "UNKNOWN";
		break;
	case INTERIOR:
		temp =  "Interior (INTERIOR)";
		break;
	case RESERVOIR:
		temp =  "Reservoir (RESERVOIR)";
		break;
	case JUNCTION:
		//temp =  "JUNCTION";
		temp =  "Junction (JUNCTION)";
		break;
	case ASSEMBLY:
		temp =  "Assembly (ASSEMBLY)";
		break;
	case ANECHOIC:
		temp = "Anechoic (ANECHOIC)";
		break;
	case END_ENVIRONMENT:
		temp = "End Environment (END_ENVIRONMENT)";
		break;
	case OPEN:
		temp =  "Open (OPEN)";
		break;
	case CLOSED:
		temp =  "Closed (CLOSED)";
		break;
	case NOZZLE:
		temp =  "Nozzle (NOZZLE)";
		break;
	case INFLOWEND:
		temp =  "Inflow (INFLOWEND)";
		break;
	case TRANSMISSIVE:
		//temp =  "TRANSMISSIVE";
		temp =  "Transmissive (TRANSMISSIVE)";
		break;
	case CYLINDER:
		temp =  "Cylinder (CYLINDER)";
		break;
	//case EXHVALVE:
	//	temp =  "Exhaust Valve (EXHVALVE)";
	//	break;
	//case INTVALVE:
	//	temp =  "Intake Valve (INTVALVE)";
	//	break;
	case VALVE:
		temp = "Valve (VALVE)";
		break;
	case SUDDEN_LEFT:
		temp =  "Sudden Expansion/Contraction, left side (SUDDEN_LEFT)";
		break;
	case SUDDEN_RIGHT:
		temp =  "Sudden Expansion/Contraction, right side (SUDDEN_RIGHT)";
		break;
	case APLDEV_LEFT:
		//temp =  "APL Device, left side (APLDEV_LEFT)";
		temp =  "APL Device";
		break;
	case APLDEV_RIGHT:
		temp =  "APL Device, right side (APLDEV_RIGHT)";
		break;
	case TURB_INLET_OR_OUTLET:
		temp =  "Turbine inlet or outlet (TURB_INLET_OR_OUTLET)";
		break;
	case TURBINLET:
		//temp =  "Turbine";
		temp =  "Turbine inlet (TURBINLET)";
		break;
	case TURBOUTLET:
		temp =  "Turbine outlet (TURBOUTLET)";
		break;
	case COMPINLET:
		temp =  "Compressor inlet (COMPINLET)";
		break;
	case COMPOUTLET:
		temp =  "Compressor outlet (COMPOUTLET)";
		break;
	case VOLUMEVALVE:
		temp =  "Volume valve (VOLUMEVALVE)";
		break;
	default:
		temp =  "UNDEFINED";
		break;
	}
	return temp;
}

char* StringToChar(std::string sz)
{
	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

std::string IntToString(int num)
{
	using namespace std;
	ostringstream myStream; //creates an ostringstream object
	myStream << num << flush;

	/*
	 * outputs the number into the string stream and then flushes
	 * the buffer (makes sure the output is put into the stream)
	 */
	return(myStream.str()); //returns the string form of the stringstream object
}

std::string CharToString(char* charstar)
{
	using namespace std;
	ostringstream myStream; //creates an ostringstream object
	myStream << charstar << flush;

  /*
   * outputs the number into the string stream and then flushes
   * the buffer (makes sure the output is put into the stream)
   */
	return(myStream.str()); //returns the string form of the stringstream object
}

char* Line(char* title_str, char* line_char)
{
	bool OVERLINE = false;
	int i;
	std::string s;
	s = title_str;
	
	std::string sz;

	if(OVERLINE)
	{
		for(i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
		sz += "\n"; // Add a newline character
	}
//	sz += s;
//	sz += "\n"; // Add a newline character
	for(i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += "\n"; // Add a newline character

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(char* title_str, char* line_char)
{
	bool OVERLINE = false;
	int i;
	std::string s;
	s = title_str;
	
	std::string sz;

	if(OVERLINE)
	{
		for(i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
		sz += "\n"; // Add a newline character
	}
	sz += s;
	sz += "\n"; // Add a newline character
	for(i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += "\n"; // Add a newline character

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(char* title_str, char* line_char, bool OVERLINE)
{
	int i;
	std::string s;
	s = title_str;
	
	std::string sz;

	if(OVERLINE)
	{
		for(i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
		sz += "\n"; // Add a newline character
	}
	sz += s;
	sz += "\n"; // Add a newline character
	for(i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += "\n"; // Add a newline character

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(char* title_str, char* line_char, bool OVERLINE, int line_length, int indent_tabs)
{
	int i;
	std::string s;
	s = title_str;

	std::string sz;

	if (OVERLINE)
	{
		for (i = 0; i < line_length; ++i) sz += line_char; // Add a line character for the length specified
		sz += "\n"; // Add a newline character
	}
	for (i = 0; i < indent_tabs; ++i) sz += "\t"; // Add tabs for the indent
	sz += s;
	sz += "\n"; // Add a newline character
	for (i = 0; i < line_length; ++i) sz += line_char; // Add a line character for the length specified
	sz += "\n"; // Add a newline character

	char* szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(char* title_str, char* line_char, bool OVERLINE, char* prefix)
{
	int i;
	std::string s;
	s = title_str;
	
	std::string sz;
	sz = prefix;

	// Top
	if(OVERLINE)
	{
		for(i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
		sz += "\n"; // Add a newline character
		sz += prefix;
	}
	
	// Middle
	sz += s;
	sz += "\n"; // Add a newline character
	
	// Bottom
	sz += prefix;
	for(i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += "\n"; // Add a newline character

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(char* title_str, char* line_char, int id)
{
	std::string s;
	s = title_str;

	s += " [";
	//s += "[";
	if(id>=10) s += int(id/10) + 48;
	s += (id - int(id/10)*10) + 48; 
	s += "]";
	
	std::string sz;
	sz = s;

	sz += "\n"; // Add a newline character
	for(int i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += "\n"; // Add a newline character

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(char* title_str, char* line_char, char* prefix)
{
	std::string s;
	s = title_str;
	
	std::string sz;
	sz = prefix;
	sz += s;
	sz += "\n"; // Add a newline character
	sz += prefix; // Add the prefix again
	for(int i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += "\n"; // Add a newline character

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(char* title_str, char* line_char, char* prefix, char* postfix)
{
	int i;
	std::string s_title;
	std::string s_postfix;
	s_title = title_str;
	s_postfix = postfix;
	for(i=0; i<int(s_postfix.length()); ++i) if(s_postfix[i] == '_') s_postfix[i] = ' ';
	
	std::string sz;
	sz = prefix;
	sz += s_title;
	sz += " - ";
	sz += s_postfix;
	sz += "\n";

	// Count the number of underline characters required
	sz += prefix; // Add the prefix again
	for(i=0; i<int(s_title.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += line_char;
	sz += line_char;
	sz += line_char;
	for(i=0; i<int(s_postfix.length()); ++i) sz += line_char; // Add a line character for each of the postfix characters
	sz += "\n"; // Add a newline character

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(char* title_str, char* line_char, int id, char* prefix)
{
	std::string s;
	s = title_str;

	s += " [";
	//s += "[";
	if(id>=10) s += int(id/10) + 48;
	s += (id - int(id/10)*10) + 48; 
	s += "]";
	
	std::string sz;
	sz = prefix + s;

	sz += "\n"; // Add a newline character
	sz += prefix; // Add the prefix again
	for(int i=0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += "\n"; // Add a newline character

	char *szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(string title_str, char* line_char)
{
	bool OVERLINE = false;
	int i;
	std::string s;
	s = title_str;

	std::string sz;

	if (OVERLINE)
	{
		for (i = 0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
		sz += "\n"; // Add a newline character
	}
	sz += s;
	sz += "\n"; // Add a newline character
	for (i = 0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += "\n"; // Add a newline character

	char* szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

char* Underline(string title_str, char* line_char, bool OVERLINE)
{
	int i;
	std::string s;
	s = title_str;

	std::string sz;

	if (OVERLINE)
	{
		for (i = 0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
		sz += "\n"; // Add a newline character
	}
	sz += s;
	sz += "\n"; // Add a newline character
	for (i = 0; i<int(s.length()); ++i) sz += line_char; // Add a line character for each of the title characters
	sz += "\n"; // Add a newline character

	char* szz;
	szz = new char[sz.length() + 1];
	strcpy(szz, sz.c_str()); // Put the title in the string
	return szz;
}

double AminusB(CProperties* pPpt, double A, double B)
{
	if (fabs(A - B) < pPpt->ZERO_TOL) return 0;
	else return A - B;
}

void Interpolate(double xLower, double yLower, double xUpper, double yUpper, double xInterp, double& rYInterp)
{
	rYInterp = yLower +
				(((yUpper - yLower)/(xUpper - xLower))	// Gradient
				*(xInterp - xLower));					// Difference
	;
}

double TotalPressurePa(CProperties* pPpt, double ps_Pa /*Pa*/, double Ts, double u)
// ====================================================================================================
// Returns total pressure in Pa using either 
//(isentropic) or incompressible relations
// ====================================================================================================
{
	if(pPpt->COMPRESSIBLE){
		double k = pPpt->gammaAir(Ts);
		return ps_Pa*pow(1 + ((k - 1)/2)*pow(u/sqrt(k*pPpt->R_air*Ts)/*Mach no.*/, 2), k/(k - 1));
	}
	else return ps_Pa + 0.5*(ps_Pa/(pPpt->R_air*Ts))/*rho*/*pow(u,2); // Incompressible flow
}

double TotalPressureBar(CProperties* pPpt, double ps_bar/*bar*/, double Ts, double u)
// ====================================================================================================
// Returns total pressure in bar using either compressible (isentropic) or incompressible relations
// ====================================================================================================
{
	if(pPpt->COMPRESSIBLE){
		double k = pPpt->gammaAir(Ts);
		return ps_bar*pow(1 + ((k - 1)/2)*pow(u/sqrt(k*pPpt->R_air*Ts)/*Mach no.*/, 2), k/(k - 1));
	}
	else return ps_bar + ( (0.5*( (ps_bar*1e5)/(pPpt->R_air*Ts))/*rho*/*pow(u,2)) /1e5); 
	// Incompressible flow (must convert second term to and from Pa due to addition)
}

double TotalTemperature(CProperties* pPpt, double Ts, double u)
// ====================================================================================================
// Returns total temperature in K
// ====================================================================================================
{
	return Ts*(1 + ((pPpt->gammaAir(Ts)-1)/2)*pow(u/(sqrt(pPpt->gammaAir(Ts)*pPpt->R_air*Ts))/*Mach no.*/,2));
}

double StaticTemperature(CProperties* pPpt, double T0, double u)
// ====================================================================================================
// Returns static temperature in K
// ====================================================================================================
{
	double Ts, Ts_old, tol;
	Ts = T0; // Initial estimate
	tol = 1e-6;
//pPpt->Out("\n"); 
//pPpt->Out("u = "); pPpt->Out(u); pPpt->Out("\n");
//pPpt->Out("T0 = "); pPpt->Out(T0); pPpt->Out("\n");
	do{
		Ts_old = Ts;
		Ts = T0/(1 + ((pPpt->gammaAir(Ts)-1)/2)*pow(u/(sqrt(pPpt->gammaAir(Ts)*pPpt->R_air*Ts)),2));
//pPpt->Out("Ts = "); pPpt->Out(Ts); pPpt->Out("\n"); 
	}while(fabs(Ts - Ts_old)>tol);
//pPpt->Out("\n");
	return Ts;
}

void PrecisionInformation()
{
	cout << "epsilon(float) = " << numeric_limits<float>::epsilon( ) << endl; //output: 2.22045e-016
	cout << "epsilon(double) = "<< numeric_limits<double>::epsilon( ) << endl; //output: 2.22045e-016
	//if ( ((d1*10)-(20.0/d2)) <= numeric_limits<double>::epsilon())
	//{
	//  do_equal();
	//}
	cout << "epsilon(long double) = " << numeric_limits<long double>::epsilon( ) << endl; //output: 2.22045e-016
	cout << endl;

	cout << "digits10(float) = " << numeric_limits<float>::digits10 << endl;//output: 15
	cout << "digits10(double) = " << numeric_limits<double>::digits10 << endl;//output: 15
	cout << "digits10(long double) = " << numeric_limits<long double>::digits10 << endl; // 18
	cout << endl;

	double d=0;
	if(numeric_limits<double>::has_quiet_NaN)
	 d=numeric_limits<double>::quiet_NaN();
	else if (numeric_limits<double>::has_signaling_NaN)
	 d=numeric_limits<double>::signaling_NaN();
	else cerr<<"NaN for double isn't supported";

	float f=0;
	if(numeric_limits<float>::has_infinity)
	 f=numeric_limits<float>::infinity();
	else cerr<<"infinity for float isn't supported";

	double one = 1.1234567890123400000;
	double two = 1.1234567890123500000;

	int temp_precision = cout.precision(); cout << setprecision(25);
	cout << one << endl;
	cout << two << endl;
	cout << one - two << endl;

	cout << setprecision(temp_precision);
	exit(1);
}
char* InflowOrOutflow(int flowDir){
    switch(flowDir){
case INFLOW:
    return "INFLOW";
    break;
case OUTFLOW:
    return "OUTFLOW";
    break;
case NOFLOW:
    return "NOFLOW";
    break;
default:
    return "UNKNOWN";
    break;
    }
}
int InflowOrOutflowNumber(int flowDir){
    switch(flowDir){
case INFLOW:
    return 1;
    break;
case OUTFLOW:
    return -1;
    break;
case NOFLOW:
    return 0;
    break;
default:
    return 999;
    break;
    }
}