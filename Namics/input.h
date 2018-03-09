#ifndef INPUTxH
#define INPUTxH
#include "namics.h"

class Input {
public:
	Input(string); 

~Input();

	string name; 
	ifstream in_file;

	
	std::string In_buffer;
	string string_value;
	bool Input_error; 
	string filename; 
	
	std::vector<string> KEYS;
	std::vector<string> SysList;
	std::vector<string> MolList;
	std::vector<string> MonList;
	std::vector<string> LatList;
	std::vector<string> AliasList; 
	std::vector<string> NewtonList;
	std::vector<string> OutputList; 
	std::vector<string> EngineList;
	std::vector<string> MesodynList;
	std::vector<string> VarList;
	std::vector<std::string> elems;
	 
	void PrintList(std::vector<std::string>);
	std::vector<std::string>& split( std::string , char, std::vector<std::string>&); 
	bool IsDigit(string &);
	int Get_int(string, int );
	bool Get_int(string, int &, const std::string &);
	bool Get_int(string, int &, int, int, const std::string &);
	string Get_string(string, const string &);
	bool Get_string(string, string &, const  std::string & );
	bool Get_string(string , string &, std::vector<std::string>&, const std::string &);
	Real Get_Real(string, Real );
	bool Get_Real(string, Real &, const std::string &);
	bool Get_Real(string, Real &, Real, Real, const std::string &);
	bool Get_bool(string, bool );
	bool Get_bool(string, bool &, const std::string &);
	bool TestNum(std::vector<std::string> &, string ,int, int, int);
	int GetNumStarts(void); 
	bool CheckParameters(string, string,int,std::vector<std::string> &, std::vector<std::string> &,std::vector<std::string> &);
	bool LoadItems(string, std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string> &);
	bool CheckInput(void); 
	bool InSet(std::vector<std::string> &, string);
	bool InSet(vector<int> &, int);
	bool InSet(vector<int> &, int&, int);
	bool InSet(std::vector<std::string> &, int&, string);
	bool ReadFile(string,string &); 
	bool ArePair(char ,char ); 
	bool EvenBrackets(string,vector<int> &, vector<int>&); 
	bool EvenSquareBrackets(string,vector<int> &, vector<int>&); 
	bool MakeLists(int);

};
#endif
