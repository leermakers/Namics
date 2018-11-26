#ifndef INPUTxH
#define INPUTxH
#include "namics.h"
#include "output_info.h"
#include <functional>

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
	OutputInfo output_info;

	std::vector<string> KEYS;
	std::vector<string> SysList;
	std::vector<string> MolList;
	std::vector<string> MonList;
	std::vector<string> LatList;
	std::vector<string> AliasList;
	std::vector<string> NewtonList;
	std::vector<string> OutputList;
	std::vector<string> MesodynList;
	std::vector<string> ClengList;
	std::vector<string> TengList;
	std::vector<string> VarList;
	std::vector<std::string> elems;
	std::vector<string> StateList;
	std::vector<string> ReactionList;


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

private:
	void parseOutputInfo();
};

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

class Reader {
public:
  Reader(string);
  ~Reader();

  enum error {
    ERROR_FILE_FORMAT,
    ERROR_EXTENSION,
  };

  enum file_format {
    NONE,
    VTK,
    PRO,
  };

  size_t dimensions;
  size_t component_no;
  size_t MX;
  size_t MY;
  size_t MZ;

  file_format filetype;
  vector<string> headers;

  void get_rho(vector<double>&);
  void get_rho(vector< vector<double> >&);

  template <typename T>
  inline T val(vector<T>& v, int x, int y, int z) {
      return v[x * JX + y * JY + z * JZ];
  }

  template <typename T>
  inline T* val_ptr(vector<T>& v, int x, int y, int z) {
    return &v[x * JX + y * JY + z * JZ];
  }

  vector<double> rho;
  vector< vector<double> > multicomponent_rho;

//	void skip_bounds(function<void(int, int, int)>);

  static vector<string> tokenize(string, char);

	vector<double> with_bounds(vector<double>);

private:
  size_t JX;
  size_t JY;
  size_t JZ;


  int init_rho_fromvtk(string);
  vector<string> init_rho_frompro(string);
	void adjust_pro_indexing();
};

#endif
