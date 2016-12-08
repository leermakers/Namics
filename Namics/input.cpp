class Input {
public:
	Input(string); 

~Input();

	string name; 
	ifstream in_file;
	
	std::string In_buffer;
	string string_value;
	bool Input_error; 

	
	std::vector<string> KEYS;
	std::vector<string> SysList;
	std::vector<string> MolList;
	std::vector<string> MonList;
	std::vector<string> LatList;
	std::vector<string> NewtonList;
	std::vector<string> OutputList; 
	std::vector<string> EngineList;
	std::vector<std::string>  elems;

	 
	void PrintList(std::vector<std::string> &);
	std::vector<std::string>& split(const std::string &, char, std::vector<std::string>&); 
	bool IsDigit(string &);
	bool Get_int(string, int &);
	bool Get_int(string, int &,const std::string &);
	bool Get_int(string, int &, int, int, const std::string &);
	bool Get_string(string, string &,  const std::string);
	bool Get_string(string, string &, std::vector<std::string>&, const std::string);
	bool Get_double(string, double &);
	bool Get_double(string, double &, const std::string &);
	bool Get_double(string, double &, double, double, const std::string &);
	bool Get_bool(string, bool &,const std::string &);
	bool TestNum(std::vector<std::string> &, string ,int, int, int);
	int GetNumStarts(void); 
	bool CheckParameters(string, string name,std::vector<std::string> &, std::vector<std::string> &,std::vector<std::string> &);
	bool LoadItems(string, std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string> &);
	bool CheckInput(void); 

};
Input::Input(string name_) {
	name=name_;
	KEYS.push_back("start"); KEYS.push_back("sys");KEYS.push_back("mol"); KEYS.push_back("mon"); 	
 	KEYS.push_back("lat"); KEYS.push_back("newton"); KEYS.push_back("engine"); KEYS.push_back("output");
	in_file.open(name.c_str()); Input_error=false; 
	if (in_file.is_open()) {
		int line_nr=0;
		bool add=true; 
		std:: string In_line;
		while (in_file) {
			line_nr++; 
			std::getline(in_file,In_line); 
			In_line.erase(std::remove(In_line.begin(), In_line.end(), ' '), In_line.end());
			if (In_line.length()==0) add = false;
			if (In_line.length()>2) {if (In_line.substr(0,2) == "//") {add = false;}}
			if (add) {In_buffer.append(SSTR(line_nr)).append(":").append(In_line).append("#");} else {add = true;} 
		}
		in_file.close();
		split(In_buffer,'#',elems);
		if (!CheckInput()) Input_error=true;
	} else {cout <<  "Inputfile " << name << " is not found. " << endl; Input_error=true; }

}
Input::~Input() {
}

void Input::PrintList(std::vector<std::string> &LIST) {
	int length=LIST.size(); 
	int i=0;
	while (i<length) {cout << LIST[i] << " ; "; i++; } 
}

std::vector<std::string>& Input::split(const std::string &s, char delim, std::vector<std::string>&elems){
	bool add=true; 
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss,item,delim)) {
		item.erase(std::remove(item.begin(), item.end(), ' '), item.end());
		std::size_t pos = item.find("//");
		if (add) {elems.push_back(item.substr(0,pos));}
		add=true;
	}
	return elems;
}

bool Input:: IsDigit(string &s) {
	return (s=="0" || s=="1" ||s=="2" || s=="3" || s=="4" || s=="5" ||s=="6" || s=="7" || s=="8" || s=="9"); 
}

bool Input:: Get_int(string s, int &ss) {
	bool success = false;
	stringstream string_s;
        string_s << s;
	string_s >> ss; 	
	string sub; if (s.length()>0) sub = s.substr(0,1);
	string sub2; if (s.length()>1) sub2= s.substr(1,1); 
	if (sub=="-") { if (s.length()>1) {if (IsDigit(sub2)) {success=true;} }}
 	else {if (s.length()>0) {if (IsDigit(sub)) { success = true;}}} 
	return success;
}

bool Input:: Get_int(string s, int &ss,const std::string &error) {
	bool success = false;
	stringstream string_s;
        string_s << s;
	string_s >> ss; 	
	string sub; if (s.length()>0) sub = s.substr(0,1);
	string sub2; if (s.length()>1) sub2= s.substr(1,1); 
	if (sub=="-") { if (s.length()>1) {if (IsDigit(sub2)) {success=true;} }}
 	else {if (s.length()>0) {if (IsDigit(sub)) { success = true;}}} 
	if (!success) cout << error << endl;
	return success;
}

bool Input:: Get_int(string s, int &ss, int low, int high, const std::string &error) {
	bool success = false;
	stringstream string_s;
        string_s << s;
	string_s >> ss; 	
	string sub; if (s.length()>0) sub = s.substr(0,1);
	string sub2; if (s.length()>1) sub2= s.substr(1,1); 
	if (sub=="-") { if (s.length()>1) {if (IsDigit(sub2)) {success=true;} }}
 	else {if (s.length()>0) {if (IsDigit(sub)) { success = true;}}} 
	if (!success) {cout << error << endl;
	} else {
		success=false;
		if (ss<low || ss>high) {cout << "Value out of range: " << error << endl ; } else success=true;
	}	
	return success;
}

bool Input:: Get_string(string s, string &ss,  const std::string error) {
	bool success = false;
	if (s.length() > 0) { ss=s; success = true;}
	if (!success) cout << error << endl;
	return success; 
}
bool Input:: Get_string(string s, string &ss, std::vector<std::string>&S, const std::string error) {
	bool success = false;
	if (s.length() > 0) { ss=s; success = true;}
	if (!success) {cout << error << endl;} 
	else {
		int i=0; 
		success = false;
		int length = S.size();
		if (i<length) {
			if (S[i]==ss) success=true;  
			i++;
		}
		if (!success) {cout << "Value of " << ss << " is not recognised. Select from: " << endl; 
			PrintList(S); cout << endl;
		} 
	}
	return success; 
}

bool Input:: Get_double(string s, double &ss) {
	bool success=false;
	stringstream string_s;
	string_s << s ;
	string_s >> ss;
	string sub; if (s.length()>0) sub = s.substr(0,1);
	string sub2; if (s.length()>1) sub2= s.substr(1,1); 
	if (sub=="-") { if (s.length()>1) {if (IsDigit(sub2)) {success=true;} }}
 	else {if (s.length()>0) {if (IsDigit(sub)) { success = true;}}} 
	return success;
}
bool Input:: Get_double(string s, double &ss, const std::string &error) {
	bool success=false;
	stringstream string_s;
	string_s << s ;
	string_s >> ss;
	string sub; if (s.length()>0) sub = s.substr(0,1);
	string sub2; if (s.length()>1) sub2= s.substr(1,1); 
	if (sub=="-") { if (s.length()>1) {if (IsDigit(sub2)) {success=true;} }}
 	else {if (s.length()>0) {if (IsDigit(sub)) { success = true;}}} 
	if (!success) cout << error << endl;
	return success;
}
bool Input:: Get_double(string s, double &ss, double low, double high, const std::string &error) {
	bool success=false;
	stringstream string_s;
	string_s << s ;
	string_s >> ss;
	string sub; if (s.length()>0) sub = s.substr(0,1);
	string sub2; if (s.length()>1) sub2= s.substr(1,1); 
	if (sub=="-") { if (s.length()>1) {if (IsDigit(sub2)) {success=true;} }}
 	else {if (s.length()>0) {if (IsDigit(sub)) { success = true;}}} 
	if (!success) {cout << error << endl;
	} else {
		success=false;
		if (ss<low || ss>high) {cout << "Value out of range: " << error << endl ; } else success=true;
	}	
	return success;
}

bool Input:: Get_bool(string s, bool &ss,const std::string &error) {
	bool success=false;
	std::string string_s;
	string_s=s;
	if (string_s =="true" ||  string_s =="True" || string_s =="TRUE") {ss=true; success=true;} else 
	if (string_s =="false" ||  string_s =="False" || string_s =="FALSE") ss=false; else success=false; 
	if (!success) cout << error << endl;
	return success;
}

bool Input:: TestNum(std::vector<std::string> &S, string c,int num_low, int num_high, int UptoStartNumber ) {
	bool InList=false;
	int i=0;
	int length = elems.size();
	int n_starts=0;
	while (i<length) {  
		std::vector<std::string> set;
	       split(elems[i],':',set);
		if (set[1]=="start") n_starts++; 
		if (c==set[1] && n_starts<UptoStartNumber ){
			InList=false;
			int S_length=S.size();
			for (int k=0; k<S_length; k++) if (set[2]==S[k]) InList=true;
			if (!InList) S.push_back(set[2]);
		} 
		i++; 
	}
	int number=S.size();
	if (number>num_low-1 && number<num_high+1) {return true;}
	return false;
}

int Input:: GetNumStarts() {
	int length = elems.size();
	int number=0;
	for (int i=0; i<length; i++) {
		vector<std::string> set;
		split(elems[i],':',set);
		if (set[1] == "start") number++;
	}
	return number;
}
bool Input:: CheckParameters(string keyword, string name,std::vector<std::string> &Standard, std::vector<std::string> &Input,std::vector<std::string> &Input_values) {
	bool success=true;
	bool prop_found; 
	int length = elems.size();
	int S_length = Standard.size();
	int I_length;
	string parameter; 
	int i=0;
	int j; 
	while (i<length){
		vector<std::string> set;
		split(elems[i],':',set);
		if (set[1] == keyword && set[2]== name) {
			parameter=set[3];
			j=0; prop_found=false; 
			while (j<S_length && !prop_found) {
				if (Standard[j]==parameter) prop_found=true;
				j++; 
			}
			if (!prop_found) {success=false; cout <<"In line " << set[0] << " "  << keyword << " property '" << parameter << "' is unknown. Select from: "<< endl; 
				for (int k=0; k<S_length; k++) cout << Standard[k] << endl; 
			} else {
				j=0; I_length = Input.size(); prop_found=false;
				while (j<I_length &&!prop_found) {
					if (Input[j]==parameter) prop_found = true; 
					j++;
				}
				if (prop_found) {success=false; cout <<"In line " << set[0] << " " << keyword << " property '" << parameter << "' is already defined. "<< endl; }
				else {Input.push_back(parameter); Input_values.push_back(set[4]);}
			}
		}
		i++;
	}
	return success; 
}


bool Input:: LoadItems(string template_,std::vector<std::string> &Out_name, std::vector<std::string> &Prop, std::vector<std::string> &Value) {
	int length = elems.size();
	int key_length = KEYS.size();
	int name_length;
	bool key_found,name_found; 
	int i=0; 
 	int j,k; 
	while (i<length){
		vector<std::string> set;
		split(elems[i],':',set);
		if (set[1] == template_) {
			j=0; key_found=false;
			while (j<key_length) {
				if (KEYS[j] == set[2]) {
					key_found=true;
					switch (j) {
						case 1:
							if (SysList[0]!=set[3]) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl; 
								PrintList(SysList); name_found=false;} 
							break;
						case 2:
							k=0; name_length=MolList.size(); name_found=false; 
							while (k<name_length) {
								if (MolList[k]==set[3]) name_found=true; 
								k++;
							}
							if (!name_found) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: " << endl;
								PrintList(MolList);
							}
							break;
						case 3:
							k=0; name_length=MonList.size(); name_found=false; 
							while (k<name_length) {
								if (MonList[k]==set[3]) name_found=true; 
								k++;
							}
							if (!name_found) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(MonList); 
							}
							break;
						case 4:
							if (LatList[0]!=set[3]) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(LatList); name_found=false;
							}
							break;
						case 5:
							if (NewtonList[0]!=set[3]) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(NewtonList); name_found=false;
							}
							break;
						case 6:
							if (EngineList[0]!=set[3]) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(EngineList); name_found=false;
							}
							break;
						default:
							key_found=false; 
						}					
					
				} 
				j++;
			} 
			if (!key_found) {cout<< "In line " << set[0] << " the keyword '" << set[2] << "' not recognized. Choose keywords from: " << endl; 
				for (int k=1; k<7; k++) cout << KEYS[k] << endl; return false; 
			} 
			if (!name_found) {cout << endl; return false; 
			}   
			Out_name.push_back(set[2]); Prop.push_back(set[3]); Value.push_back(set[4]);	
		} //end temp
		i++;
	} //end i; 
	return true; 
}

bool Input:: CheckInput(void) {
	bool success=true;
	bool Keywordfound; 
	int key_length=KEYS.size();
	int length = elems.size();  
	int j; 
	
	string word;
	
	if (length==0) {cout << "inputfile is empty " << endl; success=false; } 
	int i=0;
	while (i<length) { 
		vector<std::string> set;
		split(elems[i],':',set); 
		if (set.size() !=5) {if (set[1]!="start") {cout << " Line number " << set[0] << " does not contain 4 items" << endl; success=false; return false;}}
		i++;
	}
	i=0;
	while (success && i<length) {
		vector<std::string> set;
		split(elems[i],':',set); 
		if (set[1]=="output" && set[2]=="in") {cout << "In line " << set[0] << "file extension 'in' is not allowed. This would overwride the inputfile."<< endl; success=false;}
		if (set[1]=="output" && set[3]=="template") {
			word=set[4];
			j=0; Keywordfound=false; key_length =KEYS.size();
 
			while (j<key_length) { 
				if (word==KEYS[j]) {Keywordfound=true; }
				j++; 
			}
			if (!Keywordfound) {
				KEYS.push_back(word);
				key_length++;
			} else {cout << " In line " << set[0] << "template name '" << word << "' is already taken or is a standard Keyword." <<endl; success=false;}
		}
		i++;
	} 	 
	i=0;
	while (i<length) {	
		vector<std::string> set;
		split(elems[i],':',set); 
		word=set[1];
		//word.erase(std::remove(word.begin(), word.end(), ' '), word.end());
		j=0; Keywordfound=false;
		while (j<key_length) { 
			if (word==KEYS[j]) {Keywordfound=true; }
			j++;
		}
		if (!Keywordfound) {cout << word << " is not valid keyword in line " << set[0] << endl;
			cout << "select one of the following:" << endl; 
			for( int k=0; k<key_length; k++) cout << KEYS[k] << endl;  success=false;}
		i++; 
	}

	if (!TestNum(SysList,"sys",0,1,1)) {cout << "There can be no more than 1 sys name in the input" << endl; }
	if (!TestNum(LatList,"lat",1,1,1)) {cout << "There must be exactly one lat name in the input" << endl; success=false;}
	if (!TestNum(NewtonList,"newton",0,1,1)) {cout << "There can be no more than 1 newton name in input" << endl; success=false;}
	if (!TestNum(MonList,"mon",1,1000,1)) {cout << "There must be at least one mon name in input" << endl; success=false;}
	if (!TestNum(MolList,"mol",1,1000,1)) {cout << "There must be at least one mol name in input" << endl; success=false;}
	if (!TestNum(OutputList,"output",1,1000,1)) {cout << "No output defined! " << endl;}
	if (!TestNum(EngineList,"engine",1,1,1)) {cout << "You need to define exactly one engine in the input " << endl; success=false;}
	
	return success;
}

