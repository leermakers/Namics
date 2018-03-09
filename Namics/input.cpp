#include "input.h"
Input::Input(string name_) {
	name=name_;
	KEYS.push_back("start"); 
	KEYS.push_back("sys");KEYS.push_back("mol"); KEYS.push_back("mon"); KEYS.push_back("alias"); 	
 	KEYS.push_back("lat"); KEYS.push_back("newton"); KEYS.push_back("engine"); KEYS.push_back("mesodyn"); KEYS.push_back("output"); 
	KEYS.push_back("var");    
	

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
			if (add) elems.push_back(SSTR(line_nr).append(":").append(In_line)); else add=true;
		}
		in_file.close();
		if (!CheckInput()) Input_error=true;
	} else {cout <<  "Inputfile " << name << " is not found. " << endl; Input_error=true; }

}
Input::~Input() {
}

bool Input::ArePair(char opening,char closing){
	if (opening == '(' && closing == ')') return true;
	else if (opening == '[' && closing == ']') return true;
	//else if (opening == '{' && closing == '}') return true;
	return false; 
}

bool Input::EvenSquareBrackets(string exp,vector<int> &open, vector<int> &close) {
	vector <char> S;
	int length = exp.size(); 
	for (int i=0; i<length; i++) {
		if (exp[i] == '[' ) {S.push_back(exp[i]); open.push_back(i);} 
		else if (exp[i] == ']') { close.push_back(i);
			if (S.size()==0 || !ArePair(S[S.size()-1],exp[i])) return false;
			else
			S.pop_back();
		}
	}
	return S.size()==0 ? true:false; 
}

bool Input::EvenBrackets(string exp,vector<int> &open, vector<int> &close) {
	vector <char> S;
	int length = exp.size(); 
	for (int i=0; i<length; i++) {
		if (exp[i] == '(' ) {S.push_back(exp[i]); open.push_back(i);} 
		else if (exp[i] == ')') { close.push_back(i);
			if (S.size()==0 || !ArePair(S[S.size()-1],exp[i])) return false;
			else
			S.pop_back();
		}
	}
	return S.size()==0 ? true:false; 
}

bool Input::ReadFile(string fname, string &In_buffer) {
	ifstream this_file;
	bool success=true; 
	bool add;  
	this_file.open(fname.c_str());
	std:: string In_line; 
	if (this_file.is_open()) {
		while (this_file) {
			add=true;
			std::getline(this_file,In_line);
			In_line.erase(std::remove(In_line.begin(), In_line.end(), ' '), In_line.end());
			if (In_line.length()==0) add = false;
			if (In_line.length()>2) {if (In_line.substr(0,2) == "//") {add = false;}}
			if (add) {In_buffer.append(In_line).append("#");}; 
		}	
		this_file.close();
		if (In_buffer.size()==0) {cout << "File " + fname + " is empty " << endl; success=false; }
	} else {cout <<  "Inputfile " << fname << " is not found. " << endl; success=false; }
	return success;
}


void Input::PrintList(std::vector<std::string> LIST) {
	int length=LIST.size(); 
	int i=0;
	while (i<length) {cout << LIST[i] << " ; "; i++; } 
}

std::vector<std::string>& Input::split(std::string s, char delim, std::vector<std::string>&elems){
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

int Input:: Get_int(string s, int ss) {
	bool success = false;
	int sss; 
	stringstream string_s;
        string_s << s;
	string_s >> sss; 	
	string sub; if (s.length()>0) sub = s.substr(0,1);
	string sub2; if (s.length()>1) sub2= s.substr(1,1); 
	if (sub=="-") { if (s.length()>1) {if (IsDigit(sub2)) {success=true;} }}
 	else {if (s.length()>0) {if (IsDigit(sub)) { success = true;}}} 
	if (!success) sss=ss;
	return  sss; 
}

bool Input:: Get_int(string s, int &ss, const std::string &error) {
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

string Input:: Get_string(string s, const string &ss) {
	if (s.length() > 0) { return s;} else {return ss;} 

}
bool Input:: Get_string(string s, string &ss, const  std::string &error) {
	bool success = false;
	if (s.length() > 0) { ss=s; success = true;}
	if (!success) cout << error << endl;
	return success; 
}
bool Input:: Get_string(string s, string &ss, std::vector<std::string>&S, const std::string &error) {
	bool success = false; 
	if (s.length() > 0) { ss=s; success = true;}
	if (!success) {cout << error << endl;} 
	else {	success=InSet(S,ss);  
		if (!success) {cout << error << " value '" << ss << "' is not allowed. Select from: " << endl; 
			PrintList(S); cout << endl;
		} 
	}
	return success; 
}

Real Input:: Get_Real(string s, Real ss) {
	bool success=false;
	Real sss;
	stringstream string_s;
	string_s << s ;
	string_s >> sss;
	string sub; if (s.length()>0) sub = s.substr(0,1);
	string sub2; if (s.length()>1) sub2= s.substr(1,1); 
	if (sub=="-") { if (s.length()>1) {if (IsDigit(sub2)) {success=true;} }}
 	else {if (s.length()>0) {if (IsDigit(sub)) { success = true;}}} 
	if (!success) sss=ss;
	return sss; 
}
bool Input:: Get_Real(string s, Real &ss, const std::string &error) {
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
bool Input:: Get_Real(string s, Real &ss, Real low, Real high, const  std::string &error) {
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

bool Input:: Get_bool(string s, bool ss) {
	bool success=true;
	bool sss; 
	if (s =="true" || s =="True" || s =="TRUE") {sss=true; success=true;} else 
	if (s =="false" ||  s =="False" || s =="FALSE") sss=false; else success=false; 
	if (!success) sss=ss;
	return sss; 
}
bool Input:: Get_bool(string s, bool &ss, const std::string &error) {
	bool success=false;
	if (s =="true" ||  s =="True" || s =="TRUE") {ss=true; success=true;} else 
	if (s =="false" ||  s =="False" || s =="FALSE") ss=false; else success=false; 
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
		if ( (i==length-1)&& ( set[1]!="start") ) number++;
	}
	return number;
}

bool Input:: InSet(std::vector<std::string> &Standard, string keyword){
	bool success=false; 
	int S_length = Standard.size();
	int i=0;
	while (i<S_length && !success) {
		if (Standard[i] == keyword) success=true; 
		i++;
	}
	return success; 
}
bool Input:: InSet(std::vector<std::string> &Standard, int &pos, string keyword){
	bool success=false;
	int S_length = Standard.size();
	int i=0;
	while (i<S_length && !success) {
		if (Standard[i] == keyword) {success=true; pos=i;} 
		i++;
	}
	return success; 
}
bool Input:: InSet(vector<int> &Standard, int keyword){
	bool success=false; 
	int S_length = Standard.size();
	int i=0;
	while (i<S_length && !success) {
		if (Standard[i] == keyword) success=true; 
		i++;
	}
	return success; 
}
bool Input:: InSet(vector<int> &Standard, int &pos, int keyword){
	bool success=false; 
	int S_length = Standard.size();
	int i=0;
	while (i<S_length && !success) {
		if (Standard[i] == keyword) {success=true; pos=i;} 
		i++;
	}
	return success; 
}

bool Input:: CheckParameters(string keyword, string name,int start, std::vector<std::string> &Standard, std::vector<std::string> &Input,std::vector<std::string> &Input_values) {
	bool success=true;
	bool prop_found; 
	int length = elems.size();
	int S_length = Standard.size();
	int I_length;
	string parameter; 
	int n_start=0;
	int n_found=0;
	int i=0;
	int j; 
	while (i<length && n_start<start){
		vector<std::string> set;
		split(elems[i],':',set);
		if (set[1]=="start") {
			n_start++; 
			int k=0;
			int k_length=Input.size();  
			while (k<k_length) { //remove doubles and keep last value
				int l=k+1;
				int l_length=Input.size();  
				while (l<l_length) {
					if (Input[k]==Input[l]) {
						Input_values[k]=Input_values[l];
						Input.erase(Input.begin()+l);
						Input_values.erase(Input_values.begin()+l);
						l_length--;
						k_length--;
						l--;
					}
					l++;	
				}
				k++; 	
			}	
		}
		else {
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
					j=0; I_length = Input.size(); prop_found=false; n_found=0;
					while (j<I_length) {
						if (Input[j]==parameter) {prop_found = true; n_found++;}
						j++;
					}
					if (prop_found && n_found>1 && n_start==0) {success=false; cout <<n_start<<" "  << start << endl;  cout <<"In line " << set[0] << " " << keyword << " property '" << parameter << "' is already defined. "<< endl; }
					else {
						if (prop_found && n_found>1) {success=false; cout <<"After 'start' " << n_start << ", in line " << set[0] << " " << keyword << " property '" << parameter << "' is already defined. "<< endl; }
						else {Input.push_back(parameter); Input_values.push_back(set[4]);}
					}
				}
			}
		}
		i++;
	}
	return success; 
}


bool Input:: LoadItems(string template_,std::vector<std::string> &Out_key, std::vector<std::string> &Out_name, std::vector<std::string> &Out_prop) {
if (debug) cout <<"LoadItems in Input " << endl; 
	bool success=true; 
	int length = elems.size();
	bool wild_monlist=false;
	bool wild_mollist=false;
	bool wild_aliaslist=false;
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
					switch (j-1) {
						case 0:
							if (set[3]=="*") set[3]=SysList[0];
							//if (set[1]=="var" || set[1]=="search") name_found=true;
							if (SysList[0]!=set[3] && !name_found) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl; 
								PrintList(SysList); name_found=false;} 
							break;
						case 1:
							name_found=false;
							if (set[3]=="*") { name_found = true; wild_mollist=true;}
							//if (set[1]=="var" || set[1]=="search") name_found=true;
							k=0; name_length=MolList.size();  
							while (k<name_length && !name_found) {
								if (MolList[k]==set[3]) name_found=true; 
								k++;
							}
							if (!name_found) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: " << endl;
								PrintList(MolList);
							}
							break;
						case 2:
							name_found=false;
							if (set[3]=="*") {name_found=true; wild_monlist=true;} 
							//if (set[1]=="var") name_found=true;
							k=0; name_length=MonList.size();  
							while (k<name_length && !name_found) {
								if (MonList[k]==set[3]) name_found=true; 
								k++;
							}
							if (!name_found) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(MonList); 
							}
							break;
						case 3:
							name_found=false;
							k=0; name_length=AliasList.size(); 
							if (set[3]=="*") {name_found=true; wild_aliaslist=true;}
                                                        //if (set[1]=="var") name_found=true;

							while (k<name_length && !name_found) {
								if (AliasList[k]==set[3]) name_found=true; 
								k++;
							}
							if (!name_found) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(AliasList); 
							}
							break;
						case 4:
							if (set[3]=="*") set[3]=LatList[0];
							if (LatList[0]!=set[3]) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;

							//if (LatList[0]!=set[3] && set[1]!="var") {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(LatList); name_found=false;
							}
							break;
						case 5:
							if (set[3]=="*") set[3]=NewtonList[0];
							if (NewtonList[0]!=set[3]) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;

							//if (NewtonList[0]!=set[3] && set[1]!="var") {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(NewtonList); name_found=false;
							}
							break;
						case 6:
							if (set[3]=="*") set[3]=EngineList[0];
							if (EngineList[0]!=set[3]) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
							//if (EngineList[0]!=set[3] && set[1]!= "var") {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(EngineList); name_found=false;
							}
							break;
						case 7:
							k=0; name_length=OutputList.size(); name_found=false;
                                                        //if (set[1]=="var") name_found=true; 
							while (k<name_length) {
								if (OutputList[k]==set[3]) name_found=true; 
								k++;
							}
							if (!name_found) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(OutputList); 
							}
							break;
						case 8:
							name_found=false;
							k=0; name_length=VarList.size();  
							while (k<name_length && !name_found) {
								if (VarList[k]==set[3]) name_found=true; 
								k++;
							}
							if (!name_found) {cout << "In line " << set[0] << " name '" << set[3] << "' not recognised. Select from: "<< endl;
								PrintList(MonList); 
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
			if (!name_found) {return false;}   
			if (wild_aliaslist) {
				int length =AliasList.size();
				for (int i=0; i<length; i++) {Out_key.push_back(set[2]); Out_name.push_back(AliasList[i]); Out_prop.push_back(set[4]); }
			}
			if (wild_mollist) {
				int length =MolList.size();
				for (int i=0; i<length; i++) {Out_key.push_back(set[2]); Out_name.push_back(MolList[i]); Out_prop.push_back(set[4]); }
			}
			if (wild_monlist) {
				int length =MonList.size();
				for (int i=0; i<length; i++) {Out_key.push_back(set[2]); Out_name.push_back(MonList[i]); Out_prop.push_back(set[4]); }
			}
			if (!(wild_monlist || wild_mollist ||wild_aliaslist)) {Out_key.push_back(set[2]); Out_name.push_back(set[3]); Out_prop.push_back(set[4]);}
			wild_aliaslist=false;
			wild_mollist=false;
			wild_monlist=false;
			
			if (set[1]=="vtk" && Out_key.size()>1) {
				cout << "vtk output can have only one entry: the following entries were found:" << endl; success =false; 
				int length =Out_key.size();
				for (int i=0; i<length; i++) {cout << set[1] << " : " << Out_key[i] << " : " << Out_name[i] << " : " << Out_prop[i] << endl; }
			} 
		} //end temp
		i++;
	} //end i; 
	return success; 
}

bool Input:: CheckInput(void) {
	bool success=true;
	bool Keywordfound; 
	int key_length=KEYS.size();
	int length = elems.size();  
	int j;
	bool last_start=false; 
	
	string word;
	
	if (length==0) {cout << "inputfile is empty " << endl; success=false; } 
	int i=0;
	while (i<length) { 
		last_start=false;
		vector<std::string> set;
		split(elems[i],':',set); 
		if (set.size() !=5) {if (set[1]!="start") {
			cout <<elems[i] << endl; 
			cout << " Line number " << set[0] << " does not contain 4 items" << endl; success=false; return false;}
		}
		if (set[1]=="start") last_start=true;
		i++;
	}
	if (!last_start) {
		elems.push_back("0:start"); length++;
	}
	i=0;
	while (success && i<length) {
		vector<std::string> set;
		split(elems[i],':',set); 

		if (set[1]=="output") {
			vector<string> options; 
			options.push_back("ana"); options.push_back("vtk"); options.push_back("kal"); options.push_back("pro");
			string option;
			if (!Get_string(set[2],option,options,"Value for output extension '" + set[2] + "' not allowed. ")) {success=false;} 
			else {
				word=set[2];
				
				j=0; Keywordfound=false; key_length =KEYS.size();
 				if (word=="ana") Keywordfound=true; 
				while (j<key_length) { 
					if (word==KEYS[j]) {Keywordfound=true; }
					j++; 
				}
				if (!Keywordfound) {
					KEYS.push_back(word);
					key_length++;
				}
			}
			
		}

/*
                if (set[1]=="engine") {
                        vector<string> options;
                        options.push_back("sfbox"); 
                        string option;
                        if (!Get_string(set[4],option,options,"Value for Engine brand '" + set[4] + "' not allowed. ")) {success=false;}
                        else {
                                word=set[4];

                                j=0; Keywordfound=false; key_length =KEYS.size();
                                while (j<key_length) {
                                        if (word==KEYS[j]) {Keywordfound=true; }
                                        j++;
                                }
                                if (!Keywordfound) {
                                        KEYS.push_back(word);
                                        key_length++;
                                }
				

                        }

                }
*/


/*
			if (set[1]=="search"){
			vector<string> options;
			options.push_back("sys"); options.push_back("mol");
			string option;
			if (!Get_string(set[2],option,options, "Cannot search in ' " +set[2]+"'. Choose from options.")) {success=false;}
			else { 
				cout << "checking input" << endl; 
				key_length=MolList.size(); bool molfound=false;
				cout << key_length << endl;
				for (int j=0; j<key_length; j++) {cout << set[3]<<":" <<MolList[j] << endl; if(set[3]==MolList[j]) {cout << "Molecule '"+MolList[j]+"' value of '"+set[4]+"' will be searched to find target" << endl; molfound=true;}}
				if (set[2]==options[0] && !molfound) {cout << "check1" << endl; if (set[3]!="GrandPotential") {cout << "In 'sys' I can only search for GrandPotential(case sensitive) and not '"+set[3]+"'. Execution stopped." << endl; success=false; return 0;} }
				else if (set[2]==options[1] && !molfound){if(set[3]!="phibulk"){cout << "In 'mol' I can only serach for 'phibulk' and not '"+set[3]+"'. Executions stopped." << endl; success=false; return 0;}}
			}
		}
*/


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
	MakeLists(1);
	
	return success;
}

bool Input::MakeLists(int start) {
	bool success=true;
	SysList.clear();
	LatList.clear();
	NewtonList.clear();
	MonList.clear();
	MolList.clear();
	OutputList.clear();
	EngineList.clear();
	MesodynList.clear();
	VarList.clear();

	if (!TestNum(SysList,"sys",0,1,start)) {cout << "There can be no more than 1 'sys name' in the input" << endl; }
	if (SysList.size()==0) SysList.push_back("noname"); 
	if (!TestNum(LatList,"lat",1,1,start)) {cout << "There must be exactly one 'lat name' in the input" << endl; success=false;}
	if (!TestNum(NewtonList,"newton",0,1,start)) {cout << "There can be no more than 1 'newton name' in input" << endl; success=false;}
	if (NewtonList.size()==0) NewtonList.push_back("noname");
	if (!TestNum(NewtonList,"newton",0,1,start)) {cout << "There can be no more than 1 'newton name' in input" << endl; success=false;}
	if (!TestNum(MonList,"mon",1,1000,start)) {cout << "There must be at least one 'mon name' in input" << endl; success=false;}
	TestNum(AliasList,"alias",0,1000,start); 
	if (AliasList.size()==0) AliasList.push_back("noname"); 
	if (!TestNum(MolList,"mol",1,1000,start)) {cout << "There must be at least one 'mol name' in input" << endl; success=false;}
	if (!TestNum(OutputList,"output",1,1000,start)) {cout << "No output defined! " << endl;}
	if (!TestNum(EngineList,"engine",0,1,start)) {cout << "There can be no more than 1 'engine' name in the input " << endl; success=false;}
	if (!TestNum(MesodynList,"mesodyn",0,1,start)) {cout << "There can be no more than 1 'engine' name in the input " << endl; success=false;}
	if (EngineList.size()==0) EngineList.push_back("noname");
	if (!TestNum(VarList,"var",0,10,start)) 
	if (VarList.size()==0) VarList.push_back("noname");
	return success;
}

