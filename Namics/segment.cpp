class Segment {
public:
	Segment(vector<Input*>,string,int,int);

~Segment();

	string name;
	vector<Input*> In; 
	//vector<Lattice*> Lat;
	int n_seg; 
	int seg_nr;  
	double epsilon;
	double valence; 
	string freedom;
		

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES; 
	bool CheckInput(); 
	void PutChiKEY(string); 
	string GetValue(string); 

};
Segment::Segment(vector<Input*> In_, string name_,int segnr,int N_seg) {
	In=In_, name=name_; n_seg=N_seg; seg_nr=segnr; 
	KEYS.push_back("freedom"); 
	KEYS.push_back("valence");
	KEYS.push_back("epsilon");
	KEYS.push_back("tagged");
 
}
Segment::~Segment() {
}
void Segment::PutChiKEY(string new_name) {
	if(name != new_name) KEYS.push_back("chi-" + new_name); 
}


bool Segment::CheckInput() {
	return In[0]->CheckParameters("mon",name,KEYS,PARAMETERS,VALUES);
}

string Segment::GetValue(string parameter) {
	int length = PARAMETERS.size(); 
	int i=0;
	while (i<length) {
		if (PARAMETERS[i]==parameter) { return VALUES[i];}		
		i++;
	} 
	return " "; 
}
 
