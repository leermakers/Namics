#include "mesodyn.h"
#include <random>

Mesodyn::Mesodyn(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,vector<System*> Sys_,vector<Newton*> New_, string name_) {
	In=In_; name=name_;   Lat=Lat_; Mol=Mol_; Seg = Seg_; Sys=Sys_; New=New_;
	KEYS.push_back("timesteps");
	KEYS.push_back("timebetweensaves");
}
Mesodyn::~Mesodyn() {
}

void Mesodyn::AllocateMemory() {
	if (debug) cout <<"nothing to allocate in Mesodyn" << endl;
}


void Mesodyn::PutParameter(string new_param) {
	KEYS.push_back(new_param);
}

bool Mesodyn::CheckInput(int start) {
if (debug) cout <<"Check Mesodyn" << endl;
	bool success=true;
//right now all the checks for engine are done in input.cpp. could be move back here.
	success=In[0]->CheckParameters("mesodyn",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		vector<string> options;
		if (GetValue("timesteps").size()>0) {
                   success=In[0]->Get_int(GetValue("timesteps"),timesteps, 1, 10000, "The number of timesteps should be between 1 and 10000");
		}
	cout <<"timesteps is " << timesteps << endl;

	}
	return success;
}

/* Generates a vector of length count, contianing gaussian noise of given mean, standard deviation.
	 Noise is stored in vector<Real> Mesodyn::noise
	 Possible errors: What if count > sizeof(unsinged long)?
	 Called by langevinFlux()
*/
void Mesodyn::gaussianNoise(Real mean, Real stdev, unsigned long count) {

	random_device generator;

	//Mersenne Twister 19937 bit state size PRNG
  mt19937 prng(generator());

  normal_distribution<> dist(mean, stdev);

	this->noise.resize(count);

  for (unsigned int i = 0; i < count; ++i) {
  	this->noise[i] = dist(prng);
 	}

/* Debugging code (output value in all elements):
	for (auto const &element: mesodyn.thisNoise)
				std::cout << element << ' ';
*/
}

void Mesodyn::langevinFlux() {
	//Flux + gaussianNoise(Real, Real, unsigned long);
}

string Mesodyn::GetValue(string parameter){
	int i=0;
	int length = PARAMETERS.size();
	while (i<length) {
		if (parameter==PARAMETERS[i]) {
			return VALUES[i];
		}
		i++;
	}
	return "" ;
}

void Mesodyn::push(string s, Real X) {
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void Mesodyn::push(string s, int X) {
	ints.push_back(s);
	ints_value.push_back(X);
}
void Mesodyn::push(string s, bool X) {
	bools.push_back(s);
	bools_value.push_back(X);
}
void Mesodyn::push(string s, string X) {
	strings.push_back(s);
	strings_value.push_back(X);
}
void Mesodyn::PushOutput() {
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
}
Real* Mesodyn::GetPointer(string s) {
	//vector<string> sub;
	//nothing yet
	return NULL;
}

int Mesodyn::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
	int i=0;
	int length = ints.size();
	while (i<length) {
		if (prop==ints[i]) {
			int_result=ints_value[i];
			return 1;
		}
		i++;
	}
	i=0;
	length = Reals.size();
	while (i<length) {
		if (prop==Reals[i]) {
			Real_result=Reals_value[i];
			return 2;
		}
		i++;
	}
	i=0;
	length = bools.size();
	while (i<length) {
		if (prop==bools[i]) {
			if (bools_value[i]) string_result="true"; else string_result="false";
			return 3;
		}
		i++;
	}
	i=0;
	length = strings.size();
	while (i<length) {
		if (prop==strings[i]) {
			string_result=strings_value[i];
			return 3;
		}
		i++;
	}
	return 0;
}
