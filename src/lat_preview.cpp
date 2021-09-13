#include "lat_preview.h"
Lat_preview::Lat_preview(vector<Input*> In_,string name_) {
if (debug) cout <<"Lattice constructor" << endl;
	In=In_; name=name_;
	KEYS.push_back("gradients");
	KEYS.push_back("n_layers");
	KEYS.push_back("offset_first_layer");
	KEYS.push_back("geometry");
	KEYS.push_back("n_layers_x");
	KEYS.push_back("n_layers_y");
	KEYS.push_back("n_layers_z");
	KEYS.push_back("lowerbound");
	KEYS.push_back("upperbound");
	KEYS.push_back("lowerbound_x");
	KEYS.push_back("upperbound_x");
	KEYS.push_back("lowerbound_y");
	KEYS.push_back("upperbound_y");
	KEYS.push_back("lowerbound_z");
	KEYS.push_back("upperbound_z");
 	KEYS.push_back("bondlength");
	KEYS.push_back("ignore_site_fraction");
	KEYS.push_back("fcc_site_fraction");
  	KEYS.push_back("lattice_type");
	KEYS.push_back("stencil_full");
	KEYS.push_back("FJC_choices");
	KEYS.push_back("b/l");
	//KEYS.push_back("Markov");
	KEYS.push_back("k_stiff");
}


Lat_preview::~Lat_preview() {
if (debug) cout <<"lat_preview destructor " << endl;

}


bool Lat_preview::CheckInput(int start) {
if (debug) cout <<"CheckInput in lattice " << endl;
	bool success;
	string Value;

	success = In[0]->CheckParameters("lat",name,start,KEYS,PARAMETERS,VALUES);
	if (!success) return success;

	vector<string> options;

	options.clear();
	options.push_back("spherical");
	options.push_back("cylindrical");
	options.push_back("flat");
	options.push_back("planar");

	if (GetValue("geometry").size()>0) {
		if (!In[0]->Get_string(GetValue("geometry"),geometry,options,"In lattice input for 'geometry' not recognized."))
		success=false;
	} else
	geometry = "planar";
	if (geometry=="flat") geometry="planar";


	gradients=In[0]->Get_int(GetValue("gradients"),1);
	if (gradients<0||gradients>3) {cout << "value of gradients out of bounds 1..3; default value '1' is used instead " << endl; gradients=1;}
	return success;
}


string Lat_preview::GetValue(string parameter){
if (debug) cout << "GetValue in lattice " << endl;
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



