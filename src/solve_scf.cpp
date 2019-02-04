#include "solve_scf.h"
#include <iostream>

Solve_scf::Solve_scf(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_, vector<State*> Sta_, vector<Reaction*> Rea_, vector<Molecule*> Mol_,vector<System*> Sys_,vector<Variate*>Var_,string name_) {
	In=In_; name=name_; Sys=Sys_; Seg=Seg_; Lat=Lat_; Mol=Mol_;Var=Var_;  Sta=Sta_; Rea=Rea_;
if(debug) cout <<"Constructor in Solve_scf " << endl;
	KEYS.push_back("method");
	KEYS.push_back("gradient_type");
	KEYS.push_back("e_info"); KEYS.push_back("s_info");KEYS.push_back("i_info");KEYS.push_back("t_info");
	KEYS.push_back("iterationlimit" ); KEYS.push_back("tolerance");
	KEYS.push_back("stop_criterion");
	KEYS.push_back("deltamin");KEYS.push_back("deltamax");
	KEYS.push_back("linesearchlimit");
	//KEYS.push_back("samehessian");
	KEYS.push_back("max_accuracy_for_hessian_scaling");
	KEYS.push_back("n_iterations_for_hessian");
	KEYS.push_back("small_alpha");
	KEYS.push_back("target_function");
	KEYS.push_back("max_n_small_alpha");
	KEYS.push_back("min_accuracy_for_hessian");
	KEYS.push_back("max_fr_reverse_direction");
	KEYS.push_back("print_hessian_at_it");
	KEYS.push_back("super_e_info");
	KEYS.push_back("super_s_info");
	KEYS.push_back("super_i_info");
	KEYS.push_back("super_tolerance");
	KEYS.push_back("super_iterationlimit");
	KEYS.push_back("m");
	KEYS.push_back("super_deltamax");
	max_g = false; // compute g based on max error
	rescue_status = NONE;
}

Solve_scf::~Solve_scf() {
	DeAllocateMemory();
}

void Solve_scf :: DeAllocateMemory(){
if (debug) cout <<"Destructor in Solve " << endl;

#ifdef CUDA
	cudaFree(xx);
	cudaFree(x0);
	cudaFree(g);
	cudaFree(xR);
	cudaFree(x_x0);
#else
	free(xx);
#endif

}

void Solve_scf::AllocateMemory() {
if(debug) cout <<"AllocateMemeory in Solve " << endl;
	int M=Lat[0]->M;
	if (mesodyn) {
		iv = Sys[0]->SysMolMonList.size()*M;
	} else {
		iv = (Sys[0]->ItMonList.size() + Sys[0]->ItStateList.size())* M;
	}
	if (Sys[0]->charged) iv += M;
	if (SCF_method=="Picard") iv += M;
	
	if (Sys[0]->beta_set) iv++; 
#ifdef CUDA
	xx  = (Real*)AllOnDev(iv);
	x0  = (Real*)AllOnDev(iv);
	g   = (Real*)AllOnDev(iv);
	xR  = (Real*)AllOnDev(m*iv);
	x_x0= (Real*)AllOnDev(m*iv);
#else
	xx=(Real*) malloc(iv*sizeof(Real));
#endif
	Zero(xx,iv);
	Sys[0]->AllocateMemory();
}

bool Solve_scf::PrepareForCalculations() {
if (debug) cout <<"PrepareForCalculations in Solve " << endl;
	bool success=true;
	return success;
}

bool Solve_scf::CheckInput(int start_) { start=start_;
if(debug) cout <<"CheckInput in Solve " << endl;
	pseudohessian =false;
	deltamin =0.1;
	s_info=false;
	e_info=false;
	t_info=false;
	i_info=1;
	hessian =false;
	bool success=true;
	control=proceed;
	mesodyn=false;
	string value;
	solver=PSEUDOHESSIAN;
	SCF_method="pseudohessian";
	gradient=classical;
	residual=1;
	m=10;
	success=In[0]->CheckParameters("newton",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		iterationlimit=In[0]->Get_int(GetValue("iterationlimit"),1000);
		if (iterationlimit < 0 || iterationlimit>1e6) {iterationlimit = 1000;}

		e_info=In[0]->Get_bool(GetValue("e_info"),true); value_e_info=e_info;
		s_info=In[0]->Get_bool(GetValue("s_info"),false); value_s_info =s_info;
		t_info=In[0]->Get_bool(GetValue("t_info"),false);
		i_info=In[0]->Get_int(GetValue("i_info"),1);
		if (i_info == 0) {
		// We cannot divide by zero (see modulus statements in sfnewton), but this will probably be what the user means.
		cerr << "WARNING: i_info cannot be zero ! Defaulting to iterationlimit + 1."<< endl;
		i_info = iterationlimit+1;
		}
		value_i_info=i_info;

		super_e_info=In[0]->Get_bool(GetValue("super_e_info"),e_info);
		super_s_info=In[0]->Get_bool(GetValue("super_s_info"),s_info);
		super_i_info=In[0]->Get_bool(GetValue("super_i_info"),i_info);
		super_iterationlimit=In[0]->Get_int(GetValue("super_iterationlimit"),iterationlimit/10);

		if (GetValue("target_function").size() > 0) {
			string target;
      target = In[0]->Get_string(GetValue("target_function"), target);
			using namespace std::placeholders;
			if ( target.find("log") != string::npos  ) target_function = bind(&Solve_scf::gradient_log, this, _1, _2, _3, _4, _5);
			else if ( target.find("quotient") != string::npos  ) target_function = bind(&Solve_scf::gradient_quotient, this, _1, _2, _3, _4, _5);
			else if ( target.find("minus") != string::npos  ) target_function = bind(&Solve_scf::gradient_minus, this, _1, _2, _3, _4, _5);
			else {
				cerr << "Target function not found, please choose from log, quotient or minus. Defaulting to minus." << endl;
			}
		} else {
			using namespace std::placeholders;
			target_function = bind(&Solve_scf::gradient_minus, this, _1, _2, _3, _4, _5);
		}

		deltamax=In[0]->Get_Real(GetValue("deltamax"),0.1);
		super_deltamax=In[0]->Get_Real(GetValue("super_deltamax"),0.5);
		if (deltamax < 0 || deltamax>100) {deltamax = 0.1;  cout << "Value of deltamax out of range 0..100, and value set to default value 0.1" <<endl; }
		deltamin=0; super_deltamin=0;
		deltamin=In[0]->Get_Real(GetValue("deltamin"),deltamin);
		if (deltamin < 0 || deltamin>100) {deltamin = deltamax/100000;  cout << "Value of deltamin out of range 0..100, and value set to default value deltamax/100000" <<endl; }
		tolerance=In[0]->Get_Real(GetValue("tolerance"),1e-7);
		super_tolerance=In[0]->Get_Real(GetValue("super_tolerance"),tolerance*10);
		if (tolerance < 1e-12 ||tolerance>10) {tolerance = 1e-5;  cout << "Value of tolerance out of range 1e-12..10 Value set to default value 1e-5" <<endl; }
		if (GetValue("method").size()==0) {SCF_method="pseudohessian";} else {
			vector<string>method_options;
			method_options.push_back("DIIS");
			method_options.push_back("Picard"); //can be included again when adjusted for charges and guess
			method_options.push_back("pseudohessian");
			method_options.push_back("hessian");
			method_options.push_back("conjugate_gradient");
			if (!In[0]->Get_string(GetValue("method"),SCF_method,method_options,"In 'solve_scf' the entry for 'method' not recognized: choose from:")) success=false;
		}
		if (SCF_method=="hessian" || SCF_method=="pseudohessian") {
			if (SCF_method=="hessian") {pseudohessian=false; hessian=true; solver=HESSIAN;} else { pseudohessian=true; hessian=false; solver=PSEUDOHESSIAN;}
			samehessian=false; //In[0]->Get_bool(GetValue("samehessian"),false);
			max_accuracy_for_hessian_scaling=In[0]->Get_Real(GetValue("max_accuracy_for_hessian_scaling"),0.1);
			if (max_accuracy_for_hessian_scaling<1e-7 || max_accuracy_for_hessian_scaling>1) {
				cout <<"max_accuracy_for_hessian_scaling is out of range: 1e-7...1; default value 0.1 is used instead" << endl;
				max_accuracy_for_hessian_scaling=0.1;
			}
			minAccuracyForHessian=In[0]->Get_Real(GetValue("min_accuracy_for_hessian"),0);
			if (minAccuracyForHessian<0 ||minAccuracyForHessian>1) {
				cout <<"min_accuracy_for_hessian is out of range: 0...0.1; default value 0 is used instead (no hessian computation)" << endl;
				minAccuracyForHessian=0;
			}
			maxFrReverseDirection =In[0]->Get_Real(GetValue("max_fr_reverse_direction"),0.4);
			if (maxFrReverseDirection <0.1 ||maxFrReverseDirection >0.5) {
				cout <<"max_fr_reverse_direction is out of range: 0.1...0.5; default value 0.4 is used instead" << endl;
				maxFrReverseDirection =0.4;
			}

			n_iterations_for_hessian=In[0]->Get_int("n_iterations_for_hessian",100);
			if (n_iterations_for_hessian<1 ||n_iterations_for_hessian>1000) {
				cout <<" n_iterations_for_hessian setting is out of range: 1, ..., 1000; hessian evaluations will not be done " << endl;
				n_iterations_for_hessian=iterationlimit+100;
			}
			maxNumSmallAlpha=In[0]->Get_int("max_n_small_alpha",50);
			if (maxNumSmallAlpha<10 ||maxNumSmallAlpha>1000) {
				cout <<" max_n_small_alpha is out of range: 10, ..., 100;  max_n_small_alpha is set to default: 50 " << endl;
				maxNumSmallAlpha=50;
			}

			deltamin=In[0]->Get_Real(GetValue("delta_min"),0);
			if (deltamin <0 || deltamin>deltamax) {
				cout <<"delta_min is out of range; 0, ..., " << deltamax << "; delta_min value set to 0 " << endl;
				deltamin=0;
			}
			smallAlpha=In[0]->Get_Real(GetValue("small_alpha"),0.00001);
			if (smallAlpha <0 || smallAlpha>1) {
				cout <<"small_alpha is out of range; 0, ..., 1; small_alpha value set to default: 1e-5 " << endl;
				smallAlpha=0.00001;
			}
		}
		if (SCF_method=="DIIS") {
			solver=diis;
			m=In[0]->Get_int(GetValue("m"),10);
			if (m < 0 ||m>100) {m=10;  cout << "Value of 'm' out of range 0..100, value set to default value 10" <<endl; }
		}
		if (SCF_method=="Picard") {
			solver= PICARD;
			gradient=Picard;
		}
		if (SCF_method=="conjugate_gradient") {
			solver= conjugate_gradient;
			linesearchlimit=In[0]->Get_int(GetValue("linesearchlimit"),linesearchlimit);
		}
		if (GetValue("gradient_type").size()==0) {gradient=classical;} else {
			vector<string>gradient_options;
			gradient_options.push_back("classical");
			gradient_options.push_back("mesodyn");
			gradient_options.push_back("Picard");
			if (!In[0]->Get_string(GetValue("gradient_type"),gradients,gradient_options,"In 'solve_scf' the entry for 'gradient_type' not recognized: choose from:")) success=false;
			if (gradients=="classical") gradient=classical;
			if (gradients=="mesodyn") {gradient = MESODYN; mesodyn=true;}
			if (gradients=="Picard")  gradient=Picard;
		}

		StoreFileGuess=In[0]->Get_string(GetValue("store_guess"),"");
		ReadFileGuess=In[0]->Get_string(GetValue("read_guess"),"");
		if (GetValue("stop_criterion").size() > 0) {
			vector<string>options;
			options.push_back("norm_of_g");
			options.push_back("max_of_element_of_|g|");
			if (!In[0]->Get_string(GetValue("stop_criterion"),stop_criterion,options,"In newton the stop_criterion setting was not recognised")) {success=false; };
			if(GetValue("stop_criterion") == options[1]) {
				max_g = true;
			}
		}
	}
	return success;
}


void Solve_scf::PutParameter(string new_param) {
if(debug) cout <<"PutParameter in Solve " << endl;
	KEYS.push_back(new_param);
}

string Solve_scf::GetValue(string parameter){
if(debug) cout <<"GetValue " + parameter + " in  Solve " << endl;
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

void Solve_scf::push(string s, Real X) {
if(debug) cout <<"push (Real) in  Solve " << endl;
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void Solve_scf::push(string s, int X) {
if(debug) cout <<"push (int) in  Solve " << endl;
	ints.push_back(s);
	ints_value.push_back(X);
}
void Solve_scf::push(string s, bool X) {
if(debug) cout <<"push (bool) in  Solve " << endl;
	bools.push_back(s);
	bools_value.push_back(X);
}
void Solve_scf::push(string s, string X) {
if(debug) cout <<"push (string) in  Solve " << endl;
	strings.push_back(s);
	strings_value.push_back(X);
}
void Solve_scf::PushOutput() {
if(debug) cout <<"PushOutput in  Solve " << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();
	push("method",SCF_method);
	push("m",m);
	push("delta_max",deltamax);
	push("residual",residual);
	push("tolerance",tolerance);
	push("iterations",iterations);
	push("iterationlimit",iterationlimit);
	push("stop_criterion",stop_criterion);
	if (pseudohessian || hessian) {
		//push("same_hessian",samehessian);
		push("linesearchlimit",linesearchlimit);
		push("max_accuracy_for_hessian_scaling",max_accuracy_for_hessian_scaling);
		push("n_iteratons_for_hessian",n_iterations_for_hessian);
		push("small_alpha",smallAlpha);
		push("max_n_small_alpha",maxNumSmallAlpha);
		push("min_accuracy_for_hessian",minAccuracyForHessian);
	}
	Lat[0]->PushOutput();
	int length = In[0]->MonList.size();
	for (int i=0; i<length; i++) Seg[i]->PushOutput();
	length = In[0]->MolList.size();
	for (int i=0; i<length; i++){
		int al_length=Mol[i]->MolAlList.size();
		for (int k=0; k<al_length; k++) {
			Mol[i]->Al[k]->PushOutput();
		}
		Mol[i]->PushOutput();
	}
	length = In[0]->StateList.size();
	for (int i=0; i<length; i++) Sta[i]->PushOutput();
	length = In[0]->ReactionList.size();
	for (int i=0; i<length; i++) Rea[i]->PushOutput();
	Sys[0]->PushOutput();
}

int Solve_scf::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
if(debug) cout <<"GetValue (long) in  Solve " << endl;
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

void Solve_scf::Copy(Real* x, Real* X, int MX, int MY, int MZ, int fjc_old) {
	int mx=Lat[0]->MX;
	int my=Lat[0]->MY;
	int mz=Lat[0]->MZ;
	int jx=Lat[0]->JX;
	int jy=Lat[0]->JY;
	int i,j,k;
	int pos_i,pos_o;
	int JX=(MY+2)*(MZ+2);
	int JY=(MZ+2);
	Real Xvalue;
	int fjc=Lat[0]->fjc;

	switch (Lat[0]->gradients) {
		case 1:
			if (fjc==1 and fjc_old==1) {
			if (MY>0||MZ>0) {
				cout <<" Copy from more than one gradient to one gradient: (i) =(1,i) or (1,1,i) is used "<< endl;
			}
			if (MZ>0) { pos_i=JX+JY; pos_o=MZ+2;} else {if (MY>0) {pos_i=JX; pos_o=MY+2; } else { pos_i=0; pos_o=MX+2; } }
			for (i=0; i<mx+2; i++)  if (i<pos_o) x[i]=X[pos_i+i];
			} else {
				for (i=0; i<mx+2; i++) {
					Xvalue=0; pos_i=0; pos_o=MX+2;
					if (i<pos_o) {
						for (j=0; j<fjc_old; j++) Xvalue+=X[i*fjc_old+j];
						Xvalue/=fjc_old;
						for (j=0; j<fjc; j++) {if (fjc!=fjc_old) x[i*fjc+j]=Xvalue; else  x[i*fjc+j]= X[i*fjc+j];}
					} else for (j=0; j<fjc; j++) x[i*fjc+j]=0;
				}
			}
			break;
		case 2:
			if (MY==0) {
				cout <<" Copy from one-gradient to two gradients: one-gradient (x,i)=(i) is used for all x " << endl;
				for (i=0; i<mx+2; i++)
				for (j=0; j<my+2; j++) if (j<MX+2) x[i*jx+j]=X[j];
			} else {
				if (MZ>0) {
					cout <<" Copy from three gradients to two gradients: (i,j)=(1,i,j) is used " <<endl;
					JX=(MY+2)*(MZ+2);
					JY=(MZ+2);
					for (i=0; i<mx+2; i++)
					for (j=0; j<my+2; j++) if (i<MY+2 && j<MZ+2) x[i*jx+j]=X[JX+i*JY+j];
				} else {
					JX=(MY+2);
					for (i=0; i<mx+2; i++)
					for (j=0; j<my+2; j++) if (i<MX+2 && j<MY+2) x[i*jx+j]=X[i*JX+j];
				}
			}
			break;
		case 3:
			if (MY==0) {
				cout <<"Copy from one gradient to three gradients: (x,y,i) = (i) is used for all x,y " << endl;
				for (i=0; i<mx+2; i++)
				for (j=0; j<my+2; j++)
				for (k=0; k<mz+2; k++) if (k<MX+2) x[i*jx+j*jy+k]=X[k];
			} else {
				if (MZ==0) {
					cout <<"Copy form two gradients to three: (x,i,j) = (i,j) for all x " << endl;
					JX=(MY+2);
					for (i=0; i<mx+2; i++)
					for (j=0; j<my+2; j++)
					for (k=0; k<mz+2; k++) if (j<MX+2 && k<MY+2) x[i*jx+j*jy+k]=X[j*JX+k];
				} else {
					JX=(MY+2)*(MZ+2);
					JY=(MZ+2);
					for (i=0; i<mx+2; i++)
					for (j=0; j<my+2; j++)
					for (k=0; k<mz+2; k++) if (i<MX+2 && j<MY+2 && k<MZ+2) x[i*jx+j*jy+k]=X[i*JX+j*JY+k];
				}
			}
			break;
		default:
			break;
	}
}

bool Solve_scf::Guess(Real *X, string METHOD, vector<string> MONLIST, vector<string> STATELIST, bool CHARGED, int MX, int MY, int MZ,int fjc_old){
	if (debug) cout << "Guess in Solve" << endl;
	int M=Lat[0]->M;
	bool success=true;

	if (start ==1 && Sys[0]->GuessType != "")  {
		Lat[0]->GenerateGuess(xx,Sys[0]->CalculationType,Sys[0]->GuessType,Seg[Sys[0]->MonA]->guess_u,Seg[Sys[0]->MonB]->guess_u);
	} else {
		int m;
		if (MZ>0) {m=(MX+2)*(MY+2)*(MZ+2); } else { if (MY>0) { m=(MX+2)*(MY+2); } else {  m=(MX+2);}}
		if (fjc_old>1) m*=fjc_old;
		int length_old_mon=MONLIST.size();
		int length_old_state=STATELIST.size();
		int length_new_mon=Sys[0]->ItMonList.size();
		int length_new_state=Sys[0]->ItStateList.size();
		for (int i = 0; i<length_old_mon; i++) {
			for (int j=0; j<length_new_mon; j++) {
				if (MONLIST[i]==Seg[Sys[0]->ItMonList[j]]->name) {
					Copy(xx+M*j,X+i*m,MX,MY,MZ,fjc_old);
				}
			}
		}
		for (int i = 0; i<length_old_state; i++) {
			for (int j=0; j<length_new_state; j++) {
				if (STATELIST[i]==Sta[Sys[0]->ItStateList[j]]->name) {
					Copy(xx+M*(j+length_new_mon),X+(i+length_old_mon)*m,MX,MY,MZ,fjc_old);
				}
			}
		}

		if (CHARGED && Sys[0]->charged) {cout <<"both charged" << endl;
			Copy(xx+(length_new_mon+length_new_state)*M,X+(length_old_mon+length_old_state)*m,MX,MY,MZ,fjc_old);
		}
	}
	return success;
}

bool Solve_scf::Solve(bool report_errors_) { //going SCF here
if(debug) cout <<"Solve in  Solve_scf " << endl;
	bool success=true;
	bool report_errors=report_errors_;
	int niv = In[0]->ReactionList.size();
	if (niv>0) {
		int i_solver=0;
		if (solver==HESSIAN) i_solver=1;
		if (solver==PSEUDOHESSIAN) i_solver=2;
		if (solver==diis) i_solver=3;
		bool ee_info, ss_info;
		if (e_info) ee_info=true; else ee_info=false; e_info=false;
		if (s_info) ss_info=true; else ss_info=false; s_info=false;
		gradient = WEAK;
		control= super;
		Real* yy=(Real*) malloc((iv)*sizeof(Real)); Cp(yy,xx,iv);
		SIGN=(int*) malloc((niv)*sizeof(int)); for (int i=0; i<niv; i++) SIGN[i]=1.0;
		pseudohessian=false;  hessian =true;
		Zero(xx,niv);
		success=iterate(xx,niv,100,1e-8,0.5,deltamin,true);
		if (!success) cout <<"iteration for alphabuk values for internal states failed. Check eqns. " << endl;
		e_info=ee_info;
		s_info=ss_info;
//int n_states=In[0]->StateList.size();
//for (int i=0; i<n_states; i++) cout <<Seg[Sta[i]->mon_nr]->state_alphabulk[Sta[i]->state_nr] << endl;
		if (i_solver==1) solver=HESSIAN;
		if (i_solver==2) {solver=PSEUDOHESSIAN; pseudohessian=true;}
		if (i_solver==3) solver=diis;
		gradient = classical;
		control = proceed;
		Cp(xx,yy,iv);
/*
		int n_segments=In[0]->MonList.size();
		for (int i=0; i<n_segments; i++) {
			int ns=Seg[i]->ns;
			if (ns>1) {
				cout <<"segment " << In[0]->MonList[i] << endl;
				for (int k=0; k<ns; k++) cout << "state " << k << ":" <<Seg[i]->state_alphabulk[k] << endl;
			}
		}
*/
		free(yy);
		free(SIGN);
	}
	//gradient=classical;
	//control=proceed;
	switch(solver) {
		case HESSIAN:
			success=iterate(xx,iv,iterationlimit,tolerance,deltamax,deltamin,true);
		break;
		case PSEUDOHESSIAN:
			success=iterate(xx,iv,iterationlimit,tolerance,deltamax,deltamin,true);
		break;
		case PICARD:
			success=iterate_Picard(xx,iv,iterationlimit,tolerance,deltamax);
		break;
		case diis:
			success=iterate_DIIS(xx,iv,m,iterationlimit,tolerance,deltamax);
		break;
		default:
			cout <<"Solve is lost" << endl; success=false;
		break;
	}
	success=Sys[0]->CheckResults(report_errors);
	return success;
}

bool Solve_scf::attempt_DIIS_rescue() {
	cout << "Attempting rescue!" << endl;
	switch (rescue_status) {
		case NONE:
			cout << "Zeroing iteration variables." << endl;
			Zero(xx,iv);
			rescue_status = ZERO;
			break;
		case ZERO:
			cout << "Adjusting memory depth." << endl;
			m *= 0.5;
			cout << "Zeroing iteration variables." << endl;
			Zero(xx,iv);
			rescue_status = M;
			break;
		case M:
			cout << "Decreasing delta_max." << endl;
			deltamax *= 0.1;
			cout << "Zeroing iteration variables." << endl;
			Zero(xx,iv);
			rescue_status = DELTA_MAX;
			break;
		case DELTA_MAX:
			cerr << "Exhausted all rescue options. Crash is imminent, exiting." << endl;
			exit(0);
			break;
	}
	return true;
}

bool Solve_scf::SolveMesodyn(function< void(Real*, size_t) > alpha_callback, function< Real*() > flux_callback) {
	if(debug) cout <<"Solve (mesodyn) in  Solve_scf " << endl;
	//iv should have been set at AllocateMemory.
	mesodyn_flux = flux_callback;
	mesodyn_load_alpha = alpha_callback;
  mesodyn =true;
	gradient=MESODYN;

	bool success=true;
	Real old_deltamax = deltamax;

	switch (solver) {
		case diis:
			gradient=MESODYN;
			try {
				success=iterate_DIIS(xx,iv,m,iterationlimit,tolerance,deltamax);
			}
			catch (int error) {
				if (error == -1)
				{
					cerr << "Detected GN not larger than 0." << endl;
					attempt_DIIS_rescue();
					cout << "Restarting iteration." << endl;
					SolveMesodyn(alpha_callback, flux_callback);
				}
				if (error == -2)
				{
					cerr << "Detected nan in U in Ax." << endl;
					attempt_DIIS_rescue();
					cout << "Restarting iteration." << endl;
					SolveMesodyn(alpha_callback, flux_callback);
				}
				if (error == -3)
				{
					cerr << "Detected nan in svdcmp." << endl;
					attempt_DIIS_rescue();
					cout << "Restarting iteration." << endl;
					SolveMesodyn(alpha_callback, flux_callback);
				}
				if (error == -4)
				{
					cerr << "Detected negative phibulk." << endl;
					attempt_DIIS_rescue();
					cout << "Restarting iteration." << endl;
					SolveMesodyn(alpha_callback, flux_callback);
				}
			}
			if (success == false) {
				cerr << "Detected failure to converge" << endl;
				attempt_DIIS_rescue();
				cout << "Restarting iteration." << endl;
				SolveMesodyn(alpha_callback, flux_callback);
			}	else {
				//m = old_m;
				deltamax = old_deltamax;
				rescue_status = NONE;
			}
			break;
		case PSEUDOHESSIAN:
			success=iterate(xx,iv,iterationlimit,tolerance,deltamax,deltamin,true);
		break;
		case conjugate_gradient:
			SFNewton::conjugate_gradient(xx,iv,iterationlimit,tolerance);
		break;
		default:
			success = false; cout <<" in SolveMesodyn the iteration method is unknown. " << endl;
		break;
	}

	/*if (Sys[0]->charged) {
		Sys[0]->DoElectrostatics(alpha+sysmon_length*M,xx+sysmon_length*M);
		Lat[0]->UpdateEE(Sys[0]->EE,Sys[0]->psi,Sys[0]->eps);
	}*/


/*		if (Sys[0]->charged){
			YplusisCtimesX(alpha+i*M,Sys[0]->EE,Seg[Sys[0]->SysMolMonList[i]]->epsilon,M);
			if (Seg[Sys[0]->SysMolMonList[i]]->valence !=0)
			YplusisCtimesX(alpha+i*M,Sys[0]->psi,-1.0*Seg[Sys[0]->SysMolMonList[i]]->valence,M);
		}*/

	Sys[0]->CheckResults(false);
	return success;
}


bool Solve_scf::SuperIterate(int search, int target,int ets,int etm, int bm) {
if(debug) cout <<"SuperIteration in  Solve_scf " << endl;
	Real* x = (Real*) malloc(sizeof(Real)); H_Zero(x,1);
	bool success=true;
	value_search=search;				//cp superiteration cnontrols to solve. These can be used in residuals.
	value_target=target;
	value_ets=ets;
	value_etm=etm;
	value_bm =bm;
	if (ets==-1 && etm==-1 && bm ==-1)  {                  //pick up initial guess;
		x[0] =Var[search]->GetValue();
	} else {
		if (ets>-1) x[0] = Var[ets]->GetValue();
		if (etm>-1) x[0] = Var[etm]->GetValue();
		if (bm>-1) { x[0] = Var[bm]->GetValue(); super_tolerance *=10;}
	}
	cout <<"Your guess for X: " << x[0] << endl; 
	gradient=custum; 				//the proper gradient is used
	control=super;				//this is for inneriteration
	e_info=super_e_info;
	s_info=super_s_info;
	i_info=super_i_info;
	tolerance=super_tolerance;
	cout << "super tolerance and super deltamax: " << super_tolerance << "\t" << super_deltamax << endl;
	solver=diis;
	    if (ets==-1 && etm==-1 && bm ==-1) success=iterate_RF(x,1,iterationlimit,super_tolerance,super_deltamax,"Regula-Falsi search: ");
	    if (ets>-1) success=iterate_RF(x,1,iterationlimit,super_tolerance,super_deltamax,"Regula-Falsi Eq-to_solvent search: ");
	    if (etm>-1) success=iterate_RF(x,1,iterationlimit,super_tolerance,super_deltamax,"Regula-Falsi Eq-to_mu search: ");
	    if (bm>-1) success=iterate_RF(x,1,iterationlimit,super_tolerance,super_deltamax,"Regula-Falsi balance-membrane search: ");
	//    success=iterate_DIIS(x,1,m,iterationlimit,super_tolerance,super_deltamax);
	//success=iterate(x,1,super_iterationlimit,super_tolerance,super_deltamax,deltamin,false);	//iterate is called with just one iteration variable
	if (bm>-1) super_tolerance /=10;
	cout <<"My guess for X: " << x[0] << endl;
	return success;
}


void Solve_scf::residuals(Real* x, Real* g){
 if (debug) cout <<"residuals in Solve_scf " << endl;
	int M=Lat[0]->M;
	Real chi;
	//Real valence;
	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size(); //also frozen segments
	int i,k,xi;

	switch(gradient) {
		case WEAK:
			if (debug) cout <<"Residuals for weak iteration " << endl;

			xi=0;
			Zero(g,In[0]->ReactionList.size());
			//for (i=0; i<In[0]->ReactionList.size(); i++) {
			//	x[i]=0.5*(1.0+tanh(x[i]));
			//}

			for (i=0; i<sysmon_length; i++) {
				xi=Seg[Sys[0]->SysMonList[i]]->PutAlpha(x,xi);
			}

			for (size_t i = 0; i<In[0]->ReactionList.size(); i++) {
				g[i]=SIGN[i]*Rea[i]->Residual_value();
			}

		break;
		case MESODYN:
		{
			if (debug) cout << "Residuals for mesodyn in Solve_scf " << endl;
			ComputePhis();
			#ifdef CUDA
			Real* temp_alpha = (Real*)AllOnDev(M);
			#else
			Real* temp_alpha = (Real*)malloc(M*sizeof(Real));
			#endif
			for (size_t i = 0; i < Sys[0]->SysMolMonList.size() ; i++) {
					Cp(temp_alpha, &xx[i*M] , M);
				for (int k=0; k<mon_length; k++) {
					chi = Sys[0]->CHI[Sys[0]->SysMolMonList[i]*mon_length+k];
					if (chi!=0)
						PutAlpha(temp_alpha,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
				}
				mesodyn_load_alpha(temp_alpha, i);
			}

			RHO = mesodyn_flux();

			#if defined(PAR_MESODYN) || ! defined(CUDA)
			Cp(g,RHO,iv);
			#else
			TransferDataToDevice(RHO, g, iv);
			#endif
			

			size_t k = 0;
			for (size_t i = 0 ; i < In[0]->MolList.size() ; ++i) {
				for (size_t a = 0 ; a < Mol[i]->MolMonList.size(); ++a) {
					target_function(g, k, M, i, a);
					Lat[0]->remove_bounds(g+k*M);
					Times(g+k*M,g+k*M,Sys[0]->KSAM,M);
					k++;
				}
			}

			#ifdef CUDA
			cudaFree(temp_alpha);
			#else
			free(temp_alpha);
			#endif
		}
		break;
		case custum:
			if (debug) cout <<"Residuals in custum mode in Solve_scf " << endl;
			if (value_ets==-1 && value_etm==-1 && value_bm==-1) 			//guess from newton is stored in place.
				Var[value_search]->PutValue(x[0]);
			else {
				if (value_ets>-1) Var[value_ets]->PutValue(x[0]);
				if (value_etm>-1) Var[value_etm]->PutValue(x[0]);
				if (value_bm>-1) Var[value_bm]->PutValue(x[0]);
			}

			if (value_ets==-1|| value_search<0 || value_etm>-1 || value_bm>-1) {
				control=proceed;					//prepare for solving scf eqns. (both for inneriteration and residuals)
				gradient=classical;
				if (SCF_method=="hessian") {hessian=true; pseudohessian=false; solver=HESSIAN;}
				if (SCF_method=="pseudohessian") {hessian=false; pseudohessian=true; solver=PSEUDOHESSIAN;}
				if (SCF_method=="DIIS") {solver=diis;}
				if (SCF_method=="Picard") {solver=PICARD;}
				//e_info=value_e_info;
				s_info=value_s_info;
				e_info=false;
				i_info=value_i_info;
				//tolerance = value_tolerance;
				Solve(false);						//find scf solution
				control=super;
				gradient=custum;
				tolerance=super_tolerance;
				e_info=super_e_info;
				i_info=super_i_info;
				s_info=super_s_info;
			} else {
				old_value_bm=value_bm;
				old_value_ets=value_ets;				//prepare of next level of super-iteration.
				old_value_etm=value_etm;
				SuperIterate(value_search,value_target,-1,-1,-1);	//go to superiteration with new ets and etm values
				value_ets = old_value_ets;				//reset conditions so that old iteration can continue
				value_etm=old_value_etm;
				value_bm=old_value_bm;
			}

			if (value_ets==-1 && value_etm==-1 && value_bm==-1) {			//get value for gradient.
				g[0]=Var[value_target]->GetError();
			} else {
				if (value_ets>-1) g[0]=Var[value_ets]->GetError();
				if (value_etm>-1) g[0]=Var[value_etm]->GetError();
				if (value_bm>-1) g[0]=Var[value_bm]->GetError();
			}
		break;
		case Picard:
		{
			if (debug) cout <<"Residuals in Picard mode in Solve_scf " << endl;
			int jump=sysmon_length;
			if (Sys[0]->charged) jump++;
			Cp(alpha,xx+jump*M,M);
			ComputePhis();
			if (Sys[0]->charged) {
				Sys[0]->DoElectrostatics(g+sysmon_length*M,xx+sysmon_length*M);
				Lat[0]->UpdateEE(Sys[0]->EE,Sys[0]->psi);
				Lat[0]->set_bounds(Sys[0]->psi);
				Lat[0]->UpdatePsi(g+sysmon_length*M,Sys[0]->psi,Sys[0]->q,Sys[0]->eps,Sys[0]->psiMask);
				Lat[0]->remove_bounds(g+sysmon_length*M);
			}
			YisAplusC(g+jump*M,Sys[0]->phitot,-1.0,M);
			for (i=0; i<sysmon_length; i++) {
				Cp(g+i*M,xx+i*M,M);
				for (k=0; k<mon_length; k++) {
                       		chi= -1.0*Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k];  //The minus sign here is to change the sign of x! just a trick due to properties of PutAlpha where a minus sing is implemented....

					if (chi!=0) PutAlpha(g+i*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
				}
				if (Sys[0]->charged){
					YplusisCtimesX(g+i*M,Sys[0]->EE,Seg[Sys[0]->SysMonList[i]]->epsilon,M);
					if (Seg[Sys[0]->SysMonList[i]]->valence !=0)
					YplusisCtimesX(g+i*M,Sys[0]->psi,-1.0*Seg[Sys[0]->SysMonList[i]]->valence,M);
				}
				Lat[0]->remove_bounds(g+i*M);
				Times(g+i*M,g+i*M,Sys[0]->KSAM,M);
			}
		break;
		}
		default:
			if (debug) cout <<"Residuals in scf mode in Solve_scf " << endl;
			int itmonlistlength=Sys[0]->ItMonList.size();
			int state_length = In[0]->StateList.size();
			int itstatelistlength=Sys[0]->ItStateList.size();
			
			if (Sys[0]->beta_set) {
				Sys[0]->BETA=xx[iv-1]; 
			}

 			ComputePhis();

			Cp(g,xx,iv); Zero(alpha,M);
			for (i=0; i<itmonlistlength; i++) {
				for (k=0; k<mon_length; k++) {
					if (Seg[k]->ns<2) {
						chi =Seg[Sys[0]->ItMonList[i]]->chi[k];
						if (chi!=0)
							PutAlpha(g+i*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
					}
				}
			 	for (k=0; k<state_length; k++) {
					chi =Seg[Sys[0]->ItMonList[i]]->chi[mon_length+k];
					if (chi!=0) {
						PutAlpha(g+i*M,Sys[0]->phitot,Seg[Sta[k]->mon_nr]->phi_side + Sta[k]->state_nr*M,chi,Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr],M);
					}
				}
			}
			for (i=0; i<itmonlistlength; i++) Add(alpha,g+i*M,M);
			for (i=0; i<itstatelistlength; i++) {
				for (k=0; k<mon_length; k++) {
					if (Seg[k]->ns<2) {
						chi =Sta[Sys[0]->ItStateList[i]]->chi[k];
						if (chi!=0)
//cout <<"for segment k " << k << " chi " << chi << endl;
							PutAlpha(g+(itmonlistlength+i)*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
					}
				}
				for (k=0; k<state_length; k++) {
					chi =Sta[Sys[0]->ItStateList[i]]->chi[mon_length+k];
					if (chi!=0)
						PutAlpha(g+(itmonlistlength+i)*M,Sys[0]->phitot,Seg[Sta[k]->mon_nr]->phi_side + Sta[k]->state_nr*M,chi,Seg[Sta[k]->mon_nr]->state_phibulk[Sta[k]->state_nr],M);

				}
			}
			for (i=0; i<itstatelistlength; i++) Add(alpha,g+(itmonlistlength+i)*M,M);
			Norm(alpha,1.0/(itmonlistlength+itstatelistlength),M);
			for (i=0; i<itmonlistlength; i++) {
				AddG(g+i*M,Sys[0]->phitot,alpha,M);
				Lat[0]->remove_bounds(g+i*M);
				Times(g+i*M,g+i*M,Sys[0]->KSAM,M);
			}
			for (i=0; i<itstatelistlength; i++) {
				AddG(g+(itmonlistlength+i)*M,Sys[0]->phitot,alpha,M);
				Lat[0]->remove_bounds(g+(itmonlistlength+i)*M);
				Times(g+(itmonlistlength+i)*M,g+(itmonlistlength+i)*M,Sys[0]->KSAM,M);
			}
			if (Sys[0]->charged) {
				Cp(g+(itmonlistlength+itstatelistlength)*M,xx+(itmonlistlength+itstatelistlength)*M,M);
				Sys[0]->DoElectrostatics(g+(itmonlistlength+itstatelistlength)*M,xx+(itmonlistlength+itstatelistlength)*M);
				Lat[0]->set_bounds(Sys[0]->psi);
				Lat[0]->UpdatePsi(g+(itmonlistlength+itstatelistlength)*M,Sys[0]->psi,Sys[0]->q,Sys[0]->eps,Sys[0]->psiMask);
				Lat[0]->remove_bounds(g+(itmonlistlength+itstatelistlength)*M);
			}
			if (Sys[0]->beta_set){
			   g[iv-1]=(Mol[Sys[0]->MolBeta]->phitot[Mol[Sys[0]->MolBeta]->beta]-Mol[Sys[0]->solvent]->phitot[Mol[Sys[0]->solvent]->beta]); 
			}
		break;
	}
}

void Solve_scf::gradient_log(Real* g, int k, int M, int i, int j) {
	//Target function: ln(g/phi) < tolerance
	for (int z = 0 ; z < M ; ++z) {
		Real frac = (g+k*M)[z]/(Mol[i]->phi+j*M)[z];
		(g+k*M)[z] = log( frac );
	}
}

void Solve_scf::gradient_quotient(Real* g, int k, int M, int i, int j) {
	//Target function: g/phi-1 < tolerance
	for (int z = 0 ; z < M ; ++z) {
		Real frac = (g+k*M)[z]/(Mol[i]->phi+j*M)[z];
		(g+k*M)[z] = frac - 1;
	}
}

void Solve_scf::gradient_minus(Real* g, int k, int M, int i, int j) {
	//Target function: g - phi < tolerance
		YplusisCtimesX(g+k*M,Mol[i]->phi+j*M,-1.0,M);
}

void Solve_scf::inneriteration(Real* x, Real* g, float* h, Real accuracy, Real& deltamax, Real ALPHA, int nvar) {
if(debug) cout <<"inneriteration in Solve_scf " << endl;
	residual=accuracy; //hoping this is not creating problems with the use of residual...
	switch(control) {
		case super:
			for (int i=0; i<nvar; i++) { //this is to control the sing in the WEAK iteration.
				if (h[i+i*nvar]<0) {SIGN[i] = -1;  }
			}

		break;
		default:
			if (iterations > 0) samehessian = false;
			if (reset_pseudohessian) {reset_pseudohessian=false; pseudohessian = true;}

			if (accuracy < minAccuracySoFar && iterations > 0 && accuracy == fabs(accuracy) ) {
				minAccuracySoFar = accuracy;
			}

			if (accuracy > minAccuracySoFar*resetHessianCriterion && accuracy == fabs(accuracy) ) {
				if (s_info) {
					cout << accuracy << '\t' << minAccuracySoFar << '\t' << resetHessianCriterion << endl;
					cout << "walking backwards: newton reset" << endl;
				}
				resethessian(h,g,x,nvar);
				minAccuracySoFar *=1.5;

				if (deltamax >0.005) deltamax *=0.9;
				//if (!reverseDirection) {reverseDirection = new int[reverseDirectionRange]; H_Zero(reverseDirection,reverseDirectionRange); }
				numIterationsSinceHessian = 0;
			}

			if (ALPHA < smallAlpha) smallAlphaCount++; else smallAlphaCount = 0;

			if (smallAlphaCount == maxNumSmallAlpha) {
				smallAlphaCount = 0;
				//if (!reverseDirection) {reverseDirection=new Real[reverseDirectionRange];H_Zero(reverseDirection,reverseDirectionRange); }
				if (!s_info) {
					cout << "too many small alphas: newton reset" << endl;
				}
				resethessian(h,g,x,nvar);
				if (deltamax >0.005) deltamax *=0.9;
				numIterationsSinceHessian = 0;
			}

			if (!newtondirection && pseudohessian) {
				reverseDirection[iterations%reverseDirectionRange] = 1;
			} else {
				reverseDirection[iterations%reverseDirectionRange] = 0;
			}

			numReverseDirection = 0;
			for (int i=0; i<reverseDirectionRange; i++) {
				if (reverseDirection[i] == 1)
					numReverseDirection++;
			}

			numIterationsSinceHessian++;
			Real frReverseDirection = Real(numReverseDirection)/reverseDirectionRange;
			if ((frReverseDirection > maxFrReverseDirection && pseudohessian && accuracy < minAccuracyForHessian)) {
				cout <<"Bad convergence (reverse direction), computing full hessian..." << endl;
				pseudohessian = false; reset_pseudohessian =true;
				//if (!reverseDirection) {reverseDirection=new Real[reverseDirectionRange]; H_Zero(reverseDirection,reverseDirectionRange); }
				numIterationsSinceHessian = 0;
			} else if ((numIterationsSinceHessian >= n_iterations_for_hessian &&
						iterations > 0 && accuracy < minAccuracyForHessian && minimum < minAccuracyForHessian)) {
				cout << "Still no solution, computing full hessian..." << endl;
				pseudohessian = false; reset_pseudohessian =true;
				numIterationsSinceHessian = 0;
			}

		break;
	}
}

void Solve_scf::ComputePhis() {
if(debug) cout <<"ComputPhis in  Solve_scf " << endl;
	PutU();
	Sys[0]->PrepareForCalculations();
	Sys[0]->ComputePhis();
}

bool Solve_scf::PutU() {
if(debug) cout <<"PutU in  Solve " << endl;
	int M=Lat[0]->M;
	Real valence;
	Real *u;
	if (SCF_method == "Picard") {cout << " Picard not implemented properly " << endl; }
	bool success=true;
	alpha=Sys[0]->alpha;

	if (Sys[0]->charged) {
		Cp(Sys[0]->psi,xx+iv-M,M);
		Lat[0]->UpdateEE(Sys[0]->EE,Sys[0]->psi);
	}

	int itmonlistlength=Sys[0]->ItMonList.size();
	int itstatelistlength=Sys[0]->ItStateList.size();
	int monlistlength =In[0]->MonList.size();
	int statelistlength=In[0]->StateList.size();
	int k=0;
	for (int i=0; i<itmonlistlength; i++) {
		int IM=Sys[0]->ItMonList[i];
		u=Seg[IM]->u;
		Cp(u,xx+k*M,M);
		if (Sys[0]->charged){
			YplusisCtimesX(u,Sys[0]->EE,-1.0*Seg[IM]->epsilon,M);
			valence=Seg[IM]->valence;
			if (valence !=0)
				YplusisCtimesX(u,Sys[0]->psi,valence,M);
		}
		for (int j=0; j<monlistlength; j++) {
			if (Seg[j]->seg_nr_of_copy==IM && Seg[j]->ns<2) {
				u=Seg[j]->u;
				Cp(u,xx+k*M,M);
				if (Sys[0]->charged){
					YplusisCtimesX(u,Sys[0]->EE,-1.0*Seg[j]->epsilon,M);
					valence=Seg[j]->valence;
					if (valence !=0)
						YplusisCtimesX(u,Sys[0]->psi,valence,M);
				}
			}

		}
		for (int j=0; j<statelistlength; j++) {
			if (Sta[j]->seg_nr_of_copy==IM) {
				u=Seg[Sta[j]->mon_nr]->u+Sta[j]->state_nr*M;
				Cp(u,xx+k*M,M);
				if (Sys[0]->charged){
					YplusisCtimesX(u,Sys[0]->EE,-1.0*Seg[Sta[j]->mon_nr]->epsilon,M);
					valence=Sta[j]->valence;
					if (valence !=0)
						YplusisCtimesX(u,Sys[0]->psi,valence,M);

				}
			}
		}
		k++;
	}

	for (int i=0; i<itstatelistlength; i++) {
		int IS=Sys[0]->ItStateList[i];
		u=Seg[Sta[IS]->mon_nr]->u+(Sta[IS]->state_nr)*M;
		Cp(u,xx+k*M,M);
		if (Sys[0]->charged){
			YplusisCtimesX(u,Sys[0]->EE,-1.0*Seg[Sta[IS]->mon_nr]->epsilon,M);
			valence=Sta[IS]->valence;
			if (valence !=0)
				YplusisCtimesX(u,Sys[0]->psi,valence,M);
		}
		for (int j=0; j<statelistlength; j++) {
			if (Sta[j]->state_nr_of_copy==IS) {
				u=Seg[Sta[j]->mon_nr]->u+Sta[j]->state_nr*M;
				Cp(u,xx+k*M,M);
				if (Sys[0]->charged){
					YplusisCtimesX(u,Sys[0]->EE,-1.0*Seg[Sta[j]->mon_nr]->epsilon,M);
					valence=Sta[j]->valence;
					if (valence !=0)
						YplusisCtimesX(u,Sys[0]->psi,valence,M);
				}
			}
		}
		k++;
	}
	return success;
}

/*---------------------------------------------to be saved for a while--------------------

				if (Sys[0]->charged){
					YplusisCtimesX(g+i*M,Sys[0]->EE,Seg[Sys[0]->ItMonList[i]]->epsilon,M);
					valence=Seg[Sys[0]->ItMonList[i]]->valence;
					if (valence !=0)
						YplusisCtimesX(g+i*M,Sys[0]->psi,-1.0*valence,M);
				}



	if (In[0]->MesodynList.size()>0) {

		int i=0; int k=0;
		int length = In[0]->MolList.size();
		while (i<length) {
			int j=0;
			int LENGTH=Mol[i]->MolMonList.size();
			while (j<LENGTH) {
				Cp(Mol[i]->u+j*M,xx+k*M,M);
				k++; j++;
			}
			i++;
		}

		cout <<"Daniel: PutU in sf_solve is modified for Mesodyn. contact frans in case of problems."<<endl;
//This code must be modified in case of mon's with internal states.
		int k=0;
		int length = In[0]->MolList.size();
		for (int i=0; i<length; i++) {
			int LENGTH=Mol[i]->MolMonList.size();
			for (int j=0; j<LENGTH; j++) {
				Cp(Seg[Mol[i]->MolMonList[j]]->u,xx+k*M,M);
				Seg[Mol[i]->MolMonList[j]]->DoBoltzmann();
				k++;
			}
			Mol[i]->CpBoltzmann();
		}
	} else {
*/
