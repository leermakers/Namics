#include "newttool.h"
#include "newton.h"  

Newton::Newton(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,vector<System*> Sys_,string name_) {
	In=In_; name=name_; Sys=Sys_; Seg=Seg_; Lat=Lat_; Mol=Mol_;
if(debug) cout <<"Constructor in Newton " << endl;  
	KEYS.push_back("method"); KEYS.push_back("e_info"); KEYS.push_back("s_info");
	KEYS.push_back("delta_max"); KEYS.push_back("m"); KEYS.push_back("i_info");
	KEYS.push_back("iterationlimit" ); KEYS.push_back("tolerance"); KEYS.push_back("store_guess"); KEYS.push_back("read_guess"); 
	KEYS.push_back("stop_criterion"); 
	KEYS.push_back("delta_min");
	KEYS.push_back("hessian_width"); 
	KEYS.push_back("linesearchlimit");
	//KEYS.push_back("samehessian");  
	KEYS.push_back("max_accuracy_for_hessian_scaling");
	KEYS.push_back("n_iterations_for_hessian");
	KEYS.push_back("small_alpha");    
	KEYS.push_back("max_n_small_alpha");
	KEYS.push_back("min_accuracy_for_hessian");
	KEYS.push_back("max_fr_reverse_direction");
}

Newton::~Newton() {
if(debug) cout <<"Destructor in Newton " << endl; 
	free(Aij);
	free(Ci);
	free(Apij); 
	if (h) delete [] h;
#ifdef CUDA
	cudaFree(xx);
	cudaFree(x0);
	cudaFree(g);
	cudaFree(xR);
	cudaFree(x_x0);
#else
	free(xx);
	free(x0);
	free(g);
	free(xR);
	free(x_x0);
#endif
}

void Newton::AllocateMemory() {
if(debug) cout <<"AllocateMemeory in Newton " << endl; 
	int M=Lat[0]->M;
	iv = Sys[0]->SysMonList.size() * M;	
	if (method=="DIIS-ext") iv +=M;
	Aij =(double*) malloc(m*m*sizeof(double)); H_Zero(Aij,m*m);
	Ci =(double*) malloc(m*sizeof(double)); H_Zero(Ci,m);
	Apij =(double*) malloc(m*m*sizeof(double)); H_Zero(Apij,m*m);
#ifdef CUDA
	xx = (double*)AllOnDev(iv); 
	x0 = (double*)AllOnDev(iv);
	g= (double*)AllOnDev(iv);
	xR= (double*)AllOnDev(m*iv);
	x_x0= (double*)AllOnDev(m*iv);
#else
	xx =(double*) malloc(iv*sizeof(double)); 
	x0 =(double*) malloc(iv*sizeof(double)); 
	g =(double*) malloc(iv*sizeof(double)); 
	xR =(double*) malloc(m*iv*sizeof(double)); 
	x_x0 =(double*) malloc(m*iv*sizeof(double)); 
#endif
	Zero(xx,iv);
	Zero(x0,iv);
	Zero(g,iv);
	Zero(xR,m*iv);
	Zero(x_x0,m*iv);
if (debug){
	double test; Sum(test,xx,iv); cout <<"Sum of xx after putting to zero" << test << endl; 
}
	Sys[0]->AllocateMemory();
}

bool Newton::PrepareForCalculations() {
if (debug) cout <<"PrepareForCalculations in Newton " << endl; 
	bool success=true;
	return success; 
}

bool Newton::CheckInput(int start_) { start=start_;
if(debug) cout <<"CheckInput in Newton " << endl; 
	bool success=true;
	string value;
	success=In[0]->CheckParameters("newton",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		e_info=In[0]->Get_bool(GetValue("e_info"),true); 
		s_info=In[0]->Get_bool(GetValue("s_info"),false); 
		i_info=In[0]->Get_int(GetValue("i_info"),1);
		iterationlimit=In[0]->Get_int(GetValue("iterationlimit"),1000);
		if (iterationlimit < 0 || iterationlimit>1e6) {iterationlimit = 1000; 
		cout << "Value of 'iterationlimit' out of range 0..1e6, and value set to default value 1000" <<endl; } 

		delta_max=In[0]->Get_double(GetValue("delta_max"),0.1); 
		if (delta_max < 0 || delta_max>100) {delta_max = 0.1;  cout << "Value of delta_max out of range 0..100, and value set to default value 0.1" <<endl; } 
		tolerance=In[0]->Get_double(GetValue("tolerance"),1e-5);
		if (tolerance < 1e-12 ||tolerance>10) {tolerance = 1e-5;  cout << "Value of tolerance out of range 1e-12..10 Value set to default value 1e-5" <<endl; }  
		if (GetValue("method").size()==0) {method="DIIS";} else {
			vector<string>method_options; 
			method_options.push_back("DIIS");
			method_options.push_back("DIIS-ext"); 
			method_options.push_back("Picard");
			method_options.push_back("pseudohessian"); 
			if (!In[0]->Get_string(GetValue("method"),method,method_options,"In 'newton' the entry for 'method' not recognized: choose from:")) success=false;
		}
		m=In[0]->Get_int(GetValue("m"),10); 
		if (method=="pseudohessian") {
			pseudohessian=true;
			hessianwidth=In[0]->Get_int(GetValue("hessianwidth"),0);	
			if (hessianwidth < 0) {cout <<"hessianwidth should have a positive value " << endl; hessianwidth=0;}
			int niv=Sys[0]->SysMonList.size()*Lat[0]->M;
			if ( (hessianwidth*2-1) > niv ) {
				cout << "hessianwidth is limited to the number of iteration variables: hessianwidth value is ignored. " << endl;
				hessianwidth=0; 
			}
			samehessian=false; //In[0]->Get_bool(GetValue("samehessian"),false);
			linesearchlimit=In[0]->Get_int(GetValue("hessianwidth"),20);	
			if (linesearchlimit<0 || linesearchlimit >100) {
				cout <<"linesearchlimit is out of range: 1..100; default value of 20 is used instead" << endl;
				linesearchlimit=20;
			}
			max_accuracy_for_hessian_scaling=In[0]->Get_double(GetValue("max_accuracy_for_hessian_scaling"),0.1);
			if (max_accuracy_for_hessian_scaling<1e-7 || max_accuracy_for_hessian_scaling>1) {
				cout <<"max_accuracy_for_hessian_scaling is out of range: 1e-7...1; default value 0.1 is used instead" << endl; 
				max_accuracy_for_hessian_scaling=0.1; 
			}
			minAccuracyForHessian=In[0]->Get_double(GetValue("min_accuracy_for_hessian"),0);
			if (minAccuracyForHessian<0 ||minAccuracyForHessian>0.1) {
				cout <<"min_accuracy_for_hessian is out of range: 0...0.1; default value 0 is used instead (no hessian computation)" << endl; 
				minAccuracyForHessian=0; 
			}
			maxFrReverseDirection =In[0]->Get_double(GetValue("max_fr_reverse_direction"),0.4);
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

			delta_min=In[0]->Get_double(GetValue("delta_min"),0);
			if (delta_min <0 || delta_min>delta_max) {
				cout <<"delta_min is out of range; 0, ..., " << delta_max << "; delta_min value set to 0 " << endl; 
				delta_min=0; 
			}
			smallAlpha=In[0]->Get_double(GetValue("small_alpha"),0.00001);
			if (smallAlpha <0 || smallAlpha>1) {
				cout <<"small_alpha is out of range; 0, ..., 1; small_alpha value set to default: 1e-5 " << endl; 
				smallAlpha=0.00001; 
			}


		}		
		if (m < 0 ||m>100) {m=10;  cout << "Value of 'm' out of range 0..100, value set to default value 10" <<endl; }
		StoreFileGuess=In[0]->Get_string(GetValue("store_guess"),"");		
		ReadFileGuess=In[0]->Get_string(GetValue("read_guess"),"");
		if (GetValue("stop_criterion").size() > 0) {
			vector<string>options;
			options.push_back("norm_of_g");
			options.push_back("max_of_element_of_|g|");
			if (!In[0]->Get_string(GetValue("stop_criterion"),stop_criterion,options,"In newton the stop_criterion setting was not recognised")) {success=false; };
		}
	}
	return success; 
}


void Newton::PutParameter(string new_param) {
if(debug) cout <<"PutParameter in Newton " << endl; 
	KEYS.push_back(new_param); 
}

string Newton::GetValue(string parameter){
if(debug) cout <<"GetValue " + parameter + " in  Newton " << endl;
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

void Newton::push(string s, double X) {
if(debug) cout <<"push (double) in  Newton " << endl;
	doubles.push_back(s);
	doubles_value.push_back(X); 
}
void Newton::push(string s, int X) {
if(debug) cout <<"push (int) in  Newton " << endl;
	ints.push_back(s);
	ints_value.push_back(X); 
}
void Newton::push(string s, bool X) {
if(debug) cout <<"push (bool) in  Newton " << endl;
	bools.push_back(s);
	bools_value.push_back(X); 
}
void Newton::push(string s, string X) {
if(debug) cout <<"push (string) in  Newton " << endl;
	strings.push_back(s);
	strings_value.push_back(X); 	
}
void Newton::PushOutput() {
if(debug) cout <<"PushOutput in  Newton " << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	doubles.clear();
	doubles_value.clear();
	ints.clear();
	ints_value.clear();  
	push("method",method);
	push("m",m);
	push("delta_max",delta_max);
	push("residual",residual);
	push("tolerance",tolerance);
	push("iterations",it);
	push("iterationlimit",iterationlimit);
	push("stop_criterion",stop_criterion);
	if (pseudohessian) {
		push("hessian_width",hessianwidth);
		//push("same_hessian",samehessian);
		push("linesearchlimit",linesearchlimit);
		push("max_accuracy_for_hessian_scaling",max_accuracy_for_hessian_scaling);
		push("n_iteratons_for_hessian",n_iterations_for_hessian);
		push("small_alpha",smallAlpha);
		push("max_n_small_alpha",maxNumSmallAlpha); 
		push("min_accuracy_for_hessian",minAccuracyForHessian); 
	}
}

int Newton::GetValue(string prop,int &int_result,double &double_result,string &string_result){
if(debug) cout <<"GetValue (long) in  Newton " << endl;
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
	length = doubles.size();
	while (i<length) {
		if (prop==doubles[i]) { 
			double_result=doubles_value[i];
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

void Newton::inneriteration(double*h, double *g, double *x,double accuracy, int nvar, int m) {

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
		resethessian(h,g,x,nvar,m);
		minAccuracySoFar *=1.5;

		if (delta_max >0.005) delta_max *=0.9;
		if (!reverseDirection) {reverseDirection = new double[reverseDirectionRange]; H_Zero(reverseDirection,reverseDirectionRange); } 
		numIterationsSinceHessian = 0;
	}

	if (ALPHA < smallAlpha) smallAlphaCount++; else smallAlphaCount = 0;

	if (smallAlphaCount == maxNumSmallAlpha) {
		smallAlphaCount = 0;
		if (!reverseDirection) {reverseDirection=new double[reverseDirectionRange];H_Zero(reverseDirection,reverseDirectionRange); }
		if (!s_info) {
			cout << "too many small alphas: newton reset" << endl;
		}
		resethessian(h,g,x,nvar,m);
		if (delta_max >0.005) delta_max *=0.9;
		numIterationsSinceHessian = 0;
	}

	if (!newtondirection && pseudohessian) {
		reverseDirection[iterations%reverseDirectionRange + 1] = 1;
	} else {
		reverseDirection[iterations%reverseDirectionRange + 1] = 0;
	}

	numReverseDirection = 0;
	for (int i=1; i<=reverseDirectionRange; i++) {
		if (reverseDirection[i] > 0.5) numReverseDirection++;
	}

	numIterationsSinceHessian++;
	double frReverseDirection = double(numReverseDirection)/reverseDirectionRange;
	if ((frReverseDirection > maxFrReverseDirection && pseudohessian && accuracy < minAccuracyForHessian)) {
		cout <<"Bad convergence (reverse direction), computing full hessian..." << endl;
		pseudohessian = false; reset_pseudohessian =true;
		if (!reverseDirection) {reverseDirection=new double[reverseDirectionRange]; H_Zero(reverseDirection,reverseDirectionRange); }
		numIterationsSinceHessian = 0;
	} else if ((numIterationsSinceHessian >= n_iterations_for_hessian &&
				iterations > 0 && accuracy < minAccuracyForHessian && minimum < minAccuracyForHessian)) {
		cout << "Still no solution, computing full hessian..." << endl;
		pseudohessian = false; reset_pseudohessian =true;
		numIterationsSinceHessian = 0;
	}
}

void Newton::direction(double *h, double *p, double *g, double *g0, double *x, int nvar, int m, double alpha){
	newtondirection = true;
	gausa(h,p,g,nvar,m);
	gausb(h,p,nvar,m); 
	if (ignore_newton_direction) {
		newtondirection = true;
	} else {
		newtondirection = signdeterminant(h,nvar,m)>0;
	}
	if ( !newtondirection ) {
		for (int i=0; i<nvar; i++) {
			p[i] *= -1;
		}
		if ( e_info) cout << "*";
	}
}

void Newton::direction(double *h, double *p, double *g, double *g0, double *x, int nvar, double alpha){

	newtondirection = true;
	gausa(h,p,g,nvar); 
	gausb(h,p,nvar); 
	if (ignore_newton_direction) {
		newtondirection = true;
	} else {
		newtondirection = signdeterminant(h,nvar)>0;
	}
	if ( !newtondirection ) {
		for (int i=0; i<nvar; i++) {
			p[i] *= -1;
		}
		if ( e_info) cout << "*";
	}
}

void Newton::startderivatives(double *h, double *g, double *x, int nvar, int m){
	double diagonal = 1+norm2(g,nvar); 
	H_Zero(h,nvar*(2*m-1));
	double* ha=&h[-1];
	double* hai; int i0; 
	for (int i=1; i<=nvar; i++) {
		hai=&ha[(i-1)*nvar];
		if (i<m) i0=i; else i0=m;
		hai[i0] = diagonal; 
	}
}

void Newton::startderivatives(double *h, double *g, double *x, int nvar){
	double diagonal = 1+norm2(g,nvar); 
	H_Zero(h,nvar*nvar);
	for (int i=0; i<nvar; i++) {
		h[i+nvar*i] = diagonal;  
	}
}

void Newton::resethessian(double *h,double* g,double* x,int nvar, int m){
	trouble = 0;		
	startderivatives(h,g,x,nvar,m);
	resetiteration = iterations;
}

void Newton::resethessian(double *h,double *g,double *x,int nvar){
	trouble = 0;		
	startderivatives(h,g,x,nvar);
	resetiteration = iterations;
}

void Newton::newhessian(double *h, double *g, double *g0, double *x, double *p, int nvar, int m) { 
	double dmin=0,sum=0,theta=0,php=0,dg=0,gg=0,g2=0,py=0,y2=0;
	dmin = 1/pow(2.0,nbits); // alternative: DBL_EPSILON or DBL_MIN
	double *ha = &h[-1];
	double *pa = &p[-1]; 
	int i0,nn;
	double *hai; 

	if (!pseudohessian) findhessian(h,g,x,nvar,m); 
	else {
		if (!samehessian && ALPHA!=0 && iterations!=0) {
			double *y = new double[nvar];
			double *hp = new double[nvar];
			py = php = y2 = gg = g2 = 0;
			for (int i=0; i<nvar; i++) { //because i runs from 0 unit nvar-1 we do not need ga, ya
				y[i] = dg = g[i]-g0[i];
				py += p[i]*dg;
				y2 += pow(y[i],2);
				gg += g[i]*g0[i];
				g2 += pow(g[i],2);
			}
	
			if ( !newtondirection ) {
				multiply(hp,1,h,p,nvar,m); 
			} else {
				for (int i=0; i<nvar; i++) hp[i] = -g0[i];
			}

			Dot(php,p,hp,nvar);
			theta = py/(10*dmin+ALPHA*php);
			if ( nvar>=1 && theta>0 && iterations==resetiteration+1 && accuracy > max_accuracy_for_hessian_scaling) {
				if (e_info ) {
					cout << "hessian scaling: " << theta << endl;
				}
				ALPHA *= theta;
				py /= theta;
				php /= theta;
				for (int i=1; i<=nvar; i++) {
					if (i>m) i0=i-m; else i0=0;
					nn=m+i-1; if (nn>nvar) nn=nvar; 
					hai=&ha[(i-1)*nvar];
					pa[i] /= theta;
					for (int j=i0+1; j<=nn; j++)
					hai[j-i0] *= theta; 				
				}
			}
			if (nvar>=1) trustfactor *= (4/(pow(theta-1,2)+1)+0.5);
			if ( nvar>1 ) {
				sum = ALPHA*pow(norm2(p,nvar),2);
				theta = fabs(py/(ALPHA*php));
				if ( theta<.01 ) sum /= 0.8;
				else if ( theta>100 ) sum *= theta/50;
				for (int i=0; i<nvar; i++) y[i] -= ALPHA*hp[i]; //because loop goes from i=0; classical c++
				updatpos(h,y,p,nvar,m,1.0/sum); 
				trouble -= signdeterminant(h,nvar,m); 
				if ( trouble<0 ) trouble = 0; else if ( trouble>=3 ) {
					resethessian(h,g,x,nvar,m);				}
			} else if ( nvar>=1 && py>0 ) {
				trouble = 0;
				theta = py>0.2*ALPHA*php ? 1 : 0.8*ALPHA*php/(ALPHA*php-py);
				if ( theta<1 ) {
					py = 0;
					for (int i=0; i<nvar; i++) {
						y[i] = theta*y[i]+(1-theta)*ALPHA*hp[i];
						py += p[i]*y[i];
					}
				}
				updatpos(h,y,y,nvar,m,1.0/(ALPHA*py)); 
				updateneg(h,hp,nvar,m,-1.0/php);	 
			}

			delete [] y;
			delete [] hp;
		} else if ( !samehessian ) resethessian(h,g,x,nvar,m);
	}
}


void Newton::newhessian(double *h, double *g, double *g0, double *x, double *p, int nvar) { 
	double dmin=0,sum=0,theta=0,php=0,dg=0,gg=0,g2=0,py=0,y2=0;
	dmin = 1/pow(2.0,nbits); // alternative: DBL_EPSILON or DBL_MIN


	if (!pseudohessian){ 
		findhessian(h,g,x,nvar); 
	} else {
		if (!samehessian && ALPHA!=0 && iterations!=0) {
			double *y = new double[nvar];
			double *hp = new double[nvar];
			py = php = y2 = gg = g2 = 0;
			for (int i=0; i<nvar; i++) {
				y[i] = dg = g[i]-g0[i];
				py += p[i]*dg;
				y2 += pow(y[i],2);
				gg += g[i]*g0[i];
				g2 += pow(g[i],2);
			}
	
			if ( !newtondirection ) {
				multiply(hp,1,h,p,nvar); 
			} else {
				for (int i=0; i<nvar; i++) hp[i] = -g0[i];
			}

			Dot(php,p,hp,nvar);
			theta = py/(10*dmin+ALPHA*php);
			if ( nvar>=1 && theta>0 && iterations==resetiteration+1 && accuracy > max_accuracy_for_hessian_scaling) {
				if (e_info ) {
					cout << "hessian scaling: " << theta << endl;
				}
				ALPHA *= theta;
				py /= theta;
				php /= theta;
				for (int i=0; i<nvar; i++) {
					p[i] /= theta;
					h[i+nvar*i] *= theta; 
				}
			}
			if (nvar>=1) trustfactor *= (4/(pow(theta-1,2)+1)+0.5);
			if ( nvar>1 ) {
				sum = ALPHA*pow(norm2(p,nvar),2);
				theta = fabs(py/(ALPHA*php));
				if ( theta<.01 ) sum /= 0.8;
				else if ( theta>100 ) sum *= theta/50;
				for (int i=0; i<nvar; i++) y[i] -= ALPHA*hp[i];

				updatpos(h,y,p,nvar,1.0/sum); 
				trouble -= signdeterminant(h,nvar);  
				if ( trouble<0 ) trouble = 0; else if ( trouble>=3 ) {
					resethessian(h,g,x,nvar);
				}
			} else if ( nvar>=1 && py>0 ) {
				trouble = 0;
				theta = py>0.2*ALPHA*php ? 1
				: 0.8*ALPHA*php/(ALPHA*php-py);
				if ( theta<1 ) {
					py = 0;
					for (int i=0; i<nvar; i++) {
						y[i] = theta*y[i]+(1-theta)*ALPHA*hp[i];
						py += p[i]*y[i];
					}
				}
				updatpos(h,y,y,nvar,1.0/(ALPHA*py)); 
				updateneg(h,hp,nvar,-1.0/php);	 
			}

			delete [] y;
			delete [] hp;
		} else if ( !samehessian ) resethessian(h,g,x,nvar);
	}
}

void Newton::numhessian(double *h,double *g, double *x,int nvar, int m) {
	double dmax2=0,dmax3=0,di=0, hji;
	int i0,j0, nn;
	double *hai,*haj;
	double* g1 = new double[nvar];
	double* d = new double[nvar];
	double* x0 = new double[nvar];
	double* ga = &g[-1];
	double* ha = &h[-1];
	double* g1a = &g[-1];
	double* da = &d[-1];
	double* xa = &x[-1];
	double* x0a= &x0[-1];
	dmax2 = pow(2.0,nbits/2); //alternative 2*pow(DBL_EPSILON,-0.5)?
	dmax3 = pow(2.0,nbits/3); //alternative 2*pow(DBL_EPSILON,-1.0/3)?

	int mm1=2*m-1; if (mm1>nvar) mm1=nvar; 
	for (int ia=1; ia<=mm1; ia +=1) {
		for (int i=ia; i <=nvar; i += mm1) {
			if (i>m) i0=i-m; else i0=0;
			x0a[i]=xa[i];
			hai=&h[(i-1)*nvar];
			di = (1/(dmax3*dmax3*fabs(hai[i-i0])+dmax3+fabs(ga[i])) +1/dmax2)*(1+fabs(xa[i]));
				if (di<delta_min) di=delta_min;
				xa[i]=xa[i]+di; da[i]=di;
		}
		ComputeG(g1); 
		for (int i=ia; i<=nvar; i+=mm1) {
			if (i>m) i0=i-m; else i0=0; 
			nn=m+i-1; if (nn>nvar) nn=nvar;
			xa[i]=x0[i];
			for (int j=i0+1; j<=nn; j++) {
				if (j>m) j0=j-m; else j0=0;
				hji = (g1a[i]-ga[i])/da[i];
				
				haj=&ha[(j-1)*nvar];
				haj[i-j0] = hji; 
			}
		}
	}
	delete [] g1;
	delete [] d;
	delete [] x0; 
	ComputeG(g);
}

void Newton::numhessian(double* h,double* g, double* x,int nvar) {
	double dmax2=0,dmax3=0,di=0;
	double *g1, *xt;
	g1 = new double[nvar];
	xt = new double[nvar];
	dmax2 = pow(2.0,nbits/2); //alternative 2*pow(DBL_EPSILON,-0.5)?
	dmax3 = pow(2.0,nbits/3); //alternative 2*pow(DBL_EPSILON,-1.0/3)?
	for (int i=0; i<nvar; i++) {
		xt[i] = x[i];
		di = (1/(dmax3*dmax3*fabs(h[i+nvar*i])+dmax3+fabs(g[i]))  
			+1/dmax2)*(1+fabs(x[i]));
		if ( di<delta_min ) {
			di = delta_min;
		}
		x[i] += di;
		ComputeG(g1);
		x[i] = xt[i];
		for (int j=0; j<nvar; j++ ) {
			h[j+nvar*i] = (g1[j]-g[j])/di; //aanpassen voor m
		}
	}
	delete [] g1;
	delete [] xt;
	ComputeG(g);
}


void Newton::decomposition(double *h, int nvar, int m, int &trouble){
	int ntr=0;
	decompos(h,nvar,m,ntr); 
	if (e_info) {
		if (iterations==0) {
			if (ntr==0) {
				cout << "Not too bad.";
			} else {
				cout << " sign might be wrong.";
			}
		} else if (ntr>0 && trouble==0) {
			if (ntr==1) {
				cout << "Wait a sec.";
			} else {
				cout << "some TROUBLES appear.";
			}
		} else if (ntr>trouble+1) {
			for (int i=1; i<= ntr; i++) {
				cout << "O";
			}
			cout << "H!";
		} else if (trouble>0) {
			if (ntr==0) {
				cout << "Here we go.";
			} else if (ntr<trouble-1) {
				cout << "Hold on.";
			} else if (ntr < trouble) {
				cout << "There is some hope for you.";
			} else if (ntr == trouble) {
				cout << "no Progress.";
			} else if (ntr > 4) {
				cout << "We won't fix it.";
			} else {
				cout << "More troubles.";
			}
		}
		if (iterations==0 || trouble>0 || ntr>0) {
			cout <<  endl;
		}
	}
	trouble = ntr;

}
void Newton::decomposition(double *h,int nvar, int &trouble){
	int ntr=0;
	decompos(h,nvar,ntr);

	if (e_info) {
		if (iterations==0) {
			if (ntr==0) {
				cout << "Not too bad.";
			} else {
				cout << " sign might be wrong.";
			}
		} else if (ntr>0 && trouble==0) {
			if (ntr==1) {
				cout << "Wait a sec.";
			} else {
				cout << "some TROUBLES appear.";
			}
		} else if (ntr>trouble+1) {
			for (int i=1; i<= ntr; i++) {
				cout << "O";
			}
			cout << "H!";
		} else if (trouble>0) {
			if (ntr==0) {
				cout << "Here we go.";
			} else if (ntr<trouble-1) {
				cout << "Hold on.";
			} else if (ntr < trouble) {
				cout << "There is some hope for you.";
			} else if (ntr == trouble) {
				cout << "no Progress.";
			} else if (ntr > 4) {
				cout << "We won't fix it.";
			} else {
				cout << "More troubles.";
			}
		}
		if (iterations==0 || trouble>0 || ntr>0) {
			cout <<  endl;
		}
	}
	trouble = ntr;

}

void Newton::findhessian(double *h, double *g, double *x, int nvar, int m) {
	if ( !samehessian ) {
		if ( iterations==0 ) resethessian(h,g,x,nvar,m);
		numhessian(h,g,x,nvar,m); // passes through residuals so check pseudohessian
		if (!pseudohessian) decomposition(h,nvar,m,trouble);
	}
}

void Newton::findhessian(double *h, double *g, double *x,int nvar) {
	if ( !samehessian ) {
		if ( iterations==0 ) resethessian(h,g,x,nvar); 
		numhessian(h,g,x,nvar); // passes through residuals so check pseudohessian
		if (!pseudohessian) decomposition(h,nvar,trouble);
	}
}

void Newton::newdirection(double *h, double *p, double *p0, double *g, double *g0, double *x, int nvar, int m, double alphabound) {
	inneriteration(h,g,x,accuracy,nvar,m);
	memcpy(p0, p, sizeof(*p0)*nvar);
	if (m>0) direction(h,p,g,g0,x,nvar,m, ALPHA); else direction(h,p,g,g0,x,nvar,ALPHA);
	accuracy = residue(g,p,x,nvar,ALPHA);
}

void Newton::newtrustregion(double *g, double *g0, double *p, double *p0, int nvar){
	double normp0 = norm2(p0,nvar);

	if ( normp0>0 && trustregion>2*ALPHA*normp0 ) {
		trustregion = 2*ALPHA*normp0;
	}
	trustregion *= trustfactor;
	trustfactor = 1.0;
	if ( trustregion>delta_max ) trustregion = delta_max;
	if ( trustregion<delta_min ) trustregion = delta_min;
}

double Newton::linesearch(double *g, double *g0, double *p, double *x, double *x0, int nvar, int m, double alphabound) {
	double newalpha = alphabound<1 ? alphabound : 1;
	newalpha = zero(g,g0,p,x,x0,nvar,m,newalpha);
	return newalpha;
}

double Newton::zero(double *g, double *g0, double *p, double *x, double *x0, int nvar, int m, double newalpha) {
	double alpha=newalpha;
	bool valid, timedep;
	lineiterations++;
	if ( (lineiterations==5)) {
		memcpy(x, x0, sizeof(*x)*nvar);
		ComputeG(g);
		valid = true;
		timedep = false;
		for (int i=0; i<nvar && valid && !timedep; i++) {
			if ( g[i]!=g0[i] && !timedep) {
				cout <<"[NEWTON:ERROR?: your functions are time dependent!]"<< endl; 
				timedep = true;
			} else if (!finite(g[i]) && valid) {
				cout <<"invalid numbers in gradient, reversing search direction"<<endl; ;
				valid = false;
				alpha *= -1; // reverse search direction
			}
		}
	}

	for (int i=0; i<nvar; i++) x[i] = x0[i]+alpha*p[i];
	valid = true;

	ComputeG(g);
	for (int i=0; i<nvar && valid; i++) {
		if (!finite(g[i])) {
			valid = false;
			cout <<"invalid numbers in gradient"<<endl;
				g[i] = 1;
		}
	}
	minimum = newfunction(g,x,nvar);
	return alpha;

}

double Newton::stepchange(double *g, double *g0, double *p, double *p0, double *x, double *x0, int nvar, int m, double &alpha){
	double change, crit;
	change = crit = linecriterion(g,g0,p,p0,nvar);
	while ( crit<0.35 && lineiterations<linesearchlimit ) {
		alpha /= 4;
		zero(g,g0,p,x,x0,nvar,m,alpha);
		crit = linecriterion(g,g0,p,p0,nvar);
		change = 1;
	}
	return change;
}

double Newton::linecriterion(double *g, double *g0, double *p, double *p0,int nvar) { 
	double normg,gg0;
	normg = norm2(g0,nvar);
	Dot(gg0,g,g0,nvar);
	gg0=gg0/normg/normg;
	normg = pow(norm2(g,nvar)/normg,2);
	if ( (gg0>1 || normg>1) && normg-gg0*fabs(gg0)<0.2 ) {
		normg = 1.5*normg;
	}
	if (gg0<0 && normg<2) {
		return 1;
	} else if (normg>10) {
		return .01;
	} else {
		return 0.4*(1+0.75*normg)/(normg-gg0*fabs(gg0)+0.1);
	}
}

void Newton::iterate(double *x,int nvar) {
	nbits=52;
	it = iterations=0; lineiterations=0; numIterationsSinceHessian = 0; 
	if (nvar<1) return;
	double alphamax=0; alphabound=0; alphaMax=delta_max; alphaMin=delta_min; 
	ALPHA=1;
	minAccuracySoFar = 1e30; reverseDirectionRange = 50; trouble=0; 
	resetHessianCriterion = 1e5; reset_pseudohessian = false; 
	trustregion=delta_max;
	trustfactor =1;
	int nh=0;
	if (hessianwidth*2-1>nvar || hessianwidth<1) nh=nvar; else nh=2*hessianwidth-1;  
	p = new double[nvar]; H_Zero(p,nvar); 
	g0 = new double[nvar]; H_Zero(g0,nvar); 
	p0 = new double[nvar]; H_Zero(p0,nvar); 
	h = new double [nvar*nh]; H_Zero(h,nvar*nh);  //consider to put hessian in Engine to enable samehessian (hessian-transfer)
	if (hessianwidth>0) newhessian(h,g,g0,x,p,nvar,hessianwidth); else newhessian(h,g,g0,x,p,nvar);

	if (e_info) {cout <<"NEWTON has been notified."<< endl; 
		cout << "Your guess:"; 
	}
	ComputeG(g); 
	minimum = newfunction(g,x,nvar);
	newdirection(h,p,p0,g,g0,x,nvar,hessianwidth,alphabound);
	normg=sqrt(minimum);
	while ((tolerance < accuracy || tolerance*10<normg) && iterations<iterationlimit && accuracy == fabs(accuracy) ) {
		if (e_info && it%i_info == 0){
			printf("it = %i \t E = %f \t |g| = %f \t alpha = %1e \n",it,accuracy,normg,ALPHA);
		}
		it++; iterations=it; lineiterations=0;
		newtrustregion(g,g0,p,p0,nvar);
		alphabound = alphamax = trustregion/(norm2(p,nvar)+1/pow(2.0,nbits));
		memcpy(x0, x, sizeof(*x0)*nvar);
		memcpy(g0, g, sizeof(*g0)*nvar);
		ALPHA = linesearch(g,g0,p,x,x0,nvar,hessianwidth,alphabound);
		trustfactor *= stepchange(g,g0,p,p0,x,x0,nvar,hessianwidth,ALPHA); // alpha is modified as well!
		trustfactor *= ALPHA/alphabound;
		newdirection(h,p,p0,g,g0,x,nvar,hessianwidth,alphabound);
		normg=sqrt(minimum); 
	}
	free(p); free(g0); free(p0); 
	free(h); 
}

bool Newton::PutU() {
if(debug) cout <<"PutU in  Newton " << endl;
	int M=Lat[0]->M;
	bool success=true;
	int sysmon_length = Sys[0]->SysMonList.size();
	for (int i=0; i<sysmon_length; i++) {
		double *u=Seg[Sys[0]->SysMonList[i]]->u;
		alpha=Sys[0]->alpha; 
		Cp(u,xx+i*M,M); 
		if (method == "Picard") Add(u,alpha,M);
		if (method == "DIIS-ext") Add(u,xx+sysmon_length*M,M); 
	}
	return success;
}

void Newton::Ax(double* A, double* X, int N){//From Ax_B; below B is not used: it is assumed to contain a row of unities.
if(debug) cout <<"Ax in  Newton " << endl;
	double* U = new double[N*N];
	double* S = new double[N];
	double* VT = new double[N*N];
	integer MM = (integer)N, NN = (integer)N;
	integer LDA=MM, LDU=MM, LDVT=NN, INFO, LWORK;
	int lwork;
	double WKOPT;
	double* WORK;
	char JOBU='S'; //'S' is nodig om alleen de eerste N colommen in U te schrijven.
	char JOBVT='A';

	LWORK = -1; //grootte hulpgeheugen aanvragen
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, &WKOPT, &LWORK, &INFO );
	lwork = (int)WKOPT;
	WORK = (double*)malloc( lwork*sizeof(double) );
	LWORK = (integer)lwork; //nu uitrekenen.
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT,WORK, &LWORK, &INFO );
	if (INFO >0) { cout <<"error in Ax " << endl; 
	};
	free(WORK);
	for (int i=0; i<N; i++) X[i]=0;
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += U[i*N + j];//*B[j];
	for (int i=0; i<N; i++) {S[i] = X[i]/S[i]; X[i]=0;} //S is use decause it is no longer needed.
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += VT[i*N + j]*S[j];
	delete U;
	delete S;
	delete VT;
}
void Newton::DIIS(double* xx, double* x_x0, double* xR, double* Aij, double* Apij,double* Ci, int k, int m, int iv) {
if(debug) cout <<"DIIS in  Newton " << endl;
	double normC=0; int posi;
	if (k_diis>m) { k_diis =m;
		for (int i=1; i<m; i++) for (int j=1; j<m; j++)
		Aij[m*(i-1)+j-1]=Aij[m*i+j]; //remove oldest elements
	}
	for (int i=0; i<k_diis; i++) {posi = k-k_diis+1+i; if (posi<0) posi +=m;
		double Dvalue; Dot(Dvalue,x_x0+posi*iv, x_x0+k*iv,iv);
		Aij[i+m*(k_diis-1)] = Aij[k_diis-1+m*i] = Dvalue; }
		// write to (compressed) matrix Apij
	for (int i=0; i<k_diis; i++) for (int j=0; j<k_diis; j++) {
		Apij[j+k_diis*i] = Aij[j+m*i];
	}
	Ax(Apij,Ci,k_diis);
	for (int i=0; i<k_diis; i++) normC +=Ci[i];
	for (int i=0; i<k_diis; i++) {Ci[i] =Ci[i]/normC; }
	Zero(xx,iv);
	posi = k-k_diis+1; if (posi<0) posi +=m;

	YplusisCtimesX(xx,xR+posi*iv,Ci[0],iv); //pv = Ci[0]*xR[0];
	for (int i=1; i<k_diis; i++) {
		posi = k-k_diis+1+i; if (posi<0) posi +=m;
		YplusisCtimesX(xx,xR+posi*iv,Ci[i],iv);
	}
}


bool Newton::Solve(void) {
if(debug) cout <<"Solve in  Newton " << endl;
	bool success=true;
	if (start==1) {
		if (ReadFileGuess!="") {
			Lat[0]->ReadGuess(ReadFileGuess,xx);
		} else {
			Lat[0]->GenerateGuess(xx,Sys[0]->CalculationType,Sys[0]->GuessType,Seg[Sys[0]->MonA]->guess_u,Seg[Sys[0]->MonB]->guess_u);
		}
	}
	if (method=="Picard") success=Iterate_Picard(); else success=Iterate_DIIS(); 
	Sys[0]->CheckResults();
	if (StoreFileGuess!="") Lat[0]->StoreGuess(StoreFileGuess,xx); 
	return success; 
}

void Newton::Message(int it, int iterationlimit,double residual, double tolerance) {
if(debug) cout <<"Messiage in  Newton " << endl;
	if (it < iterationlimit/10) cout <<"That was easy." << endl;
	if (it > iterationlimit/10 && it < iterationlimit ) cout <<"That will do." << endl;
	if (it <2 && iterationlimit >1 ) cout <<"You hit the nail on the head." << endl; 
	if (residual > tolerance) { cout << "Iterations failed." << endl;
		if (residual < tolerance/10) cout <<"I almost made it..." << endl;
	}
}

bool Newton::Iterate_Picard() {
if(debug) cout <<"Iterate_Picard in  Newton " << endl;
	int M=Lat[0]->M;
	double chi; 
	alpha=Sys[0]->alpha;	
	bool success=true;
	
	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size();	
	it=0;
	cout <<"Picard has been notified" << endl; 
	residual=1;
	while (residual > tolerance && it < iterationlimit) {
		Cp(x0,xx,iv); 
		ComputePhis(); 
		Zero(xx,iv);
		for (int i=0; i<sysmon_length; i++) for (int k=0; k<mon_length; k++) { 
                        chi= -1.0*Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k];  //The minus sign here is to change the sign of x! just a trick due to properties of PutAlpha where a minus sing is implemented....
			if (chi!=0) PutAlpha(xx+i*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
		}
		for (int i=0; i<sysmon_length; i++) {
			Lat[0]->remove_bounds(xx+i*M);
			Times(xx+i*M,xx+i*M,Sys[0]->KSAM,M);
		}

		Picard(xx,x0,delta_max,iv); 
		YisAminB(g,xx,x0,iv); 
		for (int i=0; i<sysmon_length; i++) Times(g+i*M,g+i*M,Sys[0]->KSAM,M);  

		Dot(residual,g,g,iv);
		UpdateAlpha(alpha, Sys[0]->phitot, delta_max, M);
		YisAplusC(g,Sys[0]->phitot,-1.0,M); 
		Lat[0]->remove_bounds(g); 

		double result; Dot(result,g,g,M);
		residual=residual+result;
		residual=sqrt(residual); 
		if(it%i_info == 0){
			printf("it = %i g = %1e \n",it,residual);
		}
		it++; 
	}
	Message(it,iterationlimit,residual,tolerance); 
	return success; 
}

bool Newton::Iterate_DIIS() {
if(debug) cout <<"Iterate_DIIS in  Newton " << endl;
	Zero(x0,iv);
	it=0; k_diis=1; 
	int k=0;
	if (method=="DIIS-ext") ComputeG_ext(); else ComputeG(g);
	YplusisCtimesX(xx,g,delta_max,iv);
	YisAminB(x_x0,xx,x0,iv);
	Cp(xR,xx,iv);
	Dot(residual,g,g,iv);
	residual=sqrt(residual); 
	printf("DIIS has been notified\n");
	printf("Your guess = %1e \n",residual);
	while (residual > tolerance && it < iterationlimit) {
		it++;
		Cp(x0,xx,iv);
		if (method=="DIIS-ext") ComputeG_ext(); else ComputeG(g);
		k=it % m; k_diis++; //plek voor laatste opslag
		YplusisCtimesX(xx,g,-delta_max,iv);
		Cp(xR+k*iv,xx,iv); YisAminB(x_x0+k*iv,xx,x0,iv);
		DIIS(xx,x_x0,xR,Aij,Apij,Ci,k,m,iv);
		Dot(residual,g,g,iv);
		residual=sqrt(residual);
		if(it%i_info == 0){
			printf("it = %i g = %1e \n",it,residual);
		}
	}

	Message(it,iterationlimit,residual,tolerance); 
	return it<iterationlimit+1;
} 

void Newton::ComputePhis() {
if(debug) cout <<"ComputPhis in  Newton " << endl;
	PutU();
	Sys[0]->PrepareForCalculations();
	Sys[0]->ComputePhis();
}

void Newton::ComputeG(double* g){ 
if(debug) cout <<"ComputeG in  Newton " << endl;
	int M=Lat[0]->M;
	double chi; 
	ComputePhis();

	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size();	

	Cp(g,xx,iv); Zero(alpha,M);
	for (int i=0; i<sysmon_length; i++) {
		for (int k=0; k<mon_length; k++) { 
                        chi= Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k];  
			if (chi!=0) {
				PutAlpha(g+i*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
			}

		}
	}
	for (int i=0; i<sysmon_length; i++) {
		Add(alpha,g+i*M,M);
	}

	Norm(alpha,1.0/sysmon_length,M);
	for (int i=0; i<sysmon_length; i++) {
		AddG(g+i*M,Sys[0]->phitot,alpha,M);
		Lat[0]->remove_bounds(g+i*M);
		//Times(g+i*M,g+i*M,Sys[0]->KSAM,M);
	} 
}

void Newton::ComputeG_ext(){ 
if(debug) cout <<"CompueG_est() in  Newton " << endl;
	int M=Lat[0]->M;
	alpha=Sys[0]->alpha;
	double chi; 
	ComputePhis();
	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size();	

	Cp(g,xx,iv); 
	for (int i=0; i<sysmon_length; i++) { 
		for (int k=0; k<mon_length; k++){ 
                        chi= Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k];  
			if (chi!=0) PutAlpha(g+i*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
		}
	}
	for (int i=0; i<sysmon_length; i++){
		Lat[0]->remove_bounds(g+i*M);
		Times(g+i*M,g+i*M,Sys[0]->KSAM,M);
	}
	YisAplusC(g+sysmon_length*M,Sys[0]->phitot, -1.0, M);
	Lat[0]->remove_bounds(g+sysmon_length*M);
	Norm(g+sysmon_length*M,-1,M);//you can improve ......
	//Norm(g,-1,sysmon_length*M);
	
}
