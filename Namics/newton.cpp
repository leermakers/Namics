#include "newttool.h"
#include "newton.h"

Newton::Newton(vector<Input*> In_,vector<Lattice*> Lat_,vector<Segment*> Seg_,vector<Molecule*> Mol_,vector<System*> Sys_,vector<Variate*>Var_,string name_) {
	In=In_; name=name_; Sys=Sys_; Seg=Seg_; Lat=Lat_; Mol=Mol_;Var=Var_;
if(debug) cout <<"Constructor in Newton " << endl;
	KEYS.push_back("method");  KEYS.push_back("m");
	KEYS.push_back("e_info"); KEYS.push_back("s_info");KEYS.push_back("i_info");

	KEYS.push_back("iterationlimit" ); KEYS.push_back("tolerance");
	KEYS.push_back("stop_criterion");
	KEYS.push_back("deltamin");KEYS.push_back("deltamax");
	KEYS.push_back("linesearchlimit");
	//KEYS.push_back("samehessian");
	KEYS.push_back("max_accuracy_for_hessian_scaling");
	KEYS.push_back("n_iterations_for_hessian");
	KEYS.push_back("small_alpha");
	KEYS.push_back("max_n_small_alpha");
	KEYS.push_back("min_accuracy_for_hessian");
	KEYS.push_back("max_fr_reverse_direction");
	//KEYS.push_back("print_hessian_at_it");
	KEYS.push_back("super_e_info");
	KEYS.push_back("super_s_info");
	KEYS.push_back("super_tolerance");
	KEYS.push_back("super_iterationlimit");
}

Newton::~Newton() {
	DeAllocateMemory();
}

void Newton :: DeAllocateMemory(){
if(debug) cout <<"Destructor in Newton " << endl;
	free(Aij);
	free(Ci);
	free(Apij);
	free(mask);
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
	if (method == "DIIS-mesodyn") {
		iv = Sys[0]->SysMolMonList.size()*M;
	} else iv = Sys[0]->SysMonList.size() * M;
	if (Sys[0]->charged) iv +=M;
	if (method=="DIIS-ext") iv +=M;
	Aij  =(Real*) malloc(m*m*sizeof(Real)); H_Zero(Aij,m*m);
	Ci   =(Real*) malloc(m*sizeof(Real)); H_Zero(Ci,m);
	Apij =(Real*) malloc(m*m*sizeof(Real)); H_Zero(Apij,m*m);
	mask= (int*) malloc(iv*sizeof(int));
#ifdef CUDA
	xx  = (Real*)AllOnDev(iv);
	x0  = (Real*)AllOnDev(iv);
	g   = (Real*)AllOnDev(iv);
	xR  = (Real*)AllOnDev(m*iv);
	x_x0= (Real*)AllOnDev(m*iv);
#else
	xx   =(Real*) malloc(iv*sizeof(Real));
	x0   =(Real*) malloc(iv*sizeof(Real));
	g    =(Real*) malloc(iv*sizeof(Real));
	xR   =(Real*) malloc(m*iv*sizeof(Real));
	x_x0 =(Real*) malloc(m*iv*sizeof(Real));
#endif
	Zero(xx,iv);
	Zero(x0,iv);
	Zero(g,iv);
	Zero(xR,m*iv);
	Zero(x_x0,m*iv);
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
		if (iterationlimit < 0 || iterationlimit>1e6) {iterationlimit = 1000;}

		super_e_info=In[0]->Get_bool(GetValue("super_e_info"),e_info);
		super_s_info=In[0]->Get_bool(GetValue("super_s_info"),s_info);
		super_iterationlimit=In[0]->Get_int(GetValue("super_iterationlimit"),iterationlimit/10);

		delta_max=In[0]->Get_Real(GetValue("deltamax"),0.1);
		if (delta_max < 0 || delta_max>100) {delta_max = 0.1;  cout << "Value of deltamax out of range 0..100, and value set to default value 0.1" <<endl; }
		delta_min=delta_max/100000;
		delta_min=In[0]->Get_Real(GetValue("deltamin"),delta_min);
		if (delta_min < 0 || delta_min>100) {delta_min = delta_max/100000;  cout << "Value of deltamin out of range 0..100, and value set to default value deltamax/100000" <<endl; }

		tolerance=In[0]->Get_Real(GetValue("tolerance"),1e-7);
		super_tolerance=In[0]->Get_Real(GetValue("super_tolerance"),tolerance*10);
		if (tolerance < 1e-12 ||tolerance>10) {tolerance = 1e-5;  cout << "Value of tolerance out of range 1e-12..10 Value set to default value 1e-5" <<endl; }
		if (GetValue("method").size()==0) {method="DIIS";} else {
			vector<string>method_options;
			method_options.push_back("DIIS");
			method_options.push_back("DIIS-mesodyn");
			//method_options.push_back("DIIS-ext"); //can be included again when adjusted for charges and guess
			//method_options.push_back("Picard"); //can be included again when adjusted for charges and guess
			method_options.push_back("pseudohessian");
			method_options.push_back("hessian");
			if (!In[0]->Get_string(GetValue("method"),method,method_options,"In 'newton' the entry for 'method' not recognized: choose from:")) success=false;
		}
		m=In[0]->Get_int(GetValue("m"),10);
		if (method=="pseudohessian" || method == "hessian") {
			print_hessian_at_it=-1;
			if (GetValue("print_hessian_at_it").size()>0) {
				print_hessian_at_it=In[0]->Get_int(GetValue("print_hessian_at_it"),-123);
				if (print_hessian_at_it <0) cout <<"print_hessian_at_it did not find positive integer value. Input ignored " << endl;
			}
			if (method=="hessian") {pseudohessian=false; hessian=true;} else { pseudohessian=true; hessian=false; }
			samehessian=false; //In[0]->Get_bool(GetValue("samehessian"),false);
			linesearchlimit=In[0]->Get_int(GetValue("hessianwidth"),20);
			if (linesearchlimit<0 || linesearchlimit >100) {
				cout <<"linesearchlimit is out of range: 1..100; default value of 20 is used instead" << endl;
				linesearchlimit=20;
			}
			max_accuracy_for_hessian_scaling=In[0]->Get_Real(GetValue("max_accuracy_for_hessian_scaling"),0.1);
			if (max_accuracy_for_hessian_scaling<1e-7 || max_accuracy_for_hessian_scaling>1) {
				cout <<"max_accuracy_for_hessian_scaling is out of range: 1e-7...1; default value 0.1 is used instead" << endl;
				max_accuracy_for_hessian_scaling=0.1;
			}
			minAccuracyForHessian=In[0]->Get_Real(GetValue("min_accuracy_for_hessian"),0);
			if (minAccuracyForHessian<0 ||minAccuracyForHessian>0.1) {
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

			delta_min=In[0]->Get_Real(GetValue("delta_min"),0);
			if (delta_min <0 || delta_min>delta_max) {
				cout <<"delta_min is out of range; 0, ..., " << delta_max << "; delta_min value set to 0 " << endl;
				delta_min=0;
			}
			smallAlpha=In[0]->Get_Real(GetValue("small_alpha"),0.00001);
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

void Newton::push(string s, Real X) {
if(debug) cout <<"push (Real) in  Newton " << endl;
	Reals.push_back(s);
	Reals_value.push_back(X);
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
	Reals.clear();
	Reals_value.clear();
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
	if (pseudohessian||hessian) {
		//push("same_hessian",samehessian);
		push("linesearchlimit",linesearchlimit);
		push("max_accuracy_for_hessian_scaling",max_accuracy_for_hessian_scaling);
		push("n_iteratons_for_hessian",n_iterations_for_hessian);
		push("small_alpha",smallAlpha);
		push("max_n_small_alpha",maxNumSmallAlpha);
		push("min_accuracy_for_hessian",minAccuracyForHessian);
	}
}

int Newton::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
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

void Newton::inneriteration(float*h, Real *g, Real *x,Real accuracy, int nvar) {
if(debug) cout <<"inneriteration in Newton " << endl;

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

		if (delta_max >0.005) delta_max *=0.9;
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
		if (delta_max >0.005) delta_max *=0.9;
		numIterationsSinceHessian = 0;
	}

	if (!newtondirection && pseudohessian) {
		reverseDirection[iterations%reverseDirectionRange] = 1;
	} else {
		reverseDirection[iterations%reverseDirectionRange] = 0;
	}

	numReverseDirection = 0;
	for (int i=1; i<=reverseDirectionRange; i++) {
		if (reverseDirection[i] > 0.5) numReverseDirection++;
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
}

void Newton::multiply(Real *v,Real alpha, float *h, Real *w, int nvar) {
if(debug) cout <<"multiply in Newton" << endl;
	int i=0,i1=0,j=0;
	Real sum=0;
	Real *x = new Real[nvar];
	for (i=0; i<nvar; i++) {
		sum = 0;
		i1 = i-1;
		for (j=i+1; j<nvar; j++) {
			sum += w[j] * h[i+nvar*j];
		}
		x[i] = (sum+w[i])*h[i+nvar*i];
		sum = 0;
		for (j=0; j<=i1; j++) {
			sum += x[j] * h[i+nvar*j];
		}
		v[i] = alpha*(sum+x[i]);
	}
	delete [] x;
}

Real Newton::norm2(Real *x, int nvar) {
if(debug) cout <<"norm2 in Newton" << endl;

	Real sum=0;
	for (int i=0; i<nvar; i++) sum += pow(x[i],2);
	return sqrt(sum);
}

int Newton::signdeterminant(float *h,int nvar) {
if(debug) cout <<"signdeterminant in Newton" << endl;
	int sign=1;
	for (int i=0; i<nvar; i++) {
		if ( h[i+i*nvar]<0 ) {
			sign = -sign;
		}
	}
	return sign;
}

void Newton::updateneg(float *l,Real *w, int nvar, Real alpha) {
if(debug) cout <<"updateneg in Newton" << endl;
	int i=0,i1=0,j=0;
	Real dmin=0,sum=0,b=0,d=0,p=0,lji=0,t=0;
	dmin = 1.0/pow(2.0,54);
	alpha = sqrt(-alpha);
	for (i=0; i<nvar; i++) {
		i1 = i-1;
		sum = 0;
		for (j=0;j<=i1; j++) {
			sum += l[i+nvar*j]*w[j];
		}
		w[i] = alpha*w[i]-sum;
		t += (w[i]/l[i+nvar*i])*w[i];
	}
	t = 1-t;
	if ( t<dmin ) t = dmin;
	for (i=nvar-1; i>=0; i--) {
		p = w[i];
		d = l[i+nvar*i];
		b = d*t;
		t += (p/d)*p;
		l[i+nvar*i] = b/t;
		b = -p/b;
		for (j=i+1; j<nvar; j++) {
			lji = l[j+nvar*i];
			l[j+nvar*i] = lji+b*w[j];
			w[j] += p*lji;
		}
	}
}

void Newton::decompos(float *h, int nvar, int &ntr) {
if(debug) cout <<"decompos in Newton" << endl;
	int i,j,k;//itr,ntr;
	Real sum,lsum,usum,phi,phitr,c,l;
	float *ha,*hai,*haj;
	ha = &h[-1];
	phitr = FLT_MAX;
	ntr = 0;
	i = 0;
	while (i++<nvar) {
		hai = &ha[(i-1)*nvar];
		sum = 0;
		j = 0;
		while (j++<i-1) {
			haj = &ha[(j-1)*nvar];
			c = haj[i];
			l = c/haj[j];
			haj[i] = l;
			c = hai[j];
			hai[j] = c/haj[j];
			sum += l*c;
		}
		phi = hai[i] - sum;
		hai[i] = phi;
		if (phi<0) ntr++;
		if (phi<phitr) phitr = phi;
		j = i;
		while (j++<nvar) {
			haj = &ha[(j-1)*nvar];
			lsum = 0;
			usum = 0;
			k = 0;
			while (k++<i-1) {
				lsum += ha[(k-1)*nvar+j]*hai[k];
				usum += ha[(k-1)*nvar+i]*haj[k];
			}
			hai[j] -= lsum;
			haj[i] -= usum;
		}
	}
}

void Newton::updatpos(float *l, Real *w, Real *v, int nvar, Real alpha) {
if(debug) cout <<"updatepos in Newton" << endl;
	int i,j;
	Real b,c,d;
	Real vai,waj,vaj;
	float *lai,*laj;
	Real * wa = &w[-1];
	Real * va = &v[-1];
	i = 0;
	while (i++<nvar) {
		vai = va[i];
		lai = &l[-1 + (i-1)*nvar];
		d = lai[i];
		b = d+(alpha*wa[i])*vai;
		lai[i] = b;
		d /= b;
		c = vai*alpha/b;
		b = wa[i]*alpha/b;
		alpha *= d;
		j = i;
		while (j++<nvar) {
			waj = wa[j];
			wa[j] -= wa[i]*lai[j];
			lai[j] *= d;
			lai[j] += c*waj;
			laj = &l[-1 + (j-1)*nvar];
			vaj = va[j];
			va[j] -= vai*laj[i];
			laj[i] *= d;
			laj[i] += b*vaj;
		}
	}
}

void Newton::gausa(float *l, Real *dup, Real *g, int nvar) {
if(debug) cout <<"gausa in Newton" << endl;
	int i,j;
	Real*dupa,sum;
	Real *ga;
	float *lai;

	dupa = &dup[-1];
	ga = &g[-1];

	i = 0;
	while (i++<nvar) {
		sum = 0;
		lai = &l[i-1];
		j = 0;
		while (j++<i-1) {
			sum += lai[(j-1)*nvar]*dupa[j];
		}
		dupa[i] = - ga[i] - sum;
	}
}

void Newton::gausb(float *du, Real *p, int nvar) {
if(debug) cout <<"gausb in Newton " << endl;
	int i,j;
	Real *pa,sum;
	float *duai;
	pa = &p[-1];
	i = nvar+1;
	while (i-- > 1) {
		sum = 0;
		duai = &du[i-1];
		j = i;
		while (j++<nvar) {
			sum += duai[(j-1)*nvar]*pa[j];
		}
		pa[i] = pa[i]/duai[(i-1)*nvar] - sum;
	}
}

Real Newton::residue(Real *g, Real *p, Real *x, int nvar, Real alpha) {
if(debug) cout <<"residue in Newton " << endl;
	return sqrt(norm2(p,nvar)*norm2(g,nvar)/(1+norm2(x,nvar)));
}

Real Newton::linecriterion(Real *g, Real *g0, Real *p, Real *p0, int nvar) {
if(debug) cout <<"linecriterion in Newton " << endl;
	Real normg,gg0;
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

Real Newton::newfunction(Real *g, Real *x, int nvar) {
if(debug) cout <<"newfunction in Newton " << endl;
	return pow(norm2(g,nvar),2);
}

void Newton::direction(float *h, Real *p, Real *g, Real *g0, Real *x, int nvar, Real alpha){
if(debug) cout <<"direction in Newton " << endl;

	newtondirection = true;
	newhessian(h,g,g0,x,p,nvar);
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

void Newton::startderivatives(float *h, Real *g, Real *x, int nvar){
if(debug) cout <<"startderivatives in Newton" << endl;
	float diagonal = 1+norm2(g,nvar);
	H_Zero(h,nvar*nvar);
	for (int i=0; i<nvar; i++) {
		h[i+nvar*i] = diagonal;
	}
}

void Newton::resethessian(float *h,Real *g,Real *x,int nvar){
if(debug) cout <<"resethessian in Newton" << endl;
	trouble = 0;
	startderivatives(h,g,x,nvar);
	resetiteration = iterations;
}

void Newton::newhessian(float *h, Real *g, Real *g0, Real *x, Real *p, int nvar) {
if(debug) cout <<"newhessian in Newton" << endl;
	Real dmin=0,sum=0,theta=0,php=0,dg=0,gg=0,g2=0,py=0,y2=0;
	dmin = 1/pow(2.0,nbits); // alternative: DBL_EPSILON or DBL_MIN
	if (!pseudohessian){
		findhessian(h,g,x,nvar);
	} else {
		if (!samehessian && ALPHA!=0 && iterations!=0) {
			Real *y = new Real[nvar];
			Real *hp = new Real[nvar];
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
				if ( trouble<0 ) trouble = 0; else if ( trouble>=3 ) resethessian(h,g,x,nvar);
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
				updatpos(h,y,y,nvar,1.0/(ALPHA*py));
				updateneg(h,hp,nvar,-1.0/php);
			}

			delete [] y;
			delete [] hp;
		} else if ( !samehessian ) resethessian(h,g,x,nvar);
	}
}
void print_hessian(float* h, int nvar) {
	cout <<"hessian: " << endl;
	for (int i=0; i<nvar; i++)
	for (int j=0; j<nvar; j++)
		if (h[i*nvar+j]!=0) cout <<"i " << i << " j " << j << " h= " << h[i*nvar+j] << endl;
}

void Newton::numhessian(float* h,Real* g, Real* x, int nvar) {
if(debug) cout <<"numhessian in Newton" << endl;
	Real dmax2=0,dmax3=0,di=0;
	Real *g1;
	g1 = new Real[iv];
	Real xt;
	dmax2 = pow(2.0,nbits/2); //alternative 2*pow(DBL_EPSILON,-0.5)?
	dmax3 = pow(2.0,nbits/3); //alternative 2*pow(DBL_EPSILON,-1.0/3)?
	for (int i=0; i<nvar; i++) {
		xt = x[i];
		di = (1/(dmax3*dmax3*fabs(h[i+nvar*i])+dmax3+fabs(g[i]))
			+1/dmax2)*(1+fabs(x[i]));
		if ( di<delta_min ) {
			di = delta_min;
		}
		x[i] += di;
		COMPUTEG(x,g1,nvar);
		x[i] = xt;
		for (int j=0; j<nvar; j++ ) {
			h[j+nvar*i] = (g1[j]-g[j])/di;
		}
	}
	delete [] g1;
	COMPUTEG(x,g,nvar);
}

void Newton::decomposition(float *h,int nvar, int &trouble){
if(debug) cout <<"decomposition in Newton" << endl;
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

void Newton::findhessian(float *h, Real *g, Real *x,int nvar) {
if(debug) cout <<"findhessian in Newton" << endl;
	if ( !samehessian ) {
		if ( iterations==0 ) resethessian(h,g,x,nvar);
		numhessian(h,g,x,nvar); // passes through residuals so check pseudohessian
		if (!pseudohessian) {
			decomposition(h,nvar,trouble);
		}
	}
}

void Newton::newdirection(float *h, Real *p, Real *p0, Real *g, Real *g0, Real *x, int nvar, Real alphabound) {
if(debug) cout <<"newdirection in Newton" << endl;
	inneriteration(h,g,x,accuracy,nvar);
	memcpy(p0, p, sizeof(*p0)*nvar);
	direction(h,p,g,g0,x,nvar,ALPHA);
	accuracy = residue(g,p,x,nvar,ALPHA);
}

void Newton::newtrustregion(Real *g, Real *g0, Real *p, Real *p0, int nvar){
if(debug) cout <<"newtrustregion in Newton" << endl;
	Real normp0 = norm2(p0,nvar);

	if ( normp0>0 && trustregion>2*ALPHA*normp0 ) {
		trustregion = 2*ALPHA*normp0;
	}
	trustregion *= trustfactor;
	trustfactor = 1.0;
	if ( trustregion>delta_max ) trustregion = delta_max;
	if ( trustregion<delta_min ) trustregion = delta_min;
}

Real Newton::linesearch(Real *g, Real *g0, Real *p, Real *x, Real *x0, int nvar, Real alphabound) {
if(debug) cout <<"linesearch in Newton" << endl;
	Real newalpha = alphabound<1 ? alphabound : 1;
	newalpha = zero(g,g0,p,x,x0,nvar,newalpha);
	return newalpha;
}

Real Newton::zero(Real *g, Real *g0, Real *p, Real *x, Real *x0, int nvar, Real newalpha) {
if(debug) cout <<"zero in Newton " << endl;
	Real alpha=newalpha;
	bool valid, timedep;
	lineiterations++;
	if ( (lineiterations==5)) {
		memcpy(x, x0, sizeof(*x)*nvar);
		COMPUTEG(x,g,nvar);
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

	COMPUTEG(x,g,nvar);
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

Real Newton::stepchange(Real *g, Real *g0, Real *p, Real *p0, Real *x, Real *x0, int nvar, Real &alpha){
if(debug) cout <<"stepchange in Newton" << endl;
	Real change, crit;
	change = crit = linecriterion(g,g0,p,p0,nvar);
	while ( crit<0.35 && lineiterations<linesearchlimit ) {
		alpha /= 4;
		zero(g,g0,p,x,x0,nvar,alpha);
		crit = linecriterion(g,g0,p,p0,nvar);
		change = 1;
	}
	return change;
}

void Newton::COMPUTEG(Real *x, Real *g, int nvar) {
	int pos=nvar;
	for (int i=iv-1; i>=0; i--) {
		if (mask[i]==1) {
			pos--;
			x[i]=x[pos];
		} else x[i]=0;
	}
	ComputeG(g);
	pos=0;
	for (int i=0; i<iv; i++) {
		if (mask[i]==1) {x[pos]=x[i]; g[pos]=g[i];pos++;}
	}
}

void Newton::ResetX(Real* x,int nvar) {
	int pos=nvar;
	for (int i=iv-1; i>=0; i--) {
		if (mask[i]==1) {
			pos--;
			x[i]=x[pos];
		} else x[i]=0;
	}
}

/*	Iteration scheme used by pseudohessian
*/
void Newton::iterate(Real *x,int nvar) {
if(debug) cout <<"iterate in Newton" << endl;
	nbits=52;
	ignore_newton_direction = false;
	it = iterations=0; lineiterations=0; numIterationsSinceHessian = 0;
	if (nvar<1) return;
	Real alphamax=0; alphabound=0; alphaMax=delta_max; alphaMin=delta_min;
	ALPHA=1;
	minAccuracySoFar = 1e30; reverseDirectionRange = 50; trouble=0;
	resetiteration=0;
	reverseDirection = new int [reverseDirectionRange]; H_Zero(reverseDirection,reverseDirectionRange);
	resetHessianCriterion = 1e5; reset_pseudohessian = false;
	trustregion=delta_max;
	trustfactor =1;
	srand (1);
	Cp(x0,x,iv);
	for (int i=0; i<nvar; i++) x[i]+=1e-10*(Real)rand() / (Real)((unsigned)RAND_MAX + 1);
	ComputeG(g);
	Cp(x,x0,iv);
	int xxx=0;
	for (int i=0; i<nvar; i++) {if (g[i]==0) mask[i]=0; else {xxx++;  mask[i]=1;}}

	nvar=xxx;
	p = new Real[nvar]; H_Zero(p,nvar);
	g0 = new Real[nvar]; H_Zero(g0,nvar);
	p0 = new Real[nvar]; H_Zero(p0,nvar);
	h = new float [nvar*nvar]; H_Zero(h,nvar*nvar);

	if (e_info) {cout <<"NEWTON has been notified."<< endl;
		cout << "Your guess:";
	}

	int pos=0;
	for (int i=0; i<iv; i++) {
		if (mask[i]==1) {g[pos]=g[i]; x[pos]=x[i];pos++;}
	}

	newhessian(h,g,g0,x,p,nvar);
	minimum = newfunction(g,x,nvar);
	newdirection(h,p,p0,g,g0,x,nvar,alphabound);
	normg=sqrt(minimum);
	if (print_hessian_at_it==0) print_hessian(h,nvar);

	while ((tolerance < accuracy || tolerance*10<normg) && iterations<iterationlimit && accuracy == fabs(accuracy) ) {
		if (e_info && it%i_info == 0){
			printf("it =  %i  E = %e |g| = %e alpha = %e \n",it,accuracy,normg,ALPHA);
		}
		it++; iterations=it; lineiterations=0;
		newtrustregion(g,g0,p,p0,nvar);
		alphabound = alphamax = trustregion/(norm2(p,nvar)+1/pow(2.0,nbits));
		memcpy(x0, x, sizeof(*x0)*nvar);
		memcpy(g0, g, sizeof(*g0)*nvar);
		ALPHA = linesearch(g,g0,p,x,x0,nvar,alphabound);
		trustfactor *= stepchange(g,g0,p,p0,x,x0,nvar,ALPHA); // alpha is modified as well!
		trustfactor *= ALPHA/alphabound;
		//if (it==1) {newhessian(h,g,g0,x,p,nvar);}
		newdirection(h,p,p0,g,g0,x,nvar,alphabound);
		normg=sqrt(minimum);
		if (print_hessian_at_it==it) print_hessian(h,nvar);
	}
	Message(e_info,s_info,it,iterationlimit,accuracy,tolerance,"");
	ResetX(xx,nvar);
	delete [] p; delete [] g0 ; delete [] p0;
	delete [] h; delete [] reverseDirection;
}

bool Newton::PutU() {
if(debug) cout <<"PutU in  Newton " << endl;
	int M=Lat[0]->M;
	bool success=true;
	int sysmon_length = Sys[0]->SysMonList.size();
	alpha=Sys[0]->alpha;
	if (method=="DIIS-mesodyn") {
		int i=0; int k=0;
		int length = In[0]->MolList.size();
		while (i<length) {
			int j=0;
			int LENGTH=Mol[i]->MolMonList.size();
			while (j<LENGTH) {Cp(Mol[i]->u+j*M,xx+k*M,M); k++; j++;}
			i++;
		}
	} else {
		for (int i=0; i<sysmon_length; i++) {
			Real *u=Seg[Sys[0]->SysMonList[i]]->u;
			Cp(u,xx+i*M,M);
			if (method == "Picard") Add(u,alpha,M);
			if (method == "DIIS-ext") Add(u,xx+sysmon_length*M,M);
		}
		int i=0;
		int length = In[0]->MolList.size();
		while (i<length) {
			int j=0;
			int LENGTH=Mol[i]->MolMonList.size();
			while (j<LENGTH) {Cp(Mol[i]->u+j*M,Seg[Mol[i]->MolMonList[j]]->u,M); j++;}
			i++;
		}
	}
	return success;
}

/*
void Newton::Ax(Real* A, Real* X, int N){//From Ax_B; below B is not used: it is assumed to contain a row of unities.
if(debug) cout <<"Ax in  Newton " << endl;
	Real* U = new Real[N*N];
	Real* S = new Real[N];
	Real* VT = new Real[N*N];
#ifdef  __CLAPACK_H
	integer MM = (integer)N, NN = (integer)N;
	integer LDA=MM, LDU=MM, LDVT=NN, INFO, LWORK;
#else
        int MM = (int)N, NN = (int)N;
	int LDA=MM, LDU=MM, LDVT=NN, INFO, LWORK;
#endif

	int lwork;
	Real WKOPT;
	Real* WORK;
	char JOBU='S'; //'S' is nodig om alleen de eerste N colommen in U te schrijven.
	char JOBVT='A';

	LWORK = -1; //grootte hulpgeheugen aanvragen
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT, &WKOPT, &LWORK, &INFO );
	lwork = (int)WKOPT;
	WORK = (Real*)malloc( lwork*sizeof(Real) );
#ifdef  __CLAPACK_H
	LWORK = (integer)lwork; //nu uitrekenen.
#else
       LWORK = (int)lwork; //nu uitrekenen.
#endif
	dgesvd_( &JOBU, &JOBVT, &MM, &NN, A, &LDA, S, U, &LDU, VT, &LDVT,WORK, &LWORK, &INFO );
	if (INFO >0) { cout <<"error in Ax " << endl;
	};
	free(WORK);
	for (int i=0; i<N; i++) X[i]=0;
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += U[i*N + j];// *B[j];
	for (int i=0; i<N; i++) {S[i] = X[i]/S[i]; X[i]=0;} //S is used decause it is no longer needed.
	for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += VT[i*N + j]*S[j];
	delete [] U;
	delete [] S;
	delete [] VT;
}
*/

void Newton::Ax(Real* A, Real* X, int N){//From Ax_B; below B is not used: it is assumed to contain a row of unities.
if(debug) cout <<"Ax in  Newton (own svdcmp) " << endl;

	Real **U = new Real*[N];
	Real **V = new Real*[N];
	Real *S = new Real[N];

	for (int i=0; i < N; i++) {
		U[i] = new Real[N];
		V[i] = new Real[N];
	}

	for (int i=0; i<N; i++)
		for (int j=0; j<N; j++) 
			U[j][i] = A[i*N + j];
  if (N > 1) {
		//old function svdcmp still exists, simply remove modern_ prefix to switch back. The new function uses vectors for safety.
  	modern_svdcmp(U, N, N, S, V);
		if (debug) cout << "SVCDMP done, continuing.." << endl;
		for (int i=0; i<N; i++) X[i]=0;
		for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += U[j][i];// *B[j];
		for (int i=0; i<N; i++) {S[i] = X[i]/S[i]; X[i]=0;} //S is use because it is no longer needed.
		for (int i=0; i<N; i++) for (int j=0; j<N; j++) X[i] += V[i][j]*S[j];
	} else {
		X[0]=1;
	}

for (int i=0; i<N; i++) {delete [] U[i]; delete [] V[i];}
	delete [] U;
	delete [] S;
	delete [] V;
}

void Newton::DIIS(Real* xx, Real* x_x0, Real* xR, Real* Aij, Real* Apij,Real* Ci, int k, int k_diis, int m, int iv) {
if(debug) cout <<"DIIS in  Newton " << endl;
	Real normC=0; int posi;
	if (k_diis>m) { k_diis =m;
		for (int i=1; i<m; i++) for (int j=1; j<m; j++)
		Aij[m*(i-1)+j-1]=Aij[m*i+j]; //remove oldest elements
	}
	for (int i=0; i<k_diis; i++) {posi = k-k_diis+1+i; if (posi<0) posi +=m;
		Real Dvalue; Dot(Dvalue,x_x0+posi*iv, x_x0+k*iv,iv);
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

string Newton::GetNewtonInfo(int &IV) {
	IV=iv;
	return method;
}

void Newton::Copy(Real* x, Real* X, int MX, int MY, int MZ) {
	int mx=Lat[0]->MX;
	int my=Lat[0]->MY;
	int mz=Lat[0]->MZ;
	int jx=Lat[0]->JX;
	int jy=Lat[0]->JY;
	int i,j,k;
	int pos_i,pos_o;
	int JX=(MY+2)*(MZ+2);
	int JY=(MZ+2);


	switch (Lat[0]->gradients) {
		case 1:
			if (MY>0||MZ>0) {
				cout <<" Copy from more than one gradient to one gradient: (i) =(1,i) or (1,1,i) is used "<< endl;
			}
			if (MZ>0) { pos_i=JX+JY; pos_o=MZ+2;} else {if (MY>0) {pos_i=JX; pos_o=MY+2; } else { pos_i=0; pos_o=MX+2; } }
			for (i=0; i<mx+2; i++)  if (i<pos_o) x[i]=X[pos_i+i];
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

bool Newton::Guess(Real *X, string METHOD, vector<string> MONLIST, bool CHARGED, int MX, int MY, int MZ){
	if (debug) cout << "Guess in Newton" << endl;
	int M=Lat[0]->M;
	bool success=true;

	if (start ==1 && Sys[0]->GuessType != "")  {
		Lat[0]->GenerateGuess(xx,Sys[0]->CalculationType,Sys[0]->GuessType,Seg[Sys[0]->MonA]->guess_u,Seg[Sys[0]->MonB]->guess_u);
	} else {
		int m;
		if (MZ>0) {m=(MX+2)*(MY+2)*(MZ+2); } else { if (MY>0) { m=(MX+2)*(MY+2); } else {  m=(MX+2);}}
		int length_old_mon=MONLIST.size();
		int length_new_mon=Sys[0]->SysMonList.size();
		for (int i = 0; i<length_old_mon; i++) {
			for (int j=0; j<length_new_mon; j++) {
				if (MONLIST[i]==Seg[Sys[0]->SysMonList[j]]->name) {
					Copy(xx+M*j,X+i*m,MX,MY,MZ);
				}
			}
		}
		if (CHARGED && Sys[0]->charged) {cout <<"both charged" << endl;  Copy(xx+length_new_mon*M,X+length_old_mon*m,MX,MY,MZ); }
	}
	return success;
}

/* 	Decides which computation method to use, depending on the method variable.
		Iterate_DIIS supports mesodyn and the classic way.
		Called by main
*/
bool Newton::Solve(bool report_errors) {
if(debug) cout <<"Solve in  Newton " << endl;
	bool success=true;
	if (method=="pseudohessian") {
		iterate(xx,iv);
	}
	else if (method=="Picard") {success=Iterate_Picard();}
	else {
		success=Iterate_DIIS();
	}
	Sys[0]->CheckResults(report_errors);
	return success;
}

//TODO: this form of function calling loses bounds checking in vectors.
bool Newton::SolveMesodyn(Real* rho) {
	if(debug) cout <<"Solve (mesodyn) in  Newton " << endl;
	int M=Lat[0]->M;
	Real chi;
	int sysmon_length = Sys[0]->SysMolMonList.size();
	int mon_length = In[0]->MonList.size(); //also frozen segments

	bool success=true;
	success=Iterate_DIIS(rho);
	Cp(alpha,xx,iv);
	if (Sys[0]->charged) {
		Sys[0]->DoElectrostatics(alpha+sysmon_length*M,xx+sysmon_length*M);
		Lat[0]->UpdateEE(Sys[0]->EE,Sys[0]->psi,Sys[0]->eps);
	}
	for (int i=0; i<sysmon_length; i++) {
		for (int k=0; k<mon_length; k++) {
                        chi= Sys[0]->CHI[Sys[0]->SysMolMonList[i]*mon_length+k];
			if (chi!=0) {
				PutAlpha(alpha+i*M,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
			}
		}
		if (Sys[0]->charged){
			YplusisCtimesX(alpha+i*M,Sys[0]->EE,Seg[Sys[0]->SysMolMonList[i]]->epsilon,M);
			if (Seg[Sys[0]->SysMolMonList[i]]->valence !=0)
			YplusisCtimesX(alpha+i*M,Sys[0]->psi,-1.0*Seg[Sys[0]->SysMolMonList[i]]->valence,M);
		}
	}


	//Sys[0]->CheckResults(report_errors);
	return success;
}

void Newton::Message(bool e_info, bool s_info, int it, int iterationlimit,Real residual, Real tolerance, string s) {
	if (debug) cout <<"Message in  Newton " << endl;
	if (it == iterationlimit) cout <<"Warning: "<<s<<"iteration not solved. Residual error= " << residual << endl;
	if (e_info || s_info) {
		cout <<" " <<s<<"problem solved" << endl;
		if (e_info) {
			if (it < iterationlimit/10) cout <<" That was easy." << endl;
			if (it > iterationlimit/10 && it < iterationlimit ) cout <<"That will do." << endl;
			if (it <2 && iterationlimit >1 ) cout <<" You hit the nail on the head." << endl;
			if (residual > tolerance) { cout << " Iterations failed." << endl;
				if (residual < tolerance/10) cout <<" I almost made it..." << endl;
			}
		}
		if (s_info) cout <<it << " iterations used to reach residual " << residual << endl;
	}
}


bool Newton::Iterate_Picard() {
if(debug) cout <<"Iterate_Picard in  Newton " << endl;
	int M=Lat[0]->M;
	Real chi;
	alpha=Sys[0]->alpha;
	bool success=true;

	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size();
	it=0;
	if (e_info) cout <<"Picard has been notified" << endl;
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

		Real result; Dot(result,g,g,M);
		residual=residual+result;
		residual=sqrt(residual);
		if(it%i_info == 0){
			printf("it = %i g = %1e \n",it,residual);
		}
		it++;
	}
	Message(e_info,s_info,it,iterationlimit,residual,tolerance,"");
	return success;
}

/*	Direct inversion in iterative subspace. Called from Solve(bool)
*/
bool Newton::Iterate_DIIS() {
if(debug) cout <<"Iterate_DIIS in  Newton " << endl;
	Zero(x0,iv);
	it=0; k_diis=1;
	int k=0;
	// computeG_ext() has been ommented in CheckInput():
	// if (method=="DIIS-ext") ComputeG_ext();

	ComputeG(g); // Or fall back to the classical method.


	YplusisCtimesX(xx,g,delta_max,iv);
	YisAminB(x_x0,xx,x0,iv);
	Cp(xR,xx,iv);
	Dot(residual,g,g,iv);
	residual=sqrt(residual);
	if (e_info) printf("DIIS has been notified\n");
	if (e_info) printf("Your guess = %1e \n",residual);
	while (residual > tolerance && it < iterationlimit) {
		it++;
		Cp(x0,xx,iv);
		ComputeG(g);
		k=it % m; k_diis++; //plek voor laatste opslag
		YplusisCtimesX(xx,g,-delta_max,iv);
		Cp(xR+k*iv,xx,iv); YisAminB(x_x0+k*iv,xx,x0,iv);
		DIIS(xx,x_x0,xR,Aij,Apij,Ci,k,k_diis,m,iv);
		Dot(residual,g,g,iv);
		residual=sqrt(residual);
		if(e_info && it%i_info == 0){
			printf("it = %i g = %1e \n",it,residual);
		}
	}

	Message(e_info,s_info,it,iterationlimit,residual,tolerance,"");
	return it<iterationlimit+1;
}

bool Newton::Iterate_DIIS(Real* rho) {
if(debug) cout <<"Iterate_DIIS for mesodyn in Newton " << endl;
	Zero(x0,iv);
	it=0; k_diis=1;
	int k=0;
	ComputeG_mesodyn(rho);
	YplusisCtimesX(xx,g,delta_max,iv);
	YisAminB(x_x0,xx,x0,iv);
	Cp(xR,xx,iv);
	Dot(residual,g,g,iv);
	residual=sqrt(residual);
	if (e_info) printf("DIIS has been notified\n");
	if (e_info) printf("Your guess = %1e \n",residual);
	while (residual > tolerance && it < iterationlimit) {
		it++;
		Cp(x0,xx,iv);
		ComputeG_mesodyn(rho);
		k=it % m; k_diis++; //plek voor laatste opslag
		YplusisCtimesX(xx,g,-delta_max,iv);
		Cp(xR+k*iv,xx,iv);
		YisAminB(x_x0+k*iv,xx,x0,iv);
		DIIS(xx,x_x0,xR,Aij,Apij,Ci,k,k_diis,m,iv);
		Dot(residual,g,g,iv);
		residual=sqrt(residual);
		if(e_info && it%i_info == 0){
			printf("it = %i g = %1e \n",it,residual);
		}
	}
	Message(e_info,s_info,it,iterationlimit,residual,tolerance,"");
	return it<iterationlimit+1;
}

/* 	The entry point that handles DIIS.
		Called by main if success = false for pseudohessian, Iterate_Picard or iterate_DIIS.
*/
bool Newton::SuperIterate(int search, int target,int ets,int etm) {
if(debug) cout <<"SuperIteration in  Newton " << endl;
	int iv=1;
	int m=10;
	Real* x = (Real*)malloc(iv*sizeof(Real));
	if (ets==-1 && etm==-1) x[0] =Var[search]->GetValue(); else {if (ets>-1) x[0] = Var[ets]->GetValue(); else x[0]=Var[etm]->GetValue();}
	Real* x_x0 = (Real*)malloc(m*sizeof(Real));
	Real* xR = (Real*)malloc(m*sizeof(Real));
	Real* g = (Real*)malloc(iv*sizeof(Real));
	Real* x0 = (Real*)malloc(iv*sizeof(Real)); x0[0]=0;
	Real* Aij = (Real*)malloc(m*m*sizeof(Real));
	Real* Apij = (Real*)malloc(m*m*sizeof(Real));
	Real* Ci = (Real*)malloc(m*sizeof(Real));
	Real delta_max=0.5;
	Real residual;
	Real tol=super_tolerance;
	int iterationlimit=super_iterationlimit;
	int it=0;
	int k_diis=1;
	int k=0;
	if (ets==-1 ||  search<0 || etm>-1) Solve(false); else SuperIterate(search,target,-1,-1);
	if (ets==-1 && etm==-1) g[0]=Var[target]->GetError(); else {if (ets>-1) g[0]=Var[ets]->GetError(); else g[0]=Var[etm]->GetError();}
	YplusisCtimesX(x,g,delta_max,iv);
	YisAminB(x_x0,x,x0,iv);
	Cp(xR,x,iv);
	Dot(residual,g,g,iv);
	residual=sqrt(residual);
	if (super_e_info) printf("Super DIIS has been notified\n");
	if (super_e_info) printf("Your guess = %1e \n",residual);
	while (residual > tol && it < iterationlimit) {
		it++;
		Cp(x0,x,iv);
		if (ets==-1 && etm==-1) Var[search]->PutValue(x[0]); else {if (ets>-1) Var[ets]->PutValue(x[0]); else Var[etm]->PutValue(x[0]);}
		if (ets==-1||search<0 || etm>-1) Solve(false); else SuperIterate(search,target,-1,-1);
		if (ets==-1 && etm==-1) g[0]=Var[target]->GetError(); else {if (ets>-1) g[0]=Var[ets]->GetError(); else g[0]=Var[etm]->GetError();}
		k=it % m; k_diis++; //plek voor laatste opslag
		YplusisCtimesX(x,g,-delta_max,iv);
		Cp(xR+k,x,1); YisAminB(x_x0+k,x,x0,iv);
		DIIS(x,x_x0,xR,Aij,Apij,Ci,k,k_diis,m,iv);
		Dot(residual,g,g,iv);
		residual=sqrt(residual);
		if(super_e_info){
			printf("super-it = %i g = %1e \n",it,residual);
		}
	}
	Message(super_e_info,super_s_info,it,iterationlimit,residual,tol,"super");

	free(x); free(x_x0); free(xR); free(g); free(x0); free(Aij); free(Apij); free(Ci);
	return it<iterationlimit+1;
}


void Newton::ComputePhis() {
if(debug) cout <<"ComputPhis in  Newton " << endl;
	PutU();
	Sys[0]->PrepareForCalculations();
	Sys[0]->ComputePhis();
}

void Newton::ComputeG(Real* g){
 cout <<"ComputeG in Newton " << endl;
	int M=Lat[0]->M;
	Real chi;
	int sysmon_length = Sys[0]->SysMonList.size();
	int mon_length = In[0]->MonList.size(); //also frozen segments
	ComputePhis();

	if (Sys[0]->charged) {
		Sys[0]->DoElectrostatics(g+sysmon_length*M,xx+sysmon_length*M);
		Lat[0]->UpdateEE(Sys[0]->EE,Sys[0]->psi,Sys[0]->eps);
	}

	Cp(g,xx,iv); Zero(alpha,M);
	for (int i=0; i<sysmon_length; i++) {
		for (int k=0; k<mon_length; k++) {
                        chi= Sys[0]->CHI[Sys[0]->SysMonList[i]*mon_length+k];
			if (chi!=0) {
				PutAlpha(g+i*M,Sys[0]->phitot,Seg[k]->phi_side,chi,Seg[k]->phibulk,M);
			}
		}
		if (Sys[0]->charged){
			YplusisCtimesX(g+i*M,Sys[0]->EE,Seg[Sys[0]->SysMonList[i]]->epsilon,M);
			if (Seg[Sys[0]->SysMonList[i]]->valence !=0)
			YplusisCtimesX(g+i*M,Sys[0]->psi,-1.0*Seg[Sys[0]->SysMonList[i]]->valence,M);
		}
	}
	for (int i=0; i<sysmon_length; i++) Add(alpha,g+i*M,M);
	Norm(alpha,1.0/sysmon_length,M);
	for (int i=0; i<sysmon_length; i++) {
		AddG(g+i*M,Sys[0]->phitot,alpha,M);
		Lat[0]->remove_bounds(g+i*M);
		//Times(g+i*M,g+i*M,Sys[0]->KSAM,M);
	}

	if (Sys[0]->charged) {
		Lat[0]->set_bounds(Sys[0]->psi);
		Lat[0]->UpdatePsi(g+sysmon_length*M,Sys[0]->psi,Sys[0]->q,Sys[0]->eps,Sys[0]->psiMask);
		//Lat[0]->remove_bounds(g+sysmon_length*M);


	//	Cp(g+sysmon_length*M,Sys[0]->psi,M);
	//	Lat[0]->UpdatePsi(Sys[0]->psi,Sys[0]->q,Sys[0]->eps);
	//	Lat[0]->set_bounds(Sys[0]->psi);
	//	//Lat[0]->UpdateEE(Sys[0]->EE,Sys[0]->psi,Sys[0]->eps);
	//	//YisAminB(g+sysmon_length*M,g+sysmon_length*M,Sys[0]->psi,M);
	//	YplusisCtimesX(g+sysmon_length*M,Sys[0]->psi,-1.0,M);
	}
}

void Newton::ComputeG_ext(){
if(debug) cout <<"CompueG_ext() in  Newton " << endl;
	int M=Lat[0]->M;
	alpha=Sys[0]->alpha;
	Real chi;
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


/*  Entry point for the mesodyn iteration. Contains the target function.
		Called by Solve(Real* )
*/
void Newton::ComputeG_mesodyn(Real* rho) {
	if (debug) cout << "ComputeG_mesodyn in  Newton " << endl;
	int M = Lat[0]->M;
	//int sysmolmon_length = Sys[0]->SysMolMonList.size();
	ComputePhis();
	Cp(g,rho,iv);
	int i=0; int k=0;
	int length = In[0]->MolList.size();
	while (i<length) {
		int j=0;
		int LENGTH=Mol[i]->MolMonList.size();
		while (j<LENGTH) {YplusisCtimesX(g+k*M,Mol[i]->phi+j*M,-1.0,M); k++; j++;}
		i++;
	}
}
