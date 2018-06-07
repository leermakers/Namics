#include <iostream>
#include "sfnewton.h"

SFNewton::class SFNewton () {

/*      class for
        unconstrained minimization,
        systems of nonlinear algebraic equations,
        curve fitting.

References:

 (1)    "The State of the Art in Numerical Analysis"
        Edited by D.Jacobs. (1977) Academic Press.
 (2)    "Numerical Methods for Unconstrained Optimization"
        Edited by W.Murray. (1972) Academic Press.
 (3)    "Numerical Methods for Constrained Optimization"
        Edited by P.E.Gill and W.Murray. (1974) Academic Press.
 (4)    "Numerical Methods for Non-linear Algebraic Equations"
        Edited by P.Rabinowitz. (1970) Gordon and Breach.
 (5)    "Minimization Subject to Bounds on the Variables"
        P.E.Gill and W.Murray. NPL Report NAC 72 (1976).
 (6)    "Safeguarded Steplength Algorithms for Optimization Using
        Descent Methods" P.E.Gill and W.Murray. NPL Report NAC 37
        (1974).
 (7)    "The Implementation of Two Modified Newton Methods for
        Unconstrained Optimization"
        P.E.Gill, W.Murray, and S.M.Picken. NPL Rpt. NAC 24  (1972)
 (8)    "The Implementation of Two Quasi-Newton Methods for
        Unconstrained Optimization"
        P.E.Gill, W.Murray, and R.A.Pitfield. NPL Rpt NAC 11 (1972)
 (9)    "Conjugate Gradient Methods with Inexact Searches"
        D.F.Shanno. Math. of Operations Res. 3 (1978) 244-256

Author:
Jan Scheutjens (1947-1992), Wageningen Agricultural University, NL.

C Copyright (1980) (1981-1989) Wageningen Agricultural University, NL.

C++ translation:
Peter Barneveld, Wageningen Agricultural University, NL.
Addaptation for namics:
Frans Leermakers, Wageningen Agricultural University, NL.

C Copyright (1993) Wageningen Agricultural University, NL.

 *NO PART OF THIS WORK MAY BE REPRODUCED, EITHER ELECTRONICALLY OF OTHERWISE*

*/
	pseudohessian=samehessian=d_info=e_info=i_info=g_info=h_info=s_info=v_info=x_info=false;
	max_accuracy_for_hessian_scaling = 0.1;

}

SFNewton::~SFNewton() {
}

void SFNewton::residuals(Real*,Real*);
void SFNewton::inneriteration(float*,Real*,Real*,Real,int);
void SFNewton::iterate(Real*,int);

void SFNewton::multiply(Real *v,Real alpha, float *h, Real *w, int nvar) {
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

Real SFNewton::norm2(Real *x, int nvar) {
if(debug) cout <<"norm2 in Newton" << endl;

	Real sum=0;
	for (int i=0; i<nvar; i++) sum += pow(x[i],2);
	return sqrt(sum);
}

int SFNewton::signdeterminant(float *h,int nvar) {
if(debug) cout <<"signdeterminant in Newton" << endl;
	int sign=1;
	for (int i=0; i<nvar; i++) {
		if ( h[i+i*nvar]<0 ) {
			sign = -sign;
		}
	}
	return sign;
}

void SFNewton::updateneg(float *l,Real *w, int nvar, Real alpha) {
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

void SFNewton::decompos(float *h, int nvar, int &ntr) {
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

void SFNewton::updatpos(float *l, Real *w, Real *v, int nvar, Real alpha) {
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

void SFNewton::gausa(float *l, Real *dup, Real *g, int nvar) {
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

void SFNewton::gausb(float *du, Real *p, int nvar) {
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

Real SFNewton::residue(Real *g, Real *p, Real *x, int nvar, Real alpha) {
if(debug) cout <<"residue in Newton " << endl;
	return sqrt(norm2(p,nvar)*norm2(g,nvar)/(1+norm2(x,nvar)));
}

Real SFNewton::linecriterion(Real *g, Real *g0, Real *p, Real *p0, int nvar) {
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

Real SFNewton::newfunction(Real *g, Real *x, int nvar) {
if(debug) cout <<"newfunction in Newton " << endl;
	return pow(norm2(g,nvar),2);
}

void SFNewton::direction(float *h, Real *p, Real *g, Real *g0, Real *x, int nvar, Real alpha){
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
		if ( e_info && t_info) cout << "*";
	}
}

void SFNewton::startderivatives(float *h, Real *g, Real *x, int nvar){
if(debug) cout <<"startderivatives in Newton" << endl;
	float diagonal = 1+norm2(g,nvar);
	H_Zero(h,nvar*nvar);
	for (int i=0; i<nvar; i++) {
		h[i+nvar*i] = diagonal;
	}
}

void SFNewton::resethessian(float *h,Real *g,Real *x,int nvar){
if(debug) cout <<"resethessian in Newton" << endl;
	trouble = 0;
	startderivatives(h,g,x,nvar);
	resetiteration = iterations;
}

void SFNewton::newhessian(float *h, Real *g, Real *g0, Real *x, Real *p, int nvar) {
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
void SFprint_hessian(float* h, int nvar) {
	cout <<"hessian: " << endl;
	for (int i=0; i<nvar; i++)
	for (int j=0; j<nvar; j++)
		if (h[i*nvar+j]!=0) cout <<"i " << i << " j " << j << " h= " << h[i*nvar+j] << endl;
}

void SFNewton::numhessian(float* h,Real* g, Real* x, int nvar) {
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

void SFNewton::decomposition(float *h,int nvar, int &trouble){
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

void SFNewton::findhessian(float *h, Real *g, Real *x,int nvar) {
if(debug) cout <<"findhessian in Newton" << endl;
	if ( !samehessian ) {
		if ( iterations==0 ) resethessian(h,g,x,nvar);
		numhessian(h,g,x,nvar); // passes through residuals so check pseudohessian
		if (!pseudohessian) {
			decomposition(h,nvar,trouble);
		}
	}
}

void SFNewton::newdirection(float *h, Real *p, Real *p0, Real *g, Real *g0, Real *x, int nvar, Real alphabound) {
if(debug) cout <<"newdirection in Newton" << endl;
	inneriteration(h,g,x,accuracy,nvar);
	memcpy(p0, p, sizeof(*p0)*nvar);
	direction(h,p,g,g0,x,nvar,ALPHA);
	accuracy = residue(g,p,x,nvar,ALPHA);
}

void SFNewton::newtrustregion(Real *g, Real *g0, Real *p, Real *p0, int nvar){
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

Real SFNewton::linesearch(Real *g, Real *g0, Real *p, Real *x, Real *x0, int nvar, Real alphabound) {
if(debug) cout <<"linesearch in Newton" << endl;
	Real newalpha = alphabound<1 ? alphabound : 1;
	newalpha = zero(g,g0,p,x,x0,nvar,newalpha);
	return newalpha;
}

Real SFNewton::zero(Real *g, Real *g0, Real *p, Real *x, Real *x0, int nvar, Real newalpha) {
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

Real SFNewton::stepchange(Real *g, Real *g0, Real *p, Real *p0, Real *x, Real *x0, int nvar, Real &alpha){
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

void SFNewton::COMPUTEG(Real *x, Real *g, int nvar) {
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

void SFNewton::ResetX(Real* x,int nvar) {
	int pos=nvar;
	for (int i=iv-1; i>=0; i--) {
		if (mask[i]==1) {
			pos--;
			x[i]=x[pos];
		} else x[i]=0;
	}
}

void SFNewton::iterate(Real *x,int nvar) {
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
	printf("it =  %i  E = %e |g| = %e alpha = %e \n",it,accuracy,normg,ALPHA);
	Message(e_info,s_info,it,iterationlimit,accuracy,tolerance,"");
	ResetX(xx,nvar);
	delete [] p; delete [] g0 ; delete [] p0;
	delete [] h; delete [] reverseDirection;
}

string SFNewton::GetNewtonInfo(int &IV) {
	IV=iv;
	return method;
}

void SFNewton::Message(bool e_info, bool s_info, int it, int iterationlimit,Real residual, Real tolerance, string s) {
	if (debug) cout <<"Message in  Newton " << endl;
	if (it == iterationlimit) cout <<"Warning: "<<s<<"iteration not solved. Residual error= " << residual << endl;
	if (e_info || s_info) {
		
		cout <<s<<"Problem solved." << endl;
		if (e_info) {
			if (it < iterationlimit/10) cout <<"That was easy." << endl;
			if (it > iterationlimit/10 && it < iterationlimit ) cout <<"That will do." << endl;
			if (it <2 && iterationlimit >1 ) cout <<"You hit the nail on the head." << endl;
			if (residual > tolerance) { cout << " Iterations failed." << endl;
				if (residual < tolerance/10) cout <<"... I almost made it..." << endl;
			}
		}
		if (s_info) cout <<it << " iterations used to reach residual " << residual <<"."<< endl;
	}
}


