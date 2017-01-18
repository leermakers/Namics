#include "newttool.h"

double innerproduct(const double *const x, const double *const y,int n) {
	int i=0;
	double sum=0;
	for (i=0; i<n; i++) {
		sum += x[i]*y[i];
	}
	return sum;
}

void multiply(double *const v,double alpha,const float *const h,
		const double *const w, int n) {
	//COMMENT: computes v=(l)(d)(u)*w, ldu are stored in h;
	int i=0,i1=0,j=0;
	double sum=0;
	double *x = new double[n];
	for (i=0; i<n; i++) {
		sum = 0;
		i1 = i-1;
		for (j=i+1; j<n; j++) {
			sum += w[j] * h[i+n*j];
		}
		x[i] = (sum+w[i])*h[i+n*i];
		sum = 0;
		for (j=0; j<=i1; j++) {
			sum += x[j] * h[i+n*j];
		}
		v[i] = alpha*(sum+x[i]);
	}
	delete [] x;
}

// inproduct x*x met een beveiliging....
/*double norm2(const double *const x,int n) {
	int i=0;
	double sum=0,max=0;
	max = DBL_MIN;
	for (i=0; i<n; i++) {
		if ( max < fabs(x[i]) ) {
			sum = 1+sum*pow(max/x[i],2);
			max = fabs(x[i]);
		} else {
			sum += pow(x[i]/max,2);
		}
	}
	return max*sqrt(sum);
}*/
double norm2(const double *const x, int n) {
	int i=0;
	double sum=0;
	for (i=0; i<n; i++) {
		sum += pow(x[i],2);
	}
	return sqrt(sum);
}

// teken van de determinant: als teken negatief geen newtonstap,
// in Newton wordt dan een stap de andere kant op gedaan
int signdeterminant(const float *const h,int n) {
	int i=0,sign=1;
	for (i=0; i<n; i++) {
		if ( h[i+i*n]<0 ) {
			sign = -sign;
		}
	}
	return sign;
}

// zelfde als updatepos: decompositie is LU decompositie: matrix inversie
void updateneg(float *const l,double *const w,int n,double alpha_) {
	//COMMENT See (1) pages 62-63;
	int i=0,i1=0,j=0;
	double dmin=0,sum=0,b=0,d=0,p=0,lji=0,t=0;
	dmin = 1.0/pow(2.0,54);
	alpha_ = sqrt(-alpha_);
	for (i=0; i<=n; i++) {
		i1 = i-1;
		sum = 0;
		for (j=0;j<=i1; j++) {
			sum += l[i+n*j]*w[j];
		}
		w[i] = alpha_*w[i]-sum;
		t += (w[i]/l[i+n*i])*w[i];
	}
	t = 1-t;
	if ( t<dmin ) {
		t = dmin;
	}
	for (i=n-1; i>=0; i--) {
		p = w[i];
		d = l[i+n*i];
		b = d*t;
		t += (p/d)*p;
		l[i+n*i] = b/t;
		b = -p/b;
		for (j=i+1; j<n; j++) {
			lji = l[j+n*i];
			l[j+n*i] = lji+b*w[j];
			w[j] += p*lji;
		}
	}
}

void decompos(float * const h, const int nvar, int ntr) {
//float *h; 	matrix
//int nvar;	dimension of h
//int m;
//int trouble;
//int itrouble;
	int i,j,k;//itr,ntr;
	float sum,lsum,usum,phi,phitr,c,l;
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
		if (phi<0) {
			ntr++;
		}	
		if (phi<phitr) {
			phitr = phi;
		}
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
	return;
}


void updatpos(float * const l, double * const w,
	double * const v,const int nvar, double alpha) {
//float *l; 	lower matrix
//double *w;	update vector
//double *v;	update vector
//int nvar;	dimension of l, w and v
//double alpha;	scaling factor
	int i,j;
	double b,c,d,*wa,*va,vai,waj,vaj;
	float *lai,*laj;
	wa = &w[-1];
	va = &v[-1];

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
	return;
}

void gausa(const float *const l, double *const dup, const double *const g,int nvar) {
//float *l; 	lower matrix
//double *dup;	update vector
//double *g;	update vector
//int nvar;	dimension of l, w and v
	int i,j;
	double*dupa,sum;
	const double *ga;
	const float *lai;

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
	return;
}

void gausb(const float *const du,double *const p, const int nvar) {
//double *du;   lower matrix
//double *p;	update vector
//int nvar;	dimension of l, w and v
	int i,j;
	double *pa,sum;
	const float *duai;

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
	return;
}
