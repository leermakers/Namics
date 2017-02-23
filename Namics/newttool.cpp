#include "newttool.h"
#include "tools.h"
#include <stdio.h>
void multiply(double *v, double alpha, float *h, double *w, int nvar, int m) {
	int i0,i1,nn;
	double sum=0;
	double *va, *wa, *xa;
	float* hai;
	float* ha = &h[-1];
	va = &v[-1];
	wa = &w[-1];
	double *x = new double[nvar];
	xa = &x[-1];
	for (int i=1; i<=nvar; i++) {
               sum = 0;
               if (i>m) i0=i-m; else i0=0;
               i1 = i-1;
               nn=i+m-1; if (nn>nvar) nn=nvar;
               hai=&ha[(i-1)*nvar];
               for (int j=i+1; j<=nn; j++) {
                       sum+= wa[j] * hai[j-i0];
               }
               xa[i] = (sum+wa[i])*hai[i-i0];
               sum = 0;
               for (int j=i0+1; j<=i1; j++) {
                       sum += xa[j] * hai[j-i0];
               }
               va[i] = alpha*(sum+xa[i]);
       }
       delete [] x;
}

void multiply(double *v,double alpha, float *h, double *w, int nvar) {
	int i=0,i1=0,j=0;
	double sum=0;
	double *x = new double[nvar];
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

double norm2(double *x, int nvar) {
	double sum=0;
	for (int i=0; i<nvar; i++) sum += pow(x[i],2);
	return sqrt(sum);
}

int signdeterminant(float *h, int nvar, int m) {
	int sign=1; int i0;
	float *ha, *hai;
	ha = &h[-1];
	for (int i=1; i<=nvar; i++) {
		if (i<m) i0=i; else i0=m;
		hai=&ha[(i-1)*nvar];
		if (hai[i0]<0) sign=-sign;
	}
	return sign;
}

int signdeterminant(float *h,int nvar) {
	int sign=1;
	for (int i=0; i<nvar; i++) {
		if ( h[i+i*nvar]<0 ) {
			sign = -sign;
		}
	}
	return sign;
}

void updateneg(float *l, double *w, int nvar, int m, double alpha) {
	int i1=0,i0=0,nn=0,j0=0;
	double dmin=0,sum=0,b=0,d=0,p=0,lji=0,t=0;
	float *lai,*laj;
	float *la = &l[-1];
	double *wa = &w[-1]; 
	dmin = 1.0/pow(2.0,54);
	alpha = sqrt(-alpha); //is it sure that alpha is negative in the argument?
	for (int i=1; i<=nvar; i++) {
		if (i>m) i0=i-m; else i0=0;
		i1 = i-1;
		sum = 0;
		lai = &la[(i-1)*nvar];
		for (int j=i0+1; j<=i1; j++) sum += lai[j-i0]*wa[j];
		wa[i] = alpha*wa[i]-sum;
		t += (wa[i]/lai[i-i0])*wa[i];
	}
	t = 1-t;
	if ( t<dmin ) t = dmin;
	for (int i=nvar; i>=1; i--) {
		if (i>m) i0=i-m; else i0=0;
		nn=i+m-1;
		if (nn>nvar) nn=nvar;
		p = wa[i];
		lai=&la[(i-1)*nvar];
		d = lai[i-i0];
		b = d*t;
		t += (p/d)*p;
		d=b/t; lai[i-i0]=d; 
		b = -p/b;
		for (int j=i+1; j<=nn; j++) {
			if (j>m) j0=j-m; else j0=0;
			laj = &la[(j-1)*nvar]; 
			lji = laj[i-j0];
			laj[i-j0] = lji+b*wa[j];
			wa[j] += p*lji;
		}
	}
}

void updateneg(float *l,double *w, int nvar, double alpha) {
	int i=0,i1=0,j=0;
	double dmin=0,sum=0,b=0,d=0,p=0,lji=0,t=0;
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

void decompos(float *h, int nvar, int m, int &ntr) {
	int i0,i1,nn,j0,k0;
	double sum,lsum,usum,phi,phitr,c,l;
	float *ha,*hai,*haj,*hak;
	ha = &h[-1];
	phitr = FLT_MAX;
	ntr = 0;
	for (int i=1; i<=nvar; i++) {
		hai = &ha[(i-1)*nvar];
		sum = 0;
		i0=i-m; if (i0<0) i0=0; 
		i1=i-1;
		nn=m+i-1; if (nn>nvar) nn=nvar;
		for (int j=i0+1; j<=nn; j++) {
			j0=j-m; if (j0<0) j0=0;
			haj=&ha[(j-1)*nvar];
			if (j<i) {
				c=hai[j-i0]; l=c/haj[j-j0]; hai[j-i0]=l;
				c=haj[i-j0];
				haj[i-j0]=c/haj[j-j0];
				sum=sum+l*c;
			} else {
				if (j==i) {
					phi=hai[i-i0]-sum;
					hai[i-i0]=phi;
					if (phi<0) ntr++;
					if (phi<phitr) phitr = phi; 
				} else {
					lsum=usum=0;
					for (int k=j0+1; k<=i1; k++) {
						k0=k-m; if (k0<0) k0=0;
						hak=&ha[(k-1)*nvar];
						lsum+=haj[k-j0]*hak[i-k0];
						usum+=hai[k-i0]*hak[j-k0];
					}	
					haj[i-j0]-=lsum;
					hai[j-i0]-=usum;
				}
			}
		}
	}
}

void decompos(float *h, int nvar, int &ntr) {
	int i,j,k;//itr,ntr;
	double sum,lsum,usum,phi,phitr,c,l;
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

void updatpos(float *l, double *w, double *v, int nvar, int m, double alpha) {
	double b,c,d,p,q,wj,vj,lj;
	int i0=0,j0=0,nn=0;
	float *lai,*laj;
	float *la = &l[-1];
	double *wa = &w[-1];
	double *va = &v[-1];
	for (int i=1; i<=nvar; i++) {
		if (i>m) i0=i-m; else i0=0;
		nn=i+m-1; if (nn>nvar) nn=nvar;
		p=wa[i]; 
		q=va[i];
		lai= &la[(i-1)*nvar]; 
		d = lai[i-i0];
		b = d+(alpha*p)*q;
		lai[i-i0] = b;
		d /= b;
		c = q*alpha/b;
		b = p*alpha/b;
		alpha *= d;
		for (int j=i+1; j<=nn; j++) {
			if (j>m) j0=j-m; else j0=0;
			wj = wa[j];
			laj=&la[(j-1)*nvar]; lj=laj[i-j0];
			wa[j]=wj-p*lj;
			laj[i-j0]=lj*d+c*wj; 
			vj=va[j]; lj=lai[j-i0]; va[j]=vj-q*lj;
			lai[j-i0] = lj*d+b*vj;
		}
	}
}

void updatpos(float *l, double *w, double *v, int nvar, double alpha) {
	int i,j;
	double b,c,d;
	double vai,waj,vaj;
	float *lai,*laj;
	double * wa = &w[-1];
	double * va = &v[-1];
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

void gausa(float *l, double *dup, double *g, int nvar, int m) {
	int i0=0,i1=0;
	double sum;
	double *dupa;
	float *lai;

	dupa = &dup[-1];
	double* ga = &g[-1];
	float* la = &l[-1];

	for (int i=1; i<=nvar; i++) {
		if (i>m) i0=i-m; else i0=0;
		sum=0;
		i1=i-1; 
		lai = &la[(i-1)*nvar]; 
		for (int j=i0+1; j<=i1; j++){
			sum += lai[j-i0]*dupa[j];
		}
		dupa[i] = - ga[i] - sum;
	}
}

void gausa(float *l, double *dup, double *g, int nvar) {
	int i,j;
	double sum;
	double*dupa;
	double *ga;
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

void gausb(float *du, double *p, int nvar, int m) {
	int i0=0,im=0,nn=0;
	double sum;
	float *duai;
	double *pa = &p[-1];
	float *dua= &du[-1]; 
	for (int i=nvar; i>=1; i--) {
		nn=m+i-1; if (nn>nvar) nn=nvar;
		sum = 0;
		duai = &dua[(i-1)*nvar];
		for (int j=i+1; i<=nn; i++) {
			if (i>m) i0=j-i+m; else i0=j;
			sum += duai[i0]*pa[j];
		}
		if (i>m) im=m; else im=i;
		pa[i] = pa[i]/duai[im] - sum;
	}
}

void gausb(float *du, double *p, int nvar) {
	int i,j;
	double sum;
	double *pa;
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

double residue(double *g, double *p, double *x, int nvar, double alpha) {
	return sqrt(norm2(p,nvar)*norm2(g,nvar)/(1+norm2(x,nvar)));
}

double linecriterion(double *g, double *g0, double *p, double *p0, int nvar) {
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

double newfunction(double *g, double *x, int nvar) {
	return pow(norm2(g,nvar),2);
}





