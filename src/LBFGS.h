// Addopted from LBFGSsolverby:
// Copyright (C) 2016-2021 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef LBFGS_H
#define LBFGS_H

#include <Eigen/Core>
#include "LBFGSpp/Param.h"
#include "LBFGSpp/BFGSMat.h"
#include <iostream>

namespace LBFGSpp {

template < typename Scalar>
class LBFGSSolver
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Map<Vector> MapVec;

    LBFGSParam<Scalar>& m_param;  // Parameters to control the LBFGS algorithm
    BFGSMat<Scalar>           m_bfgs;   // Approximation to the Hessian matrix
    Vector                    m_xp;     // Old x
    Vector                    m_grad;   // New gradient
    Vector                    m_gradp;  // Old gradient
    Vector                    m_drt;    // Moving direction
    Scalar alphabound 	  = Scalar(0);
    Scalar trustregion    = 1;
    Scalar trustfactor    = 1;
    Scalar ALPHA 	  = Scalar(1);


    inline void reset(int n) {
        int m = m_param.m;
        m_bfgs.reset(n, m);
        m_xp.resize(n);
        m_grad.resize(n);
        m_gradp.resize(n);
        m_drt.resize(n);
    }


	Scalar stepchange(Vector g, Vector g0, int nvar) {
		Scalar normg,normg0,gg0;
		normg0= g0.norm();
		gg0=g.dot(g0);
		gg0=gg0/normg0/normg0;
		normg=g.norm()/normg0;
		normg=normg*normg;
		if ((gg0>1||normg>1) && normg-gg0*fabs(gg0) < 0.2) normg=1.5*normg;
		if (gg0<0&& normg<2) return 1;
		else if (normg>10) return 0.01;
		else { //most of the time....
			return 0.4*(1+0.75*normg)/(normg-gg0*fabs(gg0)+0.1);
		}
	}

    void newtrustregion(Scalar normp,Scalar delta_max,Scalar delta_min) {
		if (normp>0 &&trustregion>2*ALPHA*normp) trustregion = 2*ALPHA*normp;
		trustregion *=trustfactor;
		trustfactor =1;
		if (trustregion > delta_max) trustregion = delta_max;
		if (trustregion < delta_min) trustregion = delta_min;
    }

public:
    LBFGSSolver(LBFGSParam<Scalar>& param) :
       m_param(param)
    {
        m_param.check_param();
    }

    template <typename Foo>
    inline int minimize(Foo& residuals, Vector& x, Scalar& gnorm)
    {

	 Scalar pnorm;
	 Scalar delta_max=m_param.delta_max;
	 Scalar delta_min=m_param.delta_min;
	 trustregion=delta_max;
        const int n = x.size();
        reset(n);

        gnorm = residuals(x, m_grad);

        if(gnorm <= m_param.epsilon) {
		cout <<"You hit the nail on its head" <<endl;
	   	return 1;
	 }

        m_drt.noalias() = -m_grad;
        int k = 1;
        for( ; ; )
        {
            m_xp.noalias() = x;  //x0=x;
            m_gradp.noalias() = m_grad; //g0=g;

	     pnorm=m_drt.norm();
	     newtrustregion(pnorm,delta_max,delta_min);
	     alphabound = trustregion/pnorm;
	     if (alphabound>1) ALPHA=1; else ALPHA=alphabound;

	     x.noalias() = m_xp + ALPHA * m_drt;
	     gnorm=residuals(x, m_grad);

	     trustfactor *=stepchange(m_grad,m_gradp,n);
	     trustfactor *=ALPHA/alphabound;
            if (ALPHA < m_param.epsilon) ALPHA *=2;
	     if (m_param.e_info) {
		if (k==1) cout <<"Your Guess: "<<gnorm << endl;
		if (k%m_param.i_info==0 && k>1) {
			cout << "i = " << k << " |g| = " << gnorm << " alpha = " << ALPHA << endl;
		}
	     }

            if(gnorm <= m_param.epsilon  || k >= m_param.max_iterations){
                return k;
            }

            m_bfgs.add_correction(x - m_xp, m_grad - m_gradp);
            m_bfgs.apply_Hv(m_grad, -Scalar(1), m_drt);

            k++;
        }

        return k;
    }

};

} // namespace LBFGSpp

#endif // LBFGS_H
