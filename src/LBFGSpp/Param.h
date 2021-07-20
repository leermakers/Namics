// Copyright (C) 2016-2021 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef PARAM_H
#define PARAM_H

#include <Eigen/Core>
#include <stdexcept>  // std::invalid_argument


namespace LBFGSpp {

template <typename Scalar = double>
class LBFGSParam
{
public:
    ///
    /// The number of corrections to approximate the inverse Hessian matrix.
    /// The L-BFGS routine stores the computation results of previous \ref m
    /// iterations to approximate the inverse Hessian matrix of the current
    /// iteration. This parameter controls the size of the limited memories
    /// (corrections). The default value is \c 6. Values less than \c 3 are
    /// not recommended. Large values will result in excessive computing time.
    ///
    int    m;
    ///
    /// Absolute tolerance for convergence test.
    /// This parameter determines the absolute accuracy \f$\epsilon_{abs}\f$
    /// with which the solution is to be found. A minimization terminates when
    /// \f$||g|| < \max\{\epsilon_{abs}, \epsilon_{rel}||x||\}\f$,
    /// where \f$||\cdot||\f$ denotes the Euclidean (L2) norm. The default value is
    /// \c 1e-5.
    ///
    Scalar epsilon;
    ///
    ///
    /// The maximum number of iterations.
    /// The optimization process is terminated when the iteration count
    /// exceeds this parameter. Setting this parameter to zero continues an
    /// optimization process until a convergence or error. The default value
    /// is \c 0.
    ///
    int    max_iterations;
    Scalar delta_max;
    Scalar delta_min;
    bool e_info;
    int i_info;

public:
    LBFGSParam()
    {
        m              = 6;
        epsilon        = Scalar(1e-5);
        max_iterations = 100;
	 delta_max = Scalar(0.5);
	 delta_min = Scalar(1e-4);
	 e_info = true;
	 i_info =1;  	 
    }
    
    inline void check_param() const
    {
        if (m <= 0)
           throw std::invalid_argument("'m' must be positive");
        if(epsilon < 0)
            throw std::invalid_argument("'epsilon' must be non-negative");
        if(max_iterations < 0)
            throw std::invalid_argument("'max_iterations' must be non-negative");
	 if (delta_max <0) 
		throw std::invalid_argument("delta_max must be positive");
        if (delta_min <0) 
		throw std::invalid_argument("delta_min must be positive");
    }
};


} // namespace LBFGSpp

#endif // PARAM_H
