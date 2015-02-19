/*! \file rho.h
    \brief 電子密度ρ(r)を求めるクラスの実装

    Copyright © 2015 @dc1394 All Rights Reserved.
*/

#include "rho.h"
#include <boost/assert.hpp>

namespace schrac {
    // #region コンストラクタ

    Rho::Rho(std::vector<double> const & r_mesh, std::vector<double> const & rho) :
        acc_(gsl_interp_accel_alloc(), gsl_interp_accel_deleter),
        spline_(gsl_spline_alloc(gsl_interp_cspline, r_mesh.size()), gsl_spline_deleter),
        rho_(rho)
    {
        BOOST_ASSERT(r_mesh.size() == rho.size());
        
        gsl_spline_init(spline_.get(), r_mesh.data(), rho_.data(), r_mesh.size());
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数

    double Rho::operator()(double r) const
    {
        return gsl_spline_eval(spline_.get(), r, acc_.get());
    }

    // #endregion publicメンバ関数
}
