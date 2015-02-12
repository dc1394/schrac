/*! \file vhartree.h
    \brief Hartreeポテンシャルを求めるクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "vhartree.h"
#include <array>
#include <boost/numeric/odeint.hpp>     // for boost::numeric::odeint

namespace schrac {
    // #region コンストラクタ
    
    Vhartree::Vhartree(std::shared_ptr<Data> const & pdata, std::vector<double> const & rho, std::vector<double> const & r_mesh) :
        acc(gsl_interp_accel_alloc()),
        spline(gsl_spline_alloc(gsl_interp_cspline, pdata->grid_num_ + 1))
    {
        using namespace boost::numeric::odeint;

        std::array<double, 2> state = { 0.0, 1.0 };

        
        integrate_adaptive()
            stepper,
            [this](myarray const & f, myarray & dfdx, double x) { return derivs(f, dfdx, x); },
            state,
            pdiffdata_->x_i_[0],
            pdiffdata_->x_i_[pdiffdata_->mp_i_] - pdiffdata_->dx_,
            -pdiffdata_->dx_);
    }

}