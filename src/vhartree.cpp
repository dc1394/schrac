/*! \file vhartree.h
    \brief Hartreeポテンシャルを求めるクラスの実装

    Copyright © 2015 @dc1394 All Rights Reserved.
*/

#include "vhartree.h"
#include <array>
#include <boost/assert.hpp>

namespace schrac {
    // #region コンストラクタ
    
    Vhartree::Vhartree(std::vector<double> const & r_mesh) :
        PVhart([this]{ return vhart_; }, nullptr),
        acc_(gsl_interp_accel_alloc(), gsl_interp_accel_deleter),
        r_mesh_(r_mesh),
        spline_(gsl_spline_alloc(gsl_interp_cspline, r_mesh.size()), gsl_spline_deleter)
    {
        vhart_.reserve(r_mesh.size());
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数

    double Vhartree::operator()(double r) const
    {
        return gsl_spline_eval(spline_.get(), r, acc_.get());
    }

    void Vhartree::vhart_init()
    {
        gsl_spline_init(spline_.get(), r_mesh_.data(), vhart_.data(), r_mesh_.size());
    }

    // #endregion publicメンバ関数
}