/*! \file vhartree.h
    \brief Hartreeポテンシャルを求めるクラスの実装

    Copyright © 2015 @dc1394 All Rights Reserved.
*/

#include "vhartree.h"
#include <array>
#include <boost/assert.hpp>     // for BOOST_ASSERT
#include <boost/cast.hpp>       // for boost::numeric_cast

namespace schrac {
    // #region コンストラクタ
    
    Vhartree::Vhartree(std::vector<double> const & r_mesh) :
        Vhart([this]{ return vhart_; }, [this](std::vector<double> const & v) { return vhart_ = v; }),
        acc_(gsl_interp_accel_alloc(), gsl_interp_accel_deleter),
        r_mesh_(r_mesh),
        spline_(gsl_spline_alloc(gsl_interp_cspline, r_mesh.size()), gsl_spline_deleter)
    {
        vhart_.reserve(r_mesh.size());
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数

    double Vhartree::dvhartree_dr(double r) const
    {
        return gsl_spline_eval_deriv(spline_.get(), r, acc_.get());
    }

    void Vhartree::set_vhartree_boundary_condition(double Z)
    {
        auto const size = boost::numeric_cast<std::int32_t>(vhart_.size());
        for (auto i = 0; i < size; i++) {
            vhart_[i] += (Z / r_mesh_.back() - vhart_.back());
        }
    }

    void Vhartree::vhart_init()
    {
        gsl_spline_init(spline_.get(), r_mesh_.data(), vhart_.data(), r_mesh_.size());
    }

    double Vhartree::vhartree(double r) const
    {
        return gsl_spline_eval(spline_.get(), r, acc_.get());
    }

    // #endregion publicメンバ関数
}