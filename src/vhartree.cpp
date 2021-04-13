/*! \file vhartree.h
    \brief Hartreeポテンシャルを求めるクラスの実装

    Copyright © 2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "vhartree.h"
#include <array>                // for std::array
#include <boost/assert.hpp>     // for BOOST_ASSERT

namespace schrac {
    // #region コンストラクタ
    
    Vhartree::Vhartree(std::vector<double> const & r_mesh) :
        Vhart([this]{ return std::cref(vhart_); }, [this](std::vector<double> const & v) { return vhart_ = v; }),
        acc_(gsl_interp_accel_alloc(), gsl_interp_accel_free),
        r_mesh_(r_mesh),
        spline_(gsl_spline_alloc(gsl_interp_cspline, r_mesh.size()), gsl_spline_free)
    {
        vhart_.reserve(r_mesh.size());
    }

    Vhartree::Vhartree(Vhartree const & rhs) :
        Vhartree(rhs.r_mesh_)
    {
        vhart_ = rhs.vhart_;
        gsl_spline_init(spline_.get(), r_mesh_.data(), vhart_.data(), r_mesh_.size());
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数

    double Vhartree::dvhartree_dr(double r) const
    {
        return gsl_spline_eval_deriv(spline_.get(), r, acc_.get());
    }

    void Vhartree::set_vhartree_boundary_condition(double Z)
    {
        for (auto && v : vhart_) {
            v += (Z / r_mesh_.back() - vhart_.back());
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
