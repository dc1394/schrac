﻿/*! \file rho.h
    \brief 電子密度ρ(r)を求めるクラスの実装

    Copyright © 2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "diffdata.h"
#include "rho.h"
#include <cmath>    // for std::pow

namespace schrac {
    // #region コンストラクタ

    Rho::Rho(std::shared_ptr<DiffData> const & pdiffdata) :
        PRho([this] { return std::cref(rho_); }, nullptr),
        acc_(gsl_interp_accel_alloc(), gsl_interp_accel_free),
        pdiffdata_(pdiffdata),
        spline_(gsl_spline_alloc(gsl_interp_cspline, pdiffdata->r_mesh_.size()), gsl_spline_free)
    {
        auto const & pdata = pdiffdata_->pdata_;
        if (pdata->chemical_symbol_ == Data::Chemical_Symbol[0]) {
            rho_.resize(pdata->grid_num_ + 1);

            return;
        }
        else if (pdata->rho0_c_ == std::nullopt || pdata->rho0_alpha_ == std::nullopt) {
            // デフォルト値を代入
            auto const w = 2.0;
            pdata->rho0_c_ = std::make_optional<double>(std::pow(w, 4) / 16.0);
            pdata->rho0_alpha_ = std::make_optional<double>(0.5 * w);
            *pdata->rho0_c_ *= (pdata->Z_ / w);
        }

        rho_.reserve(pdata->grid_num_ + 1);
        for (auto i = 0; i <= pdata->grid_num_; i++) {
            rho_.push_back(*pdata->rho0_c_ * std::exp(-*pdata->rho0_alpha_ * pdiffdata->r_mesh_[i]));
        }
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数

    void Rho::init()
    {
        gsl_spline_init(spline_.get(), pdiffdata_->r_mesh_.data(), rho_.data(), pdiffdata_->r_mesh_.size());
    }

    double Rho::operator()(double r) const
    {
        return gsl_spline_eval(spline_.get(), r, acc_.get());
    }

    void Rho::rhomix(dvector const & newrho)
    {
        auto const & pdata = pdiffdata_->pdata_;
        for (auto i = 0; i <= pdata->grid_num_; i++) {
            rho_[i] = (1.0 - pdata->scf_mixing_weight_) * rho_[i] +
                pdata->scf_mixing_weight_ * newrho[i];
        }
    }

    // #endregion publicメンバ関数
}
