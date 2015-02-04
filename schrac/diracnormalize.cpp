/*! \file schnormalize.cpp
    \brief Sch方程式を解いて得られた波動関数を正規化するクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "diracnormalize.h"

namespace schrac {
    // #region publicメンバ関数

    void SchNormalize::evaluate()
    {
        auto const mp_im1 = pdiffdata_->mp_i_ - 1;

        auto const & li(pdiffdata_->li_);
        auto const & lo(pdiffdata_->lo_);
        auto const & mi(pdiffdata_->mi_);
        auto const & mo(pdiffdata_->mo_);
        auto const & r_mesh_i(pdiffdata_->r_mesh_i_);
        auto const & r_mesh_o(pdiffdata_->r_mesh_o_);

        auto const mpval = pdiffsolver_->getMPval();
        auto const ratio = (std::get<0>(mpval))[0] / (std::get<0>(mpval))[1];

        auto const mp_o = pdiffdata_->mp_o_;
        
        rf_.reserve(pdata_->grid_num_);
        pf_.reserve(pdata_->grid_num_);

        for (auto i = 0; i <= mp_o; i++) {
            auto const h = 1.0 / (2.0 / Data::al + Data::al * pdiffdata_->E_ - Data::al * pdiffdata_->vr_o_[i]);
            auto const dG = std::pow(
                r_mesh_i[i],
                static_cast<double>(pdata_->l_ * (pdata_->l_ + 1)) * lo[i] + mo[i]);
            rf_.push_back(h * (dG + pdata_->kappa_ * schrac::pow(r_mesh_[i], pdata_->l_) * lo[i]));
        	pf_.push_back(r_mesh_[i] * rf_[i]);
        }        

        for (auto i = mp_im1; i >= 0; i--) {
        	r_mesh_.push_back(r_mesh_i[i]);
        	rf_.push_back(schrac::pow(r_mesh_.back(), pdata_->l_) * ratio * li[i]);
        	pf_.push_back(r_mesh_.back() * rf_.back());
        }

        auto const n = 1.0 / std::sqrt(simpson(pf_));
        for (auto i = 0; i < pdata_->grid_num_; i++) {
            rf_[i] *= n;
            pf_[i] *= n;
        }
    }

    // #endregion publicメンバ関数
}
