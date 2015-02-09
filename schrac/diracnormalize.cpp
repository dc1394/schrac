/*! \file diracnormalize.cpp
    \brief Dirac方程式を解いて得られた波動関数を正規化するクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "diracnormalize.h"
#include "simpson.h"
#include <utility>          // for std::move

namespace schrac {
    // #region publicメンバ関数

    void DiracNormalize::evaluate()
    {
        auto const mp_im1 = pdiffdata_->mp_i_ - 1;

        auto const & li(pdiffdata_->li_);
        auto const & lo(pdiffdata_->lo_);
        auto const & mi(pdiffdata_->mi_);
        auto const & mo(pdiffdata_->mo_);
        auto const & r_mesh_o(pdiffdata_->r_mesh_o_);

        auto const mpval = pdiffsolver_->getMPval();
        auto const ratio = (std::get<0>(mpval))[0] / (std::get<0>(mpval))[1];

        auto const mp_o = pdiffdata_->mp_o_;
        
        rf_.reserve(pdata_->grid_num_);
        pf_large_.reserve(pdata_->grid_num_);
        pf_small_.reserve(pdata_->grid_num_);

        for (auto i = 0; i <= mp_o; i++) {
            rf_.push_back(std::pow(r_mesh_[i], pdata_->l_) * lo[i]);

            auto const h = 1.0 / (2.0 / Data::al + Data::al * pdiffdata_->E_ - Data::al * pdiffdata_->vr_o_[i]);
            auto const dG = std::pow(
                r_mesh_[i],
                static_cast<double>(pdata_->l_ * (pdata_->l_ + 1)) * lo[i] + mo[i]);

            pf_large_.push_back(r_mesh_[i] * rf_[i]);
            pf_small_.push_back(h * (dG + pdata_->kappa_ * std::pow(r_mesh_[i], pdata_->l_) * lo[i]));
        }        

        for (auto i = mp_im1; i >= 0; i--) {
        	r_mesh_.push_back(r_mesh_o[i]);
        	rf_.push_back(std::pow(r_mesh_.back(), pdata_->l_) * ratio * li[i]);

            auto const h = 1.0 / (2.0 / Data::al + Data::al * pdiffdata_->E_ - Data::al * pdiffdata_->vr_o_[i]);
            auto const dG = std::pow(
                r_mesh_.back(),
                static_cast<double>(pdata_->l_ * (pdata_->l_ + 1)) * li[i] + mi[i]);

            pf_large_.push_back(r_mesh_.back() * rf_[i]);
            pf_small_.push_back(h * (dG + pdata_->kappa_ * std::pow(r_mesh_.back(), pdata_->l_) * li[i]));
        }

        normalize();
    }

    Normalize<DiracNormalize>::myhash DiracNormalize::getresult() const
    {
        myhash hash;
        hash.insert(std::make_pair("Mesh (r)", std::move(r_mesh_)));
        hash.insert(std::make_pair("Eigen function", std::move(rf_)));
        hash.insert(std::make_pair("Eigen function large (multiply r)", std::move(pf_large_)));
        hash.insert(std::make_pair("Eigen function small (multiply r)", std::move(pf_small_)));

        return std::move(hash);
    }

    void DiracNormalize::normalize()
    {
        Simpson simpson(pdiffdata_->dx_);
        auto const n = 1.0 / std::sqrt(simpson(pf_large_, r_mesh_) + simpson(pf_small_, r_mesh_));
        for (auto i = 0; i < pdata_->grid_num_; i++) {
            rf_[i] *= n;
            pf_large_[i] *= n;
            pf_small_[i] *= n;
        }
    }

    // #endregion publicメンバ関数
}
