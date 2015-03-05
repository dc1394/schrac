/*! \file diracnormalize.cpp
    \brief Dirac方程式を解いて得られた波動関数を正規化するクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "diracnormalize.h"
#include "simpson.h"
#include <utility>          // for std::move

namespace schrac {
    // #region publicメンバ関数

    void DiracNormalize::evaluate(boost::optional<std::vector<double>> const & prho)
    {
        auto const mp_im1 = pdiffdata_->mp_i_ - 1;

        auto const & li(pdiffdata_->li_);
        auto const & lo(pdiffdata_->lo_);
        auto const & mi(pdiffdata_->mi_);
        auto const & mo(pdiffdata_->mo_);
        auto const & r_mesh_i(pdiffdata_->r_mesh_i_);

        auto const mpval = pdiffsolver_->getMPval();
        auto const ratio = (std::get<0>(mpval))[0] / (std::get<0>(mpval))[1];

        auto const mp_o = pdiffdata_->mp_o_;
        
        rf_.reserve(pdata_->grid_num_ + 1);
        pf_large_.reserve(pdata_->grid_num_ + 1);
        pf_small_.reserve(pdata_->grid_num_ + 1);

        for (auto i = 0; i <= mp_o; i++) {
            rf_.push_back(std::pow(pdiffdata_->r_mesh_[i], pdata_->l_) * lo[i]);
            pf_large_.push_back(pdiffdata_->r_mesh_[i] * rf_[i]);

            auto const h = 1.0 /
                (2.0 / Data::al + Data::al * pdiffdata_->E_ - Data::al * pdiffsolver_->V_(pdiffdata_->r_mesh_[i]));
            auto const dG = std::pow(
                pdiffdata_->r_mesh_[i],
                static_cast<double>(pdata_->l_ * (pdata_->l_ + 1)) * lo[i] + mo[i]);

            pf_small_.push_back(h * (dG + pdata_->kappa_ * std::pow(pdiffdata_->r_mesh_[i], pdata_->l_) * lo[i]));
        }        

        for (auto i = mp_im1; i >= 0; i--) {
            rf_.push_back(std::pow(r_mesh_i[i], pdata_->l_) * ratio * li[i]);
            pf_large_.push_back(r_mesh_i[i] * rf_.back());
            
            auto const h = 1.0 / (2.0 / Data::al + Data::al * pdiffdata_->E_ - Data::al * pdiffsolver_->V_(r_mesh_i[i]));
            auto const dG = ratio * std::pow(r_mesh_i[i], static_cast<double>(pdata_->l_)) *
                (static_cast<double>(pdata_->l_ + 1) * li[i] + mi[i]);

            pf_small_.push_back(h * (dG + pdata_->kappa_ * std::pow(r_mesh_i[i], pdata_->l_) * ratio * li[i]));
        }

        normalize();

        if (prho) {
            rho_.reserve(pdata_->grid_num_ + 1);
            for (auto i = 0; i <= pdata_->grid_num_; i++) {
                rho_.push_back(sqr(pdiffdata_->r_mesh_[i]) * (*prho)[i]);
            }
        }
        else {
            rho_.reserve(pdata_->grid_num_ + 1);
            for (auto i = 0; i <= pdata_->grid_num_; i++) {
                rho_.push_back(sqr(pf_large_[i]) + sqr(pf_small_[i]));
            }
        }
    }

    Normalize<DiracNormalize>::mymap DiracNormalize::getresult() const
    {
        Normalize<DiracNormalize>::mymap wf;
        wf["1 Mesh (r)"] = pdiffdata_->r_mesh_;
        wf["2 Eigen function"] = std::move(rf_);
        wf["3 Rho (multiply 4 * pi * r ** 2)"] = std::move(rho_);
        wf["4 Eigen function large (multiply r)"] = std::move(pf_large_);
        wf["5 Eigen function small (multiply r)"] = std::move(pf_small_);

        return std::move(wf);
    }

    void DiracNormalize::normalize()
    {
        Simpson simpson(pdiffdata_->dx_);
        auto const n = 1.0 / 
            std::sqrt(simpson(pf_large_, pdiffdata_->r_mesh_) + simpson(pf_small_, pdiffdata_->r_mesh_));
        for (auto i = 0; i <= pdata_->grid_num_; i++) {
            rf_[i] *= n;
            pf_large_[i] *= n;
            pf_small_[i] *= n;
        }
    }

    // #endregion publicメンバ関数
}

