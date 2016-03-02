/*! \file schnormalize.cpp
    \brief Sch方程式を解いて得られた波動関数を正規化するクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "schnormalize.h"
#include "simpson.h"
#include <utility>                              // for std::move

namespace schrac {
    // #region publicメンバ関数

    void SchNormalize::evaluate()
    {
        auto const mp_im1 = pdiffdata_->mp_i_ - 1;

        auto const & li(pdiffdata_->li_);
        auto const & lo(pdiffdata_->lo_);
        auto const & r_mesh_i(pdiffdata_->r_mesh_i_);

        auto const mpval = pdiffsolver_->getMPval();
        auto const ratio = (std::get<0>(mpval))[0] / (std::get<0>(mpval))[1];

        auto const mp_o = pdiffdata_->mp_o_;
        
        L_.reserve(pdata_->grid_num_ + 1);
        L_.assign(lo.begin(), lo.end());

        rf_.reserve(pdata_->grid_num_ + 1);
        pf_.reserve(pdata_->grid_num_ + 1);

        M_.reserve(pdata_->grid_num_ + 1);
        M_.assign(pdiffdata_->mo_.begin(), pdiffdata_->mo_.end());

        for (auto i = 0; i <= mp_o; i++) {
            rf_.push_back(std::pow(pdiffdata_->r_mesh_[i], pdata_->l_) * lo[i]);
            pf_.push_back(pdiffdata_->r_mesh_[i] * rf_[i]);
        }        

        for (auto i = mp_im1; i >= 0; i--) {
            L_.push_back(ratio * li[i]);
            M_.push_back(ratio * pdiffdata_->mi_[i]);

        	rf_.push_back(std::pow(r_mesh_i[i], pdata_->l_) * L_.back());
        	pf_.push_back(r_mesh_i[i] * rf_.back());
        }

        normalize();

        rho_.reserve(pdata_->grid_num_ + 1);
        for (auto const v : pf_) {
            rho_.push_back(sqr(v));
        }
    }

    Normalize<SchNormalize>::mymap SchNormalize::getresult() const
    {
        Normalize<SchNormalize>::mymap wf;
        wf["1 Mesh (r)"] = pdiffdata_->r_mesh_;
        wf["2 Eigen function"] = std::move(rf_);
        wf["3 Rho (mutiplied 4 * pi * r ** 2)"] = std::move(rho_);
        wf["4 Eigen function (mutiplied r)"] = std::move(pf_);

        return wf;
    }

    void SchNormalize::normalize()
    {
        Simpson simpson(pdiffdata_->dx_);
        auto const n = 1.0 / std::sqrt(simpson(pf_, pdiffdata_->r_mesh_));
        for (auto i = 0; i <= pdata_->grid_num_; i++) {
            rf_[i] *= n;
            pf_[i] *= n;
        }
    }

    // #endregion publicメンバ関数
}
