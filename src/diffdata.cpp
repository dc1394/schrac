/*! \file diffdata.cpp
    \brief 微分方程式のデータを集めた構造体の実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "diffdata.h"
#include <boost/cast.hpp>   // for boost::numeric_cast

namespace schrac {
    // #region コンストラクタ

    DiffData::DiffData(std::shared_ptr<Data> const & pdata) :
        node_(pdata->n_ - pdata->l_ - 1),
        pdata_(pdata),
        thisnode_(0),
        Z_(pdata->Z_)
    {
        auto const grid_num = pdata_->grid_num_;

        mp_o_ = boost::numeric_cast<std::int32_t>(std::round(static_cast<double>(grid_num) * pdata_->mat_po_ratio_));
        mp_i_ = grid_num - mp_o_;
        
        auto const osize = boost::numeric_cast<dvector::size_type>(mp_o_ + 1);
        auto const isize = boost::numeric_cast<dvector::size_type>(mp_i_ + 1);

        dx_ = (pdata_->xmax_ - pdata_->xmin_) / static_cast<double>(grid_num - 1);

        // メモリ確保
        x_o_.resize(osize);
        x_i_.resize(isize);
        r_mesh_i_.resize(isize);
        lo_.reserve(osize);
        li_.reserve(isize);
        mo_.reserve(osize);
        mi_.reserve(isize);

        auto const len = grid_num - boost::numeric_cast<std::int32_t>(isize);

        for (auto i = 0; i <= mp_o_; i++) {
            auto const x = pdata_->xmin_ + static_cast<double>(i) * dx_;
            x_o_[i] = x;
        }

        for (auto i = grid_num; i > len; i--) {
            auto const x = pdata_->xmin_ + static_cast<double>(i) * dx_;
            x_i_[grid_num - i] = x;
            r_mesh_i_[grid_num - i] = std::exp(x);
        }
    }
}

