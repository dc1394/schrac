/*! \file diffdata.cpp
    \brief 微分方程式のデータを集めた構造体の実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "DiffData.h"
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

		mp_o_ = boost::numeric_cast<std::int32_t>(std::round(static_cast<double>(grid_num - 1) * pdata_->mat_po_ratio_));
		mp_i_ = grid_num - mp_o_ - 1;
		
        auto const osize = boost::numeric_cast<dvector::size_type>(mp_o_ + 1);
        auto const isize = boost::numeric_cast<dvector::size_type>(mp_i_ + 1);

		DX_ = (pdata_->xmax_ - pdata_->xmin_) / static_cast<double>(grid_num - 1);

		// メモリ確保
		x_o_.resize(osize);
		x_i_.resize(isize);
		r_o_.resize(osize);
		r_i_.resize(isize);
		vr_o_.resize(osize);
		//VP_I.resize(isize);
		lo_.reserve(osize);
		li_.reserve(isize);
		mo_.reserve(osize);
		mi_.reserve(isize);

		auto const len = boost::numeric_cast<std::int32_t>(grid_num - isize);
		//if (pdata_->ompthread_) {
			//#pragma omp parallel
			{
				//#pragma omp for nowait
				for (auto i = 0; i < osize; i++) {
					auto const x = pdata_->xmin_ + static_cast<double>(i) * DX_;
					x_o_[i] = x;
					r_o_[i] = std::exp(x);
					vr_o_[i] = fnc_V(x);
				}
				//#pragma omp for nowait
				for (auto i = boost::numeric_cast<std::int32_t>(grid_num - 1); i >= len; i--) {
					auto const x = pdata_->xmin_ + static_cast<double>(i) * DX_;
					x_i_[grid_num - 1 - i] = x;
					r_i_[grid_num - 1 - i] = std::exp(x);
					//VP_I[grid_num - 1 - i] = fnc_V(x);
				}
			}
		/*} else {
			for (int i = 0; i < osize; i++) {
				const double x = pdata_->xmin + static_cast<const double>(i) * DX;
				XV_O[i] = x;
				RV_O[i] = std::exp(x);
				VP_O[i] = fnc_V(x);
			}
			for (int i = boost::numeric_cast<const int>(grid_num - 1); i >= len; i--) {
				const double x = pdata_->xmin + static_cast<const double>(i) * DX;
				XV_I[grid_num - 1 - i] = x;
				RV_I[grid_num - 1 - i] = std::exp(x);
				VP_I[grid_num - 1 - i] = fnc_V(x);
			}
		}*/
	}

    // #region メンバ関数

    double DiffData::fnc_V(double x) const
    {
        return -Z_ * std::exp(-x);
    }

    void DiffData::node_count(std::int32_t i, dvector const & WF)
    {
        if (WF[i] * WF[i - 1] < 0.0) {
            thisnode_++;
        }
    }

    // #endregion メンバ関数
}
