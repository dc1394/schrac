/*! \file diffdata.cpp
    \brief 微分方程式のデータを集めた構造体の実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "DiffData.h"
#include <boost/cast.hpp>   // for boost::numeric_cast

namespace schrac {
	// #region コンストラクタ

    DiffData::DiffData(double E, std::shared_ptr<Data> const & pdata, double TINY) :
        TINY_(TINY),
        E_(E),
        node_(pdata->n_ - pdata->l_ - 1),
        pdata_(pdata),
		thisnode_(0),
        Z_(pdata->Z_)
	{
		auto const grid_num = pdata_->grid_num_;

		MP_O_ = boost::numeric_cast<std::int32_t>(std::round(static_cast<double>(grid_num - 1) * pdata_->mat_po_ratio_));
		MP_I_ = grid_num - MP_O_ - 1;
		
        auto const osize = boost::numeric_cast<dvector::size_type>(MP_O_ + 1);
        auto const isize = boost::numeric_cast<dvector::size_type>(MP_I_ + 1);

		DX_ = (pdata_->xmax_ - pdata_->xmin_) / static_cast<double>(grid_num - 1);

		// メモリ確保
		XV_O.resize(osize);
		XV_I.resize(isize);
		RV_O_.resize(osize);
		RV_I_.resize(isize);
		VP_O_.resize(osize);
		//VP_I.resize(isize);
		LO_.resize(osize);
		LI_.resize(isize);
		MO_.resize(osize);
		MI_.resize(isize);

		auto const len = boost::numeric_cast<std::int32_t>(grid_num - isize);
		//if (pdata_->ompthread_) {
			//#pragma omp parallel
			{
				//#pragma omp for nowait
				for (auto i = 0; i < osize; i++) {
					const double x = pdata_->xmin_ + static_cast<double>(i) * DX_;
					XV_O[i] = x;
					RV_O_[i] = std::exp(x);
					VP_O_[i] = fnc_V(x);
				}
				//#pragma omp for nowait
				for (auto i = boost::numeric_cast<std::int32_t>(grid_num - 1); i >= len; i--) {
					auto const x = pdata_->xmin_ + static_cast<double>(i) * DX_;
					XV_I[grid_num - 1 - i] = x;
					RV_I_[grid_num - 1 - i] = std::exp(x);
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
