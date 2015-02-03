/*! \file normalize.h
    \brief 得られた波動関数を正規化するクラスの宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
*/

#ifndef _NORMALIZE_H_
#define _NORMALIZE_H_

#pragma once

#include "eigenvaluesearch.h"
#include <tuple>                // for std::tuple
#include <unordered_map>

namespace schrac {
	class Normalize final {
        // #region 型エイリアス

        using large_small_wf_hash = std::unordered_map<std::string, dvector>;

        using r_rf_pf_tuple = std::tuple<dvector, dvector, dvector>;

        using r_rfls_pfls_tuple = std::tuple<dvector, large_small_wf_hash, large_small_wf_hash>;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdiff 微分方程式のデータオブジェクト
        */
        Normalize(std::shared_ptr<Diff> const & pdiff);

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~Normalize()
        {
        }

        // #endregion コンストラクタ・デストラクタ

        // #region プロパティ

        //template <typename T>

        // #region メンバ関数

        //! A public member function (const).
        /*!
        L[0] = LO(rMP), L[1] = LI(rMP), M[0] = M0(rMP), M[1] = MI(rMP)を代入し、
        LとMのstd::pairを返す
        \return LとMのstd::pair
        */

	private:
		const shared_ptr<Data> pdata_;
		const shared_ptr<DiffData> pdiffdata_;

		ldvector RV;
		ldvector XV;
		ldvector RF;
		ldvector PF;

		long double phi(std::size_t n) const;
		void WF_coalesce(const shared_ptr<Diff> & pdiff);
		void WF_coalesce_omp(const shared_ptr<Diff> & pdiff);
		long double simpson() const;
		long double simpson_omp() const;

	public:
				void operator()();
		const WF_Normalize::d3tup getptup() const
		{ return make_tuple(RV, RF, PF); }
	};
}

#endif // _NORMALIZE_H_
