#ifndef _DIFFDATA_H_
#define _DIFFDATA_H_

#pragma once

#include "data.h"
#include <memory>
#include <vector>

namespace schrac {
    using dvectpr = std::vector < double > ;

	struct DiffData {
        // #region コンストラクタ・デストラクタ

        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param E エネルギー固有値 
            \param pdata データオブジェクト
            \param TINY 絶対値が小さい方の閾値
        */
        DiffData(double E, std::shared_ptr<Data> const & pdata, double TINY);

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~DiffData()
        {
        }

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ変数

    private:
        //!  A private static member variable (constant expression).
        /*!
            行列Aのサイズ
        */
		static std::size_t constexpr AVECSIZE = 3;
		
        //!  A private static member variable (constant expression).
        /*!
            行列Bのサイズ
        */
        static std::size_t constexpr BVECSIZE = 5;
        
        //!  A private member variable (constant).
        /*!
            節の数
        */
        std::int32_t const node_;

        //!  A private static member variable (constant expression).
        /*!
            行列Bのサイズ
        */
        double E_;

		const std::shared_ptr<const Data> pdata_;

		
		int thisnode;

		std::size_t MP_O;
		std::size_t MP_I;
		int OSIZE;
		int ISIZE;

		const double TINY_;

		double Z;
		double DX;
		


		array<double, AVECSIZE> V_A;
		array<double, BVECSIZE> V_B;

		ldvector XV_I;
		ldvector XV_O;
		ldvector RV_I;
		ldvector RV_O;
		ldvector VP_I;
		ldvector VP_O;
		ldvector LO;
		ldvector MO;
		ldvector LI;
		ldvector MI;

		
	};

    double fnc_V(double x);
    void node_count(int i, const ldvector & WF);

	inline double DiffData::fnc_V(double x) const
	{
		return - Z * std::exp(- x);
	}

	inline void DiffData::node_count(int i, const ldvector & WF)
	{
		if (WF[i] * WF[i - 1] < 0.0)
			thisnode++;
	}
}

#endif	// _DIFFDATA_H_
