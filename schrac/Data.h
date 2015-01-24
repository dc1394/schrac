/*! \file Data.h
    \brief インプットファイルの各種データの構造体の宣言
    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _DATA_H_
#define _DATA_H_

#include <array>
#include "ci_string.h"

namespace schrac {
    //! A class.
    /*!
        インプットファイルの各種データの構造体
    */
	struct final Data {
        // #region メンバ変数

        //!  A public variable (constant).
        /*!
            「ALPHA」の文字列
        */
        static auto const ci_string ALPHA = "ALPHA";

        //!  A public variable (constant).
        /*!
            元素記号の配列
        */
        static auto const AtomSymbol = { "H", "He" };

        //!  A public variable (constant).
        /*!
            「BETA」の文字列
        */
        static auto const ci_string BETA = "BETA";

        //!  A public variable (constant expression).
        /*!
            光速（原子単位系）
        */
		static auto constexpr c = 137.035999;

        //!  A public variable (constant expression).
        /*!
            原子単位系での微細構造定数(= 1/c)
        */
		static auto constexpr al = 1.0 / c;

        //!  A public variable (constant expression).
        /*!
            微細構造定数の2乗の1/2
        */
        static auto constexpr al2half = 0.5 * al * al;

        //!  A public variable.
        /*!
            原子名
        */
        std::string atom;
        
		int const threadnum_;

		unsigned int Z;
		unsigned int n;
		unsigned int l;


		std::string orbital;

		ci_string spin_orbital;

		long double j_;									// 全角運動量
		long double kappa;								// κ

		enum eq_type { SCH, SDIRAC, DIRAC } eqtype;
		static const eq_type def_eq_type = SDIRAC;

		long double xmin;
		static const long double XMIN_DEFAULT;
		long double xmax;
		static const long double XMAX_DEFAULT;
		std::size_t grid_num;
		static const std::size_t GRID_NUM_DEFAULT = 20000;

		long double eps;
		static const long double EPS_DEFAULT;
		enum solve_type {
			MOD_EULER, RUNGE_KUTTA, RK_ADAPSTEP, BULIRSCH_STOER } stype;

		static const solve_type def_solve_type = BULIRSCH_STOER;

		boost::optional<long double> search_lowerE;

		static const std::size_t NUM_OF_PARTITION_DEFAULT = 300;
		std::size_t num_of_partiton;

		static const long double MAT_PO_RATIO_DEFAULT;
		long double mat_po_ratio;

		explicit Data(int ompthread)
		 :	ompthread_(ompthread) {}
		‾Data() {}
	};
}

#endif	// _DATA_H_
