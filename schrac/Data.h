/*! \file Data.h
    \brief インプットファイルの各種データの構造体の宣言
    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _DATA_H_
#define _DATA_H_

#include <array>
#include <cstdint>
#include <boost/optional.hpp>
#include "ci_string.h"

namespace schrac {
    //! A class.
    /*!
        インプットファイルの各種データの構造体
    */
	struct final Data {
        // #region 列挙型
        
        //!  A enumerated type
        /*!
            解く方程式のタイプを表す列挙型
        */
        enum class eq_type {
            // Schrödinger方程式
            SCH,
            // スカラ相対論補正
            SDIRAC,
            // Dirac方程式
            DIRAC
        };

        //!  A enumerated type
        /*!
            微分方程式のソルバーの種類を表す列挙型
        */
		enum class solver_type {
            // Adams Bashforth Moulton法
            ADAMS_BASHFORTH_MOULTON, 
            // Bulirsch-Stoer法
            BULIRSCH_STOER,
            // 誤差がコントロールされたRunge-Kutta法
			CONTROLLED_RUNGE_KUTTA
        };

        // #endregion 列挙型

        // #region メンバ変数
        
        //!  A public static member variable (constant expression).
        /*!
            光速（原子単位系）
        */
		static auto constexpr c = 137.035999;

        //!  A public static membver variable (constant expression).
        /*!
            原子単位系での微細構造定数(= 1/c)
        */
		static auto constexpr al = 1.0 / c;

        //!  A public static member variable (constant expression).
        /*!
            微細構造定数の2乗の1/2
        */
        static auto constexpr al2half = 0.5 * al * al;
        
        //!  A public static member variable (constant expression).
        /*!
            微分方程式を解くときの許容誤差のデフォルト値
        */
        static auto constexpr EPS_DEFAULT = 1.0E-15;

        //!  A public static member variable (constant expression).
        /*!
            微分方程式を解くときのメッシュの数のデフォルト値
        */
        static auto constexpr GRID_NUM_DEFAULT = 20000;
        
        //!  A public static member variable (constant expression).
        /*!
            マッチングポイント（xmin〜xmaxまでの比率で表す）
        */
        static auto constexpr MAT_PO_RATIO_DEFAULT = 0.67;
        
        //!  A public static member variable (constant expression).
        /*!
            search_lowerE_から0までをいくつに分割して検索するかのデフォルトの値
        */
        static auto constexpr NUM_OF_PARTITION_DEFAULT = 300;
		
        //!  A public static member variable (constant expression).
        /*!
            微分方程式を解くときのメッシュの最小値のデフォルト値
        */
        static auto constexpr XMAX_DEFAULT = 5.0;

        //!  A public static member variable (constant expression).
        /*!
            微分方程式を解くときのメッシュの最小値のデフォルト値
        */
        static auto constexpr XMIN_DEFAULT = -7.0;
		
        //!  A public static member variable (constant).
        /*!
            「ALPHA」の文字列
        */
        static auto const ci_string ALPHA = "ALPHA";

        //!  A public static member variable (constant).
        /*!
            元素記号の配列
        */
        static auto const AtomSymbol = { "H", "He" };

        //!  A public static member variable (constant).
        /*!
            「BETA」の文字列
        */
        static auto const ci_string BETA = "BETA";
        
        //!  A public member variable.
        /*!
            原子名
        */
        std::string atom_;
        
        //!  A public member variable.
        /*!
            微分方程式を解くときの許容誤差
        */
        auto eps_ = EPS_DEFAULT;

        //!  A public member variable.
        /*!
            解く方程式のタイプ
        */
        auto eq_type_ = eq_type::SDIRAC;
        
        //!  A public member variable.
        /*!
            微分方程式を解くときのメッシュの数
        */
        auto grid_num_ = GRID_NUM_DEFAULT;

        //!  A public member variable.
        /*!
            全角運動量
        */
        double j_;

		//!  A public member variable.
        /*!
            量子数κ
        */
        double kappa_;

        //!  A public member variable.
        /*!
            方位量子数
        */
        std::uint8_t l_;
        
        //!  A public member variable.
        /*!
            マッチングポイント（xmin〜xmaxまでの比率で表す）
        */
        auto mat_po_ratio_ = Data::MAT_PO_RATIO_DEFAULT;

        //!  A public member variable.
        /*!
            主量子数
        */
        std::uint8_t n_;
        
        //!  A public member variable.
        /*!
            search_lowerE_から0までをいくつに分割して検索するか    
        */
        auto num_of_partiton_ = Data::NUM_OF_PARTITION_DEFAULT;
        
        //!  A public member variable.
        /*!
            計算対象の軌道
        */
        std::string orbital_;
        
        //!  A public member variable.
        /*!
            固有値検索を始める値
        */
        boost::optional<double> search_lowerE_;

        //!  A public member variable.
        /*!
            使用する微分方程式のソルバー
        */
		auto solve_type = Data::solve_type::BULIRSCH_STOER;

        //!  A public member variable.
        /*!
            計算対象のスピン軌道
        */
        ci_string spin_orbital_;

        //!  A public member variable.
        /*!
            並列計算をするときのスレッド数
        */
		auto threadnum_ = boost::none;

        //!  A public member variable.
        /*!
            微分方程式を解くときのメッシュの最大値
        */
        auto xmax_ = Data::XMAX_DEFAULT;

        //!  A public member variable.
        /*!
            微分方程式を解くときのメッシュの最小値
        */
        auto xmin_ = Data::XMIN_DEFAULT;

        //!  A public member variable.
        /*!
            原子番号
        */
        std::uint8_t Z_;

        // #endregion メンバ変数
    };
}

#endif	// _DATA_H_

