/*! \file data.h
    \brief インプットファイルの各種データの構造体の宣言
    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD-2 License.
*/

#ifndef _DATA_H_
#define _DATA_H_

#pragma once

#include "ci_string.h"
#include <array>                // for std::array
#include <cstdint>              // for std::int32_t, std::uint8_t
#include <boost/optional.hpp>   // for boost::optional

namespace schrac {
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
            SCFの収束判定条件の値のデフォルト値
        */
        static auto constexpr SCF_CRITERION_DEFAULT = 5.0;

        //!  A public static member variable (constant expression).
        /*!
            SCFの最大ループ回数のデフォルト値
        */
        static auto constexpr SCF_MAXITER_DEFAULT = 40;

        //!  A public static member variable (constant expression).
        /*!
            電子密度を合成するときの重みのデフォルト値
        */
        static auto constexpr SCF_MIXING_WEIGHT_DEFAULT = 0.3;

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
		

    //! A struct.
    /*!
        インプットファイルの各種データの構造体
    */
	struct Data final {
        // #region 列挙型
        
        //!  A enumerated type
        /*!
            解く方程式のタイプを表す列挙型
        */
        enum class Eq_type {
            // Schrödinger方程式
            SCH,
            // スカラ相対論補正
            SDIRAC,
            // Dirac方程式
            DIRAC
        };

        //!  A enumerated type
        /*!
            微分方程式の解法の種類を表す列挙型
        */
        enum class Solver_type {
            // Adams Bashforth Moulton法
            ADAMS_BASHFORTH_MOULTON,
            // Bulirsch-Stoer法
            BULIRSCH_STOER,
            // コントロールされたRunge-Kutta法
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
        
        //!  A public static member variable (constant).
        /*!
            「ALPHA」の文字列
        */
        static ci_string const ALPHA;

        //!  A public static member variable (constant).
        /*!
            「BETA」の文字列
        */
        static ci_string const BETA;

        //!  A public static member variable (constant).
        /*!
            元素記号の配列
        */
        static std::array<std::string, 2> const Chemical_Symbol;
        
        //!  A public member variable.
        /*!
            原子名
        */
        std::string chemical_symbol_;
        
        //!  A public member variable.
        /*!
            微分方程式を解くときの許容誤差
        */
        double eps_ = EPS_DEFAULT;

        //!  A public member variable.
        /*!
            解く方程式のタイプ
        */
        Data::Eq_type eq_type_ = Data::Eq_type::SCH;
        
        //!  A public member variable.
        /*!
            微分方程式を解くときのメッシュの数
        */
        std::int32_t grid_num_ = GRID_NUM_DEFAULT;

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
        double mat_po_ratio_ = MAT_PO_RATIO_DEFAULT;

        //!  A public member variable.
        /*!
            主量子数
        */
        std::uint8_t n_;
        
        //!  A public member variable.
        /*!
            search_lowerE_から0までをいくつに分割して検索するか    
        */
        std::int32_t num_of_partition_ = NUM_OF_PARTITION_DEFAULT;
        
        //!  A public member variable.
        /*!
            計算対象の軌道
        */
        std::string orbital_;
        
        //!  A public member variable.
        /*!
            密度の初期値ρ0(r)のための係数c（ρ0(r) = c * exp(- alpha * r)
        */
        boost::optional<double> rho0_c_;
        
        //!  A public member variable.
        /*!
            密度の初期値ρ0(r)のための係数alpha（ρ0(r) = c * exp(- alpha * r)
        */
        boost::optional<double> rho0_alpha_;

        //!  A public member variable.
        /*!
            SCFの収束判定条件の値
        */
        double scf_criterion_ = SCF_CRITERION_DEFAULT;

        //!  A public member variable.
        /*!
            SCFの最大ループ回数
        */
        std::int32_t scf_maxiter_ = SCF_MAXITER_DEFAULT;

        //!  A public member variable.
        /*!
            電子密度を合成するときの重み
        */
        double scf_mixing_weight_ = SCF_MIXING_WEIGHT_DEFAULT;

        //!  A public member variable.
        /*!
            固有値探索を始める値
        */
        boost::optional<double> search_lowerE_;

        //!  A public member variable.
        /*!
            使用する微分方程式のソルバー
        */
        Data::Solver_type solver_type_ = Data::Solver_type::BULIRSCH_STOER;

        //!  A public member variable.
        /*!
            計算対象のスピン軌道
        */
        ci_string spin_orbital_;

        //!  A public member variable.
        /*!
            TBBを使用するかどうか
        */
        bool usetbb_ = false;

        //!  A public member variable.
        /*!
            微分方程式を解くときのメッシュの最大値
        */
        double xmax_ = XMAX_DEFAULT;

        //!  A public member variable.
        /*!
            微分方程式を解くときのメッシュの最小値
        */
        double xmin_ = XMIN_DEFAULT;

        //!  A public member variable.
        /*!
            原子核の電荷
        */
        double Z_;

        // #endregion メンバ変数
    };
}

#endif	// _DATA_H_

