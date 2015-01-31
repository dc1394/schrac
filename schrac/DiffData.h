﻿/*! \file diffdata.h
    \brief 微分方程式のデータを集めた構造体の宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _DIFFDATA_H_
#define _DIFFDATA_H_

#pragma once

#include "data.h"
#include <memory>   // for std::shared_ptr
#include <vector>   // for std::vector

namespace schrac {
    // #region 型エイリアス

    using dvector = std::vector < double > ;

    // #endregion 型エイリアス

    //! A struct.
    /*!
        微分方程式のデータを集めた構造体
    */
	struct DiffData final {
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

        // #region メンバ関数

        //!  A public member function (const).
        /*!
            xのメッシュにおけるV、すなわちV(x)を計算する
            \param x xの値
            \return V(x)の値
        */
        double fnc_V(double x) const;
        
        //!  A public member function.
        /*!
            波動関数の節の数を調べる
            \param i メッシュのインデックス
            \param WF 波動関数φ
        */
        void node_count(std::int32_t i, dvector const & WF);

        // #endregion メンバ関数

        // #region メンバ変数
            
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
        
        //!  A public member variable (constant).
        /*!
            節の数
        */
        std::int32_t const node_;

        //!  A public member variable (constant).
        /*!
            データオブジェクト
        */
        std::shared_ptr<Data> const pdata_;
        
        //!  A public member variable (constant).
        /*!
            絶対値が小さい方の閾値
        */
        double const TINY_;
        
        //!  A public member variable (constant).
        /*!
            原子核の電荷
        */
        double const Z_;

        //!  A public member variable (constant).
        /*!
            微分方程式を解くときのメッシュの間隔
        */
        double DX_;

        //!  A public member variable.
        /*!
            エネルギー固有値
        */
        double E_;

        //!  A public member variable.
        /*!
            無限遠に近い点から解いた関数Lの数表が格納された可変長配列
        */
        dvector LI_;

        //!  A public member variable.
        /*!
            原点に近い点から解いた関数Lの数表が格納された可変長配列
        */
        dvector LO_;

        //!  A public member variable.
        /*!
            無限遠に近い点から関数Mの数表が格納された可変長配列
        */
        dvector MI_;
        
        //!  A public member variable.
        /*!
            原点に近い点から関数Mの数表が格納された可変長配列
        */
        dvector MO_;

        //!  A public member variable.
        /*!
            無限遠に近い点から数えたマッチングポイント
        */
        std::int32_t MP_I_;
        
        //!  A public member variable.
        /*!
            原点から近い方から数えたマッチングポイント
        */
        std::int32_t MP_O_;

        //!  A public member variable.
        /*!
            無限遠に近い点からrのメッシュが格納された可変長配列
        */
        dvector RV_I_;

        //!  A public member variable.
        /*!
            原点に近い点からrのメッシュが格納された可変長配列
        */
        dvector RV_O_;

        //!  A public member variable.
        /*!
            今回微分方程式を解くことによって得た節の数
        */
        std::int32_t thisnode_;

        //!  A public member variable.
        /*!
            行列A
        */
        std::array<double, AVECSIZE> V_A_;

        //!  A public member variable.
        /*!
            行列B
        */
        std::array<double, BVECSIZE> V_B_;
        
        //!  A public member variable.
        /*!
            原点に近い点から微分方程式を解くときに使うポテンシャルVの可変長配列（Rのメッシュ）
        */
        dvector VP_O_;

        //!  A public member variable.
        /*!
            無限遠に近い点からxのメッシュが格納された可変長配列
        */
		dvector XV_I_;

        //!  A public member variable.
        /*!
            原点に近い点からxのメッシュが格納された可変長配列
        */
        dvector XV_O_;

		//dvector VP_I;
		

        //int OSIZE;
        //int ISIZE;

        // #endregion メンバ変数
	};
}

#endif	// _DIFFDATA_H_