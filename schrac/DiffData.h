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
            \param pdata データオブジェクト
        */
        DiffData(std::shared_ptr<Data> const & pdata);

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
        
        // #endregion メンバ関数

        // #region メンバ変数
        
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
            無限遠に近い点から解いた関数Lの数表
        */
        dvector li_;

        //!  A public member variable.
        /*!
            原点に近い点から解いた関数Lの数表
        */
        dvector lo_;

        //!  A public member variable.
        /*!
            無限遠に近い点から解いた関数Mの数表
        */
        dvector mi_;
        
        //!  A public member variable.
        /*!
            原点に近い点から解いた関数Mの数表
        */
        dvector mo_;

        //!  A public member variable.
        /*!
            無限遠に近い点から数えたマッチングポイント
        */
        std::int32_t mp_i_;
        
        //!  A public member variable.
        /*!
            原点から近い方から数えたマッチングポイント
        */
        std::int32_t mp_o_;

        //!  A public member variable.
        /*!
            無限遠に近い点からのrのメッシュ
        */
        dvector r_i_;

        //!  A public member variable.
        /*!
            原点に近い点からのrのメッシュ
        */
        dvector r_o_;

        //!  A public member variable.
        /*!
            今回微分方程式を解くことによって得た節の数
        */
        std::int32_t thisnode_;
        
        //!  A public member variable.
        /*!
            原点に近い点から微分方程式を解くときに使うポテンシャルV(r)
        */
        dvector vr_o_;

        //!  A public member variable.
        /*!
            無限遠に近い点からのxのメッシュ
        */
		dvector x_i_;

        //!  A public member variable.
        /*!
            原点に近い点からのxのメッシュ
        */
        dvector x_o_;

		//dvector VP_I;
		

        //int OSIZE;
        //int ISIZE;

        // #endregion メンバ変数
	};
}

#endif	// _DIFFDATA_H_
