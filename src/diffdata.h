/*! \file diffdata.h
    \brief 微分方程式のデータを集めた構造体の宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
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
            デフォルトデストラクタ
        */
        ~DiffData() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ変数

        //!  A public member variable (constant).
        /*!
            ノードの数
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
        double dx_;

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
            rのメッシュ
        */
        dvector r_mesh_;

        //!  A public member variable.
        /*!
            無限遠に近い点からのrのメッシュ
        */
        dvector r_mesh_i_;

        //!  A public member variable.
        /*!
            今回微分方程式を解くことによって得たノードの数
        */
        std::int32_t thisnode_;
        
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

        // #endregion メンバ変数
	};
}

#endif	// _DIFFDATA_H_

