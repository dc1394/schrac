﻿/*! \file simpson.h
    \brief std::vectorに格納された関数を、Simpsonの法則で積分するクラスの宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _SIMPSON_H_
#define _SIMPSON_H_

#pragma once

#include <cstdint>  // for std::int32_t
#include <vector>   // for std::vector

namespace schrac {
    //! A class.
    /*!
        std::vectorに格納された関数を、Simpsonの法則で積分するクラス
    */
    class Simpson final {
        // #region 型エイリアス

        using dvector = std::vector<double>;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param dx メッシュの間隔
        */
        Simpson(double dx) : dx_(dx)
        {
        }

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Simpson() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //! A public member function (const).
        /*!
            f(r) * f(r) * rをシンプソンの公式で積分する
            \param f 関数f(r)のstd::vector
            \param r rのメッシュが格納されたstd::vector
            \return 積分した値
        */
        double operator()(dvector const & f, dvector const & r) const;

        //! A public member function (const).
        /*!
            f(r) * g(r) * r ** nをシンプソンの公式で積分する
            \param f 関数f(r)のstd::vector
            \param g 関数g(r)のstd::vector
            \param r rのメッシュが格納されたstd::vector
            \param n r ** nのnの値
            \return 積分した値
        */
        double operator()(dvector const & f, dvector const & g, dvector const & r, std::int32_t n) const;

        // #endregion メンバ関数
        
        // #region メンバ変数

    private:
        //! A private member variable (const).
        /*!
            メッシュの間隔
        */
        double const dx_;

        // #endregion メンバ変数
        
        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Simpson() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        Simpson(Simpson const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Simpson & operator=(Simpson const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif  // _SIMPSON_H_
