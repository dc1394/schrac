/*! \file sdiracnormalize.h
    \brief Scalar Dirac方程式を解いて得られた波動関数を正規化するクラスの宣言

    Copyright c 2015 @dc1394 All Rights Reserved.
*/

#ifndef _SDIRACNORMALIZE_H_
#define _SDIRACNORMALIZE_H_

#pragma once

#include "normalize.h"

namespace schrac {
    //! A class.
    /*!
        Scalar Dirac方程式を解いて得られた波動関数を正規化するクラス
    */
    class SDiracNormalize final : public Normalize<SDiracNormalize> {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdiffsolver 微分方程式のデータオブジェクト
        */
        SDiracNormalize(std::shared_ptr<DiffSolver> const & pdiffsolver) :
            Normalize<SDiracNormalize>(pdiffsolver)
        {
        }

        //! A destructor.
        /*!
        デフォルトデストラクタ
        */
        ~SDiracNormalize() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //! A public member function.
        /*!
        波動関数を求める
        */
        void evaluate();

    private:
        //! A public member function.
        /*!
        波動関数を正規化する
        */
        void normalize();

    public:
        //! A public member function.
        /*!
        求めた結果を返す
        */
        Normalize<SDiracNormalize>::mymap getresult() const;

        // #endregion メンバ関数

        // #region メンバ変数

        //! A private member variable.
        /*!
        固有関数
        */
        dvector rf_;

        //! A private member variable.
        /*!
        角度方向のrをかけた固有関数のlarge成分
        */
        dvector pf_large_;

        //! A private member variable.
        /*!
        角度方向のrをかけた固有関数のsmall成分
        */
        dvector pf_small_;

        // #endregion メンバ変数

    private:
        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        SDiracNormalize() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        SDiracNormalize(SDiracNormalize const &) = delete;

        //! A private member function (deleted).
        /*!
        operator=()の宣言（禁止）
        \param コピー元のオブジェクト（未使用）
        \return コピー元のオブジェクト
        */
        SDiracNormalize & operator=(SDiracNormalize const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif // _SCHNORMALIZE_H_
