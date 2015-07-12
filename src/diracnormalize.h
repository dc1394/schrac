/*! \file diracnormalize.h
    \brief Dirac方程式を解いて得られた波動関数を正規化するクラスの宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _DIRACNORMALIZE_H_
#define _DIRACNORMALIZE_H_

#pragma once

#include "normalize.h"

namespace schrac {
    //! A class.
    /*!
        Dirac方程式を解いて得られた波動関数を正規化するクラス
    */
    class DiracNormalize final : public Normalize<DiracNormalize> {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdiffsolver 微分方程式のデータオブジェクト
        */
        DiracNormalize(std::shared_ptr<DiffSolver> const & pdiffsolver) :
            Normalize<DiracNormalize>(pdiffsolver)
        {
        }

        //! A destructor.
        /*!
        デフォルトデストラクタ
        */
        ~DiracNormalize() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function.
        /*!
            波動関数と4πr ** 2のかかった形の電子密度を求める
        */
        void evaluate();

    public:
        //! A public member function.
        /*!
            求めた結果を返す
            \return メッシュと波動関数が格納されたmap
        */
        Normalize<DiracNormalize>::mymap getresult() const;

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function.
        /*!
            波動関数を正規化する
        */
        void normalize();

        // #endregion privateメンバ関数

        // #region メンバ変数

    private:
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

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        DiracNormalize() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        DiracNormalize(DiracNormalize const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        DiracNormalize & operator=(DiracNormalize const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif // _DIRACNORMALIZE_H_
