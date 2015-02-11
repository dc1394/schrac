/*! \file diracnormalize.h
    \brief Dirac方程式を解いて得られた波動関数を正規化するクラスの宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
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
        何もしないデストラクタ
        */
        ~DiracNormalize()
        {
        }

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
        Normalize<DiracNormalize>::myhash getresult() const;

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

#endif // _SCHNORMALIZE_H_
