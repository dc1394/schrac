/*! \file schnormalize.h
    \brief Sch方程式を解いて得られた波動関数を正規化するクラスの宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
*/

#ifndef _SCHNORMALIZE_H_
#define _SCHNORMALIZE_H_

#pragma once

#include "normalize.h"
#include <tuple>                                // for std::tuple
#include <boost/optional.hpp>                   // for boost::optional

namespace schrac {
    //! A class.
    /*!
        Sch方程式を解いて得られた波動関数を正規化するクラス
    */
    class SchNormalize final : public Normalize<SchNormalize> {
        // #region 型エイリアス

    public:
        using r_rf_pf_tuple = std::tuple<dvector, dvector, dvector>;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdiffsolver 微分方程式のデータオブジェクト
        */
        SchNormalize(std::shared_ptr<DiffSolver> const & pdiffsolver) :
            Normalize<SchNormalize>(pdiffsolver)
        {
        }

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~SchNormalize() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function.
        /*!
            波動関数と4πr ** 2のかかった形の電子密度を求める
            \param prho 4πr ** 2のかかった形の電子密度が格納された可変長配列へのboost::optional
        */
        void evaluate(boost::optional<std::vector<double>> const & prho);
        
        //! A public member function.
        /*!
            求めた結果を返す
        */
        Normalize<SchNormalize>::mymap getresult() const;

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
        
        //! A private member variable.
        /*!
            角度方向のrをかけた固有関数
        */
        dvector pf_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        SchNormalize() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        SchNormalize(SchNormalize const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        SchNormalize & operator=(SchNormalize const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif // _SCHNORMALIZE_H_
