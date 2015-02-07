/*! \file potential_energy.h
    \brief ポテンシャルエネルギーを計算するクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _POTENTIAL_ENERGY_H_
#define _POTENTIAL_ENERGY_H_

#pragma once

#include "diffdata.h"
#include "simpson.h"
#include <memory>

namespace schrac {
    //! A class.
    /*!
        ポテンシャルエネルギーを計算するクラス
    */
    class Potential_Energy final {
        // #region 型エイリアス

        using dvector = std::vector < double > ;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdata データオブジェクト
            */
        Potential_Energy(std::shared_ptr<DiffData> pdiffdata, dvector const & pf) :
            pdiffdata_(pdiffdata), pf_(pf)
        {
        }

        //! A destructor.
        /*!
            何もしないデストラクタ
            */
        ~Potential_Energy()
        {
        }

        // #region メンバ関数

        //! A public member function (const).
        /*!
            ポテンシャルエネルギーを計算する
            \param L 関数L(x)のstd::vector
            \param M 関数L(x)のstd::vector
            \param r rのメッシュが格納されたstd::vector
            \return ポテンシャルエネルギー
        */
        double operator()(dvector const & L, dvector const & M, dvector const & r) const;

        // #endregion メンバ関数

        // #region メンバ変数

    private:
        //! A private member variable (const).
        /*!
            データオブジェクト
        */
        std::shared_ptr<DiffData> const pdiffdata_;

        //! A private member variable (const).
        /*!
            波動関数PF
        */
        dvector const pf_;

        // #endregion メンバ変数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Potential_Energy() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        Potential_Energy(Potential_Energy const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Potential_Energy & operator=(Potential_Energy const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };

    double Potential_Energy::operator()(Potential_Energy::dvector const & L, Potential_Energy::dvector const & M, Potential_Energy::dvector const & r) const
    {
        return Simpson(pdiffdata_->dx_)(pf_, pf_, r, 2);
    }
}

#endif  // _POTENTIAL_ENERGY_H_
