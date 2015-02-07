/*! \file kinetic_energy.h
    \brief 運動エネルギーを計算するクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _KINETIC_ENERGY_
#define _KINETIC_ENERGY_

#include "diffdata.h"
#include <memory>

#pragma once

namespace schrac {
    //! A class.
    /*!
        運動エネルギーを計算するクラス
    */
	class Kinetic_Energy final {
        // #region 型エイリアス

        using dvector = std::vector<double>;

        // #endregion 型エイリアス

		// #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdata データオブジェクト
        */
        Kinetic_Energy(std::shared_ptr<DiffData> pdiffdata) :
            pdiffdata_(pdiffdata)
        {
        }

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~Kinetic_Energy()
        {
        }

        // #region メンバ関数

        //! A public member function (const).
        /*!
            運動エネルギーを計算する
            \param L 関数L(x)のstd::vector
            \param M 関数L(x)のstd::vector
            \param r rのメッシュが格納されたstd::vector
            \return 運動エネルギー
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

        // #endregion メンバ変数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Kinetic_Energy() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        Kinetic_Energy(Kinetic_Energy const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Kinetic_Energy & operator=(Kinetic_Energy const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif  // _KINETIC_ENERGY_H_
