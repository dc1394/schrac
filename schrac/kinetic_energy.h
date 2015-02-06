/*! \file kinetic_energy.h
    \brief 運動エネルギーを計算するクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _KINETIC_ENERGY_
#define _KINETIC_ENERGY_

#include "data.h"
#include "simpson.h"
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
        Kinetic_Energy(std::shared_ptr<Data> pdata) :
            pdata_(pdata)
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
            \param f 関数f(r)のstd::vector
            \param r rのメッシュが格納されたstd::vector
            \return 運動エネルギー
        */
        void operator()(dvector const & L, dvector const & M, dvector const & r) const;

        // #endregion メンバ関数

        // #region メンバ変数

    private:
        //! A private member variable (const).
        /*!
            データオブジェクト
        */
        std::shared_ptr<Data> const pdata_;

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
