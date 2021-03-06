﻿/*! \file rho.h
    \brief 電子密度ρ(r)を求めるクラスの宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "data.h"
#include "property.h"
#include <memory>           // for std::unique_ptr 
#include <vector>           // for std::vector
#include <gsl/gsl_spline.h> // for gsl_interp_accel, gsl_interp_accel_free, gsl_spline, gsl_spline_free

namespace schrac {
    //! A class.
    /*!
        電子密度ρ(r)を求めるクラス
    */
    class Rho final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdiffdata 微分方程式のデータオブジェクト
        */
        Rho(std::shared_ptr<DiffData> const & pdiffdata);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Rho() = default;

        // #region メンバ関数

        //!  A public member function.
        /*!
            現在の電子密度ρ(r)で初期化する
        */
        void init();

        //!  A public member function (const).
        /*!
            電子密度ρ(r)を返す
            \param r rの値
            \return ρ(r)の値
        */
        double operator()(double r) const;

        //!  A public member function (const).
        /*!
            新しい電子密度ρnew(r)と、電子密度ρ(r)を一次混合する
            \param rhonew 新しい電子密度ρnew(r)
        */
        void rhomix(dvector const & rhonew);

        // #endregion メンバ関数

        // #region プロパティ

    public:
        //! A property.
        /*!
            電子密度が格納された可変長配列へのプロパティ
        */
        Property< std::vector<double> > const PRho;

        // #endregion プロパティ

        // #region メンバ変数

    private:
        //! A private member variable.
        /*!
            gsl_interp_accelへのスマートポインタ
        */
        std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)> const acc_;

        //!  A private member variable.
        /*!
            データオブジェクト
        */
        std::shared_ptr<DiffData> pdiffdata_;

        //!  A private member variable.
        /*!
            密度ρ(r)
        */
        std::vector<double> rho_;

        //! A private member variable.
        /*!
            gsl_splineへのスマートポインタ
        */
        std::unique_ptr<gsl_spline, decltype(&gsl_spline_free)> const spline_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Rho() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        Rho(Rho const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Rho & operator=(Rho const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}
