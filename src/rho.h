/*! \file rho.h
    \brief 電子密度ρ(r)を求めるクラスの宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
*/

#include "data.h"
#include "deleter.h"
#include "property.h"
#include <memory>
#include <vector>

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
            何もしないデストラクタ
        */
        ~Rho()
        {
        }

        // #region メンバ関数

        //!  A public member function (const).
        /*!
            電子密度ρ(r)を返す
            \param r rの値
            \return ρ(r)の値
        */
        double operator()(double r) const;

        //!  A public member function (const).
        /*!
            電子密度ρ(r)を返す
            \param r rの値
            \return ρ(r)の値
        */
        void rhomix();

        // #endregion メンバ関数

        // #region プロパティ

    public:
        //! A property.
        /*!
            電子密度が格納された可変長配列へのプロパティ
        */
        Property<std::vector<double>> const PRho;

        // #endregion プロパティ

        // #region メンバ変数

    private:
        //! A private member variable.
        /*!
            gsl_interp_accelへのスマートポインタ
        */
        std::unique_ptr<gsl_interp_accel, decltype(gsl_interp_accel_deleter)> const acc_;

        //!  A private member variable.
        /*!
            密度ρ(r)
        */
        std::shared_ptr<Data> pdata_;

        //!  A private member variable.
        /*!
            密度ρ(r)
        */
        std::vector<double> rho_;

        //! A private member variable.
        /*!
            gsl_interp_typeへのスマートポインタ
        */
        std::unique_ptr<gsl_spline, decltype(gsl_spline_deleter)> const spline_;

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