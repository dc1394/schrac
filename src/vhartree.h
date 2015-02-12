/*! \file vhartree.h
    \brief Hartreeポテンシャルを求めるクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _VHARTREE_H_
#define _VHARTREE_H_

#pragma once

#include "data.h"
#include <functional>
#include <memory>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace schrac {
    //! A function.
    /*!
        gsl_interp_accelへのポインタを解放するクラス
    */
    auto const gsl_interp_accel_deleter = [](gsl_interp_accel * acc) {
        gsl_interp_accel_free(acc);
    };

    auto const gsl_spline_deleter = [](gsl_spline * spline) {
        gsl_spline_free(spline);
    };

    //! A class.
    /*!
        Hartreeポテンシャルを求めるクラス
    */
    class Vhartree final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param rho 電子密度
            \param r_mesh rのメッシュ
        */
        Vhartree(std::shared_ptr<Data> const & pdata, std::vector<double> const & rho, std::vector<double> const & r_mesh);

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~Vhartree()
        {
        }

        // #region メンバ関数

    private:
        bool addrho();

        //!  A private member function (const).
        /*!
            Hartreeポテンシャルの関数オブジェクト
        */
        double operator()(double x) const;

        template <typename Stepper>
        //!  A private member function.
        /*!
            Hartreeポテンシャルのメッシュを生成する
        */
        void make_vhartree_mesh();

        // #endregion メンバ関数

        // #region メンバ変数

        //! A private member variable.
        /*!
            gsl_interp_accelへのスマートポインタ
        */
        std::unique_ptr<gsl_interp_accel, decltype(gsl_interp_accel_deleter)> const acc;

        //! A private member variable.
        /*!
            gsl_interp_typeへのスマートポインタ
        */
        std::unique_ptr<gsl_spline, decltype(gsl_spline_deleter)> const spline;

        //!  A private member variable.
        /*!
            Hartreeポテンシャルの関数オブジェクト
        */
        std::function<double (double)> vhartree_;
        
        // #endregion メンバ変数
    };

}

#endif  // _VHARTREE_H_
