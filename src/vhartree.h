/*! \file vhartree.h
    \brief Hartreeポテンシャルを求めるクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _VHARTREE_H_
#define _VHARTREE_H_

#pragma once

#include "deleter.h"
#include "diffdata.h"
#include "property.h"
#include <memory>       // for std::unique_ptr

namespace schrac {
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
            \param r_mesh_ rのメッシュ
        */
        Vhartree(std::vector<double> const & r_mesh);

        //! A private copy constructor.
        /*!
            コピーコンストラクタ
        */
        Vhartree(Vhartree const & rhs);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Vhartree() = default;

        // #region メンバ関数
        
        //!  A public member function (const).
        /*!
            Hartreeポテンシャルの微分値を返す
            \param r 極座標のr
            \return Hartreeポテンシャルの微分値
        */
        double dvhartree_dr(double r) const;

        //!  A public member function.
        /*!
            Hartreeポテンシャルが境界条件を満たすようにセットする
            \param Z 原子核の電荷
        */
        void set_vhartree_boundary_condition(double Z);

        //!  A public member function.
        /*!
            Hartreeポテンシャルを初期化する
        */
        void vhart_init();

        //!  A public member function (const).
        /*!
            Hartreeポテンシャルの値を返す
            \param r 極座標のr
            \return Hartreeポテンシャルの値
        */
        double vhartree(double r) const;

        // #endregion メンバ関数

        // #region プロパティ

    public:
        //! A property.
        /*!
            Hartreeポテンシャルが格納された可変長配列へのプロパティ
        */
        Property<std::vector<double>> Vhart;

        // #endregion プロパティ

        // #region メンバ変数

    private:
        //! A private member variable.
        /*!
            gsl_interp_accelへのスマートポインタ
        */
        std::unique_ptr<gsl_interp_accel, decltype(gsl_interp_accel_deleter)> const acc_;
        
        //! A private member variable.
        /*!
            rのメッシュが格納された可変長配列
        */
        std::vector<double> const r_mesh_;

        //! A private member variable.
        /*!
            gsl_interp_typeへのスマートポインタ
        */
        std::unique_ptr<gsl_spline, decltype(gsl_spline_deleter)> const spline_;

        //! A private member variable.
        /*!
            Hartreeポテンシャルが格納された可変長配列
        */
        std::vector<double> vhart_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

    private:
        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Vhartree() = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Vhartree & operator=(Vhartree const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif  // _VHARTREE_H_
