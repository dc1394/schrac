﻿/*! \file diffsolver.h
    \brief 微分方程式を解くクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _DIFFSOLVER_H_
#define _DIFFSOLVER_H_

#pragma once

#include "diffdata.h"
#include "property.h"
#include "rho.h"
#include "solvelinearequ.h"
#include "vhartree.h"

namespace schrac {
    // #region 型エイリアス

    using myarray = std::array < double, 2 >;

    // #endregion 型エイリアス

    //! A class.
    /*!
        微分方程式を解くクラス
    */
	class DiffSolver final {
        // #region 型エイリアス

    private:
        using mypair = std::pair < myarray, myarray > ;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:

        //! A constructor.
        /*!
            コンストラクタ
            \param pdata データオブジェクト
            \param pdiffdata 微分方程式のデータオブジェクト
        */
        DiffSolver(std::shared_ptr<Data> const & pdata, std::shared_ptr<DiffData> const & pdiffdata);

        //! A constructor.
        /*!
            コンストラクタ
            \param pdata データオブジェクト
            \param pdiffdata 微分方程式のデータオブジェクト
            \param rho 密度ρのオブジェクト
            \param pvh Hartreeポテンシャルオブジェクト
        */
        DiffSolver(std::shared_ptr<Data> const & pdata, std::shared_ptr<DiffData> const & pdiffdata, std::shared_ptr<Rho> const & prho, std::shared_ptr<Vhartree> const & pvh);

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~DiffSolver()
        {
        }

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数
        
        //! A public member function (const).
        /*!
            L[0] = LO(rMP), L[1] = LI(rMP), M[0] = M0(rMP), M[1] = MI(rMP)を代入し、
            LとMのstd::pairを返す            
            \return LとMのstd::pair
        */
        mypair getMPval() const;

        //! A public member function.
        /*!
            指定したエネルギー固有値Eで初期化を行う
            \param E 指定したエネルギー固有値
        */
        void initialize(double E);

        //! A public member function.
        /*!
            原点に近い点と無限遠に近い点から、それぞれ微分方程式を解く
            \return それぞれの微分方程式が正常に解けたかどうか
        */
        void solve_diff_equ();
        
        //! A public member function.
        /*!
            （境界値の考慮されていない）Poisson方程式を解く
        */
        void solve_poisson();

        //! A public member function (const).
        /*!
            ポテンシャルV(r)の値を返す
            \param r rの値
            \return ポテンシャルV(r)の値
        */
        double V(double r) const;

        // #endregion publicメンバ関数 

        // #region privateメンバ関数

    private:
        //! A private member function.
        /*!
            V(r)の級数展開の係数amを求める
        */
        void am_evaluate();

        //! A private member function.
        /*!
            L(r)の級数展開の係数bmを求める
        */
        void bm_evaluate();

        //! A private member function (const).
        /*!
            微分方程式の式を定義する
            \param f f[0] = L, f[1] = M
            \param dfdx dfdx[0] = dL / dx, dfdx[1] = dM / dx 
            \param x xの値
        */
        void derivs(myarray const & f, myarray & dfdx, double x) const;
        
        //! A private member function (const).
        /*!
            the radial differential equation with a full relativistic treatment

            d = alpha^2 / 2 * r /MQ * dV / dr
            dM / dx = - (2NL + 1 + d)M + 2r^2 * MQ * (V - ep) * L - d * (NL + 1 + kappa) * L
            \param L L(x)の値
            \param M M(x)の値
            \param x xの値
            \return dM / dxの値
        */
        double dM_dx_dirac(double L, double M, double x) const;

        //! A private member function (const).
        /*!
            the usual radial differential equation without relativistic corrections 

            dM/dx = -(2NL + 1)M + 2r^2(V - ep)L
            \param L L(x)の値
            \param M M(x)の値
            \param x xの値
            \return dM / dxの値
        */
        double dM_dx_sch(double L, double M, double x) const;
        
        //! A private member function (const).
        /*!
            the scalar relativistic radial differential equation
            ref: D.D.Koelling and B.N.Harmon,
                 J.Phys.C: Solid State Phys. 10, 3107 (1977).

		    d = alpha^2 / 2 * r / MQ * dV / dr
		    dM / dx = -(2 * NL + 1 + d)M +(2 * r^2 * MQ(V - ep) - d * NL)L
            \param L L(x)の値
            \param M M(x)の値
            \param x xの値
            \return dM / dxの値
        */
        double dM_dx_sdirac(double L, double M, double x) const;

        //! A private member function (const).
        /*!
            ポテンシャルの微分V'(r)の値を返す
            \param r rの値
            \return ポテンシャルの微分V'(r)の値
        */
        double dV_dr(double r) const;

        //! A private member function.
        /*!
            lo_とmo_を初期化する
        */
        void init_lm_o();

        //!  A private member function.
        /*!
            L(x)のノードの数をカウントする
            \param L L(x)の格納されたstd::vector
        */
        void node_count(dvector const & L);

        //! A private member function.
        /*!
            li_とmi_の初期値を求める
            \return li_とmi_の初期値
        */
        myarray req_lm_i_init_val();

        //! A private member function.
        /*!
            lo_とmo_の初期値を求める
            \return lo_とmo_の初期値
        */
        myarray req_lm_o_init_val();

        //!  A private member function.
        /*!
            poisson方程式の初期値を求める
            \return poisson方程式の初期値
        */
        myarray req_poisson_init_val();

        template <typename Stepper>
        //! A private member function.
        /*!
            無限遠に近い点から、微分方程式を解く
            \param stepper 微分方程式のソルバーのアルゴリズム
        */
        void solve_diff_equ_i(Stepper const & stepper);

        template <typename Stepper>
        //! A private member function.
        /*!
            原点に近い点から、微分方程式を解く
            \param stepper 微分方程式のソルバーのアルゴリズム
        */
        void solve_diff_equ_o(Stepper const & stepper);

        template <typename Stepper>
        //!  A private member function.
        /*!
            （境界値の考慮されていない）Poisson方程式を解く
        */
        void solve_poisson_run(Stepper const & stepper);

        // #endregion privateメンバ関数

        // #region プロパティ

    public:
        //! A property.
        /*!
            微分方程式のデータオブジェクトへのプロパティ
        */
        Property<std::shared_ptr<DiffData>> const PDiffData;

        // #endregion プロパティ

        // #region メンバ変数

    private:
        //!  A private static member variable (constant expression).
        /*!
            L(r)の級数展開の係数bmの最大値
        */
        static std::size_t constexpr BMMAX = 5;

        //!  A private static member variable (constant expression).
        /*!
            ポテンシャルVの最小値
        */
        static auto constexpr MINVALUE = 1.0E-200;

        //!  A private member variable.
        /*!
            V(r)の級数展開の係数am
        */
        std::array<double, AMMAX> am;

        //!  A private member variable.
        /*!
            L(r)の級数展開の係数bm
        */
        std::array<double, DiffSolver::BMMAX> bm;
        
    public:
        //! A private member variable.
        /*!
            エネルギー固有値
        */
        double E_;

        //!  A private member variable (constant).
        /*!
            データオブジェクト
        */
        std::shared_ptr<Data> const pdata_;

        //!  A private member variable (constant).
        /*!
            微分方程式のデータオブジェクト
        */
		std::shared_ptr<DiffData> const pdiffdata_;

        //!  A private member variable (constant).
        /*!
            関数ρ(r)
        */
        std::shared_ptr<Rho> const prho_;

        //!  A private member variable.
        /*!
            Hartreeポテンシャルオブジェクト
        */
        std::shared_ptr<Vhartree> pvh_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        DiffSolver() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        DiffSolver(DiffSolver const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        DiffSolver & operator=(DiffSolver const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};
    
    // #region 非メンバ関数
    
	template <typename T>
    //! A template function.
    /*!
        x ** 2を計算する
        \param x xの値
        \return x ** 2の値
    */
	T sqr(T x)
	{
		return x * x;
	}
    
    // #endregion 非メンバ関数
}

#endif	// _DIFFSOLVER_H_
