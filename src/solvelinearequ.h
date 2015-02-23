/*! \file solvelinearequ.h
    \brief 連立一次方程式を解く非メンバ関数の宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _SOLVELINEAREQU_H_
#define _SOLVELINEAREQU_H_

#pragma once

#include <array>
#include <gsl/gsl_linalg.h>	// for gsl_linalg

namespace schrac {
	//!  A static variable (constant expression).
    /*!
        級数展開の係数am_の最大値
    */
    static std::size_t constexpr AMMAX = 3;

    //! #region 型エイリアス
    /*!
    */
    using myvector = std::array < double, AMMAX > ;

    //! A function.
    /*!
        連立一次方程式を解く
        \param a 連立一次方程式Ax = bにおける左辺の行列A
        \param b 連立一次方程式Ax = bにおける右辺のベクトルb
        \return 方程式の解ベクトル
    */
    myvector solve_linear_equ(std::array<double, AMMAX * AMMAX> & a, myvector & b);
}

#endif  // _SOLVELINEAREQU_H_
