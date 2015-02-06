/*! \file pow.h
    \brief x ** nを計算するtemplate関数powの宣言と実装

    Copyright © 2015 @dc1394 All Rights Reserved.
*/

#ifndef _POW_H_
#define _POW_H_

#pragma once

#include <cstdint>  // for std::int32_t

namespace schrac {
	template <typename T>
    //! A template function.
    /*!
        x ** nを計算する
        \param x xの値
        \param n nの値
        \return x ** nの値 
    */
    T pow(T x, std::int32_t n)
    {
        T p = x, y = 1.0;

        while (true) {
            if (n & 1) {
                y *= p;
            }

            n >>= 1;

            if (!n) {
                return y;
            }

            p *= p;
        }
    }
}

#endif	// _POW_H_
