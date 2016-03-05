/*! \file fastarenaobject.h
    \brief 指定された型の指定された要素数のメモリを確保するクラス

    Copyright ©  2014 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _FASTARENAOBJECT_H_
#define _FASTARENAOBJECT_H_

#pragma once

#include "arraiedallocator.h"
#include <ccomplex>

namespace checkpoint {
    //! A template class.
    /*!
        指定された型の指定された要素数のメモリを確保するクラス
        \param TTypeSize 収納する型のサイズ
        \param TnumArray 収納する要素の数
    */
    template <size_t TTypeSize, size_t TNumArray = 1>
    struct FastArenaObject final
    {
        // サイズは絶対０より大きくなくちゃダメ
        BOOST_STATIC_ASSERT(TNumArray > 0);

        // #region メンバ関数

        //! A public member function.
        /*!
            operator newの宣言と実装
            \param 未使用
        */
        static void * operator new(std::size_t) {
            return ArraiedAllocator<TTypeSize, TNumArray>::GetAllocator().Alloc();
        }

        //! A public member function.
        /*!
            operator deleteの宣言と実装
            \param p 解放するメモリの先頭アドレス
        */
        static void operator delete(void * p) {
            ArraiedAllocator<TTypeSize, TNumArray>::GetAllocator().Free(p);
        }

    private:
        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        FastArenaObject() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        FastArenaObject(FastArenaObject const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト（未使用）
        */
        FastArenaObject & operator=(FastArenaObject const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif // _FASTARENAOBJECT_H_
