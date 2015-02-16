/*! \file Functional.h
    \brief std::function<double (double)>の代わりになるクラス

    Copyright ©  2014 @dc1394 All Rights Reserved.
*/
#ifndef _FUNCTONAL_H_
#define _FUNCTONAL_H_

#pragma once

namespace myfunctional {
    //! A template class.
    /*!
        std::function<double (double)>の代わりになるtemplate class
    */
    template <typename FUNCTYPE>
    class Functional final
    {
    public:
        // #region コンストラクタ

        //! A constructor.
        /*!
        \param func operator()で呼び出す関数
        */
        Functional(FUNCTYPE const & func) : func_(func) {}

        // #endregion コンストラクタ

        // #region メンバ関数

        //! A public member function.
        /*!
            operator()の宣言と実装
            関数f(x)の値を返す
            \param x xの値
            \return f(x)の値
        */
        double operator()(double x) const
        {
            return func_(x);
        }

        // #endregion メンバ関数

    private:
        // #region メンバ変数
        
        //! A private const variable (reference).
        /*!
            operator()で呼び出す関数
        */
        FUNCTYPE const & func_;

        // #endregion メンバ変数
    };

    //! A template function（非メンバ関数）.
    /*!
        Function<FUNCTYPE>を生成する
        \param func 格納する関数
        \return 生成されたFunction<FUNCTYPE>
    */
    template <typename FUNCTYPE>
    Functional<FUNCTYPE> make_functional(FUNCTYPE const & func)
    {
        return Functional<FUNCTYPE>(func);
    }
}

#endif  // _FUNCTIONAL_H_
