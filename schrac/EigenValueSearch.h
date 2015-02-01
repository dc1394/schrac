﻿/*! \file eigenvaluesearch.h
    \brief エネルギー固有値検索を行うクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _EIGENVALUESEARCH_H_
#define _EIGENVALUESEARCH_H_

#pragma once

#include "diff.h"
#include "readinputfile.h"

namespace schrac {
    //! A class.
    /*!
        エネルギー固有値検索を行うクラス
    */
	class EigenValueSearch final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param arg インプットファイル名とTBBを使用するかどうかのstd::pair
        */
        explicit EigenValueSearch(std::pair<std::string, bool> const & arg);

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~EigenValueSearch()
        {
        }

        // #endregion コンストラクタ・デストラクタ

        // #region プロパティ

        //! A property.
        /*!
            データオブジェクトを得る
            \return データオブジェクト
        */
        Property<std::shared_ptr<Data>> const PData;

        //! A property.
        /*!
            微分方程式オブジェクトを得る
            \return 微分方程式オブジェクト
        */
        Property<std::shared_ptr<Diff>> const PDiff;

        // #endregion プロパティ

        // #region メンバ関数

    public:
        //! A public member function.
        /*!
            固有値を検索する
            \return 固有値が見つかったかどうか
        */
        bool search();

    private:
        //! A private member function.
        /*!
            Brent法で、関数Dの根を計算する
            \return 根が見つかったかどうか
        */
        bool brent();

        //! A private member function (const).
        /*!
            現在のループをメッセージで報告する
        */
        void info() const;

        //! A private member function (const).
        /*!
            固有値が見つかったことをメッセージで報告する
        */
        void info(double E) const;

        //! A private member function (const).
        /*!
            固有値が見つかったことをメッセージで報告する
            \param b 関数Dの引数b
            \param fb 関数Dの引数fb
        */
        void info(double b, double fb) const;

        //! A private member function (const).
        /*!
            状態の初期化を行う
        */
        void initialize();

        //! A private member function (const).
        /*!
            解く微分方程式についてメッセージを表示する
        */
        void msg() const;

        //! A private member function.
        /*!
            固有値をおおざっぱに検索する
            \return 固有値が見つかったかどうか
        */
        bool rough_search();

        //! A private member function (const).
        /*!
            表示する浮動小数点の桁を設定する
        */
        void setoutstream() const;
        
        // #endregion メンバ関数

        // #region メンバ変数

        //! A private member variable (constant expression).
        /*!
            エネルギー固有値探索の最大のループ回数
        */
		static constexpr auto EVALSEARCHMAX = 1000;

        //! A private member variable (constant expression).
        /*!
            閾値（絶対値の大きい方）
        */
        static constexpr auto HUGE = 1.0E+7;
        
    public:
        //! A public static member variable.
        /*!
            固有関数のノードが一致しているかどうか
        */
        static bool nodeok;
           
    private:
        //! A private member variable.
        /*!
            エネルギー固有値探索の幅
        */
        double DE_;

        //! A private member variable.
        /*!
            エネルギー固有値
        */
        double E_;

        //! A private member variable.
        /*!
            Brent法におけるエネルギー固有値の大きい方
        */
        double Emax_;

        //! A private member variable.
        /*!
            Brent法におけるエネルギー固有値の小さい方
        */
        double Emin_;

        //! A private member variable.
        /*!
            エネルギー固有値の近似値
        */
        double Eapprox_;

        //! A private member variable.
        /*!
            関数Dの古い値
        */
        double Dold;
		
        //! A private member variable.
        /*!
            エネルギー固有値探索のループ回数
        */
        std::int32_t loop_;
        
        //! A private member variable.
        /*!
            インプットファイルのデータオブジェクト
        */
		std::shared_ptr<Data> pdata_;

        //! A private member variable.
        /*!
            微分方程式オブジェクト
        */
		std::shared_ptr<Diff> pdiff_;

        //! A private member variable.
        /*!
            微分方程式データのオブジェクト
        */
		std::shared_ptr<DiffData> pdiffdata_;
        
        // #endregion メンバ変数

    private:
        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
        デフォルトコンストラクタ（禁止）
        */
        EigenValueSearch() = delete;

        //! A private copy constructor (deleted).
        /*!
        コピーコンストラクタ（禁止）
        */
        EigenValueSearch(EigenValueSearch const &) = delete;

        //! A private member function (deleted).
        /*!
        operator=()の宣言（禁止）
        \param コピー元のオブジェクト（未使用）
        \return コピー元のオブジェクト
        */
        EigenValueSearch & operator=(EigenValueSearch const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};

    // #region 非メンバ関数

    //! A function.
    /*!
        関数Dの値を求める
        \param E エネルギー固有値
        \param params Diffオブジェクトのポインタを無理矢理Void *にキャスト
        \return 関数Dの値
    */
    double func_D(double E, void * params);

    //! A function.
    /*!
        対象の方程式がDirac方程式の場合に、エネルギー固有値の下限を概算する
        \param pdata データオブジェクト
        \return 大体のエネルギー固有値
    */
    double Eexact_dirac(std::shared_ptr<Data> const & pdata);
    
    //! A function.
    /*!
        対象の方程式がSch方程式の場合に、エネルギー固有値の下限を概算する
        \param pdata データオブジェクト
        \return 大体のエネルギー固有値
    */
    double Eexact_sch(std::shared_ptr<Data> const & pdata);

    //! A function.
    /*!
        対象の方程式がscalar Dirac方程式の場合に、エネルギー固有値の下限を概算する
        \param pdata データオブジェクト
        \return 大体のエネルギー固有値
    */
    double Eexact_sdirac(std::shared_ptr<Data> const & pdata);

	template <typename T>
    //! A function (template function).
    /*!
        bが正の値の場合にaの絶対値を、bが負の値の場合はaの絶対値に-をかけた値を返す
        \param a 対象の値
        \param b 正負を判断するための値
        \return bが正の値の場合はaの絶対値、bが負の値の場合はaの絶対値に-をかけた値
    */
	T sign(T a, T b)
	{
		return (b >= 0.0) ? std::fabs(a) : - std::fabs(a);
	}

    // #endregion 非メンバ関数
}

#endif // _EIGENVALUESEARCH_H_
