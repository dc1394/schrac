/*! \file eigenvaluesearch.h
    \brief エネルギー固有値検索を行うクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _EIGENVALUESEARCH_H_
#define _EIGENVALUESEARCH_H_

#pragma once

#include "diffsolver.h"

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
        EigenValueSearch(std::shared_ptr<Data> const & pdata, std::shared_ptr<DiffData> const & pdiffdata, std::shared_ptr<Rho> const & prho, std::shared_ptr<Vhartree> const & pvh);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~EigenValueSearch() = default;

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
        Property<std::shared_ptr<DiffSolver>> const PDiffSolver;

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
            \param D 関数Dの値
            \param E 関数Dの引数E
        */
        void info(double D, double E) const;
        
        //! A private member function.
        /*!
            状態の初期化を行う
            \param prho Rhoオブジェクトへのスマートポインタ
        */
        void initialize(std::shared_ptr<Rho> const & prho);
        
        //! A private member function.
        /*!
            固有値を大まかに検索する
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
        static constexpr auto EVALSEARCHMAX = 10000;

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
        std::shared_ptr<Data> const pdata_;

        //! A private member variable.
        /*!
            微分方程式オブジェクト
        */
        std::shared_ptr<DiffSolver> pdiffsolver_;

        //! A private member variable.
        /*!
            微分方程式データのオブジェクト
        */
        std::shared_ptr<DiffData> const pdiffdata_;

        //!  A private member variable.
        /*!
            Hartreeポテンシャルオブジェクト
        */
        std::shared_ptr<Vhartree> pvh_;
        
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
        \param params DiffSolverオブジェクトのポインタを無理矢理Void *にキャスト
        \return 関数Dの値
    */
    double func_D(double E, void * params);

    //! A function.
    /*!
        対象の方程式がDirac方程式の場合に、エネルギー固有値の近似値を求める
        \param pdata データオブジェクト
        \return 大体のエネルギー固有値
    */
    double Eapprox_dirac(std::shared_ptr<Data> const & pdata);
    
    //! A function.
    /*!
        対象の方程式がSch方程式の場合に、エネルギー固有値の近似値を求める
        \param pdata データオブジェクト
        \return 大体のエネルギー固有値
    */
    double Eapprox_sch(std::shared_ptr<Data> const & pdata);

    // #endregion 非メンバ関数
}

#endif // _EIGENVALUESEARCH_H_
