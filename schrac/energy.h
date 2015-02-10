/*! \file energy.h
    \brief エネルギーを計算するクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _ENERGY_H_
#define _ENERGY_H_

#pragma once

#include "diffdata.h"
#include "simpson.h"
#include <memory>

namespace schrac {
    //! A class.
    /*!
        エネルギーを計算するクラス
    */
    class Energy final {
        // #region 型エイリアス

        using dvector = std::vector < double > ;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdiffdata 微分方程式のデータオブジェクト
            \param pf 規格化された波動関数
        */
        Energy(std::shared_ptr<DiffData> const & pdiffdata, dvector const & rf, dvector const & r);

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~Energy()
        {
        }

        // #region メンバ関数

        //! A public member function (const).
        /*!
            固有エネルギーを表示する
        */
        void eigenvalue() const;

        //! A public member function (const).
        /*!
            運動エネルギーを表示する
        */
        void kinetic_energy() const;

        //! A public member function (const).
        /*!
            ポテンシャルエネルギーを表示する
        */
        void potential_energy() const;

        // #endregion メンバ関数

        // #region メンバ変数

    private:
        //! A private member variable (const).
        /*!
            データオブジェクト
        */
        std::shared_ptr<DiffData> const pdiffdata_;

        //! A private member variable (const).
        /*!
            規格化された波動関数
        */
        dvector const rf_;
        
        //! A private member variable (const).
        /*!
            ポテンシャルエネルギー
        */
        double const potential_energy_;

        // #endregion メンバ変数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Energy() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        Energy(Energy const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Energy & operator=(Energy const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };        
}

#endif  // _ENERGY_H_
