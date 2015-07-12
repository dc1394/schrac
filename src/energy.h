/*! \file energy.h
    \brief エネルギーを計算するクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
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
            \param Z 原子核の電荷
        */
        Energy(std::shared_ptr<DiffData> const & pdiffdata, dvector const & r, dvector const & rf, double Z);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Energy() = default;

        // #region メンバ関数
        
        //! A public member function (const).
        /*!
            エネルギーを表示する
        */
        void express_energy(boost::optional<double> const & ehartree) const;

    private:
        //! A private member function (const).
        /*!
            Coulombエネルギーを表示する
        */
        void coulomb_energy() const;

        //! A private member function (const).
        /*!
            固有エネルギーを表示する
        */
        void eigenvalue() const;

        //! A private member function (const).
        /*!
            ポテンシャルエネルギーを表示する
        */
        void hartree_energy(double ehartree) const;
        
        //! A private member function (const).
        /*!
            運動エネルギーを表示する
        */
        void kinetic_energy() const;

        //! A private member function (const).
        /*!
            ポテンシャルエネルギーを表示する
        */
        void potential_energy() const;

        //! A private member function (const).
        /*!
            全エネルギーを表示する
        */
        void total_energy(boost::optional<double> const & ehartree) const;

        // #endregion メンバ関数

        // #region メンバ変数
    
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
        double const potcoulomb_energy_;

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
