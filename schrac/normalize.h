/*! \file normalize.h
    \brief 得られた波動関数を正規化するクラスの宣言と実装

    Copyright © 2015 @dc1394 All Rights Reserved.
*/

#ifndef _NORMALIZE_H_
#define _NORMALIZE_H_

#pragma once

#include "eigenvaluesearch.h"
#include <unordered_map>        // for std::unordered_map

namespace schrac {
    template <typename Derived>
    //! A template class.
    /*! 
        得られた波動関数を正規化するクラス
    */
	class Normalize {
        // #region 型エイリアス

    public:
        using myhash = std::unordered_map<std::string, dvector>;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdiffsolver 微分方程式のデータオブジェクト
        */
        Normalize(std::shared_ptr<DiffSolver> const & pdiffsolver);

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~Normalize()
        {
        }

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

    public:
        //! A public member function.
        /*!
            波動関数を求める
        */
        void base_evaluate();

        //! A public member function.
        /*!
            結果を返す
        */
        myhash base_getresult() const;

    private:
        //! A public member function.
        /*!
            波動関数を正規化する
        */
        void base_normalize();

        // #endregion メンバ関数

        // #region メンバ変数

    protected:
        //! A protected member variable.
        /*!
            データオブジェクト
        */
		std::shared_ptr<Data> const pdata_;

        //! A protected member variable.
        /*!
            微分方程式のデータオブジェクト
        */
        std::shared_ptr<DiffData> const pdiffdata_;

        //! A protected member variable.
        /*!
            微分方程式オブジェクト
        */
        std::shared_ptr<DiffSolver> const pdiffsolver_;

        //! A protected member variable.
        /*!
            rのメッシュが格納されたstd::vector
        */
        dvector r_mesh_;
        
        // #endregion メンバ変数 

    private:
        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
        デフォルトコンストラクタ（禁止）
        */
        Normalize() = delete;

        //! A private copy constructor (deleted).
        /*!
        コピーコンストラクタ（禁止）
        */
        Normalize(Normalize const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Normalize & operator=(Normalize const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};

    // #region コンストラクタの実装

    template <typename Derived>
    Normalize<Derived>::Normalize(std::shared_ptr<DiffSolver> const & pdiffsolver) :
        pdata_(pdiffsolver->PDiffData()->pdata_),
        pdiffdata_(pdiffsolver->PDiffData),
        pdiffsolver_(pdiffsolver)
    {
        r_mesh_.reserve(pdata_->grid_num_);

        auto const & rmesh_o(pdiffdata_->r_mesh_o_);

        r_mesh_.assign(rmesh_o.begin(), rmesh_o.end());
    }

    // #endregion コンストラクタの実装

    // #region protectedメンバ関数の実装

    template <typename Derived>
    void Normalize<Derived>::base_evaluate()
    {
        static_cast<Derived &>(*this).evaluate();
    }

    template <typename Derived>
    Normalize<Derived>::myhash Normalize<Derived>::base_getresult() const
    {
        return static_cast<Derived &>(*this).getresult();
    }

    template <typename Derived>
    void Normalize<Derived>::base_normalize()
    {
        static_cast<Derived &>(*this).normalize();
    }

    // #region protectedメンバ関数の実装
}

#endif // _NORMALIZE_H_
