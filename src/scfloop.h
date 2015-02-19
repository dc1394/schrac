﻿/*! \file scfloop.h
    \brief SCFを行うクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _SCFLOOP_H_
#define _SCFLOOP_H_

#pragma once

#include "diffsolver.h"

namespace schrac {
	class ScfLoop final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param arg インプットファイル名とTBBを使用するかどうかのstd::pair
        */
        explicit ScfLoop(std::pair<std::string, bool> const & arg);
        
        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~ScfLoop()
        {
        }

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        bool check_converge();
                
        //! A public member function (const).
        /*!
            対象の原子と解く方程式についてメッセージを表示する
        */
        void message() const;
        
        //! A public member function.
        /*!
            実際にSCFを実行する    
        */
        void operator()();

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function.
        /*!
            状態の初期化を行う
        */
        void initialize();
        
        //! A private member function.
        /*!
            Hartreeポテンシャルを生成する
        */
        void make_vhartree();

        //! A private member function.
        /*!
            与えられた固有値、密度及びHartreeポテンシャルから全エネルギーを求める
            \param eigen 固有値
            \param rho 密度
            \param vhartree Hartreeポテンシャル
        */
        double req_energy(double eigen, dvector const & rho, dvector const & vhartree) const;

        //! A private member function.
        /*!
            残差ノルムを求める
            \param newrho 新しい密度ρnew
            \param oldrho 密度ρ
            \return 残差
        */
        double req_normrd(dvector const & newrho, dvector const & oldrho) const;

        //! A private member function.
        /*!
            RF(r)から、新しい密度ρnew(r)を求める
            \param rf RF(r)
            \return 密度ρnew(r)
        */
        dvector req_newrho(dvector const & rf) const;

        // #endregion privateメンバ関数

        // #region プロパティ

    public:
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
        Property<std::shared_ptr<DiffData>> const PDiffData;

        // #endregion プロパティ

        // #region メンバ変数

    private:
        //!  A private member variable (constant).
        /*!
            データオブジェクト
        */
        std::shared_ptr<Data> pdata_;

        //!  A private member variable (constant).
        /*!
            微分方程式のデータオブジェクト
        */
        std::shared_ptr<DiffData> pdiffdata_;
        
        //!  A private member variable (constant).
        /*!
            微分方程式のソルバー
        */
        std::shared_ptr<DiffSolver> pdiffsolver_;

        //!  A private member variable (constant).
        /*!
            関数ρ(r)
        */
        std::shared_ptr<Rho> prho_;

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
        ScfLoop() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        ScfLoop(ScfLoop const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        ScfLoop & operator=(ScfLoop const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};

	/*inline void SCFLoop::setExp() const
	{
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout << std::setprecision(
			boost::numeric_cast<const std::streamsize>(std::fabs(std::log10(pdata_->eps))) - 2);
	}*/
}

#endif  // _SCFLOOP_H_