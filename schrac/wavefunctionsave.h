/*! \file wavefunction.h
    \brief 得られた波動関数をファイルに書き出すクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _WAVEFUNCTION_H_
#define _WAVEFUNCTION_H_

#pragma once

#include "data.h"
#include <memory>
#include <unordered_map>

namespace schrac {
    //! A class.
    /*!
        得られた波動関数をファイルに書き出すクラス
    */
	class WaveFunctionSave final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param pdata データオブジェクト
            \param hash 波動関数が格納されたハッシュ
        */
        WaveFunctionSave(std::unordered_map<std::string, std::vector<double>> const & myhash, std::shared_ptr<Data> const & pdata);

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~WaveFunctionSave()
        {
        }

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

	    std::string make_filename() const;

	public:
		bool operator()();

        // #endregion メンバ関数

    private:
        // #region メンバ変数
        
        //!  A private member variable.
        /*!
            波動関数が格納されたhash
        */
        std::unordered_map<std::string, std::vector<double>> myhash_;

        //!  A private member variable.
        /*!
            データオブジェクト
        */
        std::shared_ptr<Data> const pdata_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        WaveFunctionSave() = delete;

        //! A private copy constructor (deleted).
        /*!
        コピーコンストラクタ（禁止）
        */
        WaveFunctionSave(WaveFunctionSave const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        WaveFunctionSave & operator=(WaveFunctionSave const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};
}

#endif  // _WAVEFUNCTION_H_
