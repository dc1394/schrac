/*! \file wavefunction.h
    \brief 得られた波動関数をファイルに書き出すクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#ifndef _WAVEFUNCTION_H_
#define _WAVEFUNCTION_H_

#pragma once

#include "data.h"
#include <memory>                       // for std::shared_ptr
#include <utility>                      // for std::pair
#include <vector>                       // for std::vector
#include <boost/container/flat_map.hpp> // for boost::container::flat_map
#include <boost/optional.hpp>           // for boost::optional

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
            \param wf 波動関数が格納されたハッシュ
        */
        WaveFunctionSave(boost::container::flat_map<std::string, std::vector<double>> const & wf, std::shared_ptr<Data> const & pdata);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~WaveFunctionSave() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数
    
    private:
        //!  A private member function.
        /*!
            スピン軌道を表す文字列を返す関数
            \return スピン軌道を表す文字列が格納されたboost::optional
        */
        boost::optional<std::string> get_spin_orbital() const;

        //!  A private member function.
        /*!
            ファイル名を生成する関数
            \return 波動関数のファイル名と電子密度のファイル名のstd::pair
        */
	    std::pair<std::string, std::string> make_filename() const;

	public:
        //!  A public member function.
        /*!
            実際にファイル出力処理を実行する関数
            \return ファイル出力処理が成功したかどうか
        */
		bool operator()();

        // #endregion メンバ関数

    private:
        // #region メンバ変数
        
        //!  A private member variable.
        /*!
            波動関数が格納されたboost::container::flat_map
        */
        boost::container::flat_map<std::string, std::vector<double>> wf_;

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
