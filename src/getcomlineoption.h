/*! \file getcomlineoption.h
    \brief コマンドラインオプションの解析を行うクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _GETCOMLINEOPTION_H_
#define _GETCOMLINEOPTION_H_

#pragma once

#include <cstdint>  // for std::int32_t
#include <string>   // for std::string
#include <utility>  // for std::pair

namespace schrac {
    //! A class.
    /*!
        コマンドラインオプションを解析するクラス    
    */
    class GetComLineOption final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            何もしないコンストラクタ
        */
        GetComLineOption() = default;

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~GetComLineOption() = default;
        
        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //! A public member function.
        /*!
            コマンドラインオプションを解析する
            \param argc コマンドライン引数の数
            \param argv コマンドライン引数
            \return 解析に成功したら0、失敗したら-1
        */
        std::int32_t getopt(int argc, char * const argv[]);

        //! A public member function (constant).
        /*!
            インプットファイル名とTBBを使用するかどうかを、std::pairで返す
            \return インプットファイル名とTBBを使用するかどうかのstd::pair
        */
        std::pair<std::string, bool> getpairdata() const;

        // #endregion メンバ関数

    private:
        // #region メンバ変数

        //!  A private static member variable (constant).
        /*!
            デフォルトのインプットファイル名
        */
        static std::string const DEFINPNAME;
        
        //!  A private member variable.
        /*!
            デフォルトのインプットファイル名
        */
        std::string inpname_;

        //!  A private member variable.
        /*!
            TBBを使用するかどうか    
        */
        bool usetbb_ = false;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        GetComLineOption(GetComLineOption const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        GetComLineOption & operator=(GetComLineOption const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif  // _GETCOMLINEOPTION_H_
