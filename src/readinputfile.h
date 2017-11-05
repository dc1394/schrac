/*! \file readinputfile.h
    \brief インプットファイルの読み込みを行うクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _READINPUTFILE_H_
#define _READINPUTFILE_H_

#pragma once

#include "ci_string.h"
#include "data.h"
#include "property.h"
#include <fstream>                  // for std::ifstream
#include <memory>                   // for std::shared_ptr
#include <optional>					// for std::optional
#include <vector>                   // for std::vector
#include <boost/lexical_cast.hpp>   // for boost::lexical_cast

namespace schrac {
    //! A class.
    /*!
        インプットファイルを読み込み、Dataクラスのオブジェクトに格納するクラス
    */
    class ReadInputFile final {
        // #region 型エイリアス

        using strvec = std::vector < ci_string > ;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param arg インプットファイル名と、TBBを使用するかどうかのstd::pair
        */
        explicit ReadInputFile(std::pair<std::string, bool> const & arg);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~ReadInputFile() = default;

        // #region コンストラクタ・デストラクタ
        
        // #region メンバ関数

    public:
        //! A public member function.
        /*!
            ファイルを読み込む
        */
        void readFile();

    private:
        //! A private member function (const).
        /*!
            エラーを表示する
            \param s エラーのトークン
            \return 読み込みが成功したかどうか
        */
        void errorMessage(ci_string const & s) const;

        //! A private member function (const).
        /*!
            エラーを表示する
            \param line エラーのある行
            \param s1 エラーのトークン
            \param s2 エラーのトークン2
        */
        void errorMessage(std::int32_t line, ci_string const & s, ci_string const & s2) const;

        //! A private member function (const).
        /*!
            解析対象の文字列をトークンに分割する
            \param article 解析対象の文字列
            \return 関数が成功したかどうかと、トークンのstd::pair
        */
        std::pair<std::int32_t, std::optional<ReadInputFile::strvec>> getToken(ci_string const & article);

        //! A private member function.
        /*!
            原子に関するデータを読み込む
            \return 読み込みが成功したかどうか
        */
        bool readAtom();

        //! A private member function.
        /*!
            文字列を解析して、データとして読み込んで返す
            \param article 解析対象の文字列
            \return 読みこんだ文字列データ
        */
        std::optional<ci_string> readData(ci_string const & article);
        
        //! A private member function.
        /*!
            データを読み込む
            \param article 解析対象の文字列
            \param def デフォルトの文字列
            \return 読みこんだ文字列
        */
        std::optional<ci_string> readData(ci_string const & article, ci_string const & def);

        template <typename T>
        //! A private member function (template function).
        /*!
            データを読み込む
            \param article 解析対象の文字列
            \param default_value デフォルト値
            \return 読みこんだ文字列
        */
        std::optional<T> readData(ci_string const & article, T const & default_value);
        
        //! A private member function.
        /*!
            データを読み込む
            \param article 解析対象の文字列
            \return 読みこんだ文字列（読み込みに失敗したならstd::nullopt）
        */
        std::optional<ci_string> readDataAuto(ci_string const & article);

        //! A private member function.
        /*!
            対象の方程式を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readEq();

        //! A private member function.
        /*!
            SCFの一次混合の重みを読み込む
            \return 読み込みが成功したかどうか
        */
        bool readScfMixingWeight();

        //! A private member function.
        /*!
            微分方程式の解法を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readType();

        template <typename T>
        //! A private member function.
        /*!
            対象の要素の値をその行から読み込む
            \param article 要素名
            \param default_value デフォルトの値
            \param value 読み込んだ値
        */
        void readValue(ci_string const & article, T const & default_value, T & value);

        template <typename T>
        //! A private member function.
        /*!
            対象の要素の値をその行から読み込む
            \param article 要素名
            \param value 読み込んだ値
            \return 読み込みが成功したかどうか
        */
        bool readValueAuto(ci_string const & article, std::optional<T> & value);

        // #endregion メンバ関数

        // #region プロパティ

    public:
        //! A property.
        /*!
            読み込んだデータを返す
            \return 読み込んだデータ
        */
        Property<std::shared_ptr<Data>> const PData;

        // #endregion プロパティ

        // #region メンバ変数

        //! A private member variable (constant expression).
        /*!
            バッファサイズ
        */
        static std::streamsize constexpr BUFSIZE = 1024;

        //! A private member variable (constant).
        /*!
            「chemical.symbol」の文字列
        */
        static const ci_string CHEMICAL_SYMBOL;

        //! A private member variable (constant).
        /*!
            デフォルトの「eq.type」の文字列
        */
        static const ci_string EQ_TYPE_DEFAULT;

        //! A private member variable (constant).
        /*!
            「eq.type」の文字列
        */
        static const ci_string EQ_TYPE;

        //! A private member variable (constant).
        /*!
            方程式の種類の文字列の配列
        */
        static const std::array<ci_string, 4> EQ_TYPE_ARRAY;

        //! A private member variable (constant).
        /*!
            「orbital」の文字列
        */
        static const ci_string ORBITAL;

        //! A private member variable (constant).
        /*!
            微分方程式の数値解法の文字列の配列
        */
        static const std::array<ci_string, 4> SOLVER_TYPE_ARRAY;

        //! A private member variable (constant).
        /*!
            デフォルトの微分方程式の数値解法
        */
        static const ci_string SOLVER_TYPE_DEFAULT;

        //! A private member variable (constant).
        /*!
            「spin_orbital」の文字列
        */
        static const ci_string SPIN_ORBITAL;

        //! A private member variable (constant).
        /*!
        デフォルトの「spin_orbital」の文字列
        */
        static const ci_string SPIN_ORBITAL_DEFAULT;

        //! A private member variable.
        /*!
            ファイル読み込み用のストリーム
        */
        std::ifstream ifs_;
        
        //! A private member variable.
        /*!
            現在の行数
        */
        std::size_t lineindex_;

        //! A private member variable.
        /*!
            インプットファイルから読み込んだデータ
        */
        std::shared_ptr<Data> pdata_;
            
        // #endregion メンバ変数
        
        // #region 禁止されたコンストラクタ・メンバ関数

    private:
        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        ReadInputFile() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        ReadInputFile(ReadInputFile const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        ReadInputFile & operator=(ReadInputFile const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };

    template <typename T>
    std::optional<T> ReadInputFile::readData(ci_string const & article, T const & default_value)
    {
        // グリッドを読み込む
        for (; true; lineindex_++) {
            auto const ret = getToken(article);

            switch (std::get<0>(ret))
            {
            case -1:
                return std::nullopt;
                break;

            case 0:
            {
                auto const tokens = *(std::get<1>(ret));
                auto itr(++tokens.begin());

                lineindex_++;
                // 読み込んだトークンの数をはかる
                switch (tokens.size()) {
                case 1:
                    // デフォルト値を返す
                    return std::make_optional<T>(default_value);
                    break;

                case 2:
                    if (*(++itr) == "DEFAULT") {
                        // デフォルト値を返す
                        return std::make_optional<T>(default_value);
                    }
                    else {
                        try {
                            return std::make_optional<T>(boost::lexical_cast<T>(itr->c_str()));
                        }
                        catch (boost::bad_lexical_cast const &) {
                            errorMessage(lineindex_ - 1, article, *itr);
                            return std::nullopt;
                        }
                    }

                default:
                {
                    auto val = *itr;

                    if (val == "DEFAULT" || val[0] == '#') {
                        return std::make_optional<T>(default_value);
                    }
                    else if ((*(++itr))[0] != '#') {
                        errorMessage(lineindex_ - 1, article, *itr);
                        return std::nullopt;
                    }

                    try {
                        return std::make_optional<T>(boost::lexical_cast<T>(val.c_str()));
                    }
                    catch (boost::bad_lexical_cast const &) {
                        errorMessage(lineindex_ - 1, article, val);
                        return std::nullopt;
                    }
                }
                }
            }
                break;

            case 1:
                break;

            default:
                BOOST_ASSERT(!"何かがおかしい!");
                break;
            }
        }
    }

    template <typename T>
    void ReadInputFile::readValue(ci_string const & article, T const & default_value, T & value)
    {
        if (auto const p = readData(article, default_value)) {
            value = *p;
        }
        else {
            throw std::runtime_error("インプットファイルが異常です");
        }
    }

    template <typename T>
    bool ReadInputFile::readValueAuto(ci_string const & article, std::optional<T> & value)
    {
        if (auto const val = readDataAuto(article)) {
            if (!val->empty()) {
                try {
                    std::size_t idx;
                    auto const v = std::stod(val->c_str(), &idx);
                    if (idx != val->length()) {
                        throw std::invalid_argument("");
                    }
                    value = std::make_optional<double>(v);
                }
                catch (std::invalid_argument const &) {
                    errorMessage(lineindex_ - 1, article, *val);
                    return false;
                }
            }
            else {
                value = std::nullopt;
            }
        }
        else {
            return false;
        }

        return true;
    }
}

#endif  // _READINPUTFILE_H_
