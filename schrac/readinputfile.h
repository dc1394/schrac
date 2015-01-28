#ifndef _READINPUTFILE_H_
#define _READINPUTFILE_H_

#include "ci_string.h"
#include "data.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>

namespace schrac {
    //! A class.
    /*!
        インプットファイルを読み込み、Dataクラスのオブジェクトに格納するクラス    
    */
	class ReadInputFile final {
        // #region 型エイリアス

        using strvec = std::vector<ci_string>;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
        */
        ReadInputFile(std::pair<std::string, bool> const & arg);

        //! A destructor.
        /*!
        何もしないデストラクタ
        */
        ~ReadInputFile()
        {
        }

        // #region コンストラクタ・デストラクタ
        
        // #region メンバ関数

    public:
        //! A public member function.
        /*!
            読み込んだデータを返す
            \return 読み込んだデータ
        */
        std::shared_ptr<Data> && getpData();
        
        //! A public member function.
        /*!
            ファイルを読み込む
            \return 読み込みが成功したかどうか
        */
        void readFile();

    private:

        //! A private member function (const).
        /*!
            エラーを表示する
            \param s エラーのトークン
            \return 読み込みが成功したかどうか
        */
        void errMsg(ci_string const & s) const;

        //! A private member function (const).
        /*!
            エラーを表示する
            \param line エラーのある行
            \param s1 エラーのトークン
            \param s2 エラーのトークン2
        */
        void errMsg(std::int32_t line, ci_string const & s, ci_string const & s2) const;

        //! A private member function (const).
        /*!
            解析対象の文字列をトークンに分割する
            \param article 解析対象の文字列
            \return 関数が成功したかどうかと、トークンのstd::pair
        */
        std::pair<std::int32_t, boost::optional<ReadInputFile::strvec>> getToken(ci_string const & article);

        //! A private member function.
        /*!
            原子名を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readAtom();

        //! A private member function.
        /*!
            文字列を解析して、データとして読み込んで返す
            \param article 解析対象の文字列
            \return 読みこんだ文字列データ
        */
        boost::optional<ci_string> readData(ci_string const &);
        
        //! A private member function.
        /*!
            データを読み込む
            \param article 解析対象の文字列
            \param def デフォルトの文字列
            \return 読みこんだ文字列
        */
        boost::optional<ci_string> readData(ci_string const & article, ci_string const & def);
        
        //! A private member function (template function).
        /*!
            データを読み込む
            \param article 解析対象の文字列
            \param def_val デフォルト値
            \return 読みこんだ文字列
        */
        template <typename T>
        boost::optional<T> readData(ci_string const & article, T const & def_val);
        
        //! A private member function.
        /*!
            データを読み込む
            \param article 解析対象の文字列
            \return 読みこんだ文字列（読み込みに失敗したならboost::none）
        */
        boost::optional<ci_string> readDataAuto(ci_string const & article);

        //! A private member function.
        /*!
            許容誤差を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readEps();

        //! A private member function.
        /*!
            対象の方程式を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readEq();

        //! A private member function.
        /*!
            メッシュの設定を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readGrid();

        //! A private member function.
        /*!
            固有値の検索を始める値を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readLowerE();

        //! A private member function.
        /*!
            原子番号の値を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readNumofp();

        //! A private member function.
        /*!
            微分方程式の解法を読み込む
            \return 読み込みが成功したかどうか
        */
        bool readType();

        //! A private member function.
        /*!
            マッチングポイントを読み込む
            \return 読み込みが成功したかどうか
        */
        bool readRatio();

        // #endregion メンバ関数

        // #region メンバ変数

        //! A private member variable (constant expression).
        /*!
            バッファサイズ
        */
        static constexpr std::streamsize BUFSIZE = 1024;

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
    boost::optional<T> ReadInputFile::readData(ci_string const & article, T const & def_val)
    {
        // グリッドを読み込む
        for (; true; lineindex_++) {
            auto const ret = getToken(article);

            switch (std::get<0>(ret))
            {
            case -1:
                return nullptr;
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
                    return boost::optional<T>(def_val);
                    break;

                case 2:
                    if (*(++itr) == "DEFAULT") {
                        // デフォルト値を返す
                        return boost::optional<T>(def_val);
                    }
                    else {
                        try {
                            return boost::optional<T>(boost::lexical_cast<T>(itr->c_str()));
                        }
                        catch (boost::bad_lexical_cast const &) {
                            errMsg(lineindex_ - 1, article, *itr);
                            return boost::none;
                        }
                    }

                default:
                {
                    auto val = *itr;

                    if (val == "DEFAULT" || val[0] == '#') {
                        return boost::optional<T>(def_val);
                    }
                    else if ((*(++itr))[0] != '#') {
                        errMsg(lineindex_ - 1, article, *itr);
                        return boost::none;
                    }

                    try {
                        return boost::optional<T>(boost::lexical_cast<T>(val.c_str()));
                    }
                    catch (boost::bad_lexical_cast const &) {
                        errMsg(lineindex_ - 1, article, val);
                        return boost::none;
                    }
                    break;
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
}

#endif	// _READINPUTFILE_H_
