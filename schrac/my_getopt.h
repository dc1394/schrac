#ifndef _MY_GETOPT_H_
#define _MY_GETOPT_H_

#include <cstdint>
#include <string>
#include <utility>

namespace schrac {
    //! A class.
    /*!
        コマンドラインオプションを解析するクラス    
    */
	class my_getopt final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            何もしないコンストラクタ
        */
        my_getopt()
        {
        }

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~my_getopt()
        {
        }
        
        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //! A public member function.
        /*!
            コマンドラインオプションを解析する関数の宣言
            \param argc コマンドライン引数の数
            \param argv コマンドライン引数
            \return 解析に成功したら0、失敗したら-1
        */
        std::int32_t getopt(int argc, char * const argv[]);

        //! A public member function (constant).
        /*!
            \return インプットファイルと並列計算するときのスレッド数のstd::pair
        */
        std::pair<std::string, bool> getpairdata() const;

        // #endregion メンバ関数

    private:
        // #region メンバ変数

        //!  A public static member variable (constant).
        /*!
            デフォルトのインプットファイル名
        */
        static std::string const DEFINPNAME;
        
        //!  A public member variable.
        /*!
            デフォルトのインプットファイル名
        */
        std::string inpname_;

        //!  A public member variable.
        /*!
            TBBを使用するかどうか    
        */
        bool usetbb_ = false;

        // #endregion メンバ変数
	};
}

#endif	// _MY_GETOPT_H_

