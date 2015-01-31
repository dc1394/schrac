#ifndef _EIGENVALUESEARCH_H_
#define _EIGENVALUESEARCH_H_

#include "Diff.h"
#include "ReadInputFile.h"

namespace schrac {
	class EigenValueSearch final {
        // #region コンストラクタ・デストラクタ

    public:

        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param arg インプットファイル名とTBBを使用するかどうかのstd::pair
        */
        explicit EigenValueSearch(std::pair<std::string, bool> const & arg);

        //! A destructor.
        /*!
            何もしないデストラクタ
        */
        ~EigenValueSearch()
        {
        }

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

    public:
        //! A public member function.
        /*!
            データオブジェクトを得る
            \return データオブジェクト
        */
        const std::shared_ptr<Data> & getpData() const;
        
        //! A public member function.
        /*!
            微分方程式オブジェクトを得る
            \return 微分方程式オブジェクト
        */
        const std::shared_ptr<Diff> & getpDiff() const;

        //! A public member function.
        /*!
            固有値を検索する
            \return 固有値が見つかったかどうか
        */
        bool search();

    private:
        //! A private member function.
        /*!
            関数Dの値を返す
            \return 関数Dの値（微分方程式が正常に解けなかったときはboost::none）
        */
        boost::optional<double> fnc_D();

        //! A private member function (const).
        /*!
            現在のループをメッセージで報告する
        */
        void info() const;

        //! A private member function (const).
        /*!
            固有値が見つかったことをメッセージで報告する
        */
        void info(double E) const;

        //! A private member function (const).
        /*!
            固有値が見つかったことをメッセージで報告する
            \param b 関数Dの引数b
            \param fb 関数Dの引数fb
        */
        void info(double b, double fb) const;

        //! A private member function (const).
        /*!
            状態の初期化を行う
        */
        void init();

        //! A private member function (const).
        /*!
            解く微分方程式についてメッセージを表示する
        */
        void msg() const;

        //! A private member function.
        /*!
            固有値をおおざっぱに検索する
            \return 固有値が見つかったかどうか
        */
        bool rough_search();

        //! A private member function (const).
        /*!
            表示する浮動小数点の桁を設定する
        */
        void setoutstream() const;

        // #endregion メンバ関数

        // #region メンバ変数

        //! A private member variable (constant expression).
        /*!
            エネルギー固有値探索の最大のループ回数
        */
		static constexpr auto EVALSEARCHMAX = 1000;

        //! A private member variable (constant expression).
        /*!
            閾値（絶対値の大きい方）
        */
        static constexpr auto HUGE = 1.0E+7;

        //! A private member variable (constant expression).
        /*!
            閾値（絶対値の小さい方）
        */
		static constexpr auto TINY = 1.0E-30;
        
        //! A private member variable (constant).
        /*!
            許容誤差
        */
		double eps_;

        //! A private member variable (constant).
        /*!
            許容誤差（eps * 10.0）
        */
        double tol_;

        //! A private member variable.
        /*!
            Brent法における関数Dの大きい方
        */
        double Dmax;

        //! A private member variable.
        /*!
            Brent法における関数Dの小さい方
        */
        double Dmin;

        //! A private member variable.
        /*!
            前のループと今のエネルギー固有値の差
        */
        double DE;

        //! A private member variable.
        /*!
            エネルギー固有値
        */
        double E;

        //! A private member variable.
        /*!
            Brent法におけるエネルギー固有値の大きい方
        */
        double Emax;

        //! A private member variable.
        /*!
            Brent法におけるエネルギー固有値の小さい方
        */
        double Emin;

        //! A private member variable.
        /*!
            エネルギー固有値の大体の値
        */
        double Erough_exact_;

        //! A private member variable.
        /*!
            関数Dの古い値
        */
        double Dold;
		
        //! A private member variable.
        /*!
            エネルギー固有値探索のループ回数
        */
        std::int32_t loop_;

        //! A private member variable.
        /*!
            固有関数のノードが一致しているかどうか
        */
        bool noden_;

        //! A private member variable.
        /*!
            インプットファイルのデータオブジェクト
        */
		std::shared_ptr<Data> pdata_;

        //! A private member variable.
        /*!
            微分方程式オブジェクト
        */
		std::shared_ptr<Diff> pdiff_;

        //! A private member variable.
        /*!
            微分方程式データのオブジェクト
        */
		std::shared_ptr<DiffData> pdiffdata_;
        
        // #endregion メンバ関数
		
		bool zbrent();

    private:
        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
        デフォルトコンストラクタ（禁止）
        */
        EigenValueSearch() = delete;

        //! A private copy constructor (deleted).
        /*!
        コピーコンストラクタ（禁止）
        */
        EigenValueSearch(EigenValueSearch const &) = delete;

        //! A private member function (deleted).
        /*!
        operator=()の宣言（禁止）
        \param コピー元のオブジェクト（未使用）
        \return コピー元のオブジェクト
        */
        EigenValueSearch & operator=(EigenValueSearch const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
	};

    // #region 非メンバ関数

    //! A function.
    /*!
        対象の方程式がDirac方程式の場合に、エネルギー固有値の下限を概算する
        \param pdata データオブジェクト
        \return 大体のエネルギー固有値
    */
    double Eexact_dirac(std::shared_ptr<Data> const & pdata);
    
    //! A function.
    /*!
        対象の方程式がSch方程式の場合に、エネルギー固有値の下限を概算する
        \param pdata データオブジェクト
        \return 大体のエネルギー固有値
    */
    double Eexact_sch(std::shared_ptr<Data> const & pdata);

    //! A function.
    /*!
        対象の方程式がscalar Dirac方程式の場合に、エネルギー固有値の下限を概算する
        \param pdata データオブジェクト
        \return 大体のエネルギー固有値
    */
    double Eexact_sdirac(std::shared_ptr<Data> const & pdata);

	template <typename T>
    //! A function (template function).
    /*!
        bが正の値の場合にaの絶対値を、bが負の値の場合はaの絶対値に-をかけた値を返す
        \param a 対象の値
        \param b 正負を判断するための値
        \return bが正の値の場合はaの絶対値、bが負の値の場合はaの絶対値に-をかけた値
    */
	T sign(T a, T b)
	{
		return (b >= 0.0) ? std::fabs(a) : - std::fabs(a);
	}

    // #endregion 非メンバ関数
}

#endif // _EIGENVALUESEARCH_H_