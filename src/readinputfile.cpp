/*! \file readinputfile.cpp
    \brief インプットファイルの読み込みを行うクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "readinputfile.h"
#include <iostream>                     // for std::cerr
#include <stdexcept>                    // for std::runtime_error
#include <boost/algorithm/string.hpp>   // for boost::algorithm
#include <boost/assert.hpp>             // for BPOOST_ASSERT
#include <boost/cast.hpp>               // for boost::numeric_cast
#include <boost/range/algorithm.hpp>    // for boost::find

namespace schrac {
    // #region staticメンバ変数

    ci_string const ReadInputFile::CHEMICAL_SYMBOL = "chemical.symbol";
	ci_string const ReadInputFile::EQ_TYPE_DEFAULT = "sch";
	ci_string const ReadInputFile::EQ_TYPE = "eq.type";
    std::array<ci_string, 4> const ReadInputFile::EQ_TYPE_ARRAY =
    {
        ci_string("sch"),
        ci_string("sdirac"),
        ci_string("dirac")
    };
    ci_string const ReadInputFile::ORBITAL = "orbital";
    const std::array<ci_string, 4> ReadInputFile::SOLVER_TYPE_ARRAY =
    {
        ci_string("adams_bashforth_moulton"),
        ci_string("bulirsch_stoer"),
        ci_string("controlled_runge_kutta")
    };
    ci_string const ReadInputFile::SOLVER_TYPE_DEFAULT = "controlled_runge_kutta";
	ci_string const ReadInputFile::SPIN_ORBITAL = "spin.orbital";
	ci_string const ReadInputFile::SPIN_ORBITAL_DEFAULT = "alpha";

    // #endregion staticメンバ変数

    // #region コンストラクタ

    ReadInputFile::ReadInputFile(std::pair<std::string, bool> const & arg) :
        PData([this]() { return std::cref(pdata_); }, nullptr), 
        ifs_(std::get<0>(arg).c_str()),
        lineindex_(1),
        pdata_(std::make_shared<Data>())
    {
        pdata_->usetbb_ = std::get<1>(arg);
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数
        
    void ReadInputFile::readFile()
    {
        if (!ifs_.is_open())
            throw std::runtime_error("インプットファイルが開けませんでした");

        auto const errorendfunc = []() {
            throw std::runtime_error("インプットファイルが異常です");
        };

        if (!readAtom()) {
            errorendfunc();
        }

        if (!readEq()) {
            errorendfunc();
        }

        // グリッドの最小値を読み込む
        readValue("grid.xmin", XMIN_DEFAULT, pdata_->xmin_);

        // グリッドの最大値を読み込む
        readValue("grid.xmax", XMAX_DEFAULT, pdata_->xmax_);

        // グリッドのサイズを読み込む
        readValue("grid.num", GRID_NUM_DEFAULT, pdata_->grid_num_);

        // 許容誤差を読み込む
        readValue("eps", EPS_DEFAULT, pdata_->eps_);

        // 解く方程式のタイプを読み込む
        if (!readSolverType()) {
            errorendfunc();
        }

        // 固有値探索をはじめる値を読み込む
        if (!readValueAuto("search.LowerE", pdata_->search_lowerE_)) {
            errorendfunc();
        }

        // 固有値検索の間隔を読み込む
        readValue("num.of.partition", NUM_OF_PARTITION_DEFAULT, pdata_->num_of_partition_);

        // マッチングポイントを読み込む
        readValue("matching.point.ratio", MAT_PO_RATIO_DEFAULT, pdata_->mat_po_ratio_);

        // 密度の初期値ρ0(r)のための係数cを読み込む
        if (!readValueAuto("rho0.c", pdata_->rho0_c_)) {
            errorendfunc();
        }

        // 密度の初期値ρ0(r)のための係数alphaを読み込む
        if (!readValueAuto("rho0.alpha", pdata_->rho0_alpha_)) {
            errorendfunc();
        }

        // SCFの最大ループ回数を読み込む
        readValue("scf.maxIter", SCF_MAXITER_DEFAULT, pdata_->scf_maxiter_);

        // SCFの一次混合の重みを読み込む
        if (!readScfMixingWeight()) {
            errorendfunc();
        }
        
        // SCFの収束判定条件の値を読み込む
        readValue("scf.criterion", SCF_CRITERION_DEFAULT, pdata_->scf_criterion_);
    }
    
    // #endregion publicメンバ関数

    // #region privateメンバ関数

    void ReadInputFile::errorMessage(ci_string const & s) const
    {
        std::cerr << "インプットファイルに" << s.c_str() << "の行が見つかりませんでした。\n";
    }

    void ReadInputFile::errorMessage(std::int32_t line, ci_string const & s1, ci_string const & s2) const
    {
        std::cerr << "インプットファイルの[" << s1.c_str() << "]の行が正しくありません。\n";
        std::cerr << line << "行目, 未知のトークン:" << s2.c_str() << std::endl;
    }

    std::pair<std::int32_t, std::optional<ReadInputFile::strvec>> ReadInputFile::getToken(ci_string const & article)
    {
        using namespace boost::algorithm;

        std::array<char, BUFSIZE> buf;
        ifs_.getline(buf.data(), BUFSIZE);
        ci_string const line(buf.data());

        // もし一文字も読めなかったら
        if (!ifs_.gcount()) {
            errorMessage(article);
            return std::make_pair(-1, std::nullopt);
        }

        // 読み込んだ行が空、あるいはコメント行でないなら
        if (!line.empty() && (line[0] != '#')) {
            // トークン分割
            strvec tokens;
            split(tokens, line, is_any_of(" \t"), token_compress_on);
            
            auto const itr(tokens.begin());

            if (*itr != article) {
                errorMessage(lineindex_, article, *itr);
                return std::make_pair(-1, std::nullopt);
            }

            return std::make_pair(0, std::make_optional<strvec>(std::move(tokens)));
        }
        else {
            return std::make_pair(1, std::nullopt);
        }
    }

    bool ReadInputFile::readAtom()
    {
        // 原子の種類を読み込む
        auto const chemsym(readData(ReadInputFile::CHEMICAL_SYMBOL));
        if (!chemsym) {
            return false;
        }

        try {
            auto const itr(boost::find(Data::Chemical_Symbol, chemsym->c_str()));
            if (itr == Data::Chemical_Symbol.end()) {
                throw std::invalid_argument("");
            }

            pdata_->chemical_symbol_ = *itr;
            pdata_->Z_ = static_cast<double>(std::distance(Data::Chemical_Symbol.begin(), itr)) + 1.0;
        }
        catch (std::invalid_argument const &) {
            errorMessage(lineindex_ - 1, ReadInputFile::CHEMICAL_SYMBOL, *chemsym);
            return false;
        }

        // 軌道を読み込む
        auto const porbital(readData(ReadInputFile::ORBITAL));
        if (!porbital) {
            return false;
        }

        auto const orbital(*porbital);
        if (orbital.length() != 2) {
            errorMessage(lineindex_ - 1, ReadInputFile::ORBITAL, orbital);
            return false;
        }

        if (!std::isdigit(orbital[0])) {
            errorMessage(lineindex_ - 1, ReadInputFile::ORBITAL, orbital);
            return false;
        }
        pdata_->orbital_ = orbital[0];
        pdata_->n_ = boost::numeric_cast<std::uint8_t>(orbital[0] - '0');

        switch (orbital[1]) {
        case 's':
            pdata_->l_ = 0;
            pdata_->orbital_ += 's';
            break;

        case 'p':
            pdata_->l_ = 1;
            pdata_->orbital_ += 'p';
            break;

        case 'd':
            pdata_->l_ = 2;
            pdata_->orbital_ += 'd';
            break;

        case 'f':
            pdata_->l_ = 3;
            pdata_->orbital_ += 'f';
            break;

        case 'g':
            pdata_->l_ = 4;
            pdata_->orbital_ += 'g';
            break;

        default:
            errorMessage(lineindex_ - 1, ReadInputFile::ORBITAL, orbital);
            return false;
            break;
        }

        if (pdata_->n_ - pdata_->l_ < 1) {
            std::cerr << "量子数の指定が異常です。\n";
            return false;
        }

        // スピン軌道を読み込む
        auto const pspin_orbital(readData(ReadInputFile::SPIN_ORBITAL, ReadInputFile::SPIN_ORBITAL_DEFAULT));
        if (!pspin_orbital) {
            return false;
        }

        pdata_->spin_orbital_ = *pspin_orbital;
        if (pdata_->spin_orbital_ != Data::ALPHA && pdata_->spin_orbital_ != Data::BETA) {
            errorMessage(lineindex_ - 1, ReadInputFile::SPIN_ORBITAL, pdata_->spin_orbital_);
            return false;
        }

        if (!pdata_->l_) {
            // s軌道は特別なケース
            // j = l_ + 1/2
            pdata_->j_ = 0.5;
            pdata_->kappa_ = -1.0;
        }
        else if (pdata_->spin_orbital_ == Data::ALPHA) {
            // j = l_ + 1/2に対して
            pdata_->j_ = static_cast<double>(pdata_->l_) + 0.5;
            pdata_->kappa_ = -static_cast<double>(pdata_->l_) - 1.0;
        }
        else {
            // j = l_ - 1/2に対して
            pdata_->j_ = static_cast<double>(pdata_->l_) - 0.5;
            pdata_->kappa_ = static_cast<double>(pdata_->l_);
        }

        return true;
    }

    std::optional<ci_string> ReadInputFile::readData(ci_string const & article)
    {
        for (; true; lineindex_++) {
            auto const ret(getToken(article));

            switch (std::get<0>(ret))
            {
            case -1:
                return std::nullopt;
                break;

            case 0:
            {
                auto const tokens(*(std::get<1>(ret)));

                // 読み込んだトークンの数がもし2個以外だったら
                if (tokens.size() != 2 || tokens[1].empty()) {
                    std::cerr << "インプットファイル" << lineindex_ << "行目の、[" << article.c_str() << "]の行が正しくありません。\n";
                    return std::nullopt;
                }

                ++lineindex_;
                return *(++tokens.begin());
            }
            case 1:
                break;

            default:
                BOOST_ASSERT(!"何かがおかしい!");
                break;
            }
        }
    }

    std::optional<ci_string> ReadInputFile::readData(ci_string const & article, ci_string const & def)
    {
        // グリッドを読み込む
        for (; true; lineindex_++) {
            auto const ret(getToken(article));

            switch (std::get<0>(ret))
            {
            case -1:
                return std::nullopt;
                break;

            case 0:
            {
                auto const tokens(*(std::get<1>(ret)));
                auto itr(++tokens.begin());
                ++lineindex_;

                // 読み込んだトークンの数をはかる
                switch (tokens.size()) {
                case 1:
                    // デフォルト値を返す
                    return def;
                    break;

                case 2:
                    return *itr == "DEFAULT" ? def : *itr;
                    break;

                default:
                    {
                        auto val(*itr);

                        if (val == "DEFAULT" || val[0] == '#') {
                            // デフォルト値を返す
                            return def;
                        }
                        else if ((*(++itr))[0] != '#') {
                            errorMessage(lineindex_ - 1, article, *itr);
                            return std::nullopt;
                        }

                        return val;
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

    std::optional<ci_string> ReadInputFile::readDataAuto(ci_string const & article)
    {
        for (; true; lineindex_++) {
            auto const ret(getToken(article));

            switch (std::get<0>(ret))
            {
            case -1:
                return std::nullopt;
                break;

            case 0:
            {
                auto const tokens(*(std::get<1>(ret)));
                auto itr(++tokens.begin());
                ++lineindex_;

                // 読み込んだトークンの数をはかる
                switch (tokens.size()) {
                case 1:
                    return std::nullopt;
                    break;

                case 2:
                    return (*itr == "DEFAULT" || *itr == "AUTO") ?
                        std::make_optional<ci_string>(ci_string()) : std::make_optional<ci_string>(*itr);
                    break;

                default:
                    {
                        auto val = *itr;

                        if (val == "DEFAULT" || val == "AUTO" || val[0] == '#') {
                            // デフォルト値を返す
                            return std::make_optional<ci_string>(ci_string());
                        } else if ((*(++itr))[0] != '#') {
                            errorMessage(lineindex_ - 1, article, *itr);

                            // エラー
                            return std::nullopt;
                        }

                        return std::make_optional<ci_string>(std::move(val));
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
    
    bool ReadInputFile::readEq()
    {
        auto const peqtype(readData(ReadInputFile::EQ_TYPE, ReadInputFile::EQ_TYPE_DEFAULT));
        if (!peqtype) {
            return false;
        }

        auto const eqtype(*peqtype);
        auto const itr(boost::find(ReadInputFile::EQ_TYPE_ARRAY, eqtype));
        
        if (itr == ReadInputFile::EQ_TYPE_ARRAY.end()) {
            errorMessage(lineindex_ - 1, ReadInputFile::EQ_TYPE, eqtype);
            return false;
        }
        else if (itr != ReadInputFile::EQ_TYPE_ARRAY.begin()) {
            pdata_->eq_type_ = boost::numeric_cast<Data::Eq_type>(std::distance(ReadInputFile::EQ_TYPE_ARRAY.begin(), itr));
        }

        if (pdata_->eq_type_ == Data::Eq_type::SDIRAC) {
            // scalar Dirac方程式の場合
            // j = l_ + 1/2
            pdata_->j_ = 0.5;
            pdata_->kappa_ = -1.0;
        }

        return true;
    }

    bool ReadInputFile::readScfMixingWeight()
    {
        readValue("scf.Mixing.Weight", SCF_MIXING_WEIGHT_DEFAULT, pdata_->scf_mixing_weight_);
        if (pdata_->scf_mixing_weight_ <= 0.0 || pdata_->scf_mixing_weight_ > 1.0) {
            std::cerr << "インプットファイルの[scf.Mixing.Weight]の行が正しくありません。\n";
            return false;
        }
        return true;
    }

    bool ReadInputFile::readSolverType()
    {
        auto const psolvetype(readData("solver.type", ReadInputFile::SOLVER_TYPE_DEFAULT));
        if (!psolvetype) {
            return false;
        }

        auto const solvertype = *psolvetype;
        auto const itr(boost::find(ReadInputFile::SOLVER_TYPE_ARRAY, solvertype));
        if (itr == ReadInputFile::SOLVER_TYPE_ARRAY.end()) {
            errorMessage(lineindex_ - 1, "solver.type", solvertype);
            return false;
        } else {
            pdata_->solver_type_ = boost::numeric_cast<Data::Solver_type>(
                std::distance(ReadInputFile::SOLVER_TYPE_ARRAY.begin(), itr));
        }

        return true;
    }

    // #endregion privateメンバ関数
}
