#include "readinputfile.h"
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include <boost/cast.hpp>
#include <boost/range/algorithm.hpp>

namespace schrac {
    using namespace boost::algorithm;

    const ci_string ReadInputFile::CHEMICAL_SYMBOL = "chemical.symbol";
    const ci_string ReadInputFile::EQ_TYPE_DEFAULT = "sch";
    const ci_string ReadInputFile::EQ_TYPE = "eq.type";
    const std::array<ci_string, 4> ReadInputFile::EQ_TYPE_ARRAY =
    {
        ci_string("sch"),
        ci_string("sdirac"),
        ci_string("dirac")
    };
    const ci_string ReadInputFile::ORBITAL = "orbital";
    const std::array<ci_string, 4> ReadInputFile::SOLVER_TYPE_ARRAY =
    {
        ci_string("adams_bashforth_moulton"),
        ci_string("bulirsch_stoer"),
        ci_string("controlled_runge_kutta")
    };
    const ci_string ReadInputFile::SOLVER_TYPE_DEFAULT = "bulirsch_stoer";
    const ci_string ReadInputFile::SPIN_ORBITAL = "spin.orbital";
    
    ReadInputFile::ReadInputFile(std::pair<std::string, bool> const & arg) :
        ifs_(std::get<0>(arg).c_str()),
        lineindex_(1),
        pdata_(std::make_shared<Data>())
    {
        pdata_->usetbb_ = std::get<1>(arg);
    }

    std::shared_ptr<Data> && ReadInputFile::getpData()
    {
        return std::move(pdata_);
    }

	void ReadInputFile::readFile()
	{
		if (!ifs_.is_open())
			throw std::runtime_error("インプットファイルが開けませんでした");

		if (!readAtom())
			throw std::runtime_error("インプットファイルが異常です");

		if (!readEq())
			throw std::runtime_error("インプットファイルが異常です");

		if (!readGrid())
			throw std::runtime_error("インプットファイルが異常です");

		if (!readEps())
			throw std::runtime_error("インプットファイルが異常です");

		if (!readType())
			throw std::runtime_error("インプットファイルが異常です");

		if (!readLowerE())
			throw std::runtime_error("インプットファイルが異常です");

		if (!readNumofp())
			throw std::runtime_error("インプットファイルが異常です");

		if (!readRatio())
			throw std::runtime_error("インプットファイルが異常です");
	}
    
    void ReadInputFile::errMsg(ci_string const & s) const
    {
        std::cerr << "インプットファイルに" << s << "の行が見つかりませんでした" << std::endl;
    }

    void ReadInputFile::errMsg(std::int32_t line, ci_string const & s1, ci_string const & s2) const
    {
        std::cerr << "インプットファイルの[" << s1 << "]の行が正しくありません" << std::endl;
        std::cerr << line << "行目, 未知のトークン:" << s2 << std::endl;
    }

    std::pair<std::int32_t, boost::optional<ReadInputFile::strvec>> ReadInputFile::getToken(ci_string const & article)
    {
        std::array<char, BUFSIZE> buf;
        ifs_.getline(buf.data(), BUFSIZE);
        ci_string const line(buf.data());

        // もし一文字も読めなかったら
        if (!ifs_.gcount()) {
            errMsg(article);
            return std::make_pair(-1, boost::none);
        }

        // 読み込んだ行が空、あるいはコメント行でないなら
        if (!line.empty() && (line[0] != '#')) {
            // トークン分割
            strvec tokens;
            split(tokens, line, is_any_of(" \t"), token_compress_on);
            
            auto itr(tokens.begin());

            if (*itr != article) {
                errMsg(lineindex_, article, *itr);
                return std::make_pair(-1, boost::none);
            }

            return std::make_pair(0, boost::optional<strvec>(std::move(tokens)));
        }
        else {
            return std::make_pair(1, boost::none);
        }
    }

    boost::optional<ci_string> ReadInputFile::readData(ci_string const & article)
	{
		for (; true; lineindex_++) {
            auto const ret = getToken(article);

            switch (std::get<0>(ret))
            {
            case -1:
                return boost::none;
                break;

            case 0:
            {
                auto const tokens = *(std::get<1>(ret));

                // 読み込んだトークンの数がもし2個以外だったら
                if (tokens.size() != 2 || tokens[1].empty()) {
                    std::cerr << "インプットファイル" << lineindex_ << "行目の、[" << article << "]の行が正しくありません" << std::endl;
                    return boost::none;
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

    boost::optional<ci_string> ReadInputFile::readData(ci_string const & article, ci_string const & def)
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
                ++lineindex_;

                // 読み込んだトークンの数をはかる
                switch (tokens.size()) {
                case 1:
                    // デフォルト値を返す
                    return std::move(def);
                    break;

                case 2:
                    return *itr == "DEFAULT" ? std::move(def) : *itr;
                    break;

                default:
                    {
                        auto const val = *itr;

                        if (val == "DEFAULT" || val[0] == '#') {
                            // デフォルト値を返す
                            return std::move(def);
                        }
                        else if ((*(++itr))[0] != '#') {
                            errMsg(lineindex_ - 1, article, *itr);
                            return nullptr;
                        }

                        return std::move(val);
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

    boost::optional<ci_string> ReadInputFile::readDataAuto(ci_string const & article)
	{
		for (; true; lineindex_++) {
            auto const ret = getToken(article);

            switch (std::get<0>(ret))
            {
            case -1:
                return boost::none;
                break;

            case 0:
            {
                auto const tokens = *(std::get<1>(ret));
                auto itr(++tokens.begin());
                ++lineindex_;

                // 読み込んだトークンの数をはかる
                switch (tokens.size()) {
                case 1:
                    return boost::none;
                    break;

                case 2:
                    return (*itr == "DEFAULT" || *itr == "AUTO") ?
                        boost::optional<ci_string>(ci_string()) : boost::optional<ci_string>(*itr);
                    break;

                default:
                    {
                        auto val = *itr;

                        if (val == "DEFAULT" || val == "AUTO" || val[0] == '#') {
                            // デフォルト値を返す
                            return boost::optional<ci_string>(ci_string());
                        } else if ((*(++itr))[0] != '#') {
                            errMsg(lineindex_ - 1, article, *itr);

                            // エラー
                            return boost::none;
                        }

                        return boost::optional<ci_string>(std::move(val));
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
	
	bool ReadInputFile::readAtom()
	{
		// 原子の種類を読み込む
        auto const chemsym(readData(ReadInputFile::CHEMICAL_SYMBOL));
        if (!chemsym) {
            return false;
        }

		try {
            auto const itr = boost::find(Data::Chemical_Symbol, chemsym->c_str());
            if (itr == Data::Chemical_Symbol.end()) {
                throw std::invalid_argument("");
            }

            pdata_->chemical_symbol_ = *itr;
            pdata_->Z_ = static_cast<double>(std::distance(Data::Chemical_Symbol.begin(), itr)) + 1.0;
		} catch (std::invalid_argument const &) {
            errMsg(lineindex_ - 1, ReadInputFile::CHEMICAL_SYMBOL, *chemsym);
			return false;
		}

		// 軌道を読み込む
        auto const porbital(readData(ReadInputFile::ORBITAL));
		if (!porbital) {
			return false;
		}
        
        auto const orbital(*porbital);
        if (orbital.length() != 2) {
            errMsg(lineindex_ - 1, ReadInputFile::ORBITAL, orbital);
			return false;
		}

		if (!std::isdigit(orbital[0])) {
            errMsg(lineindex_ - 1, ReadInputFile::ORBITAL, orbital);
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
                errMsg(lineindex_ - 1, ReadInputFile::ORBITAL, orbital);
				return false;
			break;
		}

		if (pdata_->n_ - pdata_->l_ < 1) {
			std::cerr << "量子数の指定が異常です" << std::endl;
			return false;
		}

		// スピン軌道を読み込む
        auto const pspin_orbital(readData(ReadInputFile::SPIN_ORBITAL));
		if (!pspin_orbital) {
			return false;
		}
        
        pdata_->spin_orbital_ = *pspin_orbital;
        if (pdata_->spin_orbital_ != Data::ALPHA && pdata_->spin_orbital_ != Data::BETA) {
            errMsg(lineindex_ - 1, ReadInputFile::SPIN_ORBITAL, pdata_->spin_orbital_);
			return false;
		}

		if (!pdata_->l_) {
            // s軌道は特別なケース
			// j = l_ + 1/2
			pdata_->j_ = 0.5;
			pdata_->kappa_ = - 1.0;
		} else if (pdata_->spin_orbital_ == Data::ALPHA) {
            // j = l_ + 1/2に対して
			pdata_->j_ = static_cast<long double>(pdata_->l_) + 0.5;
			pdata_->kappa_ = - static_cast<long double>(pdata_->l_) - 1.0;
		} else {
            // j = l_ - 1/2に対して
			pdata_->j_ = static_cast<long double>(pdata_->l_) - 0.5;
			pdata_->kappa_ = static_cast<long double>(pdata_->l_);
		}

		return true;
	}
	
    bool ReadInputFile::readEps()
    {
        if (auto const peps = readData<long double>("eps", pdata_->EPS_DEFAULT)) {
            pdata_->eps_ = *peps;
        }
        else {
            return false;
        }

        return true;
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
            errMsg(lineindex_ - 1, ReadInputFile::EQ_TYPE, eqtype);
			return false;
        }
        else if (itr != ReadInputFile::EQ_TYPE_ARRAY.begin()) {
            pdata_->eq_type_ = boost::numeric_cast<const Data::Eq_type>(std::distance(ReadInputFile::EQ_TYPE_ARRAY.begin(), itr));
		}

		return true;
	}

    bool ReadInputFile::readGrid()
    {
        if (auto const pxmin = readData<double>("grid.xmin", pdata_->XMIN_DEFAULT)) {
            pdata_->xmin_ = *pxmin;
        }
        else {
            return false;
        }

        if (auto const pxmax = readData<double>("grid.xmax", pdata_->XMAX_DEFAULT)) {
            pdata_->xmax_ = *pxmax;
        }
        else {
            return false;
        }

        if (auto const pgrid_num = readData<std::size_t>("grid.num", Data::GRID_NUM_DEFAULT)) {
			pdata_->grid_num_ = *pgrid_num;
		} else {
			return false;
		}

        return true;
	}

    bool ReadInputFile::readLowerE()
    {
        if (auto const val = readDataAuto("search.LowerE")) {
            if (!val->empty()) {
                try {
                    std::size_t idx;
                    auto const v = std::stod(val->c_str(), &idx);
                    if (idx != val->length()) {
                        throw std::invalid_argument("");
                    }
                    pdata_->search_lowerE_ = boost::optional<double>(v);
                }
                catch (std::invalid_argument const &) {
                    errMsg(lineindex_ - 1, "search.LowerE", *val);
                    return false;
                }
            }
            else {
                pdata_->search_lowerE_ = boost::none;
            }
        }
        else {
            return false;
        }

        return true;
    }

    bool ReadInputFile::readNumofp()
    {
        if (auto const pnumofp = readData<std::size_t>("num.of.partition", Data::NUM_OF_PARTITION_DEFAULT)) {
            pdata_->num_of_partiton_ = *pnumofp;
        }
        else {
            return false;
        }

        return true;
    }

    bool ReadInputFile::readRatio()
    {
        if (auto const pmpr = readData<double>("matching.point.ratio", pdata_->MAT_PO_RATIO_DEFAULT)) {
            pdata_->mat_po_ratio_ = *pmpr;
        }
        else {
            return false;
        }

        return true;
    }

	bool ReadInputFile::readType()
	{
        auto const psolvetype(readData("solver.type", ReadInputFile::SOLVER_TYPE_DEFAULT));
        if (!psolvetype) {
            return false;
        }

        auto const solvertype = *psolvetype;
		auto const itr = boost::find(ReadInputFile::SOLVER_TYPE_ARRAY, solvertype);
        if (itr == ReadInputFile::SOLVER_TYPE_ARRAY.end()) {
			errMsg(lineindex_ - 1, "solver.type", solvertype);
			return false;
		} else {
            pdata_->solver_type_ = boost::numeric_cast<Data::Solver_type>(std::distance(ReadInputFile::SOLVER_TYPE_ARRAY.begin(), itr));
		}

		return true;
	}    
}
