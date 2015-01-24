#include "ReadInpFile.h"

namespace HydroSchDirac {
	void ReadInpFile::readFile()
	{
		if (!ifs.is_open())
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

#if (_MSC_VER >= 1600)
	ci_string ReadInpFile::readData(const char * const article)
#else
	const ci_string ReadInpFile::readData(const char * const article)
#endif
	{
		for (; true; i_++) {
			char buf[BUFSIZE];
			ifs.getline(buf, BUFSIZE);
			const ci_string line(buf);
 
			// もし一文字も読めなかったら
			if (!ifs.gcount()) {
				errMsg(article);
				return ci_string();
			}

			// 読み込んだ行が空、あるいはコメント行でないなら
			if (!line.empty() && (line[0] != '#')) {
				// トークン分割
				strvec tokens;
				split(tokens, line, is_any_of(" \t"), token_compress_on);
				
				// 読み込んだトークンの数がもし2個以外だったら
				if (tokens.size() != 2) {
					std::cerr << "インプットファイルの" << article << "の行が正しくありません" << std::endl;
					return ci_string();
				}

				strvec::const_iterator citr(tokens.begin());
				if (*citr == article) {
					++citr;

					i_++;
#if (_MSC_VER >= 1600)
					return std::move(*citr);
#else
					return *citr;
#endif
				} else {
					errMsg(article, citr->c_str());					
					return ci_string();
				}
			}
		}
	}

#if (_MSC_VER >= 1600)
	ci_string ReadInpFile::readData(const char * const article, ci_string && def)
#else
	const ci_string ReadInpFile::readData(const char * const article, const ci_string & def)
#endif
	{
		// グリッドを読み込む
		for (; true; i_++) {
			char buf[BUFSIZE];
			ifs.getline(buf, BUFSIZE);
			const ci_string line(buf);
 
			// もし一文字も読めなかったら
			if (!ifs.gcount()) {
				errMsg(article);
				return ci_string();
			}

			// 読み込んだ行が空、あるいはコメント行でないなら
			if (!line.empty() && (line[0] != '#')) {
				// トークン分割
				strvec tokens;
				split(tokens, line, is_any_of(" \t"), token_compress_on);

				strvec::const_iterator citr(tokens.begin());

				if (*citr != article) {
					errMsg(article, *citr);
					return ci_string();
				}

				// 読み込んだトークンの数をはかる
				const strvec::size_type size = tokens.size();

				ci_string val;
				switch (size) {
					case 1:
#if (_MSC_VER >= 1600)
						return std::move(def);				// デフォルト
#else
						return def;
#endif
					break;

					case 2:
						++citr;
						if (*citr == "DEFAULT") {
							i_++;
#if (_MSC_VER >= 1600)
							return std::move(def);				// デフォルト
#else
							return def;
#endif
						} else {
							i_++;
#if (_MSC_VER >= 1600)
							return std::move(*citr);
#else
							return *citr;
#endif
						}
					break;

					default:
						++citr;
						val = *citr;
						if (val == "DEFAULT" || val[0] == '#') {
							i_++;
#if (_MSC_VER >= 1600)
							return std::move(def);				// デフォルト
#else
							return def;
#endif
						}
						
						++citr;

						if ((*citr)[0] != '#') {
							errMsg(article, *citr);
							return ci_string();
						}
					
						i_++;
#if (_MSC_VER >= 1600)
						return std::move(val);
#else
						return val;
#endif
					break;
				}
			}
		}
	}

	const boost::optional<const ci_string> ReadInpFile::readDataAuto(const char * const article)
	{
		for (; true; i_++) {
			char buf[BUFSIZE];
			ifs.getline(buf, BUFSIZE);
			const ci_string line(buf);
 
			// もし一文字も読めなかったら
			if (!ifs.gcount()) {
				errMsg(article);
				return boost::none;
			}

			// 読み込んだ行が空、あるいはコメント行でないなら
			if (!line.empty() && (line[0] != '#')) {
				// トークン分割
				strvec tokens;
				split(tokens, line, is_any_of(" \t"), token_compress_on);
				strvec::const_iterator citr(tokens.begin());

				if (*citr != article) {
					errMsg(article, *citr);
					return boost::none;
				}

				// 読み込んだトークンの数をはかる
				const strvec::size_type size = tokens.size();

				ci_string val;
				switch (size) {
					case 1:
						return boost::optional<const ci_string>(ci_string());
					break;
				
					case 2:
						++citr;
						if (*citr == "DEFAULT" || *citr == "AUTO") {
							return boost::optional<const ci_string>(ci_string());			// デフォルト
						} else {
							return boost::optional<const ci_string>(*citr);
						}
					break;
				
					default:
						++citr;
						val = *citr;
						if (val == "DEFAULT" || val == "AUTO" || val[0] == '#') {
							return boost::optional<const ci_string>(ci_string());			// デフォルト
						}

						++citr;

						if ((*citr)[0] != '#') {
							errMsg(article, *citr);
							return boost::none;
						} 

						return boost::optional<const ci_string>(val);
					break;
				}
			}
		}
	}

	template <typename T>
	const boost::optional<const T> ReadInpFile::readData(const char * const article, T def_val)
	{
		// グリッドを読み込む
		for (; true; i_++) {
			char buf[BUFSIZE];
			ifs.getline(buf, BUFSIZE);
			const ci_string line(buf);
 
			// もし一文字も読めなかったら
			if (!ifs.gcount()) {
				errMsg(article);
				return boost::none;
			}

			// 読み込んだ行が空、あるいはコメント行でないなら
			if (!line.empty() && (line[0] != '#')) {
				// トークン分割
				strvec tokens;
				split(tokens, line, is_any_of(" \t"), token_compress_on);

				strvec::const_iterator citr(tokens.begin());

				if (*citr != article) {
					errMsg(article, *citr);
					return boost::none;
				}

				// 読み込んだトークンの数をはかる
				const strvec::size_type size = tokens.size();

				if (size == 1) {
					return boost::optional<const T>(def_val);					// デフォルト
				} else if (size == 2) {
					++citr;
					if (*citr == "DEFAULT") {
						i_++;
						return boost::optional<const T>(def_val);				// デフォルト
					} else {
						try {
							i_++;
							return boost::optional<const T>(boost::lexical_cast<T>(citr->c_str()));
						} catch (const boost::bad_lexical_cast &) {
							errMsg(article, *citr);
							return boost::none;
						}
					}
				} else {
					++citr;
					const ci_string val(*citr);
					if (val == "DEFAULT" || val[0] == '#') {
						i_++;
						return boost::optional<const T>(def_val);
					}
					
					++citr;

					if ((*citr)[0] != '#') {
						errMsg(article, *citr);
						return boost::none;
					}

					try {
						i_++;
						return boost::optional<const T>(boost::lexical_cast<T>(val.c_str()));
					} catch (const boost::bad_lexical_cast &) {
						errMsg(article, val);
						return boost::none;
					}
				}
			}
		}
	}
	
	bool ReadInpFile::readAtom()
	{
		// 原子の種類を読み込む
		const ci_string atomnum(readData("atomic.number"));
		if (atomnum.empty())
			return false;

		try {
			const int tmp = boost::lexical_cast<int>(atomnum.c_str());
			if (tmp <= 0)
				throw boost::bad_lexical_cast();
			
			pdata_->Z = boost::numeric_cast<const unsigned int>(tmp);
			if (pdata_->Z > Data::atomName.size())
				throw boost::bad_lexical_cast();
		} catch (const boost::bad_lexical_cast &) {
			errMsg("atomic.number", atomnum);
			return false;
		}

		const unsigned int charge = pdata_->Z - 1;			
		pdata_->atom += Data::atomName[charge];
		
		switch (charge) {
			case 0:
			break;

			case 1:
				pdata_->atom += '+';
			break;

			default:
				pdata_->atom += boost::lexical_cast<std::string>(charge) + '+';
			break;
		}
		
		// 軌道を読み込む
		const ci_string orbtmp(readData("orbital"));
		if (orbtmp.empty()) {
			return false;
		} else if (orbtmp.length() != 2) {
			errMsg("orbital", orbtmp);
			return false;
		}

		if (!std::isdigit(orbtmp[0])) {
			errMsg("orbital", orbtmp);
			return false;
		}
		pdata_->orbital = orbtmp[0];
		pdata_->n = boost::numeric_cast<const unsigned int>(orbtmp[0] - '0');
					
		switch (orbtmp[1]) {
			case 's':
				pdata_->l = 0;
				pdata_->orbital += 's';
			break;

			case 'p':
				pdata_->l = 1;
				pdata_->orbital += 'p';
			break;

			case 'd':
				pdata_->l = 2;
				pdata_->orbital += 'd';
			break;

			case 'f':
				pdata_->l = 3;
				pdata_->orbital += 'f';
			break;

			case 'g':
				pdata_->l = 4;
				pdata_->orbital += 'g';
			break;

			default:
				errMsg("orbital", orbtmp);
				return false;
			break;
		}

		if (pdata_->n - pdata_->l < 1) {
			std::cerr << "量子数の指定が異常です" << std::endl;
			return false;
		}

		// スピン軌道を読み込む
		pdata_->spin_orbital = readData("spin.orbital");
		if (pdata_->spin_orbital.empty()) {
			return false;
		} else if (pdata_->spin_orbital != Data::ALPHA && pdata_->spin_orbital != Data::BETA) {
			errMsg("spin.orbital", pdata_->spin_orbital);
			return false;
		}

		if (!pdata_->l) {									// s軌道は特別なケース
			// j = l + 1/2
			pdata_->j_ = 0.5;
			pdata_->kappa = - 1.0;
		} else if (pdata_->spin_orbital == Data::ALPHA) {	// j = l + 1/2に対して
			pdata_->j_ = static_cast<const long double>(pdata_->l) + 0.5;
			pdata_->kappa = - static_cast<const long double>(pdata_->l) - 1.0;
		} else {											// j = l - 1/2に対して
			pdata_->j_ = static_cast<const long double>(pdata_->l) - 0.5;
			pdata_->kappa = static_cast<const long double>(pdata_->l);
		}

		return true;
	}
	
	bool ReadInpFile::readEq()
	{
		const ci_string eqtype(readData("eq.type", ci_string("sdirac")));

		if (eqtype.empty())
			return false;

		const array<const ci_string, 4> eqtypeary =
			{ ci_string("Default"), ci_string("sch"),
			  ci_string("sdirac"), ci_string("dirac") };

		int i;
		array<const ci_string, 4>::const_iterator citr(eqtypeary.begin());
		const array<const ci_string, 4>::const_iterator citr_end(eqtypeary.end());
		for (i = 0; citr != citr_end; ++citr, i++) {
			if (*citr == eqtype)
				break;
		}

		if (citr == eqtypeary.end()) {
			errMsg("eq.type", eqtype);
			return false;
		} else if (!i) {
			pdata_->eqtype = Data::def_eq_type;
		} else {
			pdata_->eqtype = boost::numeric_cast<const Data::eq_type>(i - 1);
		}

		return true;

	}

	bool ReadInpFile::readGrid()
	{
		const boost::optional<const long double &> pxmin(
			readData<long double>("grid.xmin", pdata_->XMIN_DEFAULT));
		if (pxmin)
			pdata_->xmin = *pxmin;
		else
			return false;

		const boost::optional<const long double &> pxmax(
			readData<long double>("grid.xmax", pdata_->XMAX_DEFAULT));
		if (pxmax)
			pdata_->xmax = *pxmax;
		else
			return false;

		const boost::optional<const std::size_t &> p(
			readData<std::size_t>("grid.num", Data::GRID_NUM_DEFAULT));
		if (p) {
			pdata_->grid_num = *p;

			return true;
		} else {
			return false;
		}
	}

	bool ReadInpFile::readEps()
	{
		const boost::optional<const long double &> p(
			readData<long double>("eps", pdata_->EPS_DEFAULT));
		if (p) {
			pdata_->eps = *p;

			return true;
		} else {
			return false;
		}
	}

	bool ReadInpFile::readType()
	{
		const ci_string solvetype(readData("solve.type", ci_string("Bulirsch_Stoer")));

		if (solvetype.empty())
			return false;

		const array<const ci_string, 5> stypeary =
			{ ci_string("Default"), ci_string("Mod_Euler"),
			  ci_string("Runge_Kutta"), ci_string("RK_AdapStep"),
			  ci_string("Bulirsch_Stoer") };

		int i;
		array<const ci_string, 5>::const_iterator citr(stypeary.begin());
		const array<const ci_string, 5>::const_iterator citr_end(stypeary.end());
		for (i = 0; citr != citr_end; ++citr, i++) {
			if (*citr == solvetype)
				break;
		}

		if (citr == stypeary.end()) {
			errMsg("solve.type", solvetype);
			return false;
		} else if (!i) {
			pdata_->stype = Data::def_solve_type;
		} else {
			pdata_->stype = boost::numeric_cast<const Data::solve_type>(i - 1);
		}

		return true;
	}

	bool ReadInpFile::readLowerE()
	{
		const boost::optional<const ci_string &> val(readDataAuto("search.LowerE"));

		if (val) {
			if (!val->empty()) {
				try {
					pdata_->search_lowerE = boost::optional<const long double>(
						boost::lexical_cast<long double>(val->c_str()));
				} catch (const boost::bad_lexical_cast &) {
					errMsg("search.LowerE", *val);
					return false;
				}
			} else {
				pdata_->search_lowerE = boost::none;
			}
		} else {
			return false;
		}
		
		i_++;
		return true;
	}

	bool ReadInpFile::readNumofp()
	{
		const boost::optional<std::size_t> p(
			readData<std::size_t>("num.of.partition",
								  Data::NUM_OF_PARTITION_DEFAULT));
		if (p) {
			pdata_->num_of_partiton = *p;

			return true;
		} else {
			return false;
		}
	}

	bool ReadInpFile::readRatio()
	{
		const boost::optional<const long double &> p(
			readData<long double>("matching.point.ratio",
								  pdata_->MAT_PO_RATIO_DEFAULT));
		if (p) {
			pdata_->mat_po_ratio = *p;

			return true;
		} else {
			return false;
		}
	}
}
