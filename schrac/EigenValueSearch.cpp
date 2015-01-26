#include "EigenValueSearch.h"

namespace HydroSchDirac {
	const long double EigenValueSearch::TINY = 1.0E-30;
	const long double EigenValueSearch::HUGE = 1.0E+7;

	// constructor
	EigenValueSearch::EigenValueSearch(const tuple<const std::string, const int> & arg)
	 :	i_(1), bNoden(false)
	{
#ifdef _OPENMP
		ReadInpFile rif(arg);			// ファイルを読み込む
		rif.readFile();
		pdata_ = rif.getpData();
		
		msg();

		const int nthread = pdata_->ompthread_;
		if (nthread) {
			// OpenMPの動的調整有効
			omp_set_dynamic(true);
			// OpenMPで使用するスレッドを指定
			omp_set_num_threads(nthread);
		}

		EPS = pdata_->eps;
		TOL = EPS * 10.0;

		init();
		setExp();
#else
		BOOST_STATIC_ASSERT(false);
#endif
	}

	/*	private method	*/
	
	void EigenValueSearch::init()
	{
		switch (pdata_->eqtype) {
			case Data::SCH:
				Eexact = Eexact_sch(pdata_);
			break;

			case Data::SDIRAC:
				Eexact = Eexact_sdirac(pdata_);
			break;

			case Data::DIRAC:
				Eexact = Eexact_dirac(pdata_);
			break;

			default:
				BOOST_ASSERT("何かがおかしい！！");
			break;
		}

		if (pdata_->search_lowerE) {
			E = *(pdata_->search_lowerE);
			DE = - E / static_cast<const long double>(pdata_->num_of_partiton);
		} else {
			E = Eexact;
			DE = - E / static_cast<const long double>(pdata_->num_of_partiton);
			E -= 3.0 * DE;
		}

		switch (pdata_->stype) {
			case Data::MOD_EULER:
#ifdef _DEBUG
				pdiff_ = make_shared<ModEuler>(pdata_, E, TINY);
#else
				throw std::runtime_error("「solve.type：Mod_Euler」はデバッグ用のオプションです！");
#endif
			break;

			case Data::RUNGE_KUTTA:
				pdiff_ = make_shared<RungeKutta>(pdata_, E, TINY);
			break;

			case Data::RK_ADAPSTEP:
				pdiff_ = make_shared<RK_AdapStep>(pdata_, E, TINY);
			break;

			case Data::BULIRSCH_STOER:
				pdiff_ = make_shared<Bulirsch_Stoer>(pdata_, E, TINY);
			break;

			default:
				BOOST_ASSERT(!"何かがおかしい！！");
			break;
		}

		pdiffdata_ = pdiff_->getpDiffData();
	}

	const boost::optional<const long double> EigenValueSearch::fnc_D()
	{
		if (!pdiff_->solve_diff_equ()) {
			std::cerr << "微分方程式が正常に解けませんでした。" << std::endl;
			return boost::none;
		}

		bNoden = (pdiffdata_->node == pdiffdata_->thisnode);

		array<long double, 2> L, M;
		tie(L, M) = pdiff_->getMPval();
		const long double ratio = L[0] / L[1];
		return boost::optional<const long double>(M[0] - ratio * M[1]);
	}

	bool EigenValueSearch::rough_search()
	{
		const boost::optional<const long double> pDold(fnc_D());
		if (pDold)
			Dold = *pDold;
		else
			return false;

		info();

		i_++;

		for (; i_ < EVALSEARCHMAX; i_++) {
   			E += DE;

			if (E > 0.0)
				return false;

			pdiff_->Initialize(E);

			const boost::optional<const long double> pDnew(fnc_D());
			long double Dnew;
			if (pDnew)
				Dnew = *pDnew;
			else
				return false;
	 
			if (Dnew * Dold < 0.0) {
				Emax = E;
   				Emin = E - DE;
				Dmax = Dnew; 
				Dmin = Dold;

				break;
   			} else {
				Dold = Dnew;
   			}

			info();
		}

		return (i_ != EVALSEARCHMAX);
	}

	bool EigenValueSearch::zbrent()
	{
		long double a = Emin, b = Emax, c = Emax, d = 0.0, e = 0.0;
		long double fa = Dmin, fb = Dmax;
		long double fc = fb;

		if (fa * fb > 0.0) {
			BOOST_ASSERT(!"Root must be bracketed in zbrent");
		}

		for (; i_ < EVALSEARCHMAX; i_++) {
			info(b, fb);

			if (std::fabs(fb) > HUGE) {
				std::cerr << "関数Dに極があります。極を無視して探索を続けます..." << std::endl;
				return false;
			}

			if (fb * fc > 0.0) {
				c = a;														// a, b, cの名前を付け替えて区間幅調整
				fc = fa;
				e = d = b - a;
			}
			if (std::fabs(fc) < std::fabs(fb)) {
				a = b;
				b = c;
				c = a;
				fa = fb;
				fb = fc;
				fc = fa;
			}

			const long double tol1 = 2.0 * EPS * std::fabs(b) + 0.5 * TOL;	// 収束のチェック
			const long double xm = 0.5 * (c - b);

			if (std::fabs(xm) <= tol1 || std::fabs(fb) < TINY) {
				E = b;
				return true;
			}

			if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
				const long double s = fb / fa;								// 逆二乗補間を試みる
				long double q, r, p;

				if (std::fabs(a - c) < TINY) {
					p = 2.0 * xm * s;
					q = 1.0 - s;
				} else {
					q = fa / fc;
					r = fb / fc;
					p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
					q = (q - 1.0) * (r - 1.0) * (s - 1.0);
				}
				if (p > 0.0)
					q = - q;												// 区間内かどうかチェック

				p = std::fabs(p);
				const long double min1 = 3.0 * xm * q - std::fabs(tol1 * q);
				const long double min2 = std::fabs(e * q);
				const long double min = min1 < min2 ? min1 : min2;

				if (2.0 * p < min) {
					e = d;													// 補間値を採用
					d = p / q;
				} else {
					d = xm;													// 補間失敗、二分法を使う
					e = d;
				}
			} else {														// 区間幅の減少が遅すぎるので二分法を使う
				d = xm;
				e = d;
			}

			a = b;															// 区間幅の最良値をaに移す
			fa = fb;

			if (fabs(d) > tol1)												// 新しい根の候補を計算
				b += d;
			else
				b += sign<long double>(tol1, xm);
			
			pdiff_->Initialize(b);
			const boost::optional<const long double> pfb(fnc_D());
			if (pfb)
				fb = *pfb;
			else
				throw std::runtime_error("");
		}

		throw std::runtime_error("Maximum number of iterations exceeded in zbrent");
	}

	void EigenValueSearch::msg() const
	{
		std::cout << pdata_->atom;

		switch (pdata_->Z) {
			case 1:
				std::cout << "原子の";
			break;

			default:
				std::cout << "イオンの";
			break;
		}

		std::cout << pdata_->orbital << "軌道";

		if (pdata_->eqtype == Data::DIRAC && pdata_->spin_orbital == Data::ALPHA)
				std::cout << "スピン上向きの";
		else if (pdata_->eqtype == Data::DIRAC && pdata_->spin_orbital == Data::BETA)
				std::cout << "スピン下向きの";

		std::cout << "波動関数とエネルギー固有値を計算します。\n" << std::endl;

	}
	
	/*	public method	*/

	bool EigenValueSearch::search()
	{
		for (; i_ < EVALSEARCHMAX; i_++) {
			if (!rough_search())
				return false;

			bool b;
			try {
				b = zbrent();
			} catch (const std::runtime_error & e) {
				std::cerr << e.what() << std::endl;

				return false;
			}

			if (b && bNoden) {
				info(E);
				return true;
			} else {
				E += DE;
				if (E > 0.0)
					return false;
				pdiff_->Initialize(E);
			}
		}

		return false;
	}
}
