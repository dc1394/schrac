#include "EigenValueSearch.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <boost/cast.hpp>

namespace schrac {
	EigenValueSearch::EigenValueSearch(std::pair<std::string, bool> const & arg)
	 :	loop_(1), noden_(false)
	{
		ReadInputFile rif(arg);			// ファイルを読み込む
		rif.readFile();
		pdata_ = rif.getpData();
		
		msg();

        eps_ = pdata_->eps_;
		tol_ = eps_ * 10.0;

		init();
		setoutstream();
	}

    /* public method */

    const std::shared_ptr<Data> & EigenValueSearch::getpData() const
    {
        return pdata_;
    }

    const std::shared_ptr<Diff> & EigenValueSearch::getpDiff() const
    {
        return pdiff_;
    }

    bool EigenValueSearch::search()
    {
        for (; loop_ < EigenValueSearch::EVALSEARCHMAX; loop_++) {
            if (!rough_search()) {
                return false;
            }

            /*bool b;
            try {
                b = zbrent();
            }
            catch (const std::runtime_error & e) {
                std::cerr << e.what() << std::endl;
                return false;
            }

            if (b && noden_) {
                info(E);
                return true;
            }
            else {
                E += DE;

                if (E > 0.0) {
                    return false;
                }

                pdiff_->Initialize(E);
            }*/
        }

        return false;
    }

	/*	private method	*/
	
    boost::optional<double> EigenValueSearch::fnc_D()
    {
        if (!pdiff_->solve_diff_equ()) {
            std::cerr << "微分方程式が正常に解けませんでした。" << std::endl;
            return boost::none;
        }

        noden_ = pdiffdata_->node_ == pdiffdata_->thisnode_;

        std::array<double, 2> L, M;
        std::tie(L, M) = pdiff_->getMPval();
        return boost::optional<double>(M[0] - (L[0] / L[1]) * M[1]);
    }

    inline void EigenValueSearch::info() const
    {
        std::cout << "i = " << loop_ << ", D = " << Dold << ", node = "
            << pdiffdata_->thisnode_;

        if (noden_) {
            std::cout << " (OK)" << std::endl;
        } 
        else {
            std::cout << " (NG)" << std::endl;
        }
    }

    void EigenValueSearch::info(double E) const
    {
        std::cout << "ノード数が一致する固有値を発見しました！" << std::endl;
        //std::cout << "E(厳密)   = " << Erough_exact_ << " (Hartree)" << std::endl;
        std::cout << "E(計算値) = " << E << " (Hartree)" << std::endl;

    }

    inline void EigenValueSearch::info(double b, double fb) const
    {
        std::cout << "i = " << loop_ << ", D = "
            << fb << ", E = " << b
            << ", node = "
            << pdiffdata_->thisnode_;

        if (noden_) {
            std::cout << " (OK)" << std::endl;
        } 
        else {
            std::cout << " (NG)" << std::endl;
        }
    }

	void EigenValueSearch::init()
	{
		switch (pdata_->eq_type_) {
        case Data::Eq_type::SCH:
				Erough_exact_ = Eexact_sch(pdata_);
			break;

        case Data::Eq_type::SDIRAC:
				Erough_exact_ = Eexact_sdirac(pdata_);
			break;

        case Data::Eq_type::DIRAC:
				Erough_exact_ = Eexact_dirac(pdata_);
			break;

        default:
				BOOST_ASSERT("何かがおかしい！！");
			break;
		}

		if (pdata_->search_lowerE_) {
			E = *pdata_->search_lowerE_;
			DE = - E / static_cast<double>(pdata_->num_of_partition_);
		} else {
			E = Erough_exact_;
			DE = - E / static_cast<const double>(pdata_->num_of_partition_);
			E -= 3.0 * DE;
		}

		switch (pdata_->solver_type_) {
        case Data::Solver_type::BULIRSCH_STOER:
            pdiff_ = std::make_shared<Diff>(E, pdata_, TINY);
			break;

			//case Data::RK_ADAPSTEP:
			//	pdiff_ = make_shared<RK_AdapStep>(pdata_, E, TINY);
			//break;

			//case Data::BULIRSCH_STOER:
			//	pdiff_ = make_shared<Bulirsch_Stoer>(pdata_, E, TINY);
			//break;

			default:
				BOOST_ASSERT(!"何かがおかしい！！");
			break;
		}

		pdiffdata_ = pdiff_->getpDiffData();
	}

    void EigenValueSearch::msg() const
    {
        std::cout << pdata_->chemical_symbol_;

        switch (static_cast<std::int32_t>(pdata_->Z_)) {
        case 1:
            std::cout << "原子の";
            break;

        default:
            std::cout << "イオンの";
            break;
        }

        std::cout << pdata_->orbital_ << "軌道";

        if (pdata_->eq_type_ == Data::Eq_type::DIRAC && pdata_->spin_orbital_ == Data::ALPHA) {
            std::cout << "スピン上向きの";
        }
        else if (pdata_->eq_type_ == Data::Eq_type::DIRAC && pdata_->spin_orbital_ == Data::BETA) {
            std::cout << "スピン下向きの";
        }

        std::cout << "波動関数とエネルギー固有値を計算します。\n" << std::endl;
    }

	bool EigenValueSearch::rough_search()
	{
		const boost::optional<const double> pDold(fnc_D());
		if (pDold)
			Dold = *pDold;
		else
			return false;

		info();

		loop_++;

        for (; loop_ < EVALSEARCHMAX; loop_++) {
   			E += DE;

			if (E > 0.0)
				return false;

			pdiff_->Initialize(E);

			auto const pDnew = fnc_D();
			double Dnew;
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

		return loop_ != EVALSEARCHMAX;
	}


    void EigenValueSearch::setoutstream() const
    {
        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        std::cout << std::setprecision(
            boost::numeric_cast<std::streamsize>(std::fabs(std::log10(pdata_->eps_))) - 2);
    }

	//bool EigenValueSearch::zbrent()
	//{
	//	double a = Emin, b = Emax, c = Emax, d = 0.0, e = 0.0;
	//	double fa = Dmin, fb = Dmax;
	//	double fc = fb;

	//	if (fa * fb > 0.0) {
	//		BOOST_ASSERT(!"Root must be bracketed in zbrent");
	//	}

	//	for (; loop_ < EVALSEARCHMAX; loop_++) {
	//		info(b, fb);

	//		if (std::fabs(fb) > HUGE) {
	//			std::cerr << "関数Dに極があります。極を無視して探索を続けます..." << std::endl;
	//			return false;
	//		}

	//		if (fb * fc > 0.0) {
	//			c = a;														// a, b, cの名前を付け替えて区間幅調整
	//			fc = fa;
	//			e = d = b - a;
	//		}
	//		if (std::fabs(fc) < std::fabs(fb)) {
	//			a = b;
	//			b = c;
	//			c = a;
	//			fa = fb;
	//			fb = fc;
	//			fc = fa;
	//		}

	//		const double tol1 = 2.0 * EPS * std::fabs(b) + 0.5 * TOL;	// 収束のチェック
	//		const double xm = 0.5 * (c - b);

	//		if (std::fabs(xm) <= tol1 || std::fabs(fb) < TINY) {
	//			E = b;
	//			return true;
	//		}

	//		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
	//			const double s = fb / fa;								// 逆二乗補間を試みる
	//			double q, r, p;

	//			if (std::fabs(a - c) < TINY) {
	//				p = 2.0 * xm * s;
	//				q = 1.0 - s;
	//			} else {
	//				q = fa / fc;
	//				r = fb / fc;
	//				p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
	//				q = (q - 1.0) * (r - 1.0) * (s - 1.0);
	//			}
	//			if (p > 0.0)
	//				q = - q;												// 区間内かどうかチェック

	//			p = std::fabs(p);
	//			const double min1 = 3.0 * xm * q - std::fabs(tol1 * q);
	//			const double min2 = std::fabs(e * q);
	//			const double min = min1 < min2 ? min1 : min2;

	//			if (2.0 * p < min) {
	//				e = d;													// 補間値を採用
	//				d = p / q;
	//			} else {
	//				d = xm;													// 補間失敗、二分法を使う
	//				e = d;
	//			}
	//		} else {														// 区間幅の減少が遅すぎるので二分法を使う
	//			d = xm;
	//			e = d;
	//		}

	//		a = b;															// 区間幅の最良値をaに移す
	//		fa = fb;

	//		if (fabs(d) > tol1)												// 新しい根の候補を計算
	//			b += d;
	//		else
	//			b += sign<double>(tol1, xm);
	//		
	//		pdiff_->Initialize(b);
	//		const boost::optional<const double> pfb(fnc_D());
	//		if (pfb)
	//			fb = *pfb;
	//		else
	//			throw std::runtime_error("");
	//	}

	//	throw std::runtime_error("Maximum number of iterations exceeded in zbrent");
	//}
	
    // #endregion privateメンバ関数 

    // #region 非メンバ関数

    double Eexact_dirac(std::shared_ptr<Data> const & pdata)
    {
        auto const nr = static_cast<double>(pdata->n_) - pdata->j_ - 0.5;
        auto const lambda = std::sqrt(sqr(pdata->kappa_) -
            sqr(pdata->Z_ / Data::c));
        auto const denominator = std::sqrt(sqr(nr + lambda) +
            sqr(pdata->Z_ / Data::c));

        return ((nr + lambda) / denominator - 1.0) * sqr(Data::c);
    }

    double Eexact_sch(std::shared_ptr<Data> const & pdata)
    {
        return -sqr(pdata->Z_ /
            static_cast<double>(pdata->n_)) / 2.0;
    }

    double Eexact_sdirac(std::shared_ptr<Data> const & pdata)
    {
        auto const nr = static_cast<double>(pdata->n_) - 1.0;
        auto const lambda = std::sqrt(1.0 - sqr(pdata->Z_ / Data::c));
        auto const denominator = std::sqrt(sqr(nr + lambda) +
            sqr(pdata->Z_ / Data::c));

        return ((nr + lambda) / denominator - 1.0) * sqr(Data::c);
    }

    // #endregion 非メンバ関数
}
