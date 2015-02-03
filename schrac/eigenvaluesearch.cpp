﻿/*! \file eigenvaluesearch.cpp
    \brief エネルギー固有値検索を行うクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include <conio.h>
#include "eigenvaluesearch.h"
#include <cmath>                // for std::fabs, std::log10
#include <iomanip>              // for std::setprecision
#include <iostream>             // for std::cot, std::cerr
#include <tuple>                // for std::tie
#include <boost/assert.hpp>     // for BOOST_ASSERT
#include <boost/cast.hpp>       // boost::numeric_cast
#include <gsl/gsl_errno.h>      // for GSL_SUCCESS
#include <gsl/gsl_roots.h>

namespace schrac {
    // #region staticメンバ変数

    bool EigenValueSearch::nodeok = false;

    // #endregion staticメンバ変数 

    // #region コンストラクタ

	EigenValueSearch::EigenValueSearch(std::pair<std::string, bool> const & arg) :
        PData([this]() { return pdata_; }, nullptr),
        PDiff([this]() { return pdiff_; }, nullptr),
        loop_(1)
	{
		ReadInputFile rif(arg);			// ファイルを読み込む
		rif.readFile();
		pdata_ = rif.PData;
		
		msg();

		initialize();
		setoutstream();
	}

    // #endregion コンストラクタ

    // #region publicメンバ関数

    bool EigenValueSearch::search()
    {
        for (; loop_ < EigenValueSearch::EVALSEARCHMAX; loop_++) {
            if (!rough_search()) {
                return false;
            }

            bool b = false;
            try {
                b = brent();
            }
            catch (const std::runtime_error & e) {
                std::cerr << e.what() << std::endl;
                return false;
            }

            if (b && EigenValueSearch::nodeok) {
                info(E_);

                return true;
            }
            else {
                E_ += DE_;

                if (E_ > 0.0) {
                    return false;
                }
            }
        }

        return false;
    }

    // #endregion publicメンバ関数

	// #region privateメンバ関数

    bool EigenValueSearch::brent()
    {
        auto const fsolver_deleter = [](gsl_root_fsolver * s)
        {
            gsl_root_fsolver_free(s);
        };
        std::unique_ptr<gsl_root_fsolver, decltype(fsolver_deleter)> s(
            gsl_root_fsolver_alloc(gsl_root_fsolver_brent),
            fsolver_deleter);

        gsl_function F;

        F.function = &func_D;
        F.params = reinterpret_cast<void *>(pdiff_.get());

        auto x_lo = Emin_, x_hi = Emax_;
        gsl_root_fsolver_set(s.get(), &F, x_lo, x_hi);
        for (; loop_ < EVALSEARCHMAX; loop_++) {
            auto status = gsl_root_fsolver_iterate(s.get());
            E_ = gsl_root_fsolver_root(s.get());
            x_lo = gsl_root_fsolver_x_lower(s.get());
            x_hi = gsl_root_fsolver_x_upper(s.get());
            status = gsl_root_test_interval(x_lo, x_hi, 0.0, pdata_->eps_);

            if (status == GSL_SUCCESS) {
                break;
            }
        }

        return loop_ != EVALSEARCHMAX;
    }

    void EigenValueSearch::info() const
    {
        std::cout << "i = " << loop_ << ", D = " << Dold << ", node = "
            << pdiffdata_->thisnode_;

        if (EigenValueSearch::nodeok) {
            std::cout << " (OK)" << std::endl;
        } 
        else {
            std::cout << " (NG)" << std::endl;
        }
    }

    void EigenValueSearch::info(double E) const
    {
        std::cout << "ノード数が一致する固有値を発見しました！" << std::endl;
        std::cout << "E(計算値) = " << E << " (Hartree)" << std::endl;
    }

    void EigenValueSearch::info(double b, double fb) const
    {
        std::cout << "i = " << loop_ << ", D = "
            << fb << ", E_ = " << b
            << ", node = "
            << pdiffdata_->thisnode_;

        if (EigenValueSearch::nodeok) {
            std::cout << " (OK)" << std::endl;
        } 
        else {
            std::cout << " (NG)" << std::endl;
        }
    }

	void EigenValueSearch::initialize()
	{
		switch (pdata_->eq_type_) {
        case Data::Eq_type::SCH:
				Eapprox_ = Eapprox_sch(pdata_);
			break;

        case Data::Eq_type::SDIRAC:
        case Data::Eq_type::DIRAC:
				Eapprox_ = Eapprox_dirac(pdata_);
			break;

        default:
				BOOST_ASSERT("何かがおかしい！！");
			break;
		}

		if (pdata_->search_lowerE_) {
			E_ = *pdata_->search_lowerE_;
			DE_ = - E_ / static_cast<double>(pdata_->num_of_partition_);
		} else {
			E_ = Eapprox_;
			DE_ = - E_ / static_cast<const double>(pdata_->num_of_partition_);
			E_ -= 3.0 * DE_;
		}

        pdiff_ = std::make_shared<Diff>(pdata_);
        pdiffdata_ = pdiff_->PDiffData;
	}

    void EigenValueSearch::msg() const
    {
        std::cout << pdata_->chemical_symbol_
            << "原子の"
            << pdata_->orbital_
            << "軌道";

        if (pdata_->eq_type_ == Data::Eq_type::DIRAC && pdata_->spin_orbital_ == Data::ALPHA) {
            std::cout << "スピン上向きの";
        }
        else if (pdata_->eq_type_ == Data::Eq_type::DIRAC && pdata_->spin_orbital_ == Data::BETA) {
            std::cout << "スピン下向きの";
        }

        std::cout << "の波動関数とエネルギー固有値を計算します。\n" << std::endl;
    }

	bool EigenValueSearch::rough_search()
	{
        auto pdiff = reinterpret_cast<void *>(pdiff_.get());
        Dold = func_D(E_, pdiff);

		info();

		loop_++;

        auto const E0 = E_;
        for (; loop_ < EVALSEARCHMAX; loop_++) {
   			E_ += DE_;

            if (E_ > 0.0) {
                return false;
            }

			auto const Dnew = func_D(E_, pdiff);
	 
			if (Dnew * Dold < 0.0) {
				Emax_ = E_;
   				Emin_ = E_ - DE_;

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
	
    // #endregion privateメンバ関数 

    // #region 非メンバ関数

    double func_D(double E_, void * params)
    {
        auto pdiff = reinterpret_cast<Diff *>(params);
        pdiff->initialize(E_);
        pdiff->solve_diff_equ();

        EigenValueSearch::nodeok = pdiff->PDiffData()->node_ == pdiff->PDiffData()->thisnode_;

        myarray L, M;
        std::tie(L, M) = pdiff->getMPval();
        return M[0] - (L[0] / L[1]) * M[1];
    }

    double Eapprox_dirac(std::shared_ptr<Data> const & pdata)
    {
        auto const nr = static_cast<double>(pdata->n_) - pdata->j_ - 0.5;
        auto const lambda = std::sqrt(sqr(pdata->kappa_) - sqr(pdata->Z_ / Data::c));
        auto const denominator = std::sqrt(sqr(nr + lambda) + sqr(pdata->Z_ / Data::c));

        return ((nr + lambda) / denominator - 1.0) * sqr(Data::c);
    }

    double Eapprox_sch(std::shared_ptr<Data> const & pdata)
    {
        return -sqr(pdata->Z_ / static_cast<double>(pdata->n_)) / 2.0;
    }

    // #endregion 非メンバ関数
}