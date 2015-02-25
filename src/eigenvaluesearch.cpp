/*! \file eigenvaluesearch.cpp
    \brief エネルギー固有値検索を行うクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "eigenvaluesearch.h"
#include <cmath>                // for std::fabs, std::log10
#include <iomanip>              // for std::setprecision
#include <iostream>             // for std::cot, std::cerr
#include <tuple>                // for std::tie
#include <boost/assert.hpp>     // for BOOST_ASSERT
#include <boost/cast.hpp>       // boost::numeric_cast
#include <gsl/gsl_errno.h>      // for GSL_SUCCESS
#include <gsl/gsl_roots.h>      // gsl_root_fsolver

namespace schrac {
    // #region staticメンバ変数

    bool EigenValueSearch::nodeok = false;

    // #endregion staticメンバ変数 

    // #region コンストラクタ

    EigenValueSearch::EigenValueSearch(std::shared_ptr<Data> const & pdata, std::shared_ptr<DiffData> const & pdiffdata, std::shared_ptr<Rho> const & prho, std::shared_ptr<Vhartree> const & pvh) :
        PData([this]() { return pdata_; }, nullptr),
        PDiffSolver([this]() { return pdiffsolver_; }, nullptr),
        loop_(1),
        pdata_(pdata),
        pdiffdata_(pdiffdata),
        pvh_(pvh)
    {
        initialize(prho);
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

            if (brent() && EigenValueSearch::nodeok) {
                return true;
            }
            else {
                pdiffsolver_->E_ += DE_;

                if (pdiffsolver_->E_ > 0.0) {
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
        F.params = reinterpret_cast<void *>(pdiffsolver_.get());

        gsl_root_fsolver_set(s.get(), &F, Emin_, Emax_);

        for (; loop_ < EVALSEARCHMAX; loop_++) {
            auto const ret = gsl_root_fsolver_iterate(s.get());

            switch (ret)
            {
            case GSL_EBADFUNC:
                std::cerr << "the iteration encountered a singular point where"
                    << "the function or its derivative evaluated to Inf or NaN.\n";
                    return false;
                break;

            case GSL_EZERODIV:
                std::cerr << "the derivative of the function vanished at the iteration point,"
                    << "preventing the algorithm from continuing without a division by zero.\n";
                    return false;
                break;

            default:
                break;
            }

            if (pdata_->chemical_symbol_ == Data::Chemical_Symbol[0]) {
                info(func_D(pdiffsolver_->E_, F.params), pdiffsolver_->E_);
            }

            pdiffsolver_->E_ = gsl_root_fsolver_root(s.get());
            
            auto const status = gsl_root_test_interval(
                gsl_root_fsolver_x_lower(s.get()),
                gsl_root_fsolver_x_upper(s.get()),
                0.0,
                pdata_->eps_);

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
        
    void EigenValueSearch::info(double D, double E) const
    {
        std::cout << "i = " << loop_ << ", D = "
            << D << ", E = " << E
            << ", node = "
            << pdiffdata_->thisnode_;

        if (EigenValueSearch::nodeok) {
            std::cout << " (OK)" << std::endl;
        } 
        else {
            std::cout << " (NG)" << std::endl;
        }
    }

    void EigenValueSearch::initialize(std::shared_ptr<Rho> const & prho)
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

        pdiffsolver_ = std::make_shared<DiffSolver>(pdata_, pdiffdata_, prho, pvh_);

        if (pdata_->search_lowerE_) {
            pdiffsolver_->E_ = *pdata_->search_lowerE_;
            DE_ = - pdiffsolver_->E_ / static_cast<double>(pdata_->num_of_partition_);
        } else {
            pdiffsolver_->E_ = Eapprox_;
            DE_ = - pdiffsolver_->E_ / static_cast<const double>(pdata_->num_of_partition_);
            pdiffsolver_->E_ -= 3.0 * DE_;
        }
    }

    bool EigenValueSearch::rough_search()
    {
        auto pdiffsolver = reinterpret_cast<void *>(pdiffsolver_.get());
        Dold = func_D(pdiffsolver_->E_, pdiffsolver);

        if (pdata_->chemical_symbol_ == Data::Chemical_Symbol[0]) {
            info();
        }

        ++loop_;

        for (; loop_ < EVALSEARCHMAX; loop_++) {
            pdiffsolver_->E_ += DE_;

            if (pdiffsolver_->E_ > 0.0) {
                return false;
            }

            auto const Dnew = func_D(pdiffsolver_->E_, pdiffsolver);
     
            if (Dnew * Dold < 0.0) {
                Emax_ = pdiffsolver_->E_;
                Emin_ = pdiffsolver_->E_ - DE_;

                break;
            } else {
                Dold = Dnew;
            }

            if (pdata_->chemical_symbol_ == Data::Chemical_Symbol[0]) {
                info();
            }
        }

        return loop_ != EVALSEARCHMAX;
    }

    void EigenValueSearch::setoutstream() const
    {
        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        std::cout << std::setprecision(
            boost::numeric_cast<std::streamsize>(std::fabs(std::log10(pdata_->eps_))));
    }
    
    // #endregion privateメンバ関数 

    // #region 非メンバ関数

    double func_D(double E, void * params)
    {
        auto pdiffsolver = reinterpret_cast<DiffSolver *>(params);
        pdiffsolver->initialize(E);
        pdiffsolver->solve_diff_equ();

        EigenValueSearch::nodeok = 
            pdiffsolver->PDiffData()->node_ == pdiffsolver->PDiffData()->thisnode_;

        myarray L, M;
        std::tie(L, M) = pdiffsolver->getMPval();
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
