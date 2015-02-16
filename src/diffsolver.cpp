/*! \file diffsolver.cpp
    \brief 微分方程式を解くクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "diffsolver.h"
#include <algorithm>                    // for std::copy
#include <stdexcept>                    // for std::runtime_error
#include <boost/cast.hpp>               // for boost::numeric_cast
#include <boost/numeric/odeint.hpp>     // for boost::numeric::odeint
#include <tbb/parallel_invoke.h>        // for tbb::parallel_invoke

namespace schrac {
    using namespace boost::numeric::odeint;

    // #region 型エイリアス

    using error_stepper_type = runge_kutta_dopri5< myarray >;
    using controlled_stepper_type = controlled_runge_kutta< error_stepper_type >;

    // #endregion 型エイリアス

	// #region コンストラクタ

    DiffSolver::DiffSolver(std::shared_ptr<Data> const & pdata, std::shared_ptr<DiffData> const & pdiffdata, std::shared_ptr<Rho> const & prho, std::shared_ptr<Vhartree> const & pvh) :
        PDiffData([this]() { return pdiffdata_; }, nullptr),
        pdata_(pdata),
        pdiffdata_(pdiffdata),
        prho_(prho),
        pvh_(pvh)
	{
        am_evaluate();
	}

    // #endregion コンストラクタ

    // #region publicメンバ関数

    DiffSolver::mypair DiffSolver::getMPval() const
    {
        myarray L, M;

        L[0] = pdiffdata_->lo_[pdiffdata_->mp_o_];
        L[1] = pdiffdata_->li_[pdiffdata_->mp_i_];
        M[0] = pdiffdata_->mo_[pdiffdata_->mp_o_];
        M[1] = pdiffdata_->mi_[pdiffdata_->mp_i_];

        return std::make_pair(L, M);
    }

    void DiffSolver::initialize(double E)
    {
        pdiffdata_->E_ = E;			// エネルギーを代入
        pdiffdata_->thisnode_ = 0;	// ノード数初期化
        bm_evaluate();				// bmを求める

        // 初期値の作成
        // 原点から解く方の初期化
        init_lm_o();

        // 無限遠から解く方の初期化
        init_lm_i();
    }

    void DiffSolver::solve_diff_equ()
    {
        if (pdata_->usetbb_) {
            switch (pdata_->solver_type_) {
            case Data::Solver_type::ADAMS_BASHFORTH_MOULTON:
                tbb::parallel_invoke(
                    [this]{ solve_diff_equ_o(adams_bashforth_moulton< 2, myarray >()); },
                    [this]{ solve_diff_equ_i(adams_bashforth_moulton< 2, myarray >()); });
                break;

            case Data::Solver_type::BULIRSCH_STOER:
                tbb::parallel_invoke(
                    [this]{ solve_diff_equ_o(bulirsch_stoer < myarray >()); },
                    [this]{ solve_diff_equ_i(bulirsch_stoer < myarray >()); });
                break;

            case Data::Solver_type::CONTROLLED_RUNGE_KUTTA:
                tbb::parallel_invoke(
                    [this]{ solve_diff_equ_o(controlled_stepper_type()); },
                    [this]{ solve_diff_equ_i(controlled_stepper_type()); });
                break;

            default:
                BOOST_ASSERT(!"何かがおかしい！");
                break;
            }
        }
        else {
            switch (pdata_->solver_type_) {
            case Data::Solver_type::ADAMS_BASHFORTH_MOULTON:
                solve_diff_equ_o(adams_bashforth_moulton< 2, myarray >());
                solve_diff_equ_i(adams_bashforth_moulton< 2, myarray >());
                break;

            case Data::Solver_type::BULIRSCH_STOER:
                solve_diff_equ_o(bulirsch_stoer < myarray >());
                solve_diff_equ_i(bulirsch_stoer < myarray >());
                break;

            case Data::Solver_type::CONTROLLED_RUNGE_KUTTA:
                solve_diff_equ_o(controlled_stepper_type());
                solve_diff_equ_i(controlled_stepper_type());
                break;

            default:
                BOOST_ASSERT(!"何かがおかしい！");
                break;
            }
        }
    }

    void DiffSolver::solve_poisson()
    {
        switch (pdata_->solver_type_) {
        case Data::Solver_type::ADAMS_BASHFORTH_MOULTON:
            solve_poisson_run(adams_bashforth_moulton< 2, myarray >()); 
            break;

        case Data::Solver_type::BULIRSCH_STOER:
            solve_poisson_run(bulirsch_stoer < myarray >());
            break;

        case Data::Solver_type::CONTROLLED_RUNGE_KUTTA:
            solve_poisson_run(controlled_stepper_type());
            break;

        default:
            BOOST_ASSERT(!"何かがおかしい！");
            break;
        }
    }

    double DiffSolver::V(double r) const
    {
        return -pdiffdata_->Z_ / r;
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数

    void DiffSolver::am_evaluate()
	{
        std::array<double, DiffSolver::AMMAX * DiffSolver::AMMAX> a;
		myvector b;

		for (std::size_t i = 0; i < DiffSolver::AMMAX; i++) {
			auto rtmp = 1.0;

			for (std::size_t j = 0; j < DiffSolver::AMMAX; j++) {
                a[DiffSolver::AMMAX * i + j] = rtmp;
				rtmp *= pdiffdata_->r_mesh_[i];
			}

			b[i] = V(pdata_->xmin_ + static_cast<double>(i) * pdiffdata_->dx_);    
		}
            
        am = solve_linear_equ(std::move(a), std::move(b));
	}

    void DiffSolver::bm_evaluate()
	{
		bm[0] = 1.0;
		bm[1] = 0.0;
		bm[2] = (am[0] - pdiffdata_->E_) / static_cast<double>(2 * pdata_->l_ + 3) * bm[0];
		bm[3] = am[1] / static_cast<double>(3 * pdata_->l_ + 6) * bm[0];
		bm[4] = (am[0] * bm[2] + am[2] * bm[0] - pdiffdata_->E_ * bm[2]) / static_cast<double>(4 * pdata_->l_ + 10);
	}

    void DiffSolver::derivs(myarray const & f, myarray & dfdx, double x) const
    {
        auto const dL_dx = [](double M) { return M; };

        // dL / dx = M
        dfdx[0] = dL_dx(f[1]);
                
        switch (pdata_->eq_type_) {
        case Data::Eq_type::DIRAC:
            // dM / dx 
            dfdx[1] = dM_dx_dirac(f[0], f[1], x);
            break;

        case Data::Eq_type::SCH:
            // dM / dx 
            dfdx[1] = dM_dx_sch(f[0], f[1], x);
            break;

        case Data::Eq_type::SDIRAC:
            // dM / dx 
            dfdx[1] = dM_dx_sdirac(f[0], f[1], x);
            break;

        default:
            BOOST_ASSERT(!"何かがおかしい！！");
            break;
        }
    }

    double DiffSolver::dM_dx_dirac(double L, double M, double x) const
    {
        auto const r = std::exp(x);

        auto const mass = 1.0 + Data::al2half * (pdiffdata_->E_ - V(r));
        auto const d = Data::al2half * r / mass * dV_dr(r);
        auto const l = static_cast<double>(pdata_->l_);

        // dependence on all angular momentum
        auto const d1 = -(2.0 * l + 1.0 + d) * M;
        auto const d2 = (2.0 * sqr(r) * mass * (V(r) - pdiffdata_->E_) -
            d * (l + 1.0 + pdata_->kappa_)) * L;

        return d1 + d2;
    }

    double DiffSolver::dM_dx_sch(double L, double M, double x) const
    {
        auto const r = std::exp(x);

        return -(2.0 * static_cast<double>(pdata_->l_) + 1.0) * M +
               2.0 * sqr(r) * (V(r) - pdiffdata_->E_) * L;
    }

    double DiffSolver::dM_dx_sdirac(double L, double M, double x) const
    {
        auto const r = std::exp(x);

        auto const mass = 1.0 + Data::al2half * (pdiffdata_->E_ - V(r));
        auto const d = Data::al2half * r / mass * dV_dr(r);
        auto const l = static_cast<double>(pdata_->l_);

        // scaler treatment
        auto const d1 = -(2.0 * l + 1.0 + d) * M;
        auto const d2 = (2.0 * sqr(r) * mass * (V(r) - pdiffdata_->E_) - d * l) * L;
        
        return d1 + d2;
    }

    double DiffSolver::dV_dr(double r) const
    {
        return pdiffdata_->Z_ / (r * r);
    }

    void DiffSolver::init_lm_i()
    {
        pdiffdata_->li_.clear();
        pdiffdata_->mi_.clear();

        auto const rmax = pdiffdata_->r_mesh_i_[0];
        auto const a = std::sqrt(-2.0 * pdiffdata_->E_);
        auto const d = std::exp(-a * rmax);
        
        pdiffdata_->li_.push_back(d / std::pow(rmax, pdata_->l_ + 1));

        if (pdiffdata_->li_[0] < DiffSolver::MINVALUE) {
            pdiffdata_->li_[0] = DiffSolver::MINVALUE;
        }

        pdiffdata_->mi_.push_back(-pdiffdata_->li_[0] *
            (a + static_cast<double>(pdata_->l_ + 1) / rmax));
        
        if (std::fabs(pdiffdata_->mi_[0]) < DiffSolver::MINVALUE) {
            pdiffdata_->mi_[0] = -DiffSolver::MINVALUE;
        }
    }

	void DiffSolver::init_lm_o()
	{
        pdiffdata_->lo_.clear();
        pdiffdata_->mo_.clear();

		pdiffdata_->lo_.push_back(bm[DiffSolver::BMMAX - 1]);
		pdiffdata_->mo_.push_back(4.0 * bm[DiffSolver::BMMAX - 1]);

        auto const cnt = static_cast<std::int32_t>(DiffSolver::BMMAX - 2);
		for (auto i = cnt; i >= 0; i--) {
			pdiffdata_->lo_[0] *= pdiffdata_->r_mesh_[0];
			pdiffdata_->lo_[0] += bm[i];
		}

		for (auto i = cnt; i > 0; i--) {
			pdiffdata_->mo_[0] *= pdiffdata_->r_mesh_[0];
			pdiffdata_->mo_[0] += static_cast<double>(i) * bm[i];
		}
		pdiffdata_->mo_[0] *= pdiffdata_->r_mesh_[0];
	}
    
    void DiffSolver::node_count(dvector const & L)
    {
        if (L.size() > 1 && (L.back() * *(++L.rbegin()) < 0.0)) {
            pdiffdata_->thisnode_++;
        }
    }

    template <typename Stepper>
    void DiffSolver::solve_poisson_run(Stepper const & stepper)
    {
        myarray state = { 0.0, 1.0 };
                
        integrate_adaptive(
            stepper,
            [this](myarray const & f, myarray & dfdx, double r) {
                dfdx[0] = f[1];
                dfdx[1] = -r * (*prho_)(r);
            },
            state,
            0.0,
            pdiffdata_->r_mesh_[0],
            pdiffdata_->r_mesh_[0]);

        auto const loop = boost::numeric_cast<std::int32_t>(pdiffdata_->r_mesh_.size() - 1);
        for (auto i = 0; i < loop; i++) {
            integrate_adaptive(
                stepper,
                [this](myarray const & f, myarray & dfdx, double r) {
                dfdx[0] = f[1];
                dfdx[1] = -r * (*prho_)(r);
            },
            state,
            pdiffdata_->r_mesh_[i],
            pdiffdata_->r_mesh_[i + 1],
            pdiffdata_->r_mesh_[i + 1] - pdiffdata_->r_mesh_[i]);

            pvh_->PVhart().push_back(state[0]);
        }
    }

    template <typename Stepper>
    void DiffSolver::solve_diff_equ_i(Stepper const & stepper)
    {
        myarray state = { pdiffdata_->li_[0], pdiffdata_->mi_[0] };

        pdiffdata_->li_.pop_back();
        pdiffdata_->mi_.pop_back();

        integrate_const(
            stepper,
            [this](myarray const & f, myarray & dfdx, double x) { return derivs(f, dfdx, x); },
            state,
            pdiffdata_->x_i_[0],
            pdiffdata_->x_i_[pdiffdata_->mp_i_] - pdiffdata_->dx_,
            - pdiffdata_->dx_,
            [this](myarray const & f, double const x)
        {
            pdiffdata_->li_.push_back(f[0]);
            pdiffdata_->mi_.push_back(f[1]);
            node_count(pdiffdata_->li_);
        });
    }

    template <typename Stepper>
    void DiffSolver::solve_diff_equ_o(Stepper const & stepper)
    {
        myarray state = { pdiffdata_->lo_[0], pdiffdata_->mo_[0] };

        pdiffdata_->lo_.pop_back();
        pdiffdata_->mo_.pop_back();

        integrate_const(
            stepper,
            [this](myarray const & f, myarray & dfdx, double x) { return derivs(f, dfdx, x); },
            state,
            pdiffdata_->x_o_[0],
            pdiffdata_->x_o_[pdiffdata_->mp_o_],
            pdiffdata_->dx_,
            [this](myarray const & f, double const x)
        {
            pdiffdata_->lo_.push_back(f[0]);
            pdiffdata_->mo_.push_back(f[1]);
            node_count(pdiffdata_->lo_);
        });
    }

    // #endregion privateメンバ関数

    // #region 非メンバ関数

    DiffSolver::myvector solve_linear_equ(std::array<double, DiffSolver::AMMAX * DiffSolver::AMMAX> a, DiffSolver::myvector b)
    {
        // save original handler, install new handler
        auto old_handler = gsl_set_error_handler(
            [](char const * reason, char const * file, std::int32_t line, std::int32_t)
        {
            auto const str = std::string(reason) + "\nFile: " + file + "\nline: " + std::to_string(line);
            throw std::runtime_error(str);
        });

        auto m = gsl_matrix_view_array(a.data(), DiffSolver::AMMAX, DiffSolver::AMMAX);
        auto const v = gsl_vector_view_array(b.data(), DiffSolver::AMMAX);

        auto const gsl_vector_deleter = [](gsl_vector * p)
        {
            gsl_vector_free(p);
        };
        std::unique_ptr<gsl_vector, decltype(gsl_vector_deleter)> x(
            gsl_vector_alloc(DiffSolver::AMMAX),
            gsl_vector_deleter);

        auto const gsl_perm_deleter = [](gsl_permutation * p)
        {
            gsl_permutation_free(p);
        };
        std::unique_ptr<gsl_permutation, decltype(gsl_perm_deleter)> p(
            gsl_permutation_alloc(DiffSolver::AMMAX),
            gsl_perm_deleter);

        std::int32_t s;
        gsl_linalg_LU_decomp(&m.matrix, p.get(), &s);
        gsl_linalg_LU_solve(&m.matrix, p.get(), &v.vector, x.get());
        DiffSolver::myvector solution{ 0.0, 0.0, 0.0 };
        std::copy(x->data, x->data + DiffSolver::AMMAX, solution.begin());

        // restore original handler
        gsl_set_error_handler(old_handler);

        return std::move(solution);
    }
    
    // #endregion 非メンバ関数

    // #region templateメンバ関数の実体化

    template void DiffSolver::solve_poisson_run<adams_bashforth_moulton< 2, myarray > >(adams_bashforth_moulton< 2, myarray > const & stepper);
    template void DiffSolver::solve_poisson_run<bulirsch_stoer < myarray > >(bulirsch_stoer < myarray > const & stepper);
    template void DiffSolver::solve_poisson_run<controlled_stepper_type>(controlled_stepper_type const & stepper);
    template void DiffSolver::solve_diff_equ_i<adams_bashforth_moulton< 2, myarray > >(adams_bashforth_moulton< 2, myarray > const & stepper);
    template void DiffSolver::solve_diff_equ_i<bulirsch_stoer < myarray > >(bulirsch_stoer < myarray > const & stepper);
    template void DiffSolver::solve_diff_equ_i<controlled_stepper_type>(controlled_stepper_type const & stepper);
    template void DiffSolver::solve_diff_equ_o<adams_bashforth_moulton< 2, myarray > >(adams_bashforth_moulton< 2, myarray > const & stepper);
    template void DiffSolver::solve_diff_equ_o<bulirsch_stoer < myarray > >(bulirsch_stoer < myarray > const & stepper);
    template void DiffSolver::solve_diff_equ_o<controlled_stepper_type>(controlled_stepper_type const & stepper);

    // #endregion templateメンバ関数の実体化
}
