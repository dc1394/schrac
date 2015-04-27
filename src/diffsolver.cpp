/*! \file diffsolver.cpp
    \brief 微分方程式を解くクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD-2 License.
*/

#include "diffsolver.h"
#include <algorithm>                    // for std::copy
#include <stdexcept>                    // for std::runtime_error
#include <utility>                      // for std::move
#include <boost/cast.hpp>               // for boost::numeric_cast
#include <boost/numeric/odeint.hpp>     // for boost::numeric::odeint
#include <tbb/parallel_invoke.h>        // for tbb::parallel_invoke

namespace schrac {
    using namespace boost::numeric::odeint;

    // #region 型エイリアス

    using error_stepper_type = runge_kutta_dopri5< myarray >;

    // #endregion 型エイリアス

    // #region コンストラクタ

    DiffSolver::DiffSolver(std::shared_ptr<Data> const & pdata, std::shared_ptr<DiffData> const & pdiffdata) :
        DiffSolver(pdata, pdiffdata, nullptr, nullptr)
    {
    }

    DiffSolver::DiffSolver(std::shared_ptr<Data> const & pdata, std::shared_ptr<DiffData> const & pdiffdata, std::shared_ptr<Rho> const & prho, std::shared_ptr<Vhartree> const & pvh) :
        PDiffData([this]() { return pdiffdata_; }, nullptr),
        pdata_(pdata),
        pdiffdata_(pdiffdata),
        prho_(prho),
        pvh_(pvh),
        pvh2_(pvh ? std::make_shared<Vhartree>(*pvh_) : nullptr)
    {
        if (pdata_->chemical_symbol_ == Data::Chemical_Symbol[0]) {
            V_ = [this](double r)
            {
                return -pdiffdata_->Z_ / r;
            };
                        
            V2_ = V_;

            dV_dr_ = [this](double r)
            {
                return pdiffdata_->Z_ / (r * r);
            };

            dV_dr2_ = dV_dr_;
        } else {
            V_ = [this](double r)
            {
                return -pdiffdata_->Z_ / r + pvh_->vhartree(r);
            };

            V2_ = [this](double r)
            {
                return -pdiffdata_->Z_ / r + pvh2_->vhartree(r);
            };

            dV_dr_ = [this](double r)
            {
                return pdiffdata_->Z_ / (r * r) + pvh_->dvhartree_dr(r);
            };

            dV_dr2_ = [this](double r)
            {
                return pdiffdata_->Z_ / (r * r) + pvh2_->dvhartree_dr(r);
            };
        };
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
        pdiffdata_->E_ = E;         // エネルギーを代入
        pdiffdata_->thisnode_ = 0;  // ノード数初期化
        am_evaluate();              // am_を求める
        bm_evaluate();              // bm_を求める

        pdiffdata_->li_.clear();
        pdiffdata_->mi_.clear();
        pdiffdata_->lo_.clear();
        pdiffdata_->mo_.clear();
    }

    void DiffSolver::solve_diff_equ()
    {
        if (pdata_->usetbb_) {
            switch (pdata_->solver_type_) {
            case Data::Solver_type::ADAMS_BASHFORTH_MOULTON:
                tbb::parallel_invoke(
                    [this]{ solve_diff_equ_o(adams_bashforth_moulton< 2, myarray >(), V_, dV_dr_); },
                    [this]{ solve_diff_equ_i(adams_bashforth_moulton< 2, myarray >(), V2_, dV_dr2_); });
                break;

            case Data::Solver_type::BULIRSCH_STOER:
                tbb::parallel_invoke(
                    [this]{ solve_diff_equ_o(bulirsch_stoer < myarray >(pdata_->eps_, pdata_->eps_), V_, dV_dr_); },
					[this]{ solve_diff_equ_i(bulirsch_stoer < myarray >(pdata_->eps_, pdata_->eps_), V2_, dV_dr2_); });
                break;

            case Data::Solver_type::CONTROLLED_RUNGE_KUTTA:
                tbb::parallel_invoke(
					[this]{ solve_diff_equ_o(make_controlled(pdata_->eps_, pdata_->eps_, error_stepper_type()), V_, dV_dr_); },
					[this]{ solve_diff_equ_i(make_controlled(pdata_->eps_, pdata_->eps_, error_stepper_type()), V2_, dV_dr2_); });
                break;

            default:
                BOOST_ASSERT(!"何かがおかしい！");
                break;
            }
        }
        else {
            switch (pdata_->solver_type_) {
            case Data::Solver_type::ADAMS_BASHFORTH_MOULTON:
                solve_diff_equ_o(adams_bashforth_moulton< 2, myarray >(), V_, dV_dr_);
                solve_diff_equ_i(adams_bashforth_moulton< 2, myarray >(), V_, dV_dr_);
                break;

            case Data::Solver_type::BULIRSCH_STOER:
				solve_diff_equ_o(bulirsch_stoer < myarray >(pdata_->eps_, pdata_->eps_), V_, dV_dr_);
				solve_diff_equ_i(bulirsch_stoer < myarray >(pdata_->eps_, pdata_->eps_), V_, dV_dr_);
                break;

            case Data::Solver_type::CONTROLLED_RUNGE_KUTTA:
				solve_diff_equ_o(make_controlled(pdata_->eps_, pdata_->eps_, error_stepper_type()), V_, dV_dr_);
				solve_diff_equ_i(make_controlled(pdata_->eps_, pdata_->eps_, error_stepper_type()), V_, dV_dr_);
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
			solve_poisson_run(bulirsch_stoer < myarray >(pdata_->eps_, pdata_->eps_));
            break;

        case Data::Solver_type::CONTROLLED_RUNGE_KUTTA:
			solve_poisson_run(make_controlled(pdata_->eps_, pdata_->eps_, error_stepper_type()));
            break;

        default:
            BOOST_ASSERT(!"何かがおかしい！");
            break;
        }
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数

    void DiffSolver::am_evaluate()
    {
        std::array<double, AMMAX * AMMAX> a;
        myvector b;

        for (std::size_t i = 0; i < AMMAX; i++) {
            auto rtmp = 1.0;

            for (std::size_t j = 0; j < AMMAX; j++) {
                a[AMMAX * i + j] = rtmp;
                rtmp *= pdiffdata_->r_mesh_[i];
            }

            b[i] = V_(std::exp(pdata_->xmin_ + static_cast<double>(i) * pdiffdata_->dx_));    
        }
            
        am_ = solve_linear_equ(a, b);
    }

    void DiffSolver::bm_evaluate()
    {
        bm_[0] = 1.0;
        bm_[1] = 0.0;
        bm_[2] = (am_[0] - pdiffdata_->E_) / static_cast<double>(2 * pdata_->l_ + 3) * bm_[0];
        bm_[3] = am_[1] / static_cast<double>(3 * pdata_->l_ + 6) * bm_[0];
        bm_[4] = (am_[0] * bm_[2] + am_[2] * bm_[0] - pdiffdata_->E_ * bm_[2]) / static_cast<double>(4 * pdata_->l_ + 10);
    }

    void DiffSolver::derivs(myarray const & f, myarray & dfdx, double x, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr) const
    {
        auto const dL_dx = [](double M) { return M; };

        // dL / dx = M
        dfdx[0] = dL_dx(f[1]);
                
        switch (pdata_->eq_type_) {
        case Data::Eq_type::DIRAC:
            // dM / dx 
            dfdx[1] = dM_dx_dirac(f[0], f[1], x, V, dV_dr);
            break;

        case Data::Eq_type::SCH:
            // dM / dx 
            dfdx[1] = dM_dx_sch(f[0], f[1], x, V);
            break;

        case Data::Eq_type::SDIRAC:
            // dM / dx 
            dfdx[1] = dM_dx_sdirac(f[0], f[1], x, V, dV_dr);
            break;

        default:
            BOOST_ASSERT(!"何かがおかしい！！");
            break;
        }
    }

    double DiffSolver::dM_dx_dirac(double L, double M, double x, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr) const
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

    double DiffSolver::dM_dx_sch(double L, double M, double x, std::function<double(double)> const & V) const
    {
        auto const r = std::exp(x);

        return -(2.0 * static_cast<double>(pdata_->l_) + 1.0) * M +
               2.0 * sqr(r) * (V(r) - pdiffdata_->E_) * L;
    }

    double DiffSolver::dM_dx_sdirac(double L, double M, double x, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr) const
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
    
    void DiffSolver::node_count(dvector const & L)
    {
        if (L.size() > 1 && (L.back() * *(++L.rbegin()) < 0.0)) {
            pdiffdata_->thisnode_++;
        }
    }

    myarray DiffSolver::req_lm_i_init_val()
    {
        auto const rmax = pdiffdata_->r_mesh_i_[0];
        auto const a = std::sqrt(-2.0 * pdiffdata_->E_);
        auto const d = std::exp(-a * rmax);

        myarray state;
        state[0] = d / std::pow(rmax, pdata_->l_ + 1);

        if (state[0] < DiffSolver::MINVALUE) {
            state[0] = DiffSolver::MINVALUE;
        }

        state[1] = -state[0] * (a + static_cast<double>(pdata_->l_ + 1) / rmax);

        if (std::fabs(state[1]) < DiffSolver::MINVALUE) {
            state[1] = -DiffSolver::MINVALUE;
        }

        return std::move(state);
    }

    myarray DiffSolver::req_lm_o_init_val()
    {
        myarray state;
        state[0] = bm_[DiffSolver::BMMAX - 1];
        state[1] = 4.0 * bm_[DiffSolver::BMMAX - 1];

        auto const cnt = static_cast<std::int32_t>(DiffSolver::BMMAX - 2);
        for (auto i = cnt; i >= 0; i--) {
            state[0] *= pdiffdata_->r_mesh_[0];
            state[0] += bm_[i];
        }

        for (auto i = cnt; i > 0; i--) {
            state[1] *= pdiffdata_->r_mesh_[0];
            state[1] += static_cast<double>(i) * bm_[i];
        }
        state[1] *= pdiffdata_->r_mesh_[0];

        return std::move(state);
    }
    
    myarray DiffSolver::req_poisson_init_val()
    {
        std::array<double, AMMAX * AMMAX> a;
        myvector b;

        for (std::size_t i = 0; i < AMMAX; i++) {
            auto rtmp = pdiffdata_->r_mesh_[i];

            for (std::size_t j = 0; j < AMMAX; j++) {
                a[AMMAX * i + j] = rtmp;
                rtmp *= pdiffdata_->r_mesh_[i];
            }

            b[i] = - pdiffdata_->r_mesh_[i] * (*prho_)(pdiffdata_->r_mesh_[i]);
        }

        auto const bn = solve_linear_equ(a, b);
        
        myarray state{ 0.0, 0.0 };
        auto const r0 = pdiffdata_->r_mesh_[0];

		state[0] = ((0.2 * bn[2] * r0 + bn[1] / 3.0) * 0.5 * r0 + bn[0] / 3.0) * 0.5 * r0;
		state[1] = ((0.9 * bn[2] * r0 + bn[1]) * r0 + bn[0]) / 6.0;

        return std::move(state);
    }

    template <typename Stepper>
    void DiffSolver::solve_poisson_run(Stepper const & stepper)
    {
        auto state = req_poisson_init_val();
        auto const loop = boost::numeric_cast<std::int32_t>(pdiffdata_->r_mesh_.size() - 1);

        std::vector<double> vhart;
        vhart.reserve(pdiffdata_->r_mesh_.size());
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

            vhart.push_back(state[0] / pdiffdata_->r_mesh_[i]);
         
        }

        vhart.push_back(state[0] / pdiffdata_->r_mesh_.back());
        pvh_->Vhart(vhart);
    }

    template <typename Stepper>
    void DiffSolver::solve_diff_equ_i(Stepper const & stepper, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr)
    {
        myarray state = req_lm_i_init_val();

        integrate_const(
            stepper,
            [this, &V, &dV_dr](myarray const & f, myarray & dfdx, double x) { return derivs(f, dfdx, x, V, dV_dr); },
            state,
            pdiffdata_->x_i_[0],
            pdiffdata_->x_i_[pdiffdata_->mp_i_] - pdiffdata_->dx_,
            - pdiffdata_->dx_,
            [this](myarray const & f, double const)
        {
            pdiffdata_->li_.push_back(f[0]);
            pdiffdata_->mi_.push_back(f[1]);
            node_count(pdiffdata_->li_);
        });
    }

    template <typename Stepper>
    void DiffSolver::solve_diff_equ_o(Stepper const & stepper, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr)
    {
        auto state = req_lm_o_init_val();

        integrate_const(
            stepper,
            [this, &V, &dV_dr](myarray const & f, myarray & dfdx, double x) { return derivs(f, dfdx, x, V, dV_dr); },
            state,
            pdiffdata_->x_o_[0],
            pdiffdata_->x_o_[pdiffdata_->mp_o_],
            pdiffdata_->dx_,
            [this](myarray const & f, double const)
        {
            pdiffdata_->lo_.push_back(f[0]);
            pdiffdata_->mo_.push_back(f[1]);
            node_count(pdiffdata_->lo_);
        });

        if (pdiffdata_->lo_.size() != static_cast<std::vector<double>::size_type>(pdiffdata_->mp_o_ + 1)) {
            pdiffdata_->lo_.pop_back();
            pdiffdata_->mo_.pop_back();

            integrate_const(
                stepper,
                [this, &V, &dV_dr](myarray const & f, myarray & dfdx, double x) { return derivs(f, dfdx, x, V, dV_dr); },
                state,
                pdiffdata_->x_o_[pdiffdata_->mp_o_],
                pdiffdata_->x_o_[pdiffdata_->mp_o_] + pdiffdata_->dx_,
                pdiffdata_->dx_,
                [this](myarray const & f, double const)
            {
                pdiffdata_->lo_.push_back(f[0]);
                pdiffdata_->mo_.push_back(f[1]);
                node_count(pdiffdata_->lo_);
            });
        }
    }

    // #endregion privateメンバ関数

    // #region templateメンバ関数の実体化

    template void DiffSolver::solve_poisson_run<adams_bashforth_moulton< 2, myarray > >(adams_bashforth_moulton< 2, myarray > const & stepper);
    template void DiffSolver::solve_poisson_run<bulirsch_stoer < myarray > >(bulirsch_stoer < myarray > const & stepper);
	template void DiffSolver::solve_poisson_run<controlled_stepper_type>(controlled_stepper_type const & stepper);
    template void DiffSolver::solve_diff_equ_i<adams_bashforth_moulton< 2, myarray > >(adams_bashforth_moulton< 2, myarray > const & stepper, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr);
    template void DiffSolver::solve_diff_equ_i<bulirsch_stoer < myarray > >(bulirsch_stoer < myarray > const & stepper, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr);
	template void DiffSolver::solve_diff_equ_i<controlled_stepper_type>(controlled_stepper_type const & stepper, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr);
    template void DiffSolver::solve_diff_equ_o<adams_bashforth_moulton< 2, myarray > >(adams_bashforth_moulton< 2, myarray > const & stepper, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr);
    template void DiffSolver::solve_diff_equ_o<bulirsch_stoer < myarray > >(bulirsch_stoer < myarray > const & stepper, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr);
	template void DiffSolver::solve_diff_equ_o<controlled_stepper_type>(controlled_stepper_type const & stepper, std::function<double(double)> const & V, std::function<double(double)> const & dV_dr);

    // #endregion templateメンバ関数の実体化
}
