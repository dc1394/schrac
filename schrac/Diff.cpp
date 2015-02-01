#include "Diff.h"
#include <algorithm>
#include <stdexcept>
#include <boost/numeric/odeint.hpp>
#include <gsl/gsl_linalg.h>

namespace schrac {
	// #region コンストラクタ

    Diff::Diff(std::shared_ptr<Data> const & pdata) :
        pdata_(pdata),
        pdiffdata_(std::make_shared<DiffData>(pdata)),
        PDiffData([this]() { return pdiffdata_; }, nullptr)
	{		
        am_evaluate();
	}

    // #endregion コンストラクタ

    // #region publicメンバ関数

    Diff::mypair Diff::getMPval() const
    {
        myarray L, M;

        L[0] = pdiffdata_->lo_[pdiffdata_->mp_o_];
        L[1] = pdiffdata_->li_[pdiffdata_->mp_i_];
        M[0] = pdiffdata_->mo_[pdiffdata_->mp_o_];
        M[1] = pdiffdata_->mi_[pdiffdata_->mp_i_];

        return std::make_pair(L, M);
    }

    void Diff::initialize(double E)
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

    bool Diff::solve_diff_equ()
    {
        //if (pdata_->ompthread_) {
        /*	#pragma omp parallel sections
        {
        #pragma omp section
        {
        solve_diff_equ_O();
        }

        #pragma omp section
        {
        solve_diff_equ_I();
        }

        }
        } else {*/
        solve_diff_equ_o();
        solve_diff_equ_i();
        //}

        return true;
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数

    void Diff::am_evaluate()
	{
        std::array<double, Diff::AMMAX * Diff::AMMAX> a;
		myvector b;

		for (std::size_t i = 0; i < Diff::AMMAX; i++) {
			auto rtmp = 1.0;

			for (std::size_t j = 0; j < Diff::AMMAX; j++) {
                a[Diff::AMMAX * i + j] = rtmp;
				rtmp *= pdiffdata_->r_o_[i];
			}

			b[i] = pdiffdata_->vr_o_[i];
		}
            
        am = solve_linear_equ(std::move(a), std::move(b));
	}

    void Diff::bm_evaluate()
	{
		bm[0] = 1.0;
		bm[1] = 0.0;
		bm[2] = (am[0] - pdiffdata_->E_) / static_cast<double>(2 * pdata_->l_ + 3) * bm[0];
		bm[3] = am[1] / static_cast<double>(3 * pdata_->l_ + 6) * bm[0];
		bm[4] = (am[0] * bm[2] + am[2] * bm[0] - pdiffdata_->E_ * bm[2]) / static_cast<double>(4 * pdata_->l_ + 10);
	}

    void Diff::derivs(myarray const & f, myarray & dfdx, double x) const
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

    double Diff::dM_dx_dirac(double L, double M, double x) const
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

    double Diff::dM_dx_sch(double L, double M, double x) const
    {
        auto const r = std::exp(x);

        return -(2.0 * static_cast<double>(pdata_->l_) + 1.0) * M +
               2.0 * sqr(r) * (V(r) - pdiffdata_->E_) * L;
    }

    double Diff::dM_dx_sdirac(double L, double M, double x) const
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

    double Diff::dV_dr(double r) const
    {
        return pdiffdata_->Z_ / (r * r);
    }

    void Diff::init_lm_i()
    {
        pdiffdata_->li_.clear();
        pdiffdata_->mi_.clear();

        auto const rmax = pdiffdata_->r_i_[0];
        auto const rmaxm = pdiffdata_->r_i_[1];
        auto const a = std::sqrt(-2.0 * pdiffdata_->E_);
        auto const d = std::exp(-a * rmax);
        auto const dm = std::exp(-a * rmaxm);

        pdiffdata_->li_.push_back(d / schrac::pow(rmax, pdata_->l_ + 1));
        pdiffdata_->li_.push_back(dm / schrac::pow(rmaxm, pdata_->l_ + 1));
        if (pdiffdata_->li_[0] < Diff::MINV) {
            pdiffdata_->li_[0] = Diff::MINV;
            pdiffdata_->li_[1] = Diff::MINV;
        }

        pdiffdata_->mi_.push_back(-pdiffdata_->li_[0] *
            (a * rmax + static_cast<double>(pdata_->l_ + 1)));
        pdiffdata_->mi_.push_back(-pdiffdata_->li_[1] *
            (a * rmaxm + static_cast<double>(pdata_->l_ + 1)));
        if (std::fabs(pdiffdata_->mi_[0]) < Diff::MINV) {
            pdiffdata_->mi_[0] = -Diff::MINV;
            pdiffdata_->mi_[1] = -Diff::MINV;
        }
    }

	void Diff::init_lm_o()
	{
        pdiffdata_->lo_.clear();
        pdiffdata_->mo_.clear();

		pdiffdata_->lo_.push_back(bm[Diff::BMMAX - 1]);
		pdiffdata_->mo_.push_back(4.0 * bm[Diff::BMMAX - 1]);

		for (int i = Diff::BMMAX - 2; i >= 0; i--) {
			pdiffdata_->lo_[0] *= pdiffdata_->r_o_[0];
			pdiffdata_->lo_[0] += bm[i];
		}

		for (int i = Diff::BMMAX - 2; i > 0; i--) {
			pdiffdata_->mo_[0] *= pdiffdata_->r_o_[0];
			pdiffdata_->mo_[0] += static_cast<double>(i) * bm[i];
		}
		pdiffdata_->mo_[0] *= pdiffdata_->r_o_[0];
	}
    
    void Diff::node_count(std::int32_t i, dvector const & wf)
    {
        if (WF[i] * WF[i - 1] < 0.0) {
            thisnode_++;
        }
    }

    bool Diff::solve_diff_equ_i()
    {
        using namespace boost::numeric::odeint;

        bulirsch_stoer<myarray> stepper(pdata_->eps_, pdata_->eps_);
        myarray initial_val = { pdiffdata_->li_[1], pdiffdata_->mi_[1] };
        pdiffdata_->li_.pop_back();
        pdiffdata_->mi_.pop_back();
        std::int32_t index = 1;
        integrate_const(
            stepper,
            [this](myarray const & f, myarray & dfdx, double x) { return derivs(f, dfdx, x); },
            initial_val,
            pdiffdata_->x_i_[1],
            pdiffdata_->x_i_[pdiffdata_->mp_i_],
            - pdiffdata_->DX_,
            [this, &index](myarray const & f, double const x)
        {
            pdiffdata_->li_.push_back(f[0]);
            pdiffdata_->mi_.push_back(f[1]);
            auto const size = pdiffdata_->li_.size();
            pdiffdata_->node_count(size - 1, pdiffdata_->li_);
        });

        return true;
    }

    bool Diff::solve_diff_equ_o()
    {
        using namespace boost::numeric::odeint;

        bulirsch_stoer<myarray> stepper(pdata_->eps_, pdata_->eps_);
        myarray initial_val = { pdiffdata_->lo_[0], pdiffdata_->mo_[0] };
        pdiffdata_->lo_.pop_back();
        pdiffdata_->mo_.pop_back();
        integrate_const(
            stepper,
            [this](myarray const & f, myarray & dfdx, double x) { return derivs(f, dfdx, x); },
            initial_val,
            pdiffdata_->x_o_[0],
            pdiffdata_->x_o_[pdiffdata_->mp_o_],
            pdiffdata_->DX_,
            [this](myarray const & f, double const x)
            {
                pdiffdata_->lo_.push_back(f[0]);
                pdiffdata_->mo_.push_back(f[1]);
                auto const size = pdiffdata_->lo_.size();
                if (size > 1) {
                    pdiffdata_->node_count(size - 1, pdiffdata_->lo_);
                }
            });

        return true;
    }

    double Diff::V(double r) const
    {
        return -pdiffdata_->Z_ / r;
    }

    // #endregion privateメンバ関数

    // #region 非メンバ関数

    Diff::myvector solve_linear_equ(std::array<double, Diff::AMMAX * Diff::AMMAX> a, Diff::myvector b)
    {
        // save original handler, install new handler
        auto old_handler = gsl_set_error_handler(
            [](const char * reason, const char * file, std::int32_t line, std::int32_t)
        {
            auto const str = std::string(reason) + "\nFile: " + file + "\nline: " + std::to_string(line);
            throw std::runtime_error(str);
        });

        auto m = gsl_matrix_view_array(a.data(), Diff::AMMAX, Diff::AMMAX);
        auto const v = gsl_vector_view_array(b.data(), Diff::AMMAX);

        auto const gsl_vector_deleter = [](gsl_vector * p)
        {
            gsl_vector_free(p);
        };
        std::unique_ptr<gsl_vector, decltype(gsl_vector_deleter)> x(
            gsl_vector_alloc(Diff::AMMAX),
            gsl_vector_deleter);

        auto const gsl_perm_deleter = [](gsl_permutation * p)
        {
            gsl_permutation_free(p);
        };
        std::unique_ptr<gsl_permutation, decltype(gsl_perm_deleter)> p(
            gsl_permutation_alloc(Diff::AMMAX),
            gsl_perm_deleter);

        std::int32_t s;
        gsl_linalg_LU_decomp(&m.matrix, p.get(), &s);
        gsl_linalg_LU_solve(&m.matrix, p.get(), &v.vector, x.get());
        Diff::myvector solution{};
        std::copy(x->data, x->data + Diff::AMMAX, solution.begin());

        // restore original handler
        gsl_set_error_handler(old_handler);

        return std::move(solution);
    }

    // #endregion 非メンバ関数
}
