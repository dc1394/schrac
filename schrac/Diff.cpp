#include "Diff.h"
#include <algorithm>
#include <stdexcept>
#include <boost/numeric/odeint.hpp>
#include <gsl/gsl_linalg.h>

namespace schrac {
	// constructor
    Diff::Diff(double E, std::shared_ptr<Data> const & pdata, double TINY)
	 :	pdata_(pdata), pdiffdata_(std::make_shared<DiffData>(E, pdata, TINY))
	{
		switch (pdata_->eq_type_) {
        case Data::Eq_type::SCH:
		    dM_dx = &dM_dx_sch;
			break;

        case Data::Eq_type::SDIRAC:
			dM_dx = &dM_dx_sdirac;
			break;

        case Data::Eq_type::DIRAC:
			dM_dx = &dM_dx_dirac;
			break;

			default:
			BOOST_ASSERT(!"何かがおかしい！！");
			break;
		}

        a_init();
		Initialize(E);
	}

	void Diff::a_init()
	{
        std::array<double, DiffData::AVECSIZE * DiffData::AVECSIZE> a;
		myarray b;

		for (std::size_t i = 0; i < DiffData::AVECSIZE; i++) {
			auto rtmp = 1.0;

			for (std::size_t j = 0; j < DiffData::AVECSIZE; j++) {
                a[DiffData::AVECSIZE * i + j] = rtmp;
				rtmp *= pdiffdata_->RV_O_[i];
			}

			b[i] = pdiffdata_->VP_O_[i];
		}
            
        pdiffdata_->V_A_ = solve_linear_equ(std::move(a), std::move(b));
	}

	void Diff::b_init()
	{
		pdiffdata_->V_B_[0] = 1.0;
		pdiffdata_->V_B_[1] = 0.0;
		pdiffdata_->V_B_[2] = (pdiffdata_->V_A_[0] - pdiffdata_->E_) /
							 static_cast<double>(2 * pdata_->l_ + 3) * pdiffdata_->V_B_[0];
		pdiffdata_->V_B_[3] = pdiffdata_->V_A_[1] / static_cast<double>(3 * pdata_->l_ + 6) *
							 pdiffdata_->V_B_[0];
		pdiffdata_->V_B_[4] = (pdiffdata_->V_A_[0] * pdiffdata_->V_B_[2] + pdiffdata_->V_A_[2] * pdiffdata_->V_B_[0] -
							  pdiffdata_->E_ * pdiffdata_->V_B_[2]) /
							  static_cast<double>(4 * pdata_->l_ + 10);
	}

	void Diff::Initialize(double E)
	{
		pdiffdata_->E_ = E;			// エネルギーを代入
		pdiffdata_->thisnode_ = 0;	// ノード数初期化
		b_init();					// ベクトルBを作成

		// 初期値の作成
		// 原点から解く方の初期化
		init_LM_O();
		// 無限遠から解く方の初期化
		init_LM_I();
	}

	void Diff::init_LM_O()
	{
		pdiffdata_->LO_[0] = pdiffdata_->V_B_[DiffData::BVECSIZE - 1];
		pdiffdata_->MO_[0] = 4.0 * pdiffdata_->V_B_[DiffData::BVECSIZE - 1];

		for (int i = DiffData::BVECSIZE - 2; i >= 0; i--) {
			pdiffdata_->LO_[0] *= pdiffdata_->RV_O_[0];
			pdiffdata_->LO_[0] += pdiffdata_->V_B_[i];
		}

		for (int i = DiffData::BVECSIZE - 2; i > 0; i--) {
			pdiffdata_->MO_[0] *= pdiffdata_->RV_O_[0];
			pdiffdata_->MO_[0] += static_cast<double>(i) * pdiffdata_->V_B_[i];
		}
		pdiffdata_->MO_[0] *= pdiffdata_->RV_O_[0];
	}

	void Diff::init_LM_I()
	{
		auto const rmax = pdiffdata_->RV_I_[0];
		auto const rmaxm = pdiffdata_->RV_I_[1];
		auto const a = std::sqrt(- 2.0 * pdiffdata_->E_);
		auto const d = std::exp(- a * rmax);
		auto const dm = std::exp(- a * rmaxm);

		pdiffdata_->LI_[0] = d / schrac::pow(rmax, pdata_->l_ + 1);
		pdiffdata_->LI_[1] = dm / schrac::pow(rmaxm, pdata_->l_ + 1);
		if (pdiffdata_->LI_[0] < MINV) {
			pdiffdata_->LI_[0] = MINV;
			pdiffdata_->LI_[1] = MINV;
		}

		pdiffdata_->MI_[0] = - pdiffdata_->LI_[0] *  
						  (a * rmax + static_cast<double>(pdata_->l_ + 1));
		pdiffdata_->MI_[1] = - pdiffdata_->LI_[1] *  
						  (a * rmaxm + static_cast<double>(pdata_->l_ + 1));
		if (std::fabs(pdiffdata_->MI_[0]) < MINV) {
			pdiffdata_->MI_[0] = - MINV;
			pdiffdata_->MI_[1] = - MINV;
		}
	}

    Diff::myarray Diff::solve_linear_equ(std::array<double, DiffData::AVECSIZE * DiffData::AVECSIZE> a, myarray b) const
    {
        // save original handler, install new handler
        auto old_handler = gsl_set_error_handler(
            [](const char * reason, const char * file, std::int32_t line, std::int32_t)
        {
            auto const str = std::string(reason) + "\nFile: " + file + "\nline: " + std::to_string(line);
            throw std::runtime_error(str);
        });

        auto m = gsl_matrix_view_array(a.data(), DiffData::AVECSIZE, DiffData::AVECSIZE);
        auto const v = gsl_vector_view_array(b.data(), DiffData::AVECSIZE);

        auto const gsl_vector_deleter = [](gsl_vector * p)
        {
            gsl_vector_free(p);
        };
        std::unique_ptr<gsl_vector, decltype(gsl_vector_deleter)> x(
            gsl_vector_alloc(DiffData::AVECSIZE),
            gsl_vector_deleter);

        auto const gsl_perm_deleter = [](gsl_permutation * p)
        {
            gsl_permutation_free(p);
        };
        std::unique_ptr<gsl_permutation, decltype(gsl_perm_deleter)> p(
            gsl_permutation_alloc(DiffData::AVECSIZE),
            gsl_perm_deleter);

        std::int32_t s;
        gsl_linalg_LU_decomp(&m.matrix, p.get(), &s);
        gsl_linalg_LU_solve(&m.matrix, p.get(), &v.vector, x.get());
        myarray solution{};
        std::copy(x->data, x->data + DiffData::AVECSIZE, solution.begin());

        // restore original handler
        gsl_set_error_handler(old_handler);

        return std::move(solution);
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
			solve_diff_equ_O();
			//solve_diff_equ_I();
		//}

		return true;
	}

    bool Diff::solve_diff_equ_O()
    {
        using namespace boost::numeric::odeint;

        bulirsch_stoer<std::array<double, 2>> stepper(pdata_->eps_, pdata_->eps_);
        std::array<double, 2> initial_val = { pdiffdata_->MO_[0], pdiffdata_->LO_[0] };
        integrate_const(
            stepper,
            [this](std::array<double, 2> const & f, std::array<double, 2> & dfdx, double x) { return derivs(f, dfdx, x); },
            initial_val,
            pdiffdata_->XV_O_[0],
            pdiffdata_->XV_O_[pdiffdata_->MP_O_],
            pdiffdata_->DX_);

        return true;
    }

    bool Diff::solve_diff_equ_I()
    {
        return true;
    }

	/****************************************************
		the usual radial differential equation without 
        relativistic corrections 

        dM/dx = -(2NL + 1)M + 2r^2(V - ep)L
        dL/dx = M
	*****************************************************/
	double dM_dx_sch(double x, double L, double M,
						  const std::shared_ptr<DiffData> & pdiffdata)
	{
		double r = std::exp(x);

		return - (2.0 * static_cast<double>(pdiffdata->pdata_->l_) + 1.0) * M +
			      2.0 * sqr(r) * (fnc_V(r, pdiffdata) - pdiffdata->E_) * L;
	}

	/***********************************************************
		the scalar relativistic radial differential equation
        ref: D.D.Koelling and B.N.Harmon,
             J.Phys.C: Solid State Phys. 10, 3107 (1977).

		d = alpha^2/2*r/MQ*dV/dr
		dM/dx = -(2*NL + 1 + d)M
				+(2*r^2*MQ(V - ep) - d*NL)L
	************************************************************/
	double dM_dx_sdirac(double x, double L, double M,
							 const std::shared_ptr<DiffData> & pdiffdata)
	{
		double r = std::exp(x);

		double mass = 1.0 + Data::al2half * (pdiffdata->E_ - fnc_V(r, pdiffdata));
		double d = Data::al2half * r / mass * dV_dr(r, pdiffdata);
		double l = static_cast<double>(pdiffdata->pdata_->l_);

		// scaler treatment
		double d1 = - (2.0 * l + 1.0 + d) * M;
		double d2 = (2.0 * sqr(r) * mass * (fnc_V(r, pdiffdata) - pdiffdata->E_) -
								d * l) * L;
		return d1 + d2;
	}

	/****************************************************
		the radial differential equation with a full
        relativistic treatment

		d = alpha^2/2*r/MQ*dV/dr
        dM/dx = - (2NL + 1 + d)M 
                + 2r^2*MQ*(V - ep)*L
                - d*(NL+1+kappa)*L
        dL/dx = M
	*****************************************************/
	double dM_dx_dirac(double x, double L, double M,
							const std::shared_ptr<DiffData> & pdiffdata)
	{
		double r = std::exp(x);

		double mass = 1.0 + Data::al2half * (pdiffdata->E_ - fnc_V(r, pdiffdata));
		double d = Data::al2half * r / mass * dV_dr(r, pdiffdata);
		const std::shared_ptr<const Data> & pdata(pdiffdata->pdata_);
		double l = static_cast<double>(pdata->l_);

		// dependence on all angular momentum
		double d1 = - (2.0 * l + 1.0 + d) * M;
		double d2 = (2.0 * sqr(r) * mass * (fnc_V(r, pdiffdata) - pdiffdata->E_) -
								d * (l + 1.0 + pdata->kappa_)) * L;
		return d1 + d2;
	}


    void Diff::derivs(std::array<double, 2> const & f, std::array<double, 2> & dfdx, double x) const
    {
        // M = dL / dx
        dfdx[0] = dL_dx(f[1]);

        // L = dM / dx
        dfdx[1] = Diff::dM_dx(x, f[0], f[1], pdiffdata_);
    }

    double fnc_V(double r, const std::shared_ptr<DiffData> & pdiffdata)
    {
        return -pdiffdata->Z_ / r;
    }

    double dV_dr(double r, const std::shared_ptr<DiffData> & pdiffdata)
    {
        return pdiffdata->Z_ / (r * r);
    }

    double dL_dx(double M)
    {
        return M;
    }

    Diff::mytuple Diff::getMPval() const
    {
        std::array<double, 2> L, M;

        L[0] = pdiffdata_->LO_[pdiffdata_->MP_O_];
        L[1] = pdiffdata_->LI_[pdiffdata_->MP_I_];
        M[0] = pdiffdata_->MO_[pdiffdata_->MP_O_];
        M[1] = pdiffdata_->MI_[pdiffdata_->MP_I_];

        return std::make_pair(L, M);
    }
}
