#include "Diff.h"
#include <boost/numeric/odeint.hpp>

namespace schrac {
	const double Diff::MINV = 1.0E-200;
	
	std::function<double(double, double, double,
			const std::shared_ptr<DiffData> &)> Diff::dM_dx;

	// constructor
    Diff::Diff(double E, std::shared_ptr<Data> const & pdata, double TINY)
	 :	pdata_(pdata), pdiffdata_(std::make_shared<DiffData>(pdata, E, TINY))
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

		if (!a_init())
			throw std::runtime_error("初期値の作成に失敗しました");
		Initialize(E);
	}

	bool Diff::a_init()
	{
		std::array<std::array<double, DiffData::AVECSIZE>, DiffData::AVECSIZE> a;
		std::array<double, DiffData::AVECSIZE> y;

		for (std::size_t i = 0; i < DiffData::AVECSIZE; i++) {
			double rtmp = 1.0;

			for (std::size_t j = 0; j < DiffData::AVECSIZE; j++) {
				a[i][j] = rtmp;
				rtmp *= pdiffdata_->RV_O_[i];
			}

			y[i] = pdiffdata_->VP_O_[i];
		}

		const boost::optional<const std::array<double, DiffData::AVECSIZE> > pV_A(
			 S_gausswp(a, y));
		if (pV_A)
			pdiffdata_->V_A_ = *pV_A;
		else
			return false;

		return true;
	}

	void Diff::b_init()
	{
		pdiffdata_->V_B_[0] = 1.0;
		pdiffdata_->V_B_[1] = 0.0;
		pdiffdata_->V_B_[2] = (pdiffdata_->V_A_[0] - pdiffdata_->E_) /
							 static_cast<const double>(2 * pdata_->l_ + 3) * pdiffdata_->V_B_[0];
		pdiffdata_->V_B_[3] = pdiffdata_->V_A_[1] / static_cast<const double>(3 * pdata_->l + 6) *
							 pdiffdata_->V_B_[0];
		pdiffdata_->V_B_[4] = (pdiffdata_->V_A_[0] * pdiffdata_->V_B_[2] + pdiffdata_->V_A_[2] * pdiffdata_->V_B_[0] -
							  pdiffdata_->E_ * pdiffdata_->V_B_[2]) /
							  static_cast<const double>(4 * pdata_->l_ + 10);
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
			pdiffdata_->MO_[0] += static_cast<const double>(i) * pdiffdata_->V_B_[i];
		}
		pdiffdata_->MO_[0] *= pdiffdata_->RV_O_[0];
	}

	void Diff::init_LM_I()
	{
		const double rmax = pdiffdata_->RV_I_[0];
		const double rmaxm = pdiffdata_->RV_I_[1];
		const double a = std::sqrt(- 2.0 * pdiffdata_->E_);
		const double d = std::exp(- a * rmax);
		const double dm = std::exp(- a * rmaxm);

		pdiffdata_->LI_[0] = d / schrac::pow(rmax, pdata_->l_ + 1);
		pdiffdata_->LI_[1] = dm / schrac::pow(rmaxm, pdata_->l_ + 1);
		if (pdiffdata_->LI_[0] < MINV) {
			pdiffdata_->LI_[0] = MINV;
			pdiffdata_->LI_[1] = MINV;
		}

		pdiffdata_->MI_[0] = - pdiffdata_->LI_[0] *  
						  (a * rmax + static_cast<const double>(pdata_->l_ + 1));
		pdiffdata_->MI_[1] = - pdiffdata_->LI_[1] *  
						  (a * rmaxm + static_cast<const double>(pdata_->l_ + 1));
		if (std::fabs(pdiffdata_->MI_[0]) < MINV) {
			pdiffdata_->MI_[0] = - MINV;
			pdiffdata_->MI_[1] = - MINV;
		}
	}

	const boost::optional<const std::array<double, DiffData::AVECSIZE> >
			Diff::S_gausswp(std::array<std::array<double, DiffData::AVECSIZE>, DiffData::AVECSIZE> & a,
								std::array<double, DiffData::AVECSIZE> & b) const
	{
		// defnition
		unsigned int j = 0, pk;
		unsigned int p[DiffData::AVECSIZE];
		std::array<double, DiffData::AVECSIZE> x;

		for (unsigned int k = 0; k < DiffData::AVECSIZE; k++) {
    		// initial pivotting number
			p[k] = k;
		}

		for (unsigned int k = 0; k < DiffData::AVECSIZE - 1; k++) {
			pk = p[k];
			double max = 0.0;
			for (unsigned int i = k; i < DiffData::AVECSIZE; i++) {
				const double aik = a[p[i]][k];
				if (std::fabs(aik) > max) {
            		max = std::fabs(aik);
					j = i;
				}
			}

			if (max < pdiffdata_->TINY_) {
				//std::cerr << "gausswp, Failed! the system is singular." << std::endl;
				return boost::none;
			}

			// exchange of pivot
			if (j != k) {
				p[k] = p[j];
				p[j] = pk;
				pk = p[k]; 
			}

			// forward elimanation
			double pivot = 1.0 / a[pk][k];

			for (unsigned int j = k + 1; j < DiffData::AVECSIZE; j++) {
				int pj = p[j]; 
				double multi = a[pj][k] * pivot;
				if (std::fabs(multi) > pdiffdata_->TINY_) {
					a[pj][k] = multi;

					for (unsigned int i = k + 1; i < DiffData::AVECSIZE; i++)
                		a[pj][i] -= multi * a[pk][i];

					b[pj] -= multi * b[pk];
				} else {
            		a[pj][k] = 0.0;
				}
			}
		}

		if (std::fabs(a[p[DiffData::AVECSIZE - 1]][DiffData::AVECSIZE - 1]) <
			pdiffdata_->TINY_) {
			//std::cerr << "gausswp, Failed! the system is singular." << std::endl;
			return boost::none;
		}

		// backward substitution
		x[DiffData::AVECSIZE - 1] = b[p[DiffData::AVECSIZE - 1]] /
									a[p[DiffData::AVECSIZE - 1]][DiffData::AVECSIZE - 1];
		for (int k = DiffData::AVECSIZE - 2; k >= 0; k--) {
			pk = p[k];
			double sum = b[pk];

			for (unsigned int i = k + 1; i < DiffData::AVECSIZE; i++)
				sum -= a[pk][i] * x[i];

				x[k] = sum / a[pk][k];
		}
	    
		return boost::optional<const std::array<double, DiffData::AVECSIZE> >(x);
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
			solve_diff_equ_I();
		//}

		return true;
	}

    bool Diff::solve_diff_equ_O()
    {
        using namespace boost::numeric::odeint;

        bulirsch_stoer<std::array<double, 2>> stepper(pdata_->eps_, pdata_->eps_);
        std::array<double, 2> initial_val = { pdiffdata_->LO_[0], pdiffdata_->MO_[0] }
        integrate_const(stepper, Diff::derivs, initial_val, , , dx_);
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
		const double r = std::exp(x);

		return - (2.0 * static_cast<const double>(pdiffdata->pdata_->l_) + 1.0) * M +
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
		const double r = std::exp(x);

		const double mass = 1.0 + Data::al2half * (pdiffdata->E_ - fnc_V(r, pdiffdata));
		const double d = Data::al2half * r / mass * dV_dr(r, pdiffdata);
		const double l = static_cast<const double>(pdiffdata->pdata_->l_);

		// scaler treatment
		const double d1 = - (2.0 * l + 1.0 + d) * M;
		const double d2 = (2.0 * sqr(r) * mass * (fnc_V(r, pdiffdata) - pdiffdata->E_) -
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
		const double r = std::exp(x);

		const double mass = 1.0 + Data::al2half * (pdiffdata->E_ - fnc_V(r, pdiffdata));
		const double d = Data::al2half * r / mass * dV_dr(r, pdiffdata);
		const std::shared_ptr<const Data> & pdata(pdiffdata->pdata_);
		const double l = static_cast<const double>(pdata->l_);

		// dependence on all angular momentum
		const double d1 = - (2.0 * l + 1.0 + d) * M;
		const double d2 = (2.0 * sqr(r) * mass * (fnc_V(r, pdiffdata) - pdiffdata->E_) -
								d * (l + 1.0 + pdata->kappa_)) * L;
		return d1 + d2;
	}
}
