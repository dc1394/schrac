#include "Diff.h"

namespace HydroSchDirac {
	const long double Diff::MINV = 1.0E-200;
	
	function<long double(long double, long double, long double,
			const shared_ptr<DiffData> &)> Diff::dM_dx;

	// constructor
	Diff::Diff(const shared_ptr<const Data> & pdata, long double E, long double TINY)
	 :	pdata_(pdata), pdiffdata_(make_shared<DiffData>(pdata, E, TINY))
	{
		switch (pdata_->eqtype) {
			case Data::SCH:
				dM_dx = &dM_dx_sch;
			break;

			case Data::SDIRAC: 
				dM_dx = &dM_dx_sdirac;
			break;

			case Data::DIRAC:
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
		array<array<long double, DiffData::AVECSIZE>, DiffData::AVECSIZE> a;
		array<long double, DiffData::AVECSIZE> y;

		for (std::size_t i = 0; i < DiffData::AVECSIZE; i++) {
			long double rtmp = 1.0;

			for (std::size_t j = 0; j < DiffData::AVECSIZE; j++) {
				a[i][j] = rtmp;
				rtmp *= pdiffdata_->RV_O[i];
			}

			y[i] = pdiffdata_->VP_O[i];
		}

		const boost::optional<const array<long double, DiffData::AVECSIZE> > pV_A(
			 S_gausswp(a, y));
		if (pV_A)
			pdiffdata_->V_A = *pV_A;
		else
			return false;

		return true;
	}

	void Diff::b_init()
	{
		pdiffdata_->V_B[0] = 1.0;
		pdiffdata_->V_B[1] = 0.0;
		pdiffdata_->V_B[2] = (pdiffdata_->V_A[0] - pdiffdata_->E_) /
							 static_cast<const long double>(2 * pdata_->l + 3) * pdiffdata_->V_B[0];
		pdiffdata_->V_B[3] = pdiffdata_->V_A[1] / static_cast<const long double>(3 * pdata_->l + 6) *
							 pdiffdata_->V_B[0];
		pdiffdata_->V_B[4] = (pdiffdata_->V_A[0] * pdiffdata_->V_B[2] + pdiffdata_->V_A[2] * pdiffdata_->V_B[0] -
							  pdiffdata_->E_ * pdiffdata_->V_B[2]) /
							  static_cast<const long double>(4 * pdata_->l + 10);
	}

	void Diff::Initialize(long double E)
	{
		pdiffdata_->E_ = E;			// エネルギーを代入
		pdiffdata_->thisnode = 0;	// ノード数初期化
		b_init();					// ベクトルBを作成

		// 初期値の作成
		// 原点から解く方の初期化
		init_LM_O();
		// 無限遠から解く方の初期化
		init_LM_I();
	}

	void Diff::init_LM_O()
	{
		pdiffdata_->LO[0] = pdiffdata_->V_B[DiffData::BVECSIZE - 1];
		pdiffdata_->MO[0] = 4.0 * pdiffdata_->V_B[DiffData::BVECSIZE - 1];

		for (int i = DiffData::BVECSIZE - 2; i >= 0; i--) {
			pdiffdata_->LO[0] *= pdiffdata_->RV_O[0];
			pdiffdata_->LO[0] += pdiffdata_->V_B[i];
		}

		for (int i = DiffData::BVECSIZE - 2; i > 0; i--) {
			pdiffdata_->MO[0] *= pdiffdata_->RV_O[0];
			pdiffdata_->MO[0] += static_cast<const long double>(i) * pdiffdata_->V_B[i];
		}
		pdiffdata_->MO[0] *= pdiffdata_->RV_O[0];
	}

	void Diff::init_LM_I()
	{
		const long double rmax = pdiffdata_->RV_I[0];
		const long double rmaxm = pdiffdata_->RV_I[1];
		const long double a = std::sqrt(- 2.0 * pdiffdata_->E_);
		const long double d = std::exp(- a * rmax);
		const long double dm = std::exp(- a * rmaxm);

		pdiffdata_->LI[0] = d / HydroSchDirac::pow(rmax, pdata_->l + 1);
		pdiffdata_->LI[1] = dm / HydroSchDirac::pow(rmaxm, pdata_->l + 1);
		if (pdiffdata_->LI[0] < MINV) {
			pdiffdata_->LI[0] = MINV;
			pdiffdata_->LI[1] = MINV;
		}

		pdiffdata_->MI[0] = - pdiffdata_->LI[0] *  
						  (a * rmax + static_cast<const long double>(pdata_->l + 1));
		pdiffdata_->MI[1] = - pdiffdata_->LI[1] *  
						  (a * rmaxm + static_cast<const long double>(pdata_->l + 1));
		if (std::fabs(pdiffdata_->MI[0]) < MINV) {
			pdiffdata_->MI[0] = - MINV;
			pdiffdata_->MI[1] = - MINV;
		}
	}

	const boost::optional<const array<long double, DiffData::AVECSIZE> >
			Diff::S_gausswp(array<array<long double, DiffData::AVECSIZE>, DiffData::AVECSIZE> & a,
								array<long double, DiffData::AVECSIZE> & b) const
	{
		// defnition
		unsigned int j = 0, pk;
		unsigned int p[DiffData::AVECSIZE];
		array<long double, DiffData::AVECSIZE> x;

		for (unsigned int k = 0; k < DiffData::AVECSIZE; k++) {
    		// initial pivotting number
			p[k] = k;
		}

		for (unsigned int k = 0; k < DiffData::AVECSIZE - 1; k++) {
			pk = p[k];
			long double max = 0.0;
			for (unsigned int i = k; i < DiffData::AVECSIZE; i++) {
				const long double aik = a[p[i]][k];
				if (std::fabs(aik) > max) {
            		max = std::fabs(aik);
					j = i;
				}
			}

			if (max < pdiffdata_->TINY_) {
				std::cerr << "gausswp, Failed! the system is singular." << std::endl;
				return boost::none;
			}

			// exchange of pivot
			if (j != k) {
				p[k] = p[j];
				p[j] = pk;
				pk = p[k]; 
			}

			// forward elimanation
			long double pivot = 1.0 / a[pk][k];

			for (unsigned int j = k + 1; j < DiffData::AVECSIZE; j++) {
				int pj = p[j]; 
				long double multi = a[pj][k] * pivot;
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
			std::cerr << "gausswp, Failed! the system is singular." << std::endl;
			return boost::none;
		}

		// backward substitution
		x[DiffData::AVECSIZE - 1] = b[p[DiffData::AVECSIZE - 1]] /
									a[p[DiffData::AVECSIZE - 1]][DiffData::AVECSIZE - 1];
		for (int k = DiffData::AVECSIZE - 2; k >= 0; k--) {
			pk = p[k];
			long double sum = b[pk];

			for (unsigned int i = k + 1; i < DiffData::AVECSIZE; i++)
				sum -= a[pk][i] * x[i];

				x[k] = sum / a[pk][k];
		}
	    
		return boost::optional<const array<long double, DiffData::AVECSIZE> >(x);
	}

	bool Diff::solve_diff_equ()
	{
#ifdef _OPENMP
		if (pdata_->ompthread_) {
			#pragma omp parallel sections
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
		} else {
			solve_diff_equ_O();
			solve_diff_equ_I();
		}

		return true;
#else
		BOOST_STATIC_ASSERT(false);

		return false;
#endif
	}

	/****************************************************
		the usual radial differential equation without 
        relativistic corrections 

        dM/dx = -(2NL + 1)M + 2r^2(V - ep)L
        dL/dx = M
	*****************************************************/
	long double dM_dx_sch(long double x, long double L, long double M,
						  const shared_ptr<DiffData> & pdiffdata)
	{
		const long double r = std::exp(x);

		return - (2.0 * static_cast<const long double>(pdiffdata->pdata_->l) + 1.0) * M +
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
	long double dM_dx_sdirac(long double x, long double L, long double M,
							 const shared_ptr<DiffData> & pdiffdata)
	{
		const long double r = std::exp(x);

		const long double mass = 1.0 + Data::al2half * (pdiffdata->E_ - fnc_V(r, pdiffdata));
		const long double d = Data::al2half * r / mass * dV_dr(r, pdiffdata);
		const long double l = static_cast<const long double>(pdiffdata->pdata_->l);

		// scaler treatment
		const long double d1 = - (2.0 * l + 1.0 + d) * M;
		const long double d2 = (2.0 * sqr(r) * mass * (fnc_V(r, pdiffdata) - pdiffdata->E_) -
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
	long double dM_dx_dirac(long double x, long double L, long double M,
							const shared_ptr<DiffData> & pdiffdata)
	{
		const long double r = std::exp(x);

		const long double mass = 1.0 + Data::al2half * (pdiffdata->E_ - fnc_V(r, pdiffdata));
		const long double d = Data::al2half * r / mass * dV_dr(r, pdiffdata);
		const shared_ptr<const Data> & pdata(pdiffdata->pdata_);
		const long double l = static_cast<const long double>(pdata->l);

		// dependence on all angular momentum
		const long double d1 = - (2.0 * l + 1.0 + d) * M;
		const long double d2 = (2.0 * sqr(r) * mass * (fnc_V(r, pdiffdata) - pdiffdata->E_) -
								d * (l + 1.0 + pdata->kappa)) * L;
		return d1 + d2;
	}

	template <typename T>
	T pow(T x, unsigned int n)
	{
		T p = x, y = 1.0;

		while (true) {
			if (n & 1)
				y *= p;

			n >>= 1;

			if (!n)
				return y;

			p *= p;
		}
	}
}
