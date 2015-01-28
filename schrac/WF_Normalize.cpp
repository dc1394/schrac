#include "WF_Normalize.h"

namespace schrac {
	// constructor
	WF_Normalize::WF_Normalize(const shared_ptr<Diff> & pdiff)
	 :	pdata_(pdiff->getpDiffData()->pdata_),
		pdiffdata_(pdiff->getpDiffData()),
		RV(pdata_->grid_num, 0.0), XV(pdata_->grid_num, 0.0),
		RF(pdata_->grid_num, 0.0), PF(pdata_->grid_num, 0.0)
	{
		if (pdata_->ompthread_)
			WF_coalesce_omp(pdiff);
		else
			WF_coalesce(pdiff);
	}

	/*	private method	*/
	
	void WF_Normalize::WF_coalesce(const shared_ptr<Diff> & pdiff)
	{
		const std::size_t GRID_NUM = pdata_->grid_num;
		const int MP_O = boost::numeric_cast<const int>(pdiffdata_->MP_O);
		const int MP_I = boost::numeric_cast<const int>(pdiffdata_->MP_I - 1);

		const ldvector & LI(pdiffdata_->LI);
		const ldvector & LO(pdiffdata_->LO);	
		const ldvector & MI(pdiffdata_->MI);
		const ldvector & MO(pdiffdata_->MO);
		const ldvector & XV_I(pdiffdata_->XV_I);
		const ldvector & XV_O(pdiffdata_->XV_O);
		const ldvector & RV_I(pdiffdata_->RV_I);
		const ldvector & RV_O(pdiffdata_->RV_O);
		const ldvector & VP_I(pdiffdata_->VP_I);

		const Diff::mytuple tup(pdiff->getMPval());
		const long double ratio = (get<0>(tup))[0] / (get<0>(tup))[1];

		RV.assign(RV_O.begin(), RV_O.end());
		XV.assign(XV_O.begin(), XV_O.end());
		for (int i = 0; i <= MP_O; i++) {
			RF[i] = schrac::pow(RV[i], pdata_->l) * LO[i];
			PF[i] = RV[i] * RF[i];
		}

		RV.resize(GRID_NUM, 0.0);
		XV.resize(GRID_NUM, 0.0);

		for (int i = MP_I; i >= 0; i--) {
			const int j = MP_O + MP_I - i + 1;
			RV[j] = RV_I[i];
			XV[j] = XV_I[i];
			RF[j] = schrac::pow(RV[j], pdata_->l) * ratio * LI[i];
			PF[j] = RV[j] * RF[j];
			
		}
	}

	void WF_Normalize::WF_coalesce_omp(const shared_ptr<Diff> & pdiff)
	{
#ifdef _OPENMP
		const int MP_O = boost::numeric_cast<const int>(pdiffdata_->MP_O);
		const int MP_I = boost::numeric_cast<const int>(pdiffdata_->MP_I - 1);

		const ldvector & LI(pdiffdata_->LI);
		const ldvector & LO(pdiffdata_->LO);	
		const ldvector & MI(pdiffdata_->MI);
		const ldvector & MO(pdiffdata_->MO);
		const ldvector & XV_I(pdiffdata_->XV_I);
		const ldvector & XV_O(pdiffdata_->XV_O);
		const ldvector & RV_I(pdiffdata_->RV_I);
		const ldvector & RV_O(pdiffdata_->RV_O);
		const ldvector & VP_I(pdiffdata_->VP_I);

		const Diff::mytuple tup(pdiff->getMPval());
		const long double ratio = (get<0>(tup))[0] / (get<0>(tup))[1];
		
		#pragma omp parallel
		{
			#pragma omp for nowait
			for (int i = 0; i <= MP_O; i++) {
				RV[i] = RV_O[i];
				XV[i] = XV_O[i];
				RF[i] = schrac::pow(RV[i], pdata_->l) * LO[i];
				PF[i] = RV[i] * RF[i];
			}

			#pragma omp for nowait
			for (int i = MP_I; i >= 0; i--) {
				const int j = MP_O + MP_I - i + 1;
				RV[j] = RV_I[i];
				XV[j] = XV_I[i];
				RF[j] = schrac::pow(RV[j], pdata_->l) * ratio * LI[i];
				PF[j] = RV[j] * RF[j];
			}
		}
#else
		BOOST_STATIC_ASSERT(false);
#endif
	}

	long double WF_Normalize::simpson() const
	{
		long double sum = 0.0;
		const int max = boost::numeric_cast<const int>(pdata_->grid_num - 2);

		for (int i = 0; i < max; i += 2) {
			const long double f0 = PF[i] * PF[i] * RV[i];
			const long double f1 = PF[i + 1] * PF[i + 1] * RV[i + 1];
			const long double f2 = PF[i + 2] * PF[i + 2] * RV[i + 2];
			sum += (f0 + 4.0 * f1 + f2);
		}

		return sum * pdiffdata_->DX / 3.0;
	}

	long double WF_Normalize::simpson_omp() const
	{
#ifdef _OPENMP
		volatile long double sum = 0.0;

		const int max = boost::numeric_cast<const int>(pdata_->grid_num - 2);

		#pragma omp parallel for reduction(+:sum)
		for (int i = 0; i < max; i += 2) {
			const long double f0 = PF[i] * PF[i] * RV[i];
			const long double f1 = PF[i + 1] * PF[i + 1] * RV[i + 1];
			const long double f2 = PF[i + 2] * PF[i + 2] * RV[i + 2];
			sum += (f0 + 4.0 * f1 + f2);
		}

		return sum * pdiffdata_->DX / 3.0;
#else
		BOOST_STATIC_ASSERT(false);

		return 0.0;
#endif
	}

	/*	public method	*/

	void WF_Normalize::operator()()
	{
#ifdef _OPENMP
		const int grid_num = boost::numeric_cast<const int>(pdata_->grid_num);

		if (pdata_->ompthread_) {
			const long double a = 1.0 / std::sqrt(simpson_omp());
			#pragma omp parallel for
			for (int i = 0; i < grid_num; i++) {
				RF[i] *= a;
				PF[i] *= a;
			}
		} else {
			const long double a = 1.0 / std::sqrt(simpson());
			for (int i = 0; i < grid_num; i++) {
				RF[i] *= a;
				PF[i] *= a;
			}
		}
#else
		BOOST_STATIC_ASSERT(false);
#endif
	}
}
