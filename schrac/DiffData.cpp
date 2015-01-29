#include "DiffData.h"

namespace schrac {
	// constructor
	DiffData::DiffData(const shared_ptr<const Data> & pdata, double E, double TINY)
	 :	pdata_(pdata), node(pdata_->n - pdata_->l - 1),
		thisnode(0), TINY_(TINY), E_(E)
	{
#ifdef _OPENMP
		Z = static_cast<const double>(pdata_->Z);

		const std::size_t GRID_NUM = pdata_->grid_num;

		MP_O = round(static_cast<const double>(GRID_NUM - 1) * pdata_->mat_po_ratio);
		MP_I = GRID_NUM - MP_O - 1;
		
		OSIZE = boost::numeric_cast<const int>(MP_O + 1);
		ISIZE = boost::numeric_cast<const int>(MP_I + 1);

		DX = (pdata_->xmax - pdata_->xmin) / static_cast<const double>(GRID_NUM - 1);

		// ÉÅÉÇÉäämï€
		XV_O.resize(OSIZE);
		XV_I.resize(ISIZE);
		RV_O.resize(OSIZE);
		RV_I.resize(ISIZE);
		VP_O.resize(OSIZE);
		VP_I.resize(ISIZE);
		LO.resize(OSIZE);
		LI.resize(ISIZE);
		MO.resize(OSIZE);
		MI.resize(ISIZE);

		const int len = boost::numeric_cast<const int>(GRID_NUM - ISIZE);
		if (pdata_->ompthread_) {
			#pragma omp parallel
			{
				#pragma omp for nowait
				for (int i = 0; i < OSIZE; i++) {
					const double x = pdata_->xmin + static_cast<const double>(i) * DX;
					XV_O[i] = x;
					RV_O[i] = std::exp(x);
					VP_O[i] = fnc_V(x);
				}
				#pragma omp for nowait
				for (int i = boost::numeric_cast<const int>(GRID_NUM - 1); i >= len; i--) {
					const double x = pdata_->xmin + static_cast<const double>(i) * DX;
					XV_I[GRID_NUM - 1 - i] = x;
					RV_I[GRID_NUM - 1 - i] = std::exp(x);
					VP_I[GRID_NUM - 1 - i] = fnc_V(x);
				}
			}
		} else {
			for (int i = 0; i < OSIZE; i++) {
				const double x = pdata_->xmin + static_cast<const double>(i) * DX;
				XV_O[i] = x;
				RV_O[i] = std::exp(x);
				VP_O[i] = fnc_V(x);
			}
			for (int i = boost::numeric_cast<const int>(GRID_NUM - 1); i >= len; i--) {
				const double x = pdata_->xmin + static_cast<const double>(i) * DX;
				XV_I[GRID_NUM - 1 - i] = x;
				RV_I[GRID_NUM - 1 - i] = std::exp(x);
				VP_I[GRID_NUM - 1 - i] = fnc_V(x);
			}
		}
#else
		BOOST_STATIC_ASSERT(false);
#endif
	}
}
