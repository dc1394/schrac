#include "Diff.h"

namespace HydroSchDirac {
	struct AdapStepHelper : 
		private boost::noncopyable {
			static const int OO = 0;
			static const int INFINITY = 1;

			static const int IMAX = 11;
			static const int NUSE = 7;

			static const long double PGROW;
			static const long double PSHRNK;
			static const long double FCOR;
			static const long double SAFETY;
			static const long double ERRCON;

			static const std::size_t NVAR = 2;
			typedef array<long double, NVAR> darray;

			const shared_ptr<DiffData> pdiffdata_;

			const unsigned int MAXSTEP;
			const int MP;

			const long double H1;
			const long double HMIN;
			const long double EPS;
			const long double TINY;

			const long double x1_;
			const long double x2_;

			const boost::optional<const ldvector &> xp_;
			array<boost::optional<ldvector &>, NVAR> yp_;

			array<array<long double, NUSE>, NVAR> d_;
			array<long double, IMAX> x_;
		
			int i_;
			
			darray y_;
			darray dydx_;

			long double xx_;
			long double hdid_;
			long double hnext_;

			AdapStepHelper(const shared_ptr<DiffData> & pdiffdata, const boost::mpl::int_<OO> &);
			AdapStepHelper(const shared_ptr<DiffData> & pdiffdata, const boost::mpl::int_<INFINITY> &);
	};

	inline AdapStepHelper::AdapStepHelper(const shared_ptr<DiffData> & pdiffdata,
										  const boost::mpl::int_<OO> &)
	 :	pdiffdata_(pdiffdata),
		MAXSTEP(static_cast<const unsigned int>(pdiffdata_->pdata_->grid_num) * 100),
		MP(boost::numeric_cast<const int>(pdiffdata_->MP_O)),
		H1(pdiffdata_->DX), HMIN(pdiffdata_->TINY_),
		EPS(pdiffdata_->pdata_->eps), TINY(pdiffdata_->TINY_),
		x1_(pdiffdata_->XV_O[DiffData::AVECSIZE - 1]),
		x2_(pdiffdata_->XV_O[pdiffdata_->MP_O]), 
		xp_(boost::optional<const ldvector &>(pdiffdata_->XV_O)),
		i_(DiffData::AVECSIZE), xx_(x1_)
	{
		yp_[0] = boost::optional<ldvector &>(pdiffdata_->LO);
		yp_[1] = boost::optional<ldvector &>(pdiffdata_->MO);

		y_[0] = pdiffdata_->LO[DiffData::AVECSIZE - 1];
		y_[1] = pdiffdata_->MO[DiffData::AVECSIZE - 1];
	}

	inline AdapStepHelper::AdapStepHelper(const shared_ptr<DiffData> & pdiffdata,
										  const boost::mpl::int_<INFINITY> &)
	 :	pdiffdata_(pdiffdata),
		MAXSTEP(static_cast<const unsigned int>(pdiffdata_->pdata_->grid_num) * 100),
		MP(boost::numeric_cast<const int>(pdiffdata_->MP_I)),
		H1(pdiffdata_->DX), HMIN(pdiffdata_->TINY_),
		EPS(pdiffdata_->pdata_->eps), TINY(pdiffdata_->TINY_),
		x1_(pdiffdata_->XV_I[0]),
		x2_(pdiffdata_->XV_I[pdiffdata_->MP_I]),
		xp_(boost::optional<const ldvector &>(pdiffdata_->XV_I)),
		i_(2), xx_(x1_)
	{
		yp_[0] = boost::optional<ldvector &>(pdiffdata_->LI);
		yp_[1] = boost::optional<ldvector &>(pdiffdata_->MI);

		y_[0] = pdiffdata_->LI[1];
		y_[1] = pdiffdata_->MI[1];
	}
}

