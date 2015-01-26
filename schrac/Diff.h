#ifndef _DIFF_H_
#define _DIFF_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include "DiffData.h"

namespace HydroSchDirac {
	long double fnc_V(long double r, const shared_ptr<DiffData> & pdiffdata);
	long double dV_dr(long double r, const shared_ptr<DiffData> & pdiffdata);
	long double dL_dx(long double M);
	long double dM_dx_sch(long double x, long double L, long double M,
					 const shared_ptr<DiffData> & pdiffdata);
	long double dM_dx_sdirac(long double x, long double L, long double M,
						const shared_ptr<DiffData> & pdiffdata);
	long double dM_dx_dirac(long double x, long double L, long double M,
					   const shared_ptr<DiffData> & pdiffdata);

	template <typename T> T pow(T x, unsigned int n);
	template <typename T> T sqr(T x);

	class Diff :
		private boost::noncopyable {
			static const long double MINV;

			bool a_init();
			void b_init();

			void init_LM_O();
			void init_LM_I();

			const boost::optional<const array<long double, DiffData::AVECSIZE> >
				S_gausswp(array<array<long double, DiffData::AVECSIZE>, DiffData::AVECSIZE> & a,
						  array<long double, DiffData::AVECSIZE> & b) const;

			virtual bool solve_diff_equ_O() = 0;
			virtual bool solve_diff_equ_I() = 0;

	protected:
		const shared_ptr<const Data> pdata_;
		const shared_ptr<DiffData> pdiffdata_;

		Diff(const shared_ptr<const Data> & pdata, long double E, long double TINY);

	public:
		static function<long double(long double, long double, long double,
			const shared_ptr<DiffData> &)> dM_dx;
		typedef tuple<const array<long double, 2>, const array<long double, 2> > mytuple;

		void Initialize(long double E);
		const shared_ptr<DiffData> & getpDiffData() const
		{ return pdiffdata_; }
		virtual bool solve_diff_equ();
		virtual ~Diff() {}
		const mytuple getMPval() const;
	};

	inline long double fnc_V(long double r, const shared_ptr<DiffData> & pdiffdata)
	{
		return - pdiffdata->Z / r;
	}

	inline long double dV_dr(long double r, const shared_ptr<DiffData> & pdiffdata)
	{
		return pdiffdata->Z / (r * r);
	}
		
	inline long double dL_dx(long double M)
	{
		return M;
	}

	inline const Diff::mytuple Diff::getMPval() const
	{
		array<long double, 2> L, M;

		L[0] = pdiffdata_->LO[pdiffdata_->MP_O];
		L[1] = pdiffdata_->LI[pdiffdata_->MP_I];
		M[0] = pdiffdata_->MO[pdiffdata_->MP_O];
		M[1] = pdiffdata_->MI[pdiffdata_->MP_I];

		return make_tuple(L, M);
	}
	
	template <typename T>
	inline T sqr(T x)
	{
		return x * x;
	}
}

#endif	// _DIFF_H_
