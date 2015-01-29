#ifndef _DIFF_H_
#define _DIFF_H_

#pragma once

#include "DiffData.h"

namespace schrac {
	double fnc_V(double r, const shared_ptr<DiffData> & pdiffdata);
	double dV_dr(double r, const shared_ptr<DiffData> & pdiffdata);
	double dL_dx(double M);
	double dM_dx_sch(double x, double L, double M,
					 const shared_ptr<DiffData> & pdiffdata);
	double dM_dx_sdirac(double x, double L, double M,
						const shared_ptr<DiffData> & pdiffdata);
	double dM_dx_dirac(double x, double L, double M,
					   const shared_ptr<DiffData> & pdiffdata);
    
    class Diff final {
			static const double MINV;

			bool a_init();
			void b_init();

			void init_LM_O();
			void init_LM_I();

			const boost::optional<const array<double, DiffData::AVECSIZE> >
				S_gausswp(array<array<double, DiffData::AVECSIZE>, DiffData::AVECSIZE> & a,
						  array<double, DiffData::AVECSIZE> & b) const;

			virtual bool solve_diff_equ_O() = 0;
			virtual bool solve_diff_equ_I() = 0;

	protected:
		const shared_ptr<const Data> pdata_;
		const shared_ptr<DiffData> pdiffdata_;

		Diff(const shared_ptr<const Data> & pdata, double E, double TINY);

	public:
		static function<double(double, double, double,
			const shared_ptr<DiffData> &)> dM_dx;
		typedef tuple<const array<double, 2>, const array<double, 2> > mytuple;

		void Initialize(double E);
		const shared_ptr<DiffData> & getpDiffData() const
		{ return pdiffdata_; }
		virtual bool solve_diff_equ();
		virtual ~Diff() {}
		const mytuple getMPval() const;
	};

	inline double fnc_V(double r, const shared_ptr<DiffData> & pdiffdata)
	{
		return - pdiffdata->Z / r;
	}

	inline double dV_dr(double r, const shared_ptr<DiffData> & pdiffdata)
	{
		return pdiffdata->Z / (r * r);
	}
		
	inline double dL_dx(double M)
	{
		return M;
	}

	inline const Diff::mytuple Diff::getMPval() const
	{
		array<double, 2> L, M;

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
