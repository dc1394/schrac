#ifndef _ADAPTIVE_STEP_H_
#define _ADAPTIVE_STEP_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include "AdapStepHelper.h"

namespace HydroSchDirac {
	class Adaptive_Step : public Diff {
		virtual bool odeint(const shared_ptr<AdapStepHelper> & pasa) const = 0;
		virtual bool solve_diff_equ_O();
		virtual bool solve_diff_equ_I();

	protected:
		
		void RungeKutta();
		bool rkqc(long double htry,
				  const AdapStepHelper::darray & yscal,
				  const shared_ptr<AdapStepHelper> & pasa) const;
		void rk4(long double x, long double h,
				 const AdapStepHelper::darray & y,
				 const AdapStepHelper::darray & dydx,
				 AdapStepHelper::darray & yout) const;
		void derivs(long double x,
					const AdapStepHelper::darray & y,
					AdapStepHelper::darray & dydx) const;
	public:
		Adaptive_Step(const shared_ptr<const Data> & pdata,
					  long double E, long double TINY)
		 :	Diff(pdata, E, TINY) {}
		virtual ~Adaptive_Step() {}
		virtual bool solve_diff_equ();
	};

	inline void Adaptive_Step::derivs(long double x,
									  const AdapStepHelper::darray & y,
									  AdapStepHelper::darray & dydx) const
	{
		dydx[0] = dL_dx(y[1]);
		dydx[1] = Diff::dM_dx(x, y[0], y[1], pdiffdata_);
	}
}

#endif	// _ADAPTIVE_STEP_H_
