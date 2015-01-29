#ifndef _BULIRSCH_STOER_H_
#define _BULIRSCH_STOER_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include "Adaptive_Step.h"

namespace schrac {
	int div2(int n);

	class Bulirsch_Stoer : public Adaptive_Step {
		static const int THRESHOLD = 3;

		static const double SHRINK;
		static const double GROW;

		static const array<const int, AdapStepHelper::IMAX> nseq;

		virtual bool odeint(const shared_ptr<AdapStepHelper> & pasa) const;

		bool bsstep(double htry, const AdapStepHelper::darray & yscal,
					const shared_ptr<AdapStepHelper> & pasa) const;
		void rzextr(std::size_t iest, double xest,
					const AdapStepHelper::darray & yest,
					AdapStepHelper::darray & dy,
					const shared_ptr<AdapStepHelper> & pasa) const;
		void mmid(double xs, double htot,
				  int nstep,
				  AdapStepHelper::darray & y,
				  AdapStepHelper::darray & dydx,
				  AdapStepHelper::darray & yout,
				  const shared_ptr<AdapStepHelper> & pasa) const;

	public:
		Bulirsch_Stoer(const shared_ptr<const Data> & pdata,
					   double E, double TINY)
		 :	Adaptive_Step(pdata, E, TINY) {}
		virtual ~Bulirsch_Stoer() {}
	};

	inline int div2(int n)
	{
		n += (n >> 31);

		return n >> 1;
	}
}

#endif	// _BULIRSCH_STOER_H_
