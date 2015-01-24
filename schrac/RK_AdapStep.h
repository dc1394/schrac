#ifndef _RK_AdapStep_H_
#define _RK_AdapStep_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include "Adaptive_Step.h"

namespace HydroSchDirac {
	class RK_AdapStep : public Adaptive_Step {
		virtual bool odeint(const shared_ptr<AdapStepHelper> & pasa) const;

	public:
		RK_AdapStep(const shared_ptr<const Data> & pdata,
				   long double E, long double TINY)
		 :	Adaptive_Step(pdata, E, TINY) {}
		virtual ~RK_AdapStep() {}
	};
}
#endif	// _RK_AdapStep_H_
