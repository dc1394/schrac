#include "Diff.h"

namespace HydroSchDirac {
	class RungeKutta : public Diff {
		virtual bool solve_diff_equ_O();
		virtual bool solve_diff_equ_I();

	public:
		RungeKutta(const shared_ptr<const Data> & pdata,
				   long double E, long double TINY) : Diff(pdata, E, TINY) {}
		virtual ~RungeKutta() {}
	};
}
