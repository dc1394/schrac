#include "Diff.h"

#ifdef _DEBUG
namespace schrac {
	class ModEuler : public Diff {
		virtual bool solve_diff_equ_O();
		virtual bool solve_diff_equ_I();

	public:
		ModEuler(const shared_ptr<const Data> & pdata, long double E, long double TINY)
		 :	Diff(pdata, E, TINY)	{}
		virtual ~ModEuler() {}
	};
}
#endif
