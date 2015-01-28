#include "EigenValueSearch.h"

namespace schrac {
	class WF_Normalize
		: private boost::noncopyable {
	public:
		typedef tuple<const ldvector, const ldvector, const ldvector> d3tup;

	private:
		const shared_ptr<const Data> pdata_;
		const shared_ptr<const DiffData> pdiffdata_;

		ldvector RV;
		ldvector XV;
		ldvector RF;
		ldvector PF;

		long double phi(std::size_t n) const;
		void WF_coalesce(const shared_ptr<Diff> & pdiff);
		void WF_coalesce_omp(const shared_ptr<Diff> & pdiff);
		long double simpson() const;
		long double simpson_omp() const;

	public:
		explicit WF_Normalize(const shared_ptr<Diff> & pdiff);
		void operator()();
		const WF_Normalize::d3tup getptup() const
		{ return make_tuple(RV, RF, PF); }
	};
}

