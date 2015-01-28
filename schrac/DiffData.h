#ifndef _DIFFDATA_H_
#define _DIFFDATA_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <memory>

#include "data.h"

namespace schrac {
	std::size_t round(long double val);

	struct DiffData {
		static const std::size_t AVECSIZE = 3;
		static const std::size_t BVECSIZE = 5;

		const std::shared_ptr<const Data> pdata_;

		const int node;
		int thisnode;

		std::size_t MP_O;
		std::size_t MP_I;
		int OSIZE;
		int ISIZE;

		const long double TINY_;

		long double Z;
		long double DX;
		long double E_;


		array<long double, AVECSIZE> V_A;
		array<long double, BVECSIZE> V_B;

		ldvector XV_I;
		ldvector XV_O;
		ldvector RV_I;
		ldvector RV_O;
		ldvector VP_I;
		ldvector VP_O;
		ldvector LO;
		ldvector MO;
		ldvector LI;
		ldvector MI;

		DiffData(const shared_ptr<const data> & pdata, long double E, long double TINY);
		long double fnc_V(long double x) const;
		void node_count(int i, const ldvector & WF);
	};

	inline long double DiffData::fnc_V(long double x) const
	{
		return - Z * std::exp(- x);
	}

	inline void DiffData::node_count(int i, const ldvector & WF)
	{
		if (WF[i] * WF[i - 1] < 0.0)
			thisnode++;
	}

	inline std::size_t round(long double val)
	{
		std::ostringstream oss;

		// ***** ŠÛ‚ß‚ðs‚¢•¶Žš—ñ‚É•ÏŠ· *****
		oss << boost::format("%.0f") % val;

		// ***** •¶Žš—ñ‚©‚ç”’l‚ÉÄ•ÏŠ·‚µ‚Ä•Ô‚· *****
		return boost::lexical_cast<std::size_t>(oss.str());
	}
}

#endif	// _DIFFDATA_H_
