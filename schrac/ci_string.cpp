#include "ci_string.h"

namespace HydroSchDirac {
#ifndef _MSC_VER
	// Assumption: assume that both p1 and p2 are of size n
	int memIcmp(char const* p1, char const* p2, size_t n)
	{
		typedef std::pair<char const*, char const*> ci_diff_pair;
#ifdef __GXX_EXPERIMENTAL_CXX0X__
		ci_diff_pair p = std::mismatch(p1, p1 + n, p2, [](char left, char right) {
			return std::toupper(left) == std::toupper(right);
		});
#else
		ci_diff_pair p = std::mismatch(p1, p1 + n, p2, ci_equal_to());
#endif
		// both characters match exactly (case insensitive)
		if (p.first == p1 + n && p.second == p2 + n) {
			return 0;
		}

		return *(p.first) < *(p.second) ? -1 : 1;
	}
#endif

	int ci_char_traits::compare(const char * const s1, const char * const s2, std::size_t n)
	{
#ifdef _MSC_VER
		return _memicmp(s1, s2, n);
#else
		return memIcmp(s1, s2, n);
#endif
	}

	std::ostream & operator<<(std::ostream & os, const ci_string & s)
	{
		os << s.c_str();

		return os;
	}

}