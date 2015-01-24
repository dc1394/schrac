#ifndef _CI_STRING_H_
#define _CI_STRING_H

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <string>
#include <cstring>
#include <cctype>

#ifndef _MSC_VER
	#include <utility>
	#include <functional>
#endif

namespace HydroSchDirac {
#if !defined(_MSC_VER) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
	// Should I have just derived this class from
	// binary_function instead of equal_to?
	struct ci_equal_to : public std::equal_to<char> {
		bool operator()(char left, char right) {
			return std::toupper(left) == std::toupper(right);
		}
	};
#endif

	struct ci_char_traits : public std::char_traits<char> {
		static bool eq(char lhs, char rhs) {
			return std::toupper(lhs) == std::toupper(rhs);
		}
		static bool lt(char lhs, char rhs) {
			return std::toupper(lhs) < std::toupper(rhs);
		}
		static int compare(const char * const s1, const char * const s2, std::size_t n);
	};

	typedef std::basic_string<char, ci_char_traits> ci_string;

	std::ostream & operator<<(std::ostream & os, const ci_string & s);
}

#endif	// _CISTRING_H_
