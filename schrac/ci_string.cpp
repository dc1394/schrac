#include "ci_string.h"
#include <cstring>

#ifndef _MSC_VER
    #include <algorithm>    // std::mismatch
#endif

namespace schrac {
    //! A public static member function.
    /*!
        2個のバッファー (大文字と小文字を区別しない) を比較する関数の実装
        \param s1 一つ目の文字列
        \param s2 二つ目の文字列
        \param 文字列の長さ
        \return std::strcmpの戻り値と同様
    */
	int ci_char_traits::compare(char const * s1, char const * s2, std::size_t n)
	{
#ifdef _MSC_VER
		return _memicmp(s1, s2, n);
#else
		return memIcmp(s1, s2, n);
#endif
	}

#ifndef _MSC_VER
    //! A function.
    /*!
        2個のバッファー (大文字と小文字を区別しない) を比較する関数の実装
        \param s1 一つ目の文字列
        \param s2 二つ目の文字列
        \param 文字列の長さ
        \return std::strcmpの戻り値と同様
    */
    int memIcmp(char const * s1, char const * s2, size_t n)
    {
        typedef std::pair<char const*, char const*> ci_diff_pair;

        auto p = std::mismatch(p1, p1 + n, p2, [](char left, char right) {
            return std::toupper(left) == std::toupper(right);
        });

        // both characters match exactly (case insensitive)
        if (p.first == p1 + n && p.second == p2 + n) {
            return 0;
        }

        return *(p.first) < *(p.second) ? -1 : 1;
    }
#endif

    //! A function.
    /*!
        ci_stringに対するoperator<<の実装
        \param os 対象のstd::ostream
        \param s 対象の（大文字小文字を区別しない）文字列
        \return 結果の参照
    */
	std::ostream & operator<<(std::ostream & os, const ci_string & s)
	{
		os << s.c_str();

		return os;
	}
}
