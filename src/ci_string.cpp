/*! \file ci_string.cpp
    \brief 大文字小文字を区別しない文字列クラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD-2 License.
    This software is released under the BSD-2 License.
*/    This software is released under the BSD-2 License.*/

#include "ci_string.h"
#include <cstring>          // for _memicmp, std::toupper

#ifndef _MSC_VER
    #include <algorithm>    // std::mismatch
#endif

namespace schrac {
    // #region staticメンバ関数

    int ci_char_traits::compare(char const * s1, char const * s2, std::size_t n)
    {
#ifdef _MSC_VER
        return _memicmp(s1, s2, n);
#else
        return memIcmp(s1, s2, n);
#endif
    }
    
    // #region staticメンバ関数

    // #region 非メンバ関数

#ifndef _MSC_VER
    int memIcmp(char const * s1, char const * s2, std::size_t n)
    {
        using ci_diff_pair = std::pair<char const *, char const *>;

        auto const p = std::mismatch(s1, s1 + n, s2, [](char left, char right) {
            return std::toupper(left) == std::toupper(right);
        });

        // both characters match exactly (case insensitive)
        if (p.first == s1 + n && p.second == s2 + n) {
            return 0;
        }

        return *(p.first) < *(p.second) ? -1 : 1;
    }
#endif

    std::ostream & operator<<(std::ostream & os, const ci_string & s)
    {
        os << s.c_str();

        return os;
    }

    // #endregion 非メンバ関数
}
