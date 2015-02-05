#include "simpson.h"
#include <cstdint>          // for std::int32_t
#include <boost/cast.hpp>   // for boost::numeric_cast

namespace schrac {
    double Simpson::operator()(Simpson::dvector const & f, Simpson::dvector const & r) const
    {
        auto sum = 0.0;
        auto const max = boost::numeric_cast<std::int32_t>(f.size() - 2);
        for (auto i = 0; i < max; i += 2) {
            auto const f0 = f[i] * f[i] * r[i];
            auto const f1 = f[i + 1] * f[i + 1] * r[i + 1];
            auto const f2 = f[i + 2] * f[i + 2] * r[i + 2];
            sum += (f0 + 4.0 * f1 + f2);
        }

        return sum * dx_ / 3.0;
    }

    double Simpson::operator()(Simpson::dvector const & f, Simpson::dvector const & g, Simpson::dvector const & r, std::int32_t n) const
    {
        auto sum = 0.0;
        auto const max = boost::numeric_cast<std::int32_t>(f.size() - 2);
        for (auto i = 0; i < max; i += 2) {
            auto const f0 = f[i] * g[i] * schrac::pow(r[i], n);
            auto const f1 = f[i + 1] * g[i + 1] * schrac::pow(r[i + 1], n);
            auto const f2 = f[i + 2] * g[i + 2] * schrac::pow(r[i + 2], n);
            sum += (f0 + 4.0 * f1 + f2);
        }

        return sum * dx_ / 3.0;
    }
}
