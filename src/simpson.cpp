#include "simpson.h"
#include <boost/cast.hpp>   // for boost::numeric_cast

namespace schrac {
    double Simpson::operator()(Simpson::dvector const & f, Simpson::dvector const & r) const
    {
        return (*this)(f, f, r, 1);
    }

    double Simpson::operator()(Simpson::dvector const & f, Simpson::dvector const & g, Simpson::dvector const & r, std::int32_t n) const
    {
        auto sum = 0.0;
        auto const max = boost::numeric_cast<std::int32_t>(f.size() - 2);
        for (auto i = 0; i < max; i += 2) {
            auto const f0 = f[i] * g[i] * std::pow(r[i], n);
            auto const f1 = f[i + 1] * g[i + 1] * std::pow(r[i + 1], n);
            auto const f2 = f[i + 2] * g[i + 2] * std::pow(r[i + 2], n);
            sum += (f0 + 4.0 * f1 + f2);
        }
        
        return sum * dx_ / 3.0;
    }
}
