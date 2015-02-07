/*! \file kinetic_energy.cpp
    \brief 運動エネルギーを計算するクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "kinetic_energy.h"
#include "simpson.h"
#include <boost/math/constants/constants.hpp>   // for boost::math::constants::pi

namespace schrac {
    double Kinetic_Energy::operator()(Kinetic_Energy::dvector const & L, Kinetic_Energy::dvector const & M, Kinetic_Energy::dvector const & r) const
    {
        Simpson simpson(pdiffdata_->dx_);

        auto const pi = boost::math::constants::pi<double>();

        double xxx = -4.0 * pi * simpson(L, M, r, 3);
        xxx += 4.0 * pi * simpson(L, M, r, 2);
        xxx += 2.0 * pi * simpson(M, M, r, 3);

        return xxx;
    }
}
