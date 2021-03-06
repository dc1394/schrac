﻿/*! \file normalize.cpp
    \brief 実際に波動関数を正規化する関数の実装

    Copyright © 2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "normalization.h"
#include "schnormalize.h"
#include "diracnormalize.h"

namespace schrac {
    boost::container::flat_map<std::string, dvector> nomalization(std::shared_ptr<DiffSolver> const & pdiffsolver)
    {
        switch (pdiffsolver->pdata_->eq_type_) {
        case Data::Eq_type::DIRAC:
        case Data::Eq_type::SDIRAC:
        {
            DiracNormalize dn(pdiffsolver);
            dn.evaluate();
            return dn.getresult();
        }
            break;

        case Data::Eq_type::SCH:
        {
            SchNormalize sn(pdiffsolver);
            sn.evaluate();
            return sn.getresult();
        }
            break;

        default:
            BOOST_ASSERT(!"何かがおかしい！");
            return boost::container::flat_map<std::string, dvector>();
            break;
        }
    }
}
