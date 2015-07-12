/*! \file normalization.h
    \brief 実際に波動関数を正規化する関数の宣言

    Copyright © 2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _NORMALIZATION_H_
#define _NORMALIZATION_H_

#include "diffsolver.h"
#include <boost/container/flat_map.hpp>

#pragma once

namespace schrac {
    //! A function.
    /*!
        波動関数の正規化を行う関数
        \param pdiffsolver 微分方程式のデータオブジェクト
        \return メッシュと波動関数が格納されたmap
    */
    boost::container::flat_map<std::string, dvector> nomalization(std::shared_ptr<DiffSolver> const & pdiffsolver);
}

#endif	// _NORMALIZATION_H_
