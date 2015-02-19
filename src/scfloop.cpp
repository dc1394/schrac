/*! \file scfloop.cpp
    \brief SCFを行うクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "eigenvaluesearch.h"
#include "normalization.h"
#include "readinputfile.h"
#include "scfloop.h"
#include "simpson.h"
#include <iostream>         // for std::cout
#include <boost/math/constants/constants.hpp>

namespace schrac {
    // #region コンストラクタ

    ScfLoop::ScfLoop(std::pair<std::string, bool> const & arg) :
        PData([this]{ return pdata_; }, nullptr),
        PDiffData([this]{ return pdiffdata_; }, nullptr)
    {
        ReadInputFile rif(arg);			// ファイルを読み込む
        rif.readFile();
        pdata_ = rif.PData;

        initialize();
        
        if (pdata_->chemical_symbol_ == Data::Chemical_Symbol[0]) {
            pdiffsolver_ = std::make_shared<DiffSolver>(pdata_, pdiffdata_);
        }
        else {
            prho_ = std::make_shared<Rho>(pdiffdata_);
            pvh_ = std::make_shared<Vhartree>(pdiffdata_->r_mesh_);
            pdiffsolver_ = std::make_shared<DiffSolver>(pdata_, pdiffdata_, prho_, pvh_);
        }
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数

    void ScfLoop::message() const
    {
        std::cout << pdata_->chemical_symbol_
            << "原子の"
            << pdata_->orbital_
            << "軌道";

        if (pdata_->eq_type_ == Data::Eq_type::DIRAC && pdata_->spin_orbital_ == Data::ALPHA) {
            std::cout << "スピン上向きの";
        }
        else if (pdata_->eq_type_ == Data::Eq_type::DIRAC && pdata_->spin_orbital_ == Data::BETA) {
            std::cout << "スピン下向きの";
        }

        std::cout << "の波動関数と固有値を計算します。\n" << std::endl;
    }

    void ScfLoop::operator()()
    {
        if (pdata_->chemical_symbol_ != Data::Chemical_Symbol[0]) {
            make_vhartree();
        }

        EigenValueSearch evs(pdata_, pdiffdata_, prho_, pvh_);

        //cp.checkpoint("入力ファイル読み込み処理", __LINE__);

        if (!evs.search()) {
            std::cerr << "固有値が見つかりませんでした。終了します。" << std::endl;
            //goexit();
            //return EXIT_FAILURE;
        }

        auto const pdsol = nomalization(evs.PDiffSolver);
        auto const newrho = req_newrho(pdsol.at("2 Eigen function"));
        auto const res = req_normrd(newrho, prho_->PRho);

        std::cout << "NormRD = " << res;
        std::cout << ", Energy = " << req_energy(pdiffdata_->E_, newrho, pvh_->Vhart) << std::endl;
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数
    
    void ScfLoop::initialize()
    {        
        pdiffdata_ = std::make_shared<DiffData>(pdata_);

        pdiffdata_->r_mesh_.reserve(pdata_->grid_num_ + 1);
        for (auto i = 0; i <= pdata_->grid_num_; i++) {
            pdiffdata_->r_mesh_.push_back(std::exp(pdata_->xmin_ + static_cast<double>(i) * pdiffdata_->dx_));
        }
    }

    void ScfLoop::make_vhartree()
    {
        pdiffsolver_->solve_poisson();
        pvh_->set_vhartree_boundary_condition(pdata_->Z_);
        pvh_->vhart_init();
    }

    double ScfLoop::req_energy(double eigen, dvector const & rho, dvector const & vhartree) const
    {
        dvector u2;
        u2.reserve(pdata_->grid_num_ + 1);
        for (auto i = 0; i <= pdata_->grid_num_; i++) {
            u2.push_back(sqr(pdiffdata_->r_mesh_[i]) * rho[i]);
        }

        Simpson simpson(pdiffdata_->dx_);
        return 2.0 * eigen - simpson(vhartree, u2, pdiffdata_->r_mesh_, 1);
    }

    double ScfLoop::req_normrd(dvector const & newrho, dvector const & oldrho) const
    {
        BOOST_ASSERT(newrho.size() == oldrho.size());

        dvector residual;
        residual.reserve(newrho.size());

        for (auto i = 0; i <= pdata_->grid_num_; i++) {
            residual.push_back(newrho[i] - oldrho[i]);
        }

        Simpson const simpson(pdiffdata_->dx_);
        return 4.0 * boost::math::constants::pi<double>() *
            simpson(residual, residual, pdiffdata_->r_mesh_, 3);
    }

    dvector ScfLoop::req_newrho(dvector const & rf) const
    {
        dvector newrho;
        newrho.reserve(pdata_->grid_num_ + 1);
        for (auto i = 0; i <= pdata_->grid_num_; i++) {
            newrho.push_back(sqr(rf[i]));
        }

        return std::move(newrho);
    }

    // #endregion privateメンバ関数
}
