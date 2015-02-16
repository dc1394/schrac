/*! \file scfloop.cpp
    \brief SCFを行うクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "eigenvaluesearch.h"
#include "readinputfile.h"
#include "scfloop.h"
#include <iostream>         // for std::cout

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
        make_rho();
        
        prho_ = std::make_shared<Rho>(pdiffdata_->r_mesh_, pdiffdata_->rho_);
        pvh_ = std::make_shared<Vhartree>(pdiffdata_->r_mesh_);
        pdiffsolver_ = std::make_shared<DiffSolver>(pdata_, pdiffdata_, prho_, pvh_);
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
        make_vhartree();
        EigenValueSearch evs(pdata_, pdiffdata_, prho_, pvh_);

        //cp.checkpoint("入力ファイル読み込み処理", __LINE__);

        if (!evs.search()) {
            std::cerr << "固有値が見つかりませんでした。終了します。" << std::endl;
            //goexit();
            //return EXIT_FAILURE;
        }

        //auto const pdsol = nomalization(evs.PDiffSolver);
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

    void ScfLoop::make_rho()
    {
        if (pdata_->chemical_symbol_ == Data::Chemical_Symbol[0])
        {
            pdiffdata_->rho_.resize(pdata_->grid_num_ + 1, 0.0);
            return;
        } else if (pdata_->rho0_c_ == boost::none || pdata_->rho0_alpha_ == boost::none) {
            // デフォルト値を代入
            auto const w = 50.0;
            pdata_->rho0_c_ = boost::optional<double>(std::pow(w, 4) / 16.0);
            pdata_->rho0_alpha_ = boost::optional<double>(0.5 * w);
            *pdata_->rho0_c_ *= pdata_->Z_ / w;
        }

        pdiffdata_->rho_.reserve(pdata_->grid_num_ + 1);
        for (auto i = 0; i <= pdata_->grid_num_; i++) {
            pdiffdata_->rho_.push_back(*pdata_->rho0_c_ * std::exp(- *pdata_->rho0_alpha_ * pdiffdata_->r_mesh_[i]));
        }
    }

    void ScfLoop::make_vhartree()
    {
        pdiffsolver_->solve_poisson();
        int i = 1;
    }

    // #endregion privateメンバ関数

}
