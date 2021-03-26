/*! \file schracmain.cpp
    \brief メインファイル

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "checkpoint/checkpoint.h"
#include "energy.h"
#include "getcomlineoption.h"
#include "goexit.h"
#include "normalization.h"
#include "scfloop.h"
#include "wavefunctionsave.h"
#include <cstdlib>                              // for EXIT_FAILURE, EXIT_SUCCESS
#include <iostream>                             // for std::cerr
#include <optional>								// for std::optional
#include <boost/format.hpp>                     // for boost::format

int main(int argc, char * argv[])
{
    using namespace schrac;
    checkpoint::CheckPoint cp;

    cp.checkpoint("処理開始", __LINE__);

    GetComLineOption mg;
    switch (mg.getopt(argc, argv)) {
    case -1:
        goexit();

        return EXIT_FAILURE;
        break;

    case 0:
        break;

    case 1:
        goexit();

        return EXIT_SUCCESS;
        break;

    default:
        BOOST_ASSERT(!"何かがおかしい！");
        break;
    }

    cp.checkpoint("コマンドラインオプション解析処理", __LINE__);

    try {
        ScfLoop sl(mg.getpairdata());

        cp.checkpoint("初期化処理", __LINE__);

        auto [pdiffdata, wavefunctions] = sl();

        cp.checkpoint("微分方程式の積分と固有値探索処理及び規格化処理", __LINE__);

        Energy(
            pdiffdata,
            wavefunctions.at("1 Mesh (r)"),
            wavefunctions.at("2 Eigen function"),
            pdiffdata->pdata_->Z_).express_energy(sl.PEhartree);

        cp.checkpoint("エネルギー出力処理", __LINE__);

		WaveFunctionSave wfs(wavefunctions, pdiffdata->pdata_);
		wfs();

		cp.checkpoint("ファイル書き込み処理", __LINE__);
    }
    catch (std::runtime_error const & e) {
        std::cerr << e.what() << std::endl;
        goexit();
        return EXIT_FAILURE;
    }
        
    cp.checkpoint_print();
    cp.totalpassageoftime();

    checkpoint::usedmem();

    goexit();

    return EXIT_SUCCESS;
}

