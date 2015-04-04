/*! \file schracmain.cpp
    \brief メインファイル

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD-2 License.
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
#include <boost//format.hpp>                    // for boost::format
#include <boost/optional.hpp>                   // for boost::optional
#include <boost/utility/in_place_factory.hpp>   // for boost::in_place

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

    std::shared_ptr<DiffData> pdiffdata;
    ScfLoop::mymap wavefunctions;
    try {
        ScfLoop sl(mg.getpairdata());

        cp.checkpoint("初期化処理", __LINE__);

        std::tie(pdiffdata, wavefunctions) = sl();

        cp.checkpoint("微分方程式の積分と固有値探索処理及び正規化処理", __LINE__);

        Energy(
            pdiffdata,
            wavefunctions.at("1 Mesh (r)"),
            wavefunctions.at("2 Eigen function"),
            pdiffdata->pdata_->Z_).express_energy(sl.PEhartree);

        cp.checkpoint("エネルギー出力処理", __LINE__);
    }
    catch (std::runtime_error const & e) {
        std::cerr << e.what() << std::endl;
        goexit();
        return EXIT_FAILURE;
    }
    
    WaveFunctionSave wfs(wavefunctions, pdiffdata->pdata_);
    wfs();

    cp.checkpoint("ファイル書き込み処理", __LINE__);

    cp.checkpoint_print();
    cp.totalpassageoftime();

    checkpoint::usedmem();

    goexit();

    return EXIT_SUCCESS;
}

