/*! \file schracmain.cpp
    \brief メインファイル

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "energy.h"
#include "getcomlineoption.h"
#include "normalization.h"
#include "scfloop.h"
#include "wavefunctionsave.h"

#if defined(_WIN32) || defined(_WIN64)
    #include <conio.h>                          // for _getch
#endif

#include <cstdlib>                              // for EXIT_FAILURE, EXIT_SUCCESS
#include <iostream>                             // for std::cerr
#include <boost/optional.hpp>                   // for boost::optional
#include <boost/utility/in_place_factory.hpp>   // for boost::in_place
//#include "ChkPoint.h"

namespace schrac {
    void goexit();
}

int main(int argc, char * argv[])
{
    using namespace schrac;
    //	CheckPoint::ChkPoint cp("処理開始", __LINE__);

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

    //cp.checkpoint("コマンドラインオプション解析処理", __LINE__);

    ScfLoop sl(mg.getpairdata());
    sl();
    
    //Energy en(evs.PDiffSolver()->PDiffData, pdsol.at("1 Mesh (r)"), pdsol.at("2 Eigen function"));
    //en.express_energy();
    //WaveFunctionSave wfs(pdsol, evs.PData);
    //wfs();

    goexit();
}
	/*cp.checkpoint("微分方程式の積分と固有値探索処理", __LINE__);

	schrac::WF_Normalize wfn(pevs->getpDiff());

	cp.checkpoint("波動関数のコピー処理", __LINE__);

	wfn();
	
	cp.checkpoint("波動関数の正規化処理", __LINE__);
	
	schrac::WF_Save wfs;
	if (!wfs(pevs->getpData(), wfn.getptup())) {
		std::cerr << "ファイルを作成できませんでした。終了します。" << std::endl;
		return EXIT_FAILURE;
	}
	cp.checkpoint("ファイル書き込み処理", __LINE__);

	pevs = boost::none;

	cp.checkpoint("終了処理", __LINE__);

	std::cout << "正常に終了しました。" << std::endl << std::endl;
	schrac::showomp(go);

	cp.checkpoint_print();
	double cpu, real;
	tie(cpu, real) = cp.totalpassageoftime();

	std::cout << "総処理時間： CPU:" << boost::format("%.4f") % cpu
			  << " (msec), 実時間:" << boost::format("%.4f") % real
			  << " (msec), 並列化効率:" << boost::format("%.4f") % (cpu / real)
			  << "（倍）" << std::endl << std::endl;

#if defined(_WIN32) || defined(_WIN64)
	CheckPoint::usedmem();
#endif

	std::cout << "終了するには何かキーを押してください..." << std::endl;

#if defined(_WIN32) || defined(_WIN64)
	::_getch();
#else
	std::getchar();
#endif

	return EXIT_SUCCESS;
}*/

namespace schrac {
    void goexit()
    {
        std::cout << "終了するには何かキーを押してください..." << std::endl;

#if defined(_WIN32) || defined(_WIN64)
        ::_getch();
#else
        std::getchar();
#endif
    }
}
