#include "my_getopt.h"
//#include "WF_Save.h"
//#include "ChkPoint.h"

/*namespace HydroSchDirac {
	void showomp(const HydroSchDirac::My_getOpt & go);
}*/

int main(int argc, char * argv[])
{
}
/*#if (_MSC_VER >= 1600)
	using std::tie;
#else
	using boost::tie;
#endif

	CheckPoint::ChkPoint cp("処理開始", __LINE__);
	
	HydroSchDirac::My_getOpt go;
	switch (go.getOpt(argc, argv)) {
		case -1:
			return EXIT_FAILURE;
		break;

		case 1:
			return EXIT_SUCCESS;
		break;

		case 0:
		break;

		default:
			BOOST_ASSERT(!"何かがおかしい！！");
		break;
	}

	cp.checkpoint("コマンドラインオプション解析処理", __LINE__);

	boost::optional<HydroSchDirac::EigenValueSearch> pevs;

	try {
		pevs = boost::in_place(go.getpData());
	} catch (const std::runtime_error & e) {
		std::cerr << e.what() << std::endl;

		return EXIT_FAILURE;
	}

	cp.checkpoint("入力ファイル読み込み処理", __LINE__);

	if (!pevs->search()) {
		std::cerr << "固有値が見つかりませんでした。終了します。" << std::endl;
		return EXIT_FAILURE;
	}
	
	cp.checkpoint("微分方程式の積分と固有値探索処理", __LINE__);

	HydroSchDirac::WF_Normalize wfn(pevs->getpDiff());

	cp.checkpoint("波動関数のコピー処理", __LINE__);

	wfn();
	
	cp.checkpoint("波動関数の正規化処理", __LINE__);
	
	HydroSchDirac::WF_Save wfs;
	if (!wfs(pevs->getpData(), wfn.getptup())) {
		std::cerr << "ファイルを作成できませんでした。終了します。" << std::endl;
		return EXIT_FAILURE;
	}
	cp.checkpoint("ファイル書き込み処理", __LINE__);

	pevs = boost::none;

	cp.checkpoint("終了処理", __LINE__);

	std::cout << "正常に終了しました。" << std::endl << std::endl;
	HydroSchDirac::showomp(go);

	cp.checkpoint_print();
	long double cpu, real;
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
}

namespace HydroSchDirac {
	void showomp(const HydroSchDirac::My_getOpt & go)
	{
		std::cout << "OpenMP: ";

		const int nthread = go.getOmpThread();
		if (nthread) {
			std::cout << "使用" << std::endl;
			std::cout << "使用スレッド数: " << nthread
					  << '\n' << std::endl;
		} else {
			std::cout << "未使用\n" << std::endl;
		}
	}
}*/
