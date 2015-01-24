#include "My_getopt.h"

namespace HydroSchDirac {
#ifdef _OPENMP
	const std::string My_getOpt::DEFINPNAME("input.inp");
	
	int My_getOpt::getOpt(int argc, char * const argv[])
	{
		using namespace boost::program_options;

		options_description opt("オプション");
		// 引数の書式を定義
		opt.add_options()
			("help,h", "ヘルプを表示")
			("inputfile,I", value<std::string>()->default_value(DEFINPNAME), "インプットファイル名")
			("openmp,O", value<int>()->implicit_value(defompthread),
			 "OpenMPで使用するスレッド数（0かオプション省略で未使用、デフォルトはプロセッサ数）");

		// 引数の書式に従って実際に指定されたコマンドライン引数を解析
		variables_map vm;
		try {
			store(parse_command_line(argc, argv, opt), vm);
		} catch (const std::exception & e) {
			std::cerr << e.what()
					  << ". コマンドライン引数が異常です。終了します。" << std::endl;

			return -1;
		}
		notify(vm);

		// ヘルプ表示指定がある場合、ヘルプ表示して終了
		if (vm.count("help")) {
			std::cout << opt << std::endl;

			return 1;
		}

		// インプットファイル名指定がある場合
		if (vm.count("inputfile")) {
			inpname_ = vm["inputfile"].as<std::string>();
		}

		// OpenMP指定がある場合
		if (vm.count("openmp")) {
			ompthread_ = vm["openmp"].as<int>();
			if (ompthread_ < 0) {
				std::cerr << "使用スレッド数の指定が異常です。"
						  << "コマンドライン引数が異常です。終了します。" << std::endl;

				return -1;
			}
		} else {
			ompthread_ = 0;
		}

		return 0;
	}
#else
	BOOST_STATIC_ASSERT(false);
#endif
}
