#include "my_getopt.h"
#include <iostream>
#include <boost/program_options.hpp>

namespace schrac {
    std::string const my_getopt::DEFINPNAME = "input.inp";        
        
    //! A public member function.
    /*!
        コマンドラインオプションを解析する関数の実装
        \param argc コマンドライン引数の数
        \param argv コマンドライン引数
        \return 解析に成功したら0、失敗したら-1、ヘルプ表示の場合は1
    */
    std::int32_t my_getopt::getopt(int argc, char * const argv[])
	{
		using namespace boost::program_options;

		options_description opt("オプション");
		// 引数の書式を定義
		opt.add_options()
			("help,h", "ヘルプを表示")
			("inputfile,I", value<std::string>()->default_value(my_getopt::DEFINPNAME), "インプットファイル名")
			("tbb,T", value<bool>()->implicit_value(false),
			 "Threading Building Blocksを使用するかどうか。デフォルトでは使用しない。");

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

		// TBB指定がある場合
		if (vm.count("tbb")) {
			usetbb_ = vm["tbb"].as<bool>();
		}

		return 0;
	}

    //! A public member function (constant).
    /*!
        \return インプットファイル名とTBBを使用するかどうかのstd::pair
    */
    std::pair<std::string, bool> my_getopt::getpairdata() const
	{
        return std::make_pair(inpname_, usetbb_);
    }
}

