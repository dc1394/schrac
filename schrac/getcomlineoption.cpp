/*! \file getcomlineoption.cpp
    \brief コマンドラインオプションの解析を行うクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/
#include "getcomlineoption.h"
#include <iostream>                     // for std::cerr, std::cout
#include <boost/program_options.hpp>    // for boost::program_options

namespace schrac {
    // #region staticメンバ変数
    
    std::string const GetComLineOption::DEFINPNAME = "input.inp";        
        
    // #endregion staticメンバ変数

    // #region publicメンバ関数

    std::int32_t GetComLineOption::getopt(int argc, char * const argv[])
	{
		using namespace boost::program_options;

        // オプションの設計
		options_description opt("option");

		// 引数の書式を定義
        opt.add_options()
            ("help,h", "ヘルプを表示")
            ("inputfile,I", value<std::string>()->default_value(GetComLineOption::DEFINPNAME), "インプットファイル名")
			("tbb,T", value<bool>()->implicit_value(false),
			 "TBBを使用して並列計算を行うかどうか（デフォルトはTBBを使用しない）");

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
    
    std::pair<std::string, bool> GetComLineOption::getpairdata() const
	{
        return std::make_pair(inpname_, usetbb_);
    }

    // #endregion publicメンバ関数
}

