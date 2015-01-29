#include "getcomlineoption.h"
#include <iostream>
#include <boost/program_options.hpp>

namespace schrac {
    std::string const GetComLineOption::DEFINPNAME = "input.inp";        
        
    //! A public member function.
    /*!
        コマンドラインオプションを解析する関数の実装
        \param argc コマンドライン引数の数
        \param argv コマンドライン引数
        \return 解析に成功したら0、失敗したら-1、ヘルプ表示の場合は1
    */
    std::int32_t GetComLineOption::getopt(int argc, char * const argv[])
	{
		using namespace boost::program_options;

        // オプションの設計
		options_description opt("option");

		// 引数の書式を定義
        opt.add_options()
            ("help,h", "Show help")
            ("inputfile,I", value<std::string>()->default_value(GetComLineOption::DEFINPNAME), "Specifying the input file name")
			("tbb,T", value<bool>()->implicit_value(false),
			 "Whether or not to use the TBB (default: is not used)");

		// 引数の書式に従って実際に指定されたコマンドライン引数を解析
		variables_map vm;
		try {
			store(parse_command_line(argc, argv, opt), vm);
		} catch (const std::exception & e) {
			std::cerr << e.what()
					  << ". command line argument is valid!" << std::endl;

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
    std::pair<std::string, bool> GetComLineOption::getpairdata() const
	{
        return std::make_pair(inpname_, usetbb_);
    }
}

