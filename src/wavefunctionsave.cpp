/*! \file wavefunction.cpp
    \brief 得られた波動関数をファイルに書き出すクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "wavefunctionsave.h"
#include <cstdlib>              // for std::fclose, std::fprintf
#include <iostream>             // for std::cerr
#include <utility>              // for std::move
#include <boost/cast.hpp>       // for boost::numeric_cast

namespace schrac {
    WaveFunctionSave::WaveFunctionSave(boost::container::flat_map<std::string, std::vector<double>> const & wf, std::shared_ptr<Data> const & pdata) :
        wf_(wf),
        pdata_(pdata)
    {
    }

    std::string WaveFunctionSave::make_filename() const
    {
        std::string filename("wavefunction_");

        filename += pdata_->chemical_symbol_ + '_';
        filename += pdata_->orbital_.c_str();

        switch (pdata_->eq_type_) {
        case Data::Eq_type::DIRAC:
            filename += '_';
            if (pdata_->spin_orbital_ == Data::ALPHA) {
                filename += "alpha";
            }
            else {
                filename += "beta";
            }
            break;

        case Data::Eq_type::SCH:
        case Data::Eq_type::SDIRAC:
            break;

        default:
            BOOST_ASSERT(!"何かがおかしい！");
        }

        filename += ".csv";

        return std::move(filename);
    }

    bool WaveFunctionSave::operator()()
    {
        auto const filename(make_filename());

        auto const fcloser = [](FILE * fp)
        {
            if (fp) {
                std::fclose(fp);
            }
        };

        std::unique_ptr<FILE, decltype(fcloser)> fp(
            std::fopen(filename.c_str(), "w"),
            fcloser);

        if (!fp) {
            std::cerr << "波動関数のファイルが作成できませんでした。" << std::endl;
            return false;
        }
        
        std::unique_ptr<FILE, decltype(fcloser)> fp2(
            std::fopen(filename.c_str(), "w"),
            fcloser);

        if (!fp) {
            std::cerr << "波動関数のファイルが作成できませんでした。" << std::endl;
            return false;
        }

        auto const end(wf_.end());
        for (auto itr(wf_.begin()); itr != end; ++itr) {
            std::fprintf(fp.get(), "%s,", itr->first.c_str());
        }
        std::fputs("\n", fp.get());

        auto const size = boost::numeric_cast<std::int32_t>(wf_.begin()->second.size());
        for (auto i = 0; i < size; i++) {
            for (auto itr(wf_.begin()); itr != end; ++itr) {
                std::fprintf(fp.get(), "%.15f,", itr->second[i]);
            }
            std::fputs("\n", fp.get());
        }

        std::cout << '\n' << filename << "に波動関数を書き込みました。" << std::endl;

        return true;
    }
}
