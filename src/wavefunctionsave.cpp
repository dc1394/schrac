/*! \file wavefunction.cpp
    \brief 得られた波動関数をファイルに書き出すクラスの実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/

#include "wavefunctionsave.h"
#include <cstdlib>              // for std::fclose, std::fprintf
#include <iostream>             // for std::cerr
#include <tuple>                // for std::tie
#include <boost/cast.hpp>       // for boost::numeric_cast

namespace schrac {
    WaveFunctionSave::WaveFunctionSave(boost::container::flat_map<std::string, std::vector<double>> const & wf, std::shared_ptr<Data> const & pdata) :
        wf_(wf),
        pdata_(pdata)
    {
    }
    
    boost::optional<std::string> WaveFunctionSave::get_spin_orbital() const
    {
        switch (pdata_->eq_type_) {
        case Data::Eq_type::DIRAC:
            if (pdata_->spin_orbital_ == Data::ALPHA) {
                return boost::optional<std::string>("alpha");
            }
            else {
                return boost::optional<std::string>("beta");
            }
            break;

        case Data::Eq_type::SCH:
        case Data::Eq_type::SDIRAC:
            return boost::none;
            break;

        default:
            BOOST_ASSERT(!"何かがおかしい！");
            return boost::none;
        }
    }

    std::pair<std::string, std::string> WaveFunctionSave::make_filename() const
    {
        std::string wffilename("wavefunction_");
        std::string rhofilename("wf_");

        auto filename = pdata_->chemical_symbol_ + '_';
        filename += pdata_->orbital_.c_str();

        auto const pspin_orbital(get_spin_orbital());
        if (pspin_orbital) {
            filename += '_' + *pspin_orbital;
        }

        filename += ".csv";

        return std::make_pair(wffilename + filename, rhofilename + filename);
    }

    bool WaveFunctionSave::operator()()
    {
        std::string wffilename, rhofilename;
        std::tie(wffilename, rhofilename) = make_filename();

        auto const fcloser = [](FILE * fp)
        {
            if (fp) {
                std::fclose(fp);
            }
        };

        std::unique_ptr<FILE, decltype(fcloser)> wffp(
            std::fopen(wffilename.c_str(), "w"),
            fcloser);

        if (!wffp) {
            std::cerr << "波動関数のファイルが作成できませんでした。" << std::endl;
            return false;
        }
        
        std::unique_ptr<FILE, decltype(fcloser)> rhofp(
            std::fopen(rhofilename.c_str(), "w"),
            fcloser);

        if (!rhofp) {
            std::cerr << "波動関数のファイルが作成できませんでした。" << std::endl;
            return false;
        }

        auto const end(wf_.end());
        for (auto itr(wf_.begin()); itr != end; ++itr) {
            std::fprintf(wffp.get(), "%s,", itr->first.c_str());
        }
        std::fputs("\n", wffp.get());
        
        auto const size = boost::numeric_cast<std::int32_t>(wf_.begin()->second.size());
        for (auto i = 0; i < size; i++) {
            for (auto itr(wf_.begin()); itr != end; ++itr) {
                std::fprintf(wffp.get(), "%.15f,", itr->second[i]);
            }
            std::fputs("\n", wffp.get());

            std::fprintf(rhofp.get(), "%.15f,", wf_["1 Mesh (r)"][i]);
            std::fprintf(rhofp.get(), "%.15f\n", wf_["2 Eigen function"][i]);
        }

        std::cout << wffilename << " に波動関数を、 " << rhofilename << " に電子密度を書き込みました。" << std::endl;

        return true;
    }
}
