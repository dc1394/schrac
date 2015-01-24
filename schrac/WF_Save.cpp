#include "WF_Save.h"

namespace HydroSchDirac {
#if (_MSC_VER >= 1600)
	std::string WF_Save::make_filename(const shared_ptr<const Data> & pdata) const
#else
	const std::string WF_Save::make_filename(const shared_ptr<const Data> & pdata) const
#endif
	{
		std::string filename("wavefunction_");

		filename += pdata->atom + '_';		
		filename += pdata->orbital.c_str();
		
		if (pdata->eqtype == Data::DIRAC) {
			filename += '_';
			if (pdata->spin_orbital == Data::ALPHA)
				filename += "Alpha";
			else
				filename += "Beta";
		}
		
		filename += ".csv";

#if (_MSC_VER >= 1600)
		return std::move(filename);
#else
		return filename;
#endif
	}

	bool WF_Save::operator()(const shared_ptr<const Data> & pdata, const WF_Normalize::d3tup & tup) const
	{
		const std::string filename(make_filename(pdata));
		const unique_ptr<FILE, const FILEDeleter> fp(std::fopen(filename.c_str(), "w"), FILEDeleter());
		if (fp == NULL) {
			std::cerr << "ファイルが開けません。" << std::endl;
			return false;
		}

		const ldvector & RV(get<0>(tup));
		const ldvector & RF(get<1>(tup));
		const ldvector & PF(get<2>(tup));

		const std::size_t size = RV.size();
		for (std::size_t i = 0; i < size; i++) {
			fprintf(fp.get(), "%.15f,%.15f,%.15f\n", RV[i], RF[i], PF[i]);
		}

		std::cout << '\n' << filename << "に波動関数を書き込みました。" << std::endl;

		return true;
	}
}
