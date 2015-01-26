#ifndef _READINWFILE_H_
#define _READINWFILE_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include "Data.h"

namespace HydroSchDirac {
	using boost::algorithm::split;
	using boost::algorithm::is_any_of;
	using boost::algorithm::token_compress_on;

	class ReadInpFile :
		private boost::noncopyable {
			static const std::streamsize BUFSIZE = 1024;					// バッファサイズ(適当）
			typedef std::vector<ci_string> strvec;							// 何回も使うのでtypedef
		
			std::size_t i_;
			std::ifstream ifs;
			shared_ptr<Data> pdata_;

#if (_MSC_VER >= 1600)
			ci_string readData(const char * const article);
			ci_string readData(const char * const article, ci_string && def);
#else
			const ci_string readData(const char * const article);
			const ci_string readData(const char * const article, const ci_string & def);
#endif
			const boost::optional<const ci_string> readDataAuto(const char * const article);
			template <typename T>
			const boost::optional<const T> readData(const char * const article, T def_val);

			bool readAtom();
			bool readEq();
			bool readGrid();
			bool readEps();
			bool readType();
			bool readLowerE();
			bool readNumofp();
			bool readRatio();

			void errMsg(const char * const s) const;
			void errMsg(const char * const s1, const ci_string & s2) const;

	public:
		ReadInpFile(const tuple<const std::string, const int> & arg);
		void readFile();

#if (_MSC_VER >= 1600)
		shared_ptr<Data> && getpData()
		{ return std::move(pdata_); }
#else
		const shared_ptr<Data> & getpData() const
		{ return pdata_; }
#endif
	};

	inline ReadInpFile::ReadInpFile(const tuple<const std::string, const int> & arg)
	 :	i_(1), ifs(get<0>(arg).c_str()),
		pdata_(make_shared<Data>(get<1>(arg)))
	{
	}

	inline void ReadInpFile::errMsg(const char * const s) const
	{
		std::cerr << "インプットファイルに" << s << "の行が見つかりませんでした" << std::endl;
	}

	inline void ReadInpFile::errMsg(const char * const s1, const ci_string & s2) const
	{
		std::cerr << "インプットファイルの" << s1 << "の行が正しくありません" << std::endl;
		std::cerr << i_ << "行目, 未知のトークン:" << s2 << std::endl;
	}
}

#endif	// _READINWFILE_H_
