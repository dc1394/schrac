#ifndef _MY_GETOPT_H_
#define _MY_GETOPT_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#endif

#include <boost/program_options.hpp>
#include "stdHydroSchDirac.h"

namespace HydroSchDirac {
	class My_getOpt :
		private boost::noncopyable {
			const int defompthread;
			static const std::string DEFINPNAME;
			int ompthread_;
			std::string inpname_;
	public:
		My_getOpt();
		int getOpt(int argc, char * const argv[]);
		const tuple<const std::string, const int> getpData() const
		{ return make_tuple(inpname_, ompthread_); }
		int getOmpThread() const
		{ return ompthread_; }
	};

	inline My_getOpt::My_getOpt()
	 :	defompthread(omp_get_num_procs()),
		inpname_(DEFINPNAME)
	{
	}
}

#endif	// _MY_GETOPT_H_
