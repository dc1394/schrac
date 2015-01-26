#ifndef _STDHYDROSCHDIRAC_H_
#define _STDHYDROSCHDIRAC_H_

#if defined(_MSC_VER) && (_MSC_VER >= 1020)
#pragma once
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <cstdio>
#include <boost/cast.hpp>
#include <boost/assert.hpp>
#include <boost/format.hpp>
#include <boost/mpl/int.hpp>
#include <boost/optional.hpp>
#include <boost/noncopyable.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/static_assert.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/utility/in_place_factory.hpp>

#ifdef _DEBUG
	#include <vector>
#else
	#include <tbb/concurrent_vector.h>
#endif 

#if (_MSC_VER >= 1600)
#include <array>
#include <tuple>
#include <memory>
#include <utility>
#include <functional>
#else
#include <boost/array.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/interprocess/smart_ptr/unique_ptr.hpp>
#endif

#if defined(_WIN32) || defined(_WIN64)
#include <conio.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#else
BOOST_STATIC_ASSERT(false);
#endif

namespace HydroSchDirac {
#ifdef _DEBUG
	typedef std::vector<long double> ldvector;
#else
	typedef tbb::concurrent_vector<long double> ldvector;
#endif 

#if (_MSC_VER >= 1600)
	using std::get;
	using std::tie;
	using std::tuple;
	using std::make_tuple;
	using std::array;
	using std::function;
	using std::shared_ptr;
	using std::make_shared;
	using std::unique_ptr;
#else
	using boost::get;
	using boost::tie;
	using boost::tuples::tuple;
	using boost::make_tuple;
	using boost::array;
	using boost::function;
	using boost::shared_ptr;
	using boost::make_shared;
	using boost::interprocess::unique_ptr;
#endif
}

#endif	//  _STDHYDROSCHDIRAC_H_
