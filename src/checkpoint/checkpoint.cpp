/*! \file checkpoint.cpp
    \brief 時間計測のためのクラスの実装

    Copyright ©  2014 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "checkpoint.h"
#include <iostream>             // for std::cout
#include <optional>				// for std::optional
#include <system_error>         // for std::system_category
#include <boost/assert.hpp>     // for boost::assert
#include <boost/cast.hpp>       // for boost::numeric_cast
#include <boost/format.hpp>     // for boost::format

#ifdef _WIN32
    #include <Windows.h>        // for GetCurrentProcess
    #include <Psapi.h>          // for GetProcessMemoryInfo

    #pragma comment(lib, "Psapi.Lib")
#else
    #include <errno.h>          // for errno 
    #include <sys/time.h>       // for struct timeval
    #include <sys/resource.h>   // for getrusage
#endif

namespace checkpoint {
    CheckPoint::CheckPoint()
        : cfp(
            reinterpret_cast<CheckPoint::CheckPointFastImpl *>(
                FastArenaObject<sizeof(CheckPoint::CheckPointFastImpl)>::operator new(0)))
    {
    }

    CheckPoint::~CheckPoint()
    {
    }

    void CheckPoint::checkpoint(char const * action, std::int32_t line)
    {
        BOOST_ASSERT(cfp->cur < static_cast<std::int32_t>(CheckPoint::CheckPointFastImpl::N));

        auto const p = cfp->points.begin() + cfp->cur;

        p->action = action;
        p->line = line;
        p->realtime = std::chrono::high_resolution_clock::now();

        cfp->cur++;
    }
    
    void CheckPoint::checkpoint_print() const
    {
        using namespace std::chrono;

        std::optional<high_resolution_clock::time_point> prevreal(std::nullopt);

        auto itr = cfp->points.begin();
        for (auto i = 0; i < cfp->cur; ++i, ++itr) {
            if (prevreal) {
                auto const realtime(duration_cast<duration<double, std::milli>>(itr->realtime - *prevreal));
                std::cout << itr->action
                          << boost::format(" elapsed time = %.4f (msec)\n") % realtime.count();
            }

            prevreal = std::optional<high_resolution_clock::time_point>(itr->realtime);
        }
    }

    void CheckPoint::totalpassageoftime() const
    {
        using namespace std::chrono;

        auto const realtime = duration_cast<duration<double, std::milli>>(
            cfp->points[cfp->cur - 1].realtime - cfp->points[0].realtime);

        std::cout << boost::format("Total elapsed time = %.4f (msec)") % realtime.count() << std::endl;
    }

    // #region 非メンバ関数

#ifdef _WIN32
    void usedmem()
    {
        PROCESS_MEMORY_COUNTERS memInfo = { 0 };
        
        if (!::GetProcessMemoryInfo(::GetCurrentProcess(), &memInfo, sizeof(memInfo))) {
            throw std::system_error(std::error_code(::GetLastError(), std::system_category()));
        }

        std::cout << "Used Memory Size: "
                  << boost::numeric_cast<std::uint32_t>(memInfo.PeakWorkingSetSize >> 10)
                  << "(kB)"
                  << std::endl; 
    }
#else
    void usedmem()
    {
        struct rusage r;

        if (getrusage(RUSAGE_SELF, &r)) {
            throw std::system_error(errno, std::system_category());
        }

        std::cout << "Used Memory Size: "
                  << boost::numeric_cast<std::uint32_t>(r.ru_maxrss)
                  << "(kB)"
                  << std::endl;
    }
#endif

    // #endregion 非メンバ関数
}
