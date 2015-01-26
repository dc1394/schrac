#ifndef _CHKPOINT_H_
#define _CHKPOINT_H_

#include <utility>
#include <boost/noncopyable.hpp>

#if (_MSC_VER >= 1600)
	#include <memory>
#else
	#include <boost/checked_delete.hpp>
	#include <boost/interprocess/smart_ptr/unique_ptr.hpp>
#endif

namespace CheckPoint {
#if (_MSC_VER >= 1600)
	using std::unique_ptr;
#else
	using boost::interprocess::unique_ptr;
#endif

#ifdef _WIN32
	void usedmem();
#endif

	class ChkPoint :
		private boost::noncopyable {
			struct ChkPointFastImpl;
			struct Timestamp;

			template <typename T>
			struct fastpimpl_deleter {
				void operator()(T * const p) const {
					FastArenaObject<sizeof(ChkPoint::ChkPointFastImpl)>::
						operator delete(reinterpret_cast<void *>(p));
				}
			};

			const unique_ptr<ChkPointFastImpl, const fastpimpl_deleter<ChkPointFastImpl> > cfp;
	
	public:
		typedef std::pair<const double, const double> dpair;

		ChkPoint(const char * const func, int line);
		~ChkPoint();
		void checkpoint(const char * const func, int line);
#if (_MSC_VER >= 1600)
		ChkPoint::dpair ChkPoint::totalpassageoftime() const;
#else
		const ChkPoint::dpair ChkPoint::totalpassageoftime() const;
#endif
		void checkpoint_print() const;
	};
}

#endif	// _CHKPOINT_H_
