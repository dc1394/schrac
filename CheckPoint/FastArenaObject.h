#ifndef _FAST_ARENA_OBJECT_H_
#define _FAST_ARENA_OBJECT_H_

#include "ArraiedAllocator.h"

namespace CheckPoint {
	template <size_t TTypeSize, size_t TNumArray = 1>
	struct FastArenaObject :
		private boost::noncopyable {
			// サイズは絶対０より大きくなくちゃダメ
			BOOST_STATIC_ASSERT(TNumArray > 0);

			static void * operator new(std::size_t dummy) {
				return ArraiedAllocator<TTypeSize, TNumArray>::GetAllocator().Alloc();
			}
			static void operator delete(void * p) {
				ArraiedAllocator<TTypeSize, TNumArray>::GetAllocator().Free(p);
			}
	};
}

#endif // _FAST_ARENA_OBJECT_H_
