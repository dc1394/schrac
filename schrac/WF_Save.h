#include "WF_Normalize.h"

namespace schrac {
	class WF_Save
		: private boost::noncopyable {
#if (_MSC_VER >= 1600)
			std::string make_filename(const std::shared_ptr<const Data> & pdata) const;
#else
			const std::string make_filename(const shared_ptr<const Data> & pdata) const;
#endif
	public:
		WF_Save() {}
		bool operator()(const std::shared_ptr<const Data> & pdata, const WF_Normalize::d3tup & tup) const;
	};

	struct FILEDeleter {
		void operator()(FILE * const fp) const {
			std::fclose(fp);
		};
	};
}
