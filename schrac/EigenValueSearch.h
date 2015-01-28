#include "Diff.h"
#include "ReadInpFile.h"
#include "RungeKutta.h"
#include "RK_AdapStep.h"
#include "Bulirsch_Stoer.h"

#ifdef _DEBUG
	#include "ModEuler.h"
#endif

namespace schrac {
	template <typename T> T sign(T a, T b);

	long double Eexact_sch(const shared_ptr<const Data> & pdata);
	long double Eexact_sdirac(const shared_ptr<const Data> & pdata);
	long double Eexact_dirac(const shared_ptr<const Data> & pdata);

	class EigenValueSearch
		: private boost::noncopyable {
		static const int EVALSEARCHMAX = 1000;

		static const long double TINY;
		static const long double HUGE;

		long double EPS;
		long double TOL;
		
		
		shared_ptr<const Data> pdata_;
		shared_ptr<Diff> pdiff_;
		shared_ptr<DiffData> pdiffdata_;

		int i_;
		bool bNoden;

		long double Eexact;
		long double E;
		long double DE;
		long double Emax;
		long double Emin;
		long double Dold;
		long double Dmax;
		long double Dmin;

		const boost::optional<const long double> fnc_D();
		bool rough_search();
		bool zbrent();
		void setExp() const;
		void msg() const;
		void info() const;
		void info(long double E) const;
		void info(long double b, long double fb) const;
		void init();

	public:
		explicit EigenValueSearch(const tuple<const std::string, const int> & arg);
		bool search();
		const shared_ptr<const Data> & getpData() const
		{ return pdata_; }
		const shared_ptr<Diff> & getpDiff() const
		{ return pdiff_; }
	};

	inline void EigenValueSearch::setExp() const
	{
		std::cout.setf(std::ios::fixed, std::ios::floatfield);
		std::cout << std::setprecision(
			boost::numeric_cast<const std::streamsize>(std::fabs(std::log10(pdata_->eps))) - 2);
	}

	inline void EigenValueSearch::info() const
	{
		std::cout << "i = " << i_ << ", D = " << Dold << ", node = "
				  << pdiffdata_->thisnode;
		if (bNoden)
			std::cout << " (OK)" << std::endl;
		else
			std::cout << " (NG)" << std::endl;
	}

	inline void EigenValueSearch::info(long double E) const
	{
		std::cout << "ノード数が一致する固有値を発見しました！" << std::endl;
		std::cout << "E(厳密)   = " << Eexact << " (Hartree)" << std::endl;
		std::cout << "E(計算値) = " << E << " (Hartree)" << std::endl;
	}

	inline void EigenValueSearch::info(long double b, long double fb) const
	{
		std::cout << "i = " << i_ << ", D = "
				  << fb << ", E = " << b
				  << ", node = "
				  << pdiffdata_->thisnode;
		if (bNoden)
			std::cout << " (OK)" << std::endl;
		else
			std::cout << " (NG)" << std::endl;
	}

	template <typename T>
	inline T sign(T a, T b)
	{
		return (b >= 0.0) ? std::fabs(a) : - std::fabs(a);
	}

	inline long double Eexact_sch(const shared_ptr<const Data> & pdata)
	{
		return - sqr<long double>(static_cast<const long double>(pdata->Z) /
			   static_cast<const long double>(pdata->n)) / 2.0;
	}

	inline long double Eexact_sdirac(const shared_ptr<const Data> & pdata)
	{
		const long double nr = static_cast<const long double>(pdata->n) - 1;
		const long double lambda = std::sqrt(1 - sqr(static_cast<const long double>(pdata->Z) / Data::c));
		const long double denominator = std::sqrt(sqr(nr + lambda) +
			sqr(static_cast<const long double>(pdata->Z) / Data::c));

		return ((nr + lambda) / denominator - 1.0) * sqr(Data::c);
	}

	inline long double Eexact_dirac(const shared_ptr<const Data> & pdata)
	{
		const long double nr = static_cast<const long double>(pdata->n) - pdata->j_- 0.5;
		const long double lambda = std::sqrt(sqr(pdata->kappa) -
			sqr(static_cast<const long double>(pdata->Z) / Data::c));
		const long double denominator = std::sqrt(sqr(nr + lambda) +
			sqr(static_cast<const long double>(pdata->Z) / Data::c));

		return ((nr + lambda) / denominator - 1.0) * sqr(Data::c);
	}
}
