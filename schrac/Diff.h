/*! \file diff.h
    \brief 微分方程式を解くクラスの宣言

    Copyright ©  2015 @dc1394 All Rights Reserved.
*/
#ifndef _DIFF_H_
#define _DIFF_H_

#pragma once

#include "DiffData.h"

namespace schrac {
	class Diff {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param E エネルギー固有値
            \param pdata データオブジェクト
            \param TINY 絶対値が小さい方の閾値
            */
        Diff(double E, std::shared_ptr<Data> const & pdata, double TINY);

        //! A destructor.
        /*!
        何もしないデストラクタ
        */
        ~Diff()
        {
        }

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数
        
    public:
        static const double MINV;

			bool a_init();
			void b_init();

			void init_LM_O();
			void init_LM_I();

			const boost::optional<const std::array<double, DiffData::AVECSIZE> >
				S_gausswp(std::array<std::array<double, DiffData::AVECSIZE>, DiffData::AVECSIZE> & a,
						  std::array<double, DiffData::AVECSIZE> & b) const;

			virtual bool solve_diff_equ_O();
			virtual bool solve_diff_equ_I();

            void derivs(double x, std::array<double, 2> const & y, std::array<double, 2> & dydx) const;

	protected:
		const std::shared_ptr<const Data> pdata_;
		const std::shared_ptr<DiffData> pdiffdata_;

		

	public:
		static std::function<double(double, double, double,
			const std::shared_ptr<DiffData> &)> dM_dx;
		typedef std::pair<const std::array<double, 2>, const std::array<double, 2> > mytuple;

		void Initialize(double E);
		const std::shared_ptr<DiffData> & getpDiffData() const
		{ return pdiffdata_; }
		virtual bool solve_diff_equ();
		virtual ~Diff() {}
		const mytuple getMPval() const;
	};

    double fnc_V(double r, const std::shared_ptr<DiffData> & pdiffdata);
    double dV_dr(double r, const std::shared_ptr<DiffData> & pdiffdata);
    double dL_dx(double M);
    double dM_dx_sch(double x, double L, double M,
        const std::shared_ptr<DiffData> & pdiffdata);
    double dM_dx_sdirac(double x, double L, double M,
        const std::shared_ptr<DiffData> & pdiffdata);
    double dM_dx_dirac(double x, double L, double M,
        const std::shared_ptr<DiffData> & pdiffdata);

    void Diff::derivs(double x, std::array<double, 2> const & y, std::array<double, 2> & dydx) const
    {
        dydx[0] = dL_dx(y[1]);
        dydx[1] = Diff::dM_dx(x, y[0], y[1], pdiffdata_);
    }

	inline double fnc_V(double r, const std::shared_ptr<DiffData> & pdiffdata)
	{
		return - pdiffdata->Z_ / r;
	}

	inline double dV_dr(double r, const std::shared_ptr<DiffData> & pdiffdata)
	{
		return pdiffdata->Z_ / (r * r);
	}
		
	inline double dL_dx(double M)
	{
		return M;
	}

	inline const Diff::mytuple Diff::getMPval() const
	{
		std::array<double, 2> L, M;

		L[0] = pdiffdata_->LO_[pdiffdata_->MP_O_];
		L[1] = pdiffdata_->LI_[pdiffdata_->MP_I_];
		M[0] = pdiffdata_->MO_[pdiffdata_->MP_O_];
		M[1] = pdiffdata_->MI_[pdiffdata_->MP_I_];

		return std::make_pair(L, M);
	}
	
	template <typename T>
	T sqr(T x)
	{
		return x * x;
	}

    template <typename T>
    T pow(T x, std::uint32_t n)
    {
        T p = x, y = 1.0;

        while (true) {
            if (n & 1)
                y *= p;

            n >>= 1;

            if (!n)
                return y;

            p *= p;
        }
    }
}

#endif	// _DIFF_H_
