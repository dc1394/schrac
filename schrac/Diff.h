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
        // #region 型エイリアス

        using myarray = std::array < double, DiffData::AVECSIZE > ;

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
        static auto constexpr MINV = 1.0E-200;

		void a_init();
		void b_init();

		void init_LM_O();
		void init_LM_I();

        Diff::myarray solve_linear_equ(std::array<double, DiffData::AVECSIZE * DiffData::AVECSIZE> a, myarray b) const;

		bool solve_diff_equ_O();
		bool solve_diff_equ_I();

        void derivs(std::array<double, 2> const & f, std::array<double, 2> & dfdx, double x) const;

	protected:
		const std::shared_ptr<const Data> pdata_;
		const std::shared_ptr<DiffData> pdiffdata_;

		

	public:
		std::function<double(double, double, double, std::shared_ptr<DiffData> const &)> dM_dx;
		typedef std::pair<const std::array<double, 2>, const std::array<double, 2> > mytuple;

		void Initialize(double E);
		const std::shared_ptr<DiffData> & getpDiffData() const
		{ return pdiffdata_; }
		virtual bool solve_diff_equ();
		mytuple getMPval() const;
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

	template <typename T>
	T sqr(T x)
	{
		return x * x;
	}
}

#endif	// _DIFF_H_
