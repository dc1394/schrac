#include "Adaptive_Step.h"

namespace schrac {
	void Adaptive_Step::RungeKutta()
	{
		long double k1[2], k2[2], k3[2], k4[2];

		for (std::size_t i = 0; i < DiffData::AVECSIZE - 1; i++) {
			k1[0] = pdiffdata_->DX * dL_dx(pdiffdata_->MO[i]);
			k1[1] = pdiffdata_->DX * Diff::dM_dx(pdiffdata_->XV_O[i], pdiffdata_->LO[i],
												 pdiffdata_->MO[i], pdiffdata_);

			k2[0] = pdiffdata_->DX * dL_dx(pdiffdata_->MO[i] + k1[1] / 2.0);
			k2[1] = pdiffdata_->DX * Diff::dM_dx(pdiffdata_->XV_O[i] + pdiffdata_->DX / 2.0,
												 pdiffdata_->LO[i] + k1[0] / 2.0,
												 pdiffdata_->MO[i] + k1[1] / 2.0,
												 pdiffdata_);

			k3[0] = pdiffdata_->DX * dL_dx(pdiffdata_->MO[i] + k2[1] / 2.0);
			k3[1] = pdiffdata_->DX * Diff::dM_dx(pdiffdata_->XV_O[i] + pdiffdata_->DX / 2.0,
												 pdiffdata_->LO[i] + k2[0] / 2.0,
												 pdiffdata_->MO[i] + k2[1] / 2.0,
												 pdiffdata_);

			k4[0] = pdiffdata_->DX * dL_dx(pdiffdata_->MO[i] + k3[1]);
			k4[1] = pdiffdata_->DX * Diff::dM_dx(pdiffdata_->XV_O[i] + pdiffdata_->DX,
												 pdiffdata_->LO[i] + k3[0],
												 pdiffdata_->MO[i] + k3[1],
												 pdiffdata_);
		
			pdiffdata_->LO[i + 1] = pdiffdata_->LO[i] + 
									(k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0;
			pdiffdata_->MO[i + 1] = pdiffdata_->MO[i] +
									(k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0;
		}
	}

	bool Adaptive_Step::rkqc(long double htry,
							 const AdapStepHelper::darray & yscal,
							 const shared_ptr<AdapStepHelper> & pasa) const
	{
		AdapStepHelper::darray yold, dyold, ytemp;
		long double xold = pasa->xx_;						// 初期値の格納

		for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++) {
			yold[i] = pasa->y_[i];
			dyold[i] = pasa->dydx_[i];
		}

		long double h = htry;								// 刻み幅を初期の仮の値に設定

		while (true) {
			const long double hh = 0.5 * h;					// 2つの半ステップを取る

			rk4(xold, hh, yold, dyold, ytemp);

			pasa->xx_ = xold + hh;
			derivs(pasa->xx_, ytemp, pasa->dydx_);

			rk4(pasa->xx_, hh, ytemp, pasa->dydx_, pasa->y_);

			pasa->xx_ = xold + h;

			if (std::fabs(pasa->xx_ - xold) <= pasa->TINY) {
				std::cerr << "Step size too small in method rkqc" << std::endl;
				return false;
			}

			rk4(xold, h, yold, dyold, ytemp);				// 全ステップを取る
			long double errmax = 0.0;						// 精度の計算

			for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++) {
				ytemp[i] = pasa->y_[i] - ytemp[i];			// ここで評価誤差はytempに入れられる
				const long double temp = std::fabs(ytemp[i] / yscal[i]);

				if (errmax < temp)
					errmax = temp;
			}

			errmax /= pasa->EPS;							// 許容誤差に対してスケール

			if (errmax <= 1.0) {							// ステップがうまくいったので、次の誤差を計算
				pasa->hdid_ = h;
				pasa->hnext_ = 
					(errmax > AdapStepHelper::ERRCON ?
					AdapStepHelper::SAFETY * h * std::exp(
					AdapStepHelper::PGROW * std::log(errmax)) : 4.0 * h);
				break;
			}

			h = AdapStepHelper::SAFETY * h *
				std::exp(AdapStepHelper::PSHRNK * std::log(errmax));// 打ち切り誤差が大きすぎるとき、刻み幅を減少
		}

		for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++)
			pasa->y_[i] += ytemp[i] * AdapStepHelper::FCOR;			// 5次の打ち切り誤差の後始末

		return true;
	}

	void Adaptive_Step::rk4(long double x, long double h,
							 const AdapStepHelper::darray & y,
							 const AdapStepHelper::darray & dydx,
							 AdapStepHelper::darray & yout) const
	{
		AdapStepHelper::darray yt, dyt, dym;
		const long double hh = h * 0.5;
		const long double h6 = h / 6.0;
		const long double xh = x + hh;

		for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++)
			yt[i] = y[i] + hh * dydx[i];					// First step

		derivs(xh, yt, dyt);								// Second step

		for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++)
			yt[i] = y[i] + hh * dyt[i];

		derivs(xh, yt, dym);								// Third step

		for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++) {
			yt[i] = y[i] + h * dym[i];
			dym[i] += dyt[i];
		}

		derivs(x + h, yt, dyt);								// Fourth step

		for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++)	// 重みをつけて総和を取る
			yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
	}

	bool Adaptive_Step::solve_diff_equ_O()
	{
		// 最初の三項のみ通常のルンゲ・クッタ法で計算
		RungeKutta();
		boost::mpl::int_<AdapStepHelper::OO> o;
		const shared_ptr<AdapStepHelper> pasa(
			make_shared<AdapStepHelper>(pdiffdata_, o));	

		return odeint(pasa);
	}

	bool Adaptive_Step::solve_diff_equ_I()
	{
		boost::mpl::int_<AdapStepHelper::INFINITY> inf;
		const shared_ptr<AdapStepHelper> pasa(
			make_shared<AdapStepHelper>(pdiffdata_, inf));	

		return odeint(pasa);
	}

	bool Adaptive_Step::solve_diff_equ()
	{
#ifdef _OPENMP
		if (pdata_->ompthread_) {
			volatile bool bRet = true;
			bool bRet2 = false;

			#pragma omp parallel sections
			{
				#pragma omp section
				{
					bRet = solve_diff_equ_O();
				}

				#pragma omp section
				{
					if (bRet)
						bRet2 = solve_diff_equ_I();
				}
			}

			return (bRet && bRet2);
		} else {
			return (solve_diff_equ_O() && solve_diff_equ_I());
		}
#else
		BOOST_STAIC_ASSERT(false);

		return false;
#endif
	}

}
