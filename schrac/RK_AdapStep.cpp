#include "RK_AdapStep.h"

namespace HydroSchDirac {
	bool RK_AdapStep::odeint(const shared_ptr<AdapStepHelper> & pasa) const
	{
		AdapStepHelper::darray yscal;
		long double h = (pasa->x1_ < pasa->x2_) ? std::fabs(pasa->H1) : - std::fabs(pasa->H1);

		for (unsigned int nstp = 0; nstp < pasa->MAXSTEP; nstp++) {						// Take at most MAXSTP steps
			derivs(pasa->xx_, pasa->y_, pasa->dydx_);

			for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++)						// 精度を監視するためのスケーリング
				yscal[i] = std::fabs(pasa->y_[i]) + std::fabs(pasa->dydx_[i] * h) + pasa->TINY;

			if (std::fabs(pasa->xx_ - (*pasa->xp_)[pasa->i_]) < pasa->TINY) {
				for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++)
					(*pasa->yp_[i])[pasa->i_] = pasa->y_[i];

				pdiffdata_->node_count(pasa->i_, *pasa->yp_[0]);						// ノード数チェック
				pasa->i_++;
			}

			if ((pasa->x1_ < pasa->x2_) && (pasa->xx_ + h) > (*pasa->xp_)[pasa->i_] ||
				(pasa->x1_ > pasa->x2_) && (pasa->xx_ + h) < (*pasa->xp_)[pasa->i_]) {
				h = (*pasa->xp_)[pasa->i_] - pasa->xx_;
			}

			if (!rkqc(h, yscal, pasa))
				return false;

			if (std::fabs(pasa->xx_ - (*pasa->xp_)[pasa->MP]) < pasa->TINY && pasa->i_ == pasa->MP) {	// 完了したか？
				for (std::size_t i = 0; i < AdapStepHelper::NVAR; i++)
					(*pasa->yp_[i])[pasa->i_] = pasa->y_[i];

				return true;															// 正常終了
			}

			if (std::fabs(pasa->hnext_) <= pasa->HMIN) {
				std::cerr << "Step size too small in odeint" << std::endl;
				return false;
			}

			h = pasa->hnext_;
		}

		std::cerr << "Too many steps in routine odeint" << std::endl;
		return false;
	}
}