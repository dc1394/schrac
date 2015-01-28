#include "RungeKutta.h"

namespace schrac {
	bool RungeKutta::solve_diff_equ_O()
	{
		long double k1[2], k2[2], k3[2], k4[2];

		for (int i = 0; i < pdiffdata_->OSIZE - 1; i++) {
			pdiffdata_->node_count(i, pdiffdata_->LO);

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
		
			pdiffdata_->LO[i + 1] = pdiffdata_->LO[i] + (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0;
			pdiffdata_->MO[i + 1] = pdiffdata_->MO[i] + (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0;
		}

		return true;
	}

	bool RungeKutta::solve_diff_equ_I()
	{
		long double k1[2], k2[2], k3[2], k4[2];

		for (int i = 1; i < pdiffdata_->ISIZE - 1; i++) {
			pdiffdata_->node_count(i, pdiffdata_->LI);

			k1[0] = pdiffdata_->DX * dL_dx(pdiffdata_->MI[i]);
			k1[1] = pdiffdata_->DX * Diff::dM_dx(pdiffdata_->XV_I[i], pdiffdata_->LI[i],
									  pdiffdata_->MI[i], pdiffdata_);

			k2[0] = pdiffdata_->DX * dL_dx(pdiffdata_->MI[i] - k1[1] / 2.0);
			k2[1] = pdiffdata_->DX * Diff::dM_dx(pdiffdata_->XV_I[i] - pdiffdata_->DX / 2.0,
									  pdiffdata_->LI[i] - k1[0] / 2.0,
									  pdiffdata_->MI[i] - k1[1] / 2.0,
									  pdiffdata_);

			k3[0] = pdiffdata_->DX * dL_dx(pdiffdata_->MI[i] - k2[1] / 2.0);
			k3[1] = pdiffdata_->DX * Diff::dM_dx(pdiffdata_->XV_I[i] - pdiffdata_->DX / 2.0,
									  pdiffdata_->LI[i] - k2[0] / 2.0,
									  pdiffdata_->MI[i] - k2[1] / 2.0,
									  pdiffdata_);

			k4[0] = pdiffdata_->DX * dL_dx(pdiffdata_->MI[i] - k3[1]);
			k4[1] = pdiffdata_->DX * Diff::dM_dx(pdiffdata_->XV_I[i] - pdiffdata_->DX,
									  pdiffdata_->LI[i] - k3[0],
									  pdiffdata_->MI[i] - k3[1],
									  pdiffdata_);
		
			pdiffdata_->LI[i + 1] = pdiffdata_->LI[i] - (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) / 6.0;
			pdiffdata_->MI[i + 1] = pdiffdata_->MI[i] - (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) / 6.0;
		}

		return true;
	}
}