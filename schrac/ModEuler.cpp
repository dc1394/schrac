#include "ModEuler.h"

#ifdef _DEBUG
namespace HydroSchDirac {
	bool ModEuler::solve_diff_equ_O()
	{
		for (int i = 0; i < pdiffdata_->OSIZE - 1; i++) {
			pdiffdata_->node_count(i, pdiffdata_->LO);

			const long double r0 = pdiffdata_->RV_O[i];
			const long double f0 = - (2.0 * static_cast<const long double>(pdata_->l) + 1.0) * pdiffdata_->MO[i] +
							  2.0 * r0 * r0 * (pdiffdata_->VP_O[i] - pdiffdata_->E_) * pdiffdata_->LO[i];    

			const long double Mtmp = pdiffdata_->MO[i] + pdiffdata_->DX * f0;
			pdiffdata_->LO[i + 1] = pdiffdata_->LO[i] + 0.5 * pdiffdata_->DX * (pdiffdata_->MO[i] + Mtmp);

			const long double r1 = pdiffdata_->RV_O[i + 1];
			const long double f1 = - (2.0 * static_cast<const long double>(pdata_->l) + 1.0) * Mtmp +
								  2.0 * r1 * r1 * (pdiffdata_->VP_O[i + 1] - pdiffdata_->E_) * pdiffdata_->LO[i + 1];
			pdiffdata_->MO[i + 1] = pdiffdata_->MO[i] + 0.5 * pdiffdata_->DX * (f0 + f1);
		}

		return true;
	}

	bool ModEuler::solve_diff_equ_I()
	{
		for (int i = 1; i < pdiffdata_->ISIZE - 1; i++) {
			pdiffdata_->node_count(i, pdiffdata_->LI);

			const long double r0 = pdiffdata_->RV_I[i];
			const long double f0 = - (2.0 * static_cast<const long double>(pdata_->l) + 1.0) * pdiffdata_->MI[i] +
							  2.0 * r0 * r0 * (pdiffdata_->VP_I[i] - pdiffdata_->E_) * pdiffdata_->LI[i];    
			const long double Mtmp = pdiffdata_->MI[i] + pdiffdata_->DX * f0;
			pdiffdata_->LI[i + 1] = pdiffdata_->LI[i] - 0.5 * pdiffdata_->DX * (pdiffdata_->MI[i] + Mtmp);

			const long double r1 = pdiffdata_->RV_I[i + 1];
			const long double f1 = - (2.0 * static_cast<const long double>(pdata_->l) + 1.0) * Mtmp +
								  2.0 * r1 * r1 * (pdiffdata_->VP_I[i + 1] - pdiffdata_->E_) * pdiffdata_->LI[i + 1];
			pdiffdata_->MI[i + 1] = pdiffdata_->MI[i] - 0.5 * pdiffdata_->DX * (f0 + f1);
		}

		return true;
	}
}
#endif
