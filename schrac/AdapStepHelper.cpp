#include "AdapStepHelper.h"

namespace schrac {
	const long double AdapStepHelper::PGROW = - 0.2;
	const long double AdapStepHelper::PSHRNK = - 0.25;
	const long double AdapStepHelper::FCOR = 1.0 / 15.0;
	const long double AdapStepHelper::SAFETY = 0.9;
	const long double AdapStepHelper::ERRCON = std::pow(4.0 / SAFETY, 1.0 / PGROW);
}
