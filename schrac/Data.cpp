#include "Data.h"

namespace HydroSchDirac {
	const long double Data::c = 137.035999;
	const long double Data::al = 1.0 / c;
	const long double Data::al2half = 0.5 * al * al;

	const ci_string Data::ALPHA("ALPHA");
	const ci_string Data::BETA("BETA");
	const array<const std::string, 6> Data::atomName = { "H", "He", "Li", "Be", "B", "C" };

	const long double Data::XMIN_DEFAULT = - 7.0;
	const long double Data::XMAX_DEFAULT = 5.0;
	const long double Data::EPS_DEFAULT = 1.0E-15;
	const long double Data::MAT_PO_RATIO_DEFAULT = 0.67;
}
