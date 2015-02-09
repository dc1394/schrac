#include "energy.h"
#include "simpson.h"
#include <iostream>

namespace schrac {
    Energy::Energy(std::shared_ptr<DiffData> const & pdiffdata, dvector const & rf, dvector const & r) :
        pdiffdata_(pdiffdata),
        rf_(rf),
        potential_energy_(Simpson(pdiffdata->dx_)(rf, rf, r, 2))
    {
    }

    void Energy::eigenenergy() const
    {
        std::cout << "E(Eigenenergy) = " << pdiffdata_->E_ << std::endl;
    }

    void Energy::kinetic_energy() const
    {
        std::cout << "E(Kinetic Energy) = " << -0.5 * potential_energy_ << std::endl;
    }

    void Energy::potential_energy() const
    {
        std::cout << "E(Potential Energy) = " << potential_energy_ << std::endl;
    }
}
