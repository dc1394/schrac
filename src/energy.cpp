#include "energy.h"
#include "simpson.h"
#include <iostream>             // for std::cout

namespace schrac {
    // #region コンストラクタ

    Energy::Energy(std::shared_ptr<DiffData> const & pdiffdata, dvector const & r, dvector const & rf, double Z) :
        pdiffdata_(pdiffdata),
        rf_(rf),
        potcoulomb_energy_(- Z * Simpson(pdiffdata->dx_)(rf, rf, r, 2))
    {
    }

    //#endregion コンストラクタ
    
    // #region メンバ関数

    void Energy::express_energy(boost::optional<double> const & ehartree) const
    {
        kinetic_energy();
        if (ehartree) {
            coulomb_energy();
            hartree_energy(*ehartree);
        }
        else {
            potential_energy();
        }

        eigenvalue();
    }

    void Energy::coulomb_energy() const
    {
        std::cout << "E(Coulomb Energy)\t= " << potcoulomb_energy_ << std::endl;
    }

    void Energy::eigenvalue() const
    {
        std::cout << "E(Eigenvalue)\t\t= " << pdiffdata_->E_ << std::endl;
    }

    void Energy::hartree_energy(double ehartree) const
    {
        std::cout << "E(Hartree Energy)\t= " << ehartree << std::endl;
    }

    void Energy::kinetic_energy() const
    {
        std::cout << "E(Kinetic Energy)\t= " << -0.5 * potcoulomb_energy_ << std::endl;
    }

    void Energy::potential_energy() const
    {
        std::cout << "E(Potential Energy)\t= " << potcoulomb_energy_ << std::endl;
    }

    // #endregion メンバ関数
}
