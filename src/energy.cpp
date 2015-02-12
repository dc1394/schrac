#include "energy.h"
#include "simpson.h"
#include <iostream>

namespace schrac {
    // #region コンストラクタ

    Energy::Energy(std::shared_ptr<DiffData> const & pdiffdata, dvector const & rf, dvector const & r) :
        pdiffdata_(pdiffdata),
        rf_(rf),
        potential_energy_(- Simpson(pdiffdata->dx_)(rf, rf, r, 2))
    {
    }

    //#endregion コンストラクタ
    
    // #region メンバ関数

    void Energy::express_energy() const
    {
        kinetic_energy();
        potential_energy();
        eigenvalue();
    }

    void Energy::eigenvalue() const
    {
        std::cout << "E(Eigenvalue)\t\t= " << pdiffdata_->E_ << std::endl;
    }

    void Energy::kinetic_energy() const
    {
        std::cout << "E(Kinetic Energy)\t= " << -0.5 * potential_energy_ << std::endl;
    }

    void Energy::potential_energy() const
    {
        std::cout << "E(Potential Energy)\t= " << potential_energy_ << std::endl;
    }

    // #endregion メンバ関数
}
