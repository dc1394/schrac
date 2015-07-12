/*! \file solvelinearequ.cpp
    \brief 連立一次方程式を解く関数の実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#include "solvelinearequ.h"
#include <cstdint>          // for std::int32_t
#include <memory>           // for std::unique_ptr
#include <stdexcept>        // for std::runtime_error
#include <string>           // for std::to_string

namespace schrac {
    myvector solve_linear_equ(std::array<double, AMMAX * AMMAX> & a, myvector & b)
    {
        // save original handler, install new handler
        auto old_handler = gsl_set_error_handler(
            [](char const * reason, char const * file, std::int32_t line, std::int32_t)
        {
            auto const str = std::string(reason) + "\nFile: " + file + "\nline: " + std::to_string(line);
            throw std::runtime_error(str);
        });

        auto m = gsl_matrix_view_array(a.data(), AMMAX, AMMAX);
        auto const v = gsl_vector_view_array(b.data(), AMMAX);

        auto const gsl_vector_deleter = [](gsl_vector * p)
        {
            gsl_vector_free(p);
        };
        std::unique_ptr<gsl_vector, decltype(gsl_vector_deleter)> x(
            gsl_vector_alloc(AMMAX),
            gsl_vector_deleter);

        auto const gsl_perm_deleter = [](gsl_permutation * p)
        {
            gsl_permutation_free(p);
        };
        std::unique_ptr<gsl_permutation, decltype(gsl_perm_deleter)> p(
            gsl_permutation_alloc(AMMAX),
            gsl_perm_deleter);

        std::int32_t s;
        gsl_linalg_LU_decomp(&m.matrix, p.get(), &s);
        gsl_linalg_LU_solve(&m.matrix, p.get(), &v.vector, x.get());
        myvector solution{ 0.0, 0.0, 0.0 };
        std::copy(x->data, x->data + AMMAX, solution.begin());

        // restore original handler
        gsl_set_error_handler(old_handler);

        return std::move(solution);
    }
}
