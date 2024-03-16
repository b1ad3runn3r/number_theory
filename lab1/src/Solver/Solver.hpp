#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace linalg = boost::numeric::ublas;

class Solver {
private:
    linalg::matrix<int> matrix;
    linalg::vector<int> values;

    static void swap_min(linalg::matrix_range<linalg::matrix<int>> &mr);
    static void gcd(linalg::matrix_range<linalg::matrix<int>> &mr);

    std::vector<linalg::vector<int>> general_solution;
    linalg::vector<int> partial_solution;

public:
    Solver(linalg::matrix<int> &new_matrix, linalg::vector<int> &new_values)
    : matrix(new_matrix), values(new_values) {};

    [[nodiscard]] const auto &getMatrix() const noexcept { return matrix; }
    [[nodiscard]] const auto &getValues() const noexcept { return values; }
    [[nodiscard]] const auto &getGeneral() const noexcept { return general_solution;}
    [[nodiscard]] const auto &getPartial() const noexcept { return partial_solution;}


    Solver &expand();
    Solver &convert();
    Solver &solve();

};


#endif // SOLVER_HPP
