#include <iostream>
#include "Solver/Solver.hpp"
#include <boost/numeric/ublas/io.hpp>

int main() {
    size_t m, n;
    std::cin >> m >> n;
    linalg::matrix<int> matrix(m, n);

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            std::cin >> matrix(i,j);
        }
    }

    linalg::vector<int> vector(m);
    for (auto &i : vector) {
        std::cin >> i;
    }
    auto solver = Solver(matrix, vector);
    try {
        solver = solver
                .expand()
                .convert()
                .solve();
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    std::cout << solver.getMatrix() << std::endl;
    for (const auto &sol : solver.getGeneral()) {
        std::cout << sol << " + ";
    }
    std::cout << std::endl;
    std::cout << solver.getPartial() << std::endl;

    return 0;
}
