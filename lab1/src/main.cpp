#include <iostream>
#include "Solver/Solver.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <format>

int main() {
    size_t m, n;
    std::cout << "Enter number of rows: " << std::endl;
    std::cin >> m;

    std::cout << "Enter number of columns: " << std::endl;
    std::cin >> n;

    linalg::matrix<int> matrix(m, n);

    std::cout << "Enter matrix coeffs: " << std::endl;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            std::cin >> matrix(i,j);
        }
    }

    std::cout << "Enter values: " << std::endl;

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

    auto partial = solver.getPartial();
    auto general = solver.getGeneral();

    if (partial.empty()) {
        std::cout << "No solution :(" << std::endl;
        return 1;
    }

    std::cout << "Solution: " << std::endl;
    for (size_t i = 0; i < general.size(); ++i) {
        std::cout << std::format("t{} * ", i + 1) << general[i] << " + ";
    }

    std::cout << partial << std::endl;


    return 0;
}
