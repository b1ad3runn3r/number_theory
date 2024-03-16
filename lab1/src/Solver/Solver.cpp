#include "Solver.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <algorithm>
#include <ranges>
// #include <boost/numeric/ublas/io.hpp>

void Solver::swap_min(linalg::matrix_range<linalg::matrix<int>> &mr) {
    linalg::matrix_row row(mr, 0);
    try {
        size_t min_idx = 0;
        int min = std::numeric_limits<int>::max();
        for (size_t i = 0; i < row.size(); ++i) {
            if (std::abs(row[i]) < std::abs(min) && row[i] != 0) {
                min_idx = i;
                min = row[i];
            }
        }

        for (size_t i = 0; i < mr.size1(); ++i) {
            std::swap(mr(i, 0), mr(i, min_idx));
        }
    }
    catch (...) {
        throw;
    }
}

void Solver::gcd(linalg::matrix_range<linalg::matrix<int>> &mr) {
    auto q = mr(0, 0);

    for (size_t i = 1; i < mr.size2(); ++i) {
        auto k = mr(0, i);
        auto d = (k - k % q) / q;
        linalg::column(mr, i) -= d * linalg::column(mr, 0);
    }
}

Solver &Solver::expand() {
    auto m = matrix.size1();
    auto n = matrix.size2();

    matrix.resize(m + n, n + 1);

    linalg::vector<int> col(m + n, 0);
    linalg::vector_range vr(col, linalg::range(0, m));
    vr = values;

    linalg::column(matrix, n) = col;

    linalg::matrix_range sub_matrix(matrix, linalg::range(m, m + n), linalg::range(0, n ));
    sub_matrix = linalg::identity_matrix<int>(n);


    return *this;
}

Solver &Solver::convert() {
    auto m = matrix.size1();
    auto n = matrix.size2();
    for (size_t i = 0; i < matrix.size1() - matrix.size2(); ++i) {
        linalg::matrix_range sub_matrix(matrix, linalg::range(i, matrix.size1()), linalg::range(i, matrix.size2()));
        linalg::matrix_row row(matrix, 0);
        linalg::vector_range vr(row, linalg::range(1, row.size()));
        while (!std::all_of(vr.begin(), vr.end(), [](auto i) {return i == 0; })) {
            swap_min(sub_matrix);
            gcd(sub_matrix);
        }
    }

    return *this;
}

Solver &Solver::solve() {
    linalg::matrix_column col(matrix, matrix.size2() - 1);
    linalg::vector<int> b_exp(matrix.size1());
    linalg::vector_range vr(col, linalg::range(0, matrix.size1() - matrix.size2() + 1));
    vr = values;
    vr *= -1;
    for (size_t i = 0; i < matrix.size1(); ++i) {
        b_exp[i] = vr[i];
    }

    auto m = matrix.size1() - matrix.size2();
    for (size_t i = 0; i < m; ++i) {
        auto a = b_exp[i];
        auto q = matrix (i,i);

        if (q == 0 && a != 0) {
            return *this;
        }

        if (q != 0) {
            auto k = (a - a % q) / q;
            b_exp += -k * linalg::column(matrix, i);
        }
    }

    linalg::vector_range vrr(b_exp, linalg::range(0, m));
    if (std::all_of(vrr.begin(), vrr.end(), [](auto a){return a == 0;})) {
        auto t = linalg::vector_range(b_exp, linalg::range(m, b_exp.size()));
        if (std::all_of(t.begin(), t.end(), [](auto a){ return a == 0; })) {
            linalg::vector<int> f(t.size(),0);
            general_solution.push_back(f);
        }

        // par = Some(t)
        for (size_t i = m; i < matrix.size2(); ++i) {
            linalg::matrix_range mr(matrix, linalg::range(m, matrix.size1()), linalg::range(i, i + 1));
            general_solution.push_back(linalg::matrix_column(mr, 0));
        }
    }

    return *this;
}
