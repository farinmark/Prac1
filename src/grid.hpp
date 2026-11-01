#pragma once
#include <cmath>
#include <cstddef>
#include <vector>

struct Grid {
    int M, N;
    double h1, h2, ih1, ih2;
    double eps;
    int stride;
    inline int id(int i, int j) const { return i * stride + j; }
    size_t size() const { return (size_t)(M + 1) * (N + 1); }
};

// Коэффициенты разностной задачи.
struct Problem {
    std::vector<double> Q, dinv, b;
    double area_sum;
    long mixed_cells;
};

inline double v0(double x, double y) { return std::sin(x + y); }

Grid make_grid(int M, int N);

void assemble(const Grid& g, Problem& pb);
