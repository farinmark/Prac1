#include "grid.hpp"
#include "geometry.hpp"

#include <algorithm>

Grid make_grid(int M, int N) {
    Grid g;
    g.M = M; g.N = N;
    g.h1 = 10.0 / M; g.h2 = 6.0 / N;
    g.ih1 = 1.0 / (g.h1 * g.h1); g.ih2 = 1.0 / (g.h2 * g.h2);
    const double h = std::max(g.h1, g.h2);
    g.eps = h * h;
    g.stride = N + 1;
    return g;
}

void assemble(const Grid& g, Problem& pb) {
    const size_t L = g.size();
    pb.Q.assign(L, 0.0); pb.dinv.assign(L, 0.0); pb.b.assign(L, 0.0);
    double* Q = pb.Q.data(); double* Di = pb.dinv.data(); double* B = pb.b.data();

    const double h1 = g.h1, h2 = g.h2, cell = h1 * h2;
    const double qfict = 1.0 / (g.eps * g.eps);
    const double dlap = 2.0 * g.ih1 + 2.0 * g.ih2;
    double area = 0.0;
    long mixed = 0;
#pragma omp parallel for reduction(+:area,mixed) schedule(static)
    for (int i = 1; i < g.M; ++i) {
        const double x = -5.0 + i * h1;
        for (int j = 1; j < g.N; ++j) {
            const double y = -3.0 + j * h2;
            const double S = geom::cell_area(x - 0.5 * h1, x + 0.5 * h1, y - 0.5 * h2, y + 0.5 * h2);
            area += S;
            if (S > 0.0 && S < cell * (1.0 - 1e-12)) ++mixed;
            const double Qij = (cell - S) / cell * qfict;   // (h1h2 - S)/(h1h2 eps^2)
            const double Fij = S / cell;
            double bij = Fij + Qij * v0(x, y);
            if (i == 1)       bij += g.ih1 * v0(-5.0, y);
            if (i == g.M - 1) bij += g.ih1 * v0( 5.0, y);
            if (j == 1)       bij += g.ih2 * v0(x, -3.0);
            if (j == g.N - 1) bij += g.ih2 * v0(x,  3.0);
            const int k = g.id(i, j);
            Q[k] = Qij;
            Di[k] = 1.0 / (dlap + Qij);
            B[k] = bij;
        }
    }
    pb.area_sum = area;
    pb.mixed_cells = mixed;
}
