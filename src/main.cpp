#include <cmath>
#include <cstdio>
#include <vector>

#include "grid.hpp"
#include "params.hpp"
#include "pcg.hpp"
#include "util.hpp"

int main(int argc, char** argv) {
    const Params P = parse_args(argc, argv);

    const double T0 = util::wtime();
    const Grid g = make_grid(P.M, P.N);
    Problem pb;
    assemble(g, pb);
    Result res;
    res.t_setup = util::wtime() - T0;

    const double T1 = util::wtime();
    std::vector<double> w;
    pcg_solve(g, pb, P, w, res, nullptr);
    res.t_solve = util::wtime() - T1;

    res.relres_true = true_residual(g, pb, w);
    if (res.status == "CONVERGED" && !std::isfinite(res.relres_true)) {
        res.status = "NAN"; res.reason = "nan_in_solution";
    }
    res.t_total = util::wtime() - T0;

#ifdef _OPENMP
    const char* code = "omp";
#else
    const char* code = "serial";
#endif
    std::printf("code=%s M=%d N=%d threads=%d host=%s\n",
                code, g.M, g.N, util::num_threads(), util::host().c_str());
    std::printf("stop: delta=%.3e status=%s reason=%s iters=%ld\n",
                P.delta, res.status.c_str(), res.reason.c_str(), res.iters);
    std::printf("norms: step_E=%.6e relres_final=%.6e relres_true=%.6e\n",
                res.step, res.relres, res.relres_true);
    std::printf("functional: J_final=%.15e J_violations=%ld\n", res.J, res.J_viol);
    std::printf("time: setup=%.6f solve=%.6f total=%.6f (s)\n",
                res.t_setup, res.t_solve, res.t_total);
    return res.status == "CONVERGED" ? 0 : 1;
}
