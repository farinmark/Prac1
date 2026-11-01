#pragma once
#include <string>
#include <vector>

#include "grid.hpp"
#include "params.hpp"

struct Result {
    long iters = 0;
    double step = 0, relres0 = 0, relres = 0, relres_true = 0, J = 0, J0 = 0;
    long J_viol = 0, restarts = 0;
    double J_max_incr = 0;
    std::string status = "RUNNING", reason;
    double t_setup = 0, t_solve = 0, t_total = 0;
};

// Решает A w = b с w0 = 0. Остановка: ||w^{k+1}-w^k||_E < P.delta.
// Если hist != nullptr, туда пишутся строки истории
// (по 6 чисел: k, step_E, relres, J, alpha, beta) с шагом P.hist_every.
void pcg_solve(const Grid& g, const Problem& pb, const Params& P,
               std::vector<double>& w, Result& res, std::vector<double>* hist);

// Истинная относительная невязка ||b - A w|| / ||b||.
double true_residual(const Grid& g, const Problem& pb, const std::vector<double>& w);
