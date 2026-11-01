#include "pcg.hpp"

#include <cmath>

static void apply_A(const Grid& g, const std::vector<double>& Q, const std::vector<double>& u,
                    std::vector<double>& q, double* s) {
    const int M = g.M, N = g.N, st = g.stride;
    const double ih1 = g.ih1, ih2 = g.ih2;
    const double* U = u.data();
    const double* QQ = Q.data();
    double* R = q.data();
    double uq = 0.0, uu = 0.0;
#pragma omp parallel for reduction(+:uq,uu) schedule(static)
    for (int i = 1; i < M; ++i) {
        const int base = i * st;
        for (int j = 1; j < N; ++j) {
            const int k = base + j;
            const double c = U[k];
            const double v = (2.0 * c - U[k - st] - U[k + st]) * ih1
                           + (2.0 * c - U[k - 1] - U[k + 1]) * ih2
                           + QQ[k] * c;
            R[k] = v;
            uq += c * v;
            uu += c * c;
        }
    }
    s[0] += uq; s[1] += uu;
}

void pcg_solve(const Grid& g, const Problem& pb, const Params& P,
               std::vector<double>& w, Result& res, std::vector<double>* hist) {
    const size_t L = g.size();
    const int M = g.M, N = g.N, st = g.stride;
    std::vector<double> r(L, 0.0), z(L, 0.0), p(L, 0.0), q(L, 0.0);
    w.assign(L, 0.0);
    const double cellE = g.h1 * g.h2;   // ||u||_E^2 = h1 h2 sum u^2

    double bnorm2 = 0, rho = 0;
    {
        // w = 0, r = b, z = D^{-1} r, p = z
        double* R = r.data(); double* Z = z.data(); double* PP = p.data();
        const double* B = pb.b.data(); const double* Di = pb.dinv.data();
#pragma omp parallel for reduction(+:bnorm2,rho) schedule(static)
        for (int i = 1; i < M; ++i)
            for (int j = 1; j < N; ++j) {
                const int k = i * st + j;
                R[k] = B[k];
                Z[k] = Di[k] * R[k];
                PP[k] = Z[k];
                bnorm2 += B[k] * B[k];
                rho += R[k] * Z[k];
            }
    }
    const double bnorm = std::sqrt(bnorm2);
    res.relres0 = 1.0;
    res.J0 = 0.0;
    double Jprev = 0.0;
    long k = 0;
    if (!(bnorm > 0.0)) { res.status = "CONVERGED"; res.reason = "zero_rhs"; }

    while (res.status == "RUNNING") {
        if (k >= P.maxit) { res.status = "MAXIT"; res.reason = "iteration_limit"; break; }
        ++k;
        double s1[2] = {0, 0};
        apply_A(g, pb.Q, p, q, s1);
        const double pq = s1[0], pp = s1[1];
        if (!std::isfinite(pq) || !std::isfinite(pp)) { res.status = "NAN"; res.reason = "nan_in_pAp"; break; }
        if (pq <= 0.0) { res.status = "BREAKDOWN"; res.reason = "pAp_nonpositive"; break; }
        const double alpha = rho / pq;

        double rz = 0, rr = 0, bw = 0, rw = 0;
        {
            double* W = w.data(); double* R = r.data(); double* Z = z.data();
            const double* PP = p.data(); const double* QQ = q.data();
            const double* B = pb.b.data(); const double* Di = pb.dinv.data();
#pragma omp parallel for reduction(+:rz,rr,bw,rw) schedule(static)
            for (int i = 1; i < M; ++i)
                for (int j = 1; j < N; ++j) {
                    const int kk = i * st + j;
                    const double wn = W[kk] + alpha * PP[kk];
                    const double rn = R[kk] - alpha * QQ[kk];
                    const double zn = Di[kk] * rn;
                    W[kk] = wn; R[kk] = rn; Z[kk] = zn;
                    rz += rn * zn; rr += rn * rn;
                    bw += B[kk] * wn; rw += rn * wn;
                }
        }

        const double step = std::fabs(alpha) * std::sqrt(cellE * pp);
        const double relres = std::sqrt(rr) / bnorm;
        const double J = -0.5 * (bw + rw);
        res.iters = k; res.step = step; res.relres = relres; res.J = J;

        if (!std::isfinite(step) || !std::isfinite(relres) || !std::isfinite(J)) {
            res.status = "NAN"; res.reason = "nan_in_iteration"; break;
        }
        bool restart_now = false;
        const double incr = J - Jprev;
        const double tolJ = 1e-12 * (std::fabs(J) + std::fabs(bw)) + 1e-300;
        if (k > 1 && incr > tolJ) {
            ++res.J_viol;
            if (incr > res.J_max_incr) res.J_max_incr = incr;
            if (P.restart) { restart_now = true; ++res.restarts; }
        }
        Jprev = J;

        if (step < P.delta) { res.status = "CONVERGED"; res.reason = "step_E_below_delta"; }
        const double beta = restart_now ? 0.0 : rz / rho;
        if (hist && (k % P.hist_every == 0 || res.status != "RUNNING")) {
            const double row[6] = {(double)k, step, relres, J, alpha, beta};
            hist->insert(hist->end(), row, row + 6);
        }
        if (res.status != "RUNNING") break;
        rho = rz;
        {
            double* PP = p.data(); const double* Z = z.data();
#pragma omp parallel for schedule(static)
            for (int i = 1; i < M; ++i)
                for (int j = 1; j < N; ++j) {
                    const int kk = i * st + j;
                    PP[kk] = Z[kk] + beta * PP[kk];
                }
        }
    }
}

double true_residual(const Grid& g, const Problem& pb, const std::vector<double>& w) {
    const int M = g.M, N = g.N, st = g.stride;
    std::vector<double> q(g.size(), 0.0);
    double s[2] = {0, 0};
    apply_A(g, pb.Q, w, q, s);
    double tr = 0, bb = 0;
    const double* B = pb.b.data(); const double* QQ = q.data();
#pragma omp parallel for reduction(+:tr,bb) schedule(static)
    for (int i = 1; i < M; ++i)
        for (int j = 1; j < N; ++j) {
            const int kk = i * st + j;
            const double d = B[kk] - QQ[kk];
            tr += d * d;
            bb += B[kk] * B[kk];
        }
    return bb > 0 ? std::sqrt(tr / bb) : 0.0;
}
