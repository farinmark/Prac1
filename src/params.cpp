#include "params.hpp"
#include "util.hpp"

#include <cstdio>
#include <cstdlib>
#include <string>

static void usage() {
    std::printf(
        "usage: poisson_<serial|omp> [options]\n"
        "  -M <int>  -N <int>      число интервалов по x и y (по умолчанию 80 80)\n"
        "  --delta <real>          порог остановки ||w^{k+1}-w^k||_E < delta (1e-8)\n"
        "  --maxit <int>           лимит итераций (1000000)\n"
        "  --restart               перезапуск PCG при нарушении монотонности J\n"
        "  --history               напечатать историю итераций (строки hist: k step_E relres J)\n");
        
}

Params parse_args(int argc, char** argv) {
    Params p;
    for (int k = 1; k < argc; ++k) {
        const std::string a = argv[k];
        const bool has = k + 1 < argc;
        if (a == "-M" && has) p.M = std::atoi(argv[++k]);
        else if (a == "-N" && has) p.N = std::atoi(argv[++k]);
        else if (a == "--delta" && has) p.delta = std::atof(argv[++k]);
        else if (a == "--maxit" && has) p.maxit = std::atol(argv[++k]);
        else if (a == "--restart") p.restart = true;
        else if (a == "--history") p.history = true;
        else if (a == "-h" || a == "--help") { usage(); std::exit(0); }
        else { usage(); util::die("unknown option " + a); }
    }
    if (p.M < 2 || p.N < 2) util::die("M and N must be >= 2");
    return p;
}
