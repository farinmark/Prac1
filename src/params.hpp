#pragma once

struct Params {
    int M = 80, N = 80;
    double delta = 1e-8;
    long maxit = 1000000;
    int hist_every = 1;
    bool restart = false;
};

Params parse_args(int argc, char** argv);
