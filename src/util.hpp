#pragma once
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace util {

inline double wtime() {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
#endif
}

[[noreturn]] inline void die(const std::string& msg) {
    std::fprintf(stderr, "ERROR: %s\n", msg.c_str());
    std::exit(2);
}

inline std::string host() {
    const char* h = std::getenv("HOSTNAME");
    if (!h) h = std::getenv("COMPUTERNAME");
    return h ? h : "unknown";
}

inline int num_threads() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

}
