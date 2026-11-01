#pragma once
#include <cmath>

namespace geom {

struct Pt { double x, y; };

static const int D_NV = 8;
static const Pt D_POLY[D_NV] = {
    { 5.0, 0.0}, { 2.0, 3.0}, { 0.0, 1.0}, {-2.0, 3.0},
    {-5.0, 0.0}, {-2.0,-3.0}, { 0.0,-1.0}, { 2.0,-3.0}
};
static const double D_AREA = 34.0;

inline bool inD(double x, double y) {
    return std::fabs(std::fabs(x) - 2.0) + std::fabs(y) < 3.0;
}

static const int CLIP_MAX = 64;

inline int clip_half(const Pt* in, int n, Pt* out, int axis, double c, double s) {
    int m = 0;
    for (int k = 0; k < n; ++k) {
        const Pt& P = in[k];
        const Pt& R = in[(k + 1) % n];
        const double fp = s * ((axis == 0 ? P.x : P.y) - c);
        const double fr = s * ((axis == 0 ? R.x : R.y) - c);
        if (fp >= 0.0) out[m++] = P;
        if ((fp >= 0.0) != (fr >= 0.0)) {
            const double t = fp / (fp - fr);
            Pt X;
            X.x = P.x + t * (R.x - P.x);
            X.y = P.y + t * (R.y - P.y);
            if (axis == 0) X.x = c; else X.y = c;
            out[m++] = X;
        }
    }
    return m;
}

inline double poly_area(const Pt* p, int n) {
    double a = 0.0;
    for (int k = 0; k < n; ++k) {
        const Pt& P = p[k];
        const Pt& R = p[(k + 1) % n];
        a += P.x * R.y - R.x * P.y;
    }
    return 0.5 * a;
}

inline double cell_area(double x0, double x1, double y0, double y1) {
    if (x1 <= -5.0 || x0 >= 5.0 || y1 <= -3.0 || y0 >= 3.0) return 0.0;
    Pt a[CLIP_MAX], b[CLIP_MAX];
    int n = D_NV;
    for (int k = 0; k < D_NV; ++k) a[k] = D_POLY[k];
    n = clip_half(a, n, b, 0, x0, +1.0); if (n == 0) return 0.0;
    n = clip_half(b, n, a, 0, x1, -1.0); if (n == 0) return 0.0;
    n = clip_half(a, n, b, 1, y0, +1.0); if (n == 0) return 0.0;
    n = clip_half(b, n, a, 1, y1, -1.0); if (n == 0) return 0.0;
    const double s = poly_area(a, n);
    return s > 0.0 ? s : 0.0;
}

}
