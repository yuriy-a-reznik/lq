/*!
 *  \file lq.c
 *
 *  \brief LQ: a very basic lattice quantization library.
 *
 *  This library implements lattice quantization algorithms for lattices LQ_Zn, LQ_An, LQ_Dn, LQ_En and their duals.
 *  Dimensions up to 8 are supported. 
 *
 *  This includes best known lattices for minimizing covering radius:
 *
 *      1	Z	        Integer lattice
 *      2	A2* ~ A2	Hexagonal lattice
 *      3	A3* ~ D3*	BCC / FCC dual
 *      4	A4*	        A4 dual
 *      5	A5*	        A5 dual
 *      6	A6*	        A6 dual
 *      7	A7*	        A7 dual
 *      8	E8	        Gossett lattice
 *
 *  and best lattices for minimizing quantization error:
 *
 *      1	Z	        Integer lattice
 *      2	A2	        Hexagonal lattice
 *      3	A3* ~ D3*	BCC / FCC dual
 *      4	D4	        4D checkerboard lattice, Hurwitz integers
 *      5	D5*	        D5 dual
 *      6	E6*	        E6 dual
 *      7	E7*	        E7 dual
 *      8	E8	        Gossett lattice
 *
 *  In all cases, this library uses the Voronoi region scaled by a factor r >= 2 as the granular region.
 *  The resulting number of points is r^n, and the effective rate of the quantizer is n·log2(r) bits.
 *  The distortion is upper-bounded by the covering radius for points within the granular region.
 *  Points outside the granular region are clipped to the nearest point within it, with the overload flag set to 1.
 *  The overload region is unbounded.
 *
 *  Main references:
 *     [1] J.H. Conway and N.J.A. Sloane, "Fast quantizing and decoding algorithms for lattice quantizers and codes,"
 *         IEEE Trans. Inform. Theory, vol. IT-28, no. 2, pp. 227-232, March 1982.
 *     [2] J.H. Conway and N.J.A. Sloane, "A fast encoding method for lattice codes and quantizers,"
 *         IEEE Trans. Inform. Theory, vol. IT-29, no. 6, pp. 820-824, November 1983.
 *     [3] K. Takizawa, H. Yagi and T. Kawabata, "Closest point algorithms with lp norm for root lattices,"
 *         2010 IEEE International Symposium on Information Theory, Austin, TX, USA, 2010, pp. 1042-1046.
 *     [4] J.H. Conway and N.J.A. Sloane, "Sphere packings, lattices and groups," Vol. 290, Springer, 2013.
 *
 *  Implementation notes:  
 * 
 *    All library API functions operate with n-dimensional vectors. 
 *    The extensions to ambient lattice spaces are handled internally. 
 * 
 *    LQ_Zn, LQ_Dn, Dn_star, E8:  ambient = n,   no padding.
 *    LQ_An, An_star:             ambient = n+1, pad with x[n] = -sum(x[0..n-1]).
 *    E7, E7_star:                ambient = 8,   pad with x[7] = -sum(x[0..6]).
 *    E6, E6_star:                ambient = 8,   pad with x[6] = -sum(x[1..5]), x[7] = -x[0]. (constraints: total sum=0, x[0]+x[7]=0)
 * 
 *    The generator matrices G and their pseudo-inverses Ginv are pre-computed for all lattices and dimensions, and stored in static tables.
 *    The center points (a) for all lattices, dimensions, and granularity parameters (r) are also pre-computed and stored in static tables.
 * 
 *    The clipping of vectors in overload region is done as part of quantization function, returning overload flag set to 1. 
 * 
 *    Otherwise, the implementation follow closely the algorithms described in the references [1-3].
 * 
 *    Pure C89, preferring clarity and simplicity over performance optimizations.
 * 
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.02
 *  \date    April 27, 2026
 */

#include <string.h>
#include <limits.h>
#include <math.h>
#include "lq.h"

/*!**************
 *
 *  Rounding and other helper functions utilized in finding closest lattice point algorithms [1]:
 */

/*! Round to nearest integer, breaking ties towards zero (smallest absolute value). */
static double f(double x)
{
    double lo = floor(x), hi = lo + 1.0;
    double dl = x - lo, dh = hi - x;
    if (dl < dh) return lo;
    if (dh < dl) return hi;
    return (fabs(lo) <= fabs(hi)) ? lo : hi;
}

/*! Round the "wrong way", breaking ties away from zero (largest absolute value). */
static double w(double x)
{
    double lo = floor(x), hi = lo + 1.0;
    double dl = x - lo, dh = hi - x;
    if (dl < dh) return hi;
    if (dh < dl) return lo;
    return (fabs(lo) <= fabs(hi)) ? hi : lo; /* round away from zero */
}

/*! Vector version of f() */
static void fv(const double* x, int n, double* out)
{
    int i;
    for (i = 0; i < n; i++) out[i] = f(x[i]);
}

/*! The same as fv(), except that component with farthest distance from integer is rounded wrong way: */
static void gv(const double* x, int n, double* out)
{
    int i, k = 0;
    double md, d;
    /* compute f(x) */
    fv(x, n, out);
    /* find worst component: */
    md = fabs(x[0] - out[0]);
    for (i = 1; i < n; i++) {
        d = fabs(x[i] - f(x[i]));
        if (d > md) { md = d; k = i; }
    }
    /* round it wrong way: */
    out[k] = w(x[k]);
}

/*! Integer sum of rounded vector components */
static int isum(const double* v, int n)
{
    int s = 0, i;
    for (i = 0; i < n; i++)
        s += (int)(v[i] + (v[i] >= 0 ? 0.5 : -0.5));
    return s;
}

/*! Squared distance between two vectors */
static double dist2(const double* a, const double* b, int n)
{
    double s = 0, d; int i;
    for (i = 0; i < n; i++) { d = a[i] - b[i]; s += d * d; }
    return s;
}

/*! Modulo operation that always returns non-negative result */
static int mod_pos(int a, int m)
{
    int r = a % m;
    return r < 0 ? r + m : r;
}

/*! Delta structure for sorting */
typedef struct { int idx; double val; } delta_t;

/*! Comparator for delta_t: sort by val ascending, break ties by idx ascending */
static int cmp_delta(const void* a, const void* b)
{
    const delta_t* p = (const delta_t*)a;
    const delta_t* q = (const delta_t*)b;
    if (p->val < q->val) return -1;
    if (p->val > q->val) return  1;
    return (p->idx < q->idx) ? -1 : (p->idx > q->idx) ? 1 : 0;
}


/*************** 
 * 
 *  Closest-point algorithms.
 * 
 *    Implemented according to [1] & [3] (lattices E6, E7, and their duals).
 */

/*! Closest-point algorithm for Z^n lattice [1] */
static int closest_Zn(const double* x, int n, double* out)
{
    int i;

	/* check arguments */
    if (!x || !out || n < 1 || n > 8) return LQ_INVARG;

	/* round each component to nearest integer: */
    for (i = 0; i < n; i++) out[i] = f(x[i]);

    return LQ_SUCCESS;
}

/*! Closest-point algorithm for D^n lattice [1] */
static int closest_Dn(const double* x, int n, double* out)
{
    double fx[LQ_MAX_AMB], gx[LQ_MAX_AMB];

	/* check arguments */
    if (!x || !out || n < 2 || n > LQ_MAX_DIM) return LQ_INVARG;

	/* find closest point in LQ_Zn, and in LQ_Zn shifted by (0.5,...,0.5), and pick the best: */
    fv(x, n, fx);
    gv(x, n, gx);
    memcpy(out, (isum(fx, n) & 1) == 0 ? fx : gx, n * sizeof(double));

    return LQ_SUCCESS;
}

/*! Closest-point algorithm for D^n* lattice [1] */
static int closest_Dn_star(const double* x, int n, double* out)
{
    double y0[LQ_MAX_AMB], y1[LQ_MAX_AMB], tmp[LQ_MAX_AMB];
    int i;

	/* check arguments */
    if (!x || !out || n < 2 || n > LQ_MAX_DIM) return LQ_INVARG;

	/* find closest point in LQ_Dn, and in LQ_Dn shifted by (0.5,...,0.5), and pick the best: */
    fv(x, n, y0);
    for (i = 0; i < n; i++) tmp[i] = x[i] - 0.5;
    fv(tmp, n, y1);
    for (i = 0; i < n; i++) y1[i] += 0.5;
    memcpy(out, dist2(x, y0, n) <= dist2(x, y1, n) ? y0 : y1, n * sizeof(double));
    
    return LQ_SUCCESS;
}

/*! Closest-point algorithm for A^n lattice [1] */
static int closest_An(const double* x, int n, double* out)
{
    int np1, i, def, absdef;
    double s, xp[LQ_MAX_AMB];
    delta_t sd[LQ_MAX_AMB];

	/* check arguments */
    if (!x || !out || n < 1 || n > LQ_MAX_DIM) return LQ_INVARG;

	/* compute centered coordinates: xp[i] = x[i] - mean(x) */
    np1 = n + 1;
    s = 0.0;
    for (i = 0; i < np1; i++) s += x[i];
    s /= (double)np1;
    for (i = 0; i < np1; i++) xp[i] = x[i] - s;

	/* round to nearest integer, compute rounding error, and sum of rounded values (defect): */
    def = 0;
    for (i = 0; i < np1; i++) {
        out[i] = f(xp[i]);
        def += (int)(out[i] + (out[i] >= 0 ? 0.5 : -0.5));
        sd[i].idx = i;
        sd[i].val = xp[i] - out[i];
    }

    /* If defect != 0, we need to round some components the wrong way. 
	 * Sort components by rounding error, and pick the ones with largest error: */
    if (def != 0) {
        qsort(sd, (size_t)np1, sizeof(delta_t), cmp_delta);
        absdef = def > 0 ? def : -def;
        if (absdef > np1) absdef = np1;
        if (def > 0)
            for (i = 0; i < absdef; i++) out[sd[np1 - 1 - i].idx] -= 1.0;
        else
            for (i = 0; i < absdef; i++) out[sd[i].idx] += 1.0;
    }

    return LQ_SUCCESS;
}

/*! Closest-point algorithm for A^n* lattice [1] */
static int closest_An_star(const double* x, int n, double* out)
{
    int np1, i, j, cnt;
    double best, d, xr[LQ_MAX_AMB], yi[LQ_MAX_AMB], ri[LQ_MAX_AMB];

    /* check arguments */
    if (!x || !out || n < 1 || n > LQ_MAX_DIM) return LQ_INVARG;

    /* try all n+1 cosets of LQ_An in LQ_An*, find the closest point in each, and pick the best: */
    np1 = n + 1;
    best = 1e30;
    for (i = 0; i <= n; i++) {
        cnt = np1 - i;
        for (j = 0; j < cnt; j++)  ri[j] = (double)i / (double)np1;
        for (j = cnt; j < np1; j++) ri[j] = -(double)cnt / (double)np1;
        for (j = 0; j < np1; j++) xr[j] = x[j] - ri[j];
        if (closest_An(xr, n, yi) != LQ_SUCCESS) return LQ_INVARG;
        for (j = 0; j < np1; j++) yi[j] += ri[j];
        d = dist2(x, yi, np1);
        if (d < best) { best = d; memcpy(out, yi, np1 * sizeof(double)); }
    }

    return LQ_SUCCESS;
}

/*! Closest-point algorithm for E6 lattice [3] */
static int closest_E6(const double* x, int n, double* out)
{
    static const double d1[2] = { -0.5, 0.5 };
    static const double d2[6] = { -0.5,-0.5,-0.5, 0.5, 0.5, 0.5 };
    double x1[2], x2[6], k1a[2], k1b[2], k2a[6], k2b[6];
    double c0[8], c1[8], xr2[2], xr6[6];
    int i;

	/* check arguments */
    if (!x || !out || n!=6) return LQ_INVARG;

    /* 
     * E6 can be represented as the intersection of two LQ_An lattices: 
     *     E6 = {x in R^8 : sum(x)=0, x[0]+x[7]=0, and x[1..6] sum to 0} 
     * We find the closest point in each LQ_An, and then pick the best. 
     */
    x1[0] = x[0]; x1[1] = x[7];
    for (i = 0; i < 6; i++) x2[i] = x[1 + i];

    if (closest_An(x1, 1, k1a) != LQ_SUCCESS) return LQ_INVARG;
    if (closest_An(x2, 5, k2a) != LQ_SUCCESS) return LQ_INVARG;

    xr2[0] = x1[0] - d1[0]; xr2[1] = x1[1] - d1[1];
    if (closest_An(xr2, 1, k1b) != LQ_SUCCESS) return LQ_INVARG;
    k1b[0] += d1[0]; k1b[1] += d1[1];

    for (i = 0; i < 6; i++) xr6[i] = x2[i] - d2[i];
    if (closest_An(xr6, 5, k2b) != LQ_SUCCESS) return LQ_INVARG;
    for (i = 0; i < 6; i++) k2b[i] += d2[i];

    c0[0] = k1a[0]; for (i = 0; i < 6; i++) c0[1 + i] = k2a[i]; c0[7] = k1a[1];
    c1[0] = k1b[0]; for (i = 0; i < 6; i++) c1[1 + i] = k2b[i]; c1[7] = k1b[1];

    memcpy(out, dist2(x, c0, 8) <= dist2(x, c1, 8) ? c0 : c1, 8 * sizeof(double));
    return LQ_SUCCESS;
}

/*! Closest-point algorithm for E6* lattice [3] */
static int closest_E6_star(const double* x, int n, double* out)
{
    static const double r[3][8] = {
        {0, 0, 0, 0, 0, 0, 0, 0},
        {0,-2.0 / 3.0,-2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0},
        {0, 2.0 / 3.0, 2.0 / 3.0,-1.0 / 3.0,-1.0 / 3.0,-1.0 / 3.0,-1.0 / 3.0, 0}
    };
    double yi[8], xr[8], tmp[8], best, d;
    int c, j;

    /* check arguments */
    if (!x || !out || n != 6) return LQ_INVARG;

	/* try all 3 cosets of E6 in LQ_An, find the closest point in each, and pick the best: */
    best = 1e30;
    for (c = 0; c < 3; c++) {
        for (j = 0; j < 8; j++) xr[j] = x[j] - r[c][j];
        if (closest_E6(xr, 6, tmp) != LQ_SUCCESS) return LQ_INVARG;
        for (j = 0; j < 8; j++) yi[j] = tmp[j] + r[c][j];
        d = dist2(x, yi, 8);
        if (d < best) { best = d; memcpy(out, yi, 8 * sizeof(double)); }
    }

    return LQ_SUCCESS;
}

/*! Closest-point algorithm for E7 lattice */
static int closest_E7(const double* x, int n, double* out)
{
    static const double r7[8] = { -0.5,-0.5,0.5,0.5,0.5,0.5,0.5,0.5 };
    double y0[8], y1[8], xr[8];
    int i;

    /* check arguments */
    if (!x || !out || n!=7) return LQ_INVARG;

	/* find closest point in LQ_An, and in LQ_An shifted by r7, and pick the best: */
    if (closest_An(x, 7, y0) != LQ_SUCCESS) return LQ_INVARG;
    for (i = 0; i < 8; i++) xr[i] = x[i] - r7[i];
    if (closest_An(xr, 7, y1) != LQ_SUCCESS) return LQ_INVARG;
    for (i = 0; i < 8; i++) y1[i] += r7[i];
    memcpy(out, dist2(x, y0, 8) <= dist2(x, y1, 8) ? y0 : y1, 8 * sizeof(double));

    return LQ_SUCCESS;
}

/*! Closest-point algorithm for E7* lattice */
static int closest_E7_star(const double* x, int n, double* out)
{
    static const double s[4][8] = {
        {-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5},
        { 0.5, 0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5},
        { 0.5, 0.5, 0.5, 0.5,-0.5,-0.5,-0.5,-0.5},
        { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,-0.5,-0.5}
    };
    double yi[8], xr[8], tmp[8], best, d;
    int c, j;
    
    /* check arguments */
    if (!x || !out || n!=7) return LQ_INVARG;

	/* try all 4 cosets of E7 in LQ_An, find the closest point in each, and pick the best: */
    best = 1e30;
    for (c = 0; c < 4; c++) {
        for (j = 0; j < 8; j++) xr[j] = x[j] - s[c][j];
        if (closest_An_star(xr, 7, tmp) != LQ_SUCCESS) return LQ_INVARG;
        for (j = 0; j < 8; j++) yi[j] = tmp[j] + s[c][j];
        d = dist2(x, yi, 8);
        if (d < best) { best = d; memcpy(out, yi, 8 * sizeof(double)); }
    }
    return LQ_SUCCESS;
}

/*! Closest-point algorithm for E8 lattice [1] */
static int closest_E8(const double* x, int n, double* out)
{
    double y0[8], y1[8], xh[8], fx[8], gx[8];
    int i;

    /* check arguments */
    if (!x || !out || n != 8) return LQ_INVARG;

    /* find closest point in A8, and in A8 shifted by (0.5,...,0.5), and pick the best: */
    fv(x, 8, fx); gv(x, 8, gx);
    memcpy(y0, (isum(fx, 8) & 1) == 0 ? fx : gx, 8 * sizeof(double));
    for (i = 0; i < 8; i++) xh[i] = x[i] - 0.5;
    fv(xh, 8, fx); gv(xh, 8, gx);
    if ((isum(fx, 8) & 1) == 0)
        for (i = 0; i < 8; i++) y1[i] = fx[i] + 0.5;
    else
        for (i = 0; i < 8; i++) y1[i] = gx[i] + 0.5;
    memcpy(out, dist2(x, y0, 8) <= dist2(x, y1, 8) ? y0 : y1, 8 * sizeof(double));

    return LQ_SUCCESS;
}

/*! Unified closest-point function type */
typedef int (*closest_fn_t)(const double* x, int n, double* out);


/*************
 *
 * Lattice point encoding and decoding algorithms based on methods from [2, Sec.III-IV].
 *
 */

 /*!
  *  \brief Encode a point (in ambient coords, already inside a + V_r) to an index.
  *
  *  Implements an algorithm from [2, Sec.III]:
  *          lp = closest(xc);
  *          c_j = lp * col_j(Ginv);
  *          index_j = c_j mod r;
  *
  *  \param[in]  xc     point in ambient coords, already clipped into a + V_r
  *  \param[in]  Ginv   inverse generator matrix, column-major (amb x n)
  *  \param[in]  n      lattice dimension n
  *  \param[in]  amb    ambient dimension
  *  \param[in]  cfn    closest-point function
  *  \param[in]  r      Voronoi region scale
  *  \param[out] idx    output index vector (length n), each in [0, r-1].
  *
  *  \return LQ_SUCCESS or error code.
  */
static int voronoi_encode(const double* xc, const double* Ginv, int n, int amb, closest_fn_t cfn, int r, int* idx)
{
    double lp[LQ_MAX_AMB];
    double cj;
    int i, j;

	/* check arguments: */
	if (!xc || !Ginv || !idx ||!cfn) return LQ_INVARG;

	/* retrieve closest point: */
    if (cfn(xc, n, lp) != LQ_SUCCESS) return LQ_INVARG;

	/* compute c_j = lp * col_j(Ginv), and reduce modulo r: */
    for (j = 0; j < n; j++) {
        cj = 0.0;
        for (i = 0; i < amb; i++)
            cj += lp[i] * Ginv[j * amb + i];
        idx[j] = mod_pos((int)floor(cj + 0.5), r);
    }

    return LQ_SUCCESS;
}

/*!
 *  \brief Decode index to a lattice point (ambient coords) given center a.
 *
 *  Implements an algorithm from [2, Sec.IV]:
 *     x' = sum(k_i*v_i);
 *     z = (x'-a)/r;
 *     lambda = closest(z);
 *     lp = x'-r*lambda;
 *
 *  \param[in]  idx    input index vector (length = n, values in [0,r-1])
 *  \param[in]  a      center offset vector (length = amb)
 *  \param[in]  G      generator matrix (n x amb)
 *  \param[in]  n      lattice dimension
 *  \param[in]  amb    ambient dimension
 *  \param[in]  cfn    closest point function
 *  \param[in]  r      lattice parameter r
 *  \param[out] lp     output lattice point (length = amb)
 *
 *  \return LQ_SUCCESS if success, LQ_INVARG if error.
 */
static int voronoi_decode(const int* idx, const double* a, const double* G, int n, int amb, closest_fn_t cfn, int r, double* lp)
{
    double xp[LQ_MAX_AMB], z[LQ_MAX_AMB], lambda[LQ_MAX_AMB];
    double inv_r = 1.0 / (double)r;
    int i, j;

    /* check arguments: */
	if (!idx || !a || !G || !lp || !cfn) return LQ_INVARG;

	/* compute x' = sum(k_i*v_i) */
    for (j = 0; j < amb; j++) {
        xp[j] = 0.0;
        for (i = 0; i < n; i++)
            xp[j] += (double)idx[i] * G[i * amb + j];
    }
    
	/* compute z = (x' - a) / r: */
    for (j = 0; j < amb; j++)
        z[j] = (xp[j] - a[j]) * inv_r;

	/* compute lambda = closest(z) */
    if (cfn(z, n, lambda) != LQ_SUCCESS) return LQ_INVARG;

	/* compute lp = x' - r * lambda: */
    for (j = 0; j < amb; j++)
        lp[j] = xp[j] - (double)r * lambda[j];

    return LQ_SUCCESS;
}


/*************
 * 
 *  Module's constants, static data, and tables.
 * 
 *  All these tables are filled once at the initialization.
 *  They should not cause any threading issues, since they are read-only after initialization, and the initialization is idempotent.
 */

/* Initialization state: */
static int g_initialized = 0;

/*! Ambient dimension table [lattice][n]; 0 = unsupported combination. */
static const int amb_table[LQ_TYPES][LQ_MAX_DIM+1] = {
    /* LQ_Zn */      {0, 1, 2, 3, 4, 5, 6, 7, 8},
    /* LQ_An */      {0, 2, 3, 4, 5, 6, 7, 8, 9},
    /* LQ_An_star */ {0, 2, 3, 4, 5, 6, 7, 8, 9},
    /* LQ_Dn */      {0, 0, 2, 3, 4, 5, 6, 7, 8},
    /* LQ_Dn_star */ {0, 0, 2, 3, 4, 5, 6, 7, 8},
    /* LQ_En */      {0, 0, 0, 0, 0, 0, 8, 8, 8},
    /* LQ_En_star */ {0, 0, 0, 0, 0, 0, 8, 8, 8}
};

/*! Closest-point functions dispatch table [lattice][n] */
static closest_fn_t fn_table[LQ_TYPES][LQ_MAX_DIM+1] = {
    /* LQ_Zn */      {NULL, closest_Zn, closest_Zn, closest_Zn, closest_Zn, closest_Zn, closest_Zn, closest_Zn, closest_Zn},
    /* LQ_An */      {NULL, closest_An, closest_An, closest_An, closest_An, closest_An, closest_An, closest_An, closest_An},
    /* LQ_An_star */ {NULL, closest_An_star, closest_An_star, closest_An_star, closest_An_star, closest_An_star, closest_An_star, closest_An_star, closest_An_star},
    /* LQ_Dn */      {NULL, NULL, closest_Dn, closest_Dn, closest_Dn, closest_Dn, closest_Dn, closest_Dn, closest_Dn},
    /* LQ_Dn_star */ {NULL, NULL, closest_Dn_star, closest_Dn_star, closest_Dn_star, closest_Dn_star, closest_Dn_star, closest_Dn_star, closest_Dn_star},
    /* LQ_En */      {NULL, NULL, NULL, NULL, NULL, NULL, closest_E6, closest_E7, closest_E8},
    /* LQ_En_star */ {NULL, NULL, NULL, NULL, NULL, NULL, closest_E6_star, closest_E7_star, closest_E8}
};

/*! Generator matrices for all latices: */
static double G_Dn_data[7][LQ_MAX_AMB * LQ_MAX_AMB];    /* n=2..8, [n-2] */
static double Ginv_Dn_data[7][LQ_MAX_AMB * LQ_MAX_AMB];
static double G_An_data[8][LQ_MAX_AMB * LQ_MAX_AMB];    /* n=1..8, [n-1] */
static double Ginv_An_data[8][LQ_MAX_AMB * LQ_MAX_AMB];
static double G_E8_data[64];
static double Ginv_E8_data[64];
static double G_E7_data[56];
static double Ginv_E7_data[56];
static double G_E6_data[48];
static double Ginv_E6_data[48];

/*! G matrix dispatch table [lattice][n]. NULL = identity or unsupported. */
static const double* G_table[LQ_TYPES][LQ_MAX_DIM+1] = {
    /* LQ_Zn */      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
    /* LQ_An */      {NULL, G_An_data[0], G_An_data[1], G_An_data[2], G_An_data[3], G_An_data[4], G_An_data[5], G_An_data[6], G_An_data[7]},
    /* LQ_An_star */ {NULL, G_An_data[0], G_An_data[1], G_An_data[2], G_An_data[3], G_An_data[4], G_An_data[5], G_An_data[6], G_An_data[7]},
    /* LQ_Dn */      {NULL, NULL, G_Dn_data[0], G_Dn_data[1], G_Dn_data[2], G_Dn_data[3], G_Dn_data[4], G_Dn_data[5], G_Dn_data[6]},
    /* LQ_Dn_star */ {NULL, NULL, G_Dn_data[0], G_Dn_data[1], G_Dn_data[2], G_Dn_data[3], G_Dn_data[4], G_Dn_data[5], G_Dn_data[6]},
    /* LQ_En */      {NULL, NULL, NULL, NULL, NULL, NULL, G_E6_data, G_E7_data, G_E8_data},
    /* LQ_En_star */ {NULL, NULL, NULL, NULL, NULL, NULL, G_E6_data, G_E7_data, G_E8_data}
};

/*! Ginv matrix dispatch table [lattice][n]. NULL = identity or unsupported. */
static const double* Ginv_table[LQ_TYPES][LQ_MAX_DIM+1] = {
    /* LQ_Zn */      {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
    /* LQ_An */      {NULL, Ginv_An_data[0], Ginv_An_data[1], Ginv_An_data[2], Ginv_An_data[3], Ginv_An_data[4], Ginv_An_data[5], Ginv_An_data[6], Ginv_An_data[7]},
    /* LQ_An_star */ {NULL, Ginv_An_data[0], Ginv_An_data[1], Ginv_An_data[2], Ginv_An_data[3], Ginv_An_data[4], Ginv_An_data[5], Ginv_An_data[6], Ginv_An_data[7]},
    /* LQ_Dn */      {NULL, NULL, Ginv_Dn_data[0], Ginv_Dn_data[1], Ginv_Dn_data[2], Ginv_Dn_data[3], Ginv_Dn_data[4], Ginv_Dn_data[5], Ginv_Dn_data[6]},
    /* LQ_Dn_star */ {NULL, NULL, Ginv_Dn_data[0], Ginv_Dn_data[1], Ginv_Dn_data[2], Ginv_Dn_data[3], Ginv_Dn_data[4], Ginv_Dn_data[5], Ginv_Dn_data[6]},
    /* LQ_En */      {NULL, NULL, NULL, NULL, NULL, NULL, Ginv_E6_data, Ginv_E7_data, Ginv_E8_data},
    /* LQ_En_star */ {NULL, NULL, NULL, NULL, NULL, NULL, Ginv_E6_data, Ginv_E7_data, Ginv_E8_data}
};

/* Pre-computed center offsets: [type][n][r], in ambient coordinates */
static double center_a[LQ_TYPES][LQ_MAX_DIM + 1][LQ_MAX_R + 1][LQ_MAX_AMB];

/***********
 * 
 * Matrix utilities 
 *
 */

#define PIVOT_EPS   1e-15

/*!
 *  \brief Invert an n x n matrix using Gaussian elimination with partial pivoting. 
 * 
 *  \param[in]  src     The source matrix (column-major).
 *  \param[in]  n       The dimension of the matrix.
 *  \param[out] dst     The destination matrix (column-major) to store the inverse.
 * 
 *  \return LQ_SUCCESS on success, LQ_INVARG on invalid argument.
 */
static int mat_inv(const double* src, int n, double* dst)
{
    double aug[LQ_MAX_DIM][LQ_MAX_DIM * 2];
    double pivot, fac, t;
    int i, j, k, p;

    /* check arguments: */
    if (!src || !dst || n <= 0 || n > LQ_MAX_DIM) return LQ_INVARG;

	/* create augmented matrix [src | I]: */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) aug[i][j] = src[i * n + j];
        for (j = 0; j < n; j++) aug[i][n + j] = (i == j) ? 1.0 : 0.0;
    }

	/* perform Gaussian elimination with partial pivoting: */
    for (k = 0; k < n; k++) {
        for (p = k, i = k + 1; i < n; i++)
            if (fabs(aug[i][k]) > fabs(aug[p][k])) p = i;
        if (fabs(aug[p][k]) < PIVOT_EPS) return LQ_INVARG;
        if (p != k)
            for (j = 0; j < 2 * n; j++) {
                t = aug[k][j]; aug[k][j] = aug[p][j]; aug[p][j] = t;
            }
        pivot = aug[k][k];
        for (j = 0; j < 2 * n; j++) aug[k][j] /= pivot;
        for (i = 0; i < n; i++) {
            if (i == k) continue;
            fac = aug[i][k];
            for (j = 0; j < 2 * n; j++) aug[i][j] -= fac * aug[k][j];
        }
    }

	/* extract the inverse from the augmented matrix: */
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            dst[j * n + i] = aug[i][n + j];

    return LQ_SUCCESS;
}

/*!
 *  \brief Compute the pseudo-inverse of a generator matrix G (n x amb) using the formula Ginv = (G^T G)^{-1} G^T.
 * 
 *  \param[in]  G        The generator matrix (n x amb, row-major).
 *  \param[in]  nrows    The number of rows in G (lattice dimension n).
 *  \param[in]  ncols    The number of columns in G (ambient dimension amb).
 *  \param[out] Ginv     The destination matrix (amb x n, column-major) to store the pseudo-inverse.
 * 
 *  \return LQ_SUCCESS on success, LQ_INVARG on invalid argument or if G^T G is not invertible.
 */
static int mat_pinv(const double* G, int nrows, int ncols, double* Ginv)
{
    double GGt[LQ_MAX_DIM*LQ_MAX_DIM], GGt_inv_cm[LQ_MAX_DIM*LQ_MAX_DIM], GGt_inv[LQ_MAX_DIM*LQ_MAX_DIM];
    int i, j, k;

	/* check arguments: */
	if (!G || !Ginv || nrows <= 0 || nrows > LQ_MAX_DIM || ncols <= 0 || ncols > LQ_MAX_DIM) return LQ_INVARG;

	/* Compute GGt = G * G^T: */
    memset(GGt, 0, sizeof(GGt));
    for (i = 0; i < nrows; i++)
        for (j = 0; j < nrows; j++)
            for (k = 0; k < ncols; k++)
                GGt[i * nrows + j] += G[i * ncols + k] * G[j * ncols + k];

	/* Invert GGt to get (G^T G)^{-1} in column-major order: */
    if (mat_inv(GGt, nrows, GGt_inv_cm) != LQ_SUCCESS) return LQ_INVARG;

	/* Transpose GGt_inv_cm to get (G^T G)^{-1} in row-major order: */
    for (i = 0; i < nrows; i++)
        for (j = 0; j < nrows; j++)
            GGt_inv[i * nrows + j] = GGt_inv_cm[j * nrows + i];

	/* Compute Ginv = (G^T G)^{-1} G^T: */
    memset(Ginv, 0, ncols * nrows * sizeof(double));
    for (i = 0; i < ncols; i++) 
        for (j = 0; j < nrows; j++) {
            double s = 0.0;
            for (k = 0; k < nrows; k++)
                s += G[k * ncols + i] * GGt_inv[k * nrows + j];
            Ginv[j * ncols + i] = s;
        }

    return LQ_SUCCESS;
}


/***********
 * 
 * Builders for generator matrices for all supported lattices.
 * 
 */

/*! Builds the generator matrix G for LQ_Dn lattice (n x n, row-major). */
static void build_G_Dn(int n, double* G)
{
    int i, sz = n * n;
    memset(G, 0, sz * sizeof(double));
    G[0] = 2.0;
    for (i = 1; i < n; i++) {
        G[i * n + (i - 1)] = 1.0;
        G[i * n + i] = 1.0;
    }
}

/*! Builds the generator matrix G for LQ_An lattice (n x (n+1), row-major). */
static void build_G_An(int n, double* G)
{
    int amb = n + 1, i;
    memset(G, 0, n * amb * sizeof(double));
    for (i = 0; i < n; i++) {
        G[i * amb + i] = 1.0;
        G[i * amb + i + 1] = -1.0;
    }
}

/*! Builds the generator matrix G for E6 lattice (6 x 8, row-major). */
static void build_G_E6(double* G)
{
    int i;
    memset(G, 0, 48 * sizeof(double));
    G[0 * 8 + 0] = 1.0; G[0 * 8 + 7] = -1.0;
    for (i = 1; i < 6; i++) {
        G[i * 8 + i] = 1.0;
        G[i * 8 + i + 1] = -1.0;
    }
}

/*! Builds the generator matrix G for E7 lattice (7 x 8, row-major). */
static void build_G_E7(double* G)
{
    int i;
    memset(G, 0, 56 * sizeof(double));
    for (i = 0; i < 7; i++) {
        G[i * 8 + i] = 1.0;
        G[i * 8 + i + 1] = -1.0;
    }
}

/*! Builds the generator matrix G for E8 lattice (8 x 8, row-major). */
static void build_G_E8(double* G)
{
    static const double data[64] = {
         2,    0,    0,    0,    0,    0,    0,    0,
         1,    1,    0,    0,    0,    0,    0,    0,
         1,    0,    1,    0,    0,    0,    0,    0,
         1,    0,    0,    1,    0,    0,    0,    0,
         1,    0,    0,    0,    1,    0,    0,    0,
         1,    0,    0,    0,    0,    1,    0,    0,
         1,    0,    0,    0,    0,    0,    1,    0,
        -0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5
    };
    memcpy(G, data, 64 * sizeof(double));
}


/************
 * 
 * Dimension embedding helper functions 
 *
 */

/*! Embeds a point from lattice coordinates to ambient coordinates. */
static void embed(int lattice, int n, const double* x, double* xa)
{
    int i; double s;

    switch (lattice) {
    case LQ_Zn: case LQ_Dn: case LQ_Dn_star:
        memcpy(xa, x, n * sizeof(double));
        break;
    case LQ_An: case LQ_An_star:
        s = 0.0;
        for (i = 0; i < n; i++) { xa[i] = x[i]; s += x[i]; }
        xa[n] = -s;
        break;
    case LQ_En: case LQ_En_star:
        if (n == 8) {
            memcpy(xa, x, 8 * sizeof(double));
        }
        else if (n == 7) {
            s = 0.0;
            for (i = 0; i < 7; i++) { xa[i] = x[i]; s += x[i]; }
            xa[7] = -s;
        }
        else {
            xa[0] = x[0];
            s = 0.0;
            for (i = 1; i < 6; i++) { xa[i] = x[i]; s += x[i]; }
            xa[6] = -s;
            xa[7] = -x[0];
        }
        break;
    }
}

/*! Extracts the lattice coordinates from ambient coordinates (inverse of embed). */
static void extract(int n, const double* ya, double* x)
{
    memcpy(x, ya, n * sizeof(double));
}

/************
 * 
 *  Voronoi-region specific helper functions: validation, center computation, and point counting.
 *
 */

#define MAX_CENTROID_SIZE INT_MAX
#define CENTROID_ITERS 5

/*! Validates the input parameters for lattice point encoding/decoding. */
static int validate(int n, int lattice, int r)
{
    if (r < 2 || r > LQ_MAX_R) return LQ_INVARG;
    if (n < 1 || n > 8) return LQ_INVARG;
    if (lattice < 0 || lattice >= LQ_TYPES) return LQ_NOTSUP;
	if ((lattice == LQ_Dn || lattice == LQ_Dn_star) && n < 2) return LQ_NOTSUP;
    if ((lattice == LQ_En || lattice == LQ_En_star) && n < 6) return LQ_NOTSUP;
    if (!fn_table[lattice][n]) return LQ_NOTSUP;
    return LQ_SUCCESS;
}

/*! Center offset computation [2, Sec.II] */
static int idx_inc(int* idx, int n, int r)
{
    int i = n - 1;
    while (i >= 0) {
        idx[i]++;
        if (idx[i] < r) return 0;
        idx[i] = 0;
        i--;
    }
    return 1;
}

/*! 
 *  \brief Computes the number of points in the Voronoi region a + V_r, which is r^n. 
 * 
 *  \param[in]  r   Voronoi region scale.
 *  \param[in]  n   Lattice dimension.
 * 
 *  \return 0 if overflow. 
 */
static int n_points(int r, int n)
{
    int result = 1, i;
    for (i = 0; i < n; i++) {
        if (result > MAX_CENTROID_SIZE / r) 
            return 0;
        result *= r;
    }
    return result;
}

/*!
 *  \brief Compute center a for (lattice, n, r).
 * 
 *   This function implements the iterative centroid computation described in [2, Sec.II].
 */
static int compute_center(int lattice, int n, int r)
{
    const double* Gmat;
    int amb, M, iter, i, j;
    double a[LQ_MAX_AMB], centroid[LQ_MAX_AMB], lp[LQ_MAX_AMB];
    int idx[LQ_MAX_DIM];
    closest_fn_t cfn;

    /* validate parameters */
	if (validate(n, lattice, r) != LQ_SUCCESS) return LQ_INVARG;

	/* compute the number of points in the Voronoi region */
    M = n_points(r, n);
    if (M <= 0 || M > MAX_CENTROID_SIZE) return LQ_INVARG;

    /* LQ_Zn: center is zero */
    if (lattice == LQ_Zn) {
        memset(center_a[lattice][n][r], 0, LQ_MAX_AMB * sizeof(double));
        return LQ_SUCCESS;
    }

    /* get the closest-point function, generator matrix, and ambient dimension for this lattice. */
    cfn = fn_table[lattice][n];
    if (!cfn) return LQ_NOTSUP;
    Gmat = G_table[lattice][n];
    if (!Gmat) return LQ_NOTSUP;
    amb = amb_table[lattice][n];
    if (amb == 0) return LQ_NOTSUP;

    /* pick initial a: small fractional offset avoiding boundary */
    for (j = 0; j < amb; j++)
        a[j] = (double)(j + 1) / (double)(n * r + 1);

    /* enforce embedding constraints */
    if (lattice == LQ_An || lattice == LQ_An_star) {
        double s = 0.0;
        for (j = 0; j < amb; j++) s += a[j];
        s /= (double)amb;
        for (j = 0; j < amb; j++) a[j] -= s;
    }
    else 
    if (lattice == LQ_En || lattice == LQ_En_star) {
        if (n == 7) {
            double s = 0.0;
            for (j = 0; j < 7; j++) s += a[j];
            a[7] = -s;
        }
        else if (n == 6) {
            double s = 0.0;
            a[7] = -a[0];
            for (j = 1; j < 6; j++) s += a[j];
            a[6] = -s;
        }
    }

    /* iterative centroid refinement */
    for (iter = 0; iter < CENTROID_ITERS; iter++) {
        memset(centroid, 0, sizeof(centroid));
        memset(idx, 0, n * sizeof(int));

        for (i = 0; i < M; i++) {
            if (voronoi_decode(idx, a, Gmat, n, amb, cfn, r, lp) != LQ_SUCCESS)
                return LQ_INVARG;
            for (j = 0; j < amb; j++)
                centroid[j] += lp[j];
            idx_inc(idx, n, r);
        }

        for (j = 0; j < amb; j++)
            a[j] = centroid[j] / (double)M;
    }

    memcpy(center_a[lattice][n][r], a, LQ_MAX_AMB * sizeof(double));
    return LQ_SUCCESS;
}

/*! threshold for inside region check */
#define INSIDE_THRESHOLD 1e-12

/*!
 *  \brief Clip xa into the Voronoi region a + V_r, returning the clipped point and overload flag.
 * 
 *  This function implements the following clipping algorithm:
 *      z = (xa - a) / r;
 *      lambda = closest(z);            // 0 if xa is inside the Voronoi region, non-zero if outside
 *      out = xa - r * lambda;          // if xa is inside, this returns xa; if outside, this returns the closest point on the boundary of a + V_r.
 * 
 *  \param[in]  lattice  lattice type
 *  \param[in]  n        lattice dimension 
 *  \param[in]  xa       input vector (augmented to ambient dimension)
 *  \param[in]  amb      ambient dimension
 *  \param[in]  r        Voronoi region scale (>= 2).
 *  \param[in]  a        center of the Voronoi region
 *  \param[out] out      clipped output vector
 *  \param[out] overload 1 if xa is outside the Voronoi region (overload), 0 if inside.
 *
 *  \return LQ_SUCCESS, LQ_INVARG, or LQ_NOTSUP.
 */
static int voronoi_clip(int lattice, int n, const double* xa, int amb, int r, const double* a, double* out, int *overload)
{
    double z[LQ_MAX_AMB], lam[LQ_MAX_AMB];
    double inv_r = 1.0 / (double)r;
    closest_fn_t cfn;
    int i, rc;

    /* validate parameters */
    if (!xa || !a || !overload) return LQ_INVARG;
    if ((rc = validate(n, lattice, r)) != LQ_SUCCESS) return rc;

    /* get the closest-point function  */
    cfn = fn_table[lattice][n];
    if (!cfn) return LQ_NOTSUP;

	/* compute z = (xa - a) / r: */
    for (i = 0; i < amb; i++)
        z[i] = (xa[i] - a[i]) * inv_r;

	/* compute lambda = closest(z) */
    if (cfn(z, n, lam) != LQ_SUCCESS) return LQ_INVARG;

    /* check if lambda is zero (inside) or not (overload): */
    for (i = 0; i < amb; i++)
        if (fabs(lam[i]) > INSIDE_THRESHOLD) {
            /* lambda is non-zero, xa is outside the Voronoi region */
            *overload = 1;
			/* clip xa to the closest point on the boundary of a + V_r: */
            for (i = 0; i < amb; i++)
                out[i] = xa[i] - (double)r * lam[i];
            return LQ_SUCCESS;
        }

    /* no overload detected, copy input to output */
    *overload = 0;
	memcpy(out, xa, amb * sizeof(double));

    return LQ_SUCCESS;
}

/******************
 *
 *  Public API functions:
 *
 */

/*! 
 *  \brief Initialize library: build lattice generator matrices, pre-compute center offsets.
 */
int lq_init(void)
{
    int n, lattice, r, M;

    /* build G and Ginv matrices: */
    for (n = 2; n <= 8; n++) {
        build_G_Dn(n, G_Dn_data[n - 2]);
        mat_inv(G_Dn_data[n - 2], n, Ginv_Dn_data[n - 2]);
    }
    for (n = 1; n <= 8; n++) {
        build_G_An(n, G_An_data[n - 1]);
        mat_pinv(G_An_data[n - 1], n, n + 1, Ginv_An_data[n - 1]);
    }
    build_G_E8(G_E8_data);
    mat_inv(G_E8_data, 8, Ginv_E8_data);
    build_G_E7(G_E7_data);
    mat_pinv(G_E7_data, 7, 8, Ginv_E7_data);
    build_G_E6(G_E6_data);
    mat_pinv(G_E6_data, 6, 8, Ginv_E6_data);

    /* compute center offsets */
    memset(center_a, 0, sizeof(center_a));
    for (lattice = 0; lattice < LQ_TYPES; lattice++) {
        for (n = 1; n <= LQ_MAX_DIM; n++) {
            if (!fn_table[lattice][n]) 
                continue;
            for (r = 2; r <= LQ_MAX_R; r++) {
                M = n_points(r, n);
                if (M <= 0 || M > MAX_CENTROID_SIZE)
                    return LQ_INVARG;
                compute_center(lattice, n, r);
            }
        }
    }

	/* mark initialized & exit: */
    g_initialized = 1;
    return LQ_SUCCESS;
}

/*!
 *  \brief Quantize input vector x to a lattice index.
 *
 *  \param[in]  x        Input vector (length n).
 *  \param[in]  n        Dimension (1..8). For LQ_En/En_star only n=6,7,8 supported.
 *  \param[in]  lattice  Lattice type.
 *  \param[in]  r        Voronoi region scale (>= 2).
 *  \param[out] index    Output index vector (length n).
 *  \param[out] overload 0 if x is inside V_r, 1 if outside.
 *
 *  \return LQ_SUCCESS, LQ_INVARG, LQ_NOTSUP, or LQ_NOTINIT.
 */
int lq_quantize(const double* x, int n, int lattice, int r, int* index, int* overload)
{
	double xa[LQ_MAX_AMB];      /* x embedded in ambient space */
    double xc[LQ_MAX_AMB];      /* x clipped to a + V_r */  
    const double *Ginv, *a;
    closest_fn_t cfn;
    int amb, j, rc;

    /* check initialization */
    if (!g_initialized) return LQ_NOTINIT;

    /* check parameters */
    if (!x || !index || !overload) return LQ_INVARG;
    if ((rc = validate(n, lattice, r)) != LQ_SUCCESS) return rc;

    /* retrieve ambient, closest-point function & center point: */
    amb = amb_table[lattice][n];
    cfn = fn_table[lattice][n];
    a = center_a[lattice][n][r];
    if (!amb || !cfn) return LQ_NOTSUP;

    /* embed user n-dim input into ambient space */
    embed(lattice, n, x, xa);

    /* clip into a + V_r */
    if (voronoi_clip(lattice, n, xa, amb, r, a, xc, overload) != LQ_SUCCESS)
        return LQ_INVARG;

    /* LQ_Zn special case: G=I, Ginv=I */
    if (lattice == LQ_Zn) {
        double lp[LQ_MAX_AMB];
        if (cfn(xc, n, lp) != LQ_SUCCESS) return LQ_INVARG;
        for (j = 0; j < n; j++)
            index[j] = mod_pos((int)(lp[j] + (lp[j] >= 0 ? 0.5 : -0.5)), r);
        return LQ_SUCCESS;
    }

    /* retrieve inverse generator matrix */
    Ginv = Ginv_table[lattice][n];
    if (!Ginv) return LQ_NOTSUP;

    /* Encode clipped point to index */
    return voronoi_encode(xc, Ginv, n, amb, cfn, r, index);
}

/*!
 *  \brief Reconstruct vector from lattice index.
 *
 *  \param[out] x        Reconstructed vector (length n).
 *  \param[in]  n        Dimension (same as used in lq_quantize).
 *  \param[in]  lattice  Lattice type.
 *  \param[in]  r        Voronoi region scale.
 *  \param[in]  index    Index vector (length n).
 *
 *  \return LQ_SUCCESS, LQ_INVARG, LQ_NOTSUP, or LQ_NOTINIT.
 */
int lq_reconstruct(double* x, int n, int lattice, int r, const int* index)
{
    double lp[LQ_MAX_AMB];
    closest_fn_t cfn;
    const double* G, * a;
    int amb, i, j, rc;

    /* check initialization */
    if (!g_initialized) return LQ_NOTINIT;

    /* check parameters */
    if (!x || !index) return LQ_INVARG;
    if ((rc = validate(n, lattice, r)) != LQ_SUCCESS) return rc;

    /* retrieve ambient, closest-point function, & center point: */
    amb = amb_table[lattice][n];
    cfn = fn_table[lattice][n];
    a = center_a[lattice][n][r];
	if (!amb || !cfn) return LQ_NOTSUP;

    /* validate index entries */
    for (i = 0; i < n; i++)
        if (index[i] < 0 || index[i] >= r) return LQ_INVARG;

    /* LQ_Zn special case */
    if (lattice == LQ_Zn) {
        for (j = 0; j < n; j++)
            x[j] = (index[j] >= r / 2) ? (double)(index[j] - r): (double)index[j];
        return LQ_SUCCESS;
    }

	/* retrieve generator matrix */
    G = G_table[lattice][n];
    if (!G) return LQ_NOTSUP;

    /* index -> lattice point in ambient coords */
    if (voronoi_decode(index, a, G, n, amb, cfn, r, lp) != LQ_SUCCESS)
        return LQ_INVARG;

    /* extract user n-dim output */
    extract(n, lp, x);
    return LQ_SUCCESS;
}

/* lq.c -- end of file */
