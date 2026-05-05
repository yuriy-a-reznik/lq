/*!
 *  \file lq_test.c
 *
 *  \brief Test program for LQ (lattice quantization) library.
 *
 *  This program implements few basic tests for lattice quantization and reconstruction functions.
 *  It tests supported lattices (Z_n, D_n, D_n*, A_n, A_n*, E_8, E_7, E_6) with various input vectors and r values.
 *  The results are printed to the console, including quantized indices, reconstructed vectors, and round-trip verification.
 *
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.02
 *  \date    April 27, 2026
 */

#include <stdio.h>
#include <math.h>
#include "lq.h"

/*! Print floating-point vector with label. */
static void print_vector(const char *s, const double *v, int n)
{
    int i;
    printf("  %s = (", s);
    for (i = 0; i < n; i++) printf("%7.4f%s", v[i], i < n-1 ? "," : "");
    printf(")\n");
}

/*! Print integer vector with label. */
static void print_int_vector(const char *s, const int *v, int n)
{
    int i;
    printf("  %s = (", s);
    for (i = 0; i < n; i++) printf("%d%s", v[i], i < n-1 ? "," : "");
    printf(")\n");
}

/*! Compute squared error between two vectors. */
static double sq_err(const double *a, const double *b, int n)
{
    double s = 0; int i;
    for (i = 0; i < n; i++) { double d = a[i]-b[i]; s += d*d; }
    return s;
}

/*!
 *  \brief Quantize then reconstruct, verify round-trip. 
 *
 *  \param name    test name (for printing)
 *  \param x       input vector
 *  \param n       dimension of the vector
 *  \param lattice lattice type
 *  \param r       granular region scale parameter r 
 */
static void test(const char *name, const double *x, int n, int lattice, int r)
{
    double y[8], x2[8];
    int idx[8], idx2[8], ovl, ovl2, rc, match, i;

	/* Print input */
    printf("%s, n=%d, r=%d:\n", name, n, r);
    print_vector("x", x, n);

	/* Quantize */
    rc = lq_quantize(x, n, lattice, r, idx, &ovl);
    if (rc != LQ_SUCCESS) { printf("  quantize error: %d\n\n", rc); return; }
    printf("  overload=%d\n", ovl);
    print_int_vector("idx", idx, n);

	/* Reconstruct */
    rc = lq_reconstruct(y, n, lattice, r, idx);
    if (rc != LQ_SUCCESS) { printf("  reconstruct error: %d\n\n", rc); return; }
    print_vector("y", y, n);
    printf("  ||x-y||^2 = %.6f\n", sq_err(x, y, n));

    /* Verify: re-quantize reconstructed point, should get same index */
    rc = lq_quantize(y, n, lattice, r, idx2, &ovl2);
    if (rc != LQ_SUCCESS) { printf("  re-quantize error: %d\n\n", rc); return; }
    for (match = 1, i = 0; i < n; i++)
        if (idx[i] != idx2[i]) { match = 0; break; }
    printf("  round-trip: %s\n\n", match ? "PASS" : "FAIL");
}

/*!
 *  \brief Main test program. 
 */
int main(void)
{
    double x[8];
    int idx[8], ovl, rc;

    printf("Lattice Quantization Library Tests:\n\n");

    /* Test error before init */
    x[0] = 0;
    rc = lq_quantize(x, 4, LQ_Dn, 4, idx, &ovl);
    printf("Before init: rc=%d (expect %d)\n\n", rc, LQ_NOTINIT);

	/* Initialize library */
    lq_init();
    printf("Initialized.\n\n");

    /* Test unsupported lattice type: E5 */
    rc = lq_quantize(x, 5, LQ_En, 4, idx, &ovl);
    printf("E5: rc=%d (expect %d)\n\n", rc, LQ_NOTSUP);

    /* Z_4 */
    x[0]=0.6; x[1]=-1.1; x[2]=1.7; x[3]=0.1;
    test("Z4", x, 4, LQ_Zn, 4);

    /* D_4 */
    x[0]=0.6; x[1]=-1.1; x[2]=1.7; x[3]=0.1;
    test("D4", x, 4, LQ_Dn, 4);

    /* D_4* */
    x[0]=0.2; x[1]=0.5; x[2]=0.8; x[3]=0.1;
    test("D4*", x, 4, LQ_Dn_star, 4);

    /* A_2 (user gives 2 components) */
    x[0]=0.3; x[1]=-0.2;
    test("A2", x, 2, LQ_An, 4);

    /* A_3 */
    x[0]=0.4; x[1]=-0.3; x[2]=0.1;
    test("A3", x, 3, LQ_An, 4);

    /* A_2* */
    x[0]=0.3; x[1]=-0.2;
    test("A2*", x, 2, LQ_An_star, 4);

    /* E_8 */
    x[0]=0.1; x[1]=0.1; x[2]=0.8; x[3]=1.3;
    x[4]=2.2; x[5]=-0.6; x[6]=-0.7; x[7]=0.9;
    test("E8", x, 8, LQ_En, 4);

    /* E_8* (self-dual) */
    test("E8*", x, 8, LQ_En_star, 4);

    /* E_7 */
    x[0]=0.1; x[1]=-0.2; x[2]=0.3; x[3]=-0.1;
    x[4]=0.2; x[5]=-0.3; x[6]=0.1;
    test("E7", x, 7, LQ_En, 4);

    /* E_7* */
    test("E7*", x, 7, LQ_En_star, 4);

    /* E_6 */
    x[0]=0.1; x[1]=-0.2; x[2]=0.3; x[3]=-0.1;
    x[4]=0.2; x[5]=-0.3;
    test("E6", x, 6, LQ_En, 4);

    /* E_6* */
    test("E6*", x, 6, LQ_En_star, 4);

    /* Overload test */
    printf("Overload detection tests:\n");
    x[0]=10.0; x[1]=10.0; x[2]=10.0; x[3]=10.0;
    lq_quantize(x, 4, LQ_Dn, 4, idx, &ovl);
    printf("  D4 far: overload=%d (expect 1)\n", ovl);
    x[0]=0.1; x[1]=0.1; x[2]=-0.1; x[3]=-0.1;
    lq_quantize(x, 4, LQ_Dn, 4, idx, &ovl);
    printf("  D4 near: overload=%d (expect 0)\n\n", ovl);

    /* Overload clipping: quantize far point, reconstruct, verify valid */
    printf("Overload clipping tests:\n");
    x[0]=10.0; x[1]=10.0; x[2]=10.0; x[3]=10.0;
    {
        double y[4];
        int idx2[4], ovl2;
        lq_quantize(x, 4, LQ_Dn, 4, idx, &ovl);
        printf("  original overload=%d\n", ovl);
        print_int_vector("idx", idx, 4);
        lq_reconstruct(y, 4, LQ_Dn, 4, idx);
        print_vector("recon", y, 4);
        lq_quantize(y, 4, LQ_Dn, 4, idx2, &ovl2);
        printf("  recon overload=%d (expect 0)\n", ovl2);
        {
            int match = 1, i;
            for (i = 0; i < 4; i++)
                if (idx[i] != idx2[i]) { match = 0; break; }
            printf("  clip round-trip: %s\n\n", match ? "PASS" : "FAIL");
        }
    }

    /* Try different r values */
    printf("D4 with r=2,3,4:\n");
    x[0]=0.6; x[1]=-1.1; x[2]=1.7; x[3]=0.1;
    test("D4 r=2", x, 4, LQ_Dn, 2);
    test("D4 r=3", x, 4, LQ_Dn, 3);
    test("D4 r=4", x, 4, LQ_Dn, 4);

    printf("All tests done.\n");
    return 0;
}

/* lq_test.c - end of file */
