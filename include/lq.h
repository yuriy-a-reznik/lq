/*!
 *  \file lq.h
 * 
 *  \brief LQ: a very basic lattice quantization library.
 *
 *  This library implements lattice quantization algorithms for lattices Zn, An, Dn, En and their duals.
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
 *  Copyright (c) 2026 Yuriy A. Reznik
 *  Licensed under the MIT License: https://opensource.org/licenses/MIT
 *
 *  \author  Yuriy A. Reznik
 *  \version 1.02
 *  \date    April 27, 2026
 */

#ifndef LQ_H
#define LQ_H

#ifdef __cplusplus
extern "C" {
#endif

/*! Constants / limits: */
#define LQ_MAX_DIM 8		/*!< maximum supported dimension n */
#define LQ_MAX_AMB 9		/*!< maximum internal (ambient) working dimension */
#define LQ_MAX_R   5		/*!< maximum Voronoi region scale */

/*! Lattice types: */
enum {
	LQ_Zn = 0,
	LQ_An = 1,
	LQ_An_star = 2,
	LQ_Dn = 3,
	LQ_Dn_star = 4,
	LQ_En = 5,
	LQ_En_star = 6,
	LQ_TYPES = 7
};

/*! Return codes: */
enum {
	LQ_SUCCESS	= 0,
	LQ_INVARG	= 1,
	LQ_NOTSUP	= 3,
	LQ_NOTINIT	= 4
};

/*! function prototypes: */

/*!
 *  \brief Initialize the lattice quantizer. Must be called once before use.
 *  \return LQ_SUCCESS.
 */
int lq_init(void);

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
int lq_quantize(const double* x, int n, int lattice, int r, int* index, int* overload);

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
int lq_reconstruct(double* x, int n, int lattice, int r, const int* index);

#ifdef __cplusplus
}
#endif

#endif /* LQ_H */