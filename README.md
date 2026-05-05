# lq
LQ: a very basic lattice quantization library.

This library implements lattice quantization algorithms for lattices Zn, An, Dn, En and their duals.
Dimensions up to 8 are supported.
 
This includes best known lattices for minimizing covering radius:

     1	Z	        Integer lattice
     2	A2* ~ A2	Hexagonal lattice
     3	A3* ~ D3*	BCC / FCC dual
     4	A4*	        A4 dual
     5	A5*	        A5 dual
     6	A6*	        A6 dual
     7	A7*	        A7 dual
     8	A8*	        A8 dual

and best lattices for minimizing quantization error:

     1	Z	        Integer lattice
     2	A2	        Hexagonal lattice
     3	A3* ~ D3*	BCC / FCC dual
     4	D4	        4D checkerboard lattice, Hurwitz integers
     5	D5*	        D5 dual
     6	E6*	        E6 dual
     7	E7*	        E7 dual
     8	E8	        Gossett lattice
 
In all cases, this library uses the Voronoi region scaled by a factor r >= 2 as the granular region.
The resulting number of points is r^n, and the effective rate of the quantizer is n·log2(r) bits.
The distortion is upper-bounded by the covering radius for points within the granular region.
Points outside the granular region are clipped to the nearest point within it, with the overload flag set to 1.
The overload region is unbounded.

Main references:
1. J.H. Conway and N.J.A. Sloane, "Fast quantizing and decoding algorithms for lattice quantizers and codes," IEEE Trans. Inform. Theory, vol. IT-28, no. 2, pp. 227-232, March 1982.
2. J.H. Conway and N.J.A. Sloane, "A fast encoding method for lattice codes and quantizers," IEEE Trans. Inform. Theory, vol. IT-29, no. 6, pp. 820-824, November 1983.
3. K. Takizawa, H. Yagi and T. Kawabata, "Closest point algorithms with lp norm for root lattices," 2010 IEEE International Symposium on Information Theory, Austin, TX, USA, 2010, pp. 1042-1046.
4. J.H. Conway and N.J.A. Sloane, "Sphere packings, lattices and groups," Vol. 290, Springer, 2013.
