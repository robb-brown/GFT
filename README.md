# GFT
The General Fourier Family Transform (GFT) is a time-frequency transform; the Fourier, short-time Fourier and S- transforms are special cases. This is an efficient algorithm for computing the GFT, which is most interesting to use to compute the S-transform (ST). The ST is like a short time Fourier transform in that it produces a frequency vs. time (or space) spectrogram, but it is frequency adaptive like a wavelet.

This algorithm has a complexity of O(n log n) (in 1-dimension) compared to the original "fast S-transform" algorithm and subsequent optimzations such as the Discrete Orthogonal S-Transform (DOST), which are all O(N^2).

This code should be functional under Python 3.x, but is very much a work in progress. Run gft.py with 

	python -i gft.py

to see a demo of the 1-dimensional forward and inverse fast ST.

The GFT and this efficient algorithm for computing it was published in 2009 in IEEE Transactions on Signal Processing. The full text is available [here](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5184926&casa_token=-wdb__eqE3EAAAAA:XheUu232GAVUPrMsvwAdFBZH_2wyUkcpV9aPtt4G10Ay-CaH3D-Hk07XVW7xttm4XjsRFRcK_w).

A previous Python/C version of the GFT was hosted at [sourceforge.net/projects/fst-uofc/](https://sourceforge.net/projects/fst-uofc/). That project is no longer maintained.

The O(n log n) algorithm is the subject of US [patent US8458240B2](https://patents.google.com/patent/US8458240B2/en).
