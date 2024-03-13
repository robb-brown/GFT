# GFT
The General Fourier Family Transform (GFT) is a time-frequency transform; the Fourier, short-time Fourier, S- and many wavelet transforms are special cases. This is an efficient algorithm for computing the GFT, which is most interesting to use to compute the S-transform (ST). The ST is like a short time Fourier transform in that it produces a frequency vs. time (or space) spectrogram, including phase information, but it is frequency adaptive like a wavelet.

This algorithm has a complexity of O(n log n) (in 1-dimension) compared to the original "fast S-transform" algorithm and subsequent optimizations such as the Discrete Orthogonal S-Transform (DOST), which are all O(N^2).

This module implements the basic forward and inverse transform; simple dyadic window
construction; boxcar and Gaussian windows; and interpolation to create a displayable
spectrum. Customized window layouts and windows are fairly easily constructed by
following the example of the existing dyadic/Gaussian ones.

TODO:

- demo/tutorial of creating custom windows
- N dimensional transforms
- real transforms


gft.m is a self contained MatLab/Octave script file that implements the forward and inverse GFT, and includes a demo. It has only been tested in Octave.


To see a demo, import the gft package and call the runDemo() function:

	import gft
	gft.runDemo()

This requires matplotlib, which you may have to install separately:

	pip install matplotlib



The GFT and this efficient algorithm for computing it was published in 2009 in IEEE Transactions on Signal Processing. The full text is available [here](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5184926&casa_token=-wdb__eqE3EAAAAA:XheUu232GAVUPrMsvwAdFBZH_2wyUkcpV9aPtt4G10Ay-CaH3D-Hk07XVW7xttm4XjsRFRcK_w). The GFT can also be computed from the [time domain](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4649729&casa_token=1K0y20cH5_cAAAAA:vT9IxeMjzAPF1-rgM_28gYaRTgKgniPoioVUvZd3zr02TF5kwAQtrsY4S5-8W0j25H5CsBgy&tag=1). Finally, a [book chaper](https://www.intechopen.com/books/recent-advances-in-biomedical-engineering/developments-in-time-frequency-analysis-of-biomedical-signals-and-images-using-a-generalized-fourier) providing an overview of time frequency transforms and the GFT, including the link to the wavelet transform, is available.

A previous Python/C version of the GFT was hosted at [sourceforge.net/projects/fst-uofc/](https://sourceforge.net/projects/fst-uofc/). That project is no longer maintained.

The O(n log n) algorithm is the subject of US [patent US8458240B2](https://patents.google.com/patent/US8458240B2/en).
