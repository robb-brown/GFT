# GFT
The General Fourier Family Transform (GFT) is a time-frequency transform; the Fourier, short-time Fourier, S- and many wavelet transforms are special cases. This is an efficient algorithm for computing the GFT, which is most interesting to use to compute the S-transform (ST). The ST is like a short time Fourier transform in that it produces a frequency vs. time (or space) spectrogram, including phase information, but it is frequency adaptive like a wavelet.

This algorithm has a complexity of O(n log n) (in 1-dimension) compared to the original "fast S-transform" algorithm and subsequent optimizations such as the Discrete Orthogonal S-Transform (DOST), which are all O(N^2). Equally importantly, the default GFT transform is one-to-one: the GFT (by default) produces an N-point spectrum from an N-point signal. This is in contrast to e.g. the fast ST, which produces an N^2 point spectrum. The N point spectrum not only makes large transforms practical to store but also means that the transform is uniquely invertible. The spectrum can be modified and inverse transformed without ambiguity. If desired, the GFT formalism supports oversampling although this implementation does not.

This module implements the basic forward and inverse transform; simple dyadic window
construction; boxcar and Gaussian windows; and interpolation to create a displayable
spectrum. Customized window layouts and windows are fairly easily constructed by
following the example of the existing dyadic/Gaussian ones.

TODO:

- demo/tutorial of creating custom windows, especially for musical scales
- N dimensional transforms
- real transforms


gft.m is a self contained MatLab/Octave script file that implements the forward and inverse GFT, and includes a demo. It has only been tested in Octave.


To see a demo, import the gft package and call the runDemo() function:

	import gft
	gft.runDemo()

This requires matplotlib, numpy and scipy, which you may have to install separately:

	pip install matplotlib, numpy, scipy


Here is an example of a basic transform with a test signal and using the GFT object's numpy-style slicing. The slicing is mostly useful for displaying spectrograms because it will geenrally interpolate the varying-shaped GFT samples onto a regular grid. The second last example shows how you can retrieve the actual samples; note that because each frequency band has a different time-sampling rate you have to do this one frequency band at a time.

	import gft
	
	# define the signal length and create a test signal
	N = 1024
	sig = gft.signalDelta(N)
	
	# create dyadic partitions and boxcar windows
	partitions = gft.GFTPartitions.dyadicPartitions(N)
	windows = gft.GFTPartitions.boxcarWindows(N,partitions)
	
	# create a GFT object, which computes the forward gft
	# sampleRate is optional. If not provided gft object indexing (see below)
	# will be by sample # on both f and t axes. Otherwise it is by 
	# time unit and reciprocal time unit (frequency)
	# interpolation can be 'nearest' or 'linear'. Beware of interpolated phase.
	SIG = gft.GFT(signal=sig,partitions=partitions,windows=windows,sampleRate=None,interpolation='nearest')
	
	# The GFT object supports numpy-style slicing, with the first coordinate f and the
	# second coordinate t. f ranges from -Nyquist to +Nyquist, 0 is DC.
	# If not supplied, f and/or t will default to the entire spectrum represented by N 
	# interpolated samples.
	# Slicing is mostly useful for display because in most cases it will interpolate the
	# spectrum.
	
	# get the entire spectrum interpolated to NxN (don't do this if your N is too large)
	spectrum = SIG[:]
	
	# get the positive frequency half only
	positiveSpectrum = SIG[0:]
	
	# get a spectrum with a specified sampling (1000x1000; useful for large signals)
	spectrum = SIG[::N/1000.,::N/1000.]
	
	# get a subsection of the spectrum
	spectrum = SIG[0:N/2,N/4:3*N/4]
	
	# get actual (not interpolated) samples
	# first get the actual frequency sample coordinates
	fCoords = SIG.fCoordinates()
	
	# The t sampling depends on what frequency band we're looking at
	# let's use the second higheset frequency band for an N=1024 point signal
	tCoords = SIG.tCoordinates(f=fCoords[8])
	
	# get the samples
	# we transformed a unit delta so there should be one sample that is 1, all others 0
	spectrum = SIG[fCoords[8],tCoords]
	
	# get the interpolated samples
	# now we should have 8 samples equal to 1
	spectrum = SIG[fCoords[8],:]
	
	



The GFT and this efficient algorithm for computing it was published in 2009 in IEEE Transactions on Signal Processing. The full text is available [here](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5184926&casa_token=-wdb__eqE3EAAAAA:XheUu232GAVUPrMsvwAdFBZH_2wyUkcpV9aPtt4G10Ay-CaH3D-Hk07XVW7xttm4XjsRFRcK_w). The GFT can also be computed from the [time domain](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4649729&casa_token=1K0y20cH5_cAAAAA:vT9IxeMjzAPF1-rgM_28gYaRTgKgniPoioVUvZd3zr02TF5kwAQtrsY4S5-8W0j25H5CsBgy&tag=1). Finally, a [book chaper](https://www.intechopen.com/books/recent-advances-in-biomedical-engineering/developments-in-time-frequency-analysis-of-biomedical-signals-and-images-using-a-generalized-fourier) providing an overview of time frequency transforms and the GFT, including the link to the wavelet transform, is available.

A previous Python/C version of the GFT was hosted at [sourceforge.net/projects/fst-uofc/](https://sourceforge.net/projects/fst-uofc/). That project is no longer maintained.

The O(n log n) algorithm is the subject of US [patent US8458240B2](https://patents.google.com/patent/US8458240B2/en).
