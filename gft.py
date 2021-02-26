import numpy, scipy.interpolate

class GFTPartitions(object):
	""" Superclass for GFT window sets. Implements generic boxcar octave windows."""
	
	@classmethod 
	def widths(self,partitions):
		"""Takes a partitions array and returns the widths of the partitions"""
		widths = (numpy.roll(partitions,-1) - partitions); widths[-1] += max(partitions)+1
		return widths
	
	@classmethod
	def windows(self,N):
		# octave partitions
		widths = (2**numpy.arange(numpy.round(numpy.log(N)/numpy.log(2))-1)).astype(numpy.uint32)
		widths = numpy.concatenate([[1],widths,widths[::-1]])
		partitions = numpy.concatenate([[0],numpy.cumsum(widths)])
		# Boxcar windows are flat
		windows = numpy.ones(N)
		return windows,partitions
	
	

class GFT(object):

	@classmethod
	def transform(self,signal,windows,partitions):
		SIG = numpy.fft.fft(signal)
		SIG *= windows
		for p in range(len(partitions)-1):
			SIG[partitions[p]:partitions[p+1]] = numpy.fft.ifft(SIG[partitions[p]:partitions[p+1]])
		return SIG
		
	@classmethod
	def invert(self,SIGNAL,windows,partitions):
		SIG = numpy.array(SIGNAL)
		for p in range(len(partitions)-1):
			SIG[partitions[p]:partitions[p+1]] = numpy.fft.fft(SIG[partitions[p]:partitions[p+1]])
		SIG /= windows
		return numpy.fft.ifft(SIG)

	def interpolate(SIG,partitions,M=None,axes=None,kind='linear',**args):
		""" Interpolate a 1D GFT onto a grid. If axes is specified it should be a
			list or tuple consisting of two arrays, the sampling points on the time and frequency
			axes, respectively. Alternatively, M can be specified, which gives the number
			of points along each axis."""
		N = len(SIG)
		M = M if not M is None else N
		factor = N / M
		axes = [numpy.arange(0,N,factor),numpy.arange(0,N,factor)] if axes is None else axes
		newT = axes[0]		# New t axis
		widths = GFTPartitions.widths(partitions)
		spectrogram = []
		# interpolate each frequency band in time
		for p in range(len(partitions)):
			# indices of sample points, plus 3 extra on each side in case of cubic interpolation
			indices = (numpy.arange(-3,widths[p]+3))
			# time coordinates of samples
			t = indices * (N/widths[p])
			# values at sample points; indices can index multiple times into SIG[]
			f = SIG[partitions[p]:partitions[p+1]][indices % widths[p]] if p < len(partitions)-1 else SIG[partitions[p]:][indices % widths[p]]
			spectrogram.append(scipy.interpolate.interp1d(t,f,kind=kind,**args)(newT) if len(f) > 1 else f)
		
		spectrogram = numpy.array(spectrogram)
		
		# Interpolate in frequency
		newF = axes[1]		# New f axis
		indices = numpy.arange(-3,len(partitions)+3) % len(partitions)
		f = partitions[indices] + widths[indices]/2; f[0:3] -= N; f[-3:] += N
		t = spectrogram[indices]
		spectrogram = scipy.interpolate.interp1d(f,t,axis=0,kind=kind,**args)(newF)
		return spectrogram



# ****** Main ***********
from pylab import *; ion()
import PyGFT as gft
numpy.set_printoptions(precision=4,suppress=True)
N = 256
x = numpy.arange(N); sig = numpy.zeros(N,dtype=numpy.complex64); sig[N//2] = 1.0

windows,partitions = GFTPartitions.windows(N)
SIG = GFT.transform(sig,windows,partitions)
SIG1 = gft.gft1d(sig,'box')

sigR = GFT.invert(SIG,windows,partitions)

fig,ax = subplots(3,1,clear=True,num='Test')
ax[0].plot(x,sig.real,label='Original Signal',alpha=0.5)
ax[0].plot(x,sigR.real,label='Recovered Signal',alpha=0.5)
_ = [ax[1].axvline(p / (N-1) * x.max(),0,1,color='r',alpha=0.5,linestyle='--') for p in partitions]
ax[1].plot(x,abs(SIG),label='Pure Python',alpha=0.5)
#ax[1].plot(x,SIG1,label='Cython',alpha=0.5)
ax[0].legend(); ax[1].legend()

axes = [numpy.arange(0,N),numpy.arange(0,N)]
spectrogram = GFT.interpolate(SIG,partitions,axes=axes,kind='linear')
ax[2].imshow(abs(spectrogram),aspect='auto',origin='lower'); ax[2].set_title('Interpolated Spectrogram')



