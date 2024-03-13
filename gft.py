import numpy, scipy.interpolate, scipy.signal

from numpy.fft import fft, ifft, fftn, ifftn

class GFTPartitions(object):
	""" Superclass for GFT window sets. Implements generic boxcar octave windows."""
	
	@classmethod 
	def widths(self,N,partitions):
		"""Takes a partitions array and returns the widths of the partitions"""
		widths = (numpy.roll(partitions,-1) - partitions); widths[-1] += N
		return widths
	
	@classmethod
	def dyadicPartitions(self,N):
		widths = (2**numpy.arange(numpy.round(numpy.log(N)/numpy.log(2))-1)).astype(numpy.uint32)
		widths = numpy.concatenate([[1],widths,widths[::-1]])		
		partitions = numpy.concatenate([[0],numpy.cumsum(widths)])
		return partitions
	
	@classmethod
	def boxcarWindows(self,N,partitions):
		# Boxcar windows are flat
		windows = numpy.ones(N)
		return windows
	
	@classmethod
	def gaussianWindows(self,N,partitions,sigma=0.5,symmetric=False):
		windows = numpy.zeros(N)
		widths = GFTPartitions.widths(N,partitions)
		for p in range(len(partitions)):
			windows[partitions[p]:partitions[p]+widths[p]] = \
				scipy.signal.gaussian(widths[p],widths[p]*sigma,sym=symmetric)
			windows[partitions[p]:partitions[p]+widths[p]] /= sum(windows[partitions[p]:partitions[p]+widths[p]]) / widths[p] 
		return windows
	

class GFT(object):

	@classmethod
	def transform(self,signal,windows,partitions):
		SIG = numpy.fft.fft(signal)
		SIG *= windows
		for p in range(len(partitions)-1):
			SIG[partitions[p]:partitions[p+1]] = numpy.fft.ifft(SIG[partitions[p]:partitions[p+1]])
		return SIG

	@classmethod
	def transform2(self,signal,windows,partitions):
		"""Not sure what this is trying to do. Maybe working on real GFT?"""
		N = len(signal)
		widths = GFTPartitions.widths(N,partitions)
		SIG = numpy.fft.fft(signal)
		SIG *= windows
		gft = numpy.zeros_like(SIG)
		p = 0
		SIG = numpy.roll(SIG,-N//2)
		while partitions[p] + widths[p] <= N//2:
			c = numpy.concatenate([SIG[N//2-partitions[p+1]:N//2-partitions[p]],SIG[N//2+partitions[p]:N//2+partitions[p+1]]])
			gft[partitions[p]*2:partitions[p+1]*2] = numpy.fft.ifft(numpy.roll(c,widths[p]))
			p += 1
		return gft
		
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
		widths = GFTPartitions.widths(N,partitions)
		spectrogramM = []
		spectrogramP = []
		# interpolate each frequency band in time
		for p in range(len(partitions)):
			# indices of sample points, plus 3 extra on each side in case of cubic interpolation
			indices = (numpy.arange(-3,widths[p]+3))
			# time coordinates of samples
			t = indices * (N/widths[p])
			# values at sample points; indices can index multiple times into SIG[]
			a = SIG[partitions[p]:partitions[p+1]][indices % widths[p]] if p < len(partitions)-1 else SIG[partitions[p]:][indices % widths[p]]
			spectrogramM.append(scipy.interpolate.interp1d(t,abs(a),kind=kind,**args)(newT) if len(a) > 1 else abs(a))
			spectrogramP.append(scipy.interpolate.interp1d(t,a,kind=kind,**args)(newT) if len(a) > 1 else a)
		
		spectrogramM = numpy.array(spectrogramM)
		spectrogramP = numpy.array(spectrogramP)
		
		# Interpolate in frequency
		newF = axes[1]		# New f axis
		indices = numpy.arange(-3,len(partitions)+3) % len(partitions)
		f = partitions[indices] + widths[indices]/2; f[0:3] -= N; f[-3:] += N
		spectrogramM = scipy.interpolate.interp1d(f,abs(spectrogramM[indices]),axis=0,kind=kind,**args)(newF)
		#spectrogram = spectrogramM * (cos(spectrogramP) + 1j*sin(spectrogramP))
		spectrogramP = scipy.interpolate.interp1d(f,spectrogramP[indices],axis=0,kind=kind,**args)(newF)
		spectrogram = spectrogramM * (numpy.cos(numpy.angle(spectrogramP)) + 1j*numpy.sin(numpy.angle(spectrogramP)))
		return spectrogram


class GFTND(object):
	""" N dimensional GFT. Work in progress."""
	
	@classmethod
	def transform(self,signal,windows,partitions):
		pass



def demo1(sig,windowShape='boxcar',interpolation='linear'):
	"""Compute and display the GFT for a demo signal"""
	N = len(sig)
	# Create the partitions and windows
	partitions = GFTPartitions.dyadicPartitions(N)
	windows = GFTPartitions.boxcarWindows(N,partitions) if windowShape == 'boxcar' else \
							GFTPartitions.gaussianWindows(N,partitions,sigma=0.25,symmetric=False)

	# Do the GFT
	SIG = GFT.transform(sig,windows,partitions)
	
	# Do the inverse GFT
	sigR = GFT.invert(SIG,windows,partitions)

	# Interpolate the GFT to get a spectrogram
	axes = [numpy.arange(0,N),numpy.arange(0,N)]
	spectrogram = GFT.interpolate(SIG,partitions,axes=axes,kind=interpolation)

	# Plot
	import pylab; pylab.ion()
	fig,ax = pylab.subplots(5,1,clear=True,num=f'GFT Demo, {windowShape} Windows',figsize=(6,8))
	ax[0].plot(sig.real,label='Original Signal',alpha=0.5)
	ax[0].plot(sigR.real,label='Recovered Signal',alpha=0.5)
	ax[0].set_xlim(0,N)
	_ = [ax[1].axvline(p / (N-1) * N,0,1,color='r',alpha=0.5,linestyle='--') for p in partitions]
	ax[1].plot(abs(fft(sig)),label='FFT Magnitude',alpha=0.5)
	ax[1].plot(numpy.angle(fft(sig)),label='FFT Phase',alpha=0.5)
	ax[1].plot(windows,label='Windows',alpha=0.5)
	ax[1].set_xlim(0,N)
	_ = [ax[2].axvline(p / (N-1) * N,0,1,color='r',alpha=0.5,linestyle='--') for p in partitions]
	ax[2].plot(abs(SIG),label='GFT Magnitude',alpha=0.5)
	ax[2].plot(numpy.angle(SIG),label='GFT Phase',alpha=0.5)
	ax[2].set_xlim(0,N)
	ax[3].imshow(abs(spectrogram),aspect='auto',origin='lower'); ax[3].set_title('Spectrogram Magnitude')
	ax[4].imshow(numpy.angle(spectrogram),aspect='auto',origin='lower'); ax[4].set_title('Spectrogram Phase')
	ax[0].legend(); ax[1].legend(); ax[2].legend()
	ax[3].set_axis_off(); ax[4].set_axis_off()
	pylab.tight_layout()
	

def demoWindows(sig):
	"""
	Show the effect of window selection on the GFT of a demo signal.
	The width of Gaussian windows is scaled by sigma, demonstrating the tradeoff
	between frequency and time resolution. The windows are also symmetric or asymmetric;
	for even N, symmetric Gaussians place their peak between sample points while asymmetric
	always have their peak on a sample. Symmetric windows are shifted by half a pixel
	relative to asymmetric, which manifests as a phase ramp in the S spectrum.
	"""
	import pylab; pylab.ion()
	N = len(sig)
	axes = [numpy.arange(0,N),numpy.arange(0,N)]
	partitions = GFTPartitions.dyadicPartitions(N)
	fig,ax = pylab.subplots(4,4,clear=True,num='Demo Windows')
	for r,sigma in enumerate([0.1,0.25,0.5,1000]):
		for c,symmetric in enumerate([True,False]):
			windows = GFTPartitions.gaussianWindows(N,partitions,sigma=sigma,symmetric=symmetric)
			SIG = GFT.transform(sig,windows,partitions)
			spectrogram = GFT.interpolate(SIG,partitions,axes=axes,kind='linear')
			ax[r,c*2].imshow(abs(spectrogram),aspect='auto',origin='lower'); ax[r,c*2].set_title(f'σ = {sigma} {"symmetric" if symmetric else ""}')
			ax[r,c*2+1].imshow(numpy.angle(spectrogram),aspect='auto',origin='lower'); ax[r,c*2+1].set_title(f'Phase')
			ax[r,c*2].set_axis_off(); ax[r,c*2+1].set_axis_off()
	pylab.tight_layout()


def signalDelta(N):
	sig = numpy.zeros(N,dtype=numpy.complex64); sig[N//2] = 1.0
	return sig

def signalDoubleDelta(N):
	sig = numpy.zeros(N,dtype=numpy.complex64); 
	sig[N//4] = 1.0; sig[3*N//4] = 1.0
	return sig


def signalTwoTone(N):
	x = numpy.arange(0,N) / N
	sig = numpy.cos(2*pi*N/12*x)*numpy.roll(scipy.signal.gaussian(N,N/4),-N//4) + \
			numpy.cos(2*pi*N/6*x)*numpy.roll(scipy.signal.gaussian(N,N/4),N//4)
	return sig

def signalChirp(N):
	x = numpy.arange(0,N) / N
	sig = numpy.cos(2*pi*N/4*x**2)
	return sig


def runDemo():
	numpy.set_printoptions(precision=4,suppress=True)
	sig = signalDoubleDelta(N = 512)
	demo1(sig,windowShape='boxcar',interpolation='linear')
	demo1(sig,windowShape='gaussian',interpolation='linear')
	demoWindows(sig)


if __name__ == '__main__':
	runDemo()
