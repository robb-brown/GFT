import numpy, scipy.interpolate, scipy.signal, numbers

from numpy.fft import fft, ifft, fftn, ifftn

gaussian = scipy.signal.windows.gaussian
pi = numpy.pi

# TODO: indexer and maybe interpolation functions need to be adjusted so they range from
# frequency indices 0:N/2,-N/2:0 rather than -N/2:0:N/2
# change interpolator so it takes - to + indices and adds N/2?

class GFTPartitions(object):
	""" Superclass for GFT window sets. Implements generic boxcar octave windows."""
	
	@classmethod 
	def widths(self,N,partitions):
		"""Takes a partitions array and returns the widths of the partitions"""
		widths = (numpy.roll(partitions,-1) - partitions); widths[-1] += N
		return widths
	
	@classmethod
	def dyadicPartitions(self,N):
		if N < 2:
			raise Exception('Cannot create dyadic partitions for signals < 2 samples.')
		# effectively rounds N to the nearest power of 2.
		widths = (2**numpy.arange(numpy.round(numpy.log(N)/numpy.log(2))-1)).astype(numpy.uint32)
		middle = len(widths)
		widths = numpy.concatenate([[1],widths,widths[::-1]])
		# If N is not a power of two then add or subtract the extra points from the highest frequency band
		extra = int(N - 2**numpy.round(numpy.log(N)/numpy.log(2)))
		if not (extra == 0):
			widths[middle] += extra // 2
			widths[middle+1] += extra - extra//2
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
				gaussian(widths[p],widths[p]*sigma,sym=symmetric)
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
		
	@classmethod
	def interpolate(self,SIG,partitions,M=None,axes=None,kind='linear',**args):
		""" Interpolate a 1D GFT onto a grid. If axes is specified it should be a
			list or tuple consisting of two arrays, the sampling points on the time and frequency
			axes, respectively. Alternatively, M can be specified, which gives the number
			of points along each axis. Indexing is [time,frequency]"""
		N = len(SIG)
		M = M if not M is None else N
		factor = N / M
		axes = [numpy.arange(0,N,factor),numpy.arange(0,N,factor)] if axes is None else axes
		
		# Roll indices so our actual indexing of -N/2 - 0 - N/2 looks like 0 - N for the
		# interpolator
		#axes[1] -= N/2
		
		newT = axes[0]		# New t axis
		widths = GFTPartitions.widths(N,partitions)
		spectrogramM = []
		spectrogramP = []
		extra = 3 if kind=='cubic' else 1
		# interpolate each frequency band in time
		# TODO: make this more efficient by only computing the bands we need
		for p in range(len(partitions)):
			# indices of sample points, plus extra points on each side depending on order of interpolation
			indices = (numpy.arange(-extra,widths[p]+extra))
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
		indices = numpy.arange(-extra,len(partitions)+extra) % len(partitions)
		f = partitions[indices] + widths[indices]/2;
		if extra > 0:
			f[0:extra] -= N; f[-extra:] += N
		spectrogramM = scipy.interpolate.interp1d(f,abs(spectrogramM[indices]),axis=0,kind=kind,**args)(newF)
		#spectrogram = spectrogramM * (cos(spectrogramP) + 1j*sin(spectrogramP))
		spectrogramP = scipy.interpolate.interp1d(f,spectrogramP[indices],axis=0,kind=kind,**args)(newF)
		spectrogram = spectrogramM * (numpy.cos(numpy.angle(spectrogramP)) + 1j*numpy.sin(numpy.angle(spectrogramP)))
		return spectrogram
	
	def __init__(self,N=None,signal=None,partitions=None,windows=None,sampleRate=None,interpolation='nearest'):
		""" 
		If sampleRate is None, indexing will be by sample number in both t and f.
		"""
		self.N = N = N if not N is None else len(signal) if not signal is None else None
		self.signal = signal
		self.partitions = partitions
		self.windows = windows
		self.sampleRate = sampleRate
		self.interpolation = interpolation
		self.SIG = None
		if not ((self.signal is None) or (self.partitions is None) or (self.windows is None)):
			self.SIG = self.transform(signal=self.signal,windows=self.windows,partitions=self.partitions)
		if not self.partitions is None and not self.sampleRate is None:
			p = numpy.array([i if (i < N/2) else -(N/2-(i-N/2)) for i in self.partitions])
			self.partitionsf = p * 1./self.sampleRate
			
	
	def partitionFromf(self,f):
		""" Given a frequency, figure out the index of the band it's in"""
		N = self.N
		f2 = numpy.array(f).reshape((-1))
		if not self.sampleRate is None:
			f2 = f2 * self.sampleRate / N
		f2 = numpy.array([i if (i >= 0) else i+N for i in f2])
		indices = numpy.searchsorted(self.partitions, f2+0.00000001, side="left")-1
		if not hasattr(f,'__len__'):
			return indices[0]
		else:
			return indices
	
	def sampleIndexFromt(self,t):
		""" Given a time coordinate, figure out the index of the sample """
		N = self.N
		t2 = numpy.array(t).reshape((-1))
		if not self.sampleRate is None:
			t2 = t2 * self.sampleRate
		if not hasattr(t,'__len__'):
			return t2[0]
		else:
			return t2
		
	def sampleIndexFromf(self,f):
		""" Given a frequency coordinate, figure out the index of the sample """
		N = self.N
		f2 = numpy.array(f).reshape((-1))
		if not self.sampleRate is None:
			f2 = f2 * self.sampleRate / N
		if not hasattr(f,'__len__'):
			return f2[0]
		else:
			return f2
			
	
	def tCoordinates(self,f):
		""" Given a frequency, get the coordinates of the middle of each time sample"""
		N = self.N
		f = self.partitionFromf(f)
		widths = GFTPartitions.widths(self.N,self.partitions)
		indices = numpy.arange(widths[f])
		t = indices * (self.N/widths[f])
		if not self.sampleRate is None:
			 t = t / self.sampleRate
		return t
		
	def fCoordinates(self):
		""" Get the frequency (f) coordinates of the middle of each frequency sample """
		N = self.N
		if self.partitions is None:
			return None
		else:
			widths = GFTPartitions.widths(self.N,self.partitions)
			indices = numpy.arange(len(self.partitions)) % len(self.partitions)
			f = self.partitions[indices]
			f = numpy.array([i if (i < N/2) else -(N/2-(i-N/2)) for i in f])
			fc = f + widths[indices]/2
			if not self.sampleRate is None:
				fc = fc * N/self.sampleRate
			return fc
	
	
	def __getitem__(self,key):
		""" 'Slicing'. Frequency first, time second. """
		interpArgs = dict(SIG=self.SIG,partitions=self.partitions,kind=self.interpolation)
		axes = None
		N = len(self.SIG)
		# Handle sigle index
		if (not hasattr(key,'__len__')):
			key = [key,slice(None)]
		if (len(key) == 1):
			key = [key[0],slice(None)]
		
		fi = key[0]; ti = key[1]
		# Convert scalar or slice coordinates into array
		import pudb; pu.db
		if isinstance(fi,numbers.Number):
			fi = [fi]
		elif isinstance(fi,slice):
			args = {}
			sampleSize = (N/self.sampleRate) if not self.sampleRate is None else 1.
			args['start'] = fi.start if not fi.start is None else -N/2 * sampleSize
			args['stop'] = fi.stop if not fi.stop is None else N/2 * sampleSize
			args['step'] = fi.step if not fi.step is None else sampleSize				
			fi = numpy.arange(**args)
		if isinstance(ti,numbers.Number):
			ti = [ti]
		elif isinstance(ti,slice):
			args = {}
			sampleSize = (1./self.sampleRate) if not self.sampleRate is None else 1.
			maxSample = N / self.sampleRate if not self.sampleRate is None else N
			args['start'] = ti.start if not ti.start is None else 0
			args['stop'] = ti.stop if not ti.stop is None else maxSample
			args['step'] = ti.step if not ti.step is None else sampleSize
			ti = numpy.arange(**args)
		
		# convert coordinates to indices (N-based, not partitions)
		fi = (self.sampleIndexFromf(fi) + (N+1)) % N
		ti = self.sampleIndexFromt(ti)
		
		# use interpolate to generate our value(s)
		axes = [ti,fi]
		result = self.interpolate(axes=axes,**interpArgs)
		if isinstance(key[0],numbers.Number) and isinstance(key[1],numbers.Number):
			# our coordinates were both scalar so we should return a scalar
			result = result[0,0]
		return result



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
	sig = numpy.cos(2*pi*N/12*x)*numpy.roll(gaussian(N,N/4),-N//4) + \
			numpy.cos(2*pi*N/6*x)*numpy.roll(gaussian(N,N/4),N//4)
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
