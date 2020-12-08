from math import *
from numpy import *
from numpy.fft import *
from matplotlib.pylab import *
from time import *

ion()

def shift(a,nx):
	nx = -int(round(nx))
	b = concatenate((a[nx:],a[0:nx]))
	return b

def downsampleft(a,factor):
	if (factor < 2): return a
	A = fft(a)
	N = len(A)
	return ifft(concatenate((A[0:N//(2*factor)],A[N-N//(2*factor):N+1])))

def fourierKernel(signal,k):
	x = arange(0,1,1./len(signal)).astype(complex64)
	ker = e**(1j*2*pi*k*x)
	return ker

def window(signal,t,k):
	l = len(signal)
	x = arange(0,1,1./l).astype(complex64)
	win = 1*abs(k)/sqrt(2*pi)*e**(-(x-0.5)**2*abs(k)**2/2.)
	win = win/(sum(win))
	win = shift(win,-l//2)
	win = shift(win,t)
	return win


def boxedWindow(signal,t,k,width):
	l = len(signal)
	k = float(k)
	x = arange(0,1,1./l).astype(complex64)
	win = 1*abs(k)/sqrt(2*pi)*e**(-(x-0.5)**2*abs(k)**2/2.)
	sigma = l/k
	w = int(round(sigma*width))
	if (w*2) < l:
		if (t == 0): print("Width: %d" % w)
		win[0:l//2-w] = 0
		win[l//2+w:] = 0
	else:
		if (t == 0): print("Width: %d" % l)
	win = win/(sum(win))
	win = shift(win,-l//2)
	win = shift(win,t)
	return win


def plotc(sig,color='b',prefix='',**args):
	plot(sig.real,'%c-' % color,label=prefix+'Real',**args); plot(sig.imag,'%c--' % color,label=prefix+'Imaginary',**args)



class GFT(object):
	
	windowmult = 1/0.8
	windowcut = 16.0
	
	@classmethod
	def transform(self,signal,kernel=fourierKernel,window=boxedWindow):
		transform = []
		transformn = []
		sig = array(signal)
		N = len(signal)
		
		k = len(sig) 			# will get divided by 2 before we start
		
		while k > 1:
			k = k // 2
			M = len(sig)
			sigma = k*self.windowmult
			ker = kernel(sig,k)
			kern = kernel(sig,-k)
			print("k= %d  M= %d" % (k,M))
			sk = sig*ker
			skn = sig*kern
			line = []
			linen = []
			for i in range(0,M):
				win = window(sig,i,sigma,self.windowcut)				
				line.append(sum(sk*win))		
				linen.append(sum(skn*win))

			transform.append(line)
			transformn.append(linen)
			if (k <= N//4 and M > 2):
				sig = downsampleft(sig,2)
		transform.append([sum(sig)])
		return (transform,transformn)

	@classmethod		
	def invert(self,transform,kernel=fourierKernel,window=boxedWindow):
		trans = transform[0]
		transn = transform[1]
		N = len(trans[0])
		M = len(trans)
		SIG = zeros(N,complex64)
			
		dk = N//4
		k = N//2
		pSIG = k
		pT = 0
		while (pT <= M-2):
			t = array(trans[pT])
			tn = array(transn[pT])
			win = window(t,0,k*self.windowmult,self.windowcut)

			WIN = fft(win)
			print("Line %d: k=%d dk= %d, pSIG= %d to %d, %d to %d" % (pT,k,dk,pSIG,pSIG-dk,-pSIG,-pSIG+dk))
			for i in range(0,dk+1):
				kernelcorrection = kernel(t,-i)
				kernelcorrectionn = kernel(t,i)
				windowcorrection = WIN[i]
				SIG[pSIG-i] = sum(t*kernelcorrection)/windowcorrection
				SIG[-pSIG+i] = sum(tn*kernelcorrectionn)/windowcorrection
			pT += 1
			pSIG -= dk
			dk //= 2
			k //= 2
		SIG[0] = sum(trans[len(trans)-1])				#	DC		
		return SIG

	
	@classmethod
	def interpolate(self,transform):
		N = len(transform[0])
		t = zeros((N//2,N),complex64)
		end = 0
		dk = len(transform[0]) // 4
		if (dk < 1): dk = 1
		for i in range(0,len(transform)):
			start = end
			end = start+dk
			t[start:end,:] = repeat(array(transform[i]),N//len(transform[i]))

			# Modify spectral weighting
			t[start:end,:] = t[start:end,:] * (len(transform[i]))			
			
			dk //= 2
		return t



if (__name__ == '__main__'):

	rc('figure',facecolor='w')
	
	# Create a fake signal
	N = 256
	x = arange(0,1,1./N)
	sig = zeros(len(x),complex64)
	sig[0:N//2] += (sin((N/16)*2*pi*x)*1.0)[0:N//2]
	sig[N//2:N] += (cos((N/8)*2*pi*x)*1.0)[N//2:N]
	sig[2*N//16:3*N//16] += (sin((N//4)*2*pi*x)*1.0)[2*N//16:3*N//16]
	sig[N//2+N//4] = 2.0

	# Do the fast ST... it returns the postive and negative parts separately
	ST = GFT.transform(sig)
	
	# Interpolate the postive part to make a spectrogram
	STi = GFT.interpolate(ST[0])
	figure('Fast ST Spectrogram'); imshow(abs(STi))
	
	# Compute the inverse.  Inverse ST is in the frequency domain
	STInverse = GFT.invert(ST)
	# inverse FFT and shift to recover the original signal
	sigST = shift(ifft(STInverse)[::-1],1)
	figure('Original Signal and Inverted Fast ST'); clf()
	plotc(sig,prefix='Signal ',color='b',alpha=0.5); plotc(sigST,prefix='Inverted ST ',color='r',alpha=0.5)
	legend()