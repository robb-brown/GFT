from gft import *
from pylab import *

numpy.set_printoptions(precision=4,suppress=True)

success = []
ns = []
for N in arange(2,1000):
	partitions = GFTPartitions.dyadicPartitions(N)
	success.append(partitions[-1])
	ns.append(N)

print(numpy.all((numpy.array(ns) - numpy.array(success)) == 1))



N = 2000
sig = signalChirp(N)+1
partitions = GFTPartitions.dyadicPartitions(N)
windows = GFTPartitions.boxcarWindows(N,partitions)

gft = GFT(signal=sig,partitions=partitions,windows=windows,sampleRate=None)
self = gft


spectrogram = GFT.interpolate(SIG=gft.SIG,partitions=gft.partitions,kind='nearest')

s = gft[:]


demo1(sig,windowShape='boxcar',interpolation='nearest')
figure('test')
s = gft[:]; 
clf(); imshow(log(abs(s)),origin='lower'); axhline(abs(s).argmax(axis=0)[0],color='r'); 
print(s.argmax(axis=0)[0])


# cubic interpolation doesn't seem to be working right