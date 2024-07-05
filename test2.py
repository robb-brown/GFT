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
sig = signalChirp(N)+1./N
partitions = GFTPartitions.dyadicPartitions(N)
windows = GFTPartitions.boxcarWindows(N,partitions)

gft = GFT(signal=sig,partitions=partitions,windows=windows,sampleRate=N)
self = gft

s = gft[0:256,0.:1.]

figure('test',clear=True)
imshow(numpy.abs(s),origin='lower');




print(s.argmax(axis=0)[0])


# cubic interpolation doesn't seem to be working right