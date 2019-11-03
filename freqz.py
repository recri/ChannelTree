# The following is a Python/scipy snippet to evaluate the 
# coefficients for a halfband filter.  A halfband filter
# is a filter where the cutoff frequency is Fs/4 and every
# even coeffecient is zero except the center tap.
# Note: every other (even except 0) is 0, most of the coefficients
#       will be close to zero, force to zero actual

import sys
import numpy
from numpy import log10, abs, pi
import scipy
from scipy import signal
import matplotlib
import matplotlib.pyplot
import matplotlib as mpl

bands = numpy.array([0., .22, .28, .5]) # typical transition band for fs/4 halfband lowpass

h = []
name = sys.argv[1]
for x in sys.argv[2:]:
    h.append(eval(x))

N = len(h)-1  # Filter order

# Dump the coefficients for comparison and verification
print('          coeff')
print('------------------------------------')
for ii in range(N+1):
    print(' tap %2d   %-3.6f' % (ii, h[ii]))

# compute the frequency response
(w,H) = signal.freqz(h)

## ~~[Plotting]~~
# Note: the pylab functions can be used to create plots,
#       and these might be easier for beginners or more familiar
#       for Matlab users.  pylab is a wrapper around lower-level
#       MPL artist (pyplot) functions.

fig = mpl.pyplot.figure()
ax1 = fig.add_subplot(111)
ax1.plot(w, 20*log10(abs(H)), color='b')
ax1.plot(-w, 20*log10(abs(H)), color='b')
ax1.plot(w-pi, 20*log10(abs(H)), color='r')
ax1.plot(-w+pi, 20*log10(abs(H)), color='r')
#bx = bands*2*pi
#ax1.axvspan(bx[1], bx[2], facecolor='0.5', alpha=0.33)
ax1.plot(pi/2, -6, 'go')
for s in [-1,1]:
    ax1.axvline(s*pi/2, color='g', linestyle='--')
    ax1.axvline(s*pi/4, color='g', linestyle='--')
    ax1.axvline(s*pi/8, color='g', linestyle='--')
ax1.axis([-pi,pi,-124,3])
ax1.grid('on')
ax1.set_ylabel('Magnitude (dB)')
ax1.set_xlabel('Normalized Frequency (radians)')
ax1.set_title(name + ' frequency response')
fig.savefig('rsp-'+name+'.png')
