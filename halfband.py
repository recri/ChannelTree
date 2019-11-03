# The following is a Python/scipy snippet to generate the 
# coefficients for a halfband filter.  A halfband filter
# is a filter where the cutoff frequency is Fs/4 and every
# other coeffecient is zero except the cetner tap.
# Note: every other (even except 0) is 0, most of the coefficients
#       will be close to zero, force to zero actual

import numpy
from numpy import log10, abs, pi
import scipy
from scipy import signal
import matplotlib
import matplotlib.pyplot
import matplotlib as mpl

# ~~[Filter Design with Parks-McClellan Remez]~~
N = 10  # Filter order
# Filter symetric around 0.25 (where .5 is pi or Fs/2)
# bands = numpy.array([0., .22, .28, .5])    # original guess
bands = numpy.array([0., 0.125, 0.375, 0.5]) # harris recommended
h = signal.remez(N+1, bands, [1,0], [1,1])
h[abs(h) <= 1e-4] = 0.
(w,H) = signal.freqz(h)

# ~~[Filter Design with Windowed freq]~~
b = signal.firwin(N+1, 0.25)
b[abs(h) <= 1e-4] = 0.
(wb, Hb) = signal.freqz(b)

# ~~[Filter Design by fiat]
# 1/4 1/2 1/4, easy shift and add
f121 = [0.0, 0.25, 0.5, 0.25, 0.0]
(wf121, Hf121) = signal.freqz(f121)
# 3/8 1/2 3/8, less easy shift and add
f323 = [0.0, 0.375, 0.5, 0.375, 0.0]
(wf323, Hf323) = signal.freqz(f323)

# Dump the coefficients for comparison and verification
if N == 4:
    print('          remez       firwin        1/4 1/2 1/4        3/8 1/2 3/8')
    print('------------------------------------')
    for ii in range(N+1):
        print(' tap %2d   %-3.6f    %-3.6f    %-3.6f    %-3.6f' % (ii, h[ii], b[ii], f121[ii], f323[ii]))
else:
    for ii in range(N+1):
        print(' tap %2d   %-3.6f' % (ii, h[ii]))

## ~~[Plotting]~~
# Note: the pylab functions can be used to create plots,
#       and these might be easier for beginners or more familiar
#       for Matlab users.  pylab is a wrapper around lower-level
#       MPL artist (pyplot) functions.
if False:
    fig = mpl.pyplot.figure()
    ax0 = fig.add_subplot(211)
    ax0.stem(numpy.arange(len(h)), h)
    ax0.grid(True)
    ax0.set_title('Parks-McClellan (remez) Impulse Response')
    ax1 = fig.add_subplot(212)
    ax1.stem(numpy.arange(len(b)), b)
    ax1.set_title('Windowed Frequency Sampling (firwin) Impulse Response')
    ax1.grid(True)
    ax2 = fig.add_subplot(212)
    ax2.stem(numpy.arange(len(f121)), f121)
    ax2.set_title('1/4 1/2 1/4 Frequency Sampling Impulse Response')
    ax2.grid(True)
    fig.savefig('hb_imp.png')

fig = mpl.pyplot.figure()
ax1 = fig.add_subplot(111)
ax1.plot(w, 20*log10(abs(H)))
ax1.plot(w, 20*log10(abs(Hb)))
ax1.plot(w, 20*log10(abs(Hf121)))
ax1.plot(w, 20*log10(abs(Hf323)))
ax1.legend(['remez', 'firwin', '1/4 1/2 1/4', '3/8 1/2 3/8'])
bx = bands*2*pi
ax1.axvspan(bx[1], bx[2], facecolor='0.5', alpha=0.33)
ax1.plot(pi/2, -6, 'go')
ax1.axvline(pi/2, color='g', linestyle='--')
ax1.axis([0,pi,-64,3])
ax1.grid('on')
ax1.set_ylabel('Magnitude (dB)')
ax1.set_xlabel('Normalized Frequency (radians)')
ax1.set_title('Half Band Filter Frequency Response')
fig.savefig('hb_rsp.png')
