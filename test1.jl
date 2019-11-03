#!/usr/bin/julia
#
# This module provides tools for constructing and testing
# the frequency shifting stacked halfband filter decimators
# described in "An efficient channelizer tree for portable
# software defined radios" by Fred Harris & Elettra Venosa &
# Xiaofei Chen & Chris Dick
#
# Googling the paper title I discovered there are two papers
# with this name and these authors, one published in Ann.
# Telecommun. 2014 69:1-2, 99-110, and a simpler one perhaps
# as presented at the 14th International Symposium on Wireless
# Personal Multimedia Communications, WPMC 2011, Brest, France,
# October 3-7, 2011
# Both papers are associated with DOI 10.1007/s12243-013-0354-y.
#
# The prototype low pass filter coefficients are spread over
# four paths.  path0 has the center coefficient of the filter
# and is otherwise zero.  path2 is all zero. path1 and path3
# split the rest of the non-zero coefficients in the filter
# path/band band0 band1 band2 band3
# path0     1      1      1     1
# path1     1      im    -1    -im
# path2     1     -1      1    -1
# path3     1     -im    -1     im

using Plots
using FFTW
using DSP
using Printf

# Generate complex signals that cover the frequencies specified in ks
# t contains the time points in the sample
function iqsignals(ks, t)
    iqsignal(k, t) = cos.(2π * k .* t) + sin.(2π * k .* t)im
    signal = iqsignal(ks[1], t)
    for i = ks[2:end]; signal += iqsignal(i,t); end
    signal ./= length(ks)
end

# Generate a spread of signals over fs by df
iqspread(fs, df, t) = iqsignals(df:df:fs, t)

# Generate the target signals over fs
iqtargets(fs, t) = iqsignals(0:fs/4:3*fs/4, t)

# Generate the target signals over fs with labeling
iqlabels(fs, t) = iqsignals([0, fs/4,fs/4+100, 2fs/4,2fs/4+100,2fs/4+200, 3fs/4,3fs/4+100,3fs/4+200,3fs/4+300], t)

# Generate a spread of signals over fs by df excluding channel b
iqblank(fs, df, b, t) = iqsignals((b*fs/4+fs/8):df:((b+4)*fs/4-fs/8), t)

# compute a spectrogram and plot it
function spectrogram(signal, N, Ts, title)
    # Fourier Transform of signal
    F = fft(signal) |> fftshift
    freqs = fftfreq(N, 1.0/Ts) |> fftshift
    plot(freqs, abs.(F), title = title, xlim=(-6500, +6500), legend=:none)
end

# compute the coefficient combining path p into band k
# note that path 2 is always 0
function coeff(p,k)
    simplifyimag(z) = if abs(imag(z)) > 0 && abs(imag(z)) < 1e-15; Complex(real(z), 0) else z end
    simplifyreal(z) = if abs(real(z)) > 0 && abs(real(z)) < 1e-15; Complex(0, imag(z)) else z end
    simplify(z) = simplifyreal(simplifyimag(z))
    simplify(exp(im * 2 * π / 4 * p * k))
end

function main(args)

    N = 2^14 - 1                # Number of points in sampled signal
    fs = 12000                  # Sample frequency
    Ts = 1 / fs                 # Sample period
    t0 = 0                      # Start time
    tmax = t0 + N * Ts          # end time
    t = t0:Ts:tmax              # time coordinate

    # initial signal
    #signal = iqsignals(0:fs/4:3*fs/4, t)
    # signal = iqspread(fs, 100, t)
    # signal = iqtargets(fs, t)
    # signal = iqlabels(fs, t)
    # signal = iqblank(fs, 100, 0, t)
    signal = iqblank(fs, 100, 1, t)

    # coefficients of paths and bands
    coeffs = [
        [ coeff.(0,0:3) ], # path 0
        [ coeff.(1,0:3) ], # path 1
        [ coeff.(2,0:3) ], # path 2
        [ coeff.(3,0:3) ], # path 3
    ]
    #show(coeffs)

    # filter coefficients and name strings and polynomial representation
    # there's confusion about whether c[1] or c[3] belong to path 1 or 3
    # this layout centers band 0 at 0, band 1 at fs/4, band 2 at 2fs/4, and
    # band 3 at 3fs/4
    j = im
    c = [ [+1/4, +1/2, +1/4],
          [-j/4, +1/2, +j/4],
          [-1/4, +1/2, -1/4],
          [+j/4, +1/2, -j/4] ]
    n = [ "+1/4 +1/2 +1/4",
          "-j/4 +1/2 +j/4",
          "-1/4 +1/2 -1/4",
          "+j/4 +1/2 -j/4" ]
    c = [ [+3/8, +1/2, +3/8],
          [-3j/8, +1/2, +3j/8],
          [-3/8, +1/2, -3/8],
          [+3j/8, +1/2, -3j/8] ]
    n = [ "+3/8 +1/2 +3/8",
          "-j3/8 +1/2 +j3/8",
          "-3/8 +1/2 -3/8",
          "+j3/8 +1/2 -j3/8" ]
    f = [ DSP.Filters.PolynomialRatio(c[1], [1.0]),
          DSP.Filters.PolynomialRatio(c[2], [1.0]),
          DSP.Filters.PolynomialRatio(c[3], [1.0]),
          DSP.Filters.PolynomialRatio(c[4], [1.0]) ]

    # filtered signals
    s = [ filt(f[1], signal),
          filt(f[2], signal),
          filt(f[3], signal),
          filt(f[4], signal) ]

    # downsample filtered signals
    Tsd = 2*Ts
    td = t0:Tsd:tmax
    sd = [ s[1][1:2:N],
           s[2][1:2:N],
           s[3][1:2:N],
           s[4][1:2:N] ]

    # spectrogram of initial signal
    spec = spectrogram(signal, length(t), Ts, "Spectrum")

    # spectrograms of filtered signals
    f = [ spectrogram(s[1], length(t), Ts, n[1] * " /1"),
          spectrogram(s[2], length(t), Ts, n[2] * " /1"),
          spectrogram(s[3], length(t), Ts, n[3] * " /1"),
          spectrogram(s[4], length(t), Ts, n[4] * " /1") ]

    # spectrograms of filtered downsampled signals
    fd = [ spectrogram(sd[1], length(td), Tsd, n[1] * " /2"),
           spectrogram(sd[2], length(td), Tsd, n[2] * " /2"),
           spectrogram(sd[3], length(td), Tsd, n[3] * " /2"),
           spectrogram(sd[4], length(td), Tsd, n[4] * " /2") ]

    plot(spec, f[1], f[2], f[3], f[4], layout = (5,1))
    savefig("unfolded.pdf")
    plot(spec, fd[1], fd[2], fd[3], fd[4], layout = (5,1))
    savefig("folded.pdf")

end

main(ARGS)
