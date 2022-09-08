'''

# A python module for microseismic location using
# arrival-time picking based method, and
# migration style method, i.e semblance, linear stacking, semblance weighted stacking  
#
# (C) Zhengguang Zhao, 2016

'''

# !/Users/Uqer/anaconda/bin/python

import os
from math import floor
import numpy as np
from numpy import array, ones, diff, square, cumsum, zeros, linspace, append, exp,\
    savetxt, loadtxt, argmin, abs
import matplotlib.pyplot as plt
from psmodules import psarray, pssynthetic, psraytrace, pswavelet, \
    psplot, pspicker, pspdf
# from psmodules import psarray
# import psarray
# import pssynthetic
# import psraytrace
# import pswavelet
# import psplot
# import pspicker
# import pspdf


def main():
    # Generate square grid array
    geox, geoy, geoz = psarray.gridarray(25, 30000, 30000)

    # plot geophone
    psarray.plotarray(geox, geoy, tit='Grid Array', xlen=30000, ylen=30000, offset=.0)

    # geox = linspace(0, 6000, 31)
    # geoy = linspace(0, 0, 31)
    # geoz = geoy + 2.0
    nt = len(geox)

    # Define source coordinates
    sourcex = array([15000])
    sourcey = array([15000])
    sourcez = array([3000])

    # Define geological model
    zlayer = array([0, 540, 1070, 1390, 1740, 1950, 2290,
                    2630, 4000])
    # zlayer_ss103h = array([0, 10, 30, 110, 210, 315, 395, 568, 821, 923, 1177, 1270, 1331, 1383, 1593, 1656, 1746,
    #                        2155, 2301, 2488, 2701, 2863, 2974, 3114, 3200])

    # Define velocity model
    # P wave velocity
    vp = array([2100, 2500, 2950, 3300, 3700, 4200,
                4700, 5800])
    # vp_ss103h = array([880.00, 1200.00, 1732.00, 1732.00, 1732.00, 2122.00, 2193.00, 2550.00, 2522.00,
    #                    2709.00, 2793.00, 2588.00, 2596.00, 3047.00, 3170.00, 3061.00, 2537.00, 3573.00, 3655.00, 4205.00, 4345.00,
    #                    4486.00, 4629.00, 4460.00, 4843.00])
    # vs_ss103h = array([400, 545, 787.00, 787.00, 870.00, 1060.00, 1100.00, 1272.00, 1256.00, 1350.00,
    #                    1470.00, 1362.00, 1366.00, 1603.00, 1668.00, 1611.00, 1270.00, 1985.00, 2030.00, 2336.00, 2508.00, 2590.00,
    #                    2673.00, 2575.00, 2796.00])

    #  S wave velocity based on Castagna's rule
    # vs = (vp - 1360) / 1.16
    #  S wave velocity based on literature:
    #  2019 Optimal design of microseismic monitoring network: Synthetic study for the Kimberlina CO2 storage demonstration site
    vs = vp / 1.73


    # 3D passive seismic raytracing
    print("3D passive seismic raytracing example is running[Waiting...]")
    dg = 10
    # ptimes, pthetas = psraytrace.raytracing(
    #     vp_ss103h, vs_ss103h, zlayer_ss103h, dg, sourcex, sourcey, sourcez, geox, geoy, geoz)
    ptimes, pthetas = psraytrace.raytracing(
        vp, vs, zlayer, dg, sourcex, sourcey, sourcez, geox, geoy, geoz)

    print("3D passive seismic raytracing completed[OK]")
    # print(pthetas)
    # print("Trave times:")
    # print(ptimes)

    # Generate wavelet
    # oscillator wavelet
    # osciwlet, tw = pswavelet.oscillator(0.0005, 65, 0.2555, 3, 3, 1, 80, 50)
    osciwlet, tw = pswavelet.oscillator(0.002, 65, 1.022, 3, 3, 1, 80, 50)

    # Generate synthetic microseismogram
    # surface acquisition
    ns = 5000
    dt = 0.002
    # # downhole acquisition
    # ns = 600
    # dt = 0.0005

    data = zeros((ns, nt), dtype='float32')
    att = ones((nt, 1), dtype='float32')
    syndata = pssynthetic.genSynNoiseFree(
        ns, nt, osciwlet, pthetas, ptimes, dt, att)

    # nt = 62
    # Plot synthetic traces
    psplot.hseisplot(syndata, ns, nt)
    psplot.vseisplot(syndata, ns, nt)

    # Pick and plot first arrival time of microseismic events

    # pick arrival times
    pics = []
    for i in range(nt):
        tr = syndata[:, i]
        pic = pspicker.merpicker(tr, ns, 20, 600, "False")
        pics.append(pic)
    pickers = array(pics, dtype='float32')
    pickers.shape = (len(pickers), 1)

    # plot pickers
    psplot.hseispickplot(syndata, pickers, ns, nt)
    psplot.vseispickplot(syndata, pickers, ns, nt)

    # Locate microseismic event
    x1 = 1000
    x2 = 30000
    dx = 10
    nx = floor((x2 - x1) / dx)

    y1 = 1000
    y2 = 30000
    dy = 10
    ny = floor((y2 - y1) / dy)

    z1 = 100
    z2 = 3900
    dz = 10
    nz = floor((z2 - z1) / dz)

    n = 1
    t0 = 0
    sigma = 1

    minErrs = []
    for i in range(x1, x2, dx):
        for j in range(y1, y2, dy):
            for k in range(z1, z2, dz):
                sx = array([i])
                sy = array([j])
                sz = array([k])

                print('sx,sy,sz:')
                print(sx)
                print(sy)
                print(sz)
                # tps, tetas = psraytrace.raytracing(
                #     vp_ss103h, vs_ss103h, zlayer_ss103h, dg, sx, sy, sz, geox, geoy, geoz)
                tps, tetas = psraytrace.raytracing(
                    vp, vs, zlayer, dg, sx, sy, sz, geox, geoy, geoz)
                tps = tps / dt

                minErr = timediff(tps, pickers)
                minErrs.append(minErr)

    # pdfs = array(pdfs)
    minErrs = array(minErrs)
    minErrIndex = argmin(minErrs)
    print(minErrIndex)
    ind = findIndex(nx, ny, 1, minErrIndex)

    print(ind)
    print('source x: ',  range(x1, x2, dx)[ind[0]])
    print('source y: ', range(y1, y2, dy)[ind[1]])


def timediff(tp, tmp):
    """
    timediff(): microseismic event location method based on arrival time differences. 
    This method needs to pick first arrival times of microseismic event and 
    generally aims to process high signal-to-noise ratio.

    """
    tpdiff = abs(diff(tp, axis=0))
    tmpdiff = abs(diff(tmp, axis=0))

    temp = square(tpdiff - tmpdiff)
    sumErrs = cumsum(temp)
    minErr = sumErrs[len(sumErrs) - 1]
    return minErr


def time():
    """
    time(): microseismic event location method based on arrival times. 
    This method needs to pick first arrival times of microseismic event and 
    generally aims to process high signal-to-noise ratio.

    """
    print('time()')


def sws():
    """
    sws(): semblance weighted stacking for microseismic event location. 
    This method doesn't need to pick first arrival times of microseismic event and 
    generally aims to process low signal-to-noise ratio.

    """


def semblance():
    """
    linearstack(): semblance for microseismic event location. 
    This method doesn't need to pick first arrival times of microseismic event and 
    generally aims to process low signal-to-noise ratio.

    """


def linearstack():
    """
    linearstack(): linear stacking for microseismic event location. 
    This method doesn't need to pick first arrival times of microseismic event and 
    generally aims to process low signal-to-noise ratio.

    """


def findIndex(nx, ny, nz, ind):
    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                value = k * nx * ny + i * ny + j
                if value == ind:
                    return np.array([i, j, k])

# This will actually run the code if called stand-alone:
if __name__ == '__main__':
    main()

# %% single forward check
import os
from math import floor
import numpy as np
from numpy import array, ones, diff, square, cumsum, zeros, linspace, append, exp,\
    savetxt, loadtxt, argmin, abs
import matplotlib.pyplot as plt
from psmodules import psarray, pssynthetic, psraytrace, pswavelet, \
    psplot, pspicker, pspdf

def timediff(tp, tmp):
    """
    timediff(): microseismic event location method based on arrival time differences.
    This method needs to pick first arrival times of microseismic event and
    generally aims to process high signal-to-noise ratio.

    """
    tpdiff = abs(diff(tp, axis=0))
    tmpdiff = abs(diff(tmp, axis=0))

    temp = square(tpdiff - tmpdiff)
    sumErrs = cumsum(temp)
    minErr = sumErrs[len(sumErrs) - 1]
    return minErr


def time():
    """
    time(): microseismic event location method based on arrival times.
    This method needs to pick first arrival times of microseismic event and
    generally aims to process high signal-to-noise ratio.

    """
    print('time()')


def sws():
    """
    sws(): semblance weighted stacking for microseismic event location.
    This method doesn't need to pick first arrival times of microseismic event and
    generally aims to process low signal-to-noise ratio.

    """


def semblance():
    """
    linearstack(): semblance for microseismic event location.
    This method doesn't need to pick first arrival times of microseismic event and
    generally aims to process low signal-to-noise ratio.

    """


def linearstack():
    """
    linearstack(): linear stacking for microseismic event location.
    This method doesn't need to pick first arrival times of microseismic event and
    generally aims to process low signal-to-noise ratio.

    """


def findIndex(nx, ny, nz, ind):
    for k in range(nz):
        for i in range(nx):
            for j in range(ny):
                value = k * nx * ny + i * ny + j
                if value == ind:
                    return np.array([i, j, k])

# Generate square grid array
geox, geoy, geoz = psarray.gridarray(25, 30000, 30000)

# plot geophone
psarray.plotarray(geox, geoy, tit='Grid Array', xlen=30000, ylen=30000, offset=.0)

nt = len(geox)

# Define source coordinates
sourcex = array([15000])
sourcey = array([15000])
sourcez = array([3000])

# Define geological model
zlayer = array([0, 540, 1070, 1390, 1740, 1950, 2290,
                2630, 4000])


# Define velocity model
# P wave velocity
vp = array([2100, 2500, 2950, 3300, 3700, 4200,
            4700, 5800])


#  S wave velocity based on Castagna's rule
# vs = (vp - 1360) / 1.16
#  S wave velocity based on literature:
#  2019 Optimal design of microseismic monitoring network: Synthetic study for the Kimberlina CO2 storage demonstration site
vs = vp / 1.73


# 3D passive seismic raytracing
print("3D passive seismic raytracing example is running[Waiting...]")
dg = 10
# ptimes, pthetas = psraytrace.raytracing(
#     vp_ss103h, vs_ss103h, zlayer_ss103h, dg, sourcex, sourcey, sourcez, geox, geoy, geoz)
src = np.array([sourcex,sourcey,sourcez]).T
rcv = np.array([geox,geoy,geoz]).T

ptimes,_,pthetas = psraytrace.raytrace(vp, vs, zlayer, dg, src, rcv)


print("3D passive seismic raytracing completed[OK]")
# print(pthetas)
# print("Trave times:")
# print(ptimes)

# Generate wavelet
# oscillator wavelet
# osciwlet, tw = pswavelet.oscillator(0.0005, 65, 0.2555, 3, 3, 1, 80, 50)
osciwlet, tw = pswavelet.oscillator(0.002, 65, 1.022, 3, 3, 1, 80, 50)

# Generate synthetic microseismogram
# surface acquisition
ns = 8000
dt = 0.002
# # downhole acquisition
# ns = 600
# dt = 0.0005

data = zeros((ns, nt), dtype='float32')
att = ones((nt, 1), dtype='float32')
syndata = pssynthetic.genSynNoiseFree(
    ns, nt, osciwlet, pthetas, ptimes, dt, att)

# nt = 62
# Plot synthetic traces
psplot.hseisplot(syndata, ns, nt)
psplot.vseisplot(syndata, ns, nt)

# Pick and plot first arrival time of microseismic events

# pick arrival times
pics = []
for i in range(nt):
    tr = syndata[:, i]
    pic = pspicker.merpicker(tr, ns, 20, 600, "False")
    pics.append(pic)
pickers = array(pics, dtype='float32')
pickers.shape = (len(pickers), 1)

# plot pickers
psplot.hseispickplot(syndata, pickers, ns, nt)
psplot.vseispickplot(syndata, pickers, ns, nt)


# Locate microseismic event
x1 = 5000
x2 = 20000
dx = 5000
nx = floor((x2 - x1) / dx)

y1 = 5000
y2 = 20000
dy = 5000
ny = floor((y2 - y1) / dy)

z1 = 1000
z2 = 4000
dz = 1000
nz = floor((z2 - z1) / dz)

n = 1
t0 = 0
sigma = 1

minErrs = []
for i in range(x1, x2, dx):
    for j in range(y1, y2, dy):
        for k in range(z1, z2, dz):
            sx = array([i])
            sy = array([j])
            sz = array([k])

            print('sx,sy,sz:')
            print(sx)
            print(sy)
            print(sz)
            # tps, tetas = psraytrace.raytracing(
            #     vp_ss103h, vs_ss103h, zlayer_ss103h, dg, sx, sy, sz, geox, geoy, geoz)

            try_src = np.array([sx, sy, sz]).T

            tps, _, tetas = psraytrace.raytrace(vp, vs, zlayer, dg, try_src, rcv)

            # tps, tetas = psraytrace.raytracing(
            #     vp, vs, zlayer, dg, sx, sy, sz, geox, geoy, geoz)
            tps = tps / dt

            minErr = timediff(tps, pickers)
            minErrs.append(minErr)

# pdfs = array(pdfs)
minErrs = array(minErrs)
minErrIndex = argmin(minErrs)
print(minErrIndex)
ind = findIndex(nx, ny, nz, minErrIndex)

print(ind)
print('source x: ',  range(x1, x2, dx)[ind[0]])
print('source y: ', range(y1, y2, dy)[ind[1]])
print('source z: ', range(z1, z2, dz)[ind[2]])