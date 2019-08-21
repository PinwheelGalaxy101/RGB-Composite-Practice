# Python3 program to preprocess an image of Saturn

# load dark images
# use arithmetic mean of dark images as dark frame for all images
from astropy.io import fits
sd_list = [ fits.open('/Users/kawaii/Documents/obs/190619/Saturn_Dark_0.1_'+n+'.fits')[0].data[0] for n in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10'] ]
import numpy as np
sd = (1/10)*np.sum(sd_list, 0)

# load r, v, and b images and subtract dark frame
sr_list = [ fits.open('/Users/kawaii/Documents/obs/190619/Saturn_R_0.1_0'+n+'.fits')[0].data[0] - sd for n in['1', '2', '3'] ]
sv_list = [ fits.open('/Users/kawaii/Documents/obs/190619/Saturn_V_0.1_0'+n+'.fits')[0].data[0] - sd for n in['1', '2', '3'] ]
sb_list = [ fits.open('/Users/kawaii/Documents/obs/190619/Saturn_B_0.1_0'+n+'.fits')[0].data[0] - sd for n in['1', '2', '3'] ]

# load flat dark images for each filter
# use arithmetic means of flat dark images as flat dark frames
frd_list = [ fits.open('/Users/kawaii/Documents/obs/190620/flat_rdark_0.05_0'+n+'.fits')[0].data[0] for n in ['1', '2', '3', '4', '5'] ]
fvd_list = [ fits.open('/Users/kawaii/Documents/obs/190620/flat_vdark_0.08_0'+n+'.fits')[0].data[0] for n in ['1', '2', '3', '4', '5'] ]
fbd_list = [ fits.open('/Users/kawaii/Documents/obs/190620/flat_bdark_0.3_0'+n+'.fits')[0].data[0] for n in ['1', '2', '3', '4', '5'] ]
frd = (1/5)*np.sum(frd_list, 0)
fvd = (1/5)*np.sum(fvd_list, 0)
fbd = (1/5)*np.sum(fbd_list, 0)

# load flat images and subtract flat dark frames
fr_list = [ fits.open('/Users/kawaii/Documents/obs/190620/flat_r_0.05_0'+n+'.fits')[0].data[0] for n in ['1', '2', '3', '4', '5'] ]
fv_list = [ fits.open('/Users/kawaii/Documents/obs/190620/flat_v_0.08_0'+n+'.fits')[0].data[0] for n in ['1', '2', '3', '4', '5'] ]
fb_list = [ fits.open('/Users/kawaii/Documents/obs/190620/flat_b_0.3_0'+n+'.fits')[0].data[0] for n in ['1', '2', '3', '4', '5'] ]
fr = (1/5)*np.sum(fr_list, 0) - frd
fv = (1/5)*np.sum(fv_list, 0) - fvd
fb = (1/5)*np.sum(fb_list, 0) - fbd

# normalize flat images by setting median pixel to 1.0
fb = fb/np.median(fb)
fv = fv/np.median(fv)
fr = fr/np.median(fr)

# adjust r, v, and b images for sensitivity using flat images
for n in range(3):
	sr_list[n] = sr_list[n]/fr
	sv_list[n] = sv_list[n]/fv
	sb_list[n] = sb_list[n]/fb

"""
# find position of saturn in each image
# the first r image is used as reference
from scipy import signal
saturn, pos = [], []
for n in range(3):
	saturn.append(np.copy(sr_list[n][900:1150,1100:1600]))
for n in range(3):
	saturn.append(np.copy(sv_list[n][900:1150,1100:1600]))
for n in range(3):
	saturn.append(np.copy(sb_list[n][900:1150,1100:1600]))
for sat in saturn[1:]:
	corr = signal.correlate2d(sat,saturn[0], boundary='symm', mode='same')
	pos.append(np.unravel_index(np.argmax(corr, axis=None), corr.shape))
# pos is a list of eight tuples that contain
# the position of saturn in each cropped image
print(pos)
"""

# values of pos determined from previous run
pos = [(118, 186), (115, 214), (121, 275), (114, 247), (118, 265), (101, 236), (97, 242), (100, 235)]


# crop images around position of saturn
sr_list[0] = sr_list[0][900:1150,1100:1600]
for n in range(2):
	yr, xr = pos[n][0], pos[n][1]
	sr_list[n+1] = sr_list[n+1][yr+775:yr+1025,xr+850:xr+1350]
for n in range(3):
	yv, xv = pos[n+2][0], pos[n+2][1]
	sv_list[n] = sv_list[n][yv+775:yv+1025,xv+850:xv+1350]
	yb, xb = pos[n+5][0], pos[n+5][1]
	sb_list[n] = sb_list[n][yb+775:yb+1025,xb+850:xb+1350]

# combine r, v, and b images
from astropy.visualization import make_lupton_rgb
sr = np.sum(sr_list,0)
sv = np.sum(sv_list,0)
sb = np.sum(sb_list,0)
sr = sr + np.min(sr)
sv = sv + np.min(sv)
sb = sb + np.min(sb)
scale = 20/np.max([sr, sv, sb])
sr = 0.9*scale*sr
sv = 0.5*scale*sv
sb = scale*sb
saturn_image = make_lupton_rgb(sr, sv, sb)

# test plot
"""
import matplotlib.pyplot as plt
plt.figure()
plt.imshow(saturn_image)
plt.colorbar()
plt.show()
"""
"""
for n in range(3):
	plt.figure()
	plt.imshow(sr_list[n], cmap='gray')
	plt.colorbar()
	plt.show()
for n in range(3):
	plt.figure()
	plt.imshow(sv_list[n], cmap='gray')
	plt.colorbar()
	plt.show()
for n in range(3):
	plt.figure()
	plt.imshow(sb_list[n], cmap='gray')
	plt.colorbar()
	plt.show()
"""


# normalize image with zscale and plot
import matplotlib.pyplot as plt
from astropy.visualization import (ImageNormalize, ZScaleInterval)
plt.figure()
plt.imshow(saturn_image, norm=ImageNormalize(saturn_image, ZScaleInterval()))
plt.show()
