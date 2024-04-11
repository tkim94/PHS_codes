#%load_ext autoreload

#%autoreload 2
#%matplotlib inline

from matplotlib import rcParams

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from math import pi
import healpy as hp

# import classy module
from classy import Class

# Random sampling
import random

from random import seed
from random import gauss
from random import choice
from random import sample
from random import shuffle

from matplotlib import rc

import pywt
from functools import partial
import math

#rc('font',**{'family':'serif','serif':['Helvetica']})
#rc('text', usetex=True)

batch_num = int(sys.argv[1])
batch_init = int(40*(batch_num-1))
batch_end = int(40*batch_num)

# Number of maps per a package
#nclass = 1000
nclass = 500
#nclass = 1

# N_pix = 12*nside^2
#nside = 256
nside = 1024
#nside = 512

# pixel size in radian for nside.
Pixel_Size = hp.pixelfunc.nside2resol(nside)

# A number of pixels ( in total sphere )
Npixels = 12*nside**2

# The central position to generate hotspots
#                  theta    phi
vec_center = hp.ang2vec(np.pi / 2,  0  )

# A degree of window to look up the sky
Window_degree = 2.58
Wsize = (np.pi/180)*Window_degree

# Defining a square
vertex_1 = hp.ang2vec(np.pi / 2 + Wsize, 0 + Wsize)
vertex_2 = hp.ang2vec(np.pi / 2 + Wsize, 0 - Wsize)
vertex_3 = hp.ang2vec(np.pi / 2 - Wsize, 0 - Wsize)
vertex_4 = hp.ang2vec(np.pi / 2 - Wsize, 0 + Wsize)

square_array = np.array([vertex_1,vertex_2,vertex_3,vertex_4])

square_indices = hp.query_polygon(nside, square_array)

# Number of pixels for x and y coordinates
Car_Dim = int(math.sqrt(len(square_indices))) 
print(str(Car_Dim))


# Select to window
# Range in latitude. Default: [-90,90]
latra = [-Window_degree, Window_degree]
        
# Range in longitude. Default: [-180,180]
lonra = [-Window_degree, Window_degree]

for rep in range(batch_init, batch_end):
     # CMS only Maps
    Map_bkg = np.empty(( 0,  Car_Dim, Car_Dim  ))

    # Loop over each map
    for c in range(nclass):

        print("nclass = " + str(c) )


        #__________________________________________________
        # Projection scheme
        proj = hp.projector.CartesianProj(
            lonra=lonra, latra=latra,
            coord='G',
            xsize=Car_Dim, ysize=Car_Dim)
        #__________________________________________________


        # Max cl for power spectrum (class)
        ClMax = 3500

        LambdaCDM = Class()

        LambdaCDM.set({'omega_b':0.022032,
                        'omega_cdm':0.12038,
                        'h':0.67556,
                        'A_s':2.215e-9,
                        'n_s':0.9619,
                        'tau_reio':0.0925})

        LambdaCDM.set({'output':'tCl,pCl,lCl,mPk',
                        'lensing':'yes',
                        'P_k_max_1/Mpc':3.0})
        # Running Class
        LambdaCDM.compute()

        # get all C_l output
        cls = LambdaCDM.lensed_cl(ClMax)
        # To check the format of cls

        ll = cls['ell'][2:]
        clTT = cls['tt'][2:]


        # Generating a map from PSD
        Map_Combo = hp.synfast((2.7255*10**6)**2*clTT,nside)


        # Projecting the CMB images into the Cartesian coordinates.
        # CMB only
        Combo_im_bkg = proj.projmap(Map_Combo, vec2pix_func=partial(hp.vec2pix, nside))
        #print(str(Combo_im_bkg.shape))


        # Reshaping the images..
        # CMB only
        Combo_im_bkg = np.reshape(Combo_im_bkg, (1, Car_Dim, Car_Dim) )


        # Stacking images... 
        Map_bkg = np.vstack((Map_bkg, Combo_im_bkg))
        
    file_loc = '/afs/crc.nd.edu/user/t/tkim12/Work/CMB_ML/Data/Nside1024/bkg/'
    file_name = str(nclass)+'_Events_BKG_'+str(Car_Dim)+'by'+str(Car_Dim)+'_'
    np.save(file_loc+file_name+str(rep), Map_bkg)
    

#Map_bkg.tofile("text.csv",sep=',',format='%10.5f')
#plt.imshow(Map_bkg[0])
    # End of the Loop for range(nclass)
