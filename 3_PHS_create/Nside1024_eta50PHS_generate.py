from random import randint
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import sys

batch_num = int(sys.argv[1])

### Definitions ###

def angle_to_eta(theta):
    return LSS_radius*theta#/np.sqrt(4*np.pi)

def LSS_loc(r):
    loc = -eta+grid_len_r*r # Find distance from LSS (negative = inside, positive = outside)
    return loc

def dist_2spot(r1,r2,x,y):    # flag:  if separation is large, the curvature is involved so flat space approximation may fail
    delta_r = (r2-r1)*grid_len_r # in Mpc
    x_Mpc = x*grid_len_x # in Mpc
    y_Mpc = y*grid_len_x # in Mpc
    dist = np.sqrt(delta_r**2+x_Mpc**2+y_Mpc**2)
    return dist

def Get_Dis(x1,y1, x2,y2):
    return np.math.sqrt( (x1-x2)**2 + (y1-y2)**2 )

def interpolation(x, y):
    return interpolate.interp1d(x, y)

### Definitions ###

num_sigs = 500

nside = 1024
LSS_radius = 13800

# A number of pixels
Npixels = 12*nside**2

# pixel size analytic(Rad)
pix_radsq = 4*np.pi/Npixels # Solid angle
pix_rad = np.sqrt(pix_radsq) # Square pixel side length in rad

pix_mpc = angle_to_eta(pix_rad) # (Mpc) 1 pixel size ~ 6.7 Mpc

eta = 50
r_bins = 50 # No need to change
grid_len_r = eta*2/r_bins # Mpc (each grid length)
grid_len_x = pix_mpc # Mpc = grid_len_y



ang_arr = np.logspace(-3.7, -0.7, num=100)
radius_arr = angle_to_eta(ang_arr) # Mpc
pix_arr = (0.2+np.arange(0,25))* grid_len_x # Mpc, Offset is introduced to mimic center pixel
#print(radius_arr[0]) 2.69 Mpc
#print(pix_arr[0])    2.70 Mpc

SW_angle_arr = np.loadtxt('SW_stack_eta50.txt')
ISW_angle_arr = np.loadtxt('ISW_stack_eta50.txt')



# After this can be in for loop if were to create larger samples
Map_dim = 90

Signal_Map = np.empty(( 0,  Map_dim, Map_dim  ))

for trials in range(num_sigs):

    print(trials)
    
    r1 = randint(1,49) # Spot 1 depth
    #r1 = randint(20,34)  # Spot 1 depth constrained

    #r1 = 25

    for rep in range(20000000):
        r2 = randint(1,50)
        #r2 = 25

        deltax = randint(-100,100)
        deltay = randint(-100,100)

        dist = dist_2spot(r1,r2,deltax,deltay)

        if dist<= eta: # distance smaller than eta*
            print("Distance between PHS = "+str(round(dist,2))+" Mpc ")
            print(r1,r2,deltax,deltay)
            print("Last Scattering Surface Spot 1: "+str(round(LSS_loc(r1),2))+"Mpc , Spot 2: "+str(round(LSS_loc(r2),2))+" Mpc")
            break

    dist = dist_2spot(r1,r2,deltax,deltay)
    if dist > eta:
        sys.exit()
    
    
    '''
    # Alternate method v2

    def spherical_dist(r1, r2, theta):
        r1full = LSS_radius+(r1-25)*grid_len_r
        r2full = LSS_radius+(r2-25)*grid_len_r
        dist = np.sqrt( (r2full*np.cos(theta)-r1full)**2 + (r2full*np.sin(theta))**2 )
        return dist

    r1 = randint(20,34)

    for rep in range(2000):
        r2 = randint(1,50)

        theta = random.uniform(0,160/13000)
        phi = random.uniform(0,2*np.pi)

        dist = spherical_dist(r1,r2,theta)

        if dist<= eta: # distance smaller than eta*
            deltax = int(LSS_radius*np.sin(theta)*np.cos(phi)/grid_len_x)
            deltay = int(LSS_radius*np.sin(theta)*np.sin(phi)/grid_len_x)

            print(dist)
            print(r1,r2,deltax,deltay)
            print("Last Scattering Surface Spot 1:"+str(round(LSS_loc(r1),2))+"Mpc , Spot 2:"+str(round(LSS_loc(r2),2))+" Mpc")
            break

    '''       

    pix_temp1 = []
    pix_temp2 = []

    r1_value = r1 
    r2_value = r2


    for rad in range(len(pix_arr)):
        SWprofle1 = interpolation(radius_arr, SW_angle_arr[r1_value])
        ISWprofle1 = interpolation(radius_arr, ISW_angle_arr[r1_value])

        SWprofle2 = interpolation(radius_arr, SW_angle_arr[r2_value])
        ISWprofle2 = interpolation(radius_arr, ISW_angle_arr[r2_value])

        pix_temp1.append(SWprofle1(pix_arr[rad])+ISWprofle1(pix_arr[rad]))
        pix_temp2.append(SWprofle2(pix_arr[rad])+ISWprofle2(pix_arr[rad]))


    # Create hot spot profile

    SingleSpot1 = np.zeros((101, 101))
    SingleSpot2 = np.zeros((101, 101))


    nbin = len(pix_temp1)

    # Loop over each bin for a profile function (shape of a hot spot)
    for j in range(nbin):
        # A list of positions for each shell of profile
        x_pos = np.empty(0)

        y_pos = np.empty(0)

        # Delta(T) of a given pixel
        result1 = pix_temp1[j]
        result2 = pix_temp2[j]
        #print(str(result))

        # Collecting a list of positions for a given shell
        for ii in range (0,101):
            for jj in range (0,101):
                if j <= Get_Dis(51, 51, ii, jj) and Get_Dis(51, 51, ii, jj) < j+1:
                    x_pos = np.append(x_pos, ii)
                    y_pos = np.append(y_pos, jj)

        size_x_pos = len(x_pos)


        # Injecting a Delta T into the given shell
        for ii in range (0,size_x_pos):

            x_inject = int(x_pos[ii])
            y_inject = int(y_pos[ii])

            SingleSpot1[x_inject][y_inject] = SingleSpot1[x_inject][y_inject] + result1
            SingleSpot2[x_inject][y_inject] = SingleSpot2[x_inject][y_inject] + result2    


    framedim = 200
    frame = np.zeros((framedim, framedim))
    # 1st signal boundary

    shiftX = int(framedim/2)
    shiftY = int(framedim/2)

    spot2X = shiftX+deltax
    spot2Y = shiftY+deltay

    s1Xmin = shiftX-50
    s1Xmax = shiftX+51
    s1Ymin = shiftY-50
    s1Ymax = shiftY+51


    #2nd signal boundary
    s2Xmin = spot2X-50
    s2Xmax = spot2X+51
    s2Ymin = spot2Y-50
    s2Ymax = spot2Y+51

    frame[s1Xmin:s1Xmax, s1Ymin:s1Ymax] += SingleSpot1
    frame[s2Xmin:s2Xmax, s2Ymin:s2Ymax] += SingleSpot2

    frame_nonzero_Xmin, frame_nonzero_Xmax = np.nonzero(frame)[0].min(), np.nonzero(frame)[0].max()+1
    frame_nonzero_Ymin, frame_nonzero_Ymax = np.nonzero(frame)[1].min(), np.nonzero(frame)[1].max()+1

    PHS_signal = frame[frame_nonzero_Xmin:frame_nonzero_Xmax, frame_nonzero_Ymin:frame_nonzero_Ymax]

    PHS_shape = PHS_signal.shape

    print("Window Size = ["+str(PHS_signal.shape[0])+", "+str(PHS_signal.shape[1])+"]")

    framedim = 90
    frame = np.zeros((framedim, framedim))
    
    random_shiftX = randint(1,framedim-PHS_shape[0]-1)
    random_shiftY = randint(1,framedim-PHS_shape[1]-1)

    print("Random shift = ("+str(random_shiftX)+", "+str(random_shiftY)+")")

    frame[random_shiftX:random_shiftX+PHS_signal.shape[0],random_shiftY:random_shiftY+PHS_signal.shape[1]]+=PHS_signal
    
    frame = np.reshape(frame, (1, Map_dim, Map_dim) )
    
    Signal_Map = np.vstack((Signal_Map, frame))


    
file_loc = '/afs/crc.nd.edu/user/t/tkim12/Work/CMB_ML/Data/Nside1024/PHS_signal/'
file_name = str(num_sigs)+'_eta50PHS'

#print(Signal_Map.shape)

np.save(file_loc+file_name+'_g1_Sig_offLSS_'+str(batch_num), Signal_Map)
    #plt.imshow(frame)
    #plt.show()