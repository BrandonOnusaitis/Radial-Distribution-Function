#!/home/aog/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

N = 108 # number of LJ beads

nf = 100 # number of frames

L_min = np.zeros(3)-2.6 # negative box vectors (x,y,z)
L_max = np.zeros(3) + 2.6 # postive box vectors (x,y,z)

Vbox = (2.6*2)**3 # volume of the box

f = open("traj.txt", "r")

raw_traj = [] # unshaped container of trajectories

xyz = [] # container that will hold each frame 

line = f.readline()

while line != "":

    raw_traj.append(np.array(line.strip().split()).astype(float))

    line = f.readline()

f.close()

frame = []

for i in range(len(raw_traj)): 

    if (i+1)% 108 == 0 and (i != len(raw_traj)-1):

        frame.append(raw_traj[i]) 

        xyz.append(np.array(frame))

        frame = []
    
    elif i == (len(raw_traj)-1): 

        frame.append(raw_traj[i])

        xyz.append(np.array(frame))

        break

    else:

        frame.append(raw_traj[i]) 

xyz = np.array(xyz)

nbins = 200

rmax = 2.6

bin_width = rmax/nbins # the maximum length we will measure is L/2

hist = np.zeros(nbins)

test = []

for i in xyz: # loop through the frames

    dr = [] # this will be our pairwise separation container 

    for r in range(108): # loop through the rows of the frame matrix

        for c in range(r+1,108):# loop through column nums that > row num

            r_ij = i[r] - i[c] # computes pair separation

            bool_m = np.abs(r_ij) > 2.6 # creates a True False array that tests for what dimensions are greater than L/2

            r_ij[bool_m] = 5.2 - np.abs(r_ij[bool_m]) # applies minimum image convention to each dimension

            if np.linalg.norm(r_ij) > (np.linalg.norm(np.array([rmax,0,0]))): # tests whether the pair_separation greater than the allowed r val
               
                continue
            
            else: 

                dr.append(np.linalg.norm(r_ij))
    
    # now we need to loop through dr and assign its value to the bin
    
    dr = np.array(dr)

    dr = np.floor(dr/bin_width).astype(int) # gives the index of the bin that the pair separation falls in

    for b_i in dr: 

        hist[b_i] += 1 # we iterate through our dr and add 2 to the corresponding bin

# now we need to normalize our histogram and compure g(r)

g = [] # radial distribution function

r = []

for i in range(len(hist)): 

    if i < len(hist) -1: 

        volume_bin = (4/3)*np.pi*((bin_width)**3)*((i+1)**3 - i**3) # volume of spherical shell

        Prob = hist[i]/(N*(N-1)*(1/2)*nf) # probability of observing the prescribe density
       
        g.append(Prob*(Vbox/(volume_bin))) # radial distribution function assuming 1 atom in shell and vol
                                             # occupied by a single atom in the bulk
        r.append((bin_width/2)*(i+1 + i)) # append the midpoint of the bin 

rcParams['mathtext.fontset'] = 'cm'

fig,ax = plt.subplots()

ax.set_xlabel(r'$\rm{Pair \,\, Separation \,\,(\Delta R_{ij})}$', fontsize = 14)
ax.set_ylabel(r'$\rm{Radial \,\,Distribution \,\,Function,\,\, g(\Delta R_{ij})}$', fontsize =14)

ax.plot(r,g,color = 'indigo', marker = "+", linestyle= "-")

plt.show()