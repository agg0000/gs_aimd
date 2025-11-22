#!/usr/bin/env python3
'''
Random-velocities generation [DOI: 10.1063/1.3175799]
'''

import sys
import numpy as np

import softapi
from constant import evtoau, mas, amu2au

def center(x, natom, mass):
    mcenter = np.zeros(3)
    for i in range(natom):
        mcenter += x[i] * mass[i]

    tmass = np.sum(mass)
    mcenter /= tmass

    for i in range(natom):
        x[i] -= mcenter

    return x
    
def get_vector(ekin, mass, coord, n):
    ekin *= evtoau
    
    size = coord.shape
    natom = size[0]

    tmass = np.sum(mass)

    coord = center(coord, natom, mass)

    Ti = np.zeros((3, 3))
    for i in range(natom):
        Ti[0, 0] += mass[i] * (coord[i, 1] ** 2 + coord[i, 2] ** 2)
        Ti[1, 1] += mass[i] * (coord[i, 0] ** 2 + coord[i, 2] ** 2)
        Ti[2, 2] += mass[i] * (coord[i, 0] ** 2 + coord[i, 1] ** 2)

    for i in range(natom):
        Ti[0, 1] -= mass[i] * (coord[i, 0] * coord[i, 1])
        Ti[0, 2] -= mass[i] * (coord[i, 0] * coord[i, 2])
        Ti[1, 2] -= mass[i] * (coord[i, 1] * coord[i, 2])

    Ti[1, 0] = Ti[0, 1]
    Ti[2, 0] = Ti[0, 2]
    Ti[2, 1] = Ti[1, 2]

    velc_tot = np.random.normal(size = [n, natom, 3])

    mom_tot = []
    for velc in velc_tot:
        for i in range(natom):
            velc[i] /= np.sqrt(mass[i])

        velc = center(velc, natom, mass)
        p_mom = np.zeros_like(velc)
        for i in range(natom):
            p_mom[i] = velc[i] * mass[i]

        L_ang = np.zeros(3)
        for i in range(natom):
            L_ang += np.cross(coord[i], p_mom[i])

        w_ang = np.linalg.solve(Ti, L_ang)

        for i in range(natom):
            velc[i] -= np.cross(w_ang, coord[i])

        
        Ev_tot = 0.0
        for i in range(natom):
            Ev_tot += mass[i] * np.dot(velc[i], velc[i])

        a = np.sqrt(2 * ekin / Ev_tot)
        velc *= a

        mom = np.zeros_like(velc)
        for i in range(natom):
            mom[i] = velc[i] * mass[i]

        mom_tot.append(mom)

    return mom_tot

infile = sys.argv[1]
fname = infile.split('.')[0]

inp = open(infile).read().lower().splitlines()
fd = {}
for keyword in inp:
	ka = keyword.split('=')
	fd[ka[0]] = ka[1]

soft_split = fd['softuse'].split() if 'softuse' in fd else 'gaussian'
softuse = soft_split[0]

numele = int(fd['numele'])
numcycle = int(fd['numcycle'])
filename = fd['filename'] if 'filename' in fd else fname
nstate = int(fd['nstate']) if 'nstate' in fd else 1
istate = int(fd['istate']) if 'istate' in fd else 0
mcenter = int(fd['mcenter']) if 'mcenter' in fd else 0
theory = fd['theory'] if 'theory' in fd else 'td'
keyname = fd['keyname'].lower() if 'keyname' in fd else 'i'
flip = fd['flip'].lower() if 'flip' in fd else 'd'

Ekin = float(fd['ekin']) if 'ekin' in fd else 0.4
nsample = int(fd['nsample']) if 'nsample' in fd else 200

quansoft = softapi.softname(softuse)

pos0, symb = quansoft.getpos(filename, numele)
coord = np.array(pos0, dtype = float)
mass = [mas[s] for s in symb]

mom_tot = get_vector(Ekin, mass, coord, nsample)

for n in range(nsample):
    newname = "{}{}.inf".format(fname, n)

    mom0 = mom_tot[n]
    with open(newname, 'w') as newf:
        for line in inp:
            newf.write(line + '\n')

        for i, m in enumerate(mom0):
            newf.write("p{0}={0}, {1[0]:15.8f},{1[1]:15.8f},{1[2]:15.8f}\n".format(i, m))

