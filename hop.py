#!/usr/bin/env python
'''
for hop algorithm
'''
import os
import math
import random
import datetime
import numpy as np

import softapi
import writefile
from constant import evtoau, mas

def verletx(pos1, far1, mom1, symb, numele, dt):
	pos1 = np.array(pos1, dtype = float)
	far1 = np.array(far1, dtype = float)
	mom1 = np.array(mom1, dtype = float)

	for i in range(numele):
		mass = float(mas[symb[i]])
		x1 = pos1[i] + mom1[i] * dt / mass + far1[i] * dt ** 2 / (2 * mass)
		pos1[i] = x1
	pos0 = tuple(map(tuple, pos1))

	return pos0

#-------------------------------------------------------------------------------------

def verletp(mom1, far1, far2, dt):
	mom1 = np.array(mom1, dtype = float)
	far1 = np.array(far1, dtype = float)
	far2 = np.array(far2, dtype = float)

	mom1 = mom1 + dt * (far1 + far2) / 2
	mom0 = tuple(map(tuple, mom1))

	return mom0

#-------------------------------------------------------------------------------------

def beemanx(pos1, far1, far2, mom1, symb, numele, dt):
	pos1 = np.array(pos1, dtype = float)
	far1 = np.array(far1, dtype = float)
	far2 = np.array(far2, dtype = float)
	mom1 = np.array(mom1, dtype = float)

	for i in range(numele):
		mass = float(mas[symb[i]])
		x1 = pos1[i] + mom1[i] * dt / mass + (4 * far1[i] - far2[i]) * dt ** 2 / (6 * mass)
		pos1[i] = x1
	pos0 = tuple(map(tuple, pos1))

	return pos0

#-------------------------------------------------------------------------------------

def beemanp(mom1, far1, far2, far3, dt):
	mom1 = np.array(mom1, dtype = float)
	far1 = np.array(far1, dtype = float)
	far2 = np.array(far2, dtype = float)
	far3 = np.array(far3, dtype = float)

	mom1 = mom1 + dt * (2 * far1 + 5 * far2 - far3) / 6
	mom0 = tuple(map(tuple, mom1))

	return mom0

#-------------------------------------------------------------------------------------

def xcenter(pos1, symb, numele):
	pos1 = np.array(pos1, dtype = float)
	mass = np.array([mas[n] for n in symb], dtype = float)
	tmass = np.sum(mass)

	x1 = np.zeros_like(pos1)
	for i in range(numele):
		x1[i] = pos1[i] * mass[i]

	mceter = np.sum(x1, axis = 0) / tmass
	for i in range(numele):
		x1[i] = pos1[i] - mceter

	pos0 = tuple(map(tuple, x1))

	return pos0

#-------------------------------------------------------------------------------------

def pcenter(mom1, symb, numele):
	mom1 = np.array(mom1, dtype = float)
	tmom = np.sum(mom1, axis = 0)

	mass = np.array([mas[n] for n in symb], dtype = float)
	tmass = np.sum(mass)

	tv = tmom / tmass
	p1 = np.zeros_like(mom1)
	for i in range(numele):
		p1[i] = mom1[i] - mass[i] * tv

	mom0 = tuple(map(tuple, p1))

	return mom0

#-------------------------------------------------------------------------------------

def calkine(pp, symb):
	pp = np.array(pp, dtype = float)
	kine = 0
	for i in range(len(pp)):
		k = pp[i] * pp[i] / (mas[symb[i]] * 2)
		kine += np.sum(k)

	return kine

#-------------------------------------------------------------------------------------

def shopping(filename, outname, state, quansoft, deltav, intcycle, numele, jumpword, totalv, posa, fara, moma, symb, start, theory, spin):
	deltav = deltav * evtoau
	mom2 = np.array(moma[-2], dtype = float)
	pos1 = np.array(posa[-1], dtype = float)
	pos2 = np.array(posa[-2], dtype = float)
	pos3 = np.array(posa[-3], dtype = float)

	gfile1 = quansoft.getname(filename + str(intcycle))
	gfile2 = quansoft.getname(filename + str(intcycle - 1))
	gfile3 = quansoft.getname(filename + str(intcycle - 2))

	jumps = state + 1 if jumpword == 'upper' else state - 1

	keyw = 'g' if state == 0 else 'r'

	inp1 = quansoft.jpfile(gfile1, outname, symb, jumps, state, jumpword, theory, spin, False, start)
	inp2 = quansoft.jpfile(gfile2, outname, symb, jumps, state, jumpword, theory, spin, False, start)
	inp3 = quansoft.jpfile(gfile3, outname, symb, jumps, state, jumpword, theory, spin, False, start)

	oup1 = quansoft.runsoft(inp1, jumps, keyw)
	oup2 = quansoft.runsoft(inp2, jumps, keyw)
	oup3 = quansoft.runsoft(inp3, jumps, keyw)

	nfar1, _, _ = quansoft.parseout(oup1, numele, outname, jumps, [], start)
	nfar2, _, _ = quansoft.parseout(oup2, numele, outname, jumps, [], start)
	nfar3, _, _ = quansoft.parseout(oup3, numele, outname, jumps, [], start)

	nfar1 = np.array(nfar1, dtype = float)
	nfar3 = np.array(nfar3, dtype = float)
	ifar1 = np.array(fara[-1], dtype = float)
	ifar3 = np.array(fara[-3], dtype = float)

	lowerf = np.zeros_like(nfar1)
	upperf = np.zeros_like(nfar3)

	coefx = pos1 - pos3
	dx1 = pos2 - pos1
	dx3 = pos2 - pos3

	for i in range(numele):
		for j in range(3):
			if coefx[i][j]:
				lowerf[i][j] = (nfar3[i][j] * dx1[i][j] - ifar1[i][j] * dx3[i][j]) / coefx[i][j]
				upperf[i][j] = (ifar3[i][j] * dx1[i][j] - nfar1[i][j] * dx3[i][j]) / coefx[i][j]

	fn = []
	sn = []
	for i in range(numele):
		fvector = (lowerf[i] - upperf[i]) ** 2 / mas[symb[i]]
		#svector = np.sqrt(fvector)
		svector = (lowerf[i] - upperf[i]) / math.sqrt(mas[symb[i]])

		fn.append(fvector)
		sn.append(svector)

	fn = np.array(fn, dtype = float)
	sn = np.array(sn, dtype = float)
	ff = np.sum(fn)
	fs = fn / ff
	ss = sn / math.sqrt(ff)

	mp = []
	for i in range(numele):
		nn = ss[i] / math.sqrt(np.sum(fs[i]))
		dp = np.sum(mom2[i] * nn)
		pn = dp * nn

		mp.append(pn)

	mp = np.array(mp, dtype = float)

	nmom2 = np.zeros_like(mom2)

	tp = calkine(mp, symb)
	ps = 0

	asquare = 0
	bsquare = 0
	if tp >= deltav or jumpword == 'lower':
		k = math.sqrt(1 - deltav / tp) if jumpword == 'upper' else math.sqrt(1 + deltav / tp)
		nmom2 = mom2 + (k - 1) * mp

		ep = totalv[-2] + tp
	
		vx = deltav / 2
		ex = totalv[-2] + vx if jumpword == 'upper' else totalv[-2] - vx

		asquare = ff / (16 * vx ** 3)
		bsquare = (ep - ex) / (2 * vx)

		deno = math.sqrt(bsquare ** 2 + 1) if (np.sum(lowerf * upperf)) > 0 else math.sqrt(abs(bsquare ** 2 - 1))
		ps = math.exp((-math.pi / math.sqrt(16 * asquare)) * math.sqrt(2 / (bsquare + deno)))

	pr = random.random()

	p = ps if ps > pr else 0.0

	nmom2 = tuple(map(tuple, nmom2))

	return p, nmom2, asquare, bsquare

#-------------------------------------------------------------------------------------

def tchopping(filename, outname, state, quansoft, dsquare, isoc, soc, intcycle, numele, totalv, mom2, symb, start, theory, spin):
	cfile = quansoft.getname(filename + str(intcycle))
	coup = quansoft.runsoft(cfile, state, 'g')
	currf, curre, _ = quansoft.parseout(coup, numele, outname, state, [], start)
	
	tfile = quansoft.jpfile(cfile, outname, symb, isoc, state, 'trans', theory, spin, True, start)
	toup = quansoft.runsoft(tfile, state, 'g')
	tranf, trane, _ = quansoft.parseout(toup, numele, outname, state, [], start)

	currf = np.array(currf, dtype = float)
	tranf = np.array(tranf, dtype = float)
	jumpword = 'upper' if curre[0] < trane[0] else 'lower'

	deltav = curre[0] - trane[0]
	fn = []
	sn = []
	for i in range(numele):
		svector = (currf[i] - tranf[i]) / math.sqrt(mas[symb[i]])
		fvector = svector ** 2

		fn.append(fvector)
		sn.append(svector)

	fn = np.array(fn, dtype = float)
	sn = np.array(sn, dtype = float)
	ff = np.sum(fn)
	fs = fn / ff
	ss = sn / math.sqrt(ff)

	mp = []
	for i in range(numele):
		nn = ss[i] / math.sqrt(np.sum(fs[i]))
		dp = np.sum(mom2[i] * nn)
		pn = dp * nn

		mp.append(pn)

	mp = np.array(mp, dtype = float)

	nmom2 = np.zeros_like(mom2)

	tp = calkine(mp, symb)
	ps = 0

	asquare = 0
	bsquare = 0
	if tp >= deltav or jumpword == 'lower':
		k = math.sqrt(1 - deltav / tp) if jumpword == 'upper' else math.sqrt(1 + deltav / tp)
		nmom2 = mom2 + (k - 1) * mp

		ep = totalv[-1] + tp
	
		vx = soc
		ex = (curre[0] + trane[0]) / 2

		asquare = ff / (16 * vx ** 3)
		bsquare = (ep - ex) / (2 * vx)

		deno = math.sqrt(bsquare ** 2 + 1) if (np.sum(currf * tranf)) > 0 else math.sqrt(abs(bsquare ** 2 - 1))
		ps = 1 - math.exp((-math.pi / math.sqrt(16 * asquare)) * math.sqrt(2 / (bsquare + deno)))
		'''
		sqbs = math.sqrt(2 / (bsquare + math.sqrt(bsquare ** 2 + 1)))
		detp = math.pi * sqbs / ( 8 * math.sqbs(asquare))
		ps = 1 - math.exp(-2 * detp) if dsquare * detp else 1 - (math.exp(2 * dsquare * detp - 2 * detp) - 1) / (math.exp(2 * dsquare * detp) - 1)
		'''

	pr = random.random()

	p = ps if ps > pr else 0.0

	nmom2 = tuple(map(tuple, nmom2))

	return p, nmom2, asquare, bsquare

#-------------------------------------------------------------------------------------

def tphopping(filename, outname, state, quansoft, dsquare, isoc, soc, intcycle, numele, totalv, mom2, symb, start, theory, spin):
	cfile = quansoft.getname(filename + str(intcycle))
	coup = quansoft.runsoft(cfile, state, 'g')
	currf, curre, _ = quansoft.parseout(coup, numele, outname, state, [], start)
	
	tfile = quansoft.jpfile(cfile, outname, symb, isoc, state, 'trans', theory, spin, True, start)
	toup = quansoft.runsoft(tfile, state, 'g')
	tranf, trane, _ = quansoft.parseout(toup, numele, outname, state, [], start)

	currf = np.array(currf, dtype = float)
	tranf = np.array(tranf, dtype = float)
	jumpword = 'upper' if curre[0] < trane[0] else 'lower'

	deltav = curre[0] - trane[0]
	fn = []
	sn = []
	for i in range(numele):
		svector = (currf[i] + tranf[i]) / math.sqrt(mas[symb[i]])
		fvector = svector ** 2

		fn.append(fvector)
		sn.append(svector)

	fn = np.array(fn, dtype = float)
	sn = np.array(sn, dtype = float)
	ff = np.sum(fn)
	fs = fn / ff
	ss = sn / math.sqrt(ff)

	mp = []
	for i in range(numele):
		nn = ss[i] / math.sqrt(np.sum(fs[i]))
		dp = np.sum(mom2[i] * nn)
		pn = dp * nn

		mp.append(pn)

	mp = np.array(mp, dtype = float)

	nmom2 = np.zeros_like(mom2)

	tp = calkine(mp, symb)
	ps = 0

	asquare = 0
	bsquare = 0
	if tp >= deltav or jumpword == 'lower':
		k = math.sqrt(1 - deltav / tp) if jumpword == 'upper' else math.sqrt(1 + deltav / tp)
		nmom2 = mom2 + (k - 1) * mp

		ep = totalv[-1] + tp
	
		vx = soc
		ex = (curre[0] + trane[0]) / 2

		asquare = ff / (64 * vx ** 3)
		bsquare = (ep - ex) / (2 * vx)

		deno = math.sqrt(bsquare ** 2 + 1) # if (np.sum(currf * tranf)) > 0 else math.sqrt(bsquare ** 2 - 1)
		detp = math.pi * math.sqrt(2 / (bsquare + math.sqrt(bsquare ** 2 + 1))) / math.sqrt(64 * asquare)
		ps = 1 - math.exp(-2 * detp) if dsquare * detp > 6 else 1 - (math.exp(2 * dsquare * detp - 2 * detp) - 1) / (math.exp(2 * dsquare * detp) - 1)

	pr = random.random()

	p = ps if ps > pr else 0.0

	nmom2 = tuple(map(tuple, nmom2))

	return p, nmom2, asquare, bsquare

