#!/usr/bin/env python3
'''
======================================================================
gobal switch surface hopping program
base on Zhu-Nakamura theory
write by Zongzi(Linfeng Ye) E-mail:agg0000@163.com / 943109694@qq.com
======================================================================
as the main program to run the molecular dynamics
'''

import os
import re
import sys
import math
import datetime
import numpy as np

import ic
import hop
import softapi
import writefile
from constant import fstoau, evtoau, au2wavnum
from constant import molpro_np

start = datetime.datetime.now()
infile = sys.argv[1]
fname = infile.split('.')[0]

inp = open(infile).read().lower().splitlines()
fd = {}
for keyword in inp:
	ka = keyword.split('=')
	fd[ka[0]] = ka[1]

if 'fsdt' in fd:
	fsdt = float(fd['fsdt'])
	dt = fsdt / fstoau
elif 'dt' in fd:
	fsdt = 0
	dt = float(fd['dt'])
else:
	print('please input true unit of time')
	exit(1)

if 'socwm' in fd:
	socpara = float(fd['socwm']) / au2wavnum
	fd['flip'] = 'r'
elif 'socau' in fd:
	socpara = float(fd['socau'])
	fd['flip'] = 'r'
else:
	socpara = 0.0

soft_split = fd['softuse'].split() if 'softuse' in fd else 'gaussian'
softuse = soft_split[0]
if softuse == 'molpro' and len(soft_split) == 2:
    molpro_np.np = int(soft_split[1])

if softuse == 'molpro' and molpro_np.np == 1:
    molpro_np.np = int(fd['np']) if 'np' in fd else 1

numele = int(fd['numele'])
numcycle = int(fd['numcycle'])
filename = fd['filename'] if 'filename' in fd else fname
nstate = int(fd['nstate']) if 'nstate' in fd else 1
istate = int(fd['istate']) if 'istate' in fd else 0
mcenter = int(fd['mcenter']) if 'mcenter' in fd else 0
theory = fd['theory'] if 'theory' in fd else 'td'
keyname = fd['keyname'].lower() if 'keyname' in fd else 'i'
flip = fd['flip'].lower() if 'flip' in fd else 'd'

mdint = 0
if 'mdint' in fd:
    mdint = 1 if fd['mdint'] == 'beeman' else 0

totalmom = []
for word in fd:
	if word[0] == 'p':
		totalmom.append(fd[word])

quansoft = softapi.softname(softuse)
out0 = filename + ".oup"

if keyname == 'i':
	posa, fara, moma, totalv, deltav, symb, intcycle, state, spin, hopitem, transgap = ic.initp(numele, totalmom, filename, out0, dt, fsdt, numcycle, quansoft, start, flip, istate, nstate, theory, mcenter)
elif keyname == 'c':
	posa, fara, moma, totalv, deltav, symb, intcycle, state, spin, hopitem, transgap = ic.continuep(numele, out0, numcycle, nstate, start, mcenter)
else:
	writefile.werror(out0, 'keyname', start)
	print('no keyname for i or c')
	exit(1)

oldspin = spin
while intcycle < numcycle:
	intcycle += 1

	# Verlet algorithm for x
	#pos0 = hop.verletx(posa[-1], fara[-1], moma[-1], symb, numele, dt)
	if len(fara) > 1 and mdint: # Beeman algorithm for x
		pos0 = hop.beemanx(posa[-1], fara[-1], fara[-2], moma[-1], symb, numele, dt)
	else:
		pos0 = hop.verletx(posa[-1], fara[-1], moma[-1], symb, numele, dt)
	posa.append(pos0)

	newinp = quansoft.reinpotfile(filename, posa[-1], symb, intcycle, numele)
	oup0 = quansoft.runsoft(newinp, state, 'r')

	far0, totalv, genergy = quansoft.parseout(oup0, numele, out0, state, totalv, start)
	fara.append(far0)

	# Verlet algorithm for p
	#mom0 = hop.verletp(moma[-1], fara[-1], fara[-2], dt)
	if len(fara) > 2 and mdint: # Beeman algorithm for p
		mom0 = hop.beemanp(moma[-1], fara[-1], fara[-2], fara[-3], dt)
	else:
		mom0 = hop.verletp(moma[-1], fara[-1], fara[-2], dt)
	moma.append(mom0)

	kine = hop.calkine(moma[-1], symb)

	oldstate = state

	if nstate != 1:
		""" 'c': PESs cross
		    't': include triplet
		    'p': PESs parallel
		    's': only singlet
		    for example, 'tagcstate' means a include cross state tags"""
		probs = 0.0
		probtc = 0.0
		probtp = 0.0

		deltav, allpes = quansoft.gap(out0, oup0, newinp, state, state, nstate, deltav, start, flip, spin, theory, False, symb)
		hopitem += 1
		tagcstate = {} #tag for cross pes
		if flip == 'r' and '{}{}'.format(spin, state) != 'S0':
			_, tpes = quansoft.gap(out0, oup0, newinp, state, state, nstate, [], start, flip, spin, theory, True, symb)
			transgap.append(tpes)

			if spin == 'T' and hopitem >= 2:
				for i, ec in enumerate(tpes):
					if (totalv[-1] - ec) * (totalv[-2] - transgap[-2][i]) < 0:
						tagcstate[i] = totalv[-1] - ec

			elif spin == 'S' and hopitem >= 2:
				for i, ec in enumerate(tpes[1:]):
					if (totalv[-1] - ec) * (totalv[-2] - transgap[-2][i + 1]) < 0:
						tagcstate[i + 1] = totalv[-1] - ec

		if hopitem >= 3 or tagcstate:
			tagpstate = {} #tag for parallel pes
			if flip == 'r' and '{}{}'.format(spin, state) != 'S0' and not tagcstate:
				for i, ec in enumerate(tpes[1:]):
					ccoss = totalv[-1] - ec
					lcoss = totalv[-2] - transgap[-2][i + 1]
					scoss = totalv[-3] - transgap[-3][i + 1]
					if lcoss * scoss > 0 and np.all(np.abs([ccoss, lcoss, scoss]) < 0.25 * evtoau):
						tagpstate[i + 1] = [ccoss, lcoss, scoss]
					
			if tagcstate:
				soc = quansoft.getsoc(newinp, out0, state, nstate, spin, theory, intcycle, start) if not socpara else socpara
				probtca = {}
				for isoc in tagcstate:
					cursoc = soc[isoc] if not socpara else socpara
					dsquare = 1 + (2 * cursoc / tagcstate[isoc]) ** 2
					probtci, nmom2tc, asquaretc, bsquaretc = hop.tchopping(filename, out0, state, quansoft, dsquare, isoc, cursoc, intcycle, numele, totalv, moma[-1], symb, start, theory, spin)
					probtca[probtci] = isoc

				if probtca:
					probtc = max(probtca)
					tcinp = quansoft.getname(filename + str(intcycle))
					ncinp = quansoft.getname(filename + str(intcycle) + 'trans')

					tcoup = quansoft.runsoft(tcinp, state, 'g')
					ncoup = quansoft.runsoft(ncinp, state, 'g')

					newtcstate = probtca[probtc]
					oldspin = spin

			if tagpstate:
				socfile0 = quansoft.getname(filename + str(intcycle))
				socfile1 = quansoft.getname(filename + str(intcycle - 1))
				socfile2 = quansoft.getname(filename + str(intcycle - 2))

				soc0 = quansoft.getsoc(socfile0, out0, state, nstate, spin, theory, intcycle, start) if not socpara else socpara
				soc1 = quansoft.getsoc(socfile1, out0, state, nstate, spin, theory, intcycle, start) if not socpara else socpara
				soc2 = quansoft.getsoc(socfile2, out0, state, nstate, spin, theory, intcycle, start) if not socpara else socpara

				probtpa = {}
				for isoc in tagpstate:
					cursoc0 = soc0[isoc] if not socpara else socpara
					cursoc1 = soc1[isoc] if not socpara else socpara
					cursoc2 = soc2[isoc] if not socpara else socpara
					dsquare0 = 1 + (2 * cursoc0 / tagpstate[isoc][0]) ** 2
					dsquare1 = 1 + (2 * cursoc1 / tagpstate[isoc][1]) ** 2
					dsquare2 = 1 + (2 * cursoc2 / tagpstate[isoc][2]) ** 2
					if dsquare1 > dsquare0 and dsquare1 > dsquare2: # and dsquare1 > 1.15 ** 2:
						cursoc = soc1[isoc] if not socpara else socpara
						probtpi, nmom2tp, asquaretp, bsquaretp = hop.tphopping(filename, out0, state, quansoft, dsquare1, isoc, cursoc, intcycle, numele, totalv, moma[-2], symb, start, theory, spin)
						probtpa[probtpi] = isoc

				if probtpa:
					probtp = max(probtpa)
					tpinp = quansoft.getname(filename + str(intcycle - 1))
					npinp = quansoft.getname(filename + str(intcycle - 1) + 'trans')

					tpoup = quansoft.runsoft(tpinp, state, 'g')
					npoup = quansoft.runsoft(npinp, state, 'g')

					oldspin = spin
					newtpstate = probtpa[probtp]

			if state == max(range(nstate)) and '{}{}'.format(spin, state) != 'T1':
				if deltav[-2][0] < deltav[-3][0] and deltav[-2][0] < deltav[-1][0] and deltav[-2][0] < 0.2:
					probs, nmom2s, asquares, bsquares = hop.shopping(filename, out0, state, quansoft, deltav[-2][0], intcycle, numele, 'lower', totalv, posa, fara, moma, symb, start, theory, spin)

					if probs:
						newastate = state -1
						einp = quansoft.getname(filename + str(intcycle - 1))
						ninp = quansoft.getname(filename + str(intcycle - 1) + 'lower')

						eoup = quansoft.runsoft(einp, state, 'g')
						noup = quansoft.runsoft(ninp, state, 'g')

			elif state == 0 or '{}{}'.format(spin, state) == 'T1':
				if deltav[-2][0] < deltav[-3][0] and deltav[-2][0] < deltav[-1][0] and deltav[-2][0] < 0.2:
					probs, nmom2s, asquares, bsquares = hop.shopping(filename, out0, state, quansoft, deltav[-2][0], intcycle, numele, 'upper', totalv, posa, fara, moma, symb, start, theory, spin)

					if probs:
						newastate = state + 1
						einp = quansoft.getname(filename + str(intcycle - 1))
						ninp = quansoft.getname(filename + str(intcycle - 1) + 'upper')

						eoup = quansoft.runsoft(einp, state, 'g')
						noup = quansoft.runsoft(ninp, state, 'g')

			else:
				p1 = 0.0
				p2 = 0.0

				a1 = 0.0
				b1 = 0.0
				a2 = 0.0
				b2 = 0.0

				if deltav[-2][0] < deltav[-3][0] and deltav[-2][0] < deltav[-1][0] and deltav[-2][0] < 0.2:
					p1, nmoml, a1, b1 = hop.shopping(filename, out0, state, quansoft, deltav[-2][0], intcycle, numele, 'lower', totalv, posa, fara, moma, symb, start, theory, spin)
				if deltav[-2][1] < deltav[-3][1] and deltav[-2][1] < deltav[-1][1] and deltav[-2][1] < 0.2:
					p2, nmomu, a2, b2 = hop.shopping(filename, out0, state, quansoft, deltav[-2][1], intcycle, numele, 'upper', totalv, posa, fara, moma, symb, start, theory, spin)

				probs = p1 if p1 > p2 else p2
				if probs:
					einp = quansoft.getname(filename + str(intcycle - 1))
					eoup = quansoft.runsoft(einp, state, 'g')

					if probs == p1:
						newastate = state - 1
						ninp = quansoft.getname(filename + str(intcycle - 1) + 'lower')
						noup = quansoft.runsoft(ninp, state, 'g')

						nmom2s = nmoml
						asquares = a1
						bsquares = b1

					elif probs == p2:
						newastate = state + 1
						ninp = quansoft.getname(filename + str(intcycle - 1) + 'upper')
						noup = quansoft.runsoft(ninp, state, 'g')

						nmom2s = nmomu
						asquares = a2
						bsquares = b2

			prob = max([probs, probtc, probtp])
			if prob:
				oldstate == state
				if prob == probtc:
					newspin = 'S' if oldspin == 'T' else 'T'
					writefile.ncycle(out0, intcycle)
					writefile.outenergy(out0, totalv[-1], 'potential')
					writefile.outenergy(out0, kine, 'kinetic')

					writefile.spindr(out0, spin, state)
					writefile.wpes(out0, allpes[1:], genergy, oldspin)
					writefile.wpes(out0, tpes[1:], genergy, newspin)

					writefile.writecon(out0, symb, posa[-1], 'coordinate', numele)
					writefile.writecon(out0, symb, fara[-1], 'Forces', numele)
					writefile.writecon(out0, symb, moma[-1], 'momentum', numele)

					asquare, bsquare, state, spin = asquaretp, bsquaretp, newtcstate, newspin
					newinp, oldinp, oup = ncinp, tcinp, ncoup
					hopword = 'cross'

					moma[-1] = nmom2tc

				else:
					intcycle -= 1
					writefile.trunline(out0)

					del posa[-1]
					del fara[-1]
					del moma[-1]
					del totalv[-1]

					if prob == probtp:
						spin = 'S' if oldspin == 'T' else 'T'
						asquare, bsquare, state = asquaretp, bsquaretp, newtpstate
						newinp, oldinp, oup = npinp, tpinp, npoup
						moma[-1] = nmom2tp
						hopword = 'parallel'
					else:
						if theory == "td":
							kine2 = hop.calkine(moma[-2], symb)
							kine1 = hop.calkine(moma[-1], symb)
							kines = hop.calkine(nmom2s, symb)

							tot2 = kine2 + totalv[-2]
							tot1 = kine1 + totalv[-1]

							alpha = math.sqrt(1 + (tot2 - tot1) / kines)
							nmom2s = alpha * np.array(nmom2s)
							nmom2s = tuple(map(tuple, nmom2s))

						asquare, bsquare, state = asquares, bsquares, newastate
						newinp, oldinp, oup = ninp, einp, noup
						moma[-1] = nmom2s
						hopword = 'avoid'

				hopitem = -3
				farn, totalv, genergy = quansoft.parseout(oup, numele, out0, state, totalv, start)

				fara[-1] = farn
				newkine = hop.calkine(moma[-1], symb)

				newinp, oldinp = quansoft.getrename(newinp, oldinp)

				os.rename(newinp, oldinp)
				writefile.whop(out0, intcycle, state, spin, oldstate, oldspin)
				writefile.outenergy(out0, totalv[-1], 'potential')
				writefile.outenergy(out0, newkine, 'kinetic')
				writefile.writeab(out0, asquare, bsquare, hopword)

				writefile.writecon(out0, symb, fara[-1], 'newForces', numele)
				writefile.writecon(out0, symb, moma[-1], 'newmomentum', numele)
				writefile.cyclend(out0, intcycle)

				continue

	writefile.ncycle(out0, intcycle)
	writefile.outenergy(out0, totalv[-1], 'potential')
	writefile.outenergy(out0, kine, 'kinetic')

	if nstate != 1:
		ge = allpes[0] if spin == 'S' else tpes[0]
		writefile.outground(out0, ge)
		writefile.spindr(out0, spin, state)
		if spin == 'S':
			writefile.wpes(out0, allpes[1:], genergy, 'S')
			if flip == 'r':
				writefile.wpes(out0, tpes[1:], genergy, 'T')
		else:
			writefile.wpes(out0, tpes[1:], genergy, 'S')
			writefile.wpes(out0, allpes[1:], genergy, 'T')

	else:
		writefile.outground(out0, genergy)
		writefile.adiabatic(out0)

	writefile.writecon(out0, symb, posa[-1], 'coordinate', numele)
	writefile.writecon(out0, symb, fara[-1], 'Forces', numele)

	writefile.writecon(out0, symb, moma[-1], 'momentum', numele)
	writefile.cyclend(out0, intcycle)

	if nstate != 1:
		boola = np.all([num > 0 for num in deltav[-1]])
		
		if not boola:
			writefile.werror(out0, 'negative', start)
			print('negative energy')
			exit(1)

writefile.endtime(out0, 'cycle', intcycle, start)
print('cycle normal exit')
exit(0)
