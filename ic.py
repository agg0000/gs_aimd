#!/usr/bin/env python
'''
for choose init or continue
'''

import re
import copy
import numpy as np
import collections

import hop
import softapi
import writefile
from constant import autoan

def initp(numele, totalmom, filename, outname, dt, fsdt, numcycle, quansoft, start, flip, istate, nstate, theory, mcenter):
	'''
	init function
	'''

	intcycle = 0
	spin = 'S'
	posa = collections.deque(maxlen = 4)
	moma = collections.deque(maxlen = 4)
	fara = collections.deque(maxlen = 4)
	totalv = collections.deque(maxlen = 4)
	deltav = collections.deque(maxlen = 4)
	transgap = collections.deque(maxlen = 4)

	numbermom = {}
	for mom in totalmom:
		everymom = re.split(r',|\s+', mom)

		while '' in everymom:
			everymom.remove('')

		symbolm = int(everymom[0])
		numbermom[symbolm] = everymom[1: 4]

	for i in range(numele):
		if i not in numbermom:
			numbermom[i] = [0, 0, 0]

	originmom = []
	for i in range(numele):
		originmom.append(numbermom[i])

	mom0 = np.array(originmom, dtype = float)

	state = istate
	pos0, symb = quansoft.getpos(filename, numele)
	mom0 = tuple(map(tuple, mom0))

	mpos = copy.deepcopy(pos0)
	mmom = copy.deepcopy(mom0)

	pos0 = hop.xcenter(pos0, symb, numele) if mcenter else pos0
	mom0 = hop.pcenter(mom0, symb, numele) if mcenter else mom0
	
	inp0 = quansoft.reinpotfile(filename, pos0, symb, intcycle, numele)
	oup0 = quansoft.runsoft(inp0, state, 'r')
	oldstate = state

	hopitem = 1

	far0, totalv, genergy = quansoft.parseout(oup0, numele, outname, state, totalv, start)
	kine = hop.calkine(mom0, symb)

	if nstate != 1:
		deltav, allpes = quansoft.gap(outname, oup0, inp0, state, state, nstate, deltav, start, flip, spin, theory, False , symb)
		if flip == 'r':
			_, tpes = quansoft.gap(outname, oup0, inp0, state, state, nstate, [], start, flip, spin, theory, True, symb)
			transgap.append(tpes)
	
	posa.append(pos0)
	fara.append(far0)
	moma.append(mom0)

	outfile = open(outname, 'w')
	outfile.write('*'*84 + '\n')
	outfile.write(' Welcome to Global switch program(v1.1.3beta) '.center(84, '*') + '\n')
	outfile.write('*'*84 + '\n'*2)

	outfile.write('input file name is <%s> \n' %filename)

	if fsdt:
		outfile.write('step(time) is %f fs \n' %fsdt)
	else:
		outfile.write('step(time) is %f a.u. \n' %dt)

	if flip == 'r':
		outfile.write('algorithm for spin-crossover\n')

	outfile.write('the total state is %d, the initial state is %d \n' %(nstate, istate))
	outfile.write('number of atoms need to calculate is %d \n\n' %numele)
	outfile.close()

	writefile.ncycle(outname, intcycle)
	writefile.outenergy(outname, totalv[-1], 'potential')
	writefile.outenergy(outname, kine, 'kinetic')
	writefile.outground(outname, genergy)

	if nstate != 1:
		writefile.spindr(outname, spin, state)
		writefile.wpes(outname, allpes[1:], genergy, spin)
		if flip == 'r':
			writefile.wpes(outname, tpes[1:], genergy, 'T')

	else:
		writefile.adiabatic(outname)

	if mcenter == 2:
		writefile.wdashed(outname)
		writefile.masscen(outname)
		writefile.writecon(outname, symb, mpos, 'before mcenter pos', numele)
		writefile.writecon(outname, symb, mpos, 'before mcenter mom', numele)
		writefile.wdashed(outname)

	writefile.writecon(outname, symb, pos0, 'coordinate', numele)
	writefile.writecon(outname, symb, far0, 'Forces', numele)
	writefile.writecon(outname, symb, mom0, 'momentum', numele)
	writefile.cyclend(outname, intcycle)

	return posa, fara, moma, totalv, deltav, symb, intcycle, state, spin, hopitem, transgap

#------------------------------------------------------------------------------------

def continuep(numele, outname, numcycle, nstate, start, mcenter):
	'''
	continue function
	'''

	newtxt = open(outname, 'a+')
	newtxt.seek(0)
	txt = newtxt.readlines()

	ecycle = 0
	for i, line in enumerate(txt):
		if ' cycle end' in line:
			ecycle = i

	if not ecycle:
		writefile.werror(out0, 'initial', start)
		print('error initial')
		exit(1)

	txt = txt[:ecycle + 1]
	cyclinearray = []
	hoplinearray = 0
	for i, line in enumerate(txt):
		if 'Cycle' in line:
			cyclinearray.append(i)

		if 'hop' in line:
			hoplinearray = i

	spin = 'S'
	posa = collections.deque(maxlen = 4)
	fara = collections.deque(maxlen = 4)
	moma = collections.deque(maxlen = 4)
	totalv = collections.deque(maxlen = 4)
	deltav = collections.deque(maxlen = 4)
	transgap = collections.deque(maxlen = 4)

	if hoplinearray < cyclinearray[-3]:
		ii = 3
		icycle = cyclinearray[-2]
		spes = []
		tpes = []
		while 'cycle end' not in txt[cyclinearray[-2] + ii]:
			line = txt[icycle + ii]

			if 'coordinate' in line:
				coornum = ii + icycle

			if 'Forces' in line:
				forcnum = ii + icycle

			if 'momentum' in line:
				momenum = ii + icycle

			if 'potential energy' in line:
				totalv.append(float(line.split()[-4]))

			if 'ground state' in line:
				genergy = float(line.split()[-4])

			if 'S' in line and 'pes' in line:
				senergy = float(line.split()[-2])
				spes.append(senergy)

			if 'T' in line and 'pes' in line:
				tenergy = float(line.split()[-2])
				tpes.append(tenergy)

			if 'spin' in line:
				sline = line.split()[-1]
				spin = re.search('\D', sline).group()
				state = int(re.search('\d', sline).group())

			ii += 1

		print(spes)
		dev = []
		if nstate != 1:
			if spin == 'S':
				spes1 = [genergy] + spes
				if state == 0:
					dev = [spes[0] - spes1[0]]
				elif state == nstate - 1:
					dev = [spes[state] - spes1[state]]
				else:
					dev = [spes[state - 1] - spes1[state - 1], spes[state] - spes1[state]]

				allpes = np.array([genergy] + tpes, dtype = float)
				transgap.append(allpes)

			else:
				tpes1 = tpes[1:]
				if tpes1:
					if state == 1:
						dev = [tpes1[0] - tpes[0]]
					elif state == nstate - 1:
						dev = [tpes1[state - 1] - tpes[state - 1]]
					else:
						dev = [tpes1[state - 2] - tpes[state - 2], tpes1[state - 1] - tpes[state - 1]]

					allpes = np.array([genergy] + spes, dtype = float)
					transgap.append(allpes)

			deltav.append(dev)

		pos0 = []
		symb = []
		for ii in range(numele):
			fs = txt[coornum + 4 + ii].split()
			xq = fs[1:]
			symb.append(fs[0])
			pos0.append(xq)
	
		pos0 = np.array(pos0, dtype = float) / autoan
		pos0 = tuple(map(tuple, pos0))
		posa.append(pos0)

		far0 = []
		for ii in range(numele):
			f = txt[forcnum + 4 + ii].split()[1:]
			far0.append(f)

		far0 = tuple(map(tuple, far0))
		fara.append(far0)

		if nstate != 1:
			if spin == 'S':
				pesa = [genergy] + spes
			else:
				pesa = [genergy] + tpes

			ev = pesa[state]
			totalv.append(ev)

		mom0 = []
		for ii in range(numele):
			p = txt[momenum + 4 + ii].split()[1:]
			mom0.append(p)

		mom0 = tuple(map(tuple, mom0))
		moma.append(mom0)

	icycle = cyclinearray[-1]
	intcycle = int(txt[icycle].split()[1])

	ii = 3
	hoptag = 0
	spes = []
	tpes = []
	while not txt[icycle + ii].count('cycle end') or txt[icycle + ii].count('ERROR'):
		line = txt[icycle + ii]

		if 'coordinate' in line:
			coornum = ii + icycle

		if 'Forces' in line:
			forcnum = ii + icycle

		if 'momentum' in line:
			momenum = ii + icycle

		if 'hop ' in line:
			hoptag = ii + icycle
			sline = line.split()[-1]
			spin = re.search('\D', sline).group()
			state = int(re.search('\d', sline).group())

		if 'newForces' in line:
			forcnum = ii + icycle

		if 'newmomentum' in line:
			momenum = ii + icycle

		if 'potential energy' in line:
			totalv.append(float(line.split()[-4]))

		if 'ground state' in line:
			genergy = float(line.split()[-4])

		if 'S' in line and 'pes' in line:
			senergy = float(line.split()[-2])
			spes.append(senergy)

		if 'T' in line and 'pes' in line:
			tenergy = float(line.split()[-2])
			tpes.append(tenergy)

		if 'spin' in line:
			sline = line.split()[-1]
			spin = re.search('\D', sline).group()
			state = int(re.search('\d', sline).group())

		ii += 1

	dev = []
	if nstate != 1:
		if spin == 'S':
			spes1 = [genergy] + spes
			if state == 0:
				dev = [spes[0] - spes1[0]]
			elif state == nstate - 1:
				dev = [spes[state] - spes1[state]]
			else:
				dev = [spes[state - 1] - spes1[state - 1], spes[state] - spes1[state]]

			allpes = np.array([genergy] + tpes, dtype = float)
			transgap.append(allpes)

		else:
			tpes1 = tpes[1:]
			if tpes1:
				if state == 1:
					dev = [tpes1[0] - tpes[0]]
				elif state == nstate - 1:
					dev = [tpes1[state - 1] - tpes[state - 1]]
				else:
					dev = [tpes1[state - 2] - tpes[state - 2], tpes1[state - 1] - tpes[state - 1]]
			
				allpes = np.array([genergy] + spes, dtype = float)
				transgap.append(allpes)

		deltav.append(dev)

	pos0 = []
	symb = []
	for ii in range(numele):
		fs = txt[coornum + 4 + ii].split()
		xq = fs[1:]
		symb.append(fs[0])
		pos0.append(xq)

	pos0 = np.array(pos0, dtype = float) / autoan
	pos0 = tuple(map(tuple, pos0))
	posa.append(pos0)

	far0 = []
	for ii in range(numele):
		f = txt[forcnum + 4 + ii].split()[1:]
		far0.append(f)

	far0 = tuple(map(tuple, far0))
	fara.append(far0)

	mom0 = []
	for ii in range(numele):
		p = txt[momenum + 4 + ii].split()[1:]
		mom0.append(p)

	mom0 = tuple(map(tuple, mom0))
	moma.append(mom0)

	cyclinearray.append(hoplinearray)
	cyclinearray.sort()

	hopitem = len(cyclinearray) - cyclinearray.index(hoplinearray) - 3

	newtxt.write('\n' + '#'*99 + '\n')
	newtxt.write(('restart at %s' %str(start)).center(99))
	newtxt.write('\n' + '#'*99 + '\n')
	newtxt.close()

	return posa, fara, moma, totalv, deltav, symb, intcycle, state, spin, hopitem, transgap
