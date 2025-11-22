#!/usr/bin/env python
'''
set up the API for use quantum chemistry self
'''

import os
import re
import shutil
import datetime
import numpy as np

import writefile
from constant import evtoau, autoan, au2wavnum, ldaspin, pbespin
from constant import molpro_np

def softname(softuse):
	'''
	convince the quantum chemistry from the keyword
	'''
	softall = {"gaussian" : gaussian(),  "molpro" : molpro(), "dftbplus" : dftbp(),
                   "turbmole" : turbmole(),  "bdf"    : BDF(),    "oniom"  : oniom() }

	if softuse not in softall:
		print("please use useful key words")
		exit(1)

	return softall.get(softuse)

#####################################################################################
#                            class gaussian                                         #
#####################################################################################

class gaussian():
	'''
	use the gaussian soft
	'''
	def __init__(self):
		self = self

#-------------------------------------------------------------------------------------

	def getname(self, filename):
		outname = filename + '.gjf'

		return outname

#-------------------------------------------------------------------------------------

	def getrename(self, iname, oname):

		return iname, oname

#-------------------------------------------------------------------------------------

	def sample(self, filename, pos0, mom0, freq):
		runfile = filename.split('.')[0] + '.gjf'
		newinpf = filename.split('.')[0] + '0.gjf'

		oldfile = open(runfile).readlines()
		newfile = open(newinpf, 'w')

		for line in oldfile:
			newfile.write(line.lower())

		newfile.close()

		return newinpf, pos0, mom0

#-------------------------------------------------------------------------------------

	def reinpotfile(self, filename, newpos, symb, intcycle, numele):
		'''
		use the new position to create a new gjf file
		'''
		if intcycle:
			oldname = self.getname(filename + str(intcycle - 1))
			newname = self.getname(filename + str(intcycle))
		else:
			oldname = self.getname(filename)
			newname = self.getname(filename + str(intcycle))

		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		newpos = np.array(newpos, dtype = float) * autoan
		linenu = 0 # ensure the positine line number.
		linena = 0
		for i, line in enumerate(oldfile):
			if line.strip() == "":
				linena += 1

			if linena == 2:
				linenu = i + 2
				break

		for i in range(linenu): 
			newfile.write(oldfile[i])

		for i in range(linenu, linenu + numele): 
			newfile.write(' {0:<5}{1[0]:15.8f}{1[1]:15.8f}{1[2]:15.8f}\n'.format(symb[i - linenu], newpos[i - linenu]))

		for i in range(linenu + numele, len(oldfile)):
			newfile.write(oldfile[i])

		newfile.close()
		return newname

#-------------------------------------------------------------------------------------

	def getpos(self, filename, numele):
		'''
		from the first gjf file to get the first position
		'''
		runfile = filename + ".gjf"
		oldfile = open(runfile).readlines()
		filelen = len(oldfile)

		linenu = 0 # ensure the positine line number.
		linena = 0
		for i, line in enumerate(oldfile):
			if line.strip() == "":
				linena += 1

			if linena == 2:
				linenu = i + 2
				break

		numpos  = oldfile[linenu: linenu + numele]

		pos = []
		sym = []
		for i in range(numele):
			posarray = numpos[i].split()
			symbol = posarray[0].title()

			pos.append(posarray[1:])
			sym.append(symbol)

		pos = np.array(pos, dtype = float) / autoan
		sym = tuple(sym)
		pos = tuple(map(tuple, pos))

		return pos, sym
	
#-------------------------------------------------------------------------------------

	def jpfile(self, oldname, outname, symb, newstate, state, jumpword, theory, spin, trans, starttime):
		newname = self.getname(oldname.split(".")[0] + jumpword)
		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		if trans:
			oldstr = True if spin == 'S' else False
			newstr = not oldstr if trans else oldstr
	
			oldspin = 'singlet' if oldstr else 'triplet'
			newspin = 'singlet' if newstr else 'triplet'

		for line in oldfile:
			if line[0] == '#':
				if state == 0:
					line = line.strip() + ' {}(singlet,root=1,nstate=8)\n'.format(theory)

				elif newstate == 0:
					line = re.sub(' {}\(.+\)'.format(theory), '', line)
				
				elif trans:
					line = line.replace(oldspin, newspin)
					line = line.replace('root=%d' %state, 'root=%d' %newstate)

				elif state != newstate:
					line = line.replace('root=%d' %state, 'root=%d' %newstate)

			newfile.write(line)

		newfile.close()

		return newname

#-------------------------------------------------------------------------------------

	def gap(self, outname, oup, inp, newstate, state, nstate, deltav, start, flip, spin, theory, trans, symb):
		efile = oup
		if state == 0:
			inpe = self.jpfile(inp, outname, symb, 1, 0, 'upper', theory, spin, False, start)
			oupe = self.runsoft(inpe, state, 'r')
			efile = oupe

		if trans:
			inpe = self.jpfile(inp, outname, symb, newstate, state, 'trans', theory, spin, True, start)
			oupe = self.runsoft(inpe, state, 'r')
			efile = oupe

		enefile = open(efile).readlines()

		genergy = 0.0
		excitarray = [0.0]
		for line in enefile:
			if 'SCF Done' in line:
				genergy = float(line.split()[4])

			if 'Excited State' in line:
				gapenergy = line.split()[4]
				excitarray.append(gapenergy)

		excitarray = np.array(excitarray, dtype = float)
		allpes = excitarray[:nstate] * evtoau + genergy

		dev = []
		for i in range(1, nstate):
			dev.append(excitarray[i] - excitarray[i - 1])

		if not dev:
			writefile.werror(outname, 'deltav', start)
			print('no gap energy in outfile')
			exit(1)

		if not trans:
			if spin == 'S' and state == 0:
				deltav.append([dev[0]])
	
			elif spin == 'T' and state == 1:
				deltav.append([dev[1]])

			elif state == max(range(nstate)):
				deltav.append([dev[-1]])
	
			else:
				deltav.append(dev[state - 1: state + 1])

		return deltav, allpes

#-------------------------------------------------------------------------------------

	def runsoft(self, runfile, state, keysoft):
		'''
		for run the gaussian progam, and return the log file
		'''
		outname = runfile.split('.')[0] + '.log'

		if keysoft == 'r':
			os.system('g16 %s' %runfile)

		return outname

#------------------------------------------------------------------------------------

	def parseout(self, resultname, numele, outname, state, potential, starttime):
		resultfile = open(resultname).readlines()

		signf = 0
		signe = []
		for i, line in enumerate(resultfile):
			if 'Forces ' in line:
				signf = i

			if 'SCF Done' in line:
				genergy = float(line.split()[4])
				signe.append(line)

			if 'Total Energy' in line:
				signe.append(line)

		if 'Forces ' in resultfile[signf]:
			numf = resultfile[signf + 3: signf + 3 + numele]
		else:
			writefile.werror(outname, 'forces', starttime)
			print('no force in outfile')
			exit(1)

		sfar = []
		for i in range(numele):
			eachforce = numf[i].split()
			sfar.append(eachforce[2:])

		sfar = np.array(sfar, dtype = float)
		far = tuple(map(tuple, sfar))
		
		ii = signe[-1]
		if 'SCF Done' in ii:
			energyrow = ii.split()[4]

		elif 'Total Energy' in ii:
			energyrow = ii.split()[-1]

		else:
			writefile.werror(outname, 'energy', starttime)
			print('no energy in outfile')
			exit(1)

		potential.append(float(energyrow))

		return far, potential, genergy

#-------------------------------------------------------------------------------------

	def getsoc(self, filename, outname, state, nstate, spin, theory, intcycle, starttime):
		with open(filename) as openf:
			openl = openf.readlines()

		with open('gaussian.gjf', 'w') as socf:
			for line in openl:
				line = '%rwf=gaussian.rwf\n' if 'chk' in line else line
				if theory in line:
					line = re.sub('{}\(.+\)'.format(theory), '{}(50-50,nstate=8)'.format(theory), line)
					if 'gfinput' not in line.lower():
						line = line.strip() + ' 6d 10f gfinput\n'

				socf.write(line)
				
		soclog = self.runsoft('gaussian.gjf', state, 'r')

		with open('init.py', 'w') as inif:
			inif.write("import os\n")
			inif.write("import sys\n\n")
			inif.write("QM_ex_flag=False\n")
			inif.write("QM_code='gauss_tddft'\n")
			inif.write("n_s=range(1, {})\n".format(nstate + 1))
			inif.write("n_t=range(1, {})\n".format(nstate + 1))
			inif.write("n_g=['True']\n")
			inif.write("soc_scal=1.0\n")
			inif.write("cicoeff_thresh=[1.0-5]\n\n")
			
			inif.write("g09root='{}'\n".format(os.environ['g16root']))
			inif.write("sys.path.append(g09root+'/g16')\n")
			inif.write("molsoc_path='molsoc0.1.exe'\n\n")
			inif.write("qm_out=['gaussian.log','gaussian.rwf']\n")
			inif.write("soc_key=['ANG', 'Zeff', 'DIP']\n")
			inif.write("molsoc_input=['molsoc.inp', 'molsoc_basis']\n\n")

		os.system('pysoc.py')

		if os.path.exists('soc_out.dat'):
			writefile.werror(outname, 'soc', starttime)
			print('something error in calculate soc')
			exit(1)

		with open('soc_out.dat') as sof:
			sol = sof.readlines()

		soca = []
		for line in sol:
			if '{}{}'.format(oldspin, state) in line:
				soc = float(line.split()[3]) / au2wavnum
				soca.append(soc)

		return soca

#####################################################################################
#                             class molpro                                          #
#####################################################################################

class molpro():
	'''
	use the molpro soft
	'''
	def __init__(self):
		self = self

#-------------------------------------------------------------------------------------

	def getname(self, filename):
		outname = filename + '.in'

		return outname

#-------------------------------------------------------------------------------------

	def getrename(self, iname, oname):

		return iname, oname

#-------------------------------------------------------------------------------------

	def sample(self, filename, pos0, mom0, freq):
		fname = filename.split('.')[0]
		nfname = fname + str(0)

		runfile = fname + '.in'
		newinpf = nfname + '.in'
		#os.system('molpro %s' %runfile)

		oldfile = open(runfile).readlines()
		newfile = open(newinpf, 'w')

		for line in oldfile:
			if "file,2" in line:
				line0 = line.replace('wfu', 'wfu,old')
				line1 = line.replace(fname, nfname)
				line = line0 + line1
				
			newfile.write(line)
	
		newfile.close()
		
		return newinpf, pos0, mom0

#-------------------------------------------------------------------------------------

	def reinpotfile(self, filename, newpos, symb, intcycle, numele):
		'''
		use the new position to create a new in file
		'''
		if intcycle:
			oldname = self.getname(filename + str(intcycle - 1))
			newname = self.getname(filename + str(intcycle))
		else:
			oldname = self.getname(filename)
			newname = self.getname(filename + str(intcycle))

		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		newpos = np.array(newpos, dtype = float) * autoan
		
		wline = 99999
		for i, line in enumerate(oldfile):
			if 'geometry' in line:
				geoline = i + 3

			if 'file,2' in line and 'old' not in line:
				wline = i

		for i in range(geoline):
			if i == wline:
				oldfile[i] = oldfile[i].replace(oldname, newname)

			newfile.write(oldfile[i])

		for i in range(geoline, geoline + numele):
			newfile.write(' {0:<5}{1[0]:15.8f}{1[1]:15.8f}{1[2]:15.8f}\n'.format(symb[i - geoline], newpos[i - geoline]))			

		for i in range(geoline + numele, len(oldfile)):
			newfile.write(oldfile[i])

		newfile.close()
		return newname

#-------------------------------------------------------------------------------------

	def getfre(self, outfile):
		'''
		get the frequence
		'''
		freq = 0
	
		return freq

#-------------------------------------------------------------------------------------

	def getpos(self, filename, numele):
		'''
		from the first input file to get the first position
		'''
		runfile = filename + ".in"
		oldfile = open(runfile).readlines()

		for i, line in enumerate(oldfile):
			if 'geometry' in line:
				geoline = i + 3

		numpos = oldfile[geoline: geoline + numele]

		pos = []
		sym = []
		for i in range(numele):
			posarray = numpos[i].split()
			symbol = posarray[0].title()

			pos.append(posarray[1:])
			sym.append(symbol)
		
		pos = np.array(pos, dtype = float) / autoan
		sym = tuple(sym)
		pos = tuple(map(tuple, pos))

		return pos, sym

#-------------------------------------------------------------------------------------

	def jpfile(self, oldname, outname, symb, newstate, state, jumpword, theory, spin, trans, starttime):
		newname = oldname[:-3] + jumpword + oldname[-3:]
		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		if trans:
			oldstr = True if spin == 'S' else False
			newstr = not oldstr if trans else oldstr
	
			oldspin = 0 if oldstr else 2
			newspin = 0 if newstr else 2

		mi = 0
		for i, l in enumerate(oldfile):
			l = l.lower()
			if 'multi' in l:
				mi = i

		for ii, line in enumerate(oldfile):
			line = line.lower()
			if 'wf,' in line and trans and ii > mi:
				wfline = re.compile('wf,\d+,\d+')
				ne = wfline.search(line).group().split(',')[1]
				line = line.replace('wf,{},1,{}'.format(ne, oldspin), 'wf,{},1,{}'.format(ne, newspin))

			if 'cpmcscf' in line:
				line = line.replace('grad,%d.1' %(state + 1), 'grad,%d.1' %(newstate + 1))

			if 'root' in line:
				line = line.replace('root=%d' %(state + 1), 'root=%d'%(newstate + 1))

			newfile.write(line)

		newfile.close()

		return newname

#-------------------------------------------------------------------------------------

	def gap(self, outname, oup, inp, newstate, state, nstate, deltav, start, flip, spin, theory, trans, symb):
		efile = oup
		if trans:
			inpe = self.jpfile(inp, outname, symb, newstate, state, 'trans', theory, spin, True, start)
			oupe = self.runsoft(inpe, state, 'r')
			efile = oupe

		enefile = open(efile).readlines()

		excity1 = []
		excity2 = []
		for i, line in enumerate(enefile):
			if 'MCSCF STATE' in line and 'Energy' in line :
				excity1.append(line.split()[-1])

			if 'RSPT2 STATE' in line and 'Energy' in line:
				excity2.append(line.split()[-1])

		excitarray = excity2 if excity2 else excity1

		if spin == 'S' and trans or spin == 'T' and not trans:
			excitarray.pop()
			excitarray.insert(0, excitarray[0])

		allpes = np.array(excitarray, dtype = float)
		excitarray = np.array(excitarray, dtype = float) / evtoau

		dev = []
		for i in range(1, nstate):
			dev.append(excitarray[i] - excitarray[i - 1])

		if not dev:
			writefile.werror(outname, 'deltav', start)
			print('no gap energy in outfile')
			exit(1)

		if not trans:
			if spin == 'S' and state == 0:
				deltav.append([dev[0]])
	
			elif spin == 'T' and state == 1:
				deltav.append([dev[1]])

			elif state == max(range(nstate)):
				deltav.append([dev[-1]])
	
			else:
				deltav.append(dev[state - 1: state + 1])

		return deltav, allpes

#-------------------------------------------------------------------------------------

	def runsoft(self, runfile, state, keysoft):
		'''
		for run the molpro progam, and return the output file
		'''
		outname = runfile.split('.')[0] + '.out'

		os.system('molpro -n %d -t 1 %s' %(molpro_np.np, runfile))

		return outname

#------------------------------------------------------------------------------------

	def parseout(self, resultname, numele, outname, state, potential, starttime):
		resultfile = open(resultname).readlines()

		signf = 0
		genergy = 0.0
		ev = 0.0
		for i, line in enumerate(resultfile):
			if 'GRADIENT FOR ' in line:
				signf = i

			if 'STATE 1.1 Energy' in line:
				genergy = float(line.split()[-1])

			if 'STATE {:d}.1 Energy'.format(state + 1) in line:
				ev = line.split()[-1]

		if 'GRADIENT FOR ' in resultfile[signf]:
			numf = resultfile[signf + 4: signf + 4 + numele]
		else:
			writefile.werror(outname, 'forces', starttime)
			print('no force in outfile')
			exit(1)

		sfar = []
		for i in range(numele):
			eachforce = numf[i].split()
			sfar.append(eachforce[1:])

		far = np.array(sfar, dtype = float) * -1.0
		far = tuple(map(tuple, far))
		
		if ev:
			outenergy = float(ev)		
			potential.append(outenergy)
		else:
			writefile.werror(outname, 'energy', starttime)
			print('no energy in outfile')
			exit(1)

		return far, potential, genergy

#------------------------------------------------------------------------------------

	def getsoc(self, filename, outname, state, nstate, spin, theory, intcycle, starttime):
		filename = filename.split('.')[0]
		with open('{}.in'.format(filename)) as openf:
			openl = openf.readlines()

		wf = 0
		occ = 0
		closed = 0
		for line in openl:
			if 'wf,' in line:
				wfline = re.compile('wf,\d+,\d+')
				wf = wfline.search(line).group().split(',')[1]

			if 'occ' in line:
				occline = re.compile('occ,\d+')
				occ = occline.search(line).group().split(',')[1]

			if 'closed' in line:
				closedline = re.compile('closed,\d+')
				closed = closedline.search(line).group().split(',')[1]

			if wf and occ and closed:
				break

		if not wf and not occ and not closed:
			writefile.werror(outname, 'soc keyword', starttime)
			print('something error in soc input file')
			exit(1)

		if not os.path.exists('{}_soc.in'.format(filename)):
			with open('{}_soc.in'.format(filename), 'w') as socf:
				for line in openl:
					if 'file,2' in line:
						continue

					socf.write(line)

					if '}' in line:
						break

				socf.write('\nint\n{{hf;wf,{},1,0;}}\n\n'.format(wf))
				socf.write('{{multi;occ,{};closed,{};\n'.format(occ, closed))
				socf.write('wf,{},1,0;state,{};\n'.format(wf, nstate))
				socf.write('wf,{},1,2;state,{};}}\n\n'.format(wf, nstate))
				socf.write('{{ci;occ,{0};closed,{1};core,{1};wf,{2},1,0;state,{3};save,3040.1;noexc}}\n'.format(occ, closed, wf, nstate))
				socf.write('{{ci;occ,{0};closed,{1};core,{1};wf,{2},1,2;state,{3};save,3041.1;noexc}}\n\n'.format(occ, closed, wf, nstate))
				socf.write('lsint\n\n')
				socf.write('{ci;hlsmat,ls,3040.1,3041.1;print,hls=2,vls=0}\n\n')
				socf.write('---\n\n\n\n')

		keyw = 'g' if os.path.exists('{}_soc.out'.format(filename)) else 'r'
		socfile = self.runsoft('{}_soc.in'.format(filename), state, keyw)

		with open(socfile) as soctf:
			soctl = soctf.readlines()

		soci = 0
		bpob = 0
		endp = 0
		totalsoc = []
		for i, line in enumerate(soctl):
			bra = 'Bra-wavefunction' in line and '3040.1' in line
			if bra:
				ket = 'Ket-wavefunction' in soctl[i + 1] and '3041.1' in soctl[i + 1]
				soci = 1 if ket else 0

			if soci:
				if 'Breit-Pauli' in line:
					bpob = i + 2

				if bpob and line == '\n':
					endp += 1

				if endp == 2:
					totalsoc.append(soctl[bpob: i])
					break

		if not totalsoc:
			writefile.werror(outname, 'soc', starttime)
			print('something error in soc output file')
			exit(1)

		smark = str(state + 0.1) + '>' if spin == 'T' else '<' + str(state + 1.1)
		socarray = []
		for dsoc in totalsoc:
			socd = []
			for line in dsoc:
				if smark in line:
					sd = line.split()[3].replace('i', '')
					socd.append(float(sd))

			socarray.append(socd)

		socnew = np.array(socarray).T
		soca = np.zeros(socnew.shape[0])
		for i, num in enumerate(socnew):
			soca[i] = np.linalg.norm(num)

		return soca

#####################################################################################
#                              class turbmole                                       #
#####################################################################################

class turbmole():
	'''
	use the turbomole soft
	'''
	def __init__(self):
		self = self

#-------------------------------------------------------------------------------------

	def getname(self, filename):
		outname = filename

		return outname

#-------------------------------------------------------------------------------------

	def getrename(self, iname, oname):
		iname += '/control'
		oname += '/control'

		return iname, oname

#-------------------------------------------------------------------------------------

	def sample(self, filename, pos0, mom0, freq):
		runfile = filename
		newinpf = filename + '0'

		if os.path.isfile(newinpf):
			os.remove(newinpf)

		if not os.path.isdir(newinpf):
			os.mkdir(newinpf)

		if os.path.exists('{}/auxbasis'.format(runfile)):
			shutil.copy('{}/auxbasis'.format(runfile), newinpf)

		if os.path.exists('{}/control'.format(runfile)):
			shutil.copy('{}/control'.format(runfile), newinpf)

		if os.path.exists('{}/coord'.format(runfile)):
			shutil.copy('{}/coord'.format(runfile), newinpf)

		if os.path.exists('{}/basis'.format(runfile)):
			shutil.copy('{}/basis'.format(runfile), newinpf)

		if os.path.exists('{}/mos'.format(runfile)):
			shutil.copy('{}/mos'.format(runfile), newinpf)

		if os.path.exists('{}/alpha'.format(runfile)):
			shutil.copy('{}/alpha'.format(runfile), newinpf)

		if os.path.exists('{}/beta'.format(runfile)):
			shutil.copy('{}/beta'.format(runfile), newinpf)

		return newinpf, pos0, mom0

#-------------------------------------------------------------------------------------

	def reinpotfile(self, filename, newpos, symb, intcycle, numele):
		'''
		use the new position to create a new in dir
		'''
		oldname = filename + str(intcycle - 1)
		newname = filename + str(intcycle)

		if not os.path.isdir(newname):
			os.mkdir(newname)

		oldfile = open('{}/coord'.format(oldname)).readlines()
		newfile = open('{}/coord'.format(newname), 'w')

		newfile.write(oldfile[0])
		for i in range(1, numele + 1):
			newfile.write('{1[0]:15.8f}{1[1]:15.8f}{1[2]:15.8f}{0:>5}\n'.format(symb[i - 1], newpos[i - 1]))			
		
		for i in range(numele + 1, len(oldfile)):
			newfile.write(oldfile[i])

		newfile.close()

		if os.path.exists('{}/auxbasis'.format(oldname)):
			shutil.copy('{}/auxbasis'.format(oldname), newname)

		if os.path.exists('{}/control'.format(oldname)):
			shutil.copy('{}/control'.format(oldname), newname)

		if os.path.exists('{}/basis'.format(oldname)):
			shutil.copy('{}/basis'.format(oldname), newname)

		if os.path.exists('{}/mos'.format(oldname)):
			shutil.copy('{}/mos'.format(oldname), newname)

		if os.path.exists('{}/alpha'.format(runfile)):
			shutil.copy('{}/alpha'.format(runfile), newinpf)

		if os.path.exists('{}/beta'.format(runfile)):
			shutil.copy('{}/beta'.format(runfile), newinpf)

		return newname

#-------------------------------------------------------------------------------------

	def getfre(self, outfile):
		'''
		get the frequence
		'''
		freq = 0
	
		return freq

#-------------------------------------------------------------------------------------

	def getpos(self, filename, numele):
		'''
		from the first input coord file to get the first position
		'''
		oldfile = open('{}/coord'.format(filename)).readlines()

		pos = []
		sym = []
		for i in range(1, numele + 1):
			posarray = oldfile[i].split()
			pos.append(posarray[:3])
			symbol = posarray[-1].title()

			sym.append(symbol)
			
		pos = np.array(pos, dtype = float)
		sym = tuple(sym)
		pos = tuple(map(tuple, pos))

		return pos, sym

#-------------------------------------------------------------------------------------

	def jpfile(self, oldname, outname, symb, newstate, state, jumpword, theory, spin, trans, starttime):
		newname = oldname + jumpword

		if not os.path.isdir(newname):
			os.mkdir(newname)

		oldfile = open('{}/control'.format(oldname)).readlines()
		newfile = open('{}/control'.format(newname), 'w')

		if trans:
			oldstr = True if spin == 'S' else False
			newstr = not oldstr if trans else oldstr
	
			oldspin = 's' if oldstr else 't'
			newspin = 's' if newstr else 't'

		for line in oldfile:
			if 'exopt' in line:
				line = line.replace(str(state), str(newstate))

			if trans:
				if 'scfinstab' in line:
					line = line.replace('{}{}'.format(theory, oldspin), '{}{}'.format(theory, newspin))

			newfile.write(line)

		newfile.close()

		if os.path.exists('{}/auxbasis'.format(oldname)):
			shutil.copy('{}/auxbasis'.format(oldname), newname)

		if os.path.exists('{}/coord'.format(oldname)):
			shutil.copy('{}/coord'.format(oldname), newname)

		if os.path.exists('{}/basis'.format(oldname)):
			shutil.copy('{}/basis'.format(oldname), newname)

		if os.path.exists('{}/mos'.format(oldname)):
			shutil.copy('{}/mos'.format(oldname), newname)

		if os.path.exists('{}/alpha'.format(runfile)):
			shutil.copy('{}/alpha'.format(runfile), newinpf)

		if os.path.exists('{}/beta'.format(runfile)):
			shutil.copy('{}/beta'.format(runfile), newinpf)

		return newname

#-------------------------------------------------------------------------------------

	def gap(self, outname, oup, inp, newstate, state, nstate, deltav, start, flip, spin, theory, trans, symb):
		if trans:
			inpe = self.jpfile(inp, outname, symb, newstate, state, 'trans', theory, spin, True, start)
			oup = self.runsoft(inpe, state, 'r')

		efile = oup + '/out'
		enef = open(efile)
		enefile = enef.readlines()
		enef.close()

		excitline = []
		for i, line in enumerate(enefile):
			if 'a excitation' in line:
				excitline.append(i)

		if not excitline:
			os.chdir(oup)
			os.system('rm *a')

			sizestr = os.popen('grep gridsize control').read().split()[1]
			sizenum = int(re.sub(r'\D', '', sizestr))

			if sizenum < 5:
				os.system("sed -i 's/m{0}/m{1}/' control".format(sizenum, sizenum + 1))

			if state == 0:
				os.system('escf >> out')
			else:
				os.system('egrad >> out')
				
			os.chdir('..')

			enef = open(efile)
			enefile = enef.readlines()
			enef.close()

		for i, line in enumerate(enefile):
			if 'a excitation' in line:
				excitline.append(i)

		if not excitline:
			writefile.werror(outname, 'deltav', start)
			print('no gap energy in outfile')
			exit(1)

		dev = []
		excitarray = [0.0]
		for ii in excitline:
			if 'Excitation' in enefile[ii + 5]:
				gapenergy = float(enefile[ii + 7].split()[-1])
				excitarray.append(gapenergy)

		allpes = np.array(excitarray[:nstate], dtype = float) * evtoau
		for i in range(1, nstate):
			dev.append(excitarray[i] - excitarray[i - 1])

		if state == 0:
			deltav.append([dev[0]])

		elif state == max(range(nstate)):
			deltav.append([dev[-1]])

		else:
			deltav.append(dev[state - 1: state + 1])

		return deltav, allpes

#-------------------------------------------------------------------------------------

	def runsoft(self, runfile, state, keysoft):
		'''
		for run the turbomole progam, and return the output file
		'''
		outname = runfile

		if keysoft == 'r':
			os.chdir(runfile)
			os.system('ridft >> out')

			if state == 0:
				os.system('rdgrad >> out')
				os.system('escf >> out')

			else:
				os.system('egrad >> out')

			os.chdir('..')

		return outname

#------------------------------------------------------------------------------------

	def parseout(self, resultname, numele, outname, state, potential, starttime):
		gradname = resultname + '/gradient'

		if not os.path.exists(gradname) and state != 0:
			os.chdir(resultname)
			os.system('rm *a')
			if state != 0:
				os.system('egrad >> out')

			os.chdir('..')

		if not os.path.exists(gradname):
			writefile.werror(outname, 'forces', starttime)
			print('no force in outfile')
			exit(1)

		resultfile = open(gradname).readlines()

		numf = resultfile[numele + 2: numele + 2 + numele]

		sfar = []
		for i in range(numele):
			nuf = numf[i].replace('D', 'e')
			sfar.append(nuf.split())

		far = np.array(sfar, dtype = float) * -1.0
		far = tuple(map(tuple, far))

		ename = resultname + '/energy'

		if not os.path.exists(ename):
			writefile.werror(outname, 'energy', starttime)
			print('no energy in outfile')
			exit(1)

		resultfile = open(ename).readlines()

		earray = np.array(resultfile[-2].split(), dtype = float)
		
		if earray[-1] > 0:
			genergy = earray[1]
			potential.append(genergy + earray[-1])

		else:
			genergy = earray[1]
			potential.append(genergy)

		return far, potential, genergy

#-------------------------------------------------------------------------------------

	def getsoc(self, filename, outname, state, nstate, spin, theory, intcycle, starttime):
		exit('turbmole cannot calculate soc on-the-fly')

#####################################################################################
#                              class dftbplus                                       #
#####################################################################################

class dftbp():
	'''
	use the molpro soft
	'''
	def __init__(self):
		self = self

#-------------------------------------------------------------------------------------

	def getname(self, filename):
		outname = filename

		return outname

#-------------------------------------------------------------------------------------

	def getrename(self, iname, oname):
		iname += '/dftb_in.hsd'
		oname += '/dftb_in.hsd'

		return iname, oname

#-------------------------------------------------------------------------------------

	def sample(self, filename, pos0, mom0, freq):
		runfile = filename
		newinpf = filename + '0'

		if os.path.isfile(newinpf):
			os.remove(newinpf)

		if not os.path.isdir(newinpf):
			os.mkdir(newinpf)

		shutil.copy('{}/geom'.format(runfile), newinpf)
		shutil.copy('{}/dftb_in.hsd'.format(runfile), newinpf)

		return newinpf, pos0, mom0

#-------------------------------------------------------------------------------------

	def reinpotfile(self, filename, newpos, symb, intcycle, numele):
		'''
		use the new position to create a new in file
		'''
		if intcycle:
			oldname = self.getname(filename + str(intcycle - 1))
			newname = self.getname(filename + str(intcycle))
		else:
			oldname = self.getname(filename)
			newname = self.getname(filename + str(intcycle))

		if not os.path.isdir(newname):
			os.mkdir(newname)

		oldfile = open('{}/geom'.format(oldname)).readlines()
		newfile = open('{}/geom'.format(newname), 'w')

		newpos = np.array(newpos, dtype = float) * autoan

		ditele = {}
		arrele = oldfile[1].split()

		for i, n in enumerate(arrele):
			ditele[n] = i + 1

		newfile.write(oldfile[0])
		newfile.write(oldfile[1])
		for i in range(1, numele + 1):
			newfile.write(' {2:<5}{0:<3}{1[0]:15.8f}{1[1]:15.8f}{1[2]:15.8f}\n'.format(ditele[symb[i - 1]], newpos[i - 1], i))			
		
		newfile.close()

		shutil.copy('{}/dftb_in.hsd'.format(oldname), newname)

		return newname

#-------------------------------------------------------------------------------------

	def getfre(self, outfile):
		'''
		get the frequence
		'''
		freq = 0
	
		return freq

#-------------------------------------------------------------------------------------

	def getpos(self, filename, numele):
		'''
		from the first gjf file to get the first position
		'''
		oldfile = open('{}/geom'.format(filename)).readlines()

		ditele = {}
		arrele = oldfile[1].split()

		for i, n in enumerate(arrele):
			ditele[i + 1] = n

		pos = []
		sym = []
		for i in range(2, numele + 2):
			posarray = oldfile[i].split()
			pos.append(posarray[2:])
			symbol = ditele[int(posarray[1])].title()

			sym.append(symbol)
			
		pos = np.array(pos, dtype = float) / autoan
		sym = tuple(sym)
		pos = tuple(map(tuple, pos))

		return pos, sym

#-------------------------------------------------------------------------------------

	def jpfile(self, oldname, outname, symb, newstate, state, jumpword, theory, spin, trans, starttime):
		newname = oldname + jumpword

		if not os.path.isdir(newname):
			os.mkdir(newname)

		oldfile = open('{}/dftb_in.hsd'.format(oldname)).read()
		newfile = open('{}/dftb_in.hsd'.format(newname), 'w')

		oldstr = oldfile.split('\n\n')
		while '' in oldstr:
			oldstr.remove('')

		if trans:
			olds = True if spin == 'S' else False
			news = not olds if trans else olds
	
			oldspin = 'Singlet' if olds else 'Triplet'
			newspin = 'Singlet' if news else 'Triplet'

			if theory == "lda":
				condir = ldaspin
			elif theory == "pbe":
				condir = pbespin
			else:
				writefile.werror(outname, 'dftbplus theory', starttime)
				print('error theory in dftbplus')
				exit(1)
				
			conline = ' SpinConstants={\n'
			for n in set(symb):
				conline += '  {0}={{\n   {1}\n  }}\n'.format(n, condir[n])
			conline += ' }\n'
				

		if state == 0:
			exstr = 'ExcitedState={\n  Casida={\n    NrOfExcitations=12\n    StateOfInterest=1\n    Symmetry=Singlet\n  }\n}'
			oldstr.append(exstr)

		elif newstate == 0:
			for line in oldstr:
				if 'ExcitedState' in line:
					oldstr.remove(line)

		elif trans:
			for i, line in enumerate(oldstr):
				if 'Hamiltonian' in line:
					if spin == 'S':
						nket = line.rfind('}')
						oldstr[i] = line[:nket] + conline + line[nket:]
					elif spin == 'T':
						oldstr[i] = line.replace(conline, '')
						
				if 'Symmetry' in line:
					oldstr[i] = line.replace(oldspin, newspin)

		else:
			for line in oldstr:
				if 'StateOfInterest' in line:
					oldstr[i] = line.replace(str(state), str(newstate))

		newline = '\n\n'.join(oldstr)
		newfile.write(newline)

		newfile.close()

		shutil.copy('{}/geom'.format(oldname), newname)

		return newname

#-------------------------------------------------------------------------------------

	def gap(self, outname, oup, inp, newstate, state, nstate, deltav, start, flip, spin, theory, trans, symb):
		afile = oup

		if state == 0:
			inpe = self.jpfile(inp, outname, symb, 1, 0, 'upper', theory, spin, False, start)
			oupe = self.runsoft(inpe, state, 'r')
			afile = oupe

		if trans:
			inpe = self.jpfile(inp, outname, symb, newstate, state, 'trans', theory, spin, True, start)
			oupe = self.runsoft(inpe, state, 'r')
			afile = oupe

		gfile = afile + '/detailed.out'
		efile = afile + '/EXC.DAT'
		if not os.path.isfile(efile) or not os.path.isfile(gfile):
			writefile.werror(outname, 'deltav', start)
			print('no gap energy in outfile')
			exit(1)
			
		enef = open(efile)
		enefile = enef.readlines()
		enef.close()

		with open(gfile) as gf:
			gefile = gf.readlines()

		for line in gefile:
			if 'Total energy' in line:
				genergy = float(line.split()[2])

		dev = []
		excitarray = [0.0]
		for line in enefile:
			if '->' in line:
				gapenergy = float(line.split()[0])
				excitarray.append(gapenergy)

		allpes = np.array(excitarray[:nstate], dtype = float) * evtoau + genergy
		for i in range(1, nstate):
			dev.append(excitarray[i] - excitarray[i - 1])

		if state == 0:
			deltav.append([dev[0]])

		elif state == max(range(nstate)):
			deltav.append([dev[-1]])

		else:
			deltav.append(dev[state - 1: state + 1])

		return deltav, allpes

#-------------------------------------------------------------------------------------

	def runsoft(self, runfile, state, keysoft):
		'''
		for run the gaussian progam, and return the log file
		'''
		outname = runfile

		if keysoft == 'r':
			os.chdir(runfile)
			os.system('dftb+ >> out')
			os.chdir('..')

		return outname

#------------------------------------------------------------------------------------

	def parseout(self, resultname, numele, outname, state, potential, starttime):
		newname = resultname + '/detailed.out'

		if not os.path.exists(newname):
			writefile.werror(outname, 'energy', starttime)
			print('no energy in outfile')
			exit(1)

		resultfile = open(newname).readlines()

		signe = 0
		excit = 0
		signf = 0
		error = False
		for i, line in enumerate(resultfile):
			if 'Total energy' in line:
				signe = float(line.split()[2])

			if 'Excitation Energy' in line:
				excit = float(line.split()[2])

			if 'Total Forces' in line:
				signf = i

			if 'NOT converged' in line:
				error = True

		if error:
			writefile.werror(outname, 'forces', starttime)
			print('no force in outfile')
			exit(1)

		if signe:
			genergy = signe
			potential.append(signe + excit)
		else:
			writefile.werror(outname, 'energy', starttime)
			print('no energy in outfile')
			exit(1)

		if 'Total Forces' in resultfile[signf]:
			numf = resultfile[signf + 1: signf + 1 + numele]
		else:
			writefile.werror(outname, 'forces', starttime)
			print('no force in outfile')
			exit(1)

		sfar = []
		for i in range(numele):
			eachforce = numf[i].split()[1: 4]
			sfar.append(eachforce)

		far = np.array(sfar, dtype = float)
		far = tuple(map(tuple, far))

		return far, potential, genergy

#-------------------------------------------------------------------------------------

	def getsoc(self, filename, outname, state, nstate, spin, theory, intcycle, starttime):
		exit('dftbplus cannot calculate soc on-the-fly')

#####################################################################################
#                              class BDF                                            #
#####################################################################################

class BDF():
	'''
	use the BDF soft
	'''
	def __init__(self):
		self = self

#-------------------------------------------------------------------------------------

	def getname(self, filename):
		outname = filename + '.inp'

		return outname

#-------------------------------------------------------------------------------------

	def getrename(self, iname, oname):

		return iname, oname

#-------------------------------------------------------------------------------------

	def sample(self, filename, pos0, mom0, freq):
		runfile = filename.split('.')[0] + '.inp'
		newinpf = filename.split('.')[0] + '0.inp'

		oldfile = open(runfile).readlines()
		newfile = open(newinpf, 'w')

		for line in oldfile:
			newfile.write(line.lower())

		newfile.close()

		return newinpf, pos0, mom0

#-------------------------------------------------------------------------------------

	def reinpotfile(self, filename, newpos, symb, intcycle, numele):
		'''
		use the new position to create a new inp file
		'''
		if intcycle:
			oldname = self.getname(filename + str(intcycle - 1))
			newname = self.getname(filename + str(intcycle))
		else:
			oldname = self.getname(filename)
			newname = self.getname(filename + str(intcycle))

		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		newpos = np.array(newpos, dtype = float) * autoan

		geoline = 0
		for i, line in enumerate(oldfile):
			line = line.lower()
			if 'geometry' in line and 'end' not in line:
				geoline = i + 1

		for i in range(geoline): 
			newfile.write(oldfile[i].lower())

		for i in range(geoline, geoline + numele): 
			newfile.write(' {0:<5}{1[0]:15.8f}{1[1]:15.8f}{1[2]:15.8f}\n'.format(symb[i - geoline], newpos[i - geoline]))

		for i in range(geoline + numele, len(oldfile)):
			newfile.write(oldfile[i].lower())

		newfile.close()
		return newname

#-------------------------------------------------------------------------------------

	def getpos(self, filename, numele):
		'''
		from the first inp file to get the first position
		'''
		runfile = filename + ".inp"
		oldfile = open(runfile).readlines()
		filelen = len(oldfile)

		geoline = 0
		for i, line in enumerate(oldfile):
			line = line.lower()
			if 'geometry' in line and "end" not in line:
				geoline = i + 1

		numpos  = oldfile[geoline: geoline + numele]

		pos = []
		sym = []
		for i in range(numele):
			posarray = numpos[i].split()
			symbol = posarray[0].title()

			pos.append(posarray[1:])
			sym.append(symbol)

		pos = np.array(pos, dtype = float) / autoan
		sym = tuple(sym)
		pos = tuple(map(tuple, pos))

		return pos, sym
	
#-------------------------------------------------------------------------------------

	def jpfile(self, oldname, outname, symb, newstate, state, jumpword, theory, spin, trans, starttime):
		newname = self.getname(oldname.split(".")[0] + jumpword)
		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		if trans:
			oldstr = True if spin == 'S' else False
			newstr = not oldstr if trans else oldstr
	
			oldspin = 0 if oldstr else 1
			newspin = 0 if newstr else 1

		respline0 = "$resp\ngeom\n$end\n"
		respline1 = "$resp\ngeom\nimethod\n 2\nnfiles\n 1\n$end\n"

		itda = 0
		itest = 0
		icorrect = 0
		if "tda" in theory or "cis" in theory:
			itda = 1
		if "xtd" in theory:
			itest = 1
			icorrect = 1

		tdline = "$tddft\nisf\n 0\n"
		if itda:
			tdline += "itda\n 1\n"
		if itest:
			tdline += "imethod\n 2\nitest\n 1\nicorrect\n 1\n"

		tdline += "nroot\n 8\nistore\n 1\n$end\n"

		tdline += respline1

		scfline = 0
		scfend = False
		writetd = False

		tdend = False
		passtd = False
		tdisf = False

		respend = False
		resproot = False
		resp0 = False

		passline = False
		for line in oldfile:
			# ==== find scf module ====
			if "$scf" in line:
				scfend = True
			
			if "$end" in line and scfend:
				writetd = True
				scfend = False
			# =========================

			# ==== find tddft module ====
			if "$tddft" in line:
				tdend = True

			if tdend and "isf" in line:
				tdisf = True

			if "$end" in line and tdend:
				passtd = True
			# ===========================

			# ==== find resp module ====
			if "$resp" in line:
				respend = True

			if respend and "root" in line:
				resproot = True

			if "$end" in line and respend:
				resp0 = True
			# ==========================

			if state == 0:
				if writetd:
					line += tdline
					writetd = False

				if respend:
					if resp0:
						respend = False

					continue

			elif newstate == 0 and tdend:
				if passtd:
					tdend = False
				continue
			
			elif trans:
				if passline:
					passline = False
					continue

				if tdisf:
					line += "{}\n".format(newspin)
					tdisf = False
					passline = True

				if resproot:
					line += "{}\n".format(newstate)
					resproot = False
					passline = True

			elif state != newstate:
				if passline:
					passline = False
					continue

				if resproot:
					line += "{}\n".format(newstate)
					resproot = False
					passline = True

			newfile.write(line)

		newfile.close()

		return newname

#-------------------------------------------------------------------------------------

	def gap(self, outname, oup, inp, newstate, state, nstate, deltav, start, flip, spin, theory, trans, symb):
		efile = oup
		if state == 0:
			inpe = self.jpfile(inp, outname, symb, 1, 0, 'upper', theory, spin, False, start)
			oupe = self.runsoft(inpe, state, 'r')
			efile = oupe

		if trans:
			inpe = self.jpfile(inp, outname, symb, newstate, state, 'trans', theory, spin, True, start)
			oupe = self.runsoft(inpe, state, 'r')
			efile = oupe

		enefile = open(efile).readlines()

		genergy = 0.0
		for i, line in enumerate(enefile):
			if 'E_tot' in line:
				genergy = float(line.split()[-1])

			if 'ExEnergies' in line:
				exeline = i

		excitarray = [0.0]
		for line in enefile[exeline + 2: exeline + nstate + 2]:
			exene = float(line.split()[4])
			excitarray.append(exene)

		excitarray = np.array(excitarray, dtype = float)
		allpes = excitarray[:nstate] * evtoau + genergy

		dev = []
		for i in range(1, nstate):
			dev.append(excitarray[i] - excitarray[i - 1])

		if not dev:
			writefile.werror(outname, 'deltav', start)
			print('no gap energy in outfile')
			exit(1)

		if not trans:
			if spin == 'S' and state == 0:
				deltav.append([dev[0]])
	
			elif spin == 'T' and state == 1:
				deltav.append([dev[1]])

			elif state == max(range(nstate)):
				deltav.append([dev[-1]])
	
			else:
				deltav.append(dev[state - 1: state + 1])

		return deltav, allpes

#-------------------------------------------------------------------------------------

	def runsoft(self, runfile, state, keysoft):
		'''
		for run the gaussian progam, and return the log file
		'''
		outname = runfile.split('.')[0] + '.out'

		if keysoft == 'r':
			os.system('bdfdrv.py %s > %s' %(runfile, outname))

		return outname

#------------------------------------------------------------------------------------

	def parseout(self, resultname, numele, outname, state, potential, starttime):
		resultfile = open(resultname).readlines()

		signf = 0
		signe = []
		for i, line in enumerate(resultfile):
			if 'Gradient contribution from GS+ex' in line or 'Gradient contribution from Tot-egrad' in line:
				signf = i

			if 'E_tot' in line:
				genergy = float(line.split()[-1])
				signe.append(line)

			if 'EXGRAD_estate=' in line:
				signe.append(line)

		if 'Gradient contribution from GS+ex' in resultfile[signf] or 'Gradient contribution from Tot-egrad' in resultfile[signf]:
			numf = resultfile[signf + 1: signf + 1 + numele]
		else:
			writefile.werror(outname, 'forces', starttime)
			print('no force in outfile')
			exit(1)

		sfar = []
		for i in range(numele):
			eachforce = numf[i].split()
			sfar.append(eachforce[1:])

		sfar = np.array(sfar, dtype = float) * -1.0
		far = tuple(map(tuple, sfar))
		
		ii = signe[-1]
		if signe:
			energyrow = ii.split()[-1]

		else:
			writefile.werror(outname, 'energy', starttime)
			print('no energy in outfile')
			exit(1)

		potential.append(float(energyrow))

		return far, potential, genergy

#-------------------------------------------------------------------------------------

	def getsoc(self, filename, outname, state, nstate, spin, theory, intcycle, starttime):
		exit('BDF calculate soc on-the-fly not implement')

#####################################################################################
#                             class oniom                                           #
#####################################################################################

class oniom(gaussian):
	'''
	use the gaussian oniom 
	'''
	def __init__(self):
		self = self

#-------------------------------------------------------------------------------------

	def reinpotfile(self, filename, newpos, symb, intcycle, numele):
		'''
		use the new position to create a new gjf file
		'''
		oldname = filename + str(intcycle - 1) + '.gjf'
		newname = filename + str(intcycle) + '.gjf'

		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		newpos = np.array(newpos, dtype = float) * autoan
		linenu = 0 # ensure the positine line number.
		linena = 0
		for i, line in enumerate(oldfile):
			if line.strip() == "":
				linena += 1

			if linena == 2:
				linenu = i + 2
				break

		for i in range(linenu): 
			newfile.write(oldfile[i])

		for i in range(linenu, linenu + numele): 
			newline = oldfile[i].split()
			newline[2: 5] = newpos[i - linenu]
			newfile.write(' {0[0]:<32s}{0[1]:<5s}{0[2]:15.8f}{0[3]:15.8f}{0[4]:15.8f}{0[5]:>2s}\n'.format(newline))

		for i in range(linenu + numele, len(oldfile)):
			newfile.write(oldfile[i])

		newfile.close()
		return newname

#-------------------------------------------------------------------------------------

	def getpos(self, filename, numele):
		'''
		from the first gjf file to get the first position
		'''
		runfile  = filename + ".in"
		oldfile = open(runfile).readlines()
		filelen = len(oldfile)

		linenu = 0 # ensure the positine line number.
		linena = 0
		for i, line in enumerate(oldfile):
			if line.strip() == "":
				linena += 1

			if linena == 2:
				linenu = i + 2
				break

		numpos  = oldfile[linenu: linenu + numele]

		pos = []
		sym = []
		for i in range(numele):
			posarray = numpos[i].split()
			symbol = posarray[0].split('-')[0].title()

			pos.append(posarray[2: 5])
			sym.append(symbol)

		pos = np.array(pos, dtype = float) / autoan
		sym = tuple(sym)
		pos = tuple(map(tuple, pos))

		return pos, sym

#-------------------------------------------------------------------------------------

	def jpfile(self, oldname, outname, symb, newstate, state, jumpword, theory, spin, trans, starttime):
		newname = self.getname(oldname.split(".")[0] + jumpword)
		oldfile = open(oldname).readlines()
		newfile = open(newname, 'w')

		if trans:
			oldstr = True if spin == 'S' else False
			newstr = not oldstr if trans else oldstr
	
			oldspin = 'singlet' if oldstr else 'triplet'
			newspin = 'singlet' if newstr else 'triplet'

		for line in oldfile:
			if line[0] == '#':
				if state == 0:
					line = line.replace(':', ' {}(singlet,root=1,nstate=8):'.format(theory))

				elif newstate == 0:
					line = re.sub(' {}\(.+\)'.format(theory), '', line)
				
				elif trans:
					line = line.replace(oldspin, newspin)
					line = line.replace('root=%d' %state, 'root=%d' %newstate)

				elif state != newstate:
					line = line.replace('root=%d' %state, 'root=%d' %newstate)

			newfile.write(line)

		newfile.close()

		return newname

#-------------------------------------------------------------------------------------

	def parseout(self, resultname, numele, outname, state, potential, starttime):
		resultfile = open(resultname).readlines()

		level = 0
		signf = 0
		signe = []
		oniom = []
		for i, line in enumerate(resultfile):
			if 'low   level' in line:
				level = i

			if 'Integrated Forces ' in line:
				signf = i

			if 'SCF Done' in line:
				genergy = float(line.split()[4])
				signe.append(line)

			if 'Total Energy' in line:
				signe.append(line)

			if 'ONIOM: gridpoint  ' in line:
				oenergy = float(line.split()[-1])
				oniom.append(oenergy)

		fixarray = []
		fixline = resultfile[level + 1: level + 1 + numele]
		for i, nline in enumerate(fixline):
			num = nline.split()[1]
			if num == '-1':
				fixarray.append(i)

		if 'Forces ' in resultfile[signf]:
			numf = resultfile[signf + 3: signf + 3 + numele]
		else:
			writefile.werror(outname, 'forces', starttime)
			print('no force in outfile')
			exit(1)

		sfar = []
		for i in range(numele):
			eforce = np.zeros(3) if i in fixarray else numf[i].split()[2:]
			sfar.append(eforce)

		sfar = np.array(sfar, dtype = float)
		far = tuple(map(tuple, sfar))
		
		ii = signe[-1]
		if 'SCF Done' in ii:
			energyrow = ii.split()[4]

		elif 'Total Energy' in ii:
			energyrow = ii.split()[-1]

		else:
			writefile.werror(outname, 'energy', starttime)
			print('no energy in outfile')
			exit(1)

		genergy = genergy - oniom[0] + oniom[2]
		oniomrow = float(energyrow) - oniom[0] + oniom[2]
		potential.append(oniomrow)

		return far, potential, genergy

#-------------------------------------------------------------------------------------

	def getsoc(self, filename, outname, state, nstate, spin, theory, intcycle, starttime):
		exit('oniom cannot calculate soc on-the-fly')


