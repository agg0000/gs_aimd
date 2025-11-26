#!/usr/bin/env python
'''
for create output file
'''

import numpy as np
import datetime

from constant import autoan, evtoau

def writecon(outname, sym, con, group, numele):
	'''
	writecon means write content.
	because the out file has to many similarity content about force, momentum and position
	so use this function to write down the similarity content in every cycle
	'''
	con = np.array(con, dtype = float)
	con = con * autoan if group == 'coordinate' else con

	outfile = open(outname, 'a+')
	outfile.write('\n' + 'The {}'.format(group).center(63) + '\n')
	outfile.write('-'*63 + '\n')
	outfile.write(' {0:<9}{1:>16}{2:>16}{3:>16}\n'.format('Symbol', 'x', 'y', 'z'))
	outfile.write('-'*63 + '\n')

	for i in range(numele):
		outfile.write('  {0:<12}{1[0]:>16.9f}{1[1]:>16.9f}{1[2]:>16.9f}\n'.format(sym[i], con[i]))

	outfile.write('-'*63 + '\n')
	outfile.close()

	return

#------------------------------------------------------------------------------------

def werror(outname, keyword, starttime):
	'''
	for write down the error infomation
	'''

	endtime = datetime.datetime.now()
	usetime = endtime - starttime

	outfile = open(outname, 'a+')
	outfile.write('*'*99 + '\n')
	outfile.write('ERROR %s' %keyword + '\n')
	outfile.write('Grad'*24 + '\n')
	outfile.write('*'*99 + '\n')
	outfile.write('no {} in outfile\n'.format(keyword))
	outfile.write('start program at {}\n'.format(starttime))
	outfile.write('use time {}\n'.format(usetime))
	outfile.close()

	return

#------------------------------------------------------------------------------------

def endtime(outname, keyword, intcycle, starttime):
	endtime = datetime.datetime.now()
	usetime = endtime - starttime

	outfile = open(outname, 'a+')
	outfile.write('\n' + '*'*75 + '\n')
	outfile.write('cycle %d normal exit' %intcycle + '\n')
	outfile.write('Start program at  ' + str(starttime) + '\n')
	outfile.write('Use time  ' + str(usetime) + '\n')
	outfile.write('%s normal exit\n' %keyword)
	outfile.close()
	
	return

#------------------------------------------------------------------------------------

def outstate(outname, state):
	outfile = open(outname, 'a+')
	outfile.write('the calculate state is %s' %state)
	outfile.close()

	return

#------------------------------------------------------------------------------------

def outenergy(outname, energy, keyword):
	outfile = open(outname, 'a+')
	outfile.write('the total %s energy is %f a.u. %.6f (eV)\n' %(keyword, energy, energy / evtoau))
	outfile.close()

	return

#------------------------------------------------------------------------------------

def adiabatic(outname):
	outfile = open(outname, 'a+')
	outfile.write('\n  adiabatic molecular dynamic\n')
	outfile.close()

	return

#------------------------------------------------------------------------------------

def whop(outname, intcycle, state, spin, oldstate, oldspin):
	outfile = open(outname, 'a+')
	outfile.write('\nthe state at %d from %s%d hop to state %s%d \n' %(intcycle, oldspin, oldstate, spin, state))
	outfile.close()

	return

#------------------------------------------------------------------------------------

def writeab(outname, asquare, bsquare, keyword):
	outfile = open(outname, 'a+')
	outfile.write('\nusing {} hopping algorithm\n'.format(keyword))
	outfile.write('effective nonadiabatic coupling asquare is {:.9f}\n'.format(asquare))
	outfile.write('effective collision energy bsquare is {:.9f}\n'.format(bsquare))
	outfile.close()

	return

#------------------------------------------------------------------------------------

def wpes(outname, allpes, genergy, spin):
	outfile = open(outname, 'a+')
	outfile.write('\n')

	for i, energy in enumerate(allpes):
		outfile.write('the %s%d state pes energy is %f a.u.\n' %(spin, i + 1, energy))

	outfile.close()

	return

#------------------------------------------------------------------------------------

def ncycle(outname, intcycle):

	outfile = open(outname, 'a+')
	outfile.write('\n' + '='*72 + '\n')
	outfile.write('Cycle {:d}'.format(intcycle).center(72) + '\n')
	outfile.write('='*72 + '\n'*2)
	outfile.close()

	return

#------------------------------------------------------------------------------------

def outground(outname, genergy):
	outfile = open(outname, 'a+')
	outfile.write('\nthe ground state energy is %f a.u. %f (eV) \n' %(genergy, genergy / evtoau))
	outfile.close()

	return

#------------------------------------------------------------------------------------

def spindr(outname, spin, state):
	outfile = open(outname, 'a+')
	outfile.write('\nthe spin state is %s%d \n' %(spin, state))
	outfile.close()

	return

#------------------------------------------------------------------------------------

def cyclend(outname, intcycle):
	outfile = open(outname, 'a+')
	outfile.write('\n %d cycle end \n' %intcycle)
	outfile.close()

	return

#------------------------------------------------------------------------------------

def trunline(outname):
	outfile = open(outname, 'rb+')
	allfile = outfile.readlines()
	endfile = len(allfile[-1])
	outfile.seek(-(endfile + 1), 1)
	outfile.truncate()
	outfile.close()

	return

#------------------------------------------------------------------------------------

def wdashed(outname):
	outfile = open(outname, 'a+')
	outfile.write('\n' + '-'*65 + '\n')
	outfile.close()

	return

#------------------------------------------------------------------------------------

def masscen(outname):
	outfile = open(outname, 'a+')
	outfile.write('  Some information about before mcenter.\n')
	outfile.close()

	return

#------------------------------------------------------------------------------------

