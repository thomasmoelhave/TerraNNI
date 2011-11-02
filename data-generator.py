#!/usr/bin/env python
# encoding: utf-8

import sys
import os
import math
import random
from optparse import OptionParser

# Simple function to write out results in a (which contains three arrays) to file fn
def writeArray(fn, a):
	if fn:
		f = open(fn, 'w')
		for i in range(0,len(a)):
			f.write("%f %f %f %d\n" % (a[i][0],a[i][1],a[i][2],a[i][3]))
		f.close()
	else:
		for i in range(0,len(a)):
			print "%f %f %f %d\n" % (a[i][0],a[i][1],a[i][2],a[i][3])

def linearTF(t,m = 1):
	return t * m

def symetricSmooth(tf,to):
	xo = yo = 0.5
	def retFunc(x,y,t):
		return to + tf(t) * math.cos( math.pow(x-xo,2) + math.pow(y-yo,2) )
	return retFunc

def wall(xmin, xmax, tmin, tmax):
	h = 20
	b = 1
	width = 5
	def retFunc(x,y,t):
		tp = 1 - (1.0 * (tmax - t - tmin) / (tmax - tmin))
		xm = tp * (xmax - xmin) + xmin
		if ( x >= xm and x < xm + width ):
			return h
		return b
	return retFunc


def randData():
	def retFunc(x,y,t):
		return myrand(1,100,0.01)
	return retFunc


def myrand(mn, mx, g):
	m = 1
	while (g * m) != int(g*m):
		m = m*10
	return (1.0 * random.randrange(m*mn, m*mx, m*g)) / m

def fill(f,n,p,ta,fn):
	if fn:
		fl = open(fn, 'w')
	#r = []
	for j in range(0,len(ta)):
		t1 = ta[j]
		print 'time: ', t1, n
		for i in range(0,n):
			x1 = myrand(p['xmin'],p['xmax'],p['gran'])
			y1 = myrand(p['ymin'],p['ymax'],p['gran'])
			if fn:
				fl.write("%f %f %f %d\n" % (x1,y1,f(x1,y1,t1),t1))
			else:
				print "%f %f %f %d\n" % (x1,y1,f(x1,y1,t1),t1)
			#r.append( (x1,y1,f(x1,y1,t1),t1) )
	if fn:
		fl.close()
	#return r


def main():
	parser = OptionParser()
	parser.add_option("-f", "--file", action="store", type="string", dest="filename")
	parser.add_option("-n", type="int", dest="n", default=1000)
	(options, args) = parser.parse_args()

	params = {'xmin': -100, 'ymin': -100, 'xmax': 100, 'ymax': 100, 'gran': 0.001 }

	print "n: ", options.n

	# Only use ONE Of these functions for any given run
	#writeArray( options.filename, fill( symetricSmooth(linearTF,100.0), options.n, params, (0,10) ) )
	#writeArray( options.filename, fill( wall(params['xmin'], params['xmax'],0,100), options.n, params, range(0,110,10) )
	#writeArray( options.filename, fill( randData, 200000000, params, range(0,21,1) ))

	#fill( symetricSmooth(linearTF,100.0), options.n, params, (0,10), options.filename ) 
	ta= range(0,110,1)
	ta.remove(10)
	ta.remove(30)
	ta.remove(50)
	ta.remove(70)
	ta.remove(90)
	#fill( wall(params['xmin'], params['xmax'],0,100), options.n, params, ta, options.filename )
	#fill( randData(), 200000000, params, range(0,21,1), options.filename)
	#fill( randData(), 50000000, params, range(0,11,1), options.filename)  # BIG DATA
	fill( randData(), 10000000, params, range(0,11,1), options.filename) 

if __name__ == '__main__':
	main()
