#!/usr/bin/env python

import matplotlib as mpl
#mpl.use('Cairo')
import matplotlib.pyplot as plt
import numpy as np
import sys, getopt
import json
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def main(argv):

	# Read specifications for the c++ json input and output files
	ifile = ''
	ofile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print 'pyplot.py -i <inputfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'pyplot.py -i <inputfile> -o <outputfile>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			ifile = arg
		elif opt in ("-o", "--ofile"):
			ofile = arg

	# Load the data from the json files
	with open(ifile) as in_file:
		params = json.load(in_file)
	with open(ofile) as out_file:
		output = json.load(out_file)

	# Analytical solution for P and u
	goldx = np.arange(params["Mesh"]["x_min"],params["Mesh"]["x_max"]+0.001, 0.001)
	function = params["Function"]
	formulation = params["Formulation"]
	if formulation == "Mimetic":
		formulation = "Dual-Mixed"
	orderP = params["Order"]
	orderu = int(orderP) - 1
	if formulation == "Mixed":
		orderu = orderP
		orderP = orderu - 1
	elems = int(params["Mesh"]["N"]) + 1
	filetype = params["Plot"]["Filetype"]
	if function == "ExpX3pX":
		goldP = np.exp(goldx)*goldx*(1-goldx)
		if formulation != "Standard": goldu = np.exp(goldx)*(goldx-1+np.power(goldx,2))
		f = (r'$f = e^{x}x(x-3)$')
	elif function == "Two":
		goldP = goldx*(1-goldx)
		if formulation != "Standard": goldu = 2*goldx-1
		f = (r'$f = 2$')

	# Create output plots
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	mpl.rc('xtick', labelsize=12)
	mpl.rc('ytick', labelsize=12)

	plt.figure(1, figsize=(5,5))
	plt.plot(goldx, goldP, 'k')
	plt.plot(output["xgrid"],output["Pgrid"],'b.')
	plt.xlabel(r'Position', fontsize=14)
	plt.ylabel(r'$p$ value', fontsize=14)
	t1 = (r'$p$ vs $x$, ')
	t2 = (', Formulation: ') 
	t3 = ('\n Poly. degree: ')
	t4 = (r', Elements: ')
	plt.title("%s %s%s%s%s%s" % (t1,f,t3,orderP,t4,elems), fontsize=16)
	plt.grid(True)
	plt.axis([params["Mesh"]["x_min"],params["Mesh"]["x_max"],0,0.5])
	plt.tight_layout()
	plt.savefig("%s%s" % ("./plot/P.",filetype))

	if formulation != "Standard":
		plt.figure(2, figsize=(5,5))
		plt.plot(goldx, goldu, 'k')
		plt.plot(output["xgrid"],output["ugrid"],'r.')
		t1 = (r'$u$ vs $x$, ')
		plt.xlabel(r'Position', fontsize=14)
		plt.ylabel(r'$u$ value', fontsize=14)
		plt.title("%s %s%s%s%s%s" % (t1,f,t3,orderu,t4,elems), fontsize=16)
		plt.grid(True)
		plt.axis([params["Mesh"]["x_min"],params["Mesh"]["x_max"],-1,2.8])
		plt.tight_layout()
		plt.savefig("%s%s" % ("./plot/u.",filetype))


if __name__ == "__main__":
   main(sys.argv[1:])
