#!/usr/bin/env python

import matplotlib as mpl
#mpl.use('Cairo')
import matplotlib.pyplot as plt
import numpy as np
import sys, getopt
import json
import subprocess
import math
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def main(argv):

	# Read specifications for the c++ json input and output files
	ifile = ''
	ofile = ''
	cfile = 'out2d.json'
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print './plot/conv2.py -i <inputfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './plot/conv2.py -i plot/in2dp.json -o in2dc.json'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			ifile = arg
		elif opt in ("-o", "--ofile"):
			ofile = arg
	# i.e. run ./plot/conv2.py -i plot/inp.json -o inc.json from src

	# Load the data from the json and prepare the c++ input json
	with open(ifile) as in_file:
		params = json.load(in_file)
	p1 = params["Converge"]["p_min"]
	p2 = params["Converge"]["p_max"]
	h1 = params["Converge"]["h_min"]
	h2 = params["Converge"]["h_max"]
	p1h = params["Converge"]["p_minh"]
	p2h = params["Converge"]["p_maxh"]
	c2 = params["Converge"]["curve"]

	# Stores the output
	# ((curvature, formulations, {opposite parameter}, parameter), double)
	errorp = np.zeros((2,3,(p2-p1)),dtype=np.float64)
	errorh = np.zeros((2,3,(p2h-p1h),(h2-h1)),dtype=np.float64)
	nelems = np.zeros(h2-h1,int)

	# Loop through the desired tests
	for c, curve in enumerate([0.0, c2]):
		params["Mesh"]["c"] = curve
		for f, formulation in enumerate(['Standard', 'Mixed', 'Dual-Mixed']):
			print ('')
			print "p convergence 1 thru 9"
			print ("%s%s%s%s" % ('c= ',curve,' f= ',formulation))
			for p in range(p1,p2):
				params["Formulation"] = formulation
				if formulation == "Mixed":
					params["Mesh"]["el_order"] = p+1
				else:
					params["Mesh"]["el_order"] = p

				# Run the Poisson code
				with open(ofile, 'w') as outfile:
					json.dump(params, outfile)
				subprocess.call(["./Poisson2D", ofile, cfile])

				# Disregard outputs, acquire error
				with open(cfile) as c_file:
					result = json.load(c_file)
				errorp[c][f][p-p1] = result["error"]

	for c, curve in enumerate([0.0, c2]):
		params["Mesh"]["c"] = curve
		for f, formulation in enumerate(['Standard', 'Mixed', 'Dual-Mixed']):
			for p in range(p1h,p2h):
				print ('')
				print "h convergence (2^1,2^1) thru (2^4,2^4)"
				print ("%s%s%s%s%s%s" % ('c= ',curve,' f= ',formulation,' p= ', p))
				for h in range(h1,h2):
					nelems[h-h1] = int(math.pow(2,h))
					params["Mesh"]["N_x"] = nelems[h-h1]
					params["Mesh"]["N_y"] = nelems[h-h1]
					params["Formulation"] = formulation
					if formulation == "Mixed":
						params["Mesh"]["el_order"] = p+1
					else:
						params["Mesh"]["el_order"] = p

					# Run the Poisson code
					with open(ofile, 'w') as outfile:
						json.dump(params, outfile)
					subprocess.call(["./Poisson2D", ofile, cfile])

					# Disregard outputs, acquire error
					with open(cfile) as c_file:
						result = json.load(c_file)
					errorh[c][f][p-p1h][h-h1] = result["error"]

	# Generate plots
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	mpl.rc('xtick', labelsize=12)
	mpl.rc('ytick', labelsize=12)
	size = 5.0
	bfsize = 16
	fsize = 14
	lsize = 8
	markers = ["o", "s", "^"]
	colors = ["r", "m", "g", "c", "b", "k"]
	handles = []

	# p plots
	for c, curve in enumerate([0.0, c2]):
		plt.figure((c+1), figsize=(size,size))
		for f, formulation in enumerate(['standard', 'mixed', 'dual']):
			plt.semilogy(range(p1,p2), errorp[c][f][0:(p2-p1)], marker=markers[f], 
				color='k', label=r'{}'.format(formulation))
		plt.title(r'Convergence for $p$-refinement', fontsize=bfsize)
		plt.xlabel(r'Polynomial Degree $k$', fontsize=fsize)
		plt.ylabel(r'$\left\|p_h - p\right\|_2$', fontsize=fsize)
		plt.legend(loc='best', fontsize=lsize)
		plt.tight_layout()
		plt.savefig("%s%s%s" % ('plot/convergence_2dp',curve,'.pdf'))

	# h plots
	for p in range(p1h, p2h):
	    handles.append(plt.Line2D([], [], color=colors[
		              p - 1], label=r'$k = {}$'.format(p)))
	for c, curve in enumerate([0.0, c2]):
		plt.figure((c+3), figsize=(size,size))
		for f, formulation in enumerate(['standard', 'mixed', 'dual']):
			if c == 0:
				handles.append(plt.Line2D([], [], color='k', linestyle='',
					marker=markers[f], label=r'{}'.format(formulation)))
			for p in range(p1h, p2h):
				plt.loglog(nelems, errorh[c][f][p-p1h][0:(h2-h1)], marker=markers[f], 
					color=colors[p-p1h], basex=2)
		plt.title(r'Convergence for $h$-refinement', fontsize=bfsize)
		plt.xlabel(r'$N_{elements}$ per dimension', fontsize=fsize)
		plt.ylabel(r'$\left\|p_h - p\right\|_2$', fontsize=fsize)
		plt.legend(handles=handles, ncol=2, loc='best', fontsize=8)
		plt.tight_layout()
		plt.savefig("%s%s%s" % ('plot/convergence_2dh',curve,'.pdf'))

if __name__ == "__main__":
	main(sys.argv[1:])
