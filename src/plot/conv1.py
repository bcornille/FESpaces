#!/usr/bin/env python

import matplotlib as mpl
#mpl.use('Cairo')
import matplotlib.pyplot as plt
import numpy as np
import sys, getopt
import json
import subprocess
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def main(argv):

	# Read specifications for the c++ json input and output files
	ifile = ''
	ofile = ''
	cfile = 'out.json'
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print './plot/conv1sh.py -i <inputfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print './plot/conv1sh.py -i plot/inp.json -o inc.json'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			ifile = arg
		elif opt in ("-o", "--ofile"):
			ofile = arg

	# Load the data from the json and prepare the c++ input json
	with open(ifile) as in_file:
		params = json.load(in_file)
	p1 = params["Converge"]["p_min"]
	p2 = params["Converge"]["p_max"]

	# Stores the output
	errorp = np.zeros((3,(p2-p1),2),dtype=np.float64)

	# Loop through the desired tests
	for c, curve in enumerate([0.0, 0.4]):
		params["Mesh"]["c"] = curve
		for f, formulation in enumerate(['Standard', 'Mixed', 'Mimetic']):
			for p in range(p1,p2):
				params["Order"] = p
				params["Formulation"] = formulation

				# Run the Poisson code
				with open(ofile, 'w') as outfile:
					json.dump(params, outfile)
				subprocess.call(["./Poisson1D", ofile, cfile])

				# Disregard outputs, acquire error
				with open(cfile) as c_file:
					result = json.load(c_file)
				errorp[f][p-p1][c] = result["error"]
	

if __name__ == "__main__":
	main(sys.argv[1:])
