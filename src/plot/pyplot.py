#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys, getopt
import json
#from pprint import pprint 

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

	# Analytical solution for P
	goldx = np.arange(params["Mesh"]["x_min"],params["Mesh"]["x_max"]+0.001, 0.001)
	function = params["Function"]
	if function == "ExpX3pX":
		goldP = np.exp(goldx)*goldx*(1-goldx)
	elif function == "Two":
		goldP = goldx*(1-goldx)

	# Create output plot
	plt.plot(goldx, goldP, 'k')
	plt.plot(output["xgrid"],output["Pgrid"])
	title = ('P vs Position, ');
	plt.xlabel('position')
	plt.ylabel('P value')
	plt.title("%s %s" % (title,params["Formulation"]))
	plt.grid(True)
	plt.axis([params["Mesh"]["x_min"],params["Mesh"]["x_max"],0,0.5])
	plt.savefig("./plot/P.png")
	#plt.show()

if __name__ == "__main__":
   main(sys.argv[1:])
