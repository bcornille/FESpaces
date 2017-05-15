#!/usr/bin/env python

import matplotlib as mpl
#mpl.use('Cairo')
import matplotlib.pyplot as plt
import numpy as np
import sys, getopt
import json
import math as m
from matplotlib import rc
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D as ax
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

	# Input specs and analytical solution for P
	function = params["Function"]
	formulation = params["Formulation"]
	orderP = params["Mesh"]["el_order"]
	orderu = int(orderP) - 1
	if formulation == "Mixed":
		orderu = orderu + 2
	elemsx = int(params["Mesh"]["N_x"])
	elemsy = int(params["Mesh"]["N_y"])
	filetype = params["Plot"]["Filetype"]
	xmin = params["Mesh"]["x_min"]
	xmax = params["Mesh"]["x_max"]
	ymin = params["Mesh"]["y_min"]
	ymax = params["Mesh"]["y_max"]
	dx = (xmax-xmin)/1000.0
	dy = (ymax-ymin)/1000.0

	# Gold solution
	goldxf = np.arange(xmin, xmax+0.001, dx) 
	goldyf = np.arange(ymin, ymax+0.001, dy)
	XG, YG = np.meshgrid(goldxf, goldyf)
	PG = np.multiply(np.sin(m.pi*XG),np.sin(m.pi*YG));

	# Calculated solution
	X = output["xgrid"]
	Y = output["ygrid"]
	P = output["Pgrid"]

	# Plot gold
	f1 = plt.figure(1,figsize=(8,8))
	ax = f1.add_subplot(111, projection='3d')
	surf = ax.plot_surface(XG,YG,PG, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	f1.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()
	plt.savefig("%s%s" % ("./plot/GoldP.",filetype))

	# Plot calculated
	f2 = plt.figure(2,figsize=(8,8))
	ax = f2.add_subplot(111, projection='3d')
	surf = ax.plot_surface(X,Y,P, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	f2.colorbar(surf, shrink=0.5, aspect=5)
	plt.show()
	plt.savefig("%s%s" % ("./plot/CalcP.",filetype))

	if formulation != "Standard":
		U = output["uxgrid"]
		V = output["uygrid"]

		f3 = plt.figure(3,figsize=(8,8))
		ax = f3.add_subplot(111, projection='3d')
		surf = ax.plot_surface(X,Y,U, cmap=cm.coolwarm, linewidth=0, antialiased=False)
		f3.colorbar(surf, shrink=0.5, aspect=5)
		plt.show()
		plt.savefig("%s%s" % ("./plot/CalcUx.",filetype))

		f4 = plt.figure(4,figsize=(8,8))
		ax = f4.add_subplot(111, projection='3d')
		surf = ax.plot_surface(X,Y,V, cmap=cm.coolwarm, linewidth=0, antialiased=False)
		f4.colorbar(surf, shrink=0.5, aspect=5)
		plt.show()
		plt.savefig("%s%s" % ("./plot/CalcVx.",filetype))


if __name__ == "__main__":
   main(sys.argv[1:])
