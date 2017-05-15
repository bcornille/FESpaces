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
	elems = elemsx*elemsy
	filetype = params["Plot"]["Filetype"]
	stride = params["Plot"]["SPEc"]
	xmin = params["Mesh"]["x_min"]
	xmax = params["Mesh"]["x_max"]
	ymin = params["Mesh"]["y_min"]
	ymax = params["Mesh"]["y_max"]
	dx = (xmax-xmin)/1000.0
	dy = (ymax-ymin)/1000.0

	# Create output plots
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	mpl.rc('xtick', labelsize=12)
	mpl.rc('ytick', labelsize=12)
	size = 8.0
	bfsize = 16
	fsize = 14
	lsize = 8

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
	f1 = plt.figure(1,figsize=(size,size))
	ax = f1.add_subplot(111, projection='3d')
	surf = ax.plot_surface(XG,YG,PG, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	f1.colorbar(surf, shrink=0.5, aspect=5)
	ax.set_xlabel(r'x position', fontsize=fsize)
	ax.set_ylabel(r'y position', fontsize=fsize)
	ax.set_zlabel(r'$P$ value', fontsize=fsize)
	t1 = (r'$P$ vs position, Exact solution')
	t2 = (', Formulation: ') 
	t3 = ('\n Poly. degree: ')
	t4 = (r', Elements: ')
	plt.title("%s %s%s%s%s%s" % (t1,formulation,t3,orderP,t4,elems), fontsize=bfsize)
	plt.tight_layout()
	plt.savefig("%s%s" % ("./plot/GoldP.",filetype))
	plt.show()

	# Plot calculated
	f2 = plt.figure(2,figsize=(size,size))
	ax = f2.add_subplot(111, projection='3d')
	surf = ax.plot_surface(X,Y,P, cmap=cm.coolwarm, linewidth=0, antialiased=False)
	f2.colorbar(surf, shrink=0.5, aspect=5)
	ax.set_xlabel(r'x position', fontsize=fsize)
	ax.set_ylabel(r'y position', fontsize=fsize)
	ax.set_zlabel(r'$P$ value', fontsize=fsize)
	t1 = (r'$P$ vs position, ')
	t2 = (', Formulation: ') 
	t3 = ('\n Poly. degree: ')
	t4 = (r', Elements: ')
	plt.title("%s %s%s%s%s%s" % (t1,formulation,t3,orderP,t4,elems), fontsize=bfsize)
	plt.tight_layout()
	plt.savefig("%s%s" % ("./plot/CalcP.",filetype))
	plt.show()

	if formulation != "Standard":
		U = output["uxgrid"]
		V = output["uygrid"]

		f3 = plt.figure(3,figsize=(size,size))
		ax = f3.add_subplot(111, projection='3d')
		surf = ax.plot_surface(X,Y,U, cmap=cm.coolwarm, linewidth=0, antialiased=False)
		f3.colorbar(surf, shrink=0.5, aspect=5)
		ax.set_xlabel(r'x position', fontsize=fsize)
		ax.set_ylabel(r'y position', fontsize=fsize)
		ax.set_zlabel(r'$u$ value, x-component', fontsize=fsize)
		t1 = (r'$U_{x}$ vs position, ')
		t2 = (', Formulation: ') 
		t3 = ('\n Poly. degree: ')
		t4 = (r', Elements: ')
		plt.title("%s %s%s%s%s%s" % (t1,formulation,t3,orderP,t4,elems), fontsize=bfsize)
		plt.tight_layout()
		plt.savefig("%s%s" % ("./plot/CalcUx.",filetype))
		plt.show()

		f4 = plt.figure(4,figsize=(size,size))
		ax = f4.add_subplot(111, projection='3d')
		surf = ax.plot_surface(X,Y,V, cmap=cm.coolwarm, linewidth=0, antialiased=False)
		f4.colorbar(surf, shrink=0.5, aspect=5)
		ax.set_xlabel(r'x position', fontsize=fsize)
		ax.set_ylabel(r'y position', fontsize=fsize)
		ax.set_zlabel(r'$u$ value, y-component', fontsize=fsize)
		t1 = (r'$U_{y}$ vs position, ')
		t2 = (', Formulation: ') 
		t3 = ('\n Poly. degree: ')
		t4 = (r', Elements: ')
		plt.title("%s %s%s%s%s%s" % (t1,formulation,t3,orderP,t4,elems), fontsize=bfsize)
		plt.tight_layout()
		plt.savefig("%s%s" % ("./plot/CalcVx.",filetype))
		plt.show()

		Xc = (np.array(X))[::stride, ::stride]
		Yc = (np.array(Y))[::stride, ::stride]
		Uc = (np.array(U))[::stride, ::stride]
		Vc = (np.array(V))[::stride, ::stride]

		f5 = plt.figure(5,figsize=(size,size))
		Q = plt.quiver(Xc, Yc, Uc, Vc)
		plt.xlabel(r'x position', fontsize=fsize)
		plt.ylabel(r'y position', fontsize=fsize)
		t1 = (r'$U$ vector field vs position, ')
		t2 = (', Formulation: ') 
		t3 = ('\n Poly. degree: ')
		t4 = (r', Elements: ')
		plt.title("%s %s%s%s%s%s" % (t1,formulation,t3,orderP,t4,elems), fontsize=bfsize)
		plt.tight_layout()
		plt.savefig("%s%s" % ("./plot/CalcUV.",filetype))
		plt.show()



if __name__ == "__main__":
   main(sys.argv[1:])
