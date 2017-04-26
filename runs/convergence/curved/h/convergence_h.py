#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import glob
import json

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

size = 5.0
lwidth = 2
bfsize = 16
fsize = 14

handles = []
markers = ["o", "s", "^"]
colors = ["r", "m", "g", "c", "b"]
plt.figure(figsize=(size, size))

for p in range(2, 7):
    handles.append(plt.Line2D([], [], color=colors[
                   p - 2], label=r'$p = {}$'.format(p)))

for i, formulation in enumerate(['standard', 'mixed', 'dual']):
    N_temp_form = []
    err_temp_form = []
    handles.append(plt.Line2D([], [], color='k', linestyle='',
                              marker=markers[i],
                              label=r'{}'.format(formulation)))
    for p in range(2, 7):
        N_temp_p = []
        err_temp_p = []
        for file in sorted(glob.glob('./' + formulation +
                                     '_p' + str(p) + '_N*.json')):
            with open(file) as out_file:
                data_out = json.load(out_file)
            N_temp_p.append(data_out["N_el"])
            err_temp_p.append(data_out["error"])
        N = np.array(N_temp_p)
        err = np.array(err_temp_p)
        plt.loglog(N, err, basex=2, marker=markers[i], color=colors[p - 2])

plt.title(r'Convergence for $h$-refinement', fontsize=bfsize)
plt.xlabel(r'$N_{elements}$', fontsize=fsize)
plt.ylabel(r'$\left\|p_h - p\right\|_2$', fontsize=fsize)
plt.legend(handles=handles, ncol=2, loc='best', fontsize=8)
plt.tight_layout()

plt.savefig('convergence_h_curve.pdf')

# plt.show()
