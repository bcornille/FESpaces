#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import json

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

size = 5.0
lwidth = 2
bfsize = 16
fsize = 14

handles = []
markers = ["o", "s", "^"]
plt.figure(figsize=(size, size))

for i, formulation in enumerate(['standard', 'mixed', 'dual']):
    err_temp_p = []
    for p in range(1, 10):
        with open('./' + formulation + '_p' + str(p) +
                  '_N04.json') as out_file:
            data_out = json.load(out_file)
        err_temp_p.append(data_out["error"])
    err = np.array(err_temp_p)
    plt.semilogy(range(1, 10), err, marker=markers[i],
                 color='k', label=r'{}'.format(formulation))

plt.title(r'Convergence for $p$-refinement', fontsize=bfsize)
plt.xlabel(r'Polynomial Degree $k$', fontsize=fsize)
plt.ylabel(r'$\left\|p_h - p\right\|_2$', fontsize=fsize)
plt.legend(loc='best', fontsize=8)
plt.tight_layout()

plt.savefig('convergence_p.pdf')

# plt.show()
