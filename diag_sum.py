#!/usr/bin/env python2

#This takes a NxN matrix file and sums each diagonal parallel to the
#main diagonal, and normalizes each element with the one in the main
#diagonal. Not all elements are included, only those that are on a
#"perpendicular" line from the main diagonal. Point is we want to see
#how fast the (relative) element sizes dropp from the main diagonal.

import math
import sys

def get_matrix(file_name):
    "Read in all non commented lines and return columns as a matrix"
    matrix = []
    myfile = open(file_name,'r')
    for line in myfile.readlines():
        li=line.strip()               #remove leading whitespace
        if not li.startswith("#"):    #ignore comments
            tmp = [float(value) for value in line.split()]
            matrix.append(tmp)
    myfile.close()
    return matrix

def make_matrix_small(input, numberOfElementsToSkipp):
    """Remove every n elements. (to check convergence of different
    matrices)"""
    out = []
    for i in range(0, len(input), numberOfElementsToSkipp):
        out.append(input[i][::numberOfElementsToSkipp])
    return out

#======== MAIN PROGRAM ============================

if len(sys.argv) != 2 and len(sys.argv) != 3:
    sys.stderr.write("Usage:\t%s matrix-file\n" % sys.argv[0])
    sys.stderr.write("or:\t%s matrix-file every_n_value\n" % sys.argv[0])
    sys.exit()

if len(sys.argv) == 3:
    skipp_every_n_data = int(sys.argv[2])
    H = make_matrix_small(get_matrix(sys.argv[1]), skipp_every_n_data)
    print("# skipping every %g element" % skipp_every_n_data)
else:
    H = get_matrix(sys.argv[1])

N = int(len(H))
save = []

for k in range(int(N/2)):
    save.append(0.0)
    for i in range(k,N-k):
        save[k] += H[i-k][i+k]/H[i][i]
    save[k] = save[k] / (N-k)

#pprint.pprint(normalization)
print("# off diagonal \t value")
for i in range(len(save)):
    print("%f\t%g" % (1.0*i/len(save), save[i]))
