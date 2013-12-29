#!/usr/bin/env python2

#This takes a NxN matrix file and sums each diagonal from center
#diagonal (N elements), out to last diagonal (1 element), and prints
#normalized average value of element in that diagonal.
import math
import sys
import pprint

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

# pprint.pprint(make_matrix_small(H,200))
# print("\n")
# pprint.pprint(H)

H = get_matrix(sys.argv[1])
N = int(len(H))
save = []

print("# off diagonal \t value")
for k in range(0, N):
    diag_sum = 0.0
    for i in range(0, N-k):
        j = i + k
        diag_sum += math.fabs(H[i][j])
    save.append(diag_sum/(i+1))   # i+1 = number of elements processed

max_value = max(save)
for i in range(len(save)):
    print("%g\t%g" % (i, save[i]/max_value))
