#!/usr/bin/env python2

#program takes an matrix as input (writen for B-matrix) and makes it
#into an H-matrix, by summing up the elements. Purpose of this script
#is to filer out elements that are small (> 1% of main diagonal).

import sys
import numpy as np

def get_matrix(file_name):
    """Read in all non commented lines and return columns as a single
    numpy array matrix
    """
    matrix = []
    myfile = open(file_name,'r')
    for line in myfile.readlines():
        li = line.strip()             #remove leading whitespace
        if not li.startswith("#"):    #ignore comments
            tmp = [float(value) for value in line.split()]
            matrix.append(tmp)
    myfile.close()
    return np.array(matrix)           #convert list -> numpy array

def print_usage():
    msg = """Usage:\t%s B-matrix_x B-matrix_y\n """ % (sys.argv[0])
    sys.stderr.write(msg)
    sys.exit()

if len(sys.argv) != 3:
    print_usage()

B_x = get_matrix(sys.argv[1])
B_y = get_matrix(sys.argv[2])
noSamplingTimes = len(B_x)
H = np.zeros(shape=(noSamplingTimes,noSamplingTimes))

max_x_element = np.amax(B_x)
min_x_element = np.amin(B_x)
max_y_element = np.amax(B_y)
min_y_element = np.amin(B_y)

# print("maxX:%g, minX:%g" % (max_x_element,min_x_element))
# print("maxY:%g, minY:%g" % (max_y_element,min_y_element))

invM = 1.0/250000

for i in range(noSamplingTimes):
    sys.stderr.write("%g\r" % (100.0*i/noSamplingTimes))
    for j in range(noSamplingTimes):
          sumX = 0
          sumY = 0
          for k in range(i+1):           # k <= i
              for l in range(j+1):       # l <= i
                  if B_x[k,l] > 0.01* max_x_element:
                      sumX += B_x[k,l]
                  if B_y[k,l] > 0.01* max_y_element:
                      sumY += B_y[k,l]
          H[i,j] = 2 * invM * (np.power(sumX, 2) + np.power(sumY, 2))




for i in range(len(H)):
    for j in range(len(H)):
        print("%g\t" % H[i,j])
    print("\n")
