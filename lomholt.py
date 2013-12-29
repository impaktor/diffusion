#!/usr/bin/env python2
#This takes a file out_H-matrix and then uses the Lomholt procedure to
#get the slope and error. (This file is just lomholt_old2.m rewritten)

# We want to compute slope and sigma, from:
#   slope = (y * inv(H) * t) / (t * inv(H) * t) = (y * z) / (t * z)
#   sigma^2 = 1 / (t * inv(H) * t)= 1 / (t * z)
# where we defined z := inv(H) * t, and y = MSD, and "*" = dot product

import sys
import numpy as np
import bisect
import argparse

ERROR = "\033[1;31mERROR \033[0m"

def parse_arg():
    "Use argparse to parse the arguments from the command line"

    prog_descrip = """Apply lomholt correlated fit to output from
    simulation to get slope and error. Can also get as function of
    sampling points if H matrix provided"""

    parser = argparse.ArgumentParser(description=prog_descrip)

    parser.add_argument('-i','--input', help='Input file (with columns: t, MSD, stderr, pearson, z)', required=True)

    matrix_descr = """H-matrix, only needed if we want to remove sampling points in time,
    and/or check slope as function of included sampling points"""

    parser.add_argument('-m','--matrix', help=matrix_descr, required=False)
    parser.add_argument('--min-time', help='Earliest time to include in the fitting', required=False, default = 0)
    parser.add_argument('-N','--timepoints', help='Change size of H matrix to this. Must use -m as well.', required=False)

    args = parser.parse_args()
    return args


# Fore some nice colors in warning output
# http://stackoverflow.com/questions/287871/print-in-terminal-with-colors-using-python
class color:
    HEADER  = '\033[95m'
    OKBLUE  = '\033[94m'
    OKGREEN = '\033[92m'
    YELLOW  = '\033[93m'
    RED     = '\033[91m'
    ENDC    = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.YELLOW = ''
        self.RED = ''
        self.ENDC = ''

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


def remove_timepoints(first_time_to_use, time):
    """For early times tagged particle does not know its neighbors,
    remove these poins from fit."""

    #find index to insert first_time_to_use to maintain ordering in
    #time-list.
    first_index = bisect.bisect(time, float(first_time_to_use))
    return first_index


def prune_to_size_n(H, n, t, y):
    "Make H, t, y into size N(xN)."

    N, N = np.shape(H)         #get number of elements, returns tuple

    #construct a new smaller matrix h, with just n sampling points
    h = np.zeros((n,n))       #initiate new smaller matrix
    t_new = np.zeros(n)       #initiate the time-vector
    y_new = np.zeros(n)       #initiate the new MSD vector

    #construct new (smaller) H-matrix
    for i in range(len(t_new)):
        I = int(np.fix(i*N/n))
        for j in range(0,i+1):
            J = int(np.fix(j*N/n))
            h[i,j] = H[I,J]
            h[j,i] = h[i,j]
        t_new[i] = t[I]
        y_new[i] = y[I]
    z_new = np.dot(np.linalg.inv(h), t_new)
    return (t_new, y_new, z_new)


def compute_slope(t, y, z):     #compute MSD^T dot z / time^T dot z
    sigma = 1.0/np.sqrt(np.dot(t.transpose(), z))
    slope = np.dot(y, z) / np.dot(t, z)
    return (slope, sigma)


def prune_many(H, t, y):
    """Successively remove parts of matrix H to see how sigma and mu
    depend on number of sampling points (in time), see if it
    converges. Requires the correlation matrix H to be given.
    """
    R = 100                    #Number of different sampling times to try.
    N, N = np.shape(H)         #get number of elements, returns tuple
    if R > N:
        sys.stderr.write(color.RED + "Warning: " + color.ENDC +
                         "Matrix dim. (%g) < resolution (%g)\n" % (N,R))
        sys.stderr.write("Force-setting resolution to %g\n" % N)
        R = N
    step = np.fix((N/R))       #steplength to make R number of matrices

    for n in np.linspace(step, N, R):
        (t_new, y_new, z_new) = prune_to_size_n(H, n, t, y)
        (slope, sigma)        = compute_slope(t_new, y_new, z_new)
        print("%i \t %g \t %g" % (n, slope, sigma))


def main(args):
    data = get_matrix(args.input)

    s = remove_timepoints(args.min_time, data[:,0])

    #tmp = "%s %s" % (s, len(data[:,0]))
    #print >> sys.stderr, tmp

    if s >= len(data[:,0]):
        error_str= ERROR + "Cant remove more data points (%s) than there are (%s)." % (s, len(data[:,0]))
        #print(error_str, file = sys.stderr)   #Python3 syntax
        sys.stderr.write(error_str)            #Python2 and 3 syntax
        sys.exit()

    time = data[s:,0]
    MSD = data[s:,1]
    z = data[s:,4]    # (might recalculate this from H*z=t)

    if args.matrix != None:
        H = get_matrix(args.matrix)
        H = H[s:, s:]

        # TODO / XXX This Z works better (always) than the one from C++?!
        z = np.dot(np.linalg.inv(H), time)  #recompute z, will be same as data[s:,4].

        if args.timepoints != None:
            (time, MSD, z) = prune_to_size_n(H, int(args.timepoints), time, MSD)
        else:
            #If given an H-matrix, start removing sampling points, and
            #fit line, repetedly for smaller and smaller H:
            prune_many(H, time, MSD)

    (slope, sigma) = compute_slope(time, MSD, z)
    print("# mu: %g" % slope)
    print("# sigma: %g" % sigma)


if __name__ == "__main__":
    args = parse_arg()
    main(args)
