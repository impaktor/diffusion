#!/usr/bin/env python2
# Process Pigon's raw data files, to search for bugs. Read in x,y,z
# coordinates for each trajectory, sum over the squaresverage

import sys
import argparse
import numpy

def parse_arg():
    "Use argparse to parse the arguments from the command line"

    prog_descrip = """Read in raw data with every trajectory printed,
    seperated by a line with: \"SEED = 000\""""

    parser = argparse.ArgumentParser(description=prog_descrip)

    parser.add_argument('-i','--input', help='Input file (with columns: t, x,y,z,r)', required=True)
    args = parser.parse_args()
    return args


# String -> ListOfFloat
# Open file, read in first trajectory and return list of sampling times
def get_time(filename):
    myfile = open(filename,'r')

    for i in range(5):            # skipp header comments
        tmp = myfile.readline()

    time = []
    for line in myfile.readlines():
        if line.startswith("#"):  # starts with "# SEED = " actually
            break
        tmp = [float(value) for value in line.split()] #make line float
        time.append(tmp[0])       # first column is time

    return time


# List, ListOfList, Integer -> nil
# pretty print data in matrix, with first List to the left
def print_data(time, matrix, M):

    for i in range(len(matrix)):
        MSD = (matrix[i][0] + matrix[i][1] + matrix[i][2]) / M
        print("%s\t %s\t %s" % (time[i], MSD, matrix[i][3] / M))


def main(args):
    myfile = open(args.input,'r')
    for i in range(4):
        tmp = myfile.readline()

    time = get_time(args.input)
    data = []
    for i in range(len(time)):
        data.append([0,0,0,0])  #dx, dy, dz, dr

    index = 0
    M = 0
    for line in myfile.readlines():
        if line.startswith("#"):
            index = 0    # reset index, new trajectory
            M = M + 1    # count which trajectory we're doing
            sys.stderr.write("#Doing trajectory: %s\n" % M)
            continue     # skipp this ("#SEED=...") line
        line_float = [float(value) for value in line.split()]

        for j in range(1,len(line_float)):
            data[index][j-1] += line_float[j]**2
            #print("index: %s,\t j: %s =\t %s" % (index, j, data[index][j]))

        index += 1

    myfile.close()
    print_data(time, data, M)



if __name__ == "__main__":
    args = parse_arg()
    main(args)
