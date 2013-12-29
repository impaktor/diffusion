#!/usr/bin/env python2

# Run read_data on all inputfiles (filled with every single trajectory
# for one simulation), then run lomholt.py on the output (t, MSD, and
# Z) and H-matrices, and extract the slope and sigma, and output them.

import sys
import subprocess    # popen
import math

#Path to these programs, relative to where this script is run:
PATH_TO_read_data = "~/crowding/src/"
PATH_TO_lomholt = "~/crowding/src/"

# TODO: get environment for use_H_matrix

#Either compute slope from H-matrix, or just parse from
#jackknife-shrinkage-outputfile from read_data.cpp
use_H_matrix = False


ERROR = "\033[1;31mERROR \033[0m"

def run_read_data(file_to_read, M, min_time = 0):
    """get the mean data, and H-matrixes from big file with every
    single trajectory in"""

    def extract_slope_and_error(stream):
        "extract from stdout of lomholt.py"
        split_stream = stream.split("\n")
        if len(split_stream) >= 2:
            str1 = split_stream[0].split()
            str2 = split_stream[1].split()
            if str1[1] == 'mu:':
                mu = str1[-1]       #slope
            else:
                print >> sys.stderr, 'error in parsing mu: %s' % stream.split('\n')[0]
                mu = -1
            if str2[1] == 'sigma:':
                sigma = str2[-1]    #error of slope
            else:
                print >> sys.stderr, 'error in parsing sigma: %s' % stream.split('\n')[1]
                sigma = -1
        else:
            sigma = -2
            mu = -2
        return mu, sigma

    def extract_slope_and_error2(filename):
        "extract from the jackknife+shrinkage output of read_data"

        f = open(filename,'r')
        string = f.readline()
        string += f.readline()
        f.close()

        mu, sigma = extract_slope_and_error(string)
        return mu,sigma


    file_name_without_full_path = file_to_read.split("/")[-1]
    output_name = "processed_" + file_name_without_full_path
    command1 = PATH_TO_read_data + ("read_data %s %s" % (file_to_read, output_name))

    #make data of trajectories into MSD and H-matrix
    p1 = subprocess.Popen(command1, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    stdout, stderr = p1.communicate()
    if p1.poll():
        print(ERROR + "1: in popen: %s" % stdout)
    if stderr:
        print(ERROR + "1: %s" % stderr)

    new_N = math.sqrt(M/2.0) # don't use all points in the data set (reduce H-matrix)

    if use_H_matrix:         # either compute lomholt or just read result from jackknife+shrink-file
        #apply correlated chi**2 fitting on our newly generated output from command1 above
        command2 = (PATH_TO_lomholt + "lomholt.py -i %s --min-time %s -N %s" % (output_name, min_time, new_N))
        p2 = subprocess.Popen(command2, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
        stdout, stderr = p2.communicate()
        if p2.poll():
            print(ERROR + "2: in popen: %s" % stdout)
        if stderr:
            print(ERROR + "2: %s" % stderr)
        slope, error = extract_slope_and_error(stdout)
    else:
        slope, error = extract_slope_and_error2((output_name + "_shrink_jack_slope"))

    return slope, error


def extract_k_and_c_and_M(file_to_read):
    "extract the jumprate and concentration from the head of the file"

    f = open(file_to_read,'r')
    line = []
    for i in range(0,3):                             #get the first three lines
        line.append(f.readline().split('\t'))
    f.close()

    inv_k_in_parenthesis = line[1][2].split()[-1]
    tmp = inv_k_in_parenthesis
    for char in '()':
        tmp = tmp.replace(char, '')
    k = 1.0 / float(tmp)

    M = int(line[0][0].split('=')[-1])
    concentration = float(line[1][0].split('=')[-1]) #second line, first instance
    #k_tagg = float(line[2][1].split('=')[-1])        #third line, second instance
    return k, concentration, M



def asFuncOfConcentration(files_to_run):
    "Get correlation factor as function of concentration"
    data = {}

    for f in files_to_run:
        print("\033[1;32m Doing File:\033[0m\t" + f)
        k, c, M = extract_k_and_c_and_M(f)

        k_t = 1.0
        min_time = 10.0 / (min(k,k_t)*c)   #Min time to start fit from

        slope, error = run_read_data(f, M, min_time)

        one_output_line_tuple = (c, slope, error)

        if data.has_key(k):
            data[k].append(one_output_line_tuple)
        else:
            data[k] = []
            data[k].append(one_output_line_tuple)

    #print result to file [Correlation factor as function of jump rate k]
    for key,val in data.items():                  # loop over key, and values
        save_file_name = "k%s.dat" % (key)        # name and key is jump rate
        f = open(save_file_name, 'w')
        f.write("# c\t slope\t error\n")
        for i in sorted(val, key=lambda x: x[0]): # sorted: 10.0 commes after 1.0
            for j in i:
                f.write("%s\t" % j)
            f.write('\n')
        f.close()


def asFuncOfJumprate(files_to_run):
    "Get correlation factor as function of jump rate"

    data = {}

    for f in files_to_run:
        print("\033[1;32m Doing File:\033[0m\t" + f)
        k, c, M = extract_k_and_c_and_M(f)

        k_t = 1.0
        min_time = 10.0 / (min(k,k_t)*c)   #Min time to start fit from

        slope, error = run_read_data(f, M, min_time)

        one_output_line_tuple = (k, slope, error)

        if data.has_key(c):
            data[c].append(one_output_line_tuple)
        else:
            data[c] = []
            data[c].append(one_output_line_tuple)

    #print result to file [Correlation factor as function of jump rate k]
    for key,val in data.items():             # loop over key, and values
        save_file_name = "c%s.dat" % (key)   # name and key is concentration
        f = open(save_file_name, 'w')
        f.write("# k\t slope\t error\n")
        for i in sorted(val, key=lambda x: x[0]): # sorted: 10.0 commes after 1.0
            for j in i:
                f.write("%s\t" % j)
            f.write('\n')
        f.close()



def print_usage():
    print("usage:\n %s [k | c] [file1, file2, file3...]" % sys.argv[0])
    print("c\t\t if  processing files with different concentrations.")
    print("k\t\t if  processing files with different jumprates.")
    print("file1...\t files to process.")

def main(argv):
    if not len(argv) > 1:
        print_usage()
    else:
        files_to_run = argv[2::1]   # every element from 2. (skipp argv[0] & [1]) argv[1:]

        if(argv[1] == "k" or argv[1] == "c"):
        #Either get correlation as function of concentration, and each
        #output gets its name from the jump rate...
            if(argv[1] == "c"):
                print("asFuncOfConcentration")
                asFuncOfConcentration(files_to_run)
        #...or get correlation as function of jumprate k, and each
        #output gets its name from the concentration
            if(argv[1] == "k"):
                print("asFuncOfJumprate")
                asFuncOfJumprate(files_to_run)
        else:
            print_usage()

if __name__ == "__main__":
    main(sys.argv)
