INSTRUCTIONS AND COMMENTS
=========================

Compiling
---------
Just running the command "make" in the source directory defaults to running the Makefile that's placed therein. One could also run 'g++ classes.cpp prog.cpp' in the command line. Optimization flags I usually use is -O2 or -O3. If you intend to debug the code (compile with -g flag), don't use optimization!

Running
-------
For help on command flags and arguments just run the program with the -h flag. For convenience they are given below:
-h         Prints help to terminal.
-m         Use the (much) less memory. Advisable if 
           (Number of values to store) x (number of ensembles) / 2^30 ~ RAM (in Gb) 
-l         Use logarithmic spacing of the data points
-i n       the -i option activates the interaction-algorithm, with "n"
           being the interaction strength. Note that n >= 0. With n = 0
           the simulation gives the same result as without, but it will run slower.
-w target  changes the default output: "out.dat" to "target", which can include 
           a path and name up to 99 characters long.

Dimension
---------
The program accurately accounts for all dimensions larger than 1. This means that a lattice with size 100x100x0 or 100x100x1 is considered two dimensional lattice, and only moves in X and Y are made. All particles still have a Z-coordinate, but no move will be made in this direction as long as Z<2, or else we would "leak" some probability/speed of the particles. 


Jump-rates
----------
The code has been written so if anyone is so inclined they can have different jump rates in different directions, to allow simulations during some biased/potential influence on the particles. 

setting jump rates in three dimensions if only 2D lattice will ... XXX-TODO-XXX!

Coding Conventions
------------------
The code follows the following rules:
 ______________________________________________________
|DESCRIPTION              | EXAMPLE                    |
|                         |                            | 
|-------------------------|----------------------------| 
|All variables:           | camelCase, myVariable      | 
|-------------------------|----------------------------| 
|Constants:               | MY_CONSTANT                | 
|-------------------------|----------------------------| 
|Private class variables: | camelCase_                 | 
|-------------------------|----------------------------| 
|Types:                   | CamelCase, MyType, Particle| 
|(always a noun)          | Direction, etc.            | 
|-------------------------|----------------------------| 
|Functions & methods:     | myFunction(), getVariable()| 
|(always a verb)          | setVariable(), computeD()  | 
'______________________________________________________'

-"Generic variables should have the same name as their type"

-operators are preceded and followed by space, " == ", " + " etc.

-Note: underscore "_" can be used in the standard LaTeX-sense: pos_0 etc.  

use "int&" rather than "int &"


