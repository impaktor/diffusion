INSTRUCTIONS AND COMMENTS
=========================

Files
-----
nr/      Folder where all the needed header files taken from Numerical
         Recipes (W.Press). This can be replaced with a symbolic link
         to the actual folder where all NR-files are located. 
         
         


Compiling
---------
Just running the command "make" in the source directory defaults to running the Makefile that's placed therein. One could also run 'g++ classes.cpp prog.cpp' in the command line. Optimization flags I usually use is -O2 or -O3. If you intend to debug the code (compile with -g flag), don't use optimization!

Running
-------
For help on command flags and arguments just run the program with the -h flag. For convenience they are given below:
-h         Prints help to terminal.
-m         Use the (much) less memory algorithm. Advisable if 
           (Number of values to store) x (number of ensembles) / 2^30 ~ RAM (in Gb) 
-l         Use logarithmic spacing of the data points
-i n       the -i option activates the interaction-algorithm, with "n"
           being the interaction strength. Note that n >= 0. If n = 0
           the simulation gives the same result as without, but it will run slower.
-w target  changes the default output: "out.dat" to "target", which can include 
           a path and name up to 99 characters long.
-q         Run in "quiet mode", i.e. no printout to screen.

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



CHANGE LOG
==========

version 7 (current prog.cpp)
----------------------------
-Introduced the H-matrix and computeCorrelation-function, to compute
the correlation between sampling points, and use an H-matrix that
automagically removes/compensates for correlation. (Work on this
 continues in the next version, as it is not 100% complete)
-Use of svd.h evaluates the error in H-matrix due to floating point.
-renamed and split up my main (singleton) class to "Lattice"
(lattice.[h,cpp])
-made setSamplingTime take an additional optional argument for
exponential waiting time.
-made my TESTinteraction.cpp (generally called "super interaction")
into a class SuperInteraction, that inherits all the public and protected
methods/variables of the Lattice class. This makes it use the same source
code as the lattice class, plus some specific functions.  I'm keeping the
old TESTinteraction.cpp for documentation, since it is completely self-
contained, (ie. everything in one source file, except ran.h & nr3.h).
-Moved all the numerical recipes codes to their own folder.
-The three main Seeds will now be set through define-macros.

version 6 
---------
NOTE: from this version forward the output is not comparable to 
previous simulations, by use of the "diff"-command, since I now have
a new algorithm to set the sampling times! The simulation result is 
still the same, physically, but not truly identical.

-changed the calculation of average diffusion constant to omit the
tagged particle.
-activated the LowMem command flag ('-m')
-resolved the spaghetti that was my sampling time code. In doing so,
removed get_[dx,dy,dz,dr], set_dt, etc, and other garbage
-renamed the matrices dx,dy,dz,dr to store_dx, store_dy. etc
-new single row vector called dx,dy,dz,dr, to store the displacement
for an ensemble. these are then stored in store_dx, store_dy, etc.
-made a template of the getNonCommentInput function, and removed the 
two functions specific for int/double return value. (in auxiliary.cpp) 
-created the computeWaitingTime()-method, and encapsulated the move()
function in a while(time is still running) condition.
-changed the logarithmic spacing between the data points, to be more 
like the algorithm in the Gillespie-source code. 
-Calculate the pearsonCoefficient, and print it to the output as a 
fourth column. This is work in progress. Not complete. This also 
depends on a LU-decomposition and adds the dependency to ludcmp.h 


version 5:
---------
Note: this is the last version that gives output-files that are comparable 
with all previous simulations used in my master thesis, and elsewhere.  

-removed the six jump-rate vectors, and the three overloaded setJumpRate-
functions. Implemented a new Jump class, with members .x.r, .x.l, .y.r, etc.

-I do have an auxiliary.h file
 

version 4:
---------
????


-removed the Xpos, Ypos, Zpos vectors and made the class Particle, to
store the position of the particles:
(Particle pos[i].x pos[i].y was Xpos[i], Ypos[i] etc.) 



version 3
---------
Just stating what it looks like:
-all vectors normal int/double, ie. 3 for pos. + 6 for jumprate 
-Uses auxiliary.cpp (no h-file)
-has quiet-flag
-has LowMem
-interaction + (normal) hard-core interaction merged

version 2.5
-----------
-merged the two separate codes prog.cpp + progInteraction.cpp
into one single file. progInteraction.cpp was the one I previously 
used for particles with an attractive potential. (see Master thesis.)
-I ran many test to see that nothing got messed up when I merged them,
but it was straight forward. only one if-statement, basically, to 
check if we should use the interaction algorithm.   


version 2
---------
-Last version with prog.cpp and progInteraction.cpp separate. 


-Fixed serious bug in original Gillespie code for finding the right mu. 
It got stuck in an infinite loop for large simulations.

-all versions see a continuous revision of comments and naming of 
variables, as I have been consistently inconsistent in the naming convention.
(See README.txt for my new rules on how to do this)

