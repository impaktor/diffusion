TODO:
Auxiliary: isComment
Classes: windingnumber (speedgain?)
Save: Reserve size for the matrix
Lattice: short vs. int ?
Fixa main_*.cpp
no->n
search: XXX & TODO

INSTRUCTIONS AND COMMENTS
=========================

Getting started with the program
--------------------------------
Upon execution, the program will ask for input parameters, if no
input-file was specified, (use the '-r' flag). Example:
  ./prog -w output.dat -r input.dat
  ./prog -w output_42.dat < input.dat
See below for more options. (or do ./prog -h)


Running
-------
For help on command flags and arguments just run the program with the
-h flag. For convenience they are given below:

-h        Displays a help.
-m        Use the (much) less memory consuming algorithm. Advisable if
          (Number of values to store) x (number of ensembles) / 2^30 ~RAM (in Gb)
-l        Use logarithmic spacing of the data points (time sampling points)
-i M      The -i option activates the interaction-algorithm, with "M"
          being the interaction strength. Note that M >= 0. If M = 0
          the simulation gives the same result as without, but it will run slower.
-w target Write simulation output to 'target' instead of the default "out.dat"
          (can include a path and name up to 99 characters long).
-r source Read an input-file instead of asking for simulation parameters.
-q        Run in "quiet mode", i.e. less printout to screen.
-B N      Generate N different output-files with the same input parameters,
          by using bootstrap; after the simulation it generates N synthetic output-files.
-b N      Generate N different output-files with the same input parameters,
          by using Brute force, i.e. Just rerun the simulation N times.
-j        Turn on the jackknife algorithm to print out a corrected error to the screen
          and file appended with "_jackknife". can not be used with -m flag.
-t k      Set the tracer jumprate to "k", where "k" is a float.


Output-files
------------
The primary output-file, (default name is out.dat) has the following
form:

    <time>  <dr^2>  <dr^2_error>    <4th column>,

where the "4th column" is either the Pearson coefficient, which
determines the correlation between the current value and the next; or
it is "z" (H*z=t) needed to compute the MSD (<dr^2>) and error, and
does not require us to save the entire H-matrix. This is set from the
Save::computeCorrelation()-function, which can also print the entire
H-matrix to a separate file ("name"_matrix), where "name" is the name
of the primary out-file (default: out.dat), but this is primarily for
testing purpose.


Getting started with the code
-----------------------------
If the daunting task of reading the source code has been bestowed upon
you, I'd recommend limiting yourself to lattice.[h,cpp] and
main.cpp. Do not venture in to my own personal playground that is
auxiliary.cpp. Also save.[h,cpp] is redundant for understanding the
physics that is simulated, and understanding the definitions in
classes.[h,cpp] is crucial.


Running on many computers
-------------------------
 use ssh:
      $ ssh -n 'cd /to/path; nice -19 ./prog < input.dat'
 Or:  $ ssh -X nice -19 ./prog output.dat < input.dat'
Use che2 to find available computers.


Files
-----
Makefile        Needed to compile the project. Contains description to compile
                the binary *.o files, and on which files they depend upon.

Makedep.rule    Is created by running "make depend", and will read all include-
                macros and contains which files depend on which.

main.cpp        main function. Most of the code is for finding input parameters,
                and setting the jumprate and sampling time vector.

main_thread.cpp Same as main.cpp, but uses threads, and depends on libboost_thread
                and does not use nr/nr3.h, but rather a modified nr/ran_mod.h, or
                else boost will not compile. The random number sequence is the same.

nr/             Folder where all the needed header files taken from Numerical
                Recipes (W.Press). This can be replaced with a symbolic link
                to the actual folder where all NR-files are located.

save.[h,cpp]    Save routine. During simulation it stores the MSD for each
                individual "ensemble", and then computes standard deviation,
                correlation, etc, and prints results to a file.

lattice.[h,cpp] This is where the physics happens.

classes.[h,cpp] Defines the Particle class which depends on Jump class,
                which depends on the Direction class. (in use since version 8).
                Also has a nice class for printing time left on a simulation.

global.h        Some hardcoded variables that I use when testing/checking the code.

auxiliary.[h,cpp]   Personal playground with "nifty (unnecessary?) shit".
                    The most important content is the printError-function,
                    Disregard the rest.

Compiling and the Makefile
--------------------------
If the file Makedep.rule does not exist, or any file has changed its
#include-macros, run "make depend". Then doing a simple "make"
compiles it. Check the header of the Makefile for options. Do not use
optimization flags ( -O2 or -O3) if you intend to debug the code
(compile with -g flag). The Makefile can take a number of different
commands by running the following from the command line:

make           Compile the binary "prog" from main.cpp
make thread    Compile using main_thread.cpp. Gives 4 threads.
make super     Compile the binary "super" from main_super.cpp which is the
               "super interaction"-code, currently in testing stage.
make tar       Create prog.tar.gz with all needed source files
make zip       Create prog.zip with all needed source files.
make clean     Remove all binary and object files (*.o)
make tags      Create a Tags-file for Emacs. Use M-. in emacs to navigate the files.
make lines     Returns the total number of lines of code. ( :D )

The Makefile will recompile any files that have been changed and all
files that depend on it, thus speeding up the compilation, and always
making sure the binary is correctly linked. Running "make -j4 -k"
compiles on 4 threads, and keeps going even when encountering errors.


Hard-coded variables:
--------------------
To not use winding number when using periodic boundary conditions,
check the commented code in void Lattice::move(). With the use of
winding numbers we keep track on how many times the tracer has circled
the lattice.

Generally, hard-coded variables are defined in the head of the
implementation (*.cpp) files. Especially save.cpp and main.cpp.


Dimension
---------
The program accurately accounts for all dimensions larger than 1. This
means that a lattice with size 100x100x0 or 100x100x1 is considered
two dimensional lattice, and only moves in X and Y are made. All
particles still have a Z-coordinate, but no move will be made in this
direction as long as Z<2, or else we would "leak" some
probability/speed of the particles.


Jump-rates
----------
The code allows for different jump rates in different directions,
although this has never been used.


Coding Conventions
------------------
The code follows the following rules:
(see http://geosoft.no/development/cppstyle.html)
 _______________________________________________________
|DESCRIPTION               | EXAMPLE                    |
|==========================|============================|
|All variables:            | camelCase, myVariable      |
|--------------------------|----------------------------|
|Constants:                | MY_CONSTANT                |
|--------------------------|----------------------------|
|Private class variables:  | camelCase_                 |
|--------------------------|----------------------------|
|Types:                    | CamelCase, MyType, Particle|
|(always a noun, or is it?)| Direction, etc.            |
|--------------------------|----------------------------|
|Functions & methods:      | getVariable(), computeD(), |
|(always a verb, ...is it?)| setVariable()              |
'_______________________________________________________'

-"Generic variables should have the same name as their type"

-operators are preceded and followed by space, " == ", " + " etc.

-Note: underscore "_" can be used in the standard LaTeX-sense: pos_0 etc.

use "int&" rather than "int &"


Additional speed
----------------
When running huge simulations one could try commenting some of the
testing if-statements that are there to find bugs. The most used
functions are:

  Lattice::moveAndBoundaryCheck
  Lattice::vacancyCheck
  Lattice::convertMuToParticle

Also investigate convertMuToParticle for changing to binary search
method, which I'm kind of already doing. (I split the interval in 3
then in 2, I think...)

Also it is a good idea to compile with -D NDEBUG (see Makeflile) when
running real (non-testing) simulations, as this disables any
asserts(). Also, if running a simulation for long time, such as is
needed for percolation, use lowMem = true. This saves a lot of RAM,
but unfortunately we can not use bootstrap since that needs all
individual trajectories after the simulation.


CHANGE LOG
==========

NOTE: think about the k_t_ in save-class!

version 11 2011-11
----------
Finishing it up for Tobias:
-use of ?-operator in moveAndBoundaryCheck()
-redid the Makefile completely.
-New c++ standard here! Includes threads from boost library, i.e. main_thread.cpp no longer depends on the boots_thread library!
-introduced several functions as "const".
-removed Lattice:buildCluster2(). I don't know what it did there?.
-made vacancy-matrix use int (again) as short is too small for big systems.
-removed some {} from one-line if-statements
-inserted a number of assert() to not need the __LINE__ macro (in superInteraction).
-minor correction to nr/fitab and nr/LUdecomp to remove compilation warnings.
-moved global debug/testing-variables to global.h
-introduced the InputValues-struct to make the reading of command line arg. easier.
-restructured the computeCorrelation to two functions. One for pearson and one for z.
-changed name of Save::saveBinning --> computeDistribution, and moved it from save::save() to be called from main().


version 10 2010
----------
-finally made my superInteraction-code compile. (note to self: "inline" in base-class, not a good idea).
-cleaning up code, removed a bunch of no longer used functions in main.cpp, some of them are still defined in lattice.[h,cpp] though for use again.
-changed vacancy matrix from <int> (4 byte) to <short> (2 byte), only uses half as much RAM!
-added header to each column in "out.dat_txyz"-file.
-significant recoding of auxiliary.cpp. Re-wrote getNonCommentInput completely, and it now works, which it did not before.
-added a '-r' flag for reading input files (with comments!), I bet Tobias will be happy now...
-added copy constructors to Particle:: and Direction:: (even though I'm not using them), and made the void constructor initiate everything to 0, for both Particle, Jump, and Direction -classes.
-removed:
    Lattice::vacancyCheckOld
    Lattice::moveOlod
    Lattice::snapshot,
    Lattice::dumpSimulation,
    Lattice::saveCluster & Lattice::printCluster (meant to print clustersize-distribution when using interaction),
    TESTinteraction.cpp (removed entire file, old code. Not used.)
-minor changes in Lattice:checkVacancyMatrix
-removed 3 redundant '#include <...>' in lattice.cpp, and wrote an '#ifndef' in nr/ran_mod.h.

version 9
----------
-implemented "windingnumbers", for use with periodic boundary. One could have it as it is now, with a class-reachable windingnumber for just the tracer particle, alternatively one could add this to the particle class, which would be less hassle, (no need to fiddle about in Lattice::vacancyCheck), only in Lattice:MoveAndBoundaryCheck(), but would require (slightly) more memory. To turn on/off, just (un)/comment the (indicated) code in Lattice::Move.
-(23/2 -11) To make the Lattice class be able to multi-thread, setting the Ran class (NR 3 ed.) as "static" is not optimal, since we will initiate several instances of the Lattice-class, and static reaches across the class encapsulation and we don't know what happens if two threads try to access the random number generator at the same time. Therefore, we initiate a single non-static random number generator in the class, that is only reachable for the members of that instance of the class. We do this by using initialization list in the class. This makes my simulation non-comparable with previous simulations, since I now only use one seed for each Lattice class instead of three (waiting time, place, and move).
-Made the bootknife in save.cpp be able to make several fittings with different starting times.
-All previous simulations ever made had the exact same distribution of crowders, i.e. "quenched", now I WILL/CAN change this [EDIT: first I must get Nakazato to give reasonable results], and start by moving the generation of jumprates to a separate function in the main() function.
-Added the 't' flag for setting the tracer jumprate, i.e. not needing to compile 22 different binary files. Distribution, and nakazato crowding rate (k_c) is still hard coded.
-Made the save-class show remaining time when bootstrapping and computing the H-matrix. This called for some minor changes in the RemainingTime-class and also made the save-class depend on classes.[cpp,h].
-my super-interaction code from this summer was probably broken many versions ago by modifications in various classes. Now it works, and is renamed (main_super.cpp), and is directly based on main.cpp
-Code now uses the Particle-class, which was implemented some time between version 8 & 9.


version 8 (15/2 2011)
---------------------
-Renamed "prog.cpp" to more appropriate "main.cpp".
-Made a threading version of main.cpp, called main_thread.cpp. It uses the
 thread-package in the boost library. Increase in speed more than two-fold.
-Re-made the Makefile, now it is glorious.
-Changed the search-function (slightly) in Lattice::convertMuToParticle()
-Moved all the save-routines to their own class, save.[h,cpp],
 and compensated for it in main.cpp and superProg.cpp
-Continued the renaming of variables according to mentioned
 convention above. (dx -> dx_, etc.)
-Removed all if-statements in main.cpp that allowed for using
 Lattice::moveOld()-function. (I never use it)
-Converted all "char*, char name[]" types to std::strings.
-Output-files now stores which waiting time distribution was used.
-Made the binning function in save.cpp work.
-Added a bootstrapping function in the save-class, (with b flag)
-Added a Jackknife function in the save-class, and a command flag "-j".
-Added a "hybrid" bootknife function in save::computeBootstrap. It is
 hard coded and turned on/off by boolean at head of save.cpp.
-Put an additional loop outside the ensemble loop, so now I can
 run many simulations in a row (with the same parameters), to generate
 numerous output-files, which I need to check the correlation between
 sampling points. Also made two new command flags, b (use bootstrap), and
 B (use Brute force), both requires one to specify the number of re-runs.
-Introduced two global variables in save.cpp, for easy hard coding
 access to turn on/off of correlation function, and binning function.

version 7 (Nov 2010)
---------------------
-Renamed and split up my main (singleton) class to "Lattice"
(lattice.[h,cpp])
-Made my TESTinteraction.cpp (generally called "super interaction")
into a class SuperInteraction, that inherits all the public and protected
methods/variables of the Lattice class. This makes it use the same source
code as the lattice class, plus some specific functions. I'm keeping the
old TESTinteraction.cpp for documentation, since it is completely self-
contained, (ie. everything in one source file, except ran.h & nr3.h).
-Introduced the H-matrix and computeCorrelation-function, to compute
the correlation between sampling points, and use an H-matrix that
automagically removes/compensates for correlation. (Work on this
 continues in the next version, as it is not 100% complete)
-Use of svd.h evaluates the error in H-matrix due to floating point.
-made setSamplingTime take an additional optional argument for
exponential waiting time.
-Moved all the numerical recipes codes to their own folder.
-The three main Seeds will now be set through define-macros.

version 6
---------
NOTE: from this version forward the output is not comparable to
previous simulations, by use of the "diff"-command, since I now have
a new algorithm to set the sampling times! The simulation result is
still the same, physically, but not truly identical.

-Changed the calculation of average diffusion constant to omit the
tracer particle.
-Activated the LowMem command flag ('-m')
-Resolved the spaghetti that was my sampling time code. In doing so,
removed get_[dx,dy,dz,dr], set_dt, etc, and other garbage
-Renamed the matrices dx,dy,dz,dr to store_dx, store_dy. etc
-New single row vector called dx,dy,dz,dr, to store the displacement
for a single "ensemble". These are then stored in store_dx, store_dy, etc.
-Made a template of the getNonCommentInput function, and removed the
two functions specific for int/double return value. (in auxiliary.cpp)
-Created the computeWaitingTime()-method, and encapsulated the move()
function in a while(time is still running) condition.
-Changed the logarithmic spacing between the data points, to be more
like the algorithm in the Gillespie-source code.
-Calculate the pearsonCoefficient, and print it to the output as a
fourth column. This is work in progress. Not complete. This also
depends on a LU-decomposition and adds the dependency to ludcmp.h


version 5:
---------
Note: this is the last version that gives output-files that are comparable
with all previous simulations used in my master thesis, and elsewhere.

-Removed the six jump-rate vectors, and the three overloaded setJumpRate-
functions. Implemented a new Jump class, with members .x.r, .x.l, .y.r, etc.

-I do have an auxiliary.h file


version 4:
---------
????


-removed the Xpos, Ypos, Zpos vectors and made the class Particle,
to store the position of the particles:
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
-I ran many test to see that nothing got messed up when I merged
them, but it was straight forward. only one if-statement, basically,
to check if we should use the interaction algorithm.


version 2
---------
-Last version with prog.cpp and progInteraction.cpp separate.


-Fixed serious bug in original Gillespie code for finding the right mu.
It got stuck in an infinite loop for large simulations.

-all versions see a continuous revision of comments and naming of
variables, as I have been consistently inconsistent in the naming
convention. (See README.txt for my new rules on how to do this)

