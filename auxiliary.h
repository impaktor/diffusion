#ifndef AUXILARY_H
#define AUXILARY_H

#include <string>

//This file contains implementations that have nothing to do with
//physics, such as printing messages, ask for input parameters, read
//command arguments (flags), print error messages, ignore commented
//lines, etc.


//container for the default simulations values.
//modifiable by passing arguments on command line.
struct InputValues{
public:
  //use logarithmic spacing of sampling points
  bool isLogScale;

  //run in low memory mode
  bool isLowMem;

  //use nearest neighbor interaction
  bool isInteracting;

  //print simulation information to screen
  bool isQuiet;

  //if  isInteracting = true, then use:
  int interactionStrength;

  //save simulation output to files with this prefix
  std::string outputFileName;

  // file (path) to read i input values from
  std::string inputFileName;

  //number of output files/MSD-trajectories.
  //Only valid for "bootstrap" or "Brute force"
  int nOutputs;

  //use Bootstrap of Bruteforce?
  char method;

  //use jack knife method?
  bool isJackknife;

  //jumprate of tracer particle
  float jmpTracer;
};


//only these functions needs to be accessible from the "outside",
void argumentFlags(int, char**, InputValues& input);

//parameters to ask for during start (if no input-file was given)
void askUserForInputParameters(int&,int&,int&,int&,int&,int&,double&);

//get simulation parameters from input-file instead.
bool readInputFile(std::string&,int&,int&,int&,int&,int&,int&,double&);

//Prints message, and aborts the program.
void printError(std::string message);
void printError(std::string message, int line);
void printError(std::string message, std::string file, int line);

//only used by AskUserForInputParameters():
void printHelp(std::string, char**);

#endif
