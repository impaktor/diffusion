#ifndef AUXILARY_H
#define AUXILARY_H

#include <string>
#include <sstream>              // for ostringstream

//This file contains implementations that have nothing to do with
//physics, such as reading command arguments (flags), print error
//messages

// convert anything to string:
template<class T>
static std::string tostring(const T &t){
  std::ostringstream s;
  s << t;
  return s.str();
}


struct InputValues{
  //container for the default simulations values. Modifiable by
  //passing arguments on command line. The first in each std::pair
  //holds the value, and the second is whether it has been changed
  //through the command line.

public:

  std::pair<bool, bool> isLowMem;   //run in low memory mode

  bool isInteracting;               //use nearest neighbor interaction
  int interactionStrength;          //if isInteracting = true


  std::string outputFileName;   //save simulation output to files with this prefix
  std::string inputFileName;    //file (path) to read i input values from

  //number of output files/MSD-trajectories, only valid for "bootstrap"
  int nOutputs;

  bool isJackknife;   //use jackknife method?
  bool isBruteForce;  //use brute force, i.e. only use one point/trajectory
  char method;        //use Bootstrap ?
  bool isBootknife;   //use Bootknife method?

  //number of trajectories to generate when/if bootknifing.
  size_t nBootknife;

  //jumprate of tracer particle
  std::pair<float, bool> jmpTracer;

  InputValues(){
    isLowMem =            std::make_pair(false, false);
    isInteracting =       false;
    interactionStrength = 0;
    outputFileName =      "out.dat";
    inputFileName =       "input.ini";
    nOutputs =            1;
    method =              '-';
    isBootknife =         false;
    nBootknife =          1;    //what is this? XXX TODO
    isJackknife =         false;
    isBruteForce =        false;
    jmpTracer =           std::make_pair(1.0, false);
  }
};


namespace aux{
  void argumentFlags(int argc, char** argv, InputValues &in);

  bool convertToBool(const char * c, std::string);

  //Prints message, and aborts the program.
  void printHelp(char**);
}

#endif
