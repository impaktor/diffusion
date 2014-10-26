#include <iostream>    //for std::cout, among others
#include <cstdlib>     //for abort()-function
#include <cmath>       //gives us sqrt, pow, fabs,
#include <string>
#include <vector>
#include <omp.h>
//#include <thread>      //new C++ standard! same as old boost_threads
#include <map>

#include <sstream>     //for [o]stringstream

//Random number generator from numerical recipes 3 ed, but
//modified not to depend on nr3.h or else threads wont compile.
#include "nr/ran_mod.h"

#include "simpleini/SimpleIni.h" //parse ini-file
#include <stdlib.h>    //atoi

#include "auxiliary.h" //non-physics stuff. (print messages etc.)
#include "classes.h"   //various data structures /classes (Jump, Particle)
#include "lattice.h"   //this is the main class, that does the physics.
#include "save.h"      //class to save MSD, compute errors, and print to file


const std::string getWaitTime(enum waitingtime wt)
{
  switch (wt)
  {
  case LIN: return "lin";
  case EXP: return "exp";
  case POW: return "pow";
  }
  return "Error";
}



int main(int argc, char* argv[]){

  InputValues def;                     //initiate default paramter values
  aux::argumentFlags(argc, argv, def); //read in command line options


  //read input *.ini-file (from command line) and set everything that
  //wasn't set from the command line:
  //=================================================================

  const char * input_file = def.inputFileName.c_str();

  bool isUTF8 = true;

  //only one value for each key, on one line each:
  bool useMultiKey = false, useMultiLine = false;

  // load from a data file
  CSimpleIniA ini(isUTF8, useMultiKey, useMultiLine);
  SI_Error rc = ini.LoadFile(input_file);
  if (rc < 0){
    aux::printHelp(argv);
    aux::printError("invalid inputfile path: " + std::string(input_file));
  }

  //start parsing "key = value" from input-file with section "lattice"
  const char * tmpValue = ini.GetValue("lattice", "fix_boundary", NULL /* default*/);
  bool isBoundaryFix = aux::convertToBool(tmpValue, "fix_boundary");

  int X = atoi(ini.GetValue("lattice", "x", NULL));
  int Y = atoi(ini.GetValue("lattice", "y", NULL));
  int Z = atoi(ini.GetValue("lattice", "z", NULL));

  int nParticles = atoi(ini.GetValue("particles", "N", NULL));

  //jumprate of tracer particle, k_t
  tmpValue = ini.GetValue("particles", "rk", NULL);
  float jumprateTracer = (def.jmpTracer.second) ? def.jmpTracer.first : atof(tmpValue);

  tmpValue = ini.GetValue("time", "t_stop", NULL);
  double maxTime = atof(tmpValue);

  tmpValue = ini.GetValue("time", "waiting", NULL);
  waitingtime waitingTime = static_cast<waitingtime>(atoi(tmpValue));

  tmpValue = ini.GetValue("time", "N", NULL);
  int nSamplings = atoi(tmpValue);

  //use logarithmic spacing of sampling points
  tmpValue = ini.GetValue("time", "log_spacing", NULL);
  bool isLogScale = aux::convertToBool(tmpValue, "logspace");

  tmpValue = ini.GetValue("simulation", "M", NULL);
  int ensembles = atoi(tmpValue);

  float seed = atof(ini.GetValue("simulation", "seed", NULL));
  //=============================================

  //Construct vector of sampling times
  std::vector<double> samplingTimes;
  if(!isLogScale){
    double minTime = 1;
    double deltaTime = (double) (maxTime - minTime) / nSamplings;
    double timeSum = minTime; //start time
    for (int i = 0; i < nSamplings; i++){
      samplingTimes.push_back(timeSum);
      timeSum += deltaTime;
    }
  }
  else{
    //NOTE should take difference log(stop) - log(1)
    double deltaTimeLog = (double) log(maxTime)/nSamplings;
    double logTimeSum = log(1);
    while (logTimeSum < log(maxTime)){
      logTimeSum += deltaTimeLog;
      samplingTimes.push_back( exp(logTimeSum) );
    }
  }


  //save jump rate u to right and left direction
  Jump::Direction dir(jumprateTracer,jumprateTracer);

  Jump jmp;
  //same jumprate in all dimensions (x,y,z)
  jmp.x = dir;
  jmp.y = dir;
  jmp.z = dir;
  std::vector<Jump> jumpRates;
  jumpRates.push_back(jmp);

  //initiate the save-class to store, save, and print result
  //(also computes standard error and correlation estimate)
  Save save(samplingTimes, ensembles, def.isLowMem.first);

  //non-important nice-to-have simulation info.
  RemainingTime printToScreen(ensembles);

  std::vector<Lattice> lattices;   //possibility to store many lattice obj.
  int noLattice = 1;               //default to one Lattice object...
  int lattice_index = 0;           //... stored in first element of vector
# pragma omp parallel shared(save,lattices,noLattice) private(lattice_index)
  {
#   ifdef _OPENMP                      //If compiling with OpenMP...
    noLattice = omp_get_num_threads(); //change to threads" number of lattices.
    std::cerr << "Notice: I don't know if this 'brute' implementation is thread safe!" << std::endl;
#   endif

#   pragma omp single                  //only do this once in total.
    for(int i = 0; i < noLattice; ++i){
      //Initiate the Lattice class, which is what does all the physics.

      Lattice tmp(X, Y, Z, nParticles, seed * (i+1), isBoundaryFix);
      tmp.setJumpRate(jumpRates);

      //if use Interaction algorithm, with InteractStr
      if(def.isInteracting)
        tmp.setInteraction(def.interactionStrength);

      lattices.push_back(tmp);
    }

#   ifdef _OPENMP                         //if compiling with OpenMP:
    lattice_index = omp_get_thread_num(); //each thread operates on its own
#   endif                                 //lattice index in lattices-vector

//for(size_t i = 0; i < samplingTimes.size(); ++i) std::cout << "SAMPLINGTIMES\t" << samplingTimes[i] << std::endl;

#   pragma omp for
    for(int E = 0; E < ensembles; ++E){
      if(lattice_index == 0)               //only first thread prints progress
        printToScreen.printProgress(E*noLattice);

      //Store tracer position at each point for this "ensembles", needed by
      //class to compute standard error/deviation. (+for binning, etc.)
      //extract displacement coordinates of the tracer...
      //store the current ensemble values in these:
      std::vector<int> dx,dy,dz;
      std::vector<double> dr;

      for(size_t k = 0; k < samplingTimes.size(); ++k){
        std::vector<double> subSamplingTimes = samplingTimes;

//std::cout << "k\t" << k << "\t" << E << std::endl;

        // simulate to t=k
        subSamplingTimes.resize(k+1);

//for(size_t i = 0; i < subSamplingTimes.size(); ++i) std::cout << "\t" << subSamplingTimes[i]; std::cout << std::endl;

        lattices[lattice_index].setSamplingTimes(subSamplingTimes, waitingTime);
//std::cout << "subSamplingTimes.back()\t" << subSamplingTimes.back() << std::endl;

        lattices[lattice_index].place();
        lattices[lattice_index].move();

        // take last element in x, y, z
        std::vector<int> dx_sub,dy_sub,dz_sub;
        std::vector<double> dr_sub;
        lattices[lattice_index].getDisplacement(dx_sub, dy_sub, dz_sub, dr_sub);

        // push back to trajectories to save
        dx.push_back(dx_sub.back());
        dy.push_back(dy_sub.back());
        dz.push_back(dz_sub.back());
        dr.push_back(dr_sub.back());
      }

#     pragma omp critical  //only one thread may write at a time:
      {
        save.store(dx,dy,dz,dr);
      }
    }
  }

  //To print which distribution we used to head of output-file
  std::string dist[]     = {"uniform","exponential","powerlaw","nakazato"};
  std::string onOff[]    = {"Off","On"};
  std::string bound[]    = {"periodic","fix"};
  float d_eff = lattices[0].computeEffectiveDiffusionConst();
  float d_av = lattices[0].computeAverageDiffusionConst();

  //print three lines of info about simulation to head of each outfile:
  std::ostringstream print;
  print << "#E = " << ensembles <<"\t N = " << nParticles << "\t X-Y-Z: " //line1
        << X << "x" << Y << "x" << Z << "\t 2*d*D_naka: "
        << lattices[0].computeNakazato() << "  Waiting time: "
        << getWaitTime(waitingTime) << std::endl;
  print << "#Conc.: "<< (float) nParticles/(X*Y*Z) << "\t MSD_equil: "  //line2
        << "bound: " << bound[isBoundaryFix] << std::endl;
  print << "#D_eff: " << d_eff << "\t k_tagg: "                //line3
        << jumpRates[0].x.r << "\t D_av: " << d_av
        << "\t Interaction: " << onOff[def.isInteracting] << " "
        << def.interactionStrength << std::endl;

  std::string head = print.str();

  // dump raw trajectories (non-mean) to file. Do rest in script.
  save.dump(def.outputFileName, head);

  return 0;
}
