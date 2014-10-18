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


void computeJumpRates(std::vector<Jump>&, float&, int, float, float, float, int);

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

  tmpValue = ini.GetValue("particles", "jumprate_distr", NULL);
  const int jumprateDistribution  = atoi(tmpValue);

  //jumprate of tracer particle, k_t
  tmpValue = ini.GetValue("particles", "rk", NULL);
  float jumprateTracer = (def.jmpTracer.second) ? def.jmpTracer.first : atof(tmpValue);

  tmpValue = ini.GetValue("particles", "rc", NULL);
  float jumprateCrowders = atof(tmpValue);

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
  float seedJump = atof(ini.GetValue("simulation", "seed_jump", NULL));

  tmpValue = ini.GetValue("save", "print_histogram", NULL);
  bool printHistogram = aux::convertToBool(tmpValue, "print_histogram");

  tmpValue = ini.GetValue("save", "noHistogramBins", NULL);
  int noHistogramBins = atoi(tmpValue);

  //=============================================


  //When using the interaction code, make sure I have Nakazato distribution
  //with identical jumprates!
  if ((def.isInteracting && jumprateDistribution != 3) ||
      (jumprateCrowders != jumprateTracer
       && jumprateDistribution == 3 && def.isInteracting)){
    std::cout << "\nWARNING!! using interaction-code but"
              << " not identical jumprates! \a\n" << std::endl;
  }


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


  //If Brute force, store the "base name" and append number later (see for-loop)
  const std::string baseName = def.outputFileName;

  bool loopAgain = true;
  //for-loop to produce many output-files, but only if Brute force.
  for(int W = 0; W < def.nOutputs && loopAgain; W++){
    if(def.method == 'b')    //if bootstrap, we'll generate the noOutFiles
      loopAgain = false; //after the first loop run, so only loop once.

    //If we rerun the simulation many times, we must change the output file
    //name by appending numbers to them.
    if(def.nOutputs > 1 && def.method == 'B'){
      std::cout << "\nOutput-file number: " << W + 1 << " of "
                << def.nOutputs << std::endl;

      //change name of outfile to baseNameW, like: out.dat[number]
      std::stringstream temp;
      temp << baseName << W;
      def.outputFileName = temp.str();
    }

    std::vector<Jump> jumpRates;                  //returned by reference

    //store characteristic trait of chosen distribution
    float info;

    //generate jumprates for all crowding particles.
    computeJumpRates(jumpRates, info, nParticles, jumprateTracer,
                     jumprateCrowders, seedJump, jumprateDistribution);


    //initiate the save-class to store, save, and print result
    //(also computes standard error and correlation estimate)
    Save save(samplingTimes, ensembles, def.isLowMem.first);

    //non-important nice-to-have simulation info.
    RemainingTime printToScreen(ensembles);

    std::vector<Lattice> lattices;   //possibility to store many lattice obj.
    int noLattice = 1;               //default to one Lattice object...
    int lattice_index = 0;           //... stored in first element of vector
#   pragma omp parallel shared(save,lattices,noLattice) private(lattice_index)
    {
#     ifdef _OPENMP                      //If compiling with OpenMP...
      noLattice = omp_get_num_threads(); //change to threads" number of lattices.
#     endif

#     pragma omp single                  //only do this once in total.
      for(int i = 0; i < noLattice; ++i){
        //Initiate the Lattice class, which is what does all the physics.

        Lattice tmp(X, Y, Z, nParticles, seed * (i+1), isBoundaryFix);
        tmp.setSamplingTimes(samplingTimes, waitingTime);

        //Only actually needed for our "computeNakazato"-function.
        tmp.setJumpNaka(jumprateCrowders, jumprateTracer);
        tmp.setJumpRate(jumpRates);

        //if use Interaction algorithm, with InteractStr
        if(def.isInteracting)
          tmp.setInteraction(def.interactionStrength);

        lattices.push_back(tmp);
      }

#     ifdef _OPENMP                         //if compiling with OpenMP:
      lattice_index = omp_get_thread_num(); //each thread operates on its own
#     endif                                 //lattice index in lattices-vector

#    pragma omp for
      for(int E = 0; E < ensembles; ++E){
       if(lattice_index == 0)               //only first thread prints progress
         printToScreen.printProgress(E*noLattice);

       lattices[lattice_index].place();
       lattices[lattice_index].move();

      //Store tracer position at each point for this "ensembles", needed by
      //class to compute standard error/deviation. (+for binning, etc.)
      //extract displacement coordinates of the tracer...
      //store the current ensemble values in these:
      std::vector<int> dx,dy,dz;
      std::vector<double> dr;

      lattices[lattice_index].getDisplacement(dx,dy,dz,dr);

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
          << lattices[0].computeErgodicity(X,Y,Z)
          << "\t distr: " << dist[jumprateDistribution] << " (" << info
          << ")\t bound: " << bound[isBoundaryFix] << std::endl;
    print << "#D_eff: " << d_eff << "\t k_tagg: "                //line3
          << jumpRates[0].x.r << "\t D_av: " << d_av
          << "\t Interaction: " << onOff[def.isInteracting] << " "
          << def.interactionStrength << std::endl;

    std::string head = print.str();

    if(def.method != 'd')
      //calculate standard deviation, error-bars, and save to file
      save.save(def.outputFileName, head);
    else
      // dump raw trajectories (non-mean) to file. Do rest in script.
      save.dump(def.outputFileName, head);

    //prints distribution, and saves to "outputFileName" + "_histogram"
    if(printHistogram)
      save.computeDistribution(def.outputFileName, noHistogramBins);

    if(def.method == 'b')
      //bootstrap the shit out of this. Note, the bootstrapped output
      //files can be distinguished from the "real" simulation by them
      //ending with a number i: {0 < i < (noOutFiles - 1)}
      save.computeBootstrap(def.outputFileName, def.nOutputs,
                            jumpRates[0].x.r, head);

    if(def.isBootknife)
      //If true, use the hybrid "bootknife" method to compute slope,
      //and error, and write to separate file: "out.dat_bootknife"
      save.computeBootknife(def.outputFileName, def.nBootknife,
                            jumpRates[0].x.r);

    if(def.isJackknife)
      save.computeJackknife(def.outputFileName);


    // if(def.isShrinkage){
    //   //print H-matrix with shrinkage applied to it
    //   save.printShrinkage();
    // }


    // not implemented yet, XXX
    // if(def.isInvestigateConvergeance){
    //   save.investigateMuConvergeance();
    // }


  }  //END: of loop for re-running the simulation

  return 0;
}






void computeJumpRates(std::vector<Jump>& hopRate, float& info, int N, float jumpTracer,
                      float jumpCrowders, float seedJump, int jumpDistribution){

  //we store by hopRate.push_back() further down
  hopRate.clear();

  static Ran randomHopDistribution(seedJump);   //just any seed will do

  //Characteristic trait of the chosen distribution,
  //just used to print info to file/screen
  info = 0;

  for (int particle = 0; particle < N; particle ++){
    double u;
    do {
      u = randomHopDistribution.doub();
    }while (u == 0);

    int n = jumpDistribution;  //Choose distribution

    //NOTE: The jumprate is in EACH direction! Meaning, jumprate =1
    //is actually a jumprate =4 in 2D. (I think)
    if (0 <= n && n <= 3 && (u < 1 && u > 0) ){
      float lambda = 1.0;            //used in option 1 Exp-dist.
      float y_c = 1.0;               //used in option 2 Power-law
      float alpha = 0.5;             //used in option 2 Power-law

      switch(n){
      case 0:                        // 0 = uniform distribution
        u = u;
        info = 1.0;                  //Max value for uniform random number
        break;
      case 1:                        // 1 = exponential distribution
        u = (-log(u)) / lambda;
        u = 1.0/u;
        info = lambda;
        break;
      case 2:                        // 2 = power law distribution
        // p(y) = alpha/y_c * (y/y_c)^{-1-alfa} if y > y_c, || 0
        u = y_c * pow(u,-1/alpha);
        u = 1.0/u;                   //convert friction coeff. to jump rate
        info = alpha;
        break;
      case 3:                        // 3 = nakazato distribution
        if (particle == 0) u = jumpTracer;   //set juprate for the first
        else u = jumpCrowders;               //set juprate for the rest
        info = jumpTracer/jumpCrowders;
        break;
      }
    }
    else{
      if (n > 3 || n <0 )
        std::cout << "Invalid value/choice for prob.distribution ("<< n <<")"
                  << std::endl;
      else
        aux::printError("Random number for jumprate must be 0 < r < 1");
    }

    //Manually set the jumprate of the tracer particle (first one)
    if (particle == 0 ) u = jumpTracer;

    //save jump rate u to right and left direction
    Jump::Direction dir(u,u);

    Jump jump;

    //same jumprate in all dimensions (x,y,z)
    jump.x = dir;
    jump.y = dir;
    jump.z = dir;

    hopRate.push_back(jump);
  }
}
