#include <iostream>    //for std::cout, among others
#include <cstdlib>     //for abort()-function
#include <cmath>       //gives us sqrt, pow, fabs,
#include <string>
#include <vector>
#include <omp.h>
//#include <thread>      //new C++ standard! same as old boost_threads
#include <map>
#include <cassert>
#include <memory>      //for smart pointers
#include <sstream>     //for [o]stringstream

//Random number generator from numerical recipes 3 ed, but
//modified not to depend on nr3.h or else threads wont compile.
#include "nr/ran_mod.h"

#include "simpleini/SimpleIni.h" //parse ini-file
#include <stdlib.h>              //atoi

#include "auxiliary.h"           //non-physics stuff. (print messages etc.)
#include "classes.h"             //various data structures /classes (Jump, Particle)
#include "baselattice.h"         //this is the main class, that does the physics.
#include "lattices.h"            //this is the main class, that does the physics.
#include "save.h"                //class to save MSD, compute errors, and print to file


void computeJumpRates(std::vector<Jump>&, float&, int, float, float, float, int);

enum geometry {
  square1d = 1,
  square2d,
  square3d,
  honeycomb2d,
};

// Just to print info to output file.
inline std::ostream & operator<<(std::ostream & str, geometry X){
  switch (X) {
  case square1d: return str << "Square 1D";
  case square2d: return str << "Square 2D";
  case square3d: return str << "Square 3D";
  case honeycomb2d: return str << "Honeycomb 2D";
  default: return str << (int) X;
  }
}

int main(int argc, char* argv[]){

  InputValues def;                     //initiate default paramter values

  try{
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
      throw std::string("invalid inputfile path: " + std::string(input_file));
    }

    //start parsing "key = value" from input-file with section "lattice"
    const char * tmpValue = ini.GetValue("lattice", "fix_boundary", NULL /* default*/);
    bool isBoundaryFix = aux::convertToBool(tmpValue, "fix_boundary");

    geometry latticeGeometry = static_cast<geometry>(atoi(ini.GetValue("lattice", "geometry", NULL)));

    int latticeSize = atoi(ini.GetValue("lattice", "l", NULL));

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
         && jumprateDistribution == 3 && def.isInteracting))
      throw std::string("using interaction-code but not identical jumprates!");

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

    std::vector<Jump> jumpRates;            //returned by reference
    float info;                             //store trait of chosen distribution
    RemainingTime printToScreen(ensembles); //non-important nice-to-have simulation info

    //generate jumprates for all crowding particles
    computeJumpRates(jumpRates, info, nParticles, jumprateTracer,
                     jumprateCrowders, seedJump, jumprateDistribution);

    //initiate the save-class to store, save, and print result
    //(also computes standard error and correlation estimate)
    Save save(samplingTimes, ensembles, def.isLowMem.first);

    //store (pointers to) lattices, one for each thread, and they can be
    //of any derived class from BaseLattice, and set a run-time.
    std::vector<std::unique_ptr<BaseLattice> > lattices;

    int dim = 0;                                       // dimension of lattice
    int noLattice = 1;                                 //default to one Lattice object...
    int lattice_index = 0;                             //... stored in first element of vector
#   pragma omp parallel shared(save,lattices,noLattice) private(lattice_index)
    {
#     ifdef _OPENMP                      //If compiling with OpenMP...
      noLattice = omp_get_num_threads(); //change to threads" number of lattices.
#     endif

#     pragma omp single                   //only do this once in total.
      for(int i = 0; i < noLattice; ++i){ //Initiate the Lattice, (does all the physics)

        if(latticeGeometry == geometry::square1d){
          dim = 1;
          // lattices.push_back(std::make_unique<Square>(latticeSize,1,1, dim, nParticles, seed * (i+1), isBoundaryFix));
          lattices.push_back(std::unique_ptr<Square>(new Square(latticeSize,1,1, dim,
                                                                nParticles, seed * (i+1), isBoundaryFix)));
        }
        else if(latticeGeometry == geometry::square2d){
          dim = 2;
          lattices.push_back(std::unique_ptr<Square>(new Square(latticeSize,latticeSize,1, dim,
                                                                nParticles, seed * (i+1), isBoundaryFix)));
        }
        else if(latticeGeometry == geometry::square3d){
          dim = 3;
          lattices.push_back(std::unique_ptr<Square>(new Square(latticeSize,latticeSize,latticeSize, dim,
                                                                nParticles, seed * (i+1), isBoundaryFix)));
        }
        else if(latticeGeometry == geometry::honeycomb2d){
          dim = 2;
          lattices.push_back(std::unique_ptr<Honeycomb2d>(new Honeycomb2d(latticeSize, nParticles, seed * (i+1), isBoundaryFix)));
        }
        else
          throw std::string("Wrong lattice specified in input file"); // Note: must compile in serial to see this string

        lattices.back()->setSamplingTimes(samplingTimes, waitingTime);

        //Only actually needed for our "computeNakazato"-function.
        lattices.back()->setJumpNaka(jumprateCrowders, jumprateTracer);
        lattices.back()->setJumpRate(jumpRates);

        //if use Interaction algorithm, with InteractStr
        if(def.isInteracting)
          lattices.back()->setInteraction(def.interactionStrength);
      }

#     ifdef _OPENMP                         //if compiling with OpenMP:
      lattice_index = omp_get_thread_num(); //each thread operates on its own
#     endif                                 //lattice index in lattices-vector

#    pragma omp for
      for(int E = 0; E < ensembles; ++E){
        if(lattice_index == 0)              //only first thread prints progress
          printToScreen.printProgress(E*noLattice);

        //Store tracer position at each point for this "ensembles", to
        //compute standard error/deviation, etc.
        std::vector<int> dx,dy,dz;
        std::vector<double> dr;

        if (def.isBruteForce){              //only use one data point from each trajectory

          dx.reserve(samplingTimes.size()); //make room, for speed
          dy.reserve(samplingTimes.size());
          dz.reserve(samplingTimes.size());
          dr.reserve(samplingTimes.size());
          assert(dx.size() == 0);           //since using push back in for loop

          for(size_t k = 0; k < samplingTimes.size(); ++k) {

            std::vector<double>::iterator it_begin = samplingTimes.begin();
            lattices[lattice_index]->setSamplingTimes(std::vector<double>(it_begin, it_begin + k+1), waitingTime);

            lattices[lattice_index]->place();
            lattices[lattice_index]->move();

            std::vector<int> dx_sub,dy_sub,dz_sub; //save subset of x, y, z, data
            std::vector<double> dr_sub;
            lattices[lattice_index]->getDisplacement(dx_sub, dy_sub, dz_sub, dr_sub);

            dx.push_back(dx_sub.back());           //only store/use last element
            dy.push_back(dy_sub.back());
            dz.push_back(dz_sub.back());
            dr.push_back(dr_sub.back());
          }
        }
        else{
          lattices[lattice_index]->place();
          lattices[lattice_index]->move();
          lattices[lattice_index]->getDisplacement(dx, dy, dz, dr);
        }

#       pragma omp critical  //only one thread may write at a time:
        {
          save.store(dx,dy,dz,dr);
        }
      }
    }

    //To print which distribution we used to head of output-file
    std::string dist[]     = {"uniform","exponential","powerlaw","nakazato"};
    std::string onOff[]    = {"Off","On"};
    std::string bound[]    = {"periodic","fix"};
    float d_eff = lattices[0]->computeEffectiveDiffusionConst();
    float d_av = lattices[0]->computeAverageDiffusionConst();

    //print three lines of info about simulation to head of each outfile:
    std::ostringstream print;
    print << "#E = " << ensembles <<"\t N = " << nParticles << "\t latticeSize: " //line1
          << latticeSize <<"\t geometry: " << latticeGeometry << "\t 2*d*D_naka: "
          << lattices[0]->computeNakazato() << "  Waiting time: " << waitingTime << std::endl;
    print << "#Conc.: "<< (float) nParticles/pow(latticeSize,dim) << "\t MSD_equil: "  //line2
          << lattices[0]->computeErgodicity(latticeSize)
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
  }
  catch(std::string s)
  {
    std::cerr << "Error: " << s << "\n";
    return 1;
  }

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
        throw std::string("Invalid value/choice for prob.distribution ("
                          + tostring(n) + ")");
      else
        throw std::string("Random number for jumprate must be 0 < r < 1");
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
