#include <iostream>    //for cout, among others
#include <cstdlib>     //for abort()-function
#include <cmath>       //gives us sqrt, pow, fabs,
#include <string>
#include <vector>

#include <thread>      //new C++ standard! same as old boost_threads

#include <sstream>     //for [o]stringstream

//Random number generator from numerical recipes 3 ed, but
//modified not to depend on nr3.h or else threads wont compile.
#include "nr/ran_mod.h"

#include "auxiliary.h" //non-physics stuff. (print messages etc.)
#include "classes.h"   //various data structures /classes (Jump, Particle)
#include "lattice.h"   //this is the main class, that does the physics.
#include "save.h"      //class to save MSD, compute errors, and print to file

//HARD-CODED VARIABLES:

const int DISTRIBUTION  = 3;           //0 = uniform, 1 = exponential,
//                                      2 = power-law, 3 = nakazato
const float JUMPRATE_TRACER = 1;       //default value
const float JUMPRATE_CROWDER = 0.5;    //(if nakazato-distribution)
const bool FIXBOUNDARY =      0;       //1=fix wall, 0=periodic boundary
const bool EXPONENTAIL_WAITING_TIME = 1;
const bool PRINT_HISTOGRAM = false;     //MSD-distribution to out.dat_histogram

const double SEED_JUMP     = 17;       //usually: 17

const double SEED_LATTICE1 = 15;       //usually: 15,87,64,32
const double SEED_LATTICE2 = 87;
const double SEED_LATTICE3 = 64;
const double SEED_LATTICE4 = 32;


void computeJumpRates(vector<Jump>&, float&, const int,
                      const float, const float);


int main(int argc, char* argv[]){

  //DEFAULT VALUES:
  //--------------
  //All these default values can be changed by passing argument flags
  //at the command line when executing the binary file.

  //struct (from auxiliary.h) that keeps all default simulation options, that
  //are modified by passing argument flags. These are the default values:
  InputValues def;

  //use logarithmic spacing of sampling points
  def.isLogScale = false;

  //run in low memory mode
  def.isLowMem = false;

  //use nearest neighbor interaction, or bypass that algorithm.
  def.isInteracting = false;

  //print simulation information to screen
  def.isQuiet = false;

  //if  isInteracting = true, then use:
  def.interactionStrength = 0;

  //save simulation output to files with this prefix
  def.outputFileName = "out.dat";

  // file (path) to read i input values from
  def.inputFileName = "";

  //number of outputfiles/MSD-trajectories, if bootstrapping or brute force.
  def.nOutputs = 1;

  //is a char, use Bootstrap ('b') of Bruteforce ('B')?
  def.method = 'B';

  //use jack knife method?
  def.isJackknife = false;

  //jumprate of tracer particle, k_t
  def.jmpTracer = JUMPRATE_TRACER;

  //more simulation options:
  bool fixBoundaryOn = FIXBOUNDARY;  //boundary condition, false = periodic
  bool printHistogram = PRINT_HISTOGRAM;

  //change default values depending on arg. flags given on cmdline.
  argumentFlags(argc,argv, def);

  //----------END--------------




  int N;                            //Number of particles
  int X,Y,Z;                        //lattice size
  int ensemble,nSamplings;
  double maxTime;

  if (fixBoundaryOn)  cout << "# Fix boundary"     << endl;
  else                cout << "# Periodic boundary" << endl;

  //get the simultion parameters: If no inputfilename given, or could
  //not find it, run interactivley - i.e. quiz the user for parameters.
  if (def.inputFileName.empty() ||
      !readInputFile(def.inputFileName,X,Y,Z,N,ensemble,
                     nSamplings,maxTime))
    askUserForInputParameters(X,Y,Z,N,ensemble,nSamplings,maxTime);

  //Construct vector of sampling times
  vector<double> samplingTimes;
  if(!def.isLogScale){
    double minTime = 1;
    double deltaTime = (double) (maxTime - minTime) / nSamplings;
    double timeSum = minTime; //start time
    for (int i = 0; i < nSamplings; i++){
      samplingTimes.push_back(timeSum);
      timeSum = timeSum + deltaTime;
    }
  }
  else{
    //NOTE should take difference log(stop) - log(1)
    double deltaTimeLog = (double) log(maxTime)/nSamplings;
    double logTimeSum = log(1);
    while (logTimeSum < log(maxTime)){
      logTimeSum = logTimeSum + deltaTimeLog;
      samplingTimes.push_back( exp(logTimeSum) );
    }
  }


  //Initiate the Lattice class, which is what does all the physics.
  Lattice crowd1(X,Y,Z,N,SEED_LATTICE1,fixBoundaryOn);
  Lattice crowd2(X,Y,Z,N,SEED_LATTICE2,fixBoundaryOn);
  Lattice crowd3(X,Y,Z,N,SEED_LATTICE3,fixBoundaryOn);
  Lattice crowd4(X,Y,Z,N,SEED_LATTICE4,fixBoundaryOn);


  bool isExpWaitingTime = EXPONENTAIL_WAITING_TIME;
  crowd1.setSamplingTimes(samplingTimes, isExpWaitingTime);
  crowd2.setSamplingTimes(samplingTimes, isExpWaitingTime);
  crowd3.setSamplingTimes(samplingTimes, isExpWaitingTime);
  crowd4.setSamplingTimes(samplingTimes, isExpWaitingTime);

  //use Interaction algorithm, with InteractStr
  if(def.isInteracting){
    crowd1.setInteraction(def.interactionStrength);
    crowd2.setInteraction(def.interactionStrength);
    crowd3.setInteraction(def.interactionStrength);
    crowd4.setInteraction(def.interactionStrength);
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
      cout << endl << "Output-file number: " << W + 1 << " of "
           << def.nOutputs << endl;

      //change name of outfile to baseNameW, like: out.dat[number]
      stringstream temp;
      temp << baseName << W;
      def.outputFileName = temp.str();
    }


    float jumpCrowders = JUMPRATE_CROWDER;  //k_c (only if Nakazato)
    vector<Jump> jumpRates;                  //returned by reference

    //store characteristic trait of chosen distribution
    float info;

    //generate jumprates for all crowding particles.
    computeJumpRates(jumpRates, info, N, def.jmpTracer, jumpCrowders);

    //Only actually needed for our "computeNakazato"-function.
    crowd1.setJumpNaka(jumpCrowders, def.jmpTracer);
    crowd2.setJumpNaka(jumpCrowders, def.jmpTracer);
    crowd3.setJumpNaka(jumpCrowders, def.jmpTracer);
    crowd4.setJumpNaka(jumpCrowders, def.jmpTracer);


    crowd1.setJumpRate(jumpRates);
    crowd2.setJumpRate(jumpRates);
    crowd3.setJumpRate(jumpRates);
    crowd4.setJumpRate(jumpRates);


    //====================
    //When I use the interaction code, make sure I have Nakazato distribution
    //with identical jumprates!
    if ((def.isInteracting && DISTRIBUTION != 3) ||
        (JUMPRATE_CROWDER != JUMPRATE_TRACER
         && DISTRIBUTION == 3 && def.isInteracting)){
      cout << endl
           << "WARNING!! using interaction-code but not identical jumprates! \a"
           << endl << endl;
    }//====================


    //initiate the save-class to store, save, and print result
    //(also computes standard error and correlation estimate)
    Save save(samplingTimes,ensemble,def.isLowMem);

    //non-important nice-to-have simulation info.
    RemainingTime printToScreen(ensemble);

    for(int E = 0; E < ensemble;){ //E is iterated from the std::ref(E)
      if (!def.isQuiet)  printToScreen.printProgress(E);

      //When not using threads: (could as well use crowd1.generateTrajectory(E))
      //-------------------
      // crowd1.place();
      // crowd1.move();
      // ++E;

      //When using threads:
      //-------------------
      std::thread thrd1(&Lattice::generateTrajectory, &crowd1,std::ref(E));
      std::thread thrd2(&Lattice::generateTrajectory, &crowd2,std::ref(E));
      std::thread thrd3(&Lattice::generateTrajectory, &crowd3,std::ref(E));
      std::thread thrd4(&Lattice::generateTrajectory, &crowd4,std::ref(E));

      thrd1.join();
      thrd2.join();
      thrd3.join();
      thrd4.join();

      //Store tracer position at each point for this "ensemble",
      //needed by class to compute standard error/deviation. (+for
      //binning, etc.)  extract displacement coordinates of the
      //tracer...  store the current ensemble values in these:
      vector<int> dx,dy,dz;
      vector<double> dr;

      crowd1.getDisplacement(dx,dy,dz,dr);
      save.store(dx,dy,dz,dr);

      crowd2.getDisplacement(dx,dy,dz,dr);
      save.store(dx,dy,dz,dr);

      crowd3.getDisplacement(dx,dy,dz,dr);
      save.store(dx,dy,dz,dr);

      crowd4.getDisplacement(dx,dy,dz,dr);
      save.store(dx,dy,dz,dr);
    }

    //To print which distribution we used to head of output-file
    string dist[]     = {"uniform","exponential","powerlaw","nakazato"};
    string onOff[]    = {"Off","On"};
    string bound[]    = {"periodic","fix"};
    string waitTime[] = {"lin","exp"};
    float d_eff = crowd1.computeEffectiveDiffusionConst();
    float d_av = crowd1.computeAverageDiffusionConst();

    //print three lines of info about simulation to head of each outfile:
    ostringstream print;
    print << "#E = " << ensemble<<"\t N = " << N << "\t X-Y-Z: " //line1
          << X << "x" << Y << "x" << Z << "\t 2*d*D_naka: "
          << crowd1.computeNakazato() << "  Waiting time: "
          << waitTime[EXPONENTAIL_WAITING_TIME] << endl;
    print << "#Conc.: "<< (float) N/(X*Y*Z) << "\t MSD_equil: "  //line2
          << crowd1.computeErgodicity(X,Y,Z)
          << "\t distr: " << dist[DISTRIBUTION] << " (" << info
          << ")\t bound: " << bound[FIXBOUNDARY] << endl;
    print << "#D_eff: " << d_eff << "\t k_tagg: "                //line3
          << jumpRates[0].x.r << "\t D_av: " << d_av
          << "\t Interaction: " << onOff[def.isInteracting] << " "
          << def.interactionStrength << endl;

    string head = print.str();

    //calculate standard deviation, error-bars, and save to file
    save.save(def.outputFileName, head);

    //prints distribution, and saves to "outputFileName" + "_histogram"
    if(printHistogram) save.computeDistribution(def.outputFileName);

    if(def.method == 'b' && !def.isLowMem){
      //bootstrap the shit out of this. Note, the bootstrapped
      //output files can be distinguished from the "real" simulation
      //by them ending with a number i: {0 < i < (noOutFiles - 1)}

      save.setJumprate(jumpRates[0].x.r);

      save.computeBootstrap(def.outputFileName, def.nOutputs, head);
    }

    if(def.isJackknife && !def.isLowMem)
      save.computeJackknife(def.outputFileName);

  }  //END: of loop for re-running the simulation

  return 0;
}






void computeJumpRates(vector<Jump>& hopRate, float& info, const int N,
                      const float jumpTracer, const float jumpCrowders){

  //we store by hopRate.push_back() further down
  hopRate.clear();

  static Ran randomHopDistribution(SEED_JUMP);   //just any seed will do

  //Characteristic trait of the chosen distribution,
  //just used to print info to file/screen
  info = 0;

  for (int particle = 0; particle < N; particle ++){
    double u;
    do {
      u = randomHopDistribution.doub();
    }while (u == 0);

    int n = DISTRIBUTION;  //Choose distribution

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
        cout << "Invalid value/choice for prob.distribution ("<< n <<")"
             << endl;
      else
        printError("Random number for jumprate must be 0 < r < 1");
    }

    //Manually set the jumprate of the tracer particle (first one)
    if (particle == 0 ) u = jumpTracer;

    //save jump rate u to right and left direction
    Direction dir(u,u);

    Jump jump;

    //same jumprate in all dimensions (x,y,z)
    jump.x = dir;
    jump.y = dir;
    jump.z = dir;

    hopRate.push_back(jump);
  }
}
