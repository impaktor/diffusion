#include <iostream>    //for cout, among others
#include <cstdlib>     //for abort()-function
#include <cmath>       //gives us sqrt, pow, fabs,
#include <string>
#include <vector>

//Random number generator from numerical recipes 3 ed, but
//modified not to depend on nr.h since that breaks boost::thread
#include "nr/ran_mod.h"

#include "auxiliary.h" //non-physics stuff. (print messages etc.)
#include "classes.h"   //various data structures /classes (Jump, Particle)
#include "lattice.h"   //this is the main class, that does the physics.
#include "save.h"      //class to save MSD, compute errors, and print to file

#include <sstream>     //XXX temp (multiple output-files)

#include <boost/thread/thread.hpp>

//BUG! LowMem 1: k_t 1, k_c 0.5, fix bound, 400^2, N=1, t_max 2000, E=1000, values=1000
//MEN! den skriver ut perfekt output??!!

// use ssh:
//      $ ssh -n 'cd /to/path; nice -19 ./program <input.dat'
// Or:  $ ssh -X nice -19 ./program output.dat < input.dat'

//HARD-CODED VARIABLES:
#define DISTRIBUTION     0          //0 = uniform, 1 = exponential,
//                                    2 = power-law, 3 = nakazato
#define JUMPRATE_TRACER  5.5        //
#define JUMPRATE_CROWDER 0.5        //(if nakazato-distribution)
#define FIXBOUNDARY      1          //1=fix wall, 0=periodic boundary
#define EXPONENTAIL_WAITING_TIME 1

#define SEED_JUMP         17
//TODO: testa med lowMem, mångakörningar!, parametrar ovan, och nedan.

/*
  Testar nu den nya v.6 med samplingTime istallet!
  -LowMem on/off   Identical. 1% minne vid On, 50 % vid off!
*/

//-----------------------------------------------------------
//now test my InteractionRoutine, with LowMem. (and without!):
// at -i 2, c=0.7 k_t=0.25, LowMem = ON  Identical!
// at -i 3, c=0.4 k_t=0.25, LowMem = OFF Identical!
// at -i 2, c=0.5 k_t=0.25, LowMem = ON  Identical!

//now test LowMem ON (and this new merged code) for interactions off
// at k_t=0.5 c=0.5,         2Duniform,  Identical!
// at k_t=0.5 c=0.5 -i 0 (!) 2Duniform,  Identical! (got warning about not using Nakazato)


//////////////////////////////////////
//   ****************************   //
//   * M A I N  F U N C T I O N *   //
//   ****************************   //
//////////////////////////////////////

int main(int argc, char* argv[]){
  int N;                             //Number of particles
  int X,Y,Z;                         //lattice size
  int ensemble,antalPunkter;
  double maxTime;

  //DEFAULT VALUES:
  //All these default values can be changed by passing argument flags
  //at the command line when executing the binary file. This is done
  //by the function argumentFlags() below.
  bool fixBoundaryOn = FIXBOUNDARY; //boundary condition, false = periodic
  bool logscale = false;            //spacing between sampling times log/lin
  bool useLowMem = false;           //don't save memory
  bool quiet = false;               //print simulation progress to screen
  string nameOfOutFile = "out.dat"; //default output name
  int noOutFiles = 1;               //only run one simulation (1 output file)
  int interactStr = 0;              //interaction strenth
  bool interactOn = false;          //don't use the interaction algorithm
  bool jackknife = false;           //if true: use jackknife, take lots of time
  char method;                      //use bootstrap ('b') or Brute force ('B')

  //change default values depending on arg. flags given on cmdline.
  argumentFlags(argc, argv, logscale,useLowMem, interactOn, quiet,
                interactStr, nameOfOutFile, noOutFiles, method,jackknife);

  if (fixBoundaryOn)
    cout <<"# Fix boundary"<<endl;
  else
    cout <<"# Periodic boundary"<<endl;

  AskUserForInputParameters(X,Y,Z,N,ensemble,antalPunkter,maxTime);

  //Initiate the Lattice class, which is what does all the physics.
  Lattice crowd1(X,Y,Z,N,fixBoundaryOn);
  Lattice crowd2(X,Y,Z,N,fixBoundaryOn);
  Lattice crowd3(X,Y,Z,N,fixBoundaryOn);
  Lattice crowd4(X,Y,Z,N,fixBoundaryOn);

  //If Brute force, store the "base name" and append number later (see for-loop)
  const std::string baseName = nameOfOutFile;

  bool loopAgain = true;
  //for-loop to produce many output-files, but only if Brute force.
  for(int W = 0; W < noOutFiles && loopAgain; W++){
    if(method == 'b')    //if bootstrap, we'll generate the noOutFiles
      loopAgain = false; //after the first loop run, so only loop once.

    //If we rerun the simulation many times, we must change the output file
    //name. We only need to do this for Brute force; computeBootstrap() has
    //its own function to rename files, by appending numbers to file names.
    if(noOutFiles > 1 && method == 'B'){
      cout << endl << "Output-file number: " << W + 1 << " of "
           << noOutFiles << endl;

      //change name of outfile to baseNameW, like: out.dat[number]
      stringstream temp;
      temp << baseName << W;
      nameOfOutFile = temp.str();
    }

    //Construct vector of sampling times
    vector<double> samplingTimes;
    if(!logscale){
      double minTime = 1;
      double deltaTime = (double) (maxTime - minTime) / antalPunkter;
      double timeSum = minTime; //start time
      for (int i = 0; i < antalPunkter; i++){
        samplingTimes.push_back(timeSum);
        timeSum = timeSum + deltaTime;
      }
    }
    else{
      //NOTE should take difference log(stop) - log(1)
      double deltaTimeLog = (double) log(maxTime)/antalPunkter;
      double logTimeSum = log(1);
      while (logTimeSum < log(maxTime)){
        logTimeSum = logTimeSum + deltaTimeLog;
        samplingTimes.push_back( exp(logTimeSum) );
      }
    }
    bool expWaitingTime = EXPONENTAIL_WAITING_TIME;
    crowd1.setSamplingTimes(samplingTimes, expWaitingTime);
    crowd2.setSamplingTimes(samplingTimes, expWaitingTime);
    crowd3.setSamplingTimes(samplingTimes, expWaitingTime);
    crowd4.setSamplingTimes(samplingTimes, expWaitingTime);

    //use Interaction algorithm, with InteractStr
    if(interactOn){
      crowd1.setInteraction(interactStr);
      crowd2.setInteraction(interactStr);
      crowd3.setInteraction(interactStr);
      crowd4.setInteraction(interactStr);
    }

    static Ran randomHopDistribution(SEED_JUMP);   //just anny seed will do


    vector<Jump> hopRate;

    //Characteristic trait of the chosen distribution,
    //just used to print info to file/screen
    float info = 0;

    for (int particle = 0; particle < N; particle ++){
      double u;
      do {
        u = randomHopDistribution.doub();
      }while (u == 0);

      int n = DISTRIBUTION;  //Choose distribution

      float jumpTracer;

      //NOTE: The jumprate is in EACH direction! Meaning, jumprate =1
      //is actually a jumprate =4 in 2D. (I think)
      if (0 <= n && n<=3 && (u < 1 && u > 0) ){
        jumpTracer = JUMPRATE_TRACER;  //Jumprate for the tracer paricle!
        float lambda = 1.0;            //used in option 1 Exp-dist.
        float y_c = 1.0;               //used in option 2 Power-law
        float alpha = 0.5;             //used in option 2 Power-law
        float jumpCrowders = JUMPRATE_CROWDER; //used in option 3 Nakazato
        crowd1.setJumpNaka(jumpCrowders, jumpTracer);
        crowd2.setJumpNaka(jumpCrowders, jumpTracer);
        crowd3.setJumpNaka(jumpCrowders, jumpTracer);
        crowd4.setJumpNaka(jumpCrowders, jumpTracer);

        switch(n){
        case 0:      // 0 = uniform distribution
          u = u;
          info = 1.0;               //Maximum value for uniform random number
          break;
        case 1:      // 1 = exponential distribution
          u = (-log(u)) / lambda;
          u = 1.0/u;
          info = lambda;
          break;
        case 2:      // 2 = power law distribution
          // p(y) = alpha/y_c * (y/y_c)^{-1-alfa} if y > y_c, || 0
          u = y_c * pow(u,-1/alpha);
          u = 1.0/u; //convert friction coefficient to jump rate
          info = alpha;
          break;
        case 3:      // 3 = nakazato distribution
          if (particle == 0) u = jumpTracer;   //set juprate for the first
          else u = jumpCrowders;               //set juprate for the rest
          info = jumpTracer/jumpCrowders;
          break;
        }
      }
      else{
        if (n > 3 || n <0 )
          cout << "Invalid value/choise for prob.distribution ("<< n <<")" << endl;
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

    crowd1.setJumpRate(hopRate);
    crowd2.setJumpRate(hopRate);
    crowd3.setJumpRate(hopRate);
    crowd4.setJumpRate(hopRate);



    //====================
    //When I use the inteeraction code, make sure I have Nakazato distribution
    //with identical jumprates!
    if ((interactOn && DISTRIBUTION != 3) ||
        (JUMPRATE_CROWDER != JUMPRATE_TRACER
         && DISTRIBUTION == 3 && interactOn)){
      cout << endl
           << "WARNING!! using interaction-code but not identical jumprates! \a"
           << endl << endl;
    }//====================


    //initiate the save-class to store, save, and print result
    //(also computes standard error and correlation estimate)
    Save save(samplingTimes,ensemble,useLowMem);

    //non-important nice-to-have simulation info.
    RemainingTime printToScreen(ensemble);

    for(int E = 0; E < ensemble;){ //E is iterated from the ref(E)

      if (!quiet)  printToScreen.printProgress(E);

      //crowd1.place();
      //crowd1.move();

      boost::thread thrd1(&Lattice::generateTrajectory, &crowd1,boost::ref(E));
      boost::thread thrd2(&Lattice::generateTrajectory, &crowd2,boost::ref(E));
      boost::thread thrd3(&Lattice::generateTrajectory, &crowd3,boost::ref(E));
      boost::thread thrd4(&Lattice::generateTrajectory, &crowd4,boost::ref(E));

      thrd1.join();
      thrd2.join();

      if (!quiet)  printToScreen.printProgress(E);

      thrd3.join();
      thrd4.join();

      //make a print out of individual ensemble:
      //crowd.dumpSimulation(E);

      // //TEST: to get a snapshot close to the MaxTime
      // if(totalTime >= MaxTime-0.01){
      //   crowd.Snapshot();
      //   cout << "totalTime:" << totalTime << endl;
      // }


      // TODO: ej uppdaterad efter min move()-implementering!
      //  if(InteractOn){//INTERACTION
      //     This is just to save the information of cluster-size distribution
      //     crowd.saveCluster(totalTime);
      // }//BUG Augusti, gar otroligt trogt med denna!


      //Store tracer position at each point for this ensemble, needed by
      //class to compute standard error/deviation. (+for binning, etc.)
      //extract displacement coordinates of the tracer...
      //store the current ensemble values in these:
      vector<int> dx,dy,dz;
      vector<double> dr;

      crowd1.getDisplacement(dx,dy,dz,dr);

      //... and insert in our save-function:
      save.store(dx,dy,dz,dr);

      crowd2.getDisplacement(dx,dy,dz,dr);
      save.store(dx,dy,dz,dr);

      crowd3.getDisplacement(dx,dy,dz,dr);
      save.store(dx,dy,dz,dr);

      crowd4.getDisplacement(dx,dy,dz,dr);
      save.store(dx,dy,dz,dr);
    }

    if(interactOn){//INTERACTION
      //print a histogram of clusters
      //crowd.saveBinning(ensemble);
    }//BUG Augusti, gar otroligt trogt med denna!


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
          << hopRate[0].x.r << "\t D_av: " << d_av
          << "\t Interaction: " << onOff[interactOn] << " "
          << interactStr << endl;

    string head = print.str();

    //calculate standard deviation, errorbars, and save to file
    save.save(nameOfOutFile, head);

    if(method == 'b' && !useLowMem){
      //bootstrap the shit out of this. Note, the bootstrapped
      //output files can be destingushed from the "real" simulation
      //by them ending with a number i: {0 < i < (noOutFiles - 1)}
      save.computeBootstrap(nameOfOutFile, noOutFiles, head);
    }

    //Won't work if we use lowMem-algorithm.
    if (jackknife){
      if (!useLowMem ) save.computeJackknife(nameOfOutFile);
      else cout << "cant use jackknife when in low memory mode! " << endl;
    }

  }  //END: of loop for re-running the simulation


  return 0;
}





/*
  LEFT TO DO:
  ----------
  -remove various testing-loops
  -plot directly with gnuplot.
  -remove the set_interaction --> constructor! (doesn't work)
  -rename all variables using camelCase
  -kanske ha en MSD-klass, double x y z r?
  -kompilera med -Wall -pedantic
  -use __LINE__ to print line-number where error occured
  -use exceptions perhaps?
  -compare NRvector vs std::vector, for speed?
*/

//There are various testing functions in here, but they are either silenced by
//use of if-statements in combination with a boolean "test" variable set to
//false, or just commented out. A lot of work went into testing the output,
//therefore some things can seem redundant, such as printing un-squared
//displacements (<dx>~0) etc.

/*

  OUTLINE OF THE FUNCTION CALLS OF THIS PROGRAM:
  ---------------------------------------------

 With                      Without
 Interaction:              Interaction:
 -----------               ------------
 main()                     main()
  |                          |
 move()                     move()
  |                          |
 moveAndBoundaryCheck()     moveAndBoundaryCheck()
  |                          |
 vacancyCheck()             vacancyCheck()
  |                          |
 Interaction()               |
  |                          |
 CountNeighbors()            |
  |                          |
  |                          |
 [make the actual move]     [make the actual move]

*/

