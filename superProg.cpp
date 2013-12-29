#include <iostream>     //for cout, among others
#include <cstdlib>      //for abort()-function
#include <cmath>        //gives us sqrt, pow, fabs,
#include <string>
#include <vector>

#include "nr/nr3.h"     //ran.h and ludcmp.h depends on this one.
#include "nr/ran.h"     //Random number generator from numerical recipes 3 ed.

#include "auxiliary.h"  //non-physics stuff. (print messages etc.)
#include "classes.h"    //various data structures /classes (Jump, Particle)
#include "superInteraction.h"  //this is the main class, that does the physics.
#include "save.h"

// use ssh:
//      $ ssh -n 'cd /to/path; nice -19 ./program <input.dat'
// Or:  $ ssh -X nice -19 ./program output.dat < input.dat'

//HARD-CODED VARIABLES:
#define DISTRIBUTION     0    //0 = uniform, 1 = exponential,
//                              2 = power-law, 3 = nakazato
#define JUMPRATE_TRACER  0.5  //
#define JUMPRATE_CROWDER 0.5  //(if nakazato-distribution)
#define FIXBOUNDARY      0    //1 = fix wall, 0 = periodic boundary
#define EXPONENTAIL_WAITING_TIME 0 //


//////////////////////////////////////
//   ****************************   //
//   * M A I N  F U N C T I O N *   //
//   ****************************   //
//////////////////////////////////////

int main(int argc, char* argv[]){
  int N;                              //Number of particles
  int X,Y,Z;                          //lattice size
  int ensemble,antalPunkter;
  double maxTime;

  //DEFAULT VALUES:
  bool fixBoundaryOn = FIXBOUNDARY;     //boundary condition, false = periodic
  bool logscale = false;
  bool nymetod = true;                  //use move() instead of moveOld()
  bool UseLowMem = false;               //don't save memory
  bool quiet = false;                   //print remaining ensembles during simulation
  string NameOfFile = "out.dat";
  int InteractStr = 0;
  bool InteractOn = false;

  //change default values depending on arg. flags given
  argumentFlags(argc, argv, logscale,UseLowMem, InteractOn,
                quiet, InteractStr, NameOfFile);

  //TEST print result to screen:
  //-----------------------
  string OnOff[] = {"OFF","ON"};
  int tempIndex = 0;
  if (InteractOn) tempIndex = 1;
  cout << "###log = " << logscale << endl
       <<"###interaction strength " << OnOff[tempIndex]
       << " =" << InteractStr << endl
       <<"###filename = " << (string) NameOfFile << endl;
  //-----------------------

  if (fixBoundaryOn)
    cout << "# Fix boundary" << endl;
  else
    cout << "# Periodic boundary" << endl;

  AskUserForInputParameters(X,Y,Z,N,ensemble,antalPunkter,maxTime);

  //initiate my class:
  SuperInteraction crowd(X,Y,Z,N,fixBoundaryOn,InteractStr);

  static Ran tempnumb(17);                 //just anny seed will do
  vector<Jump> hopRate;

  //Characteristic trait of the chosen distribution,
  //just used to print info to file/screen
  double info = 0;

  if (nymetod){
    for (int particle = 0; particle < N; particle ++){
      double u;
      do {
        u = tempnumb.doub();
      }while (u == 0);


      int n = DISTRIBUTION;  //Choose distribution
      // 0 = uniform
      // 1 = exponential
      // 2 = power law
      // 3 = nakazato

      float jumpTracer;

      //NOTE: The jumprate is in EACH direction! Meaning,
      //jumprate =1 is actually a jumprate = 4 in 2D. (I think)
      if (0 <= n && n <= 3 && (u < 1 && u > 0) ){
        jumpTracer = JUMPRATE_TRACER;  //Jumprate for the tracer paricle!
        float lambda = 1.0;            //used in option 1 Exp-dist.
        float y_c = 1.0;               //used in option 2 Power-law
        float alpha = 0.5;             //used in option 2 Power-law
        float jumpCrowders = JUMPRATE_CROWDER; //used in option 3 Nakazato
        crowd.setJumpNaka(jumpCrowders, jumpTracer);

        switch(n){
        case 0:
          u = u;             //ie. jumprate is uniform (as I state in my thesis)
          info = 1.0;        //Maximum value for uniform random number
          break;
        case 1:
          u = (-log(u)) / lambda;
          u = 1.0/u;
          info = lambda;
          break;
        case 2:
          //powerlaw, p(y) = alfa/y_c * (y/y_c)^{-1-alfa} if y > y_c, 0 otherwise
          u = y_c * pow(u,-1/alpha);
          u=1.0/u;                   //convert friction coefficient to jump rate
          info = alpha;
          break;
        case 3:
          if (particle == 0) u = jumpTracer;   //set juprate for the first
          else u = jumpCrowders;             //set juprate for the rest
          info = jumpTracer/jumpCrowders;
          break;
        }
      }
      else{
        if (n>3 || n <0 )
          cout <<"Invalid value/choise for prob.distribution ("<<n<<")"<<endl;
        else
          cout <<"Random number for jumprate must be 0 < r < 1"<<endl;

        abort();
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

  //if using jump rates (i.e. not moveOld()-function)
  if (nymetod){
    crowd.setJumpRate(hopRate);
  }


  //make vector of sampling times
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
  crowd.setSamplingTimes(samplingTimes, expWaitingTime);

  //initiate the save-class to store, save, and print result
  //(also computes standard error and correlation estimate)
  Save save(samplingTimes,ensemble,UseLowMem);

  //store the current ensemble values in these:
  vector<int> dx,dy,dz;
  vector<double> dr;

  //non-important nice-to-have simulation info.
  RemainingTime printToScreen(ensemble);

  for(int E = 0; E < ensemble; E++){

    if (!quiet)  printToScreen.printProgress(E);

    crowd.place();

    crowd.move();

    //make a print out of individual ensemble:
    //crowd.dumpSimulation(E);


    //Store tracer position at each point for this ensemble, needed by
    //class to compute standard error/deviation. (+for binning, etc.)
    //extract displacement coordinates of the tracer...
    crowd.getDisplacement(dx,dy,dz,dr);

    //Store tracer position at each point for this ensemble, needed by
    //class to compute standard error/deviation. (+for binning, etc.)
    save.store(dx,dy,dz,dr);

  }




  //To print which distribution we used to head of output-file
  string dist[] = {"uniform","exponential","powerlaw","nakazato"};
  string onOff[] = {"Off","On"};
  string bound[] = {"periodic","fix"};
  float d_eff = crowd.computeEffectiveDiffusionConst();
  float d_av = crowd.computeAverageDiffusionConst();

  //print three lines of info about simulation to head of file:
  ostringstream print;
  print << "#E = " << ensemble<<"\t N = " << N << "\t X-Y-Z: "
        << X << "x" << Y << "x" << Z << "\t 2*d*D_naka: "
        << crowd.computeNakazato() << endl;
  print << "#Conc.: "<< (float) N/(X*Y*Z) << "\t MSD_equil: "
        << crowd.computeErgodicity(X,Y,Z)
        << "\t distr: " << dist[DISTRIBUTION] << " (" << info
        << ")\t bound: " << bound[FIXBOUNDARY] << endl;
  print << "#D_eff: " << d_eff << "\t k_tagg: "
        << hopRate[0].x.r << "\t D_av: " << d_av
        << "\t Interaction: " << onOff[InteractOn] << " "
        << InteractStr << endl;

  string head = print.str();

  save.save(NameOfFile, head);

  return 0;
}

