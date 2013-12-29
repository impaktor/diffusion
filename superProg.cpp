#include <iostream>     //for cout, among others
#include <cstdlib>      //for abort()-function
#include <cmath>        //gives us sqrt, pow, fabs, 
#include <string>
#include <vector>

#include "nr/nr3.h"     //ran.h and ludcmp.h depends on this one. 
#include "nr/ran.h"     //Random number generator from numerical recipes 3 ed.  

#include "auxiliary.h"   //non-physics stuff. (print messages etc.)
#include "classes.h"    //various data structures /classes (Jump, Particle)
#include "superInteraction.h"  //this is the main class, that does the physics.

//BUG! LowMem 1: k_t 1, k_c 0.5, fix bound, 400^2, N=1, t_max 2000, E=1000, values=1000
//MEN! den skriver ut perfekt output??!!

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
  //--------------

  char NameOfFile[99] = "out.dat";

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
  
  //kanske ha den som invariabel i Save() istallet?
  //Set the switch true/false for use of low-mem run
  crowd.setLowMem(UseLowMem);

  static Ran tempnumb(17);                 //just anny seed will do
  vector<Jump> hopRate;
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

      //Characteristic trait of the chosen distribution, 
      //just used to print info to file/screen
      double info = 0;

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
        crowd.setDist(n,info);
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
    crowd.store(); 
    
  }

  crowd.save(ensemble, NameOfFile);   // calculate std.dev. and save to file

  if (!nymetod) cout <<"NOTE: Using the old moveOld()-function"<<endl;

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
  -ta bort dx dy dz dr.
  -kompilera med -Wall -pedantic
  -use __LINE__ to print line-number where error occured
  -use exceptions perhaps?
  -compare NRvector vs std::vector, for speed?
  
*/

//There are various testing functions in here, but they are either silenced by use of if-statements in combination with a boolean "test" variable set to false, or just commented out. A lot of work went into testing the output, therefore some things can seem redundant, such as printing un-squared displacements (<dx>~0) etc. 


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


