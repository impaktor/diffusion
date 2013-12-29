//Class to show runtime information, such as remainging time
//until completion, and number of ensembles left. 

#include <iostream>        //cout 
#include <iomanip>         //setw & setprecision
#include <ctime>           //för time()
#include <time.h>

#include "classes.h"


Particle::Particle(int X, int Y, int Z){
  x = X;
  y = Y;
  z = Z;
}

Jump::Jump(double right, double left){
    r = right;
    l = left;
};

//overloading the << operator for the Particle class.
std::ostream& operator<< (std::ostream &out, Particle &particle){
  out << particle.x <<", " << particle.y <<", " << particle.z;
  return out;
}

//constructor
RemainingTime::RemainingTime(unsigned int totalEnsembles){
  startTime_ = time(NULL);
  totalEnsembles_ = totalEnsembles;
}

void RemainingTime::printProgress(unsigned int ensemble){
  time_t nowTime = time(NULL);

  //number of seconds:
  unsigned int elapsedTime = nowTime-startTime_; //works
  
  float remaining = (float) elapsedTime/(ensemble + 1) *
    (totalEnsembles_-ensemble)/60.0;
  
  float procentage = (float) ensemble/totalEnsembles_ *100;

  int ensemblesLeft = totalEnsembles_ - ensemble;

  //Delete previous line:
  std::cout<<"\r";

  //Print again!
  std::cout <<"#"<<std::setw(5)<<procentage<<" %  REMAINING ENSEMBLES:"
            <<std::setw(6)<<ensemblesLeft<<" Time (min): "<<std::fixed 
            <<std::setprecision(1)<<remaining<<std::flush;
  
  /*
If you are writing to stdout or cout or clog remember to fflush or std::flush the stream to make it output the line immediately. If you are writing to stderr or cerr then the stream is unbuffered and all output is immediate (and inefficient).
  */

}

