// test that I am sane, and the CTRW code does what it should?

#include <iostream>    //for std::cout, among others
#include <cmath>       //gives us sqrt, pow, fabs,
#include <string>
#include <vector>
#include <cstdlib>
#include <cassert>

//Random number generator from numerical recipes 3 ed, but
//modified not to depend on nr3.h or else threads wont compile.
#include "nr/nr3.h"
#include "nr/ran_mod.h"

Ran randomNumber(42);

double computeWaitingTime(){
  double r;
  do{
    r = randomNumber.doub();
  }while( r == 1 || r == 0 );

  // if 0 < alpha < 2: second mom inf. if alpha < 1, first mom inf.
  const double alpha = 0.5;
  const double a = 1;

  //one time step for crowding particle
//  double tau = a * std::pow(r, -1.0/alpha);

  // Used in Lloyd's paper, distribution: alpha/a * (1+tau/a)**(-1-alpha)
  double tau = a * (std::pow(r,-1.0/alpha) -1);

  //tau = a * (std::pow(r, -1.0/alpha) - 1.0); // correct version, but identical to above

  return tau;
}



void move(const std::vector<double> &samplingTimes, std::vector<double> &dr){

  // in orig code place() does this:
  double timeSum = 0;
  double pos = 0;
  dr.assign(samplingTimes.size(), 0);

  int i = 0;       //index of samplingTime-vector
  double tau;      //waiting time

  while(timeSum < samplingTimes.back()){

    tau = computeWaitingTime();

    //Save displacement if next time-step is beyond next sampling time
    while(timeSum <= samplingTimes[i] &&
          samplingTimes[i] < timeSum + tau  && i < samplingTimes.size()){

      //save displacement (from previous step)
      dr[i] = pos;

      i++;
    }

    double r = randomNumber.doub();

    timeSum = timeSum + tau;

    // make a move. store it.
    if(r < 0.5)
      pos = pos + 1;
    else
      pos = pos - 1;
  }
}


int main(int argc, char *argv[])
{

  // size_t N = 1e3;          // samplingpoints
  // size_t M = 1e6;          // "ensembles" aka. tranjectories
  // double maxTime = 1e6;
  if(argc != 4){
    std::cout << argv[0] << " <N> <M> <stopTime>" << std::endl;
    return 1;
  }
  size_t N = atoi(argv[1]);
  size_t M = atoi(argv[2]);
  double maxTime = (double) atoi(argv[3]);
  std::cout << "# N = " << N << " M = " << M << " stopTime = " << maxTime << std::endl;

  std::vector<double> sampling_times;
  std::vector<double> msd(N,0);
  // std::vector<std::vector<double> > trajectories;


  { // linear spacing of sampling times
    double minTime = 1;
    double deltaTime = (double) (maxTime - minTime) / N;
    double timeSum = minTime; //start time
    for(int i = 0; i < N; i++){
      sampling_times.push_back(timeSum);
      timeSum += deltaTime;
    }
  }

  // SIMULATE:
  for(size_t m = 0; m < M; ++m){
    std::vector<double> pos(N,0);
    move(sampling_times, pos);

    for(size_t i = 0; i < pos.size(); ++i)
      msd[i] += std::pow(pos[i],2);
  }

  //PROCESS:
  assert(N == sampling_times.size());
  for(size_t t = 0; t < sampling_times.size(); ++t)
    //msd[t] += std::pow(trajectories[m][t],2) / M;
    msd[t] /= M;

  //PRINT:
  for(size_t i = 0; i < msd.size(); ++i)
    std::cout << sampling_times[i] << "\t" << msd[i] << std::endl;

  return 0;
}
