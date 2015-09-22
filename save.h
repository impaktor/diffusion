#ifndef SAVE_H
#define SAVE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <cassert>
#include "classes.h"

class Save{
public:
  //Use shorter aliases instead:
  typedef std::vector<int>                   vectorI_t;
  typedef std::vector<double>                vectorD_t;
  typedef std::vector<std::vector<int> >     matrixI_t;
  typedef std::vector<std::vector<double> >  matrixD_t;

  //Constructor. lowMem defaults to false, (bool lowMem=false) or else
  // we can't use boostrap.
  Save(vectorD_t samplingTimeVector, int noEnsembles, bool lowMem);

  //Each run/trajectory/"ensemble" should be saved by calling this
  //function. dr is displacement from starting position.
  void store(const vectorD_t& dr);

  //After simulation is complete, call this to print an output-file by
  //name "fileName" with structure:
  // [time] [MSD] [MSD_error] [correlation] [Z=H^(-1)*t]
  //plus "head" which is just information about simulation parameters.
  void save(std::string fileName, std::string head, bool printCov);

  // dump raw trajectories squared to files, set by output -w
  void dump(std::string fileName, std::string head);

  //To get P(x,t) i.e. a histogram, printed to file: "fileName_histogram"
  //i.e. it appends "_histogram" to the given output file name.
  void computeDistribution(std::string fileName, int noBins);

  //compute correct error in the fitted parameter. Prints result to
  //standard out (i.e. terminal, unless piped)
  void computeJackknife(std::string outName);

  //We need k_t in the bootknife/bootstrap when we try fitting to many
  //different starting times. (since we want t*k_t to be a constant,
  //we need to know k_t)

  //Add this after save() has run, to generate "numberOfRuns" synthetic
  //output files with name: "outNameN" where N is {0 < N < (numberOfRuns-1)}
  void computeBootstrap(std::string outName, int noOfRuns, double k_t,
                        std::string head);

  //generate synthetic trajectories and fit a line to each. Don't use
  //H-matrix
  void computeBootknife(std::string outName, int noOfRuns, double k_t);

  //run shrinkage and print H-matrix with shrinkage applied to it
  void printShrinkage();

  //used from read_data.cpp
  void printSlopeJackknifeShrinkage(std::string filename, int jackgroup);


private:

  //average of displacement (MSD) squared: <X^2>, <Y^2>,...
  vectorD_t r2_mu_;

  //Store std error, in each point
  vectorD_t r2_err_;

  //needed if we use stdErrLowMem()
  vectorD_t dr4_err_;

  //Switch between the two stdErr-functions above.
  bool isLowMem_;

  //save position for each ensemble and store them all
  matrixD_t store_dr2_;

  //Number of trajectories/simulation realizations.
  int noEnsembles_;

  //Number of elements in time vector (# of sampling points)
  int noSamplingTimes_;

  //If we want to change starting time, when performing the fitting in
  //the bootstrap --> "bootknife" version. Defaults to 0.
  float minTime_;

  //times to compute the MSD. (x-axis)
  vectorD_t samplingTime_;

  //standard error, to use as errorbars.
  void computeStdErr(void);
  void computeStdErrLowMem(vectorD_t&, vectorD_t&, const vectorD_t&);


  void computeLomholt(const matrixD_t &, const vectorD_t &,
                      const vectorD_t &, double &, double &);
  void getSlopeForOne_n(bool useShrinkage, int n, int, double &slope,
                        double &sigma, double &sigmaBad);
  void investigateMuConvergeance(int, std::string, bool);

  //generates more columns to print.
  void computePearsonCoefficient(vectorD_t &pearson) const;

  void computeHmatrix(const matrixD_t&, matrixD_t&);
  //void computeHmatrix2(matrixD_t&, bool);
  void computeHmatrix2(matrixD_t&);
  void computeHmatrix3(matrixD_t& , std::string, bool, bool);
  void computeBmatrix(const matrixI_t&, matrixD_t&, bool);

  //For Schäfer 2005 paper:
  void computeShrinkage(const matrixD_t&, matrixD_t &H_star);

  double function1(const matrixD_t &trajectories, const vectorD_t &time, bool);
  double function2(const matrixD_t &trajectories, const vectorD_t &time, bool);
  double function3(const matrixD_t &trajectories, const vectorD_t &time, bool);

  void jackknife(const matrixD_t&, const vectorD_t&, int, double&, double&,
                 double (Save::*f)(const matrixD_t&, const vectorD_t&, bool), bool);

  void computeSlope(const vectorD_t&, const vectorD_t&, const vectorD_t&,
                    double& , double&);

  inline void computeError(const matrixD_t& store, const vectorD_t& msd,
                           vectorD_t& sigma);

  //used by bootknife, among other
  template<typename T>
  void computeMean(const std::vector<std::vector<T> >& store,
                   vectorD_t& mu, double exponent=1){

    //don't use the "ensemble_" or noSamplingTimes_ variables, since
    // we will be fiddling with this in the jackknife implementation.
    const int ensembles = store.size();
    const int noSamplingTimes = store[0].size();
    const double inv = 1.0 / ensembles;

    mu.assign(noSamplingTimes, 0);

#   pragma omp parallel for
    for (int j = 0; j < noSamplingTimes; j++){
      for (int i = 0; i < ensembles; i++)
        mu[j] += pow(store[i][j], exponent) * inv;
    }
  }


  template<typename Type>
  void computeVariance(const std::vector<std::vector<Type> >& trajectories,
                       matrixD_t& matrix, bool isMeanZero){
    if(samplingTime_[0] == 0)
      std::cout << "NOTE: If the first sampling-point t_0 = 0 then the \n"
                << "matrix will be singular (first row and column =0)"
                << std::endl;

    //if we use this function when bootstrapping, and might not want the
    //same number of bootstrapped synthetic simulations as ensembles:
    //(otherwise we could use the global "noEnsembles_" variable)
    const double M = trajectories.size();
    assert(M > 0);

    //might not be N == noSamplingTimes_. (jackknife unbias test)
    int N = trajectories[0].size();
    matrix.assign(N, vectorD_t(N, 0));

    //initiate class, to print progress to screen,
    RemainingTime computationProgress(N);

    if(isMeanZero){
      double invM = 1.0 / M;
      for(int i = 0; i < N; i++){
        //print remaining progress to screen:
        computationProgress.printProgress(i);

#       pragma omp parallel for
        for(int j = 0; j < N; j++){
          double sum = 0;
          for(int m = 0; m < M; m++)
            sum += (trajectories[m][i]) * (trajectories[m][j]);
          matrix[i][j] = sum * invM;
        }
      }
    }
    else{
      vectorD_t mean(N, 0);
      computeMean<Type>(trajectories, mean, 1);

      double invM = (M == 1) ? 1.0 / M : 1.0 / (M - 1);

      for(int i = 0; i < N; i++){
        computationProgress.printProgress(i);

#       pragma omp parallel for default(shared)
        for(int j = 0; j < N; j++){
          double sum = 0;
          for(int m = 0; m < M; m++)
            sum += (trajectories[m][i] - mean[i]) * (trajectories[m][j] - mean[j]);
          matrix[i][j] = sum * invM;
        }
      }
    }
  }

};

#endif
