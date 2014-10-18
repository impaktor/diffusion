#include <iostream>        //for cout, among others
#include <fstream>         //read write to files
#include <cmath>           //gives sqrt, pow, fabs,
#include <vector>          //STL vector
#include <string>          //STL string
#include <functional>      //STL function, needed for lambda-functions
#include <algorithm>       //STL for min/max element finding in histogram
#include <cassert>         //for assert()
#include <omp.h>           //openMP = multi processing, declared with #pragmas
#include "nr/nr3.h"
#include "nr/ludcmp.h"       //only depends on nr3.h
#include "nr/svd.h"          //singular value decomposition
#include "nr/ran.h"          //only needed for my bootstrap-function

#include "nr/gamma.h"        //needed by fitab.h (see NR p. 784)
#include "nr/incgammabeta.h" //needed by fitab.h
#include "nr/fitab.h"        //needed by computeSlope()

#include "classes.h"         //to print bootstrap progress to screen
#include "auxiliary.h"       //for the aux::printError()-function
#include "save.h"

//#include "nr/gaussj.h"      //needed by fitmrq.h
//#include "nr/fitmrq.h"      //not yet in use...

using std::pow;
using std::sqrt;

//number of different sampling times N to try (0=off):
const int BIAS_TEST_ITER = 0; //100

//if true: y = a + bx, if false: y = bx
const bool USE_NR_FITTING = false;

//seed for bootstrap and bootknife.
//(both will give identical virtual simulation result)
const float SEED_BOOTSTRAP = 6854;

//If true, we print each bootstrapped trajectory to separate file:
//out.dat0, out.dat1, etc, (assuming -w out.dat). Else (if "false") we
//compute the H-matrix directly, and output to out.dat_matrix
const bool BOOTSTRAP_CONSERVE_MEMORY = false;


void printHmatrix(const Save::matrixD_t&, std::string);

Save::Save(vectorD_t samplingTime,
           int ensembles, bool lowMem){

  samplingTime_ = samplingTime;

  //discard any sampling points lower than this. Currently only for
  //computeSlope (uncorrelated fit).
  minTime_ = 0;

  //number of realizations / trajectories.
  noEnsembles_ = ensembles;

  //number of elements
  noSamplingTimes_ = samplingTime_.size();

  //use the memory conserving function
  isLowMem_ = lowMem;

  //These will be our 4th and 5th column in the main output file.
  pearson_coefficient_.assign(noSamplingTimes_, 0);
  z_correlation_.assign(noSamplingTimes_, 0);

  //initiate the vectors for the summation in store()
  //only needed if we use the memory-conserving algorithm
  if(isLowMem_){
    dr4_err_.assign(noSamplingTimes_, 0);
    dx4_err_.assign(noSamplingTimes_, 0);
    dy4_err_.assign(noSamplingTimes_, 0);
    dz4_err_.assign(noSamplingTimes_, 0);

    r2_mu_.assign(noSamplingTimes_, 0);
    x2_mu_.assign(noSamplingTimes_, 0);
    y2_mu_.assign(noSamplingTimes_, 0);
    z2_mu_.assign(noSamplingTimes_, 0);
  }
  else{
    //Allocate size in memory for the matrix to store every
    //individual trajectory, will use push_back later, so don't "assign()"
    store_dx_.reserve(noEnsembles_);
    store_dy_.reserve(noEnsembles_);
    store_dz_.reserve(noEnsembles_);
    store_dr2_.reserve(noEnsembles_);
  }
}


void Save::computePearsonCoefficient(){
  //This function calculates the Pearson coefficients (see Numerical
  //recipes ch 14), to get the correlation between two sampling points
  //since if the tagged particle jumps once/sec sampling the system
  //10/s will not yield any new information, and it will mess up the
  //error estimation that assumes the sampling points are uncorrelated
  //and with a Gaussian distribution.

  //perhaps check: Physical Review D vol 34, 2469 (1986), sec. D

  if(isLowMem_)               //need data to work on
    aux::printError("Can not compute correlation, no available data.");

  if(r2_mu_.empty())          //need data to work on
    aux::printError("Must compute MSD first! (ie. r2_mu_)");

  //compute the correlation between two points, i, j, for i = j+1
  vectorD_t h(noSamplingTimes_,0);

  for(int i = 1; i < noSamplingTimes_; i++){
    for(int j = 0; j < noEnsembles_; j++)
      h[i] += 1.0 / (noEnsembles_ - 1) *
        (store_dr2_[j][i] - r2_mu_[i])*
        (store_dr2_[j][i-1] - r2_mu_[i-1]);
  }

  // pearson_coefficient_[i] = h_ij / (sigma_i*sigma_j)
  // sigma_r = r2_err_ * sqrt(ensembles);
  for(int i = 1; i < noSamplingTimes_; i++){
    if(r2_err_[i-1] != 0) //zero error at time=0, avoid division by 0.
      pearson_coefficient_[i] = (double) h[i] /
        (r2_err_[i] * r2_err_[i-1] * noEnsembles_);
  }
}



void computeZvector(const Save::matrixD_t &H,
                    const Save::vectorD_t &t, Save::vectorD_t &Z){
  //This function computes the H-matrix: H*z=t, and saves z.

  //Number of sampling points. Different if using this from jackknife
  //bias removal algorithm.
  int N = H.size();

  //In order to use NR implementations to test the matrix, and perform
  //LU-decomposition, we transform STL-matrix 'H' to NR-matrix
  NRmatrix<double> h(N, N, 0.0);

  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      h[i][j] = H[i][j];

  if(false)
    //Test the integrity of matrix. For some retarded reason it rounds
    //off the answer! i.e gives 0.0, instead of 3.5E-5 etc.
    std::cout << "conditional number: " << setprecision(5) << scientific
              << SVD(h).inv_condition() << "\tsize:\t" << H.size() << std::endl;

  //Solve Equation, with LU-Decomposition
  //-------------------------------------
  NRvector<double> b(N), z(N);

  for(int i = 0; i < N; i++)
    b[i] = t[i];

  LUdcmp LUdecomposition(h);
  LUdecomposition.solve(b,z);

  Z.assign(N, 0);

  //Convert NR vector to std::vector
  for(int i = 0; i < N; i++)
    Z[i] = z[i];
}


void Save::computeZ(std::string filename){
  //Compute H*z=t, and save z to column in out file

  //compute entire matrix, to check the correlation between every
  //value
  matrixD_t H;

  // computeHmatrix3(H, filename, true, true);
  // printHmatrix(H, filename + "_matrix3c");

  // computeHmatrix3(H, filename, false, false);
  // printHmatrix(H, filename + "_matrix3b");

  // computeHmatrix3(H, filename, true, false);
  // printHmatrix(H, filename + "_matrix3");

  // computeHmatrix2(H, false);
  // printHmatrix(H, filename + "_matrix2b");

  // computeHmatrix2(H, true);
  // printHmatrix(H, filename + "_matrix2");

  //TODO: don't use if POW waiting time / CTRW
  computeHmatrix(store_dr2_, H);
  printHmatrix(H, filename + "_matrix1");

  computeZvector(H, samplingTime_, z_correlation_);
}


void Save::store(const vectorI_t& dx, const vectorI_t& dy,
                 const vectorI_t& dz, const vectorD_t& dr){
  //This function stores all the trajectories, so we can process the data,
  //like error bars in computeStdErr().


  // Error in each point i is (in LaTeX notation):
  // (a[i])^2 = 1/N \sum_i^N (a_i - <a>)^2 = [1/N \sum_i^N a[i]^2] - [<a>]^2

  // <a>= 1/N \sum_i^N a[i]   (i.e. average value in that point)
  // N = Number of ensembles

  // The original way to do this is save the value in each point a[i] for
  // each ensemble, but that leaves us with a giant matrix if we have many
  // ensembles N, and many data points D, (N*D ~> 10^9) ~ 1Gb of mem usage
  // (10^9/2^30 ~ 1Gb). Therefore with the switch 'lowMem=true' we do the
  // left hand side of the equation above, and then we don't need a matrix
  // but just two vectors, for each dimension.


  if(isLowMem_){
    //Check that we assigned the vectors, (non-empty):
    assert(!dr4_err_.empty() && !dx4_err_.empty() &&
           !dy4_err_.empty() && !dz4_err_.empty());

    assert(!r2_mu_.empty() && !x2_mu_.empty() &&
           !y2_mu_.empty() && !z2_mu_.empty());

    //To conserve memory we do the summation directly, which might
    //cause a numerical error in the end compared to doing it the
    //more rigorous way, but all tests I've made returns identical
    //output
    for(int i = 0; i < noSamplingTimes_; i++){
      dr4_err_[i] += pow( dr[i] ,4);
      dx4_err_[i] += pow( dx[i] ,4);
      dy4_err_[i] += pow( dy[i] ,4);
      dz4_err_[i] += pow( dz[i] ,4);

      r2_mu_[i] += pow( dr[i] ,2);
      x2_mu_[i] += pow( dx[i] ,2);
      y2_mu_[i] += pow( dy[i] ,2);
      z2_mu_[i] += pow( dz[i] ,2);
    }
  }
  else{
    //if we don't use the lowMem-algorithm:
    //this gives dx[ensembles][number of values]
    store_dx_.push_back(dx);
    store_dy_.push_back(dy);
    store_dz_.push_back(dz);

// //TMP store displacement squared. TODO: have lattice return this instead. XXX
//     vectorD_t tmp(dr.size(),0);
//     for(size_t i = 0; i < dr.size(); ++i)
//       tmp[i] = pow(dr[i], 2);

//     store_dr2_.push_back(tmp);

//TMP2: Pigeon says the dr being fed into store from read_data.cpp
//might be wrong, so I try using the dx,dy,dz instead:
    vectorD_t dr_tmp(dx.size(),0);
    for(size_t i = 0; i < dx.size(); ++i)
      dr_tmp[i] = pow(dx[i], 2) + pow(dy[i], 2) + pow(dz[i], 2);

    store_dr2_.push_back(dr_tmp);
  }
}


void Save::save(std::string name, std::string head){
  //calculates the MSD from the stored positions of the tracer, &
  //prints to file.

  if(isLowMem_){
    //use the memory conserving version instead. This one does not
    //use the dx[][]-matrix, and in that way we save lots of mem.
    computeStdErrLowMem(x2_mu_, x2_err_, dx4_err_);
    computeStdErrLowMem(y2_mu_, y2_err_, dy4_err_);
    computeStdErrLowMem(z2_mu_, z2_err_, dz4_err_);
    computeStdErrLowMem(r2_mu_, r2_err_, dr4_err_);

  }
  else{

    //Compute <R^2>:
    computeMean<int>(store_dx_, x2_mu_, 2); //"2" = square trajectory
    computeMean<int>(store_dy_, y2_mu_, 2); //"2" = square trajectory
    computeMean<int>(store_dz_, z2_mu_, 2); //"2" = square trajectory
    computeMean<double>(store_dr2_, r2_mu_);

    //Calculate standard error, to generate error-bars this function
    //must come after x2_mu_...r2_mu_ have been computed.
    computeStdErr();

    computePearsonCoefficient();  //4th column (defaults to zeroes)
    computeZ(name);               //5th column (defaults to zeroes)

    printMSD(head, name); //print main results


    if(BIAS_TEST_ITER != 0){
      //For testing convergeance of slope as func. of samplingpoints.
      investigateMuConvergeance(BIAS_TEST_ITER, name + "_converge_jackknife_lomholt", false);
      investigateMuConvergeance(BIAS_TEST_ITER, name + "_converge_jackknife_shrinkage_lomholt", true);
    }
  }
}

// dump raw trajectories squared to file (with baseName)
void Save::dump(std::string baseName, std::string head){
  if(isLowMem_){
    // Not implemented, but then I should print each dr passed into
    // Save::store(..dr..) as it is called, but then Save class needs
    // to know if I'm running in "dumped" mode.
    assert(false);
  }
  else{

    for(size_t i = 0; i < store_dr2_.size(); i++){
		std::ostringstream fileName;
		fileName << baseName << "_" << i;
      std::ofstream file(fileName.str().c_str());
      file << head;
      file << "# t \t r^2 " << std::endl;

      for(int j = 0; j < noSamplingTimes_; j++)
        file << samplingTime_[j] << "\t" << store_dr2_[i][j] << std::endl;
      file.close();
      std::cout << "finished file:\t" << fileName.str() << std::endl;
    }
  }
}

void Save::printSlopeJackknifeShrinkage(std::string filename, int jackgroups){
  //To run with Pigeon's data, from read_data.cpp which removes the
  //first time points

  double sigma, slope, sigmaBad;
  int n = std::sqrt(noEnsembles_ /2.0);
  getSlopeForOne_n(true, n, jackgroups, slope, sigma, sigmaBad);  //reduce H-matrix size

  std::ofstream fileOut((filename + "_shrink_jack_slope").c_str());
  fileOut << "# mu: " << slope << "\n"
          << "# sigma: " << sigma << "\n"
          << "# crappy sigma: "<< sigmaBad << std::endl;
  fileOut << "# used shrinkgae, jackknife, and sqrt(M/2) big H-matrix" << std::endl;
  fileOut.close();
}

// void Save::printShrinkage(){
//     //test of shrinkage method for H-matrix.

//     matrixD_t H_star;
//     computeShrinkage(store_dr2_, H_star);
//     printHmatrix(H_star, filename + "_shrinkage_matrix");
// }



void Save::printMSD(std::string head, std::string name){

  std::ofstream primary(name.c_str());
  primary   << head;
  std::cout << head;

  primary << "# t \t MSD \t err \t pearson \t z" << std::endl;

  for(int i = 0; i < noSamplingTimes_; i++){
    primary << samplingTime_[i] << "\t" << r2_mu_[i] << "\t" << r2_err_[i]
            << "\t" << pearson_coefficient_[i] << "\t" << z_correlation_[i]
            << std::endl;
  }
  primary.close();
}


void Save::computeStdErr(void){
  //Calculate the standard deviation, to get error-bars.

  //Should not be using this function if lowMem = true!
  assert(!isLowMem_);

  x2_err_.assign(noSamplingTimes_, 0);
  y2_err_.assign(noSamplingTimes_, 0);
  z2_err_.assign(noSamplingTimes_, 0);
  r2_err_.assign(noSamplingTimes_, 0);

  //Optimization: avoid a[i] / b, use: c = 1/b; a[i] * c;
  double tmp = 1.0 / (1.0*noEnsembles_*(noEnsembles_ - 1));

  assert(x2_mu_.size() == noSamplingTimes_);
  assert(store_dx_.size() == noEnsembles_);

# pragma omp parallel for
  for(int i = 0; i < noSamplingTimes_; i++){
    for(int m = 0; m < noEnsembles_; m++){
      x2_err_[i] += pow (pow(store_dx_[m][i],2) - x2_mu_[i] ,2);
      y2_err_[i] += pow (pow(store_dy_[m][i],2) - y2_mu_[i] ,2);
      z2_err_[i] += pow (pow(store_dz_[m][i],2) - z2_mu_[i] ,2);
      r2_err_[i] += pow (store_dr2_[m][i] - r2_mu_[i] ,2);
    }
    x2_err_[i] = sqrt( x2_err_[i] * tmp);
    y2_err_[i] = sqrt( y2_err_[i] * tmp);
    z2_err_[i] = sqrt( z2_err_[i] * tmp);
    r2_err_[i] = sqrt( r2_err_[i] * tmp);
  }    //for E = 1 --> inf due to tmp!
}


void Save::computeStdErrLowMem(vectorD_t &i2_mu, vectorD_t &i2_err,
                               const vectorD_t &di4_err){
  //compute standard error, if using lowMem-algorithm.
  assert(isLowMem_);

  //Optimization: avoid a[i] / b, use: c = 1/b; a[i] * c;
  double tmp = 1.0 / noEnsembles_;

  //Multiply all elements in i2_mu with tmp:
  //std::for_each(i2_mu.begin(), i2_mu.end(), [&](double &n){n = n * tmp;});
  for(int i = 0; i< noSamplingTimes_; i++)
    i2_mu[i] = i2_mu[i] * tmp;

  i2_err.assign(noSamplingTimes_, 0);

  for(int i = 0; i < noSamplingTimes_; i++)
    i2_err[i] = di4_err[i] / (1.0 * noEnsembles_) - pow( i2_mu[i] ,2);

  tmp = 1.0 / (noEnsembles_ - 1.0);
  //std::for_each(i2_err.begin(), i2_err.end(), [&](double &n){n = sqrt(n * tmp);});
  for(int i = 0; i < noSamplingTimes_; i++)
    i2_err[i] = sqrt(tmp * i2_err[i]);
}


void Save::computeDistribution(std::string fileName, int noBins){
  //This function is for getting P(t,x), ie. the probability
  //distribution of finding a particle at r^2, at time t.

  if(isLowMem_){
    std::cout << "lowMem = ON, therefore no information stored!\n"
              << "can not print cluster-information" << std::endl;
  }
  else{

    class Histogram{
    private:
      vectorD_t data_;       //what we want to generate a histogram from
      size_t noBins_;        //number of bins (resolution)
      double binStep_;       //step size
      vectorD_t probAxis_, binAxis_;  //y-, and x-axis. P(x)

    public:
      Histogram(const vectorD_t& input, int number_of_bins,
                std::string name){
        data_ = input;                  //data to make histogram of
        noBins_ = number_of_bins;       //resoultion

        probAxis_.assign(noBins_, 0.0); //Probability count (Y-axis)
        binAxis_.assign(noBins_, 0.0);  //binning-axis      (X-axis)

        computeHistogram();
        printToFile(name);
      }

      void computeCumulative(std::string fileName){
        //compute the comulative distribution. This is better to deal
        //with the singularity in 1D in the origin.

        vectorD_t cumulative(noBins_, 0);

        for(size_t i = 1; i < noBins_; ++i)
          cumulative[i] = cumulative[i-1] + probAxis_[i];

        std::ofstream c(fileName.c_str());
        for(size_t i = 0; i < noBins_; i++)
          c << binAxis_[i] + binStep_/2.0
            << "\t" << cumulative[i] << std::endl;

        c.close();
        std::cout << "Printed cumulative distribution to: "
                  << fileName << std::endl;
      }

    private:
      void computeHistogram(){

        //maximuma and min value in the histogram
        double min = *std::min_element(data_.begin(), data_.end());
        double max = *std::max_element(data_.begin(), data_.end());
        binStep_ = (double) (max - min)/noBins_;

        //Build the X-axis scale
        double temp = min;
        for(size_t i = 0; i < noBins_; i++){
          binAxis_[i] = temp;
          temp += binStep_;
        }

        double value; //displacement, current value to add to histogram

        for(size_t e = 0; e < data_.size(); e++){  //loop over trajectories
          value = data_[e];

          for(size_t j = 0; j < noBins_; j++){    //loop over bin-cells
            assert(j <= probAxis_.size());

            // //why 2*PI ? r? (November 2010)
            double factor = data_.size();// * 2*3.1415*(value+binStep_)/2;
            if(binAxis_[j+1] > value && value >= binAxis_[j]){
              probAxis_[j] += 1.0 / factor;
              break;
            }
          }
        }

        //TEST, alternative way to get the MSD, by using the distribution
        //this is a good way to test it.
        double MSD = 0;
        for(size_t i = 0; i < noBins_; i++)
          MSD += binAxis_[i] * probAxis_[i];
        std::cout << "MSD calculated from distribution: " << MSD << std::endl;
      }

      void printToFile(std::string fileName){

        std::ofstream hist(fileName.c_str());

        for(size_t i = 0; i < noBins_; i++){
          if(probAxis_[i] != 0)
            hist << binAxis_[i] + binStep_/2.0 << "\t" << probAxis_[i] << std::endl;
        }

        hist.close();
        std::cout << "Printed histogram to: " << fileName << std::endl;
      }
    };

    int samplingTime = noSamplingTimes_ -1; //pick sampling time

    vectorD_t tempVectorR2(noEnsembles_, 0);
    vectorD_t tempVectorX2(noEnsembles_, 0);
    vectorD_t tempVectorY2(noEnsembles_, 0);

    for(int i = 0; i < noEnsembles_; i++){
      //(I don't know of a better way to copy the trajectories)
      tempVectorR2[i] = store_dr2_[i][samplingTime];
      tempVectorX2[i] = pow(store_dx_[i][samplingTime], 2);
      tempVectorY2[i] = pow(store_dy_[i][samplingTime], 2);
    }

    Histogram r2(tempVectorR2, noBins, fileName + "_histogram");
    Histogram x2(tempVectorX2, noBins, fileName + "_histogramX2");
    Histogram y2(tempVectorY2, noBins, fileName + "_histogramY2");

    //Also print the cumulative distribution for r^2
    r2.computeCumulative(fileName + "_cumulative");
  }
}


void Save::computeBootstrap(std::string name, int noOfRuns,
                            double k_t, std::string head){
  //Generate synthetic values, from current trajectories in store_dr2_
  //to get new values, by Quick-and-dirty Monte Carlo Bootstrapping.

  //Needs all the information for all trajectories saved in store_dr2_,
  //which will not be saved in lowMem mode.
  if(isLowMem_){
    std::cout << "Error: Can not use bootstrap when in low memory mode"
              << std::endl;
  }
  else{
    std::cout << "bootstrapping!" << std::endl;

    matrixD_t store_dr2(noEnsembles_, vectorD_t(noSamplingTimes_,0));

    Ran randomNumber(SEED_BOOTSTRAP);

    //Either print "noOfRuns" number of synthetic output-files to
    //hard drive, each file corresponds one MSD of a synthetic
    //simulation, or save all MSDs in mem, then save to one file.
    bool conserveMemory = BOOTSTRAP_CONSERVE_MEMORY;
    matrixD_t catanateOutputFiles;
    if(!conserveMemory)
      catanateOutputFiles.assign(noOfRuns,
                                 vectorD_t(noSamplingTimes_,0));

    //This part lets us choose different starting times
    // k_t * t should have these values, corresponding to errors:
    // 1.010  = 2.45
    // 1.0050 = 7.34
    // 1.0020 = 25.40
    // 1.0010 = 58.34
    // 1.0005 = 140
    const double tempArray[] = {0, 2.45, 7.34, 10, 20, 25.40, 30, 40, 50,
                                58.34, 70, 80, 90, 100, 110, 120, 130, 140};

    size_t sizeOfTempArray = (sizeof tempArray)/(sizeof tempArray[0]);

    //construct a vector with a number of starting times for the fitting
    vectorD_t startTimeForFitting;

    //now transform to time values:
    for(size_t i = 0; i < sizeOfTempArray; i++){
      float timeMark = (float) tempArray[i] / k_t;
      startTimeForFitting.push_back(timeMark);
    }

    int noReFittings = startTimeForFitting.size();

    matrixD_t mu(noOfRuns, vectorD_t(noReFittings,0));
    matrixD_t sigma_mu(noOfRuns, vectorD_t(noReFittings,0));

    //non-important nice-to-have simulation info.
    //this now makes save.cpp depend on classes.[h,cpp]
    RemainingTime progressToScreen(noOfRuns);

    //Loop to generate "noOfRuns" number of synthetic MSDs of
    //trajectories, each saved MSD corresponds to a simulation of
    //"noEnsembles_" number of individual trajectories (i.e. ensembles).
    for(int k = 0; k < noOfRuns; k++){
      progressToScreen.printProgress(k);

      //Randomly pick real trajectories into our virtual simulation
      for(int i = 0; i < noEnsembles_; i++){
        int r = (int) (noEnsembles_ * randomNumber.doub());
        store_dr2[i] = store_dr2_[r];
      }

      //Compute <R^2>, store result in r2_mu.
      vectorD_t r2_mu(noSamplingTimes_, 0);
      computeMean<double>(store_dr2, r2_mu);

      if(conserveMemory){
        //Print data to file (for processing in MatLab):
        std::stringstream tempName;
        //append run-iterator (k) to filename:
        tempName << name << k;
        std::ofstream bootOutput((tempName.str()).c_str());

        bootOutput << head;
        for(int i = 0; i < noSamplingTimes_; i++)
          bootOutput  << samplingTime_[i] << "\t" << r2_mu[i] << std::endl;
        bootOutput.close();
      }
      else{
        //store all bootstraped data in a single matrix...
        catanateOutputFiles[k] = r2_mu;
        //...and print it later, to a single file instead.
      }
    } //(end of for(k)-loop

    if(!conserveMemory){
      //Since we did not dump each MSD-trajectory to file, we
      //compute the H-matrix and dump that to a single file.

      //Compute the mean of the Mean Sq. Disp. (don't square it!)
      //use noOfRuns instead of noEnsembles_ since might not be the same
      vectorD_t meanOfMsd(noSamplingTimes_,0);
      computeMean(catanateOutputFiles, meanOfMsd, 1);

      matrixD_t H;
      //computeHmatrix1(H, catanateOutput, meanOfMsd, 1); //OLD SYNTAX TODO
      computeHmatrix(H, catanateOutputFiles);

      //print H-matrix to file, with "_matrix_bootstrap" appended:
      printHmatrix(H, name + "_matrix_bootstrap");

      //Print bootstrapped mean of mean  to file:
      //-----------------------------------------
      std::ofstream msdBootout((name + "_bootstrap").c_str());
      msdBootout << head;
      for(int i = 0; i < noSamplingTimes_; i++)
        msdBootout  << samplingTime_[i] << "\t" << meanOfMsd[i] << std::endl;
      msdBootout.close();
    }
  }
}



void Save::computeBootknife(std::string name, int noOfRuns, double k_t){
  //this is an alternative "bootknife" version of bootstrap, where we
  //refrain from using H-matrices, and just fit a line to current
  //bootstrapped "simulation" (as if uncorrelated, like in the
  //jackknife-method), and then take all my slopes and compute an
  //average, (which must be done outside this for(k)-loop).

  //Needs all the information for all trajectories saved in store_dr2_,
  //which will not be saved in lowMem mode.
  if(isLowMem_){
    std::cout << "Error: Can not use bootknife when in low memory mode"
              << std::endl;
  }
  else{
    std::cout << "bootknifing!" << std::endl;

    matrixD_t store_dr2(noEnsembles_, vectorD_t(noSamplingTimes_,0));

    Ran randomNumber(SEED_BOOTSTRAP);

    // k_t * t should have these values, corresponding to errors:
    // 1.010  = 2.45
    // 1.0050 = 7.34
    // 1.0020 = 25.40
    // 1.0010 = 58.34
    // 1.0005 = 140
    const double tempArray[] = {0, 2.45, 7.34, 10, 20, 25.40, 30, 40, 50,
                                58.34, 70, 80, 90, 100, 110, 120, 130, 140};

    size_t sizeOfTempArray = (sizeof tempArray)/(sizeof tempArray[0]);

    //construct a vector with a number of starting times for the fitting
    vectorD_t startTimeForFitting;

    //now transform to time values:
    for(size_t i = 0; i < sizeOfTempArray; i++){
      float timeMark = (float) tempArray[i] / k_t;
      startTimeForFitting.push_back(timeMark);
    }

    int noReFittings = startTimeForFitting.size();

    matrixD_t mu(noOfRuns, vectorD_t(noReFittings,0));
    matrixD_t sigma_mu(noOfRuns, vectorD_t(noReFittings,0));

    //non-important nice-to-have simulation info.
    //this now makes save.cpp depend on classes.[h,cpp]
    RemainingTime progressToScreen(noOfRuns);

    //Loop to generate "noOfRuns" number of synthetic MSDs of
    //trajectories, to each we fit a line.
    for(int k = 0; k < noOfRuns; k++){
      progressToScreen.printProgress(k);

      //Randomly pick real trajectories into our virtual simulation
      for(int i = 0; i < noEnsembles_; i++){
        int r = (int) (noEnsembles_ * randomNumber.doub());
        store_dr2[i] = store_dr2_[r];
      }

      //Compute <R^2>, store result in r2_mu.
      vectorD_t r2_mu(noSamplingTimes_, 0);
      computeMean<double>(store_dr2, r2_mu);

      //return error in each sampling-point by dr2_err, needed when
      //fitting mu. New for each bootstrapped "simulation"
      vectorD_t dr2_err;
      computeError(store_dr2, r2_mu, dr2_err);

      //Patch to modify/test minTime for fitting. get mu & sigma_mu.
      for(int i = 0; i < noReFittings; i++){
        minTime_ = startTimeForFitting[i];
        computeSlope(r2_mu, dr2_err, samplingTime_, mu[k][i], sigma_mu[k][i]);
      }
    }

    //only needs to be a vector if we fit for several different
    //starting times.
    vectorD_t mu_average(noReFittings, 0);
    vectorD_t sigma_mu_final(noReFittings, 0);
    computeMean(mu, mu_average, 1);

    for(int i = 0; i < noReFittings; i++ ){
      for(int j = 0; j < noOfRuns; j++)
        sigma_mu_final[i] += pow(mu[j][i] - mu_average[i] ,2);
      sigma_mu_final[i] = sqrt(sigma_mu_final[i] /((double) (noOfRuns-1)));
    }
    //NOTE/TODO Am I not using sigma_mu? Why do I need this?

    //Print bootstrapped/knifed mean of mean to file:
    //-----------------------------------------
    std::ofstream bootknifeOut((name + "_bootknife").c_str());

    bootknifeOut  <<  "#t_{start},  mu,  and sig_mu " << std::endl;

    for(int i = 0; i < noReFittings; i++ )
      bootknifeOut << startTimeForFitting[i] << "\t" << mu_average[i] << "\t"
                   << sigma_mu_final[i] << std::endl;

    bootknifeOut.close();
  }
}


void Save::computeLomholt(const matrixD_t &H, const vectorD_t &time,
                          const vectorD_t &dr2, double &mu, double &sigma){
  //Does what I do in python, i.e. do the Lomholt procedure to get the
  //slope (mu) and error (sigma) form the correlation matrix H.
  assert(H.size() > 0 && H[0].size() > 0);

  int N = H.size();
  vectorD_t z;
  computeZvector(H, time, z);

  double t_dot_z = 0, dr2_dot_z = 0;
  for(int i = 0; i < N; ++i){
    t_dot_z += time[i] * z[i];
    dr2_dot_z += dr2[i] * z[i];
  }

  mu = dr2_dot_z / t_dot_z;
  sigma = 1.0 / std::sqrt(t_dot_z);
}

void prune(const Save::matrixD_t &Trajectories, Save::matrixD_t &traject,
           const Save::vectorD_t &Time, Save::vectorD_t &time, int desiredSize){
  //remove elements to reduce the number of sampling-times to
  //desiredSize. Operate on TRAJECTORIES (not H-matrix), time and
  //MSD-vector.

  const double step = 1.0 * Time.size()/ desiredSize;

  time.assign(desiredSize, 0);
  traject.assign(Trajectories.size(), Save::vectorD_t(desiredSize, 0)); //same outer size (trajectories)

  for(int i = 0; i < desiredSize; ++i){
    size_t I = i*step;
    //std::cout  << "I:\t" << I << std::endl;
    assert(I <= Time.size());
    time[i] = Time[I];
    for(size_t m = 0; m < Trajectories.size(); ++m){
      traject[m][i] = Trajectories[m][I];
    }
  }
}



void Save::getSlopeForOne_n(bool useShrinkage, int n, int noGroups,
                            double &slope, double &sigma, double &sigmaBad){
  //return jackknifed slope, error, and (a bad) estimate of error in
  //slope; where we used same number of trajectories, but with fewer
  //(n) sampling points:

  matrixD_t trajectories;
  vectorD_t time;
  prune(store_dr2_, trajectories, samplingTime_, time, n);
  assert(time.size() == trajectories[0].size());

  //Remove bias from slope which we compute from Lomholt, possibly use shrinkage in jackknife
  jackknife(trajectories, time, noGroups, slope, sigmaBad, &Save::function2, useShrinkage);

  double trash; //dont need the return from jackknife

  //Remove bias from error in slope from the Lomholt/correlated fitting.
  jackknife(trajectories, time, noGroups, sigma, trash, &Save::function3, useShrinkage);
}



void Save::investigateMuConvergeance(int iterations, std::string filename,
                                     bool useShrinkage){
  //Run Jackknife on each number of removed sampling points, and in
  //the jackknife use the correlation matrix to compute the slope, for
  //subset of trajectories.

  //Store slope and error from each run with different number of
  //samplingtimes:
  vectorD_t mu(iterations, 0);
  vectorD_t sigma(iterations, 0);
  vectorD_t crappySigma(iterations, 0);  //to show that the jackknife error is crap

  //speed up in jeckknife by removing whole chunks of
  //trajectories. noGroups = how many chunks to divide all the
  //trajectories in, to remove one chunk at a time.
  const int noGroups = 100; //TODO: have it hard coded? XXX
  const double step = 1.0 * noSamplingTimes_ / iterations;

  assert(iterations <= noSamplingTimes_);
  assert(r2_mu_.size() > 0 && noGroups <= noSamplingTimes_);

  for(int i = 1; i <= iterations; ++i){
    int n = i * step;
    std::cout << "Doing: " << i << "\tof: " << iterations << " n = " << n << std::endl;
    getSlopeForOne_n(useShrinkage, n, noGroups, mu[i-1], sigma[i-1], crappySigma[i-1]);
  }

  std::ofstream unBias((filename + "").c_str());

  unBias << "#N \t mu\t sigma\t jacksigma:" << std::endl;
  for(int i = 0; i < iterations; ++i)
    unBias << (i+1)*step << "\t" << mu[i] << "\t" << sigma[i] << "\t" << crappySigma[i] << std::endl;
  unBias.close();
}



double Save::function1(const matrixD_t &trajectories, const vectorD_t &time, bool){
  //Compute slope from trajcetories, assuming all data is
  //uncorrelated. Simple diagonal chi-square fit.  last two
  //bool-arguments are needed since must have same number of arguments
  //as function2().

  vectorD_t dr2J(noSamplingTimes_, 0), sigma;
  computeMean<double>(trajectories, dr2J);
  computeError(trajectories, dr2J, sigma);

  double trash, slope;
  computeSlope(dr2J, sigma, time, slope, trash);
  return slope;
}


double Save::function2(const matrixD_t &trajectories, const vectorD_t &time,
                          bool useShrinkage){
  //Compute slope using correlated fitting, with covariance matrix.

  matrixD_t H;
  if(useShrinkage)
    computeShrinkage(trajectories, H);  //compute super good covariance
  else
    computeHmatrix(trajectories, H);    //compute standard covariance

  vectorD_t dr2;
  computeMean(trajectories, dr2);

  double sigma, slope;
  computeLomholt(H, time, dr2, slope, sigma);

  return slope;
}


double Save::function3(const matrixD_t &trajectories, const vectorD_t &time,
                          bool useShrinkage){
  //Compute slope using correlated fitting, with covariance matrix.

  matrixD_t H;
  if(useShrinkage)
    computeShrinkage(trajectories, H);  //compute super good covariance
  else
    computeHmatrix(trajectories, H);    //compute standard covariance

  vectorD_t dr2;
  computeMean(trajectories, dr2);

  double sigma, slope;
  computeLomholt(H, time, dr2, slope, sigma);

  return sigma;
}




void Save::jackknife(const matrixD_t &data, const vectorD_t &time, int g, double &X_J,
                     double &X_err_J, double (Save::*f)(const matrixD_t&, const vectorD_t&, bool),
                     bool useShrinkage){
  //http://www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife/
  //use the jackknife method to compute the error from correlated data, by
  //successively removing parts of the data (i.e. trajectories)

  //Compute variable X_J, and X_err_J from X.
  //  X_J = X * g - (g-1) * 1/g * sum_k^g X_j,k
  //  X_err_J = (g-1) * 1/g * sum_k^g (X_k_j - X)^2

  //only works if we have all trajectories. (trajectories-matrix)
  if(isLowMem_){
    std::cout << "Error: Can not use Jackknife when in low memory-mode."
              << std::endl;
  }
  else{
    std::cout << "Jackknifing!" << std::endl;

    //store all the jackknifed versions of the fitted parameter
    vectorD_t x_j(g, 0);

    //matrix of all simulation trajectories that we will "Jackknife"
    matrixD_t data_J;

    //number of trajectories in each group
    double h = 1.0 * data.size() / g;
    assert(g > 1 && g <= noEnsembles_); //g = 1 --> remove all trajectories.

    for(int i = 0; i < g; i++){
      data_J = data;

      //remove (cycklic) block of size h trajectories:
      data_J.erase(data_J.begin() + int(i*h),
                   data_J.begin() + int((i+1)*h));

      assert(data_J.size() >= (noEnsembles_ - (h + 1)));
      assert(data_J.size() <= (noEnsembles_ - (h - 1)));

      //Compute a new x_j for these specific trajectories:
      x_j[i] = (this->*f)(data_J, time, useShrinkage);
    }

    //Compute un-jackknifed fitted parameter X to data in trajectory:
    double X = (this->*f)(data, time, useShrinkage);
    X_err_J = 0;

    //compute the standard "jackknifed" error sigmaJ2 to the slope,
    //which is the end result of the jackknife method.
    double x_j_sum = 0;
    double gInv = 1.0 / g;

    for(int i = 0; i < g; i++){
      X_err_J += gInv * (g - 1.0) * pow((x_j[i] - X), 2);
      x_j_sum += gInv * x_j[i];
    }

    X_J = X*g - (g - 1) * x_j_sum;            //correction for bias.
  }
}


void Save::computeJackknife(std::string outName){
  //Public callable from main(), to remove bias from estimated
  //slope. Boldly assumes the trajectories are uncorrelated.

  double slope, sigmaJ2;
  bool foo = 0; //trash variable.
  jackknife(store_dr2_, samplingTime_, noEnsembles_, slope,
            sigmaJ2, &Save::function1, foo);

  //Print result, append "_jackknife" to the chosen output-filename
  std::ofstream printJack((outName + "_jackknife").c_str());

  printJack << "#Slope (mu), jackknifed error (sigma_J):"
            << setprecision(5) << scientific << "\n" << slope << "\n"
            << sigmaJ2 << std::endl;
  printJack.close();
}



void Save::computeSlope(const vectorD_t &dr, const vectorD_t &dr_err,
                        const vectorD_t &time, double &mu, double &sigma_mu){
  //fit a straight line to data, either by using NR: y = a + b * x
  //or using my own implementation: y = mu * x

  //since might have reduced the set for jackknife-test:
  int n = time.size();

  //Patch to remove the first values when doing the fit. Note: think
  //whether we want to use fitab or fitb when doing this.
  int minElement = 0;
  for(minElement = 0; minElement < n &&
        time[minElement] < minTime_; minElement++){}
    //just want to iterate our "minElement" to find minimum index

  if(!USE_NR_FITTING){
    //My own implementation of Ch. 15 of Numerical Recipes, where fitab.h
    //has been re-written for a y = mu*t fit. Note that Stt corresponds to
    //Sxx in NR.

    double Sty = 0, Stt = 0;

    for(int i = minElement; i < n; i++){
      Stt += 1.0 * pow(time[i],2) / pow(dr_err[i],2);
      Sty += 1.0 * time[i] * dr[i] / pow(dr_err[i],2);
    }

    //return the error in mu
    sigma_mu = sqrt(1.0 / Stt);

    //return the slope, mu (by reference)
    mu = Sty / Stt;
  }
  else{
    //If we want to use the NR code for fitting:

    //fitab wants NRvector...
    NRvector<double> sigma(n - minElement, 0.0);
    NRvector<double> X(n - minElement, 0.0);
    NRvector<double> Y(n - minElement, 0.0);

    //...so copy input STL-vectors to NRvector
    for(int i = 0; i < n - minElement; i++){
      X[i] = time[minElement + i];
      Y[i] = dr[minElement + i];
      sigma[i] = dr_err[minElement + i];
    }

    //see NR p 784, fit the function y = a + bx
    Fitab fitting(X,Y,sigma);

    //return error in mu (called 'b' in NR)
    sigma_mu = fitting.sigb;

    //return the slope (by reference)
    mu = fitting.b;
  }
}


inline void Save::computeError(const matrixD_t& store,
                               const vectorD_t& msd, vectorD_t& sigma){
  //save and return (by reference&) the error to the msd-vector
  sigma.assign(noSamplingTimes_,0);

  //don't use the "ensemble_" variable as iteration limit, since we will
  //remove trajectories (ensembles) in the jackknife implementation.
  const int ensembles = store.size();

  for(int j = 0; j < noSamplingTimes_; j++){
    for(int i = 0; i < ensembles; i++)
      sigma[j] += pow (store[i][j] - msd[j] ,2);
    sigma[j] = sqrt( sigma[j] / (1.0 * ensembles * (ensembles-1)) );
  }
}


void printHmatrix(const Save::matrixD_t& matrix,
                  std::string outputFileName){

  std::cout << "printing matrix to file: " << outputFileName;
  std::ofstream dumpMatrix(outputFileName.c_str());   //open file

  int max_i = matrix.size();     //noSamplingTimes_ if H-matrix
  int max_j = matrix[0].size();  //noSamplingTimes_ for all matrices I think

  for(int i = 0; i < max_i; i++){
    for(int j = 0; j < max_j; j++)
      dumpMatrix  << matrix[i][j] << "\t";
    dumpMatrix << std::endl;
  }
  dumpMatrix.close();
  std::cout << "\t \033[92m DONE \033[0m" << std::endl;
}



void Save::computeHmatrix(const matrixD_t& trajectories, matrixD_t& matrix){
  // std::cout << "\033[1m Computing H-matrix (method 1): \033[0m"
  //           << std::endl;
  computeVariance(trajectories, matrix, false);

  double invM = 1.0 / trajectories.size();
  size_t n = trajectories[0].size();        //dynamic size, since jackknifing

# pragma omp parallel for default(shared)
  for(size_t i = 0; i < n; i++){
    for(size_t j = 0; j < n; j++)
      matrix[i][j] = invM * matrix[i][j];
  }
}


void Save::computeHmatrix2(matrixD_t& matrix, bool isMeanZero){

  matrixD_t e_x(noSamplingTimes_, vectorD_t(noSamplingTimes_, 0));
  matrixD_t e_y(noSamplingTimes_, vectorD_t(noSamplingTimes_, 0));
  matrixD_t e_z(noSamplingTimes_, vectorD_t(noSamplingTimes_, 0));

  std::cout << "\033[1m Computing H-matrix method 2 (each dim. separately): \033[0m"
            << std::endl;

  computeVariance(store_dx_, e_x, isMeanZero);
  computeVariance(store_dy_, e_y, isMeanZero);
  computeVariance(store_dz_, e_z, isMeanZero);

  double invM = 1.0 / noEnsembles_;
# pragma omp parallel for default(shared)
  for(int i = 0; i < noSamplingTimes_; i++){
    for(int j = 0; j < noSamplingTimes_; j++)
      matrix[i][j] = invM * (pow(e_x[i][j], 2) + pow(e_y[i][j], 2) + pow(e_z[i][j], 2));
  }
}


void makeBandDiagonal(Save::matrixD_t &matrix, double factor){
  //Utility function, used in computeHmatrix3.  Find how far out from
  //the main diagonal we should go i.e. how wide the band should
  //be. Everything else is 0.

  //Find maximum element in matrix
  double max = 0;
  for(size_t i = 0; i < matrix.size(); ++i){
    double max_tmp = *std::max_element(matrix[i].begin(), matrix[i].end());
    if(max_tmp > max) max = max_tmp;
  }

  int N = matrix.size();
  double threshold = factor * max;

  int k;  // distance from main diagonal
  for(k = 0; k < int(N/2); ++k){
    double sum = 0;
    for(int i = k; i < N-k; ++i)
      sum += matrix[i-k][i+k]/matrix[i][i];
    if(sum/(N-k) < threshold)
      break;
  }

  std::cout << "k: " << k << std::endl;

# pragma omp parallel for
  for(int i = 0; i < N; ++i){
    for(int j = k+i; j < N; ++j){
      matrix[i][j] = 0;  //Matrix is symmetric.
      matrix[j][i] = 0;
    }
  }
}


void Save::computeBmatrix(const matrixI_t &store_di, matrixD_t &B_i,
                          bool isMeanZero){
  //transform normal coordinates in dimension i to B_i matrix.

  //keep displacement from previous timestep:
  matrixI_t store_du_i(noEnsembles_, vectorI_t(noSamplingTimes_, 0));

  //Compute displacement from previous step:
# pragma omp parallel for
  for(int m = 0; m < noEnsembles_; ++m){
    //First step needs special love:
    store_du_i[m][0] = store_di[m][0] - 0;

    for(int t = 1; t < noSamplingTimes_; ++t)
      store_du_i[m][t] = store_di[m][t] - store_di[m][t-1];
  }
  computeVariance(store_du_i, B_i, isMeanZero);
}


void Save::computeHmatrix3(matrixD_t& matrix, std::string filename,
                           bool isAverageZero, bool isBandDiagonal){
  std::cout << "\033[1m Computing H-matrix ";
  if(isAverageZero && !isBandDiagonal)  std::cout << "(using method 3)";
  if(!isAverageZero && !isBandDiagonal) std::cout << "(using method 3b)";
  if(isAverageZero && isBandDiagonal)   std::cout << "(using method 3c)";
  std::cout << "\033[0m" << std::endl;

  matrixD_t B_x(noSamplingTimes_, vectorD_t(noSamplingTimes_, 0));
  matrixD_t B_y(noSamplingTimes_, vectorD_t(noSamplingTimes_, 0));
  matrixD_t B_z(noSamplingTimes_, vectorD_t(noSamplingTimes_, 0));
  computeBmatrix(store_dx_, B_x, isAverageZero);
  computeBmatrix(store_dy_, B_y, isAverageZero);
  computeBmatrix(store_dz_, B_z, isAverageZero);  //no need for this one if not 3D.

  printHmatrix(B_x, filename + "_B_x");

  //initiate class, to print progress to screen,
  RemainingTime computationProgress1(noSamplingTimes_);

  if(isBandDiagonal){
    makeBandDiagonal(B_x, 0.01);
    makeBandDiagonal(B_x, 0.01);
    makeBandDiagonal(B_z, 0.01);

    printHmatrix(B_x, filename + "_B_x_band");  // for testing, XXX
  }

  double invM = 1.0 / noEnsembles_;
  matrix.assign(noSamplingTimes_, vectorD_t(noSamplingTimes_, 0));
  for(int i = 0; i < noSamplingTimes_; ++i){
    computationProgress1.printProgress(i);
#   pragma omp parallel for shared(matrix)
    for(int j = 0; j < noSamplingTimes_; ++j){
      double sumX = 0;
      double sumY = 0;
      double sumZ = 0;
      for(int k = 0; k <= i; ++k){
        for(int l = 0; l <= j; ++l){
          sumX += B_x[k][l];
          sumY += B_y[k][l];
          sumZ += B_z[k][l];
        }
      }
      matrix[i][j] = 2 * invM * (pow(sumX, 2) + pow(sumY, 2) +  pow(sumZ, 2));
    }
  }
}




//==================================================
//Here I try the Shrinkage approach of Schäfer and Strimmer 2005

double varvar(const Save::matrixD_t &trajectories, const Save::vectorD_t &r2_mu,
              const Save::matrixD_t &h, size_t i, size_t j){
  //Copmute variance of variance, see appendix in Schäfer and Strimmer 2005

  size_t m = trajectories.size();
  double factor = (m - 1.0) / m;

  double var = 0;
  for(size_t k = 0; k < m; ++k)
    var += pow((trajectories[k][i] - r2_mu[i]) * (trajectories[k][j] - r2_mu[j]) - factor * h[i][j], 2);

  return m * var / pow(m-1, 3);
}


void Save::computeShrinkage(const matrixD_t &trajectories, matrixD_t &H_star){
  //Compute a better estimator, h_star for the covariance matrix h;
  //will be better than the the unbiased one.

  size_t M = trajectories.size();
  size_t N = trajectories[0].size();
  matrixD_t h(N, vectorD_t(N, 0));

  //comute "normal" covariance. 1/(M-1) sum_k (...)_i * (...)_j
  computeVariance(trajectories, h, false);

  //compute lambda (shrinkage factor)
  double lambda = 0, denomSum = 0;

  vectorD_t r2_mu;
  computeMean<double>(trajectories, r2_mu);

  //RemainingTime printToScreen(N);                 //non-important nice-to-have simulation info
  for(size_t i = 0; i < N; ++i){
    //printToScreen.printProgress(i);
    for(size_t j = 0; j < N; ++j)
      if(i != j){
        lambda += varvar(trajectories, r2_mu, h, i, j);
        denomSum += pow(h[i][j], 2);
      }
  }

  lambda = lambda / denomSum;
  assert(denomSum > 0);

  //NOTE: if N = 1, then h = [0] --> denomSum = 0 --> lambda = -nan,
  //and after this H_star will be fucked up.

  //compute improved covariance matrix, and save to H_star:
  H_star.assign(N, vectorD_t(N, 0));
  double invM = 1.0 / M;
  for(size_t i = 0; i < N; ++i){
    for(size_t j = 0; j < N; ++j){
      if(i != j)
        H_star[i][j] = invM * h[i][j] * (1 - lambda);
      else
        H_star[i][j] = invM * (h[i][j] * (1 - lambda) + lambda * h[i][j]);
    }
  }
}





//==================================================

/*

// Here's an attempt to combine the boostrapp and bootknife functions to
// have more code in common, and then take a procedure as input argument,
// that will be the difference between the two functions.


void Save::Boots(std::string name, int noOfRuns,
                 double k_t, std::string head = ""){
  //Generate synthetic values, from current trajectories in store_dr2_
  //to get new values, by Quick-and-dirty Monte Carlo Bootstrapping.

  //Needs all the information for all trajectories saved in store_dr2_,
  //which will not be saved in lowMem mode.
  if(isLowMem_){
    std::cout << "Error: Can not use bootstrap when in low memory mode"
              << std::endl;
  }
  else{
    std::cout << "bootstrapping!" << std::endl;

    matrixD_t store_dr2(noEnsembles_, vectorD_t(noSamplingTimes_,0));

    Ran randomNumber(SEED_BOOTSTRAP);

    //This part lets us choose different starting times
    // k_t * t should have these values, corresponding to errors:
    // 1.010  = 2.45
    // 1.0050 = 7.34
    // 1.0020 = 25.40
    // 1.0010 = 58.34
    // 1.0005 = 140
    const double tempArray[] = {0, 2.45, 7.34, 10, 20, 25.40, 30, 40, 50,
                                58.34, 70, 80, 90, 100, 110, 120, 130, 140};

    size_t sizeOfTempArray = (sizeof tempArray)/(sizeof tempArray[0]);

    //construct a vector with a number of starting times for the fitting
    vectorD_t startTimeForFitting;

    //now transform to time values:
    for(size_t i = 0; i < sizeOfTempArray; i++){
      float timeMark = (float) tempArray[i] / k_t;
      startTimeForFitting.push_back(timeMark);
    }

    int noReFittings = startTimeForFitting.size();

    matrixD_t mu(noOfRuns, vectorD_t(noReFittings,0));
    matrixD_t sigma_mu(noOfRuns, vectorD_t(noReFittings,0));

    //non-important nice-to-have simulation info.
    //this now makes save.cpp depend on classes.[h,cpp]
    RemainingTime progressToScreen(noOfRuns);

    //Loop to generate "noOfRuns" number of synthetic MSDs of
    //trajectories, each saved MSD corresponds to a simulation of
    //"noEnsembles_" number of individual trajectories (i.e. ensembles).
    for(int k = 0; k < noOfRuns; k++){
      progressToScreen.printProgress(k);

      //Randomly pick real trajectories into our virtual simulation
      for(int i = 0; i < noEnsembles_; i++){
        int r = (int) (noEnsembles_ * randomNumber.doub());
        store_dr2[i] = store_dr2_[r];
      }

      //Compute <R^2>, store result in r2_mu.
      vectorD_t r2_mu(noSamplingTimes_, 0);
      computeMean<double>(store_dr2, r2_mu);


      //call this prcedure, different depending on jackknife or
      //bootstrap:
      processData();

      //Bootstrap needs:
      k;
      r2_mu;
      conserveMemory;
      name;
      catanateOutputFiles;

      //Bootknife needs:
      k;
      r2_mu;
      mu;
      sigma_mu;
      store_dr2



    }

    if(!conserveMemory){
      //Since we did not dump each MSD-trajectory to file, we
      //compute the H-matrix and dump that to a single file.

      //Compute the mean of the Mean Sq. Disp. (don't square it!)
      //use noOfRuns instead of noEnsembles_ since might not be the same
      vectorD_t meanOfMsd(noSamplingTimes_,0);
      computeMean(catanateOutputFiles, meanOfMsd, 1);

      matrixD_t H;
      //computeHmatrix1(H, catanateOutput, meanOfMsd, 1); //OLD SYNTAX TODO
      computeHmatrix(H, catanateOutputFiles);

      //print H-matrix to file, with "_matrix_bootstrap" appended:
      printHmatrix(H, name + "_matrix_bootstrap");

      //Print bootstrapped mean of mean  to file:
      //-----------------------------------------
      std::ofstream msdBootout((name + "_bootstrap").c_str());
      msdBootout << head;
      for(int i = 0; i < noSamplingTimes_; i++)
        msdBootout  << samplingTime_[i] << "\t" << meanOfMsd[i] << std::endl;
      msdBootout.close();
    }
  }
}




void Save::computeBootstrap(std::string name, int noOfRuns,
                            double k_t, std::string head = ""){

    //Either print "noOfRuns" number of synthetic output-files to
    //hard drive, each file corresponds one MSD of a synthetic
    //simulation, or save all MSDs in mem, then save to one file.
    bool conserveMemory = BOOTSTRAP_CONSERVE_MEMORY;
    matrixD_t catanateOutputFiles;
    if(!conserveMemory)
      catanateOutputFiles.assign(noOfRuns,
                                 vectorD_t(noSamplingTimes_,0));


    boot(bootstrap_procedure, );


}


void Save::bootstrap_procedure(int k){
  if(conserveMemory){
    //Print data to file (for processing in MatLab):
    std::stringstream tempName;
    //append run-iterator (k) to filename:
    tempName << name << k;
    std::ofstream bootOutput((tempName.str()).c_str());

    bootOutput << head;
    for(int i = 0; i < noSamplingTimes_; i++)
      bootOutput  << samplingTime_[i] << "\t" << r2_mu[i] << std::endl;
    bootOutput.close();
  }
  else{
    //store all bootstraped data in a single matrix...
    catanateOutputFiles[k] = r2_mu;
    //...and print it later, to a single file instead.
  }
}

void Save::jackknife_procedure(int k){

  //return error in each sampling-point by dr2_err, needed when
  //fitting mu. New for each bootstrapped "simulation"
  vectorD_t dr2_err;
  computeError(store_dr2, r2_mu, dr2_err);

  //Patch to modify/test minTime for fitting. get mu & sigma_mu.
  for(int i = 0; i < noReFittings; i++){
    minTime_ = startTimeForFitting[i];
    computeSlope(r2_mu, dr2_err, mu[k][i], sigma_mu[k][i]);
  }
}

*/
