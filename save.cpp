#include <iostream>        //for cout, among others
#include <fstream>         //read write to files
#include <cmath>           //gives sqrt, pow, fabs,
#include <vector>          //STL vector
#include <string>          //STL string
#include <algorithm>       //for min/max element finding in histogram
#include <cassert>         //for of with assert()

#include "nr/nr3.h"
#include "nr/ludcmp.h"       //only depends on nr3.h
#include "nr/svd.h"          //singular value decomposition
#include "nr/ran.h"          //only needed for my bootstrap-function

#include "nr/gamma.h"        //needed by fitab.h (see NR p. 784)
#include "nr/incgammabeta.h" //needed by fitab.h
#include "nr/fitab.h"        //needed by computeSlope()

#include "classes.h"         //to print bootstrap progress to screen
#include "auxiliary.h"       //for the printError()-function
#include "save.h"
#include "global.h"

//#include "nr/gaussj.h"      //needed by fitmrq.h
//#include "nr/fitmrq.h"      //not yet in use...

const bool USE_NR_FITTING = false;    //if true: y = a + bx, if false: y = bx
const float SEED_BOOTSTRAP = 6854;    //default seed: 6854

//If true, we print each bootstrapped trajectory to separate file:
//out.dat0, out.dat1, etc, (assuming -w out.dat). Else (if "false") we compute
//the H-matrix directly, and output to out.dat_matrix
const bool BOOTSTRAP_CONSERVE_MEMORY = false;

//If true, also use the hybrid "bootknife" method to compute
//slope, and error, and write to separate file "out.dat_bootknife"
const bool BOOTSTRAP_BOOTKNIFE = true;

//what to print to 4th column of the output file.
//Note: Mutually exclusive, only one can be true.
const bool COMPUTE_CORRELATION = false;
const bool COPMUTE_PEARSON = false;


Save::Save(std::vector<double> samplingTime,
           int ensembles, bool lowMem = false){

  samplingTime_ = samplingTime;

  //when using computeSlope, discard any sampling points lower than this.
  minTime_ = 0;

  //number of realizations / trajectories.
  noEnsembles_ = ensembles;

  //number of elements
  noSamplingTimes_ = samplingTime_.size();

  //use the memory conserving function
  isLowMem_ = lowMem;

  //initiate the vectors for the summation in store()
  //only needed if we use the memory-conserving algorithm
  if(isLowMem_){
    dr4_err_.assign(noSamplingTimes_,0);
    dx4_err_.assign(noSamplingTimes_,0);
    dy4_err_.assign(noSamplingTimes_,0);
    dz4_err_.assign(noSamplingTimes_,0);

    r2_mu_.assign(noSamplingTimes_,0);
    x2_mu_.assign(noSamplingTimes_,0);
    y2_mu_.assign(noSamplingTimes_,0);
    z2_mu_.assign(noSamplingTimes_,0);
  }
  else{
    //Allocate size in memory for the matrix to store every
    //individual trajectory, will use push_back later, so don't "assign()"
    store_dx_.reserve(noEnsembles_);
    store_dy_.reserve(noEnsembles_);
    store_dz_.reserve(noEnsembles_);
    store_dr_.reserve(noEnsembles_);
  }
}


//this is needed in the bootknife when we try fitting to many different
//starting times. (since we want t*k_t to be a constant, we need to know k_t)
void Save::setJumprate(double k_t){
  k_t_ = k_t;
}


void Save::computePearsonCoefficient(){
  //This function calculates the Pearson coefficients (see Numerical
  //recipes ch 14), to get the correlation between two sampling points
  //since if the tagged particle jumps once/sec sampling the system
  //10/s will not yield any new information, and it will mess up the
  //error estimation that assumes the sampling points are uncorrelated
  //and with a Gaussian distribution.

  //perhaps check: Physical Review D vol 34, 2469 (1986), sec. D

  if(isLowMem_)                 //need data to work on
    printError("Can not compute correlation, no available data.");

  if(r2_mu_.empty())          //need data to work on
    printError("Must compute MSD first! (ie. r2_mu_)");

  if(!fourthColumn_.empty())  //this is to where this function saves data.
    printError("Pearson Coefficient is replacing data in 4th column.");

  std::vector<double> pearsonCorrelation(noSamplingTimes_,0);

  //compute the correlation between two points, i, j, for i = j+1
  std::vector<double> h(noSamplingTimes_,0);

  for(int i = 1; i < noSamplingTimes_; i++){
    for(int j = 0; j < noEnsembles_; j++){
      h[i] += 1.0 / (noEnsembles_ - 1) *
        (pow(store_dr_[j][i],2) - r2_mu_[i])*
        (pow(store_dr_[j][i-1],2) - r2_mu_[i-1]);
    }
  }

  // pearsonCorrelation[i] = h_ij / (sigma_i*sigma_j)
  // sigma_r = r2_err_ * sqrt(ensembles);
  for(int i = 1; i < noSamplingTimes_; i++){
    if(r2_err_[i-1] != 0)
      pearsonCorrelation[i] = (double) h[i] /
        (r2_err_[i] * r2_err_[i-1] * noEnsembles_);

    //NOTE: error in first value is 0, so the correlation between value 0 & 1
    //will be inf (nan)! thus the if-statement

    //Now save the inf. for printing to the output-file.
    fourthColumn_ = pearsonCorrelation;
  }
}


void Save::computeCorrelation(){
//This function computes the H-matrix: XXX and saves z to the 4th
//column which gets printed to outputfile elsewhere.


  //compute entire matrix, to check the correlation between every value
  std::vector<std::vector<double> > H;

  //note! does not divide with the factor (1/ensembles) so I do that below
  computeHmatrix(store_dr_, r2_mu_, H);

  //In order to use NR implementations to test the matrix, and perform
  //LU-decomposition, we transform STL-matrix 'H' to NR-matrix
  NRmatrix<double> h(noSamplingTimes_, noSamplingTimes_, 0.0);

  for(int i = 0; i < noSamplingTimes_; i++){
    for(int j = 0; j < noSamplingTimes_; j++){

      //computeHmatrix only divides with (ensembles-1):
      H[i][j] = H[i][j] / (double) noEnsembles_;

      //now copy to NR matrix.
      h[i][j] = H[i][j];
    }
  }


  //Print H-matrix to file (now that we divided it with M(M-1), which is
  //done in the for-loop above (M), and in the function "computeHmatrix()"
  //where we divide it with (M-1), where M = noEnsembles_
  if(false){
    //This will collide with H-matrix computed by computeBootstrap,
    //which will also print a out.dat_matrix file, which will overwrite
    //this one, unless the "inputName" is changed.
    printHmatrix(H,"out.dat_Hmatrix_test");
  }

  if(false){
    //Now test the integrity of the matrix:
    SVD testMatrix(h);
    double conditionalNumber = testMatrix.inv_condition();
    std::cout << "conditional number is: " << setprecision(5) << scientific
              << conditionalNumber << std::endl;

    //OBSERVE: For some retarded reason it rounds off the answer to the
    //conditional number! i.e gives 0.0, instead of 3.5E-5 etc.
    //therefore, must use: << scientific <<
  }
  if(false){
    //just a small test of NR conditional number.
    NRmatrix<double> hm(2,2,0.0);
    hm[0][0] = 1.0e-16;   //epsilon
    hm[0][1] = 1.0;
    hm[1][0] = 1.0;
    hm[1][1] = 1.0;
    SVD testMatrix2(hm);
    double conditionalNumber2 = testMatrix2.inv_condition();
    std::cout << fixed << setprecision(5) << "Test-conditional number is : "
              << conditionalNumber2 << std::endl;
  }

  //Solve Equation, with LU-Decomposition
  //-------------------------------------
  NRvector<double> b(noSamplingTimes_),z(noSamplingTimes_);

  for(int i = 0; i < noSamplingTimes_; i++)
    b[i] = samplingTime_[i];

  LUdcmp LUdecomposition(h);
  LUdecomposition.solve(b,z);

  //copy z to the class-private vector, that we will print as the
  //fourth column in the output-file. Transfer a NR vector to std::
  fourthColumn_.assign(noSamplingTimes_,0);
  for (int i = 0; i < noSamplingTimes_; i++)
    fourthColumn_[i] = z[i];
}



void Save::store(const std::vector<int>& dx, const std::vector<int>& dy,
                 const std::vector<int>& dz, const std::vector<double>& dr){
  //This function stores all the trajectories, so we can process the data,
  //like error bars in computeStdErr().

  /*
    Error in each point i is (in LaTeX notation):
    (a[i])^2 = 1/N \sum_i^N (a_i - <a>)^2 = [1/N \sum_i^N a[i]^2] - [<a>]^2

    <a>= 1/N \sum_i^N a[i]   (i.e. average value in that point)
    N = Number of ensembles

    The original way to do this is save the value in each point a[i] for
    each ensemble, but that leaves us with a giant matrix if we have many
    ensembles N, and many data points D, (N*D ~> 10^9) ~ 1Gb of mem usage
    (10^9/2^30 ~ 1Gb). Therefore with the switch 'lowMem=true' we do the
    left hand side of the equation above, and then we don't need a matrix
    but just two vectors, for each dimension.
  */

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
    //if don't use the lowMem-algorithm:
    //this gives dx[ensembles][number of values]
    store_dx_.push_back(dx);
    store_dy_.push_back(dy);
    store_dz_.push_back(dz);
    store_dr_.push_back(dr);
  }
}














void Save::save(std::string name, std::string head){
  //calculates the MSD from the stored positions of the tracer,
  //and prints to file.

  //sore averaged displacement <x>, <y>, <z>, <r>
  //mostly used to check that is is ~ 0.
  std::vector<double> x_mu, y_mu, z_mu, r_mu;

  //hard-coded switch to print x,y,z, x^2,y^2,z^2
  //Not needed now (unless testing), since system is isotropic,
  bool printSecondFile = PRINT_SECOND_FILE;
  std::ofstream secondary;


  if(isLowMem_){
    //use the memory conserving version instead. This one does not
    //use the dx[][]-matrix, and in that way we save lots of mem.
    computeStdErrLowMem();
  }
  else{
    if(printSecondFile){
      //Check that <x,y,z>=0, just for fun/testing purpose.

      std::string secondFileName = name + "_txyz";
      secondary.open(secondFileName.c_str());

      //return mean, (un-squared), return through reference x_mu...
      computeMean(store_dx_,store_dy_,store_dz_, x_mu, y_mu, z_mu);
    }

    //Compute <R^2>:
    x2_mu_.assign(noSamplingTimes_,0);
    y2_mu_.assign(noSamplingTimes_,0);
    z2_mu_.assign(noSamplingTimes_,0);
    r2_mu_.assign(noSamplingTimes_,0);

    for (int j = 0; j < noSamplingTimes_; j++){
      for (int i = 0; i < noEnsembles_;i++){
        x2_mu_[j] += 1.0 * pow(store_dx_[i][j],2) / noEnsembles_;
        y2_mu_[j] += 1.0 * pow(store_dy_[i][j],2) / noEnsembles_;
        z2_mu_[j] += 1.0 * pow(store_dz_[i][j],2) / noEnsembles_;
        r2_mu_[j] += 1.0 * pow(store_dr_[i][j],2) / noEnsembles_;
      }
    }

    //Calculate standard error, to generate error-bars
    //this function must come after x2_mu_...r2_mu_ has been
    //computed.
    computeStdErr();

    if(COPMUTE_PEARSON)
      computePearsonCoefficient();
    if(COMPUTE_CORRELATION)
      computeCorrelation();
  }


  //print data to file:
  //-------------------
  std::ofstream primary;
  primary.open(name.c_str());
  primary << head;
  std::cout << head;

  if(printSecondFile){
    //print comment to file:
    secondary << "# <time>  <x_mu>  <y_mu>  <z_mu>  <time>   "
              << "<x2_mu>  <y2_mu>  <z2_mu>  "
              << "<x2_err>  <y2_err>  <z2_err>" << std::endl;
  }

  for (int i = 0; i < noSamplingTimes_; i++){
    primary << samplingTime_[i] << "\t" << r2_mu_[i] << "\t" << r2_err_[i];
    if(!fourthColumn_.empty())
      primary << "\t" << fourthColumn_[i] << std::endl;

    if(printSecondFile)
      secondary << samplingTime_[i] << "\t"
                << x_mu[i] << "\t" << y_mu[i] << "\t" << z_mu[i] <<"\t"
                << x2_mu_[i] <<"\t"<< y2_mu_[i] <<"\t"<< z2_mu_[i] <<"\t"
                << x2_err_[i]<<"\t"<<y2_err_[i]<<"\t"<<z2_err_[i]<<std::endl;
  }
  primary.close();

  //Print file with unsquared MSD for x,y,z
  if(printSecondFile) secondary.close();
}


void Save::computeStdErr(void){
  //Calculate the standard deviation, to get error-bars.

  //Should not be using this function if lowMem = true!
  assert(!isLowMem_);

  x2_err_.assign(noSamplingTimes_,0);
  y2_err_.assign(noSamplingTimes_,0);
  z2_err_.assign(noSamplingTimes_,0);
  r2_err_.assign(noSamplingTimes_,0);

  for (int j = 0; j < noSamplingTimes_; j++){
    for (int i = 0; i < noEnsembles_; i++){
      x2_err_[j] += pow (pow(store_dx_[i][j],2) - x2_mu_[j] ,2);
      y2_err_[j] += pow (pow(store_dy_[i][j],2) - y2_mu_[j] ,2);
      z2_err_[j] += pow (pow(store_dz_[i][j],2) - z2_mu_[j] ,2);
      r2_err_[j] += pow (pow(store_dr_[i][j],2) - r2_mu_[j] ,2);
    }//( "a += b" is the same as "a = a + b")
  }

  for (int j = 0; j < noSamplingTimes_; j++){
    x2_err_[j] = sqrt( x2_err_[j] / (1.0*noEnsembles_*(noEnsembles_ - 1)) );
    y2_err_[j] = sqrt( y2_err_[j] / (1.0*noEnsembles_*(noEnsembles_ - 1)) );
    z2_err_[j] = sqrt( z2_err_[j] / (1.0*noEnsembles_*(noEnsembles_ - 1)) );
    r2_err_[j] = sqrt( r2_err_[j] / (1.0*noEnsembles_*(noEnsembles_ - 1)) );
  } //for E = 1 --> inf !


  if(false){
    //Testing the Gaussian distribution of the sampling points:
    float SUM68 = 0;
    for (int j = 0; j < noSamplingTimes_; j++){    //forgot about his:
      double top = x2_mu_[j] +0.5*x2_err_[j] * sqrt(noEnsembles_-1)*2;
      double bot = x2_mu_[j] -0.5*x2_err_[j] * sqrt(noEnsembles_-1)*2;
      for (int i=0; i < noEnsembles_; i++){
        if (bot <= pow(store_dx_[i][j],2) &&
            top >= pow(store_dx_[i][j],2) ){
          bool test = false;
          if (test){ //TEST-function
            std::cout <<"lower limit: "<< bot <<std::endl;
            std::cout << "point is: "<< pow(store_dx_[i][j],2) <<std::endl;
            std::cout << "upper limit: "<< top <<std::endl;
          }
          SUM68++;
          if (test) std::cout << "SUM68 = " << SUM68 << std::endl;
        }
      }
    }
    std::cout << "68,2 % (?)="
              << 1.0 * SUM68/(noEnsembles_ * noSamplingTimes_) << std::endl;
  }

}


void Save::computeStdErrLowMem(){
  //compute standard error, if using lowMem-algorithm.

  assert(isLowMem_);

  //Divide the summed up displacement to get the average:
  for (int i = 0; i< noSamplingTimes_; i++){
    x2_mu_[i] = x2_mu_[i] / (noEnsembles_ * 1.0);
    y2_mu_[i] = y2_mu_[i] / (noEnsembles_ * 1.0);
    z2_mu_[i] = z2_mu_[i] / (noEnsembles_ * 1.0);
    r2_mu_[i] = r2_mu_[i] / (noEnsembles_ * 1.0);
  }

  x2_err_.assign(noSamplingTimes_,0);
  y2_err_.assign(noSamplingTimes_,0);
  z2_err_.assign(noSamplingTimes_,0);
  r2_err_.assign(noSamplingTimes_,0);

  //Following the notation in the comments in Lattice::store():
  // [1/N sum_i^N a_i^2] - [<a>^2], (a=dr^2)
  for(int i = 0; i < noSamplingTimes_; i++){
    x2_err_[i] = dx4_err_[i] / (1.0 * noEnsembles_) - pow( x2_mu_[i] ,2);
    y2_err_[i] = dy4_err_[i] / (1.0 * noEnsembles_) - pow( y2_mu_[i] ,2);
    z2_err_[i] = dz4_err_[i] / (1.0 * noEnsembles_) - pow( z2_mu_[i] ,2);
    r2_err_[i] = dr4_err_[i] / (1.0 * noEnsembles_) - pow( r2_mu_[i] ,2);
  }

  for(int i =0; i < noSamplingTimes_; i++){
    x2_err_[i] = sqrt(x2_err_[i] / (1.0 * (noEnsembles_ - 1) ));
    y2_err_[i] = sqrt(y2_err_[i] / (1.0 * (noEnsembles_ - 1) ));
    z2_err_[i] = sqrt(z2_err_[i] / (1.0 * (noEnsembles_ - 1) ));
    r2_err_[i] = sqrt(r2_err_[i] / (1.0 * (noEnsembles_ - 1) ));
  }
}





void Save::computeDistribution(std::string fileName){
  //This function is for getting P(t,x), ie. the probability
  //distribution of finding a particle at r^2, at time t.

  if(isLowMem_){
    std::cout << "lowMem = ON, therefore no information stored!\n"
              << "can not print cluster-information" << std::endl;
  }
  else{
    const int noBins = 300;    //number of bins
    std::vector<std::vector<float> > printToFileR2;
    std::vector<std::vector<float> > printToFileX2;
    std::vector<std::vector<float> > printToFileY2;

    //find the smallest/biggest value in store_dr_[*][t]
    //-----------------------------------------------
    int samplingTime = noSamplingTimes_ -1; //pick sampling time

    std::vector<double> tempVectorR2(noEnsembles_,0);
    std::vector<double> tempVectorX2(noEnsembles_,0);
    std::vector<double> tempVectorY2(noEnsembles_,0);

    for (int i = 0; i < noEnsembles_; i++){
      //(I don't know of a better way to do this)
      tempVectorR2[i] = pow(store_dr_[i][samplingTime],2);
      tempVectorX2[i] = pow(store_dx_[i][samplingTime],2);
      tempVectorY2[i] = pow(store_dy_[i][samplingTime],2);
    }

    double minR2 = *std::min_element(tempVectorR2.begin(), tempVectorR2.end());
    double maxR2 = *std::max_element(tempVectorR2.begin(), tempVectorR2.end());
    double minX2 = *std::min_element(tempVectorX2.begin(), tempVectorX2.end());
    double maxX2 = *std::max_element(tempVectorX2.begin(), tempVectorX2.end());
    double minY2 = *std::min_element(tempVectorY2.begin(), tempVectorY2.end());
    double maxY2 = *std::max_element(tempVectorY2.begin(), tempVectorY2.end());

    double binStepR2 = (double) (maxR2 - minR2)/noBins;
    double binStepX2 = (double) (maxX2 - minX2)/noBins;
    double binStepY2 = (double) (maxY2 - minY2)/noBins;

    // if (binStepX < 0) binStepX = -1 * binStepX;
    // if (binStepY < 0) binStepY = -1 * binStepY;

    std::cout << "binStep R2: " << binStepR2 << std::endl;
    std::cout << "binStep X2: " << binStepX2 << std::endl;
    std::cout << "binStep Y2: " << binStepY2 << std::endl;

    std::vector<float> probAxisR2(noBins,0);    //Probability count
    std::vector<float> binAxisR2(noBins,0);     //binning-axis
    std::vector<float> probAxisX2(noBins,0);    //Probability count
    std::vector<float> binAxisX2(noBins,0);     //binning-axis
    std::vector<float> probAxisY2(noBins,0);    //Probability count
    std::vector<float> binAxisY2(noBins,0);     //binning-axis

    std::cout << "min/max R2: " << minR2 << " " << maxR2 << std::endl;
    std::cout << "min/max X2: " << minX2 << " " << maxX2 << std::endl;
    std::cout << "min/max Y2: " << minY2 << " " << maxY2 << std::endl;

    //Build the X-axis scale
    float tempR2 = minR2;
    float tempX2 = minX2;
    float tempY2 = minY2;
    for(int i = 0; i < noBins; i++){
      binAxisR2[i] = tempR2;
      binAxisX2[i] = tempX2;
      binAxisY2[i] = tempY2;

      tempR2 = tempR2 + binStepR2;
      tempX2 = tempX2 + binStepX2;
      tempY2 = tempY2 + binStepY2;
    }

    //const double PI = 3.14159265;
    double r2, x2, y2;                  //displacement

    std::cout << "Calculating P(t,x) for histogram... ";

    //Analyze all sampling-times, or just one? (start with one)
    const int startTimeElement = samplingTime; // 0;
    const int endTime = startTimeElement + 1;  // noSamplingTimes_;

    for(int t = startTimeElement; t < endTime ; t++){ //loop over time
      for (int j = 0; j < noBins; j++){               //loop over bin-cells
        for(int e = 0; e < noEnsembles_; e++){          //loop over "ensembles"
          r2 = pow(store_dr_[e][t],2);
          x2 = pow(store_dx_[e][t],2);
          y2 = pow(store_dy_[e][t],2);

          //why 2*PI ? r? (November 2010)
          //if(binAxis[j+1] > r2 && r2 >= binAxis[j])
          //probAxis[j] = probAxis[j] + 1.0/(noEnsembles_ * 2 * PI * ((r2 + binStep)/2));

          if(binAxisR2[j+1] > r2 && r2 >= binAxisR2[j])
            probAxisR2[j] = probAxisR2[j] + 1.0 / noEnsembles_;

          if(binAxisX2[j+1] > x2 && x2 >= binAxisX2[j])
            probAxisX2[j] = probAxisX2[j] + 1.0 / noEnsembles_;

          if(binAxisY2[j+1] > y2 && y2 >= binAxisY2[j])
            probAxisY2[j] = probAxisY2[j] + 1.0 / noEnsembles_;

        }
      }
      printToFileR2.push_back(probAxisR2);
      printToFileX2.push_back(probAxisX2);
      printToFileY2.push_back(probAxisY2);

      probAxisR2.assign(noBins,0);  //clear and reassign vectors:
      probAxisX2.assign(noBins,0);
      probAxisY2.assign(noBins,0);

    }

    //Print to file:
    //-------------
    std::string histogramNameR2 = fileName + "_histogram";
    std::string histogramNameX2 = fileName + "_histogramX2";
    std::string histogramNameY2 = fileName + "_histogramY2";

    std::ofstream histR2(histogramNameR2.c_str());
    std::ofstream histX2(histogramNameX2.c_str());
    std::ofstream histY2(histogramNameY2.c_str());

    std::cout << " ...Printing histogram... " << std::endl;
    int noSamplingTimes = printToFileR2.size();    //or noSamplingTimes_
    for(int i = 0; i < noSamplingTimes; i++){      //loop over time
      for(int j = 0; j < noBins; j++){             //loop over bins
        if (printToFileR2[i][j] != 0){
          histR2 << binAxisR2[j] + binStepR2/2 //<< "\t" << samplingTime_[i]
                 << "\t" << printToFileR2[i][j] << std::endl;
        }

        if (printToFileX2[i][j] != 0){
          histX2 << binAxisX2[j] + binStepX2/2 //<< "\t" << samplingTime_[i]
                 << "\t" << printToFileX2[i][j] << std::endl;
        }

        if (printToFileY2[i][j] != 0){
          histY2 << binAxisY2[j] + binStepY2/2 //<< "\t" << samplingTime_[i]
                 << "\t" << printToFileY2[i][j] << std::endl;
        }
      }
      //separate data block with new (empty) line:
      histR2 << std::endl;
      histX2 << std::endl;
      histY2 << std::endl;

    }
    histR2.close();
    histX2.close();
    histY2.close();

    std::cout << "Printed output to: " << histogramNameR2 << std::endl;


    //TEST, alternative way to get the MSD, by using the distribution
    //this is a good way to test it.
    double MSD = 0;
    for (int i = 0; i < noBins; i++)
      MSD += binAxisR2[i] * printToFileR2[0][i];
    std::cout << "MSD calculated from distribution: " << MSD << std::endl;
  }
}








void Save::computeBootstrap(std::string name, int noOfRuns,
                            std::string head = ""){
  //Generate synthetic values, from current trajectories in store_dr_
  //to get new values, by Quick-and-dirty Monte Carlo Bootstrapping.
  //(TODO: actually I should re-use the previous code for this, since
  //there is an overlap, it would look better.)

  //Needs all the information for all trajectories saved in store_dr_,
  //which will not be saved in lowMem mode.
  if(isLowMem_){
    std::cout << "Error: Can not use bootstrap when in low memory mode"
              << std::endl;
  }
  else{

    std::cout << "bootstrapping!" << std::endl;

    std::vector<double> r2_mu;
    std::vector<std::vector<double> > store_dr;

    //March 2011: to check that MSD := [x - <x>]^2
    //std::vector<std::vector<int> > store_dx, store_dy, store_dz;

    //Assign store_dr (local variable) with zeros
    store_dr.assign(noEnsembles_,std::vector<double>(noSamplingTimes_,0));
    // store_dx.assign(noEnsembles_,std::vector<int>(noSamplingTimes_,0));
    // store_dy.assign(noEnsembles_,std::vector<int>(noSamplingTimes_,0));
    // store_dz.assign(noEnsembles_,std::vector<int>(noSamplingTimes_,0));

    Ran randomNumber(SEED_BOOTSTRAP);

    //Either print "noOfRuns" number of synthetic output-files to
    //hard drive, each file corresponds to the MSD of a synthetic
    //simulation, or save all MSDs in a single file. (easier for
    //matlab, but uses more memory
    bool conserveMemory = BOOTSTRAP_CONSERVE_MEMORY;
    std::vector<std::vector<double> > catanateOutputFiles;
    if(!conserveMemory)
      catanateOutputFiles.assign(noOfRuns,
                                 std::vector<double>(noSamplingTimes_,0));

    //used in bootknife-loop, this does not look good!
    std::vector<double> dr2_err;
    bool bootknife = BOOTSTRAP_BOOTKNIFE;

    //construct a vector with a number of starting times for the fitting
    std::vector<double> startTimeForFitting;

    // k_t * t should have these values, corresponding to errors:
    // 1.010  = 2.45
    // 1.0050 = 7.34
    // 1.0020 = 25.40
    // 1.0010 = 58.34
    // 1.0005 = 140
    const double tempArray[] = {0, 2.45, 7.34, 10, 20, 25.40, 30, 40, 50,
                                58.34, 70, 80, 90, 100, 110, 120, 130, 140};

    size_t sizeOfTempArray = (sizeof tempArray)/(sizeof tempArray[0]);

    //now transform to time values:
    for (size_t i = 0; i < sizeOfTempArray; i++){
      float timeMark = (float) tempArray[i] / k_t_;
      startTimeForFitting.push_back(timeMark);
    }

    int noReFittings = startTimeForFitting.size();

    std::vector<std::vector<double> > mu, sigma_mu;
    mu.assign(noOfRuns,std::vector<double>(noReFittings,0));
    sigma_mu.assign(noOfRuns,std::vector<double>(noReFittings,0));

    //non-important nice-to-have simulation info.
    //this now makes save.cpp depend on classes.[h,cpp]
    RemainingTime progressToScreen(noOfRuns);

    //Loop to generate "noOfRuns" number of synthetic MSDs of
    //trajectories, each saved MSD corresponds to a simulation of
    //"noEnsembles_" number of individual trajectories (i.e. ensembles).
    for(int k = 0; k < noOfRuns; k++){

      progressToScreen.printProgress(k);

      //Randomly pick real trajectories into our virtual simulation
      //-----------------------------------------------------------
      int r;
      for (int i = 0; i < noEnsembles_; i++){
        r = (int) (noEnsembles_ * randomNumber.doub());
        store_dr[i] = store_dr_[r];

        //March 2011: adding these, to use MSD := [x-<x>]^2 instead.
        //(i.e Not assuming <x>~0)
        // store_dx[i] = store_dx_[r];
        // store_dy[i] = store_dy_[r];
        // store_dz[i] = store_dz_[r];
      }

      //TEST: to see if using MSD = [x-<x>]^2 rather than taking <x>=0.
      //std::vector<double> dx_mu, dy_mu, dz_mu;
      //computeMean(store_dx,store_dy,store_dz, dx_mu, dy_mu, dz_mu);


      //Compute <R^2>, store result in r2_mu.
      computeMSD(store_dr, r2_mu);

      if(bootknife){
        //this is an alternative "bootknife" version of bootstrap,
        //where we refrain from using H-matrices, and just fit a line
        //to current bootstrapped "simulation" (as if uncorrelated,
        //like in the jackknife-method), and then take all my slopes
        //and compute an average, (which must be done outside this
        //for(k)-loop).

        //return error in each sampling-point by dr2_err, needed when
        //fitting mu. New for each bootstrapped "simulation"
        computeError(store_dr,r2_mu,dr2_err);

        //Patch to modify/test minTime for fitting. get mu & sigma_mu.
        for(int i = 0; i < noReFittings; i++){
          minTime_ = startTimeForFitting[i];
          computeSlope(r2_mu, dr2_err, mu[k][i], sigma_mu[k][i]);
        }
      }

      if (conserveMemory){

        //Print data to file (for processing in MatLab):
        //---------------------------------------------
        std::stringstream tempName;
        //append run-iterator (k) to filename:
        tempName << name << k;
        std::ofstream bootOutput((tempName.str()).c_str());

        bootOutput << head;
        for (int i = 0; i < noSamplingTimes_; i++){
          bootOutput  << samplingTime_[i] << "\t" << r2_mu[i] << std::endl;
        }
        bootOutput.close();
      }
      else{
        //store all bootstraped data in a single matrix...
        catanateOutputFiles[k] = r2_mu;
        //...and print it later, to a single file instead.
      }
    } //(end of for(k)-loop


    if(bootknife){
      //continue on task done in previous if(bootknife)-statement

      //only needs to be a vector if we fit for several different
      //starting times.
      std::vector<double> mu_average(noReFittings,0);
      std::vector<double> sigma_mu_final(noReFittings,0);

      for (int i = 0; i < noReFittings; i++ ){
        for (int j = 0; j < noOfRuns; j++)
          mu_average[i] += mu[j][i];
        mu_average[i] = mu_average[i] /( (double) noOfRuns);
      }

      for (int i = 0; i < noReFittings; i++ ){
        for (int j = 0; j < noOfRuns; j++)
          sigma_mu_final[i] += pow(mu[j][i] - mu_average[i] ,2);
        sigma_mu_final[i] = sqrt(sigma_mu_final[i]
                                 /((double) (noOfRuns-1)));
      }
      //NOTE Am I not using sigma_mu? Why do I need this?


      //Print bootstrapped mean of mean to file:
      //-----------------------------------------
      std::ofstream bootknifeOut;
      bootknifeOut.open((name + "_bootknife").c_str());

      bootknifeOut  <<  "#t_{start},  mu,  and sig_mu " << std::endl;

      for (int i = 0; i < noReFittings; i++ ){
        bootknifeOut << startTimeForFitting[i] << "\t" << mu_average[i] << "\t"
                     << sigma_mu_final[i] << std::endl;
      }

      bootknifeOut.close();

    }


    if(!conserveMemory){
      //Since we did not dump each MSD-trajectory to file, we
      //compute the H-matrix and dump that to a single file.

      //Compute the mean of the Mean Sq. Disp. (don't square it!)
      //use noOfRuns instead of noEnsembles_ since might not be the same
      std::vector<double> meanOfMsd(noSamplingTimes_,0);
      for (int j = 0; j < noSamplingTimes_; j++){
        for (int i = 0; i < noOfRuns; i++){
          meanOfMsd[j] = meanOfMsd[j] + 1.0 * catanateOutputFiles[i][j]
            / noOfRuns;
        }
      }

      //Save H-matrix to this:
      std::vector<std::vector<double> > H;
      computeHmatrix(catanateOutputFiles,meanOfMsd ,H, 1);

      //print H-matrix to file, with "_matrix" appended:
      printHmatrix(H,name);

      //Print bootstrapped mean of mean  to file:
      //-----------------------------------------
      std::ofstream msdBootout;
      msdBootout.open((name + "_bootstrap").c_str());

      msdBootout << head;
      for (int i = 0; i < noSamplingTimes_; i++)
        msdBootout  << samplingTime_[i] << "\t" << meanOfMsd[i] << std::endl;

      msdBootout.close();

    }
  }
}


void Save::computeJackknife(std::string outName){
  //http://www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife/
  //use the jackknife method to compute the error from correlated data,
  //by successively (TODO) removing trajectories (i.e. "ensembles").

  //only works if we have all trajectories. (store_dr_-matrix)
  if(isLowMem_){
    std::cout << "ERROR: Can not use Jackknife when in lowmem-mode."
              << std::endl;
  }
  else{

    std::cout << "Jackknifing!" << std::endl;

    //store original chi^2 fit error to the slope:
    double sigma_mu;

    //fit parameter to the full set of (un-jackknifed) data
    double slope;
    computeSlope(r2_mu_, r2_err_, slope, sigma_mu);

    //store all the jackknifed versions of the fitted parameter
    std::vector<double> slopeJ;
    slopeJ.assign(noEnsembles_,0);

    //matrix of all simulation trajectories that we will "Jackknife"
    std::vector<std::vector<double> > store_dr_J;
    std::vector<double> dr2J, sigma;
    for (int i = 0; i < noEnsembles_; i++){
      store_dr_J = store_dr_;
      store_dr_J.erase(store_dr_J.begin()+i);

      //compute mean for each time based on the noEnsembles_-1
      //trajectories and return the MSD as dr2J
      computeMSD(store_dr_J, dr2J);
      computeError(store_dr_J, dr2J, sigma);

      //compute slope from mean and error (sigma), and store as slopeJ[i]:
      double trash; //won't use this variable. (saves sig_mu here).
      computeSlope(dr2J,sigma,slopeJ[i],trash);
    }

    //compute the standard "jackknifed" error sigmaJ2 to the slope,
    //which is the end result of the jackknife method.
    double sigmaJ2 = 0;
    double slopeJtotal = 0;

    for (int i = 0; i < noEnsembles_; i++){

      sigmaJ2 = sigmaJ2 + 1.0*(noEnsembles_ - 1)/noEnsembles_ *
        pow((slopeJ[i] - slope),2);

      slopeJtotal = slopeJtotal + slopeJ[i];

    }

    //correction for bias:
    slopeJtotal = 1.0 * slopeJtotal / noEnsembles_;
    slope = slope - 1.0*(noEnsembles_ - 1) * (slopeJtotal - slope);

    //Print result:
    //-------------
    //to screen...
    std::cout << "mu:       " << setprecision(5) << scientific << slope << std::endl
              << "sigma_mu: " << sigma_mu << std::endl
              << "sigma_J:  " << sigmaJ2 << std::endl;

    //... and file. Append "_jackknife" to the chosen output-filename
    std::string outputFile = outName + "_jackknife";

    std::ofstream printJack;
    printJack.open(outputFile.c_str());

    printJack << "#Slope (mu), original err (sigma_mu), jackknifed err (sigma_J):"
              << setprecision(5) << scientific << std::endl << slope << std::endl
              << sigma_mu << std::endl
              << sigmaJ2 << std::endl;

    printJack.close();

  }
}


void Save::computeSlope(const std::vector<double>& dr,
                        const std::vector<double>& dr_err,
                        double& mu, double& sigma_mu){
  //fit a straight line to data, either by using NR: y = a + b * x
  //or using my own implementation: y = mu * x

  //patch to remove the first values when doing the fit. Note:
  //think whether we want to use fitab or fitb when doing this.
  int minElement = 0;
  for (minElement = 0; minElement < noSamplingTimes_ &&
         samplingTime_[minElement] < minTime_; minElement++){
    //just want to iterate our "minElement" to find minimum index
  }

  bool useNR_fitting = USE_NR_FITTING;

  if (!useNR_fitting){
    //My own implementation of Ch. 15 of Numerical Recipes, where fitab.h
    //has been re-written for a y = mu*t fit. Note that Stt corresponds to
    // Sxx in NR.

    double Sty = 0, Stt = 0;

    for(int i = minElement; i < noSamplingTimes_; i++){
      Stt = Stt + 1.0*pow(samplingTime_[i],2) / pow(dr_err[i],2);
      Sty = Sty + 1.0*samplingTime_[i] * dr[i] / pow(dr_err[i],2);
    }

    //return the error in mu
    sigma_mu = sqrt(1.0 / Stt);

    //return the slope, mu (by reference)
    mu = Sty / Stt;
  }
  else{
    //If we want to use the NR code for fitting:

    //fitab wants NRvector...
    NRvector<double> sigma(noSamplingTimes_ - minElement,0.0);
    NRvector<double> X(noSamplingTimes_ - minElement,0.0);
    NRvector<double> Y(noSamplingTimes_ - minElement,0.0);

    //...so copy input STL-vectors to NRvector
    for(int i = 0; i < noSamplingTimes_ - minElement; i++){
      X[i] = samplingTime_[minElement + i];
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


inline void Save::computeMSD(const std::vector<std::vector<double> >& store,
                             std::vector<double>& msd){
  assert(msd.empty());

  //save and return (by reference&) the mean square displacement:
  msd.assign(noSamplingTimes_,0);

  //don't use the "ensemble_" variable as iteration limit, since we will be
  //fiddling with this in the jackknife implementation.
  const int ensembles = store.size();

  for (int j = 0; j < noSamplingTimes_; j++){
    for (int i = 0; i < ensembles; i++){
      msd[j] = msd[j] + 1.0 * pow(store[i][j],2) / ensembles;
    }
  }
}


//Overloaded function, to test: MSD(t_i) = sum_m [x_i^(m) - <x_i>]^2
//i.e. since M != inf  --> <x> != 0.
inline void Save::computeMSD(const std::vector<std::vector<int> >& store_dx,
                             const std::vector<std::vector<int> >& store_dy,
                             const std::vector<std::vector<int> >& store_dz,
                             const std::vector<double>& dx_mu,
                             const std::vector<double>& dy_mu,
                             const std::vector<double>& dz_mu,
                             std::vector<double>& msd){

  //save and return (by reference&) the mean square displacement:
  msd.assign(noSamplingTimes_,0);

  //don't use the "ensemble_" variable as iteration limit, since we will be
  //fiddling with this in the jackknife implementation.
  const int ensembles = store_dx.size();

  for (int j = 0; j < noSamplingTimes_; j++){
    for (int i = 0; i < ensembles; i++){
      msd[j] += 1.0/(ensembles-1.0) * (
        pow( store_dx[i][j] - dx_mu[j],2) +
        pow( store_dy[i][j] - dy_mu[j],2) +
        pow( store_dz[i][j] - dz_mu[j],2));
    }
  }
}




inline void Save::computeMean(const std::vector<std::vector<int> >& store_dx,
                              const std::vector<std::vector<int> >& store_dy,
                              const std::vector<std::vector<int> >& store_dz,
                              std::vector<double>& x_mu,
                              std::vector<double>& y_mu,
                              std::vector<double>& z_mu){

  //save and return (by reference&) the mean:
  x_mu.assign(noSamplingTimes_,0);
  y_mu.assign(noSamplingTimes_,0);
  z_mu.assign(noSamplingTimes_,0);

  for (int j = 0; j < noSamplingTimes_; j++){
    for (int i = 0; i < noEnsembles_ ; i++){
      x_mu[j] = x_mu[j] + 1.0 * store_dx[i][j] / noEnsembles_;
      y_mu[j] = y_mu[j] + 1.0 * store_dy[i][j] / noEnsembles_;
      z_mu[j] = z_mu[j] + 1.0 * store_dz[i][j] / noEnsembles_;
    }
  }
}


inline void Save::computeError(const std::vector<std::vector<double> >& store,
                               const std::vector<double>& msd,
                               std::vector<double>& sigma){
  //save and return (by reference&) the error to the msd-vector
  sigma.assign(noSamplingTimes_,0);

  //don't use the "ensemble_" variable as iteration limit, since we will
  //remove trajectories (ensembles) in the jackknife implementation.
  const int ensembles = store.size();

  for (int j = 0; j < noSamplingTimes_; j++){
    for (int i = 0; i < ensembles; i++){
      sigma[j] += pow (pow(store[i][j],2) - msd[j] ,2);
    }

    sigma[j] = sqrt( sigma[j] / (1.0 * ensembles * (ensembles-1)) );
  }
}



void Save::computeHmatrix(const std::vector<std::vector<double> >& trajectories,
                          const std::vector<double>& msd,
                          std::vector<std::vector<double> >& Hmatrix,
                          double exponent){

  if(samplingTime_[0] == 0)
    std::cout << "NOTE: If the first sampling-point t_0 = 0 then the \n"
              << "H-matrix will be singular (first row and column =0)"
              << std::endl;

  std::cout << "computing H-matrix ..." << std::endl;

  //if we use this function when bootstrapping, and for some reason don't
  //want the same number of bootstrapped synthetic simulations, as ensembles:
  //(otherwise we could use the global "noEnsembles_" variable)
  int ensembles = trajectories.size();

  //set Hmatrix to 0
  std::vector<double> temp;
  temp.assign(noSamplingTimes_,0);
  Hmatrix.assign(noSamplingTimes_,temp);

  //initiate class, to print progress to screen,
  RemainingTime computationProgress(noSamplingTimes_);

  for(int i = 0; i < noSamplingTimes_; i++){

    //print remaining progress to screen:
    computationProgress.printProgress(i);

    for(int j = 0; j < noSamplingTimes_; j++){
      for(int k = 0; k < ensembles; k++){
        Hmatrix[i][j] += 1.0 / (ensembles-1) *
          ( pow( trajectories[k][i],exponent) - msd[i] ) *
          ( pow( trajectories[k][j],exponent) - msd[j] );
      }
    }
  }

  //Note: when computing the H-matrix for the mean of mean,
  //we should divide it by (M-1), but for the "original lomholt" we
  //divide with M(M-1), (M = ensembles).

  //also, don't square trajectories again when computing mean of MSD,
  //i.e. "exponent = 1 when bootstrapping"

}


void Save::printHmatrix(const std::vector<std::vector<double> >& matrix,
                        std::string name){

  //append "_matrix" to the chosen output-filename
  std::string outputFile = name + "_matrix";
  std::cout << "printing H-matrix to file: " << outputFile << std::endl;

  //construct and open file
  std::ofstream dumpMatrix(outputFile.c_str());

  for (int i = 0; i < noSamplingTimes_; i++){
    for (int j = 0; j < noSamplingTimes_; j++)
      dumpMatrix  << matrix[i][j] << "\t";

    dumpMatrix << std::endl;
  }
  dumpMatrix.close();
}
