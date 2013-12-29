#include <iostream>        //for cout, among others
#include <fstream>         //read write to files
#include <cmath>           //gives sqrt, pow, fabs,
#include <vector>          //STL vector
#include <string>          //STL string
#include <cstring>         //needed for the strcat -function
#include <algorithm>       //for min/max element finding in histogram

#include "nr/nr3.h"
#include "nr/ludcmp.h"       //only depends on nr3.h
#include "nr/svd.h"          //singular value decomposition
#include "nr/ran.h"          //only needed for my bootstrap-function

#include "nr/gamma.h"        //needed by fitab.h (see NR p. 784)
#include "nr/incgammabeta.h" //needed by fitab.h
#include "nr/fitab.h"        //needed by computeSlope()

#include "auxiliary.h"       //for the printError()-function
#include "save.h"

//#include "nr/gaussj.h"      //needed by fitmrq.h
//#include "nr/fitmrq.h"      //not yet in use...

const bool SEED_BOOTSTRAP = 6854;    //default seed: 6854

//If true, we print each bootstrapped trajectory to separate file:
//out.dat0, out.dat1, etc, (if -w out.dat). Else (if "false") we compute
//the H-matrix directley,
//Note: if "false" the computet H will be printed to out.dat_matrix, i.e.
//same name as the first "lomholt" procedure so if this is "false" and we
//use  COMPUTE_CORRELATION = true, and USE_PEARSON_INSTEAD_OF_H_MATRIX = false
//we will get two "out.dat_matrix"-files, one replacing the other.
const bool BOOTSTRAP_CONSERVE_MEMORY = false;

//If true, also use the hybrid "bootknife" method to compute
//slope, and error, and write to separate file "*_hybrid"
const bool BOOTSTRAP_HYBRID = true;

//Prints either Pearson coefficient in 4th column, or
//the entire H-matrix in a separate file:
const bool COMPUTE_CORRELATION = false;
const bool USE_PEARSON_INSTEAD_OF_H_MATRIX = false;
//NOTE TODO: if we use this to get Hmatrix and then start bootstrapping,
//and compute a H-matrix from bootstrap they will collide, and I think
//the bootstrapped will replace the first file (out.dat_matrix)

//Bin the data, and generate a histogram in a separate
//file, to get a histogram:
const bool PRINT_DISTRIBUTION  = false;

const bool PRINT_SECOND_FILE   = false;   //print <x,y,z> (un-sqrd)

Save::Save(std::vector<double> samplingTime,
           int ensembles, bool lowMem = false){

  samplingTime_ = samplingTime;

  ensembles_ = ensembles;

  //number of elements
  maxElement_ = samplingTime_.size();

  //use the memory conserving function
  lowMem_ = lowMem;

  //initiate the vectors for the summation in store()
  //only needed if we use the memory-conserving algorithm
  if(lowMem_){
    dr4_err_.assign(maxElement_,0);
    dx4_err_.assign(maxElement_,0);
    dy4_err_.assign(maxElement_,0);
    dz4_err_.assign(maxElement_,0);

    r2_mu_.assign(maxElement_,0);
    x2_mu_.assign(maxElement_,0);
    y2_mu_.assign(maxElement_,0);
    z2_mu_.assign(maxElement_,0);
  }
  else{
    //Allocate size for the matrix to store our trajectories

    /*
    //save position for each ensemble and store them all
    std::vector<int> temp;
    temp.reserve()
    std::vector< std::vector<int> > store_dy_;
    std::vector< std::vector<int> > store_dz_;
    std::vector< std::vector<double> > store_dr_;

    std::vector<double> test(10); //allocate for 10 double
    test.reserve(10);             //samma

    std::vector< std::vector<double> > store_dr_;


    */
  }

}


void Save::computeCorrelation(std::string inputName){
  //This function calculates the Pearson coefficients (see Numerical
  //recipes ch 14), to get the correlation between two sampling points
  //since if the tagged particle jumps once/sec sampling the system
  //10/s will not yeld any new information, and it will mess up the
  //error estimation that assumes the sampling points are uncorrelated
  //and with a Gaussian distribution.

  //perhaps check: Physical Review D vol 34, 2469 (1986), sec. D

  if(lowMem_)
    printError("Can not compute correlation, no available data.");

  if(r2_mu_.empty())
    printError("Must compute MSD first! (ie. r2_mu_)");

  //Either compute the Pearson correlation coefficient,
  //or use the full matrix
  if(USE_PEARSON_INSTEAD_OF_H_MATRIX){

    std::vector<double> pearsonCorrelation;
    pearsonCorrelation.assign(maxElement_,0);

    //compute the correlation between two points, i, j, for i = j+1
    std::vector<double> h;
    h.assign(maxElement_,0);

    for(int i = 1; i < maxElement_; i++){
      for(int j = 0; j < ensembles_; j++){
        h[i] += 1.0 / (ensembles_ - 1) *
          (pow(store_dr_[j][i],2) - r2_mu_[i])*
          (pow(store_dr_[j][i-1],2) - r2_mu_[i-1]);
      }
    }

    // pearsonCorrelation[i] = h_ij / (sigma_i*sigma_j)
    // sigma_r = r2_err_ * sqrt(ensembles);
    for(int i = 1; i < maxElement_; i++){
      if(r2_err_[i-1] != 0)
        pearsonCorrelation[i] = (double) h[i] /
          (r2_err_[i] * r2_err_[i-1] * ensembles_);

      //NOTE: error in first value is 0, so the correlation between value 0 & 1
      //will be inf (nan)! thus the if-statement

      //Now save the inf. for printing to the output-file.
      fourthColumn_ = pearsonCorrelation;
    }
  }
  else{
    //compute entire matrix, to check the correlation between every value
    std::vector<std::vector<double> > H;

    //note! does not devid with the factor (1/ensembles) so I do that below
    computeHmatrix(store_dr_, r2_mu_, H);

    //In order to use NR implementations to test the matrix, and perform
    //LU-decomposition, we transform STL-matrix 'H' to NR-matrix
    NRmatrix<double> h(maxElement_, maxElement_, 0.0);

    for(int i = 0; i < maxElement_; i++){
      for(int j = 0; j < maxElement_; j++){

        //computeHmatrix only devides with (ensembles-1):
        H[i][j] = H[i][j] / (double) ensembles_;

        //now copy to NR matrix.
        h[i][j] = H[i][j];
      }
    }

    //Print H-matrix to file (now that we devided it with M(M-1), which is
    //done in the for-loop above (M), and in the function "computeHmatrix()"
    //where we devide it with (M-1), where M = ensembles_
    printHmatrix(H,inputName);

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
    if(false){//TEST:  //just a small test of NR conditional number.
      NRmatrix<double> hm(2,2,0.0);
      hm[0][0] = 1.0e-16;    //epsilon
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
    NRvector<double> b(maxElement_),z(maxElement_);

    for(unsigned int i = 0; i < maxElement_; i++)
      b[i] = samplingTime_[i];

    LUdcmp LUdecomposition(h);
    LUdecomposition.solve(b,z);

    //copy z to the class-private vector, that we will print as the
    //fourth column in the output-file. Transfer a NR vector to std::
    fourthColumn_.assign(maxElement_,0);
    for (int i = 0; i < maxElement_; i++)
      fourthColumn_[i] = z[i];
    std::cout << "fourth column in " << inputName << " is z" << std::endl;

  }//end-else
}



void Save::store(const std::vector<int>& dx, const std::vector<int>& dy,
                 const std::vector<int>& dz, const std::vector<double>& dr){
  //This is needed to calculate error bars in computeStdErr().

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

  //Check that we assigned the vectors =0 if we uses the lowMem-option.
  //For lowMem=false, we do this in the Save & computeStdErr()-functions
  //since we use push_back() instead.

  if(lowMem_){
    bool error = false;
    if (dr4_err_.size() == 0) error = true;
    if (dx4_err_.size() == 0) error = true;
    if (dy4_err_.size() == 0) error = true;
    if (dz4_err_.size() == 0) error = true;

    if (r2_mu_.size() == 0) error = true;
    if (x2_mu_.size() == 0) error = true;
    if (y2_mu_.size() == 0) error = true;
    if (z2_mu_.size() == 0) error = true;

    if(error){
      printError("Incorrect assignment of vectors for stdErr-LowMem");
    }
  }

  if(!lowMem_){
    //if don't use the lowMem-algorithm:
    //this gives dx[ensembler][number of values]
    store_dx_.push_back(dx);
    store_dy_.push_back(dy);
    store_dz_.push_back(dz);
    store_dr_.push_back(dr);

  }
  else{
    //To conserve memory we do the summation directly, which might
    //cause a numerical error in the end compared to doing it the
    //more rigorous way, but all tests I've made returns identical
    //output
    for(int i = 0; i < maxElement_; i++){
      dr4_err_[i] = dr4_err_[i] + pow( dr[i] ,4);
      dx4_err_[i] = dx4_err_[i] + pow( dx[i] ,4);
      dy4_err_[i] = dy4_err_[i] + pow( dy[i] ,4);
      dz4_err_[i] = dz4_err_[i] + pow( dz[i] ,4);

      r2_mu_[i] = r2_mu_[i] + pow( dr[i] ,2);
      x2_mu_[i] = x2_mu_[i] + pow( dx[i] ,2);
      y2_mu_[i] = y2_mu_[i] + pow( dy[i] ,2);
      z2_mu_[i] = z2_mu_[i] + pow( dz[i] ,2);
    }
  }
}



void Save::save(std::string name, std::string head){
  //calculates the MSD from the stored positions of the tracer,
  //and prints to file.

  //sore averaged displacement <x>, <y>, <z>, <r>
  //mostly used to check that is is ~ 0.
  std::vector<double> x_mu, y_mu, z_mu, r_mu;

  ofstream secondary;
  bool printSecondFile = PRINT_SECOND_FILE;

  if(!lowMem_){

    //hard coded switch to print x,y,z, x^2,y^2,z^2
    //specifically. Not needed now, since system is isotropic,

    if(printSecondFile){
      //Check that <x,y,z>=0, just for testing purpose.
      //Compute <R>:

      std::string secondFileName = name + "_txyz";
      secondary.open(secondFileName.c_str());

      x_mu.assign(maxElement_,0);
      y_mu.assign(maxElement_,0);
      z_mu.assign(maxElement_,0);
      r_mu.assign(maxElement_,0);

      for (int j = 0; j < maxElement_; j++){
        for (int i = 0; i < ensembles_ ; i++){
          x_mu[j] = x_mu[j] + 1.0 * store_dx_[i][j] / ensembles_;
          y_mu[j] = y_mu[j] + 1.0 * store_dy_[i][j] / ensembles_;
          z_mu[j] = z_mu[j] + 1.0 * store_dz_[i][j] / ensembles_;
          //r_mu_[j] = r_mu_[j] + 1.0*store_dr_[i][j] / ensembles_; //not using this
        }
      }
    }//end of printSecondFile if statement

    //Compute <R^2>:
    x2_mu_.assign(maxElement_,0);
    y2_mu_.assign(maxElement_,0);
    z2_mu_.assign(maxElement_,0);
    r2_mu_.assign(maxElement_,0);

    for (int j = 0; j < maxElement_; j++){
      for (int i = 0; i < ensembles_;i++){
        x2_mu_[j] = x2_mu_[j] + 1.0 * pow(store_dx_[i][j],2) / ensembles_;
        y2_mu_[j] = y2_mu_[j] + 1.0 * pow(store_dy_[i][j],2) / ensembles_;
        z2_mu_[j] = z2_mu_[j] + 1.0 * pow(store_dz_[i][j],2) / ensembles_;
        r2_mu_[j] = r2_mu_[j] + 1.0 * pow(store_dr_[i][j],2) / ensembles_;
      }
    }

    //Calculate standard error, to generate error-bars
    //this function must come after x2_mu_...r2_mu_ has been
    //computed.
    computeStdErr();

    if(COMPUTE_CORRELATION)
      //This can be very time consuming for big matrices
      //(many ensembles and sampling points)
      computeCorrelation(name);
    else
      fourthColumn_.assign(maxElement_,0);

  }
  else{
    //use the memory conserving version instead. This one does not
    //use the dx[][]-matrix, and in that way we save lots of mem.
    computeStdErrLowMem();
  }


  if(lowMem_)  //NOTE: temporary, if lowMem is on!
    //TODO: calculate this even if lowMem_ = true!
    //does not serve any purpose at this point.
    fourthColumn_.assign(maxElement_,0);

  //print data to file:
  ofstream primary;
  primary.open(name.c_str());
  primary << head;
  std::cout << std::endl << head;

  for (int i = 0; i < maxElement_; i++){

    //Print it to file:
    primary << samplingTime_[i] << "\t" << r2_mu_[i] << "\t"
            << r2_err_[i] << "\t" << fourthColumn_[i] << std::endl;
    if(printSecondFile){
      secondary << samplingTime_[i] << "\t"
                << x_mu[i] << "\t" << y_mu[i] << "\t" << z_mu[i] <<"\t"
                << x2_mu_[i] <<"\t"<< y2_mu_[i] <<"\t"<< z2_mu_[i] <<"\t"
                << x2_err_[i]<<"\t"<<y2_err_[i]<<"\t"<<z2_err_[i]<<std::endl;
    }
  }
  primary.close();
  if(printSecondFile) secondary.close();

  if(PRINT_DISTRIBUTION)
    //generate a file with histogram, (just for testing)
    saveBinning(name);
}


void Save::computeStdErr(void){
  //Calculate the standard deviation, to get error-bars.

  //Should not be using this function if lowMem = true!
  if(lowMem_)
    printError("LowMem switch ON, but using wrong stdErr-function!");

  x2_err_.assign(maxElement_,0);
  y2_err_.assign(maxElement_,0);
  z2_err_.assign(maxElement_,0);
  r2_err_.assign(maxElement_,0);

  for (int j = 0; j < maxElement_; j++){
    for (int i = 0; i < ensembles_; i++){
      x2_err_[j] += pow (pow(store_dx_[i][j],2) - x2_mu_[j] ,2);
      y2_err_[j] += pow (pow(store_dy_[i][j],2) - y2_mu_[j] ,2);
      z2_err_[j] += pow (pow(store_dz_[i][j],2) - z2_mu_[j] ,2);
      r2_err_[j] += pow (pow(store_dr_[i][j],2) - r2_mu_[j] ,2);
    }//( "a += b" is the same as "a = a + b")
  }

  for (int j = 0; j < maxElement_  ; j++){
    x2_err_[j] = sqrt( x2_err_[j] / (1.0*ensembles_*(ensembles_ - 1)) );
    y2_err_[j] = sqrt( y2_err_[j] / (1.0*ensembles_*(ensembles_ - 1)) );
    z2_err_[j] = sqrt( z2_err_[j] / (1.0*ensembles_*(ensembles_ - 1)) );
    r2_err_[j] = sqrt( r2_err_[j] / (1.0*ensembles_*(ensembles_ - 1)) );
  } //for E = 1 --> inf !


  if(false){
    //Testing the Gaussian distribution of the sampling points:
    float SUM68 = 0;
    for (int j = 0; j < maxElement_; j++){    //forgot about his:
      double top = x2_mu_[j] +0.5*x2_err_[j] * sqrt(ensembles_-1)*2;
      double bot = x2_mu_[j] -0.5*x2_err_[j] * sqrt(ensembles_-1)*2;
      for (int i=0; i < ensembles_; i++){
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
              << 1.0 * SUM68/(ensembles_ * maxElement_) << std::endl;
  }

}


void Save::computeStdErrLowMem(void){

  if(!lowMem_)
    printError("lowMem = OFF, but using StdErrLowMem()-func!");

  //Divide the summed up displacement to get the average:
  for (int i = 0; i< maxElement_; i++){
    x2_mu_[i] = x2_mu_[i] / (ensembles_ * 1.0);
    y2_mu_[i] = y2_mu_[i] / (ensembles_ * 1.0);
    z2_mu_[i] = z2_mu_[i] / (ensembles_ * 1.0);
    r2_mu_[i] = r2_mu_[i] / (ensembles_ * 1.0);
  }

  x2_err_.assign(maxElement_,0);
  y2_err_.assign(maxElement_,0);
  z2_err_.assign(maxElement_,0);
  r2_err_.assign(maxElement_,0);

  //Following the notation in the comments in Lattice::store():
  // [1/N sum_i^N a_i^2] - [<a>^2], (a=dr^2)
  for(int i = 0; i < maxElement_; i++){
    x2_err_[i] = dx4_err_[i] / (1.0 * ensembles_) - pow( x2_mu_[i] ,2);
    y2_err_[i] = dy4_err_[i] / (1.0 * ensembles_) - pow( y2_mu_[i] ,2);
    z2_err_[i] = dz4_err_[i] / (1.0 * ensembles_) - pow( z2_mu_[i] ,2);
    r2_err_[i] = dr4_err_[i] / (1.0 * ensembles_) - pow( r2_mu_[i] ,2);
  }

  for(int i =0; i < maxElement_; i++){
    x2_err_[i] = sqrt(x2_err_[i] / (1.0 * (ensembles_ - 1) ));
    y2_err_[i] = sqrt(y2_err_[i] / (1.0 * (ensembles_ - 1) ));
    z2_err_[i] = sqrt(z2_err_[i] / (1.0 * (ensembles_ - 1) ));
    r2_err_[i] = sqrt(r2_err_[i] / (1.0 * (ensembles_ - 1) ));
  }
}


void Save::saveBinning(std::string fileName){
  //This function is for getting P(t,x), ie. the probability
  //distribution of finding a particle at r^2, at time t.

  if(!lowMem_){
    const int noBins = 300;    //number of bins
    std::vector<std::vector<float> > printToFileR2;
    std::vector<std::vector<float> > printToFileX2;
    std::vector<std::vector<float> > printToFileY2;

    //find the smallest/biggest value in store_dr_[*][t]
    //-----------------------------------------------
    int samplingTime = maxElement_ -1; //pick sampling time

    std::vector<double> tempVectorR2;
    std::vector<double> tempVectorX2;
    std::vector<double> tempVectorY2;

    tempVectorR2.assign(ensembles_,0);
    tempVectorX2.assign(ensembles_,0);
    tempVectorY2.assign(ensembles_,0);

    for (int i = 0; i < ensembles_; i++){
      //(I don't know of a better way to do this)
      tempVectorR2[i] = pow(store_dr_[i][samplingTime],2);
      tempVectorX2[i] = pow(store_dx_[i][samplingTime],2);
      tempVectorY2[i] = pow(store_dy_[i][samplingTime],2);
    }

    double minR2 = *min_element(tempVectorR2.begin(), tempVectorR2.end());
    double maxR2 = *max_element(tempVectorR2.begin(), tempVectorR2.end());
    double minX2 = *min_element(tempVectorX2.begin(), tempVectorX2.end());
    double maxX2 = *max_element(tempVectorX2.begin(), tempVectorX2.end());
    double minY2 = *min_element(tempVectorY2.begin(), tempVectorY2.end());
    double maxY2 = *max_element(tempVectorY2.begin(), tempVectorY2.end());

    double binStepR2 = (double) (maxR2 - minR2)/noBins;
    double binStepX2 = (double) (maxX2 - minX2)/noBins;
    double binStepY2 = (double) (maxY2 - minY2)/noBins;

    // if (binStepX < 0) binStepX = -1 * binStepX;
    // if (binStepY < 0) binStepY = -1 * binStepY;

    std::cout << "binStep R2: " << binStepR2 << std::endl;
    std::cout << "binStep X2: " << binStepX2 << std::endl;
    std::cout << "binStep Y2: " << binStepY2 << std::endl;

    std::vector<float> probAxisR2;    //Probability count
    std::vector<float> binAxisR2;     //binning-axis
    std::vector<float> probAxisX2;    //Probability count
    std::vector<float> binAxisX2;     //binning-axis
    std::vector<float> probAxisY2;    //Probability count
    std::vector<float> binAxisY2;     //binning-axis

    binAxisR2.assign(noBins,0);
    probAxisR2.assign(noBins,0);
    binAxisX2.assign(noBins,0);
    probAxisX2.assign(noBins,0);
    binAxisY2.assign(noBins,0);
    probAxisY2.assign(noBins,0);

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

      std::cout << "binAxisR2["<<i<<"] = " << binAxisR2[i] << std::endl;
    }

    const double PI = 3.14159265;
    double r2;                      //displacement
    double x2;                      //displacement
    double y2;                      //displacement

    std::cout << "Calculating P(t,x) for histogram... ";

    //Analyze all sampling-times, or just one? (start with one)
    const int startTimeElement = samplingTime; // 0;
    const int endTime = startTimeElement + 1;  // maxElement_;

    for(int t = startTimeElement; t < endTime ; t++){ //loop over time
      for (int j = 0; j < noBins; j++){               //loop over bin-cells
        for(int e = 0; e < ensembles_; e++){          //loop over "ensembles"
          r2 = pow(store_dr_[e][t],2);
          x2 = pow(store_dx_[e][t],2);
          y2 = pow(store_dy_[e][t],2);

          //why 2*PI ? r? (November 2010)
          //if(binAxis[j+1] > r2 && r2 >= binAxis[j])
          //probAxis[j] = probAxis[j] + 1.0/(ensembles_ * 2 * PI * ((r2 + binStep)/2));

          if(binAxisR2[j+1] > r2 && r2 >= binAxisR2[j])
            probAxisR2[j] = probAxisR2[j] + 1.0 / ensembles_;

          if(binAxisX2[j+1] > x2 && x2 >= binAxisX2[j])
            probAxisX2[j] = probAxisX2[j] + 1.0 / ensembles_;

          if(binAxisY2[j+1] > y2 && y2 >= binAxisY2[j])
            probAxisY2[j] = probAxisY2[j] + 1.0 / ensembles_;

        }
      }
      printToFileR2.push_back(probAxisR2);
      printToFileX2.push_back(probAxisX2);
      printToFileY2.push_back(probAxisY2);

      probAxisR2.clear();
      probAxisX2.clear();
      probAxisY2.clear();

      probAxisR2.assign(noBins,0);
      probAxisX2.assign(noBins,0);
      probAxisY2.assign(noBins,0);

    }

    //Print to file:
    //-------------
    ofstream histR2;
    ofstream histX2;
    ofstream histY2;

    std::string histogramNameR2 = fileName + "_histogram";
    std::string histogramNameX2 = fileName + "_histogramX2";
    std::string histogramNameY2 = fileName + "_histogramY2";

    histR2.open(histogramNameR2.c_str());
    histX2.open(histogramNameX2.c_str());
    histY2.open(histogramNameY2.c_str());

    std::cout << " ...Printing histogram... " << std::endl;
    int maxElement = printToFileR2.size();   //or maxElement_
    for(int i = 0; i < maxElement; i++){   //loop over time
      for(int j = 0; j < noBins; j++){     //loop over bins
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
    for (int i = 0; i < noBins; i++){
      MSD = MSD + binAxisR2[i] * printToFileR2[0][i];
    }
    std::cout << "New MSD value: " << MSD << std::endl;

  }
  else
    std::cout << "lowMem = ON, therefore no information stored!"
              << std::endl << "can't print cluster-information"
              << std::endl;
}


void Save::computeBootstrap(std::string name, int numberOfRuns,
                            std::string head = ""){
  //Generate synthetic values, from current trajectories in store_dr_
  //to get new values, by Quick-and-dirty Monte Carlo Bootstrapping.
  //(TODO: actually I should re-use the previous code for this, since there
  //is an overlap, it would look better.)

  //Needs all the information for all trajectories saved in store_dr_,
  //which will not be saved in lowMem mode.
  if(lowMem_){
    std::cout << "Error: Can not use bootstrap when in low memory mode"
              << std::endl;
  }
  else{

    std::cout << "bootstrapping!" << std::endl;

    std::vector<double> r2_mu, temp;
    std::vector<std::vector<double> > store_dr;

    //Assign store_dr (local variable) with zeros
    temp.assign(maxElement_,0);
    store_dr.assign(ensembles_,temp);

    Ran randomNumber(SEED_BOOTSTRAP);

    //either print "numberOfRuns" number of synthetic output-files to
    //hard drive, each file corresponds to the MSD of a synthetic simulation,
    //or save all MSDs in a single file. (easier for matlab,
    //but uses more memory
    bool conserveMemory = BOOTSTRAP_CONSERVE_MEMORY;
    std::vector<std::vector<double> > catanateOutputFiles;
    if(!conserveMemory)  catanateOutputFiles.assign(numberOfRuns,temp);


    //used in hybrid-loop, this does not look good!
    std::vector<double> mu, sigma_mu, dr2_err;
    bool hybrid = BOOTSTRAP_HYBRID;

    //store all MSD-fitted slopes (and corresponding error)
    mu.assign(numberOfRuns,0);
    sigma_mu.assign(numberOfRuns,0);


    //Loop to generate "numberOfRuns" number of synthetic MSDs of
    //trajectories, each saved MSD corresponds to a simulation of
    //"ensembles_" number of individual trajectories (i.e. ensembles).
    for(int k = 0; k < numberOfRuns; k++){

      //Randomly pick real trajectories into our virtual simulation
      //-----------------------------------------------------------
      int r;
      for (int i = 0; i < ensembles_; i++){
        r = (int) (ensembles_ * randomNumber.doub());
        store_dr[i] = store_dr_[r];
      }

      //Compute <R^2>:
      computeMSD(store_dr,r2_mu);


      if(hybrid){
        //this is an alternative "hybrid" version of bootstrap, where we
        //refrain from using H-matrices, and just fit a line to current
        //bootstrapped "simulation" (as if uncorrelated, like in the
        //jackknife-method), and then take all my slopes and compute an
        //average, (which must be done outside this for(k)-loop).


        //return error in each samplingpoint by dr2_err, needed
        //when fitting mu. New for each bootstrapped "simulation"
        computeError(store_dr,r2_mu,dr2_err);

        //get mu and sigma_mu
        computeSlope(r2_mu, dr2_err, mu[k], sigma_mu[k]);

        //TEST
        //std::cout << setprecision(5) << scientific << mu[k] << std::endl;

      }


      if (conserveMemory){

        //Print data to file (for processing in MatLab):
        //---------------------------------------------
        ofstream bootOutput;
        stringstream tempName;
        //append run-iterator (k) to filename:
        tempName << name << k;
        std::string nameOfFile = tempName.str();
        bootOutput.open(nameOfFile.c_str());

        bootOutput << head;
        for (int i = 0; i < maxElement_; i++){
          bootOutput  << samplingTime_[i] << "\t" << r2_mu[i] << std::endl;
        }
        bootOutput.close();
      }
      else{
        //store all bootstraped data in a single matrix...
        catanateOutputFiles[k] = r2_mu;
        //...and print it later, to a single file instead.
      }
    } //(end of for-loop (with k-iterator))


    if(hybrid){
      //continue on task done in previous if(hybrid)-statement

      double mu_average = 0;
      double sigma_mu_final = 0;

      for (int i = 0; i < numberOfRuns; i++){
        mu_average += mu[i];
      }

      mu_average = mu_average /( (double) numberOfRuns);

      for (int i = 0; i < numberOfRuns; i++){
        sigma_mu_final += pow(mu[i] - mu_average ,2);
      }

      sigma_mu_final = sqrt( sigma_mu_final /((double) (numberOfRuns-1) ));

      std::cout <<  "mu: " << mu_average << std::endl
                <<  "sig_mu: " << sigma_mu_final << std::endl;


      //Print bootstrapped mean of mean  to file:
      //-----------------------------------------
      ofstream hybridOut;
      std::string nameOfFile = name + "_hybrid";
      hybridOut.open(nameOfFile.c_str());

      hybridOut  <<  " # mu, and sig_mu " << std::endl
                 << mu_average << std::endl
                 << sigma_mu_final << std::endl;
      hybridOut.close();

    }


    if(!conserveMemory){
      //Since we did not dump each MSD-trajectory to file, we
      //compute the H-matrix and dump that to a single file.

      //Compute the mean of the Mean Sq. Disp. (don't square it!)
      //use numberOfRuns instead of ensembles_ since might not be the same
      std::vector<double> meanOfMsd;
      meanOfMsd.assign(maxElement_,0);
      for (int j = 0; j < maxElement_; j++){
        for (int i = 0; i < numberOfRuns; i++){
          meanOfMsd[j] = meanOfMsd[j] + 1.0 * catanateOutputFiles[i][j]
            / numberOfRuns;
        }
      }

      //Save H-matrix to this:
      std::vector<std::vector<double> > H;
      computeHmatrix(catanateOutputFiles,meanOfMsd ,H, 1);

      //print H-matrix to file, with "_matrix" appended:
      printHmatrix(H,name);

      //Print bootstrapped mean of mean  to file:
      //-----------------------------------------
      ofstream msdBootout;
      std::string nameOfFile = name + "_bootstrap";
      msdBootout.open(nameOfFile.c_str());

      msdBootout << head;
      for (int i = 0; i < maxElement_; i++){
       msdBootout  << samplingTime_[i] << "\t" << meanOfMsd[i] << std::endl;
      }
      msdBootout.close();

    }

  }
}


void Save::computeJackknife(std::string outName){
  //http://www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife/
  //use the jackknife method to compute the error from correlated data,
  //by sucsessivley (TODO) removing trajectories (i.e. "ensembles").

  //only works if we have all trajectories. (store_dr_-matrix)
  if(lowMem_){
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
    slopeJ.assign(ensembles_,0);

    //matrix of all simulation trajectories that we will "Jackknife"
    std::vector<std::vector<double> > store_dr_J;
    std::vector<double> dr2J, sigma;
    for (int i = 0; i < ensembles_; i++){
      store_dr_J = store_dr_;
      store_dr_J.erase(store_dr_J.begin()+i);

      //compute mean for each time based on the ensembles_-1
      //trajectories and return the MSD as dr2J
      computeMSD(store_dr_J, dr2J);
      computeError(store_dr_J, dr2J, sigma);
      //compute slope from mean and error (sigma), and store as slopeJ[i]:
      computeSlope(dr2J,sigma,slopeJ[i]);
    }

    //compute the standard "jackknifed" error sigmaJ2 to the slope,
    //which is the end result of the jackknife method.
    double sigmaJ2 = 0;
    double slopeJtotal = 0;

    for (int i = 0; i < ensembles_; i++){

      sigmaJ2 = sigmaJ2 + 1.0*(ensembles_ - 1)/ensembles_ *
        pow((slopeJ[i] - slope),2);

      slopeJtotal = slopeJtotal + slopeJ[i];

    }

    //correction for bias:
    slopeJtotal = 1.0 * slopeJtotal / ensembles_;
    slope = slope - 1.0*(ensembles_ - 1) * (slopeJtotal - slope);

    //Print result:
    //-------------
    //to screen...
    std::cout << "mu:       " << setprecision(5) << scientific << slope << std::endl
              << "sigma_mu: " << sigma_mu << std::endl
              << "sigma_J:  " << sigmaJ2 << std::endl;

    //... and file. Append "_jackknife" to the chosen output-filename
    std::string outputFile = outName + "_jackknife";

    ofstream printJack;
    printJack.open(outputFile.c_str());

    printJack << "#Slope (mu), original err (sigma_mu), jackknifed err (sigma_J):"
              << setprecision(5) << scientific << std::endl << slope << std::endl
              << sigma_mu << std::endl
              << sigmaJ2 << std::endl;

    printJack.close();

  }
}


//Next two functions are overloaded versions of the same thing.
//(we are only interested in the error sigma_mu for the final fitting,
//not for each M-1 fittings to get muJ.

inline void Save::computeSlope(const std::vector<double>& dr,
                               const std::vector<double>& dr_err,
                               double& mu, double& sigma_mu){

  computeSlopeInner(dr, dr_err, mu, sigma_mu);

}

inline void Save::computeSlope(const std::vector<double>& dr,
                               const std::vector<double>& dr_err,
                               double& mu){

  //if no sigma_mu provided, make one up, (it will not be used)
  double sigma_mu;

  computeSlopeInner(dr, dr_err, mu, sigma_mu);
}


void Save::computeSlopeInner(const std::vector<double>& dr,
                             const std::vector<double>& dr_err,
                             double& mu, double& sigma_mu){
  //fit a straight line to data, either by using NR: y = a + b * x
  //or using my own implementation: y = mu * x
  //(pass by reference to increase speed)

  bool useNR_fitting = false;

  if (!useNR_fitting){
    //My own implementation of Ch. 15 of Numerical Recipes, where fitab.h
    //has been re-written for a y = mu*t fit. Note that Stt corresponds to
    // Sxx in NR.

    double Sty = 0;
    double Stt = 0;

    for(int i = 0; i < maxElement_; i++){
      Stt = Stt + 1.0*pow(samplingTime_[i],2) / pow(dr_err[i],2);
      Sty = Sty + 1.0*samplingTime_[i] * dr[i] / pow(dr_err[i],2);
    }

    //return the error in mu
    sigma_mu = sqrt(1.0 / Stt);

    //return the slope (mu) by reference (see argument list)
    mu = Sty / Stt;

  }
  else{
    //If we want to use the NR code for fitting:

    //fitab wants NRvector...
    NRvector<double> sigma(maxElement_,0.0);
    NRvector<double> X(maxElement_,0.0);
    NRvector<double> Y(maxElement_,0.0);

    //...so copy input STL-vectors to NRvector
    for(int i = 0; i < maxElement_; i++){
      X[i] = samplingTime_[i];
      Y[i] = dr[i];
      sigma[i] = dr_err[i];
    }

    //see NR p 784, fit the function y = a + bx
    Fitab fitting(X,Y,sigma);

    //Output to the equation y = a + bx
    // std::cout << setprecision(10) << scientific
    //           <<"a "      << fitting.a      << std::endl
    //           <<"b "      << fitting.b      << std::endl
    //           <<"sigmA "  << fitting.siga   << std::endl
    //           <<"sigmB "  << fitting.sigb   << std::endl
    //           <<"chi2 "   << fitting.chi2   << std::endl
    //           <<"q "      << fitting.q      << std::endl
    //           <<"sigdat " << fitting.sigdat << std::endl
    //           << std::endl;

    //return error in mu (called 'b' in NR)
    sigma_mu = fitting.sigb;

    //return the slope by reference (see argument list)
    mu = fitting.b;
  }

}


inline void Save::computeMSD(const std::vector<std::vector<double> >& store,
                              std::vector<double>& msd){
  //save and return (by reference&) the mean square displacement:
  msd.assign(maxElement_,0);

  //don't use the "ensemble_" variable as iteration limit, since we will be
  //fiddling with this in the jackknife implementation.
  const int ensembles = store.size();

  for (int j = 0; j < maxElement_; j++){
    for (int i = 0; i < ensembles; i++){
      msd[j] = msd[j] + 1.0 * pow(store[i][j],2) / ensembles;
    }
  }
}


inline void Save::computeMean(const std::vector<std::vector<double> >& store,
                              std::vector<double>& mean){
  //same function as computeMSD, but don't square it.

  //save and return (by reference&) the mean:
  mean.assign(maxElement_,0);

  //don't use the "ensemble_" variable as iteration limit, since we will be
  //fiddling with this in the jackknife implementation.
  const int ensembles = store.size();

  for (int j = 0; j < maxElement_; j++){
    for (int i = 0; i < ensembles; i++){
      mean[j] = mean[j] + 1.0 * store[i][j] / ensembles;
    }
  }
}

inline void Save::computeError(const std::vector<std::vector<double> >& store,
                               const std::vector<double>& msd,
                               std::vector<double>& sigma){
  //save and return (by reference&) the error to the msd-vector
  sigma.assign(maxElement_,0);

  //don't use the "ensemble_" variable as iteration limit, since we will
  //remove trajectories (ensembles) in the jackknife implementation.
  const int ensembles = store.size();

  for (int j = 0; j < maxElement_; j++){
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

  std::cout << "\n computing H-matrix ..." << std::endl;

  //if we use this function when bootstrapping, and for some reason don't
  //want the same number of bootstrapped sythetic simulations, as ensembles:
  //(otherwivse we could use the global "ensembles_" variable)
  int ensembles = trajectories.size();

  //set Hmatrix to 0
  std::vector<double> temp;
  temp.assign(maxElement_,0);
  Hmatrix.assign(maxElement_,temp);

  for(int i = 0; i < maxElement_; i++){
    for(int j = 0; j < maxElement_; j++){
      for(int k = 0; k < ensembles; k++){

        Hmatrix[i][j] += 1.0 / (ensembles-1) *
          ( pow( trajectories[k][i],exponent) - msd[i] ) *
          ( pow( trajectories[k][j],exponent) - msd[j] );

      }
    }
  }

  //Note: when computing the H-matrix for the mean of mean,
  //we should devide it by (M-1), but for the "original lomholt" we
  //devide with M(M-1), (M = ensembles).

  //also, don't square trajectories again when computing mean of MSD,
  //i.e. "exponent = 1 when bootstrapping"

}


void Save::printHmatrix(const std::vector<std::vector<double> >& matrix,
                        std::string name){

  //append "_matrix" to the chosen output-filename
  std::string outputFile = name + "_matrix";
  std::cout << "printing H-matrix to file: " << outputFile << std::endl;

  //open file
  ofstream dumpMatrix;
  dumpMatrix.open(outputFile.c_str());

  for (int i = 0; i < maxElement_; i++){
    for (int j = 0; j < maxElement_; j++)
      dumpMatrix  << matrix[i][j] << "\t";

    dumpMatrix << std::endl;
  }
  dumpMatrix.close();
}
