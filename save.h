#ifndef SAVE_H
#define SAVE_H

//TODO: works without these, but would it be correct to exclude them?
#include <vector>
#include <string>

class Save{
private:

  //average of displacement (MSD) squared: <X^2>, <Y^2>,...
  std::vector<double> x2_mu_, y2_mu_, z2_mu_, r2_mu_;

  //Store std error, in each point
  std::vector<double> x2_err_, y2_err_ ,z2_err_ ,r2_err_;

  //needed if we use stdErrLowMem()
  std::vector<double> dr4_err_, dx4_err_, dy4_err_, dz4_err_;

  //Switch between the two stdErr-functions above.
  bool lowMem_;

  //needed only for bootknifing over many different staring times
  //defined through the dimensionless quantity k_t * t.
  float k_t_;

  //print a fourth column of inf. to the out-put file
  //could be the pearson correlation coefficient or "z"
  //depending on switch.
  std::vector<double> fourthColumn_;

  //save position for each ensemble and store them all
  std::vector< std::vector<int> > store_dx_;
  std::vector< std::vector<int> > store_dy_;
  std::vector< std::vector<int> > store_dz_;
  std::vector< std::vector<double> > store_dr_;

  int noEnsembles_;

  //Number of elements in time vector (# of sampling points)
  int noSamplingTimes_;

  //If we want to change starting time, when performing the fitting
  //in the bootstrap --> "bootknife" version. Defaults to 0.
  float minTime_;

  //times to compute the MSD. (x-axis)
  std::vector<double> samplingTime_;

  //standard error
  void computeStdErr(void);
  void computeStdErrLowMem(void);

  //compute pearson or H-matrix
  void computeCorrelation(std::string name);

  void computeHmatrix(const std::vector<std::vector<double> >& trajectories,
                      const std::vector<double>& msd,
                      std::vector<std::vector<double> >& Hmatrix, double exp=2);

  //print H-matrix to file, append "_matrix" to file name
  void printHmatrix(const std::vector<std::vector<double> >& matrix,
                    std::string name);


  //function needed for computeJackknife():
  inline void computeSlope(const std::vector<double>& msd,
                           const std::vector<double>& sigma,
                           double& mu);

  inline void computeSlope(const std::vector<double>& msd,
                           const std::vector<double>& sigma,
                           double& mu, double& sigma_mu);

  void computeSlopeInner(const std::vector<double>& dr,
                         const std::vector<double>& dr_err,
                         double& mu, double& sigma_mu);

  //used by bootknife (same as computeMSD, but don't square)
  inline void computeMean(const std::vector<std::vector<int> >& store_dx,
                          const std::vector<std::vector<int> >& store_dy,
                          const std::vector<std::vector<int> >& store_dz,
                          std::vector<double>& x_mu,std::vector<double>& y_mu,
                          std::vector<double>& z_mu);

  //Overloaded function, to test: MSD(t_i) = sum_m [x_i^(m) - <x_i>]^2
  //i.e. since M != inf  --> <x> != 0.
  inline void computeMSD(const std::vector<std::vector<int> >& store_dx,
                         const std::vector<std::vector<int> >& store_dy,
                         const std::vector<std::vector<int> >& store_dz,
                         const std::vector<double>& dx_mu,
                         const std::vector<double>& dy_mu,
                         const std::vector<double>& dz_mu,
                         std::vector<double>& msd);

  inline void computeMSD(const std::vector<std::vector<double> >& store,
                          std::vector<double>& msd);

  inline void computeError(const std::vector<std::vector<double> >& store,
                           const std::vector<double>& msd,
                           std::vector<double>& sigma);

public:
  //Constructor. (lowMem defaults to false, (bool lowMem=false))
  Save(std::vector<double> samplingTimeVector,
       int noEnsembles, bool lowMem);

  //Each single simulation run/trajectory/"ensemble" should be saved
  // by calling this function. (using reference for speed)
  void store(const std::vector<int>& dx, const std::vector<int>& dy,
             const std::vector<int>& dz, const std::vector<double>& dr);

  //After simulation is complete, call this to print an output-file by name
  //"fileName" with structure:
  //   head-head-head-head-head-head-head-head
  //   [time]   [MSD]    [MSD_error]   [correlation]
  //where "head" is simulation info string to be stored with data in file.
  void save(std::string fileName, std::string head);

  //To get P(x,t) i.e. a histogram, printed to file: "fileName_histogram"
  //i.e. it automatically appends a "_histogram" to given output file name.
  void saveBinning(std::string fileName);

  //Add this after save() has run, to generate "numberOfRuns" synthetic
  //output files with name: "outNameN" where N is {0 < N < (numberOfRuns-1)}
  void computeBootstrap(std::string outName, int numberOfRuns,
                        std::string head);

  //compute correct error in the fitted parameter. Prints result to
  //standard out (i.e. terminal, unless piped)
  void computeJackknife(std::string outName);

  //see comment on k_t_ variable above.
  void setJumprate(double);

};

#endif
