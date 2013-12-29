#ifndef LATTICE_H
#define LATTICE_H

#include <vector>

#include "auxiliary.h"      //non-physics stuff. (print messages etc.)
#include "classes.h"       //various data structures /classes (Jump, Particle)

using namespace std;

class Lattice {
private:
  bool interactionOn_;                      //use the interaction algorithm
  void interaction(int, Particle);
  vector<float> clusterDistribution_, clusterSize_; //For our save_cluster()-function
  void buildCluster2(int, vector<int>&, double);    //needed for save_cluster()-func.

protected:
  float interactionStrength_;               //only used if "InteractionOn=true"
  void checkVacancyMatrix(float);           //test-function
  void countNeighbours(Particle, vector<int>&);

  vector<Particle> pos_;                   //XYZ-position for each particle on the lattice
  Particle pos_0_;                         //Initial position for the tagged particle
  int nParticles_;                         //Total number of particles
  int latticeX_, latticeY_, latticeZ_;     //Lattice size

  vector<long double> x2_err,y2_err,z2_err,r2_err; //Store std error, in each point
  vector<double> x2_mu,y2_mu,z2_mu,r2_mu;  //average of displacement squared: <X^2>,...
  vector<double> samplingTime_;            //times to compute the MSD.
  int maxElement_;                         //number of times to sample
  vector<double> fourthColumn_;            //print a fourth column of inf. to the out-put file
  //vector<double> pearsonCorrelation_;      //correlation between two neighboring sample points.


  vector<Jump> jumpRate_;                  //jump-rate for each particle and direction
  vector <long double> partialSum_;        //cumulative sum of jump-rates

  int dim_;                                //dimension of lattice (=1 || 2 || 3)
  bool boundaryFix_;                       //"true" if fix, "false" if periodic
  float concentration_;                    //particle concentration
  double timeSum_;                         //Sum of time for each move

  //TODO------------
  vector<int> dx,dy,dz;                    //pos. of tracer particle for current ensemble 
  vector<double> dr;                       //pos. of tracer particle for current ensemble 
  
  vector< vector<int> > store_dx;          //save dx for each ensemble to compute stdErr 
  vector< vector<int> > store_dy;
  vector< vector<int> > store_dz;
  vector< vector<double> > store_dr;
  //------------

  int convertMuToParticle(int, int&);
  void computePartialSum(void);
  void moveAndBoundaryCheck(int,int);

  int vacancyCheckOld(int, Particle);     //No double Occupancy.
  int vacancyCheck(int, Particle);        //Same as above, but improved. Use this one instead.
  vector<vector<vector<int> > > vacancy_; //Needed in my more efficient "vacancyCheck()"


  void computeStdErr(int);                //calculate standard deviation
  void computeStdErrLowMem(int);          //a more memmory conservative version of the above
  bool lowMem_;                           //Switch between the two stdErr-functions above. 
  vector<double> dr4_err,dx4_err,dy4_err,dz4_err; //needed if we use stdErrLowMem()

  float jumpCrowders_, jumpTracer_;        //only used in the Nakazato()-function
  int distribution_;
  double info_;                             //save characteristic trait of particle distribution
  bool testOnOff_;                          //To print detailed information to screen
  double D_av_; 


  float nakazato(void);         
  double diffEffConstant(void);            //Calculate the effective dif. const. from the jump-rates
  float ergodicity(int,int,int);           //Calculate the theoretical MSD when ergodic

  int round(float x) {                     //Round off (used to get center of lattice)
    return int(x + 0.5);  }

  bool isExponentialWaitingTime_;
  double computeWaitingTime(void);
  int computeCorrelation(int);             //correlation between two sampling points

public:
  Lattice(int, int, int, int, bool);
  void place(void);                        //Place particles on the lattice
  void moveOld(void);                      //old obsolete move-code (no jump-rates)
  void move(void);                         //pick a particle & direction based on jump-rate
  void store(void);                        //Tell to store result
  void save(int, char[]);                  //Print simulation-run to file   
  int getDimension(void);
  void setSamplingTimes(vector<double>&, bool);
  void setInteraction(float);             //determine if we use the interaction-algorithm
  void setLowMem(bool);                   //Set the switch true/false for use of low-mem run
  void setJumpNaka(float, float);
  void setDist(int, double);

  //Set jump-rate of each particle and direction 
  void setJumpRate(vector<Jump>&);

  //Secondary not used functions
  void snapshot(void);                     //Not used. (more comments in function body)
  //  double get_dt(void);   //Size of time-step, need it for P(x,t), don't need this...

  //INTERACTION CODE:
  void saveBinning(int);                   //To get P(x,t)
  void saveCluster(double);                //To get rho(m) (cluster distribution)
  void printCluster(void);                 //Print  rho(m) (cluster distribution)

  //TEST-Code
  void dumpSimulation(int);                //dump displacement for individual run to file. 

};


#endif
