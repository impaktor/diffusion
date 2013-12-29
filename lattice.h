#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <string>

#include "nr/ran_mod.h"    //need this for the Ran-class below.
#include "auxiliary.h"     //non-physics stuff. (print messages etc.)
#include "classes.h"       //various data structures /classes (Jump, Particle)

using namespace std;

class Lattice {
public:
  Lattice(int x, int y, int z, int N, 
          double deed, bool boundaryFix);

  void place(void);                        //Place particles on the lattice
  void move(void);                         //pick a particle & direction based on jump-rate

  void generateTrajectory(int&);           //same as place() + move()

  int getDimension(void) const;
  void getDisplacement(vector<int>&, vector<int>&,
                       vector<int>&, vector<double>&) const;
  void setSamplingTimes(const vector<double>&, bool);
  void setInteraction(float);             //determine if we use the interaction-algorithm
  void setJumpNaka(float, float);
  void setDist(int, double);

  float computeNakazato(void);
  float computeErgodicity(int,int,int);   //Calculate the theoretical MSD when ergodic
  float computeAverageDiffusionConst(void);
  float computeEffectiveDiffusionConst(void);

  //Set jump-rate of each particle and direction
  void setJumpRate(const vector<Jump>&);


protected:
  Ran randomNumber;                        //generator accessible to all class members

  float interactionStrength_;              //only used if "InteractionOn=true"
  bool checkVacancyMatrix(void) const;     //test-function
  void countNeighbours(Particle, vector<int>&);

  vector<Particle> pos_;                   //XYZ-position for each particle on the lattice
  Particle pos_0_;                         //Initial position for the tagged particle
  int noParticles_;                        //Total number of particles
  int latticeX_, latticeY_, latticeZ_;     //Lattice size

  vector<double> samplingTime_;            //times to compute the MSD.
  int noSamplingTimes_;                    //number of times to sample

  vector<Jump> jumpRate_;                  //jump-rate for each particle and direction
  vector <long double> partialSum_;        //cumulative sum of jump-rates

  unsigned int dim_;                       //dimension of lattice (=1 || 2 || 3)
  bool boundaryFix_;                       //"true" if fix, "false" if periodic
  double timeSum_;                         //Sum of time for each move

  vector<int> dx_,dy_,dz_;                 //pos. of tracer particle for current ensemble
  vector<double> dr_;                      //pos. of tracer particle for current ensemble

  //count number of laps around the lattice (if periodic boundary):
  int windingNumber_0_x;
  int windingNumber_0_y;
  int windingNumber_0_z;

  void convertMuToParticle(int, int&, int&);
  void computePartialSum(void);
  void moveAndBoundaryCheck(int,int);

  int vacancyCheck(int, const Particle&); //Same as above, but improved. Use this one instead.
  vector<vector<vector<int> > > vacancy_; //Needed in my more efficient "vacancyCheck()"

  float jumpCrowders_, jumpTracer_;        //only used in the Nakazato()-function
  bool testOnOff_;                         //To print detailed information to screen

  int round(float x) {                     //Round off (used to get center of lattice)
    return int(x + 0.5);  }

  bool isExponentialWaitingTime_;
  double computeWaitingTime(void);


private:                                    //not shared with suerInteraction
  bool interactionOn_;                      //use the interaction algorithm
  void interaction(int, Particle);
  vector<float> clusterDistribution_, clusterSize_; //For our save_cluster()-function
  void buildCluster2(int, vector<int>&, double);    //needed for save_cluster()-func.
};


#endif
