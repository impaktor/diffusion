#ifndef BASELATTICE_H
#define BASELATTICE_H

#include <vector>

#include "nr/ran_mod.h"    //need this for the Ran-class below.
#include "auxiliary.h"     //non-physics stuff. (print messages etc.)
#include "classes.h"       //various data structures /classes (Particle)
#include <functional>      // for std::function

using std::vector;


class Vacancy{
protected:
  //representation of the geometry of our lattice:
  vector<vector<vector<int> > > vacancy_;

  //dimensions of lattice:
  int x_dim_, y_dim_, z_dim_;

public:
  Vacancy(int x, int y, int z){
    x_dim_ = x;   //dimension of lattice
    y_dim_ = y;
    z_dim_ = z;
  }

  //initialize the entire vacancy matrix (or what ever
  //representation)
  void initialize(const vector<Particle> &particles){

    //set all sites at vacant i.e. value = -1. A site with a
    //particle will have that particles label as value. [0 -- (N-1)]
    vacancy_.assign(x_dim_+1, vector<vector<int> >
                    (y_dim_+1, vector<int>(z_dim_+1,-1)));

    // Matrix starts at coord: (1,1,1) the 0-elements of the matrix
    //have no meaning, since our lattice starts at 1,..,x_dim_
    for(size_t n = 0; n < particles.size(); n++)
      vacancy_[ particles[n].x ][ particles[n].y ][ particles[n].z ] = n;
  }

  //check if vacant, return true or false (or label id of occupying
  //particle?)
  bool isVacant(Particle p) const{
    if(vacancy_[p.x][p.y][p.z] == -1)
      return true;
    else
      return false;
  }


  bool isVacant(int x, int y, int z) const{
    return isVacant(Particle(x,y,z));
  }

  bool isOccupied(Particle p) const{
    return !(isVacant(p));
  }

  bool isOccupied(int x, int y, int z) const{
    return !(isVacant(x, y, z));
  }

  int get_label(Particle p) const{
    return vacancy_[p.x][p.y][p.z];
  }

  int get_label(int x, int y, int z) const{
    return get_label(Particle(x,y,z));
  }

  //set site as occupied by particle with index id.
  void setAsOccupied(Particle p, int id){
    vacancy_[p.x][p.y][p.z] = id;
  }

  void setAsVacant(Particle p){
    setAsOccupied(p, -1);
  }


  bool checkOK(const vector<Particle> &p) const{
    //Just for debugging: Check that the coordinates in p-vector is
    //correct with labels in vacancy_:

    bool isCorrect = true;  //if no error

    int label;
    unsigned int count_vacant_sites = 0;

    //Loop through vacancy_ and check that the coordinates for each
    //occupied site corresponds to the right label and coordinates
    //in particles-vector.
    for (int i = 1; i <= x_dim_; i++){
      for (int j = 1; j <= y_dim_; j++)
        for (int k = 1; k <= z_dim_; k++){
          if(isOccupied(i,j,k)){
            label = get_label(i,j,k);
            if(i != p[label].x || j != p[label].y || k != p[label].z)
              isCorrect = false;
          }
          else
            count_vacant_sites++;
        }
    }

    unsigned int computeVacant = x_dim_*y_dim_*z_dim_ - p.size();
    if(count_vacant_sites != computeVacant){
      isCorrect = false;
      throw std::string( "number of vacant sites: " + tostring(count_vacant_sites));
    }

    return isCorrect;
  }


  void checkNeighbours(Particle particle, vector<int>& neighbours, bool isBoundaryFix){
    //Only needed if using some kind of interactions.  Return a vector
    //of labels of all neighbors to "particle". If "particle" is
    //alone, the vector will be empty.

    //isBoundaryFix_=0 periodic, =1 fix

    //Lambda-function, [&] = know about neighbours-vector& and bool.
    std::function<void (int, int, int)> check =
      [&](int x, int y, int z){ if(isOccupied(x, y, z))
        neighbours.push_back(get_label(x, y, z));
    };

    std::function<void (int, int, int)> checkBoundary =
      [&](int x, int y, int z){ if(isOccupied(x, y, z) && !isBoundaryFix)
        neighbours.push_back(get_label(x, y, z));
    };


    if (particle.x+1 <= x_dim_ ) //avoid going outsied the lattice
      check(particle.x, particle.y, particle.z);
    else                            //if on the edge then:
      checkBoundary(1, particle.y, particle.z);

    if (particle.x-1>0)
      //although this zero-element exists in the vacancy_ matrix
      //(i.e. no seg. fault) the first coordinate on the lattice is
      //(1,1,1). (0,0,0) is undefined
      check(particle.x - 1, particle.y, particle.z);
    else
      checkBoundary(x_dim_, particle.y, particle.z);


    if (particle.y+1 <= y_dim_)
      check(particle.x, particle.y + 1, particle.z);
    else
      checkBoundary(particle.x, 1, particle.z);

    if (particle.y-1>0)
      check(particle.x, particle.y - 1, particle.z);
    else
      checkBoundary(particle.x, y_dim_, particle.z);


    if (particle.z+1 <= z_dim_)
      check(particle.x, particle.y, particle.z + 1);
    else
      checkBoundary(particle.x, particle.y, 1);


    if (particle.z-1>0)
      check(particle.x, particle.y, particle.z - 1);
    else
      checkBoundary(particle.x, particle.y, z_dim_);
  }
};



class BaseLattice{
public:
  BaseLattice(int x, int y, int z, int dimension,
              double seed, bool boundaryFix);

  void place(void);                        //Place particles on the lattice
  void move(void);                         //pick a particle & direction based on jump-rate

  unsigned int getDimension(void) const;
  unsigned int getDirections(void) const;
  void getDisplacement(vector<int>&, vector<int>&,
                       vector<int>&, vector<double>&) const;
  void setSamplingTimes(const vector<double>&, waitingtime);
  void setInteraction(float);             //determine if we use the interaction-algorithm
  void setParticles(const vector<Particle>&); //add particles to lattice
  void setJumpNaka(float, float);
  void setDist(int, double);

  float computeNakazato(void);
  float computeErgodicity(int);            //Calculate the theoretical MSD when ergodic
  float computeAverageDiffusionConst(void);
  float computeEffectiveDiffusionConst(void);


protected:
  Ran randomNumber;                        //generator accessible to all class members

  Vacancy board_;                          //store all particle positions in a vector, for fast lookup of vacant/occupied sites.

  float interactionStrength_;              //only used if "InteractionOn=true"

  struct dir{                              // center of lattice, starting pos for tagged particle
    int x, y, z;
  } center_;

  vector<Particle> pos_;                   //XYZ-position for each particle on the lattice
  int noParticles_;                        //Total number of particles
  int latticeX_, latticeY_, latticeZ_;     //Lattice size

  vector<double> samplingTime_;            //times to compute the MSD.
  int noSamplingTimes_;                    //number of times to sample

  vector <long double> partialSum_;        //cumulative sum of jump-rates

  unsigned int dim_;                       //dimension of lattice (=1 || 2 || 3)
  unsigned int directions_;                // directions on lattice
  bool isBoundaryFix_;                     //"true" if fix, "false" if periodic
  double timeSum_;                         //Sum of time for each move

  vector<int> dx_,dy_,dz_;                 //pos. of tracer particle for current ensemble
  vector<double> dr_;                      //pos. of tracer particle for current ensemble

  //count number of laps around the lattice (if periodic boundary):
  int windingNumber_0_x;
  int windingNumber_0_y;
  int windingNumber_0_z;

  void convertMuToParticle(int, int&, int&);
  void computePartialSum();                // used to find direction to move in

  virtual void moveAndBoundaryCheck(int,int) = 0;
  virtual double distance(int &dx, int &dy, int &dz) = 0;

  int vacancyCheck(int, const Particle&);

  float jumpCrowders_, jumpTracer_;        //only used in the Nakazato()-function
  bool isTestOn_;                          //To print detailed information to screen

  int round(float x) {                     //Round off (used to get center of lattice)
    return int(x + 0.5);  }

  waitingtime tagged_waiting_time_;        // waiting time (lin,exp,pow) of the tagged particle
  double computeWaitingTime(void);

protected:                                  //not shared with suerInteraction
  bool isInteractionOn_;                    //use the interaction algorithm
  void interaction(int, Particle);
};

#endif
