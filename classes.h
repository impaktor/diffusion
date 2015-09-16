#ifndef CLASSES_H
#define CLASSES_H

#include <time.h>
#include <vector>

enum waitingtime {
  LIN = 0,
  EXP = 1,
  POW = 2,
};

// Just to print info to output file.
inline std::ostream & operator<<(std::ostream & str, waitingtime wt){
  switch (wt) {
  case LIN: return str << "lin";
  case EXP: return str << "exp";
  case POW: return str << "pow";
  default: return str << (int) wt;
  }
}


class Particle{
public:
  Particle(int, int, int);
  Particle(const Particle&);              //copy constructor
  Particle(void);                         //empty constructor
  int x;                                  //position x
  int y;
  int z;

  std::vector<double> jumprates;          //one element for each direction

  //keep track of how many laps it's made around the lattice
  //for periodic boundary condition. (only use this (temporary?) in
  //Lattice::moveAndBoundaryCheck and Lattice::move)
  // int windingNumber_x;    //can take both neg. and pos. numbers.
  // int windingNumber_y;
  // int windingNumber_z;

  //This allows us to print x,y,z with the << operator
  //like: Particle particle(1,4,2);  cout << particle <<endl;
  friend std::ostream& operator<< (std::ostream &out, Particle &particle);
  friend std::istream& operator>> (std::istream &in, Particle &particle);

  //assignment operator:
  Particle & operator=(const Particle &);

  void setPos(int, int, int);
};


//Class to show runtime information, such as remaining time
//until completion, and number of ensembles left.
class RemainingTime{
private:
  time_t startTime_;
  unsigned int totalEnsembles_;

public:
  RemainingTime(unsigned int);
  void printProgress(unsigned int);  //print progress to screen
};

#endif
