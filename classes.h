#ifndef CLASSES_H
#define CLASSES_H

#include <time.h>

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

//First three classes (Direction, Jump, Particle), will allow us to
//get the position of a Particle object as: object.x, and the jump
//rate as object.rate.x.r etc.

class Jump{
public:
  class Direction{
  public:
    Direction(double right, double left);
    Direction(const Direction&);            //copy constructor
    Direction(void);                        //empty constructor
    double r;  //right
    double l;  //left
  };

  // Jump(double xr, double xl, double yr, double yl, double zr, double zl) :
  //   x(xr,xl), y(yr,yl), z(zr,zl){ }
  Direction x; //i.e. x.r and x.l
  Direction y;
  Direction z;
};

class Particle{
public:
  Particle(int, int, int);
  Particle(const Particle&);              //copy constructor
  Particle(void);                         //empty constructor
  int x;                                  //position x
  int y;
  int z;

  //keep track of how many laps it's made around the lattice
  //for periodic boundary condition. (only use this (temporary?) in
  //Lattice::moveAndBoundaryCheck and Lattice::move)
  // int windingNumber_x;    //can take both neg. and pos. numbers.
  // int windingNumber_y;
  // int windingNumber_z;

  Jump rate;        //jump rate: rate.x.r

  //This allows us to print x,y,z with the << operator
  //like: Particle particle(1,4,2);  cout << particle <<endl;
  friend std::ostream& operator<< (std::ostream &out, Particle &particle);
  friend std::istream& operator>> (std::istream &in, Particle &particle);

  //assignment operator:
  Particle & operator=(const Particle &);

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
