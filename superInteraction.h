#ifndef SUPERINTERACTION_H
#define SUPERINTERACTION_H

#include <vector>
#include <functional>      // for std::function
#include "lattices.h"
#include "classes.h"

//the class SuperInteraction inherits all the public and protected
//methods and variables of the lattice-classs.
class SuperInteraction : public Square{
private:
  float interactionStrength_;

  void superInteractionCode(int, int);

  //sum over many particles
  double sumFrictionCoefficient(std::vector<int>) const;

  //sum over one particle
  double sumFrictionCoefficient(int) const;

  void calculateExtent(std::vector<int>&, int, int&, int&);

  void buildBlockingCluster(std::vector<int>&, std::vector<int>, int);
  void buildCluster2(int, std::vector<int>&, double);

  //TODO: import these to my lattice class, and use them there as well:
  void printPositionToScreen(std::vector<int>);
  void printPositionToScreen(int);

  void controlVacancyCheck3(int);

  void moveAndBoundaryCheck(int, int);   //overload function.

public:
  SuperInteraction(int,int,int, int, double, bool, float);
  void move(void);              //overload the lattice-class version

  //don't actually use this one, thus redefine it.
  void setInteraction(float) const;

  //for easier threading:
  void generateTrajectory(int& E);

};

#endif
