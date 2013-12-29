#ifndef SUPERINTERACTION_H
#define SUPERINTERACTION_H

#include <vector>
#include "lattice.h"
#include "classes.h"

using namespace std;

//the class SuperInteraction inherits all the public and protected 
//methods and variables of the lattice-classs.
class SuperInteraction : public Lattice{
private:
  float interactionStrength_;
  
  void superInteractionCode(int, int);

  double SumFrictionCoefficient(vector<int>);  //sum over many particles
  double SumFrictionCoefficient(int);          //sum over one particle

  void calculateExtent(vector<int>&, int, int&, int&); 

  void buildBlockingCluster(vector<int>&, vector<int>, int);
  void buildCluster2(int, vector<int>&, double);

  //TODO: import these to my lattice class, and use them there as well:
  void printPositionToScreen(vector<int>);
  void printPositionToScreen(int);

  void controlVacancyCheck3(int );

  void moveAndBoundaryCheck(int, int);   //override function.


public:
  SuperInteraction(int,int,int, int, bool, float);
  void move(void);  //override the lattice-class version
  void setInteraction(float); //don't actually use this one, thus redefine it.

};

//TODO hur funkar boundary-booleanen? 

#endif
