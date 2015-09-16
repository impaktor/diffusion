#ifndef LATTICES_H
#define LATTICES_H

#include "nr/ran_mod.h"    //need this for the Ran-class below.
#include "auxiliary.h"     //non-physics stuff. (print messages etc.)
#include "baselattice.h"   // the basic methods for a lattice


class Square: public BaseLattice{
public:

  Square(int x, int y, int z, int dimension, double seed, bool boundaryFix)
    //initiate the class by calling the base class
    : BaseLattice(x,y,z, dimension, seed, boundaryFix)
  {
  }

  // how to compute the distance on a lattice from the three new
  // coordinates in pos_[0], and the starting position in center_. Just
  // for reference, we return the displacement in x,y,z as references
  // as well, but the important one is the return of the distance "dr".
  double distance(int &dx, int &dy, int &dz){
    dx = (pos_[0].x + windingNumber_0_x * latticeX_) - center_.x;
    dy = (pos_[0].y + windingNumber_0_y * latticeY_) - center_.y;
    dz = (pos_[0].z + windingNumber_0_z * latticeZ_) - center_.z;
    return sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
  }

  void moveAndBoundaryCheck(int n, int R);
};


class Honeycomb2d: public BaseLattice{
public:

  Honeycomb2d(int squares, double seed, bool boundaryFix)
    //initiate the class by calling the base class
    : BaseLattice(squares, squares, 1, 2, seed, boundaryFix)
  {
    directions_ = 6;  // a 2D lattice, with x,y coordinates, but has 6 directions
  }

  // how to compute the distance on a lattice from three coordinates
  double distance(int &dx, int &dy, int &dz){
    dx = pos_[0].x + windingNumber_0_x * latticeX_;   // new positions in x
    dy = pos_[0].y + windingNumber_0_y * latticeY_;
    dz = 0;
    return sqrt(pow((dx - center_.x) + 0.5*(dy - center_.y), 2) + 3.0/4.0*pow(dy - center_.y, 2));
  }

  void moveAndBoundaryCheck(int n, int R);
};

#endif
