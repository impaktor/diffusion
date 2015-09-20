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
  double distance(){
    int dx = true_x - center_.x;
    int dy = true_y - center_.y;
    int dz = true_z - center_.z;
    return sqrt(pow(dx,2) + pow(dy,2) + pow(dz,2));
  }

  void moveAndBoundaryCheck(size_t n, size_t R);
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
  double distance(){
    int dx = true_x - center_.x;
    int dy = true_y - center_.y;
    return sqrt(pow(dx + 0.5*dy, 2) + 3.0/4.0*pow(dy, 2));
  }

  void moveAndBoundaryCheck(size_t n, size_t R);
};


class BCC: public BaseLattice{
public:

  BCC(int squares, double seed, bool boundaryFix)
    //initiate the class by calling the base class
    : BaseLattice(squares, squares, squares, 3, seed, boundaryFix)
  {
    directions_ = 8;  // a 3D Body-centered cubic lattice, with x,y,z coordinates, but has 8 directions
  }

  // how to compute the distance on a lattice from three coordinates
  double distance(){

    int dx = true_x - center_.x;
    int dy = true_y - center_.y;
    int dz = true_z - center_.z;

    return 0.5*sqrt(pow(  dx + dy - dz, 2) +
                    pow(  dx - dy + dz, 2) +
                    pow(- dx + dy + dz, 2));
  }

  void moveAndBoundaryCheck(size_t n, size_t R);
};


class FCC: public BaseLattice{
public:

  FCC(int squares, double seed, bool boundaryFix)
    //initiate the class by calling the base class
    : BaseLattice(squares, squares, squares, 3, seed, boundaryFix)
  {
    directions_ = 12;  // a 3D Face-centered cubic lattice, with x,y,z coordinates, but has 12 directions
  }

  // how to compute the distance on a lattice from three coordinates
  double distance(){

    int dx = true_x - center_.x;
    int dy = true_y - center_.y;
    int dz = true_z - center_.z;

    return 0.5*sqrt(pow(dx + dy, 2) +
                    pow(dx + dz, 2) +
                    pow(dy + dz, 2));
  }

  void moveAndBoundaryCheck(size_t n, size_t R);
};


#endif
