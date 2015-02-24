#include <iostream>        //for cout, among others
#include <cmath>           //gives sqrt, pow, fabs,
#include <vector>

#include "classes.h"       //various data structures /classes (Jump, Particle)
#include "lattices.h"       //this is the main class, that does the physics.


void Square::moveAndBoundaryCheck(int n, int R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant.

  //windingNumber was added in April 7th 2011, to count the number of laps
  //the particle has made around the lattice/edge.

  if (0 <= n && n < noParticles_ ){
    if (isTestOn_) std::cout << "Moving particle " << n << std::endl;
    // n = particle to move, R = direction to move, such that:
    //      R= 0 = X-right,      1 = X-left
    //         2 = Y-right,      3 = Y-left
    //         4 = Z-up,         5 = Z-down
    // boundary: boundary=0 = periodic, boundary=1 = fix boundary

    //Store in case we want to restore the previous configuration
    Particle old = pos_[n];

    //If I need to restore previous position across a boundary
    //I need to know this, (see further down)
    bool changingWindingNumber = false;

    if ( R >= 0 && R < 6 ){
      switch(R){
      case 0: pos_[n].x = pos_[n].x +1;
        if (pos_[n].x == latticeX_ +1){
            pos_[n].x = isBoundaryFix_ ? old.x : 1;

            //if the tagged particle...
            if (n==0 && !isBoundaryFix_){
              windingNumber_0_x++;
              changingWindingNumber = true;
            }
        }
        break;
      case 1: pos_[n].x = pos_[n].x -1; //move particle n to the left
        if (pos_[n].x == 0){
          pos_[n].x = isBoundaryFix_ ? old.x : latticeX_;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_x--;
            changingWindingNumber = true;
          }
        }
        break;                          //(first coordinate is 1, 0=outside)
      case 2: pos_[n].y = pos_[n].y +1;
        if (pos_[n].y == latticeY_ +1){
          pos_[n].y = isBoundaryFix_ ? old.y : 1;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_y++;
            changingWindingNumber = true;
          }
        }
        break;
      case 3: pos_[n].y = pos_[n].y -1;
        if (pos_[n].y == 0 ){
          pos_[n].y = isBoundaryFix_ ? old.y : latticeY_;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_y--;
            changingWindingNumber = true;
          }
        }
        break;
      case 4: pos_[n].z = pos_[n].z +1;
        if (pos_[n].z == latticeZ_+1){
          pos_[n].z = isBoundaryFix_ ? old.z : 1;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_z++;
            changingWindingNumber = true;
          }
        }
        break;
      case 5: pos_[n].z = pos_[n].z-1;
        if (pos_[n].z == 0){
          pos_[n].z = isBoundaryFix_ ? old.z : latticeZ_;

          if (n==0 && isBoundaryFix_ == 0){
              windingNumber_0_z--;
              changingWindingNumber = true;
          }
        }
        break;
      }

      //check that the new site is vacant if not, move the particle
      //back. (done by the vacancyCheck...)  (and if using
      //windingnumber, reset it.)

      bool resetWindingNumber = vacancyCheck(n,old);

      //now reset windingnumber if I'm colliding at the boundary:
      if (resetWindingNumber && changingWindingNumber){
        switch(R) {
        case 0:
          windingNumber_0_x--;
          break;
        case 1:
          windingNumber_0_x++;
          break;
        case 2:
          windingNumber_0_y--;
          break;
        case 3:
          windingNumber_0_y++;
          break;
        case 4:
          windingNumber_0_z--;
          break;
        case 5:
          windingNumber_0_z++;
          break;
        }
      }

    }
    else{
      std::cout << "There is a maximum of " << 2 * dim_
                << " positive directions to move in" << std::endl;
      abort();
    }
  }
  else{
    std::cout << "Invalid particle label, n=" << n
              << ", must fulfill 0 <= n < N" << std::endl;
    abort();
  }
}


void Honeycomb2d::moveAndBoundaryCheck(int n, int R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant.

  //windingNumber was added in April 7th 2011, to count the number of laps
  //the particle has made around the lattice/edge.

  if (0 <= n && n < noParticles_ ){
    if (isTestOn_) std::cout << "Moving particle " << n << std::endl;
    // Honeycomb geometry:
    // n = particle to move, R = direction to move, such that:
    //      R= 0 = X-right,      1 = X-left
    //         2 = Y-right,      3 = Y-left
    //         4 = diagonal up,  5 = diagonal down
    // boundary: boundary=0 = periodic, boundary=1 = fix boundary

    //Store in case we want to restore the previous configuration
    Particle old = pos_[n];

    //If I need to restore previous position across a boundary
    //I need to know this, (see further down)
    bool changingWindingNumber_x = false;
    bool changingWindingNumber_y = false;

    if ( R >= 0 && R < 6 ){
      switch(R){
      case 0: pos_[n].x = pos_[n].x +1;
        if (pos_[n].x == latticeX_ +1){
            pos_[n].x = isBoundaryFix_ ? old.x : 1;

            //if the tagged particle...
            if (n==0 && !isBoundaryFix_){
              windingNumber_0_x++;
              changingWindingNumber_x = true;
            }
        }
        break;
      case 1: pos_[n].x = pos_[n].x -1; //move particle n to the left
        if (pos_[n].x == 0){
          pos_[n].x = isBoundaryFix_ ? old.x : latticeX_;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_x--;
            changingWindingNumber_x = true;
          }
        }
        break;                          //(first coordinate is 1, 0=outside)
      case 2: pos_[n].y = pos_[n].y +1;
        if (pos_[n].y == latticeY_ +1){
          pos_[n].y = isBoundaryFix_ ? old.y : 1;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_y++;
            changingWindingNumber_y = true;
          }
        }
        break;
      case 3: pos_[n].y = pos_[n].y -1;
        if (pos_[n].y == 0 ){
          pos_[n].y = isBoundaryFix_ ? old.y : latticeY_;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_y--;
            changingWindingNumber_y = true;
          }
        }
        break;
      case 4:
        pos_[n].x = pos_[n].x -1;  // Honeycomb special: move diagonal up left
        pos_[n].y = pos_[n].y +1;
        if (pos_[n].x == 0){
          pos_[n].x = isBoundaryFix_ ? old.x : latticeX_;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_x--;
            changingWindingNumber_x = true;
          }
        }
        if (pos_[n].y == latticeY_ +1){
          pos_[n].y = isBoundaryFix_ ? old.y : 1;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_y++;
            changingWindingNumber_y = true;
          }
        }
        break;
      case 5:
        pos_[n].x = pos_[n].x +1;  // Honeycomb special: move diagonal down right
        pos_[n].y = pos_[n].y -1;
        if (pos_[n].x == latticeX_ +1){
          pos_[n].x = isBoundaryFix_ ? old.x : 1;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_x++;
            changingWindingNumber_x = true;
          }
        }
        if (pos_[n].y == 0 ){
          pos_[n].y = isBoundaryFix_ ? old.y : latticeY_;

          if (n==0 && !isBoundaryFix_){
            windingNumber_0_y--;
            changingWindingNumber_y = true;
          }
        }
        break;
      }

      //check that the new site is vacant if not, move the particle
      //back. (done by the vacancyCheck...)  (and if using
      //windingnumber, reset it.)

      bool resetWindingNumber = vacancyCheck(n,old);

      //now reset windingnumber if I'm colliding at the boundary:
      if (resetWindingNumber && (changingWindingNumber_x || changingWindingNumber_y)){
        switch(R) {
        case 0:
          windingNumber_0_x--;
          break;
        case 1:
          windingNumber_0_x++;
          break;
        case 2:
          windingNumber_0_y--;
          break;
        case 3:
          windingNumber_0_y++;
          break;
        case 4:
          if(changingWindingNumber_x)
            windingNumber_0_x++;
          if(changingWindingNumber_y)
            windingNumber_0_y--;
          break;
        case 5:
          if(changingWindingNumber_x)
            windingNumber_0_x--;
          if(changingWindingNumber_y)
            windingNumber_0_y++;
          break;
        }
      }

    }
    else{
      std::cout << "There is a maximum of " << 2 * dim_
                << " positive directions to move in" << std::endl;
      abort();
    }
  }
  else{
    std::cout << "Invalid particle label, n=" << n
              << ", must fulfill 0 <= n < N" << std::endl;
    abort();
  }
}
