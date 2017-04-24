#include <iostream>        //for cout, among others
#include <cmath>           //gives sqrt, pow, fabs,
#include <vector>

#include "classes.h"       //various data structures /classes (Jump, Particle)
#include "lattices.h"       //this is the main class, that does the physics.


void Square::moveAndBoundaryCheck(size_t n, size_t R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant.

  // Square geometry:
  // n = particle to move, R = direction to move, such that:
  //      R= 0 = X-right,      1 = X-left
  //         2 = Y-right,      3 = Y-left
  //         4 = Z-pos,        5 = Z-in
  // boundary: boundary=0 = periodic, boundary=1 = fix boundary

  if (n > noParticles_ -1)
    throw std::string("Invalid particle label, n=" + tostring(n) +
                      ", must fulfill 0 <= n < N");

  if (R > directions_ -1)
    throw std::string("There is a maximum of " + tostring(directions_) +
                      " positive directions to move in");

  if (isTestOn_) std::cout << "Moving particle " << n << std::endl;

  // if particle is outside and fix boundary, we know we should move it back
  bool isOutside = false;

  //Store in case we want to restore the previous configuration
  Particle old = pos_[n];

  // Move particle n on the lattice:
  switch(R){
  case 0:
    pos_[n].x = pos_[n].x +1;      // jump right in lattice system

    if(pos_[n].x == latticeX_ +1){ // if jumping outside the lattice:
      isOutside = true;            // if fix boundary, reset particle (later)
      if(!isBoundaryFix_) pos_[n].x = 1; // if periodic boundary, warp around lattice
    }
    break;
  case 1:
    pos_[n].x = pos_[n].x -1;      //move particle n to the left

    if(pos_[n].x == 0){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].x = latticeX_;
    }
    break;
  case 2:
    pos_[n].y = pos_[n].y +1;

    if(pos_[n].y == latticeY_ +1){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].y = 1;
    }
    break;
  case 3:
    pos_[n].y = pos_[n].y -1;

    if(pos_[n].y == 0){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].y = latticeY_;
    }
    break;
  case 4:
    pos_[n].z = pos_[n].z +1;

    if(pos_[n].z == latticeZ_ +1){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].z = 1;
    }
    break;
  case 5:
    pos_[n].z = pos_[n].z -1;

    if(pos_[n].z == 0){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].z = latticeZ_;
    }
    break;
  default:
    throw std::string("move in direction" + tostring(R) + "not implemented");
  }

  if(isOutside && isBoundaryFix_)
    pos_[n] = old;

  //check that the new site is vacant if not, move the particle
  //back. (done by the vacancyCheck...)
  bool collided = vacancyCheck(n,old);

  //If n is tagged particle, then also move it in its second,
  //"infinite" coordinate system, - if it didn't collied with any
  //particle on our lattice above.
  if (n == 0 && !collided && !(isBoundaryFix_ && isOutside)){
    switch(R) {
    case 0:
      true_x++;
      break;
    case 1:
      true_x--;
      break;
    case 2:
      true_y++;
      break;
    case 3:
      true_y--;
      break;
    case 4:
      true_z++;
      break;
    case 5:
      true_z--;
      break;
    default:
      throw std::string("move in direction" + tostring(R) + "not implemented");
    }
  }
}


void Honeycomb2d::moveAndBoundaryCheck(size_t n, size_t R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant.

  // Honeycomb geometry:
  // n = particle to move, R = direction to move, such that:
  //      R= 0 = X-right,      1 = X-left
  //         2 = Y-right,      3 = Y-left
  //         4 = diagonal up,  5 = diagonal down
  // boundary: boundary=0 = periodic, boundary=1 = fix boundary

  if (n > noParticles_ -1)
    throw std::string("Invalid particle label, n=" + tostring(n) +
                      ", must fulfill 0 <= n < N");

  if (R > directions_ -1)
    throw std::string("There is a maximum of " + tostring(directions_) +
                      " positive directions to move in");

  if (isTestOn_) std::cout << "Moving particle " << n << std::endl;

  // if particle is outside and fix boundary, we know we should move it back
  bool isOutside = false;

  //Store in case we want to restore the previous configuration
  Particle old = pos_[n];

    // Move particle n on the lattice:
  switch(R){
  case 0:
    pos_[n].x = pos_[n].x +1;

    if(pos_[n].x == latticeX_ +1){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].x = 1;
    }
    break;
  case 1:
    pos_[n].x = pos_[n].x -1;

    if(pos_[n].x == 0){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].x = latticeX_;
    }
    break;
  case 2:
    pos_[n].y = pos_[n].y +1;

    if(pos_[n].y == latticeY_ +1){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].y = 1;
    }
    break;
  case 3:
    pos_[n].y = pos_[n].y -1;

    if(pos_[n].y == 0){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].y = latticeY_;
    }
    break;
  case 4:
    pos_[n].x = pos_[n].x -1;  // Honeycomb special: move diagonal up left
    pos_[n].y = pos_[n].y +1;

    if(pos_[n].x == 0 || pos_[n].y == latticeY_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == 0)            pos_[n].x = latticeX_;
        if (pos_[n].y == latticeY_ +1) pos_[n].y = 1;
      }
    }
    break;
  case 5:
    pos_[n].x = pos_[n].x +1;  // Honeycomb special: move diagonal down right
    pos_[n].y = pos_[n].y -1;

    if(pos_[n].x == latticeX_ +1 || pos_[n].y == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == latticeX_ +1) pos_[n].x = 1;
        if (pos_[n].y == 0)            pos_[n].y = latticeY_;
      }
    }
    break;
  default:
    throw std::string("move in direction" + tostring(R) + "not implemented");
  }

  if(isOutside && isBoundaryFix_)
    pos_[n] = old;

  //check that the new site is vacant if not, move the particle
  //back. (done by the vacancyCheck...)
  bool collided = vacancyCheck(n,old);

  //If n is tagged particle, then also move it in its second,
  //"infinite" coordinate system, - if it didn't collied with any
  //particle on our lattice above.
  if (n == 0 && !collided && !(isBoundaryFix_ && isOutside)){
    switch(R) {
    case 0:
      true_x++;
      break;
    case 1:
      true_x--;
      break;
    case 2:
      true_y++;
      break;
    case 3:
      true_y--;
      break;
    case 4:
      true_x--;
      true_y++;
      break;
    case 5:
      true_x++;
      true_y--;
      break;
    default:
      throw std::string("move in direction" + tostring(R) + "not implemented");
    }
  }
}


void BCC::moveAndBoundaryCheck(size_t n, size_t R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant.

  // Body centered cubic geometry:
  // n = particle to move, R = direction to move
  // boundary: boundary=0 = periodic, boundary=1 = fix boundary

  if (n > noParticles_ -1)
    throw std::string("Invalid particle label, n=" + tostring(n) +
                      ", must fulfill 0 <= n < N");

  if (R > directions_ -1)
    throw std::string("There is a maximum of " + tostring(directions_) +
                      " positive directions to move in");

  if (isTestOn_) std::cout << "Moving particle " << n << std::endl;

  // if particle is outside and fix boundary, we know we should move it back
  bool isOutside = false;

  //Store in case we want to restore the previous configuration
  Particle old = pos_[n];

  // Move particle n on the lattice:
  switch(R){
  case 0:
    pos_[n].x = pos_[n].x +1;
    pos_[n].y = pos_[n].y +1;
    pos_[n].z = pos_[n].z +1;

    if (pos_[n].x == latticeX_ +1 || pos_[n].y == latticeY_ +1 || pos_[n].z == latticeZ_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == latticeX_ +1) pos_[n].x = 1;
        if (pos_[n].y == latticeY_ +1) pos_[n].y = 1;
        if (pos_[n].z == latticeZ_ +1) pos_[n].z = 1;
      }
    }
    break;
  case 1:
    pos_[n].x = pos_[n].x -1;
    pos_[n].y = pos_[n].y +1;
    pos_[n].z = pos_[n].z +1;

    if (pos_[n].x == 0 || pos_[n].y == latticeY_ +1 || pos_[n].z == latticeZ_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == 0)            pos_[n].x = latticeX_;
        if (pos_[n].y == latticeY_ +1) pos_[n].y = 1;
        if (pos_[n].z == latticeZ_ +1) pos_[n].z = 1;
      }
    }
    break;
  case 2:
    pos_[n].x = pos_[n].x +1;
    pos_[n].y = pos_[n].y -1;
    pos_[n].z = pos_[n].z +1;

    if (pos_[n].x == latticeX_ +1 || pos_[n].y == 0 || pos_[n].z == latticeZ_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == latticeX_ +1) pos_[n].x = 1;
        if (pos_[n].y == 0)            pos_[n].y = latticeY_;
        if (pos_[n].z == latticeZ_ +1) pos_[n].z = 1;
      }
    }
    break;
  case 3:
    pos_[n].x = pos_[n].x +1;
    pos_[n].y = pos_[n].y +1;
    pos_[n].z = pos_[n].z -1;

    if (pos_[n].x == latticeX_ +1 || pos_[n].y == latticeY_ +1 || pos_[n].z == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == latticeX_ +1) pos_[n].x = 1;
        if (pos_[n].y == latticeY_ +1) pos_[n].y = 1;
        if (pos_[n].z == 0)            pos_[n].z = latticeZ_;
      }
    }
    break;
  case 4:
    pos_[n].x = pos_[n].x +1;
    pos_[n].y = pos_[n].y -1;
    pos_[n].z = pos_[n].z -1;

    if (pos_[n].x == latticeX_ +1 || pos_[n].y == 0 || pos_[n].z == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == latticeX_ +1) pos_[n].x = 1;
        if (pos_[n].y == 0)            pos_[n].y = latticeY_;
        if (pos_[n].z == 0)            pos_[n].z = latticeZ_;
      }
    }
    break;
  case 5:
    pos_[n].x = pos_[n].x -1;
    pos_[n].y = pos_[n].y +1;
    pos_[n].z = pos_[n].z -1;

    if (pos_[n].x == 0 || pos_[n].y == latticeY_ +1 || pos_[n].z == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == 0)            pos_[n].x = latticeX_;
        if (pos_[n].y == latticeY_ +1) pos_[n].y = 1;
        if (pos_[n].z == 0)            pos_[n].z = latticeZ_;
      }
    }
    break;
  case 6:
    pos_[n].x = pos_[n].x -1;
    pos_[n].y = pos_[n].y -1;
    pos_[n].z = pos_[n].z +1;

    if (pos_[n].x == 0 || pos_[n].y == 0 || pos_[n].z == latticeZ_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == 0)            pos_[n].x = latticeX_;
        if (pos_[n].y == 0)            pos_[n].y = latticeY_;
        if (pos_[n].z == latticeZ_ +1) pos_[n].z = 1;
      }
    }
    break;
  case 7:
    pos_[n].x = pos_[n].x -1;
    pos_[n].y = pos_[n].y -1;
    pos_[n].z = pos_[n].z -1;

    if (pos_[n].x == 0 || pos_[n].y == 0 || pos_[n].z == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == 0) pos_[n].x = latticeX_;
        if (pos_[n].y == 0) pos_[n].y = latticeY_;
        if (pos_[n].z == 0) pos_[n].z = latticeZ_;
      }
    }
    break;
  default:
    throw std::string("move in direction" + tostring(R) + "not implemented");
  }

  if(isOutside && isBoundaryFix_)
    pos_[n] = old;

  //check that the new site is vacant if not, move the particle
  //back. (done by the vacancyCheck...)
  bool collided = vacancyCheck(n,old);

  //If n is tagged particle, then also move it in its second,
  //"infinite" coordinate system, - if it didn't collied with any
  //particle on our lattice above.
  if (n == 0 && !collided && !(isBoundaryFix_ && isOutside)){
    switch(R) {
    case 0:
      true_x++;
      true_y++;
      true_z++;
      break;
    case 1:
      true_x--;
      true_y++;
      true_z++;
      break;
    case 2:
      true_x++;
      true_y--;
      true_z++;
      break;
    case 3:
      true_x++;
      true_y++;
      true_z--;
      break;
    case 4:
      true_x++;
      true_y--;
      true_z--;
      break;
    case 5:
      true_x--;
      true_y++;
      true_z--;
      break;
    case 6:
      true_x--;
      true_y--;
      true_z++;
      break;
    case 7:
      true_x--;
      true_y--;
      true_z--;
      break;
    default:
      throw std::string("move in direction" + tostring(R) + "not implemented");
    }
  }
}

void FCC::moveAndBoundaryCheck(size_t n, size_t R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant.

  // Face centered cubic geometry:
  // n = particle to move, R = direction to move
  // boundary: boundary=0 = periodic, boundary=1 = fix boundary

  if (n > noParticles_ -1)
    throw std::string("Invalid particle label, n=" + tostring(n) +
                      ", must fulfill 0 <= n < N");

  if (R > directions_ -1)
    throw std::string("There is a maximum of " + tostring(directions_) +
                      " positive directions to move in");

  if (isTestOn_) std::cout << "Moving particle " << n << std::endl;

  // if particle is outside and fix boundary, we know we should move it back
  bool isOutside = false;

  //Store in case we want to restore the previous configuration
  Particle old = pos_[n];

  // Move particle n on the lattice:
  switch(R){
  case 0:
    pos_[n].x = pos_[n].x +1;
    pos_[n].y = pos_[n].y +1;
    //pos_[n].z = pos_[n].z +0;

    if (pos_[n].x == latticeX_ +1 || pos_[n].y == latticeY_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == latticeX_ +1) pos_[n].x = 1;
        if (pos_[n].y == latticeY_ +1) pos_[n].y = 1;
      }
    }
    break;
  case 1:
    pos_[n].x = pos_[n].x -1;
    pos_[n].y = pos_[n].y +1;
    //pos_[n].z = pos_[n].z +0;

    if (pos_[n].x == 0 || pos_[n].y == latticeY_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == 0)            pos_[n].x = latticeX_;
        if (pos_[n].y == latticeY_ +1) pos_[n].y = 1;
      }
    }
    break;
  case 2:
    pos_[n].x = pos_[n].x +1;
    pos_[n].y = pos_[n].y -1;
    //pos_[n].z = pos_[n].z +0;

    if (pos_[n].x == latticeX_ +1 || pos_[n].y == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == latticeX_ +1) pos_[n].x = 1;
        if (pos_[n].y == 0)            pos_[n].y = latticeY_;
      }
    }
    break;
  case 3:
    pos_[n].x = pos_[n].x -1;
    pos_[n].y = pos_[n].y -1;
    //pos_[n].z = pos_[n].z +0;

    if (pos_[n].x == 0 || pos_[n].y == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == 0) pos_[n].x = latticeX_;
        if (pos_[n].y == 0) pos_[n].y = latticeY_;
      }
    }
    break;
  case 4:
    pos_[n].x = pos_[n].x +1;
    //pos_[n].y = pos_[n].y +0;
    pos_[n].z = pos_[n].z +1;

    if (pos_[n].x == latticeX_ +1 || pos_[n].z == latticeZ_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == latticeX_ +1) pos_[n].x = 1;
        if (pos_[n].z == latticeZ_ +1) pos_[n].z = 1;
      }
    }
    break;
  case 5:
    pos_[n].x = pos_[n].x +1;
    //pos_[n].y = pos_[n].y +0;
    pos_[n].z = pos_[n].z -1;

    if (pos_[n].x == latticeX_ +1 || pos_[n].z == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == latticeX_ +1) pos_[n].x = 1;
        if (pos_[n].z == 0)            pos_[n].z = latticeZ_;
      }
    }
    break;
  case 6:
    pos_[n].x = pos_[n].x -1;
    //pos_[n].y = pos_[n].y +0;
    pos_[n].z = pos_[n].z +1;

    if (pos_[n].x == 0 || pos_[n].z == latticeZ_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == 0)            pos_[n].x = latticeX_;
        if (pos_[n].z == latticeZ_ +1) pos_[n].z = 1;
      }
    }
    break;
  case 7:
    pos_[n].x = pos_[n].x -1;
    //pos_[n].y = pos_[n].y +0;
    pos_[n].z = pos_[n].z -1;

    if (pos_[n].x == 0 || pos_[n].z == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].x == 0) pos_[n].x = latticeX_;
        if (pos_[n].z == 0) pos_[n].z = latticeZ_;
      }
    }
    break;
  case 8:
    //pos_[n].x = pos_[n].x +0;
    pos_[n].y = pos_[n].y +1;
    pos_[n].z = pos_[n].z +1;

    if (pos_[n].y == latticeY_ +1 || pos_[n].z == latticeZ_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].y == latticeY_ +1) pos_[n].y = 1;
        if (pos_[n].z == latticeZ_ +1) pos_[n].z = 1;
      }
    }
    break;
  case 9:
    //pos_[n].x = pos_[n].x +0;
    pos_[n].y = pos_[n].y +1;
    pos_[n].z = pos_[n].z -1;

    if (pos_[n].y == latticeY_ +1 || pos_[n].z == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].y == latticeY_ +1) pos_[n].y = 1;
        if (pos_[n].z == 0)            pos_[n].z = latticeZ_;
      }
    }
    break;
  case 10:
    //pos_[n].x = pos_[n].x +0;
    pos_[n].y = pos_[n].y -1;
    pos_[n].z = pos_[n].z +1;

    if (pos_[n].y == 0 || pos_[n].z == latticeZ_ +1){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].y == 0)            pos_[n].y = latticeY_;
        if (pos_[n].z == latticeZ_ +1) pos_[n].z = 1;
      }
    }
    break;
  case 11:
    //pos_[n].x = pos_[n].x +0;
    pos_[n].y = pos_[n].y -1;
    pos_[n].z = pos_[n].z -1;

    if (pos_[n].y == 0 || pos_[n].z == 0){
      isOutside = true;
      if(!isBoundaryFix_){
        if (pos_[n].y == 0) pos_[n].y = latticeY_;
        if (pos_[n].z == 0) pos_[n].z = latticeZ_;
      }
    }
    break;
  default:
    throw std::string("move in direction" + tostring(R) + "not implemented");
  }

  if(isOutside && isBoundaryFix_)
    pos_[n] = old;

  //check that the new site is vacant if not, move the particle
  //back. (done by the vacancyCheck...)
  bool collided = vacancyCheck(n,old);

  //If n is tagged particle, then also move it in its second,
  //"infinite" coordinate system, - if it didn't collied with any
  //particle on our lattice above.
  if (n == 0 && !collided && !(isBoundaryFix_ && isOutside)){
    switch(R) {
    case 0:
      true_x++;
      true_y++;
      break;
    case 1:
      true_x--;
      true_y++;
      break;
    case 2:
      true_x++;
      true_y--;
      break;
    case 3:
      true_x--;
      true_y--;
      break;
    case 4:
      true_x++;
      true_z++;
      break;
    case 5:
      true_x++;
      true_z--;
      break;
    case 6:
      true_x--;
      true_z++;
      break;
    case 7:
      true_x--;
      true_z--;
      break;
    case 8:
      true_y++;
      true_z++;
      break;
    case 9:
      true_y++;
      true_z--;
      break;
    case 10:
      true_y--;
      true_z++;
      break;
    case 11:
      true_y--;
      true_z--;
      break;
    default:
      throw std::string("move in direction" + tostring(R) + "not implemented");
    }
  }
}


void Special::moveAndBoundaryCheck(size_t n, size_t R){
  if (n > noParticles_ -1)
    throw std::string("Invalid particle label, n=" + tostring(n) +
                      ", must fulfill 0 <= n < N");

  if (R > directions_ -1)
    throw std::string("There is a maximum of " + tostring(directions_) +
                      " positive directions to move in");

  if (isTestOn_) std::cout << "Moving particle " << n << std::endl;

  // if particle is outside and fix boundary, we know we should move it back
  bool isOutside = false;

  //Store in case we want to restore the previous configuration
  Particle old = pos_[n];

    // Move particle n on the lattice:
  switch(R){
  case 0:
    pos_[n].x = pos_[n].x +1;

    if(pos_[n].x == latticeX_ +1){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].x = 1;
    }
    break;
  case 1:
    pos_[n].x = pos_[n].x -1;

    if(pos_[n].x == 0){
      isOutside = true;
      if(!isBoundaryFix_) pos_[n].x = latticeX_;
    }
    break;
  default:
    throw std::string("move in direction" + tostring(R) + "not implemented");
  }

  if(isOutside && isBoundaryFix_)
    pos_[n] = old;

  //check that the new site is vacant if not, move the particle
  //back. (done by the vacancyCheck...)
  bool collided = vacancyCheck(n,old);

  //If n is tagged particle, then also move it in its second,
  //"infinite" coordinate system, - if it didn't collied with any
  //particle on our lattice above.
  if (n == 0 && !collided && !(isBoundaryFix_ && isOutside)){
    switch(R) {
    case 0:
      true_x++;
      break;
    case 1:
      true_x--;
      break;
    default:
      throw std::string("move in direction" + tostring(R) + "not implemented");
    }
  }
}
