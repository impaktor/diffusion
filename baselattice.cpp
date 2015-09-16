#include <iostream>        //for cout, among others
#include <cstdlib>         //for abort()-function
#include <cmath>           //gives sqrt, pow, fabs,
#include <functional>      // for std::function
#include <vector>

#include "auxiliary.h"     //for tostring() function
#include "classes.h"       //various data structures /classes (Jump, Particle)
#include "baselattice.h"   //this is the main class, that does the physics.

using std::vector;

//CPU consuming! Use only when bug-squashing!
const bool CHECK_VACANCY_ON  = false;

//boolean, prints simulation details to screen,
const bool TEST_1            = false;


// =====================
// C O N S T R U C T O R
// =====================

BaseLattice::BaseLattice(int x, int y, int z, int dimension,
                         double seed, bool boundaryFix)
  //initiate the random number generator with "seed" just given above
  : randomNumber(seed), board_(x, y, z){

  dim_ = dimension;
  directions_ = dim_ * 2;

  //number of lattice sites
  latticeX_ = x;
  latticeY_ = dim_ >= 2 ? y : 1;
  latticeZ_ = dim_ == 3 ? z : 1;

  isBoundaryFix_ = boundaryFix;

  center_.x = round(latticeX_/2.0);
  center_.y = round(latticeY_/2.0);
  center_.z = round(latticeZ_/2.0);

  //Hardcoded Switches:
  //------------------
  //boolean, for print details to screen.
  isTestOn_ = TEST_1;

  //By default: don't use the interaction algorithm.
  interactionStrength_ = 0;
  isInteractionOn_ = false;
}


void BaseLattice::place(void){
  //This function will be called once each ensemble to
  //place the N particles (randomly), and reset the time

  timeSum_ = 0;                   //Time starts at zero each run

  if(pos_.size() < 1)
    throw std::string("No particles to place");

  //reset this counter:
  windingNumber_0_x = 0;
  windingNumber_0_y = 0;
  windingNumber_0_z = 0;

  int totalSites = latticeX_ * latticeY_ * latticeZ_;
  int leftToPlace = noParticles_;

  //place tagged particle at lattice center
  pos_[0].setPos(center_.x, center_.y, center_.z);
  leftToPlace--;                  //one less to place
  totalSites--;

  double R;                               //store random number here
  int i, j, k, n = 1;                     //First crowder at n = 1 element
  for (k=1; k <= latticeZ_; k++){         //place particles
    for (j=1; j <= latticeY_; j++){
      for (i=1; i <= latticeX_; i++){

        //don't place crowders on the tagged particle
        if (i != pos_[0].x || j != pos_[0].y || k != pos_[0].z){

          //using NR to generate 0 < R < 1 ;
          R = randomNumber.doub();
          if (R <= leftToPlace*1.0 / totalSites){
            pos_[n].x = i;
            pos_[n].y = j;
            pos_[n].z = k;
            n++;

            leftToPlace--;
          }

          totalSites--;
        }
      }
    }
  }

  //Also use a matrix representation of the lattice, for speed (fast
  //lookup).
  board_.initialize(pos_);
}


// ---------------------
// PUBLIC: SET FUNCTIONS:
// ---------------------

void BaseLattice::setSamplingTimes(const vector<double>& samplingTime,
                               waitingtime tagged_waiting_time = waitingtime::LIN){

  tagged_waiting_time_ = tagged_waiting_time;

  noSamplingTimes_ = (int) samplingTime.size();

  if (noSamplingTimes_ < 1)
    throw std::string("No sampling times set!");

  //check that it is monotonously increasing:
  bool isOrdered = true;

  if (samplingTime[0] < 0)
    isOrdered = false;

  for (int i = 1; i < noSamplingTimes_; i++){
    if (samplingTime[i-1] > samplingTime[i]){
      isOrdered = false;
      break;
    }
  }

  if (isOrdered){
    samplingTime_ = samplingTime;

    //assign the displacement container. I only need to do this
    // once/program run, as the values are then (re)set in move().
    dx_.assign(noSamplingTimes_,0);
    dy_.assign(noSamplingTimes_,0);
    dz_.assign(noSamplingTimes_,0);
    dr_.assign(noSamplingTimes_,0);
  }
  else
    throw std::string("Sampling times are not ordered");
}


void BaseLattice::setJumpNaka(const float jumpCrowders,
                              const float jumpTracer){
  jumpCrowders_ = jumpCrowders;
  jumpTracer_ = jumpTracer;
}

void BaseLattice::setParticles(const vector<Particle>& particles){
  if(particles.size() == 0)
    throw std::string("received vector of 0 particles");

  noParticles_ = particles.size();
  pos_ = particles;

  //Use these vectors to build the partial sum
  computePartialSum();
}

void BaseLattice::setInteraction(const float interactionStrength){
  interactionStrength_ = interactionStrength;
  isInteractionOn_ = true;
}


// ---------------------
// PUBLIC: GET FUNCTIONS:
// ---------------------

//could be useful to have, but not necessary...
unsigned int BaseLattice::getDimension(void) const{
  return dim_;
}

unsigned int BaseLattice::getDirections(void) const{
  return directions_;
}

void BaseLattice::getDisplacement(vector<int>& dx, vector<int>& dy,
                              vector<int>& dz, vector<double>& dr) const{
  //Get the displacement of the tracer particle
  //from the current "ensemble"
  dx = dx_;
  dy = dy_;
  dz = dz_;
  dr = dr_;
}


//check if site is vacant or not
int BaseLattice::vacancyCheck(int n, const Particle& oldPos){
  int returnvalue = 0;   //unused parameter to register collisions

  if (n > noParticles_ || 0 > n)
    throw std::string("accessing invalid particle");

  if (board_.isOccupied(pos_[n])){ //if occupied
    if (isTestOn_)
      std::cout << "\033[91mOccupied:\033[0m " << pos_[n] << "\n";
    pos_[n] = oldPos;                                 //then move back
    returnvalue = 1;
  }
  else{
    //if no interaction att all:

    //if vacant, update the board
    if(!isInteractionOn_){
       board_.setAsVacant(oldPos);         //old position is now vacant
       board_.setAsOccupied(pos_[n], n);   //new position is now occupied

       if (isTestOn_)
         std::cout << "\033[92mVacant\033[0m, " << n <<": " << pos_[n] << "\n";
    }
    //if we use the interaction algorithm:
    else
      interaction(n, oldPos);
  }

  if(isTestOn_)
    std::cout << "now particle: "<< n << " is at pos: "
              << pos_[n] << "\n" << std::endl;

  return returnvalue;
  //(the return value could be used to have a
  //    do-move-while(vacancyCheck == 1)
  //but if so, we need to move this function to be directly in the
  //Move_()-func.)
}


double BaseLattice::computeWaitingTime(void){

  //one time step for crowding particle
  double tau = 0;

  double r;

  if(tagged_waiting_time_ == waitingtime::EXP){
    do{
      r = randomNumber.doub();
    }while( r == 1 || r == 0 );

    tau = log(1.0 / r) / partialSum_.back();
  }
  else if(tagged_waiting_time_ == waitingtime::POW){

    do{
      r = randomNumber.doub();
    }while( r == 1 || r == 0 );

    // if 0 < alpha < 2: second mom inf. if alpha < 1, first mom inf.
    const double alpha = 0.5;
    const double a = 1;
    tau = a * std::pow(r, -1.0/alpha);
    //tau = a * (std::pow(r, -1.0/alpha) - 1.0); // correct version, but identical to above
  }
  else if(tagged_waiting_time_ == waitingtime::LIN)
    tau = 1.0 / partialSum_.back();

  return tau;
}


void BaseLattice::move(){

  double tau;      //waiting time
  int i = 0;       //index of samplingTime-vector
  timeSum_ = 0;

  while(timeSum_ < samplingTime_[noSamplingTimes_-1]){ // "-1" since we start on 0.
    //The time-sampling was completely rewritten in late August 2010
    //to resolve the spaghetti that was the previous version. This
    //follows the gillespie_exclusion2.cpp-implementation closely, to
    //make it easier on those already familiar with that code.

    tau = computeWaitingTime();

    //Save displacement if next time-step is beyond next sampling time
    while(timeSum_ <= samplingTime_[i] &&
          samplingTime_[i] < timeSum_ + tau  && i < noSamplingTimes_){

      //save displacement (from previous step), depends on geometry of lattice
      dr_[i] = distance(dx_[i], dy_[i], dz_[i]);

      i++;
    }

    double r2;          //stores random number

    int mu_guess;       // must be integer. (index of vector)
    int mu_left = 0;
    int mu_right = directions_ * noParticles_;
    double p_left = (double) partialSum_.front(); //value in first element (=0)
    double p_right = (double) partialSum_.back(); //value in last element (=dim_*2*N)

    do{
      r2 = randomNumber.doub();
    }while( r2 == 1 || r2 == 0 );

    double p_rand = (double) r2 * partialSum_.back();
    bool loopAgain = true;

    //NOTE: Every 10^15 run of move() will lead to a situation where
    //we actually have p_rand==p[mu_guess], and changing the interval
    //below in the if statements, from "..<= .. <.." to "..< .. <=.."
    //will result in an infinite loop every 10^15 turns since "double"
    //stores 52 bit, 2^52 ~= 10^15.  This bug has now been fixed!

    //finds the mu
    do{
      mu_guess = (int) ( (p_rand - p_left) * ( mu_right - mu_left )
                         / (p_right - p_left) + mu_left );

      if ((double) partialSum_[mu_guess] <= p_rand && p_rand <
          (double) partialSum_[mu_guess+1])
        loopAgain = false;

      else{
        if (p_rand < (double) partialSum_[mu_guess]){
          mu_right = mu_guess;
          p_right = partialSum_[mu_guess];
        }
        else{//if (r2*partialSum_.back() > partialSum_[mu_guess+1]){
          mu_left = mu_guess + 1;
          p_left = partialSum_[mu_guess+1];
        }
      }
    }while(loopAgain);

    int mu = mu_guess;
    int n,r;

    //transforms the index mu to which particle to move (n) and in
    //which direction (r, 0<= r <= 5 (in 3D))
    convertMuToParticle(mu,n,r);

    //moves the particle n in direction r, according to boundary
    //conditions. If new site is occupied, move back to previous pos.
    moveAndBoundaryCheck(n,r);

    timeSum_ = timeSum_ + tau;
  }
}

void BaseLattice::computePartialSum(){
  long double tmp = 0;

  // delete previous run
  partialSum_.clear();

  // first element must be zero
  partialSum_.push_back(tmp);

  if ((int) pos_.size() != noParticles_)
    throw std::string("number of particles, != N");

  //build partialSum-vector
  for(size_t j = 0; j < pos_[0].jumprates.size(); ++j){
    for(size_t i = 0; i < pos_.size(); ++i){
      tmp += pos_[i].jumprates[j];
      partialSum_.push_back(tmp);
    }
  }

  //check for errors: It must increase monotonically:
  for (size_t i = 1; i < pos_.size(); i++)
    if (partialSum_[i] < partialSum_[i-1])
      throw std::string("jump rate cumul. sum non-monotonically increasing");
}

// ----------------------
// JUMPRATE SPECIFIC CODE
// ----------------------



void BaseLattice::convertMuToParticle(const int mu, int& n, int& direction){

  //if no error or inconsistency...
  if (0 <= mu && mu <= (int) partialSum_.size() -1){
    //... then start checking the intervals:

    if ( 0 <= mu && mu <= 2 * noParticles_ - 1 ){
      if ( mu <= noParticles_ - 1 ){
        n = mu;
        direction = 0;
      }
      //if ( noParticles_ <= mu){
      else{
        n = mu - noParticles_;
        direction = 1;
      }
    }
    else{
      if ( 2 * noParticles_ <= mu && mu <= 4 * noParticles_ - 1){
        if ( mu <= 3 * noParticles_-1){
          n = mu - 2 * noParticles_;
          direction = 2;
        }
        //if ( 3 * noParticles_ <= mu ){
        else{
          n = mu - 3 * noParticles_;
          direction = 3;
        }
      }//if ( 4 * noParticles_ <= mu && mu <= 6 * noParticles_ - 1){
      else{
        if ( mu <= 5 * noParticles_ - 1){
          n = mu -4 * noParticles_;
          direction = 4;
        }
        //if ( 5 * noParticles_ <= mu ){
        else{
          n = mu -5 * noParticles_;
          direction = 5;
        }
      }
    }
  }
  else
    throw std::string("mu=" + tostring(mu) + "or partialSum_=" +
                      tostring(partialSum_.size()) + "not correct");
}


//called from main, to print out fun/useful fact to head of data-file.
float BaseLattice::computeEffectiveDiffusionConst(void){
  float D_eff;

  if(!partialSum_.empty()){

    //NOTE: if the lattice constant: a != 1, it must be included
    //D_i=k_i * pow(a,2)
    double sum_temp = 0;

    for(size_t j = 0; j < pos_[0].jumprates.size(); ++j){
      for(size_t i=0; i < pos_.size(); i++){
        sum_temp += 1.0/pos_[i].jumprates[j];
      }
    }

    sum_temp = (double) sum_temp / (2*noParticles_*dim_);
    D_eff = 1.0 / sum_temp;

    return D_eff;
  }
  else{
    //When we don't use jumpratios we get the diffusion constant to:
    //TODO: stämmer detta? den skall ju bli <r^2>=2*D*dim_*t ! -->
    //svar: med denna D_eff får vi f(x)=2*dim_*D_eff*x, vilket är
    //inkonsekvent med ovanstående Deff som redan är kompenserad för
    //det... tror jag... TSTAR NU: N=1, uniform distrubution...
    //RESULTAT: nej, även med ovanst. D_eff utr. måste vi ha
    //f(t)=2*dim_*D_eff*t Detta är det rimmliga, eftersom vid N>1 vill
    //vi ha den "riktiga" D_eff, dvs. som det är nu.  Note, assumed
    //lattice constant a = 1.
    D_eff = 1.0 / (2 * dim_);
    return D_eff;
  }
}


//called from main, to print out fun/useful fact to head of data-file.
float BaseLattice::computeAverageDiffusionConst(void){
  //Calculate Average:

  float D_av = 0;

  if(!partialSum_.empty()){

    int startingParticle = 1;  //don't include the tagged
    for(size_t j = 0; j < pos_[0].jumprates.size(); ++j){
      for(size_t i = startingParticle; i < pos_.size(); i++){
        D_av += pos_[i].jumprates[j];
      }
    }
    D_av = D_av * 1.0/(noParticles_ * 2 * dim_);

    return D_av;
  }
  else{
    //Note, assumed lattice constant a = 1.
    D_av  = 1.0 / (2*dim_);
    return D_av;
  }
}

//called from main, to print out fun/useful fact to head of data-file.
float BaseLattice::computeNakazato(void){

  float alpha;

  if (dim_ == 1) alpha = 1;
  if (dim_ == 2) alpha = 0.363;
  if (dim_ == 3) alpha = 0.209;

  //jump-ratio for the different particles...
  double gamma = jumpTracer_ / jumpCrowders_;

  float c = (float) noParticles_ / (latticeX_ * latticeY_ * latticeZ_);

  double nominator   = (gamma * (1-c)+1) * (1-alpha);
  double denominator = gamma*(1-c) + 1-alpha*(1 + gamma*(1-3*c));
  double f = nominator / denominator;

  //lattice constant
  float a = 1;

  float D = jumpTracer_ * (1-c) * f * pow(a,2);

  //due to laziness, we multiply the "true" diffusion constant here
  //with 2*dimension, and return that. Keep that in mind when plotting.
  D = D * dim_ * 2;

  return D;
}


//called from main, to print out fun/useful fact to head of data-file.
float BaseLattice::computeErgodicity(int squares){
  //check what the MSD should be when the system reaches equilibrium.
  //(only correct if N=1). This function has two independent parts:
  //either do it numerically or analytically.

  float sum = 0;

  //use these to "walk" through the lattice to measure the MSD
  //contribution from each square
  int x, y, z;

  int x_size = squares;
  int y_size = dim_ >= 2 ? squares : 1;
  int z_size = dim_ == 3 ? squares : 1;

  for (int i=0; i < x_size; i++)
    for (int j=0; j < y_size; j++)
      for (int k=0; k < z_size; k++){
        x = i + 1;               //First square has coord.(1,1,1)
        y = j + 1;
        z = k + 1;
        sum += pow(center_.x-x,2)+pow(center_.y-y,2)+pow(center_.z-z,2);
      }
  sum = (float) sum / (x_size*y_size*z_size);


  //ANALYTIC EXPRESSION: this is an alternate way of calculating
  //<R^2>, when the system becomes ergodic.
  bool testAnalytic = false;
  if (testAnalytic){
    int a = center_.x, b = center_.y, c = center_.z;
    float xContrib = ((x_size + 1) * (2 * x_size + 1) *
                           1.0/6 - a * (x_size+1) + a * a);
    float yContrib = ((y_size + 1) * (2 * y_size + 1) *
                           1.0/6 - b * (y_size+1) + b * b);
    float zContrib = ((z_size + 1) * (2 * z_size + 1) *
                           1.0/6 - c * (z_size + 1) + c * c);
    std::cout << "Analytic result for <R^2> (ergodicity): "
              << xContrib + yContrib + zContrib << std::endl;
  }

  return sum;
}



// -------------------------
// INTERACTION SPECIFIC CODE
// -------------------------

void BaseLattice::interaction(int n, Particle oldPos){
  //This is an updated version of the progInteraction.cpp I used in my
  //thesis It does the exact same thing, but now we can wrap around
  //the boundary also...
  //The output-files are identical between the two different version,
  //when I use fixed boundary conditions, which was what I used in my
  //master thesis.

  vector<int> NearestNeighbours;
  board_.checkNeighbours(oldPos, NearestNeighbours, isBoundaryFix_);

  int counter = NearestNeighbours.size();

  double r0;

  static Ran slump(56); // NOTE: THIS IS NOT BE THREAD SAFE!

#ifdef _OPENMP //If compiling with OPenMP:
  std::cout << "WARNING: INTERACTION CODE + THREADS?" << std::endl;
#endif

  do{
    r0=slump.doub();
  }while( r0 == 1 || r0 == 0 );

  float DeltaV = interactionStrength_;
  if (r0 < exp(-1.0*counter*DeltaV) ){        //accept move:
    board_.setAsVacant(oldPos);                //old pos. now vacant
    board_.setAsOccupied(pos_[n], n);          //new pos. now occupied
  }
  else{                                        //reject move, move back:
    pos_[n] = oldPos;                          //board_ unchanged
  }
}
