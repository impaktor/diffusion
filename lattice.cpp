#include <iostream>        //for cout, among others
#include <fstream>         //read write to files
#include <cstdlib>         //for abort()-function
#include <cmath>           //gives sqrt, pow, fabs,
#include <sstream>         //behövs kanske för "ostringstream" i save()-func.
#include <string>
#include <vector>

#include "nr/nr3.h"        //ran.h, ludcmp.h, and svd.h depends on this one.
#include "nr/ran.h"        //Random number generator from numerical recipes 3ed.

#include "classes.h"       //various data structures /classes (Jump, Particle)
#include "lattice.h"       //this is the main class, that does the physics.

using namespace std;

//HARDCODED VARIABLES:
//-------------------

//boolean, prints details to screen
const bool TEST_1            = 0;

//CPU consuming, use only when bug-squashing!
const bool CHECK_VACANCY_ON  = 0;

// =====================
// C O N S T R U C T O R
// =====================

Lattice::Lattice(int xSquare,int ySquare,int zSquare,
                 int particleNumber, double seed, bool boundary)
  //initiate the random number generator with "seed" just given
  : randomNumber(seed){

  if (xSquare >= ySquare && ySquare  >= zSquare && zSquare >= 1){
    if (particleNumber <= xSquare*ySquare*zSquare){

      //number of lattice sites
      latticeX_ = xSquare;
      latticeY_ = ySquare;
      latticeZ_ = zSquare;

      nParticles_ = particleNumber;

      boundaryFix_ = boundary;

      //place tagged particle at lattice center
      pos_0_.x = round(1.0*latticeX_/2);
      pos_0_.y = round(1.0*latticeY_/2);
      pos_0_.z = round(1.0*latticeZ_/2);

      //Hardcoded Switches:
      //------------------
      //boolean, for print details to screen.
      testOnOff_ = TEST_1;

      //Define the dimension of the lattice
      if (latticeY_<=1 && latticeZ_<=1) dim_ = 1;
      if (latticeY_>1 && latticeZ_<=1)  dim_ = 2;
      if (latticeY_>1 && latticeZ_>1)   dim_ = 3;

      //By default: don't use the interaction algorithm.
      interactionStrength_ = 0;
      interactionOn_ = false;

      //--------------------------
      //To get the Cluster Distribution (see saveCluster()):
      // Cluster_Distribution_ = number of clusters with "m" particles
      // Cluster_Size_ = number of particles in cluster (just our x-axis)
      int Maximus_Clustirus = nParticles_;
      //Maximal cluster size to include in distribution
      clusterDistribution_.assign(Maximus_Clustirus,0); //assign zeros...
      for (int i=1; i <= Maximus_Clustirus ; i++){
        clusterSize_.push_back(i);
      }//--------------------------

    }
    else {
      printError("Number of particles must be <= number of squares");
    }
  }
  else{
    printError("X,Y,Z must be 1 or bigger, and X >= Y >= Z");
  }
}


void Lattice::place(void){
  //This function will be called once each ensemble to
  //place the N particles (randomly), and reset the time

  int totalSites = latticeX_ * latticeY_ * latticeZ_;
  int leftToPlace = nParticles_;

  //assign zeros
  Particle nil(0,0,0);
  pos_.assign(nParticles_,nil);

  timeSum_ = 0;                   //Time starts at zero each run

  pos_[0] = pos_0_;               //place tracer particle at center
  leftToPlace--;                  //one less to place
  totalSites--;

  double R;                               //store random number here
  int i, j, k, n = 1;                     //First crowder at n = 1 element
  for (k=1; k <= latticeZ_; k++){         //place particles
    for (j=1; j <= latticeY_; j++){
      for (i=1; i <= latticeX_; i++){

        //don't place crowders on the tagged particle
        if (i != pos_0_.x || j != pos_0_.y || k != pos_0_.z){

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
  //TEST:
  if (testOnOff_)
    for (int q = 0; q < nParticles_; q++)
      cout <<" PLACED \n particle n = " << q << pos_[q] << endl;


  //IF we use the new vacancyCheck()-function:
  //-----------------------------------
  //initiate the vacancy-vector with "-1" (meaning vacant)
  //Must reset this matrix for each ensemble.
  vector< vector< vector<int> > >
    Vacancy(latticeX_+1, vector< vector<int> >
            (latticeY_+1, vector<int>(latticeZ_+1,-1)));
  vacancy_ = Vacancy;
  //Mark occupied sites with particle label for that site:
  //note, no change will ever occur to the 0-elements of the matrix,
  //since out lattice starts at 1,..,latticeX_, but matrix is 0,..,latticeX_
  for (n=0; n < nParticles_; n++){
    vacancy_[ pos_[n].x ][ pos_[n].y ][ pos_[n].z ] = n;
  }
}


void Lattice::dumpSimulation(int index){
  //This is a test method to simply "dump" the trajectories
  //for a single run to file. It is used to check correlations
  //between different runs.

  //save to file-name dump[i].dat
  stringstream buffer;
  buffer << "dump" << index << ".dat";

  ofstream dump;
  dump.open(buffer.str().c_str());

  for (int i = 0; i < maxElement_; i++){
    dump  << samplingTime_[i] << "\t" << dr_[i] << endl;
  }
  dump.close();
}


void Lattice::saveCluster(double Time){
  //This function calculates how many particles are in
  //the cluster that the tagged particle is part of

  //ha denna funktionen sa att den kors varje tidssteg for sista ensemblen,
  //Den skall veta nar den ar pa sista tidssteget, tror jag...

  //Cluster_Distribution defined in the constructor.
  vector<int> Cluster;
  //TODO kolla Boundary! 0 eller 1?

  //BUG augusti, detta stammer val inte?
  //ne, men jag har kommenterat ut hela den saken, sa vardet pa
  //Exponent spellar ingen roll!
  double Exponent = 0; //ie. rand > 0 always include the neighbor.

  // checkVacancyMatrix(1.01); //OK hit

  if(Time>20){  //The system is not in equilibrium from start
    //checkVacancyMatrix(10.1);
    countNeighbours(pos_[0], Cluster);
    //checkVacancyMatrix(10.2);
    buildCluster2(0, Cluster, Exponent); //Skriver denna till nagot?
    //checkVacancyMatrix(10.3); //OK hit
    int M = Cluster.size();
    //First element is # of 1-clusters (ie. single particles), etc.
    clusterDistribution_[M-1] = clusterDistribution_[M-1] + 1.0/nParticles_;
    // CheckVacancyMatrix(10.4);  //BUG! men inte OK hit!????

  }
  //CheckVacancyMatrix(10.5);

  //TODO
  //-normalisering? dela med N*Ensemble?
  //integrerar man den skall den bli 1 !
  //hur ofta skall man kolla clusterdistribution?
  //-knovergens. I borjan skall det vara en forskjutning mot 1
  //skall sedan stabilisera sig.
}


void Lattice::printCluster(void){
  //This function prints the information gathered by
  //the saveCluster()-function, and prints it to file

  //  CheckVacancyMatrix(9);

  if (clusterSize_.size() != clusterDistribution_.size())
    printError("cluster size incorrect, in printCluster()");

  ofstream clust;
  char NameOfFiles[] = "cluster.dat";
  clust.open(NameOfFiles);
  for(unsigned int j = 0; j < clusterDistribution_.size(); j++){
    clust << clusterSize_[j] <<"\t"<< clusterDistribution_[j] << endl;
  }
  clust.close();
  cout << "Pinted " << NameOfFiles << endl;
}


void Lattice::checkVacancyMatrix(float ErrorCode){
  //This is just a small function to check that the vacancy matrix is correct
  //ie. vacant sites are marked as "-1" and occupied are marked with the
  //correct particle label.
  //ErrorCode = code to output to easier find the place in the code

  if(CHECK_VACANCY_ON){

    bool error = false;
    int label;
    int count_vacant_sites = 0;

    for (int i=1; i <= latticeX_; i++){
      //cout<<"i: "<<i<<endl;
      for (int j=1; j <= latticeY_; j++){
        //cout<<"j: "<<j<<endl;
        for (int k=1; k <= latticeZ_; k++){
          //cout<<"k: "<<k<<endl;
          label = vacancy_[i][j][k];
          //cout <<"label: "<< label <<endl;
          if(label != -1){
            if(i != pos_[label].x || j != pos_[label].y || k != pos_[label].z){
              error = true;
            }
          }
          else
            count_vacant_sites++;
        }
      }
    }
    //Easier to find exactly where (above/below) error = true
    if(error) cout << "Does not compute..." << endl;

    if (count_vacant_sites != latticeX_*latticeY_*latticeZ_-nParticles_){
      cout << "Error in vacant sites" << endl
           << "# of vacant sites: " << count_vacant_sites << endl;
      error = true;
    }

    if(error){
      cout << "Vacancy Matrix is not correct, at err. code site: "
           << ErrorCode << endl;
      abort();
    }
  }// end of if(CHECK_VACANCY_ON)
}


// double Lattice::get_dt(void){
//   //Needed for the newly added binning-function
//   //to produce a histogram of P(x,t)
//   //(and for the snapshot-function, see "main()")
//   return partialSum_.back();
// }


// ---------------------
// PUBLIC: SET FUNCTIONS:
// ---------------------

void Lattice::setSamplingTimes(const vector<double>& samplingTime,
                               bool isExpWaitingTime = false){ //XXX set this in header I think?

  isExponentialWaitingTime_ = isExpWaitingTime;

  maxElement_ = (int) samplingTime.size();

  if (maxElement_ < 1)
    printError("No sampling times set!");

  //check that is is monotonously increasing:
  bool isOrdered = false;
  for (int i = 1; i < maxElement_; i++){
    if (samplingTime[i-1] < samplingTime[i])
      isOrdered = true;
  }

  if (samplingTime[0] < 0)
    isOrdered = false;


  if (isOrdered){
    samplingTime_ = samplingTime;

    //assign the displacement container. I only need to do this
    // once/program run, as the values are then (re)set in move().
    dx_.assign(maxElement_,0);
    dy_.assign(maxElement_,0);
    dz_.assign(maxElement_,0);
    dr_.assign(maxElement_,0);

  }
  else
    printError("Sampling times are not ordered.");
}


void Lattice::setJumpNaka(float jumpCrowders, float jumpTracer){
  jumpCrowders_ = jumpCrowders;
  jumpTracer_ = jumpTracer;
}

void Lattice::setInteraction(float interactionStrength){
  interactionStrength_ = interactionStrength;
  interactionOn_ = true;
}


void Lattice::generateTrajectory(int& E){
  //This is just a wrapper for easier passing to the thread-library.

  place();
  move();
  //getDisplacement();

  //increase the iterator that counts ensembles:
  ++E;

}

// ---------------------
// PUBLIC: GET FUNCTIONS:
// ---------------------

//could be useful to have, but not necessary...
int Lattice::getDimension(void){
  return dim_;
}

void Lattice::getDisplacement(vector<int>& dx, vector<int>& dy,
                              vector<int>& dz, vector<double>& dr){
  //Get the displacement of the tracer particle
  //from the current "ensemble"
  dx = dx_;
  dy = dy_;
  dz = dz_;
  dr = dr_;
}


void Lattice::moveOld(){
  //This code is obsolete, as it just works for placing identical
  //particles, but I keep it for sentimental reasons.
  //Better safe than sorry.

  int r,n;

  //Choose particle at random, 0 <= n < N
  do{
    n = randomNumber.doub() * nParticles_;
  }while(n == nParticles_);


  do{
    r = randomNumber.doub() * (2 * dim_);
  }while( r == dim_ * 2 );           //gives 0<=R<2*dim


  //Do the actual move, according to the boundary rules
  moveAndBoundaryCheck(n,r);

  timeSum_ = timeSum_ + 1.0 / nParticles_;  //Calculate the time

  if (testOnOff_)
    for (int i = 0; i < nParticles_; i++)
      cout <<" MOVED \n particle n = " << i <<": "<< pos_[i] << endl;
}


void Lattice::moveAndBoundaryCheck(int n, int R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant.
  bool boundary = boundaryFix_;

  if (0 <= n && n < nParticles_ ){
    if (testOnOff_) cout << "Moving particle " << n << endl;
    // n = particle to move, R = direction to move, such that:
    //      R= 0 = X-right,      1 = X-left
    //         2 = Y-right,      3 = Y-left
    //         4 = Z-up,         5 = Z-down
    // boundary: boundary=0 = periodic, boundary=1 = fix boundary

    //Store in case we want to restore the previous configuration
    Particle old = pos_[n];

    if ( R >= 0 && R < 6 ){
      switch(R){
      case 0: pos_[n].x=pos_[n].x +1;
        if (pos_[n].x == latticeX_ +1 && boundaryFix_ == 0) pos_[n].x = 1;
        if (pos_[n].x == latticeX_ +1 && boundaryFix_ == 1) pos_[n].x = old.x;
        break;
      case 1: pos_[n].x = pos_[n].x -1; //move particle n to the left
        if (pos_[n].x == 0 && boundaryFix_ == 0) pos_[n].x = latticeX_;//periodic
        if (pos_[n].x == 0 && boundaryFix_ == 1) pos_[n].x = old.x;    //fix
        break;                          //(first coordinate is 1, 0=outside)
      case 2: pos_[n].y = pos_[n].y +1;
        if (pos_[n].y == latticeY_ +1 && boundaryFix_ == 0) pos_[n].y = 1;
        if (pos_[n].y == latticeY_ +1 && boundaryFix_ == 1) pos_[n].y = old.y;
        break;
      case 3: pos_[n].y = pos_[n].y -1;
        if (pos_[n].y == 0 && boundaryFix_ == 0) pos_[n].y = latticeY_;
        if (pos_[n].y == 0 && boundaryFix_ == 1) pos_[n].y = old.y;
        break;
      case 4: pos_[n].z = pos_[n].z +1;
        if (pos_[n].z == latticeZ_+1 && boundaryFix_ == 0) pos_[n].z = 1;
        if (pos_[n].z == latticeZ_+1 && boundaryFix_ == 1) pos_[n].z = old.z;
        break;
      case 5: pos_[n].z = pos_[n].z-1;
        if (pos_[n].z == 0 && boundaryFix_ == 0) pos_[n].z = latticeZ_;
        if (pos_[n].z == 0 && boundaryFix_ == 1) pos_[n].z = old.z;
        break;
      }

      //check that the new site is vacant if not,
      //move the particle back. (done by the vacancyCheck...)

      //vacancyCheckOld(n,temp);  //Old code (slow)
      vacancyCheck(n,old);        //New improved, faster version.

    }
    else{
      cout << "There is a maximum of " << 2 * dim_
           << " positive directions to move in" << endl;
      abort();
    }
  }
  else{
    cout << "Invalid particle label, n=" << n
         << ", must fulfill 0 <= n < N" << endl;
    abort();
  }
}



int Lattice::vacancyCheckOld(int n, Particle previous){
  //Old and obsolete version. It doesn't need a 3D matrix as the
  //new version does, vacancyCheck(...) but this one is slow.
  //Keep it here to check that the output is exactly identical
  //(compare output-files with the "diff" command)

  if (n < nParticles_ && 0 <= n){

    int returnvalue = 0;    //unused parameter to register collisions
    for (int i=0; i<nParticles_; i++){
      if (pos_[i].x == pos_[n].x && pos_[i].y == pos_[n].y &&
          pos_[i].z == pos_[n].z && i!=n){
        if (testOnOff_) cout << "Occupied: " << pos_[n] << endl;

        //if two particles on same square: undo move
        pos_[n] = previous;
        returnvalue = 1;
        return returnvalue;     //End function here
      }
      else
        if (testOnOff_) cout <<"OK"<<endl;
    }
    if (testOnOff_) cout << "Vacant, particle "<< n
                         << ": " << pos_[n] << endl;
    return returnvalue;
    //(the returnvalue could be used to have a
    //    do-move-while(vacancyCheck == 1)
    //but if so, we need to move this function to be directly
    //in the Move_()-func.)
  }
  else{
    printError("accessing invalid particle");
  }
}



int Lattice::vacancyCheck(int n, const Particle& oldPos){
  //NOTE: This is a faster version than vacancyCheckOld(), since that
  //one loops through all particles for each "vacancy check". However
  //this uses a 3D matrix, and there will be some "double banking"
  //in the sense that the position is now stored in this matrix and
  //in XYZ-vectors. (use reference ('&') for speed)

  if (n < nParticles_ && 0 <= n){

    int returnvalue = 0;   //unused parameter to register collisions
    if (vacancy_[pos_[n].x][pos_[n].y][pos_[n].z] != -1){ //if occupied
      if (testOnOff_) cout << "Occupied: " << pos_[n] << endl;
      pos_[n] = oldPos;                                 //then move back
      returnvalue = 1;

      if (testOnOff_){
        for (int k = 1; k <= latticeZ_; k++){
          //print y in reverse direction, (better in the terminal)
          for (int j = latticeY_; j > 0; j--){
            cout << j << ": ";
            for (int i = 1; i <= latticeX_; i++){
              cout << vacancy_[i][j][k] << " ";
            }
            cout << endl;
          }
          cout << "looping for new Z" << endl;
        }
      }
    }
    else{
      //if no interaction att all:


      //if vacant, update vacancy_-matrix
      if(!interactionOn_){
        //old position is now vacant
        vacancy_[oldPos.x][oldPos.y][oldPos.z] = -1;

        //new position is now occupied:
        vacancy_[pos_[n].x][pos_[n].y][pos_[n].z] = n;

        if (testOnOff_) cout << "OK, vacant" << endl;
      }
      //if we use the interaction algorithm:
      else
        interaction(n,oldPos);
    }
    if (testOnOff_) cout << "Vacant, particle " << n <<": "
                         << pos_[n] << endl;
    return returnvalue;
    //(the return value could be used to have a
    //    do-move-while(vacancyCheck == 1)
    //but if so, we need to move this function to be directly in the
    //Move_()-func.)
  }
  else{
    printError("accessing invalid particle");
  }

}


inline double Lattice::computeWaitingTime(void){
  //compute waiting time:

  //one time step
  double tau;

  double r;

  if(isExponentialWaitingTime_){
    do{
      r = randomNumber.doub();
    }while( r == 1 || r == 0 );

    tau = log(1.0 / r) / partialSum_.back();
  }
  else //linear waiting time
    tau = 1.0 / partialSum_.back();

  return tau;
}


void Lattice::move(){

  double tau;      //waiting time
  int i = 0;       //index of samplingTime-vector
  timeSum_ = 0;

  while(timeSum_ < samplingTime_[maxElement_-1]){ // "-1" since we start on 0.
    //The time-sampling was completely rewritten in late August 2010
    //to resolve the spaghetti that was the previous version. This
    //follows the gillespie_exclusion2.cpp-implementation closely, to
    //make it easier on those already familiar with that source code.

    tau = computeWaitingTime();

    //Save displacement if next time-step is beyond next sampling time
    while(timeSum_ <= samplingTime_[i] &&
          samplingTime_[i] < timeSum_ + tau  && i < maxElement_){

      //save displacement (from previous step)
      dx_[i] = pos_[0].x - pos_0_.x;
      dy_[i] = pos_[0].y - pos_0_.y;
      dz_[i] = pos_[0].z - pos_0_.z;
      dr_[i] = sqrt( pow(dx_[i],2) + pow(dy_[i],2) + pow(dz_[i],2) );

      i++;
    }

    /*
    //TEST (this is similar to what my spagetti code used:)
    if(timeSum_ <= samplingTime_[i] &&  samplingTime_[i] < timeSum_ + tau &&
       i < maxElement_){

      //save displacement (from previous step)
      dx_[i] = pos_[0].x - pos_0.x;
      dy_[i] = pos_[0].y - pos_0.y;
      dz_[i] = pos_[0].z - pos_0.z;
      dr_[i] = sqrt( pow(dx_[i],2) + pow(dy_[i],2) + pow(dz_[i],2) );

      i++;
    }
    */

    double r2;          //store random number here

    int mu_guess;       // must be integer. (index of vector)
    int mu_left = 0;
    int mu_right = dim_ * 2 * nParticles_;
    double p_left = (double) partialSum_.front(); //value in first element (=0)
    double p_right = (double) partialSum_.back(); //value in last element (=dim_*2*N)

    do{
      r2 = randomNumber.doub();
    }while( r2 == 1 || r2 == 0 );

    double p_rand = (double) r2 * partialSum_.back();
    bool loopAgain = true;

    //NOTE: when comparing a long double against a double, C++ will fill
    //in the blanks with zeroes. Every 10^15 run of move() will lead to
    //a situation where we actually have p_rand==p[mu_guess], and
    //changing the interval below in the if statements,
    //from "..<= .. <.." to "..< .. <=.." will result in an infinite
    //loop every 10^15 turns since "double" stores 52 bit, 2^52 ~= 10^15.
    //This bug has now been fixed!

    //finds the mu
    do{
      mu_guess = (int) ( (p_rand - p_left) * ( mu_right - mu_left )
                         /(p_right - p_left) + mu_left );

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

    //transforms the index mu to which particle to move (n)
    //and in which direction (r, 0<= r <= 5 (in 3D))
    convertMuToParticle(mu,n,r);

    //moves the particle n in direction r, according to boundary
    //conditions. If new site is occupied, move back to previous pos.
    moveAndBoundaryCheck(n,r);

    timeSum_ = timeSum_ + tau;

  }
}







// ----------------------
// JUMPRATE SPECIFIC CODE
// ----------------------


void Lattice::setJumpRate(const vector<Jump>& jumpRate){
  if ( (int) jumpRate.size() == nParticles_){

    // remove jumprates from previous runs
    // but now we are quenched so not needed (I think) TODO (is this right?)
    jumpRate_.clear();
    jumpRate_ = jumpRate;

    //Use these vectors to build the partial sum
    computePartialSum();

  }
  else{
    printError("number of jumprates must equal number of particles");
  }
}


void Lattice::computePartialSum(void){
  long double tmp = 0;
  int i;

  // delete previous run
  partialSum_.clear();

  // first element must be zero
  partialSum_.push_back(tmp);

  //check the jump-rate-vectors
  if ( (int) jumpRate_.size() != nParticles_ )
    printError("number of jumprates, != N");

  //build partialSum-vector
  if(dim_>=1){
    for(i=0; i < nParticles_ ; i++){
      tmp = tmp + jumpRate_[i].x.r;
      partialSum_.push_back(tmp);  //index 1 ... N
    }
    for (i=0; i < nParticles_ ; i++){
      tmp = tmp + jumpRate_[i].x.l;
      partialSum_.push_back(tmp);  //index N+1 ... 2N
    }
  }
  if (dim_ >=2 ){
    for (i=0; i < nParticles_; i++ ){
      tmp = tmp+jumpRate_[i].y.r;
      partialSum_.push_back(tmp);   //index 2N+1 ... 3N
    }
    for (i=0; i < nParticles_; i++ ){
      tmp = tmp+jumpRate_[i].y.l;
      partialSum_.push_back(tmp);   //index 3N+1 ... 4N
    }
  }
  if(dim_==3){
    for (i=0; i < nParticles_; i++ ){
      tmp = tmp+jumpRate_[i].z.r;
      partialSum_.push_back(tmp);   //index 4N+1 ... 5N
    }
    for (i=0; i < nParticles_; i++ ){
      tmp = tmp+jumpRate_[i].z.l;
      partialSum_.push_back(tmp);   //index 5N+1 ... 6N
    }
  }

  //TEST: It must increase monotonically:
  for (int i = 1; i < nParticles_; i++){
    if (partialSum_[i] < partialSum_[i-1] ){
      cout <<"Warning! The cumulative sum of jump rates is not "<<endl
           <<"a monotonically increasing sum!"<<endl;
      abort();
    }
  }
}


void Lattice::convertMuToParticle(int mu, int& n, int& direction){
  //Perhaps using binary search here would speed things up?
  //Well, I'll settle for replacing some if --> else, see below

  //if no error or inconsistency...
  if (0 <= mu && mu <= (int) partialSum_.size() -1){
    //... then start checking the intervals:

    if ( 0 <= mu && mu <= 2 * nParticles_ - 1 ){
      if ( mu <= nParticles_ - 1 ){
        n = mu;
        direction = 0;
      }
      //if ( nParticles_ <= mu){
      else{
        n = mu - nParticles_;
        direction = 1;
      }
    }
    else{
      if ( 2 * nParticles_ <= mu && mu <= 4 * nParticles_ - 1){
        if ( mu <= 3 * nParticles_-1){
          n = mu - 2 * nParticles_;
          direction = 2;
        }
        //if ( 3 * nParticles_ <= mu ){
        else{
          n = mu - 3 * nParticles_;
          direction = 3;
        }
      }//if ( 4 * nParticles_ <= mu && mu <= 6 * nParticles_ - 1){
      else{
        if ( mu <= 5 * nParticles_ - 1){
          n = mu -4 * nParticles_;
          direction = 4;
        }
        //if ( 5 * nParticles_ <= mu ){
        else{
          n = mu -5 * nParticles_;
          direction = 5;
        }
      }
    }
  }
  else{
    cout << "mu =  "<< mu <<endl;
    cout << "partialSum_.size is: "<<partialSum_.size()<<endl;
    cout << "Error, in mu or partialSum_"<<endl;
    abort();
  }
}


float Lattice::computeEffectiveDiffusionConst(void){
  float D_eff;

  if(!partialSum_.empty()){

    //NOTE: if the lattice constant: a != 1, it must be included
    //D_i=k_i * pow(a,2)
    double sum_temp = 0;

    if (dim_ >= 1){
      for (int i=0; i < nParticles_; i++){
        sum_temp = sum_temp + 1.0/jumpRate_[i].x.l;
        sum_temp = sum_temp + 1.0/jumpRate_[i].x.l;
      }
    }

    if (dim_ >= 2){
      for (int i=0; i < nParticles_; i++){
        sum_temp = sum_temp + 1.0/jumpRate_[i].y.r;
        sum_temp = sum_temp + 1.0/jumpRate_[i].y.l;
      }
    }

    if (dim_ == 3){
      for (int i=0; i < nParticles_; i++){
        sum_temp = sum_temp + 1.0/jumpRate_[i].z.r;
        sum_temp = sum_temp + 1.0/jumpRate_[i].z.l;
      }
    }

    sum_temp = (double) sum_temp / (2*nParticles_*dim_);
    D_eff = 1.0 / sum_temp;

    return D_eff;
  }
  else{
    //When we don't use jumpratios we get the diffusion constant to:
    //TODO: stämmer detta? den skall ju bli <r^2>=2*D*dim_*t ! -->
    //svar: med denna D_eff får vi f(x)=2*dim_*D_eff*x, vilket är
    //inkonsekvent med ovanstående Deff som redan är kompenserad för
    //det... tror jag... TSTAR NU: N=1, uniform distrubution...
    //RESULTAT: nej, även med ovanst. D_eff utr. måste vi ha f(t)=2*dim_*D_eff*t
    //Detta är det rimmliga, eftersom vid N>1 vill vi ha den "riktiga" D_eff,
    //dvs. som det är nu.
    //Note, assumed lattice constant a = 1.
    D_eff = 1.0 / (2 * dim_);
    return D_eff;
  }
}



float Lattice::computeAverageDiffusionConst(void){
  //Calculate Average as well:

  float D_av = 0;

  if(!partialSum_.empty()){

    int startingParticle = 1;  //don't include the tagged
    for (int i = startingParticle; i < nParticles_; i++){
      D_av = D_av + jumpRate_[i].x.r;
      D_av = D_av + jumpRate_[i].x.l;
    }
    if(dim_ > 1){
      for (int i = startingParticle; i < nParticles_; i++){
        D_av = D_av + jumpRate_[i].y.r;
        D_av = D_av + jumpRate_[i].y.l;
      }
    }
    if (dim_>2){
      for (int i = startingParticle; i < nParticles_; i++){
        D_av = D_av + jumpRate_[i].z.r;
        D_av = D_av + jumpRate_[i].z.l;
      }
    }
    D_av = D_av * 1.0/(nParticles_ * 2 * dim_);

    return D_av;
  }
  else{
    //When we don't use jumpratios we get the diffusion constant to:
    //TODO: stämmer detta? den skall ju bli <r^2>=2*D*dim_*t ! -->
    //svar: med denna D_eff får vi f(x)=2*dim_*D_eff*x, vilket är
    //inkonsekvent med ovanstående Deff som redan är kompenserad för
    //det... tror jag... TSTAR NU: N=1, uniform distrubution...
    //RESULTAT: nej, även med ovanst. D_eff utr. måste vi ha f(t)=2*dim_*D_eff*t
    //Detta är det rimmliga, eftersom vid N>1 vill vi ha den "riktiga" D_eff,
    //dvs. som det är nu.
    //Note, assumed lattice constant a = 1.
    D_av  = 1.0 / (2*dim_);
    return D_av;
  }
}


float Lattice::computeNakazato(void){

  float alpha;

  if (dim_ == 1) alpha = 1;
  if (dim_ == 2) alpha = 0.363;
  if (dim_ == 3) alpha = 0.209;

  //jump-ratio for the different particles...
  double gamma = jumpTracer_ / jumpCrowders_;

  float c = (float) nParticles_ / (latticeX_ * latticeY_ * latticeZ_);

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



float Lattice::computeErgodicity(int x_rutor, int y_rutor, int z_rutor){
  //check what the MSD should be when the system reaches equilibrium.
  //(only correct if N=1). This function has two independent parts.
  //One where we just "step through" the lattice and add the
  //displacement from the center, and another were we actually
  //calculate it analytically.

  float sum = 0;

  //use these to "walk" through the lattice to measure the MSD
  //contribution from each square
  int x,y,z;


  //STEP THROUGH THE LATTICE
  //I re-map the starting coordinates, in case I'd
  //like to move this function out from this class...
  int mittX = pos_0_.x;
  int mittY = pos_0_.y;
  int mittZ = pos_0_.z;

  for (int i=0; i < x_rutor; i++){
    for (int j=0; j < y_rutor; j++){
      for (int k=0; k < z_rutor; k++){
        x = i + 1;                 //First square have coord.(1,1,1)
        y = j + 1;
        z = k + 1;
        sum = sum + pow(mittX-x,2)+pow(mittY-y,2)+pow(mittZ-z,2);
      }
    }
  }
  sum = (float) sum / (x_rutor*y_rutor*z_rutor);


  //ANALYTIC EXPRESSION
  //this is an alternate way of calculating <R^2>,
  //when the system becomes ergodic.
  bool testAnalytic = false;
  if (testAnalytic){
    int a,b,c;
    a = mittX;
    b = mittY;
    c = mittZ;
    float xContribution = ((x_rutor + 1) * (2 * x_rutor + 1) *
                           1.0/6 - a * (x_rutor+1) + a * a);
    float yContribution = ((y_rutor + 1) * (2 * y_rutor + 1) *
                           1.0/6 - b * (y_rutor+1) + b * b);
    float zContribution = ((z_rutor + 1) * (2 * z_rutor + 1) *
                           1.0/6 - c * (z_rutor + 1) + c * c);
    float SumAnalytic = xContribution + yContribution + zContribution;
    cout <<"Analytic result for <R^2> (ergodicity): "
         << SumAnalytic << endl;
  }

  return sum;
}













// -------------------------
// INTERACTION SPECIFIC CODE
// -------------------------

void Lattice::interaction(int n, Particle oldPos){
  vector<int> NearestNeighbours;
  //vector contains the particle labels of the nearest neighbors to
  //particle n (not sure that we actually need this at this point)

  //This is an updated version of the progInteraction.cpp I used in my
  //thesis It does the exact same thing, but now we can wrap around
  //the boundary also...
  //The output-files are identical between the two different version,
  //when I use fixed boundary conditions, which was what I used in my
  // thesis.
  countNeighbours(oldPos, NearestNeighbours);

  int counter = NearestNeighbours.size();

  double r0;
  static Ran slump(56);
  do{
    r0=slump.doub();
  }while( r0 == 1 || r0 == 0 );

  float DeltaV = interactionStrength_;
  if (r0 < exp(-1.0*counter*DeltaV) ){
    //accept move
    vacancy_[oldPos.x][oldPos.y][oldPos.z] = -1;   //old pos. now vacant
    vacancy_[pos_[n].x][pos_[n].y][pos_[n].z] = n; //new pos. now occupied
  }
  else{
    //reject move, move back.
    //vacancy_-matrix unchanged
    pos_[n] = oldPos;
  }
}




//IMPORTED FROM TESTInteractionCode() 12/7 -2010
//NOTE! BUT! "Exponent" changed to allow 100% probability of being included.
void Lattice::buildCluster2(int n, vector<int>& NearestNeighbours,
                              double Exponent){
  vector<int> unique;
  vector<int> NewParticleAdded;
  double rand;

  //Build a vector with particles that form the Cluster:
  //-------------------------------------------------
  //each neighbor has a probability given by the Boltzmann factor:
  //    exp(-1.0*Exponent)
  //to belong to the cluster.

  //Question: What happens when we're at the boundary?
  //Answer:   That's dealt with by the countNeighbors()-function

  //n included with 100 % probability:
  unique.push_back(n);

  bool CheckAgain;   //true if a new particle was added to the cluster

  do{
    CheckAgain = false;

    int limit = NearestNeighbours.size();
    for(int i = 0; i< limit; i++){
      int m = NearestNeighbours[i];
      countNeighbours(pos_[m], NearestNeighbours);
      //NearestNeighbours-vector will store more particles each
      //iteration, but only loop over the first initial set (ie. i < limit)
    }

    for(unsigned int j = 0; j < NearestNeighbours.size(); j++){
      bool SingleValued=true;

      //Check if it's stored in unique[] already:
      for(unsigned int k = 0; k < unique.size(); k++){
        if (NearestNeighbours[j] == unique[k])
          SingleValued = false;
      }
      if(SingleValued){
        //Add the new particle, with Boltzmann's consent...
        // rand = randomNumber.doub();
        //if (rand >=  exp(-1.0*Exponent)){
        unique.push_back(NearestNeighbours[j]);
        NewParticleAdded.push_back(NearestNeighbours[j]);
        //New particles added to the cluster, run the loop again
        //for particles in NewParticleAdded to check the neighbors
        // of this new particle
        CheckAgain = true;
        //}
      }
    }

    NearestNeighbours = NewParticleAdded;
    NewParticleAdded.clear();

  }while(CheckAgain);

  //Now save the final vector (cluster) and return it as a reference:
  NearestNeighbours = unique;


  if(false){ //Inactivated, for now
    //TEST:
    cout << "kollar Nearest Neighbour.size(): "<<NearestNeighbours.size()<<endl;
    for(unsigned int j=0; j < NearestNeighbours.size(); j++){
      int m = NearestNeighbours[j];
      cout <<"m = " << m << " (" << pos_[m].x << "," << pos_[m].y
           << "," << pos_[m].z << ")"<< endl;
    }
  }

  //TEST
  for (int i=0; i < nParticles_; i++){
    if (pos_[i].x==0 || pos_[i].y==0 || pos_[i].z==0)
      printError("Outside of lattice");
  }

  //TEST---------(FÖR ATT KOLLA DOUBLE COUNTING)
  bool IsUniqueTrue = true;
  for(unsigned int i = 0; i < unique.size(); i++){
    for(unsigned int j = 0; j < unique.size(); j++){
      if (unique[i] == unique[j] && i != j)
        IsUniqueTrue = false;
    }
  }
  if(!IsUniqueTrue)
    printError("Unique is not unique! Error *235897*");
  //-------------
}




void Lattice::countNeighbours(Particle particle, vector<int>& Count){
  //Each element in the vector "Count" is a particle, and the
  //value of the element is the particle label.
  //If the particle is "alone" the vector Count will be empty.
  //boundaryFix_=0 periodic, =1 fix

  //TODO invert the boolean condition, not that pretty perhaps:
  bool periodic;

  if (!boundaryFix_) periodic = true;
  else  periodic = false;

  if (particle.x+1 <= latticeX_ ){ //avoid falling off the edge of the world
    if (vacancy_[particle.x+1][particle.y][particle.z] != -1){
      Count.push_back( vacancy_[particle.x+1][particle.y][particle.z] );
    }
  }
  else{
    if (vacancy_[1][particle.y][particle.z] != -1 && periodic){
      Count.push_back( vacancy_[1][particle.y][particle.z] );
    }
  }

  if (particle.x-1>0){
    //although this zero-element exists in the vacancy_ matrix (i.e. no seg. fault)
    //the first coordinate on the lattice is (1,1,1). (0,0,0) is undefined
    if (vacancy_[particle.x-1][particle.y][particle.z] != -1){
      Count.push_back( vacancy_[particle.x-1][particle.y][particle.z] );
    }
  }
  else{
    if (vacancy_[latticeX_][particle.y][particle.z] != -1 && periodic){
      Count.push_back( vacancy_[latticeX_][particle.y][particle.z] );
    }
  }



  if (particle.y+1 <= latticeY_){
    if (vacancy_[particle.x][particle.y+1][particle.z] != -1){
      Count.push_back( vacancy_[particle.x][particle.y+1][particle.z] );
    }
  }
  else{
    if (vacancy_[particle.x][1][particle.z] != -1 && periodic){
      Count.push_back( vacancy_[particle.x][1][particle.z] );
    }
  }

  if (particle.y-1>0){
    if (vacancy_[particle.x][particle.y-1][particle.z] != -1){
      Count.push_back( vacancy_[particle.x][particle.y-1][particle.z] );
    }
  }
  else{
    if (vacancy_[particle.x][latticeY_][particle.z] != -1 && periodic){
      Count.push_back( vacancy_[particle.x][latticeY_][particle.z] );
    }
  }


  if (particle.z+1 <= latticeZ_){
    if (vacancy_[particle.x][particle.y][particle.z+1] != -1){
      Count.push_back( vacancy_[particle.x][particle.y][particle.z+1] );
    }
  }
  else{
    if (vacancy_[particle.x][particle.y][1] != -1 && periodic){
      Count.push_back( vacancy_[particle.x][particle.y][1] );
    }
  }

  if (particle.z-1>0){
    if (vacancy_[particle.x][particle.y][particle.z-1] != -1){
      Count.push_back( vacancy_[particle.x][particle.y][particle.z-1] );
    }
  }
  else{
    if (vacancy_[particle.x][particle.y][latticeZ_] != -1 && periodic){
      Count.push_back( vacancy_[particle.x][particle.y][latticeZ_] );
    }
  }
}







// ----------------------------------
// NOT COMPLETE YET. WORK IN PROGRESS
// ----------------------------------


void Lattice::snapshot(void){
  //NOTE! This is only for 2D, with isotropic jumprate, and the same
  //system side-length! Its purpose is to print a snapshot of the
  //system, to see the spatial distribution of different particles
  //as a heat map. Only used this function briefly.

  vector< vector<double> > matris;

  vector<double> temp;
  temp.assign(latticeY_,0);
  matris.assign(latticeX_,temp);

  for(int i=0; i < nParticles_; i++){
    matris[ pos_[i].x-1 ][ pos_[i].y-1 ] = jumpRate_[i].x.r; //TODO BUG!
  }

  ofstream heat;
  char NameOfFiles[]="heatmap.dat";
  heat.open(NameOfFiles);

  //true: plotar endast partiklarna, false: plottar aven tomma rutor (=0)
  bool points=true;

  if(points){
    //This does not give 0 where there are no particles:
    for (int i=0; i < nParticles_; i++)
      heat  << pos_[i].x<<" \t"<< pos_[i].y <<"\t"<< jumpRate_[i].x.r<< endl; //TODO BUG
  }
  else{
    //This gives zeroes for all empty sites:
    for(int i=0; i<latticeX_; i++){
      for(int j=0; j<latticeY_; j++){
        heat << i <<"\t"<< j <<"\t"<< matris[i][j]<<endl;
      }
    }
  }

  heat.close();
  cout<<"Pinted Heatmap"<<endl;
}

