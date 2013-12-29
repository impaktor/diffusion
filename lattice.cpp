#include <iostream>        //for cout, among others
#include <fstream>         //read write to files
#include <cstdlib>         //for abort()-function
#include <cmath>           //gives sqrt, pow, fabs,
#include <algorithm>       //needed for min_element in hist/bin-function
#include <sstream>         //behövs kanske för "ostringstream" i save()-func.
#include <string>
#include <vector>

#include "nr/nr3.h"        //ran.h, ludcmp.h, and svd.h depends on this one.
#include "nr/ran.h"        //Random number generator from numerical recipes 3ed.
#include "nr/ludcmp.h"     //only depends on nr3.h, used in computePearson 
#include "nr/svd.h"        //singlular value decomposition

#include "classes.h"       //various data structures /classes (Jump, Particle)
#include "lattice.h"       //this is the main class, that does the physics. 

using namespace std;

//-------------------
//HARDCODED VARIABLES:
//-------------------
#define TEST_1            0  //boolean, prints details to screen | (constructor)
#define CHECK_VACANCY_ON  1  //CPU consuming, use only when bug-squashing!

#define SEED_WAITING_TIME 87 //My seeds (are usually 87,1,8)
#define SEED_PLACE        1
#define SEED_MOVE         8

int Lattice::computeCorrelation(int ensembles){
  //This function calculates the Pearson coefficients (see Numerical
  //recipies ch 14), to get the correlation between two sampling points
  //since if the tagged particle jumps once/sec sampling the system 
  //10/s will not yeald any new information, and it will mess up the 
  //error estimation that assumes the sampling points are uncorrelated 
  //and with a gaussian distribution.

  //perhaps check: Physical Review D vol 34, 2469 (1986), sec. D  
  
  if(lowMem_)
    printError("Can not compute correlation, no available data.");

  if(r2_mu.empty())
    printError("Must compute MSD first! (ie. r2_mu)");

  //Either compute the Pearson correlation coefficient, or use the full matrix
  if(false){
    
    vector<double> pearsonCorrelation;
    pearsonCorrelation.assign(maxElement_,0);

    //compute the correlation between two points, i, j, for i = j+1
    vector<double> h;
    h.assign(maxElement_,0);

    for(int i = 1; i < maxElement_; i++){
      for(int j = 0; j < ensembles; j++){
        h[i] += 1.0 / (ensembles-1) *
          (pow(store_dr[j][i],2) - r2_mu[i])*
          (pow(store_dr[j][i-1],2) - r2_mu[i-1]);
      }
    }
  
    // pearsonCorrelation[i] = h_ij / (sigma_i*sigma_j)
    // sigma_r = r2_err * sqrt(ensembles);
    for(int i = 1; i < maxElement_; i++){
      if(r2_err[i-1] != 0)
        pearsonCorrelation[i] = h[i] / (r2_err[i]*r2_err[i-1] * ensembles);
      
      //NOTE: error in first value is 0, so the correlation between value 0 & 1
      //will be inf (nan)! thus the if-statement
    
      //Now save the inf. for printing to the output-file.
      fourthColumn_ = pearsonCorrelation;
    }
  }
  else{
    //compute entire matrix, to check the correlation between every value

    if(samplingTime_[0] == 0)
      cout << "NOTE: If the first sampling-point t_0 = 0 then the \n"
           << "H-matrix will be singular (first row and column =0)"
           << endl;

    NRmatrix<double> H(maxElement_, maxElement_, 0.0);
    
    for(int i = 0; i < maxElement_; i++){
      for(int j = 0; j < maxElement_; j++){
        for(int m = 0; m < ensembles; m++){
          H[i][j] += 1.0 / (ensembles*(ensembles-1)) *
            ( pow(store_dr[m][i],2) - r2_mu[i] ) *
            ( pow(store_dr[m][j],2) - r2_mu[j] );
        }
      }
    }
    



    //Now test the integrety of the matrix:
    SVD testMatrix(H);
    double conditionalNumber = testMatrix.inv_condition();
    cout << "conditional number is: " << setprecision(5) << scientific
         << conditionalNumber << endl;

    //OBSERVE: For some retarded reason it rounds off the answer to the 
    //conditional number! i.e gives 0.0, instead of 3.5E-5 etc. 
    //therefore, must use: << scientific << 
 
    if(false){//TEST:
      NRmatrix<double> hm(2,2,0.0);
      hm[0][0] = 1.0e-16;
      hm[0][1] = 1.0;
      hm[1][0] = 1.0;
      hm[1][1] = 1.0;
      SVD testMatrix2(hm);
      double conditionalNumber2 = testMatrix2.inv_condition();
      cout << fixed << setprecision(5) << "Test-conditional number is : "
           << conditionalNumber2 << endl;
    }
    
    //TEST, print matrix to file:
    //TODO: have it print the same name as the one set by the argument
    //with an suffix *_matrix.dat
    ofstream dump;
    dump.open("out_H-matrix.dat");
    for (int i = 0; i < maxElement_; i++){
      for (int j = 0; j < maxElement_; j++)
        dump  <<H[i][j] << "\t";
      
      dump << endl; // << endl;
    }
    dump.close();


    NRvector<double> b(maxElement_),z(maxElement_);
    
    for(unsigned int i = 0; i < maxElement_; i++)
      b[i] = samplingTime_[i];

    LUdcmp LUdecomposition(H);
    LUdecomposition.solve(b,z);

    //copy z to the class-private vector, that we will print as the 
    //fourth column in the output-file. Transfer a NR vector to std::
    fourthColumn_.assign(maxElement_,0);
    for (int i = 0; i < maxElement_; i++)
      fourthColumn_[i] = z[i];
    
  }//end-else
  
  return 0;

}


void Lattice::dumpSimulation(int index){
  //This is a test method to simply "dump" the trajectories
  //for a single run to file. It is used to check correlations
  //between different runs. 
  
  ofstream dump;
  dump.open("dump.dat");
  
  for (int i = 0; i < maxElement_; i++){
    dump  << samplingTime_[i] << "\t" << dr[i] << endl;
  }
  dump.close();
  
  char buffer[50];
  sprintf(buffer,"mv dump.dat dump%d.dat",index);

  cout << "Moving dump.dat to dump" << index << ".dat." << endl;
  system(buffer);
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
    printError("gaaa");

  ofstream clust;
  char NameOfFiles[]="cluster.dat";
  clust.open(NameOfFiles);
  for(unsigned int j=0; j < clusterDistribution_.size(); j++){
    clust << clusterSize_[j] <<"\t"<< clusterDistribution_[j]<<endl;
  }
  clust.close();
  cout<<"Pinted "<<NameOfFiles<<endl;
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
    if(error) cout <<"Does not compute..."<<endl;

    if (count_vacant_sites != latticeX_*latticeY_*latticeZ_-nParticles_){
      cout << "Error in vacant sites" <<endl
           <<"# of vacant sites: " <<count_vacant_sites <<endl;
      error = true;
    }

    if(error){
      cout << "Vacancy Matrix is not correct, at err. code site: "
           << ErrorCode<<endl;
      abort();
    }
  }
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

void Lattice::setSamplingTimes(vector<double>& samplingTime, 
                               bool isExpWaitingTime = false){

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

    //assign the displacement container.
    dx.assign(maxElement_,0);
    dy.assign(maxElement_,0);
    dz.assign(maxElement_,0);
    dr.assign(maxElement_,0);
 
    //initiate the vectors for the summation in store()
    //only needed if we use the memory-conserving algorithm
    if(lowMem_){
      dr4_err.assign(maxElement_,0);
      dx4_err.assign(maxElement_,0);
      dy4_err.assign(maxElement_,0);
      dz4_err.assign(maxElement_,0);
      
      r2_mu.assign(maxElement_,0);
      x2_mu.assign(maxElement_,0);
      y2_mu.assign(maxElement_,0);
      z2_mu.assign(maxElement_,0);
    }
    
  }
  else
    printError("Sampling times is not ordered.");
}


void Lattice::setJumpNaka(float jumpCrowders, float jumpTracer){
  jumpCrowders_ = jumpCrowders;
  jumpTracer_ = jumpTracer;
}

void Lattice::setDist(int distribution, double info){
  distribution_ = distribution;
  info_ = info;
}

void Lattice::setInteraction(float interactionStrength){
  interactionStrength_ = interactionStrength; 
  interactionOn_ = true; 
}

void Lattice::setLowMem(bool lowMem){
  //Set the switch true/false for use of low-mem run
  lowMem_ = lowMem;
}



// ---------------------
// PUBLIC: GET FUNCTIONS:
// ---------------------

int Lattice::getDimension(void){ //could be useful to have, but not necessary...
  return dim_;
}


// -----------
// OTHER STUFF
// -----------

float Lattice::nakazato(void){
  
  float alpha;

  if (dim_ == 1) alpha = 1;
  if (dim_ == 2) alpha = 0.363;
  if (dim_ == 3) alpha = 0.209;

  //jump-ratio for the different particles... 
  double gamma = jumpTracer_ / jumpCrowders_;

  float c = concentration_; 

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



// ===================== TODO: move to top of file
// C O N S T R U C T O R
// =====================

Lattice::Lattice(int xSquare,int ySquare,int zSquare, 
                     int particleNumber, bool boundary){

  if (xSquare >= ySquare && ySquare  >= zSquare && zSquare >= 1){
    if (particleNumber <= xSquare*ySquare*zSquare){

      latticeX_ = xSquare;            //number of lattice sites
      latticeY_ = ySquare;
      latticeZ_ = zSquare;

      nParticles_ = particleNumber;

      boundaryFix_ = boundary;

      pos_0_.x = round(1.0*latticeX_/2);//place tagged particle at lattice center
      pos_0_.y = round(1.0*latticeY_/2); 
      pos_0_.z = round(1.0*latticeZ_/2); 
      
      //TODO: what if Z = 0
      concentration_ = (float) nParticles_ / (latticeX_*latticeY_*latticeZ_);

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

      //By default: don't use the memory conserving stdErrLowMem()
      //use the set_LowMem()-function in the main()-func to change this
      lowMem_ = false;      

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



float Lattice::ergodicity(int x_rutor, int y_rutor, int z_rutor){ 
  //check what the MSD should be when the system reaches equilibrium.
  //(only correct if N=1). This function has two independent parts. 
  //One where we just "step through" the lattice and add the 
  //displacement from the center, and another were we actually 
  //calculate it analytically. 
  
  float summa = 0; 
  
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
        summa = summa + pow(mittX-x,2)+pow(mittY-y,2)+pow(mittZ-z,2);
      }
    }
  }
  summa = (float) summa / (x_rutor*y_rutor*z_rutor);


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

  return summa;
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
  
  double R;
  static Ran randomNumber1(SEED_PLACE);    //just any seed will do.
  int i, j, k, n = 1;             //First crowder at n = 1 element   
  for (k=1; k <= latticeZ_; k++){    //place particles
    for (j=1; j <= latticeY_; j++){
      for (i=1; i <= latticeX_; i++){

        //don't place crowders on the tagged particle
        if (i != pos_0_.x || j != pos_0_.y || k != pos_0_.z){ 

          //Using NR to generate 0 < R < 1 ;
          R=randomNumber1.doub();              
          if (R<=leftToPlace*1.0/totalSites){
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



void Lattice::moveOld(){ 
  //This code is obsolete, as it just works for placing identical
  //particles, but I keep it for sentimental reasons.
  //Better safe than sorry.
  
  static Ran randomNumber2(2);    //any seed will do. 
  int r,n; 
  
  //Choose particle at random, 0 <= n < N
  do{
    n = randomNumber2.doub() * nParticles_; 
  }while(n == nParticles_);                  
      

  do{
    r=randomNumber2.doub() * (2 * dim_); 
  }while( r == dim_ * 2 );           //gives 0<=R<2*dim
      

  //Do the actual move, according to the boundary rules
  moveAndBoundaryCheck(n,r);

  timeSum_ = timeSum_ + 1.0 / nParticles_;  //Calculate the time

  if (testOnOff_)
    for (int i=0; i < nParticles_; i++)
      cout <<" MOVED \n particle n = " << i <<": "<< pos_[i] << endl;
}


void Lattice::moveAndBoundaryCheck(int n, int R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant. 
  bool boundary = boundaryFix_;
  
  if (0 <= n && n < nParticles_ ){
    if (testOnOff_) cout <<"Moving particle "<< n <<endl;
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
        if (pos_[n].x == latticeX_ +1 && boundary == 0) pos_[n].x = 1;
        if (pos_[n].x == latticeX_ +1 && boundary == 1) pos_[n].x = old.x;
        break;
      case 1: pos_[n].x = pos_[n].x -1; //move particle n to the left
        if (pos_[n].x == 0 && boundary == 0) pos_[n].x = latticeX_;//periodic
        if (pos_[n].x == 0 && boundary == 1) pos_[n].x = old.x;//fix 
        break;  //(first coordinate is 1, 0=outside) 
      case 2: pos_[n].y = pos_[n].y +1;
        if (pos_[n].y == latticeY_ +1 && boundary == 0) pos_[n].y = 1;
        if (pos_[n].y == latticeY_ +1 && boundary == 1) pos_[n].y = old.y;
        break;
      case 3: pos_[n].y = pos_[n].y -1;
        if (pos_[n].y == 0 && boundary == 0) pos_[n].y = latticeY_;
        if (pos_[n].y == 0 && boundary == 1) pos_[n].y = old.y;
        break;
      case 4: pos_[n].z = pos_[n].z +1;
        if (pos_[n].z == latticeZ_+1 && boundary == 0) pos_[n].z = 1;
        if (pos_[n].z == latticeZ_+1 && boundary == 1) pos_[n].z = old.z;
        break;
      case 5: pos_[n].z = pos_[n].z-1;
        if (pos_[n].z == 0 && boundary == 0) pos_[n].z = latticeZ_;
        if (pos_[n].z == 0 && boundary == 1) pos_[n].z = old.z;
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
    if (testOnOff_) cout <<"Vacant, particle "<< n<<": "<<pos_[n]<<endl;
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



int Lattice::vacancyCheck(int n, Particle oldPos){
  //NOTE: This is a faster version than vacancyCheckOld(), since that 
  //one loops through all particles for each "vacancy check". However 
  //this uses a 3D matrix, and there will be some "double banking" 
  //in the sense that the position is now stored in this matrix and 
  //in XYZ-vectors.

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
        
        if (testOnOff_) cout <<"OK, vacant"<<endl;
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

  double tau;

  static Ran randomNumb2(SEED_WAITING_TIME);
  double r3;
  
  if(isExponentialWaitingTime_){
    do{
      r3 = randomNumb2.doub();
    }while( r3 == 1 || r3 == 0 );
    
    tau = log(1.0/r3)/partialSum_.back();
  }
  else //linear waiting time
    tau = 1.0/partialSum_.back(); 
  
  return tau;
}


void Lattice::move(){

  double tau;      //waiting time
  int i = 0;      //index of samplingTime-vector
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
      dx[i] = pos_[0].x - pos_0_.x;
      dy[i] = pos_[0].y - pos_0_.y;
      dz[i] = pos_[0].z - pos_0_.z;
      dr[i] = sqrt( pow(dx[i],2) + pow(dy[i],2) + pow(dz[i],2) ); 
      
      i++;
    }

    /*  
    //TEST (this is similar to what my spagetti code used:)
    if(timeSum_ <= samplingTime_[i] &&  samplingTime_[i] < timeSum_ + tau && 
       i < maxElement_){

      //save displacement (from previous step)
      dx[i] = pos_[0].x - pos_0.x;
      dy[i] = pos_[0].y - pos_0.y;
      dz[i] = pos_[0].z - pos_0.z;
      dr[i] = sqrt( pow(dx[i],2) + pow(dy[i],2) + pow(dz[i],2) ); 
      
      i++;
    }
    */
    
   
    //TODO remove! -------------------------------------------
    //pick a particle & direction based on jump-rate, and move.
    if (partialSum_.empty())
      printError("Partial-sum vector has not been initiated");

    static Ran randomNumb(SEED_MOVE);
    double r2;

    int mu_guess;                       // must be integer. (index of vector)
    int mu_left = 0;
    int mu_right = dim_*2*nParticles_;
    double p_left = (double) partialSum_.front(); //value in first element (=0)
    double p_right = (double) partialSum_.back(); //value in last element (=dim_*2*N)
 
    do{
      r2 = randomNumb.doub();
    }while( r2 == 1 || r2 == 0 );

    double p_rand = (double) r2*partialSum_.back();
    bool LoopAgain=true;

    //NOTE: when comparing a long double against a double, C++ will fill
    //in the blanks with zeroes. Every 10^15 run of move() will lead to
    //a situation where we actually have p_rand==p[mu_guess], and
    //changing the interval below in the if statements,
    //from "..<= .. <.." to "..< .. <=.." will result in an infinite 
    //loop every 10^15 turns since "double" stores 52 bit, 2^52 = 10^15.
    //This bug has now been fixed!

    //finds the mu
    do{ 
      mu_guess = (int) ( (p_rand - p_left) * ( mu_right - mu_left )
                         /(p_right - p_left) + mu_left );

      if ((double) partialSum_[mu_guess] <= p_rand && p_rand < 
          (double) partialSum_[mu_guess+1])
        LoopAgain = false;

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
    }while(LoopAgain);

    int mu = mu_guess;
    int n,r;        

  
    //transforms the index mu to which particle to move 
    //(index n, by reference) and in which direction,by returning 
    //an integer 0<= r <= 5 (in 3D)
    r = convertMuToParticle(mu,n);
  
    //moves the particle n in direction r, according to boundary 
    //conditions. If new site is occupied, move back to previous pos. 
    moveAndBoundaryCheck(n,r);






    timeSum_ = timeSum_ + tau;

  }
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
  static Ran RandomNUM(42);
  double randomNumber;

  //Build a vector with particles that form the Cluster:
  //-------------------------------------------------
  //each neighbor has a probability given by the Boltzmann factor:
  //    exp(-1.0*Exponent)
  //to belong to the cluster.

  //TODO! What happens when we're at the boundary?
  //Answer: that's dealt with by the countNeighbors()-function

  //n included with 100 % probability:
  unique.push_back(n);

  bool CheckAgain;   //true if a new particle was added to the cluster

  do{
    CheckAgain=false;

    int limit = NearestNeighbours.size();
    for(int i=0; i< limit; i++){
      int m = NearestNeighbours[i];
      countNeighbours(pos_[m], NearestNeighbours);
      //NearestNeighbours-vector will store more particles each 
      //iteration, but only loop over the first initial set (ie. i < limit)
    }

    for(unsigned int j=0; j < NearestNeighbours.size(); j++){
      bool SingleValued=true;
 
      //Check if it's stored in unique[] already:
      for(unsigned int k=0; k < unique.size(); k++){
        if (NearestNeighbours[j] == unique[k])
          SingleValued=false; 
      }
      if(SingleValued){
        //Add the new particle, with Boltzmann's consent...
        // randomNumber = RandomNUM.doub();
        //if (randomNumber >=  exp(-1.0*Exponent)){
        unique.push_back(NearestNeighbours[j]);
        NewParticleAdded.push_back(NearestNeighbours[j]);
        //New particles added to the cluster, run the loop again
        //for particles in NewParticleAdded to check the neighbors
        // of this new particle
        CheckAgain=true;
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



// ---------------------
// SAVE AND STORE OUTPUT
// ---------------------

void Lattice::store(void){
  //This is needed to calculate error bars in computeStdErr().
  
  /*
    Error in each point i is (in LaTeX notation):
    (a[i])^2 = 1/N \sum_i^N (a_i - <a>)^2 = [1/N \sum_i^N a[i]^2] - [<a>]^2

    <a>= 1/N \sum_i^N a[i]   (i.e. average value in that point) 
    N = Number of ensembles

    The original way to do this is save the value in each point a[i] for 
    each ensemble, but that leaves us with a giant matrix if we have many
    ensembles N, and many data points D, (N*D ~> 10^9) ~ 1Gb of mem usage
    (10^9/2^30 ~ 1Gb). Therefore with the switch 'lowMem=true' we do the 
    left hand side of the equation above, and then we don't need a matrix
    but just two vectors, for each dimension. 
  */

  //Check that we assigned the vectors =0 if we uses the lowMem-option.
  //For lowMem=false, we do this in the Save & computeStdErr()-functions
  //since we use push_back() instead.
  if(lowMem_){
    bool error = false;
    if (dr4_err.size() == 0) error = true;
    if (dx4_err.size() == 0) error = true;
    if (dy4_err.size() == 0) error = true;
    if (dz4_err.size() == 0) error = true;

    if (r2_mu.size() == 0) error = true;
    if (x2_mu.size() == 0) error = true;
    if (y2_mu.size() == 0) error = true;
    if (z2_mu.size() == 0) error = true;
   
    if(error){
      cout<<"Incorrect assignment of vectors for stdErr-LowMem"
          <<endl<<"Must be assigned by setLowMem()-function!"<<endl;
      abort();
    }
  }

  if(lowMem_ == false){
    //if don't use the lowMem-algorithm:
    //this gives dx[ensembler][number of values] 
    store_dx.push_back(dx);   
    store_dy.push_back(dy);
    store_dz.push_back(dz);
    store_dr.push_back(dr);
    
    //TODO, nollas dessa efter varje store?
    /*
      nej, initieras i setSamplingTime(), men anvander inte push back!
      dvs borde vara OK.

     */
  }
  else{
    //To conserve memory we do the summation directly, which might
    //cause a numerical error in the end compared to doing it the 
    //more rigorous way, but all tests I've made returns identical 
    //output
    for(int i = 0; i < maxElement_; i++){
      dr4_err[i] = dr4_err[i] + pow( dr[i] ,4);
      dx4_err[i] = dx4_err[i] + pow( dx[i] ,4);
      dy4_err[i] = dy4_err[i] + pow( dy[i] ,4);
      dz4_err[i] = dz4_err[i] + pow( dz[i] ,4);

      r2_mu[i] = r2_mu[i] + pow( dr[i] ,2);
      x2_mu[i] = x2_mu[i] + pow( dx[i] ,2);
      y2_mu[i] = y2_mu[i] + pow( dy[i] ,2);
      z2_mu[i] = z2_mu[i] + pow( dz[i] ,2);
    }
  }
}



void Lattice::save(int Ensemble, char name[]){
  //calculates the MSD from the stored positions of the tracer, 
  //and prints to file. 

  //print the information to file:.
  char* Primary = name;

  ofstream primary, secondary;
  primary.open(Primary);

  //hardcoded switch to print x,y,z, x^2,y^2,z^2
  //specifically. Not needed now, since system is isotropic,
  bool xyz = false;
  if(xyz) secondary.open("txyz.dat");
  
  //sore averaged displacement <x>, <y>, <z>, <r>
  //mostley used to check that is is ~ 0.
  vector<double> x_mu, y_mu, z_mu, r_mu;

  if(!lowMem_){
    //Compute <R>:
    x_mu.assign(maxElement_,0);
    y_mu.assign(maxElement_,0);
    z_mu.assign(maxElement_,0);
    r_mu.assign(maxElement_,0);

    for (int j=0; j < maxElement_; j++){  
      for (int i = 0; i < Ensemble ; i++){
        x_mu[j] = x_mu[j] + 1.0*store_dx[i][j] / Ensemble;
        y_mu[j] = y_mu[j] + 1.0*store_dy[i][j] / Ensemble;
        z_mu[j] = z_mu[j] + 1.0*store_dz[i][j] / Ensemble;
        //r_mu[j] = r_mu[j] + 1.0*store_dr[i][j] / Ensemble; //not using this
      } 
    }

    //Compute <R^2>:
    x2_mu.assign(maxElement_,0);
    y2_mu.assign(maxElement_,0);
    z2_mu.assign(maxElement_,0);
    r2_mu.assign(maxElement_,0);

    for (int j = 0; j < maxElement_; j++){  
      for (int i = 0; i < Ensemble;i++){
        x2_mu[j] = x2_mu[j] + 1.0 * pow(store_dx[i][j],2) / Ensemble;
        y2_mu[j] = y2_mu[j] + 1.0 * pow(store_dy[i][j],2) / Ensemble;
        z2_mu[j] = z2_mu[j] + 1.0 * pow(store_dz[i][j],2) / Ensemble;
        r2_mu[j] = r2_mu[j] + 1.0 * pow(store_dr[i][j],2) / Ensemble;
      }
    }

    //Calculate standard error, to generate error-bars
    //this function must come after x2_mu...r2_mu has been
    //computed.
    computeStdErr(Ensemble);

    computeCorrelation(Ensemble);

  }
  else{
    //use the memory conserving version instead. This one does not
    //use the dx[][]-matrix, and in that way we save lots of mem.
    computeStdErrLowMem(Ensemble);
  }


  //To print which distribution we used to head of output-file
  string dist[] = {"uniform","exponential","powerlaw","nakazato"};
  string onOff[] = {"Off","On"};
  string bound[] = {"periodic","fix"};
  float d_eff = diffEffConstant();  //calculate effective diffusion constant
  float d_av = D_av_;

  //TEST
  //cout <<"# Lower limit: "<<(latticeX_ * latticeX_) / (2*N*N*d_eff) <<endl;
  //cout <<"# Upper limit: "<<(latticeX_ * latticeX_) / (2*d_eff) <<endl;

  //print three lines of info about simulation to head of file:
  ostringstream print;
  print << "#E = "<<Ensemble<<"\t N = "<<nParticles_<<"\t X-Y-Z: "
        << latticeX_<<"x"<<latticeY_<<"x"<<latticeZ_
        << "\t 2*d*D_naka: "<<nakazato()<<endl;
  print << "#Conc.: "<<concentration_<<"\t MSD_equil: "
        << ergodicity(latticeX_,latticeY_,latticeZ_)
        << "\t distr: "<<dist[distribution_]<<" ("<<info_
        << ")\t bound: "<<bound[boundaryFix_]<<endl;
  print << "#D_eff: "<<d_eff<<"\t k_tagg: "
        << jumpRate_[0].x.r<<"\t D_av: "<<d_av
        << "\t Interaction: "<<onOff[interactionOn_]<<" "
        << interactionStrength_<<endl; 
  
  primary << print.str(); //convert to string and send
  cout << endl << print.str();


  if(lowMem_)  //NOTE: temporary, if lowMem is on!
    //TODO: calculate this even if lowMem_ = true!
    //does not serve any purpose at this point. 
    fourthColumn_.assign(maxElement_,0);
  
  //print data to file:
  for (int i = 0; i < maxElement_; i++){
    
    //Maybe add this:
    //cout << setprecision (5) << r2_mu[i] << ...

    //Print it to file:
    primary  << samplingTime_[i] << "\t" << r2_mu[i] << "\t" 
             << r2_err[i] << "\t" << fourthColumn_[i] << endl;
    if(xyz){
      secondary << samplingTime_[i] << "\t" 
                << x_mu[i] << "\t" << y_mu[i] << "\t" << z_mu[i] <<"\t"
               << x2_mu[i] <<"\t"<< y2_mu[i] <<"\t"<< z2_mu[i] <<"\t"
               << x2_err[i]<<"\t"<<y2_err[i]<<"\t"<<z2_err[i]<<endl;
    }
  }
  primary.close();
  if(xyz) secondary.close();
}


void Lattice::computeStdErr(int Ensemble){
  //Calculate the standard deviation, to get error-bars.
  
  //Should not be using this function if lowMem = true!
  if(lowMem_)
    printError("LowMem switch ON, but using wrong stdErr-function!");
  
  x2_err.assign(maxElement_,0);
  y2_err.assign(maxElement_,0);
  z2_err.assign(maxElement_,0);
  r2_err.assign(maxElement_,0);


  for (int j=0; j < maxElement_; j++){
    for (int i=0; i<Ensemble; i++){
      x2_err[j] += pow (pow(store_dx[i][j],2) - x2_mu[j] ,2);
      y2_err[j] += pow (pow(store_dy[i][j],2) - y2_mu[j] ,2);
      z2_err[j] += pow (pow(store_dz[i][j],2) - z2_mu[j] ,2);
      r2_err[j] += pow (pow(store_dr[i][j],2) - r2_mu[j] ,2);
    }//( "a += b" is the same as "a = a + b")
  }

  for (int j=0; j < maxElement_  ;j++){
    x2_err[j] = sqrt( x2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
    y2_err[j] = sqrt( y2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
    z2_err[j] = sqrt( z2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
    r2_err[j] = sqrt( r2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
  } //for E=1 --> inf !

  
  //Testing the Gaussian distribution of the points:
  //--------------------------------------------------
  float SUM68 = 0;
  for (int j = 0; j < maxElement_; j++){    //forgot about his: 
    double top = x2_mu[j] +0.5*x2_err[j] * sqrt(Ensemble-1)*2;
    double bot = x2_mu[j] -0.5*x2_err[j] * sqrt(Ensemble-1)*2;
    for (int i=0; i < Ensemble; i++){
      if (bot <= pow(store_dx[i][j],2) && 
          top >= pow(store_dx[i][j],2) ){
        bool test = false;
        if (test){ //TEST-function
          cout <<"lower limit: "<< bot <<endl;
          cout << "point is: "<< pow(store_dx[i][j],2) <<endl;
          cout << "upper limit: "<< top <<endl;
        }
        SUM68++;
        if (test) cout <<"SUM68 = "<<SUM68<<endl;
      }
    }
  }
  cout <<"68,2 % (?)="<<1.0*SUM68/(Ensemble*maxElement_)<<endl;

  //--------------------------------------------------
}


void Lattice::computeStdErrLowMem(int Ensemble){

  if(!lowMem_)
    printError("lowMem = OFF, but using stdErrLowMem()-func!");

  //Divide the summed up displacement to get the average:
  for (int i=0; i< maxElement_; i++){
    x2_mu[i] = x2_mu[i] / (Ensemble * 1.0);
    y2_mu[i] = y2_mu[i] / (Ensemble * 1.0);
    z2_mu[i] = z2_mu[i] / (Ensemble * 1.0);
    r2_mu[i] = r2_mu[i] / (Ensemble * 1.0);
  }

  x2_err.assign(maxElement_,0);
  y2_err.assign(maxElement_,0);
  z2_err.assign(maxElement_,0);
  r2_err.assign(maxElement_,0);

  //Following the notation in the comments in Lattice::store():
  // [1/N sum_i^N a_i^2] - [<a>^2], (a=dr^2)
  for(int i =0; i < maxElement_; i++){
    x2_err[i] = dx4_err[i] / (1.0 * Ensemble) - pow( x2_mu[i] ,2);
    y2_err[i] = dy4_err[i] / (1.0 * Ensemble) - pow( y2_mu[i] ,2);
    z2_err[i] = dz4_err[i] / (1.0 * Ensemble) - pow( z2_mu[i] ,2);
    r2_err[i] = dr4_err[i] / (1.0 * Ensemble) - pow( r2_mu[i] ,2);
  }

  for(int i =0; i < maxElement_; i++){
    x2_err[i] = sqrt(x2_err[i] / (1.0 * (Ensemble - 1) ));
    y2_err[i] = sqrt(y2_err[i] / (1.0 * (Ensemble - 1) ));
    z2_err[i] = sqrt(z2_err[i] / (1.0 * (Ensemble - 1) ));
    r2_err[i] = sqrt(r2_err[i] / (1.0 * (Ensemble - 1) ));
  }
}




// ----------------------
// JUMPRATE SPECIFIC CODE
// ----------------------


void Lattice::setJumpRate(vector<Jump>& jumpRate){
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
 

int Lattice::convertMuToParticle(int mu, int& n){

  if (0 <= mu && mu <= (int) partialSum_.size() -1){

    if ( 0 <= mu && mu <=2*nParticles_-1 ){
      if ( mu<=nParticles_-1 ){               
        n = mu;
        return 0;
      }
      if ( nParticles_ <= mu){               
        n = mu - nParticles_;
        return 1;
      }
    }

    if ( 2*nParticles_ <= mu && mu <= 4*nParticles_-1){
      if ( mu <= 3*nParticles_-1){
        n = mu -2*nParticles_;
        return 2;
      }
      if ( 3*nParticles_ <= mu ){
        n = mu -3*nParticles_;
        return 3;
      }
    }

    if (4*nParticles_ <= mu && mu <= 6*nParticles_-1){
      if ( mu <= 5*nParticles_-1){
        n = mu -4*nParticles_;
        return 4;
      }
      if ( 5*nParticles_ <= mu ){
        n = mu -5*nParticles_;
        return 5;
      }
    }
  }
  else{
    cout << "mu =  "<< mu <<endl;
    cout << "partialSum_.size is: "<<partialSum_.size()<<endl;
    cout << "Error, in mu or partialSum_"<<endl; 
    abort();
  }
  
  return -1;
}

double Lattice::diffEffConstant(void){
  double D_eff;
  
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

    if (dim_>=2){
      for (int i=0; i < nParticles_; i++){
        sum_temp = sum_temp + 1.0/jumpRate_[i].y.r;
        sum_temp = sum_temp + 1.0/jumpRate_[i].y.l;
      }
    }

    if (dim_==3){
      for (int i=0; i < nParticles_; i++){
        sum_temp = sum_temp + 1.0/jumpRate_[i].z.r;
        sum_temp = sum_temp + 1.0/jumpRate_[i].z.l;
      }
    }

    sum_temp= (double) sum_temp / (2*nParticles_*dim_); 
    D_eff = 1.0 / sum_temp;

    //Calculate Average as well:
    D_av_ = 0;
    int startingParticle = 1;  //don't include the tagged
    for (int i = startingParticle; i < nParticles_; i++){
      D_av_ = D_av_ + jumpRate_[i].x.r;
      D_av_ = D_av_ + jumpRate_[i].x.l;
    }
    if(dim_ > 1){
      for (int i = startingParticle; i < nParticles_; i++){
        D_av_ = D_av_ + jumpRate_[i].y.r;
        D_av_ = D_av_ + jumpRate_[i].y.l;
      }
    }
    if (dim_>2){
      for (int i = startingParticle; i < nParticles_; i++){
        D_av_ = D_av_ + jumpRate_[i].z.r;
        D_av_ = D_av_ + jumpRate_[i].z.l;
      }
    }
    D_av_ = D_av_ * 1.0/(nParticles_ * 2 * dim_);


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
    D_av_  = 1.0 / (2*dim_);
    D_eff = 1.0 / (2*dim_);
    return D_eff;
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


void Lattice::saveBinning(int Ensemble){
  //This function is for getting P(t,x), ie. the probability
  //distribution of finding a particle at x, at time t.
  if(!lowMem_){
    int No_bins = 300;  //number of bins
    vector<vector<float> > PrintToFile; 
  
    //-------------------
    // TODO: cutoff vid +49? skall vara -50...+50!
    // TODO: System not in equilibrium at t=0, sample anyway?

    double Min = (double) -latticeX_/2;  //max displacement from starting position (at least for dim=1) 
    double Max = (double)  latticeX_/2;
    //double Min = *min_element(dx[e].begin(),dx[e].end());  //Gives the smallest value in dx[e]
    //double Max = *max_element(dx[e].begin(),dx[e].end());  //BUG2 use ensemble averaged min/max

    double bin_step = (double) (Max - Min)/No_bins;
    vector<float> Yaxis;     //Probability count
    vector<float> Xaxis;     //binning-axis
  
    // if(bin_step<1){
    //   cout<<"Cannot have bins smaller than 1"<<endl
    //       <<"increase lattice, or decrease bin_steps"<<endl;
    //   //Must abort, or else, the result will be very strange...

    //   abort();
    // }

    Xaxis.assign(No_bins,0);
    Yaxis.assign(No_bins,0);

    float temp = Min;
    for(int i = 0; i < No_bins; i++){
      Xaxis[i] = temp;
      temp = temp + bin_step;
    }

    double PI =3.14159265;
    double R;
    cout <<"Calculating P(t,x)... ";
    for(int k = 0; k < maxElement_; k++ ){ // OBS "k" is time! (element pos is time)
      // dx[0].size = t.size() 
      for (int j=0; j < No_bins ; j++){
        for(int e = 0; e < Ensemble; e++){
          //dr=dx[e][k];
          R = sqrt(pow(store_dx[e][k],2) + 
                   pow(store_dy[e][k],2) + pow(store_dz[e][k],2));
          //R = dx[e][k] + dy[e][k] + dz[e][k];
          if(Xaxis[j+1] > R && R >= Xaxis[j]){
            Yaxis[j] = Yaxis[j] + 1.0/(Ensemble*2*PI*((R+bin_step)/2));
          }
        }
      } 
      PrintToFile.push_back(Yaxis);
      Yaxis.clear();
      Yaxis.assign(No_bins,0);
    }

    //Print to file:
    //-------------
    ofstream hist;
    char NameOfFiles[]="histogram.dat";
    hist.open(NameOfFiles);
    cout<<" ...Printing... "<<endl;  
    for(int i = 0; i < maxElement_; i++){
      for(int j=0; j < No_bins; j++){
        hist << Xaxis[j] + bin_step/2 << "\t" << samplingTime_[i] 
             <<"\t"<< PrintToFile[i][j]<<endl;
      }
      hist <<endl; //separate data block with new (empty) line
    }
    hist.close();
    cout<<"Pinted output to: "<<NameOfFiles<<endl;
  
  }
  else
    cout<<"lowMem=ON, therefore no information stored!"
        <<endl<<"can't print cluster-information"<<endl;
}

