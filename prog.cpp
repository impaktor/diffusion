#include <iostream>        //for cout, among others
#include <fstream>         //read write to files
#include <cstdlib>         //for abort()-function
#include <cmath>           //gives us sqrt, pow, fabs, 
#include <algorithm>       //needed for min_element in hist/bin-function
#include <sstream>         //behövs kanske för "ostringstream" i save()-func.
#include <string>
#include <vector>

#include "nr3.h"           //Random number generator from numerical recipes 3 ed.
#include "ran.h"  

#include "auxilary.h"      //non-physics stuff. (print messages etc.)
#include "classes.h"       //my own classes (Jump, Particle)

//BUG! LowMem 1: k_t 1, k_c 0.5, fix bound, 400^2, N=1, t_max 2000, E=1000, values=1000
//MEN! den skriver ut perfekt output??!!

// use ssh:
//      $ ssh -n 'cd /to/path; nice -19 ./program <input.dat'
// Or:  $ ssh -X nice -19 ./program output.dat < input.dat'

//HARDCODED VARIABLES:                                             SETS THEM IN:
#define DISTRIBUTION     0    //0=uniform, 1=exponential,         | (main)
//                              2=power-law, 3=nakazato           | (main)
#define JUMPRATE_TRACER  0.01//                                   | (main)
#define JUMPRATE_CROWDER 0.5  //(if nakazato-distribution)        | (main)
#define FIXBOUNDARY      1    //1=fix wall, 0=periodic boundary   | (main ->class)
#define TEST_1           0    //boolean, prints details to screen | (constructor)
#define LOWMEM           0    //0=false / 1=true                  | ()
#define EXPONENTAIL_WAITING_TIME 0 //                             | (Move2())


//now test my InteractionRoutine, with LowMem. (and without!):
// at -i 2, c=0.7 k_t=0.25, LowMem=ON  Identical!
// at -i 3, c=0.4 k_t=0.25, LowMem=OFF Identical!
// at -i 2, c=0.5 k_t=0.25, LowMem=ON  Identical!

//now test LowMem ON (and this new merged code) for interactions off
// at k_t=0.5 c=0.5,         2Duniform,  Identical!
// at k_t=0.5 c=0.5 -i 0 (!) 2Duniform,  Identical! (got warning about not using Nakazato)

using namespace std;

class Partiklar {
  vector<Particle> pos;                    //XYZ-position for each particle on the lattice
  Particle pos_0;                          //Initial position for the tagged particle
  int nParticles_;                         //Total number of particles
  int latticeX_, latticeY_, latticeZ_;     //Lattice size

  vector<long double> x2_err,y2_err,z2_err,r2_err; //Store std error, in each point
  vector<double> x2_mu,y2_mu,z2_mu,r2_mu;  //average of displacement squared: <X^2>,... 
  vector<double> x_mu, y_mu, z_mu, r_mu;   //average of displacement: <X>, <Y>...
  vector<double> t;                        //store the time to be matched to the stored pos.

  vector<Jump> jumpRate_;                   //jump-rate for each particle and direction
  vector <long double> partialSum_;        //cumulative sum of jump-rates

  int dim_;                                //dimension of lattice (=1 || 2 || 3)
  bool boundaryFix_;                       //"true" if fix, "false" if periodic
  float concentration_;                    //particle concentration
  double timeSum_;                         //Sum of time for each move

  //TODO------------
  //vector< vector<Particle> > msdTracer;    //pos. of tracer particle. needed for stdrerr
  vector< vector<int> > dx,dy,dz;            //pos. of tracer particle. needed for stdrerr
  vector< vector<double> > dr;               //pos. of tracer particle. needed for stdrerr
  //------------

  int convertMuToParticle(int, int&);
  void computePartialSum(void);
  void moveAndBoundaryCheck(int,int);

  int vacancyCheckOld(int, Particle);     //No double Occupancy.
  int vacancyCheck(int, Particle);        //Same as above, but improved. Use this one instead.
  vector<vector<vector<int> > > vacancy_; //Needed in my more efficient "vacancyCheck()"


  void stdErr(int);                        //calculate standard deviation
  void stdErrLowMem(int);                  //a more memmory conservative version of the above
  bool LowMem_;                             //Switch between the two stdErr-functions above. 
  vector<double> dr4_err,dx4_err,dy4_err,dz4_err; //needed if we use stdErrLowMem()

  float jumpCrowders_, jumpTracer_;        //only used in the Nakazato()-function
  int distribution_;
  double info_;                             //save characteristic trait of particle distribution
  bool testOnOff_;                          //To print detailed information to screen
  double D_av_; 

  bool interactionOn_;                      //use the interaction algorithm
  float interactionStrength_;               //Interaction strength to use, if "InteractionOn=true"

  float nakazato(void);         
  double diffEffConstant(void);            //Calculate the effective dif. const. from the jump-rates
  float ergodicity(int,int,int);           //Calculate the theoretical MSD when ergodic

  int round(float x) {                     //Round off (used to get center of lattice)
    return int(x + 0.5);  }

  //INTERACTION CODE ***
  void interaction(int, Particle);
  void countNeighbours(Particle, vector<int>&);
  vector<float> clusterDistribution, clusterSize;  //For our save_cluster()-function
  void buildCluster2(int, vector<int>&, double);     //needed for save_cluster()-func.
  void checkVacancyMatrix(float);                    //Just a test-function

public:
  Partiklar(int, int, int, int, bool);
  void place(void);                        //Place particles on the lattice
  void moveOld(void);                      //old obsolete move-code (no jump-rates)
  void move(void);                        //pick a particle & direction based on jump-rate
  void store(vector<int>&, vector<int>&, vector<int>&, vector<double>&);
  void save(int, char[]);                  //Print simulation-run to file   
  void set_dt(void);                       //needed to initiate time-vector
  double get_dR(void);
  double get_time(void);
  int get_dX(void);
  int get_dY(void);
  int get_dZ(void);
  int getDimension(void);
  void set_jumpNaka(float, float);
  void set_dist(int, double);
  void set_interaction(float);             //determine if we use the interaction-algorithm
  void set_LowMem(bool, int);              //Set the switch true/false for use of low-mem run

  //Set jump-rate of each particle and direction 
  void setJumpRate(vector<Jump>);

  //Secondary not used functions
  void snapshot(void);                     //Not used. (more comments in function body)
  //  double get_dt(void);   //Size of time-step, need it for P(x,t), don't need this...


  //INTERACTION CODE:
  void saveBinning(int);                     //To get P(x,t)
  void saveCluster(double);                  //To get rho(m) (cluster distribution)
  void printCluster(void);                   //Print  rho(m) (cluster distribution)

};


void Partiklar::saveCluster(double Time){
  //This funciton calculates how many particles are in 
  //the cluster that the tagged particle is part of

  //ha denna funktionen sa att den kors varje tidssteg for sista ensemblen,
  //Den skall veta nar den ar pa sista tidssteget, tror jag...

  //Cluster_Distribution defined in the constructor.
  vector<int> Cluster;
  //TODO kolla Boundary! 0 eller 1?

  //BUG augusti, detta stammer val inte?
  //ne, men jag har kommenterat ut hela den saken, sa vardet pa 
  //Exponent spellar ingen roll!
  double Exponent = 0; //ie. rand > 0 allways include the neighbor.  

  // checkVacancyMatrix(1.01); //OK hit

  if(Time>20){  //The system is not in equilibrium from start
    //checkVacancyMatrix(10.1);
    countNeighbours(pos[0], Cluster);
    //checkVacancyMatrix(10.2);
    buildCluster2(0, Cluster, Exponent); //Skriver denna till nagot?
    //checkVacancyMatrix(10.3); //OK hit
    int M = Cluster.size();
    //First element is # of 1-clusters (ie. single particles), etc.
    clusterDistribution[M-1] = clusterDistribution[M-1] + 1.0/nParticles_;
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

void Partiklar::printCluster(void){
  //This function prints the information gathered by 
  //the saveCluster()-function, and prints it to file

  //  CheckVacancyMatrix(9);

  if (clusterSize.size() != clusterDistribution.size())
    printError("gaaa");

  ofstream clust;
  char NameOfFiles[]="cluster.dat";
  clust.open(NameOfFiles);
  for(unsigned int j=0; j < clusterDistribution.size(); j++){
    clust << clusterSize[j] <<"\t"<< clusterDistribution[j]<<endl;
  }
  clust.close();
  cout<<"Pinted "<<NameOfFiles<<endl;
}


void Partiklar::saveBinning(int Ensemble){
  //This function is for getting P(t,x), ie. the probability
  //distribution of finding a particle at x, at time t.
  if(!LowMem_){
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
    for(unsigned int k = 0; k < t.size() ; k++ ){ // OBS "k" is time! (element pos is time)
      // dx[0].size = t.size() 
      for (int j=0; j < No_bins ; j++){
        for(int e = 0; e < Ensemble; e++){
          //dr=dx[e][k];
          R = sqrt(pow(dx[e][k],2) + pow(dy[e][k],2) + pow(dz[e][k],2));
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
    for(unsigned int i=0; i <t.size(); i++){
      for(int j=0; j < No_bins; j++){
        hist << Xaxis[j] + bin_step/2 <<"\t"<< t[i] <<"\t"<< PrintToFile[i][j]<<endl;
      }
      hist <<endl; //separate data block with new (empty) line
    }
    hist.close();
    cout<<"Pinted output to: "<<NameOfFiles<<endl;
  
  }
  else
    cout<<"LowMem=ON, therefore no information stored!"
        <<endl<<"can't print cluster-information"<<endl;
}



//IMPORTED FROM TESTInteractionCode() 12/7 -2010
//NOTE! BUT! "Exponent" changed to allow 100% probability of being included.
void Partiklar::buildCluster2(int n, vector<int>& NearestNeighbours, 
                              double Exponent){
  vector<int> unique;
  vector<int> NewParticleAdded;
  static Ran RandomNUM(42);
  double randomNumber;

  //Build a vector with particles that form the Cluster:
  //-------------------------------------------------
  //each neighbor has a probability given by the Boltzman factor:
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
      countNeighbours(pos[m], NearestNeighbours);
      //NearestNeighbours-vector will store more particles each 
      //iteration, but only loop over the first initial set (ie. i < limit)
    }

    for(unsigned int j=0; j < NearestNeighbours.size(); j++){
      bool SingleValued=true;
 
      //Check if it's stored in unique[] allready:
      for(unsigned int k=0; k < unique.size(); k++){
        if (NearestNeighbours[j] == unique[k])
          SingleValued=false; 
      }
      if(SingleValued){
        //Add the new particle, with Boltzmans consent...
        // randomNumber = RandomNUM.doub();
        //if (randomNumber >=  exp(-1.0*Exponent)){
        unique.push_back(NearestNeighbours[j]);
        NewParticleAdded.push_back(NearestNeighbours[j]);
        //New particels added to the cluster, run the loop again
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


  if(false){ //Inaktiverad just nu
    //TEST:
    cout << "kollar Nearest Neighbour.size(): "<<NearestNeighbours.size()<<endl;
    for(unsigned int j=0; j < NearestNeighbours.size(); j++){
      int m = NearestNeighbours[j];
      cout <<"m = "<<m<<" ("<<pos[m].x<<","<<pos[m].y<<","<<pos[m].z<<")"<<endl;    
    }
  }

  //TEST
  for (int i=0; i < nParticles_; i++){
    if (pos[i].x==0 || pos[i].y==0 || pos[i].z==0)
      printError("Outside of lattice");
  }

  //TEST---------(FÖR ATT KOLLA DOUBLE COUNTING)
  bool IsUniqueTrue=true;
  for(unsigned int i=0; i<unique.size();i++){
    for(unsigned int j=0; j<unique.size();j++){
      if (unique[i]==unique[j] && i != j)
        IsUniqueTrue=false;
    }
  }
  if(!IsUniqueTrue)
    printError("Unique is not unique! Error *235897*");
  //-------------
}



void Partiklar::checkVacancyMatrix(float ErrorCode){
  //This is just a small function to check that the vacancy matrix is correct
  //ie. vacant sites are marked as "-1" and occupied are marked with the 
  //correct particle label.
  //ErrorCode = code to output to easier find the place in the code

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
          if(i != pos[label].x || j != pos[label].y || k != pos[label].z){
            error = true;
          }
        }
        else
          count_vacant_sites++;
      }
    }
  }
  //Easier to find exactley where (above/below) error = true
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


// double Partiklar::get_dt(void){
//   //Needed for the newly added binning-function
//   //to produce a histogram of P(x,t)
//   //(and for the snapshot-function, see "main()")
//   return partialSum_.back();
// }

void Partiklar::snapshot(void){
  //NOTE! This is only for 2D, with isotropic jumprate,
  //and the same system sidelength!
  //Its purpose is to print a snapshot of the system, 
  //to see the spatial distribution of different particles
  //as a heat map. Only used this function briefly. 

  vector< vector<double> > matris;

  vector<double> temp;
  temp.assign(latticeY_,0);
  matris.assign(latticeX_,temp);

  for(int i=0; i < nParticles_; i++){
    matris[ pos[i].x-1 ][ pos[i].y-1 ] = jumpRate_[i].x.r; //TODO BUG!
  }

  ofstream heat;
  char NameOfFiles[]="heatmap.dat";
  heat.open(NameOfFiles);
  
  bool points=true; //true: plotar endast partiklarna, false: plottar aven tomma rutor (=0)
  if(points){
    //This does not give 0 where there are no particles:
    for (int i=0; i < nParticles_; i++)
      heat  << pos[i].x<<" \t"<< pos[i].y <<"\t"<< jumpRate_[i].x.r<< endl; //TODO BUG
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

void Partiklar::set_dt(void){
  //This will store the first time value:
  if ( t.empty() ){
    t.push_back(timeSum_);
  }
  //this will store the remaining values of the first enasmble
  else{
    if ( t.back() < timeSum_ ){
      t.push_back(timeSum_);
      // cout <<"tidsfunktionen: "<<t.back()<< " timeSum_: "<<timeSum_<<endl;
    }
  }
} //NOTE funkar inte om det är olika jumprates för olika ensambles

int Partiklar::getDimension(void){ //could be useful to have, but not neccesary...
  return dim_;
}

int Partiklar::get_dX(void){
  //current position minus starting position.
  return pos[0].x-pos_0.x;
}

int Partiklar::get_dY(void){
  return pos[0].y-pos_0.y;
}

int Partiklar::get_dZ(void){
  return pos[0].z-pos_0.z;
}

double Partiklar::get_dR(void){
  double r =sqrt( pow(pos[0].x-pos_0.x,2) + pow(pos[0].y-pos_0.y,2) + pow(pos[0].z-pos_0.z,2) ); 
  return r;
}

double Partiklar::get_time(void){
  return timeSum_;
}

void Partiklar::set_jumpNaka(float jumpCrowders, float jumpTracer){
  jumpCrowders_ = jumpCrowders;
  jumpTracer_ = jumpTracer;
}

void Partiklar::set_dist(int distribution, double info){
  distribution_ = distribution;
  info_ = info;
}

void Partiklar::set_interaction(float interactionStrength){
  interactionStrength_ = interactionStrength; 
  interactionOn_ = true; 
}

void Partiklar::set_LowMem(bool LowMem, int NumberOfPoints){
  //Set the switch true/false for use of low-mem run
  LowMem_ = LowMem;

  //initiate the vectors for the summation in store()
  if(LowMem){
    dr4_err.assign(NumberOfPoints,0);
    dx4_err.assign(NumberOfPoints,0);
    dy4_err.assign(NumberOfPoints,0);
    dz4_err.assign(NumberOfPoints,0);
    
    r2_mu.assign(NumberOfPoints,0);
    x2_mu.assign(NumberOfPoints,0);
    y2_mu.assign(NumberOfPoints,0);
    z2_mu.assign(NumberOfPoints,0);
  }

}


float Partiklar::nakazato(void){
  
  float alpha;

  if (dim_ == 1) alpha = 1;
  if (dim_ == 2) alpha = 0.363;
  if (dim_ == 3) alpha = 0.209;

  //jumpratio for the different particles... 
  //OBS 1 svart, resten vita!
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


////////////////////////////
//       CONSTRUCTOR      //
////////////////////////////

Partiklar::Partiklar(int xSquare,int ySquare,int zSquare, 
                     int particleNumber, bool boundary){

  if (xSquare >= ySquare && ySquare  >= zSquare && zSquare >= 1){
    if (particleNumber <= xSquare*ySquare*zSquare){

      latticeX_ = xSquare;            //number of lattice sites
      latticeY_ = ySquare;
      latticeZ_ = zSquare;

      nParticles_ = particleNumber;

      boundaryFix_ = boundary;

      pos_0.x = round(1.0*latticeX_/2); //center of lattice, where we place our tagged particle
      pos_0.y = round(1.0*latticeY_/2); 
      pos_0.z = round(1.0*latticeZ_/2); 
      
      //TODO: what if Z = 0
      concentration_=(float) nParticles_ / (latticeX_ * latticeY_ * latticeZ_); 

      //Hardcoded Switches:
      //------------------
      //boolean, for print details to screen. 
      testOnOff_=TEST_1;

      //Define the dimention of the lattice
      if (latticeY_<=1 && latticeZ_<=1) dim_ = 1; 
      if (latticeY_>1 && latticeZ_<=1)  dim_ = 2;
      if (latticeY_>1 && latticeZ_>1)   dim_ = 3;

      //By default: don't use the interaction algorithm. 
      interactionStrength_ = 0;
      interactionOn_ = false;   

      //By default: don't use the memmory conserving stdErrLowMem()
      //use the set_LowMem()-function in the main()-func to change this
      LowMem_ = false;      

      //--------------------------
      //To get the Cluster Distribution (see saveCluster()):
      // Cluster_Distribution = number of clusters with "m" particles
      // Cluster_Size = number of particles in cluster (just our x-axis)
      int Maximus_Clustirus = nParticles_; 
      //Maximal cluster size to include in distribution
      clusterDistribution.assign(Maximus_Clustirus,0); //assign zeros...
      for (int i=1; i <= Maximus_Clustirus ; i++){
        clusterSize.push_back(i);
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


float Partiklar::ergodicity(int x_rutor, int y_rutor, int z_rutor){ 
  //check what the MSD should be when the system reaches equilibrium.
  //(only correct if N=1). This function has two independent parts. 
  //One where we just "step through" the lattice and add the 
  //displacement from the center, and another were we actually 
  //calculate it analytically. 
  
  float summa = 0; 
  
  //use these to "walk" through the lattice to meassure the MSD
  //contribution from each square
  int x,y,z;
                   

  //STEP THROUGH THE LATTICE
  //I re-map the starting coordinates, in case I'd 
  //like to move this function out from this class...
  int mittX = pos_0.x;
  int mittY = pos_0.y;
  int mittZ = pos_0.z;

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


void Partiklar::place(void){
  //This function will be called once each ensemble to
  //place the N particles (randomly), and reset the time
  
  int totalSites = latticeX_ * latticeY_ * latticeZ_;
  int leftToPlace = nParticles_; 

  //assign zeros
  Particle nil(0,0,0);                 
  pos.assign(nParticles_,nil); 
                                   
  timeSum_ = 0;                   //Time starts at zero each run
                                   
  pos[0] = pos_0;                 //place tracer particle at center
  leftToPlace--;                  //one less to place
  totalSites--;                    
  
  double R;
  static Ran randomNumber1(1);    //just any seed will do.

  int i, j, k, n = 1;             //First crowder at n = 1 element   
  for (k=1; k <= latticeZ_; k++){    //place particles
    for (j=1; j <= latticeY_; j++){
      for (i=1; i <= latticeX_; i++){

        //don't place crowders on the tagged particle
        if (i != pos_0.x || j != pos_0.y || k != pos_0.z){ 

          //Using NR to generate 0 < R < 1 ;
          R=randomNumber1.doub();              
          if (R<=leftToPlace*1.0/totalSites){
            pos[n].x = i;
            pos[n].y = j;
            pos[n].z = k;
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
    for (int q=0; q < nParticles_; q++)
      cout <<" PLACED \n particle n = " << q << pos[q] << endl;


  //IF we use the new vacancyCheck()-function:
  //-----------------------------------
  //initiate the vacancy-vector with "-1" (meaning vacant) 
  //Must reset this matrix for each ensemble.
  vector< vector< vector<int> > >
    Vacancy(latticeX_+1, vector< vector<int> >
            (latticeY_+1, vector<int>(latticeZ_+1,-1)));
  vacancy_=Vacancy;
  //Mark occupied sites with particle label for that site:
  //note, no change will ever occur to the 0-elements of the matrix, 
  //since out lattice starts at 1,..,latticeX_, but matrix is 0,..,latticeX_
  for (n=0; n < nParticles_; n++){
    vacancy_[ pos[n].x ][ pos[n].y ][ pos[n].z ]=n;
  }
}



void Partiklar::moveOld(){ 
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
      

  //Do the actual move, acording to the boundary rules
  moveAndBoundaryCheck(n,r);

  timeSum_ = timeSum_ + 1.0 / nParticles_;  //Calculate the time

  if (testOnOff_)
    for (int i=0; i < nParticles_; i++)
      cout <<" MOVED \n particle n = " << i <<": "<< pos[i] << endl;
}


void Partiklar::moveAndBoundaryCheck(int n, int R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant. 
  bool boundary = boundaryFix_;
  
  if (0 <= n && n < nParticles_ ){
    if (testOnOff_) cout <<"Moving particle "<< n <<endl;
    // n = particle to move, R = direction to move, such that: 
    //      R= 0 = X-right,      1 = X-left
    //         2 = Y-right,      3 = Y-left
    //         4 = Z-up,         5 = Z-down
    // bondary: boundary=0 = periodic, boundary=1 = fix boundary

    //Store in case we want to restore the previous configuration             
    Particle old = pos[n];
    
    if ( R >= 0 && R < 6 ){
      switch(R){                                 
      case 0: pos[n].x=pos[n].x +1;
        if (pos[n].x == latticeX_ +1 && boundary == 0) pos[n].x = 1;
        if (pos[n].x == latticeX_ +1 && boundary == 1) pos[n].x = old.x;
        break;
      case 1: pos[n].x = pos[n].x -1; //move particle n to the left
        if (pos[n].x == 0 && boundary == 0) pos[n].x = latticeX_;//periodic
        if (pos[n].x == 0 && boundary == 1) pos[n].x = old.x;//fix 
        break;  //(first coordinate is 1, 0=outside) 
      case 2: pos[n].y = pos[n].y +1;
        if (pos[n].y == latticeY_ +1 && boundary == 0) pos[n].y = 1;
        if (pos[n].y == latticeY_ +1 && boundary == 1) pos[n].y = old.y;
        break;
      case 3: pos[n].y = pos[n].y -1;
        if (pos[n].y == 0 && boundary == 0) pos[n].y = latticeY_;
        if (pos[n].y == 0 && boundary == 1) pos[n].y = old.y;
        break;
      case 4: pos[n].z = pos[n].z +1;
        if (pos[n].z == latticeZ_+1 && boundary == 0) pos[n].z = 1;
        if (pos[n].z == latticeZ_+1 && boundary == 1) pos[n].z = old.z;
        break;
      case 5: pos[n].z = pos[n].z-1;
        if (pos[n].z == 0 && boundary == 0) pos[n].z = latticeZ_;
        if (pos[n].z == 0 && boundary == 1) pos[n].z = old.z;
        break;
      }

      //check that the new site is vacant if not,
      //move the particle back. (done by the vacancyCheck...)
      
      //vacancyCheckOld(n,temp);  //Old code (slow)                                        
      vacancyCheck(n,old);        //New improved, faster version.

    }
    else{
      cout <<"There is a maximum of "<<2*dim_
           <<" positive directions to move in"<<endl;
      abort();
    }
  }
  else{
    cout <<"Invalid particle label, n="<<n
         <<", must fulfill 0 <= n < N"<<endl;
    abort();
  }
}



int Partiklar::vacancyCheckOld(int n, Particle previous){
  //Old and obsolete version. It doesn't need a 3D matrix as the
  //new version does, vacancyCheck(...) but this one is slow.
  //Keep it here to check that the output is exactly identical
  //(compare output-files with the "diff" command)
  
  if (n < nParticles_ && 0 <= n){

    int returnvalue = 0;    //unused parameter to register collisions 
    for (int i=0; i<nParticles_; i++){
      if (pos[i].x == pos[n].x && pos[i].y == pos[n].y && 
          pos[i].z == pos[n].z && i!=n){ 
        if (testOnOff_) cout << "Occupied: " << pos[n] << endl;

        //if two particles on same square: undo move
        pos[n] = previous;                                      
        returnvalue = 1;
        return returnvalue;     //End function here
      }
      else
        if (testOnOff_) cout <<"OK"<<endl;
    }
    if (testOnOff_) cout <<"Vacant, particle "<< n<<": "<<pos[n]<<endl;
    return returnvalue;
    //(the returnvalue could be used to have a 
    //    do-move-while(vacancyCheck == 1) 
    //but if so, we need to move this function to be directly 
    //in the Move_()-func.)
  }
  else{
    printError("Error, accessing invalid particle");
  }
}



int Partiklar::vacancyCheck(int n, Particle oldPos){
  //NOTE: This is a faster version than vacancyCheckOld(), since that 
  //one loops through all particles for each "vacancy check". However 
  //this uses a 3D matrix, and there will be some "double banking" 
  //in the sence that the postion is now stored in this matrix and 
  //in XYZ-vectors.

  if (n < nParticles_ && 0 <= n){

    int returnvalue = 0;   //unused parameter to register collisions
    if (vacancy_[pos[n].x][pos[n].y][pos[n].z] != -1){  //if occupied
      if (testOnOff_) cout <<"Occupied: " << pos[n] << endl;
      pos[n] = oldPos;                              //then move back          
      returnvalue = 1;

      if (testOnOff_){
        for (int k=1; k <= latticeZ_; k++){
          //print y in reverse direction, (better in the terminal)
          for (int j=latticeY_ ; j > 0 ; j--){
            cout <<j<<": ";
            for (int i=1; i <= latticeX_; i++){
              cout <<vacancy_[i][j][k]<<" ";
            }
            cout << endl;
          }
          cout <<"looping for new Z"<<endl;
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
        vacancy_[pos[n].x][pos[n].y][pos[n].z] = n; 
        
        if (testOnOff_) cout <<"OK, vacant"<<endl;
      }
      //if we use the interaction algorithm:
      else
        interaction(n,oldPos);
    }
    if (testOnOff_) cout <<"Vacant, particle "<<n<<": "<<pos[n]<<endl;
    return returnvalue;
    //(the return value could be used to have a 
    //    do-move-while(vacancyCheck == 1)
    //but if so, we need to move this function to be directly in the
    //Move_()-func.)
  }
  else{
    printError("Error, accessing invalid particle");
  }
}



void Partiklar::move(){
  
  //pick a particle & direction based on jump-rate, and move.
  if (partialSum_.empty())
    printError("Partial-sum vector has not been initiated");
  
  static Ran randomNumb(8);
  double r2;

  //  for (int i=0; i < N; i++){      //make N moves

  int mu_guess;                       // must be integer. (index of vector)
  int mu_left = 0;
  int mu_right = dim_*2*nParticles_;
  double p_left = (double) partialSum_.front();  //value in first element ( =0 )
  double p_right = (double) partialSum_.back();  //value in last element ( =dim_*2*N)
 
  do{
    r2 = randomNumb.doub();
  }while( r2 == 1 || r2 == 0 );

  double p_rand = (double) r2*partialSum_.back();
  bool LoopAgain=true;

  //NOTE: when comparing a long double against a double, C++ will fill
  //in the blanks with zeroes. Every 10^15 run of move() will lead to
  //a situation where we actually have p_rand==p[mu_guess], and
  //changing the intervall below in the if statements,
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

  bool expWaitTime = EXPONENTAIL_WAITING_TIME;
  static Ran randomNumb2(87);

  //XXXX
  //Use linear or exponential waiting time
  if(expWaitTime){
    double r3;
    do{
      r3 = randomNumb2.doub();
    }while( r3 == 1 || r3 == 0 );
    
    timeSum_ = timeSum_ + log(1.0/r3)/partialSum_.back();
  }
  else{
    timeSum_ = timeSum_ + 1.0/partialSum_.back();
  }

  //} //move N times end
}

void Partiklar::interaction(int n, Particle oldPos){
  vector<int> NearestNeighbours;
  //vector contains the particle labels of the nearest neighbors to
  //particle n (not sure that we actually need this at this point)

  //This is an updated version of the progInteraction.cpp I used in my
  //thesis It does the exact same thing, but now we can wrap around 
  //the boundary also...
  //The out-putfiles are identical between the two different version,
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
    vacancy_[oldPos.x][oldPos.y][oldPos.z] = -1;//old pos. now vacant
    vacancy_[pos[n].x][pos[n].y][pos[n].z] = n; //new pos. now occupied
  }
  else{
    //reject move, move back.
    //vacancy_-matrix unchanged
    pos[n]=oldPos;
  }  
}



void Partiklar::countNeighbours(Particle particle, vector<int>& Count){
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



void Partiklar::store(vector<int>& x_displace, vector<int>& y_displace,
                      vector<int>& z_displace, vector<double>& r_displace){
  //This is needed to calculate error bars (stdErr())
  
  /*
    Error in each point i is (in LaTeX notation):
    (a[i])^2 = 1/N \sum_i^N (a_i - <a>)^2 = [1/N \sum_i^N a[i]^2] - [<a>]^2

    <a>= 1/N \sum_i^N a[i]   (i.e. average value in that point) 
    N = Number of ensembles

    The original way to do this is save the value in each point a[i] for 
    each ensemble, but that leaves us with a giant matrix if we have many
    ensembles N, and many data points D, (N*D ~> 10^9) ~ 1Gb of mem usage
    (10^9/2^30 ~ 1Gb). Therefore with the switch 'LowMem=true' we do the 
    left hand side of the equation above, and then we don't need a matrix
    but just two vectors, for each dimension. 
  */

  //Check that we assigned the vectors =0 if we uses the LowMem-option.
  //For LowMem=false, we do this in the Save & stdErr()-functions
  //since we use push_back() instead.
  if(LowMem_){
    bool error=false;
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
          <<endl<<"Must be assigned by set_LowMem()-function!"<<endl;
      abort();
    }
  }

  if(LowMem_ == false){
    //if don't use the LowMem-algorithm:
    //this gives dx[ensembler][number of values] 
    dx.push_back(x_displace);   
    dy.push_back(y_displace);
    dz.push_back(z_displace);
    dr.push_back(r_displace);
  }
  else{
    //To conserve memmory we do the summation directley, which might
    //cause a nummerical error in the end compared to doing it the 
    //more rigorous way. Compare to make sure it works.
    for(int i=0; i < r_displace.size(); i++){
      dr4_err[i] = dr4_err[i] + pow( r_displace[i] ,4);
      dx4_err[i] = dx4_err[i] + pow( x_displace[i] ,4);
      dy4_err[i] = dy4_err[i] + pow( y_displace[i] ,4);
      dz4_err[i] = dz4_err[i] + pow( z_displace[i] ,4);

      r2_mu[i] = r2_mu[i] + pow( r_displace[i] ,2);
      x2_mu[i] = x2_mu[i] + pow( x_displace[i] ,2);
      y2_mu[i] = y2_mu[i] + pow( y_displace[i] ,2);
      z2_mu[i] = z2_mu[i] + pow( z_displace[i] ,2);
    }
  }
  
}


void Partiklar::save(int Ensemble, char name[]){
  //Now process the simulation data, and print it to
  //file "name".
  
  char* Primary = name;

  int numberOfValues = t.size();
  ofstream primary,secondary;
  primary.open(Primary);

  //hardcoded switch to print x,y,z, x^2,y^2,z^2
  //specifically. Not needed now, since system is isotropic,
  bool xyz=false;
  if(xyz) secondary.open("txyz.dat");

  if(!LowMem_){
    //Compute <R>:
    x_mu.assign(numberOfValues,0);
    y_mu.assign(numberOfValues,0);
    z_mu.assign(numberOfValues,0);
    r_mu.assign(numberOfValues,0);

    for (int j=0; j<numberOfValues; j++){  
      for (int i=0; i< Ensemble ; i++){
        x_mu[j]=x_mu[j] + 1.0*dx[i][j]/Ensemble;
        y_mu[j]=y_mu[j] + 1.0*dy[i][j]/Ensemble;
        z_mu[j]=z_mu[j] + 1.0*dz[i][j]/Ensemble;
        //  r_mu[j]=r_mu[j] + 1.0*dr[i][j]/Ensemble;  //I'm not using this.
      } 
    }

    //Compute <R^2>:
    x2_mu.assign(numberOfValues,0);
    y2_mu.assign(numberOfValues,0);
    z2_mu.assign(numberOfValues,0);
    r2_mu.assign(numberOfValues,0);

    for (int j=0; j<numberOfValues; j++){  
      for (int i=0; i<Ensemble ;i++){
        x2_mu[j]=x2_mu[j] + 1.0*pow(dx[i][j],2)/Ensemble;
        y2_mu[j]=y2_mu[j] + 1.0*pow(dy[i][j],2)/Ensemble;
        z2_mu[j]=z2_mu[j] + 1.0*pow(dz[i][j],2)/Ensemble;
        r2_mu[j]=r2_mu[j] + 1.0*pow(dr[i][j],2)/Ensemble;
      }
    }

    //Calculate standard error, to generate error-bars
    //this function must come after x2_mu...r2_mu has been
    //computed.
    stdErr(Ensemble);
  }
  else{
    //use the memmory conserving version instead. This one does not
    //use the dx[][]-matrix, and in that way we save lots of mem.
    stdErrLowMem(Ensemble);
  }



  //To print which distribution we used to head of output-file
  string dist[]={"uniform","exponential","powerlaw","nakazato"};
  string onOff[]={"Off","On"};
  string bound[]={"periodic","fix"};
  float d_eff= diffEffConstant();  //calculate effective diffusion constant
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


  //print data to file:
  for (int i=0; i < numberOfValues; i++){

    //Maybe add this:
    //cout << setprecision (5) << r2_mu[i] << ...

    //Print it to file:
    primary  << t[i]<<"\t"<<r2_mu[i]<<"\t" <<r2_err[i]<<endl;  //skiter i "r"
    if(xyz){
      secondary<< t[i]<<"\t"<< x_mu[i] <<"\t"<< y_mu[i] <<"\t"<< z_mu[i] <<"\t"
               << x2_mu[i] <<"\t"<< y2_mu[i] <<"\t"<< z2_mu[i] <<"\t"
               << x2_err[i]<<"\t"<<y2_err[i]<<"\t"<<z2_err[i]<<endl;
    }
  }
  primary.close();
  if(xyz) secondary.close();
}


void Partiklar::stdErr(int Ensemble){
  //Calculate the standard deviation, to get error-bars.
  int numberOfValues = t.size();

  //Should not be using this function if LowMem=true!
  if(LowMem_)
    printError("LowMem switch ON, but using wrong stdErr-function!");

  x2_err.assign(numberOfValues,0);
  y2_err.assign(numberOfValues,0);
  z2_err.assign(numberOfValues,0);
  r2_err.assign(numberOfValues,0);


  for (int j=0; j < numberOfValues  ;j++){
    for (int i=0; i<Ensemble ; i++){
      x2_err[j]=x2_err[j] + pow( pow(dx[i][j],2) - x2_mu[j] ,2);
      y2_err[j]=y2_err[j] + pow( pow(dy[i][j],2) - y2_mu[j] ,2);
      z2_err[j]=z2_err[j] + pow( pow(dz[i][j],2) - z2_mu[j] ,2);
      r2_err[j]=r2_err[j] + pow( pow(dr[i][j],2) - r2_mu[j] ,2);
    }
  }

  for (int j=0; j < numberOfValues  ;j++){
    x2_err[j]=sqrt( x2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
    y2_err[j]=sqrt( y2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
    z2_err[j]=sqrt( z2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
    r2_err[j]=sqrt( r2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
  }     //Då E=1 --> inf !

  
  //Testing the Gaussian distribution of the points:
  //--------------------------------------------------------------------
  float SUM68=0;
  for (int j=0; j < numberOfValues; j++){    //forgot about his: 
    double top = x2_mu[j] +0.5*x2_err[j] *   sqrt(Ensemble-1)*2;
    double bot = x2_mu[j] -0.5*x2_err[j] *   sqrt(Ensemble-1)*2;
    for (int i=0; i < Ensemble; i++){
      if (bot <= pow(dx[i][j],2) &&  top >= pow(dx[i][j],2) ){
        bool test = false;
        if (test){ //TEST-function
          cout <<"lower limit: "<< bot <<endl;
          cout << "point is: "<< pow(dx[i][j],2) <<endl;
          cout << "upper limit: "<< top <<endl;
        }
        SUM68++;
        if (test) cout <<"SUM68 = "<<SUM68<<endl;
      }
    }
  }
  cout <<"68,2 % (?)="<<1.0*SUM68/(Ensemble*numberOfValues)<<endl;

  //--------------------------------------------------------------------
}


void Partiklar::stdErrLowMem(int Ensemble){
  int numberOfValues=t.size();

  if(!LowMem_)
    printError("LowMem=OFF, but using stdErrLowMem()-func!");

  //Devide the summed up displacement to get the average:
  for (int i=0; i< numberOfValues; i++){
    x2_mu[i]=x2_mu[i]/(Ensemble*1.0);
    y2_mu[i]=y2_mu[i]/(Ensemble*1.0);
    z2_mu[i]=z2_mu[i]/(Ensemble*1.0);
    r2_mu[i]=r2_mu[i]/(Ensemble*1.0);
  }

  x2_err.assign(numberOfValues,0);
  y2_err.assign(numberOfValues,0);
  z2_err.assign(numberOfValues,0);
  r2_err.assign(numberOfValues,0);

  //Following the notation in the comments in Partiklar::store():
  // [1/N sum_i^N a_i^2] - [<a>^2], (a=dr^2)
  for(int i =0; i < numberOfValues; i++){
    x2_err[i] = dx4_err[i]/(1.0*Ensemble) - pow( x2_mu[i] ,2);
    y2_err[i] = dy4_err[i]/(1.0*Ensemble) - pow( y2_mu[i] ,2);
    z2_err[i] = dz4_err[i]/(1.0*Ensemble) - pow( z2_mu[i] ,2);
    r2_err[i] = dr4_err[i]/(1.0*Ensemble) - pow( r2_mu[i] ,2);
  }

  for(int i =0; i < numberOfValues; i++){
    x2_err[i] = sqrt(x2_err[i]/(1.0*(Ensemble-1)));
    y2_err[i] = sqrt(y2_err[i]/(1.0*(Ensemble-1)));
    z2_err[i] = sqrt(z2_err[i]/(1.0*(Ensemble-1)));
    r2_err[i] = sqrt(r2_err[i]/(1.0*(Ensemble-1)));
  }
}



//////////////////////////////////////
//   ****************************   //
//   * M A I N  F U N C T I O N *   //
//   ****************************   //
//////////////////////////////////////

int main(int argc, char* argv[]){
  int N;                              //Number of particles
  int X,Y,Z;                          //lattice size
  int ensemble,antal_punkter;
  double MaxTime;

  //DEFAULT VALUES:  
  bool fixBoundaryOn = FIXBOUNDARY;     //boundary condition, false = periodic
  bool logscale = false;
  bool nymetod = true;                  //use move() instead of moveOld()
  bool UseLowMem = LOWMEM;              //TODO! Patch this to the -m flag!
  bool quiet = false;                   //print remaining ensembles during simulation
  //--------------

  char NameOfFile[99] = "out.dat";

  int InteractStr = 0;
  bool InteractOn = false;

  //change default values depending on arg. flags given
  argumentFlags(argc, argv, logscale, InteractOn,
                quiet, InteractStr, NameOfFile);

  //TEST print result to screen:
  //-----------------------
  string OnOff[]={"OFF","ON"};
  int tempIndex=0;
  if (InteractOn) tempIndex=1;
  cout<<"###log = "<<logscale<<endl
      <<"###interaction strength "<<OnOff[tempIndex]<<" ="<<InteractStr<<endl
      <<"###filename = "<<(string) NameOfFile<<endl;
  //-----------------------


  if (fixBoundaryOn) 
    cout <<"# Fix boundary"<<endl;  
  else 
    cout <<"# Periodic boundary"<<endl;


  AskUserForInputParameters(X,Y,Z,N,ensemble,antal_punkter,MaxTime);


  Partiklar crowd(X,Y,Z,N,fixBoundaryOn);
  
  
  if(InteractOn)                           //use Interaction algorithm
    crowd.set_interaction(InteractStr);    //with "InteractStr" strength.


  //kanske ha den som invariabel i Save() istallet?
  //Set the switch true/false for use of low-mem run
  crowd.set_LowMem(UseLowMem, antal_punkter);

  static Ran tempnumb(17);                 //just anny seed will do
  vector<Jump> hopRate;
  if (nymetod){
    for (int particle = 0; particle < N; particle ++){
      double u;
      do { 
        u = tempnumb.doub();   
      }while (u == 0);


      int n = DISTRIBUTION;  //Choose distribution
      // 0 = uniform
      // 1 = exponential
      // 2 = power law
      // 3 = nakazato

      //Characteristic trait of the chosen distribution, 
      //just used to print info to file/screen
      double Info=0;

      float jumpTracer;
      
      //NOTE: The jumprate is in EACH direction! Meaning, jumprate =1 is actually
      //a jumprate =4 in 2D. (I think)
      if (0 <= n && n<=3 && (u < 1 && u > 0) ){
        jumpTracer = JUMPRATE_TRACER; //Jumprate for the tracer paricle!
        float lambda = 1.0;              //used in option 1 Exp-dist.
        float y_c = 1.0;               //used in option 2 Power-law
        float alpha = 0.5;               //used in option 2 Power-law
        float jumpCrowders = JUMPRATE_CROWDER; //used in option 3 Nakazato
        crowd.set_jumpNaka(jumpCrowders, jumpTracer); 

        switch(n){
        case 0:           
          u = u;                    //ie. jumprate is uniform (as I state in my thesis)
          Info = 1.0;               //Maximum value for uniform random number
          break;
        case 1: 
          u = (-log(u)) / lambda;
          u = 1.0/u;
          Info = lambda;
          break;
        case 2:
          //powerlaw, p(y) = alfa/y_c * (y/y_c)^{-1-alfa} if y > y_c, 0 otherwise
          u = y_c * pow(u,-1/alpha);  
          u=1.0/u;                   //convert friction coefficient to jump rate
          Info = alpha;
          break;
        case 3:
          if (particle == 0) u = jumpTracer;   //set juprate for the first
          else u = jumpCrowders;             //set juprate for the rest
          Info = jumpTracer/jumpCrowders;
          break;
        }
        crowd.set_dist(n,Info);
      }
      else{
        if (n>3 || n <0 ) 
          cout <<"Invalid value/choise for prob.distribution ("<<n<<")"<<endl;
        else
          cout <<"Random number for jumprate must be 0 < r < 1"<<endl;
          
        abort();
      }

      //Manually set the jumprate of the tracer particle (first one)
      if (particle == 0 ) u = jumpTracer;
      
      //save jump rate u to right and left direction
      Direction dir(u,u);
      
      Jump jump;

      //same jumprate in all dimensions (x,y,z)
      jump.x = dir;
      jump.y = dir;
      jump.z = dir;

      hopRate.push_back(jump);
    }
  }

  //if using jump rates (i.e. not moveOld()-function)
  if (nymetod){
    crowd.setJumpRate(hopRate);
  }


  //double dt = crowd.get_dt();  //Only used this when I used my heatmap-function

  //if logaritmic spaceing between data points
  double delta_t = 1.0*MaxTime/antal_punkter;
  double delta_log_t = 1.0*log(MaxTime)/antal_punkter;    
  
  //TEST:
  //====================
  //When I use the inteeraction code, make sure I have Nakazato distribution
  //with identical jumprates!
  if ((InteractOn && DISTRIBUTION!=3) || 
      (JUMPRATE_CROWDER!=JUMPRATE_TRACER && DISTRIBUTION==3 && InteractOn)){
    cout<<endl<<"WARNING!! USING INTERACTION-CODE BUT NOT IDENTICAL JUMPRATES! \a"
        <<endl<<endl;
    //TODO! Whis isn't it working!
    // char answer;
    // cout <<"Abort simulation? y/n:"<<endl;
    // cin >> answer;
    // cin.ignore();
    // answer=getchar();
    // answer=getchar();
    // answer=getchar();
    // if (answer == 'y') 
    //   abort();
  }
  //====================

  //non-important nice-to-have simulation info.
  RemainingTime printToScreen(ensemble);  

  for(int E=0; E < ensemble; E++){

    if (!quiet)
      printToScreen.printProgress(E);
                  
    vector<int> x_disp, y_disp, z_disp;
    vector<double> r_disp;
    crowd.place();

    double t_output=0;
    int n=1;
    double totalTime;
    for (totalTime=0; totalTime <= MaxTime; totalTime = crowd.get_time()){

      if (totalTime >= t_output ){ //&& totalTime < t_output_next){
        crowd.set_dt();                            //Save the time
        x_disp.push_back(crowd.get_dX());          //get displacement
        y_disp.push_back(crowd.get_dY());
        z_disp.push_back(crowd.get_dZ());
        r_disp.push_back(crowd.get_dR());

        if(logscale){
          t_output = exp(delta_log_t *1.0*n);
        }
        else{  //t_output = t_output + delta_t;
          t_output = 1.0* n * delta_t;
        }
        
        n++;
      } //run this loop "antal_punkter"-times

      //(nymetod) ? crowd.move() : crowd.moveOld();
      crowd.move(); //faster I would guess... (moveOld() is obsolete)

      // //TEST: to get a snapshot close to the MaxTime
      // if(totalTime>=MaxTime-0.01){
      //   crowd.Snapshot();
      //   cout<<"totalTime:"<<totalTime<<endl;
      // }


      if(InteractOn){//INTERACTION
        //This is just to save the information of clustersize distribution
        //crowd.saveCluster(totalTime);
      }//BUG Augusti, gar otroligt trogt med denna!
    
    }
    //TODO: sista värdena skrivs aldrig ut, loopen slutar före utskriften...
    //...det är give or take, för vill att den börjar i nollan.
    

    //Store tracer position at each point for each ensemble,
    //inorder to get errorbars (and to bin and get a histogram of P(x,t))
    crowd.store(x_disp, y_disp, z_disp, r_disp); 

  }

  if(InteractOn){//INTERACTION
    //crowd.saveBinning(ensemble);        //print a histogram of clusters 
  }//BUG Augusti, gar otroligt trogt med denna!


  crowd.save(ensemble, NameOfFile);   // calculate std.dev. and save to file

  if (!nymetod) cout <<"NOTE: Using the old moveOld()-function"<<endl;

  return 0;
}



//////////////////////////////////
//         N Y  K O D           //
//////////////////////////////////


void Partiklar::setJumpRate(vector<Jump> jumpRate){
  if (jumpRate.size() == nParticles_){

    // remove jumprates from previous runs 
    // but now we are quenched so not needed (I think) TODO (is this right?)
    jumpRate_.clear();
    jumpRate_ = jumpRate;

    //Use these vectors to build the partial sum
    computePartialSum();         

  }
  else{
    //printError("Error, number of jumprates must equal number of particles");
    cout <<"Error, number of jumprates must equal number of particles"
         <<endl;
    abort();
  }
}


void Partiklar::computePartialSum(){
  long double tmp = 0;
  int i;

  // delete previous run
  partialSum_.clear();         

  // first element must be zero
  partialSum_.push_back(tmp);
                               
  //check the jump-rate-vectors  
  if ( jumpRate_.size()!=nParticles_ )
    printError("Error in number of jumprates, != N");
  
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
  for (int i = 1; i < partialSum_.size(); i++){
    if (partialSum_[i] < partialSum_[i-1] ){
      cout <<"Warning! The cumulative sum of jump rates is not "<<endl
           <<"a monotonically increasing sum!"<<endl;
      abort();
    }
  }
}
 

int Partiklar::convertMuToParticle(int mu, int& n){

  if (0 <= mu && mu <= partialSum_.size() -1){

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

double Partiklar::diffEffConstant(void){
  double D_eff;
  
  if(!partialSum_.empty()){ //Calculate If we use jumpratios (as opposed to moveOld())

    //NOTE: if the lattice constant: a!=1, it must be included
    //D_i=k_i * pow(a,2)
    double sum_temp=0;
  
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
    D_av_=0;
    for (int i=0; i < nParticles_; i++){
      D_av_ = D_av_ + jumpRate_[i].x.r;
      D_av_ = D_av_ + jumpRate_[i].x.l;
    }
    if(dim_ > 1){
      for (int i=0; i < nParticles_; i++){
        D_av_ = D_av_ + jumpRate_[i].y.r;
        D_av_ = D_av_ + jumpRate_[i].y.l;
      }
    }
    if (dim_>2){
      for (int i=0; i < nParticles_; i++){
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
    //svar: med denna D_eff får vi f(x)=2*dim_*D_eff*x, vilket är inkonsekvent med ovanstående 
    //Deff som redan är kompenserad för det... tror jag... TSTAR NU: N=1, uniform distrubution...
    //RESULTAT: nej, även med ovanst. D_eff utr. måste vi ha f(t)=2*dim_*D_eff*t
    //Detta är det rimmliga, eftersom vid N>1 vill vi ha den "riktiga" D_eff, dvs. som det är nu.
    //Note, assumed lattice constant a=1.
    D_av_  = 1.0 / (2*dim_);
    D_eff = 1.0 / (2*dim_);
    return D_eff;
  }
}


/*
  LEFT TO DO:
  ----------
  -remove various testing-loops
  -plot directly with gnuplot.
  -remove the set_interaction --> constructor! (doesn't work)
  -clean up main-fuction get-set-calls
  -rename all variables using camelCase
  -activate the m-flag!
  -remove the set_time, +3-4 get-functions...
  -exponential waiting time!
  -kanske ha en MSD-klass, double x y z r?
  -ta bort dx dy dz dr.
  -kompilera med -Wall -pedantic
*/

//There are various testing functions in here, but they are either silenced by use of if-statements in combination with a boolean "test" variable set to false, or just commented out. A lot of work went into testing the output, therefore some things can seem redundant, such as printing un-squared displacements (<dx>~0) etc. 


/*

  OUTLINE OF THE FUNCTION CALLS OF THIS PROGRAM:
  ---------------------------------------------

 With                      Without
 Interaction:              Interaction:
 -----------               ------------
 main()                     main()                  
  |                          |                     
 move()                     move()                 
  |                          |                     
 moveAndBoundaryCheck()     moveAndBoundaryCheck()  
  |                          |                     
 vacancyCheck()             vacancyCheck()          
  |                          |                     
 Interaction()               |                     
  |                          |                     
 CountNeighbors()            |                     
  |                          |                     
  |                          |                     
 [make the actual move]     [make the actual move]

*/
