#include <sstream>         //behövs kanske för "ostringstream" i save()-func.
#include <string>          //behövs kanske inte?
#include <iostream>        //for cout, among others
#include <fstream>         //read write to files
#include <cstdlib>         //for abort()-function
#include <vector>          //for the vector class
#include <cmath>           //gives us sqrt, pow, fabs, 
#include <algorithm>       //needed for min_element in hist/bin-function

#include "nr3.h"           //Random number generator from numerical recipes 3 ed.
#include "ran.h"  

#include "auxilary.cpp"    //non-physics functions used in main() for parameter input/output
#include "classes.h"       //my own classes


// use ssh:
//      $ ssh -n 'cd /to/path; nice -19 ./program <input.dat'
// Or:  $ ssh -X nice -19 ./program output.dat < input.dat'

//HARDCODED VARIABLES:                                             SETS THEM IN:
#define DISTRIBUTION     0    //0=uniform, 1=exponential,         | (main)
//                              2=power-law, 3=nakazato           | (main)
#define JUMPRATE_TRACER  1.5  //                                  | (main)
#define JUMPRATE_CROWDER 0.5  //(if nakazato-distribution)        | (main)
#define FIXBOUNDARY      1    //1=fix wall, 0=periodic boundary   | (main ->class)
#define TEST_1           0    //boolean, prints details to screen | (constructor)
#define LOWMEM           1    //0=false / 1=true                  | ()

//gjort kommentarer till Main():s början.

//now test my InteractionRoutine, with LowMem. (and without!):
// at -i 2, c=0.7 k_t=0.25, LowMem=ON  Identical!
// at -i 3, c=0.4 k_t=0.25, LowMem=OFF Identical!
// at -i 2, c=0.5 k_t=0.25, LowMem=ON  Identical!

//now test LowMem ON (and this new merged code) for interactions off
// at k_t=0.5 c=0.5,         2Duniform,  Identical!
// at k_t=0.5 c=0.5 -i 0 (!) 2Duniform,  Identical! (got warning about not using Nakazato)

using namespace std;

class Partiklar {
  vector<int> Xpos,Ypos,Zpos;              //XYZ-position for each particle on the lattice

  vector<long double> x2_err,y2_err,z2_err,r2_err; //Store std error, in each point
  vector<double> x2_mu,y2_mu,z2_mu,r2_mu;  //average of displacement squared: <X^2>,... 
  vector<double> x_mu, y_mu, z_mu, r_mu;   //average of displacement: <X>, <Y>...
  vector<double> t;                        //store the time to be matched to the stored pos.

  vector<double> jmpRateXright,jmpRateXleft,jmpRateYright,
    jmpRateYleft,jmpRateZup,jmpRateZdown;  //jump-rate for each particle and direction
  vector <long double> partial_sum;        //cumulative sum of jump-rates

  int RutorX,RutorY,RutorZ;                //Lattice size
  int N;                                   //Total number of particles
  int dim;                                 //dimension of lattice (=1 || 2 || 3)
  bool BoundaryFix;                        //"true" if fix, "false" if periodic
  float density;                           //particle density
  int X_0, Y_0, Z_0;                       //Initial position for the tagged particle
  double time_sum;                         //Sum of time for each move

  vector< vector<int> > dx, dy, dz;        //pos. of tracer particle. needed for stdrerr
  vector< vector<double> > dr;             //same as line above

  int ConvertMuToParticle(int mu, int& n);
  void PartialSum(void);
  void MoveAndBoundaryCheck(int,int);
  int VacancyCheckOld(int, int,int,int);   //No double Occupancy.
  int VacancyCheck(int, int,int,int);      //Same as above, but improved. Use this one instead.
  vector< vector< vector<int> > > vacancy; //Needed in my more efficient "vacancyCheck2()"


  void StdErr(int);                        //calculate standard deviation
  void StdErrLowMem(int);                  //a more memmory conservative version of the above
  bool LowMem;                             //Switch between the two StdErr-functions above. 
  vector<double> dr4_err,dx4_err,dy4_err,dz4_err; //needed if we use StdErrLowMem()

  float jump_crowders, jump_tracer;        //jumprate for the (Nakazato) particles
  int distribution;
  double Info;                             //save characteristic trait of particle distribution
  bool scaling;
  bool testOnOff;                          //To print detailed information to screen
  double D_av; 

  bool InteractionOn;                      //use the interaction algorithm
  float InteractionStrength;               //Interaction strength to use, if "InteractionOn=true"

  float Nakazato(bool);         
  double DiffEffConstant(void);            //Calculate the effective dif. const. from the jump-rates
  float Ergodicity(int,int,int, bool);     //Calculate the theoretical MSD when ergodic

  int round(float x) {                     //Round off
    return int(x + 0.5);  }

  //INTERACTION CODE
  void Interaction(int, int, int, int);
  void CountNeighbours(int, int, int, vector<int>&);
  vector<float> Cluster_Distribution, Cluster_Size;  //For our save_cluster()-function
  void BuildCluster2(int, vector<int>&, double);     //needed for save_cluster()-func.
  void CheckVacancyMatrix(float);                    //Just a test-function


public:
  Partiklar(int, int, int, int, bool);
  void Place(void);                        //Place particles on the lattice
  void MoveOld(void);                      //old obsolete move-code (no jump-rates)
  void Move2(void);                        //pick a particle & direction based on jump-rate
  void Lagra(vector<int>&, vector<int>&, vector<int>&, vector<double>&);
  void Save(int, char[]);                  //Print simulation-run to file   
  void set_dt(void);                       //needed to initiate time-vector
  double get_dR(void);
  double get_time(void);
  int get_dX(void);
  int get_dY(void);
  int get_dZ(void);
  int getDimension(void);
  void set_jumpNaka(float, float);
  void set_scaling(bool);
  void set_dist(int, double);
  void set_interaction(float);             //determine if we use the interaction-algorithm
  void set_LowMem(bool, int);              //Set the switch true/false for use of low-mem run

  //Set jump-rate of each particle and direction, depending on if dim=1,2,3. 
  void setJumpRate(vector<double>&, vector<double>&, vector<double>&,
                   vector<double>&, vector<double>&, vector<double>&);
  void setJumpRate(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
  void setJumpRate(vector<double>&, vector<double>&);

  //Secondary not used functions
  void Snapshot(void);                     //Not used. (more comments in function body)
  //  double get_dt(void);   //Size of time-step, need it for P(x,t), don't need this...


  //INTERACTION CODE:
  void save_binning(int);                     //To get P(x,t)
  void save_cluster(double);               //To get rho(m) (cluster distribution)
  void print_cluster(void);                   //Print  rho(m) (cluster distribution)


};


void Partiklar::save_cluster(double Time){
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

  // CheckVacancyMatrix(1.01); //OK hit

  if(Time>20){  //The system is not in equilibrium from start
    //CheckVacancyMatrix(10.1);
    CountNeighbours(Xpos[0], Ypos[0], Zpos[0], Cluster);
    //CheckVacancyMatrix(10.2);
    BuildCluster2(0, Cluster, Exponent); //Skriver denna till nagot?
    //CheckVacancyMatrix(10.3); //OK hit
    int M = Cluster.size();
    //First element is # of 1-clusters (ie. single particles), etc.
    Cluster_Distribution[M-1] = Cluster_Distribution[M-1] + 1.0/N;
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

void Partiklar::print_cluster(void){
  //This function prints the information gathered by 
  //the save_cluster()-function, and prints it to file

  //  CheckVacancyMatrix(9);

  if (Cluster_Size.size() != Cluster_Distribution.size()){
    cout<<"gaaa "<<endl;
    abort();
  }

  ofstream clust;
  char NameOfFiles[]="cluster.dat";
  clust.open(NameOfFiles);
  for(int j=0; j < Cluster_Distribution.size(); j++){
    clust << Cluster_Size[j] <<"\t"<< Cluster_Distribution[j]<<endl;
  }
  clust.close();
  cout<<"Pinted "<<NameOfFiles<<endl;
}


void Partiklar::save_binning(int Ensemble){
  //This function is for getting P(t,x), ie. the probability
  //distribution of finding a particle at x, at time t.
  if(!LowMem){
    int No_bins = 300;  //number of bins
    vector<vector<float> > PrintToFile; 
  
    //-------------------
    // TODO: cutoff vid +49? skall vara -50...+50!
    // TODO: System not in equilibrium at t=0, sample anyway?

    double Min = (double) -RutorX/2;  //max displacement from starting position (at least for dim=1) 
    double Max = (double)  RutorX/2;
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
    for(int k = 0; k < t.size() ; k++ ){ // OBS "k" is time! (element pos is time)
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
    for(int i=0; i <t.size(); i++){
      for(int j=0; j < No_bins; j++){
        hist << Xaxis[j] + bin_step/2 <<"\t"<< t[i] <<"\t"<< PrintToFile[i][j]<<endl;
      }
      hist <<endl; //separate data block with new (empty) line
    }
    hist.close();
    cout<<"Pinted output to: "<<NameOfFiles<<endl;
  
    //if (scaling){
    // t[i] = t[i] * pow(density,2) * DiffEffConstant();
    //}
  }
  else
    cout<<"LowMem=ON, therefore no information stored!"
        <<endl<<"can't print cluster-information"<<endl;
}



//IMPORTED FROM TESTInteractionCode() 12/7 -2010
//NOTE! BUT! "Exponent" changed to allow 100% probability of being included.
void Partiklar::BuildCluster2(int n, vector<int>& NearestNeighbours, 
                              double Exponent){
  vector<int> unique;
  vector<int> NewParticleAdded;
  static Ran RandomNUM(42);
  double random_number;

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
      int m=NearestNeighbours[i];
      CountNeighbours(Xpos[m],Ypos[m],Zpos[m], NearestNeighbours);
      //NearestNeighbours-vector will store more particles each 
      //iteration, but only loop over the first initial set (ie. i < limit)
    }

    for(int j=0; j < NearestNeighbours.size(); j++){
      bool SingleValued=true;
 
      //Check if it's stored in unique[] allready:
      for(int k=0; k < unique.size(); k++){
        if (NearestNeighbours[j] == unique[k])
          SingleValued=false; 
      }
      if(SingleValued){
        //Add the new particle, with Boltzmans consent...
        // random_number = RandomNUM.doub();
        //if (random_number >=  exp(-1.0*Exponent)){
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
    for(int j=0; j < NearestNeighbours.size(); j++){
      int m= NearestNeighbours[j];
      cout <<"m = "<<m<<" ("<<Xpos[m]<<","<<Ypos[m]<<","<<Zpos[m]<<")"<<endl;    
    }
  }

  //TEST
  for (int i=0; i < N; i++){
    if (Xpos[i]==0 || Ypos[i]==0 || Zpos[i]==0){
      cout <<"Utanför brädet!"<<endl;
      abort();
    }
  }

  //TEST---------(FÖR ATT KOLLA DOUBLE COUNTING)
  bool IsUniqueTrue=true;
  for(int i=0; i<unique.size();i++){
    for(int j=0; j<unique.size();j++){
      if (unique[i]==unique[j] && i != j)
        IsUniqueTrue=false;
    }
  }
  if(!IsUniqueTrue){
    cout <<"Unique is not unique! Error *235897*"<<endl;
    abort();
  }//-------------
}



void Partiklar::CheckVacancyMatrix(float ErrorCode){
  //This is just a small function to check that the vacancy matrix is correct
  //ie. vacant sites are marked as "-1" and occupied are marked with the 
  //correct particle label.
  //ErrorCode = code to output to easier find the place in the code

  bool error = false;
  int label;
  int count_vacant_sites = 0;

  for (int i=1; i <= RutorX; i++){
    //cout<<"i: "<<i<<endl;
    for (int j=1; j <= RutorY; j++){
      //cout<<"j: "<<j<<endl;
      for (int k=1; k <= RutorZ; k++){ 
        //cout<<"k: "<<k<<endl;
        label = vacancy[i][j][k];
        //cout <<"label: "<< label <<endl;
        if(label != -1){
          if(i != Xpos[label] || j != Ypos[label] || k != Zpos[label]){
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

  if (count_vacant_sites != RutorX*RutorY*RutorZ-N){
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
//   return partial_sum.back();
// }

void Partiklar::Snapshot(void){
  //NOTE! This is only for 2D, with isotropic jumprate,
  //and the same system sidelength!
  //Its purpose is to print a snapshot of the system, 
  //to see the spatial distribution of different particles
  //as a heat map. Only used this function briefly. 

  vector< vector<double> > matris;

  vector<double> temp;
  temp.assign(RutorY,0);
  matris.assign(RutorX,temp);

  for(int i=0; i<Xpos.size(); i++){
    matris[ Xpos[i]-1 ][ Ypos[i]-1 ] = jmpRateXright[i];
  }

  ofstream heat;
  char NameOfFiles[]="heatmap.dat";
  heat.open(NameOfFiles);
  
  bool points=true; //true: plotar endast partiklarna, false: plottar aven tomma rutor (=0)
  if(points){
    //This does not give 0 where there are no particles:
    for (int i=0; i<Xpos.size(); i++)
      heat  << Xpos[i]<<" \t"<< Ypos[i] <<"\t"<< jmpRateXright[i]<< endl;
  }
  else{
    //This gives zeroes for all empty sites:
    for(int i=0; i<RutorX; i++){
      for(int j=0; j<RutorY; j++){
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
    t.push_back(time_sum);
  }
  //this will store the remaining values of the first enasmble
  else{
    if ( t.back() < time_sum ){
      t.push_back(time_sum);
      // cout <<"tidsfunktionen: "<<t.back()<< " time_sum: "<<time_sum<<endl;
    }
  }
} //NOTE funkar inte om det är olika jumprates för olika ensambles

int Partiklar::getDimension(void){ //could be useful to have, but not neccesary...
  return dim;
}

int Partiklar::get_dX(void){
  //current position minus starting position.
  return Xpos[0]-X_0;
}

int Partiklar::get_dY(void){
  return Ypos[0]-Y_0;
}

int Partiklar::get_dZ(void){
  return Zpos[0]-Z_0;
}

double Partiklar::get_dR(void){
  double r =sqrt( pow(Xpos[0]-X_0,2) + pow(Ypos[0]-Y_0,2) + pow(Zpos[0]-Z_0,2) ); 
  return r;
}

double Partiklar::get_time(void){
  return time_sum;
}

void Partiklar::set_jumpNaka(float Jump_Crowders, float Jump_Tracer){
  jump_crowders = Jump_Crowders;
  jump_tracer = Jump_Tracer;
}

void Partiklar::set_dist(int choise, double Characteristic){
  distribution = choise;
  Info = Characteristic;
}
void Partiklar::set_scaling(bool invalue){
  scaling = invalue;
}
void Partiklar::set_interaction(float InputInteract){
  InteractionStrength=InputInteract; 
  InteractionOn=true; 
}
void Partiklar::set_LowMem(bool LowMemOn, int NumberOfPoints){
  //Set the switch true/false for use of low-mem run
  LowMem=LowMemOn;

  //initiate the vectors for the summation in Lagra()
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



float Partiklar::Nakazato(bool rescale){
  float alpha;

  if (dim==1) alpha = 1;
  if (dim==2) alpha = 0.363;
  if (dim==3) alpha = 0.209;

  //jumpratio for the different particles... 
  //OBS 1 svart, resten vita!
  double gamma = jump_tracer / jump_crowders;

  float c = density; 

  double nominator   = (gamma * (1-c)+1)*(1-alpha);
  double denominator = gamma*(1-c) + 1-alpha*(1 + gamma*(1-3*c));
  double f = nominator / denominator;
  float a = 1;            //lattice constant (might become a private variable soon)

  float D = jump_tracer * (1-c) * f * pow(a,2);

  // Compensate for the rescaling of the axis
  if (rescale) D = D / DiffEffConstant(); 

  //due to laziness, we multiply the "true" diffusion constant here
  //with 2*dimension, and return that. Keep that in mind when plotting.
  D = D * dim * 2;

  return D;
}


////////////////////////////
//       CONSTRUCTOR      //
////////////////////////////

Partiklar::Partiklar(int Xsquare,int Ysquare,int Zsquare, 
                     int ParticleNumber, bool boundary){

  if (Xsquare >= Ysquare && Ysquare  >= Zsquare && Zsquare >= 1){
    if (ParticleNumber <= Xsquare*Ysquare*Zsquare){

      RutorX=Xsquare;            //number of lattice sites
      RutorY=Ysquare;
      RutorZ=Zsquare;

      N=ParticleNumber;

      BoundaryFix=boundary;

      X_0=round(1.0*RutorX/2); //center of lattice, where we place our tagged particle
      Y_0=round(1.0*RutorY/2); 
      Z_0=round(1.0*RutorZ/2); 
      
      density=(float) N / (RutorX*RutorY*RutorZ); 

      //Hardcoded Switches:
      testOnOff=TEST_1;           //boolean, for print details to screen. 

      //Define the dimention of the lattice
      if (RutorY<=1 && RutorZ<=1) dim=1; 
      if (RutorY>1 && RutorZ<=1) dim=2;
      if (RutorY>1 && RutorZ>1) dim=3;

      //By default: don't use the interaction algorithm. 
      InteractionStrength=0;
      InteractionOn=false;   

      //By default: don't use the memmory conserving StdErrLowMem()
      //use the set_LowMem()-function in the main()-func to change this
      LowMem=false;      

      //--------------------------
      //To get the Cluster Distribution (see save_cluster()):
      // Cluster_Distribution = number of clusters with "m" particles
      // Cluster_Size = number of particles in cluster (just our x-axis)
      int Maximus_Clustirus =N; 
      //Maximal cluster size to include in distribution
      Cluster_Distribution.assign(Maximus_Clustirus,0); //assign zeros...
      for (int i=1; i <= Maximus_Clustirus ; i++){
        Cluster_Size.push_back(i);
      }//--------------------------
 

     
    }
    else {
      cout <<"Number of particles must be <= number of squares"<<endl;
      abort();
    }
  }
  else{
    cout<<"X,Y,Z must be 1 or bigger, and X >= Y >= Z"<<endl; 
    abort();
  }
}


float Partiklar::Ergodicity(int x_rutor, int y_rutor, int z_rutor, bool scaling){ 
  //check what the MSD should be when the system reaches equilibrium.
  //(only correct if N=1). This function has two independent parts. One where we just
  //"step through" the lattice and add the displacement from the center, and another
  //were we actually calculate it analytically. 
  float summa=0; 
  int x,y,z;                            //koord som vi "går" med genom alla rutor
                                        //för att mäta avståndet till mitten med

  //STEP THROUGH THE LATTICE
  //I re-map the starting coordinates, in case I'd 
  //like to move this function out from this class...
  int mittX=X_0;
  int mittY=Y_0;
  int mittZ=Z_0;

  for (int i=0; i < x_rutor; i++){
    for (int j=0; j < y_rutor; j++){
      for (int k=0; k < z_rutor; k++){
	x=i+1;                          //First square have coord.(1,1,1)
	y=j+1;
	z=k+1;
	summa = summa + pow(mittX-x,2)+pow(mittY-y,2)+pow(mittZ-z,2);
      }
    }
  }
  summa = (float) summa / (x_rutor*y_rutor*z_rutor);


  //ANALYTIC EXPRESSION
  //this is just an alternate way of calculating <R^2>, 
  //when the system becomes ergodic.
  bool testAnalytic = false;
  if (testAnalytic){
    int a,b,c;
    a = mittX;  
    b = mittY;	
    c = mittZ;	
    float xContribution = ((x_rutor+1)*(2*x_rutor+1)*1.0/6 -a*(x_rutor+1) + a*a); 
    float yContribution = ((y_rutor+1)*(2*y_rutor+1)*1.0/6 -b*(y_rutor+1) + b*b); 
    float zContribution = ((z_rutor+1)*(2*z_rutor+1)*1.0/6 -c*(z_rutor+1) + c*c); 
    float SumAnalytic = (xContribution + yContribution + zContribution);
    cout <<"Analytic result for <R^2> (ergodicity): "<< SumAnalytic << endl;
  }

  if (scaling) summa = summa * pow(density,2);

  return summa;
}


void Partiklar::Place(void){
  //This function will be called once each ensemble to
  //place the N particles (randomly), and reset the time
  
  int RutorTotalt=RutorX*RutorY*RutorZ;  //number of vacant sites
  int AttPlacera=N;			 //number of not-yet-placed particles
					
  Xpos.assign(N,0);			 //coordinates for our N particles
  Ypos.assign(N,0);			
  Zpos.assign(N,0);			
					
  time_sum=0;				 //Time starts at zero each run
					
  Xpos[0]=X_0;				 //place tracer particle at center
  Ypos[0]=Y_0;				
  Zpos[0]=Z_0;				
  AttPlacera--;                          //one less to place
  RutorTotalt--;			 
  
  double R;
  static Ran randomNumber1(1);          //just any seed will do.

  int i,j,k,n=1;                        //First crowder at n = 1 element    
  for (k=1; k<=RutorZ; k++){            //place particles
    for (j=1; j<=RutorY; j++){
      for (i=1; i<=RutorX; i++){
        if (i!=X_0 || j!=Y_0 || k!=Z_0){ //don't place crowders on the tagged particle
          R=randomNumber1.doub();                    //Using NR to generate 0 < R < 1 ;
	  if (R<=AttPlacera*1.0/RutorTotalt){
	    Xpos[n]=i;
	    Ypos[n]=j;
	    Zpos[n]=k;
	    n++;
	    AttPlacera--;			    //One particle less to place
	  }  
	  //    cout << "AttPlacera: "<<AttPlacera<<endl;
	  RutorTotalt--;                            //One site less available
	}
      }
    }
  }
  //TEST:
  if (testOnOff)
    for (int q=0; q < N; q++)
      cout <<" PLACED \n particle n="<<q<<" X="<<Xpos[q]<<" Y="<<Ypos[q]<<" Z="<<Zpos[q]<<endl;


  //IF we use the new VacancyCheck()-function:
  //-----------------------------------
  //initiate the vacancy-vector with "-1" (meaning vacant) 
  //Must reset this matrix for each ensemble.
  vector< vector< vector<int> > >
    Vacancy(RutorX+1, vector< vector<int> >(RutorY+1, vector<int>(RutorZ+1,-1)));
  vacancy=Vacancy;
  //Mark occupied sites with particle label for that site:
  //note, no change will ever occur to the 0-elements of the matrix, 
  //since out lattice starts at 1,...,RutorX, but matrix is 0,...,RutorX
  for (n=0; n<N; n++){
    vacancy[ Xpos[n] ][ Ypos[n] ][ Zpos[n] ]=n;
  }
}



void Partiklar::MoveOld(){ //BoundaryFix=1 -> fix, BoundaryFix=0 -> periodic
  //This code is obsolete, as it just works for placing identical
  //particles, but I keep it for sentimental reasons.
  //Better safe than sorry.
  
  static Ran randomNumber2(2);        //any seed will do. 
  int r,n; 
  
  do{
    n=randomNumber2.doub() * N; 
  }while(n==N);                     //Choose particle at random, 0<= n < N.
      

  do{
    r=randomNumber2.doub() * (2 * dim); 
  }while( r== dim * 2 );            //gives 0<=R<2*dim
      

  //Do the actual move, acording to the boundary rules
  MoveAndBoundaryCheck(n,r);

  time_sum = time_sum + 1.0/N;    //Calculate the time

  //TEST:
  if (testOnOff)
    for (int q=0; q < N; q++)
      cout <<" MOVED \n particle n="
           <<q<<" X="<<Xpos[q]<<" Y="<<Ypos[q]<<" Z="<<Zpos[q]<<endl;
}


void Partiklar::MoveAndBoundaryCheck(int n, int R){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant. 
  bool boundary=BoundaryFix;
  
  if (0 <= n && n < N ){
    if (testOnOff) cout <<"Moving particle "<< n <<endl;
    // n = particle to move, R = direction to move, such that: 
    //      R= 0 = X-right,      1 = X-left
    //         2 = Y-right,      3 = Y-left
    //         4 = Z-up,         5 = Z-down
    // bondary: boundary=0 = periodic, boundary=1 = fix boundary

    int Xtemp=Xpos[n];       //Store in case we want to restore
    int Ytemp=Ypos[n];       //the previous configuration             
    int Ztemp=Zpos[n];

    if ( R >= 0 && R < 6 ){
      switch(R){                                 
      case 0: Xpos[n]=Xpos[n]+1;
        if (Xpos[n]==RutorX+1 && boundary==0) Xpos[n]=1;
        if (Xpos[n]==RutorX+1 && boundary==1) Xpos[n]=Xtemp;
        break;
      case 1: Xpos[n]=Xpos[n]-1;                            //move particle n to the left
        if (Xpos[n]==0 && boundary==0) Xpos[n]=RutorX;      //periodic boundary
        if (Xpos[n]==0 && boundary==1) Xpos[n]=Xtemp;       //fix boundary, let particle stay unmoved
        break;                                        //(first coordinate is 1, 0=outside) 
      case 2: Ypos[n]=Ypos[n]+1;
        if (Ypos[n]==RutorY+1 && boundary==0) Ypos[n]=1;
        if (Ypos[n]==RutorY+1 && boundary==1) Ypos[n]=Ytemp;
        break;
      case 3: Ypos[n]=Ypos[n]-1;
        if (Ypos[n]==0 && boundary==0) Ypos[n]=RutorY;
        if (Ypos[n]==0 && boundary==1) Ypos[n]=Ytemp;
        break;
      case 4: Zpos[n]=Zpos[n]+1;
        if (Zpos[n]==RutorZ+1 && boundary==0) Zpos[n]=1;
        if (Zpos[n]==RutorZ+1 && boundary==1) Zpos[n]=Ztemp;
        break;
      case 5: Zpos[n]=Zpos[n]-1;
        if (Zpos[n]==0 && boundary==0) Zpos[n]=RutorZ;
        if (Zpos[n]==0 && boundary==1) Zpos[n]=Ztemp;
        break;
      }

      //check that the new site is vacant if not,
      //move the particle back. (done by the VacancyCheck...)
      
      //VacancyCheckOld(n,Xtemp,Ytemp,Ztemp);  //Old code (slow)                                        
      VacancyCheck(n,Xtemp,Ytemp,Ztemp);   //New improved, faster version.

    }
    else{
      cout <<"There is a maximum of "<<2*dim<<" positive directions to move in"<<endl;
      abort();
    }
  }
  else{
    cout <<"Invalid particle label, n="<<n<<", must fulfill 0 <= n < N"<<endl;
    abort();
  }
}



int Partiklar::VacancyCheckOld(int n, int Xold, int Yold, int Zold){
  //Old and obsolete version. It doesn't need a 3D matrix as the
  //new version does, VacancyCheck(...) but this one is slow.
  //Keep it here to check that the output is exactly identical
  //(compare output-files with the "diff" command)
  
  if (n < N && 0 <= n){

    int returnvalue=0;    //unused parameter to register collisions 
    for (int i=0; i<N; i++){
      if (Xpos[i]==Xpos[n] && Ypos[i]==Ypos[n] && Zpos[i]==Zpos[n] && i!=n){ 
        if (testOnOff) cout <<"Occupied ("<<Xpos[n]<<","<<Ypos[n]<<","<<Zpos[n]<<")"<<endl;
        //if two particles on same square: move back to previous position
        Xpos[n]=Xold;                                        
        Ypos[n]=Yold;
        Zpos[n]=Zold;
        returnvalue=1;
        return returnvalue;     //End function here
      }
      else
        if (testOnOff) cout <<"OK"<<endl;
    }
    if (testOnOff) cout <<"Vacant, particle "<<n<<" X = "<<Xpos[n]<<", Y = "<<Ypos[n]<<", Z = " <<Zpos[n]<<endl;
    return returnvalue;
    //(the returnvalue could be used to have a do-move-while(VacancyCheck==1) but
    //if so, we need to move this function to be directly in the Move_()-func.)
  }
  else{
    cout <<"Error, accessing invalid particle "<<endl;
    abort();
  }
}



int Partiklar::VacancyCheck(int n, int Xold, int Yold, int Zold){
  //NOTE: This is a faster version than VacancyCheckOld(), since that one
  //loops through all particles for each "vacancy check"
  //However this uses a 3D matrix, and there will be some "double banking" 
  //in the sence that the postion is now stored in this matrix and in XYZ-vectors.

  if (n < N && 0 <= n){

    int returnvalue=0;                 //unused parameter to register collisions
    if (vacancy[Xpos[n]][Ypos[n]][Zpos[n]] != -1){  //if occupied
      if (testOnOff) cout <<"Occupied ("<<Xpos[n]<<","<<Ypos[n]<<","<<Zpos[n]<<")"<<endl;
      Xpos[n]=Xold;                              //then move back          
      Ypos[n]=Yold;
      Zpos[n]=Zold;
      returnvalue=1;

      if (testOnOff){
        for (int k=0; k<RutorZ+1; k++){
          for (int j=0; j<RutorY+1; j++){
            cout <<j<<": ";
            for (int i=0; i<RutorX+1; i++){
              cout <<vacancy[i][j][k]<<" ";
            }
            cout << endl;
          }
          cout <<"looping for new Z"<<endl;
        }
      }
    }
    else{

      //if no interaction att all:
      if(!InteractionOn){                            //if vacant, update vacancy-matrix 
        vacancy[Xold][Yold][Zold]=-1;             //old position is now vacant
        vacancy[Xpos[n]][Ypos[n]][Zpos[n]]=n;     //new position is now occupied
        if (testOnOff) cout <<"OK, vacant"<<endl;
      }
      //if we use the interaction algorithm:
      else
        Interaction(n,Xold,Yold,Zold);
    }
    //if (testOnOff) cout <<"Vacant, particle "<<n<<" X = "<<Xpos[n]<<", Y = "<<Ypos[n]<<", Z = " <<Zpos[n]<<endl;
    return returnvalue;
    //(the return value could be used to have a do-move-while(VacancyCheck==1) but
    //if so, we need to move this function to be directly in the Move_()-func.)
  }
  else{
    cout <<"Error, accessing invalid particle "<<endl;
    abort();
  }
}



void Partiklar::Move2(){
  //pick a particle & direction based on jump-rate, and move.
  if (partial_sum.empty()){
    cout <<"Partial-sum vector has not been initiated"<<endl;
    abort();
  }
  
  static Ran randomNumb(8);
  double r2;

  //  for (int i=0; i < N; i++){          //make N moves

  int mu_guess;                       // must be integer. (index of vector)
  int mu_left = 0;
  int mu_right = dim*2*N;
  double p_left = (double) partial_sum.front();  //value in the first element ( =0 )
  double p_right = (double) partial_sum.back();  //value in the last element ( =dim*2*N)
 
  do{
    r2 = randomNumb.doub();
  }while( r2==1 || r2==0 );

  double p_rand = (double) r2*partial_sum.back();
  bool LoopAgain=true;

  //NOTE: when comparing a long double against a double, C++ will fill in the blanks
  //with zeroes. Every 10^15 run of Move2() will lead to a situation where we actually
  //have p_rand==p[mu_guess], and changing the intervall below in the if statements,
  //from "..<= .. <.." to "..< .. <=.." will result in an infinite loop every 10^15
  //turns since "double" stores 52 bit, 2^52 = 10^15. This bug has now been fixed!

  //finds the mu
  do{ 
    mu_guess = (int) ( (p_rand-p_left)*( mu_right-mu_left )
                       /(p_right - p_left) + mu_left );

    if ((double) partial_sum[mu_guess] <= p_rand && p_rand < (double) partial_sum[mu_guess+1])
      LoopAgain =false;

    else{
      if (p_rand < (double) partial_sum[mu_guess]){
        mu_right = mu_guess;
        p_right = partial_sum[mu_guess];
      }
      else{//if (r2*partial_sum.back() > partial_sum[mu_guess+1]){
        mu_left=mu_guess+1;
        p_left=partial_sum[mu_guess+1];
      }
    }
  }while(LoopAgain);

  int mu = mu_guess;
  int n,r;        

  
  //transforms the index mu to which particle to move (index n, by reference)
  //and in which direction,by returning an integer 0<= r <= 5 (in 3D)
  r=ConvertMuToParticle(mu,n);
   
  MoveAndBoundaryCheck(n,r);
  //moves the particle n in direction r, according to boundary conditions. 
  //if the new site is occupied, move back to previous position

  time_sum = time_sum + 1.0/partial_sum.back();   //Calculate the time

  //} //move N times end
}

void Partiklar::Interaction(int n, int Xold, int Yold, int Zold){
  vector<int> NearestNeighbours;
  //vector contains the particle labels of the nearest neighbors
  //to particle n (not sure that we actually need this at this point)

  //This is an updated version of the progInteraction.cpp I used in my thesis
  //It does the exact same thing, but now we can wrap around the boundary also...
  //The out-putfiles are identical between the two different version, when I use 
  //fixed boundary conditions, which was what I used in my thesis. 
  CountNeighbours(Xold, Yold, Zold, NearestNeighbours);

  int counter = NearestNeighbours.size();

  double r0;
  static Ran slump(56);
  do{
    r0=slump.doub();
  }while( r0==1 || r0==0 );

  float DeltaV = InteractionStrength;          //specified in constructor
  if (r0 < exp(-1.0*counter*DeltaV) ){
    //accept move
    vacancy[Xold][Yold][Zold]=-1;              //old position is now vacant
    vacancy[Xpos[n]][Ypos[n]][Zpos[n]]=n;      //new position is now occupied    
  }
  else{
    //reject move, move back.
    //vacancy-matrix unchanged
    Xpos[n]=Xold; 
    Ypos[n]=Yold;
    Zpos[n]=Zold;
  }  
}



void Partiklar::CountNeighbours(int X, int Y, int Z, vector<int>& Count){
  //Each element in the vector "Count" is a particle, and the 
  //value of the element is the particle label. 
  //If the particle is "alone" the vector Count will be empty.
  //BoundaryFix=0 periodic, =1 fix
  
  //TODO invert the boolean condition, not that pretty perhaps:
  bool periodic;
  
  if (!BoundaryFix) periodic = true;
  else  periodic = false;

  if (X+1 <= RutorX ){ //avoid falling off the edge of the world
    if (vacancy[X+1][Y][Z] != -1){
      Count.push_back( vacancy[X+1][Y][Z] );
    }
  }
  else{
    if (vacancy[1][Y][Z] != -1 && periodic){
      Count.push_back( vacancy[1][Y][Z] );
    }
  }

  if (X-1>0){
    //although this zero-element exists in the vacancy matrix (i.e. no seg. fault)
    //the first coordinate on the lattice is (1,1,1). (0,0,0) is undefined
    if (vacancy[X-1][Y][Z] != -1){
      Count.push_back( vacancy[X-1][Y][Z] );
    }
  }  
  else{
    if (vacancy[RutorX][Y][Z] != -1 && periodic){
      Count.push_back( vacancy[RutorX][Y][Z] );
    }
  }



  if (Y+1 <= RutorY){
    if (vacancy[X][Y+1][Z] != -1){
      Count.push_back( vacancy[X][Y+1][Z] );
    }
  }
  else{
    if (vacancy[X][1][Z] != -1 && periodic){
      Count.push_back( vacancy[X][1][Z] );
    }
  }

  if (Y-1>0){
    if (vacancy[X][Y-1][Z] != -1){
      Count.push_back( vacancy[X][Y-1][Z] );
    }
  }
  else{
    if (vacancy[X][RutorY][Z] != -1 && periodic){
      Count.push_back( vacancy[X][RutorY][Z] );
    }
  }


  if (Z+1 <= RutorZ){
    if (vacancy[X][Y][Z+1] != -1){
      Count.push_back( vacancy[X][Y][Z+1] );
    }
  }
  else{
    if (vacancy[X][Y][1] != -1 && periodic){
      Count.push_back( vacancy[X][Y][1] );
    }
  }

  if (Z-1>0){
    if (vacancy[X][Y][Z-1] != -1){
      Count.push_back( vacancy[X][Y][Z-1] );
    }
  }
  else{
    if (vacancy[X][Y][RutorZ] != -1 && periodic){
      Count.push_back( vacancy[X][Y][RutorZ] );
    }
  } 
}



void Partiklar::Lagra(vector<int>& x_displace, vector<int>& y_displace,
                      vector<int>& z_displace, vector<double>& r_displace){
  //This is needed to calculate error bars (StdErr())
  
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
  //For LowMem=false, we do this in the Save & StdErr()-functions
  //since we use push_back() instead.
  if(LowMem){
    bool error=false;
    if (dr4_err.size()==0) error=true;
    if (dx4_err.size()==0) error=true;
    if (dy4_err.size()==0) error=true;
    if (dz4_err.size()==0) error=true;
    if (dr4_err.size()==0) error=true;

    if (r2_mu.size()==0) error=true;
    if (x2_mu.size()==0) error=true;
    if (y2_mu.size()==0) error=true;
    if (z2_mu.size()==0) error=true;
   
    if(error){
      cout<<"Incorrect assignment of vectors for StdErr-LowMem"
          <<endl<<"Must be assigned by set_LowMem()-function!"<<endl;
      abort();
    }
  }

  if(!LowMem){ //if don't use the LowMem-algorithm:
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


void Partiklar::Save(int Ensemble, char name[]){
  //Now process the simulation data, and print it to
  //file "name".
  
  char* Primary=name;

  int numberOfValues=t.size();
  ofstream primary,secondary;
  primary.open(Primary);

  //hardcoded switch to print x,y,z, x^2,y^2,z^2
  //specifically. Not needed now, since system is isotropic,
  bool xyz=false;
  if(xyz) secondary.open("txyz.dat");

  if(!LowMem){
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
    StdErr(Ensemble);
  }
  else{
    //use the memmory conserving version instead. This one does not
    //use the dx[][]-matrix, and in that way we save lots of mem.
    StdErrLowMem(Ensemble);
  }



  //To print which distribution we used to head of output-file
  string dist[]={"uniform","exponential","powerlaw","nakazato"};
  string onOff[]={"Off","On"};
  string bound[]={"periodic","fix"};
  float d_eff= DiffEffConstant();  //calculate effective diffusion constant
  float d_av = D_av;

  //TEST
  //cout <<"# Lower limit: "<<(RutorX * RutorX) / (2*N*N*d_eff) <<endl;
  //cout <<"# Upper limit: "<<(RutorX * RutorX) / (2*d_eff) <<endl;

  //print three lines of info about simulation to head of file:
  ostringstream print;
  print <<"#E = "<<Ensemble<<"\t N = "<<N<<"\t X-Y-Z: "
        <<RutorX<<"x"<<RutorY<<"x"<<RutorZ
        <<"\t 2*d*D_naka: "<<Nakazato(scaling)<<endl;
  print <<"#Density: "<<density<<"\t MSD_equil: "
        <<Ergodicity(RutorX,RutorY,RutorZ,scaling)
        <<"\t distr: "<<dist[distribution]<<" ("<<Info
        <<")\t bound: "<<bound[BoundaryFix]<<endl;
  print <<"#D_eff: "<<d_eff<<"\t k_tagg: "
        <<jmpRateXright[0]<<"\t D_av: "<<d_av
        <<"\t Interaction: "<<onOff[InteractionOn]<<" "
        <<InteractionStrength<<endl; 
  
  primary <<print.str(); //convert to string and send
  cout <<print.str();


  //print data to file:
  for (int i=0; i<t.size(); i++){

    //rescale the axis of the primary output (dimensionless) the
    //"scaling"-boolean is rather stupid as I've now learned to
    //manipulate the data in Gnuplot instead. Don't think too much 
    //about it, since this has the value "false", mostly. 
    if (scaling){
      t[i] = t[i] * pow(density,2) * DiffEffConstant();            
      r2_mu[i] = r2_mu[i] * pow(density,2);
      r2_err[i]= r2_err[i]* pow(density,2);
    }

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


void Partiklar::StdErr(int Ensemble){
  //Calculate the standard deviation, to get error-bars.
  int numberOfValues=t.size();

  if(LowMem){//Should not be using this function if LowMem=true!
    cout <<"LowMem switch ON, but using wrong StdErr-function!"
         <<endl;
    abort();
  }

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
  for (int j=0; j<numberOfValues; j++){            //forgot about his: 
    double top = x2_mu[j] +0.5*x2_err[j] *       sqrt(Ensemble-1)*2;
    double bot = x2_mu[j] -0.5*x2_err[j] *       sqrt(Ensemble-1)*2;
    for (int i=0; i<Ensemble; i++){
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


void Partiklar::StdErrLowMem(int Ensemble){
  int numberOfValues=t.size();

  if(!LowMem){
    cout<<"LowMem=OFF, but using StdErrLowMem()-func!"<<endl;
    abort();
  }

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

  //Following the notation in the comments in Partiklar::Lagra():
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
  bool fixBoundaryOn=FIXBOUNDARY;     //boundary condition, false = periodic
  bool logscale=false;
  bool nymetod=true;                  //use Move2() instead of Move()
  bool rescale = false;               //do this in gnuplot/matlab instead
  bool UseLowMem = LOWMEM;            //TODO! Patch this to the -m flag!
  bool quiet=false;                   //print remaining ensembles during simulation
  //--------------

  char NameOfFile[99]="out.dat";

  int InteractStr=0;
  bool InteractOn=false;

  //change default values depending on arg. flags given
  argumentFlags(argc, argv, logscale, InteractOn,
                quiet, InteractStr, NameOfFile);

  //TEST print result to screen:
  //-----------------------
  string OnOff[]={"OFF","ON"};
  int tempIndex=0;
  if (InteractOn) tempIndex=1;
  cout<<"###log="<<logscale<<endl
      <<"###interaction strength "<<OnOff[tempIndex]<<" ="<<InteractStr<<endl
      <<"###filename="<<(string) NameOfFile<<endl;
  //-----------------------


  if (fixBoundaryOn) 
    cout <<"# Fix boundary"<<endl;  
  else 
    cout <<"# Periodic boundary"<<endl;


  AskUserForInputParameters(X,Y,Z,N,ensemble,antal_punkter,MaxTime);


  Partiklar crowd(X,Y,Z,N,fixBoundaryOn);
  
  
  if(InteractOn)                           //use Interaction algorithm
    crowd.set_interaction(InteractStr);    //with "InteractStr" strength.

  crowd.set_scaling(rescale);              //if we rescale the axis

  //kanske ha den som invariabel i Save() istallet?
  //Set the switch true/false for use of low-mem run
  crowd.set_LowMem(UseLowMem, antal_punkter);

  static Ran tempnumb(17);                 //just anny seed will do
  vector<double> jXr, jYr, jZu, jXl, jYl, jZd; 
  if (nymetod){
    for (int particle = 0; particle < N; particle ++){
      double u;
      do { 
        u=tempnumb.doub();   
      }while (u==0);


      int n = DISTRIBUTION;  //Choose distribution
      // 0 = uniform
      // 1 = exponential
      // 2 = power law
      // 3 = nakazato

      //Characteristic trait of the chosen distribution, 
      //just used to print info to file/screen
      double Info=0;

      float jump_tracer;
      
      //NOTE: The jumprate is in EACH direction! Meaning, jumprate =1 is actually
      //a jumprate =4 in 2D. (I think)
      if (0 <= n && n<=3 && (u < 1 && u > 0) ){
        jump_tracer = JUMPRATE_TRACER; //Jumprate for the tracer paricle!
        float lambda=1.0;              //used in option 1 Exp-dist.
        float y_c = 1.0;               //used in option 2 Power-law
        float alpha=0.5;               //used in option 2 Power-law
        float jump_crowders=JUMPRATE_CROWDER; //used in option 3 Nakazato
        crowd.set_jumpNaka(jump_crowders, jump_tracer); 

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
          if (particle==0) u = jump_tracer;   //set juprate for the first
          else u = jump_crowders;             //set juprate for the rest
          Info = jump_tracer/jump_crowders;
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
      if (particle == 0 ) u = jump_tracer;
      
      jXr.push_back(u);
      jXl.push_back(u);
      jYr.push_back(u);
      jYl.push_back(u);
      jZu.push_back(u);
      jZd.push_back(u); 

    }
  }

  
  //jumprate -START-------------------
  if (nymetod){
    //Sätter jumprate enl. tilldelat ovan:
    int dimension=crowd.getDimension();
    if (dimension==1)  crowd.setJumpRate(jXr, jXl);                  
    if (dimension==2)  crowd.setJumpRate(jXr, jXl, jYr, jYl);        
    if (dimension==3)  crowd.setJumpRate(jXr, jXl, jYr, jYl, jZu, jZd); 
  }
  //jumprate -END-------------------


  //double dt = crowd.get_dt();  //Only used this when I used my heatmap-function

  //if logaritmic spaceing between data points
  double delta_t=  1.0*MaxTime/antal_punkter;
  double delta_log_t= 1.0*log(MaxTime)/antal_punkter;    
  
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
    crowd.Place();

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

      //(nymetod) ? crowd.Move2() : crowd.MoveOld();
      crowd.Move2(); //faster I would guess... (MoveOld() is obsolete)

      // //TEST: to get a snapshot close to the MaxTime
      // if(totalTime>=MaxTime-0.01){
      //   crowd.Snapshot();
      //   cout<<"totalTime:"<<totalTime<<endl;
      // }


      if(InteractOn){//INTERACTION
        //This is just to save the information of clustersize distribution
        //crowd.save_cluster(totalTime);
      }//BUG Augusti, gar otroligt trogt med denna!
    
    }
    //TODO: sista värdena skrivs aldrig ut, loopen slutar före utskriften...
    //...det är give or take, för vill att den börjar i nollan.
    

    //Store tracer position at each point for each ensemble,
    //inorder to get errorbars (and to bin and get a histogram of P(x,t))
    crowd.Lagra(x_disp, y_disp, z_disp, r_disp); 

  }

  if(InteractOn){//INTERACTION
    //crowd.save_binning(ensemble);        //print a histogram of clusters 
  }//BUG Augusti, gar otroligt trogt med denna!


  crowd.Save(ensemble, NameOfFile);   // calculate std.dev. and save to file

  if (!nymetod) cout <<"NOTE: Using the old MoveOld()-function"<<endl;

  return 0;
}



//////////////////////////////////
//         N Y  K O D           //
//////////////////////////////////

//---------------------------------------------------------------------------
// Överlagrad set-jumprate-funktion
//---------------------------------------------------------------------------
void Partiklar::setJumpRate(vector<double>& setJmpXright, vector<double>& setJmpXleft,
                            vector<double>& setJmpYright, vector<double>& setJmpYleft, 
                            vector<double>& setJmpZup,    vector<double>& setJmpZdown){
  if (dim!=3){
    cout <<"Error: The number of jump rate-directions and dimension of lattice does not correspond"<<endl;
    abort();
  }

  if (setJmpXright.size()==N && setJmpXleft.size()==N &&
      setJmpYright.size()==N && setJmpYleft.size()==N &&
      setJmpZup.size()==N    && setJmpZdown.size()==N){

    jmpRateXright.clear();  // remove jumprates from previous runs 
    jmpRateXleft.clear();
    jmpRateYright.clear();
    jmpRateYleft.clear();
    jmpRateZup.clear();
    jmpRateZdown.clear();
    
    jmpRateXright = setJmpXright;
    jmpRateXleft  = setJmpXleft;
    jmpRateYright = setJmpYright;
    jmpRateYleft  = setJmpYleft;
    jmpRateZup    = setJmpZup;
    jmpRateZdown  = setJmpZdown;

    PartialSum();         //Use these vectors to build the partial sum

  }
  else{
    cout <<"Error, number of jumprates must equal number of particles"<<endl;
    abort();
  }
}

//in 2D:
void Partiklar::setJumpRate(vector<double>& setJmpXright, vector<double>& setJmpXleft,
                            vector<double>& setJmpYright,vector<double>& setJmpYleft){
  if (dim!=2){
    cout <<"Error: The number of jump rate-directions and dimension of lattice does not correspond"<<endl;
    abort();
  }

  if (setJmpXright.size()==N && setJmpXleft.size()==N &&
      setJmpYright.size()==N && setJmpYleft.size()==N){

    jmpRateXright.clear();  // remove jumprates from previous runs 
    jmpRateXleft.clear();
    jmpRateYright.clear();
    jmpRateYleft.clear();
    
    jmpRateXright = setJmpXright;
    jmpRateXleft  = setJmpXleft;
    jmpRateYright = setJmpYright;
    jmpRateYleft  = setJmpYleft;

    PartialSum();         //Use these vectors to build the partial sum

    //TEST
    //for (int i=0; i<N; i++)  cout <<"n="<<i<<"   "<<jmpRateXright[i]<<" "<<jmpRateYleft[i]<<endl;
  

  }
  else{
    cout <<"Error, number of jumprates must equal number of particles"<<endl;
    abort();
  }
}

//in 1D
void Partiklar::setJumpRate(vector<double>& setJmpXright, vector<double>& setJmpXleft){
  if (dim!=1){
    cout <<"Error: The number of jump rate-directions and dimension of lattice does not correspond"<<endl;
    abort();
  } 
  
  if (setJmpXright.size()==N && setJmpXleft.size()==N){
    
    jmpRateXright.clear();  // remove jumprates from previous runs 
    jmpRateXleft.clear();
    
    jmpRateXright = setJmpXright;
    jmpRateXleft  = setJmpXleft;
    
    PartialSum();         //Use these vectors to build the partial sum
  }
  else{
    cout <<"Error, number of jumprates must equal number of particles"<<endl;
    abort();
  }
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

void Partiklar::PartialSum(){
  long double tmp=0;
  int i;

  partial_sum.clear();          // delete previous run

  partial_sum.push_back(tmp);   // index 0, 
                                // first element must be zero

  //kontrollerar jump-rate-vektorerna  
  if ( jmpRateXright.size()!=N || jmpRateXleft.size()!=N ){
    cout <<"Error in X, number of jumprates != N"<<endl;
    abort();
  }
  if ( (jmpRateYright.size()!=N || jmpRateYleft.size()!=N) && dim>1){
    cout <<"Error in Y, number of jumprates != N"<<endl;
    abort(); 
  }
  if ( (jmpRateZup.size()!=N || jmpRateZdown.size()!=N) && dim==3){
    cout <<"Error in Z, number of jumprates != N"<<endl;
    abort();
  }

  //bygger partial_sum-vektorn
  if(dim>=1){
    for(i=0; i < N ; i++){
      tmp=tmp+jmpRateXright[i];
      partial_sum.push_back(tmp);  //index 1 ... N
    }
    for (i=0; i < N ; i++){ 
      tmp=tmp+jmpRateXleft[i];
      partial_sum.push_back(tmp);  //index N+1 ... 2N
    }
  }
  if (dim >=2 ){
    for (i=0; i < N; i++ ){
      tmp=tmp+jmpRateYright[i];
      partial_sum.push_back(tmp);   //index 2N+1 ... 3N
    }
    for (i=0; i < N; i++ ){
      tmp=tmp+jmpRateYleft[i];
      partial_sum.push_back(tmp);   //index 3N+1 ... 4N
    }
  }
  if(dim==3){
    for (i=0; i < N; i++ ){
      tmp=tmp+jmpRateZup[i];
      partial_sum.push_back(tmp);   //index 4N+1 ... 5N
    }
    for (i=0; i < N; i++ ){
      tmp=tmp+jmpRateZdown[i];
      partial_sum.push_back(tmp);   //index 5N+1 ... 6N
    }
  }

  //TEST-Print the sum:
  // for (int j=0; j<partial_sum.size(); j++)
  //   cout <<partial_sum[j]<<"\t";
  // cout <<endl;

  //TEST: It must increase monotonically:
  for (int i = 1; i < partial_sum.size(); i++){
    if (partial_sum[i] < partial_sum[i-1] ){
      cout <<"Warning! The cumulative sum of jump rates is not "<<endl
           <<"a monotonically increasing sum!"<<endl;
      abort();
    }
  }
}
 

int Partiklar::ConvertMuToParticle(int mu, int& n){

  if (0 <= mu && mu <= partial_sum.size() -1){

    if ( 0 <= mu && mu <=2*N-1 ){
      if ( mu<=N-1 ){               
        n = mu;
        return 0;
      }
      if ( N <= mu){               
        n = mu - N;
        return 1;
      }
    }

    if ( 2*N <= mu && mu <= 4*N-1){
      if ( mu <= 3*N-1){
        n = mu -2*N;
        return 2;
      }
      if ( 3*N <= mu ){
        n = mu -3*N;
        return 3;
      }
    }

    if (4*N <= mu && mu <= 6*N-1){
      if ( mu <= 5*N-1){
        n = mu -4*N;
        return 4;
      }
      if ( 5*N <= mu ){
        n = mu -5*N;
        return 5;
      }
    }
  }
  else{
    cout << "mu =  "<< mu <<endl;
    cout << "Partial_sum.size is: "<<partial_sum.size()<<endl;
    cout << "Error, in mu or partial_sum"<<endl; 
    abort();
  }
  
  return -1;
}

double Partiklar::DiffEffConstant(void){
  double D_eff;
  
  if(!partial_sum.empty()){ //Calculate If we use jumpratios (as opposed to MoveOld())

    //NOTE: if the lattice constant: a!=1, it must be included
    //D_i=k_i * pow(a,2)
    double sum_temp=0;
  
    if (dim >= 1){
      for (int i=0; i < jmpRateXright.size(); i++){
        sum_temp =sum_temp + 1.0/jmpRateXright[i];
        sum_temp =sum_temp + 1.0/jmpRateXleft[i];
      }
    }

    if (dim>=2){
      for (int i=0; i < jmpRateYright.size(); i++){
        sum_temp =sum_temp + 1.0/jmpRateYright[i];
        sum_temp =sum_temp + 1.0/jmpRateYleft[i];
      }
    }

    if (dim==3){
      for (int i=0; i < jmpRateZup.size(); i++){
        sum_temp =sum_temp + 1.0/jmpRateZup[i];
        sum_temp =sum_temp + 1.0/jmpRateZdown[i];
      }
    }

    sum_temp= (double) sum_temp / (2*N*dim); 
    D_eff = 1.0/sum_temp;

    //Calculate Average as well:
    D_av=0;
    for (int i=0; i < jmpRateXright.size(); i++){
      D_av =D_av + jmpRateXright[i];
      D_av =D_av + jmpRateXleft[i];
    }
    if(dim > 1){
      for (int i=0; i < jmpRateYright.size(); i++){
        D_av =D_av + jmpRateYright[i];
        D_av =D_av + jmpRateYleft[i];
      }
    }
    if (dim>2){
      for (int i=0; i < jmpRateZup.size(); i++){
        D_av =D_av + jmpRateZup[i];
        D_av =D_av + jmpRateZdown[i];
      }
    }
    D_av = D_av * 1.0/(N*2*dim);


    return D_eff;

  }
  else{
    //When we don't use jumpratios we get the diffusion constant to:
    //TODO: stämmer detta? den skall ju bli <r^2>=2*D*dim*t ! -->
    //svar: med denna D_eff får vi f(x)=2*dim*D_eff*x, vilket är inkonsekvent med ovanstående 
    //Deff som redan är kompenserad för det... tror jag... TSTAR NU: N=1, uniform distrubution...
    //RESULTAT: nej, även med ovanst. D_eff utr. måste vi ha f(t)=2*dim*D_eff*t
    //Detta är det rimmliga, eftersom vid N>1 vill vi ha den "riktiga" D_eff, dvs. som det är nu.
    //Note, assumed lattice constant a=1.
    D_av  = 1.0 / (2*dim);
    D_eff = 1.0 / (2*dim);
    return D_eff;
  }
}


/*

  LEFT TO DO:
  ----------
  -remove the scaling boolean. No need for it any more.
  -remove various testing-loops
  -plot directly with gnuplot.
  -remove the set_interaction --> constructor! (doesn't work)
  -implement a particle class, with members x,y,z particle.x(n) (I think)
  -clean up main-fuction get-set-calls
  -rename all variables using camelCase
  -activate the m-flag!
  -remove the set_time, +3-4 get-functions...
*/

//There are various testing functions in here, but they are either silenced by use of if-statements in combination with a boolean "test" variable set to false, or just commented out. A lot of work wen into testing the output, therefore some things can seem redundant, such as printing un-squared displacements (<dx>~0) etc. 


/*

  OUTLINE OF THE FUNCTION CALLS OF THIS PROGRAM:
  ---------------------------------------------

  With                       Without
  Interaction:               Interaction:
  -----------                ------------
  main()                     main()                  
  |                          |                     
  Move2()                    Move2()                 
  |                          |                     
  MoveAndBoundaryCheck()     MoveAndBoundaryCheck()  
  |                          |                     
  VacancyCheck()             VacancyCheck()          
  |                          |                     
  Interaction()                |                     
  |                          |                     
  CountNeighbors()             |                     
  |                          |                     
  |                          |                     
  [make the actual move]     [make the actual move]

*/
