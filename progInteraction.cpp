#include <iostream>        //för cout m.m
#include <fstream>         //läsa&skriva till filer
#include <cstdlib>         //för abort()-funktionen.
#include <ctime>           //för att använda systime som seed
#include <vector>          //för vektorklassen
#include <cmath>           //för pow-, sqrt-funktionerna
#include "nr3.h"           //Random number generator from nummerical recipes 3 ed.
#include "ran.h"           //obs, inkludera inte denna vid g++-kompileringen

using namespace std;

class Partiklar {
  vector<int> Xpos,Ypos,Zpos;                 //XYZ-position for each particle
  vector<double> t;
  vector<long double> x2_err,y2_err,z2_err,r2_err; //Sore the final error, in each point
  vector<double> x2_mu,y2_mu,z2_mu,r2_mu;     //medelvärdet av förf. i kvad: <X^2>, <Y^2>... 
  vector<double> x_mu, y_mu, z_mu, r_mu;      //medelvärde: <X>, <Y>...

  vector<double> jmpRateXright,jmpRateXleft,jmpRateYright
    ,jmpRateYleft,jmpRateZup,jmpRateZdown;
  vector <long double> partial_sum;           //inehåller summan av hopprates

  //This will store the original jumprate of each particle. It is only needed/used in
  //the interaction2()-function.
  vector<double>  JmpRateStorageXright, JmpRateStorageYright,
    JmpRateStorageZup, JmpRateStorageXleft, JmpRateStorageYleft, JmpRateStorageZdown;

  int RutorX,RutorY,RutorZ;                   //gittrets storlek
  int N, dim;                                 //antal partiklar och dimention (=1 || 2 || 3)
  float density;                              //particle density
  int X_0, Y_0, Z_0;                          //Initial positonion for the fist particle
  double time_sum;                            //Sum of time for each move
  int ConvertMuToParticle(int mu, int& n);
  void PartialSum(void);
  void MoveAndBoundaryCheck(int,int,bool);
  int VacancyCheck(int, int,int,int);         //No double Occupancy.
  int VacancyCheck2(int, int,int,int, bool);  //Same as above, but faster.

  void Interaction(int, int, int, int, bool);
  void Interaction2(int, int, int, int);

  void StdErr(int, int);                      //beräknar standard deviation

  vector< vector<int> > dx, dy, dz;           //alla pos. för vit part. behövs för sdrerr
  vector< vector<double> > dr;

  vector < vector<vector<int> > > vacancy;      //Needed in my more efficient "vacancyCheck2()"
  //  vector < vector<vector< <vector<int> > > > > vacancyStore; //Needed for "save_cluster()" 


  int seed1, seed2; 
  float gamma_white, gamma_black;             //jumprate for the Nakazato-particles
  int distribution;
  double Info;  
  int round(float x) {                        //Round off
    return int(x + 0.5);  }
  bool scaling;
  bool testOnOff;                             //för teastutskrift, för felsökning
  double D_av; //TEST/TMP...

  float Nakazato(bool);          //TODO: test//temp
  double DiffEffConstant(void);               //Calculate the effective dif. const. from the jumprates
  float Ergodicity(int,int,int, bool);        //Kontroll: ergodicitetens teoretiska värde
  void CountNeighbours(int, int, int, vector<int>&, bool );

  vector<float> Cluster_Distribution, Cluster_Size;  //For our save_cluster()-function
  void BuildCluster2(int, vector<int>&, double, bool); //needed for save_cluster()-func.
  void CheckVacancyMatrix(float);                      //Just a test-function


public:
  Partiklar(int, int, int, int, bool);
  void Place(void);
  void Move1(bool);
  void Move2(bool);
  void Lagra(vector<int>&, vector<int>&, vector<int>&, vector<double>&);
  void Save(int, char[]);

  void save_binning(int);                     //To get P(x,t)
  void save_cluster(double, bool);               //To get rho(m) (cluster distribution)
  void print_cluster(void);                   //Print  rho(m) (cluster distribution)

  void set_dt(void);                          //för initiering av tids-vektorn
  int get_dX(void);
  int get_dY(void);
  int get_dZ(void);
  void set_gamma(float, float);
  void set_scaling(bool);
  void set_dist(int, double);
  double get_dR(void);
  double get_time(void);
  int getDimension(void);

  void setJumpRate(vector<double>&, vector<double>&, vector<double>&,
                   vector<double>&, vector<double>&, vector<double>&);
  void setJumpRate(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
  void setJumpRate(vector<double>&, vector<double>&);
};

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



void Partiklar::save_cluster(double Time, bool Boundary){
  //This funciton calculates how many particles are in 
  //the cluster that the tagged particle is part of

  //ha denna funktionen sa att den kors varje tidssteg for sista ensemblen,
  //Den skall veta nar den ar pa sista tidssteget, tror jag...

  //Cluster_Distribution defined in the constructor.

  vector<int> Cluster;
  //TODO kolla Boundary! 0 eller 1?

  double Exponent = 0; //ie. rand > 0 allways include the neighbor.  

  // CheckVacancyMatrix(1.01); //OK hit

  if(Time>20){  //The system is not in equilibrium from start
    //CheckVacancyMatrix(10.1);
    CountNeighbours(Xpos[0], Ypos[0], Zpos[0], Cluster, Boundary);
    //CheckVacancyMatrix(10.2);
    BuildCluster2(0, Cluster, Exponent, Boundary); //Skriver denna till nagot?
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

//---------------------------------------------------------------
//IMPORTED FROM TESTInteractionCode() 12/7 -2010
//NOTE! BUT! "Exponent" changed to allow 100% probability of being included.
void Partiklar::BuildCluster2(int n, vector<int>& NearestNeighbours, 
                              double Exponent, bool boundary){
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

  bool testinglite;  //test-variable to find potential errors in the code
  bool CheckAgain;   //true if a new particle was added to the cluster

  do{
    CheckAgain=false;

    int limit = NearestNeighbours.size();
    for(int i=0; i< limit; i++){
      int m=NearestNeighbours[i];
      CountNeighbours(Xpos[m],Ypos[m],Zpos[m], NearestNeighbours, boundary);
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
} //BUG funkar inte om det är olika jumprates för olika ensambles

int Partiklar::getDimension(void){ //could be useful to have, but not neccesary...
  return dim;
}

int Partiklar::get_dX(void){
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

void Partiklar::set_gamma(float Gamma_White, float Gamma_Black){
 gamma_white = Gamma_White;
 gamma_black = Gamma_Black;
}

void Partiklar::set_dist(int choise, double Characteristic){
  distribution = choise;
  Info = Characteristic;
}
void Partiklar::set_scaling(bool invalue){
  scaling = invalue;
}


float Partiklar::Nakazato(bool rescale){
  //skall här försöka få in uttrycket från pappret...
  float alpha;

  if (dim==1) alpha = 1;
  if (dim==2) alpha = 0.363;
  if (dim==3) alpha = 0.209;

  //jumpratio for the different particles... 
  //OBS 1 svart, resten vita!
  double gamma = gamma_black / gamma_white;

  float c = density; 

  double nominator   = (gamma * (1-c)+1)*(1-alpha);
  double denominator = gamma*(1-c) + 1-alpha*(1 + gamma*(1-3*c));
  double f = nominator / denominator;
  float a = 1;            //lattice constant (might become a private variable soon)

  float D = gamma_black * (1-c) * f * pow(a,2);

  // Compensate for the rescaling of the axis
  if (rescale) D = D / DiffEffConstant(); 

  //due to laziness, we multiply the "true" difusion constant
  //with 2*dimenstion, and return that. Keep that in min when plotting.
  D = D * dim * 2;

  return D;
}


////////////////////////////
//       CONSTRUCTOR      //
////////////////////////////

Partiklar::Partiklar(int Xsquare,int Ysquare,int Zsquare, 
                     int ParticleNumber, bool testYesNo){

  if (Xsquare >= Ysquare && Ysquare  >= Zsquare && Zsquare >= 1){
    if (ParticleNumber <= Xsquare*Ysquare*Zsquare){
      seed1=1;
      seed2=2;                       //kanske ändra hur dessa får sina värden!

      RutorX=Xsquare;                //Antal rutor
      RutorY=Ysquare;
      RutorZ=Zsquare;
      N=ParticleNumber;

      X_0=round(1.0*RutorX/2);   //Mitten på gittret, där vi sätter "vit" partikel
      Y_0=round(1.0*RutorY/2);   
      Z_0=round(1.0*RutorZ/2);
      
      density=(float) N / (RutorX*RutorY*RutorZ); 

      //Hardcoded Switches:
      testOnOff=testYesNo;           //för testutskrift, för felsökning

      //Define the dimention of the lattice
      if (RutorY<=1 && RutorZ<=1) dim=1; 
      if (RutorY>1 && RutorZ<=1) dim=2;
      if (RutorY>1 && RutorZ>1) dim=3;
      

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


//kontrollerar vad ergodiciteten borde bli
float Partiklar::Ergodicity(int x_rutor, int y_rutor, int z_rutor, bool scaling){ 
  float summa=0; 
  int x,y,z;                            //koord som vi "går" med genom alla rutor
                                        //för att mäta avståndet till mitten med

  //I remapp the starting coordinates, incase I'd 
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


  //TEST: START - ANALYTIC EXPRESSION
  //this is just an alternate way of calculating <R^2>, 
  //when the system becomes erogodic.
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
  int RutorTotalt=RutorX*RutorY*RutorZ;  //Antal lediga rutor
  int AttPlacera=N;			 //antal ej utplacerade partiklar
					
  Xpos.assign(N,0);			 //X-koordinaten för alla N partiklar
  Ypos.assign(N,0);			
  Zpos.assign(N,0);			
					
  time_sum=0;				 //Time starts at zero each run
					
  Xpos[0]=X_0;				 //place black particle
  Ypos[0]=Y_0;				
  Zpos[0]=Z_0;				
  AttPlacera--;				 // (en mindre eftersom placerar ut "vit" först)  
  RutorTotalt--;			 // en ledig ruta mindre

  double R;
  static Ran randomNumber1(seed1);

  int i,j,k,n=1;                        //First white particle at n = 1 element    
  for (k=1; k<=RutorZ; k++){            //place particles
    for (j=1; j<=RutorY; j++){
      for (i=1; i<=RutorX; i++){
        if (i!=X_0 || j!=Y_0 || k!=Z_0){ //don't place white on black particle
          R=randomNumber1.doub();                    //Using NR to generate 0 < R < 1 ;
	  if (R<=AttPlacera*1.0/RutorTotalt){
	    Xpos[n]=i;
	    Ypos[n]=j;
	    Zpos[n]=k;
	    n++;
	    AttPlacera--;			    //En partikel mindre att placera
	  }  
	  //    cout << "AttPlacera: "<<AttPlacera<<endl;
	  RutorTotalt--;                            //En ruta mindre att använda
	}
      }
    }
  }
  //TEST:
  if (testOnOff)
    for (int q=0; q < N; q++)
      cout <<" PLACERADE \n partikel n="<<q<<" X="<<Xpos[q]<<" Y="<<Ypos[q]<<" Z="<<Zpos[q]<<endl;


  //IF We use VacancyCheck2()-function:
  //-----------------------------------
  //initiate the vacancy-vector with zeros. 
  //Must reset this matrix for each ensamble!
  vector< vector< vector<int> > >
    Vacancy(RutorX+1, vector< vector<int> >(RutorY+1, vector<int>(RutorZ+1,-1)));
  vacancy=Vacancy;
  //Mark occupied sites as 1:
  //note, no change will ever occur to the 0-elements.
  for (n=0; n<N; n++){
    vacancy[ Xpos[n] ][ Ypos[n] ][ Zpos[n] ]=n;
  }

}



void Partiklar::Move1(bool val){ //val=1 -> fix, val=0 -> periodiskgräns

  //  for (int i=0; i < N; i++){          //make N moves
  
    static Ran randomNumber2(seed2);
    int r,n; 
  
    do{
      n=randomNumber2.doub() * N; 
    }while(n==N);                     //Choose particle at random, 0<= n < N.
      

    do{ //funkar
      r=randomNumber2.doub() * (2 * dim); 
    }while( r== dim * 2 );            //gives 0<=R<2*dim
      

    MoveAndBoundaryCheck(n,r, val);
    //Do the actual move, acording to the boundary rules


    time_sum = time_sum + 1.0/N;    //Calculate the timey
 

    //TEST:
    if (testOnOff)
      for (int q=0; q < N; q++)
        cout <<" FLYTTADE \n partikel n="<<q<<" X="<<Xpos[q]<<" Y="<<Ypos[q]<<" Z="<<Zpos[q]<<endl;
  
    //  } //make N moves end.
}


void Partiklar::MoveAndBoundaryCheck(int n, int R, bool boundary){
  if (0 <= n && n < N ){
    if (testOnOff) cout <<"Moving particle "<< n <<endl;
    // n = particle to move, R = direction to move, such that: 
    //      R= 0 = X-right,      1 = X-left
    //         2 = Y-right,      3 = Y-left
    //         4 = Z-up,         5 = Z-down
    // bondary: boundary=0 = periodic, boundary=1 = fix boundary

    //TODO: skulle kunna lägga in kontroll för dim vs. R.

    int Xtemp=Xpos[n];       //Store in case we want to restore
    int Ytemp=Ypos[n];       //the previous configuration             
    int Ztemp=Zpos[n];

    if ( R >= 0 && R < 6 ){
      switch(R){                                 
      case 0: Xpos[n]=Xpos[n]+1;
        if (Xpos[n]==RutorX+1 && boundary==0) Xpos[n]=1;
        if (Xpos[n]==RutorX+1 && boundary==1) Xpos[n]=Xtemp;
        break;
      case 1: Xpos[n]=Xpos[n]-1;                            //move paricele n to the left
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

      //check that the new site is vacant
      //if not, move the particle back.

      //VacancyCheck(n,Xtemp,Ytemp,Ztemp);  //Old code
                                         
      VacancyCheck2(n,Xtemp,Ytemp,Ztemp,boundary);   //New improved, faster version.

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



int Partiklar::VacancyCheck(int n, int Xold, int Yold, int Zold){
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
        return returnvalue;
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



int Partiklar::VacancyCheck2(int n, int Xold, int Yold, int Zold, bool Boundary){
  //NOTE: This is a faster version than the old VacancyCheck(), since that one
  //loops through all particles for each "vacancy check"
  //However this uses a 3D matrix, and there will be some "double banking" 
  //in the sence that the postion is now stored in this matrix and in XYZ-vectors.

  if (n < N && 0 <= n){

    int returnvalue=0;                 //unused parameter to register collisions
    if (vacancy[Xpos[n]][Ypos[n]][Zpos[n]] != -1){  //if occupied
      // if (testOnOff) cout <<"Occupied ("<<Xpos[n]<<","<<Ypos[n]<<","<<Zpos[n]<<")"<<endl;
      Xpos[n]=Xold;                              //then move back          
      Ypos[n]=Yold;
      Zpos[n]=Zold;
      returnvalue=1;

      // if (testOnOff){
      //   for (int k=0; k<RutorZ+1; k++){
      //     for (int j=0; j<RutorY+1; j++){
      //       cout <<j<<": ";
      //       for (int i=0; i<RutorX+1; i++){
      //         cout <<vacancy[i][j][k]<<" ";
      //       }
      //       cout << endl;
      //     }
      //     cout <<"nytt Z:a"<<endl;
      //   }
      // }

      return returnvalue;
    }
    else{
      Interaction(n,Xold,Yold,Zold, Boundary);
      //if (testOnOff) cout <<"OK"<<endl;
      return returnvalue;
      //(the returnvalue could be used to have a do-move-while(VacancyCheck==1) but
      //if so, we need to move this function to be directly in the Move_()-func.)
    }
  }
 else{
   cout <<"Error, accessing invalid particle."<<endl;
   abort();
 }
}


void Partiklar::Interaction(int n, int Xold, int Yold, int Zold, bool boundary){
  vector<int> NearestNeighbours;
  //vector contains the particle labels of the nearest neighbors
  //to particle n (not sure that we actually need this at this point)

  //This is an updated version of the progInteraction.cpp I used in my thesis
  //It does the exact same thing, but now we can wrap around the boundary also...
  //The out-putfiles are identical between the two different version, when I use 
  //fixed boundary conditions, which was what I used in my thesis. 
  CountNeighbours(Xold, Yold, Zold, NearestNeighbours, boundary);

  int counter = NearestNeighbours.size();

  double r0;
  static Ran slump(56);
  do{
    r0=slump.doub();
  }while( r0==1 || r0==0 );

  float DeltaV = 3;  //Set interaction strength/energy!
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



void Partiklar::CountNeighbours(int X, int Y, int Z, vector<int>& Count, bool boundary){
  //Each element in the vector "Count" is a particle, and the 
  //value of the element is the particle label. 
  //If the particle is "alone" the vector Count will be empty.
  //boundary=0 periodic, =1 fix
  
  bool periodic;
  
  if (!boundary) periodic = true;
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


//Denna borde behöva ersätta vaccancyCheck tror jag....
void Partiklar::Interaction2(int n, int Xold, int Yold, int Zold){
  //This alternative Interaction-function geves rise to 6 new vectors 
  //to store the old jumprate of the particle in question.

  //vectors containing the labels for the particles at the old, and new site
  vector<int> clusterParticles, OldclusterParticles;
  for(int i=0; i < N; i++){
    //store the particles belonging to the cluster
    if (Xpos[i]==Xpos[n] && Ypos[i]==Ypos[n] && Zpos[i]==Zpos[n])
      clusterParticles.push_back(i); 
    if (Xpos[i]==Xold && Ypos[i]==Yold && Zpos[i]==Zold)
      OldclusterParticles.push_back(i);
  }

  //Code for moving FROM a cluster
  //----------------------------
  //Brake free:
  bool something=false;
  if(something){
    if(OldclusterParticles.size()==0)        //the site had single occupancy
      vacancy[Xold][Yold][Zold]=-1;          //set old position as vacant

    vacancy[Xpos[n]][Ypos[n]][Zpos[n]]=n;    //new position is now occupied (again?)
    //IDÉ: ha att vacancy sparar antalet partiklar på varje kluster istället?

    //update/restore the jumprate for n:
    jmpRateXright[n] = JmpRateStorageXright[n];
    jmpRateYright[n] = JmpRateStorageYright[n];
    jmpRateZup   [n] = JmpRateStorageZup   [n];
    jmpRateXleft [n] = JmpRateStorageXleft [n];
    jmpRateYleft [n] = JmpRateStorageYleft [n];
    jmpRateZdown [n] = JmpRateStorageZdown [n];
  

    //TODO: update the jumprate for the remaining cluster

  }//Move all
  else{
    //move all other particles in cluster to the new site
    for (int i=0; i < clusterParticles.size(); i++){
      Xpos[ OldclusterParticles[i] ] = Xpos[n];
      Ypos[ OldclusterParticles[i] ] = Ypos[n];
      Zpos[ OldclusterParticles[i] ] = Zpos[n];
      vacancy[Xold][Yold][Zold]=-1;         //old position is now vacant
      vacancy[Xpos[n]][Ypos[n]][Zpos[n]]=n; //new position is now occupied 
    }
    //BUG: och om det är fler partiklar på nya platsen?
    //TODO: Kan inte flytta klustret om det är ett annat kluster på nya siten


  }


  //Code for moving TO a cluster
  //----------------------------
  //REJECT:
  //denna kanske inte skall finnas:
  //sannolikhet att flytta tillbaka partikeln:
  bool someCondition=false;
  if (someCondition){
    Xpos[n]=Xold; 
    Ypos[n]=Yold;
    Zpos[n]=Zold;
  }
  //ACCEPT:
  //Set/Update the new jumprate of the cluster
  double SUMxr, SUMyr, SUMzu, SUMxl, SUMyl, SUMzd;
  SUMxr=SUMyr=SUMzu=SUMxl=SUMyl=SUMzd=0;
  for(int i=0; i<clusterParticles.size(); i++){
    SUMxr=SUMxr + 1.0/JmpRateStorageXright[ clusterParticles[i] ];
    SUMyr=SUMyr + 1.0/JmpRateStorageYright[ clusterParticles[i] ];
    SUMzu=SUMzu + 1.0/JmpRateStorageZup   [ clusterParticles[i] ];
    SUMxl=SUMxl + 1.0/JmpRateStorageXleft [ clusterParticles[i] ];
    SUMyl=SUMyl + 1.0/JmpRateStorageYleft [ clusterParticles[i] ];
    SUMzd=SUMzd + 1.0/JmpRateStorageZdown [ clusterParticles[i] ];
  }
  //don't forget to add the particle n. (or include it in the clusterPart() to begin with...)
  SUMxr=SUMxr + 1.0/JmpRateStorageXright[n];
  SUMyr=SUMyr + 1.0/JmpRateStorageYright[n];
  SUMzu=SUMzu + 1.0/JmpRateStorageZup   [n];
  SUMxl=SUMxl + 1.0/JmpRateStorageXleft [n];
  SUMyl=SUMyl + 1.0/JmpRateStorageYleft [n];
  SUMzd=SUMzd + 1.0/JmpRateStorageZdown [n];

  SUMxr = 1.0/SUMxr;
  SUMyr = 1.0/SUMyr;
  SUMzu = 1.0/SUMzu;
  SUMxl = 1.0/SUMxl;
  SUMyl = 1.0/SUMyl;
  SUMzd = 1.0/SUMzd;

  //NOTE: To increase speed implement a matrix that contains 
  //which particles belong to which cluster. Must contain the 
  //particle labels asociated with each cluster/site.

  for (int i=0; i < clusterParticles.size(); i++){
    //Give all cluster particles a new reduced common jumprate 
    jmpRateXright[ clusterParticles[i] ] =SUMxr;
    jmpRateYright[ clusterParticles[i] ] =SUMyr;
    jmpRateZup  [  clusterParticles[i] ] =SUMzu;
    jmpRateXleft[  clusterParticles[i] ] =SUMxl;
    jmpRateYleft[  clusterParticles[i] ] =SUMyl;
    jmpRateZdown[  clusterParticles[i] ] =SUMzd;
  }





}




void Partiklar::Move2(bool boundfix){
  if (partial_sum.empty()){
    cout <<"Partial-sum vector has not been initiated"<<endl;
    abort();
  }

  CheckVacancyMatrix(11); //BUG! Juli
  
  int seed4=8;
  static Ran randomNumb(seed4);
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
    //with zeroes. Every 10^8 run of Move2() will lead to a situation where we actually have
    //p_rand==p[mu_guess], and changing the intervall below in the if statements, from 
    // "..<= .. <.." to "..< .. <=.." will result in an infinite loop every 10^8 turns 
    //(if 32 bit)


    do{   //finds the mu
      mu_guess = (int) ( (p_rand-p_left)*( mu_right-mu_left )
                         /(p_right - p_left) + mu_left );

      if ((double) partial_sum[mu_guess] <= p_rand && p_rand < (double) partial_sum[mu_guess+1])
        LoopAgain =false;

      else{
        if (p_rand < (double) partial_sum[mu_guess]){
          mu_right = mu_guess;
          p_right = partial_sum[mu_guess];
          //if (p_rand == (double) partial_sum[mu_guess]){
          //NOTE! This if-statement is an addition to the original gillespie
          //algorithm that will get stuch in an inf. loop every 10^8 run of Move2()
          //without it! ...Unless we interchange the boundary to closed - open.
          //LoopAgain=false;
          //mu_guess = mu_guess -1;
          //}
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
  
    //  cout <<"r = " <<r<<"\t \t n = "<<n<<"\t \t mu = "<<mu<<endl;  //TEST
  
    MoveAndBoundaryCheck(n,r,boundfix);
    //moves the particle n in direction r, acording to boundary conditions. 
    //if the new site is occupied, move back to previous position

    // //TEST
    // double r3;
    // do{
    //   r3=randomNumb.doub();
    // }while( r3==1 || r3==0 );
    // double temp = log(1.0/r3);
    // //TEST-SLUT
    //    temp=1.0;

    time_sum = time_sum + 1.0/partial_sum.back();   //Calculate the time

    //} //move N times end
}


void Partiklar::Lagra(vector<int>& x_displace, vector<int>& y_displace,
                      vector<int>& z_displace, vector<double>& r_displace){
  dx.push_back(x_displace);   //this gives dx[ensembler][number of values] 
  dy.push_back(y_displace);
  dz.push_back(z_displace);
  dr.push_back(r_displace);
}


void Partiklar::Save(int Ensemble, char name[]){
  char* Primary=name;

  int AntalPunkter=t.size();
  ofstream primary,secondary;
  primary.open(Primary);
  secondary.open("txyz.dat");


  //Compute <R>:
  x_mu.assign(AntalPunkter,0);
  y_mu.assign(AntalPunkter,0);
  z_mu.assign(AntalPunkter,0);
  r_mu.assign(AntalPunkter,0);

  for (int j=0; j<AntalPunkter; j++){  
    for (int i=0; i< Ensemble ; i++){
      x_mu[j]=x_mu[j] + 1.0*dx[i][j]/Ensemble;
      y_mu[j]=y_mu[j] + 1.0*dy[i][j]/Ensemble;
      z_mu[j]=z_mu[j] + 1.0*dz[i][j]/Ensemble;
      //  r_mu[j]=r_mu[j] + 1.0*dr[i][j]/Ensemble;  //använder inte denna...
 } 
  }

  //Compute <R^2>:
  x2_mu.assign(AntalPunkter,0);
  y2_mu.assign(AntalPunkter,0);
  z2_mu.assign(AntalPunkter,0);
  r2_mu.assign(AntalPunkter,0);

  for (int j=0; j<AntalPunkter; j++){  
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
  StdErr(Ensemble, AntalPunkter);


  //To print which distribution we used
  string dist;
  if (distribution == 0)    dist="uniform";
  if (distribution == 1)    dist="exponential";
  if (distribution == 2)    dist="powerlaw";
  if (distribution == 3)    dist="nakazako";

  float d_eff= DiffEffConstant();
  float d_av = D_av;


  //TEST -felsök gränserna i exp-funktionen!
  cout <<"# Lower limit: "<<(RutorX * RutorX) / (2*N*N*d_eff) <<endl;
  cout <<"# Upper limit: "<<(RutorX * RutorX) / (2*d_eff) <<endl;

 
  //print info about simulation to head of file:
  primary <<"#E = "<<Ensemble<<"\t N = "<<N<<"\t X-Y-Z: "
          <<RutorX<<"x"<<RutorY<<"x"<<RutorZ<<endl;
  primary <<"#Density: "<<density<<"\t MSD when erogodic: "
          <<Ergodicity(RutorX,RutorY,RutorZ,scaling)
          <<"\t distr: "<<dist<<" ("<<Info<<")"<<endl;
  primary <<"#D_eff: "<<d_eff<<"\t 2*d*D_nakazato: "
          <<Nakazato(scaling)<<"\t D_av: "<<d_av<<endl; 

    //print info about simulation to screen:
  cout <<"#E = "<<Ensemble<<"\t N = "<<N<<"\t X-Y-Z: "
          <<RutorX<<"x"<<RutorY<<"x"<<RutorZ<<endl;
  cout <<"#Density: "<<density<<"\t MSD when erogodic: "
       <<Ergodicity(RutorX,RutorY,RutorZ,scaling)
       <<"\t distr: "<<dist<<" ("<<Info<<")"<<endl;
  cout <<"#D_eff: "<<d_eff<<"\t 2*d*D_nakazato: "
       <<Nakazato(scaling)<<"\t D_av: "<<d_av<<endl; 

  
  for (int i=0; i<t.size(); i++){

    //rescale the axis of the primary output (dimensionsless)
    if (scaling){
      t[i] = t[i] * pow(density,2) * DiffEffConstant();            
      r2_mu[i] = r2_mu[i] * pow(density,2);
      r2_err[i]= r2_err[i]* pow(density,2);
    }

    //Print it to file:
    primary  << t[i]<<"\t"<<r2_mu[i]<<"\t" <<r2_err[i]<<endl;  //skiter i "r"
    secondary<< t[i]<<"\t"<< x_mu[i] <<"\t"<< y_mu[i] <<"\t"<< z_mu[i] <<"\t"
	     << x2_mu[i] <<"\t"<< y2_mu[i] <<"\t"<< z2_mu[i] <<"\t"
	     << x2_err[i]<<"\t"<<y2_err[i]<<"\t"<<z2_err[i]<<endl;
  }
  primary.close();
  secondary.close();

}

void Partiklar::StdErr(int Ensemble, int AntalPunkter){
  bool test=0; //testvariabel

  if (test){ //TEST
    for (int j=0; j < AntalPunkter ;j++){
      for (int i=0; i < Ensemble ; i++)
	cout <<dx[i][j]<<" ";
      cout <<endl;
    }
  }

  x2_err.assign(AntalPunkter,0);
  y2_err.assign(AntalPunkter,0);
  z2_err.assign(AntalPunkter,0);
  r2_err.assign(AntalPunkter,0);


  for (int j=0; j < AntalPunkter  ;j++){
    for (int i=0; i<Ensemble ; i++){
      x2_err[j]=x2_err[j] + pow( pow(dx[i][j],2) - x2_mu[j] ,2);
      y2_err[j]=y2_err[j] + pow( pow(dy[i][j],2) - y2_mu[j] ,2);
      z2_err[j]=z2_err[j] + pow( pow(dz[i][j],2) - z2_mu[j] ,2);
      r2_err[j]=r2_err[j] + pow( pow(dr[i][j],2) - r2_mu[j] ,2);
    }
  }

  for (int j=0; j < AntalPunkter  ;j++){
     x2_err[j]=sqrt( x2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
     y2_err[j]=sqrt( y2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
     z2_err[j]=sqrt( z2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
     r2_err[j]=sqrt( r2_err[j] / (1.0*Ensemble*(Ensemble-1)) );
  }     //Då E=1 --> inf !

  
  //testar Gaussfördelningen av punkterna:
  //--------------------------------------------------------------------
  float SUM68=0;
  for (int j=0; j<AntalPunkter; j++){            //korrektion, tänkt fel!
    double top = x2_mu[j] +0.5*x2_err[j] *       sqrt(Ensemble-1)*2;
    double bot = x2_mu[j] -0.5*x2_err[j] *       sqrt(Ensemble-1)*2;
    for (int i=0; i<Ensemble; i++){
      if (bot <= pow(dx[i][j],2) &&  top >= pow(dx[i][j],2) ){
	if (test){	
	  cout <<"undre gräns "<< bot <<endl;
	  cout << "punkten är "<< pow(dx[i][j],2) <<endl;
	  cout << "övre gräns "<< top <<endl;
	}
	SUM68++;
	if (test) cout <<"SUM68 = "<<SUM68<<endl;
      }
    }
  }
  cout <<"68,2 % (?)="<<1.0*SUM68/(Ensemble*AntalPunkter)<<endl;

  //--------------------------------------------------------------------
}


//void Distribution(double&, int, float, float);

//////////////////////////////////////
//   ****************************   //
//   * M A I N  F U N K T I O N *   //
//   ****************************   //
//////////////////////////////////////

int main(int argc, char* argv[]){
  int N;                                //Antal paritklar
  int X,Y,Z;                            //gitterstorlek
  int ensemble,AntalPunkter;
  double MaxTime;

  //DEFAULT VALUES:  
  bool fixBoundaryOn=true;                 //boundary condition, false = periodic
  bool testJaNej = (argc==4) ? true : false;   //med 5 argument -> testutskrift
  bool logscale=false; 
  bool nymetod=true; 
  bool rescale = false;
  //--------------

  char* NameOfFiles;
  char DefaultNameOfFile[]="out.dat";
  NameOfFiles=DefaultNameOfFile;

  if (argc>=2){
    NameOfFiles=argv[1];  //manually set the output-filename
    cout<<"Output-file changed to: "<<NameOfFiles<<endl;
  }  

  if (argc>=3){           //if the second arg=1 -> logscale=on
    int SecondArg = atoi(argv[2]);
    if (SecondArg==1 || SecondArg==0){
      logscale = SecondArg;
    }
    else
      cout <<SecondArg<<
        " is wrong, must be: 0=logscale OFF, or 1=logscale ON"<<endl;
  }


  if (testJaNej) cout <<"Testutskrift på"<<endl;
  if (fixBoundaryOn) 
    cout <<"# Fixt gränsvilkor"<<endl;  
  else 
    cout <<"# Periodiskt gränsvilkor"<<endl;

  cout <<"# Gittrets storlek:"<<endl;
  cout <<"# X: ";                       
  cin >>X;
  cout<<"# Y: ";
  cin >>Y;
  cout<<"# Z: ";
  cin >>Z;

  cout <<"# antal partiklar: ";           // vit + svarta
  cin >> N;
 
  cout <<"# Antal rutor per partikel: "<<1.0*(X*Y*Z)/N<<endl; 
  cout <<"# Antal simuleringar: ";
  cin >>ensemble;
  cout <<"# Antal datapunkter: ";
  cin >> AntalPunkter;
  cout <<"# Max-tid: ";
  cin >>MaxTime;
 
  Partiklar crowd(X,Y,Z,N,testJaNej);            //X,Y,Z-squares, N-particles
  crowd.set_scaling(rescale);                    //if we rescale the axis

  int seed3=17; //17 //time(0);
  static Ran tempnumb(seed3);
  vector<double> jXr, jYr, jZu, jXl, jYl, jZd; 
  if (nymetod){
    for (int particle = 0; particle < N; particle ++){
      double u;               
      do {                     
        u=tempnumb.doub();   
      }while (u==0);

      //okila jumprates på svart/vit (nakazato):
      //----------------------

      int n = 3;  // <---< VAL AV METOD <----<
      // 0 = uniform
      // 1 = exponential
      // 2 = power law
      // 3 = black & white
      //Characteristic Trait of the choosen distribution, 
      //just used to print info to file/screen
      double Info=0; 
      float gamma_black;
      
      if (0 <= n && n<=3 && (u < 1 && u > 0) ){
        gamma_black = 0.25;    //Jumprate for the tracer particle!
        float lambda=1;        //used in option 1 Exp-dist.
        float y_c = 1;         //used in option 2 Power-law
        float alpha=0.5;       //used in option 2 Power-law
        float gamma_white = 0.25; //used in option 3 Nakazato
        crowd.set_gamma(gamma_white, gamma_black); 

        switch(n){
        case 0:            //Vanlig uniform,, Funkar 1D-crowding
          u = u;           // (men lite dålig vid låga t)
          Info = 1.0;      //Maximum value for uniform random number
          // cout <<"uniform"<<endl;
          break;
        case 1: 
          u = (-log(u)) / lambda;   //funkar inte alls...eller?
          Info = lambda;
          //cout <<"exp"<<endl;
          break;
        case 2:
          //powerlaw, p(y) = alfa/y_c * (y/y_c)^{-1-alfa} if y > y_c, 0 otherwise
          u = y_c * pow(u,-1/alpha);   //Funkar 1D-crowding skapligt (y_c=1, alfa=1)...
          Info = alpha;
          //cout <<"pow"<<endl;
          break;
        case 3:
          if (particle==0) u = gamma_black; //set juprate for the first
          else u = gamma_white;             //set juprate for the rest
          Info = gamma_black/gamma_white;
          //cout <<"naka"<<endl;
          break;
        }
        crowd.set_dist(n,Info);
      }
      else {
        if (n>3 || n <0 ) 
          cout <<"Invalid value/choise for prob.distribution ("<<n<<")"<<endl;
        else cout <<"Random number for jumprate must be 0 < r < 1"<<endl;
        abort();
      }
    
      //sätter jumprate för första till ;
      if (particle == 0 ) u = gamma_black;

      
      jXr.push_back(u);
      jXl.push_back(u);
      jYr.push_back(u);
      jYl.push_back(u);
      jZu.push_back(u);
      jZd.push_back(u); 

      //good to see if I get insane powerlaw-tails, (1*e+8)
      //cout <<"distribution: "<<u<<endl;
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


  double delta_t=  1.0*MaxTime/AntalPunkter;
  double delta_log_t= 1.0*log(MaxTime)/AntalPunkter;    
    
  for(int E=0; E < ensemble; E++){
    
    cout << "#---> REMAINING ENSEMBLES: "<<ensemble-E<<endl;

    vector<int> x_disp, y_disp, z_disp;
    vector<double> r_disp;
    crowd.Place();

    double t_output_next=MaxTime; //bara tillfällig variabel...
    double t_output=0;
    double punkt=0;
    int n=1;
    double totalTime;
    for (totalTime=0; totalTime <= MaxTime; totalTime = crowd.get_time()){
      //cout <<"totaltime: "<<totalTime<<endl;

      if (totalTime >= t_output ){ //&& totalTime < t_output_next){
        //cout <<"t_output: "<<t_output<<endl;
        crowd.set_dt();                            //Save the time
        x_disp.push_back(crowd.get_dX());          //get displacement
        y_disp.push_back(crowd.get_dY());
        z_disp.push_back(crowd.get_dZ());
        r_disp.push_back(crowd.get_dR());

        if(logscale){
          //TODO, ta bort första punkten 0 0 0 , funkar inte med LOG!
          // t_output = exp(punkt);
          // punkt = punkt + delta_log_t;
          t_output = exp(delta_log_t *1.0*n);
          //cout << "check "<<endl;
          //          t_output_next = exp(delta_log_t *(n+1));
        }
        else{  //t_output = t_output + delta_t;
          t_output = 1.0* n * delta_t;
          //  t_output_next =  1.0* (n+1) * delta_t;
        }
        
        n++;
      } //Denna loop körs AntalPunkter gånger...

      //(nymetod) ? crowd.Move2(fixBoundaryOn) : crowd.Move1(fixBoundaryOn);
      crowd.Move2(fixBoundaryOn); //faster I would guess... (Move1() is obsolete)

      //This is just to save the information of clustersize distribution
      crowd.save_cluster(totalTime, fixBoundaryOn);


    }
    //TODO: sista värdena skrivs aldrig ut, loopen slutar före utskriften...
    //...det är give or take, för vill att den börjar i nollan.
    crowd.Lagra(x_disp, y_disp, z_disp, r_disp); //Lagra hela ensemblens värden   0 
  }

  crowd.save_binning(ensemble);        //print a histogram
  crowd.Save(ensemble, NameOfFiles);   // calculate std.dev. and save to file

  //to print the inf. on cluster distribution:
  crowd.print_cluster(); 

  if (!nymetod) cout <<"NOTE: Using the old move1()-function"<<endl;
  
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

  //TEST-Börja-----skriv ut hela summan:

  // for (int j=0; j<partial_sum.size(); j++)
  //   cout <<partial_sum[j]<<"\t";
  // cout <<endl;

  //TEST-Slut------

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
}

double Partiklar::DiffEffConstant(void){
  double D_eff, sum_temp=0;
  
  if(!partial_sum.empty()){ //Calculate If we use jumpratios (as opposed to Move1())

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


    //TEST: (har den som en lite grisig variabel...)
    D_av=0;
    for (int i=0; i < jmpRateXright.size(); i++){
      D_av =D_av + jmpRateXright[i];
      D_av =D_av + jmpRateXleft[i];
    }
    //WORKS FOR 1D! NOT FOR MORE!
    //NEED TO PUT IT IN THE LOOPS ABOVE, IN THE CASE...
    D_av = D_av * 1.0/(N*2);
    //TEST-END

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

    D_eff = 1.0 / (2*dim);
    return D_eff;
  }
}


//BUG/TODO: Om man tilldelar nya jumprates varje ensemble, så får man olika värden på tiden för samma element i <X>...
//BUG Exp-fördelning, stämmer inte med analytiskt 1D.


/*
void printGnuplot(char name[], int ToGnuPlot, int ensemble, int distribution){
  bool logscale = true;
  ofstream plot;

  double density = Density();
  double D_eff=DiffEffConstant();
  double D_naka=Nakazato();
  double ergo = Ergodicity(RutorX,RutorY,RutorZ);

  plot.open(name);
  //Lite information:
  plot <<"# E = "<<ensemble<<"\t N = "<<N<<"\t X-Y-Z: "
       <<RutorX<<"x"<<RutorY<<"x"<<RutorZ<<"\t rho="<<density<<endl;
  //TODO: skriv ut vilken fördelning det är på jumprates...

  if (writeToFile){
    //set size 0.8, 0.8
    plot<<"set term post color linewidth 1.5"<<endl;
    plot<<"set out \""<<name<<".ps\"  "<<endl;
  }
  if(logscale) plot<<"set logscale"<<endl;

  plot <<"set ylabel \"<r^2>\""<<endl;
  plot <<"set xlabel \"t\""<<endl;
  plot <<"set legend bottom"<<endl;


  //TODO, manick som skriver vad det är för fördelning. 
  if (distribution==)

  if (dim==1) plot <<"f(x)=((1-"<<density<<")/"<<density<<")*sqrt(4*"<<D_eff<<"*x/pi)"<<endl;
  
  plot <<"e(x)="<<ergo<<endl;

  plot <<"n(x)= "<<D_naka<<"*x"<<endl;



  //om vi bill ha kvar plotten på skärmen, vill vi inte stänga filen!
  if (ToGnuPlot!=2) plot<<"quit"<<endl;
  plot.close();
  if (ToGnuPlot!=0) system("gnuplot plot.gp");

  //till att börja med: skriv alltd plot.gp-fil
  //toGnuPlot: 0= ej plotta till skärm eller fil
  //toGnuPlot: 1= plotta till ps-fil
  //toGnuPlot: 2= plotta till skärm

}

*/

