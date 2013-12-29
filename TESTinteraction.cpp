#include <iostream>        //för cout m.m
#include <fstream>         //läsa&skriva till filer
#include <cstdlib>         //för abort()-funktionen.
#include <ctime>           //för att använda systime som seed
#include <vector>          //för vektorklassen
#include <cmath>           //för pow-, sqrt-funktionerna
#include "nr3.h"           //Random number generator from nummerical recipes 3 ed.
#include "ran.h"           //obs, inkludera inte denna vid g++-kompileringen

//Main-function usage: ./program [outputN.dat] [0/1]
//where the two arguments are optional. The default is to produce an output to "out.dat", with values with linear spacing (logscale off =0). If you want logscale=ON, but still use default "out.dat", you must still specify it: ./program out.dat 1. 

// use ssh:
// $ ssh -n 'cd /to/path; nice -19 ./program <input.dat > output1.dat'

//some hardcoded stuff:
#define INTERACTION_ENERGY 0    //Interaction energy, > 0 ! (must be non-negative)
#define VERBOSE 1               //print detailed information about SuperInteraction()
#define CHECK_VACANCY_ON 0      //turn on/off the time-consuming CheckVacancyMatrix()

#define DISTRIBUTION     0    //0=uniform, 1=exponential,         | (main)
//                              2=power-law, 3=nakazato           | (main)
#define JUMPRATE_TRACER  1.0    //                                  | (main)
#define JUMPRATE_CROWDER 1.0  //(if nakazato-distribution)        | (main)
#define FIXBOUNDARY      0    //1=fix wall, 0=periodic boundary   | (main ->class)

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

  int RutorX,RutorY,RutorZ;                   //gittrets storlek
  int N, dim;                                 //antal partiklar och dimention (=1 || 2 || 3)
  float density;                              //particle density
  int X_0, Y_0, Z_0;                          //Initial positonion for the fist particle
  double time_sum;                            //Sum of time for each move
  int ConvertMuToParticle(int mu, int& n);
  void PartialSum(void);
  void MoveAndBoundaryCheck(int,int,bool);
  int VacancyCheck(int, int,int,int);         //No double Occupancy.
  int VacancyCheck2(int, int,int,int);        //Same as above, but improved. Use this one instead.
  void StdErr(int, int);                      //beräknar standard deviation
  void BuildCluster(int,  vector<int>&, int, int, int);
  vector< vector<int> > dx, dy, dz;           //alla pos. för vit part. behövs för sdrerr
  vector< vector<double> > dr;

  vector< vector< vector<int> > > vacancy;    //Needed in my more efficient "vacancyCheck2()"

  int seed1, seed2;

  void SuperInteractionCode(int, int, bool);
  void CountNeighbours(int, int, int, vector<int>&, bool);
  void CalculateExtent(vector<int>&, int, int&, int&); 
  void BuildBlockingCluster(vector<int>&, vector<int>, int, bool);
  void BuildCluster2(int, vector<int>&, double, bool);




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

  //Bugfinding and Testing functions
  //--------------------------------
  float Ergodicity(int,int,int, bool);        //Korntroll: ergodicitetens teoretiska värde
  void CheckVacancyMatrix(float);             //Simple controll function, to find bugs
  void ControlVacancyCheck3(int);             //Simple controll function, to find bugs
  void PrintPositionToScreen(int);            //Print position of particle label (int), or vector
  void PrintPositionToScreen(vector<int>);    //containing particle labels, to screen. (Test-function)


public:
  Partiklar(int, int, int, int, bool);
  void Place(void);
  void Move1(bool);
  void Move2(bool);
  void Lagra(vector<int>&, vector<int>&, vector<int>&, vector<double>&);
  void Save(int, char[]);
  void set_dt(void);                           //för initiering av tids-vektorn
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


void Partiklar::CheckVacancyMatrix(float ErrorCode){
  //This is just a small function to check that the vacancy matrix is correct
  //ie. vacant sites are marked as "-1" and occupied are marked with the 
  //correct particle label.
  //ErrorCode = code to output to easier find the place in the code
  
  if(CHECK_VACANCY_ON){
    bool error = false;
    int label;
    int count_vacant_sites = 0;


    //Check that occupied site is correct
    //and the number of vacant sites is correct
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
}

void Partiklar::PrintPositionToScreen(int particle){
//Print position of particle label (int), or vector
//containing particle labels, to screen. (Test-function)
  int m;
  m=particle;
  cout<<"Particle "<<m<<" X:Y:Z "
      <<Xpos[m]<<":"<<Ypos[m]<<":"<<Zpos[m]<<endl;
}          

void Partiklar::PrintPositionToScreen(vector<int> particle){
//Print position of particle label (int), or vector
//containing particle labels, to screen. (Test-function)
  int m;
  for(int i=0; i < particle.size(); i++){
    m=particle[i];
    cout<<"Particle "<<m<<" X:Y:Z "
        <<Xpos[m]<<":"<<Ypos[m]<<":"<<Zpos[m]<<endl;
  }
  cout<<"Total size: "<<particle.size()<<endl;
}



void Partiklar::ControlVacancyCheck3(int n){
  //This is just a control/bug finding algorithm, based on VacancyCheck2()
  //It's run before setting the new site as occupied, to make sure it was vacant 
  //before the move...

  if (n < N && 0 <= n){
    if (vacancy[Xpos[n]][Ypos[n]][Zpos[n]] != -1){  //if occupied
      cout <<"Occupied ("<<Xpos[n]<<","<<Ypos[n]<<","<<Zpos[n]<<")"<<endl;
      cout <<"Two particles on same site."<<endl;
      abort();
    }
  }
  else{
    cout <<"Error, accessing invalid particle "<<endl;
    abort();
  }
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

  //due to laziness, we multiply the "true" difusion constant here
  //with 2*dimenstion, and return that. Keep that in mind when plotting.
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
      cout <<" PLACERADE \n partikel n="<<q
           <<" X="<<Xpos[q]<<" Y="<<Ypos[q]<<" Z="<<Zpos[q]<<endl;


  //IF We use VacancyCheck2()-function:
  //-----------------------------------
  //initiate the vacancy-vector with "-1". (since
  // occupied site = paritcle label >= 0).
  //Must reset this matrix for each ensamble!
  vector< vector< vector<int> > >
    Vacancy(RutorX+1, vector< vector<int> >(RutorY+1, vector<int>(RutorZ+1,-1)));
  vacancy=Vacancy;
  //Mark occupied sites as n:
  //note, no change will ever occur to the 0-elements. (vacancy[0][0][0])
  for (n=0; n<N; n++){
    vacancy[ Xpos[n] ][ Ypos[n] ][ Zpos[n] ]=n;
  }

  CheckVacancyMatrix(316846546); //TEST that it's correct!

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
        break;                                              //(first coordinate is 1, 0=outside) 
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
      //NOTE: Used to have a Vacancy-check function here...
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
        return returnvalue;        //Exit this method/function
      }
      else
        if (testOnOff) cout <<"OK"<<endl;
    }
    if (testOnOff) cout <<"Vacant, particle "<<n<<" X = "<<Xpos[n]<<", Y = "
                        <<Ypos[n]<<", Z = " <<Zpos[n]<<endl;
    return returnvalue;
    //(the returnvalue could be used to have a do-move-while(VacancyCheck==1) but
    //if so, we need to move this function to be directly in the Move_()-func.)
  }
  else{
    cout <<"Error, accessing invalid particle "<<endl;
    abort();
  }
}



int Partiklar::VacancyCheck2(int n, int Xold, int Yold, int Zold){
  //NOTE: This is a faster version than the old VacancyCheck(), since that one
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
      returnvalue=1;                             //Exit method/function

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
    else{                                       //if vacant, update vacancy-matrix 
      vacancy[Xold][Yold][Zold]=-1;             //old position is now vacant
      vacancy[Xpos[n]][Ypos[n]][Zpos[n]]=n;     //new position is now occupied
      if (testOnOff) cout <<"OK"<<endl;
    }
    if (testOnOff) cout <<"Vacant, particle "<<n<<
                     " X = "<<Xpos[n]<<", Y = "<<Ypos[n]<<", Z = " <<Zpos[n]<<endl;
    return returnvalue;
    //(the returnvalue could be used to have a do-move-while(VacancyCheck==1) but
    //if so, we need to move this function to be directly in the Move_()-func.)
  }
  else{
    cout <<"Error, accessing invalid particle "<<endl;
    abort();
  }
}

//-----------------------------------------*----------------------------------------------
//----------------------------------------*-*---------------------------------------------
//---------------------------------------BEGIN--------------------------------------------


void Partiklar::SuperInteractionCode(int n, int direction, bool boundary ){
  //boundary==0 periodic, boundary==1 fix
  //Note: this function sould not be called from MoveAndBoundaryCheck, but
  //rather from the function that determines n, namley Move2()
  
  //TODO Flytta in CheckVacancyMatrix in i verbose-satserna, när koden funkar

  //This is a switch to print the events to the screen, 
  //to check it's actually doing what we want it to do.
  bool verbose = VERBOSE;

  if(verbose){
    string strDir[]={"X+1","X-1","Y+1","Y-1","Z+1","Z-1"};
    cout <<"-------------------------"<<endl
         <<"New run of SuperInteractionCode() to move particle n="
         <<n<<", in direction: "<<strDir[direction]<<endl;
    PrintPositionToScreen(n);
    cout <<"-------------------------"<<endl<<endl;

    //print all occupied sites to screen:
    cout<<"Lattice configuration is:"<<endl;
    for(int i=0; i<Xpos.size();i++){
      cout<<"particle "<<i<<": "<<Xpos[i]<<", "<<Ypos[i];
      if(dim==3) cout<<", "<<Zpos[i]; cout<<endl;
    }
    cout<<endl;

    //TEST
    //the purpose is just to give each run a unique label,
    //it helps when debugging. 
    static Ran Random2(69865);
    double Test_Random_number;
    Test_Random_number=Random2.doub();
    cout<<"random label for this run: "<<(double)Test_Random_number<<endl;

    if(Test_Random_number<=0.883078 && Test_Random_number>=0.883076) {
      cout<<"Bug number 2";
    }

    if(Test_Random_number<=0.127939 && Test_Random_number>=0.127935) {
      cout<<"Bug number 3";
    }

    if(Test_Random_number<=0.842700 && Test_Random_number>=0.842690) {
      cout<<"Bug number 4";

      //0.690199


    }

  }
  
  CheckVacancyMatrix(0.10);

  vector<int> Cluster;            //store number of particles in Cluster here. 
  CountNeighbours(Xpos[n],Ypos[n],Zpos[n], Cluster, boundary);
  double V=INTERACTION_ENERGY;    //interaction energy /(k_B*T) (V > 0)

  static Ran Random(56);
  double Rand_numb;

  if  (Cluster.size() == 0){
    //accept move
    vacancy[Xpos[n]][Ypos[n]][Zpos[n]]=-1;      //set as vacant
    MoveAndBoundaryCheck(n,direction,boundary); //move
    ControlVacancyCheck3(n);                    //TEST Function, checks that vacancy(n)..=-1
    vacancy[Xpos[n]][Ypos[n]][Zpos[n]]=n;       //new position occupied

    CheckVacancyMatrix(0.11); 

    if(verbose){
      cout<<"Moved single particle (no cluster/neighbors) n to "<<endl;
      PrintPositionToScreen(n);
      cout <<endl;
    }

  }
  else{
    //start with adding particle n to the cluster:
    //(or else it might not be part of the cluster,
    // due to Boltzmann prob, which would be wrong!)
    
    //    Cluster.push_back(n);
    //    Gör detta i BuildCluster2() istället

    CheckVacancyMatrix(5874);  //TEST OK
    BuildCluster2(n,Cluster,V,boundary);

    //TEST---------
    bool testinglite=true;
    for(int i=0; i<Cluster.size();i++){
      if (Cluster[i]==n)
        testinglite=false;
    }
    if(testinglite){
      cout <<"n not part of Cluster! Error 1"<<endl;
      abort();
    }//-------------


    if(verbose){
      cout<<"Particle "<<n<<" part of cluster:"<<endl;
      PrintPositionToScreen(Cluster);
      cout <<endl;
    }

    //Calculate the friction constant (k^(-1))
    double SUM_Cluster=0;
    double SUM_Particle=0;

    for (int i=0; i<Cluster.size(); i++){
      int m = Cluster[i]; //gives correct particle-labels
      SUM_Cluster = SUM_Cluster + 1.0/jmpRateXright[m];
      SUM_Cluster = SUM_Cluster + 1.0/jmpRateXleft [m]; 
    }
    SUM_Particle = SUM_Particle + 1.0/jmpRateXleft [n]; 
    SUM_Particle = SUM_Particle + 1.0/jmpRateXright[n];

    if(dim>=2){
      for (int i=0; i<Cluster.size(); i++){
        int m = Cluster[i]; //gives correct particle-labels
        SUM_Cluster = SUM_Cluster + 1.0/jmpRateYright[m]; 
        SUM_Cluster = SUM_Cluster + 1.0/jmpRateYleft [m]; 
      }
      SUM_Particle = SUM_Particle + 1.0/jmpRateYleft [n]; 
      SUM_Particle = SUM_Particle + 1.0/jmpRateYright[n]; 
    }
    if(dim==3){
      for (int i=0; i<Cluster.size(); i++){
        int m = Cluster[i]; //gives correct particle-labels
        SUM_Cluster = SUM_Cluster + 1.0/jmpRateZup   [m]; 
        SUM_Cluster = SUM_Cluster + 1.0/jmpRateZdown [m]; 
      }
      SUM_Particle = SUM_Particle + 1.0/jmpRateZup   [n]; 
      SUM_Particle = SUM_Particle + 1.0/jmpRateZdown [n];
    }

    testinglite=true; //TEST
    for(int i=0; i<Cluster.size();i++){
      if (Cluster[i]==n)
        testinglite=false;
    }
    if(testinglite){
      cout <<"n not part of Cluster!"<<endl;
      abort();
    }

    double P = SUM_Particle / ( SUM_Cluster );
    //Note, SUM_Particle is part of SUM_Cluster, since n is part
    //of Cluster-vector. This is correct, since then P < 1.
    //or P = 1 if n is not part of a cluster. 

    if(1<P){ //TEST
      cout <<"SUM_Particle"<<SUM_Particle<<endl;
      cout <<"SUM_Cluster"<<SUM_Cluster<<endl;
      cout <<"P fel.P="<< P <<endl;
      abort();
    }

    Rand_numb = Random.doub();
    if (Rand_numb > P){
      //Don't move the cluster if (P<r):
      //Nothing needs to be done, vacancy-matrix unchanged. If 
      //cluster.size()=1, then P=1, and we will always move the particle.
      if(verbose) cout<<"Cluster to heavy, particle "
                      <<n<<" stays"<<endl<<endl;
    }
    else{
      //move whole Cluster.

      if(verbose) cout<<"Move the whole cluster, does it collide? ";

      CheckVacancyMatrix(645601); //OK hit! 

      //We do a trial-and-error move with the whole cluster, since
      //it's the easiest way to check for obstructions/crowders.
      //store the old positions, temporarily.
      vector<int> XoldCluster, YoldCluster, ZoldCluster;
      for (int i=0; i<Cluster.size(); i++){
        int c = Cluster[i];
        XoldCluster.push_back(Xpos[c]);
        YoldCluster.push_back(Ypos[c]);
        ZoldCluster.push_back(Zpos[c]);
        //This loop includes particle n.

        //set their old positions as Vacant, perhaps just temporarily:
        vacancy[Xpos[c]][Ypos[c]][Zpos[c]]= -1;

        //Move all particles in the cluster:
        MoveAndBoundaryCheck(c,direction,boundary);
      }

      //Check their new positions against double occupancy. Be sure to exclude their old
      //coordinates from the vacancy-matrix. (this is already done by seting them as vacant)

      //CheckVacancyMatrix(680001); //Klipper har vilket ar ratt...

      vector<int> Blocking;
      bool occupied = false;
      for (int i=0; i<Cluster.size(); i++){
        int c = Cluster[i];
        //Now all particles in Cluster have new positions
        //but their new sites are not yet marked as occupied (=-1)
        if (vacancy[Xpos[c]][Ypos[c]][Zpos[c]] != -1){
          //store which particles are blocking the cluster;
          Blocking.push_back(vacancy[Xpos[c]][Ypos[c]][Zpos[c]]);
          occupied = true;
        }
      }

      //TEST
      for (int i=0; i < Blocking.size(); i++ ){
        for (int j=0; j < Cluster.size(); j++ ){
          if(Cluster[j]==Blocking[i]){
            cout << "Error code 000: Same particle in both vectors!"<<endl;
            abort();
          }
        }
      }

      if(verbose){
        if (occupied) cout <<"YES"<<endl<<endl;
        else cout << "NO" <<endl<<endl;
      }

      if (!occupied){
        //no double occupancies, accept move, uppdate vacancy-matrix:

        if(Blocking.size()!=0){
          cout <<"Blocking fel skall vara 0.";
          abort();
        }

        for (int i=0; i<Cluster.size(); i++){
          int c = Cluster[i];
          //Now all particles in Cluster have new positions, including n
          vacancy[Xpos[c]][Ypos[c]][Zpos[c]] = c; 
        } 

        if(verbose){
          cout<<"No collisions, moved cluster to:"<<endl;
          PrintPositionToScreen(Cluster);
          cout <<endl;
        }
          
        CheckVacancyMatrix(35648); //BUG when fix boundary!

      }
      else{ //Ouch, our cluster is colliding with some shit...
        //TODO: tänk på att vacancy är rätt märklig nu, med XoldCluster == -1 osv...


        //this is an addition made in July: For fucks sake! I need to restor the vacancy-matrix
        //now, since CountNearestneighbors() and BuildCluster...() depends on it being correct!
        //Reset the "trial-and-error" move (for now), ie. move the cluster back! Pronto!
        for(int i=0; i<Cluster.size(); i++){
          int c = Cluster[i];
          Xpos[c] = XoldCluster[i];
          Ypos[c] = YoldCluster[i];
          Zpos[c] = ZoldCluster[i];

          vacancy[Xpos[c]][Ypos[c]][Zpos[c]] = c;
        }
        CheckVacancyMatrix(9080); //This must be correct now!
        //End of July addition (but must compensate for this further down when we actually move)
        //NOTE: With the trial-and-error, we got what we wanted, namley the Blocking[]-vector
        
        //Kor denna innuti BuildBlockingCluster istallet!
        //CalculateExtent(Cluster,direction,min,max);

        if(Blocking.size()==0){ //TEST
          cout <<"Blocking fel skall INTE vara 0 err2.";
          abort();
        }

        //Construct a vector that contains all particles that are blocking the path.
        BuildBlockingCluster(Blocking,Cluster,direction,boundary);

        if(Blocking.size()==0){
          cout <<"Blocking fel skall INTE vara 0. err3"<<endl;
          abort();
        }

        if(verbose){
          cout<<"The blocking cluster is:"<<endl;
          PrintPositionToScreen(Blocking);
          cout <<endl;
        }

        double SUM_Blocking=0;

        for (int i=0; i < Blocking.size(); i++){
          int b = Blocking[i]; 
          SUM_Blocking = SUM_Blocking + 1.0/jmpRateXright[b];
          SUM_Blocking = SUM_Blocking + 1.0/jmpRateXleft [b]; 
        }

        if (dim>=2){
          for (int i=0; i < Blocking.size(); i++){
            int b = Blocking[i]; 
            SUM_Blocking = SUM_Blocking + 1.0/jmpRateYright[b]; 
            SUM_Blocking = SUM_Blocking + 1.0/jmpRateYleft [b]; 
          } 

        }
        if (dim==3){
          for (int i=0; i < Blocking.size(); i++){
            int b = Blocking[i]; 
            SUM_Blocking = SUM_Blocking + 1.0/jmpRateZup   [b]; 
            SUM_Blocking = SUM_Blocking + 1.0/jmpRateZdown [b]; 
          }

        }

        double P = (double) SUM_Cluster / (SUM_Blocking + SUM_Cluster);
        
        static Ran something(681);
        double slump;
        slump = something.doub();
        if(slump <= P ){
          /*Move Cluster and All the Blocking particles. But beware 
            when we move the cluster and update the vacancy matrix we 
            must make sure it works for a situation where for instance
            everything gets shifted, around the boundary:
            ---------------------
            | 1 | 2 | 3 | 4 | 5 |
            ---------------------
            if particle 4 is the single particle cluster that pushes
            at the blocking cluster [1 2 3 5], then we must get 
            (for X+1 move):
            ---------------------
            | 5 | 1 | 2 | 3 | 4 |
            ---------------------
            and have the vecancy-matrix be correct!
           */
          
          CheckVacancyMatrix(111156);

          if(Blocking.size()==0){//TEST
            cout <<"Blocking fel skall INTE vara 0. err5";
            abort();
          }


          //Set all old positions as vacant:
          for (int i=0; i < Blocking.size(); i++){
            int b = Blocking[i];
            vacancy[Xpos[b]][Ypos[b]][Zpos[b]]=-1;
          }
          for (int i=0; i < Cluster.size(); i++){
            int c = Cluster[i];
            vacancy[Xpos[c]][Ypos[c]][Zpos[c]]=-1;
          }

          //Move all particles, to get new positions in Xpos, Ypos,..
          for (int i=0; i < Blocking.size(); i++){
            int b = Blocking[i];
            MoveAndBoundaryCheck(b,direction,boundary);
          }
          for (int i=0; i < Cluster.size(); i++){
            int c = Cluster[i];
            MoveAndBoundaryCheck(c,direction,boundary);
          }

          //Update the vacancy-matrix
          for (int i=0; i < Blocking.size(); i++){
            int b = Blocking[i];
            vacancy[Xpos[b]][Ypos[b]][Zpos[b]]=b;
          }
          for (int i=0; i < Cluster.size(); i++){
            int c = Cluster[i];
            vacancy[Xpos[c]][Ypos[c]][Zpos[c]] = c;
          }
          
          if(verbose){
            cout <<"Move the cluster and the blocking particles."<<endl
                 <<"Blocking Particles:"<<endl;
            PrintPositionToScreen(Blocking);
            cout <<"Original Cluster:"<<endl;
            PrintPositionToScreen(Cluster);
            cout <<endl;
          }

          //check for errors
          CheckVacancyMatrix(32229);
        }
        else{
          //Move everyting back, uppdate vacancy Matrix.
        
          //SINCE THE JULY ADDITION, WE HAVE ALREADY MOVED BACK THE 
          //TRIAL-AND-ERROR MOVE, SO WE NEED NOT DO ANYTHING...
          
          // for (int i=0; i < Cluster.size(); i++){
          //   int c = Cluster[i];
          
          //   //Keep in mind: label stored in Cluster[i]
          //   //has its old position stored in X/Y/ZoldCluster[i]
          //   Xpos[c] = XoldCluster[i];
          //   Ypos[c] = YoldCluster[i];
          //   Zpos[c] = ZoldCluster[i];          
          
          //   vacancy[Xpos[c]][Ypos[c]][Zpos[c]] = c;
          // }

          if(verbose){
            cout <<"The blocking cluster is too heavy, no move is made"
                 <<endl<<"the cluster is unmoved:"<<endl;
            PrintPositionToScreen(Cluster);
            cout <<endl;
          }

          CheckVacancyMatrix(0.39);
        
        }
        CheckVacancyMatrix(0.40);
        
      }//END -- colliding Cluster
      CheckVacancyMatrix(0.41);
      
    }//if move whole cluster-END.
    CheckVacancyMatrix(0.42); 
    
  }//if nearestNeighbor != 0 END
  CheckVacancyMatrix(0.43); 
}


void Partiklar::CalculateExtent(vector<int>& Cluster,int Direction, 
                                int& Min, int& Max){
  //Function to calculate the maximum and minimum coordinates of the cluster,
  //perpendicular to the jumping-direction. Only Works for 2D so far...
  //(and 1D? hmm, ingen utbredning vinkelrät mot flyttriktning...)

  //Do I need to consider boundary periodic/fix? I don't think so, since
  //other functions... .... hmmmm

  //BUG: Potential Bug, when moving a wrapped cluster across the edge...

  int Dir;

  if ( Direction == 0 || Direction == 1) Dir = 0; //Move in X
  if ( Direction == 2 || Direction == 3) Dir = 1; //Move in Y
  if ( Direction == 4 || Direction == 5) Dir = 2; //Move in Z

  switch(Dir){

  case 0: //perpendicular to X
    Max = 1;
    Min = RutorY;
    for (int i=0; i < Cluster.size(); i++){
      int m = Cluster[i];

      //Max2 = 1;
      //Min2 = RutorZ;

      if (Ypos[m] < Min ) Min = Ypos[m];
      if (Ypos[m] > Max ) Max = Ypos[m];

      //if (Zpos[m] < Min2 ) Min2 = Zpos[m];
      //if (Zpos[m] > Max2 ) Max2 = Zpos[m];
    }
    break;


  case 1: //perpendicular to Y
    Max = 1;
    Min = RutorX;
    for (int i=0; i < Cluster.size(); i++){
      int m = Cluster[i];
    
      //Max2 = 1;
      //Min2 = RutorZ;

      if (Xpos[m] < Min ) Min = Xpos[m];
      if (Xpos[m] > Max ) Max = Xpos[m];

      //if (Zpos[m] < Min2 ) Min2 = Zpos[m];
      //if (Zpos[m] > Max2 ) Max2 = Zpos[m];
    }
    break;

  case 2: //perpendicular to Z
    cout <<"Move perpendicular to Z, error in 2D!"<<endl;
    //this code is not yet operational in 3D!
    abort();
    Max = 1;
    Min = RutorX;
    for (int i=0; i < Cluster.size(); i++){
      int m = Cluster[i];

      //Max2 = 1;
      //Min2 = RutorY;

      if (Xpos[m] < Min ) Min = Xpos[m];
      if (Xpos[m] > Max ) Max = Xpos[m];

      //if (Ypos[m] < Min2 ) Min2 = Ypos[m];
      //if (Ypos[m] > Max2 ) Max2 = Ypos[m];
    }
    break;
  }//switch end

  //Check if Cluster is wraped around the periodic boundary
  //-------------------------------------------------------
  bool wrapped=false;
  if (Max==RutorY && Min==1 && Dir==0) wrapped = true; //move perpendicular to X
  if (Max==RutorX && Min==1 && Dir==1) wrapped = true; //move perpendicular to X
  //if (Max==RutorZ && Min==1 && Dir==3) wrapped = true; //not working, but Dir!=3 now...

    if(wrapped){
    //We define a convention:
    //If the cluster is wrapped around the periodic boundary:
    //its borders will be 1<partOne<MaxX && MinX<partTwo<RutorX
    //which means that MaxX < MinX, which is how we later will 
    //know that the cluster is devided in two.

    int MaxX,MaxY,MinX,MinY;
    bool step=true;

    if(Dir==0){ //If move in X-direction
      //Scan cluster from 1 --> RutorY and stop
      // for-loop once step = flase,
      for(int i=1; step && i<= RutorY; i++){ 
        if (RutorY<i){
          cout << "Ändra gräns < error kod 44"<<endl;
          abort();
        }

        step=false; //stop for-loop if no cluster-particles at this coordinate (=i)
        for(int j=0; j<Cluster.size() && !step; j++){  
          //having "&& !step" in the condition above is just to skip rest
          //of this time-consuming loop over the cluster particles. once we know
          //the cluster has at least one particle here, we dont need to run it more   
          if(Ypos[Cluster[j]]==i){
            //OK, allow another step upwards in Y
            //to check if the cluster extends this far
            step=true;
            MaxY=i;
          }
        }
      }

      step=true;  //kan nog ta bort denna senare...(när allt funkar)
      //Scan cluster from RutorY --> 1
      for(int i=RutorY; step && i>0; i--){
        if (i<1){
          cout << "Ändra gräns <= error kod 45"<<endl;
          abort();
        }
        step=false;
        for(int j=0; j<Cluster.size() && !step; j++){
          if(Ypos[Cluster[j]]==i){
            //OK, allow another step downwards in Y
            //to check if the cluster extends this far
            step=true;  
            MinY=i;
          }
        }
      }
      if (MinY==1 && MaxY==RutorY) 
      //This is quite possible, since it can be "twisted" like a vortex or a swirl
        cout <<"Cluster spanns entire lattice!"<<endl;


      Max=MaxY;
      Min=MinY;
    }//"Move in X-direction" -END
   

    if(Dir==1){ //If move in Y-direction

      //Scan cluster from 1 --> RutorX
      for(int i=1; step && i<= RutorX; i++){ //Osäker på: <= eller < ???
        if (RutorX<i){
          cout << "Ändra gräns < error kod 42"<<endl;
          abort();
        }
        step=false;
        for(int j=0; j<Cluster.size() && !step; j++){
          if(Xpos[Cluster[j]]==i){
            //OK, allow another step to the right in X
            //to check if the cluster extends this far
            step=true;  
            MaxX=i;
          }
        }
      }

      step=true;  //kan nog ta bort denna senare...(när allt funkar)
      //Scan cluster from RutorX --> 1
      for(int i=RutorX; step && i>0; i--){
        step=false;
        if (i<1){
          cout << "Ändra gräns <= error kod 43"<<endl;
          abort();
        }

        for(int j=0; j<Cluster.size() && !step; j++){
          if(Xpos[Cluster[j]]==i){
            //OK, allow another step to the left in X
            //to check if the cluster extends this far
            step=true;  
            MinX=i;
          }
        }
      }
      if (MinX==1 && MaxX==RutorX) 
        cout <<"Cluster spanns entire lattice!"<<endl;
      Max=MaxX;
      Min=MinX;   
    }//"Move in Y-direction" -END

    
  }//if "maybe wrapped around periodic"-END
}

 
void Partiklar::BuildBlockingCluster(vector<int>& BuildBlockingVector, vector<int> Cluster,
                                     int direction, bool FixBound){
  //"Cluster" is the moving cluster that's colliding with BuildBlockingVector
  //"BuildBlockingVector" contains at least one of the blocking particles. 
  //We add the rest of the blocking particles to this, and store in unique

  vector<int> unique; 


  //TEST to see if Cluster and BuildBlockingVector contains the same particles
  //if so --> check that vacancy matrix hasn't be (incorrectley) changed
  for (int i=0; i < BuildBlockingVector.size(); i++ ){
    for (int j=0; j < Cluster.size(); j++ ){
      if(Cluster[j]==BuildBlockingVector[i]){
        cout << "Error: Same particle in both vectors!"<<endl;
          abort();
      }
    }
  }

  //First check the min,max coordinates for the Cluster
  //perpendicular to its direction of movement: 
  int Min, Max;  //these are given as references in the func. below:
  CalculateExtent(Cluster, direction, Min, Max);

  //TEST:
  CheckVacancyMatrix(39103);

  if(dim==3){
    cout <<"BuildBlockingCluster does not work in 3D (yet)"<<endl;
    abort();
  }

  //if we added a particle to the vector, run again to check its neighbors
  bool CheckAgain;

  unsigned int counter=0;  //debugging variable, see how many times the loop runs
  do{
    counter++;
    CheckAgain=false;

    int limit = BuildBlockingVector.size();
    for(int i=0; i< limit; i++){
      int m = BuildBlockingVector[i];
      CountNeighbours(Xpos[m],Ypos[m],Zpos[m], BuildBlockingVector, FixBound);
      //BuildBlockingVector will store more particles each 
      //iteration, but only loop over the old particles ( i.e. i < limit)
    }

    for(int j = 0; j < BuildBlockingVector.size(); j++){
      bool SingleValued = true;
      int s = BuildBlockingVector[j];

      //Check if it's stored in unique[] already:
      for(int k = 0; k < unique.size(); k++){
        if (s == unique[k])
          SingleValued = false;
        //NOTE: Never enters this loop on the first run (unique.size=0)
      }
      //Check that particle s is not part of "Cluster"
      if(SingleValued){ //no need to run it if it's already "false"
        for(int i=0; i<Cluster.size();i++){
          if(s==Cluster[i])  SingleValued=false;
        }
      }

      if(SingleValued){//BUG: After this Block unique.size == 0
            
        //Check if Min & Max are infact wrapped around the boundary!
        bool normal=true;                    
        if(Max < Min) normal=false;  //according to definition of Max/Min in CalculateExtent()

        if (normal && (direction == 0 || direction == 1)){ //Move in X
          if (Ypos[s] <= Max && Min <= Ypos[s]){
            unique.push_back(BuildBlockingVector[j]);
            //New particels added to the cluster, run the loop again, 
            //to check the neighbors of this new particle
            CheckAgain=true;
          }
        }
        if (!normal && (direction == 0 || direction == 1)){ //Move in X
          if ( (1 <= Ypos[s] && Ypos[s] <=Max ) || ( Min <= Ypos[s] && Ypos[s] <= RutorY)){
            unique.push_back(BuildBlockingVector[j]);
            //New particels added to the cluster, run the loop again, 
            //to check the neighbors of this new particle
            CheckAgain=true;
          }
        }
        if (normal && (direction == 2 || direction == 3)){ //Move in Y
          if (Xpos[s] <= Max && Min <= Xpos[s]){
            unique.push_back(BuildBlockingVector[j]);
            //New particels added to the cluster, run the loop again, 
            //to check the neighbors of this new particle
            CheckAgain=true;    
          }
        }
        if (!normal && (direction == 2 || direction == 3)){ //Move in Y
          if ( (1 <= Xpos[s] && Xpos[s] <=Max ) || ( Min <= Xpos[s] && Xpos[s] <= RutorX)){
            unique.push_back(BuildBlockingVector[j]);
            //New particels added to the cluster, run the loop again, 
            //to check the neighbors of this new particle
            CheckAgain=true;
          }
        } 
        if ( direction == 4 || direction == 5){ //Move in Z
          cout <<"I'm in 2D!, err-code: 769"<<endl;
          //this code is not fully functional in 3D....
          abort();
        }

        /*
        //We have now limited the Blocking cluster on the sides, But not backwards!
        //Exclude any particle that's a part of the Cluster-vector, IF a new
        //particle was added
        if(CheckAgain){
          for(int i=0; i<Cluster.size();i++){
            //cout <<endl<<"Cluster.size(): "<<Cluster.size()<<endl;
            //cout <<"i="<<i<<endl;
            //cout << "unique.back: "<<unique.back()<<endl;
            //cout << "unique.size(): "<<unique.size()<<endl<<endl;
            if(unique.back()==Cluster[i]){
              unique.pop_back(); //delete the added particle
              //cout << "unique.size() (after pop_back): "<<unique.size()<<endl;
              
              //Since we removed the particle, we should not check the neighbors
              //of it in a new run!
              CheckAgain=false;  
            }
          }
        
          //TEST
          if(unique.size()==0){
            cout << "Unique.size() err code 1214: "<<unique.size()<<endl;
            abort();       
          }
          }*/

        //TEST
        // if(unique.size()==0){ //BUG! eller?
        //   cout << "unique.size err code 1215: "<<unique.size()<<endl
        //        << "CheckAgain-value: "<< CheckAgain << endl;
        //     abort();       
        // }
      }
  
      //TEST
      // if(unique.size()==0){
      //   cout << "Unique.size() err 1216:"<<unique.size()<<endl;
      //   abort();
      // }


    }

    //TEST
    if(unique.size()==0){
      cout << "Unique.size() err 1217:"<<unique.size()<<endl;
      abort();       
    }
    //This code works even if BuildBlockingVector-vector contains 
    //the same particle multiple times, new and/or old. I think... 


    //TEST: 
    if (unique.size() > BuildBlockingVector.size() ){
      cout <<"Fel i unique tror jag (655654)"<<endl;
      abort();
    }
    
    BuildBlockingVector = unique;
    
  }while(CheckAgain);

  /*
  //TEST:
  cout << "kollar BuildBlockingVector.size():"<<BuildBlockingVector.size()<<endl;
  for(int j=0; j < BuildBlockingVector.size(); j++){
  int m= BuildBlockingVector[j];
  cout <<"n = "<<m<<" ("<<Xpos[m]<<","<<Ypos[m]<<","<<Zpos[m]<<")"<<endl;    
  }

  //TEST
  for (int i=0; i < N; i++){
  if (Xpos[i]==0 || Ypos[i]==0 || Zpos[i]==0){
  cout <<"Utanför brädet!"<<endl;
  abort();
  }
  }
  */

  //TEST to see if Cluster and BuildBlockingVector contains the same particles
  //if so --> check that vacancy matrix hasn't be (incorrectley) changed
  for (int i=0; i < BuildBlockingVector.size(); i++ ){
    for (int j=0; j < Cluster.size(); j++ ){
      if(Cluster[j]==BuildBlockingVector[i]){
        cout << "Error: Same particle in both vectors (err. location: 2)"<<endl;
          abort();
      }
    }
  }

}

void Partiklar::BuildCluster2(int n, vector<int>& BuildClusterVector, 
                              double Exponent, bool boundary){
  //Build a vector with particles that form the Cluster:
  //-------------------------------------------------
  //each neighbor has a probability given by the Boltzman factor:
  //    exp(-1.0*Exponent)
  //to belong to the cluster.

  //cout<<"Interaction Strength: "<< Exponent << endl;

  vector<int> unique;
  vector<int> NewParticleAdded;
  static Ran RandomNUM(42);
  double random_number;

  //n included with 100 % probability:
  unique.push_back(n);

  bool testinglite;  //test-variable to find potential errors in the code
  bool CheckAgain;   //true if a new particle was added to the cluster

  do{
    CheckAgain=false;

    int limit = BuildClusterVector.size();
    for(int i=0; i< limit; i++){
      int m=BuildClusterVector[i];
      CountNeighbours(Xpos[m],Ypos[m],Zpos[m], BuildClusterVector, boundary);
      //BuildClusterVector-vector will store more particles each 
      //iteration, but only loop over the first initial set (ie. i < limit)
    }

    for(int j=0; j < BuildClusterVector.size(); j++){
      bool SingleValued=true;
 
      //Check if it's stored in unique[] allready:
      for(int k=0; k < unique.size(); k++){
        if (BuildClusterVector[j] == unique[k])
          SingleValued=false; 
      }
      if(SingleValued){
        //Add the new particle, with Boltzmans consent...
        random_number = RandomNUM.doub();
        if (random_number > exp(-1.0*Exponent)){
          unique.push_back(BuildClusterVector[j]);
          NewParticleAdded.push_back(BuildClusterVector[j]);
          //New particels added to the cluster, run the loop again
          //for particles in NewParticleAdded to check the neighbors
          // of this new particle
          CheckAgain=true;
        }
      }
    }
     
    BuildClusterVector = NewParticleAdded;
    NewParticleAdded.clear();

  }while(CheckAgain);

  //Now save the final vector (cluster) and return it as a reference:
  BuildClusterVector = unique;


  if(false){ //Inaktiverad just nu
    //TEST:
    cout << "kollar Nearest Neighbour.size(): "<<BuildClusterVector.size()<<endl;
    for(int j=0; j < BuildClusterVector.size(); j++){
      int m= BuildClusterVector[j];
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


void Partiklar::CountNeighbours(int X, int Y, int Z, vector<int>& Count, bool boundary){
  //Each element in the vector "Count" is a particle, and the 
  //value of the element is the particle label. 
  //If the particle is "alone" the vector Count will be empty.
  //boundary=0 periodic, =1 fix
  
  bool periodic;
  
  if (!boundary) periodic = true;
  else  periodic = false;

  CheckVacancyMatrix(7898); //TEST Just to check it's correct.

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


//---------------------------------------E-N-D------------------------------------------
//----------------------------------------*-*-------------------------------------------
//-----------------------------------------*--------------------------------------------


void Partiklar::Move2(bool boundfix){
  CheckVacancyMatrix(981);

  if (partial_sum.empty()){
    cout <<"Partial-sum vector has not been initiated"<<endl;
    abort();
  }
  
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

  CheckVacancyMatrix(984);

  SuperInteractionCode(n,r,boundfix); //particle, direction, boundary.

  CheckVacancyMatrix(985);
  
  //    MoveAndBoundaryCheck(n,r,boundfix);
  //moves the particle n in direction r, acording to boundary conditions. 
  //if the new site is occupied, move back to previous position

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
  if (distribution == 3)    dist="nakazato";

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
  system("rm txyz.dat"); //vill inte ha denna filen pa min dator...


}


void Partiklar::StdErr(int Ensemble, int AntalPunkter){
  bool test = 0; //testvariabel

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
  bool fixBoundaryOn=FIXBOUNDARY;        //boundary condition, false = periodic
  bool testJaNej = false;
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
      cout <<"Logscale: "<<SecondArg<<endl;
    }
    else cout <<SecondArg<<
           " is wrong, must be: 0=logscale OFF, or 1=logscale ON"<<endl;
  }

  if (argc>=4){           //if the third arg=1 -> rescale=on
    int ThirdArg = atoi(argv[3]);
    if (ThirdArg==1 || ThirdArg==0){
      rescale = ThirdArg;
      cout <<"Rescale: "<<ThirdArg<<endl;
    }
    else cout <<ThirdArg<<
           " is wrong, must be: 0=normal Time, or 1=rescale ON"<<endl;
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

      int n = DISTRIBUTION;
      // 0 = uniform
      // 1 = exponential
      // 2 = power law
      // 3 = black & white
      //Characteristic Trait of the choosen distribution, 
      //just used to print info to file/screen
      double Info=0; 
      float gamma_black;
      
      //NOTE: The jumprate is in EACH direction! Meaning, jumprate =1 is actually
      //a jumprate =4 in 2D. (I think, test this NOW)
      if (0 <= n && n<=3 && (u < 1 && u > 0) ){
        gamma_black = JUMPRATE_TRACER;       //Jumprate for the tracer particle!
        float lambda=1.0;                    //used in option 1 Exp-dist.
        float y_c = 1.0;                     //used in option 2 Power-law
        float alpha=1.0;                     //used in option 2 Power-law
        float gamma_white = JUMPRATE_TRACER; //used in option 3 Nakazato
        crowd.set_gamma(gamma_white, gamma_black); 

        switch(n){
        case 0:           
          u = u;          
          Info = 1.0;      //Maximum value for uniform random number
          break;
        case 1: 
          u = (-log(u)) / lambda;   //funkar inte alls...eller?
          Info = lambda;
          break;
        case 2:
          //powerlaw, p(y) = alfa/y_c * (y/y_c)^{-1-alfa} if y > y_c, 0 otherwise
          u = y_c * pow(u,-1/alpha);  
          Info = alpha;
          break;
        case 3:
          if (particle==0) u = gamma_black; //set juprate for the first
          else u = gamma_white;             //set juprate for the rest
          Info = gamma_black/gamma_white;
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
  
  float start_time, now_time;  //just for fun, calculate time remaining on simulation...
  start_time = time(0); 


  for(int E=0; E < ensemble; E++){

    //Some simulation info:
    now_time=time(0);
    float time_remain =((now_time - start_time)/(E+1) * (ensemble - E))/60;
    cout <<"#"<<setw(5)<< 1.0*E/ensemble * 100<<" %  REMAINING ENSEMBLES:"
         <<setw(5)<<ensemble-E<<" Time (min): " <<(int) time_remain<<endl;


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
          //          t_output_next =  1.0* (n+1) * delta_t;
        }
        
        n++;
      } //Denna loop körs AntalPunkter gånger...

      (nymetod) ? crowd.Move2(fixBoundaryOn) : crowd.Move1(fixBoundaryOn);

    }
    //TODO: sista värdena skrivs aldrig ut, loopen slutar före utskriften...
    //...det är give or take, för vill att den börjar i nollan.
    crowd.Lagra(x_disp, y_disp, z_disp, r_disp); //Lagra hela ensemblens värden   0 
  }

  
  crowd.Save(ensemble, NameOfFiles);   // calculate std.dev. and save to file

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

  //TEST-Börja-----skriv ut hela summan:

  // for (int j=0; j<partial_sum.size(); j++)
  //   cout <<partial_sum[j]<<"\t";
  // cout <<endl;

  //TEST-Slut------

  //TEST
  if (partial_sum.back() < 0){
    cout<<"partial_sum too big"<<endl
        <<"...(that's what she said)"<<endl;
    abort();
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

/*
  plot "0-40000-4000.dat" using 1:2 with line title "1D",\
  0.282/sqrt(x) title "1/sqrt(t)",\
  "0-200-200-4000-(100).dat" using 1:2 with line title "2D",\
  0.19*log(x)/x title "ln(t)/t"


  plot f(x) title "dens=0.5, D=1/2", 'out.dat' title "100x1x1, N=50"
  quit


  plot f(x), '50x50-ergodicitet.dat' title "50x50x1, N=1"

  f(x)=((1-0.2)/0.2)*sqrt(4*0.258362*x/pi)*(0.2*0.2)
  plot '1dZeroAndOne.dat' with errorbars, f(x)
*/

