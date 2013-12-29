#include <iostream>        //for cout, among others
#include <fstream>         //read write to files
#include <cstdlib>         //for abort()-function
#include <ctime>           //to use system time
#include <vector>          //for the vector class
#include <cmath>           //gives us sqrt, pow, fabs, 
#include <algorithm>       //needed for min_element in hist/bin-function
#include "nr3.h"           //Random number generator from numerical recipes 3 ed.
#include "ran.h"           

//Main-function usage: ./program [outputN.dat] [0/1]
//where the two arguments are optional. The default is to produce an output to "out.dat", with values with linear spacing (logscale off =0). If you want logscale=ON, but still use default "out.dat", you must still specify it: ./program out.dat 1. 

// use ssh:
// $ ssh -n 'cd /to/path; nice -19 ./program <input.dat > output1.dat'
// Or in this more recent version:
// $ ssh -X nice -19 ./program output.dat < input.dat'

#define DISTRIBUTION     2    //0=uniform, 1=exponential, 2=power law, 3=nakazato
#define JUMPRATE_TRACER  0.25
#define JUMPRATE_CROWDER 1.0  //(if nakazato-distribution)
#define TEST_1           0    //boolean, prints details to screen

using namespace std;

class Partiklar {
  vector<int> Xpos,Ypos,Zpos;              //XYZ-position for each particle
  vector<double> t;
  vector<long double> x2_err,y2_err,z2_err,r2_err; //Store the final error, in each point
  vector<double> x2_mu,y2_mu,z2_mu,r2_mu;  //average of displacement squared: <X^2>, <Y^2>... 
  vector<double> x_mu, y_mu, z_mu, r_mu;   //average of displacement: <X>, <Y>...

  vector<double> jmpRateXright,jmpRateXleft,jmpRateYright,
    jmpRateYleft,jmpRateZup,jmpRateZdown;  //jump-rate for each particle and direction
  vector <long double> partial_sum;        //cumulative sum of jump-rates

  int RutorX,RutorY,RutorZ;                //Lattice size
  int N;                                   //Total number of particles
  int dim;                                 //dimension of lattice (=1 || 2 || 3)
  float density;                           //particle density
  int X_0, Y_0, Z_0;                       //Initial position for the fist particle
  double time_sum;                         //Sum of time for each move
  int ConvertMuToParticle(int mu, int& n);
  void PartialSum(void);
  void MoveAndBoundaryCheck(int,int,bool);
  int VacancyCheck(int, int,int,int);      //No double Occupancy.
  int VacancyCheck2(int, int,int,int);     //Same as above, but improved. Use this one instead.
  void StdErr(int, int);                   //calculate standard deviation

  vector< vector<int> > dx, dy, dz;        //pos. of tracer particle. needed for stdrerr
  vector< vector<double> > dr;             //same as line above

  vector< vector< vector<int> > > vacancy; //Needed in my more efficient "vacancyCheck2()"

  int seed1, seed2;

  float jump_crowders, jump_tracer;        //jumprate for the (Nakazato) particles
  int distribution;
  double Info;                             //save caracheristic trait of particle distribution
  int round(float x) {                     //Round off
    return int(x + 0.5);  }
  bool scaling;
  bool testOnOff;                          //To print detailed information to screen
  double D_av; 

  float Nakazato(bool);         
  double DiffEffConstant(void);            //Calculate the effective dif. const. from the jump-rates
  float Ergodicity(int,int,int, bool);     //Calculate the theoretical MSD when ergodic

public:
  Partiklar(int, int, int, int);
  void Place(void);                        //Place particles on the lattice
  void Move1(bool);                        //old obsolete move-code (no jump-rates)
  void Move2(bool);                        //pick a particle & direction based on jump-rate
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

  //Set jump-rate of each particle and direction, depending on if dim=1,2,3. 
  void setJumpRate(vector<double>&, vector<double>&, vector<double>&,
                   vector<double>&, vector<double>&, vector<double>&);
  void setJumpRate(vector<double>&, vector<double>&, vector<double>&, vector<double>&);
  void setJumpRate(vector<double>&, vector<double>&);

  //Secondary not used functions
  void Snapshot(void);                     //Not used. (more comments in function body)
  //  double get_dt(void);   //Size of time-step, need it for P(x,t), don't need this...

};


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
                     int ParticleNumber){

  if (Xsquare >= Ysquare && Ysquare  >= Zsquare && Zsquare >= 1){
    if (ParticleNumber <= Xsquare*Ysquare*Zsquare){
      seed1=1;
      seed2=2;                   //kanske ändra hur dessa får sina värden!

      RutorX=Xsquare;            //number of lattice sites
      RutorY=Ysquare;
      RutorZ=Zsquare;
      N=ParticleNumber;

      X_0=round(1.0*RutorX/2); 
      Y_0=round(1.0*RutorY/2); //center of lattice, where we place our tagged particle
      Z_0=round(1.0*RutorZ/2); 
      
      density=(float) N / (RutorX*RutorY*RutorZ); 

      //Hardcoded Switches:
      testOnOff=TEST_1;           //boolean, for print details to screen. 

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
  static Ran randomNumber1(seed1);

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


  //IF We use VacancyCheck2()-function:
  //-----------------------------------
  //initiate the vacancy-vector with zeros. 
  //Must reset this matrix for each ensamble!
  vector< vector< vector<int> > >
    Vacancy(RutorX+1, vector< vector<int> >(RutorY+1, vector<int>(RutorZ+1,-1)));
  vacancy=Vacancy;
  //Mark occupied sites as n:
  //note, no change will ever occur to the 0-elements.
  for (n=0; n<N; n++){
    vacancy[ Xpos[n] ][ Ypos[n] ][ Zpos[n] ]=n;
  }

}



void Partiklar::Move1(bool val){ //val=1 -> fix, val=0 -> periodic
  //This code is obsolete, as it just works for placing identical
  //particles, but I keep it for sentimental reasons.
  //Better safe than sorry.
  
    static Ran randomNumber2(seed2);
    int r,n; 
  
    do{
      n=randomNumber2.doub() * N; 
    }while(n==N);                     //Choose particle at random, 0<= n < N.
      

    do{
      r=randomNumber2.doub() * (2 * dim); 
    }while( r== dim * 2 );            //gives 0<=R<2*dim
      

    //Do the actual move, acording to the boundary rules
    MoveAndBoundaryCheck(n,r, val);

    time_sum = time_sum + 1.0/N;    //Calculate the time

    //TEST:
    if (testOnOff)
      for (int q=0; q < N; q++)
        cout <<" MOVED \n particle n="
             <<q<<" X="<<Xpos[q]<<" Y="<<Ypos[q]<<" Z="<<Zpos[q]<<endl;
}


void Partiklar::MoveAndBoundaryCheck(int n, int R, bool boundary){
  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant. 
  
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
      
      //VacancyCheck(n,Xtemp,Ytemp,Ztemp);  //Old code (slow)                                        
      VacancyCheck2(n,Xtemp,Ytemp,Ztemp);   //New improved, faster version.

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
  //Old and obsolete version. It doesn't need a 3D matrix as the
  //new version does, VacancyCheck2(...) but this one is slow.
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

      return returnvalue;
    }
    else{                                       //if vacant, update vacancy-matrix 
      vacancy[Xold][Yold][Zold]=-1;             //old position is now vacant
      vacancy[Xpos[n]][Ypos[n]][Zpos[n]]=n;     //new position is now occupied
      if (testOnOff) cout <<"OK, vacant"<<endl;
    }
    //if (testOnOff) cout <<"Vacant, particle "<<n<<" X = "<<Xpos[n]<<", Y = "<<Ypos[n]<<", Z = " <<Zpos[n]<<endl;
    return returnvalue;
    //(the returnvalue could be used to have a do-move-while(VacancyCheck==1) but
    //if so, we need to move this function to be directly in the Move_()-func.)
  }
 else{
   cout <<"Error, accessing invalid particle "<<endl;
   abort();
 }
}



void Partiklar::Move2(bool boundfix){
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
    //turns since double stores 52 bit, 2^52 = 10^15. 

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
  //Now process the simulation data, and print it to
  //file. 
  
  char* Primary=name;

  int AntalPunkter=t.size();
  ofstream primary,secondary;
  primary.open(Primary);

  //hardcoded switch to print x,y,z, x^2,y^2,z^2
  //specifically. Not needed now, since isotropic,
  bool xyz=false;
  if(xyz) secondary.open("txyz.dat");

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
  string dist[]={"uniform","exponential","powerlaw","nakazato"};
  float d_eff= DiffEffConstant();
  float d_av = D_av;


  //TEST -felsök gränserna i exp-funktionen!
  cout <<"# Lower limit: "<<(RutorX * RutorX) / (2*N*N*d_eff) <<endl;
  cout <<"# Upper limit: "<<(RutorX * RutorX) / (2*d_eff) <<endl;

 
  //print info about simulation to head of file:
  primary <<"#E = "<<Ensemble<<"\t N = "<<N<<"\t X-Y-Z: "
          <<RutorX<<"x"<<RutorY<<"x"<<RutorZ
          <<"\t 2*d*D_nakazato: "<<Nakazato(scaling)<<endl;
  primary <<"#Density: "<<density<<"\t MSD when erogodic: "
          <<Ergodicity(RutorX,RutorY,RutorZ,scaling)
          <<"\t distr: "<<dist[distribution]<<" ("<<Info<<")"<<endl;
  primary <<"#D_eff: "<<d_eff<<"\t k_tracer: "
          <<jmpRateXright[0]<<"\t D_av: "<<d_av<<endl; 

  //print info about simulation to screen:
  cout <<"#E = "<<Ensemble<<"\t N = "<<N<<"\t X-Y-Z: "
       <<RutorX<<"x"<<RutorY<<"x"<<RutorZ
       <<"\t 2*d*D_nakazato: "<<Nakazato(scaling)<<endl;
  cout <<"#Density: "<<density<<"\t MSD when erogodic: "
       <<Ergodicity(RutorX,RutorY,RutorZ,scaling)
       <<"\t distr: "<<dist[distribution]<<" ("<<Info<<")"<<endl;
  cout <<"#D_eff: "<<d_eff<<"\t k_tracer: "
       <<jmpRateXright[0]<<"\t D_av: "<<d_av<<endl; 
  
  for (int i=0; i<t.size(); i++){

    //rescale the axis of the primary output (dimensionsless)
    if (scaling){
      t[i] = t[i] * pow(density,2) * DiffEffConstant();            
      r2_mu[i] = r2_mu[i] * pow(density,2);
      r2_err[i]= r2_err[i]* pow(density,2);
    }

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


void Partiklar::StdErr(int Ensemble, int AntalPunkter){

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

  
  //Testing the Gaussian distribution of the points:
  //--------------------------------------------------------------------
  float SUM68=0;
  for (int j=0; j<AntalPunkter; j++){            //forgot about his: 
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
  cout <<"68,2 % (?)="<<1.0*SUM68/(Ensemble*AntalPunkter)<<endl;

  //--------------------------------------------------------------------
}


//////////////////////////////////////
//   ****************************   //
//   * M A I N  F U N K T I O N *   //
//   ****************************   //
//////////////////////////////////////

int main(int argc, char* argv[]){
  int N;                              //Number of particles
  int X,Y,Z;                          //lattice size
  int ensemble,AntalPunkter;
  double MaxTime;

  //DEFAULT VALUES:  
  bool fixBoundaryOn=true;            //boundary condition, false = periodic
  bool logscale=false;
  bool nymetod=true;                  //use Move2() instead of Move()
  bool rescale = false;               //do this in gnuplot/matlab instead
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

  if (fixBoundaryOn) 
    cout <<"# Fix boundary"<<endl;  
  else 
    cout <<"# Periodic boundary"<<endl;

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
  cout << endl; 

  Partiklar crowd(X,Y,Z,N);                   //X,Y,Z-squares, N-particles
  crowd.set_scaling(rescale);                 //if we rescale the axis
  
  static Ran tempnumb(17);                    //just anny seed will do
  vector<double> jXr, jYr, jZu, jXl, jYl, jZd; 
  if (nymetod){
    for (int particle = 0; particle < N; particle ++){
      double u;
      do { 
        u=tempnumb.doub();   
      }while (u==0);


      int n = DISTRIBUTION;  // <---< VAL AV METOD <----<
      // 0 = uniform
      // 1 = exponential
      // 2 = power law
      // 3 = nakazato

      //Characteristic trait of the choosen distribution, 
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

      //sätter jumprate för första till ;
      if (particle == 0 ) u = jump_tracer;

      
      jXr.push_back(u);
      jXl.push_back(u);
      jYr.push_back(u);
      jYl.push_back(u);
      jZu.push_back(u);
      jZd.push_back(u); 

      //TEST:
      //    cout <<"u: "<< u << endl;

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
  double delta_t=  1.0*MaxTime/AntalPunkter;
  double delta_log_t= 1.0*log(MaxTime)/AntalPunkter;    
  
  float start_time, now_time;  //just for fun, calculate time remaining on simulation...
  start_time = time(0);


  for(int E=0; E < ensemble; E++){

    //Some simulation info (like time remaining). (fabs=absolute value)
    now_time=time(0);
    float time_remain =(fabs(now_time - start_time)/(E+1) * (ensemble - E))/60;
    cout <<"#"<<setw(5)<< 1.0*E/ensemble * 100<<" %  REMAINING ENSEMBLES:"
         <<setw(6)<<ensemble-E<<" Time (min): " <<(int) time_remain<<endl;

    vector<int> x_disp, y_disp, z_disp;
    vector<double> r_disp;
    crowd.Place();

    double t_output=0;
    double punkt=0;
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
          //TODO, ta bort första punkten 0 0 0 , funkar inte med LOG!
          // t_output = exp(punkt);
          // punkt = punkt + delta_log_t;
          t_output = exp(delta_log_t *1.0*n);
        }
        else{  //t_output = t_output + delta_t;
          t_output = 1.0* n * delta_t;
        }
        
        n++;
      } //run this loop "AntalPunkter"-times

      //(nymetod) ? crowd.Move2(fixBoundaryOn) : crowd.Move1(fixBoundaryOn);
      crowd.Move2(fixBoundaryOn); //faster I would guess... (Move1() is obsolete)

      // //TEST: to get a snapshot close to the MaxTime
      // if(totalTime>=MaxTime-0.01){
      //   crowd.Snapshot();
      //   cout<<"totalTime:"<<totalTime<<endl;
      // }
    }
    //TODO: sista värdena skrivs aldrig ut, loopen slutar före utskriften...
    //...det är give or take, för vill att den börjar i nollan.
    

    //Store tracer position at each point for each ensemble,
    //inorder to get errorbars (and to bin and get a histogram of P(x,t))
    crowd.Lagra(x_disp, y_disp, z_disp, r_disp); 
    //NOTE: this method uses push.back()! --> slow and heavy on
    //memmory if large run. (many ensembles & points) 
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
       <<RutorX<<"x"<<RutorY<<"x"<<RutorZ<<"\t concentration="<<density<<endl;
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

