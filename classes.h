#ifndef CLASSES_H
#define CLASSES_H

//#include <iostream>        //for cout, among others
#include<time.h>

//First three classes (Direction, Jump, Particle), will allow us to 
//get the position of a Particle object as: object.x, and the jump 
//rate as object.rate.x.r etc.

class Direction{
public:
  Direction(double right,double left);
  Direction(void){};
  double r;  //right
  double l;  //left
};

class Jump{
public:
  Direction x; //i.e. x.r and x.l
  Direction y;
  Direction z;
};

class Particle{
public:
  Particle(int,int,int);
  Particle(void){};
  int x;            //position x
  int y;
  int z;
  
  Jump rate;        //jump rate: rate.x.r

  //This allows us to print x,y,z with the << operator
  //like: Particle particle(1,4,2);  cout << particle <<endl; 
  friend std::ostream& operator<< (std::ostream &out, Particle &particle); 
  friend std::istream& operator>> (std::istream &in, Particle &particle); 
};





//Class to show runtime information, such as remainging time
//until completion, and number of ensembles left. 
class RemainingTime{
private:
  time_t startTime_;
  unsigned int totalEnsembles_;

public:
  RemainingTime(unsigned int);
  void printProgress(unsigned int);  //print progress to screen
};


#endif




/*
#include <iostream.h>
#include <string.h>
class Base{
  char strVal[100];
public:
  Base(){
    strcpy(strVal,"");
  }
  Base(char *val){ 
    strcpy(strVal,val);
  }
  ~Base(){
    *strVal = '\0';
  }
  friend istream& operator >>(istream &is,Base &obj);
  friend ostream& operator <<(ostream &os,const Base &obj);
};

istream& operator >>(istream &is,Base &obj)
{
  is>>strVal;
  return is;
}

ostream& operator <<(ostream &os,const Base &obj)
{
  os<<obj.strVal;
  return os;
}

void main()
{
  Base b;
  cin>>b;
  cout<<"Printing the value\n";
  cout<<b<<endl;
}
*/

