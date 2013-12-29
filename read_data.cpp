#include <iostream>    //for std::cout, among others
#include <fstream>
#include <cstdlib>     //for abort()-function
#include <cmath>       //gives us sqrt, pow, fabs,
#include <string>
#include <vector>
#include <sstream>     //for [o]stringstream
#include <limits>
#include <stdlib.h>    //atoi
#include <cstdio>      //sscanf
#include <algorithm>   //std::min

#include "auxiliary.h" //non-physics stuff. (print messages etc.)
#include "save.h"      //class to save MSD, compute errors, and print to file

const bool debug = false;

//TODO: Hur är formatet på hans filer? samma som på square? är det dr
//eller dr2? (skall vara dr). Har han "#SEED = 103902" mellan varje
//körning?


template<typename T>
T parse(const std::string &s_in, std::string str){
  //Go through str until we find the matching str, then remove " = "
  //and return next thing as val;
  T val;

  std::stringstream ss;
  ss.str(s_in);

  std::string s;
  do{
    ss >> s;
    //std::cout << "throw away:\t" << s << std::endl;
  }while(!(s == str || ss.bad() || ss.eof()));

  //std::cerr << "Error when parsing header (\"ensembles\")" << std::endl;

  int pos = ss.tellg();    //get current position of point in stream
  ss >> s;                 //remove " = "
  if(s != "="){            //maybe there's no "=" (as in the case of "distr : nakazato (n)"
    std::stringstream tmp(s);
    ss << " " << tmp;     //put s back. (ss.peek() only works on single char forward)
    ss.seekg(pos);        //put point back.
  }

  ss >> val;              // get sought value
  return val;
}


void getSimulationInfo(const std::string &head, int &M, double &r_t,
                       double &r_c, double &conc){
  //get the simulation info from the comments.

  M    = parse<int>(head, "#E");
  conc = parse<double>(head, "#Concentration");
  r_t  = parse<double>(head, "k_tagg");

  std::string foo = parse<std::string>(head, "nakazato");
  foo = foo.substr(1, foo.size()-2);          //remove parenthesis
  std::stringstream ss(foo);
  ss >> r_c;                                  //convert to number
  r_c = 1.0 / r_c;                            //inverse back
}


std::string readHeader(std::istream &input){
  //Remove first four lines of comments marking the start of each new
  //trajectory from input stream.

  std::string tmp_string, str;

  //Skipp first four lines
  for(int i = 1; i <= 4; ++i){
    std::getline(input, tmp_string);
    str += tmp_string + '\n';
  }

  //ignore stream until next end of line
  if((input >> tmp_string) && tmp_string == "#SEED"){
    input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  else{
    std::cerr << "archangel!" << std::endl;
    abort();
  }

  return str;   //successfully removed header (comments)
}


bool readOneTrajectory(std::istream &stream, Save::vectorD_t &dt, Save::vectorI_t &dx,
                       Save::vectorI_t &dy, Save::vectorI_t &dz, Save::vectorD_t &dr,
                       double minTime){
  //Start reading from stream of five columns (double, int, int, int,
  //double) and stop at first instance of a '#' and return true.

  std::string s;

  while(std::getline(stream, s)){
    std::stringstream ss(s);
    int number;
    double numberDouble;

    //if reached the next block of data
    if(ss.peek() == '#')
      return true;  //done with this trajectory, ready for next
    else{
      ss >> numberDouble;
      if(debug) std::cout << "t " <<  numberDouble << std::endl;
      if(numberDouble < minTime)
        continue;
      dt.push_back(numberDouble);

      ss >> number;
      if(debug) std::cout << "x " <<  number << std::endl;
      dx.push_back(number);

      ss >> number;
      if(debug) std::cout << "y " <<  number << std::endl;
      dy.push_back(number);

      ss >> number;
      if(debug) std::cout << "z " <<  number << std::endl;
      dz.push_back(number);

      ss >> numberDouble;
      if(debug) std::cout << "r " <<  numberDouble << std::endl;
      dr.push_back(numberDouble);
    }
  }
  return false;     //reached end of stream?
}


void readFile(std::string path, Save &save, double minTime){
  //Read in all the trajectories of displacements for all ensembles

  std::ifstream input(path.c_str());
  std::string s;

  if(!input)
    throw std::string("Could not load file: " + path);

  if(input.is_open()){

    bool read = true;
    int i = 1;
    readHeader(input);   //skipp header
    while(read){
      if(debug) std::cout << "Reading trajectory: " << i << std::endl;

      Save::vectorI_t dx, dy, dz;
      Save::vectorD_t dr, dt;
      read = readOneTrajectory(input, dt, dx, dy, dz, dr, minTime);
      save.store(dx, dy, dz, dr);

      // // when testing, dont take for ever: XXX
      // i++;
      // if(i > 1000)
      //   break;
    }
  }
  else std::cerr << "No file: " << path << std::endl;
  input.close();
}



int findActualNumTrajectories(std::string path){
  //Since the "Ensembel = " in the header is wrong sometimes (Pigon
  //uses cat on several outputs)

  std::ifstream input(path.c_str());
  std::string s;

  size_t M = 0;  // number of trajectories

  if(!input)
    throw std::string("Could not load file: " + path);

  while(std::getline(ínput, s)){
    std::stringstream ss(s);

    if(ss.peek() == '#')
      M++

  }


  if(input.is_open()){

    bool read = true;
    int i = 1;
    readHeader(input);   //skipp header
    while(read){
      Save::vectorI_t dx, dy, dz;
      Save::vectorD_t dr, dt;
      read = readOneTrajectory(input, dt, dx, dy, dz, dr, minTime);
      save.store(dx, dy, dz, dr);

      // // when testing, dont take for ever: XXX
      // i++;
      // if(i > 1000)
      //   break;
    }
  }
  else std::cerr << "No file: " << path << std::endl;
  input.close();

}



void readTime(std::string path, Save::vectorD_t &time, double minTime){
  //get the time vector from the first trajectory of file at "path".

  std::ifstream input(path.c_str());
  std::string s;

  if(input.is_open()){
    readHeader(input);            //throw away

    Save::vectorI_t dx, dy, dz;   //throw away
    Save::vectorD_t dr;
    readOneTrajectory(input, time, dx, dy, dz, dr, minTime);

    std::cout << "getting sampling times. " << std::endl;
  }
  else std::cerr << "No file: " << path << std::endl;

  input.close();
}



int main(int argc, char* argv[]){

  if(argc == 3){

    bool isLowMem = false;

    std::string header;
    std::ifstream input(argv[1]);
    if(input.is_open())
      header = readHeader(input);
    else
      std::cerr << "No file: " << argv[1] << std::endl;
    input.close();

    int M;
    double r_t, r_c, c;
    getSimulationInfo(header, M, r_t, r_c, c);

    // TODO/XXX Extremly ugly fix: Pigeon runs 8000*3 for 3D lattices,
    // and uses cat:
    if(M == 8000) M = 3*8000;

    std::cout << "\nM = " << M << "\tr_t = " << r_t << "\tr_c = " << r_c
              << "\tc = " << c << "\n" << std::endl;

    //compute from which sampling time to start
    double minTime = 10.0 / (std::min(r_c,r_t)*c);

    Save::vectorD_t time;
    readTime(argv[1], time, minTime);

    //M = 1000;  //TMP TODO XXX using this when debugging

    //initiate the save-class to store, save, and print result
    //(also computes standard error and correlation estimate)
    Save save(time, M, isLowMem);
    readFile(argv[1], save, minTime);

    std::string outputName(argv[2]);
    save.save(outputName, header);

    int jackknifeGroups = 100;
    save.printSlopeJackknifeShrinkage(outputName, jackknifeGroups);

    return 0;
  }
  else
    std::cerr << "Specify file to read, and name of output:\n"
              << argv[0] << " [input] [output]" << std::endl;
}
