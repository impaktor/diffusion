#include <iostream>        //for cout, among others
#include <cstdlib>         //for abort()-function
#include <string>
#include <fstream>
#include <sstream>
#include <getopt.h>

#include "auxiliary.h"

//Read in specifications from command flags, return value through references.
//See: http://www.gnu.org/software/libc/manual/html_node/Getopt.html
//for information on the getopt()-function, which is what this entire
//function uses.
void aux::argumentFlags(int argc, char** argv, InputValues &in){

  bool error = false;

  int c;
  opterr = 0;

  while((c = getopt(argc,argv,"hmdw:r:i:jBb:k:t:")) != -1)
    switch (c){
    case 'm':
      in.isLowMem = std::make_pair(true, true);
      break;
    case 'h':
      printHelp(argv);
      break;
    case 'w':
      in.outputFileName = std::string(optarg);        //converts c-array (char*) to std::string
      std::cout << "# Default output file changed to: "
                << in.outputFileName << std::endl;
      break;
    case 'r':
      in.inputFileName = std::string(optarg);         //converts c-array (char*) to std::string
      std::cout << "# Default input file changed to: "
                << in.inputFileName << std::endl;
      break;
    case 'i':
      in.isInteracting = true;
      in.interactionStrength = atoi(optarg);
      if (in.interactionStrength < 0) error = true;
      break;
    case 'j':
      in.isJackknife = true;
      break;
    case 'k':
      in.isBootknife = true;
      in.nBootknife = atoi(optarg);
      break;
    case 'b':
      in.method = 'b';
      in.nOutputs = atoi(optarg);
      break;
    case 'B':
      in.isBruteForce = true;
      break;
    case 'd':
      in.method = 'd';
      break;
    case 't':
      in.jmpTracer = std::make_pair(atof(optarg), true);
      break;
    case '?':
      error = true;
      if(optopt =='w' || optopt =='i' || optopt =='b' || optopt =='t' || optopt == 'k')
        //fprintf(stderr, "Option -%c requires an argument. \n",optopt);
        std::cout << "Option -" << (char) optopt << " requires an argument" << std::endl;
      else if(isprint(optopt))
        //fprintf(stderr,"Unknown option '-%c'. \n",optopt);
        std::cout << "unknown option: -"<< (char) optopt << std::endl;
      else
        //fprintf(stderr,"Unknown option character '\\x%x'.",optopt);
        std::cout << "unknown option character " << optopt << std::endl;
      //return 1;
      break;
    default:
      std::cout << "Non-critical error in pogram (bad programming)"<< std::endl;
      abort();
    }

  //print error message if we entered case '?' above
  if(error){
    std::cout << argv[0] << ": undefined entry." << std::endl
         << "Try \"" << argv[0] << " -h\" for help." << std::endl;
    exit(1);
  }
}



bool aux::convertToBool(const char * c, std::string where){
  //Convert c-string starting with either"true ..." or "false ..." to bool.

  bool val = 0;

  std::string str, str0(c);
  std::stringstream ss(str0);
  ss >> str;

  if(str == "true")
    val = true;
  else if(str == "false")
    val = false;
  else
    printError("Could not parse:" + str + " in: " + where);

  return val;
}



void aux::printHelp(char** argv){
  std::cout << "USAGE" << std::endl << "-----" << std::endl;
  std::cout << argv[0] <<" [OPTIONS]" << std::endl << std::endl;

  std::cout << "OPTIONS" << std::endl <<"-------" << std::endl;
  std::cout << "-i M\t Set interactions to ON with interaction strength \"M\".\n\t Note: M >= 0, since it will be used in the Boltzmann step as \"exp(-M)\"."<<std::endl;
  std::cout << "-m\t Use the low memory consuming algorithm for computing std err." << std::endl;
  std::cout << "-w NAME\t Set output-file to 'NAME'." << std::endl;
  std::cout << "-r NAME\t Set input-file to 'NAME'." << std::endl;
  std::cout << "-j\t Turn on the jackknife algorithm. The -m flag can not be used in this mode." << std::endl;
  std::cout << "-k N\t Turn on the bootknife algorithm with N reshuffled MSD's." << std::endl;
  std::cout << "\t Can not be used with the -m flag." << std::endl;
  std::cout << "-b N\t Generate the N number of outfiles (NAME0, NAME1,...) using bootstrap." << std::endl;
  std::cout << "-B\t Bruteforce: each point uncorrelated" << std::endl;
  std::cout << "-t K\t Set tracer jumprate to K." << std::endl;
  std::cout << "-d\t Dump all raw trajectories to path specified by -w. No post processing." << std::endl;

  std::cout << std::endl << "Arguments can be given in any order, and can be grouped (like \"-li 4\"). If no flags are set the first argument following \"" <<argv[0]<<"\" will be taken as the output-file name.\n Note, with '-i 0' the interaction strength will be 0, but the interaction code will be turned ON, resulting in a slower simulation run, therefore the default circumvents the interaction calculations entirely."<<std::endl<<std::endl;

  std::cout << "EXAMPLES" << std::endl << "--------" << std::endl;

  std::cout << "\t" << argv[0] << " 5.5/uniform04.dat" << std::endl
            << "\t" << argv[0] <<" -w 5.5/uniform04.dat" << std::endl
            << "Puts the output in the folder \"5.5\" with the name \"uniform04.dat\"."
            << std::endl;
  exit(0);
}





std::string ERROR = "\033[1;31mERROR \033[0m";

void aux::printError(std::string message){
  std::cerr <<"\n" + ERROR + message << std::endl;
  abort();
}


//use __FILE__ and __LINE__ macro to get line and filename;
void aux::printError(std::string message, int line){
  std::cerr <<"\n Error occured at line: " << line
            << std::endl << message << std::endl;
  abort();
}

void aux::printError(std::string message, std::string file, int line){
  std::cerr <<"\n Error occured in file " << file << " at line: " << line
            << std::endl << message << std::endl;
  abort();
}
