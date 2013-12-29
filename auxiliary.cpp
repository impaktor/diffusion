#include <iostream>        //for cout, among others
#include <cstdlib>         //for abort()-function
#include <string>
#include <fstream>
#include <vector>
#include <sstream>


#include "auxiliary.h"

//Read in specifications from command flags, return value through references.
//See: http://www.gnu.org/software/libc/manual/html_node/Getopt.html
//for information on the getopt()-function, which is what this entire
//function uses.
void argumentFlags(int argc, char** argv, InputValues& in){

  char* writeFileName = NULL;        //output file-name
  char* readFileName = NULL;         //output file-name

  bool changedOutputFileName = false; //if output file-name changed
  bool changedInputFileName = false;  //if input file-name changed

  bool error = false;

  //make sure we only set -b or -B, not both!
  int usingMultipleOutputs = 0;

  int c;
  opterr = 0;

  while((c = getopt(argc,argv,"lhqmw:r:i:jB:b:t:")) != -1)
    switch (c){
    case 'l':
      in.isLogScale = true;
      break;
    case 'm':
      in.isLowMem = true;
      break;
    case 'h':
      printHelp(in.outputFileName, argv);
      break;
    case 'q':
      in.isQuiet = true;
      break;
    case 'w':
      writeFileName = optarg;
      changedOutputFileName = true;
      break;
    case 'r':
      readFileName = optarg;
      changedInputFileName = true;
      break;
    case 'i':
      in.interactionStrength = atoi(optarg);
      in.isInteracting = true;
      if (in.interactionStrength < 0) error = true;
      break;
    case 'j':
      in.isJackknife = true;
      break;
    case 'b':
      in.method = 'b';
      in.nOutputs = atoi(optarg);
      ++usingMultipleOutputs;
      break;
    case 'B':
      in.method = 'B';
      in.nOutputs = atoi(optarg);
      ++usingMultipleOutputs;
      break;
    case 't':
      in.jmpTracer = atof(optarg);
      break;
    case '?':
      error = true;
      if (optopt =='w' || optopt =='i' || optopt =='B'
          || optopt =='b' || optopt =='t')
        //fprintf(stderr, "Option -%c requires an argument. \n",optopt);
        std::cout << "Option -" << (char) optopt << " requires an argument"
                  << std::endl;
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

  if (usingMultipleOutputs > 1){
    std::cout << "  Correct use is: " << argv[0] << " -[b/B] N" << std::endl
              << "  where N = number of simulation re-runs "<< std::endl
              << "  ie. number of output-files, and [b/B] is either:\n"
              << "    -b   bootstrap (fast)" << std::endl
              << "    -B   Brute force (roughly N times slower)" << std::endl;
    error = true;
  }

  //print error message if we entered case '?' above
  if(error){
    std::cout << argv[0] << ": undefined entry." << std::endl
         << "Try \"" << argv[0] << " -h\" for help." << std::endl;
    exit(1);
  }


  //if none of the above specified arguments were passed:
  //take the first as new writeFileName, and print the rest as
  //errors
  int index;
  bool moreThanOneInvalidArgument = (changedOutputFileName) ? true : false;

  for(index = optind; index < argc; index++){
    //only run first iteration & if we didn't change filename (with w|f)
    if(!moreThanOneInvalidArgument){
      writeFileName = argv[index];
      changedOutputFileName = true;
    }

    //run for all invalid arguments except the first one
    if(moreThanOneInvalidArgument)
      std::cout << "Non-option argument " << argv[index] << std::endl;

    moreThanOneInvalidArgument = true;
  }

  if(changedOutputFileName){
    in.outputFileName = writeFileName; //converts c-array to std::string
    std::cout << "# Default output file changed to: "
              << in.outputFileName << std::endl;
  }

  if(changedInputFileName){
    in.inputFileName = readFileName;   //converts c-array to std::string
    std::cout << "# using input file: "
              << in.inputFileName << std::endl;
  }

}



bool isComment(const std::string& line){

  //if line contains a '#' (anywhere)...
  unsigned int elementPos = line.find_first_of("#",0);

  if (elementPos != std::string::npos){

    //...check if it is preceded by ' ', and/or '\t', else not a comment
    std::string toPrecedComments = " \t";
    if (line.find_first_not_of(toPrecedComments) == elementPos)
      return true;
  }

  return false;
}


template<typename T>
T convertString(const std::string& string){
  //this is not optimal, since it does not return the stat of the
  //operation...

  std::istringstream stream(string);

  T returnValue;

  stream >> returnValue;

  return returnValue;

}


bool readInputFile(std::string& fileToRead, int& xdim, int& ydim,
                   int& zdim, int& particleNumber, int& ensembles,
                   int& numberOfValues, double& stopTime){

  std::vector<double> store;

  std::ifstream readFile;
  readFile.open(fileToRead.c_str());

  if (!readFile){ //failed to open file.
    std::cout << "# No input-file found." << std::endl;
    return false;
  }
  else{

    std::string readString;

    while(std::getline(readFile,readString)){
      if (!isComment(readString)){

        double temp = convertString<double>(readString);

        store.push_back(temp);
      }
    }

    xdim =            (int) store[0];
    ydim =            (int) store[1];
    zdim =            (int) store[2];
    particleNumber =  (int) store[3];
    ensembles =       (int) store[4];
    numberOfValues =  (int) store[5];
    stopTime =              store[6];

    std::cout << "# Using input from file: " << fileToRead << std::endl;

    return true;
  }
}


template<typename T>
T getNonCommentInput(){

  std::string s;

  //get non-comment line:
  //(useful when starting the program with: './prog < input.dat')
  do{
    getline(std::cin,s);
  }while(isComment(s));

  return convertString<T>(s);

}


//(return simulation parameters through reference)
void askUserForInputParameters(int& xdim, int& ydim, int& zdim,
                               int& particleNumber, int& ensembles,
                               int& numberOfValues, double& stopTime){

  std::cout << "# Specify number of lattice sites in:" << std::endl;
  std::cout << "# X: ";
  xdim = getNonCommentInput<int>(); //ignore lines starting with '#'

  std::cout << "# Y: ";
  ydim = getNonCommentInput<int>(); //ignore lines starting with '#'

  std::cout << "# Z: ";
  zdim = getNonCommentInput<int>(); //ignore lines starting with '#'

  std::cout << "# Total number of particles: ";
  particleNumber = getNonCommentInput<int>();

  float info = 1.0*(xdim*ydim*zdim)/particleNumber;

  std::cout << "# Number of vacant sites/particle: " << info << std::endl;
  std::cout << "# Number of ensembles: ";
  ensembles = getNonCommentInput<int>();
  std::cout << "# Number of data-values to save: ";
  numberOfValues = getNonCommentInput<int>();
  std::cout << "# Stop-time: ";
  stopTime = getNonCommentInput<double>();
  std::cout << std::endl << std::endl;
}



void printHelp(std::string File, char** argv){
  std::string defaultFileName = File;
  std::cout << "USAGE" << std::endl << "-----" << std::endl;
  std::cout << argv[0] <<" [OPTIONS]" << std::endl << std::endl;

  std::cout << "OPTIONS" << std::endl <<"-------" << std::endl;
  std::cout << "-i M\t Set interactions to ON with interaction strength \"M\".\n\t Note: M >= 0, since it will be used in the Boltzmann step as \"exp(-M)\"."<<std::endl;
  std::cout << "-l \t Set log-spaced data values to ON." << std::endl;
  std::cout << "-m\t Use the low memory consuming algorithm for computing std err." << std::endl;
  std::cout << "-w NAME\t Set output-file to 'NAME'." << std::endl;
  std::cout << "-r NAME\t Set input-file to 'NAME'." << std::endl;
  std::cout << "NAME\t Set first non-valid options-flag as name of output-file.\n \t(shorter/faster version than using the '-w' option)"<<std::endl;
  std::cout << "-q\t Set \"quiet mode\"=true. I.e print less information to the screen." << std::endl;
  std::cout << "-j\t Turn on the jackknife algorithm. The -m flag can not be used in this mode." << std::endl;
  std::cout << "-b N\t Generate the N number of outfiles (NAME0, NAME1,...) using bootstrap." << std::endl;
  std::cout << "-B N\t Generate the N number of outfiles (NAME0, NAME1,...) using Bruteforce." << std::endl;
  std::cout << "-t K\t Set tracer jumprate to K. Default = 1" << std::endl;

  std::cout << std::endl << "Arguments can be given in any order, and can be grouped (like \"-li 4\"). If no flags are set the first argument following \"" <<argv[0]<<"\" will be taken as the output-file name.\n Note, with '-i 0' the interaction strength will be 0, but the interaction code will be turned ON, resulting in a slower simulation run, therefore the default circumvents the interaction calculations entirely."<<std::endl<<std::endl;

  std::cout << "EXAMPLES" << std::endl << "--------" << std::endl;

  std::cout << "\t" << argv[0] << " 5.5/uniform04.dat" << std::endl
            << "\t" << argv[0] <<" -w 5.5/uniform04.dat" << std::endl
            << "Puts the output in the folder \"5.5\" with the name \"uniform04.dat\"."
            << std::endl;

  std::cout << "\t" << argv[0] << " -B 1000 -w brute.dat" << std::endl
            << "Use brute force, save to brute.dat0 ... brute.dat999.\n";
  std::cout << std::endl;

  exit(0);
}


void printError(std::string message){
  std::cerr <<"\n Error: " << message << std::endl;
  abort();
}

//use __FILE__ and __LINE__ macro to get line and filename;
void printError(std::string message, int line){
  std::cerr <<"\n Error occured at line: " << line
       << std::endl << message << std::endl;
  abort();
}

void printError(std::string message, std::string file, int line){
  std::cerr <<"\n Error occured in file " << file << " at line: " << line
       << std::endl << message << std::endl;
  abort();
}

