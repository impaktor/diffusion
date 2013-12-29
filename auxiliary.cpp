#include <iostream>        //for cout, among others
#include <cstdlib>         //for abort()-function
#include <cstring>         //strcpy, strchr
#include <string>

#include "auxiliary.h"

//Read in specifications from command flags, return value through references.
//See: http://www.gnu.org/software/libc/manual/html_node/Getopt.html
//for information on the getopt()-function, which is what this entire
//function uses.
void argumentFlags(int argc, char** argv, bool& logarithm, bool& lowMem,
                   bool& interactionON, bool& quiet, int& interactionStrength,
                   std::string& returnFileName, int& noRuns, char& nMethod,
                   bool& jackknife, float& tracer){

  char* fileName = NULL;
  bool changedFileName = false;
  bool error = false;

  //make sure we only set -b or -B, not both!
  int usingMultipleOutputs = 0;

  int c;
  opterr = 0;

  while((c = getopt(argc,argv,"lhqmw:i:jB:b:t:")) != -1)
    switch (c){
    case 'l':
      logarithm = true;
      break;
    case 'm':
      lowMem = true;
      break;
    case 'h':
      printHelp(returnFileName,argv);
      break;
    case 'q':
      quiet = true;
      break;
    case 'w':
      fileName = optarg;
      changedFileName = true;
      break;
    case 'i':
      interactionStrength = atoi(optarg);
      interactionON = true;
      if (interactionStrength < 0) error = true;
      break;
    case 'j':
      jackknife = true;
      break;
    case 'b':
      nMethod = 'b';
      noRuns = atoi(optarg);
      ++usingMultipleOutputs;
      break;
    case 'B':
      nMethod = 'B';
      noRuns = atoi(optarg);
      ++usingMultipleOutputs;
      break;
    case 't':
      tracer = atof(optarg);
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
  //take the first as new fileName, and print the rest as
  //errors
  int index;
  bool more_than_one_invalid_argument = (changedFileName) ? true : false;

  for(index = optind; index < argc; index++){
    //only run first iteration & if we didn't change filename (with w|f)
    if(!more_than_one_invalid_argument){
      fileName = argv[index];
      changedFileName = true;
    }

    //run for all invalid arguments except the first one
    if(more_than_one_invalid_argument)
      std::cout << "Non-option argument " << argv[index] << std::endl;

    more_than_one_invalid_argument = true;
  }

  if(changedFileName){
    returnFileName = fileName;   //converts c-array to std::string
    std::cout << "# Default output file changed to: "
              << returnFileName << std::endl;
  }

}


template <class T>
T getNonCommentInput(void){
  //small function to allow comments in the input-file.
  //this will ignore any line containing "#", since it
  //only wants numbers.

  char invalue[999];

  for(bool repeat = true; repeat;){
    std::cin >> invalue;
    if(strchr(invalue,'#') == NULL)
      repeat = false;
  }

  return (T) atof(invalue);
}

//(return simulation parameters through reference)
void AskUserForInputParameters(int& xdim, int& ydim, int& zdim,
                               int& particleNumber, int& ensembles,
                               int& numberOfValues, double& stopTime){

  std::cout << "# Specify number of lattice sites in:" << std::endl;
  std::cout << "# X: ";
  xdim = getNonCommentInput<int>(); //ignore lines containing '#'
  std::cout << "# Y: ";
  ydim = getNonCommentInput<int>(); //ignore lines containing '#'
  std::cout << "# Z: ";
  zdim = getNonCommentInput<int>(); //ignore lines containing '#'

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
  std::cout << "\t" << argv[0] <<" [OPTIONS]" << std::endl << std::endl;

  std::cout << "OPTIONS" << std::endl <<"-------" << std::endl;
  std::cout << "-i M\t Set interactions to ON with interaction strength \"M\".\n\t Note: M >= 0, since it will be used in the Boltzmann step as \"exp(-M)\"."<<std::endl;
  std::cout << "-l \t Set log-spacing between data values to be written to\n\t output-file to ON." << std::endl;
  std::cout << "-m\t Use the low memory consuming algorithm for computing std err." << std::endl;
  std::cout << "-w NAME\t Set output-file to 'NAME'." << std::endl;
  std::cout << "NAME\t Set first non-valid options-flag as name of output-file.\n \t(shorter/faster version than using the '-w' option)"<<std::endl;
  std::cout << "-q\t Set \"quiet mode\"=true. I.e print less information to the screen." << std::endl;
  std::cout << "-j\t Turn on the jackknife algorithm. The -m flag can not be used in this mode." << std::endl;
  std::cout << "-b\t Generate the N number of outfiles (NAME0, NAME1,...) using bootstrap." << std::endl;
  std::cout << "-B\t Generate the N number of outfiles (NAME0, NAME1,...) using Bruteforce." << std::endl;

  std::cout << std::endl << "Arguments can be given in any order. If no flags are set the first argument following \""<<argv[0]<<"\" will be taken as the output-file name. With no arguments given the default is log-spacing of values to output-file \"" << defaultFileName <<"\" without any interactions at all. Note, with '-i 0' the interaction strength will be 0, but the interaction code will be turned ON, resulting in a slower simulation run, therefore the default circumvents the interaction calculations entirely."<<std::endl<<std::endl;

  std::cout << "EXAMPLES" << std::endl << "--------" << std::endl;

  std::cout << "\t" << argv[0] << " 5.5/uniform04.dat" << std::endl
       << "\t" << argv[0] <<" -w 5.5/uniform04.dat" << std::endl
       << "Puts the output in the folder \"5.5\" with the name \"uniform04.dat\"."
       << std::endl;

  std::cout <<"\t"<< argv[0] <<" -li 4" << std::endl
      <<"Use log-spacing, and interaction strength DeltaPhi/(K_b*T)=4."
      <<std::endl;

    std::cout << "\t" << argv[0] <<" -B 1000 -w brute.dat" << std::endl
       << "Generate 1000 output-files: brute.dat0, brute.dat1, etc, using brute force method"
       << std::endl;

  std::cout <<"\t" << argv[0] <<" -li 4 powerlaw.dat" << std::endl
      <<"Since the first argument without flag will be interpreted as \"-w arg.\" this is the same as the previous example, with the additional \"-w powerlaw.dat\" argunment. Note: Be careful when using this, to make sure you are not giving the name as argument to the \"-i\" flag (like \"-i name.out\" would not give an error)."<<std::endl;

  std::cout <<"\t" << argv[0] <<" -li 4 powerlaw01.dat < input01.dat" << std::endl
      <<"Same as previous example but read from the input01.dat-file instead of from the keyboard. The file input01.dat has the same form as what you would type at the keyboard after execution begins." <<std::endl;
  std::cout << std::endl;

  exit(0);

}

//====================================================================
/*#include <fstream>
  #include <iostream>

  using namespace std;

  int main ( int argc, char *argv[] )
  {
  if ( argc != 2 ) // argc should be 2 for correct execution
  // We print argv[0] assuming it is the program name
  std::cout<<"usage: "<< argv[0] <<" <filename>\n";
  else {
  // We assume argv[1] is a filename to open
  ifstream the_file ( argv[1] );
  // Always check to see if file opening succeeded
  if ( !the_file.is_open() )
  std::cout<<"Could not open file\n";
  else {
  char x;
  // the_file.get ( x ) returns false if the end of the file
  //  is reached or an error occurs
  while ( the_file.get ( x ) )
  std::cout<< x;
  }
  // the_file is closed implicitly here
  }
  }


*/
//====================================================================

/*
  void printGnuplot(char name[], int ToGnuPlot, int ensemble, int distribution){
  bool logscale = true;
  ofstream plot;

  double density = Density();
  double D_eff=DiffEffConstant();
  double D_naka=Nakazato();
  double ergo = Ergodicity(latticeX_,latticeY_,latticeZ_);

  plot.open(name);
  //Lite information:
  plot <<"# E = "<<ensemble<<"\t N = "<<N<<"\t X-Y-Z: "
  <<latticeX_<<"x"<<latticeY_<<"x"<<latticeZ_<<"\t concentration="<<density<<std::endl;
  //TODO: skriv ut vilken fördelning det är på jumprates...

  if (writeToFile){
  //set size 0.8, 0.8
  plot<<"set term post color linewidth 1.5"<<std::endl;
  plot<<"set out \""<<name<<".ps\"  "<<std::endl;
  }
  if(logscale) plot<<"set logscale"<<std::endl;

  plot <<"set ylabel \"<r^2>\""<<std::endl;
  plot <<"set xlabel \"t\""<<std::endl;
  plot <<"set legend bottom"<<std::endl;


  //TODO, manick som skriver vad det är för fördelning.
  if (distribution==)

  if (dim==1) plot <<"f(x)=((1-"<<density<<")/"<<density<<")*sqrt(4*"<<D_eff<<"*x/pi)"<<std::endl;

  plot <<"e(x)="<<ergo<<std::endl;

  plot <<"n(x)= "<<D_naka<<"*x"<<std::endl;



  //om vi bill ha kvar plotten på skärmen, vill vi inte stänga filen!
  if (ToGnuPlot!=2) plot<<"quit"<<std::endl;
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



void printError(std::string message){
  std::cout <<"\n Error: " << message << std::endl;
  abort();
}

//use __FILE__ and __LINE__ macro to get line and filename;
void printError(std::string message, int line){
  std::cout <<"\n Error occured at line: " << line
       << std::endl << message << std::endl;
  abort();
}

void printError(std::string message, std::string file, int line){
  std::cout <<"\n Error occured in file " << file << " at line: " << line
       << std::endl << message << std::endl;
  abort();
}

