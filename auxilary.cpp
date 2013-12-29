#include <iostream>        //for cout, among others
#include <cstdlib>         //for abort()-function
#include <cstring>         //strcpy, strchr 
#include <vector>
#include <string>

#include "auxilary.h"


using namespace std;

void argumentFlags(int argc, char** argv, bool& logarithm, bool& interactionON,
                   bool& quiet, int& interaction_strength, char ReturnFileName[]){
  
  char* file_name = ReturnFileName;
  char* default_file_name = file_name;

  bool changed_file_name = false;
  bool error = false;

  //TODO: flytta in som arg.
  bool low_mem = false;

  int c;
  opterr=0;

  while((c = getopt(argc,argv,"lhmqw:o:i:"))!=-1)
    switch (c){
    case 'l':
      logarithm=true;
      break;
    case 'm':
      low_mem=true;
      break;
    case 'h':
      printHelp(default_file_name,argv);
      break;
    case 'q':
      quiet=true;
      break;
    case 'w':
      file_name = optarg;
      changed_file_name=true;
      break;
    case 'o':
      file_name = optarg;
      changed_file_name=true;
      break;
    case 'i':
      interaction_strength = atoi(optarg);
      interactionON=true;
      if (interaction_strength<0) error=true;
      break;
    case '?':
      error=true;
      if (optopt=='w' || optopt=='o'|| optopt=='i')
        fprintf(stderr, "Option -%c requires an argument. \n",optopt);
        //cout<<"option "<<(char) optopt<<" requires an argument"<<endl;
      else if(isprint(optopt))
        fprintf(stderr,"Unknown option '-%c'. \n",optopt);
      //cout<<"unknown option: "<<(char) optopt<<endl;
      else
        fprintf(stderr,"Unknown option character '\\x%x'.",optopt);
      //cout<<"unknown option character "<<optopt<<endl;
      //return 1;
      break;
    default:
      cout<<"Non-critical error in pogram (bad programming)"<<endl;
      abort();
    }

  //print error message if we entered case '?' above
  if(error){
    cout<<argv[0]<<": undefined entry."<<endl
        <<"Try \""<<argv[0]<<" -h\" for help."<<endl;
    exit(1);
  }

  //if non of the above specified arguments were passed:
  //take the first as new file_name, and print the rest as 
  //errors
  int index;
  bool more_than_one_invalid_argument = (changed_file_name) ? true : false;

  for(index=optind; index < argc; index++){
    //only run first iteration & if we didn't change filename (with w|f) 
    if(!more_than_one_invalid_argument){
      file_name=argv[index];
      changed_file_name=true;
    }

    //run for all invalid arguments except the first one
    if(more_than_one_invalid_argument)
      printf("Non-option argument %s\n",argv[index]);

    more_than_one_invalid_argument=true;
  }

  if(changed_file_name)
    cout<<"### Default output file changed to: "<<file_name<<endl;
  else
    file_name=default_file_name;

  //cout<<endl<<"#TEST001: filename="<<(string) file_name<<endl;

  //ReturnFileName=file_name;

  strcpy(ReturnFileName,file_name);

}


void AskUserForInputParameters(int& xdim, int& ydim, int& zdim,
                               int& particleNumber, int& ensembles, 
                               int& number_of_values, double& stopTime){
  cout <<"# Specify number of lattice sites in:"<<endl;
  cout <<"# X: ";
  xdim = getNonCommentIntInput(); //ignore lines containing '#'
  cout<<"# Y: ";
  ydim = getNonCommentIntInput(); //ignore lines containing '#'
  cout<<"# Z: ";
  zdim = getNonCommentIntInput(); //ignore lines containing '#'

  cout <<"# Total number of particles: ";
  particleNumber = getNonCommentIntInput();

  float info = 1.0*(xdim*ydim*zdim)/particleNumber;
 
  cout <<"# Number of vacant sites/particle: "<<info<<endl;
  cout <<"# Number of ensembles: ";
  ensembles = getNonCommentIntInput();
  cout <<"# Number of data-values to save: ";
  number_of_values=getNonCommentIntInput();
  cout <<"# Stop-time: ";
  stopTime = getNonCommentDoubInput();
  cout << endl<<endl;
  
}


void printHelp(char* File, char** argv){
  char* default_file_name=File;
    cout<<"USAGE"<<endl<< "-----"<<endl;
    cout<<"\t"<<argv[0]<<" [OPTIONS]"<<endl<<endl;
    
    cout<<"OPTIONS"<<endl<<"-------"<<endl;
    cout<<"-i n\t Set interactions to ON with interaction strength \"n\".\n\t Note: n>=0, since it will be used in the Boltzmann step as \"exp(-n)\"."<<endl;
    cout<<"-l \t Set log-spacing between data values to be written to\n\t output-file to ON."<<endl;
    cout<<"-w NAME\t Set output-file to 'NAME'."<<endl;
    cout<<"NAME\t Set first non-valid options-flag as name of output-file.\n (shorter/faster version than using the '-w' option)"<<endl;

    cout<<endl<<"Arguments can be given in any order. If no flags are set the first argument following \""<<argv[0]<<"\" will be taken as the output-file name. With no arguments given the default is log-spacing of values to output-file \""<<default_file_name<<"\" without any interactions at all. Note, with '-i 0' the interaction strength will be 0, but the interaction code will be turned ON, resulting in a slower simulation run, therefore the default circumvents the interaction calculations entirely."<<endl<<endl;

    cout<<"EXAMPLES"<<endl<<"--------"<<endl;
    
    cout<<"\t"<<argv[0]<<" 5.5/uniform04.dat"<<endl
        <<"\t"<<argv[0]<<" -w 5.5/uniform04.dat"<<endl
        <<"Puts the output in the folder \"5,5\" with the name \"uniform04.dat\"."
        <<endl;
    
    cout<<"\t"<<argv[0]<<" -li 4"<<endl
        <<"Use log-spacing, and interaction strength DeltaPhi/(K_b*T)=4."
        <<endl;

    cout<<"\t"<<argv[0]<<" -li 4 powerlaw.dat"<<endl
        <<"Since the first argument without flag will be interpreted as \"-w arg.\" this is the same as the previous example, with the additional \"-w powerlaw.dat\" argunment. Note: Be careful when using this, to make sure you are not giving the name as argument to the \"-i\" flag (like \"-i name.out\" would not give an error)."<<endl;

    cout<<"\t"<<argv[0]<<" -li 4 powerlaw01.dat <<input01.dat"<<endl
        <<"Same as previous example but read from the input01.dat-file instead of from the keyboard. The file input01.dat has the same form as what you would type at the keyboard after execution begins." <<endl;
    cout<<endl;
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
    cout<<"usage: "<< argv[0] <<" <filename>\n";
  else {
    // We assume argv[1] is a filename to open
    ifstream the_file ( argv[1] );
    // Always check to see if file opening succeeded
    if ( !the_file.is_open() )
      cout<<"Could not open file\n";
    else {
      char x;
      // the_file.get ( x ) returns false if the end of the file
      //  is reached or an error occurs
      while ( the_file.get ( x ) )
        cout<< x;
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


int getNonCommentIntInput(void){
  //small function to allow comments in the input-file.
  //this will ignore any line containing "#", since it
  //only wants numbers.

  char invalue[999];

  for(bool repeat=true; repeat;){
    cin>>invalue;
    if(strchr(invalue,'#') == NULL)
      repeat=false;
  }

  return atoi(invalue);
}

double getNonCommentDoubInput(void){
  //small function to allow comments in the input-file.
  //this will ignore any line containing "#", since it
  //only wants numbers.

  char invalue[999];

  for(bool repeat=true; repeat;){
    cin>>invalue;
    if(strchr(invalue,'#') == NULL)
      repeat=false;
  }

  return atof(invalue);
}


void printError(string message){
  cout << message << endl;
  abort();
}

