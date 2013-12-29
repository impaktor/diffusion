#ifndef AUXILARY_H
#define AUXILARY_H

#include <string>

//This file contains implementations that have nothing to do with
//physics, such as printing messages, ask for input parameters, read
//command arguments (flags), print error messages, ignore commented
//lines, etc.

//only these functions needs to be accessible from the "outside",
void argumentFlags(int, char**, bool&, bool&, bool&, bool&, int&,
                   std::string&, std::string&, int&, char&, bool&, float&);

void askUserForInputParameters(int&, int&, int&, int&,
                               int&, int&, double&);

void readInputFile(std::string&, int&, int&, int&, int&,
                   int&, int&, double&);

//Prints message, and aborts the program.
void printError(std::string message);
void printError(std::string message, int line);
void printError(std::string message, std::string file, int line);

//only used by AskUserForInputParameters():
void printHelp(std::string, char**);

#endif
