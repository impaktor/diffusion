#ifndef AUXILARY_H
#define AUXILARY_H

#include <string>

//This file contains implementations that have nothing to do with 
//physics, such as printing messages, ask for input parameters, read
//command arguments (flags), print error messages, ignore commented 
//lines, etc.

//only these functions needs to be accessible from the "outside",
void argumentFlags(int, char**, bool&, bool&, bool&, bool&, int&, char[]);

void AskUserForInputParameters(int&, int&, int&, int&, int&, int&, double&);

void printError(std::string message);
void printError(std::string message, int line);
void printError(std::string message, std::string file, int line);

//only used by AskUserForInputParameters():
void printHelp(char*, char**);

#endif
