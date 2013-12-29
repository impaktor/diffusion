#ifndef AUXILARY_H
#define AUXILARY_H

#include <string>

//This file contains implementations that have nothing to do with 
//physics, such as print messages, ask for input parameters, read
//command arguments (flags), print error messages.

void printHelp(char*, char**);

void argumentFlags(int, char**, bool&, bool&, bool&, int&, char[]);

void AskUserForInputParameters(int&, int&, int&, int&, int&, int&, double&);

int getNonCommentIntInput(void);

double getNonCommentDoubInput(void);

void printError(std::string);

#endif
