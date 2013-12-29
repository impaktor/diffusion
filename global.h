#ifndef GLOBAL_H
#define GLOBAL_H
//This is for constants used when killing bugs or checking a program run.


//boolean, prints simulation details to screen, 
//sets the private variable Lattice::testOnOff_
const bool TEST_1            = 0;

//CPU consuming! Use only when bug-squashing!
//Is used in Lattice::checkVacancyMatrix().
const bool CHECK_VACANCY_ON  = 0;

//used in superInteraction.cpp
const bool VERBOSE = 0;


//Used in save.cpp 
//Print <x>,<y>,<z> (un-squared) to out.dat_txyz (used for testing.)
const bool PRINT_SECOND_FILE = true;

#endif
