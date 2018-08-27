#ifndef INSTANTON_SOLVER_MAIN_PROJECT_SPECIFIC
#define INSTANTON_SOLVER_MAIN_PROJECT_SPECIFIC

//Fix some problems with boost containers on GNU:
#define FUSION_MAX_VECTOR_SIZE 13

#define STATIC_EXPORT

//Switches to determine which part of the code will be used.
#define USING_POTENTIAL_TEST1//Our logarithmic test potential.
#define USING_POTENTIAL_BRANCHINA//For testing Branchina's potential.
#define USING_POTENTIAL_NEW_PHYSICS
//#define DE_SITTER_INSTANTONS //When we are considering only the de-Sitter
//case.
#define ANTI_DE_SITTER_FLAT_INSTANTONS //When we want to look at instantons
//with a flat space
//false vacuum or anti-de-Sitter false vacuum.
#define DE_SITTER_INSTANTONS //Indicates that we are using the de-Sitter
    //solver.
#define USING_SM_POTENTIAL//For studying the Standard Model.

//Defines the arbitrary precision arithmetic we are using.
#include "multi_precision_definitions.h"

#endif
