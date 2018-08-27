#ifndef FLAT_H
#define FLAT_H

//Header file for AdS/Flat space-time instanton finder backend. This
//code attempts to find an instanton for a potential, including gravitational
//backreaction, when the false vacuum is exactly flat (and the resulting
//space-time inside the bubble locally anti-de-sitter.


//#include "project_specific.h"//For exporting this dll to other
//programs, eg MATLAB.
#include "instanton_solver.h"//For linking in the ode solving library.
#include <iostream>

//Core code needed. Make sure to include the code versions
//so that the compiler is able to select and compile
//as and when it is necessary.
#include "bcadjust_code.h" //includes "bcadjust.h"
#include "events_functions_code.h"
#include "ode_solver_generic_code.h"
#include "potentials_code.h"
#include "ode_rhs_code.h"

//------------------------------------------------------------------------------
//Declaration of the principle function we will use for flat space instanton,
//fixed background approximation:
//typedef boost::array< multi , 3 > FlatSolTypeAction;
template< class value_type , class time_type , class solution_type,
          class solution_type_action >
DLL_EXPORT void odeSolveFlatFixedBackground
    (value_type false_vacuum , value_type true_vacuum , value_type barrier ,
     potential< value_type >& V , time_type chimax , int odeSolverToUse,
     value_type RelTol, value_type AbsTol, value_type stepError, value_type xi,
     value_type lowerBound , value_type upperBound ,
     solution_grid< time_type , solution_type_action , value_type >& solOut,
     value_type& DSout,value_type precision,std::ostream& outStream,
     bool useSimpleBC = false,
     time_type epsilon_step = std::numeric_limits<time_type>::epsilon());
//------------------------------------------------------------------------------
//Solve for a single value:
template< class value_type , class time_type , class solution_type,
          class solution_type_action >
DLL_EXPORT void odeSolveFlatFixedBackgroundSingle
    (value_type y0,value_type false_vacuum , value_type true_vacuum ,
     value_type barrier , potential< value_type >& V , time_type chimax ,
     int odeSolverToUse, value_type RelTol, value_type AbsTol,
     value_type stepError, value_type xi, value_type lowerBound ,
     value_type upperBound ,
     solution_grid< time_type , solution_type_action , value_type >& solOut,
     value_type& DSout,value_type precision,std::ostream& outStream);
//------------------------------------------------------------------------------
#endif //FLAT_H
