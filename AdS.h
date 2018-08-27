#ifndef ADS_H
#define ADS_H

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
//Declaration of the principle function we will use for AdS/Flat space
//instantons, including backreaction.
//typedef boost::array< multi , 5 > AdSFlatSolTypeAction;
template< class value_type , class time_type >
DLL_EXPORT void odeSolveAdSFlat
    (value_type false_vacuum , value_type true_vacuum , value_type barrier ,
     potential< value_type >& V , time_type chimax , int odeSolverToUse,
     value_type RelTol, value_type AbsTol, value_type stepError, value_type xi,
     value_type lowerBound , value_type upperBound ,
     solution_grid< time_type , std::vector<value_type> , value_type >& solOut,
     value_type& DSout,value_type precision,std::ostream& outStream,
     bool track_scale_factor = true,bool compute_H0_linearisation = false);
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolveAdSFlatSingle
    (value_type y0, value_type false_vacuum , value_type true_vacuum ,
     value_type barrier , potential< value_type >& V , time_type chimax ,
     int odeSolverToUse, value_type RelTol, value_type AbsTol,
     value_type stepError, value_type xi, value_type lowerBound ,
     value_type upperBound ,
     solution_grid< time_type , std::vector<value_type> , value_type >& solOut,
     value_type& DSout,value_type precision,std::ostream& outStream,
     bool track_scale_factor = true);
//------------------------------------------------------------------------------
#endif // ADS_H
