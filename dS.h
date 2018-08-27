#ifndef DS_H_INCLUDED
#define DS_H_INCLUDED

//Header file for dS space-time instanton finder backend. This
//code attempts to find an instanton for a potential, including gravitational
//backreaction, when the false vacuum is de Sitter.


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
DLL_EXPORT void odeSolve_dS(value_type false_vacuum , value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutLHS,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutRHS,
                     value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     std::vector<value_type>& other_return_values,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,bool graded_precision,
                     value_type a_scale,value_type overshoot_lower,
                     value_type overshoot_upper,bool userSuppliedBounds,
                     value_type LHS_bound_lower,value_type LHS_bound_upper,
                     value_type RHS_bound_lower,value_type RHS_bound_upper,
                     bool reverse_shooting,bool use_analytic_a,
                     bool fixed_background);
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dSsingle(value_type y0, value_type false_vacuum ,
                           value_type true_vacuum , value_type barrier ,
                           potential< value_type >& V , time_type chimax ,
                           int odeSolverToUse,
                           value_type RelTol, value_type AbsTol,
                           value_type stepError, value_type xi,
                           solution_grid< time_type , std::vector<value_type> ,
                           value_type >& solOut,value_type& DSout,
                           value_type precision,value_type h,value_type W0,
                           std::ostream& outStream,
                           bool track_scale_factor = true,
                           bool linear_BCs = false,
                           bool fixed_background = false);
//------------------------------------------------------------------------------
//Function to find the bounce for a single side (ie, don't attempt to find
//the other side, just do overshoot/undershoot on one side
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS_one_sided(value_type false_vacuum ,
                     value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     value_type& y_lower,value_type& y_upper,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,
                     bool graded_precision,
                     value_type a_scale,
                     value_type overshoot_lower,
                     value_type overshoot_upper,
                     bool userSuppliedBounds,
                     value_type LHS_bound_lower,
                     value_type LHS_bound_upper,
                     bool reverse_shooting,bool fixed_background);
//------------------------------------------------------------------------------
//Function which checks if we cross a particular threshold and integrates
//differently if so:
template<class time_type,class solution_type,class value_type>
int ode_solve_trans(odeRHS< solution_type, time_type, value_type >& odefun,
              solution_type& y0, std::vector< time_type >& tspan,
              events_function< solution_type, time_type, value_type >& events,
              int odeSolver,
              solution_grid< time_type, solution_type, value_type >& sol,
              value_type RelTol, value_type AbsTol,time_type initStepMaxFrac,
              value_type thresh,value_type y0_low,value_type y0_high,
              potential<value_type>& V,value_type V0,bool include_action,
              value_type xi,value_type h,bool true_vacuum_side,
              bool ap0termination,value_type fv);
//------------------------------------------------------------------------------
//Overloaded versions of odeSolve_dS. Used because we want to specify some
//default parameters to have their default values equal to a previous value.
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS(value_type false_vacuum , value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutLHS,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutRHS,
                     value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     std::vector<value_type>& other_return_values,
                     std::ostream& outStream,bool track_scale_factor = true,
                     bool linear_BCs = false);
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS(value_type false_vacuum , value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutLHS,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutRHS,
                     value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     std::vector<value_type>& other_return_values,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,bool graded_precision = false,
                     value_type a_scale = value_type(1.0));
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS(value_type false_vacuum , value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutLHS,
                     solution_grid< time_type , std::vector< value_type > ,
                     value_type >& solOutRHS,
                     value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     std::vector<value_type>& other_return_values,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,bool graded_precision,
                     value_type a_scale,value_type overshoot_lower,
                     value_type overshoot_upper,bool userSuppliedBounds =false);
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS_one_sided(value_type false_vacuum ,
                     value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     value_type& y_lower,value_type& y_upper,
                     std::ostream& outStream,bool track_scale_factor = true,
                     bool linear_BCs = false);
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS_one_sided(value_type false_vacuum ,
                     value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     value_type& y_lower,value_type& y_upper,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,
                     bool graded_precision = false,
                     value_type a_scale = value_type(1.0));
//------------------------------------------------------------------------------
template< class value_type , class time_type >
DLL_EXPORT void odeSolve_dS_one_sided(value_type false_vacuum ,
                     value_type true_vacuum ,
                     value_type barrier , potential< value_type >& V ,
                     time_type chimax_suggestion , int odeSolverToUse,
                     value_type RelTol, value_type AbsTol, value_type stepError,
                     value_type xi,
                     solution_grid< time_type , std::vector<value_type> ,
                     value_type >& solOut,value_type& DSout,
                     value_type precision,value_type h,value_type W0,
                     value_type& y_lower,value_type& y_upper,
                     std::ostream& outStream,bool track_scale_factor,
                     bool linear_BCs,value_type action_reltol,
                     value_type action_abstol,bool graded_precision,
                     value_type a_scale,value_type overshoot_lower,
                     value_type overshoot_upper,bool userSuppliedBounds= false);
//------------------------------------------------------------------------------
#endif // DS_H_INCLUDED
