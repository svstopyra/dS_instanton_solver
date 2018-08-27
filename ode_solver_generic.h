#ifndef ODE_SOLVER_GENERIC_INCLUDED
#define ODE_SOLVER_GENERIC_INCLUDED

//DLL stuff - determines whether the header is being used by the dll itself
//(in which case we are exporting) or
//an application which references it (in which case we are importing):
#ifdef DLL_EXPORT
#define DLL_EXPORT __declspec(dllexport)
#else
#ifdef STATIC_EXPORT
#define DLL_EXPORT //If we specifically do not want to use either of the dll
//options, and just want to compile the functions normally.
#else
#define DLL_EXPORT __declspec(dllimport)
#endif
#endif

//Core code definitions:
#include "events_functions.h"
#include <boost/numeric/odeint.hpp>

//Needed libraries:
#include <vector>
//#include "ode_rhs_class_definition.h"
//#include "event_function_definition.h"

//------------------------------------------------------------------------------
//Class to store solution and its events:
template<class time_type,class solution_type,class value_type>
class solution_grid
{
public:
    std::vector< time_type > T; //Time grid
    std::vector< solution_type > Y; //Solution grid
    std::vector< time_type > TE; // Time of detected events
    std::vector< solution_type > YE; //Value at detected events
    std::vector< int > IE; //Index of the events function which fired.
    std::vector< int > eventLocations;//Locations of the events in the grids
        //T and Y, so we can identify them later.
    int nSteps;
    bool terminated_early;
	int nProvisional;//Everything past this point is provisional data, and may
        //be deleted if the solver detects an error
	int nProvisionalEvents;
    //Default constructor:
	solution_grid();
    //Member function to erase all data, but retain the struct, effectively
    //refreshing
    //it to how it started.
	DLL_EXPORT void erase();
	DLL_EXPORT void refresh();//Deletes all the provisional data if we
        //decide we don't need it, but leaves the
        //non-provisional data intact.
    value_type S_correction;//Correction to the action that is needed by some
        //solution methods.
};
//------------------------------------------------------------------------------
template<class time_type, class solution_type, class value_type>
solution_grid< time_type, solution_type, value_type >::solution_grid()
{
	nSteps = 0;
	terminated_early = false;
	nProvisional = 0;
	nProvisionalEvents = 0;
}
//------------------------------------------------------------------------------
//Template function which computes a solution given the input starting values,
//a particular rhs of the ode, a tspan.
template< class time_type , class solution_type , class value_type >
DLL_EXPORT int ode_solve(odeRHS< solution_type , time_type , value_type >&,
                         solution_type&, std::vector< time_type >& ,
                         events_function< solution_type, time_type,
                                        value_type >& ,
                         int ,solution_grid< time_type ,solution_type,
                                            value_type >& ,
                        value_type ,value_type,time_type,
                        bool refining_events = true);//value_type RelTol
                        //= value_type(1e-6),
                        //value_type AbsTol = value_type(1e-6))
//------------------------------------------------------------------------------
//Version which switches to an analytic form for a(t) when a certain threshold
//is passed.
template<class time_type,class solution_type,class value_type>
int ode_solve_trans(odeRHS< solution_type, time_type, value_type >& odefun,
              solution_type& y0, std::vector< time_type >& tspan,
              events_function< solution_type, time_type, value_type >& events,
              int odeSolver,
              solution_grid< time_type, solution_type, value_type >& sol,
              value_type RelTol, value_type AbsTol,time_type initStepMaxFrac);
//------------------------------------------------------------------------------
#endif
