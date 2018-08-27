#ifndef ODE_SOLVER_GENERIC_CODE_H
#define ODE_SOLVER_GENERIC_CODE_H

//Needed libraries
#include <cstdlib>
#include <boost/array.hpp>

#include "ode_solver_generic.h"
//#include "event_detection.h"
//#include "event_function_definition.h"
#include "find_zero.h"

//------------------------------------------------------------------------------
//Similar definitions for events functions:
template< class solution_type, class time_type, class value_type >
void events_function< solution_type, time_type, value_type >::record_event
(value_type* value, bool* isterminal, int* direction)
{


	for(int i = 0; i < nEvents; i++)
	{
		values[i] = value[i];
		isterminals[i] = (isterminal[i] && terminate_on) && dynamicSwitch[i];
		directions[i] = direction[i];
		//std::cout << "\nvalues[" << i << "] = " << values[i];
	}
}
//------------------------------------------------------------------------------
//Code for event detection:
//eventLogger constructor
template<class time_type, class solution_type>
eventLogger< time_type, solution_type >::eventLogger
    (std::vector< solution_type >& ye_before,
     std::vector< solution_type >& ye_after,
	std::vector< time_type >& te_before, std::vector< time_type >& te_after,
	std::vector< int >& ie, int& events_found)
	: YE_before(ye_before), YE_after(ye_after), TE_before(te_before),
	TE_after(te_after), IE(ie), nEventsFound(events_found)
{
	//Start with no events. Vectors are defined as above, starting as empty.
	nEventsFound = 0;
}
//------------------------------------------------------------------------------
//output_observer_with_events constructor:
template< class time_type, class solution_type, class value_type >
output_observer_with_events< time_type, solution_type, value_type >
::output_observer_with_events
(int& steps_in, std::vector< solution_type >& Y_input,
	std::vector< time_type >& T_input,
	events_function< solution_type, time_type, value_type >& event_function,
	eventLogger< time_type, solution_type >& Event_Log, bool& terminate_early)
	: nSteps(steps_in), Y(Y_input), T(T_input), event_log(Event_Log),
	events(event_function),
	terminate(terminate_early)
{
	use_termination = true;
	//refining_mode = false;
	//Make sure there is room to store information about the events:
	/*
	for(int k = 0;k < this->events.values.size();k++)
    {
        std::cout << "\nevents.values[" << k << "] = " << this->events.values[k];
    }
    */
    this->valueLast = events.values;
    /*
    for(int k = 0;k < this->valueLast.size();k++)
    {
        std::cout << "\nvalueLast[" << k << "] = " << this->valueLast[k];
    }
    */
    this->isterminalLast = events.isterminals;
    this->directionLast = events.directions;
}
//------------------------------------------------------------------------------
//observer function. Returns whether execution should terminate early
//(false by default)
template< class time_type, class solution_type, class value_type >
bool output_observer_with_events< time_type, solution_type, value_type >
::operator()(const solution_type &x, const time_type t)
{
	bool terminate_early = false;
	/*
	//Check whether or not the file is open, and if not, open it:
	if(!(*output_file).is_open())
	{
	(*output_file).open(filename_.c_str(), ios::app);
	}
	//Proceed with output in csv format, so MATLAB can understand it:
	(*output_file) << t << "," << x[0] << ',' << x[1] << ',' << x[2] << endl;
	*/
	//output directly to MATLAB, by storing data first:

	//int old_elements_y = nSteps*Components;

	//Evaluate the events function and check each event type for new events:
	events(t, x);
    Y.push_back(x);
    T.push_back(t);
    nSteps++;
    /*
    if(events.dynamicSwitch[0])//Only switched on for the final solution.
    {
        std::cout << "\nYlast = ";
        for(int i = 0;i < int(Ylast.size());i++)
        {
            std::cout << "\nYlast[" << i << "] = " << Ylast[i];
        }
        std::cout << "\nvalueLast = ";
        for(int i = 0;i < int(valueLast.size());i++)
        {
            std::cout << "\nvalueLast[" << i << "] = " << valueLast[i];
        }
        std::cout << "\nY(t) = ";
        for(int i = 0;i < int(x.size());i++)
        {
            std::cout << "\nYlast[" << i << "] = " << x[i];
        }
        std::cout << "\nvalues = ";
        for(int i = 0;i < int(events.values.size());i++)
        {
            std::cout << "\nvalues[" << i << "] = " << events.values[i];
        }
    }*/
    //std::cout << "\nI have called the events function and saved the data "
    //          << nSteps << " times";

	value_type _0p0 = value_type(0.0);
	if (nSteps > 0)
	{
		for (int i = 0; i < events.nEvents; i++)
		{
			//Check for events:
			bool cond_decrease = (valueLast[i] > _0p0 && events.values[i] <_0p0)
                                 && directionLast[i] == -1;
			bool cond_increase = (valueLast[i] < _0p0 && events.values[i] >_0p0)
                                 && directionLast[i] == 1;
			bool cond_either = ((valueLast[i] > _0p0 && events.values[i] < _0p0)
                                || (valueLast[i] < _0p0 &&
                                events.values[i] > _0p0))
                                 && directionLast[i] == 0;
			if ((cond_decrease || cond_increase || cond_either))
			{

			    std::cout << "\nFound event!!!";
			    std::cout << "\nIE[i] = " << i + 1 << " valueLast[i] = " <<
                        valueLast[i] << " events.values[i]  = "
                        << events.values[i]
                        << "\n y = ";
                for(int k = 0;k < x.size();k++)
                {
                    std::cout << "\ny[" << k << "] = " << x[k];
                }
				//found an event! Record data about its location,
				//solution value and type of event:
				//mexPrintf("New Event!\n Ylast: %.10f\n Ynext: %.10f\n Tlast
                //%.10f\n Tnext %.10f\n",Ylast[0],x[0],Tlast,t);
				event_log.YE_before.push_back(Ylast);
				event_log.YE_after.push_back(x);
				event_log.TE_before.push_back(Tlast);
				event_log.TE_after.push_back(t);
				event_log.IE.push_back(i + 1);
				event_log.nEventsFound++;
				//Is the event terminal?
				//std::cout << "\nisterminals[" << i << "] = " << events.isterminals[i];
                //std::cout << "\ndynamicSwitch[" << i << "] = " << events.dynamicSwitch[i];
				terminate_early = (events.isterminals[i]) && use_termination;
				//Add the event to the list of solutions, and record it so we
				//can later improve it:
				nSteps++;
				event_log.EventLocations.push_back(nSteps);
				//Y.push_back(Ylast);
				//T.push_back(Tlast);
			}
		}
		if (!terminate_early)
        {
            this->valueLast = events.values;
            this->isterminalLast = events.isterminals;
            this->directionLast = events.directions;
        }
	}
	/*
	for (int i = 0; i < events.nEvents; i++)
	{
		//Store the result of the events function evaluation:
		valueLast[i] = events.values[i];
		isterminalLast[i] = events.isterminals[i];
		directionLast[i] = events.directions[i];
	}*/

    //Fill new rows with data, IF the event is not terminal:
    /*
    if (!terminate_early)
    {
        Y.push_back(x);
        T.push_back(t);
        nSteps++;
    }
    */
    /*
    std::cout << "\nSuccessfully called events function for: ";
    std::cout << "\nt = " << t;
    for(int i = 0;i < x.size();i++)
    {
        std::cout << "\ny[" << i << "] = " << x[i];
    }*/
	Ylast = x;
	Tlast = t;
	terminate = terminate || terminate_early;
	return terminate_early;
}
//------------------------------------------------------------------------------
//Constructor for event_refiner class:
template< class time_type, class solution_type, class value_type,
          class Stepper >
event_refiner< time_type, solution_type, value_type, Stepper >
::event_refiner
    (int nEvent, time_type t_lower, solution_type& Y_lower,
     odeRHS< solution_type, time_type, value_type >& ode_to_solve,
     events_function< solution_type, time_type, value_type >& Events,
     Stepper stepper_to_use)
     : Y(Y_lower), ode(ode_to_solve), events(Events), stepper(stepper_to_use)
{
	T = t_lower;
	event = nEvent;
}
//------------------------------------------------------------------------------
template< class time_type, class solution_type, class value_type,
          class Stepper >
value_type event_refiner< time_type, solution_type, value_type, Stepper >
::operator()(time_type t)
{
	odeCaller<solution_type, time_type, value_type > odeCall(&ode);
	//Integrate up to t:
	solution_type Y_guess = Y_init;
	time_type T_guess = T_init;
	/*
	std::cout.precision(50);
	std::cout << "\nRefining event IE = " << event << ".\n";
	for(int i = 0;i < int(Y_low.size());i++)
    {
        std::cout << "\nY_low[" << i << "] = " << Y_low[i];
    }*/
	boost::numeric::odeint::integrate_adaptive< Stepper,
                                                odeCaller<solution_type,
                                                          time_type,
                                                          value_type >,
                                                solution_type, time_type,
                                                null_obs<time_type,
                                                solution_type> >
		(stepper, odeCall, Y_guess, T_guess, t, time_type(0.1)*(t - T_guess),
         null_obs_event);
	//evaluate the events function:
	/*std::cout << "\nY(t) = ";
	for(int i = 0;i < int(Y_low.size());i++)
    {
        std::cout << "\nY_low[" << i << "] = " << Y_low[i];
    }*/
	events(t, Y_guess);
	//std::cout << "\nevents.values = " << events.values[event - 1];
	//T = t;
	//Record the value/time of the event:
	Y = Y_guess;
	T = T_guess;
	return events.values[event - 1];
}
//------------------------------------------------------------------------------
template< class time_type, class function, class value_type >
time_type bisection(time_type& lower, time_type& upper, function& f,
                    value_type precision)
{
	//Attempts to locate the zero of some function, known to lie between
	//lower and upper, to a given precision

	//time is some time variables, y a function return variable and
	//function f a class with an operator().
	time_type guess = (upper + lower) / 2;
	value_type flower = f(lower);
	value_type fupper = f(upper);
	value_type fguess = f(guess);
	//mexPrintf("Refine zero:\n");
	//threshold. Added this line to prevent early exits with bad precision. We
	// need to know the PRECISE location of the zero for our work:
	value_type epsilonThresh = std::numeric_limits< value_type >::epsilon();
	while ((fguess > precision || fguess < -precision)
            && (upper - lower) > epsilonThresh)
	{
		//mexPrintf("%.15f\n %.15f\n %.10f\n %.10f\n",fguess,guess,lower,upper);
		if (fguess*fupper > 0)
		{
			upper = guess;
		}
		else if (fguess*flower > 0)
		{
			lower = guess;
		}
		guess = (upper + lower) / 2;
		fguess = f(guess);
	}
	return guess;
}
//------------------------------------------------------------------------------
//Function to refine events
template<class time_type, class solution_type,
class value_type , class Stepper >
	refined_events<time_type, solution_type>
	refineEvents(eventLogger< time_type, solution_type >& eventLog,
		output_observer_with_events< time_type,
		solution_type, value_type >& record_data,
		odeRHS< solution_type, time_type, value_type >& ode_rhs,
		value_type precision, Stepper stepper)
{
	//First switch off event termination:
	record_data.use_termination = false;

	//Iterate through events and refine their locations:
	refined_events< time_type, solution_type > events_final(eventLog.IE);
	typedef event_refiner< time_type, solution_type, value_type, Stepper > er;
	for (int i = 0; i < eventLog.nEventsFound; i++)
	{
		time_type lower = eventLog.TE_before[i];
		time_type upper = eventLog.TE_after[i];
		solution_type solution = eventLog.YE_before[i];
		er refiner(eventLog.IE[i], lower, solution, ode_rhs,
                   record_data.events, stepper);
		//Refine the event's location:
		events_final.TE.push_back(bisection< time_type, er, value_type >
                                  (lower, upper, refiner, precision));
		events_final.YE.push_back(solution);//Event refiner modifies
                                            //solution to be
											//at the event in question.
	}
	//Switch termination back on again, in case we want to continue using it.
	//Otherwise, users may not realise this has been turned off and forget to
	//switch it back on again.
	record_data.use_termination = true;
	return events_final;
}
//------------------------------------------------------------------------------
//Code for ode solving:
//Member function definitions for solution_grid:
//Member function to erase all data, but retain the struct, effectiveyl
//refreshing it to how it started.
template<class time_type, class solution_type, class value_type>
void solution_grid< time_type, solution_type, value_type >::erase()
{
	T.clear();
	Y.clear();
	TE.clear();
	YE.clear();
	IE.clear();
	eventLocations.clear();
	nSteps = 0;
	terminated_early = false;
}
//------------------------------------------------------------------------------
//Member function to delete all provisional data from the solution grid:
template<class time_type, class solution_type, class value_type>
void solution_grid< time_type, solution_type, value_type >::refresh()
{
	T.erase(T.begin() + nProvisional, T.end());
	Y.erase(Y.begin() + nProvisional, Y.end());
	TE.erase(TE.begin() + nProvisionalEvents, TE.end());
	YE.erase(YE.begin() + nProvisionalEvents, YE.end());
	IE.erase(IE.begin() + nProvisionalEvents, IE.end());
	eventLocations.erase(eventLocations.begin() + nProvisionalEvents, IE.end());
}
//------------------------------------------------------------------------------
//Code for the templated ode_solve function:
template<class time_type, class solution_type, class value_type >
int ode_solve(odeRHS< solution_type, time_type, value_type >& odefun,
              solution_type& y0, std::vector< time_type >& tspan,
              events_function< solution_type, time_type, value_type >& events,
              int odeSolver,
              solution_grid< time_type, solution_type, value_type >& sol,
              value_type RelTol, value_type AbsTol,time_type initStepMaxFrac,
              bool refining_events)
              //value_type RelTol = value_type(1e-6),
              //value_type AbsTol = value_type(1e-6))
{
	/*
	*TEMPLATES:
	*time_type - type of variable used to measure time. e.g., double.
	*value_type - type of variable the components of the solution are
	*constructed from, eg, doubles.
	*solution_type - object which the solution is, eg, a vector of value_types,
	*or matrix of value_types.
	*odeRHS - a class with an operator() with arguments:
	*      (const solution_type &y , solution_type &dydt , time_type t)
	*      y - solution at time t
	*      t - time
	*      dydt - derivative at time t (ie, the ode we are solving).
	*      No return value is required - should be void.
	*event_handler - a class with an operator(), with arguments:
	*      (const time_type t,const solution_type& y)
	*      should return an event_return<value_type> struct.
	*/
	/*
	*INPUTS:
	*odefun - class with the operator() which evaluates the rhs of the ode.
	*y0 - initial condition
	*tspan - times at which to evaluate the solution. tspan[0] is the start
	*time, tspan[end] the end time.
	*  anything in between will be evaluated, but additional points
	*corresponding to those selected by the
	*  adaptive step size will also be included.
	*events - events function to evaluate at each step, searching for events
	*as we go and possibly terminating execution.
	*odeSolver - which method we should use to solve the ode.
	*  1 - Dormand-Prince45 method.
	*  2 -
	*  default - Dormand-Prince45 method.
	*sol - struct to which we output the data of the solution.
	*  Should be constructed and passed by reference to avoid having to keep
	* copying it.
	*
	*/


	//Wrapper to indirectly call the ode_solver, since the boost integrator
	//complains when we try to compile
	//it with a virtual class, due to attempting to instantiate it, for some
	//reason.
	//time_type initStepMaxFrac = time_type(0.01)
	odeCaller<solution_type, time_type, value_type> oderhs_caller(&odefun);

	//VERIFY INPUT:
	int success = 1;//Determines whether the program succeeded or not.
					/*
					*0 - failed for unknown reasons
					*1 - succeeded
					*-1 - failed because tspan had fewer than two elements.
					*-2 - failed because tspan was not monotonic.
					*-3 - integration failed because smaller step-sizes could
                            not be used.
					*-4 - integration failed by exceeding maximum number of
                            iterations.
					*-5 - integration failed due to attempted division by zero
					*-6 - integration failed due to encountering NaN value.
					*/
	bool acceptable_tspan = true;
	int nTimeSteps = int(tspan.size());
	if (tspan.size() < 2)
	{
		//Invalid tspan
		//std::cout << "tspan requires at least two elements.\n";
		acceptable_tspan = false;
		success = -1;
		throw "tspan requires at least two elements.\n";
	}
	else
	{
		//Check monotonic:
		bool increasing = tspan[0] < tspan[1];
		bool decreasing = tspan[0] > tspan[1];
		if (!(increasing || decreasing))
		{
			//Invalid tspan:
			success = -2;
			throw "tspan requires increasing or decreasing elements.\n";
		}
		//Check tspan is monotonically increasing:
		else if (increasing && int(tspan.size()) > 2)
		{
			for (int i = 1; i < tspan.size(); i++)
			{
				acceptable_tspan = (acceptable_tspan)
                                    && (tspan[i - 1] < tspan[i]);
				if (!acceptable_tspan)
				{
					success = -2;
					throw "Monotonically increasing tspan expected..\n";
				}
			}
		}
		//check tspan is monotonically decreasing:
		else if (decreasing && int(tspan.size()) > 2)
		{
			for (int i = 1; i < tspan.size(); i++)
			{
				acceptable_tspan = (acceptable_tspan)
                                    && (tspan[i - 1] > tspan[i]);
				if (!acceptable_tspan)
				{
					success = -2;
					throw "Monotonically decreasing tspan expected.\n";
				}
			}
		}
	}

	//Now perform the integration.
	solution_type y = y0;
	time_type t0 = tspan[0];

	//Setup event logging. We first get the rough position, then refine after/
	//solving the equation.
	//Event Data:
	std::vector< solution_type > YE_before; //solution just before an event.
	std::vector< time_type > TE_before; //time just before event.
	std::vector< solution_type > YE_after; //solution just after an event.
	std::vector< time_type > TE_after; //time just after event
									   //Indices of events are are stored in
                                    //sol.IE
	int nEventsFound; //total # of events we found.
					  //Test the events function to find out how many types of
					  //events there are:
	int nNoOfEvents = events.nEvents;
	//Create a log to store the events:
	eventLogger< time_type, solution_type > event_log(YE_before, YE_after,
		TE_before, TE_after, sol.IE, nEventsFound);


	//Construct an observer to deal with recording data. Pass the above
	//variables by reference so they are available out of scope:
	typedef output_observer_with_events< time_type, solution_type, value_type >
            oowe;
    events(tspan[0],y0);//Setup initial values to compare with.
	oowe record_data(sol.nSteps, sol.Y, sol.T, events, event_log,
                     sol.terminated_early);
    //Give the event detector something to compare with!
    record_data.Ylast = y0;
    record_data.Tlast = tspan[0];


	//typedefs for the possible ode solvers:
	//Cash-Karp method:
	typedef boost::numeric::odeint::runge_kutta_cash_karp54
        < solution_type, value_type, solution_type, time_type >
            error_stepper_type_cash_karp54;
	typedef boost::numeric::odeint::controlled_runge_kutta
        < error_stepper_type_cash_karp54 >
            controlled_stepper_type_cash_karp54;
	//Dormand-Prince:
	typedef boost::numeric::odeint::runge_kutta_dopri5
        < solution_type, value_type, solution_type, time_type >
            error_stepper_type_dorpri45;
	typedef boost::numeric::odeint::controlled_runge_kutta
        < error_stepper_type_dorpri45 > controlled_stepper_type_dorpri45;
	//Fehlberg78:
	typedef boost::numeric::odeint::runge_kutta_fehlberg78
        < solution_type, value_type, solution_type, time_type >
            error_stepper_type_fehlberg78;
	typedef boost::numeric::odeint::controlled_runge_kutta
        < error_stepper_type_fehlberg78 > controlled_stepper_type_fehlberg78;
	//Bulirsch Stoer:
	typedef boost::numeric::odeint::bulirsch_stoer
        < solution_type, value_type, solution_type, time_type >
        controlled_stepper_type_bulsto;

	//Instantiate all the solvers. Need to do this outside the switch scope,
	//otherwise we can't access them later when we want to do error location
	//refining.
	controlled_stepper_type_cash_karp54 controlled_stepper_cash_karp54 =
		boost::numeric::odeint::make_controlled
		< error_stepper_type_cash_karp54 >(AbsTol, RelTol);
	controlled_stepper_type_dorpri45 controlled_stepper_dorpri45 =
		boost::numeric::odeint::make_controlled
		< error_stepper_type_dorpri45 >(AbsTol, RelTol);
	controlled_stepper_type_fehlberg78 controlled_stepper_fehlberg78 =
		boost::numeric::odeint::make_controlled
		< error_stepper_type_fehlberg78 >(AbsTol, RelTol);
	controlled_stepper_type_bulsto controlled_stepper_bulsto(AbsTol, RelTol);

	//Actual integration step. Loop so that we get all the requested points in
	//tspan, including the adaptive selected points.

	//Make sure we include the initial point:
	sol.T.push_back(tspan[0]);
	sol.Y.push_back(y0);


	for (int i = 0; i < nTimeSteps - 1; i++)
	{

		//Setup initial conditions for next run:
		//time_type t_init = tspan[i];
		time_type t_init = i == 0 ? tspan[0] : sol.T[sol.T.size() - 1];
		time_type t_end = tspan[i + 1];
		value_type initialStepMax = initStepMaxFrac*(t_end - t_init);
		y = (i == 0 ? y0 : sol.Y[sol.Y.size() - 1]);
		/*std::cout << "\nt = " << t_init;
		for(int j = 0;j < y.size();j++)
        {
            std::cout << "\ny[" << j << "] = " << y[j];
        }*/
        //Reset the events function so that we don't accidentally repeat an
        //event:
        if(i > 0)
        {
            events(sol.T[sol.Y.size() - 2],sol.Y[sol.Y.size() - 2]);
            record_data.Ylast = sol.Y[sol.Y.size() - 2];
            record_data.Tlast = sol.T[sol.Y.size() - 2];
            record_data.valueLast = events.values;
            record_data.isterminalLast = events.isterminals;
            record_data.directionLast = events.directions;
        }

		do
		{
			success = 1;
			try {
				switch (odeSolver)
				{
				case 1:
				{
					//Cash-Karp algorithm:
					boost::numeric::odeint::integrate_adaptive
					< controlled_stepper_type_cash_karp54,
                      odeCaller<solution_type, time_type, value_type>,
                      solution_type, time_type, oowe >
                        (controlled_stepper_cash_karp54, oderhs_caller, y,
                         t_init, t_end,initialStepMax, record_data);
					break;
				}
				case 2:
				{
					//Dormand-Prince algorithm:
					boost::numeric::odeint::integrate_adaptive
					< controlled_stepper_type_dorpri45,
                      odeCaller<solution_type, time_type, value_type>,
                      solution_type, time_type, oowe >
						(controlled_stepper_dorpri45, oderhs_caller, y,
                         t_init, t_end,initialStepMax, record_data);
					break;
				}
				case 3:
				{
					//Runge-Kutta Fehlberg78 algorithm:
					//std::cout << "\nUsing RK Fehlberg78";
					/*std::cout << "1: Fine up to here.";
					std::cout << "\nt_init = " << t_init;
					std::cout << "\nt_end = " << t_end;
					for(int i = 0;i < y.size();i++)
                    {
                        std::cout << "\ny[" << i << "] = " << y[i];
                    }
                    std::cout << "\n";*/
					boost::numeric::odeint::integrate_adaptive
					< controlled_stepper_type_fehlberg78,
                      odeCaller<solution_type, time_type, value_type>,
                      solution_type, time_type, oowe >
						(controlled_stepper_fehlberg78, oderhs_caller, y,
                         t_init, t_end,initialStepMax, record_data);
					break;
				}
				case 4:
				{
					//Bulirsch-Stoer algorithm:
					boost::numeric::odeint::integrate_adaptive
					< controlled_stepper_type_bulsto,
                      odeCaller<solution_type, time_type, value_type>,
                      solution_type, time_type, oowe >
                        (controlled_stepper_bulsto, oderhs_caller, y, t_init,
                         t_end, initialStepMax, record_data);
					break;
				}
				default:
				{
					//Cash-Karp algorithm:
					boost::numeric::odeint::integrate_adaptive
					< controlled_stepper_type_cash_karp54,
                      odeCaller<solution_type, time_type, value_type>,
                      solution_type, time_type, oowe >
						(controlled_stepper_cash_karp54, oderhs_caller, y,
                         t_init, t_end, initialStepMax, record_data);
					break;
				}
				}
			}
			catch (boost::numeric::odeint::step_adjustment_error sae)
			{
				//Took too many steps while adjusting the step-size.
				success = -3;
				break;
			}
			catch (boost::numeric::odeint::no_progress_error npe)
			{
				//Too many iterations during integration.
				std::cout << "Error occurred during integration: \n"
                          << npe.what() << std::endl;
				success = -4;
				break;
			}
			catch (int e)
			{
				switch (e)
				{
				case 0:
					success = -5;
					std::cout << "Attempted division by zero error occurred "
                              << "during integration.\n";
					break;
				case 1:
					success = -6;
					std::cout << "NaN value encountered during integration.\n";
					break;
				case 2:
					success = -7;
					break;
				default:
					success = 0;
					std::cout << "Unknown integer exception occurred during"
                              << " integration.\n";
				}
				break;
			}
			catch (char* c)
			{
				success = 0;
				std::cout << "Unexpected integration failure. Reason: \n"
                          << c << std::endl;
			}
			catch (std::string c)
			{
				success = 0;
				std::cout << "Unexpected integration failure. Reason: \n"
                          << c << std::endl;
			}
			catch (...)
			{
				//Failed, but not sure why.
				std::cout << "An unknown error occurred. Integration failed.\n";
				success = 0;
				break;
			}
            if (record_data.terminate)
            {
                break;
            }
		}
		while (success == -7);
		//Terminate the loop if the output observer finds that we hit
		//a terminal event.

		//Indicate that the data now contained in is no lognger provisional
		//and can be kept.
		sol.nProvisional = int(sol.T.size());
		sol.nProvisionalEvents = int(sol.TE.size());
	}
	//Now refine the events, assuming we passed through the integration phase
	//unscathed:
	//First switch off event termination:
	record_data.use_termination = false;
	//Iterate through events and refine their locations:
	//Store results in sol.TE, sol.YE.
	typedef event_refiner< time_type, solution_type, value_type,
                           controlled_stepper_type_cash_karp54 > erCashKarp54;
	typedef event_refiner< time_type, solution_type, value_type,
                           controlled_stepper_type_dorpri45 > erDorPri45;
	typedef event_refiner< time_type, solution_type, value_type,
                           controlled_stepper_type_fehlberg78 > erFehlberg78;
	typedef event_refiner< time_type, solution_type, value_type,
                           controlled_stepper_type_bulsto > erBulSto54;
    std::cout << "\nSolution Complete with " << record_data.nSteps << " steps."
              <<" Refine events.";
	if (success == 1 && nEventsFound > 0)
	{
		time_type lower;
		time_type upper;
		solution_type solution;
		for (int i = 0; i < nEventsFound; i++)
		{
		    std::cout << "\nRefine event " << i;
			try {
				switch (odeSolver)
				{
				case 1:
				{
					//Define an event refiner to handle this event:
					lower = TE_before[i];
					upper = TE_after[i];
					solution = YE_before[i];
					erCashKarp54 refiner
                        (sol.IE[0], lower, solution, odefun, events,
                         controlled_stepper_cash_karp54);
					refiner.event = sol.IE[i];
					refiner.T_init = lower;
					refiner.Y_init = solution;
					//Refine the event's location:
					/*sol.TE.push_back(bisection< time_type, erCashKarp54,
                                                value_type >
                        (lower, upper, refiner, AbsTol));*/
                    if(refining_events)
                    {
                        sol.TE.push_back(find_zero_brent(refiner,lower,upper));
                        sol.YE.push_back(refiner.Y);
                    }
                    else
                    {
                        //Record the solution immediately following the event:
                        sol.TE.push_back(TE_after[i]);
                        sol.YE.push_back(YE_after[i]);
                    }
					break;
				}
				case 2:
				{
					//Define an event refiner to handle this event:
					lower = TE_before[i];
					upper = TE_after[i];
					solution = YE_before[i];
					erDorPri45 refiner(sol.IE[0], lower, solution, odefun,
                                       events, controlled_stepper_dorpri45);
					refiner.event = sol.IE[i];
					refiner.T_init = lower;
					refiner.Y_init = solution;
					//Refine the event's location:
					/*sol.TE.push_back(bisection< time_type, erDorPri45,
                                                value_type >
                        (lower, upper, refiner, AbsTol));*/
                    if(refining_events)
                    {
                        sol.TE.push_back(find_zero_brent(refiner,lower,upper));
                        sol.YE.push_back(refiner.Y);
                    }
                    else
                    {
                        //Record the solution immediately following the event:
                        sol.TE.push_back(TE_after[i]);
                        sol.YE.push_back(YE_after[i]);
                    }
					break;
				}
				case 3:
				{
					//Define an event refiner to handle this event:
					try
					{
					lower = TE_before[i];
					upper = TE_after[i];
					solution = YE_before[i];
					erFehlberg78 refiner(sol.IE[0], lower, solution, odefun,
                                         events, controlled_stepper_fehlberg78);
					refiner.event = sol.IE[i];
					refiner.T_init = lower;
					refiner.Y_init = solution;
					//Refine the event's location:
					/*sol.TE.push_back(bisection< time_type, erFehlberg78,
                                                value_type >
                        (lower, upper, refiner, AbsTol));*/

                    if(refining_events)
                    {
                        sol.TE.push_back(find_zero_brent(refiner,lower,upper));
                        sol.YE.push_back(refiner.Y);
                    }
                    else
                    {
                        //Record the solution immediately following the event:
                        sol.TE.push_back(TE_after[i]);
                        sol.YE.push_back(YE_after[i]);
                    }
					}
					catch(...)
					{
					    std::cout << "\nWarning - error occurred for event "
                                  << i << ", IE = " <<  sol.IE[i]
                                  << " during event refining. Event maybe "
                                  << "inaccurate.";
                        continue;
					}
					break;
				}
				case 4:
				{
					//Define an event refiner to handle this event:
					lower = TE_before[i];
					upper = TE_after[i];
					solution = YE_before[i];
					erBulSto54 refiner(sol.IE[0], lower, solution, odefun,
                                       events, controlled_stepper_bulsto);
					refiner.event = sol.IE[i];
					refiner.T_init = lower;
					refiner.Y_init = solution;
					//Refine the event's location:
					/*sol.TE.push_back(bisection< time_type, erBulSto54,
                                                value_type >
                        (lower, upper, refiner, AbsTol));*/
                    if(refining_events)
                    {
                        sol.TE.push_back(find_zero_brent(refiner,lower,upper));
                        sol.YE.push_back(refiner.Y);
                    }
                    else
                    {
                        //Record the solution immediately following the event:
                        sol.TE.push_back(TE_after[i]);
                        sol.YE.push_back(YE_after[i]);
                    }
					break;
				}
				default:
				{
					//Define an event refiner to handle this event:
					lower = TE_before[i];
					upper = TE_after[i];
					solution = YE_before[i];
					erCashKarp54 refiner
                        (sol.IE[0], lower, solution, odefun,events,
                         controlled_stepper_cash_karp54);
					refiner.event = sol.IE[i];
					refiner.T_init = lower;
					refiner.Y_init = solution;
					//Refine the event's location:
					/*sol.TE.push_back(bisection< time_type, erCashKarp54,
                                                value_type >
                        (lower, upper, refiner, AbsTol));*/
                    if(refining_events)
                    {
                        sol.TE.push_back(find_zero_brent(refiner,lower,upper));
                        sol.YE.push_back(refiner.Y);
                    }
                    else
                    {
                        //Record the solution immediately following the event:
                        sol.TE.push_back(TE_after[i]);
                        sol.YE.push_back(YE_after[i]);
                    }
					break;
				}
				}
				//Update the event locations in the solution grid:
				//sol.Y[event_log.EventLocations[i] - 1] = sol.YE[i];
				//sol.T[event_log.EventLocations[i] - 1] = sol.TE[i];
				//sol.eventLocations.push_back(event_log.EventLocations[i] - 1);
			}
			catch (boost::numeric::odeint::step_adjustment_error sae)
			{
				//Took too many steps while adjusting the step-size.
				std::cout << "Error occured during event refining: \n"
                          << sae.what() << std::endl;
				success = -3;
				break;
			}
			catch (boost::numeric::odeint::no_progress_error npe)
			{
				//Too many iterations during integration.
				std::cout << "Error occured during event refining: \n"
                          << npe.what() << std::endl;
				success = -4;
				break;
			}
			catch (int e)
			{
				switch (e)
				{
				case 0:
					success = -5;
					std::cout << "Attempted division by zero error occured "
					          << "during event refining.\n";
					break;
				case 1:
					success = -6;
					std::cout << "NaN value encountered during event "
                              << "refining.\n";
					break;
				default:
					success = 0;
					std::cout << "Unknown integer exception occured "
                              << "during event refining.\n";
				}
				break;
			}
			catch (...)
			{
				//Failed, but not sure why.
				std::cout << "An unknown error occured during event refining."
                          << " Integration failed.\n";
				success = 0;
				break;
			}
		}
	}
	return success;
}
//------------------------------------------------------------------------------
#endif //ODE_SOLVER_GENERIC_H
