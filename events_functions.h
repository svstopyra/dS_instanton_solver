#ifndef EVENTS_FUNCTIONS_H
#define EVENTS_FUNCTIONS_H
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

//------------------------------------------------------------------------------
//==============================================================================


//Header file to store events functions. Contains no code
//(see events_functions.cpp)

//#include "event_detection.h"
#include "ode_rhs.h"
//Definition of a potential:
#include "potentials.h"

//Needed libraries
#include <boost/math/constants/constants.hpp>
//==============================================================================
//Definition of the events functions
template< class solution_type , class time_type, class value_type >
class events_function
{
protected:
    const int N;
public:
    std::vector< value_type > values;
	std::vector< bool > isterminals;
	std::vector< int > directions;
	const int nEvents;
	std::vector< bool > dynamicSwitch;//allows us to switch specific terminal
        //events on/off.
    bool terminate_on;//Whether events are terminal. Allows us to easily
        //switch off terminal
    //events if need be without changing the events function.
    value_type AbsTol;
    value_type RelTol;//Used to check whether we really crossed a zero or just
        //seemed to because of a fluctuation in a near-event solution.
        //Basically, only trigger the event IF the function modified by the
        //precision in either direction would still trigger it. If not
        //specified, then zero is assumed. Note that it is up to the derived
        //classes to implement this feature as they decide how events are
        //computed.

    //Constructors:
	events_function(int N);
	events_function(int N, value_type RELTOL,value_type ABSTOL);
	//events_function();
	//events_function(value_type RELTOL,value_type ABSTOL);
	//Member functions:
	DLL_EXPORT virtual void operator()
        (const time_type t,const solution_type& y) =0;
	DLL_EXPORT void record_event
        (value_type* value, bool* isterminal, int* direction);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
events_function< solution_type, time_type, value_type >::events_function(int n)
//events_function< solution_type, time_type, value_type >::events_function()
 : nEvents(n), N(n), values(n, value_type(0.0)),dynamicSwitch(n,true),
 isterminals(n,true), directions(n,0)
{

    std::vector<value_type> value (n, value_type(0.0));
    values = value;
    std::vector<bool> dynSwitch (n, true);
    dynamicSwitch = dynSwitch;
    std::vector<bool> isterminal (n, true);
    isterminals = isterminal;
    std::vector<int> direction (n,0);
    directions = direction;

	//Default constructor.
	terminate_on = true;
	RelTol = value_type(0.0);
	AbsTol = value_type(0.0);

	/*
	for(int k = 0;k < this->values.size();k++)
    {
        std::cout << "\nevents.values[" << k << "] = " << this->values[k];
    }
    */
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
events_function< solution_type, time_type, value_type >::events_function
    (int n, value_type RELTOL,value_type ABSTOL)
    //(value_type RELTOL,value_type ABSTOL)
     : nEvents(n), N(n), values(n, value_type(0.0)),dynamicSwitch(n,true),
    isterminals(n,true), directions(n,0)
{
    //Constructor with tolerances:

    std::vector<value_type> value (n, value_type(0.0));
    values = value;
    std::vector<bool> dynSwitch (n, true);
    dynamicSwitch = dynSwitch;
    std::vector<bool> isterminal (n, true);
    isterminals = isterminal;
    std::vector<int> direction (n,0);
    directions = direction;

	terminate_on = true;
	RelTol = RELTOL;
	AbsTol = ABSTOL;
	/*
	for(int k = 0;k < this->values.size();k++)
    {
        std::cout << "\nevents.values[" << k << "] = " << this->values[k];
    }
    */
}
//------------------------------------------------------------------------------
//==============================================================================

//==============================================================================
template<class time_type,class solution_type>
struct eventLogger
{
    std::vector< solution_type >& YE_before;
    std::vector< time_type >& TE_before;
    std::vector< solution_type >& YE_after;
    std::vector< time_type >& TE_after;
    std::vector< int >& IE;
    std::vector< int > EventLocations; //Position of the event in the list of
                                       //solution points.
    int& nEventsFound;
    //Constructor:
	DLL_EXPORT eventLogger(std::vector< solution_type >& ye_before,
                           std::vector< solution_type >& ye_after,
                           std::vector< time_type >& te_before,
                           std::vector< time_type >& te_after,
                           std::vector< int >& ie,int& events_found);
};
//------------------------------------------------------------------------------
//Class for the output observer, including event detection:
template< class time_type, class solution_type, class value_type >
class output_observer_with_events
{
    public:
	//Addresses of vectors to store data:
    std::vector< solution_type >& Y;
    std::vector< time_type >& T;
    //Track the number of steps:
    int& nSteps;
    //Store last time and value. Useful for event detection:
    time_type Tlast;
    solution_type Ylast;
    //Events function data:
	//Pointer to events function class:
	events_function< solution_type, time_type, value_type >& events;
	std::vector< value_type > valueLast;
    std::vector< bool > isterminalLast;
    std::vector< int > directionLast;
    //Store events:
    eventLogger< time_type , solution_type >& event_log;
    //Switch for event termination. Set this to false manually to disable
    //event logging.
    bool use_termination;
    //refining_mode - indicates that we are not solving for the solution, we
    //so much as trying to find better locations of the events associated
    //to the solution. Off by default:
    //bool refining_mode;
    //Record when we terminated early:
    bool& terminate;
	//Constructor, including a pointer to the std::function object supplied by
	//the user
	DLL_EXPORT output_observer_with_events
        (int& steps_in,std::vector< solution_type >& Y_input,
         std::vector< time_type >& T_input,
         events_function< solution_type, time_type,value_type >& event_function,
         eventLogger< time_type , solution_type >& Event_Log,
         bool& terminate_early);

    //observer function. Returns whether execution should terminate early
    //(false by default)
	DLL_EXPORT bool operator()( const solution_type &x , const time_type t );
};
//------------------------------------------------------------------------------
template<class time_type,class solution_type>
class null_obs
{
public:
    bool operator()(solution_type& Y, const time_type t)
    {
        //Never triggers termination:
        return false;
    }
};
//------------------------------------------------------------------------------
//function to evaluate events function at a given time. Used for bisection
//search. Also remembers our last guess for more efficient integrating.
template< class time_type , class solution_type , class value_type ,
          class Stepper >
class event_refiner
{
public:
    //Index IE of event we wish to return with this function.
    int event;
    //Time close to the event at which the solution is known:
    time_type T;
    //Solution at T:
    solution_type& Y;
    time_type T_init;//Initial guess at the event's location.
    //ode we are integrating:
    solution_type Y_init;//Value of the solution at T_init
	odeRHS< solution_type, time_type, value_type >& ode;
    //events function to evaluate:
	events_function< solution_type, time_type, value_type >& events;
    //Stepper to use:
    Stepper stepper;
    null_obs<time_type,solution_type> null_obs_event;
    //constructor:
	DLL_EXPORT event_refiner(int nEvent,time_type t_lower,
                             solution_type& Y_lower,
                             odeRHS< solution_type, time_type,
                                     value_type >& ode_to_solve,
                             events_function< solution_type,
                                              time_type, value_type >& Events,
                             Stepper stepper_to_use);
	DLL_EXPORT value_type operator()(time_type t);
};
//------------------------------------------------------------------------------
template< class time_type , class function , class value_type >
time_type bisection(time_type& ,time_type& ,function& ,value_type );
//------------------------------------------------------------------------------
template<class time_type,class solution_type>
struct refined_events
{
	std::vector< time_type > TE;
	std::vector< solution_type > YE;
	std::vector< int > IE;
	//Constructor:
	DLL_EXPORT refined_events(std::vector<int> ie): IE(ie){}
};
//------------------------------------------------------------------------------
//Function to refine events
template<class time_type,class solution_type, class value_type, class Stepper >
DLL_EXPORT refined_events<time_type,solution_type> refineEvents
        (eventLogger< time_type , solution_type >& ,
        output_observer_with_events< time_type ,
        solution_type , value_type >& ,
        odeRHS< solution_type, time_type, value_type >&,
        value_type ,Stepper );
//------------------------------------------------------------------------------
//==============================================================================

//==============================================================================
#ifdef DE_SITTER_INSTANTONS
template< class solution_type , class time_type ,class value_type >
class testEvents1: public events_function< solution_type , time_type ,
                                           value_type>
{
	//const static int N = 3;//Number of individual events.
public:
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class testEvents2: public events_function< solution_type , time_type ,
                                           value_type>
{
	//const static int N = 3;//Number of individual events.
public:
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class eventsMethod0: public events_function< solution_type , time_type ,
                                             value_type>
{
private:
	//const static int N = 5;//Number of individual events.
public:
    value_type y0lowerBound;
	value_type y0upperBound;
    eventsMethod0(value_type y0low,value_type y0high);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};

//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsMethod0< solution_type, time_type, value_type >::eventsMethod0
(value_type y0low, value_type y0high)
 : events_function< solution_type, time_type, value_type >(5)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
}
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class events_delta_a: public events_function< solution_type , time_type ,
                                             value_type>
{
private:
	//const static int N = 5;//Number of individual events.
public:
    value_type y0lowerBound;
	value_type y0upperBound;
	value_type H0;
    events_delta_a(value_type y0low,value_type y0high,value_type h0);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
events_delta_a< solution_type, time_type, value_type >::events_delta_a
(value_type y0low, value_type y0high,value_type h0)
 : events_function< solution_type, time_type, value_type >(5)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	H0 = h0;
}
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class eventsSwitch: public events_function< solution_type , time_type ,
                                             value_type>
{
private:
	//const static int N = 5;//Number of individual events.
	value_type thresh;//Threshold to trigger switch
	potential<value_type>& V;
	value_type V0;
public:
    value_type y0lowerBound;
	value_type y0upperBound;
    eventsSwitch(value_type y0low, value_type y0high,potential<value_type>& W,
                 value_type THRESH,value_type W0);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsSwitch< solution_type, time_type, value_type >::eventsSwitch
(value_type y0low, value_type y0high,potential<value_type>& W,
 value_type THRESH,value_type W0)
 : events_function< solution_type, time_type, value_type >(6), V(W)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	thresh = THRESH;
	V0 = W0;
}
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class eventsSwitch_analytic: public events_function< solution_type , time_type ,
                                             value_type>
{
private:
	//const static int N = 5;//Number of individual events.
	value_type thresh;//Threshold to trigger switch
	value_type H0;
	value_type V0;
	value_type h;
	potential<value_type>& V;
public:
    value_type y0lowerBound;
	value_type y0upperBound;
	value_type phi;
    eventsSwitch_analytic(value_type y0low, value_type y0high,
                          potential<value_type>& W,value_type THRESH,
                          value_type W0,value_type H);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsSwitch_analytic< solution_type, time_type, value_type >
::eventsSwitch_analytic
(value_type y0low, value_type y0high,potential<value_type>& W,
 value_type THRESH,value_type W0,value_type H)
 : events_function< solution_type, time_type, value_type >(6), V(W)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	thresh = THRESH;
	V0 = W0;
	h = H;
	const value_type _3p0 = value_type(3.0);
	const value_type _0p0 = value_type(0.0);
	H0 = h*sqrt(V0/_3p0);
	phi = _0p0;
}
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class eventsGravFriction: public events_function< solution_type , time_type ,
                                             value_type>
{
private:
	//const static int N = 5;//Number of individual events.
public:
    value_type y0lowerBound;
	value_type y0upperBound;
	eventsGravFriction(value_type y0low,value_type y0high);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsGravFriction< solution_type, time_type, value_type >::eventsGravFriction
(value_type y0low, value_type y0high)
 : events_function< solution_type, time_type, value_type >(5)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
}
//------------------------------------------------------------------------------

template< class solution_type , class time_type ,class value_type >
class eventsMethod1: public events_function< solution_type , time_type ,
                                             value_type >
{
private:
    value_type y0lowerBound;
    value_type y0upperBound;
	//const static int N = 5;//Number of individual events.
public:
    bool NRstage;
	eventsMethod1(value_type y0low,value_type y0high,bool nrstage = false);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsMethod1< solution_type, time_type, value_type >::eventsMethod1
(value_type y0low, value_type y0high, bool nrstage)
 : events_function< solution_type, time_type, value_type >(5)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	NRstage = nrstage;
}
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class eventsMethod2: public events_function< solution_type , time_type ,
                                             value_type >
{
private:
    value_type y0lowerBound;
    value_type y0upperBound;
	//const static int N = 5;//Number of individual events.
public:
    bool NRstage;
	eventsMethod2(value_type y0low,value_type y0high,bool nrstage = false);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsMethod2< solution_type, time_type, value_type >::eventsMethod2
(value_type y0low, value_type y0high, bool nrstage)
 : events_function< solution_type, time_type, value_type >(5)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	NRstage = nrstage;
}
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class eventsMethod3: public events_function< solution_type , time_type ,
                                             value_type >
{
private:
    value_type y0lowerBound;
    value_type y0upperBound;
	//const static int N = 5;//Number of individual events.
public:
    bool NRstage;
	eventsMethod3(value_type y0low,value_type y0high,bool nrstage = false);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsMethod3< solution_type, time_type, value_type >::eventsMethod3
(value_type y0low, value_type y0high, bool nrstage)
 : events_function< solution_type, time_type, value_type >(5)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	NRstage = nrstage;
}
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class eventsMethod4: public events_function< solution_type , time_type ,
                                             value_type >
{
private:
    value_type y0lowerBound;
    value_type y0upperBound;
	//const static int N = 4;//Number of individual events.
public:
    bool NRstage;
	eventsMethod4(value_type y0low,value_type y0high,bool nrstage = false);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsMethod4< solution_type, time_type, value_type >::eventsMethod4
(value_type y0low, value_type y0high, bool nrstage)
 : events_function< solution_type, time_type, value_type >(4)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	NRstage = nrstage;
}
//------------------------------------------------------------------------------
#endif
//==============================================================================



//==============================================================================
#ifdef ANTI_DE_SITTER_FLAT_INSTANTONS
//Events function including backreaction, for AdS-Flat case:
template< class solution_type , class time_type ,class value_type >
class eventsAdSFlat: public events_function< solution_type , time_type ,
                                             value_type >
{
private:
    value_type y0lowerBound;
    value_type y0upperBound;
	value_type false_vacuum;
	std::ostream& outStream;
	//const static int N = 6;//Number of individual events.
public:
	//Check for when the solution passes through the following point:
	value_type value_to_check;
	eventsAdSFlat(value_type y0low,value_type y0high, value_type fv);
	eventsAdSFlat(value_type y0low,value_type y0high, value_type fv ,
                  std::ostream& os);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsAdSFlat< solution_type, time_type, value_type >::eventsAdSFlat
(value_type y0low, value_type y0high, value_type fv)
 : events_function< solution_type, time_type, value_type >(6) ,
   outStream(std::cout)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	false_vacuum = fv;
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsAdSFlat< solution_type, time_type, value_type >::eventsAdSFlat
(value_type y0low, value_type y0high, value_type fv , std::ostream& os)
 : events_function< solution_type, time_type, value_type >(6) , outStream(os)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	false_vacuum = fv;
}
//------------------------------------------------------------------------------
//Events function for flat space, fixed background:
template< class solution_type, class time_type, class value_type >
class eventsFlat : public events_function< solution_type, time_type,
                                           value_type >
{
private:
	value_type y0lowerBound;
	value_type y0upperBound;
	value_type false_vacuum;
	value_type true_vacuum;
	value_type barrier;
	std::ostream& outStream;
	//const static int N = 5;//Number of individual events.
public:
	//Check for when the solution passes through the following point:
	value_type value_to_check;
	eventsFlat(value_type y0low, value_type y0high, value_type fv,
               value_type BARRIER,value_type tv);
	eventsFlat(value_type y0low, value_type y0high, value_type fv,
               value_type BARRIER,value_type tv,std::ostream& os);
	eventsFlat(value_type y0low, value_type y0high, value_type fv,
               value_type BARRIER,value_type tv,value_type RELTOL,
               value_type ABSTOL);
	DLL_EXPORT void operator()(const time_type t, const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsFlat< solution_type, time_type, value_type >::eventsFlat
    (value_type y0low,value_type y0high, value_type fv,value_type BARRIER,
     value_type tv)
      : events_function< solution_type, time_type, value_type >(5),
        outStream(std::cout)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	false_vacuum = fv;
	barrier = BARRIER;
	true_vacuum = tv;
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsFlat< solution_type, time_type, value_type >::eventsFlat
    (value_type y0low, value_type y0high, value_type fv,value_type BARRIER,
     value_type tv , std::ostream& os)
      : events_function< solution_type, time_type, value_type >(5) ,
        outStream(os)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	false_vacuum = fv;
	barrier = BARRIER;
	true_vacuum = tv;
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
eventsFlat< solution_type, time_type, value_type >::eventsFlat
    (value_type y0low, value_type y0high, value_type fv,value_type BARRIER,
     value_type tv,value_type RELTOL,value_type ABSTOL)
      : events_function< solution_type, time_type, value_type >
      (5 , RELTOL,ABSTOL), outStream(std::cout)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	false_vacuum = fv;
	barrier = BARRIER;
	true_vacuum = tv;
}
//------------------------------------------------------------------------------
//Events function for AdS space, fixed background:
template< class solution_type, class time_type, class value_type >
class eventsAdS : public events_function< solution_type, time_type, value_type >
{
private:
	value_type y0lowerBound;
	value_type y0upperBound;
	value_type false_vacuum;
	//const static int N = 5;//Number of individual events.
public:
	//Check for when the solution passes through the following point:
	value_type value_to_check;
	eventsAdS(value_type y0low, value_type y0high, value_type fv);
	DLL_EXPORT void operator()(const time_type t, const solution_type& y);
};
template< class solution_type, class time_type, class value_type >
eventsAdS< solution_type, time_type, value_type >::eventsAdS
(value_type y0low, value_type y0high, value_type fv)
 : events_function< solution_type, time_type, value_type >(5)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	false_vacuum = fv;
}
//------------------------------------------------------------------------------
template< class solution_type , class time_type ,class value_type >
class events_dS_fixed: public events_function< solution_type , time_type ,
                                             value_type>
{
private:
	//const static int N = 5;//Number of individual events.
public:
    value_type y0lowerBound;
	value_type y0upperBound;
	value_type H0;
    events_dS_fixed(value_type y0low,value_type y0high,value_type h0);
	DLL_EXPORT void operator()(const time_type t,const solution_type& y);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
events_dS_fixed< solution_type, time_type, value_type >::events_dS_fixed
(value_type y0low, value_type y0high,value_type h0)
 : events_function< solution_type, time_type, value_type >(5)
{
	y0lowerBound = y0low;
	y0upperBound = y0high;
	H0 = h0;
}
//------------------------------------------------------------------------------
#endif
//==============================================================================

#endif // EVENTS_FUNCTIONS_H
