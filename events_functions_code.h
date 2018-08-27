#ifndef EVENTSFUNCTIONS_CODE
#define EVENTSFUNCTIONS_CODE

//#include "multi_precision_definitions.h"
//#include "event_function_definition.h"

//Core solver definitions:
#include "events_functions.h"

//Needed libraries
#include <boost/math/constants/constants.hpp>


//==============================================================================
#ifdef DE_SITTER_INSTANTONS
template< class solution_type, class time_type, class value_type >
void testEvents1< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type value[] = { y[0] - value_type(100.0),y[0] - value_type(1e-5),y[1]
                         - value_type(10.0) };
	bool isterminal[] = { false,false,true };
	int direction[] = { 0,0,0 };

	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void testEvents2< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type value[] ={y[0] +value_type(10.0),y[0],y[1] + value_type(1000.0)};
	bool isterminal[] = { false,false,false };
	int direction[] = { 0,0,0 };

	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void eventsMethod0< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type value[] = {y[3],y[1],y[2],y0upperBound-y[0],y[0]-y0lowerBound};
	bool isterminal[] = { true,true,true,true,true};
	int direction[] = { 0,0,-1,0,0};
	//std::cout << "\nt = " << t;
	//std::cout << "\ny[0] = " << y[0];
	/*for(int i = 0;i < y.size();i++)
    {
        std::cout << "\ny[" << i << "] = " << y[i];
    }*/

	//Record events:
	this->record_event(value, isterminal, direction);
	/*values.assign(value,value + this->N);
    isterminals.assign(isterminal,isterminal + this->N);
    directions.assign(direction,direction + this->N);*/
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void events_delta_a< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type arg = this->H0*value_type(t);
	value_type a = sin(arg)/this->H0 + y[2];
	value_type ap = cos(arg) + y[3];
	value_type value[] = {ap,y[1],a,y0upperBound-y[0],y[0]-y0lowerBound};
	bool isterminal[] = { true,true,true,true,true};
	int direction[] = { 0,0,-1,0,0};
	//std::cout << "\nt = " << t;
	//std::cout << "\ny[0] = " << y[0];
	/*for(int i = 0;i < y.size();i++)
    {
        std::cout << "\ny[" << i << "] = " << y[i];
    }*/

	//Record events:
	this->record_event(value, isterminal, direction);
	/*values.assign(value,value + this->N);
    isterminals.assign(isterminal,isterminal + this->N);
    directions.assign(direction,direction + this->N);*/
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void events_dS_fixed< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type arg = this->H0*value_type(t);
	value_type a = sin(arg)/this->H0;
	value_type ap = cos(arg);
	value_type value[] = {ap,y[1],a,y0upperBound-y[0],y[0]-y0lowerBound};
	bool isterminal[] = { true,true,true,true,true};
	int direction[] = { 0,0,-1,0,0};
	//std::cout << "\nt = " << t;
	//std::cout << "\ny[0] = " << y[0];
	/*for(int i = 0;i < y.size();i++)
    {
        std::cout << "\ny[" << i << "] = " << y[i];
    }*/

	//Record events:
	this->record_event(value, isterminal, direction);
	/*values.assign(value,value + this->N);
    isterminals.assign(isterminal,isterminal + this->N);
    directions.assign(direction,direction + this->N);*/
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void eventsSwitch< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type r = (abs(y[1]*y[1]) + abs(V(y[0])))/this->V0;
	value_type value[] = {y[3],y[1],y[2],y0upperBound-y[0],
                          y[0]-y0lowerBound,r - this->thresh};
	bool isterminal[] = { true,true,true,true,true ,true};
	int direction[] = { 0,0,-1,0,0 ,0};
	//std::cout << "\nt = " << t;
	//std::cout << "\ny[0] = " << y[0];
	/*for(int i = 0;i < y.size();i++)
    {
        std::cout << "\ny[" << i << "] = " << y[i];
    }*/

	//Record events:
	this->record_event(value, isterminal, direction);
	/*values.assign(value,value + this->N);
    isterminals.assign(isterminal,isterminal + this->N);
    directions.assign(direction,direction + this->N);*/
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void eventsSwitch_analytic< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type r = (abs(y[1]*y[1]) + abs(V(y[0])))/this->V0;
	value_type arg = this->H0*value_type(t) + phi;
	value_type value[] = {cos(arg),y[1],sin(arg)/this->H0,y0upperBound-y[0],
                          y[0]-y0lowerBound,r - this->thresh};
	bool isterminal[] = { true,true,true,true,true,true};
	int direction[] = { 0,0,-1,0,0 ,0};
	//std::cout << "\nt = " << t;
	//std::cout << "\ny[0] = " << y[0];
	/*for(int i = 0;i < y.size();i++)
    {
        std::cout << "\ny[" << i << "] = " << y[i];
    }*/

	//Record events:
	this->record_event(value, isterminal, direction);
	/*values.assign(value,value + this->N);
    isterminals.assign(isterminal,isterminal + this->N);
    directions.assign(direction,direction + this->N);*/
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void eventsGravFriction< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	const value_type _1p0 = value_type(1.0);
	const value_type epsilon = std::numeric_limits<value_type>::epsilon();
	value_type value[] = {y[3],y[1],y[2],y0upperBound-y[0],
                            y[0]-y0lowerBound};
	bool isterminal[] = { false,true,true,true,true };
	int direction[] = { 0,1,-1,0,0 };

	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void eventsMethod1< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type value[] = {y[3],y[1],y[2],y0upperBound-y[0],y[0]-y0lowerBound};
	bool isterminal[] = { false,true,true,true,true };
	if (NRstage)
	{
		//NR stage wants to terminate when we find a' = 0,
		//but we don't want to do this for overshoot/undershoot
		//stage.
		isterminal[0] = true;
		//Disable other termination reasons:
		isterminal[1] = false;
		isterminal[2] = false;
		isterminal[3] = false;
		isterminal[4] = false;
	}
	int direction[] = { 0,0,-1,0,0 };

	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void eventsMethod2< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type value[] = {t - time_type(0.5),y[1],y[3],y0upperBound - y[0],
                          y[0] - y0lowerBound };
	bool isterminal[] = { false,true,false,true,true };
	if (NRstage)
	{
		//NR stage wants to terminate when we find a' = 0,
		//but we don't want to do this for overshoot/undershoot
		//stage.
		isterminal[2] = true;
		//Disable other termination reasons:
		isterminal[1] = false;
		isterminal[3] = false;
		isterminal[4] = false;
	}
	int direction[] = { 0,0,0,0,0 };

	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void eventsMethod3< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type a0p = value_type(PI*cos(PI*t));
	value_type value[] = {t - time_type(0.5),y[1],a0p + y[3],
                          y0upperBound - y[0],y[0] - y0lowerBound };
	bool isterminal[] = { false,true,false,true,true };
	if (NRstage)
	{
		//NR stage wants to terminate when we find a' = 0,
		//but we don't want to do this for overshoot/undershoot
		//stage.
		isterminal[2] = true;
		//Disable other termination reasons:
		isterminal[1] = false;
		isterminal[3] = false;
		isterminal[4] = false;
	}
	int direction[] = { 0,0,0,0,0 };

	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void eventsMethod4< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type value[] = { t - PI / (time_type(2.0)),y[1],y0upperBound - y[0],
                          y[0] - y0lowerBound };
	bool isterminal[] = { false,true,true,true };
	if (NRstage)
	{
		//NR stage wants to terminate when we find a' = 0,
		//but we don't want to do this for overshoot/undershoot
		//stage.
		isterminal[0] = true;
		//Disable other termination reasons:
		isterminal[1] = false;
		isterminal[2] = false;
		isterminal[3] = false;
	}
	int direction[] = { 0,0,0,0 };

	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
#endif
//==============================================================================


//==============================================================================
#ifdef ANTI_DE_SITTER_FLAT_INSTANTONS
//Events functions for anti-de-Sitter or flat space true vacuum:
template< class solution_type, class time_type, class value_type >
void eventsAdSFlat< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
    //Define events:
	value_type value[] = { y[0] - false_vacuum,y[1],y[2],
                          y[0] - value_to_check,y0upperBound - y[0],
                          y[0] - y0lowerBound };
	bool isterminal[] = { true,true,false,false,true,true };
	int direction[] = { 0,0,0,0,0,0 };

	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
//Events functions for flat space fixed background
template< class solution_type, class time_type, class value_type >
void eventsFlat< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{


	//Define events:
	value_type value[] = { y[0] - false_vacuum,y[1],y[0] - value_to_check,
                          y0upperBound - y[0],y[0] - y0lowerBound };
	bool isterminal[] = { true,(true_vacuum - barrier)*(y[0] - barrier) < 0,
                         false,true,true };
        //Choose this strange termination condition because undershoots can
        //only occur on the opposite side of the barrier to the tv.
        //Any other zeros are thus spurious and may be due to truncation
        //errors. We continue integrating until we are sure
	int direction[] = { 0,0,0,0,0 };
	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
//Events functions for AdS fixed background:
template< class solution_type, class time_type, class value_type >
void eventsAdS< solution_type, time_type, value_type >::operator()
(const time_type t, const solution_type& y)
{
	//Define events:
	value_type value[] = {y[0] - false_vacuum,y[1],y[0] - value_to_check,
                          y0upperBound - y[0],y[0] - y0lowerBound };
	bool isterminal[] = { true,true,false,true,true };
	int direction[] = { 0,0,0,0,0 };

	//Record events:
	//this->record_event(value, isterminal, direction);
	this->values.assign(value,value + this->N);
    this->isterminals.assign(isterminal,isterminal + this->N);
    this->directions.assign(direction,direction + this->N);
}
//------------------------------------------------------------------------------
#endif
//==============================================================================
#endif
