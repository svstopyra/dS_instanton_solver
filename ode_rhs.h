#ifndef ODERHS_HEADER
#define ODERHS_HEADER

//DLL stuff - determines whether the header is being used by the dll itself (in
//which case we are exporting) or
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



//Contains only headers for the ode functions, so that we can separate them from
// the source code and compile it separately.
//Note - the separate source code file, ode_rhs.cpp only works if some of the
//needed specialisations are instantiated. This header
//otherwise contains no code.

//Needs to know what a potential is:
#include "potentials.h"
//#include "events_functions.h"
//Definition of an ode_rhs class:
//#include "ode_rhs_class_definition.h"

//Needed libraries:
#include <boost/math/constants/constants.hpp>//for pi

//==============================================================================
//==============================================================================
//Definition of the class:
//Base class for all odes:
template< class solution_type , class time_type , class value_type >
class odeRHS
{
public:
    virtual DLL_EXPORT void operator()(const solution_type& y,solution_type&
                                       dydt,const time_type& t) =0;
};
//------------------------------------------------------------------------------
//Base class for odes which store their own linearisation.
template< class solution_type , class time_type , class value_type >
class odeRHSwithLinear : public odeRHS< solution_type, time_type, value_type >
{
public:
	virtual DLL_EXPORT void linearised(const solution_type& Y,
                                solution_type& dydt, const time_type& t) = 0;
};
//------------------------------------------------------------------------------
//Because we never know what ode we are solving at compile-time (it is decided
//at run time) we need a wrapper that we can supply to the ode template. This
//will call the ode function indirectly, and simply stores a pointer to the
//relevant ode function. This has a slight overhead in that we have to use two
//function calls every time, but it should be ok.
//NB - if we know which ode we are solving at compile time, we can simply supply
// that to the ode function and not bother with this class. Only use it if we
//want to decide at run time.
template< class solution_type , class time_type,class value_type >
class odeCaller
{
private:
    odeRHS< solution_type , time_type , value_type >* odePointer;
public:
    //Constructor:
	odeCaller(odeRHS< solution_type, time_type, value_type >* pointer);
    //Function to call the odeRHS function:
	void operator()(const solution_type& y, solution_type& dydt,
                    const time_type t);
};
//------------------------------------------------------------------------------
//Class to call the Jacobian of the ode, rather than the ode itself. Essentially
//a wrapper for the oderhs class
template< class solution_type, class time_type, class value_type >
class jacobian : public odeRHS< solution_type , time_type , value_type >
{
private:
	odeRHSwithLinear< solution_type, time_type, value_type >& odeToCall;
public:

	//Constructor:
	jacobian(odeRHSwithLinear< solution_type, time_type, value_type >&
             ODETOCALL);
	//Operator:
	DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
};
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
jacobian< solution_type, time_type, value_type >::jacobian(
        odeRHSwithLinear< solution_type, time_type, value_type >& ODETOCALL)
         : odeToCall(ODETOCALL) {}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
odeCaller< solution_type, time_type, value_type >::odeCaller(
        odeRHS< solution_type, time_type, value_type >* pointer)
{
	odePointer = pointer;
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void odeCaller< solution_type, time_type, value_type >::operator()(
    const solution_type& y, solution_type& dydt, const time_type t)
{
	(*odePointer)(y, dydt, t);
}
//==============================================================================



//==============================================================================
#ifdef DE_SITTER_INSTANTONS
//Test ode 1 (simple exponential)
template< class solution_type , class time_type , class value_type >
class odeTest1RHS : public odeRHS< solution_type , time_type , value_type >
{
public:
	DLL_EXPORT void operator()(const solution_type& y,solution_type& dydt,
                               const time_type& t);
};//odeTest1RHS
//------------------------------------------------------------------------------
//Test ode 2 (precision test: sqrt(t) - 1/t )
template< class solution_type , class time_type , class value_type >
class odeTest2RHS : public odeRHS< solution_type , time_type , value_type >
{
public:
	DLL_EXPORT void operator()(const solution_type& y,solution_type& dydt,
                               const time_type& t);
};//odeTest2RHS
//------------------------------------------------------------------------------
//Gravitational instantons, using method 1, including linearisation.
template< class solution_type , class time_type , class value_type >
class ode1 : public odeRHSwithLinear< solution_type , time_type , value_type >
{
    //Overshoot/undershoot but without the action - useful for speed
    //optimisation!
private:
	potential< value_type >& V;
    value_type xi;
    value_type h;
    value_type W0;
public:
	ode1(potential< value_type >& v,value_type XI,value_type H,value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
	DLL_EXPORT void linearised(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
};//ode1
//------------------------------------------------------------------------------
//Constructor:
template< class solution_type, class time_type, class value_type >
ode1< solution_type, time_type, value_type >::ode1
    (potential< value_type >& v, value_type XI, value_type H,value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//ode1
//------------------------------------------------------------------------------
//Method 1, with action included.
template< class solution_type , class time_type , class value_type >
class ode1action : public odeRHS< solution_type , time_type , value_type >
{
private:
	potential< value_type >& V;
    value_type xi;
    value_type h;
    value_type W0;
public:
	ode1action(potential< value_type >& v,value_type XI,value_type H,value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
};//ode1action
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
ode1action< solution_type, time_type, value_type >::ode1action
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
    : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//ode1action
//------------------------------------------------------------------------------
//Gravitational instantons, using method 2, including linearisation.
template< class solution_type , class time_type , class value_type >
class ode2 : public odeRHSwithLinear< solution_type , time_type , value_type >
{
private:
	potential< value_type >& V;
    value_type xi;
    value_type h;
    value_type W0;
public:
	ode2(potential< value_type >& v,value_type XI,value_type H,value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
	DLL_EXPORT void linearised(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
};//ode2
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
ode2< solution_type, time_type, value_type >::ode2
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//ode2
//------------------------------------------------------------------------------
//Method 2, with action included.
template< class solution_type , class time_type , class value_type >
class ode2action : public odeRHS< solution_type , time_type , value_type >
{
private:
	potential< value_type >& V;
    value_type xi;
    value_type h;
    value_type W0;
public:
	ode2action(potential< value_type >& v,value_type XI,value_type H,
               value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
};//ode2action
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
ode2action< solution_type, time_type, value_type >::ode2action
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//ode2action
//------------------------------------------------------------------------------
//Gravitational instantons, using method 3, including linearisation.
template< class solution_type , class time_type , class value_type >
class ode3 : public odeRHSwithLinear< solution_type , time_type , value_type >
{
private:
	potential< value_type >& V;
    value_type xi;
    value_type h;
    value_type W0;
public:
	ode3(potential< value_type >& v,value_type XI,value_type H,value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
	DLL_EXPORT void linearised(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
};//ode3
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
ode3< solution_type, time_type, value_type >::ode3
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//ode3
//------------------------------------------------------------------------------
//Method 3, with action included.
template< class solution_type , class time_type , class value_type >
class ode3action : public odeRHS< solution_type , time_type , value_type >
{
private:
	potential< value_type >& V;
    value_type xi;
    value_type h;
    value_type W0;
public:
	ode3action(potential< value_type >& v,value_type XI,value_type H,
               value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
};//ode3action
//------------------------------------------------------------------------------
//ode3action constructor:
template< class solution_type, class time_type, class value_type >
ode3action< solution_type, time_type, value_type >::ode3action
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//ode3action
//------------------------------------------------------------------------------
//Gravitational instantons, using method 4, including linearisation.
template< class solution_type , class time_type , class value_type >
class ode4 : public odeRHSwithLinear< solution_type , time_type , value_type >
{
private:
    potential< value_type >& V;
    value_type xi;
    value_type h;
    value_type W0;
public:
	ode4(potential< value_type >& v,value_type XI,value_type H,value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
	DLL_EXPORT void linearised(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
};//ode4
//------------------------------------------------------------------------------
//ode4 constructor:
template< class solution_type, class time_type, class value_type >
ode4< solution_type, time_type, value_type >::ode4
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//ode4
//------------------------------------------------------------------------------
//Method 4, with action included.
template< class solution_type , class time_type , class value_type >
class ode4action : public odeRHS< solution_type , time_type , value_type >
{
private:
	potential< value_type >& V;
    value_type xi;
    value_type h;
    value_type W0;
public:
	ode4action(potential< value_type >& v,value_type XI,value_type H,
               value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y,solution_type& dydt,
                               const time_type& t);
};//ode4action
//------------------------------------------------------------------------------
  //ode4action constructor:
template< class solution_type, class time_type, class value_type >
ode4action< solution_type, time_type, value_type >::ode4action
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//ode4action
//------------------------------------------------------------------------------
#endif


#ifdef ANTI_DE_SITTER_FLAT_INSTANTONS
//With backreacttion - AdSFlat:
template< class solution_type, class time_type, class value_type >
class odeAdSFlat : public odeRHSwithLinear< solution_type,time_type,value_type >
{
private:
	potential< value_type >& V;
	value_type xi;
	value_type h;
	value_type W0;
public:
	value_type yscale;
	odeAdSFlat(potential< value_type >& v, value_type XI, value_type H,
               value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
	DLL_EXPORT void linearised(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
};
//------------------------------------------------------------------------------
//Constructor:
template< class solution_type, class time_type, class value_type >
odeAdSFlat< solution_type, time_type, value_type >::odeAdSFlat
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//odeAdSFlat
//------------------------------------------------------------------------------
//With backreacttion - AdSFlat:
template< class solution_type, class time_type, class value_type >
class odeAdSFlatAction : public odeRHS< solution_type, time_type, value_type >
{
private:
	potential< value_type >& V;
	value_type xi;
	value_type h;
	value_type W0;
public:
	value_type yscale;
	odeAdSFlatAction(potential< value_type >& v, value_type XI, value_type H,
                      value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
};
//------------------------------------------------------------------------------
//Constructor:
template< class solution_type, class time_type, class value_type >
odeAdSFlatAction< solution_type, time_type, value_type >::odeAdSFlatAction
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
	/*
	if (W0 > 0)
	{
		throw "'W0' must be negative semi-definite for AdS-Flat instantons.\n";
	}*/
}//odeAdSFlatAction
//------------------------------------------------------------------------------
//With backreacttion - AdSFlat:
template< class solution_type, class time_type, class value_type >
class odeAdSFlatAction_delta_a
 : public odeRHS< solution_type, time_type, value_type >
{
private:
	potential< value_type >& V;
	value_type xi;
	value_type h;
	value_type W0;
public:
	value_type yscale;
	const value_type _1p0;
	odeAdSFlatAction_delta_a(potential< value_type >& v, value_type XI,
                             value_type H,value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
};
//------------------------------------------------------------------------------
//Constructor:
template< class solution_type, class time_type, class value_type >
odeAdSFlatAction_delta_a< solution_type, time_type, value_type >
::odeAdSFlatAction_delta_a
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v), _1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
	/*
	if (W0 > 0)
	{
		throw "'W0' must be negative semi-definite for AdS-Flat instantons.\n";
	}*/
}//odeAdSFlatAction
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
class odeAdSFlat_delta_a
 : public odeRHS< solution_type,time_type,value_type >
{
private:
	potential< value_type >& V;
	value_type xi;
	value_type h;
	value_type W0;
public:
	value_type yscale;
	const value_type _1p0;

	odeAdSFlat_delta_a(potential< value_type >& v, value_type XI, value_type H,
               value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
};
//------------------------------------------------------------------------------
//Constructor:
template< class solution_type, class time_type, class value_type >
odeAdSFlat_delta_a< solution_type, time_type, value_type >::odeAdSFlat_delta_a
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v), _1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
}//odeAdSFlat
//------------------------------------------------------------------------------
//Fixed background, for flat space.
template< class solution_type, class time_type, class value_type >
class odeFlat : public odeRHSwithLinear< solution_type, time_type, value_type >
{
private:
	potential< value_type >& V;
	value_type xi;
	value_type h;
	value_type W0;
public:
	odeFlat(potential< value_type >& v, value_type XI, value_type H,
            value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
	DLL_EXPORT void linearised(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
};//ode4action
//------------------------------------------------------------------------------
//ode4action constructor:
template< class solution_type, class time_type, class value_type >
odeFlat< solution_type, time_type, value_type >::odeFlat
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//odeFlat
//------------------------------------------------------------------------------
//Fixed background, for flat space, with action:
template< class solution_type, class time_type, class value_type >
class odeFlatAction : public odeRHS< solution_type, time_type, value_type >
{
private:
	potential< value_type >& V;
	value_type xi;
	value_type h;
	value_type W0;
public:
	odeFlatAction(potential< value_type >& v, value_type XI, value_type H,
                  value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
};//ode4action
//------------------------------------------------------------------------------
//ode4action constructor:
template< class solution_type, class time_type, class value_type >
odeFlatAction< solution_type, time_type, value_type >::odeFlatAction
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//odeFlatAction
//------------------------------------------------------------------------------
//Fixed background, for AdS space.
template< class solution_type, class time_type, class value_type >
class odeAdS : public odeRHSwithLinear< solution_type, time_type, value_type >
{
private:
	potential< value_type >& V;
	value_type xi;
	value_type h;
	value_type W0;
public:
	odeAdS(potential< value_type >& v, value_type XI, value_type H,
           value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
	DLL_EXPORT void linearised(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
};//ode4action
//------------------------------------------------------------------------------
//ode4action constructor:
template< class solution_type, class time_type, class value_type >
odeAdS< solution_type, time_type, value_type >::odeAdS
    (potential< value_type >& v, value_type XI, value_type H, value_type w0)
     : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//odeAdS
//------------------------------------------------------------------------------
 //Fixed background, for AdS space, with action:
template< class solution_type, class time_type, class value_type >
class odeAdSAction : public odeRHS< solution_type, time_type, value_type >
{
private:
	potential< value_type >& V;
	value_type xi;
	value_type h;
	value_type W0;
public:
	odeAdSAction(potential< value_type >& v, value_type XI, value_type H,
                 value_type w0);
	DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                               const time_type& t);
};//ode4action
//------------------------------------------------------------------------------
//ode4action constructor:
template< class solution_type, class time_type, class value_type >
odeAdSAction< solution_type, time_type, value_type >::odeAdSAction
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v)
{
	xi = XI;
	h = H;
	W0 = w0;
}//odeFlatAction
//------------------------------------------------------------------------------
//For the gravitational case, but not tracking a and a' directly, instead
//tracking the ratio a'/a.
template< class solution_type , class time_type, class value_type >
class odeGravFriction : public odeRHS< solution_type, time_type, value_type >
{
private:
    potential< value_type>& V;
    value_type xi;
    value_type h;
    value_type W0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _1p0;
    const value_type _0p0;
public:
    odeGravFriction(potential< value_type >& v, value_type XI, value_type H,
                 value_type w0);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
};
//------------------------------------------------------------------------------
//odeGravFriction constructor:
template< class solution_type , class time_type, class value_type >
odeGravFriction< solution_type, time_type, value_type >::odeGravFriction
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_1p0(1.0),_0p0(0.0)
{
	xi = XI;
	h = H;
	W0 = w0;
}//odeGravFriction
//------------------------------------------------------------------------------
//For the gravitational case, but not tracking a and a' directly, instead
//tracking the ratio a'/a.
template< class solution_type , class time_type, class value_type >
class odeGravFrictionAction : public odeRHS< solution_type, time_type, value_type >
{
private:
    potential< value_type>& V;
    value_type xi;
    value_type h;
    value_type W0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _0p0;
    const value_type _1p0;
public:
    odeGravFrictionAction(potential< value_type >& v, value_type XI,
                          value_type H,value_type w0);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
};
//------------------------------------------------------------------------------
//odeGravFriction constructor:
template< class solution_type , class time_type, class value_type >
odeGravFrictionAction< solution_type, time_type,
                       value_type >::odeGravFrictionAction
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_0p0(0.0),_1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
}//odeGravFriction
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//For the gravitational case, but not tracking a and a' directly, instead
//tracking the ratio a'/a.
template< class solution_type , class time_type, class value_type >
class ode_dS : public odeRHS< solution_type, time_type, value_type >
{
private:
    potential< value_type>& V;
    value_type xi;
    value_type h;
    value_type W0;
    value_type H0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _0p0;
    const value_type _1p0;
public:
    ode_dS(potential< value_type >& v, value_type XI,
                          value_type H,value_type w0);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
};
//------------------------------------------------------------------------------
//ode_dS_delta_a constructor:
template< class solution_type , class time_type, class value_type >
ode_dS< solution_type, time_type, value_type >::ode_dS
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_0p0(0.0),_1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
	H0 = h*sqrt(W0/_3p0);
}//ode_dS
//------------------------------------------------------------------------------
//For the gravitational case, but not tracking a and a' directly, instead
//tracking the ratio a'/a.
template< class solution_type , class time_type, class value_type >
class ode_dS_delta_a : public odeRHS< solution_type, time_type, value_type >
{
private:
    potential< value_type>& V;
    value_type xi;
    value_type h;
    value_type W0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _0p0;
    const value_type _1p0;
    value_type H0;
public:
    ode_dS_delta_a(potential< value_type >& v, value_type XI,
                          value_type H,value_type w0);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
};
//------------------------------------------------------------------------------
//ode_dS_delta_a constructor:
template< class solution_type , class time_type, class value_type >
ode_dS_delta_a< solution_type, time_type, value_type >::ode_dS_delta_a
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_0p0(0.0),_1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
	H0 = sqrt(w0/value_type(3.0));
}//ode_dS
//------------------------------------------------------------------------------
//For the gravitational case, but not tracking a and a' directly, instead
//tracking the ratio a'/a.
template< class solution_type , class time_type, class value_type >
class ode_dS_action : public odeRHS< solution_type, time_type, value_type >
{
private:
    potential< value_type>& V;
    value_type xi;
    value_type h;
    value_type W0;
    value_type H0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _0p0;
    const value_type _1p0;
public:
    ode_dS_action(potential< value_type >& v, value_type XI,
                          value_type H,value_type w0);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
};
//------------------------------------------------------------------------------
//ode_dS constructor:
template< class solution_type , class time_type, class value_type >
ode_dS_action< solution_type, time_type, value_type >::ode_dS_action
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_0p0(0.0),_1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
	H0 = h*sqrt(W0/_3p0);
}//ode_dS
//------------------------------------------------------------------------------
//For the gravitational case, but not tracking a and a' directly, instead
//tracking the ratio a'/a.
template< class solution_type , class time_type, class value_type >
class ode_dS_action_delta_a
 : public odeRHS< solution_type, time_type, value_type >
{
private:
    potential< value_type>& V;
    value_type xi;
    value_type h;
    value_type W0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _0p0;
    const value_type _1p0;
    value_type H0;
public:
    ode_dS_action_delta_a(potential< value_type >& v, value_type XI,
                          value_type H,value_type w0);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
};
//------------------------------------------------------------------------------
//ode_dS constructor:
template< class solution_type , class time_type, class value_type >
ode_dS_action_delta_a< solution_type, time_type, value_type >
::ode_dS_action_delta_a
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_0p0(0.0),_1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
	H0 = sqrt(w0/value_type(3.0));
}//ode_dS
//------------------------------------------------------------------------------
//de Sitter case, with the option to switch to an anlytic form for a(t) when
//a certain threshold is passed.
template< class solution_type , class time_type, class value_type >
class ode_dS_switch : public odeRHS< solution_type, time_type, value_type >
{
private:
    potential< value_type>& V;
    value_type xi;//Non-minimal coupling
    value_type h;//ratio of scale to plank mass.
    value_type W0;//Cosmological constant energy density
    value_type H0;//Hubble rate
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _0p0;
    const value_type _1p0;
public:
    ode_dS_switch(potential< value_type >& v, value_type XI,
                          value_type H,value_type w0);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
    value_type phi;//Phase offset used for analytic form of a(t).
};
//------------------------------------------------------------------------------
//ode_dS constructor:
template< class solution_type , class time_type, class value_type >
ode_dS_switch< solution_type, time_type, value_type >::ode_dS_switch
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_0p0(0.0),_1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
	H0 = h*sqrt(W0/this->_3p0);//Hubble rate.
	phi = this->_0p0;//Phase offset is initially zero (needs to be set by
                                                       //matching conditions).
}//ode_dS_switch
//------------------------------------------------------------------------------
//Ode switch, including action:
template< class solution_type , class time_type, class value_type >
class ode_dS_action_switch : public odeRHS< solution_type, time_type,
                                            value_type >
{
private:
    potential< value_type>& V;
    value_type xi;
    value_type h;
    value_type W0;
    value_type H0;//Hubble rate
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _0p0;
    const value_type _1p0;
public:
    ode_dS_action_switch(potential< value_type >& v, value_type XI,
                          value_type H,value_type w0,value_type FV);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
    value_type phi;//Phase offset used for analytic form of a(t).
    value_type fv;
};
//------------------------------------------------------------------------------
//ode_dS_action_switch constructor:
template< class solution_type , class time_type, class value_type >
ode_dS_action_switch< solution_type, time_type,
                      value_type >::ode_dS_action_switch
(potential< value_type >& v, value_type XI, value_type H, value_type w0,
 value_type FV) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_0p0(0.0),_1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
	H0 = h*sqrt(W0/value_type(3.0));//Hubble rate.
	std::cout << "\nH0 = " << H0;
	phi = value_type(0.0);//Phase offset is initially zero (needs to be set by
                                                       //matching conditions).
    fv = FV;
}//ode_dS
//------------------------------------------------------------------------------
//Fixed de Sitter background, no action:
template< class solution_type , class time_type, class value_type >
class ode_dS_fixed
 : public odeRHS< solution_type, time_type, value_type >
{
private:
    potential< value_type>& V;
    value_type xi;
    value_type h;
    value_type W0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _0p0;
    const value_type _1p0;
    value_type H0;
public:
    ode_dS_fixed(potential< value_type >& v, value_type XI,
                          value_type H,value_type w0);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
};
//------------------------------------------------------------------------------
//ode_dS_fixed constructor:
template< class solution_type , class time_type, class value_type >
ode_dS_fixed< solution_type, time_type, value_type >
::ode_dS_fixed
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_0p0(0.0),_1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
	H0 = sqrt(w0/value_type(3.0));
}//ode_dS
//------------------------------------------------------------------------------
//Fixed de Sitter background
template< class solution_type , class time_type, class value_type >
class ode_dS_fixed_action
 : public odeRHS< solution_type, time_type, value_type >
{
private:
    potential< value_type>& V;
    value_type xi;
    value_type h;
    value_type W0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _0p0;
    const value_type _1p0;
    value_type H0;
public:
    ode_dS_fixed_action(potential< value_type >& v, value_type XI,
                          value_type H,value_type w0);
    DLL_EXPORT void operator()(const solution_type& Y, solution_type& dydt,
                   const time_type& t);
};
//------------------------------------------------------------------------------
//ode_dS_fixed constructor:
template< class solution_type , class time_type, class value_type >
ode_dS_fixed_action< solution_type, time_type, value_type >
::ode_dS_fixed_action
(potential< value_type >& v, value_type XI, value_type H, value_type w0) : V(v),
_4p0(4.0),_6p0(6.0),_2p0(2.0),_3p0(3.0),_0p0(0.0),_1p0(1.0)
{
	xi = XI;
	h = H;
	W0 = w0;
	H0 = sqrt(w0/value_type(3.0));
}//ode_dS
//------------------------------------------------------------------------------
#endif


#endif
