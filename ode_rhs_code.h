#ifndef ODE_RHS_CODE_H
#define ODE_RHS_CODE_H

//Core code definitions:
#include "ode_rhs.h"

//Needed libraries:
#include <boost/array.hpp>//For boost arrays (needed for ode_types)
#include <boost/math/special_functions.hpp>

//------------------------------------------------------------------------------
#ifdef DE_SITTER_INSTANTONS
//odeTest1RHS code:
template< class solution_type, class time_type, class value_type >
void odeTest1RHS< solution_type, time_type, value_type >::operator()
    (const solution_type& y, solution_type& dydt, const time_type& t)
{
	//Solve y'' - y' -2y =0, which is known to be y = Ae^2x + Be^-x.
	//y(0) = A + B, y'(0) = 2A - B so A = (y(0) + y'(0))/3,
	//B = (2*y(0) - y'(0))/3
	dydt[0] = y[1];
	dydt[1] = y[1] + 2 * y[0];
}//odeTest1RHS::operator()
//------------------------------------------------------------------------------
 //odeTest2RHS code:
template< class solution_type, class time_type, class value_type >
void odeTest2RHS< solution_type, time_type, value_type >::operator()
(const solution_type& y, solution_type& dydt, const time_type& t)
{
	if (t == time_type(0.0))
	{
		throw "odeTest2RHS attempted division by zero.";
	}
	dydt[0] = value_type(3.0) / (value_type(2.0)*(t*t))
		+ y[0] / (value_type(2.0)*t);
	dydt[1] = value_type(3.0) / (value_type(2.0)*(t*t))
		+ y[1] / (value_type(2.0)*t);
}//odeTest2RHS::operator()
//------------------------------------------------------------------------------
 //ode1 code:
 //Main ode1 function:
template< class solution_type, class time_type, class value_type >
void ode1< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (a == value_type(0.0) || (value_type(1.0) - h*h*xi*y*y) ==
        value_type(0.0) || value_type(1.0) - h*h*xi*(value_type(1.0) -
        value_type(6.0)*xi)*y*y == value_type(0.0))
	{
		throw "ode1 attempted division by zero.";
	}
	value_type ypp = -value_type(3.0)*(ap / a)*yp + (yp*yp)*
        ((xi*h*h*y*(value_type(1.0) - value_type(6.0)*xi)) /
         (value_type(1.0) - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ Wp*((value_type(1.0) - h*h*xi*y*y) / (value_type(1.0) -
        h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ h*h*(W + W0)*((value_type(4.0)*xi*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y));
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = ap;
	dydt[3] = -((h*h*a) / (value_type(3.0)*(value_type(1.0) - h*h*xi*y*y)))*
            (yp*yp + W0 + W	- value_type(3.0)*xi*(yp*yp + y*ypp +
            (ap / a)*y*yp));
}//ode1::operator()
//------------------------------------------------------------------------------
 //Linearised ode1 function:
template< class solution_type, class time_type, class value_type >
void ode1< solution_type, time_type, value_type >::linearised
    (const solution_type& Y, solution_type& dydt, const time_type& t)
{
	//Returns a linearised version of the ode. This only makes sense if we
	//are using an extended solution type where the first few variables are y,
	//and the latter variables the linearised fluctuations about y.

	//Call the ode to get the first few dydts (non-linear part, as
	//the linear part depends on these)
	(*this)(Y, dydt, t);

	//Define the variables:
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	//Non-linear part:
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type W = V(y);
	value_type Wp = V.d(y);
	value_type ypp = -value_type(3.0)*(ap / a)*yp + (yp*yp)*((xi*h*h*y
        *(value_type(1.0) - value_type(6.0)*xi)) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ Wp*((value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ h*h*(W + W0)*((value_type(4.0)*xi*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y));
	//Linear Fluctuations:
	value_type Dy = Y[4];
	value_type Dyp = Y[5];
	value_type Da = Y[6];
	value_type Dap = Y[7];
	value_type Wpp = V.d2(y);
	value_type Wy0 = W;
	//Return the derivative of the linearised fluctuations:
	dydt[4] = Dyp;
	dydt[5] = ((value_type(3.0)*(((value_type(1.0)) / (a*a))*(ap*(Da*yp))))
              + ((value_type(-3.0)*(((value_type(1.0)) / (a))*(Dap*yp)))
              + ((Dyp*((value_type(-3.0)*(((value_type(1.0)) / (a))*ap))
              + (value_type(2.0)*(h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*(y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*yp))))))))
              + (Dy*((value_type(8.0)*(h*h*h*h*((W0 + Wy0)*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*xi*(y*y*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
              - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))))))
              + ((value_type(2.0)*(h*h*(Wp*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*(y*((value_type(1.0)
              - (h*h*(xi*y*y)))*((value_type(1.0)) / ((value_type(1.0)
              - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
              - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y)))))))))))))
              + ((value_type(4.0)*(h*h*((W0 + Wy0)*(xi*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))))
              + ((value_type(2.0)*(h*h*(Wp*(xi*(y*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y)))))))))))
              + ((Wpp*((value_type(1.0) - (h*h*(xi*y*y)))*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))
              + ((value_type(2.0)*(h*h*h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*xi*(y*y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
              - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*yp*yp))))))
              + (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*yp*yp))))))))))))));
	dydt[6] = Dap;
	dydt[7] = ((Dap*(h*h*(xi*(y*(((value_type(1.0)) / ((value_type(1.0)
              - (h*h*(xi*y*y)))))*yp))))) + ((((value_type(-1.0))
              / (value_type(3.0)))*(a*(Dyp*(h*h*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))))*((value_type(2.0)*yp)
              + (value_type(-3.0)*(xi*((((value_type(1.0))
              / (a))*(ap*y)) + (value_type(2.0)*yp))))))))))
              + ((Dy*((((value_type(-1.0))
              / (value_type(3.0)))*(a*(h*h*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))))*(Wp
              + (value_type(-3.0)*(xi*((((value_type(1.0)) / (a))*(ap*yp))
              + ypp)))))))) + (((value_type(-2.0))
              / (value_type(3.0)))*(a*(h*h*h*h*(xi*(y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))*(value_type(1.0)
              - (h*h*(xi*y*y)))))*(W + (W0 + (yp*yp
              + (value_type(-3.0)*(xi*((((value_type(1.0))
              / (a))*(ap*(y*yp))) + (yp*yp + (y*ypp))))))))))))))))
              + (Da*((((value_type(-1.0))
              / (value_type(3.0)))*(h*h*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))))*(W + (W0 + (yp*yp
              + (value_type(-3.0)*(xi*((((value_type(1.0))
              / (a))*(ap*(y*yp))) + (yp*yp + (y*ypp)))))))))))
              - (((value_type(1.0)) / (a))*(ap*(h*h*(xi*(y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))))*yp)))))))))));
}//ode1::linearised
//------------------------------------------------------------------------------
 //ode1action with code:

 //Main ode1action:
template< class solution_type, class time_type, class value_type >
void ode1action< solution_type, time_type, value_type >::operator()
    (const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (a == value_type(0.0) || (value_type(1.0) - h*h*xi*y*y)
        == value_type(0.0) || value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y == value_type(0.0))
	{
		throw "ode1 attempted division by zero.";
	}
	value_type ypp = -value_type(3.0)*(ap / a)*yp
        + (yp*yp)*((xi*h*h*y*(value_type(1.0) - value_type(6.0)*xi))
        / (value_type(1.0) - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ Wp*((value_type(1.0) - h*h*xi*y*y)
        / (value_type(1.0) - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ h*h*(W + W0)*((value_type(4.0)*xi*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y));
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = ap;
	dydt[3] = -((h*h*a) / (value_type(3.0)*(value_type(1.0)
              - h*h*xi*y*y)))*(yp*yp + W0 + W
              - value_type(3.0)*xi*(yp*yp + y*ypp + (ap / a)*y*yp));

	//Action:
	value_type a0 = value_type(sin(t));
	value_type a0p = value_type(cos(t));
	dydt[4] = (value_type(2.0)*(PI*PI*((value_type(3.0)
              *(a0*a0*a0*((value_type(1.0)) / (h*h))))
              + (a*a*a*((((value_type(3.0)*(xi*((Wp*y) + yp*yp)))
              + (value_type(6.0)*(h*h*(xi*xi*(y*y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*((value_type(2.0)*(W
              + W0)) + ((((value_type(1.0)) / (value_type(2.0)))*yp*yp)
              + (value_type(-3.0)*(xi*((Wp*y) + yp*yp))))))))))) - W0) - W)))));
}//ode1action::operator()
//------------------------------------------------------------------------------
 //ode2 code:
 //Constructor ode2
 //Main ode2 function:
template< class solution_type, class time_type, class value_type >
void ode2< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type xmax = Y[4];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (a == value_type(0.0) || (value_type(1.0) - h*h*xi*y*y)
        == value_type(0.0)
		|| value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y == value_type(0.0))
	{
		throw "ode2 attempted division by zero.";
	}
	value_type ypp = -value_type(3.0)*(ap / a)*yp +
		(yp*yp)*((xi*h*h*(value_type(1.0) - value_type(6.0)*xi)*y)
         / (value_type(1.0) - h*h*xi*(value_type(1.0)
         - value_type(6.0)*xi)*y*y)) + (xmax*xmax)*Wp*(
			(value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
         - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
         + (xmax*xmax)*h*h*((value_type(4.0)*xi*y*(W + W0))
         / (value_type(1.0) - h*h*xi*(value_type(1.0)
         - value_type(6.0)*xi)*y*y));
	//Ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = ap;
	dydt[3] = -((h*h*a) / (value_type(3.0)*(1 - h*h*xi*y*y)))*(yp*yp
              + (xmax*xmax)*(W + W0)
              - value_type(3.0)*xi*(yp*yp + y*ypp + (ap / a)*y*yp));
	dydt[4] = 0;
}//ode2::operator()
//------------------------------------------------------------------------------
 //Linearised ode2 function:
template< class solution_type, class time_type, class value_type >
void ode2< solution_type, time_type, value_type >::linearised
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	//Returns a linearised version of the ode. This only makes sense
	//if we are using an
	//extended solution type where the first few variables are y, and the
	//latter variables the linearised fluctuations about y.

	//Call the ode to get the first few dydts (non-linear part, as
	//the linear part depends on these)
	(*this)(Y, dydt, t);

	//Define the variables:
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	//Non-linear part:
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type xmax = Y[4];
	value_type W = V(y);
	value_type Wp = V.d(y);
	value_type ypp = -value_type(3.0)*(ap / a)*yp
        + (yp*yp)*((xi*h*h*(value_type(1.0) - value_type(6.0)*xi)*y)
        / (value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y)) + (xmax*xmax)*Wp
        *((value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
        + (xmax*xmax)*h*h*((value_type(4.0)*xi*y*(W + W0))
        / (value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y));
	//Linear Fluctuations:
	value_type Dy = Y[5];
	value_type Dyp = Y[6];
	value_type Da = Y[7];
	value_type Dap = Y[8];
	value_type Dxmax = Y[9];
	value_type Wpp = V.d2(y);
	value_type Wy0 = W;

	//Return the derative of the linearised fluctuations:
	dydt[5] = Dyp;
	dydt[6] = ((Dxmax*((value_type(8.0)*(h*h*(Wy0*(xi*(xmax*(y
              *((value_type(1.0)) / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y)))))))))))) + (value_type(2.0)
              *(Wp*(xmax*((value_type(1.0) - (h*h*(xi*y*y)))*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))))))
              + ((value_type(3.0)*(((value_type(1.0)) / (a*a))*(ap*(Da*yp))))
              + ((value_type(-3.0)*(((value_type(1.0)) / (a))*(Dap*yp)))
              + ((Dyp*((value_type(-3.0)*(((value_type(1.0)) / (a))*ap))
              + (value_type(2.0)*(h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*(y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*yp))))))))
              + (Dy*((value_type(8.0)*(h*h*h*h*(Wy0*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*xi*(xmax*xmax*(y*y*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
              - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y)))))))))))))
              + ((value_type(2.0)*(h*h*(Wp*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*(xmax*xmax*(y*((value_type(1.0)
              - (h*h*(xi*y*y)))*((value_type(1.0)) / ((value_type(1.0)
              - (h*h*((value_type(1.0) + (value_type(-6.0)*xi))*(xi*y*y))))
              *(value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))))))))
              + ((value_type(4.0)*(h*h*(Wy0*(xi*(xmax*xmax*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y)))))))))))
              + ((value_type(2.0)*(h*h*(Wp*(xi*(xmax*xmax*(y*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))))))
              + ((Wpp*(xmax*xmax*((value_type(1.0)
              - (h*h*(xi*y*y)))*((value_type(1.0)) / ((value_type(1.0)
              - (h*h*((value_type(1.0) + (value_type(-6.0)*xi))*(xi*y*y)))))))))
              + ((value_type(2.0)*(h*h*h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*xi*(y*y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
              - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*yp*yp))))))
              + (h*h*((value_type(1.0) + (value_type(-6.0)*xi))
              *(xi*(((value_type(1.0)) / ((value_type(1.0)
              - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*yp*yp)))))))))))))));
	dydt[7] = Dap;
	dydt[8] = ((((value_type(-2.0)) / (value_type(3.0)))
              *(a*(Dxmax*(h*h*(W*(xmax*((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))))))))))
              + ((Dap*(h*h*(xi*(y*(((value_type(1.0)) / ((value_type(1.0)
              - (h*h*(xi*y*y)))))*yp))))) + ((((value_type(-1.0))
              / (value_type(3.0)))*(a*(Dyp*(h*h*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))))*((value_type(2.0)*yp)
              + (value_type(-3.0)*(xi*((((value_type(1.0)) / (a))*(ap*y))
              + (value_type(2.0)*yp)))))))))) + ((Dy*((((value_type(-1.0))
              / (value_type(3.0)))*(a*(h*h*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))))*((Wp*xmax*xmax)
              + (value_type(-3.0)*(xi*((((value_type(1.0)) / (a))*(ap*yp))
              + ypp)))))))) + (((value_type(-2.0)) / (value_type(3.0)))
              *(a*(h*h*h*h*(xi*(y*(((value_type(1.0)) / ((value_type(1.0)
              - (h*h*(xi*y*y)))*(value_type(1.0) - (h*h*(xi*y*y)))))
              *((W*xmax*xmax) + (yp*yp
              + (value_type(-3.0)*(xi*((((value_type(1.0))
              / (a))*(ap*(y*yp))) + (yp*yp + (y*ypp)))))))))))))))
              + (Da*((((value_type(-1.0)) / (value_type(3.0)))
              *(h*h*(((value_type(1.0)) / ((value_type(1.0)
              - (h*h*(xi*y*y)))))*((W*xmax*xmax) + (yp*yp + (value_type(-3.0)
              *(xi*((((value_type(1.0)) / (a))*(ap*(y*yp))) + (yp*yp
              + (y*ypp)))))))))) - (((value_type(1.0))
              / (a))*(ap*(h*h*(xi*(y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))))*yp))))))))))));
	dydt[9] = 0;
}//ode2::linearised


 //Constructor ode2action:



 //Main ode2action function:
template< class solution_type, class time_type, class value_type >
void ode2action< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type xmax = Y[4];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (a == value_type(0.0) || (value_type(1.0) - h*h*xi*y*y)
        == value_type(0.0)
		|| value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y == value_type(0.0))
	{
		throw "ode2 attempted division by zero.";
	}
	value_type ypp = -value_type(3.0)*(ap / a)*yp +
		(yp*yp)*((xi*h*h*(value_type(1.0) - value_type(6.0)*xi)*y)
        / (value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y)) + (xmax*xmax)*Wp*(
			(value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y)) +
		(xmax*xmax)*h*h*((value_type(4.0)*xi*y*(W + W0))
        / (value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y));
	//Ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = ap;
	dydt[3] = -((h*h*a) / (value_type(3.0)*(1 - h*h*xi*y*y)))*(yp*yp
              + (xmax*xmax)*(W + W0)
		- value_type(3.0)*xi*(yp*yp + y*ypp + (ap / a)*y*yp));
	dydt[4] = 0;
	//Action:
	value_type a0 = value_type(sin(PI*t));
	value_type a0p = value_type(PI*cos(PI*t));
	dydt[5] = (value_type(2.0)*(PI*PI*((value_type(3.0)*PI*a0*a0*a0
              / (h*h)) + (a*a*a*(xmax*(((value_type(3.0)*(xi*((Wp*y)
              + (((value_type(1.0)) / (xmax*xmax))*yp*yp))))
              + (value_type(6.0)*(h*h*(xi*xi*(y*y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*((value_type(2.0)*W)
              + ((((value_type(1.0)) / (value_type(2.0)))*(((value_type(1.0))
              / (xmax*xmax))*yp*yp)) + (value_type(-3.0)*(xi*((Wp*y)
              + (((value_type(1.0)) / (xmax*xmax))*yp*yp)))))))))))) - W))))));
}//ode2action::operator()
//------------------------------------------------------------------------------
 //ode3 functions:
 //ode3 Constructor:
 //Main ode3 function:
template< class solution_type, class time_type, class value_type >
void ode3< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type da = Y[2];
	value_type dap = Y[3];
	value_type dxmax = Y[4];
	value_type W = V(y);
	value_type Wp = V.d(y);
	value_type a0 = value_type(sin(PI*t));
	value_type a0p = value_type(PI*cos(PI*t));
	value_type a0pp = value_type(-PI*PI*sin(PI*t));
	if (a0 + da == value_type(0.0) || (1 - h*h*xi*y*y) == value_type(0.0)
		|| 1 - h*h*xi*(1 - 6 * xi)*y*y == value_type(0.0))
	{
		throw "ode3 attempted division by zero.";
	}
	value_type ypp = -3 * ((a0p + dap) / (a0 + da))*yp +
		(yp*yp)*((xi*h*h*(1 - 6 * xi)*y) / (1 - h*h*xi*(1 - 6 * xi)*y*y))
		+ ((PI + dxmax)*(PI + dxmax))*Wp*((1 - h*h*xi*y*y)
			/ (1 - h*h*xi*(1 - 6 * xi)*y*y))
		+ ((PI + dxmax)*(PI + dxmax))*h*h*((4 * xi*y*(W + W0))
			/ (1 - h*h*xi*(1 - 6 * xi)*y*y));
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = dap;
	dydt[3] = -a0pp - ((h*h*(a0 + da)) / (value_type(3.0)*(1
              - h*h*xi*y*y)))*(yp*yp + ((PI + dxmax)*(PI + dxmax))*(W + W0)
		- value_type(3.0)*xi*(yp*yp + y*ypp + ((a0p + dap) / (a0 + da))*y*yp));
	dydt[4] = 0;
}//ode3::operator()


 //Linearised ode3 function:
template< class solution_type, class time_type, class value_type >
void ode3< solution_type, time_type, value_type >::linearised
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	//Returns a linearised version of the ode. This only makes sense
	// if we are using an
	//extended solution type where the first few variables are y, and the
	//latter variables the linearised fluctuations about y.

	//Call the ode to get the first few dydts (non-linear part, as
	//the linear part depends on these)
	(*this)(Y, dydt, t);

	//Define the variables:
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	//Non-linear part:
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type da = Y[2];
	value_type dap = Y[3];
	value_type dxmax = Y[4];
	value_type W = V(y);
	value_type Wp = V.d(y);
	value_type a0 = value_type(sin(PI*t));
	value_type a0p = value_type(PI*cos(PI*t));
	value_type a0pp = value_type(-PI*PI*sin(PI*t));
	value_type ypp = -3 * ((a0p + dap) / (a0 + da))*yp +
		(yp*yp)*((xi*h*h*(1 - 6 * xi)*y) / (1 - h*h*xi*(1 - 6 * xi)*y*y))
		+ ((PI + dxmax)*(PI + dxmax))*Wp*((1 - h*h*xi*y*y)
			/ (1 - h*h*xi*(1 - 6 * xi)*y*y))
		+ ((PI + dxmax)*(PI + dxmax))*h*h*((4 * xi*y*(W + W0))
			/ (1 - h*h*xi*(1 - 6 * xi)*y*y));
	//Linear Fluctuations:
	value_type Dy = Y[5];
	value_type Dyp = Y[6];
	value_type Da = Y[7];
	value_type Dap = Y[8];
	value_type Ddxmax = Y[9];
	value_type Wpp = V.d2(y);
	value_type Wy0 = W;

	//Return the derative of the linearised fluctuations:
	dydt[5] = Dyp;
	dydt[6] = ((Ddxmax*((value_type(8.0)*((PI + dxmax)*(h*h*(Wy0
              *(xi*(y*((value_type(1.0)) / ((value_type(1.0)
              - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))))))
              + (value_type(2.0)*((PI + dxmax)*(Wp*((value_type(1.0)
              - (h*h*(xi*y*y)))*((value_type(1.0)) / ((value_type(1.0)
              - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))))))
              + ((Dyp*((value_type(-3.0)*(((value_type(1.0)) / ((a0
              + da)))*(a0p + dap))) + (value_type(2.0)*(h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*(y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*yp))))))))
              + (Dy*((value_type(8.0)*((PI + dxmax)*(PI
              + dxmax)*(h*h*h*h*(Wy0*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*xi*(y*y*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
              - (h*h*((value_type(1.0) + (value_type(-6.0)*xi))
              *(xi*y*y))))))))))))) + ((value_type(2.0)*((PI
              + dxmax)*(PI + dxmax)*(h*h*(Wp*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*(y*((value_type(1.0)
              - (h*h*(xi*y*y)))*((value_type(1.0)) / ((value_type(1.0)
              - (h*h*((value_type(1.0) + (value_type(-6.0)*xi))*(xi*y*y))))
              *(value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))))))))
              + ((value_type(4.0)*((PI + dxmax)*(PI + dxmax)*(h*h*(Wy0*(xi
              *((value_type(1.0)) / ((value_type(1.0)
              - (h*h*((value_type(1.0) + (value_type(-6.0)*xi))
              *(xi*y*y))))))))))) + ((value_type(2.0)*((PI + dxmax)*(PI
              + dxmax)*(h*h*(Wp*(xi*(y*((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))))))))
              + (((PI + dxmax)*(PI + dxmax)*(Wpp*((value_type(1.0)
              - (h*h*(xi*y*y)))*((value_type(1.0)) / ((value_type(1.0)
              - (h*h*((value_type(1.0) + (value_type(-6.0)*xi))*(xi*y*y)))))))))
              + ((value_type(2.0)*(h*h*h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*xi*(y*y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
              - (h*h*((value_type(1.0) + (value_type(-6.0)*xi))
              *(xi*y*y))))))*yp*yp)))))) + (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*((value_type(1.0)
              + (value_type(-6.0)*xi))*(xi*y*y))))))*yp*yp)))))))))))));
	dydt[7] = Dap;
	dydt[8] = ((((value_type(-2.0)) / (value_type(3.0)))*(Ddxmax*((PI
              + dxmax)*(h*h*(W*(((value_type(1.0)) / ((value_type(1.0)
              - (h*h*(xi*y*y)))))*(da + sin((PI*t)))))))))
              + ((((value_type(-1.0)) / (value_type(3.0)))
              *(Dyp*(h*h*(((value_type(1.0)) / ((value_type(1.0)
              - (h*h*(xi*y*y)))))*((da + sin((PI*t)))*((value_type(2.0)*yp)
              + (value_type(-3.0)*(xi*((value_type(2.0)*yp) + (y*((dap
              + (PI*cos((PI*t))))*((value_type(1.0)) / ((da
              + sin((PI*t)))))))))))))))) + (Dy*((((value_type(-1.0))
              / (value_type(3.0)))*(h*h*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))))*((da + sin((PI*t)))
              *(((PI + dxmax)*(PI + dxmax)*Wp) + (value_type(-3.0)*(xi*(ypp
              + (yp*((dap + (PI*cos((PI*t))))*((value_type(1.0)) / ((da
              + sin((PI*t))))))))))))))) + (((value_type(-2.0))
              / (value_type(3.0)))*(h*h*h*h*(xi*(y*(((value_type(1.0))
              / ((value_type(1.0) - (h*h*(xi*y*y)))*(value_type(1.0)
              - (h*h*(xi*y*y)))))*((da + sin((PI*t)))*(((PI
              + dxmax)*(PI + dxmax)*W) + (yp*yp + (value_type(-3.0)*(xi*(yp*yp
              + ((y*ypp) + (y*(yp*((dap + (PI*cos((PI*t))))*((value_type(1.0))
              / ((da + sin((PI*t))))))))))))))))))))))));
	dydt[9] = 0;
}//ode3::linearised
//------------------------------------------------------------------------------
 //ode3action functions:
 //ode3action main function:
template< class solution_type, class time_type, class value_type >
void ode3action< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type da = Y[2];
	value_type dap = Y[3];
	value_type dxmax = Y[4];
	value_type W = V(y);
	value_type Wp = V.d(y);
	value_type a0 = value_type(sin(PI*t));
	value_type a0p = value_type(PI*cos(PI*t));
	value_type a0pp = value_type(-PI*PI*sin(PI*t));
	if (a0 + da == value_type(0.0) || (1 - h*h*xi*y*y) == value_type(0.0)
		|| 1 - h*h*xi*(1 - 6 * xi)*y*y == value_type(0.0))
	{
		throw "ode3 attempted division by zero.";
	}
	value_type ypp = -3 * ((a0p + dap) / (a0 + da))*yp +
		(yp*yp)*((xi*h*h*(1 - 6 * xi)*y) / (1 - h*h*xi*(1 - 6 * xi)*y*y))
		+ ((PI + dxmax)*(PI + dxmax))*Wp*((1 - h*h*xi*y*y)
			/ (1 - h*h*xi*(1 - 6 * xi)*y*y))
		+ ((PI + dxmax)*(PI + dxmax))*h*h*((4 * xi*y*(W + W0))
			/ (1 - h*h*xi*(1 - 6 * xi)*y*y));
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = dap;
	dydt[3] = -a0pp - ((h*h*(a0 + da)) / (value_type(3.0)*(1
              - h*h*xi*y*y)))*(yp*yp
              + ((PI + dxmax)*(PI + dxmax))*(W + W0)
		- value_type(3.0)*xi*(yp*yp + y*ypp + ((a0p + dap) / (a0 + da))*y*yp));
	dydt[4] = 0;
	//Action:
	value_type xmax = dxmax + PI;
	dydt[5] = (value_type(2.0)*(PI*PI*((value_type(3.0)*(PI*(a0*a0*a0
             *((value_type(1.0)) / (h*h))))) + ((a0 + da)*(a0 + da)*(a0
            + da)*((PI + dxmax)*(((value_type(3.0)*(xi*((Wp*y)
            + (((value_type(1.0)) / ((PI + dxmax)*(PI + dxmax)))*yp*yp))))
            + (value_type(6.0)*(h*h*(xi*xi*(y*y*(((value_type(1.0))
            / ((value_type(1.0) - (h*h*((value_type(1.0)
            + (value_type(-6.0)*xi))*(xi*y*y))))))*((value_type(2.0)*W)
            + ((((value_type(1.0)) / (value_type(2.0)))*(((value_type(1.0))
            / ((PI + dxmax)*(PI + dxmax)))*yp*yp))
            + (value_type(-3.0)*(xi*((Wp*y) + (((value_type(1.0)) / ((PI
            + dxmax)*(PI + dxmax)))*yp*yp)))))))))))) - W))))));
}//ode3action::operator()
//------------------------------------------------------------------------------
//ode4 functions:
//Main ode4 function:
template< class solution_type, class time_type, class value_type >
void ode4< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (t == time_type(0.0) || t == time_type(PI))
	{
		throw "ode4 attempted division by zero.";
	}
	value_type a0 = value_type(sin(t));
	value_type a0p = value_type(cos(t));
	//ode:
	dydt[0] = yp;
	dydt[1] = -value_type(3.0)*yp / tan(t) + Wp + 12 * xi*y;
}//ode4::operator()


 //Linearised ode4 function:
template< class solution_type, class time_type, class value_type >
void ode4< solution_type, time_type, value_type >::linearised
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	//Returns a linearised version of the ode. This only makes
	//sense if we are using an
	//extended solution type where the first few variables are y, and the
	//latter variables the linearised fluctuations about y.

	//Call the ode to get the first few dydts (non-linear part, as
	//the linear part depends on these)
	(*this)(Y, dydt, t);

	//Define the variables:
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	//Non-linear part:
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type W = V(y);
	value_type Wp = V.d(y);
	//Linear Fluctuations:
	value_type Dy = Y[2];
	value_type Dyp = Y[3];
	value_type Wpp = V.d2(y);

	//Return the derative of the linearised fluctuations:
	dydt[2] = Dyp;
	dydt[3] = ((Dy*(Wpp + (value_type(12.0)*xi)))
            + (value_type(-3.0)*(Dyp / tan(t))));
}//ode4::linearised
//------------------------------------------------------------------------------
 //ode4action functions:
 //Main ode4action function:
template< class solution_type, class time_type, class value_type >
void ode4action< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (t == time_type(0.0) || t == time_type(PI))
	{
		throw "ode4 attempted division by zero.";
	}
	value_type a0 = value_type(sin(t));
	value_type a0p = value_type(cos(t));
	//ode:
	dydt[0] = yp;
	dydt[1] = -value_type(3.0)*yp / tan(t) + Wp + 12 * xi*y;
	//Action:
	/*
	dydt[2] = value_type(2.0)*value_type(PI)*value_type(PI)*(a0*a0*a0*
	(value_type(0.5)*yp*yp + W + 6*xi*y*y));
	*/
	dydt[2] = (value_type(2.0)*(a0*a0*a0*(PI*PI*(W +
		((value_type(6.0)*(xi*y*y)) +
			(((value_type(1.0)) / (value_type(2.0)))*yp*yp))))));

}//ode4action::operator()
//------------------------------------------------------------------------------
 //Class to call the Jacobian of the ode, rather than the ode itself.
 // Essentially a wrapper for the oderhs class
//Operator:
template< class solution_type, class time_type, class value_type >
void jacobian< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	//Note, this won't compile properly if the supplied template parameter
	//doesn't actually have a "linearised" member function.
	odeToCall.linearised(Y, dydt, t);
}
//------------------------------------------------------------------------------

#endif //DE_SITTER_INSTANTONS

//------------------------------------------------------------------------------
#ifdef ANTI_DE_SITTER_FLAT_INSTANTONS
//odeAdSFlat - with backreaction:
template< class solution_type, class time_type, class value_type >
void odeAdSFlat< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (a == value_type(0.0) || (value_type(1.0) - h*h*xi*y*y)
        == value_type(0.0) || value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y == value_type(0.0))
	{
		throw "ode1 attempted division by zero.";
	}
	value_type ypp = -value_type(3.0)*(ap / a)*yp
        + (yp*yp)*((xi*h*h*y*(value_type(1.0) - value_type(6.0)*xi))
        / (value_type(1.0) - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ Wp*((value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ h*h*(W + W0)*((value_type(4.0)*xi*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y));
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = ap;
	dydt[3] = -((h*h*a) / (value_type(3.0)*(value_type(1.0)
              - h*h*xi*y*y)))*(yp*yp + W0 + W
		- value_type(3.0)*xi*(yp*yp + y*ypp + (ap / a)*y*yp));
}//odeAdSFlat::operator()
//------------------------------------------------------------------------------
 //Linearised version of odeAdSFlat:
template< class solution_type, class time_type, class value_type >
void odeAdSFlat< solution_type, time_type, value_type >::linearised
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	//Returns a linearised version of the ode. This only makes sense if
	//we are using an
	//extended solution type where the first few variables are y, and the
	//latter variables the linearised fluctuations about y.

	//Call the ode to get the first few dydts (non-linear part, as
	//the linear part depends on these)
	(*this)(Y, dydt, t);

	//Define the variables:
	//Non-linear part:
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type W = V(y);
	value_type Wp = V.d(y);
	value_type ypp = -value_type(3.0)*(ap / a)*yp
                     + (yp*yp)*((xi*h*h*y*(value_type(1.0)
                     - value_type(6.0)*xi)) / (value_type(1.0)
                     - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
                     + Wp*((value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
                     - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
                     + h*h*(W + W0)*((value_type(4.0)*xi*y)
                     / (value_type(1.0) - h*h*xi*(value_type(1.0)
                     - value_type(6.0)*xi)*y*y));
	//Linear Fluctuations:
	value_type Dy = Y[4];
	value_type Dyp = Y[5];
	value_type Da = Y[6];
	value_type Dap = Y[7];
	value_type Wpp = V.d2(y);
	value_type Wy0 = W;

	//Return the derivative of the linearised fluctuations:
	dydt[4] = Dyp;
	dydt[5] = ((value_type(3.0)*(((value_type(1.0))
              / (a*a))*(ap*(Da*yp)))) + ((value_type(-3.0)*(((value_type(1.0))
              / (a))*(Dap*yp))) + ((Dyp*((value_type(-3.0)*(((value_type(1.0))
              / (a))*ap)) + (value_type(2.0)*(h*h*((value_type(1.0)
             + (value_type(-6.0)*xi))*(xi*(y*(((value_type(1.0))
             / ((value_type(1.0) - (h*h*((value_type(1.0)
             + (value_type(-6.0)*xi))*(xi*y*y))))))*yp))))))))
             + (Dy*((value_type(8.0)*(h*h*h*h*((W0 + Wy0)*((value_type(1.0)
             + (value_type(-6.0)*xi))*(xi*xi*(y*y*((value_type(1.0))
             / ((value_type(1.0) - (h*h*((value_type(1.0)
             + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
             - (h*h*((value_type(1.0) + (value_type(-6.0)*xi))
             *(xi*y*y)))))))))))) + ((value_type(2.0)*(h*h*(Wp
             *((value_type(1.0) + (value_type(-6.0)*xi))*(xi*(y
             *((value_type(1.0) - (h*h*(xi*y*y)))*((value_type(1.0))
             / ((value_type(1.0) - (h*h*((value_type(1.0)
             + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
             - (h*h*((value_type(1.0)
             + (value_type(-6.0)*xi))*(xi*y*y)))))))))))))
             + ((value_type(4.0)*(h*h*((W0 + Wy0)*(xi*((value_type(1.0))
             / ((value_type(1.0) - (h*h*((value_type(1.0)
             + (value_type(-6.0)*xi))*(xi*y*y))))))))))
             + ((value_type(2.0)*(h*h*(Wp*(xi*(y*((value_type(1.0))
             / ((value_type(1.0) - (h*h*((value_type(1.0)
             + (value_type(-6.0)*xi))*(xi*y*y)))))))))))
             + ((Wpp*((value_type(1.0) - (h*h*(xi*y*y)))*((value_type(1.0))
             / ((value_type(1.0) - (h*h*((value_type(1.0)
             + (value_type(-6.0)*xi))*(xi*y*y))))))))
            + ((value_type(2.0)*(h*h*h*h*((value_type(1.0)
            + (value_type(-6.0)*xi))*(value_type(1.0)
            + (value_type(-6.0)*xi))*(xi*xi*(y*y*(((value_type(1.0))
            / ((value_type(1.0) - (h*h*((value_type(1.0)
            + (value_type(-6.0)*xi))*(xi*y*y))))*(value_type(1.0)
            - (h*h*((value_type(1.0)
            + (value_type(-6.0)*xi))*(xi*y*y))))))*yp*yp))))))
            + (h*h*((value_type(1.0)
            + (value_type(-6.0)*xi))*(xi*(((value_type(1.0))
            / ((value_type(1.0) - (h*h*((value_type(1.0)
            + (value_type(-6.0)*xi))*(xi*y*y))))))*yp*yp))))))))))))));
	dydt[6] = Dap;
	dydt[7] = ((Dap*(h*h*(xi*(y*(((value_type(1.0)) / ((value_type(1.0)
                - (h*h*(xi*y*y)))))*yp))))) + ((((value_type(-1.0))
                / (value_type(3.0)))*(a*(Dyp*(h*h*(((value_type(1.0))
                / ((value_type(1.0) - (h*h*(xi*y*y)))))*((value_type(2.0)*yp)
                + (value_type(-3.0)*(xi*((((value_type(1.0)) / (a))*(ap*y))
                + (value_type(2.0)*yp)))))))))) + ((Dy*((((value_type(-1.0))
                / (value_type(3.0)))*(a*(h*h*(((value_type(1.0))
                / ((value_type(1.0) - (h*h*(xi*y*y)))))*(Wp
                + (value_type(-3.0)*(xi*((((value_type(1.0)) / (a))*(ap*yp))
                + ypp)))))))) + (((value_type(-2.0))
                / (value_type(3.0)))*(a*(h*h*h*h*(xi*(y*(((value_type(1.0))
                / ((value_type(1.0) - (h*h*(xi*y*y)))*(value_type(1.0)
                - (h*h*(xi*y*y)))))*(W + (W0 + (yp*yp
                + (value_type(-3.0)*(xi*((((value_type(1.0))
                / (a))*(ap*(y*yp))) + (yp*yp + (y*ypp))))))))))))))))
                + (Da*((((value_type(-1.0)) / (value_type(3.0)))*(h*h
                *(((value_type(1.0)) / ((value_type(1.0) - (h*h*(xi*y*y)))))
                *(W + (W0 + (yp*yp + (value_type(-3.0)
                *(xi*((((value_type(1.0)) / (a))*(ap*(y*yp))) + (yp*yp
                + (y*ypp))))))))))) - (((value_type(1.0))
                / (a))*(ap*(h*h*(xi*(y*(((value_type(1.0)) / ((value_type(1.0)
                - (h*h*(xi*y*y)))))*yp)))))))))));
}//odeAdSFlat::linearised
//------------------------------------------------------------------------------
 //odeAdSFlat - with backreaction:
template< class solution_type, class time_type, class value_type >
void odeAdSFlatAction< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (a == value_type(0.0) || (value_type(1.0) - h*h*xi*y*y)
        == value_type(0.0) || value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y == value_type(0.0))
	{
		throw "odeAdSFlatAction attempted division by zero.";
	}
	value_type ypp = -value_type(3.0)*(ap / a)*yp
        + (yp*yp)*((xi*h*h*y*(value_type(1.0) - value_type(6.0)*xi))
        / (value_type(1.0) - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ Wp*((value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ h*h*(W + W0)*((value_type(4.0)*xi*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y));
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = ap;
	dydt[3] = -((h*h*a) / (value_type(3.0)*(value_type(1.0)
              - h*h*xi*y*y)))*(yp*yp + W0 + W
              - value_type(3.0)*xi*(yp*yp + y*ypp + (ap / a)*y*yp));
	//Action:
	value_type H0 = sqrt(abs(W0) / value_type(3.0));
	value_type a0;
	if (abs(H0) < std::numeric_limits < value_type >::epsilon())
	{
		a0 = value_type(0.0);
	}
	else if(H0 < value_type(0.0))
	{
		value_type inter = H0*value_type(t);
		a0 = sinh(inter) / H0;
	}
	else
    {
        value_type inter = H0*value_type(t);
        a0 = sin(inter)/H0;
    }
	dydt[4] = (value_type(2.0)*(PI*PI*((a0*a0*a0*W0)
                + (a*a*a*(((value_type(3.0)*(xi*((Wp*y) + yp*yp)))
                + (value_type(6.0)*(h*h*(xi*xi*(y*y*(((value_type(1.0))
                / ((value_type(1.0) - (h*h*((value_type(1.0)
                + (value_type(-6.0)*xi))*(xi*y*y))))))*((value_type(2.0)*(W
                + W0)) + ((((value_type(1.0)) / (value_type(2.0)))*yp*yp)
                + (value_type(-3.0)*(xi*((Wp*y) + yp*yp))))))))))) - W)))));
}
//------------------------------------------------------------------------------
//Fixed background, flat space:
template< class solution_type, class time_type, class value_type >
void odeFlat< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	value_type y = Y[0];
	value_type yp = Y[1];
	//value_type W = V(y);
	value_type Wp = V.d(y);
	dydt[0] = yp;
	dydt[1] = -(value_type(3.0) / t)*yp + Wp;
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void odeFlat< solution_type, time_type, value_type >::linearised
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	(*this)(Y, dydt, t);//Calls the original ode to get the non-linear part.
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type Dy = Y[2];
	value_type Dyp = Y[3];
	//value_type W = V(y);
	//value_type Wp = V.d(y);
	value_type Wpp = V.d2(y);
	dydt[2] = Dyp;
	dydt[3] = ((Dy*Wpp) + (value_type(-3.0)*(Dyp*((value_type(1.0)) / t))));
}
//------------------------------------------------------------------------------
//Fixed background, flat space, with action:
template< class solution_type, class time_type, class value_type >
void odeFlatAction< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type W = V(y);
	value_type Wp = V.d(y);
	dydt[0] = yp;
	dydt[1] = -(value_type(3.0) / t)*yp + Wp;
	dydt[2] = value_type(2.0)*PI*PI*t*t*t*(((yp*yp) / value_type(2.0)) + W);
}
//------------------------------------------------------------------------------
//Fixed background, AdS space:
template< class solution_type, class time_type, class value_type >
void odeAdS< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	value_type y = Y[0];
	value_type yp = Y[1];
	//value_type W = V(y);
	value_type Wp = V.d(y);
	//value_type a0 = value_type(sinh(t));
	//value_type a0p = value_type(cosh(t));
	//ode:
	dydt[0] = yp;
	dydt[1] = -value_type(3.0)*yp / tanh(t) + Wp - 12 * xi*y;

}//odeAdS::operator()
//------------------------------------------------------------------------------
 //Linearised odeAdS function:
template< class solution_type, class time_type, class value_type >
void odeAdS< solution_type, time_type, value_type >::linearised
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	//Returns a linearised version of the ode. This only makes sense
	//if we are using an
	//extended solution type where the first few variables are y, and the
	//latter variables the linearised fluctuations about y.

	//Call the ode to get the first few dydts (non-linear part, as
	//the linear part depends on these)
	(*this)(Y, dydt, t);

	//Define the variables:
	//Non-linear part:
	value_type y = Y[0];
	value_type yp = Y[1];
	//value_type W = V(y);
	//value_type Wp = V.d(y);
	//Linear Fluctuations:
	value_type Dy = Y[2];
	value_type Dyp = Y[3];
	value_type Wpp = V.d2(y);

	//Return the derative of the linearised fluctuations:
	dydt[2] = Dyp;
	dydt[3] = ((Dy*(Wpp - (value_type(12.0)*xi)))
              + (value_type(-3.0)*(Dyp / tanh(t))));
}//odeAdS::linearised
//------------------------------------------------------------------------------
 //Fixed background, AdS space, with action:
template< class solution_type, class time_type, class value_type >
void odeAdSAction< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type W = V(y);
	value_type Wp = V.d(y);
	value_type H0 = sqrt(-W0 / 3);//AdS radius, in Planck units.
	time_type a0tt = sinh(H0*t);
	time_type a0ptt = H0*cosh(H0*t);
	value_type a0 = value_type(a0tt);
	value_type a0p = value_type(a0ptt);
	//ode:
	dydt[0] = yp;
	dydt[1] = -value_type(3.0)*yp / tanh(t) + Wp - 12 * xi*y;
	//Action:
	/*
	dydt[2] = value_type(2.0)*value_type(PI)*value_type(PI)*(a0*a0*a0*
	(value_type(0.5)*yp*yp + W + 6*xi*y*y));
	*/
	dydt[2] = (value_type(2.0)*(a0*a0*a0*(PI*PI*(W +
		(-(value_type(6.0)*(xi*y*y)) +
			(((value_type(1.0)) / (value_type(2.0)))*yp*yp))))));

}//odeAdSAction::operator()
//------------------------------------------------------------------------------
//Gravitational case, switching variables to theta = log(a + 1)
//to avoid overflow problems.
template< class solution_type, class time_type, class value_type >
void odeGravFriction< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
    using boost::math::expm1;
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type theta = Y[2];
	value_type thetap = Y[3];
	value_type yp2 = yp*yp;
	value_type y2 = y*y;
	value_type W = V(y);
	value_type Wp = V.d(y);
	value_type h2 = h*h;
	if ((this->_1p0 - h2*xi*y2)
        == this->_0p0 || this->_1p0 - h2*this->xi*(this->_1p0
        - this->_6p0*this->xi)*y2 == this->_0p0)
	{
		throw "ode1 attempted division by zero.";
	}

	value_type R = h2*(yp2 + this->_4p0*(W + W0)
                - this->_6p0*this->xi*(yp2 + y*Wp))
                /(this->_1p0 - this->xi*h2*(this->_1p0 - this->_6p0*this->xi));
        //R in units of H, the scale used.
    value_type _1_m_exp_mtheta = -expm1(-theta);
    value_type conformal = this->_1p0 - this->xi*y2*h2;
	value_type ypp = -this->_3p0*yp*thetap/_1_m_exp_mtheta
	 + Wp + this->xi*y*h2*R;
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = thetap;
	dydt[3] = -thetap*thetap - h2*_1_m_exp_mtheta*(yp2 + W + W0 -
                this->_3p0*this->xi*(yp2 + y*ypp))/(this->_3p0*conformal)
                + h2*this->xi*thetap*y*yp/conformal;
}//odeGravFriction::operator()
//------------------------------------------------------------------------------
//Gravitational case, tracking a'/a rather than a and a' separately.
//Includes action.
template< class solution_type, class time_type, class value_type >
void odeGravFrictionAction< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	using boost::math::expm1;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type theta = Y[2];
	value_type thetap = Y[3];
	value_type S = Y[4];
	value_type yp2 = yp*yp;
	value_type y2 = y*y;
	value_type W = V(y);
	value_type Wp = V.d(y);
	value_type h2 = h*h;
	if ((this->_1p0 - h2*xi*y2)
        == this->_0p0 || this->_1p0 - h2*this->xi*(this->_1p0
        - this->_6p0*this->xi)*y2 == this->_0p0)
	{
		throw "ode1 attempted division by zero.";
	}

	value_type R = h2*(yp2 + this->_4p0*(W + W0)
                - this->_6p0*this->xi*(yp2 + y*Wp))
                /(this->_1p0 - this->xi*h2*(this->_1p0 - this->_6p0*this->xi));
        //R in units of H, the scale used.
	value_type _1_m_exp_mtheta = -expm1(-theta);
    value_type conformal = this->_1p0 - this->xi*y2*h2;
	value_type ypp = -this->_3p0*yp*thetap/_1_m_exp_mtheta
	 + Wp + this->xi*y*h2*R;
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = thetap;
	dydt[3] = -thetap*thetap - h2*_1_m_exp_mtheta*(yp2 + W + W0 -
                this->_3p0*this->xi*(yp2 + y*ypp))/(this->_3p0*conformal)
                + h2*this->xi*thetap*y*yp/conformal;
	//Action:
	value_type H0 = sqrt(abs(W0) / this->_3p0);
	value_type a0;
	if (abs(H0) < std::numeric_limits < value_type >::epsilon())
	{
		a0 = value_type(0.0);
	}
	else
	{
		value_type inter = H0*value_type(t);
		a0 = W0 > this->_0p0 ? sin(inter)/H0 : sinh(inter) / H0;
	}
	value_type a = expm1(theta);
	dydt[4] = (this->_2p0*(PI*PI*((a0*a0*a0*W0)
                + (a*a*a*(((this->_3p0*(this->xi*((Wp*y) + yp*yp)))
                + (this->_6p0*(h2*(this->xi*this->xi*(y2*(((this->_1p0)
                / ((this->_1p0 - (h2*((this->_1p0
                + (-this->_6p0*this->xi))*(this->xi*y2))))))*((this->_2p0*(W
                + W0)) + ((((this->_1p0) / (this->_2p0))*yp2)
                + (-this->_3p0*(this->xi*((Wp*y) + yp2))))))))))) - W)))));
}//odeGravFrictionAction::operator()
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void ode_dS< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
    value_type y = Y[0];
    value_type yp = Y[1];
    value_type a = Y[2];
    value_type ap = Y[3];

    //Intermediates:
    value_type h2 = this->h*this->h;
    value_type y2 = y*y;
    value_type yp2 = yp*yp;
    value_type _1_m6xi = (this->_1p0 - this->_6p0*this->xi);
    value_type conformal = (this->_1p0 - this->xi*h2*_1_m6xi*y2);
    value_type factor = (this->_1p0 - this->xi*h2*y2);

    //Evaluate potential:
    value_type W = V(y);
    value_type Wp = V.d(y);

    //R without the V0 part:
   value_type R = h2*(yp2*_1_m6xi + this->_4p0*(W + W0)
                       - this->_6p0*this->xi*y*Wp)/conformal;
    //R with the V0 part:
    //value_type R0 = this->_4p0*this->W0/conformal;
    //value_type R = R0 + DR;

    value_type ypp = -this->_3p0*ap*yp/a  + Wp + this->xi*y*R;

    dydt[0] = yp;
    dydt[1] = ypp;
    dydt[2] = ap;
    /*dydt[3] = -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*(yp2
                      + y*Wp + this->xi*y2*R))/(this->_3p0*factor)
                - this->_2p0*h2*y*yp*ap*this->xi/factor;*/
    dydt[3] = -a*h2*(yp2 + W + this->W0
                     -this->_3p0*this->xi*(yp2 + y*ypp + ap*y*yp/a))
                     /(this->_3p0*factor);
    /*
    std::cout << "\n 1 - 6*xi = " << _1_m6xi;
    std::cout << "\n conformal = " << conformal;
    std::cout << "\n xi = " << xi;
    std::cout << "\n W0 = " << this->W0;
    std::cout << "\n W = " << W;
    std::cout << "\n h2 = " << h2;*/
    /*dydt[3] = -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*yp2
                     + this->_6p0*this->xi*y*yp*ap/a - this->_3p0*this->xi*y*Wp
                     -this->_3p0*this->xi*this->xi*y2*this->_6p0
                     *(this->_1p0 - ap*ap)/(a*a))/(this->_3p0*conformal);*/
    /*std::cout << "\nt = " << t;
    for(int i = 0;i < 4;i++)
    {
        std::cout << "\ny[" << i << "] = " << Y[i];
    }
    for(int i = 0;i < 4;i++)
    {
        std::cout << "\ndydt[" << i << "] = " << dydt[i];
    }*/
}//ode_dS::operator()
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void ode_dS_delta_a< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
    value_type y = Y[0];
    value_type yp = Y[1];
    value_type da = Y[2];
    value_type dap = Y[3];

    value_type arg = H0*value_type(t);
    value_type sinarg = sin(arg);
    value_type a = da + sinarg/H0;
    value_type ap = dap + cos(arg);

    //Intermediates:
    value_type h2 = this->h*this->h;
    value_type y2 = y*y;
    value_type yp2 = yp*yp;
    value_type _1_m6xi = (this->_1p0 - this->_6p0*this->xi);
    value_type conformal = (this->_1p0 - this->xi*h2*_1_m6xi*y2);
    value_type factor = (this->_1p0 - this->xi*h2*y2);

    //Evaluate potential:
    value_type W = V(y);
    value_type Wp = V.d(y);

    //R without the V0 part:
   value_type R = h2*(yp2*_1_m6xi + this->_4p0*(W + W0)
                       - this->_6p0*this->xi*y*Wp)/conformal;
    //R with the V0 part:
    //value_type R0 = this->_4p0*this->W0/conformal;
    //value_type R = R0 + DR;

    value_type ypp = -this->_3p0*ap*yp/a  + Wp + this->xi*y*R;

    dydt[0] = yp;
    dydt[1] = ypp;
    dydt[2] = dap;
    /*dydt[3] = H0*sinarg -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*(yp2
                      + y*Wp + this->xi*y2*R))/(this->_3p0*factor)
                - this->_2p0*h2*y*yp*ap*this->xi/factor;*/
    dydt[3] = H0*sinarg -a*h2*(yp2 + W + this->W0
                     -this->_3p0*this->xi*(yp2 + y*ypp + ap*y*yp/a))
                     /(this->_3p0*factor);
    /*
    std::cout << "\n 1 - 6*xi = " << _1_m6xi;
    std::cout << "\n conformal = " << conformal;
    std::cout << "\n xi = " << xi;
    std::cout << "\n W0 = " << this->W0;
    std::cout << "\n W = " << W;
    std::cout << "\n h2 = " << h2;*/
    /*dydt[3] = -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*yp2
                     + this->_6p0*this->xi*y*yp*ap/a - this->_3p0*this->xi*y*Wp
                     -this->_3p0*this->xi*this->xi*y2*this->_6p0
                     *(this->_1p0 - ap*ap)/(a*a))/(this->_3p0*conformal);*/
    /*std::cout << "\nt = " << t;
    for(int i = 0;i < 4;i++)
    {
        std::cout << "\ny[" << i << "] = " << Y[i];
    }
    for(int i = 0;i < 4;i++)
    {
        std::cout << "\ndydt[" << i << "] = " << dydt[i];
    }*/
}//ode_dS::operator()
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void ode_dS_action_delta_a< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
    using boost::math::constants::pi;
	//using boost::math::expm1;
	time_type PI = pi< time_type >();

    value_type y = Y[0];
    value_type yp = Y[1];
    value_type da = Y[2];
    value_type dap = Y[3];
    value_type S = Y[4];

    value_type arg = H0*value_type(t);
    value_type sinarg = sin(arg);
    value_type a = da + sinarg/H0;
    value_type ap = dap + cos(arg);

    //Intermediates:
    value_type h2 = this->h*this->h;
    value_type y2 = y*y;
    value_type yp2 = yp*yp;
    value_type _1_m6xi = (this->_1p0 - this->_6p0*this->xi);
    value_type conformal = (this->_1p0 - this->xi*h2*_1_m6xi*y2);
    value_type factor = (this->_1p0 - this->xi*h2*y2);

    //Evaluate potential:
    value_type W = V(y);
    value_type Wp = V.d(y);

    //R without the V0 part:
    value_type R = h2*(yp2*_1_m6xi + this->_4p0*(W + W0)
                       - this->_6p0*this->xi*y*Wp)/conformal;
    //R with the V0 part:
    //value_type R0 = this->_4p0*this->W0/conformal;
    //value_type R = R0 + DR;

    value_type ypp = -this->_3p0*ap*yp/a  + Wp + this->xi*y*R;

    dydt[0] = yp;
    dydt[1] = ypp;
    dydt[2] = dap;
    /*std::cout << "\n 1 - 6*xi = " << _1_m6xi;
    std::cout << "\n conformal = " << conformal;
    std::cout << "\n xi = " << xi;
    std::cout << "\n W0 = " << this->W0;
    std::cout << "\n W = " << W;
    std::cout << "\n h2 = " << h2;*/
    /*dydt[3] = -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*(yp2
                      + y*Wp + this->xi*y2*R))/(this->_3p0*factor)
                - this->_2p0*h2*y*yp*ap*this->xi/factor;*/
    dydt[3] = H0*sinarg -a*h2*(yp2 + W + this->W0
                     -this->_3p0*this->xi*(yp2 + y*ypp + ap*y*yp/a))
                     /(this->_3p0*factor);
    /*dydt[3] = -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*yp2
                     + this->_6p0*this->xi*y*yp*ap/a - this->_3p0*this->xi*y*Wp
                     -this->_3p0*this->xi*this->xi*y2*this->_6p0
                     *(this->_1p0 - ap*ap)/(a*a))/(this->_3p0*conformal);*/

    value_type a3 = a*a*a;
    /*dydt[4] = this->_2p0*PI*PI*a3*( - W - this->W0 + this->_3p0*this->xi
                                   *(yp2 + y*Wp)
                                   + this->_3p0*this->xi*this->xi*y2*R);*/
    value_type PI2 = PI*PI;
    //Only valid for xi = 0 at the moment. Neglects a constant correction
    //that we need to add on (but this was the same as before anyway).
    dydt[4] = -this->_2p0*PI2*a3*W;
    dydt[5] =  - this->_6p0*PI2*h2*(this->_3p0*sinarg*sinarg*da + this->_3p0*
                 H0*sinarg*da*da + H0*H0*da*da*da);
    /*std::cout << "\nt = " << t;
    for(int i = 0;i < 5;i++)
    {
        std::cout << "\ny[" << i << "] = " << Y[i];
    }
    for(int i = 0;i < 5;i++)
    {
        std::cout << "\ndydt[" << i << "] = " << dydt[i];
    }*/
}//ode_dS_action::operator()
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void ode_dS_action< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
    using boost::math::constants::pi;
	//using boost::math::expm1;
	time_type PI = pi< time_type >();

    value_type y = Y[0];
    value_type yp = Y[1];
    value_type a = Y[2];
    value_type ap = Y[3];
    value_type S = Y[4];

    //Intermediates:
    value_type h2 = this->h*this->h;
    value_type y2 = y*y;
    value_type yp2 = yp*yp;
    value_type _1_m6xi = (this->_1p0 - this->_6p0*this->xi);
    value_type conformal = (this->_1p0 - this->xi*h2*_1_m6xi*y2);
    value_type factor = (this->_1p0 - this->xi*h2*y2);

    //Evaluate potential:
    value_type W = V(y);
    value_type Wp = V.d(y);

    //R without the V0 part:
    value_type R = h2*(yp2*_1_m6xi + this->_4p0*(W + W0)
                       - this->_6p0*this->xi*y*Wp)/conformal;
    //R with the V0 part:
    //value_type R0 = this->_4p0*this->W0/conformal;
    //value_type R = R0 + DR;

    value_type ypp = -this->_3p0*ap*yp/a  + Wp + this->xi*y*R;

    dydt[0] = yp;
    dydt[1] = ypp;
    dydt[2] = ap;
    /*std::cout << "\n 1 - 6*xi = " << _1_m6xi;
    std::cout << "\n conformal = " << conformal;
    std::cout << "\n xi = " << xi;
    std::cout << "\n W0 = " << this->W0;
    std::cout << "\n W = " << W;
    std::cout << "\n h2 = " << h2;*/
    /*dydt[3] = -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*(yp2
                      + y*Wp + this->xi*y2*R))/(this->_3p0*factor)
                - this->_2p0*h2*y*yp*ap*this->xi/factor;*/
    dydt[3] = -a*h2*(yp2 + W + this->W0
                     -this->_3p0*this->xi*(yp2 + y*ypp + ap*y*yp/a))
                     /(this->_3p0*factor);
    /*dydt[3] = -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*yp2
                     + this->_6p0*this->xi*y*yp*ap/a - this->_3p0*this->xi*y*Wp
                     -this->_3p0*this->xi*this->xi*y2*this->_6p0
                     *(this->_1p0 - ap*ap)/(a*a))/(this->_3p0*conformal);*/

    value_type a3 = a*a*a;
    /*dydt[4] = this->_2p0*PI*PI*a3*( - W - this->W0 + this->_3p0*this->xi
                                   *(yp2 + y*Wp)
                                   + this->_3p0*this->xi*this->xi*y2*R);*/
    dydt[4] = this->_2p0*PI*PI*a3*(W + this->W0 + yp2/this->_2p0
                                   - factor*R/(this->_2p0*h2));
    /*std::cout << "\nt = " << t;
    for(int i = 0;i < 5;i++)
    {
        std::cout << "\ny[" << i << "] = " << Y[i];
    }
    for(int i = 0;i < 5;i++)
    {
        std::cout << "\ndydt[" << i << "] = " << dydt[i];
    }*/
}//ode_dS_action::operator()
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void ode_dS_switch< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
    value_type y = Y[0];
    value_type yp = Y[1];
    value_type arg = this->H0*value_type(t) + this->phi;
    value_type a = sin(arg)/this->H0;
    value_type ap = cos(arg);

    //Intermediates:
    value_type h2 = this->h*this->h;
    value_type y2 = y*y;
    value_type yp2 = yp*yp;
    value_type _1_m6xi = (this->_1p0 - this->_6p0*this->xi);
    value_type conformal = (this->_1p0 - this->xi*h2*_1_m6xi*y2);
    value_type factor = (this->_1p0 - this->xi*h2*y2);

    //Evaluate potential:
    value_type W = V(y);
    value_type Wp = V.d(y);

    //R without the V0 part:
   value_type R = h2*(yp2*_1_m6xi + this->_4p0*(W + W0)
                       - this->_6p0*this->xi*y*Wp)/conformal;
    //R with the V0 part:
    //value_type R0 = this->_4p0*this->W0/conformal;
    //value_type R = R0 + DR;

    value_type ypp = -this->_3p0*ap*yp/a  + Wp + this->xi*y*R;

    dydt[0] = yp;
    dydt[1] = ypp;
    //dydt[2] = ap;
    /*dydt[3] = -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*(yp2
                      + y*Wp + this->xi*y2*R))/(this->_3p0*factor)
                - this->_2p0*h2*y*yp*ap*this->xi/factor;*/
    /*dydt[3] = -a*h2*(yp2 + W + this->W0
                     -this->_3p0*this->xi*(yp2 + y*ypp + ap*y*yp/a))
                     /(this->_3p0*factor);*/

    /*
    std::cout << "\n 1 - 6*xi = " << _1_m6xi;
    std::cout << "\n conformal = " << conformal;
    std::cout << "\n xi = " << xi;
    std::cout << "\n W0 = " << this->W0;
    std::cout << "\n W = " << W;
    std::cout << "\n h2 = " << h2;*/
    /*dydt[3] = -h2*a*(yp2 + W + this->W0 - this->_3p0*this->xi*yp2
                     + this->_6p0*this->xi*y*yp*ap/a - this->_3p0*this->xi*y*Wp
                     -this->_3p0*this->xi*this->xi*y2*this->_6p0
                     *(this->_1p0 - ap*ap)/(a*a))/(this->_3p0*conformal);*/
    /*std::cout << "\nt = " << t;
    for(int i = 0;i < 4;i++)
    {
        std::cout << "\ny[" << i << "] = " << Y[i];
    }
    for(int i = 0;i < 4;i++)
    {
        std::cout << "\ndydt[" << i << "] = " << dydt[i];
    }*/
}//ode_dS::operator()
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void ode_dS_action_switch< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
    using boost::math::constants::pi;
	//using boost::math::expm1;
	time_type PI = pi< time_type >();

    value_type y = Y[0];
    value_type yp = Y[1];

    //Analytic scale factor:
    value_type arg = value_type(t)*this->H0 + this->phi;
    value_type a = sin(arg)/this->H0;
    value_type ap = cos(arg);


    value_type S = Y[2];

    //Intermediates:
    value_type h2 = this->h*this->h;
    value_type y2 = y*y;
    value_type yp2 = yp*yp;
    value_type xi2 = this->xi*this->xi;
    value_type _1_m6xi = (this->_1p0 - this->_6p0*this->xi);
    value_type conformal = (this->_1p0 - this->xi*h2*_1_m6xi*y2);
    value_type conformal0 =(this->_1p0 - this->xi*h2*_1_m6xi*this->fv*this->fv);
    value_type factor = (this->_1p0 - this->xi*h2*y2);
    const value_type _12p0 = this->_3p0*this->_4p0;

    //Evaluate potential:
    value_type W = V(y);
    value_type Wp = V.d(y);

    //R without the V0 part:
    value_type R = h2*(yp2*_1_m6xi + this->_4p0*(W + W0)
                       - this->_6p0*this->xi*y*Wp)/conformal;
    //R with the V0 part:
    //value_type R0 = this->_4p0*this->W0/conformal;
    //value_type R = R0 + DR;

    value_type ypp = -this->_3p0*ap*yp/a  + Wp + this->xi*y*R;

    dydt[0] = yp;
    dydt[1] = ypp;
    //dydt[2] = ap;
    /*dydt[3] = -a*h2*(yp2 + W + this->W0
                     -this->_3p0*this->xi*(yp2 + y*ypp + ap*y*yp/a))
                     /(this->_3p0*factor);*/

    value_type a3 = a*a*a;
    dydt[2] = this->_2p0*PI*PI*a3*( - W + this->_3p0*this->xi*(yp2 + y*Wp)
                + (this->_3p0*xi2*h2*y2/conformal)*
                (yp2*(this->_1p0 - this->_6p0*this->xi)
                + this->_4p0*W - this->_6p0*this->xi*y*Wp)
                +_12p0*W0*xi2*h2*( yp2/conformal - fv*fv/conformal0 ) );
    //std::cout << "\nS(t) = " << S;
    //std::cout << "\ndS/dt = " << dydt[2];
    /*std::cout << "\nt = " << t;
    for(int i = 0;i < 5;i++)
    {
        std::cout << "\ny[" << i << "] = " << Y[i];
    }
    for(int i = 0;i < 5;i++)
    {
        std::cout << "\ndydt[" << i << "] = " << dydt[i];
    }*/
}//ode_dS_action::operator()
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void odeAdSFlat_delta_a< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type da = Y[4];
	value_type dap = Y[5];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (a == value_type(0.0) || (value_type(1.0) - h*h*xi*y*y)
        == value_type(0.0) || value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y == value_type(0.0))
	{
		throw "ode1 attempted division by zero.";
	}
	value_type ypp = -value_type(3.0)*(ap / a)*yp
        + (yp*yp)*((xi*h*h*y*(value_type(1.0) - value_type(6.0)*xi))
        / (value_type(1.0) - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ Wp*((value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ h*h*(W + W0)*((value_type(4.0)*xi*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y));
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = ap;
	dydt[3] = -((h*h*a) / (value_type(3.0)*(value_type(1.0)
              - h*h*xi*y*y)))*(yp*yp + W0 + W
		- value_type(3.0)*xi*(yp*yp + y*ypp + (ap / a)*y*yp));
    //Difference between a and chi:
    dydt[4] = dap;
    dydt[5] = dydt[3];
}//odeAdSFlat::operator()
//------------------------------------------------------------------------------
//odeAdSFlat - with backreaction:
template< class solution_type, class time_type, class value_type >
void odeAdSFlatAction_delta_a< solution_type, time_type, value_type >
::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	using boost::math::constants::pi;
	time_type PI = pi< time_type >();
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type a = Y[2];
	value_type ap = Y[3];
	value_type da = Y[4];
	value_type dap = Y[5];
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (a == value_type(0.0) || (value_type(1.0) - h*h*xi*y*y)
        == value_type(0.0) || value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y == value_type(0.0))
	{
		throw "odeAdSFlatAction attempted division by zero.";
	}
	value_type ypp = -value_type(3.0)*(ap / a)*yp
        + (yp*yp)*((xi*h*h*y*(value_type(1.0) - value_type(6.0)*xi))
        / (value_type(1.0) - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ Wp*((value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ h*h*(W + W0)*((value_type(4.0)*xi*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y));
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
	dydt[2] = ap;
	dydt[3] = -((h*h*a) / (value_type(3.0)*(value_type(1.0)
              - h*h*xi*y*y)))*(yp*yp + W0 + W
              - value_type(3.0)*xi*(yp*yp + y*ypp + (ap / a)*y*yp));
    //Difference between a and chi:
    dydt[4] = dap;
    dydt[5] = dydt[3];
	//Action:
	value_type H0 = sqrt(abs(W0) / value_type(3.0));
	value_type a0;
	if (abs(H0) < std::numeric_limits < value_type >::epsilon())
	{
		a0 = value_type(0.0);
	}
	else if(H0 < value_type(0.0))
	{
		value_type inter = H0*value_type(t);
		a0 = sinh(inter) / H0;
	}
	else
    {
        value_type inter = H0*value_type(t);
        a0 = sin(inter)/H0;
    }
	dydt[6] = (value_type(2.0)*(PI*PI*((a0*a0*a0*W0)
                + (a*a*a*(((value_type(3.0)*(xi*((Wp*y) + yp*yp)))
                + (value_type(6.0)*(h*h*(xi*xi*(y*y*(((value_type(1.0))
                / ((value_type(1.0) - (h*h*((value_type(1.0)
                + (value_type(-6.0)*xi))*(xi*y*y))))))*((value_type(2.0)*(W
                + W0)) + ((((value_type(1.0)) / (value_type(2.0)))*yp*yp)
                + (value_type(-3.0)*(xi*((Wp*y) + yp*yp))))))))))) - W)))));
}
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void ode_dS_fixed< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
	value_type y = Y[0];
	value_type yp = Y[1];
	value_type arg = this->H0*value_type(t);
    value_type a = sin(arg)/this->H0;
	value_type ap = cos(arg);
	value_type W = V(y);
	value_type Wp = V.d(y);
	if (a == value_type(0.0) || (value_type(1.0) - h*h*xi*y*y)
        == value_type(0.0) || value_type(1.0) - h*h*xi*(value_type(1.0)
        - value_type(6.0)*xi)*y*y == value_type(0.0))
	{
		throw "ode1 attempted division by zero.";
	}
	value_type ypp = -value_type(3.0)*(ap / a)*yp
        + (yp*yp)*((xi*h*h*y*(value_type(1.0) - value_type(6.0)*xi))
        / (value_type(1.0) - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ Wp*((value_type(1.0) - h*h*xi*y*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y))
		+ h*h*(W + W0)*((value_type(4.0)*xi*y) / (value_type(1.0)
        - h*h*xi*(value_type(1.0) - value_type(6.0)*xi)*y*y));
	//ode:
	dydt[0] = yp;
	dydt[1] = ypp;
}//ode_dS_fixed::operator()
//------------------------------------------------------------------------------
template< class solution_type, class time_type, class value_type >
void ode_dS_fixed_action< solution_type, time_type, value_type >::operator()
(const solution_type& Y, solution_type& dydt, const time_type& t)
{
    using boost::math::constants::pi;
	//using boost::math::expm1;
	time_type PI = pi< time_type >();

    value_type y = Y[0];
    value_type yp = Y[1];
    value_type arg = this->H0*value_type(t);
    value_type a = sin(arg)/this->H0;
	value_type ap = cos(arg);
    value_type S = Y[2];

    //Intermediates:
    value_type h2 = this->h*this->h;
    value_type y2 = y*y;
    value_type yp2 = yp*yp;
    value_type _1_m6xi = (this->_1p0 - this->_6p0*this->xi);
    value_type conformal = (this->_1p0 - this->xi*h2*_1_m6xi*y2);
    value_type factor = (this->_1p0 - this->xi*h2*y2);

    //Evaluate potential:
    value_type W = V(y);
    value_type Wp = V.d(y);

    //R without the V0 part:
    value_type R = h2*(yp2*_1_m6xi + this->_4p0*(W + W0)
                       - this->_6p0*this->xi*y*Wp)/conformal;
    //R with the V0 part:
    //value_type R0 = this->_4p0*this->W0/conformal;
    //value_type R = R0 + DR;

    value_type ypp = -this->_3p0*ap*yp/a  + Wp + this->xi*y*R;
    value_type a3 = a*a*a;

    dydt[0] = yp;
    dydt[1] = ypp;
    dydt[2] = this->_2p0*PI*PI*a3*(yp*yp/_2p0 + W);
    /*std::cout << "\nt = " << t;
    for(int i = 0;i < 5;i++)
    {
        std::cout << "\ny[" << i << "] = " << Y[i];
    }
    for(int i = 0;i < 5;i++)
    {
        std::cout << "\ndydt[" << i << "] = " << dydt[i];
    }*/
}//ode_dS_fixed_action::operator()
//------------------------------------------------------------------------------
#endif //ANTI_DE_SITTER_FLAT_INSTANTONS

#endif //ODE_RHS_CODE_h
