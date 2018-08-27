#ifndef POTENTIALS_CODE_H
#define POTENTIALS_CODE_H

//Needed core code:
#include "potentials.h"

//Needed libraries:
#include <sstream>
//==============================================================================
//------------------------------------------------------------------------------
//Code for testPotential
#ifdef USING_POTENTIAL_TEST1
template< class value_type >
value_type testPotential< value_type >::operator()(value_type y)
{
	value_type x = yscale*y;
	return ((a*log((x*x) / (m*m))*log((x*x) / (m*m)) - b*log((x*x) / (m*m))
		+ c)*(x*x*x*x) / value_type(4.0)) / (yscale*yscale*yscale*yscale);
}
//------------------------------------------------------------------------------
template< class value_type >
void testPotential< value_type >::operator()
(value_type* y, value_type* arrayOut, int numel)
{
	//Vectorised version, for speed (could just call the above many times
	//, but that would lead to overheads calling the function each time.
	for (int i = 0; i < numel; i++)
	{
		value_type x = yscale*y[i];
		arrayOut[i] = ((a*log((x*x) / (m*m))*log((x*x) / (m*m)) - b*log((x*x)
                      / (m*m)) + c)*(x*x*x*x) / value_type(4.0))
                      / (yscale*yscale*yscale*yscale);
	}
}
//------------------------------------------------------------------------------
//First derivatives:
template< class value_type >
value_type testPotential< value_type >::d(value_type y)
{
	value_type x = yscale*y;
	return (-(value_type(1.0) / value_type(2.0))*(x*x*x)*(b - value_type(2.0)*c
            - value_type(2.0)*(a - b)*log((x*x) / (m*m))
            - value_type(2.0)*a*log((x*x) / (m*m))*log((x*x) / (m*m))))
            / (yscale*yscale*yscale);
}
//------------------------------------------------------------------------------
template< class value_type >
void testPotential< value_type >::d
(value_type* y, value_type* arrayOut, int numel)
{
	for (int i = 0; i < numel; i++)
	{
		value_type x = yscale*y[i];
		arrayOut[i] = (-(value_type(1.0) / value_type(2.0))*(x*x*x)*(b
                      - value_type(2.0)*c - value_type(2.0)*(a - b)*log((x*x)
                      / (m*m)) 	- value_type(2.0)*a*log((x*x)
                      / (m*m))*log((x*x) / (m*m)))) / (yscale*yscale*yscale);
	}
}
//------------------------------------------------------------------------------
//Second derivatives:
template< class value_type >
value_type testPotential< value_type >::d2(value_type y)
{
	value_type x = yscale*y;
	return ((value_type(1.0) / value_type(2.0))*(x*x)*(value_type(4.0)*a
            - value_type(7.0)*b + value_type(6.0)*c + value_type(2.0)
            *(value_type(7.0)*a - value_type(3.0)*b)*log((x*x) / (m*m))
            + value_type(6.0)*a*log((x*x) / (m*m))*log((x*x) / (m*m))))
            / (yscale*yscale);
}
//------------------------------------------------------------------------------
template< class value_type >
void testPotential< value_type >::d2
(value_type* y, value_type* arrayOut, int numel)
{
	for (int i = 0; i < numel; i++)
	{
		value_type x = yscale*y[i];
		arrayOut[i] = ((value_type(1.0) / value_type(2.0))*(x*x)
                      *(value_type(4.0)*a - value_type(7.0)*b
                      + value_type(6.0)*c + value_type(2.0)*(value_type(7.0)*a
                      - value_type(3.0)*b)*log((x*x) / (m*m))
                      + value_type(6.0)*a*log((x*x) / (m*m))*log((x*x)
                      / (m*m)))) / (yscale*yscale);
	}
}
//------------------------------------------------------------------------------
//Third derivatives:
template< class value_type >
value_type testPotential< value_type >::d3(value_type y)
{
	value_type x = yscale*y;
	return (x*(value_type(18.0)*a - value_type(13.0)*b + value_type(6.0)*c
            + (value_type(26.0)*a - value_type(6.0)*b)*log((x*x) / (m*m))
            + value_type(6.0)*a*log((x*x) / (m*m))*log((x*x) / (m*m))))
            / (yscale);
}
//------------------------------------------------------------------------------
template< class value_type >
void testPotential< value_type >::d3
(value_type* y, value_type* arrayOut, int numel)
{
	for (int i = 0; i < numel; i++)
	{
		value_type x = yscale*y[i];
		arrayOut[i] = (x*(value_type(18.0)*a - value_type(13.0)*b
                       + value_type(6.0)*c + (value_type(26.0)*a
                       - value_type(6.0)*b)*log((x*x) / (m*m))
                       + value_type(6.0)*a*log((x*x) / (m*m))*log((x*x)
                       / (m*m)))) / (yscale);
	}
}
//------------------------------------------------------------------------------
//Fourth derivatives:
template< class value_type >
value_type testPotential< value_type >::d4(value_type y)
{
	value_type x = yscale*y;
	return value_type(70.0)*a - value_type(25.0)*b + value_type(6.0)*c
           + (value_type(50.0)*a - value_type(6.0)*b)*log((x*x) / (m*m)) +
           value_type(6.0)*a*log((x*x) / (m*m))*log((x*x) / (m*m));
}
//------------------------------------------------------------------------------
template< class value_type >
void testPotential< value_type >::d4
(value_type* y, value_type* arrayOut, int numel)
{
	for (int i = 0; i < numel; i++)
	{
		value_type x = yscale*y[i];
		arrayOut[i] = value_type(70.0)*a - value_type(25.0)*b
                      + value_type(6.0)*c + (value_type(50.0)*a
                      - value_type(6.0)*b)*log((x*x) / (m*m))
                      +	value_type(6.0)*a*log((x*x) / (m*m))*log((x*x) / (m*m));
	}
}
//------------------------------------------------------------------------------
//==============================================================================
#endif


//==============================================================================
//Standard Model Potential:
#ifdef USING_SM_POTENTIAL
template< class value_type >
value_type SMHiggsPotentialSpline< value_type >::SinglePointSpline
(value_type& y, int m)
{
	//Construct a cubic spline approximation of the running coupling:
	//NB - argument y MUST be in units of GeV.
	int n = this->ChooseSpline(y, this->xdata, this->Ndata - 1);
	//const value_type _1p0 = value_type(1.0);
	return this->Spline(y,kdata[m],ydata[m],n);
}
//------------------------------------------------------------------------------
/*
template< class value_type >
value_type SMHiggsPotentialSpline< value_type >::SplineDerivative
(value_type& y,int m)
{
	//Derivative of the cubic used for interpolation:
	int n = this->ChooseSpline(y, this->xdata, this->Ndata - 1);
	//const value_type _2p0 = value_type(2.0);
	return this->SplineDerivative(y,m,n);
}*/
//------------------------------------------------------------------------------
//Same as the above, but assuming n is already known.
template< class value_type >
value_type SMHiggsPotentialSpline< value_type >::SinglePointSpline
(value_type& y, int m,int n)
{
	//Construct a cubic spline approximation of the running coupling:
	//NB - argument y MUST be in units of GeV.
	//Replace log y^2/M^2 with 0 if y = 0, to keep things finite, since
	//the result will be zero anyway.
	return this->Spline(y,kdata[m],ydata[m],n);
}
//------------------------------------------------------------------------------
/*
template< class value_type >
value_type SMHiggsPotentialSpline< value_type >::SplineDerivative
(value_type& y,int m,int n)
{
	//Derivative of the cubic used for interpolation:
	//const value_type _2p0 = value_type(2.0);

	//Replace log y^2/M^2 with 0 if y = 0, to keep things finite, since
	//the result will be zero anyway.
	value_type logfactor = (y != this->_0p0 ? log(y*y / (this->M*this->M))
                          : this->_0p0 );

	//const value_type _1p0 = value_type(1.0);

	value_type t = (logfactor - this->xdata[n]) /
                    (this->xdata[n + 1] - this->xdata[n]);
	value_type tp = (this->_1p0)/(this->xdata[n + 1] - this->xdata[n]);
	value_type a = this->kdata[m][n] * (this->xdata[n + 1] - this->xdata[n])
                    - (this->ydata[m][n + 1] - ydata[m][n]);
	value_type b = -this->kdata[m][n + 1] *(this->xdata[n + 1] - this->xdata[n])
                   + (this->ydata[m][n + 1] - this->ydata[m][n]);
	return -tp*this->ydata[m][n] + tp*this->ydata[m][n + 1]
            + tp*(this->_1p0 - t)*(a*(this->_1p0 - t) + b*t)
            - t*tp*(a*(this->_1p0 - t) + b*t)
            + t*(this->_1p0 - t)*(-a*tp + b*tp);
}
*/
//------------------------------------------------------------------------------
template< class value_type >
void SMHiggsPotentialSpline< value_type >::MultiPointSpline
(value_type* y, value_type* arrayOut, int numel, int nArrayToUse)
{
    switch(nArrayToUse)
    {
    case 0:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*this->yscale;//convert to GeV units
            //std::cout << "\nn = " << i;
            int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);

            arrayOut[i] = SinglePointSpline(x, 0, n)*(y[i]*y[i]*y[i]*y[i])
                          /(value_type(4.0));
        }
        break;
    case 1:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*this->yscale;
            int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
            value_type dldL = SinglePointSpline(x, 1, n);
            value_type l = SinglePointSpline(x, 0, n);
            value_type dldptp = value_type(2.0)*dldL;
            arrayOut[i] = (dldptp*(y[i]*y[i]*y[i])/(value_type(4.0)) +
            l*(y[i]*y[i]*y[i]));
        }
        break;
    case 2:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*this->yscale;
            int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
            value_type d2ldL2 = SinglePointSpline(x, 2, n);
            value_type dldL = SinglePointSpline(x, 1, n);
            value_type l = SinglePointSpline(x, 0, n);
            value_type dldptp = value_type(2.0)*dldL;
            value_type d2ldp2tp2 = -value_type(2.0)*dldL
                                   + value_type(4.0)*d2ldL2;
            arrayOut[i] = (d2ldp2tp2*(y[i]*y[i])/(value_type(4.0)) +
            dldptp*(value_type(2.0)*y[i]*y[i]) +
            value_type(3.0)*y[i]*y[i]*l);
        }
        break;
    case 3:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*this->yscale;
            int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
            value_type d3ldL3 = SinglePointSpline(x, 3, n);
            value_type d2ldL2 = SinglePointSpline(x, 2, n);
            value_type dldL = SinglePointSpline(x, 1, n);
            value_type l = SinglePointSpline(x, 0, n);
            value_type dldptp = value_type(2.0)*dldL;
            value_type d2ldp2tp2 = -value_type(2.0)*dldL
                                   + value_type(4.0)*d2ldL2;
            value_type d3ldp3tp3 = value_type(4.0)*dldL
                                   - value_type(12.0)*d2ldL2
                                   + value_type(8.0)*d3ldL3;
            arrayOut[i] = (d3ldp3tp3*(y[i])/(value_type(4.0)) +
            d2ldp2tp2*(value_type(3.0)*y[i]) +
            value_type(9.0)*y[i]*dldptp +
            value_type(6.0)*y[i]*l);
        }
        break;
    case 4:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*this->yscale;
            int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
            value_type d4ldL4 = SinglePointSpline(x, 4, n);
            value_type d3ldL3 = SinglePointSpline(x, 3, n);
            value_type d2ldL2 = SinglePointSpline(x, 2, n);
            value_type dldL = SinglePointSpline(x, 1, n);
            value_type l = SinglePointSpline(x, 0, n);
            value_type dldptp = value_type(2.0)*dldL;
            value_type d2ldp2tp2 = -value_type(2.0)*dldL
                                   + value_type(4.0)*d2ldL2;
            value_type d3ldp3tp3 = value_type(4.0)*dldL
                                   - value_type(12.0)*d2ldL2
                                   + value_type(8.0)*d3ldL3;
            value_type d4ldp4tp4 = -value_type(12.0)*dldL
                                   + value_type(44.0)*d2ldL2
                                   - value_type(48.0)*d3ldL3
                                   + value_type(16.0)*d4ldL4;
            arrayOut[i] = d4ldp4tp4/(value_type(4.0)) +
            d3ldp3tp3*(value_type(4.0)) +
            value_type(18.0)*d2ldp2tp2 +
            value_type(24.0)*dldptp +
            value_type(6.0)*l;
        }
        break;
    default:
        throw "Invalid derivative.\n";
    }
}
//------------------------------------------------------------------------------
template< class value_type >
value_type SMHiggsPotentialSpline< value_type >::operator()
(value_type y)
{
    value_type x = y*yscale;
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
	return SinglePointSpline(x, 0, n)*(y*y*y*y)/(value_type(4.0));
}
//------------------------------------------------------------------------------
template< class value_type >
void SMHiggsPotentialSpline< value_type >::operator()
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 0);
}
//------------------------------------------------------------------------------
//First derivatives:
template< class value_type >
value_type SMHiggsPotentialSpline< value_type >::d(value_type y)
{
    value_type x = y*yscale;
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    value_type dldL = SinglePointSpline(x, 1, n);
    //value_type dldL = SplineDerivative(x, 0, n);
    value_type l = SinglePointSpline(x, 0, n);
    value_type dldptp = value_type(2.0)*dldL;
	return dldptp*(y*y*y)/(value_type(4.0)) + l*(y*y*y);
}
//------------------------------------------------------------------------------
template< class value_type >
void SMHiggsPotentialSpline< value_type >::d
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 1);
}
//------------------------------------------------------------------------------
//Second derivatives:
template< class value_type >
value_type SMHiggsPotentialSpline< value_type >::d2(value_type y)
{
    value_type x = y*yscale;
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    value_type d2ldL2 = SinglePointSpline(x, 2, n);
    value_type dldL = SinglePointSpline(x, 1, n);
    value_type l = SinglePointSpline(x, 0, n);
    value_type dldptp = value_type(2.0)*dldL;
    value_type d2ldp2tp2 = -value_type(2.0)*dldL + value_type(4.0)*d2ldL2;
	return d2ldp2tp2*(y*y)/(value_type(4.0)) + dldptp*(value_type(2.0)*y*y)
            + value_type(3.0)*y*y*l;
}
//------------------------------------------------------------------------------
template< class value_type >
void SMHiggsPotentialSpline< value_type >::d2
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 2);
}
//------------------------------------------------------------------------------
//Third derivatives:
template< class value_type >
value_type SMHiggsPotentialSpline< value_type >::d3(value_type y)
{
    value_type x = y*yscale;
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    value_type d3ldL3 = SinglePointSpline(x, 3, n);
    value_type d2ldL2 = SinglePointSpline(x, 2, n);
    value_type dldL = SinglePointSpline(x, 1, n);
    value_type l = SinglePointSpline(x, 0, n);
    value_type dldptp = value_type(2.0)*dldL;
    value_type d2ldp2tp2 = -value_type(2.0)*dldL + value_type(4.0)*d2ldL2;
    value_type d3ldp3tp3 = value_type(4.0)*dldL - value_type(12.0)*d2ldL2
                           + value_type(8.0)*d3ldL3;
	return d3ldp3tp3*(y)/(value_type(4.0)) + d2ldp2tp2*(value_type(3.0)*y)
           + value_type(9.0)*y*dldptp + value_type(6.0)*y*l;
}
//------------------------------------------------------------------------------
template< class value_type >
void SMHiggsPotentialSpline< value_type >::d3
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 3);
}
//------------------------------------------------------------------------------
//Fourth derivatives:
template< class value_type >
value_type SMHiggsPotentialSpline< value_type >::d4(value_type y)
{
    value_type x = y*yscale;
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    value_type d4ldL4 = SinglePointSpline(x, 4, n);
    value_type d3ldL3 = SinglePointSpline(x, 3, n);
    value_type d2ldL2 = SinglePointSpline(x, 2, n);
    value_type dldL = SinglePointSpline(x, 1, n);
    value_type l = SinglePointSpline(x, 0, n);
    value_type dldptp = value_type(2.0)*dldL;
    value_type d2ldp2tp2 = -value_type(2.0)*dldL + value_type(4.0)*d2ldL2;
    value_type d3ldp3tp3 = value_type(4.0)*dldL - value_type(12.0)*d2ldL2
               + value_type(8.0)*d3ldL3;
    value_type d4ldp4tp4 = -value_type(12.0)*dldL + value_type(44.0)*d2ldL2
               - value_type(48.0)*d3ldL3 + value_type(16.0)*d4ldL4;
	return d4ldp4tp4/(value_type(4.0)) + d3ldp3tp3*(value_type(4.0))
               + value_type(18.0)*d2ldp2tp2 + value_type(24.0)*dldptp
               + value_type(6.0)*l;
}
//------------------------------------------------------------------------------
template< class value_type >
void SMHiggsPotentialSpline< value_type >::d4
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 4);
}
//------------------------------------------------------------------------------
template< class value_type >
int SplinePotential< value_type >::ChooseSpline
(value_type& x, value_type* tdata, int Nupper)
{
	//We would like to optimise this function as much as possible, as it will
	//be called each time the
	//potential is called, which is potentially (no pun intended!) a lot.

	//The list xdata is assumed to be ordered in ascending order, otherwise
	//this won't work.
	const value_type zero = value_type(0.0);
	value_type logx = ( x != zero ? log(x*x/(this->M*this->M)) : zero);
	//std::cout.precision(50);
	//std::cout << "\nlogx = " << logx;
	int lower = 0;
	int upper = Nupper;
	//std::cout << "\nm = 1";
	//std::cout << "\nlogx = " << logx;
	//std::cout << "\ntdata[lower] = " << tdata[lower];
	//std::cout << "\ntdata[upper] = " << tdata[upper];
	if(logx >= tdata[upper] || logx <= tdata[lower] || x == zero)
	{
	    //Might be able to recover if logx is negative - indicates that we
	    //should extrapolate to small values below the
	    //top quark mass. In other words, just use the first spline:
	    /*
	    if(x == zero || logx == zero)
s        {
            //std::cout << "\nx == 0";
            return 0;
        }
        else
        {
            if(logx > zero && logx <= tdata[lower])
            {
                //std::cout << "\nm = 2";
                return 0;
            }
            else
            {
                //std::cout << "\nm = 3";
                std::stringstream ss;
                ss << logx << "\n x = " << x;
                ss << "\nx is out of range for interpolation.\n";
                std::string erStr = ss.str();
                //throw erStr;
                throw erStr.c_str();
            }
        }
        */
        if(logx > tdata[upper])
        {
            //std::cout << "\nm = 3";
            std::stringstream ss;
            ss << logx << "\n x = " << x;
            ss << "\nx is out of range for interpolation.\n";
            std::string erStr = ss.str();
            //throw erStr;
            throw erStr.c_str();
        }
        else
        {
            return 0;
        }
	}
	//std::cout << "\nm = 4";
	/*
	std::cout << "\nr = " << floor(value_type(upper - lower)*(logx - tdata[lower])
            / (tdata[upper] - tdata[lower]));
    std::cout << "\nr_n = " << convert_type<value_type,int>(floor(value_type(upper - lower)*(logx - tdata[lower])
            / (tdata[upper] - tdata[lower])));
            */
	int n = lower + convert_type<value_type,int>(floor(value_type(upper - lower)*(logx - tdata[lower])
            / (tdata[upper] - tdata[lower])));
	int m = lower + convert_type<value_type,int>(ceil(value_type(upper - lower)*(logx - tdata[lower])
            / (tdata[upper] - tdata[lower])));
	int nSplineToUse = n;
	if(n != m)//If n == m, then we hit one of our data points by chance.
    {
        while (!(tdata[m] >= logx && tdata[n] <= logx))
        {
            if (logx > tdata[m])
            {
                //Undershot:
                lower = m;
            }
            else if (logx < tdata[n])
            {
                //Overshot:
                upper = n;
            }
            else
            {
                //Must have   tdata[n] <= logx <= tdata[m], so we found the right
                //spline:
                nSplineToUse = n;
                break;
            }
            //Otherwise, adjust m and n  and keep searching:
            n = lower + convert_type<value_type,int>(floor(value_type(upper - lower)*(logx - tdata[lower])
                / (tdata[upper] - tdata[lower])));
            m = lower + convert_type<value_type,int>(ceil(value_type(upper - lower)*(logx - tdata[lower])
                / (tdata[upper] - tdata[lower])));
            nSplineToUse = n;
        }
    }
	//std::cout << "\nm = 5";
	return nSplineToUse;
}
//------------------------------------------------------------------------------
#endif // USING_SM_POTENTIAL
//==============================================================================


//==============================================================================
//Potential used to study new physics operators.
#ifdef USING_POLY_POTENTIAL
//------------------------------------------------------------------------------
template< class value_type >
value_type polyPotential< value_type >::operator()(value_type y)
{
    value_type _1o2 = value_type(1.0)/value_type(2.0);
    value_type _1o4 = value_type(1.0)/value_type(4.0);
    value_type _1o6 = value_type(1.0)/value_type(6.0);
    value_type _1o8 = value_type(1.0)/value_type(8.0);
    return _1o2*m2*y*y + _1o4*l4*y*y*y*y + _1o6*l6*y*y*y*y*y*y/(M6*M6)
           + _1o8*l8*y*y*y*y*y*y*y*y/(M8*M8*M8*M8);
}
//------------------------------------------------------------------------------
template< class value_type >
void polyPotential< value_type >::operator()
(value_type* y, value_type* arrayOut, int numel)
{
    value_type _1o2 = value_type(1.0)/value_type(2.0);
    value_type _1o4 = value_type(1.0)/value_type(4.0);
    value_type _1o6 = value_type(1.0)/value_type(6.0);
    value_type _1o8 = value_type(1.0)/value_type(8.0);
	for(int i = 0; i < numel;i++)
    {
        arrayOut[i] = _1o2*m2*y[i]*y[i] + _1o4*l4*y[i]*y[i]*y[i]*y[i]
                      + _1o6*l6*y[i]*y[i]*y[i]*y[i]*y[i]*y[i]/(M6*M6)
                      + _1o8*l8*y[i]*y[i]*y[i]*y[i]*y[i]*y[i]*y[i]*y[i]
                      /(M8*M8*M8*M8);
    }
}
//------------------------------------------------------------------------------
//First derivatives:
template< class value_type >
value_type polyPotential< value_type >::d(value_type y)
{
	return m2*y + l4*y*y*y + l6*y*y*y*y*y/(M6*M6)
                  + l8*y*y*y*y*y*y*y/(M8*M8*M8*M8);
}
//------------------------------------------------------------------------------
template< class value_type >
void polyPotential< value_type >::d
(value_type* y, value_type* arrayOut, int numel)
{
	for(int i = 0; i < numel;i++)
    {
        arrayOut[i] = m2*y[i] + l4*y[i]*y[i]*y[i]
                      + l6*y[i]*y[i]*y[i]*y[i]*y[i]/(M6*M6)
                      + l8*y[i]*y[i]*y[i]*y[i]*y[i]*y[i]*y[i]/(M8*M8*M8*M8);
    }
}
//------------------------------------------------------------------------------
//Second derivatives:
template< class value_type >
value_type polyPotential< value_type >::d2(value_type y)
{
    value_type _3 = value_type(3.0);
    value_type _5 = value_type(5.0);
    value_type _7 = value_type(7.0);
    return m2 + _3*l4*y*y + _5*l6*y*y*y*y/(M6*M6)
           + _7*l8*y*y*y*y*y*y/(M8*M8*M8*M8);
}
//------------------------------------------------------------------------------
template< class value_type >
void polyPotential< value_type >::d2
(value_type* y, value_type* arrayOut, int numel)
{
	value_type _3 = value_type(3.0);
    value_type _5 = value_type(5.0);
    value_type _7 = value_type(7.0);
    for(int i = 0;i < numel;i++)
    {
        arrayOut[i] = m2 + _3*l4*y[i]*y[i] + _5*l6*y[i]*y[i]*y[i]*y[i]/(M6*M6)
                      + _7*l8*y[i]*y[i]*y[i]*y[i]*y[i]*y[i]/(M8*M8*M8*M8);
    }
}
//------------------------------------------------------------------------------
//Third derivatives:
template< class value_type >
value_type polyPotential< value_type >::d3(value_type y)
{
    value_type _6 = value_type(6.0);
    value_type _20 = value_type(20.0);
    value_type _42 = value_type(42.0);
    return _6*l4*y + _20*l6*y*y*y/(M6*M6) + _42*l8*y*y*y*y*y/(M8*M8*M8*M8);
}
//------------------------------------------------------------------------------
template< class value_type >
void polyPotential< value_type >::d3
(value_type* y, value_type* arrayOut, int numel)
{
	value_type _6 = value_type(6.0);
    value_type _20 = value_type(20.0);
    value_type _42 = value_type(42.0);
    for(int i = 0;i < numel;i++)
    {
        arrayOut[i] = _6*l4*y[i] + _20*l6*y[i]*y[i]*y[i]/(M6*M6)
                      + _42*l8*y[i]*y[i]*y[i]*y[i]*y[i]/(M8*M8*M8*M8);
    }
}
//------------------------------------------------------------------------------
//Fourth derivatives:
template< class value_type >
value_type polyPotential< value_type >::d4(value_type y)
{
    value_type _6 = value_type(6.0);
    value_type _60 = value_type(60.0);
    value_type _210 = value_type(210.0);
    return _6*l4 + _60*l6*y*y/(M6*M6) + _210*l8*y*y*y*y/(M8*M8*M8*M8);
}
//------------------------------------------------------------------------------
template< class value_type >
void polyPotential< value_type >::d4
(value_type* y, value_type* arrayOut, int numel)
{
	value_type _6 = value_type(6.0);
    value_type _60 = value_type(60.0);
    value_type _210 = value_type(210.0);
    for(int i = 0;i < numel;i++)
    {
        arrayOut[i] =  _6*l4 + _60*l6*y[i]*y[i]/(M6*M6)
                      + _210*l8*y[i]*y[i]*y[i]*y[i]/(M8*M8*M8*M8);
    }
}
//------------------------------------------------------------------------------
#endif//USING_POLY_POTENTIAL
//==============================================================================
//computeSharedData - computes data needed by the other functions, A,B,Ap,...
//so that they don't need to repeat each other's work. This is done once,
//before any of those functions are called:
/*
template<class value_type>
void dtcldtxi<value_type>::computeSharedData(beta_data<value_type>& x,
                                             int sharing_level)
{
    this->phicl_shared = this->mt*exp(x.g[8]/this->_2p0);
    this->mu_shared = this->mt*exp(x.g[9]/this->_2p0);
    this->dtdtcl_shared = this->dtdtcl_fun(x);
    if(sharing_level > 0)
    {
        this->d2tdtcl2_shared = this->dtdtcl_fun._1st_derivative(x);;
        this->dtcldtxi_shared = this->f(x);
    }
    else if (sharing_level > 1)
    {
        this->d3tdtcl3_shared = this->dtdtcl_fun._2nd_derivative(x);
        this->d2tcldtxi2_shared = this->d(x);
    }
    else if(sharing_level > 2)
    {
        this->d4tdtcl4_shared = this->dtdtcl_fun._3rd_derivative(x);
        this->d3tcldtxi3_shared = this->d2(x);
    }
}
*/
//------------------------------------------------------------------------------
//Member functions for dtcldtxi:
//NB - this code is NOT meant to be debugged by hand - it is generated by
//mathematica directly from the derivatives of the relevant expressions. If
//debugging is necessary, it should be regenerated, and processed using the
//appropriate post-processing programs (namely "value_type_parser.exe")
//A - dtcldtxi
/*
template<class value_type>
value_type dtcldtxi<value_type>::A(beta_data<value_type>& x)
{
    value_type tcl = x.g[8];
    value_type xi = x.g[6];
    value_type phicl = this->mt*exp(tcl/this->_2p0);
    return (exp((((this->_1p0)/(this->_2p0))*(x.txi - tcl)))
            *(this->_1p0 - (this->h*this->h*(phicl*phicl*xi))));
}
*/
//------------------------------------------------------------------------------
//B - dtcldtxi
/*
template<class value_type>
value_type dtcldtxi<value_type>::B(beta_data<value_type>& x)
{
    value_type tcl = x.g[8];
    value_type xi = x.g[6];
    value_type xip = x.dgdt[6];
    value_type phicl = this->phicl_shared;
    value_type dtdtcl = this->dtdtcl_shared;
    return sqrt(((this->_1p0 + (this->_6p0*(this->h*this->h*(phicl*phicl
            *(xi + (dtdtcl*xip))*(xi + (dtdtcl*xip))))))
                 - (this->h*this->h*(phicl*phicl*xi))));
}
//------------------------------------------------------------------------------
//Ap - dtcldtxi
template<class value_type>
value_type dtcldtxi<value_type>::Ap(beta_data<value_type>& x)
{
    value_type tcl = x.g[8];
    value_type xi = x.g[6];
    value_type xip = x.dgdt[6];
    value_type phicl = this->phicl_shared;
    value_type dtcldtxi = this->dtcldtxi_shared;
    value_type arg = (((this->_1p0)/(this->_2p0))*(x.txi - tcl));
    value_type exparg = exp(arg);
    return (((this->_m1p0)/(this->_2p0))*(exparg*(this->_m1p0 + (dtcldtxi
            + (((this->_1p0 + dtcldtxi)*(this->h*this->h*(phicl*phicl*xi)))
            + (this->_2p0*(this->h*this->h*(phicl*phicl*xip))))))));
}
//------------------------------------------------------------------------------
//Bp - dtcldtxi
template<class value_type>
value_type dtcldtxi<value_type>::Bp(beta_data<value_type>& x)
{
    value_type tcl = x.g[8];
    value_type xi = x.g[6];
    value_type xip = x.dgdt[6];
    value_type xipp = x.d2gdt2[6];
    value_type phicl = this->phicl_shared;
    value_type dtdtcl = this->dtdtcl_shared;
    value_type d2tdtcl2 = this->d2tdtcl2_shared;
    value_type dtcldtxi = this->dtcldtxi_shared;

    return (((this->_1p0)/(this->_2p0))*(this->h*this->h*(phicl*phicl
           *(((this->_1p0)/(sqrt(((this->_1p0 + (this->_6p0
           *(this->h*this->h*(phicl*phicl*(xi + (dtdtcl*xip))*(xi + (dtdtcl*xip))))))
           - (this->h*this->h*(phicl*phicl*xi))))))*((((this->_6p0*(dtcldtxi*(xi
           + (dtdtcl*xip))*(xi + (dtdtcl*xip)))) + (this->_12p0*((xi
           + (dtdtcl*xip))*(xip + (dtcldtxi*((d2tdtcl2*xip)
           + (dtdtcl*dtdtcl*xipp))))))) - xip) - (dtcldtxi*xi))))));
}
//------------------------------------------------------------------------------
//App - dtcldtxi
template<class value_type>
value_type dtcldtxi<value_type>::App(beta_data<value_type>& x)
{
    value_type tcl = x.g[8];
    value_type xi = x.g[6];
    value_type xip = x.dgdt[6];
    value_type xipp = x.d2gdt2[6];
    value_type phicl = this->phicl_shared;
    value_type dtcldtxi = this->dtcldtxi_shared;
    value_type d2tcldtxi2 = this->d2tcldtxi2_shared;
    value_type arg = (((this->_1p0)/(this->_2p0))*(x.txi - tcl));
    value_type exparg = exp(arg);

    return (((this->_m1p0)/(this->_4p0))*(exparg*((this->_m1p0 + ((this->_2p0*d2tcldtxi2)
            + (((this->_1p0 + ((this->_2p0*d2tcldtxi2)
            + ((this->_2p0*dtcldtxi) + dtcldtxi*dtcldtxi)))*(h*h*(phicl
            *phicl*xi))) + ((this->_4p0*(h*h*(phicl*phicl*xip)))
            + ((dtcldtxi*(this->_2p0
            + (this->_4p0*(h*h*(phicl*phicl*xip)))))
            + (this->_4p0*(h*h*(phicl*phicl*xipp))))))))
            - dtcldtxi*dtcldtxi)));
}
//------------------------------------------------------------------------------
//Bpp - dtcldtxi
template<class value_type>
value_type dtcldtxi<value_type>::Bpp(beta_data<value_type>& x)
{
    value_type tcl = x.g[8];
    value_type xi = x.g[6];
    value_type xip = x.dgdt[6];
    value_type xipp = x.d2gdt2[6];
    value_type xippp = x.d3dgt3[6];
    value_type phicl = this->phicl_shared;
    value_type dtdtcl = this->dtdtcl_shared;
    value_type d2tdtcl2 = this->d2tdtcl2_shared;
    value_type d3tdtcl3 = this->d3tdtcl3_shared;
    value_type dtcldtxi = this->dtcldtxi_shared;
    value_type d2tcldtxi2 = this->d2tcldtxi2_shared;
    value_type arg = ((this->_1p0 + (this->_6p0*(this->h*this->h*(phicl*phicl*(xi
                    + (dtdtcl*xip))*(xi + (dtdtcl*xip))))))
                    - (this->h*this->h*(phicl*phicl*xi)));
    value_type exponent = ((this->_3p0)/(this->_2p0));

    return (((this->_1p0)/(this->_4p0))*(((this->_1p0)
            /(pow(arg,exponent)))*((this->_2p0*(this->h*this->h*(phicl*phicl
            *(((this->_1p0 + (this->_6p0*(this->h*this->h*(phicl*phicl*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip))))))
            - (this->h*this->h*(phicl*phicl*xi)))*(((((this->_m2p0*(dtcldtxi*xip))
            + ((this->_6p0*(d2tcldtxi2*(xi + (dtdtcl*xip))*(xi
            + (dtdtcl*xip)))) + ((this->_6p0*(dtcldtxi*dtcldtxi*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip))))
            + ((this->_24p0*(dtcldtxi*((xi + (dtdtcl*xip))*(xip
            + (dtcldtxi*((d2tdtcl2*xip) + (dtdtcl*dtdtcl*xipp)))))))
            + ((this->_12p0*(xip + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp))))*(xip + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp))))) + (this->_12p0*((xi
            + (dtdtcl*xip))*((((d2tcldtxi2*d2tdtcl2)
            + (d3tdtcl3*dtcldtxi*dtcldtxi))*xip)
            + (xipp + (dtdtcl*((d2tcldtxi2*(dtdtcl*xipp))
            + (dtcldtxi*dtcldtxi*((this->_3p0*(d2tdtcl2*xipp))
            + (dtdtcl*dtdtcl*xippp)))))))))))))) - xipp)
            - (dtcldtxi*dtcldtxi*xi)) - (d2tcldtxi2*xi))))))
            - (this->h*this->h*this->h*this->h*(phicl*phicl*phicl*phicl*((dtcldtxi*xi)
            + (xip + ((this->_m6p0*(dtcldtxi*(xi + (dtdtcl*xip))*(xi
            + (dtdtcl*xip)))) + (this->_m12p0*((xi
            + (dtdtcl*xip))*(xip + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp)))))))))*((dtcldtxi*xi) + (xip
            + ((this->_m6p0*(dtcldtxi*(xi + (dtdtcl*xip))*(xi
            + (dtdtcl*xip)))) + (this->_m12p0*((xi
            + (dtdtcl*xip))*(xip + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp))))))))))))));
}
//------------------------------------------------------------------------------
//Appp - dtcldtxi
template<class value_type>
value_type dtcldtxi<value_type>::Appp(beta_data<value_type>& x)
{
    value_type tcl = x.g[8];
    value_type xi = x.g[6];
    value_type xip = x.dgdt[6];
    value_type xipp = x.d2gdt2[6];
    value_type xippp = x.d3gdt3[6];
    value_type phicl = this->phicl_shared;
    value_type dtcldtxi = this->dtcldtxi_shared;
    value_type d2tcldtxi2 = this->d2tcldtxi2_shared;
    value_type d3tcldtxi3 = this->d3tcldtxi3_shared;
    value_type arg = (((this->_1p0)/(this->_2p0))*(x.txi - tcl));
    value_type exparg = exp(arg);

    return (((this->_m1p0)/(this->_8p0))*(exparg*(this->_m1p0 + ((this->_6p0*d2tcldtxi2)
            + ((this->_4p0*d3tcldtxi3) + (dtcldtxi*dtcldtxi*dtcldtxi
            + (((this->_1p0 + ((this->_6p0*d2tcldtxi2)
            + ((this->_4p0*d3tcldtxi3) + (((this->_3p0
            + (this->_6p0*d2tcldtxi2))*dtcldtxi)
            + ((this->_3p0*dtcldtxi*dtcldtxi)
            + dtcldtxi*dtcldtxi*dtcldtxi)))))*(this->h*this->h*(phicl*phicl*xi)))
            + ((this->_6p0*(this->h*this->h*(phicl*phicl*xip)))
            + ((this->_12p0*(d2tcldtxi2*(this->h*this->h*(phicl*phicl*xip))))
            + ((dtcldtxi*dtcldtxi*(this->_m3p0
            + (this->_6p0*(this->h*this->h*(phicl*phicl*xip)))))
            + ((this->_12p0*(this->h*this->h*(phicl*phicl*xipp)))
            + ((this->_3p0*(dtcldtxi*(this->_1p0
            + ((this->_m2p0*d2tcldtxi2) + ((this->_4p0*(this->h*this->h*(phicl*phicl
            *xip))) + (this->_4p0*(this->h*this->h*(phicl*phicl*xipp))))))))
            + (this->_8p0*(this->h*this->h*(phicl*phicl*xippp)))))))))))))));
}
//------------------------------------------------------------------------------
//Bppp - dtcldtxi
template<class value_type>
value_type dtcldtxi<value_type>::Bppp(beta_data<value_type>& x)
{
    value_type tcl = x.g[8];
    value_type xi = x.g[6];
    value_type xip = x.dgdt[6];
    value_type xipp = x.d2gdt2[6];
    value_type xippp = x.d3dgt3[6];
    value_type xipppp = x.d4dgt4[6];
    value_type phicl = this->phicl_shared;
    value_type dtdtcl = this->dtdtcl_shared;
    value_type d2tdtcl2 = this->d2tdtcl2_shared;
    value_type d3tdtcl3 = this->d3tdtcl3_shared;
    value_type d4tdtcl4 = this->d3tdtcl3_shared;
    value_type dtcldtxi = this->dtcldtxi_shared;
    value_type d2tcldtxi2 = this->d2tcldtxi2_shared;
    value_type d3tcldtxi3 = this->d2tcldtxi2_shared;
    value_type arg = ((this->_1p0
            + (this->_6p0*(this->h*this->h*(phicl*phicl*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip))))))
            - (this->h*this->h*(phicl*phicl*xi)));
    value_type exponent = ((this->_5p0)/(this->_2p0));

    return (((this->_1p0)/(this->_8p0))*(((this->_1p0)/(pow(arg,exponent)))
            *((this->_m3p0*(this->h*this->h*this->h*this->h*this->h*this->h
            *(phicl*phicl*phicl*phicl*phicl*phicl*((dtcldtxi*xi)
            + (xip + ((this->_m6p0*(dtcldtxi*(xi + (dtdtcl*xip))*(xi
            + (dtdtcl*xip)))) + (this->_m12p0*((xi + (dtdtcl*xip))*(xip
            + (dtcldtxi*((d2tdtcl2*xip) + (dtdtcl*dtdtcl*xipp)))))))))
            *((dtcldtxi*xi) + (xip + ((this->_m6p0*(dtcldtxi*(xi + (dtdtcl
            *xip))*(xi + (dtdtcl*xip))))
            + (this->_m12p0*((xi + (dtdtcl*xip))*(xip
            + (dtcldtxi*((d2tdtcl2*xip) + (dtdtcl*dtdtcl*xipp)))))))))
            *((dtcldtxi*xi) + (xip + ((this->_m6p0*(dtcldtxi*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip)))) + (this->_m12p0*((xi
            + (dtdtcl*xip))*(xip + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp))))))))))))
            + ((this->_m6p0*(this->h*this->h*this->h*this->h
            *(phicl*phicl*phicl*phicl*(((this->_1p0
            + (this->_6p0*(this->h*this->h*(phicl*phicl*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip))))))
            - (this->h*this->h*(phicl*phicl*xi)))*(((((this->_6p0
            *(dtcldtxi*(xi + (dtdtcl*xip))*(xi + (dtdtcl*xip))))
            + (this->_12p0*((xi + (dtdtcl*xip))*(xip + (dtcldtxi
            *((d2tdtcl2*xip) + (dtdtcl*dtdtcl*xipp))))))) - xip)
            - (dtcldtxi*xi))*(((((this->_m2p0*(dtcldtxi*xip))
            + ((this->_6p0*(d2tcldtxi2*(xi + (dtdtcl*xip))*(xi
            + (dtdtcl*xip)))) + ((this->_6p0*(dtcldtxi*dtcldtxi*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip))))
            + ((this->_24p0*(dtcldtxi*((xi + (dtdtcl*xip))*(xip
            + (dtcldtxi*((d2tdtcl2*xip) + (dtdtcl*dtdtcl*xipp)))))))
            + ((this->_12p0*(xip + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp))))*(xip + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp))))) + (this->_12p0*((xi + (dtdtcl*xip))
            *((((d2tcldtxi2*d2tdtcl2) + (d3tdtcl3*dtcldtxi*dtcldtxi))*xip)
            + (xipp + (dtdtcl*((d2tcldtxi2*(dtdtcl*xipp))
            + (dtcldtxi*dtcldtxi*((this->_3p0*(d2tdtcl2*xipp))
            + (dtdtcl*dtdtcl*xippp)))))))))))))) - xipp)
            - (dtcldtxi*dtcldtxi*xi)) - (d2tcldtxi2*xi)))))))
            + (this->_4p0*(this->h*this->h*(phicl*phicl*(((this->_1p0
            + (this->_6p0*(this->h*this->h*(phicl*phicl*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip))))))
            - (this->h*this->h*(phicl*phicl*xi)))*((this->_1p0
            + (this->_6p0*(this->h*this->h*(phicl*phicl*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip))))))
            - (this->h*this->h*(phicl*phicl*xi)))*(((((this->_m3p0
            *(d2tcldtxi2*(dtcldtxi*xi)))
            + ((this->_m3p0*(d2tcldtxi2*xip))
            + ((this->_m3p0*(dtcldtxi*dtcldtxi*xip))
            + ((this->_6p0*(d3tcldtxi3*(xi + (dtdtcl*xip))*(xi
            + (dtdtcl*xip)))) + ((this->_18p0*(d2tcldtxi2*(dtcldtxi*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip)))))
            + ((this->_6p0*(dtcldtxi*dtcldtxi*dtcldtxi*(xi
            + (dtdtcl*xip))*(xi + (dtdtcl*xip))))
            + ((this->_m3p0*(dtcldtxi*xipp))
            + ((this->_36p0*(d2tcldtxi2*((xi + (dtdtcl*xip))*(xip
            + (dtcldtxi*((d2tdtcl2*xip) + (dtdtcl*dtdtcl*xipp)))))))
            + ((this->_36p0*(dtcldtxi*dtcldtxi*((xi + (dtdtcl*xip))*(xip
            + (dtcldtxi*((d2tdtcl2*xip) + (dtdtcl*dtdtcl*xipp)))))))
            + ((this->_36p0*(dtcldtxi*(xip + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp))))*(xip + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp)))))) + ((this->_36p0*(dtcldtxi*((xi
            + (dtdtcl*xip))*((((d2tcldtxi2*d2tdtcl2)
            + (d3tdtcl3*dtcldtxi*dtcldtxi))*xip) + (xipp
            + (dtdtcl*((d2tcldtxi2*(dtdtcl*xipp))
            + (dtcldtxi*dtcldtxi*((this->_3p0*(d2tdtcl2*xipp))
            + (dtdtcl*dtdtcl*xippp)))))))))) + ((this->_36p0*((xip
            + (dtcldtxi*((d2tdtcl2*xip)
            + (dtdtcl*dtdtcl*xipp))))*((((d2tcldtxi2*d2tdtcl2)
            + (d3tdtcl3*dtcldtxi*dtcldtxi))*xip)
            + (xipp + (dtdtcl*((d2tcldtxi2*(dtdtcl*xipp))
            + (dtcldtxi*dtcldtxi*((this->_3p0*(d2tdtcl2*xipp))
            + (dtdtcl*dtdtcl*xippp))))))))) + (this->_12p0*((xi
            + (dtdtcl*xip))*((d2tdtcl2*(d3tcldtxi3*xip))
            + ((d3tcldtxi3*(dtdtcl*dtdtcl*xipp)) + (xippp
            + ((this->_3p0*(d2tcldtxi2*(dtcldtxi*((d3tdtcl3*xip)
            + ((this->_3p0*(d2tdtcl2*(dtdtcl*xipp)))
            + (dtdtcl*dtdtcl*dtdtcl*xippp))))))
            + (dtcldtxi*dtcldtxi*dtcldtxi*((d4tdtcl4*xip)
            + ((this->_3p0*(d2tdtcl2*d2tdtcl2*xipp))
            + ((this->_4p0*(d3tdtcl3*(dtdtcl*xipp)))
            + ((this->_6p0*(d2tdtcl2*(dtdtcl*dtdtcl*xippp)))
            + (dtdtcl*dtdtcl*dtdtcl*dtdtcl*xipppp))))))))))))))))))))))))
            - xippp) - (dtcldtxi*dtcldtxi*dtcldtxi*xi))
            - (d3tcldtxi3*xi))))))))));
}
*/
//------------------------------------------------------------------------------
//dtdtcl_1loop_no_logs member functions, used to compute derivatives of
//A(x)/B(x).
//------------------------------------------------------------------------------
//computeSharedData - used to compute functions before A,B, etc... are
//evaluated, to prevent unnecessary repetition of work.
/*
template<class value_type>
void dtdtcl_1loop_no_logs<value_type>::computeSharedData
(beta_data_1_loop<value_type>& x,int sharing_level)
{
    this->phicl_shared = this->mt*exp(x.g[8]/this->_2p0);
    this->mu_shared = this->mt*exp(x.g[9]/this->_2p0);
    this->log_shared.assign(9,this->_0p0);
    for(int i = 0;i < 9;i++)
    {
        this->log_shared[i] = log(x.Mi2[i]/(this->mt*this->mt))
                                - x.g[9];
    }
    if(sharing_level > 0)
    {
        this->dtdtcl_shared = this->f(x);
    }
    else if (sharing_level > 1)
    {
        this->d2dtcl2_shared = this->d(x);
    }
    else if(sharing_level > 2)
    {
        this->d3dtcl3_shared = this->d2(x);
    }
}
*/
//------------------------------------------------------------------------------
#ifdef USING_SM_POTENTIAL
//computeMi2:
//
//Computes Mi^2(\phi) = k_i*phi^2(t) - k'_i + theta_i*R
//Where k_i, k'_i, theta_i depend on the particle contributing to
//the given loop correction, R is the Ricci scalar, and phi = Z(t)*phicl
//is the renormalised field (phicl is the field evaluated at the electroweak
//scale).
//
//This function takes as inputs vectors in which to store Mi^2 (Mi2),
// \frac{\partial Mi^2(phi)}{\partial t} (pMi2pt), and log(Mi^2(phi)/mu^2(t))
// (logfactor), as well as arrays of coefficients k_i (ki), k'_i (kpi) and
//theta_i (thetai), and the t derivatives of these coefficients.
//It also needs Z and its t derivative, Zp, and log(mu^2/mt^2) (logmu2).
//
template<class value_type>
void dtdtcl_1loop<value_type>::computeMi2
(std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
 std::vector<value_type>& pMi2pt,
 const std::vector<value_type>& ki,const std::vector<value_type>& kpi,
 const std::vector<value_type>& thetai,
 const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
 const std::vector<value_type>& dthetaidt,
 const value_type& phicl,const value_type& Z,const value_type& Zp,
 const value_type& logmu2)
{
    //Scalar field at chosen renormalisation scale:
    value_type phi = phicl*Z;
    //std::cout << "\nphicl = " << phicl;
    //std::cout << "\nZ = " << Z;
    //std::cout << "\nZp = " << Zp;
    //Anomalous dimension:
    //value_type gamma = -Zp;
    //Pre-allocate output vectors to prevent any going out of bounds:
    Mi2.assign(this->nLoops,this->_0p0);
    pMi2pt.assign(this->nLoops,this->_0p0);
    logfactor.assign(this->nLoops,this->_0p0);
    for(int i = 0;i < this->nLoops;i++)
    {
        //std::cout << "\nk[" << i << "] = " << ki[i];
        //std::cout << "\nk'[" << i << "] = " << kpi[i];
        //std::cout << "\ntheta[" << i << "] = " << thetai[i];
        //Mi^2(phi):
        Mi2[i] = ki[i]*phi*phi - kpi[i] + thetai[i]*this->R;
        //Log(Mi^2(phi)/mt^2)
        logfactor[i] = log(abs(Mi2[i])/(this->mt*this->mt)) - logmu2;
        // \frac{\partial Mi^2(\phi)}{\partial t} at constant phicl
        pMi2pt[i] = (dkidt[i] + this->_2p0*ki[i]*Zp/Z)*Z*Z*phicl*phicl
                                - dkpidt[i] + this->R*dthetaidt[i];
        //std::cout << "\nMi2[" << i << "] = " << Mi2[i];
        //std::cout << "\nlogfactor[" << i << "] = " << logfactor[i];
        //std::cout << "\npMi2pt[" << i << "] = " << pMi2pt[i];
    }
}
//Compute the derivatives of Mi2 with respect to tcl:
//Note that the coupling derivatives here are supplied as derivatives
//with respect to t, NOT tcl. The function converts them to tcl derivatives.
//For reference, t = log(mu^2/mt^2) and tcl = log(phicl^2/mt^2).
//
//Note that this function
// internally calls computeMi2
//in order to compute Mi2, pMi2pt, logfactor, so it is not
//necessary to compute those separately.
template<class value_type>
void dtdtcl_1loop<value_type>::computeMi2p
(std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
 std::vector<value_type>& Mi2p,
 std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
 const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
 const std::vector<value_type>& thetai,
 const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
 const std::vector<value_type>& dthetaidt,
 const std::vector<value_type>& d2kidt2,const std::vector<value_type>& d2kpidt2,
 const std::vector<value_type>& d2thetaidt2,
 const value_type& phicl,const value_type& Z,const value_type& Zp,
 const value_type& Zpp,const value_type& dtdtcl,const value_type& logmu2)
{
    //Compute Mi2, as most functions requiring Mi2p also need Mi2:
    this->computeMi2(Mi2,logfactor,pMi2pt,kiarray,kpi,thetai,dkidt,dkpidt,
                     dthetaidt,phicl,Z,Zp,logmu2);
    //Pre-allocate to prevent any going out of bounds:
    Mi2p.assign(this->nLoops,this->_0p0);
    pMi2ptp.assign(this->nLoops,this->_0p0);
    value_type ki,kip,kipp,thetaip,thetaipp,kpip,kpipp;
    for(int i = 0;i < Mi2.size();i++)
    {
        ki = kiarray[i];
        kip = dkidt[i];
        kipp = d2kidt2[i];
        thetaip = dthetaidt[i];
        thetaipp = d2thetaidt2[i];
        kpip = dkpidt[i];
        kpipp = d2kpidt2[i];

        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter these:
        Mi2p[i] = ((dtdtcl*(((this->R*thetaip) + (kip*(phicl*phicl*Z*Z)))
                - kpip)) + (ki*(phicl*phicl*(Z*(Z +
                (this->_2p0*(dtdtcl*Zp)))))));
        pMi2ptp[i] = (((kip + (dtdtcl*kipp))*(phicl*phicl*Z*Z))
                   + ((dtdtcl*(((this->R*thetaipp)+(this->_2p0*(ki*(phicl*phicl*Zp
                   *Zp)))) - kpipp)) + (this->_2p0*(phicl*phicl*(Z*((this->_2p0
                   *(dtdtcl*(kip*Zp)))
                    + (ki*(Zp + (dtdtcl*Zpp)))))))));
    }
}
//Same as computeMi2p, but for second derivatives wrt tcl.
//Internally calls computeMi2p
template<class value_type>
void dtdtcl_1loop<value_type>::computeMi2pp
(std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
 std::vector<value_type>& Mi2p,std::vector<value_type>& Mi2pp,
 std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
 std::vector<value_type>& pMi2ptpp,
 const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
 const std::vector<value_type>& thetai,
 const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
 const std::vector<value_type>& dthetaidt,
 const std::vector<value_type>& d2kidt2,const std::vector<value_type>& d2kpidt2,
 const std::vector<value_type>& d2thetaidt2,
 const std::vector<value_type>& d3kidt3,const std::vector<value_type>& d3kpidt3,
 const std::vector<value_type>& d3thetaidt3,
 const value_type& phicl,const value_type& Z,const value_type& Zp,
 const value_type& Zpp,const value_type& Zppp,
 const value_type& dtdtcl,const value_type& d2tdtcl2,const value_type& logmu2)
{
    this->computeMi2p(Mi2,logfactor,Mi2p,pMi2pt,pMi2ptp,kiarray,kpi,thetai,
                      dkidt,dkpidt,dthetaidt,d2kidt2,d2kpidt2,
                      d2thetaidt2,phicl,Z,Zp,Zpp,dtdtcl,logmu2);
    //Pre-allocate to prevent any going out of bounds:
    Mi2pp.assign(this->nLoops,this->_0p0);
    pMi2ptpp.assign(this->nLoops,this->_0p0);
    value_type ki,kip,kipp,kippp,thetaip,thetaipp,thetaippp,kpip,kpipp,kpippp;
    for(int i = 0;i < Mi2.size();i++)
    {
        ki = kiarray[i];
        kip = dkidt[i];
        kipp = d2kidt2[i];
        kippp = d3kidt3[i];
        thetaip = dthetaidt[i];
        thetaipp = d2thetaidt2[i];
        thetaippp = d3thetaidt3[i];
        kpip = dkpidt[i];
        kpipp = d2kpidt2[i];
        kpippp = d3kpidt3[i];


        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter these:
        Mi2pp[i] = ((((d2tdtcl2*(this->R*thetaip)) + ((dtdtcl*dtdtcl
                *(this->R*thetaipp))
                + (((((d2tdtcl2 + (this->_2p0*dtdtcl))*kip)
                + (dtdtcl*dtdtcl*kipp))*(phicl*phicl*Z*Z))
                + ((this->_4p0*(dtdtcl*dtdtcl*(kip*(phicl*phicl*(Z*Zp)))))
                + (ki*(phicl*phicl*(Z*Z + ((this->_2p0*(dtdtcl*dtdtcl*Zp*Zp))
                + (this->_2p0*(Z*((d2tdtcl2*Zp) + ((this->_2p0*(dtdtcl*Zp))
                + (dtdtcl*dtdtcl*Zpp)))))))))))))
                - (dtdtcl*dtdtcl*kpipp)) - (d2tdtcl2*kpip));
        pMi2ptpp[i] = ((((d2tdtcl2*(this->R*thetaipp)) + ((dtdtcl*dtdtcl*
                    (this->R*thetaippp)) + (((kip + ((d2tdtcl2*kipp)
                    + ((this->_2p0*(dtdtcl*kipp)) + (dtdtcl*dtdtcl*kippp))))
                    *(phicl*phicl*Z*Z)) + ((this->_6p0*(dtdtcl*dtdtcl*(kip
                    *(phicl*phicl*Zp*Zp)))) + ((this->_2p0*(ki*(phicl*phicl
                    *(Zp*((d2tdtcl2*Zp) + ((this->_2p0*(dtdtcl*Zp))
                    + (this->_3p0*(dtdtcl*dtdtcl*Zpp))))))))
                    + (this->_2p0*(phicl*phicl*(Z*((this->_3p0*(dtdtcl*dtdtcl
                    *(kipp*Zp))) + ((kip*((this->_2p0*(d2tdtcl2*Zp))
                    + ((this->_4p0*(dtdtcl*Zp)) + (this->_3p0*(dtdtcl*dtdtcl*Zpp)))))
                    + (ki*(Zp + ((d2tdtcl2*Zpp) + ((this->_2p0*(dtdtcl*Zpp))
                    + (dtdtcl*dtdtcl*Zppp)))))))))))))))
                    - (dtdtcl*dtdtcl*kpippp)) - (d2tdtcl2*kpipp));

    }
}
template<class value_type>
void dtdtcl_1loop<value_type>::computeMi2ppp
(std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
 std::vector<value_type>& Mi2p,std::vector<value_type>& Mi2pp,
 std::vector<value_type>& Mi2ppp,
 std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
 std::vector<value_type>& pMi2ptpp,std::vector<value_type>& pMi2ptppp,
 const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
 const std::vector<value_type>& thetai,
 const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
 const std::vector<value_type>& dthetaidt,
 const std::vector<value_type>& d2kidt2,const std::vector<value_type>& d2kpidt2,
 const std::vector<value_type>& d2thetaidt2,
 const std::vector<value_type>& d3kidt3,const std::vector<value_type>& d3kpidt3,
 const std::vector<value_type>& d3thetaidt3,
 const std::vector<value_type>& d4kidt4,const std::vector<value_type>& d4kpidt4,
 const std::vector<value_type>& d4thetaidt4,
 const value_type& phicl,const value_type& Z,const value_type& Zp,
 const value_type& Zpp,const value_type& Zppp,const value_type& Zpppp,
 const value_type& dtdtcl,const value_type& d2tdtcl2,const value_type& d3tdtcl3,
 const value_type& logmu2)
{
    this->computeMi2pp(Mi2,logfactor,Mi2p,Mi2pp,pMi2pt,pMi2ptp,pMi2ptpp,
                       kiarray,kpi,
                       thetai,dkidt,dkpidt,dthetaidt,d2kidt2,d2kpidt2,
                       d2thetaidt2,d3kidt3,d3kpidt3,d3thetaidt3,phicl,Z,Zp,Zpp,
                       Zppp,dtdtcl,d2tdtcl2,logmu2);

    //Pre-allocate to prevent any going out of bounds:
    Mi2ppp.assign(this->nLoops,this->_0p0);
    pMi2ptppp.assign(this->nLoops,this->_0p0);
    value_type ki,kip,kipp,kippp,thetaip,thetaipp,thetaippp,kpip,kpipp,kpippp;
    value_type kipppp,thetaipppp,kpipppp;
    for(int i = 0;i < Mi2.size();i++)
    {
        ki = kiarray[i];
        kip = dkidt[i];
        kipp = d2kidt2[i];
        kippp = d3kidt3[i];
        kipppp = d4kidt4[i];
        thetaip = dthetaidt[i];
        thetaipp = d2thetaidt2[i];
        thetaippp = d3thetaidt3[i];
        thetaipppp = d4thetaidt4[i];
        kpip = dkpidt[i];
        kpipp = d2kpidt2[i];
        kpippp = d3kpidt3[i];
        kpipppp = d4kpidt4[i];

        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter these:
        Mi2ppp[i] = ((((this->_m3p0*(d2tdtcl2*(dtdtcl*kpipp)))
                + ((d3tdtcl3*(this->R*thetaip)) + ((this->_3p0*(d2tdtcl2*(dtdtcl*
                (this->R*thetaipp)))) + ((dtdtcl*dtdtcl*dtdtcl
                *(this->R*thetaippp))
                + ((((((this->_3p0*d2tdtcl2) + (d3tdtcl3 + (this->_3p0*dtdtcl)))*kip)
                + (dtdtcl*((this->_3p0*(d2tdtcl2*kipp)) + ((this->_3p0*(dtdtcl*kipp))
                + (dtdtcl*dtdtcl*kippp)))))*(phicl*phicl*Z*Z))
                + ((this->_6p0*(dtdtcl*dtdtcl*dtdtcl*(kip*(phicl*phicl*Zp*Zp))))
                + ((this->_6p0*(dtdtcl*(phicl*phicl*(Z*((dtdtcl*dtdtcl*(kipp*Zp))
                + (kip*((this->_2p0*(d2tdtcl2*Zp)) + ((this->_2p0*(dtdtcl*Zp))
                + (dtdtcl*dtdtcl*Zpp))))))))) + (ki*(phicl*phicl*(Z*Z
                + ((this->_6p0*(dtdtcl*(Zp*((d2tdtcl2*Zp) + ((dtdtcl*Zp)
                + (dtdtcl*dtdtcl*Zpp)))))) + (this->_2p0*(Z*((((this->_3p0*d2tdtcl2)
                + d3tdtcl3)*Zp) + ((this->_3p0*(dtdtcl*dtdtcl*Zpp))
                + ((this->_3p0*(dtdtcl*(Zp + (d2tdtcl2*Zpp))))
                + (dtdtcl*dtdtcl*dtdtcl*Zppp)))))))))))))))))
                - (dtdtcl*dtdtcl*dtdtcl*kpippp)) - (d3tdtcl3*kpip));
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter these:
        pMi2ptppp[i] = ((((this->_m3p0*(d2tdtcl2*(dtdtcl*kpippp)))
                    + ((d3tdtcl3*(this->R*thetaipp))
                    + ((this->_3p0*(d2tdtcl2*(dtdtcl
                    *(this->R*thetaippp)))) + ((dtdtcl*dtdtcl*dtdtcl
                    *(this->R*thetaipppp))
                    + ((phicl*phicl*(Z*((kip*Z) + (this->_2p0*(ki*Zp)))))
                    + ((this->_3p0*(dtdtcl*(phicl*phicl*((kipp*Z*Z)
                    + ((this->_2p0*(ki*Zp*Zp)) + (this->_2p0*(Z*((this->_2p0*(kip*Zp))
                    + (ki*Zpp))))))))) + ((this->_3p0*(phicl*phicl*((((d2tdtcl2
                    *kipp) + (dtdtcl*dtdtcl*kippp))*Z*Z)
                    + ((this->_2p0*(ki*(Zp*((d2tdtcl2*Zp)
                    + (this->_3p0*(dtdtcl*dtdtcl*Zpp))))))
                    + ((kip*((this->_4p0*(d2tdtcl2*(Z*Zp)))
                    + (this->_6p0*(dtdtcl*dtdtcl*(Zp*Zp + (Z*Zpp))))))
                    + (this->_2p0*(Z*((d2tdtcl2*(ki*Zpp))
                    + (dtdtcl*dtdtcl*((this->_3p0*(kipp*Zp))
                    + (ki*Zppp)))))))))))
                    + (phicl*phicl*((d3tdtcl3*((kipp*Z*Z)
                    + ((this->_2p0*(ki*Zp*Zp)) + (this->_2p0*(Z*((this->_2p0*(kip*Zp))
                    + (ki*Zpp))))))) + ((this->_3p0*(d2tdtcl2*(dtdtcl*((kippp*Z*Z)
                    + ((this->_6p0*(ki*(Zp*Zpp))) + ((this->_6p0*(kip*(Zp*Zp
                    + (Z*Zpp)))) + (this->_2p0*(Z*((this->_3p0*(kipp*Zp))
                    + (ki*Zppp)))))))))) + (dtdtcl*dtdtcl*dtdtcl*((kipppp
                    *Z*Z) + ((this->_12p0*(kipp*Zp*Zp)) + ((this->_6p0*(ki*Zpp*Zpp))
                    + ((this->_8p0*(Zp*((kippp*Z) + ((this->_3p0*(kip*Zpp))
                    + (ki*Zppp))))) + (this->_2p0*(Z*((this->_6p0*(kipp*Zpp))
                    + ((this->_4p0*(kip*Zppp)) + (ki*Zpppp))))))))))))))))))))
                    - (dtdtcl*dtdtcl*dtdtcl*kpipppp)) - (d3tdtcl3*kpipp));

    }
}
template<class value_type>
void dtdtcl_1loop<value_type>::computeMi2pppp
(std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
 std::vector<value_type>& Mi2p,std::vector<value_type>& Mi2pp,
 std::vector<value_type>& Mi2ppp,std::vector<value_type>& Mi2pppp,
 std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
 std::vector<value_type>& pMi2ptpp,std::vector<value_type>& pMi2ptppp,
 std::vector<value_type>& pMi2ptpppp,
 const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
 const std::vector<value_type>& thetai,
 const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
 const std::vector<value_type>& dthetaidt,
 const std::vector<value_type>& d2kidt2,const std::vector<value_type>& d2kpidt2,
 const std::vector<value_type>& d2thetaidt2,
 const std::vector<value_type>& d3kidt3,const std::vector<value_type>& d3kpidt3,
 const std::vector<value_type>& d3thetaidt3,
 const std::vector<value_type>& d4kidt4,const std::vector<value_type>& d4kpidt4,
 const std::vector<value_type>& d4thetaidt4,
 const std::vector<value_type>& d5kidt5,const std::vector<value_type>& d5kpidt5,
 const std::vector<value_type>& d5thetaidt5,
 const value_type& phicl,const value_type& Z,const value_type& Zp,
 const value_type& Zpp,const value_type& Zppp,const value_type& Zpppp,
 const value_type& Zppppp,
 const value_type& dtdtcl,const value_type& d2tdtcl2,const value_type& d3tdtcl3,
 const value_type& d4tdtcl4,const value_type& logmu2)
{
    this->computeMi2ppp(Mi2,logfactor,Mi2p,Mi2pp,Mi2ppp,pMi2pt,pMi2ptp,pMi2ptpp,
                        pMi2ptppp,kiarray,kpi,
                        thetai,dkidt,dkpidt,dthetaidt,d2kidt2,d2kpidt2,
                        d2thetaidt2,d3kidt3,d3kpidt3,d3thetaidt3,d4kidt4,
                        d4kpidt4,d4thetaidt4,phicl,Z,Zp,Zpp,
                        Zppp,Zpppp,dtdtcl,d2tdtcl2,d3tdtcl3,logmu2);

    //Pre-allocate to prevent any going out of bounds:
    Mi2pppp.assign(this->nLoops,this->_0p0);
    pMi2ptpppp.assign(this->nLoops,this->_0p0);
    value_type ki,kip,kipp,kippp,thetaip,thetaipp,thetaippp,kpip,kpipp,kpippp;
    value_type kipppp,kippppp,thetaipppp,thetaippppp,kpipppp,kpippppp;
    for(int i = 0;i < this->nLoops;i++)
    {
        ki = kiarray[i];
        kip = dkidt[i];
        kipp = d2kidt2[i];
        kippp = d3kidt3[i];
        kipppp = d4kidt4[i];
        kippppp = d5kidt5[i];
        thetaip = dthetaidt[i];
        thetaipp = d2thetaidt2[i];
        thetaippp = d3thetaidt3[i];
        thetaipppp = d4thetaidt4[i];
        thetaippppp = d5thetaidt5[i];
        kpip = dkpidt[i];
        kpipp = d2kpidt2[i];
        kpippp = d3kpidt3[i];
        kpipppp = d4kpidt4[i];
        kpippppp = d5kpidt5[i];

        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter these:
        Mi2pppp[i] = ((((this->_m3p0*(d2tdtcl2*d2tdtcl2*kpipp))
                        + ((this->_m4p0*(d3tdtcl3*(dtdtcl*kpipp)))
                        + ((this->_m6p0*(d2tdtcl2*(dtdtcl*dtdtcl*kpippp)))
                        + ((d4tdtcl4*(R*thetaip))
                        + ((this->_3p0*(d2tdtcl2*d2tdtcl2*(R*thetaipp)))
                        + ((this->_4p0*(d3tdtcl3*(dtdtcl*(R*thetaipp))))
                        + ((this->_6p0*(d2tdtcl2*(dtdtcl*dtdtcl*(R*thetaippp))))
                        + ((dtdtcl*dtdtcl*dtdtcl*dtdtcl*(R*thetaipppp))
                        + ((((((this->_6p0*d2tdtcl2) + ((this->_4p0*d3tdtcl3)
                        + (d4tdtcl4 + (this->_4p0*dtdtcl))))*kip)
                        + ((this->_3p0*(d2tdtcl2*d2tdtcl2*kipp))
                        + ((this->_4p0*(((this->_3p0*d2tdtcl2) + d3tdtcl3)
                        *(dtdtcl*kipp)))
                        + ((this->_4p0*(dtdtcl*dtdtcl*dtdtcl*kippp))
                        + ((this->_6p0*(dtdtcl*dtdtcl*(kipp
                        + (d2tdtcl2*kippp))))
                        + (dtdtcl*dtdtcl*dtdtcl*dtdtcl*kipppp))))))
                        *(phicl*phicl*Z*Z))
                        + ((this->_36p0*(d2tdtcl2*(dtdtcl*dtdtcl*(kip
                        *(phicl*phicl*Zp*Zp)))))
                        + ((this->_24p0*(dtdtcl*dtdtcl*dtdtcl
                        *(kip*(phicl*phicl*Zp*Zp))))
                        + ((this->_12p0*(dtdtcl*dtdtcl*dtdtcl*dtdtcl
                        *(kipp*(phicl*phicl*Zp*Zp))))
                        + ((this->_24p0*(dtdtcl*dtdtcl*dtdtcl*dtdtcl
                        *(kip*(phicl*phicl*(Zp*Zpp)))))
                        + ((this->_4p0*(phicl*phicl*(Z*((dtdtcl*dtdtcl
                        *((this->_9p0*(d2tdtcl2*(kipp*Zp)))
                        + ((this->_6p0*(dtdtcl*(kipp*Zp)))
                        + (dtdtcl*dtdtcl*((this->_2p0*(kippp*Zp))
                        + (this->_3p0*(kipp*Zpp)))))))
                        + (kip*((this->_3p0*(d2tdtcl2*d2tdtcl2*Zp))
                        + ((this->_4p0*(((this->_3p0*d2tdtcl2)
                        + d3tdtcl3)*(dtdtcl*Zp)))
                        + ((this->_6p0*(dtdtcl*dtdtcl*dtdtcl*Zpp))
                        + ((dtdtcl*dtdtcl*((this->_6p0*Zp)
                        + (this->_9p0*(d2tdtcl2*Zpp))))
                        + (this->_2p0*(dtdtcl*dtdtcl*dtdtcl*dtdtcl
                        *Zppp)))))))))))
                        + (ki*(phicl*phicl*(Z*Z
                        + ((this->_6p0*(d2tdtcl2*d2tdtcl2*Zp*Zp))
                        + ((this->_8p0*(((this->_3p0*d2tdtcl2)
                        + d3tdtcl3)*(dtdtcl*Zp*Zp)))
                        + ((this->_24p0*(dtdtcl*dtdtcl*dtdtcl*(Zp*Zpp)))
                        + ((this->_12p0*(dtdtcl*dtdtcl*(Zp*(Zp
                        + (this->_3p0*(d2tdtcl2*Zpp))))))
                        + ((dtdtcl*dtdtcl*dtdtcl*dtdtcl*((this->_6p0*Zpp*Zpp)
                        + (this->_8p0*(Zp*Zppp))))
                        + (this->_2p0*(Z*((((this->_6p0*d2tdtcl2)
                        + ((this->_4p0*d3tdtcl3) + d4tdtcl4))*Zp)
                        + ((this->_3p0*(d2tdtcl2*d2tdtcl2*Zpp))
                        + ((this->_4p0*(dtdtcl*(Zp + (((this->_3p0*d2tdtcl2)
                        + d3tdtcl3)*Zpp))))
                        + ((this->_4p0*(dtdtcl*dtdtcl*dtdtcl*Zppp))
                        + ((this->_6p0*(dtdtcl*dtdtcl*(Zpp + (d2tdtcl2*Zppp))))
                        + (dtdtcl*dtdtcl*dtdtcl*dtdtcl
                        *Zpppp))))))))))))))))))))))))))))))
                        - (dtdtcl*dtdtcl*dtdtcl*dtdtcl*kpipppp))
                        - (d4tdtcl4*kpip));
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter these:
        pMi2ptpppp[i] = ((((this->_m3p0*(d2tdtcl2*d2tdtcl2*kpippp))
                        + ((this->_m4p0*(d3tdtcl3*(dtdtcl*kpippp)))
                        + ((this->_m6p0*(d2tdtcl2*(dtdtcl*dtdtcl*kpipppp)))
                        + ((d4tdtcl4*(R*thetaipp))
                        + ((this->_3p0*(d2tdtcl2*d2tdtcl2*(R*thetaippp)))
                        + ((this->_4p0*(d3tdtcl3*(dtdtcl*(R*thetaippp))))
                        + ((this->_6p0*(d2tdtcl2*(dtdtcl*dtdtcl
                        *(R*thetaipppp))))
                        + ((dtdtcl*dtdtcl*dtdtcl*dtdtcl*(R*thetaippppp))
                        + ((phicl*phicl*(Z*((kip*Z) + (this->_2p0*(ki*Zp)))))
                        + ((this->_4p0*(dtdtcl*(phicl*phicl*((kipp*Z*Z)
                        + ((this->_2p0*(ki*Zp*Zp))
                        + (this->_2p0*(Z*((this->_2p0*(kip*Zp))
                        + (ki*Zpp)))))))))
                        + ((this->_6p0*(phicl*phicl*((((d2tdtcl2*kipp)
                        + (dtdtcl*dtdtcl*kippp))*Z*Z)
                        + ((this->_2p0*(ki*(Zp*((d2tdtcl2*Zp)
                        + (this->_3p0*(dtdtcl*dtdtcl*Zpp))))))
                        + ((kip*((this->_4p0*(d2tdtcl2*(Z*Zp)))
                        + (this->_6p0*(dtdtcl*dtdtcl*(Zp*Zp + (Z*Zpp))))))
                        + (this->_2p0*(Z*((d2tdtcl2*(ki*Zpp))
                        + (dtdtcl*dtdtcl*((this->_3p0*(kipp*Zp))
                        + (ki*Zppp)))))))))))
                        + ((this->_4p0*(phicl*phicl*((d3tdtcl3*((kipp*Z*Z)
                        + ((this->_2p0*(ki*Zp*Zp))
                        + (this->_2p0*(Z*((this->_2p0*(kip*Zp))
                        + (ki*Zpp)))))))
                        + ((this->_3p0*(d2tdtcl2*(dtdtcl*((kippp*Z*Z)
                        + ((this->_6p0*(ki*(Zp*Zpp)))
                        + ((this->_6p0*(kip*(Zp*Zp + (Z*Zpp))))
                        + (this->_2p0*(Z*((this->_3p0*(kipp*Zp))
                        + (ki*Zppp))))))))))
                        + (dtdtcl*dtdtcl*dtdtcl*((kipppp*Z*Z)
                        + ((this->_12p0*(kipp*Zp*Zp))
                        + ((this->_6p0*(ki*Zpp*Zpp))
                        + ((this->_8p0*(Zp*((kippp*Z)
                        + ((this->_3p0*(kip*Zpp)) + (ki*Zppp)))))
                        + (this->_2p0*(Z*((this->_6p0*(kipp*Zpp))
                        + ((this->_4p0*(kip*Zppp)) + (ki*Zpppp))))))))))))))
                        + (phicl*phicl*((d4tdtcl4*(kipp*Z*Z))
                        + ((this->_3p0*(d2tdtcl2*d2tdtcl2*(kippp*Z*Z)))
                        + ((this->_18p0*(d2tdtcl2*d2tdtcl2*(kipp*(Z*Zp))))
                        + ((this->_2p0*(d4tdtcl4*(ki*Zp*Zp)))
                        + ((this->_2p0*(d4tdtcl4*(ki*(Z*Zpp))))
                        + ((this->_18p0*(d2tdtcl2*d2tdtcl2*(ki*(Zp*Zpp))))
                        + ((this->_2p0*(kip*((this->_2p0*(d4tdtcl4*(Z*Zp)))
                        + ((this->_9p0*(d2tdtcl2*d2tdtcl2*Zp*Zp))
                        + (this->_9p0*(d2tdtcl2*d2tdtcl2*(Z*Zpp)))))))
                        + ((this->_6p0*(d2tdtcl2*d2tdtcl2*(ki*(Z*Zppp))))
                        + ((this->_4p0*(d3tdtcl3*(dtdtcl*((kippp*Z*Z)
                        + ((this->_6p0*(ki*(Zp*Zpp)))
                        + ((this->_6p0*(kip*(Zp*Zp + (Z*Zpp))))
                        + (this->_2p0*(Z*((this->_3p0*(kipp*Zp))
                        + (ki*Zppp))))))))))
                        + ((this->_6p0*(d2tdtcl2*(dtdtcl*dtdtcl*((kipppp*Z*Z)
                        + ((this->_12p0*(kipp*Zp*Zp))
                        + ((this->_6p0*(ki*Zpp*Zpp))
                        + ((this->_8p0*(Zp*((kippp*Z)
                        + ((this->_3p0*(kip*Zpp)) + (ki*Zppp)))))
                        + (this->_2p0*(Z*((this->_6p0*(kipp*Zpp))
                        + ((this->_4p0*(kip*Zppp)) + (ki*Zpppp))))))))))))
                        + (dtdtcl*dtdtcl*dtdtcl*dtdtcl*((kippppp*Z*Z)
                        + ((this->_20p0*(kippp*Zp*Zp))
                        + ((this->_20p0*(kippp*(Z*Zpp)))
                        + ((this->_20p0*(kipp*(Z*Zppp)))
                        + ((this->_20p0*(ki*(Zpp*Zppp)))
                        + ((this->_10p0*(Zp*((kipppp*Z)
                        + ((this->_6p0*(kipp*Zpp))
                        + ((this->_4p0*(kip*Zppp)) + (ki*Zpppp))))))
                        + ((this->_10p0*(kip*((this->_3p0*Zpp*Zpp)
                        + (Z*Zpppp))))
                        + (this->_2p0*(ki
                        *(Z*Zppppp))))))))))))))))))))))))))))))))))
                        - (dtdtcl*dtdtcl*dtdtcl*dtdtcl*kpippppp))
                        - (d4tdtcl4*kpipp));

    }
}
//The functions below are used to compute the derivatives
//d^nt/dtcl^n. However, since dt/dtcl = A/B is a rational function,
//this gets complicated very fast. We need to compute things like
//A,A',A''...,B,B',B'',... and combine them in complicated ways.
//This is what these functions are - A returns A, obviously, Ap = A',App = A''
//etc...
//
//Most of the hard work is done outside the function, however, and
//they take various vectors with inputs. This reduces the workload
//if we have to compute, say, A, Ap, App and Appp in order to
//get d^4t/dtcl^4, because the functions can share input which
//then doesn't need to be recomputed by all of them.
//A(x) - dtdtcl_1loop_no_logs
//
//Note - derivatives are with respect to tcl. Zp, Zpp etc... should be
//derivatives with respect to tcl, not t, but the functions calling
//these will take care of that.
template<class value_type>
value_type dtdtcl_1loop<value_type>::A
(const std::vector<value_type>& log_shared,const value_type& phicl,
 const std::vector<value_type>& Mi2,const std::vector<value_type>& ki,
 const value_type& Z)
 {
    //Access shared data:
    value_type res = this->_0p0;
    value_type H2 = this->R/(this->_12p0);

    //Compute A:
    for(int i = 0;i < this->nLoops;i++)
    {
        value_type logfactor = log_shared[i];

        //std::cout << "\nlogfactor = " << logfactor;
        //std::cout << "nilist[" << i << "] = " << this->nilist[i];
        //std::cout << "cilist[" << i << "] = " << this->cilist[i];
        //std::cout << "npilist[" << i << "] = " << this->npilist[i];
        //std::cout << "\nMi2[" << i << "] = " << Mi2[i];
        //std::cout << "ki[" << i << "] = " << ki[i];
        //std::cout << "Z = " << Z;
        //std::cout << "\nphicl = " << phicl;
        //std::cout << "\nH2 = " << H2;
        res += -this->_2p0*(this->nilist[i]*Mi2[i]*ki[i]*Z*Z*phicl*phicl*
                (logfactor - this->cilist[i]
                 + this->_1p0/this->_2p0  )
                 + this->npilist[i]*H2*H2*ki[i]*Z*Z*phicl*phicl/Mi2[i]);
    }
    //std::cout << "\nA = " << res;
    return res;
}
//------------------------------------------------------------------------------
//B(x) - dtdtcl_1loop_no_logs
template<class value_type>
value_type dtdtcl_1loop<value_type>::B
(const std::vector<value_type>& log_shared,
 const std::vector<value_type>& Mi2,const std::vector<value_type>& pMi2pt)
 {
    //Access shared data:
    value_type res = this->_0p0;
    value_type H2 = this->R/(this->_12p0);
    //Compute B:
    for(int i = 0;i < this->nLoops;i++)
    {
        value_type logfactor = log_shared[i];
        res += this->nilist[i]*(this->_2p0*Mi2[i]*pMi2pt[i]
                      *(logfactor - this->cilist[i]
                 + this->_1p0/this->_2p0  ) - Mi2[i]*Mi2[i])
                 -this->npilist[i]*H2*H2
                 +this->npilist[i]*H2*H2*pMi2pt[i]/Mi2[i];
    }
    //std::cout << "\nB = " << res;
    return res;
}
//------------------------------------------------------------------------------
//Ap(x) - dtdtcl_1loop_no_logs
//Note - this function expects Zp to be a t derivative, not a tcl derivative.
//It will do the conversion itself.
//Mi2p array is expected as a tcl derivative too - computeMi2p will do this
//automatically (if supplied with t derivatives).
//kiparray should be t derivatives (not tcl). This function will convert them.
//Note for future - perhaps we should simplify this
//so that tcl derivatives are always expected...
template<class value_type>
value_type dtdtcl_1loop<value_type>::Ap
(const value_type& Z,const value_type& Zd,
 const std::vector<value_type>& kiarray,
 const std::vector<value_type>& kiparray,
 const value_type& mu,const std::vector<value_type>& log_shared,
 const value_type& phicl,
 const std::vector<value_type>& Mi2array,
 const std::vector<value_type>& Mi2parray,
 const value_type& dtdtcl)
 {
    //Access shared data:
    value_type res = this->_0p0;
    value_type mu2 = mu*mu;
    value_type mu2p = mu2*dtdtcl;
    value_type Mi2,Mi2p,ci,ni,npi;
    value_type H2 = this->R/(this->_12p0);
    //ki,kip,thetaip,kpip;

    //t derivatives of Z are supplied. We need tcl derivatives!
    value_type Zp = Zd*dtdtcl;
    //Compute Ap:
    for(int i = 0;i < this->nLoops;i++)
    {
        //Needed coefficients:
        Mi2 = Mi2array[i];
        Mi2p = Mi2parray[i];
        ci = this->cilist[i];
        ni = this->nilist[i];
        npi = this->npilist[i];
        value_type ki = kiarray[i];
        value_type kip = dtdtcl*kiparray[i];

        /*
        std::cout << "\nMi2 = " << Mi2 << " Mi2p = " << Mi2p << " ci = " << ci
                  << " ni = " << ni << " npi = " << npi << " ki = " << ki
                  << " kip = " << kip;*/

        value_type logfactor = log_shared[i];
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter this:
        res += (this->_m1p0*(((this->_1p0)/(Mi2*Mi2))*(((this->_1p0)/(mu2))
                *(phicl*phicl*(Z*(((ki*((this->_3p0
                + ((this->_m2p0*ci) + (this->_2p0*logfactor)))
                *(Mi2*Mi2*(Mi2p*(mu2*(ni*Z))))))
                + ((H2*H2*(ki*(Mi2*(mu2*(npi*(Z + (this->_2p0*Zp)))))))
                + (Mi2*Mi2*Mi2*(ni*((kip*((this->_1p0 + ((this->_m2p0*ci)
                + (this->_2p0*logfactor)))*(mu2*Z)))
                + (ki*((this->_m2p0*(mu2p*Z))
                - ((this->_m1p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*(mu2*(Z + (this->_2p0*Zp)))))))))))
                - (H2*H2*(ki*(Mi2p*(mu2*(npi*Z)))))))))));
    }
    return res;
}
//------------------------------------------------------------------------------
//App(x) - dtdtcl_1loop_no_logs
template<class value_type>
value_type dtdtcl_1loop<value_type>::App
(const value_type& Z,const value_type& Zd,const value_type& Zdd,
 const std::vector<value_type>& kiarray,
 const std::vector<value_type>& kiparray,
 const std::vector<value_type>& kipparray,
 const value_type& mu,const std::vector<value_type>& log_shared,
 const value_type& phicl,
 const std::vector<value_type>& Mi2array,
 const std::vector<value_type>& Mi2parray,
 const std::vector<value_type>& Mi2pparray,
 const value_type& dtdtcl,const value_type& d2tdtcl2)
{
    //Access shared data:
    value_type res = this->_0p0;
    //mu^2 and its (tcl) derivatives:
    value_type mu2 = mu*mu;
    value_type mu2p = mu2*dtdtcl;
    value_type mu2pp = (d2tdtcl2 + dtdtcl*dtdtcl)*mu2;
    value_type Mi2,Mi2p,Mi2pp,ci,ni,npi;
    //t derivatives of Z are supplied. We need tcl derivatives!
    value_type Zp = Zd*dtdtcl;
    value_type Zpp = Zd*d2tdtcl2 + Zdd*dtdtcl*dtdtcl;
    value_type H2 = this->R/(this->_12p0);
    //ki,kip,kipp,thetaip,kpip,thetaipp,kpipp;
    //Compute App:
    for(int i = 0;i < this->nLoops;i++)
    {
        //Needed coefficients:
        Mi2 = Mi2array[i];

        Mi2p = Mi2parray[i];
        Mi2pp = Mi2pparray[i];
        ci = this->cilist[i];
        ni = this->nilist[i];
        npi = this->npilist[i];
        //Derivatives of ki with respect to tcl (not t):
        //NB - kiparray, kipparray are derivatives wrt t, not tcl, so we
        //have to apply jacobian formulae:
        value_type ki = kiarray[i];
        value_type kip = dtdtcl*kiparray[i];
        value_type kipp = d2tdtcl2*kiparray[i] + dtdtcl*dtdtcl*kipparray[i];

        value_type logfactor = log_shared[i];
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter this:
        res += (this->_m1p0*(((this->_1p0)/(Mi2*Mi2*Mi2))
                *(((this->_1p0)/(mu2*mu2))*(phicl*phicl*(((this->_2p0*(H2*H2
                *(ki*(Mi2p*Mi2p*(mu2*mu2*(npi*Z*Z))))))
                + ((Mi2*Mi2*Mi2*(mu2*(ni*(Z*((this->_2p0*(kip*((this->_3p0
                + ((this->_m2p0*ci) + (this->_2p0*logfactor)))*(Mi2p*(mu2*Z)))))
                + (ki*((this->_m4p0*(Mi2p*(mu2p*Z))) - ((this->_m3p0
                + ((this->_2p0*ci) + (this->_m2p0*logfactor)))*(mu2
                *((((this->_2p0*Mi2p) + Mi2pp)*Z)
                + (this->_4p0*(Mi2p*Zp))))))))))))
                + ((Mi2*Mi2*(mu2*mu2*((((this->_2p0*(ki*(Mi2p*Mi2p*ni)))
                + (H2*H2*(ki*npi)))*Z*Z)
                + ((this->_2p0*(H2*H2*(ki*(npi*Zp*Zp))))
                + (this->_2p0*(H2*H2*(ki*(npi*(Z*((this->_2p0*Zp)
                + Zpp))))))))))
                + (Mi2*Mi2*Mi2*Mi2*(ni*((mu2*(Z*((this->_m4p0*(kip*(mu2p*Z)))
                - ((this->_m1p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*(mu2*((((this->_2p0*kip) + kipp)*Z)
                + (this->_4p0*(kip*Zp))))))))
                + (ki*(((this->_2p0*(mu2p*mu2p*Z*Z))
                + (this->_m2p0*(mu2*(Z*((((this->_2p0*mu2p) + mu2pp)*Z)
                + (this->_4p0*(mu2p*Zp)))))))
                - ((this->_m1p0 + ((this->_2p0*ci) + (this->_m2p0*logfactor)))
                *(mu2*mu2*(Z*Z + ((this->_2p0*Zp*Zp)
                + (this->_2p0*(Z*((this->_2p0*Zp) + Zpp)))))))))))))))
                - (H2*H2*(ki*(Mi2*(mu2*mu2*(npi*(Z*((((this->_2p0*Mi2p)
                + Mi2pp)*Z) + (this->_4p0*(Mi2p*Zp))))))))))))));
    }
    return res;
}
//------------------------------------------------------------------------------
//Appp(x) - dtdtcl_1loop_no_logs
template<class value_type>
value_type dtdtcl_1loop<value_type>::Appp
(const value_type& Z,const value_type& Zd,const value_type& Zdd,
 const value_type& Zddd,
 const std::vector<value_type>& kiarray,const std::vector<value_type>& kiparray,
 const std::vector<value_type>& kipparray,
 const std::vector<value_type>& kippparray,
 const value_type& mu,const std::vector<value_type>& log_shared,
 const value_type& phicl,
 const std::vector<value_type>& Mi2array,
 const std::vector<value_type>& Mi2parray,
 const std::vector<value_type>& Mi2pparray,
 const std::vector<value_type>& Mi2ppparray,
 const value_type& dtdtcl,const value_type& d2tdtcl2,
 const value_type& d3tdtcl3)
{
    //Access shared data:
    value_type res = this->_0p0;
    //mu^2 and its (tcl) derivatives:
    value_type mu2 = mu*mu;
    value_type mu2p = mu2*dtdtcl;
    value_type mu2pp = (d2tdtcl2 + dtdtcl*dtdtcl)*mu2;
    value_type mu2ppp = (d3tdtcl3
                         + this->_3p0*d2tdtcl2*dtdtcl
                         + dtdtcl*dtdtcl*dtdtcl)*mu2;
    value_type H2 = this->R/(this->_12p0);
    value_type Mi2,Mi2p,Mi2pp,Mi2ppp,ci,ni,npi;

    //t derivatives of Z are supplied. We need tcl derivatives!
    value_type Zp = Zd*dtdtcl;
    value_type Zpp = Zd*d2tdtcl2 + Zdd*dtdtcl*dtdtcl;
    value_type Zppp = Zd*d3tdtcl3 + this->_3p0*Zdd*d2tdtcl2*dtdtcl +
                        Zddd*dtdtcl*dtdtcl*dtdtcl;

    //ki,kip,kipp,thetaip,kpip,thetaipp,kpipp;
    //value_type kippp,kpippp,thetaippp;
    //Compute Appp:
    for(int i = 0;i < this->nLoops;i++)
    {
        //Needed coefficients:
        Mi2 = Mi2array[i];
        Mi2p = Mi2parray[i];
        Mi2pp = Mi2pparray[i];
        Mi2ppp = Mi2ppparray[i];
        ci = this->cilist[i];
        ni = this->nilist[i];
        npi = this->npilist[i];
        //Derivatives of ki with respect to tcl (not t):
        //NB - kiparray, kipparray are derivatives wrt t, not tcl, so we
        //have to apply jacobian formulae:
        value_type ki = kiarray[i];
        value_type kip = dtdtcl*kiparray[i];
        value_type kipp = d2tdtcl2*kiparray[i] + dtdtcl*dtdtcl*kipparray[i];
        value_type kippp = d3tdtcl3*kiparray[i]
                           + this->_3p0*dtdtcl*d2tdtcl2*kipparray[i]
                            + dtdtcl*dtdtcl*dtdtcl*kippparray[i];

        value_type logfactor = log_shared[i];
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter this:
        res += (this->_m1p0*(((this->_1p0)/(Mi2*Mi2*Mi2*Mi2))
                *(((this->_1p0)/(mu2*mu2*mu2))
                *(phicl*phicl*(((this->_m6p0*(H2*H2*(ki*(Mi2p*Mi2p*Mi2p
                *(mu2*mu2*mu2*(npi*Z*Z))))))
                + ((this->_6p0*(H2*H2*(ki*(Mi2*(Mi2p*(mu2*mu2*mu2*(npi*(Z
                *(((Mi2p + Mi2pp)*Z) + (this->_2p0*(Mi2p*Zp)))))))))))
                + ((Mi2*Mi2*Mi2*Mi2*(mu2*(ni*((this->_3p0*(mu2*(Z
                *((this->_m4p0*(kip*(Mi2p*(mu2p*Z))))
                - ((this->_m3p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*(mu2*((((kipp*Mi2p)
                + (kip*((this->_2p0*Mi2p) + Mi2pp)))*Z)
                + (this->_4p0*(kip*(Mi2p*Zp))))))))))
                + (ki*(((this->_6p0*(Mi2p*(mu2p*mu2p*Z*Z)))
                + (this->_m6p0*(mu2*(Z*((((Mi2pp*mu2p)
                + (Mi2p*((this->_2p0*mu2p) + mu2pp)))*Z)
                + (this->_4p0*(Mi2p*(mu2p*Zp))))))))
                - ((this->_m3p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*(mu2*mu2*((((this->_3p0*Mi2p)
                + ((this->_3p0*Mi2pp) + Mi2ppp))*Z*Z)
                + ((this->_6p0*(Mi2p*Zp*Zp)) + (this->_6p0*(Z*((Mi2pp*Zp)
                + (Mi2p*((this->_2p0*Zp) + Zpp)))))))))))))))
                + ((Mi2*Mi2*Mi2*(mu2*mu2*mu2
                *((((this->_6p0*(kip*(Mi2p*Mi2p*ni)))
                + ((this->_6p0*(ki*(Mi2p*((Mi2p + Mi2pp)*ni))))
                + (H2*H2*(ki*npi))))*Z*Z)
                + ((this->_6p0*(H2*H2*(ki*(npi*(Zp*(Zp + Zpp))))))
                + (this->_2p0*(Z*((this->_3p0*(((this->_2p0
                *(ki*(Mi2p*Mi2p*ni))) + (H2*H2*(ki*npi)))*Zp))
                + (H2*H2*(ki*(npi*((this->_3p0*Zpp) + Zppp)))))))))))
                + (Mi2*Mi2*Mi2*Mi2*Mi2*(ni*((mu2*(((this->_6p0
                *(kip*(mu2p*mu2p*Z*Z))) + (this->_m6p0*(mu2*(Z*((((kipp*mu2p)
                + (kip*((this->_2p0*mu2p) + mu2pp)))*Z)
                + (this->_4p0*(kip*(mu2p*Zp))))))))
                - ((this->_m1p0 + ((this->_2p0*ci) + (this->_m2p0*logfactor)))
                *(mu2*mu2*((((this->_3p0*kip)
                + ((this->_3p0*kipp) + kippp))*Z*Z) + ((this->_6p0*(kip*Zp*Zp))
                + (this->_6p0*(Z*((kipp*Zp)
                + (kip*((this->_2p0*Zp) + Zpp)))))))))))
                + (ki*(((this->_m4p0*(mu2p*mu2p*mu2p*Z*Z))
                + ((this->_6p0*(mu2*(mu2p*(Z*(((mu2p + mu2pp)*Z)
                + (this->_2p0*(mu2p*Zp)))))))
                + (this->_m2p0*(mu2*mu2*((((this->_3p0*mu2p)
                + ((this->_3p0*mu2pp) + mu2ppp))*Z*Z)
                + ((this->_6p0*(mu2p*Zp*Zp)) + (this->_6p0*(Z*((mu2pp*Zp)
                + (mu2p*((this->_2p0*Zp) + Zpp)))))))))))
                - ((this->_m1p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*(mu2*mu2*mu2*(Z*Z
                + ((this->_6p0*(Zp*(Zp + Zpp)))
                + (this->_2p0*(Z*((this->_3p0*Zp)
                + ((this->_3p0*Zpp) + Zppp)))))))))))))))))
                - (Mi2*Mi2*(mu2*mu2*mu2*((((this->_2p0*(ki*(Mi2p*Mi2p*Mi2p*ni)))
                + ((this->_3p0*(H2*H2*(ki*(Mi2p*npi))))
                + (H2*H2*(ki*(((this->_3p0*Mi2pp) + Mi2ppp)*npi)))))*Z*Z)
                + ((this->_6p0*(H2*H2*(ki*(Mi2p*(npi*Zp*Zp)))))
                + (this->_6p0*(H2*H2*(ki*(npi*(Z*((Mi2pp*Zp)
                + (Mi2p*((this->_2p0*Zp) + Zpp)))))))))))))))));
    }
    return res;
}
//Apppp(x) - dtdtcl_1loop_no_logs
template<class value_type>
value_type dtdtcl_1loop<value_type>::Apppp
(const value_type& Z,const value_type& Zd,const value_type& Zdd,
 const value_type& Zddd,
 const value_type& Zdddd,
 const std::vector<value_type>& kiarray,
 const std::vector<value_type>& kiparray,
 const std::vector<value_type>& kipparray,
 const std::vector<value_type>& kippparray,
 const std::vector<value_type>& kipppparray,
 const value_type& mu,const std::vector<value_type>& log_shared,
 const value_type& phicl,
 const std::vector<value_type>& Mi2array,
 const std::vector<value_type>& Mi2parray,
 const std::vector<value_type>& Mi2pparray,
 const std::vector<value_type>& Mi2ppparray,
 const std::vector<value_type>& Mi2pppparray,
 const value_type& dtdtcl,const value_type& d2tdtcl2,
 const value_type& d3tdtcl3,
 const value_type& d4tdtcl4)
{
    /*
    std::cout << "\nkiarray.size() = " << kiarray.size();
    std::cout << "\nkiparray.size() = " << kiparray.size();
    std::cout << "\nkipparray.size() = " << kipparray.size();
    std::cout << "\nkippparray.size() = " << kippparray.size();
    std::cout << "\nkipppparray.size() = " << kipppparray.size();
    std::cout << "\nlog_shared.size() = " << log_shared.size();
    std::cout << "\nMi2array.size() = " << Mi2array.size();
    std::cout << "\nMi2parray.size() = " << Mi2parray.size();
    std::cout << "\nMi2pparray.size() = " << Mi2pparray.size();
    std::cout << "\nMi2ppparray.size() = " << Mi2ppparray.size();
    std::cout << "\nMi2pppparray.size() = " << Mi2pppparray.size();
    */

    //Access shared data:
    value_type res = this->_0p0;
    //mu^2 and its (tcl) derivatives:
    value_type mu2 = mu*mu;
    value_type mu2p = mu2*dtdtcl;
    value_type mu2pp = (d2tdtcl2 + dtdtcl*dtdtcl)*mu2;
    value_type mu2ppp = (d3tdtcl3
                         + this->_3p0*d2tdtcl2*dtdtcl
                         + dtdtcl*dtdtcl*dtdtcl)*mu2;
    value_type mu2pppp = (d4tdtcl4 + this->_4p0*d3tdtcl3*dtdtcl
                         + this->_3p0*d2tdtcl2*d2tdtcl2
                         + this->_6p0*d2tdtcl2*dtdtcl*dtdtcl
                         + dtdtcl*dtdtcl*dtdtcl*dtdtcl)*mu2;
    value_type H2 = this->R/(this->_12p0);
    value_type Mi2,Mi2p,Mi2pp,Mi2ppp,Mi2pppp,ci,ni,npi;
    //ki,kip,kipp,thetaip,kpip,thetaipp,kpipp;
    //value_type kippp,kpippp,thetaippp;

    //t derivatives of Z are supplied. We need tcl derivatives!
    value_type Zp = Zd*dtdtcl;
    value_type Zpp = Zd*d2tdtcl2 + Zdd*dtdtcl*dtdtcl;
    value_type Zppp = Zd*d3tdtcl3 + this->_3p0*Zdd*d2tdtcl2*dtdtcl
                        + Zddd*dtdtcl*dtdtcl*dtdtcl;
    value_type Zpppp = Zd*d4tdtcl4
                        + this->_6p0*Zddd*d2tdtcl2*dtdtcl*dtdtcl
                        + this->_4p0*Zdd*d3tdtcl3*dtdtcl
                        + this->_3p0*Zdd*d2tdtcl2*d2tdtcl2
                        + Zdddd*dtdtcl*dtdtcl*dtdtcl*dtdtcl;

    //Compute Appp:
    for(int i = 0;i < this->nLoops;i++)
    {
        //Needed coefficients:
        Mi2 = Mi2array[i];
        Mi2p = Mi2parray[i];
        Mi2pp = Mi2pparray[i];
        Mi2ppp = Mi2ppparray[i];
        Mi2pppp = Mi2pppparray[i];
        ci = this->cilist[i];
        ni = this->nilist[i];
        npi = this->npilist[i];
        //Derivatives of ki with respect to tcl (not t):
        //NB - kiparray, kipparray are derivatives wrt t, not tcl, so we
        //have to apply jacobian formulae:
        value_type ki = kiarray[i];
        value_type kip = dtdtcl*kiparray[i];
        value_type kipp = d2tdtcl2*kiparray[i] + dtdtcl*dtdtcl*kipparray[i];
        value_type kippp = d3tdtcl3*kiparray[i]
                           + this->_3p0*dtdtcl*d2tdtcl2*kipparray[i]
                            + dtdtcl*dtdtcl*dtdtcl*kippparray[i];
        value_type kipppp = this->_3p0*kipparray[i]*d2tdtcl2*d2tdtcl2
                            +this->_6p0*dtdtcl*dtdtcl*d2tdtcl2*kippparray[i]
                            + this->_4p0*dtdtcl*d3tdtcl3*kipparray[i]
                            + dtdtcl*dtdtcl*dtdtcl*dtdtcl*kipppparray[i]
                            + kiparray[i]*d4tdtcl4;

        value_type logfactor = log_shared[i];
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter this:
        res += (this->_m1p0*(((this->_1p0)/(Mi2*Mi2*Mi2*Mi2*Mi2))
                *(((this->_1p0)/(mu2*mu2*mu2*mu2))
                *(phicl*phicl*(((this->_24p0*(H2*H2*(ki*(Mi2p*Mi2p*Mi2p*Mi2p
                *(mu2*mu2*mu2*mu2*(npi*Z*Z))))))
                + ((this->_m12p0*(H2*H2*(ki*(Mi2*(Mi2p*Mi2p*(mu2*mu2*mu2*mu2
                *(npi*(Z*((((this->_2p0*Mi2p) + (this->_3p0*Mi2pp))*Z)
                + (this->_4p0*(Mi2p*Zp)))))))))))
                + ((this->_2p0*(Mi2*Mi2*(mu2*mu2*mu2*mu2*((((this->_2p0
                *(ki*(Mi2p*Mi2p*Mi2p*Mi2p*ni)))
                + ((this->_6p0*(H2*H2*(ki*(Mi2p*Mi2p*npi))))
                + ((this->_3p0*(H2*H2*(ki*(Mi2pp*Mi2pp*npi))))
                + (this->_4p0*(H2*H2*(ki*(Mi2p*(((this->_3p0*Mi2pp)
                + Mi2ppp)*npi))))))))*Z*Z)
                + ((this->_12p0*(H2*H2*(ki*(Mi2p*Mi2p*(npi*Zp*Zp)))))
                + (this->_12p0*(H2*H2*(ki*(Mi2p*(npi*(Z*((this->_2p0*(Mi2pp*Zp))
                + (Mi2p*((this->_2p0*Zp) + Zpp))))))))))))))
                + ((Mi2*Mi2*Mi2*Mi2*Mi2*(mu2*(ni*((this->_2p0*(mu2
                *(((this->_12p0*(kip*(Mi2p*(mu2p*mu2p*Z*Z))))
                + (this->_m12p0*(mu2*(Z*((((kipp*(Mi2p*mu2p))
                + (kip*((Mi2pp*mu2p) + (Mi2p*((this->_2p0*mu2p) + mu2pp)))))*Z)
                + (this->_4p0*(kip*(Mi2p*(mu2p*Zp)))))))))
                - ((this->_m3p0 + ((this->_2p0*ci) + (this->_m2p0*logfactor)))
                *(mu2*mu2*((((this->_2p0*(((this->_3p0*kipp) + kippp)*Mi2p))
                + ((this->_3p0*(kipp*Mi2pp))
                + (this->_2p0*(kip*((this->_3p0*Mi2p) + ((this->_3p0*Mi2pp)
                + Mi2ppp))))))*Z*Z) + ((this->_12p0*(kip*(Mi2p*Zp*Zp)))
                + (this->_12p0*(Z*((kipp*(Mi2p*Zp)) + (kip*((Mi2pp*Zp)
                + (Mi2p*((this->_2p0*Zp) + Zpp))))))))))))))
                + (ki*(((this->_m16p0*(Mi2p*(mu2p*mu2p*mu2p*Z*Z)))
                + ((this->_12p0*(mu2*(mu2p*(Z*((((Mi2pp*mu2p)
                + (this->_2p0*(Mi2p*(mu2p + mu2pp))))*Z)
                + (this->_4p0*(Mi2p*(mu2p*Zp))))))))
                + (this->_m4p0*(mu2*mu2*((((this->_2p0*(((this->_3p0*Mi2pp)
                + Mi2ppp)*mu2p)) + ((this->_3p0*(Mi2pp*mu2pp))
                + (this->_2p0*(Mi2p*((this->_3p0*mu2p) + ((this->_3p0*mu2pp)
                + mu2ppp))))))*Z*Z) + ((this->_12p0*(Mi2p*(mu2p*Zp*Zp)))
                + (this->_12p0*(Z*((Mi2pp*(mu2p*Zp)) + (Mi2p*((mu2pp*Zp)
                + (mu2p*((this->_2p0*Zp) + Zpp)))))))))))))
                - ((this->_m3p0 + ((this->_2p0*ci) + (this->_m2p0*logfactor)))
                *(mu2*mu2*mu2*((((this->_4p0*Mi2p) + ((this->_6p0*Mi2pp)
                + ((this->_4p0*Mi2ppp) + Mi2pppp)))*Z*Z)
                + ((this->_12p0*(Zp*((Mi2pp*Zp)
                + (this->_2p0*(Mi2p*(Zp + Zpp))))))
                + (this->_4p0*(Z*((this->_2p0*(((this->_3p0*Mi2pp)
                + Mi2ppp)*Zp)) + ((this->_3p0*(Mi2pp*Zpp))
                + (this->_2p0*(Mi2p*((this->_3p0*Zp)
                + ((this->_3p0*Zpp) + Zppp))))))))))))))))))
                + ((Mi2*Mi2*Mi2*Mi2*(mu2*mu2*mu2*mu2*((((this->_12p0
                *(kipp*(Mi2p*Mi2p*ni)))
                + ((this->_24p0*(kip*(Mi2p*((Mi2p + Mi2pp)*ni))))
                + ((this->_2p0*(ki*(((this->_6p0*Mi2p*Mi2p)
                + ((this->_3p0*Mi2pp*Mi2pp)
                + (this->_4p0*(Mi2p*((this->_3p0*Mi2pp) + Mi2ppp)))))*ni)))
                + (H2*H2*(ki*npi)))))*Z*Z)
                + ((this->_24p0*(ki*(Mi2p*Mi2p*(ni*Zp*Zp))))
                + ((this->_12p0*(H2*H2*(ki*(npi*Zp*Zp))))
                + ((this->_24p0*(H2*H2*(ki*(npi*(Zp*Zpp)))))
                + ((this->_6p0*(H2*H2*(ki*(npi*Zpp*Zpp))))
                + ((this->_8p0*(H2*H2*(ki*(npi*(Zp*Zppp)))))
                + (this->_2p0*(Z*((this->_4p0*(((this->_6p0
                *(kip*(Mi2p*Mi2p*ni)))
                + ((this->_6p0*(ki*(Mi2p*((Mi2p + Mi2pp)*ni))))
                + (H2*H2*(ki*npi))))*Zp))
                + ((this->_6p0*(((this->_2p0*(ki*(Mi2p*Mi2p*ni)))
                + (H2*H2*(ki*npi)))*Zpp)) + (H2*H2*(ki*(npi*((this->_4p0*Zppp)
                + Zpppp))))))))))))))))
                + (Mi2*Mi2*Mi2*Mi2*Mi2*Mi2*(ni*((mu2
                *(((this->_m16p0*(kip*(mu2p*mu2p*mu2p*Z*Z)))
                + ((this->_12p0*(mu2*(mu2p*(Z*((((kipp*mu2p)
                + (this->_2p0*(kip*(mu2p + mu2pp))))*Z)
                + (this->_4p0*(kip*(mu2p*Zp))))))))
                + (this->_m4p0*(mu2*mu2*((((this->_2p0*(((this->_3p0*kipp)
                + kippp)*mu2p)) + ((this->_3p0*(kipp*mu2pp))
                + (this->_2p0*(kip*((this->_3p0*mu2p) + ((this->_3p0*mu2pp)
                + mu2ppp))))))*Z*Z) + ((this->_12p0*(kip*(mu2p*Zp*Zp)))
                + (this->_12p0*(Z*((kipp*(mu2p*Zp)) + (kip*((mu2pp*Zp)
                + (mu2p*((this->_2p0*Zp) + Zpp)))))))))))))
                - ((this->_m1p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*(mu2*mu2*mu2*((((this->_4p0*kip)
                + ((this->_6p0*kipp) + ((this->_4p0*kippp) + kipppp)))*Z*Z)
                + ((this->_12p0*(Zp*((kipp*Zp)
                + (this->_2p0*(kip*(Zp + Zpp))))))
                + (this->_4p0*(Z*((this->_2p0*(((this->_3p0*kipp) + kippp)*Zp))
                + ((this->_3p0*(kipp*Zpp)) + (this->_2p0*(kip*((this->_3p0*Zp)
                + ((this->_3p0*Zpp) + Zppp))))))))))))))
                + (ki*(((this->_12p0*(mu2p*mu2p*mu2p*mu2p*Z*Z))
                + ((this->_m8p0*(mu2*(mu2p*mu2p*(Z*((((this->_2p0*mu2p)
                + (this->_3p0*mu2pp))*Z) + (this->_4p0*(mu2p*Zp)))))))
                + ((this->_2p0*(mu2*mu2*((((this->_6p0*mu2p*mu2p)
                + ((this->_3p0*mu2pp*mu2pp)
                + (this->_4p0*(mu2p*((this->_3p0*mu2pp) + mu2ppp)))))*Z*Z)
                + ((this->_12p0*(mu2p*mu2p*Zp*Zp))
                + (this->_12p0*(mu2p*(Z*((this->_2p0*(mu2pp*Zp))
                + (mu2p*((this->_2p0*Zp) + Zpp))))))))))
                + (this->_m2p0*(mu2*mu2*mu2*((((this->_4p0*mu2p)
                + ((this->_6p0*mu2pp) + ((this->_4p0*mu2ppp) + mu2pppp)))*Z*Z)
                + ((this->_12p0*(Zp*((mu2pp*Zp)
                + (this->_2p0*(mu2p*(Zp + Zpp))))))
                + (this->_4p0*(Z*((this->_2p0*(((this->_3p0*mu2pp)
                + mu2ppp)*Zp)) + ((this->_3p0*(mu2pp*Zpp))
                + (this->_2p0*(mu2p*((this->_3p0*Zp)
                + ((this->_3p0*Zpp) + Zppp)))))))))))))))
                - ((this->_m1p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*(mu2*mu2*mu2*mu2*(Z*Z
                + ((this->_12p0*Zp*Zp) + ((this->_6p0*Zpp*Zpp)
                + ((this->_8p0*(Zp*((this->_3p0*Zpp) + Zppp)))
                + (this->_2p0*(Z*((this->_4p0*Zp) + ((this->_6p0*Zpp)
                + ((this->_4p0*Zppp) + Zpppp)))))))))))))))))))))
                - (Mi2*Mi2*Mi2*(mu2*mu2*mu2*mu2*((((this->_8p0*((ki + kip)
                *(Mi2p*Mi2p*Mi2p*ni)))
                + ((this->_12p0*(ki*(Mi2p*Mi2p*(Mi2pp*ni))))
                + ((this->_4p0*(H2*H2*(ki*(Mi2p*npi))))
                + (H2*H2*(ki*(((this->_6p0*Mi2pp)
                + ((this->_4p0*Mi2ppp) + Mi2pppp))*npi))))))*Z*Z)
                + ((this->_12p0*(H2*H2*(ki*(npi*(Zp*((Mi2pp*Zp)
                + (this->_2p0*(Mi2p*(Zp + Zpp)))))))))
                + (this->_4p0*(Z*((this->_4p0*(ki*(Mi2p*Mi2p*Mi2p*(ni*Zp))))
                + ((H2*H2*(ki*(npi*((this->_2p0*(((this->_3p0*Mi2pp)
                + Mi2ppp)*Zp)) + (this->_3p0*(Mi2pp*Zpp))))))
                + (this->_2p0*(H2*H2*(ki*(Mi2p*(npi*((this->_3p0*Zp)
                + ((this->_3p0*Zpp) + Zppp))))))))))))))))))));
    }
    return res;
}
//------------------------------------------------------------------------------
//Bp(x) - dtdtcl_1loop_no_logs
template<class value_type>
value_type dtdtcl_1loop<value_type>::Bp
(const value_type& dtdtcl,const value_type& mu,const value_type& phicl,
 const std::vector<value_type>& log_shared,
 const std::vector<value_type>& Mi2array,
 const std::vector<value_type>& Mi2parray,
 const std::vector<value_type>& pMi2ptarray,
 const std::vector<value_type>& pMi2ptparray)
{
    //Access shared data:
    value_type res = this->_0p0;
    //mu^2 and its (tcl) derivatives:
    value_type mu2 = mu*mu;
    value_type mu2p = mu2*dtdtcl;
    value_type H2 = this->R/(this->_12p0);

    //std::cout << "\nmu2 = " << mu2 << " mu2p = " << mu2p << " H2 = " << H2;
    value_type Mi2,Mi2p,ci,ni,npi;
    //ki,kip,kipp,thetaip,kpip,kpipp,thetaipp;
    value_type pMi2pt,pMi2ptp;
    //Compute Bp:
    //std::cout << "\nMi2 (";
    for(int i = 0;i < this->nLoops;i++)
    {
        //Needed coefficients:
        Mi2 = Mi2array[i];
        Mi2p = Mi2parray[i];
        pMi2pt = pMi2ptarray[i];
        pMi2ptp = pMi2ptparray[i];
        ci = this->cilist[i];
        ni = this->nilist[i];
        npi = this->npilist[i];


        /*
        std::cout << "\nMi2 = " << Mi2 << " Mi2p = " << Mi2p << " ci = " << ci
                  << " ni = " << ni
                  << " pMi2pt = " << pMi2pt
                  << " pMi2ptp = " << pMi2ptp;*/




        value_type logfactor = log_shared[i];
        //if(i != 0)
        //{
            //std::cout << ",";
        //}
        //std::cout << "\nComponents = (";
        /*
        std::cout << Mi2 << "," << Mi2p << "," << this->_1p0 << ","
                  << mu2 << "," << this->_3p0 << "," << this->_m2p0 << ","
                  << ci << "," << logfactor << "," << this->_2p0 << ","
                  << ni << "," << pMi2pt << "," << H2 << ","
                  << npi << "," << pMi2ptp << "," << this->_1p0 << ","
                  << Mi2 << "," << Mi2p << "," << this->_1p0 << ")";*/
        //std::cout << Mi2;
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter this:
        res += (((this->_1p0)/(Mi2*Mi2))*(((this->_1p0)/(mu2))*((((this->_3p0
                + ((this->_m2p0*ci) + (this->_2p0*logfactor)))
                *(Mi2*Mi2*(Mi2p*(mu2*(ni*pMi2pt)))))
                + ((H2*H2*H2*H2*(Mi2*(mu2*(npi*pMi2ptp))))
                + (Mi2*Mi2*Mi2*(ni*((this->_m2p0*(mu2p*pMi2pt))
                + (mu2*((this->_m2p0*Mi2p) + ((this->_1p0 + ((this->_m2p0*ci)
                + (this->_2p0*logfactor)))*pMi2ptp))))))))
                - (H2*H2*H2*H2*(Mi2p*(mu2*(npi*pMi2pt)))))));
        //std::cout << res;
    }
    //std::cout << ");";
    return res;
}
//------------------------------------------------------------------------------
//Bpp(x) - dtdtcl_1loop_no_logs
template<class value_type>
value_type dtdtcl_1loop<value_type>::Bpp
(const value_type& mu,const value_type& phicl,
 const value_type& dtdtcl,const value_type& d2tdtcl2,
 const std::vector<value_type>& log_shared,
 const std::vector<value_type>& Mi2array,
 const std::vector<value_type>& Mi2parray,
 const std::vector<value_type>& Mi2pparray,
 const std::vector<value_type>& pMi2ptarray,
 const std::vector<value_type>& pMi2ptparray,
 const std::vector<value_type>& pMi2ptpparray)
{
    //Access shared data:
    value_type res = this->_0p0;
    //mu^2 and its (tcl) derivatives:
    value_type mu2 = mu*mu;
    value_type mu2p = mu2*dtdtcl;
    value_type mu2pp = (d2tdtcl2 + dtdtcl*dtdtcl)*mu2;
    value_type H2 = this->R/(this->_12p0);
    value_type Mi2,Mi2p,Mi2pp,ci,ni,npi;
    //ki,kip,kipp,kippp,thetaip,kpip,kpipp;
    //value_type kpippp,thetaipp,thetaippp;
    value_type pMi2pt,pMi2ptp,pMi2ptpp;
    //Compute Bpp:
    for(int i = 0;i < this->nLoops;i++)
    {
        Mi2 = Mi2array[i];
        Mi2p = Mi2parray[i];
        Mi2pp = Mi2pparray[i];
        pMi2pt = pMi2ptarray[i];
        pMi2ptp = pMi2ptparray[i];
        pMi2ptpp = pMi2ptpparray[i];
        ci = this->cilist[i];
        ni = this->nilist[i];
        npi = this->npilist[i];

        value_type logfactor = log_shared[i];
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter this:
        res += (((this->_1p0)/(Mi2*Mi2*Mi2))*(((this->_1p0)/(mu2*mu2))
                *(((this->_2p0*(H2*H2*H2*H2*(Mi2p*Mi2p*(mu2*mu2*(npi*pMi2pt)))))
                + ((Mi2*Mi2*Mi2*(mu2*(ni*((this->_m4p0*(Mi2p*(mu2p*pMi2pt)))
                + (mu2*((this->_m2p0*Mi2p*Mi2p)
                + (((this->_3p0 + ((this->_m2p0*ci) + (this->_2p0*logfactor)))
                *(Mi2pp*pMi2pt)) + (this->_2p0*((this->_3p0 + ((this->_m2p0*ci)
                + (this->_2p0*logfactor)))*(Mi2p*pMi2ptp))))))))))
                + ((Mi2*Mi2*(mu2*mu2*((this->_2p0*(Mi2p*Mi2p*(ni*pMi2pt)))
                + (H2*H2*H2*H2*(npi*pMi2ptpp)))))
                + (Mi2*Mi2*Mi2*Mi2*(ni*((this->_2p0*((mu2p*mu2p - (mu2*mu2pp))
                *pMi2pt)) + (mu2*((this->_m4p0*(mu2p*pMi2ptp))
                + (mu2*((this->_m2p0*Mi2pp) + ((this->_1p0 + ((this->_m2p0*ci)
                + (this->_2p0*logfactor)))*pMi2ptpp)))))))))))
                - (H2*H2*H2*H2*(Mi2*(mu2*mu2*(npi*((Mi2pp*pMi2pt)
                + (this->_2p0*(Mi2p*pMi2ptp))))))))));
    }
    return res;
}
//------------------------------------------------------------------------------
//Bppp(x) - dtdtcl_1loop_no_logs
template<class value_type>
value_type dtdtcl_1loop<value_type>::Bppp
(const value_type& mu,const value_type& phicl,
 const value_type& dtdtcl,const value_type& d2tdtcl2,
 const value_type& d3tdtcl3,
 const std::vector<value_type>& log_shared,
 const std::vector<value_type>& Mi2array,
 const std::vector<value_type>& Mi2parray,
 const std::vector<value_type>& Mi2pparray,
 const std::vector<value_type>& Mi2ppparray,
 const std::vector<value_type>& pMi2ptarray,
 const std::vector<value_type>& pMi2ptparray,
 const std::vector<value_type>& pMi2ptpparray,
 const std::vector<value_type>& pMi2ptppparray)
{
    //Access shared data:
    value_type res = this->_0p0;
    //mu^2 and its (tcl) derivatives:
    value_type mu2 = mu*mu;
    value_type mu2p = mu2*dtdtcl;
    value_type mu2pp = (d2tdtcl2 + dtdtcl*dtdtcl)*mu2;
    value_type mu2ppp = (d3tdtcl3
                         + this->_3p0*d2tdtcl2*dtdtcl
                         + dtdtcl*dtdtcl*dtdtcl)*mu2;
    value_type Mi2,Mi2p,Mi2pp,Mi2ppp,ci,ni,npi;
    //ki,kip,kipp,kippp,kipppp,thetaip;
    //value_type kpip,kpipp,kpippp,kpipppp;
    //value_type thetaipp,thetaippp,thetaipppp;
    value_type pMi2pt,pMi2ptp,pMi2ptpp,pMi2ptppp;
    value_type H2 = this->R/(this->_12p0);
    //Compute Bppp:
    for(int i = 0;i < this->nLoops;i++)
    {
        Mi2 = Mi2array[i];
        Mi2p = Mi2parray[i];
        Mi2pp = Mi2pparray[i];
        Mi2ppp = Mi2ppparray[i];
        pMi2pt = pMi2ptarray[i];
        pMi2ptp = pMi2ptparray[i];
        pMi2ptpp = pMi2ptpparray[i];
        pMi2ptppp = pMi2ptppparray[i];
        ci = this->cilist[i];
        ni = this->nilist[i];
        npi = this->npilist[i];

        value_type logfactor = log_shared[i];
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter this:
        res += (((this->_1p0)/(Mi2*Mi2*Mi2*Mi2))*(((this->_1p0)/(mu2*mu2*mu2))
                *((this->_m4p0*(Mi2*Mi2*Mi2*Mi2*Mi2*(mu2p*mu2p*mu2p
                *(ni*pMi2pt))))
                + ((this->_6p0*(Mi2*Mi2*Mi2*Mi2*(mu2*(mu2p*(ni*((((Mi2p*mu2p)
                + (Mi2*mu2pp))*pMi2pt) + (Mi2*(mu2p*pMi2ptp))))))))
                + ((this->_m2p0*(Mi2*Mi2*Mi2*Mi2*(mu2*mu2*(ni*((((this->_3p0
                *(Mi2pp*mu2p)) + (Mi2*mu2ppp))*pMi2pt)
                + ((this->_3p0*(Mi2p*((mu2pp*pMi2pt)
                + (this->_2p0*(mu2p*pMi2ptp)))))
                + (this->_3p0*(Mi2*((mu2pp*pMi2ptp) + (mu2p*pMi2ptpp))))))))))
                + (mu2*mu2*mu2*((((Mi2*Mi2*(Mi2ppp*(((this->_3p0
                + ((this->_m2p0*ci) + (this->_2p0*logfactor)))*(Mi2*Mi2*ni))
                - (H2*H2*H2*H2*npi))))
                + ((this->_6p0*(Mi2*(Mi2p*(Mi2pp*((Mi2*Mi2*ni)
                + (H2*H2*H2*H2*npi))))))
                + (this->_m2p0*(Mi2p*Mi2p*Mi2p*((Mi2*Mi2*ni)
                + (this->_3p0*(H2*H2*H2*H2*npi)))))))*pMi2pt)
                + (Mi2*((this->_6p0*(Mi2p*Mi2p*(((Mi2*Mi2*ni)
                + (H2*H2*H2*H2*npi))*pMi2ptp)))
                + ((this->_m3p0*(Mi2p*((H2*H2*H2*H2*(Mi2*(npi*pMi2ptpp)))
                + (Mi2*Mi2*Mi2*(ni*((this->_2p0*Mi2pp)
                + ((this->_m3p0 + ((this->_2p0*ci) + (this->_m2p0*logfactor)))
                *pMi2ptpp))))))) + (Mi2*((this->_m3p0*(Mi2pp*((((this->_m3p0
                + ((this->_2p0*ci) + (this->_m2p0*logfactor)))*(Mi2*Mi2*ni))
                + (H2*H2*H2*H2*npi))*pMi2ptp)))
                + ((H2*H2*H2*H2*(Mi2*(npi*pMi2ptppp)))
                + (Mi2*Mi2*Mi2*(ni*((this->_m2p0*Mi2ppp)
                + ((this->_1p0 + ((this->_m2p0*ci) + (this->_2p0*logfactor)))
                *pMi2ptppp)))))))))))))))));
    }
    return res;
}
//Bpppp(x) - dtdtcl_1loop_no_logs
template<class value_type>
value_type dtdtcl_1loop<value_type>::Bpppp
(const value_type& mu,const value_type& phicl,
 const value_type& dtdtcl,const value_type& d2tdtcl2,
 const value_type& d3tdtcl3,
 const value_type& d4tdtcl4,
 const std::vector<value_type>& log_shared,
 const std::vector<value_type>& Mi2array,
 const std::vector<value_type>& Mi2parray,
 const std::vector<value_type>& Mi2pparray,
 const std::vector<value_type>& Mi2ppparray,
 const std::vector<value_type>& Mi2pppparray,
 const std::vector<value_type>& pMi2ptarray,
 const std::vector<value_type>& pMi2ptparray,
 const std::vector<value_type>& pMi2ptpparray,
 const std::vector<value_type>& pMi2ptppparray,
 const std::vector<value_type>& pMi2ptpppparray)
{
    //Access shared data:
    value_type res = this->_0p0;
    //mu^2 and its (tcl) derivatives:
    value_type mu2 = mu*mu;
    value_type mu2p = mu2*dtdtcl;
    value_type mu2pp = (d2tdtcl2 + dtdtcl*dtdtcl)*mu2;
    value_type mu2ppp = (d3tdtcl3
                         + this->_3p0*d2tdtcl2*dtdtcl
                         + dtdtcl*dtdtcl*dtdtcl)*mu2;
    value_type mu2pppp = (d4tdtcl4 + this->_4p0*d3tdtcl3*dtdtcl
                         + this->_3p0*d2tdtcl2*d2tdtcl2
                         + this->_6p0*d2tdtcl2*dtdtcl*dtdtcl
                         + dtdtcl*dtdtcl*dtdtcl*dtdtcl)*mu2;
    value_type Mi2,Mi2p,Mi2pp,Mi2ppp,Mi2pppp,ci,ni,npi;
    value_type H2 = this->R/(this->_12p0);

    //ki,kip,kipp,kippp,kipppp,thetaip;
    //value_type kpip,kpipp,kpippp,kpipppp;
    //value_type thetaipp,thetaippp,thetaipppp;
    value_type pMi2pt,pMi2ptp,pMi2ptpp,pMi2ptppp,pMi2ptpppp;
    //Compute Bppp:
    for(int i = 0;i < this->nLoops;i++)
    {
        Mi2 = Mi2array[i];
        Mi2p = Mi2parray[i];
        Mi2pp = Mi2pparray[i];
        Mi2ppp = Mi2ppparray[i];
        Mi2pppp = Mi2pppparray[i];
        pMi2pt = pMi2ptarray[i];
        pMi2ptp = pMi2ptparray[i];
        pMi2ptpp = pMi2ptpparray[i];
        pMi2ptppp = pMi2ptppparray[i];
        pMi2ptpppp = pMi2ptpppparray[i];
        ci = this->cilist[i];
        ni = this->nilist[i];
        npi = this->npilist[i];

        value_type logfactor = log_shared[i];
        //NB - machine generated code! Do not debug by hand!
        //Use the mathematica script if it is necessary to alter this:
        res += (((this->_1p0)/(Mi2*Mi2*Mi2*Mi2*Mi2))
                *(((this->_1p0)/(mu2*mu2*mu2*mu2))
                *((this->_12p0*(Mi2*Mi2*Mi2*Mi2*Mi2*Mi2*(mu2p*mu2p*mu2p*mu2p
                *(ni*pMi2pt))))
                + ((this->_m8p0*(Mi2*Mi2*Mi2*Mi2*Mi2*(mu2*(mu2p*mu2p
                *(ni*((((this->_2p0*(Mi2p*mu2p))
                + (this->_3p0*(Mi2*mu2pp)))*pMi2pt)
                + (this->_2p0*(Mi2*(mu2p*pMi2ptp)))))))))
                + ((this->_2p0*(Mi2*Mi2*Mi2*Mi2*Mi2*(mu2*mu2*(ni*((((this->_6p0
                *(Mi2pp*mu2p*mu2p)) + ((this->_3p0*(Mi2*mu2pp*mu2pp))
                + (this->_4p0*(Mi2*(mu2p*mu2ppp)))))*pMi2pt)
                + ((this->_12p0*(Mi2p*(mu2p*((mu2pp*pMi2pt)
                + (mu2p*pMi2ptp)))))
                + (this->_6p0*(Mi2*(mu2p*((this->_2p0*(mu2pp*pMi2ptp))
                + (mu2p*pMi2ptpp)))))))))))
                + ((this->_m2p0*(Mi2*Mi2*Mi2*Mi2*Mi2*(mu2*mu2*mu2*(ni
                *((this->_6p0*(Mi2pp*(mu2pp*pMi2pt))) + ((Mi2*(mu2pppp*pMi2pt))
                + ((this->_4p0*(Mi2*(mu2ppp*pMi2ptp)))
                + ((this->_4p0*(Mi2p*((mu2ppp*pMi2pt)
                + (this->_3p0*(mu2pp*pMi2ptp)))))
                + ((this->_6p0*(Mi2*(mu2pp*pMi2ptpp)))
                + (this->_4p0*(mu2p*((Mi2ppp*pMi2pt)
                + ((this->_3p0*(Mi2pp*pMi2ptp)) + ((this->_3p0*(Mi2p*pMi2ptpp))
                + (Mi2*pMi2ptppp)))))))))))))))
                + (mu2*mu2*mu2*mu2*((((this->_8p0*(Mi2*Mi2*(Mi2p*(Mi2ppp
                *((Mi2*Mi2*ni) + (H2*H2*H2*H2*npi))))))
                + ((this->_m12p0*(Mi2*(Mi2p*Mi2p*(Mi2pp*((Mi2*Mi2*ni)
                + (this->_3p0*(H2*H2*H2*H2*npi)))))))
                + ((this->_4p0*(Mi2p*Mi2p*Mi2p*Mi2p*((Mi2*Mi2*ni)
                + (this->_6p0*(H2*H2*H2*H2*npi)))))
                + (Mi2*Mi2*((Mi2*(Mi2pppp*(((this->_3p0 + ((this->_m2p0*ci)
                + (this->_2p0*logfactor)))*(Mi2*Mi2*ni)) - (H2*H2*H2*H2*npi))))
                + (this->_6p0*(Mi2pp*Mi2pp*((Mi2*Mi2*ni)
                + (H2*H2*H2*H2*npi)))))))))*pMi2pt)
                + (Mi2*((this->_m8p0*(Mi2p*Mi2p*Mi2p*(((Mi2*Mi2*ni)
                + (this->_3p0*(H2*H2*H2*H2*npi)))*pMi2ptp)))
                + ((this->_12p0*(Mi2*(Mi2p*Mi2p*(((Mi2*Mi2*ni)
                + (H2*H2*H2*H2*npi))*pMi2ptpp))))
                + ((this->_m4p0*(Mi2*(Mi2p*((this->_m6p0*(Mi2pp*(((Mi2*Mi2*ni)
                + (H2*H2*H2*H2*npi))*pMi2ptp)))
                + ((H2*H2*H2*H2*(Mi2*(npi*pMi2ptppp)))
                + (Mi2*Mi2*Mi2*(ni*((this->_2p0*Mi2ppp)
                + ((this->_m3p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*pMi2ptppp)))))))))
                + (Mi2*Mi2*((this->_m2p0*(H2*H2*H2*H2*(npi*((this->_2p0
                *(Mi2ppp*pMi2ptp)) + (this->_3p0*(Mi2pp*pMi2ptpp))))))
                + ((this->_m2p0*(Mi2*Mi2*(ni*((this->_3p0*Mi2pp*Mi2pp)
                + ((this->_2p0*((this->_m3p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*(Mi2ppp*pMi2ptp)))
                + (this->_3p0*((this->_m3p0 + ((this->_2p0*ci)
                + (this->_m2p0*logfactor)))*(Mi2pp*pMi2ptpp))))))))
                + ((H2*H2*H2*H2*(Mi2*(npi*pMi2ptpppp)))
                + (Mi2*Mi2*Mi2*(ni*((this->_m2p0*Mi2pppp) + ((this->_1p0
                + ((this->_m2p0*ci)
                + (this->_2p0*logfactor)))*pMi2ptpppp))))))))))))))))))));
    }
    return res;
}
//These functions do the hard work of generating data
//for A,Ap,B,Bp etc..., so that we can share data if
//necessary and avoid having to repeat work.
template<class value_type>
void dtdtcl_1loop<value_type>::compute_dt_dtcl_data
(const value_type& tcl, const std::vector<value_type>& g,
 const std::vector<value_type>& dgdt,std::vector<value_type>& logfactor,
 std::vector<value_type>& Mi2,std::vector<value_type>& pMi2pt,
 std::vector<value_type>& kiarray,value_type& Z,value_type& phicl)
{
    phicl = this->mt*exp(tcl/this->_2p0);
    //Various coefficients needed for log contribution:
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);
    std::vector<value_type> kiparray;
    std::vector<value_type> kpiparray;
    std::vector<value_type> thetaiparray;
    this->populateCoefficientArraysDerivatives(kiparray,kpiparray,
                                               thetaiparray,dgdt);

    Z = this->eval.Z(g);
    value_type Zp = this->eval.Z(dgdt);
    value_type logmu2 = this->eval.t(g);
    /*
    std::cout << "g[i] = (";
    for(int i = 0;i<g.size();i++)
    {
        std::cout << g[i];
        if(i < g.size()-1)
        {
            std::cout << ",";
        }
    }
    std::cout << ")";*/
    this->computeMi2(Mi2,logfactor,pMi2pt,kiarray,kpiarray,thetaiarray,
                     kiparray,kpiparray,thetaiparray,phicl,Z,Zp,logmu2);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::dt_dtcl
(const value_type& tcl, const std::vector<value_type>& g,
                       const std::vector<value_type>& dgdt)
{
    std::vector<value_type> Mi2 (this->nLoops,this->_0p0);
    std::vector<value_type> logfactor (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2pt (this->nLoops,this->_0p0);
    std::vector<value_type> kiarray;
    value_type phicl,Z,logmu2;
    this->compute_dt_dtcl_data(tcl,g,dgdt,logfactor,Mi2,pMi2pt,kiarray,Z,phicl);
    return this->dt_dtcl(tcl,g,dgdt,logfactor,Mi2,pMi2pt,kiarray,Z,phicl);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::dt_dtcl
(const value_type& tcl, const std::vector<value_type>& g,
const std::vector<value_type>& dgdt,
const std::vector<value_type>& logfactor,
const std::vector<value_type>& Mi2,
const std::vector<value_type>& pMi2pt,
const std::vector<value_type>& kiarray,
const value_type Z,const value_type phicl)
{
    value_type a = this->A(logfactor,phicl,Mi2,kiarray,Z);
    value_type b = this->B(logfactor,Mi2,pMi2pt);
    return a/b;
}
template<class value_type>
void dtdtcl_1loop<value_type>::compute_d2t_dtcl2_data
(const value_type& tcl, const std::vector<value_type>& g,
 const std::vector<value_type>& dgdt,const std::vector<value_type>& d2gdt2,
 const value_type& dtdtcl,
 std::vector<value_type>& logfactor,std::vector<value_type>& Mi2,
 std::vector<value_type>& pMi2pt,std::vector<value_type>& Mi2p,
 std::vector<value_type>& pMi2ptp,std::vector<value_type>& kiarray,
 std::vector<value_type>& kiparray,value_type& Z,value_type& Zp,
 value_type& Zpp,value_type& mu,value_type& phicl)
{
    phicl = this->mt*exp(tcl/this->_2p0);
    mu = this->mt*exp(this->eval.t(g)/this->_2p0);
    //Various coefficients needed for log contribution:
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);
    std::vector<value_type> kpiparray;
    std::vector<value_type> thetaiparray;
    this->populateCoefficientArraysDerivatives(kiparray,kpiparray,
                                               thetaiparray,dgdt);
    std::vector<value_type> kipparray;
    std::vector<value_type> kpipparray;
    std::vector<value_type> thetaipparray;
    this->populateCoefficientArraysDerivatives(kipparray,kpipparray,
                                               thetaipparray,d2gdt2);
    Z = this->eval.Z(g);
    //t derivatives of Z:
    Zp = this->eval.Z(dgdt);
    Zpp = this->eval.Z(d2gdt2);

    value_type logmu2 = this->eval.t(g);
    this->computeMi2p(Mi2,logfactor,Mi2p,pMi2pt,pMi2ptp,kiarray,
                      kpiarray,thetaiarray,kiparray,kpiparray,thetaiparray,
                      kipparray,kpipparray,thetaipparray,
                      phicl,Z,Zp,Zpp,dtdtcl,logmu2);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::d2t_dtcl2
(const value_type& tcl, const std::vector<value_type>& g,
                   const std::vector<value_type>& dgdt,
                   const std::vector<value_type>& d2gdt2,
                   const value_type& dtdtcl)
{
    value_type phicl, mu, Z, Zp, Zpp;
    std::vector<value_type> kiarray, kiparray;
    std::vector<value_type> logfactor,Mi2,pMi2pt,Mi2p,pMi2ptp;

    this->compute_d2t_dtcl2_data(tcl,g,dgdt,d2gdt2,dtdtcl,logfactor,Mi2,
                                 pMi2pt,Mi2p,pMi2ptp,kiarray,kiparray,Z,Zp,Zpp,
                                 mu,phicl);
    return this->d2t_dtcl2(tcl,g,logfactor,Mi2,pMi2pt,Mi2p,pMi2ptp,kiarray,
                    kiparray,Z,Zp,Zpp,mu,phicl,dtdtcl);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::d2t_dtcl2
(const value_type& tcl, const std::vector<value_type>& g,
                   const std::vector<value_type>& logfactor,
                   const std::vector<value_type>& Mi2,
                   const std::vector<value_type>& pMi2pt,
                   const std::vector<value_type>& Mi2p,
                   const std::vector<value_type>& pMi2ptp,
                   const std::vector<value_type>& kiarray,
                   const std::vector<value_type>& kiparray,
                   const value_type& Z,const value_type& Zp,
                   const value_type& Zpp,const value_type& mu,
                   const value_type phicl,const value_type& dtdtcl)
{
    value_type a = this->A(logfactor,phicl,Mi2,kiarray,Z);
    value_type b = this->B(logfactor,Mi2,pMi2pt);
    value_type ap = this->Ap(Z,Zp,kiarray,kiparray,mu,logfactor,
                             phicl,Mi2,Mi2p,dtdtcl);
    value_type bp = this->Bp(dtdtcl,mu,phicl,logfactor,
                              Mi2,Mi2p,pMi2pt,pMi2ptp);
    //std::cout << "\nA = " << a << " Ap = " << ap << " B = " << b << " Bp = " << bp;
    return ap/b - a*bp/(b*b);
}
template<class value_type>
void dtdtcl_1loop<value_type>::compute_d3t_dtcl3_data
(const value_type& tcl, const std::vector<value_type>& g,
 const std::vector<value_type>& dgdt,const std::vector<value_type>& d2gdt2,
 const std::vector<value_type>& d3gdt3,
 std::vector<value_type>& logfactor,std::vector<value_type>& Mi2,
 std::vector<value_type>& pMi2pt,std::vector<value_type>& Mi2p,
 std::vector<value_type>& pMi2ptp,std::vector<value_type>& Mi2pp,
 std::vector<value_type>& pMi2ptpp,std::vector<value_type>& kiarray,
 std::vector<value_type>& kiparray,
 std::vector<value_type>& kipparray,value_type& Z,
 value_type& Zp,value_type& Zpp,value_type& Zppp,
 value_type& mu,value_type& phicl,const value_type& dtdtcl,
 const value_type& d2tdtcl2)
{
    phicl = this->mt*exp(tcl/this->_2p0);
    mu = this->mt*exp(this->eval.t(g)/this->_2p0);
    //Various coefficients needed for log contribution:
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);
    std::vector<value_type> kpiparray;
    std::vector<value_type> thetaiparray;
    this->populateCoefficientArraysDerivatives(kiparray,kpiparray,
                                               thetaiparray,dgdt);
    std::vector<value_type> kpipparray;
    std::vector<value_type> thetaipparray;
    this->populateCoefficientArraysDerivatives(kipparray,kpipparray,
                                               thetaipparray,d2gdt2);
    std::vector<value_type> kippparray;
    std::vector<value_type> kpippparray;
    std::vector<value_type> thetaippparray;
    this->populateCoefficientArraysDerivatives(kippparray,kpippparray,
                                               thetaippparray,d3gdt3);

    Z = this->eval.Z(g);
    Zp = this->eval.Z(dgdt);
    Zpp = this->eval.Z(d2gdt2);
    Zppp = this->eval.Z(d3gdt3);
    value_type logmu2 = this->eval.t(g);
    this->computeMi2pp(Mi2,logfactor,Mi2p,Mi2pp,pMi2pt,pMi2ptp,pMi2ptpp,kiarray,
                      kpiarray,thetaiarray,kiparray,kpiparray,thetaiparray,
                      kipparray,kpipparray,thetaipparray,
                      kippparray,kpippparray,thetaippparray,
                      phicl,Z,Zp,Zpp,Zppp,dtdtcl,d2tdtcl2,logmu2);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::d3t_dtcl3
(const value_type& tcl, const std::vector<value_type>& g,
 const std::vector<value_type>& dgdt,const std::vector<value_type>& d2gdt2,
 const std::vector<value_type>& d3gdt3,const value_type& dtdtcl,
 const value_type& d2tdtcl2)
{
    value_type phicl,mu,Z,Zp,Zpp,Zppp,logmu2;
    std::vector<value_type> kiarray,kiparray,kipparray;

    std::vector<value_type> Mi2 (this->nLoops,this->_0p0);
    std::vector<value_type> logfactor (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2pt (this->nLoops,this->_0p0);
    std::vector<value_type> Mi2p (this->nLoops,this->_0p0);
    std::vector<value_type> Mi2pp (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2ptp (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2ptpp (this->nLoops,this->_0p0);

    this->compute_d3t_dtcl3_data(tcl,g,dgdt,d2gdt2,d3gdt3,logfactor,Mi2,pMi2pt,
                                 Mi2p,pMi2ptp,Mi2pp,pMi2ptpp,kiarray,kiparray,
                                 kipparray,Z,Zp,Zpp,Zppp,mu,phicl,dtdtcl,
                                 d2tdtcl2);

    return this->d3t_dtcl3(tcl,g,logfactor,Mi2,pMi2pt,Mi2p,pMi2ptp,Mi2pp,
                    pMi2ptpp,kiarray,kiparray,kipparray,Z,Zp,Zpp,Zppp,mu,
                    phicl,dtdtcl,d2tdtcl2);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::d3t_dtcl3
(const value_type& t, const std::vector<value_type>& g,
                   const std::vector<value_type>& logfactor,
                   const std::vector<value_type>& Mi2,
                   const std::vector<value_type>& pMi2pt,
                   const std::vector<value_type>& Mi2p,
                   const std::vector<value_type>& pMi2ptp,
                   const std::vector<value_type>& Mi2pp,
                   const std::vector<value_type>& pMi2ptpp,
                   const std::vector<value_type>& kiarray,
                   const std::vector<value_type>& kiparray,
                   const std::vector<value_type>& kipparray,
                   const value_type& Z,const value_type& Zp,
                   const value_type& Zpp,const value_type& Zppp,
                   const value_type& mu,
                   const value_type phicl,const value_type& dtdtcl,
                   const value_type& d2tdtcl2)
{
    value_type a = this->A(logfactor,phicl,Mi2,kiarray,Z);
    value_type b = this->B(logfactor,Mi2,pMi2pt);
    value_type ap = this->Ap(Z,Zp,kiarray,kiparray,mu,logfactor,
                             phicl,Mi2,Mi2p,dtdtcl);
    value_type bp = this->Bp(dtdtcl,mu,phicl,logfactor,
                              Mi2,Mi2p,pMi2pt,pMi2ptp);
    value_type app = this->App(Z,Zp,Zpp,kiarray,kiparray,kipparray,
                               mu,logfactor,phicl,Mi2,Mi2p,Mi2pp,dtdtcl,
                               d2tdtcl2);
    value_type bpp = this->Bpp(dtdtcl,d2tdtcl2,mu,phicl,logfactor,
                               Mi2,Mi2p,Mi2pp,pMi2pt,pMi2ptp,pMi2ptpp);
    return app/b - this->_2p0*ap*bp/(b*b) - a*bpp/(b*b)
            + this->_2p0*a*bp*bp/(b*b*b);
}
template<class value_type>
void dtdtcl_1loop<value_type>::compute_d4t_dtcl4_data
(const value_type& tcl, const std::vector<value_type>& g,
 const std::vector<value_type>& dgdt,const std::vector<value_type>& d2gdt2,
 const std::vector<value_type>& d3gdt3,const std::vector<value_type>& d4gdt4,
 std::vector<value_type>& logfactor,std::vector<value_type>& Mi2,
 std::vector<value_type>& pMi2pt,std::vector<value_type>& Mi2p,
 std::vector<value_type>& pMi2ptp,std::vector<value_type>& Mi2pp,
 std::vector<value_type>& pMi2ptpp,std::vector<value_type>& Mi2ppp,
 std::vector<value_type>& pMi2ptppp,
 std::vector<value_type>& kiarray,
 std::vector<value_type>& kiparray,
 std::vector<value_type>& kipparray,
 std::vector<value_type>& kippparray,value_type& Z,
 value_type& Zp,value_type& Zpp,value_type& Zppp,value_type& Zpppp,
 value_type& mu,value_type& phicl,const value_type& dtdtcl,
 const value_type& d2tdtcl2,const value_type& d3tdtcl3)
{
    phicl = this->mt*exp(tcl/this->_2p0);
    mu = this->mt*exp(this->eval.t(g)/this->_2p0);
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);
    std::vector<value_type> kpiparray;
    std::vector<value_type> thetaiparray;
    this->populateCoefficientArraysDerivatives(kiparray,kpiparray,
                                               thetaiparray,dgdt);
    std::vector<value_type> kpipparray;
    std::vector<value_type> thetaipparray;
    this->populateCoefficientArraysDerivatives(kipparray,kpipparray,
                                               thetaipparray,d2gdt2);
    std::vector<value_type> kpippparray;
    std::vector<value_type> thetaippparray;
    this->populateCoefficientArraysDerivatives(kippparray,kpippparray,
                                               thetaippparray,d3gdt3);
    std::vector<value_type> kipppparray;
    std::vector<value_type> kpipppparray;
    std::vector<value_type> thetaipppparray;
    this->populateCoefficientArraysDerivatives(kipppparray,kpipppparray,
                                               thetaipppparray,d4gdt4);

    Z = this->eval.Z(g);
    Zp = this->eval.Z(dgdt);
    Zpp = this->eval.Z(d2gdt2);
    Zppp = this->eval.Z(d3gdt3);
    Zpppp = this->eval.Z(d4gdt4);
    value_type logmu2 = this->eval.t(g);
    this->computeMi2ppp(Mi2,logfactor,Mi2p,Mi2pp,Mi2ppp,pMi2pt,pMi2ptp,
                        pMi2ptpp,pMi2ptppp,kiarray,kpiarray,thetaiarray,
                        kiparray,kpiparray,thetaiparray,
                        kipparray,kpipparray,thetaipparray,
                        kippparray,kpippparray,thetaippparray,
                        kipppparray,kpipppparray,thetaipppparray,
                        phicl,Z,Zp,Zpp,Zppp,Zpppp,dtdtcl,d2tdtcl2,d3tdtcl3,
                        logmu2);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::d4t_dtcl4
(const value_type& tcl, const std::vector<value_type>& g,
                     const std::vector<value_type>& dgdt,
                     const std::vector<value_type>& d2gdt2,
                     const std::vector<value_type>& d3gdt3,
                     const std::vector<value_type>& d4gdt4,
                     const value_type& dtdtcl,const value_type& d2tdtcl2,
                     const value_type& d3tdtcl3)
{
    value_type phicl, mu,Z,Zp,Zpp,Zppp,Zpppp,logmu2;
    std::vector<value_type> kiarray,kiparray,kipparray,kippparray;

    //Various coefficients needed for log contribution:


    std::vector<value_type> Mi2 (this->nLoops,this->_0p0);
    std::vector<value_type> logfactor (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2pt (this->nLoops,this->_0p0);
    std::vector<value_type> Mi2p (this->nLoops,this->_0p0);
    std::vector<value_type> Mi2pp (this->nLoops,this->_0p0);
    std::vector<value_type> Mi2ppp (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2ptp (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2ptpp (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2ptppp (this->nLoops,this->_0p0);

    this->compute_d4t_dtcl4_data
    (tcl,g,dgdt,d2gdt2,d3gdt3,d4gdt4,logfactor,Mi2,pMi2pt,Mi2p,pMi2ptp,Mi2pp,
     pMi2ptpp,Mi2ppp,pMi2ptppp,kiarray,kiparray,kipparray,kippparray,
     Z,Zp,Zpp,Zppp,Zpppp,mu,phicl,dtdtcl,d2tdtcl2,d3tdtcl3);

    return this->d4t_dtcl4(tcl,g,logfactor,Mi2,pMi2pt,Mi2p,pMi2ptp,Mi2pp,
                    pMi2ptpp,Mi2ppp,pMi2ptppp,kiarray,kiparray,kipparray,
                    kippparray,Z,Zp,Zpp,Zppp,mu,
                    phicl,dtdtcl,d2tdtcl2,d3tdtcl3);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::d4t_dtcl4
(const value_type& t, const std::vector<value_type>& g,
                     const std::vector<value_type>& logfactor,
                     const std::vector<value_type>& Mi2,
                     const std::vector<value_type>& pMi2pt,
                     const std::vector<value_type>& Mi2p,
                     const std::vector<value_type>& pMi2ptp,
                     const std::vector<value_type>& Mi2pp,
                     const std::vector<value_type>& pMi2ptpp,
                     const std::vector<value_type>& Mi2ppp,
                     const std::vector<value_type>& pMi2ptppp,
                     const std::vector<value_type>& kiarray,
                     const std::vector<value_type>& kiparray,
                     const std::vector<value_type>& kipparray,
                     const std::vector<value_type>& kippparray,
                     const value_type& Z,const value_type& Zp,
                     const value_type& Zpp,const value_type& Zppp,
                     const value_type& mu,
                     const value_type phicl,const value_type& dtdtcl,
                     const value_type& d2tdtcl2,const value_type& d3tdtcl3)
{
    value_type a = this->A(logfactor,phicl,Mi2,kiarray,Z);
    value_type b = this->B(logfactor,Mi2,pMi2pt);
    value_type ap = this->Ap(Z,Zp,kiarray,kiparray,mu,logfactor,
                             phicl,Mi2,Mi2p,dtdtcl);
    value_type bp = this->Bp(dtdtcl,mu,phicl,logfactor,
                              Mi2,Mi2p,pMi2pt,pMi2ptp);
    value_type app = this->App(Z,Zp,Zpp,kiarray,kiparray,kipparray,
                               mu,logfactor,phicl,Mi2,Mi2p,Mi2pp,dtdtcl,
                               d2tdtcl2);
    value_type bpp = this->Bpp(dtdtcl,d2tdtcl2,mu,phicl,logfactor,
                              Mi2,Mi2p,Mi2pp,pMi2pt,pMi2ptp,pMi2ptpp);
    value_type appp = this->Appp(Z,Zp,Zpp,Zppp,kiarray,kiparray,
                                kipparray,kippparray,mu,logfactor,phicl,
                                Mi2,Mi2p,Mi2pp,Mi2ppp,
                                dtdtcl,d2tdtcl2,d3tdtcl3);
    value_type bppp = this->Bppp(dtdtcl,d2tdtcl2,d3tdtcl3,mu,phicl,
                                logfactor,Mi2,Mi2p,Mi2pp,Mi2ppp,
                                pMi2pt,pMi2ptp,pMi2ptpp,pMi2ptppp);
    return appp/b - this->_3p0*app*bp/(b*b)- this->_3p0*ap*bpp/(b*b)
            + this->_6p0*ap*bp*bp/(b*b*b) - a*bppp/(b*b)
            + this->_6p0*a*bpp*bp/(b*b*b)
            - this->_6p0*a*bp*bp*bp/(b*b*b*b);
}
template<class value_type>
void dtdtcl_1loop<value_type>::compute_d5t_dtcl5_data
(const value_type& tcl, const std::vector<value_type>& g,
 const std::vector<value_type>& dgdt,const std::vector<value_type>& d2gdt2,
 const std::vector<value_type>& d3gdt3,const std::vector<value_type>& d4gdt4,
 const std::vector<value_type>& d5gdt5,
 std::vector<value_type>& logfactor,std::vector<value_type>& Mi2,
 std::vector<value_type>& pMi2pt,std::vector<value_type>& Mi2p,
 std::vector<value_type>& pMi2ptp,std::vector<value_type>& Mi2pp,
 std::vector<value_type>& pMi2ptpp,std::vector<value_type>& Mi2ppp,
 std::vector<value_type>& pMi2ptppp,std::vector<value_type>& Mi2pppp,
 std::vector<value_type>& pMi2ptpppp,
 std::vector<value_type>& kiarray,
 std::vector<value_type>& kiparray,
 std::vector<value_type>& kipparray,
 std::vector<value_type>& kippparray,std::vector<value_type>& kipppparray,
 value_type& Z,
 value_type& Zp,value_type& Zpp,value_type& Zppp,value_type& Zpppp,
 value_type& Zppppp,
 value_type& mu,value_type& phicl,const value_type& dtdtcl,
 const value_type& d2tdtcl2,const value_type& d3tdtcl3,
 const value_type& d4tdtcl4)
{
    phicl = this->mt*exp(tcl/this->_2p0);
    mu = this->mt*exp(this->eval.t(g)/this->_2p0);
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);
    std::vector<value_type> kpiparray;
    std::vector<value_type> thetaiparray;
    this->populateCoefficientArraysDerivatives(kiparray,kpiparray,
                                               thetaiparray,dgdt);
    std::vector<value_type> kpipparray;
    std::vector<value_type> thetaipparray;
    this->populateCoefficientArraysDerivatives(kipparray,kpipparray,
                                               thetaipparray,d2gdt2);
    std::vector<value_type> kpippparray;
    std::vector<value_type> thetaippparray;
    this->populateCoefficientArraysDerivatives(kippparray,kpippparray,
                                               thetaippparray,d3gdt3);
    //std::vector<value_type> kipppparray;
    std::vector<value_type> kpipppparray;
    std::vector<value_type> thetaipppparray;
    this->populateCoefficientArraysDerivatives(kipppparray,kpipppparray,
                                               thetaipppparray,d4gdt4);
    std::vector<value_type> kippppparray;
    std::vector<value_type> kpippppparray;
    std::vector<value_type> thetaippppparray;
    this->populateCoefficientArraysDerivatives(kippppparray,kpippppparray,
                                               thetaippppparray,d5gdt5);

    Z = this->eval.Z(g);
    Zp = this->eval.Z(dgdt);
    Zpp = this->eval.Z(d2gdt2);
    Zppp = this->eval.Z(d3gdt3);
    Zpppp = this->eval.Z(d4gdt4);
    Zpppp = this->eval.Z(d5gdt5);
    value_type logmu2 = this->eval.t(g);
    this->computeMi2pppp(Mi2,logfactor,Mi2p,Mi2pp,Mi2ppp,Mi2pppp,pMi2pt,pMi2ptp,
                        pMi2ptpp,pMi2ptppp,pMi2ptpppp,
                        kiarray,kpiarray,thetaiarray,
                        kiparray,kpiparray,thetaiparray,
                        kipparray,kpipparray,thetaipparray,
                        kippparray,kpippparray,thetaippparray,
                        kipppparray,kpipppparray,thetaipppparray,
                        kippppparray,kpippppparray,thetaippppparray,
                        phicl,Z,Zp,Zpp,Zppp,Zpppp,Zppppp,
                        dtdtcl,d2tdtcl2,d3tdtcl3,d4tdtcl4,
                        logmu2);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::d5t_dtcl5
(const value_type& tcl, const std::vector<value_type>& g,
                     const std::vector<value_type>& dgdt,
                     const std::vector<value_type>& d2gdt2,
                     const std::vector<value_type>& d3gdt3,
                     const std::vector<value_type>& d4gdt4,
                     const std::vector<value_type>& d5gdt5,
                     const value_type& dtdtcl,const value_type& d2tdtcl2,
                     const value_type& d3tdtcl3,const value_type& d4tdtcl4)
{
    value_type phicl, mu,Z,Zp,Zpp,Zppp,Zpppp,Zppppp,logmu2;
    std::vector<value_type> kiarray,kiparray,kipparray,kippparray,kipppparray;

    //Various coefficients needed for log contribution:


    std::vector<value_type> Mi2 (this->nLoops,this->_0p0);
    std::vector<value_type> logfactor (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2pt (this->nLoops,this->_0p0);
    std::vector<value_type> Mi2p (this->nLoops,this->_0p0);
    std::vector<value_type> Mi2pp (this->nLoops,this->_0p0);
    std::vector<value_type> Mi2ppp (this->nLoops,this->_0p0);
    std::vector<value_type> Mi2pppp (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2ptp (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2ptpp (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2ptppp (this->nLoops,this->_0p0);
    std::vector<value_type> pMi2ptpppp (this->nLoops,this->_0p0);

    this->compute_d5t_dtcl5_data
    (tcl,g,dgdt,d2gdt2,d3gdt3,d4gdt4,d5gdt5,logfactor,Mi2,pMi2pt,Mi2p,pMi2ptp,
     Mi2pp,
     pMi2ptpp,Mi2ppp,pMi2ptppp,Mi2pppp,pMi2ptpppp,
     kiarray,kiparray,kipparray,kippparray,kipppparray,
     Z,Zp,Zpp,Zppp,Zpppp,Zppppp,mu,phicl,dtdtcl,d2tdtcl2,d3tdtcl3,d4tdtcl4);

    return this->d5t_dtcl5(tcl,g,logfactor,Mi2,pMi2pt,Mi2p,pMi2ptp,Mi2pp,
                    pMi2ptpp,Mi2ppp,pMi2ptppp,Mi2pppp,pMi2ptpppp,
                    kiarray,kiparray,kipparray,
                    kippparray,kipppparray,Z,Zp,Zpp,Zppp,Zpppp,mu,
                    phicl,dtdtcl,d2tdtcl2,d3tdtcl3,d4tdtcl4);
}
template<class value_type>
value_type dtdtcl_1loop<value_type>::d5t_dtcl5
(const value_type& t, const std::vector<value_type>& g,
                     const std::vector<value_type>& logfactor,
                     const std::vector<value_type>& Mi2,
                     const std::vector<value_type>& pMi2pt,
                     const std::vector<value_type>& Mi2p,
                     const std::vector<value_type>& pMi2ptp,
                     const std::vector<value_type>& Mi2pp,
                     const std::vector<value_type>& pMi2ptpp,
                     const std::vector<value_type>& Mi2ppp,
                     const std::vector<value_type>& pMi2ptppp,
                     const std::vector<value_type>& Mi2pppp,
                     const std::vector<value_type>& pMi2ptpppp,
                     const std::vector<value_type>& kiarray,
                     const std::vector<value_type>& kiparray,
                     const std::vector<value_type>& kipparray,
                     const std::vector<value_type>& kippparray,
                     const std::vector<value_type>& kipppparray,
                     const value_type& Z,const value_type& Zp,
                     const value_type& Zpp,const value_type& Zppp,
                     const value_type& Zpppp,const value_type& mu,
                     const value_type phicl,const value_type& dtdtcl,
                     const value_type& d2tdtcl2,const value_type& d3tdtcl3,
                     const value_type& d4tdtcl4)
{
    value_type a = this->A(logfactor,phicl,Mi2,kiarray,Z);
    value_type b = this->B(logfactor,Mi2,pMi2pt);
    value_type ap = this->Ap(Z,Zp,kiarray,kiparray,mu,logfactor,
                             phicl,Mi2,Mi2p,dtdtcl);
    value_type bp = this->Bp(dtdtcl,mu,phicl,logfactor,
                              Mi2,Mi2p,pMi2pt,pMi2ptp);
    value_type app = this->App(Z,Zp,Zpp,kiarray,kiparray,kipparray,
                               mu,logfactor,phicl,Mi2,Mi2p,Mi2pp,dtdtcl,
                               d2tdtcl2);
    value_type bpp = this->Bpp(dtdtcl,d2tdtcl2,mu,phicl,logfactor,
                              Mi2,Mi2p,Mi2pp,pMi2pt,pMi2ptp,pMi2ptpp);
    value_type appp = this->Appp(Z,Zp,Zpp,Zppp,kiarray,kiparray,
                                kipparray,kippparray,mu,logfactor,phicl,
                                Mi2,Mi2p,Mi2pp,Mi2ppp,
                                dtdtcl,d2tdtcl2,d3tdtcl3);
    value_type bppp = this->Bppp(dtdtcl,d2tdtcl2,d3tdtcl3,mu,phicl,
                                logfactor,Mi2,Mi2p,Mi2pp,Mi2ppp,
                                pMi2pt,pMi2ptp,pMi2ptpp,pMi2ptppp);
    value_type apppp = this->Apppp(Z,Zp,Zpp,Zppp,Zpppp,kiarray,kiparray,
                                kipparray,kippparray,kipppparray,mu,logfactor,
                                phicl,Mi2,Mi2p,Mi2pp,Mi2ppp,Mi2pppp,
                                dtdtcl,d2tdtcl2,d3tdtcl3,d4tdtcl4);
    value_type bpppp = this->Bpppp(dtdtcl,d2tdtcl2,d3tdtcl3,d4tdtcl4,
                                   mu,phicl,logfactor,Mi2,Mi2p,Mi2pp,Mi2ppp,
                                   Mi2pppp,pMi2pt,pMi2ptp,pMi2ptpp,pMi2ptppp,
                                   pMi2ptpppp);
    return -this->_24p0*ap*bp*bp*bp/(b*b*b*b)
                +this->_24p0*a*bp*bp*bp*bp/(b*b*b*b*b)
                +this->_12p0*bp*bp*app/(b*b*b)
                +this->_24p0*ap*bp*bpp/(b*b*b)
                -this->_36p0*a*bp*bp*bpp/(b*b*b*b)
                -this->_6p0*app*bpp/(b*b)
                + this->_6p0*a*bpp*bpp/(b*b*b)
                - this->_4p0*bp*appp/(b*b)
                - this->_4p0*ap*bppp/(b*b)
                + this->_8p0*a*bp*bppp/(b*b*b)
                +apppp/b - a*bpppp/(b*b);
}
//------------------------------------------------------------------------------
//==============================================================================
//------------------------------------------------------------------------------
//Standard Model Potential:
//Reference list of couplings:
    /*
    BASIC:
    g[0]=lambda
    g[1]=g_1^2
    g[2]=g_2^2
    g[3]=g_3^2
    g[4]=m^2
    g[5]=y_e^2
    g[6]=y_mu^2
    g[7]=y_tau^2
    g[8]=y_u^2
    g[9]=y_d^2
    g[10]=y_c^2
    g[11]=y_s^2
    g[12]=y_t^2
    g[13]=y_b^2
    EXTENDED:
    g[14]=Z
    g[15] = log(phi^2/mt^2) = t_cl
    g[16] = log(mu^2/mt^2) = t
    GRAVITATIONAL:
    g[17]=xi
    g[18]=V_0
    g[19]=kappa
    g[20]=alpha_1
    g[21]=alpha_2
    g[22]=alpha_3

    OLD SYSTEM - DON'T USE!!!!!!
    //Second argument, [i], refers to the element of the computed grid. First
    //argument determined the coupling.
    ydata[0][i] - lambda (Higgs self coupling)
    ydata[1][i] - g1^2 (U(1) coupling)
    ydata[2][i] - g2^2 (SU(2) coupling)
    ydata[3][i] - g3^2 (SU(3) coupling)
    ydata[4][i] - yt^2 (Top quark Yukawa coupling)
    ydata[5][i] - m^2 (Higgs tachyonic mass term)
    ydata[6][i] - xi (Higgs-curvature coupling)
    ydata[7][i] - Z (Higgs field renormalisation) - technically the RHS of this
                  is a gamma function, but we solve it with the beta functions)
    //Additionally, we have to solve for the non-beta functions:
    ydata[8][i] = t_cl = log(phicl^2/mt^2)
    ydata[9][i] = t = log(mu^2/mt^2)
    //Cosmological constant, if we include it:
    ydata[10][i] = Lambda*Mp^2 (cosmological constant, in units of the Planck
                                mass)
    */
//Compute log mass^2 terms that appear in the 1-loop potential. Keep
//separate so that we can compute these once if V needs to be called
//multiple times.
template<class value_type>
void SMpotential_dS_1loop<value_type>::populateCoefficientArrays
(std::vector<value_type>& k,std::vector<value_type>& kp,
 std::vector<value_type>& theta,const std::vector<value_type>& g)
{
    this->coeff.populateCoefficientArrays(k,kp,theta,g);
}
template<class value_type>
void SMpotential_dS_1loop<value_type>::populateCoefficientArraysDerivatives
(std::vector<value_type>& k,std::vector<value_type>& kp,
 std::vector<value_type>& theta,const std::vector<value_type>& gp)
{
    this->coeff.populateCoefficientArraysDerivatives(k,kp,theta,gp);
}
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2
(value_type y,std::vector<value_type>& g,std::vector<value_type>& Mi2,
 std::vector<value_type>& logMi2)
{
    //Various coefficients needed for log contribution:
    std::vector<value_type> kiarray;
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);

    this->computeMi2(y,g,Mi2,logMi2,kiarray,kpiarray,thetaiarray);
}
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2
(value_type y,std::vector<value_type>& g,std::vector<value_type>& Mi2,
 std::vector<value_type>& logMi2,
 std::vector<value_type>& kiarray,std::vector<value_type>& kpiarray,
 std::vector<value_type>& thetaiarray)
{
    value_type Z = this->eval.Z(g);
    value_type phi = Z*y;
    Mi2.assign(this->nLoops,_0p0);
    logMi2.assign(this->nLoops,_0p0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    value_type Mi2val;
    for(int i = 0;i < this->nLoops;i++)
    {
        Mi2val = kiarray[i]*phi*phi - kpiarray[i] + thetaiarray[i]*this->R;
        Mi2[i] = Mi2val == this->_0p0 ? eps : Mi2val;
        //We have problems with limits like x^2log(x) as x->0 if we allow Mi2
        //to vanish. To cure this, set Mi2 to the lowest number
        //allowed if this happens.


        //If Mi2 vanishes,the log will diverge too, which may cause a problem!
        //However, it always appears as (Mi2)^2log(Mi2), (or H^4log(Mi2)
        //for which the limit
        //is zero. To avoid this problem, just set the log to zero
        logMi2[i] = log(abs(Mi2[i])/(this->mt*this->mt)) - this->eval.t(g);
        Mi2[i] = Mi2val;
    }
}
//1st Derivatives of Mi2. Note, derivatives MUST be d^ng/dphicl^n,
//NOT d^ng/dtcl^n, which is what we store in the potential. Don't confuse them!
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2p
(value_type y, std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& Mi2,std::vector<value_type>& logMi2,
 std::vector<value_type>& Mi2p,std::vector<value_type>& logMi2p)
{
    //Various coefficients needed for log contribution:
    std::vector<value_type> kiarray;
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);
    std::vector<value_type> kiparray;
    std::vector<value_type> kpiparray;
    std::vector<value_type> thetaiparray;
    this->populateCoefficientArraysDerivatives
            (kiparray,kpiparray,thetaiparray,gp);

    this->computeMi2p(y,g,gp,Mi2,logMi2,Mi2p,logMi2p,
                      kiarray,kpiarray,thetaiarray,
                      kiparray,kpiparray,thetaiparray);

}
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2p
(value_type y, std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& Mi2,std::vector<value_type>& logMi2,
 std::vector<value_type>& Mi2p,std::vector<value_type>& logMi2p,
 std::vector<value_type>& kiarray,std::vector<value_type>& kpiarray,
 std::vector<value_type>& thetaiarray,std::vector<value_type>& kiparray,
 std::vector<value_type>& kpiparray,std::vector<value_type>& thetaiparray)
{
    value_type Z = this->eval.Z(g);
    value_type Zp = this->eval.Z(gp);
    value_type phip = Zp*y + Z;
    value_type phi = Z*y;

    Mi2.assign(this->nLoops,_0p0);
    logMi2.assign(this->nLoops,_0p0);
    Mi2p.assign(this->nLoops,_0p0);
    logMi2p.assign(this->nLoops,_0p0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    value_type Mi2val;
    for(int i = 0;i < this->nLoops;i++)
    {
        Mi2val = kiarray[i]*phi*phi - kpiarray[i] + thetaiarray[i]*this->R;
        Mi2[i] = Mi2val == this->_0p0 ? eps : Mi2val;
        //logMi2[i] = log(Mi2[i]/(this->mt*this->mt)) - this->eval.t(g);
        logMi2[i] = log(abs(Mi2[i])/(this->mt*this->mt)) - this->eval.t(g);

        Mi2p[i] = kiparray[i]*phi*phi + this->_2p0*kiarray[i]*phi*phip
                  - kpiparray[i] + thetaiparray[i]*this->R;
        logMi2p[i] = Mi2p[i]/Mi2[i] - this->eval.t(gp);
        Mi2[i] = Mi2val;
    }
}
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2pp
(value_type y, std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,
 std::vector<value_type>& Mi2,std::vector<value_type>& logMi2,
 std::vector<value_type>& Mi2p,std::vector<value_type>& logMi2p,
 std::vector<value_type>& Mi2pp,std::vector<value_type>& logMi2pp)
{
    //Various coefficients needed for log contribution:
    std::vector<value_type> kiarray;
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);
    std::vector<value_type> kiparray;
    std::vector<value_type> kpiparray;
    std::vector<value_type> thetaiparray;
    this->populateCoefficientArraysDerivatives
            (kiparray,kpiparray,thetaiparray,gp);
    std::vector<value_type> kipparray;
    std::vector<value_type> kpipparray;
    std::vector<value_type> thetaipparray;
    this->populateCoefficientArraysDerivatives
            (kipparray,kpipparray,thetaipparray,gpp);


    this->computeMi2pp(y,g,gp,gpp,Mi2,logMi2,Mi2p,logMi2p,
                      Mi2pp,logMi2pp,
                      kiarray,kpiarray,thetaiarray,
                      kiparray,kpiparray,thetaiparray,
                      kipparray,kpipparray,thetaipparray);
}
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2pp
(value_type y, std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,
 std::vector<value_type>& Mi2,std::vector<value_type>& logMi2,
 std::vector<value_type>& Mi2p,std::vector<value_type>& logMi2p,
 std::vector<value_type>& Mi2pp,std::vector<value_type>& logMi2pp,
 std::vector<value_type>& kiarray,std::vector<value_type>& kpiarray,
 std::vector<value_type>& thetaiarray,std::vector<value_type>& kiparray,
 std::vector<value_type>& kpiparray,std::vector<value_type>& thetaiparray,
 std::vector<value_type>& kipparray,
 std::vector<value_type>& kpipparray,std::vector<value_type>& thetaipparray)
{
    value_type Z = this->eval.Z(g);
    value_type Zp = this->eval.Z(gp);
    value_type Zpp = this->eval.Z(gpp);
    value_type phi = Z*y;
    value_type phip = Zp*y + Z;
    value_type phipp = Zpp*y + this->_2p0*Zp;

    Mi2.assign(this->nLoops,_0p0);
    logMi2.assign(this->nLoops,_0p0);
    Mi2p.assign(this->nLoops,_0p0);
    logMi2p.assign(this->nLoops,_0p0);
    Mi2pp.assign(this->nLoops,_0p0);
    logMi2pp.assign(this->nLoops,_0p0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    value_type Mi2val;
    for(int i = 0;i < this->nLoops;i++)
    {
        Mi2val = kiarray[i]*phi*phi - kpiarray[i] + thetaiarray[i]*this->R;
        Mi2[i] = Mi2val == this->_0p0 ? eps : Mi2val;
        //logMi2[i] = log(Mi2[i]/(this->mt*this->mt)) - g[9];
        logMi2[i] = log(abs(Mi2[i])/(this->mt*this->mt)) - this->eval.t(g);

        Mi2p[i] = kiparray[i]*phi*phi + this->_2p0*kiarray[i]*phi*phip
                  - kpiparray[i] + thetaiparray[i]*this->R;
        logMi2p[i] = Mi2p[i]/Mi2[i] - this->eval.t(gp);
        Mi2pp[i] = this->_4p0*kiparray[i]*phi*phip
                    + this->_2p0*kiarray[i]*phip*phip + phi*phi*kipparray[i]
                    -kpipparray[i] + this->_2p0*kiarray[i]*phi*phipp
                    + this->R*thetaipparray[i];
        logMi2pp[i] = Mi2pp[i]/Mi2[i] - Mi2p[i]*Mi2p[i]/(Mi2[i]*Mi2[i])
                        - this->eval.t(gpp);
        Mi2[i] = Mi2val;
    }
}
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2ppp
(value_type y, std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,std::vector<value_type>& gppp,
 std::vector<value_type>& Mi2,std::vector<value_type>& logMi2,
 std::vector<value_type>& Mi2p,std::vector<value_type>& logMi2p,
 std::vector<value_type>& Mi2pp,std::vector<value_type>& logMi2pp,
 std::vector<value_type>& Mi2ppp,std::vector<value_type>& logMi2ppp)
{
    //Various coefficients needed for log contribution:
    std::vector<value_type> kiarray;
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);
    std::vector<value_type> kiparray;
    std::vector<value_type> kpiparray;
    std::vector<value_type> thetaiparray;
    this->populateCoefficientArraysDerivatives
            (kiparray,kpiparray,thetaiparray,gp);
    std::vector<value_type> kipparray;
    std::vector<value_type> kpipparray;
    std::vector<value_type> thetaipparray;
    this->populateCoefficientArraysDerivatives
            (kipparray,kpipparray,thetaipparray,gpp);
    std::vector<value_type> kippparray;
    std::vector<value_type> kpippparray;
    std::vector<value_type> thetaippparray;
    this->populateCoefficientArraysDerivatives
            (kippparray,kpippparray,thetaippparray,gppp);


    this->computeMi2ppp(y,g,gp,gpp,gppp,Mi2,logMi2,Mi2p,logMi2p,
                      Mi2pp,logMi2pp,Mi2ppp,logMi2ppp,
                      kiarray,kpiarray,thetaiarray,
                      kiparray,kpiparray,thetaiparray,
                      kipparray,kpipparray,thetaipparray,
                      kippparray,kpippparray,thetaippparray);
}
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2ppp
(value_type y, std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,std::vector<value_type>& gppp,
 std::vector<value_type>& Mi2,std::vector<value_type>& logMi2,
 std::vector<value_type>& Mi2p,std::vector<value_type>& logMi2p,
 std::vector<value_type>& Mi2pp,std::vector<value_type>& logMi2pp,
 std::vector<value_type>& Mi2ppp,std::vector<value_type>& logMi2ppp,
 std::vector<value_type>& kiarray,std::vector<value_type>& kpiarray,
 std::vector<value_type>& thetaiarray,std::vector<value_type>& kiparray,
 std::vector<value_type>& kpiparray,std::vector<value_type>& thetaiparray,
 std::vector<value_type>& kipparray,
 std::vector<value_type>& kpipparray,std::vector<value_type>& thetaipparray,
 std::vector<value_type>& kippparray,
 std::vector<value_type>& kpippparray,std::vector<value_type>& thetaippparray)
{
    value_type Z = this->eval.Z(g);
    value_type Zp = this->eval.Z(gp);
    value_type Zpp = this->eval.Z(gpp);
    value_type Zppp = this->eval.Z(gppp);
    value_type phi = Z*y;
    value_type phip = Zp*y + Z;
    value_type phipp = Zpp*y + this->_2p0*Zp;
    value_type phippp = Zppp*y + this->_3p0*Zpp;

    Mi2.assign(this->nLoops,_0p0);
    logMi2.assign(this->nLoops,_0p0);
    Mi2p.assign(this->nLoops,_0p0);
    logMi2p.assign(this->nLoops,_0p0);
    Mi2pp.assign(this->nLoops,_0p0);
    logMi2pp.assign(this->nLoops,_0p0);
    Mi2ppp.assign(this->nLoops,_0p0);
    logMi2ppp.assign(this->nLoops,_0p0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    value_type Mi2val;
    for(int i = 0;i < this->nLoops;i++)
    {
        Mi2val = kiarray[i]*phi*phi - kpiarray[i] + thetaiarray[i]*this->R;
        Mi2[i] = Mi2val == this->_0p0 ? eps : Mi2val;
        //logMi2[i] = log(Mi2[i]/(this->mt*this->mt)) - this->eval.t(g);
        logMi2[i] = log(abs(Mi2[i])/(this->mt*this->mt)) - this->eval.t(g);

        Mi2p[i] = kiparray[i]*phi*phi + this->_2p0*kiarray[i]*phi*phip
                  - kpiparray[i] + thetaiparray[i]*this->R;
        logMi2p[i] = Mi2p[i]/Mi2[i] - this->eval.t(gp);
        Mi2pp[i] = this->_4p0*kiparray[i]*phi*phip
                    + this->_2p0*kiarray[i]*phip*phip + phi*phi*kipparray[i]
                    -kpipparray[i] + this->_2p0*kiarray[i]*phi*phipp
                    + this->R*thetaipparray[i];
        logMi2pp[i] = Mi2pp[i]/Mi2[i] - Mi2p[i]*Mi2p[i]/(Mi2[i]*Mi2[i])
                        - this->eval.t(gpp);
        Mi2ppp[i] = this->_6p0*kiparray[i]*phip*phip
                    + this->_6p0*phi*phip*kipparray[i]
                    + this->_6p0*phi*kiparray[i]*phipp
                    + this->_6p0*kiarray[i]*phip*phipp + phi*phi*kippparray[i]
                    - kpippparray[i] + this->_2p0*kiarray[i]*phi*phippp
                    +this->R*thetaippparray[i];
        logMi2ppp[i] = Mi2ppp[i]/Mi2[i]
                    - this->_3p0*Mi2pp[i]*Mi2p[i]/(Mi2[i]*Mi2[i])
                    + this->_2p0*Mi2p[i]*Mi2p[i]*Mi2p[i]/(Mi2[i]*Mi2[i]*Mi2[i])
                        - this->eval.t(gppp);
        Mi2[i] = Mi2val;
    }
}
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2pppp
(value_type y, std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,std::vector<value_type> gppp,
 std::vector<value_type> gpppp,
 std::vector<value_type>& Mi2,std::vector<value_type>& logMi2,
 std::vector<value_type>& Mi2p,std::vector<value_type>& logMi2p,
 std::vector<value_type>& Mi2pp,std::vector<value_type>& logMi2pp,
 std::vector<value_type>& Mi2ppp,std::vector<value_type>& logMi2ppp,
 std::vector<value_type>& Mi2pppp,std::vector<value_type>& logMi2pppp)
{
    //Various coefficients needed for log contribution:
    std::vector<value_type> kiarray;
    std::vector<value_type> kpiarray;
    std::vector<value_type> thetaiarray;
    this->populateCoefficientArrays(kiarray,kpiarray,thetaiarray,g);
    std::vector<value_type> kiparray;
    std::vector<value_type> kpiparray;
    std::vector<value_type> thetaiparray;
    this->populateCoefficientArraysDerivatives
            (kiparray,kpiparray,thetaiparray,gp);
    std::vector<value_type> kipparray;
    std::vector<value_type> kpipparray;
    std::vector<value_type> thetaipparray;
    this->populateCoefficientArraysDerivatives
            (kipparray,kpipparray,thetaipparray,gpp);
    std::vector<value_type> kippparray;
    std::vector<value_type> kpippparray;
    std::vector<value_type> thetaippparray;
    this->populateCoefficientArraysDerivatives
            (kippparray,kpippparray,thetaippparray,gppp);
    std::vector<value_type> kipppparray;
    std::vector<value_type> kpipppparray;
    std::vector<value_type> thetaipppparray;
    this->populateCoefficientArraysDerivatives
            (kipppparray,kpipppparray,thetaipppparray,gppp);


    this->computeMi2pppp(y,g,gp,gpp,gppp,gpppp,Mi2,logMi2,Mi2p,logMi2p,
                      Mi2pp,logMi2pp,Mi2ppp,logMi2ppp,
                      Mi2pppp,logMi2pppp,
                      kiarray,kpiarray,thetaiarray,
                      kiparray,kpiparray,thetaiparray,
                      kipparray,kpipparray,thetaipparray,
                      kippparray,kpippparray,thetaippparray,
                      kipppparray,kpipppparray,thetaipppparray);
}
template<class value_type>
void SMpotential_dS_1loop<value_type>::computeMi2pppp
(value_type y, std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,std::vector<value_type>& gppp,
 std::vector<value_type>& gpppp,
 std::vector<value_type>& Mi2,std::vector<value_type>& logMi2,
 std::vector<value_type>& Mi2p,std::vector<value_type>& logMi2p,
 std::vector<value_type>& Mi2pp,std::vector<value_type>& logMi2pp,
 std::vector<value_type>& Mi2ppp,std::vector<value_type>& logMi2ppp,
 std::vector<value_type>& Mi2pppp,std::vector<value_type>& logMi2pppp,
 std::vector<value_type>& kiarray,std::vector<value_type>& kpiarray,
 std::vector<value_type>& thetaiarray,std::vector<value_type>& kiparray,
 std::vector<value_type>& kpiparray,std::vector<value_type>& thetaiparray,
 std::vector<value_type>& kipparray,
 std::vector<value_type>& kpipparray,std::vector<value_type>& thetaipparray,
 std::vector<value_type>& kippparray,
 std::vector<value_type>& kpippparray,std::vector<value_type>& thetaippparray,
 std::vector<value_type>& kipppparray,
 std::vector<value_type>& kpipppparray,std::vector<value_type>& thetaipppparray)
{
    value_type Z = this->eval.Z(g);
    value_type Zp = this->eval.Z(gp);
    value_type Zpp = this->eval.Z(gpp);
    value_type Zppp = this->eval.Z(gppp);
    value_type Zpppp = this->eval.Z(gpppp);
    value_type phi = Z*y;
    value_type phip = Zp*y + Z;
    value_type phipp = Zpp*y + this->_2p0*Zp;
    value_type phippp = Zppp*y + this->_3p0*Zpp;
    value_type phipppp = Zpppp*y + this->_4p0*Zppp;

    Mi2.assign(this->nLoops,_0p0);
    logMi2.assign(this->nLoops,_0p0);
    Mi2p.assign(this->nLoops,_0p0);
    logMi2p.assign(this->nLoops,_0p0);
    Mi2pp.assign(this->nLoops,_0p0);
    logMi2pp.assign(this->nLoops,_0p0);
    Mi2ppp.assign(this->nLoops,_0p0);
    logMi2ppp.assign(this->nLoops,_0p0);
    Mi2pppp.assign(this->nLoops,_0p0);
    logMi2pppp.assign(this->nLoops,_0p0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    value_type Mi2val;
    for(int i = 0;i < this->nLoops;i++)
    {
        Mi2val = kiarray[i]*phi*phi - kpiarray[i] + thetaiarray[i]*this->R;
        Mi2[i] = Mi2val == this->_0p0 ? eps : Mi2val;
        //logMi2[i] = log(Mi2[i]/(this->mt*this->mt)) - this->eval.t(g);
        logMi2[i] = log(abs(Mi2[i])/(this->mt*this->mt)) - this->eval.t(g);

        Mi2p[i] = kiparray[i]*phi*phi + this->_2p0*kiarray[i]*phi*phip
                  - kpiparray[i] + thetaiparray[i]*this->R;
        logMi2p[i] = Mi2p[i]/Mi2[i] - this->eval.t(gp);
        Mi2pp[i] = this->_4p0*kiparray[i]*phi*phip
                    + this->_2p0*kiarray[i]*phip*phip + phi*phi*kipparray[i]
                    -kpipparray[i] + this->_2p0*kiarray[i]*phi*phipp
                    + this->R*thetaipparray[i];
        logMi2pp[i] = Mi2pp[i]/Mi2[i] - Mi2p[i]*Mi2p[i]/(Mi2[i]*Mi2[i])
                        - this->eval.t(gpp);
        Mi2ppp[i] = this->_6p0*kiparray[i]*phip*phip
                    + this->_6p0*phi*phip*kipparray[i]
                    + this->_6p0*phi*kiparray[i]*phipp
                    + this->_6p0*kiarray[i]*phip*phipp + phi*phi*kippparray[i]
                    - kpippparray[i] + this->_2p0*kiarray[i]*phi*phippp
                    +this->R*thetaippparray[i];
        logMi2ppp[i] = Mi2ppp[i]/Mi2[i]
                    - this->_3p0*Mi2pp[i]*Mi2p[i]/(Mi2[i]*Mi2[i])
                    + this->_2p0*Mi2p[i]*Mi2p[i]*Mi2p[i]/(Mi2[i]*Mi2[i]*Mi2[i])
                        - this->eval.t(gppp);
        Mi2pppp[i] = this->_12p0*phip*phip*kipparray[i]
                    + this->_24p0*kiparray[i]*phip*phipp
                    + this->_12p0*phi*kipparray[i]*phipp
                    + this->_6p0*kiarray[i]*phipp*phipp
                    + this->_8p0*phi*phip*kippparray[i]
                    + this->_8p0*phi*kiparray[i]*phippp
                    + this->_8p0*kiarray[i]*phip*phippp
                    + phi*phi*kipppparray[i] - kpipppparray[i]
                    + this->_2p0*kiarray[i]*phi*phipppp
                    + this->R*thetaipppparray[i];
        logMi2pppp[i] = -this->_6p0*Mi2p[i]*Mi2p[i]*Mi2p[i]*Mi2p[i]
                        /(Mi2[i]*Mi2[i]*Mi2[i]*Mi2[i])
                        + this->_12p0*Mi2p[i]*Mi2p[i]*Mi2pp[i]
                        /(Mi2[i]*Mi2[i]*Mi2[i])
                        -this->_3p0*Mi2pp[i]*Mi2pp[i]/(Mi2[i]*Mi2[i])
                        -this->_4p0*Mi2p[i]*Mi2ppp[i]/(Mi2[i]*Mi2[i])
                        + Mi2pppp[i]/Mi2[i] - this->eval.t(gpppp);
        Mi2[i] = Mi2val;
    }
}
//Jordan frame potential, in terms of Jordan frame classical field y,
//and the couplings at the chosen scale.
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::V
(value_type y,std::vector<value_type>& g,bool include_xi,bool include_grav)
{
    //Compute log terms:
    std::vector<value_type> Mi2;
    std::vector<value_type> logMi2;
    this->computeMi2(y,g,Mi2,logMi2);

    //Compute potential:
    return this->V(y,g,logMi2,Mi2,include_xi,include_grav);
}
//Same as the above, but uses a pre-computed version of log[Mi^2/mt^2], in
//order to save work if multiple functions need this.
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::V
(value_type y,std::vector<value_type>& g,std::vector<value_type>& logMi2,
 std::vector<value_type>& Mi2v,bool include_xi,bool include_grav)
{
    //Field renormalisation:
    value_type Z = this->eval.Z(g);
    //std::cout << "\nZ = " << Z;
    //Renormalised field:
    value_type phi2 = y*y*Z*Z;
    //std::cout << "\nphi2 = " << phi2;

    //Non-minimal coupling. If we are evaluating the Einstein-frame
    //potential, then we want to exclude the xi*R*phi^2/2 term, which is
    //conveniently done by setting xi = 0. Note that this doesn't apply to
    //the log terms! These are still xi dependent, but receive this via
    //the Mi2 functions, which are computed with the correct (non-zero in
    //general) xi.
    value_type xi = include_xi ? this->eval.xi(g) : this->_0p0;
    //std::cout << "\nxi = " << xi;
    //Check whether we want to include the cosmological constant:
    //value_type V0 = include_grav ? g[10] : this->_0p0;

    //Classical potential part:
    value_type v =  -this->eval.m2(g)*phi2/this->_2p0 + xi*this->R*phi2/this->_2p0
                    + this->eval.l(g)*phi2*phi2/this->_4p0;
    /*
    std::cout << "\nm2 = " << g[5];
    std::cout << "\nlambda = " << g[0];
    std::cout << "\nV_{cl} = " << v;
    */

    if(this->include_loops)
    {
        for(int i = 0;i < this->nLoops;i++)
        {
            /*
            std::cout << "\nni = " << this->ni[i];
            std::cout << "\nMi2v = " << Mi2v[i];
            std::cout << "\nlogMi2 = " << logMi2[i];
            std::cout << "\nci = " << this->ci[i];*/
            //log mass term for loop corrections:
            //value_type Mi2 = kiarray[i]*phi2 - kpiarray[i]
            //                    + thetaiarray[i]*this->R;
            //loop corrections:
            //These should sum to zero if the scale choice was a good one.
            //Numerical errors may break this, but this should account for that:
            v += (this->ni[i]*this->one_loop_factor)*Mi2v[i]*Mi2v[i]
                        *(logMi2[i]- this->ci[i]);
        }
    }

    if(include_grav)
    {
        value_type H2 = this->R/this->_12p0;
        value_type H4 = H2*H2;
        value_type V0 = this->eval.V0(g);
        value_type kappa =this->eval.kappa(g);
        value_type a1 = this->eval.alpha1(g);
        value_type a2 = this->eval.alpha2(g);
        value_type a3 = this->eval.alpha3(g);
        value_type alpha = this->_144p0*a1 + this->_36p0*a2 + this->_24p0*a3;
        v += V0;
        if(this->include_R2)
        {
            v += this->_12p0*kappa*H2 + alpha*H4;
        }
        /*
        std::cout << "\nH2 = " << H2;
        std::cout << "\nV0 = " << V0;
        std::cout << "\nkappa = " << kappa;
        std::cout << "\nalpha1 = " << a1;
        std::cout << "\nalpha2 = " << a2;
        std::cout << "\nalpha3 = " << a3;
        std::cout << "\nalpha = " << alpha;
        std::cout << "\nV = " << v;*/
        if(this->include_loops)
        {
            for(int i = 0;i < this->nLoops;i++)
            {
                //std::cout << "\nnpi = " << this->npi[i];
                //log mass term for loop corrections:
                //value_type Mi2 = kiarray[i]*phi2 - kpiarray[i]
                //                    + thetaiarray[i]*this->R;
                //loop corrections:
                //These should sum to zero if the scale choice was a good one.
                //Numerical errors may break this, but this should account for that:
                v += this->npi[i]*this->one_loop_factor*H4*logMi2[i];
            }
        }
    }
    //std::cout << "\nV = " << v;
    return v;
}
//1st derivative:
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::Vp
(value_type y,std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& logMi2,std::vector<value_type>& Mi2v,
 std::vector<value_type>& logMi2p,std::vector<value_type>& Mi2vp,
 bool include_xi,bool include_grav)
{
    value_type lambda = this->eval.l(g);
    value_type lambdap = this->eval.l(gp);
    //value_type xi = g[6];
    //value_type xip = gp[6];
    value_type m2 = this->eval.m2(g);
    value_type m2p = this->eval.m2(gp);
    value_type Z = this->eval.Z(g);
    value_type Zp = this->eval.Z(gp);
    value_type phi = Z*y;
    value_type phip = (Zp*y + Z);

    //Non-minimal coupling. If we are evaluating the Einstein-frame
    //potential, then we want to exclude the xi*R*phi^2/2 term, which is
    //conveniently done by setting xi = 0. Note that this doesn't apply to
    //the log terms! These are still xi dependent, but receive this via
    //the Mi2 functions, which are computed with the correct (non-zero in
    //general) xi.
    value_type xi = include_xi ? this->eval.xi(g) : this->_0p0;
    value_type xip = include_xi ? this->eval.xi(gp) : this->_0p0;

    /*std::cout << "\nxi = " << xi << "\nlambda = " << lambda
              << "\nxip = " << xip << "\nlambdap = " << lambdap
              << "\nm2 = " << m2 << "\nm2p = " << m2p;*/

    value_type vp = (((this->_1p0)/(this->_4p0))*(phi*((lambdap*phi*phi*phi)
                    + ((this->_4p0*(lambda*(phi*phi*phip)))
                    + ((this->_m4p0*(phip*(m2 - (this->R*xi))))
                       + (this->_m2p0*(phi*(m2p - (this->R*xip)))))))));

    //std::cout << "\nvp = " << vp;


    if(this->include_loops)
    {
        for(int i = 0;i < this->nLoops;i++)
        {
            value_type logFactor = logMi2[i];
            value_type logFactorp = logMi2p[i];
            value_type Mi2 = Mi2v[i];
            value_type Mi2p = Mi2vp[i];
            /*std::cout << "\nlogFactor = " << logFactor << " logFactorp = "
                      << logFactorp << " Mi2 = " << Mi2 << " Mi2p = " << Mi2p;*/
            vp += (((this->_1p0)/(this->_64p0))*(Mi2*(((logFactorp*Mi2)
                    + (this->_2p0*((logFactor - this->ci[i])*Mi2p)))
                    *(this->ni[i]*((this->_1p0)/(this->PI*this->PI))))));

        }
    }
    //std::cout << "\nvp = " << vp;

    if(include_grav)
    {
        value_type H2 = this->R/this->_12p0;
        value_type H4 = H2*H2;
        value_type V0p = this->eval.V0(gp);
        value_type kappap = this->eval.kappa(gp);
        value_type a1p = this->eval.alpha1(gp);
        value_type a2p = this->eval.alpha2(gp);
        value_type a3p = this->eval.alpha3(gp);
        value_type alphap = this->_144p0*a1p + this->_36p0*a2p
                            + this->_24p0*a3p;
        vp += V0p;
        if(this->include_R2)
        {
            vp += this->_12p0*kappap*H2 + alphap*H4;
        }
        //std::cout << "\nvp = " << vp;
        if(this->include_loops)
        {
            for(int i = 0;i < this->nLoops;i++)
            {
                //log mass term for loop corrections:
                //value_type Mi2 = kiarray[i]*phi2 - kpiarray[i]
                //                    + thetaiarray[i]*this->R;
                //loop corrections:
                //These should sum to zero if the scale choice was a good one.
                //Numerical errors may break this, but this should account for that:
                vp += this->npi[i]*this->one_loop_factor*H4*logMi2p[i];
            }
        }
    }
    //std::cout << "\nvp = " << vp;
    return vp;
}
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::Vp
(value_type y,std::vector<value_type>& g,std::vector<value_type>& gp,
 bool include_xi,bool include_grav)
{
    //Compute log terms:
    std::vector<value_type> Mi2;
    std::vector<value_type> logMi2;
    std::vector<value_type> Mi2p;
    std::vector<value_type> logMi2p;
    this->computeMi2p(y,g,gp,Mi2,Mi2p,logMi2,logMi2p);

    //Compute potential:
    return this->Vp(y,g,gp,logMi2,Mi2,logMi2p,Mi2p,include_xi,include_grav);
}
//2nd derivative:
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::Vpp
(value_type y,std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,
 std::vector<value_type>& logMi2,std::vector<value_type>& Mi2v,
 std::vector<value_type>& logMi2p,std::vector<value_type>& Mi2vp,
 std::vector<value_type>& logMi2pp,std::vector<value_type>& Mi2vpp,
 bool include_xi,bool include_grav)
{
    value_type lambda = this->eval.l(g);
    value_type lambdap = this->eval.l(gp);
    value_type lambdapp = this->eval.l(gpp);
    //value_type xi = g[6];
    //value_type xip = gp[6];
    //value_type xipp = gpp[6];
    value_type m2 = this->eval.m2(g);
    value_type m2p = this->eval.m2(gp);
    value_type m2pp = this->eval.m2(gpp);
    value_type Z = this->eval.Z(g);
    value_type Zp = this->eval.Z(gp);
    value_type Zpp = this->eval.Z(gpp);
    value_type phi = Z*y;
    value_type phip = (Zp*y + Z);
    value_type phipp = Zpp*y + this->_2p0*Zp;

    //Non-minimal coupling. If we are evaluating the Einstein-frame
    //potential, then we want to exclude the xi*R*phi^2/2 term, which is
    //conveniently done by setting xi = 0. Note that this doesn't apply to
    //the log terms! These are still xi dependent, but receive this via
    //the Mi2 functions, which are computed with the correct (non-zero in
    //general) xi.
    value_type xi = include_xi ? this->eval.xi(g) : this->_0p0;
    value_type xip = include_xi ? this->eval.xi(gp) : this->_0p0;
    value_type xipp = include_xi ? this->eval.xi(gpp) : this->_0p0;

    value_type vpp = ((((this->_1p0)/(this->_4p0))*(lambdapp*phi*phi*phi*phi))
                     + ((phi*phi*phi*((this->_2p0*(lambdap*phip))
                    + (lambda*phipp))) + ((phip*phip*((this->R*xi) - m2))
                    + ((phi*((this->_m2p0*(m2p*phip))
                    + ((phipp*((this->R*xi) - m2))
                    + (this->_2p0*(phip*(this->R*xip))))))
                    + (((this->_1p0)/(this->_2p0))*(phi*phi*(((this->_6p0
                    *(lambda*phip*phip)) + (this->R*xipp)) - m2pp)))))));


    if(this->include_loops)
    {
        for(int i = 0;i < this->nLoops;i++)
        {
            value_type logFactor = logMi2[i];
            value_type logFactorp = logMi2p[i];
            value_type logFactorpp = logMi2pp[i];
            value_type Mi2 = Mi2v[i];
            value_type Mi2p = Mi2vp[i];
            value_type Mi2pp = Mi2vpp[i];
            vpp += (((this->_1p0)/(this->_64p0))*(((logFactorpp*Mi2*Mi2)
                    + ((this->_2p0*((logFactor - this->ci[i])*Mi2p*Mi2p))
                    + (this->_2p0*(Mi2*((this->_2p0*(logFactorp*Mi2p))
                    + ((logFactor - this->ci[i])*Mi2pp))))))
                    *(this->ni[i]*((this->_1p0)/(this->PI*this->PI)))));
        }
    }

    if(include_grav)
    {
        value_type H2 = this->R/this->_12p0;
        value_type H4 = H2*H2;
        value_type V0pp = this->eval.V0(gpp);
        value_type kappapp = this->eval.kappa(gpp);
        value_type a1pp = this->eval.alpha1(gpp);
        value_type a2pp = this->eval.alpha2(gpp);
        value_type a3pp = this->eval.alpha3(gpp);
        value_type alphapp = this->_144p0*a1pp + this->_36p0*a2pp
                            + this->_24p0*a3pp;
        vpp += V0pp;
        if(this->include_R2)
        {
            vpp += this->_12p0*kappapp*H2 + alphapp*H4;
        }
        if(this->include_loops)
        {
            for(int i = 0;i < this->nLoops;i++)
            {
                //log mass term for loop corrections:
                //value_type Mi2 = kiarray[i]*phi2 - kpiarray[i]
                //                    + thetaiarray[i]*this->R;
                //loop corrections:
                //These should sum to zero if the scale choice was a good one.
                //Numerical errors may break this, but this should account for that:
                vpp += this->npi[i]*this->one_loop_factor*H4*logMi2pp[i];
            }
        }
    }
    return vpp;
}
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::Vpp
(value_type y,std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,bool include_xi,bool include_grav)
{
    //Compute log terms:
    std::vector<value_type> Mi2;
    std::vector<value_type> logMi2;
    std::vector<value_type> Mi2p;
    std::vector<value_type> logMi2p;
    std::vector<value_type> Mi2pp;
    std::vector<value_type> logMi2pp;
    this->computeMi2pp(y,g,gp,gpp,Mi2,Mi2p,Mi2pp,logMi2,logMi2p,logMi2pp);

    //Compute potential:
    return this->Vpp(y,g,gp,gpp,logMi2,Mi2,logMi2p,Mi2p,logMi2pp,Mi2pp,
                     include_xi,include_grav);
}
//3rd derivative:
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::Vppp
(value_type y,std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,std::vector<value_type>& gppp,
 std::vector<value_type>& logMi2,std::vector<value_type>& Mi2v,
 std::vector<value_type>& logMi2p,std::vector<value_type>& Mi2vp,
 std::vector<value_type>& logMi2pp,std::vector<value_type>& Mi2vpp,
 std::vector<value_type>& logMi2ppp,std::vector<value_type>& Mi2vppp,
 bool include_xi,bool include_grav)
{
    value_type lambda = this->eval.l(g);
    value_type lambdap = this->eval.l(gp);
    value_type lambdapp = this->eval.l(gpp);
    value_type lambdappp = this->eval.l(gppp);
    //value_type xi = g[6];
    //value_type xip = gp[6];
    //value_type xipp = gpp[6];
    //value_type xippp = gppp[6];
    value_type m2 = this->eval.m2(g);
    value_type m2p = this->eval.m2(gp);
    value_type m2pp = this->eval.m2(gpp);
    value_type m2ppp = this->eval.m2(gppp);
    value_type Z = this->eval.Z(g);
    value_type Zp = this->eval.Z(gp);
    value_type Zpp = this->eval.Z(gpp);
    value_type Zppp = this->eval.Z(gppp);
    value_type phi = Z*y;
    value_type phip = Zp*y + Z;
    value_type phipp = Zpp*y + this->_2p0*Zp;
    value_type phippp = Zppp*y + this->_3p0*Zpp;

    //Non-minimal coupling. If we are evaluating the Einstein-frame
    //potential, then we want to exclude the xi*R*phi^2/2 term, which is
    //conveniently done by setting xi = 0. Note that this doesn't apply to
    //the log terms! These are still xi dependent, but receive this via
    //the Mi2 functions, which are computed with the correct (non-zero in
    //general) xi.
    value_type xi = include_xi ? this->eval.xi(g) : this->_0p0;
    value_type xip = include_xi ? this->eval.xi(gp) : this->_0p0;
    value_type xipp = include_xi ? this->eval.xi(gpp) : this->_0p0;
    value_type xippp = include_xi ? this->eval.xi(gppp) : this->_0p0;

    value_type vppp = ((((this->_1p0)/(this->_4p0))*(lambdappp*phi*phi*phi*phi))
                     + ((this->_m3p0*(m2p*(phip*phip + (phi*phipp))))
                     + ((phi*phi*phi*((this->_3p0*(lambdapp*phip))
                     + ((this->_3p0*(lambdap*phipp)) + (lambda*phippp))))
                     + ((this->_3p0*(phip*((phipp*((this->R*xi) - m2))
                     + (phip*(this->R*xip))))) + ((phi*(((this->_6p0
                     *(lambda*phip*phip*phip)) + ((phippp*(this->R*xi))
                     + ((this->_3p0*(phipp*(this->R*xip)))
                     + (this->_m3p0*(phip*(m2pp - (this->R*xipp)))))))
                     - (m2*phippp))) + (((this->_1p0)/(this->_2p0))*(phi*phi
                     *(((this->_18p0*(lambdap*phip*phip))
                     + ((this->_18p0*(lambda*(phip*phipp))) + (this->R*xippp)))
                     - m2ppp))))))));




    if(this->include_loops)
    {
        for(int i = 0;i < this->nLoops;i++)
        {
            value_type logFactor = logMi2[i];
            value_type logFactorp = logMi2p[i];
            value_type logFactorpp = logMi2pp[i];
            value_type logFactorppp = logMi2ppp[i];
            value_type Mi2 = Mi2v[i];
            value_type Mi2p = Mi2vp[i];
            value_type Mi2pp = Mi2vpp[i];
            value_type Mi2ppp = Mi2vppp[i];
            vppp += (((this->_1p0)/(this->_64p0))*(((logFactorppp*Mi2*Mi2)
                    + ((this->_6p0*((logFactor - this->ci[i])*(Mi2p*Mi2pp)))
                    + ((this->_6p0*(logFactorp*(Mi2p*Mi2p + (Mi2*Mi2pp))))
                    + (this->_2p0*(Mi2*((this->_3p0*(logFactorpp*Mi2p))
                    + ((logFactor - this->ci[i])*Mi2ppp)))))))
                    *(this->ni[i]*((this->_1p0)/(this->PI*this->PI)))));
        }
    }

    if(include_grav)
    {
        value_type H2 = this->R/this->_12p0;
        value_type H4 = H2*H2;
        value_type V0ppp = this->eval.V0(gppp);
        value_type kappappp = this->eval.kappa(gppp);
        value_type a1ppp = this->eval.alpha1(gppp);
        value_type a2ppp = this->eval.alpha2(gppp);
        value_type a3ppp = this->eval.alpha3(gppp);
        value_type alphappp = this->_144p0*a1ppp + this->_36p0*a2ppp
                            + this->_24p0*a3ppp;
        vppp += V0ppp;
        if(this->include_R2)
        {
            vppp += this->_12p0*kappappp*H2 + alphappp*H4;
        }
        if(this->include_loops)
        {
            for(int i = 0;i < this->nLoops;i++)
            {
                //log mass term for loop corrections:
                //value_type Mi2 = kiarray[i]*phi2 - kpiarray[i]
                //                    + thetaiarray[i]*this->R;
                //loop corrections:
                //These should sum to zero if the scale choice was a good one.
                //Numerical errors may break this, but this should account for that:
                vppp += this->npi[i]*this->one_loop_factor*H4*logMi2ppp[i];
            }
        }
    }
    return vppp;
}
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::Vppp
(value_type y,std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,std::vector<value_type>& gppp,bool include_xi,bool include_grav)
{
    //Compute log terms:
    std::vector<value_type> Mi2;
    std::vector<value_type> logMi2;
    std::vector<value_type> Mi2p;
    std::vector<value_type> logMi2p;
    std::vector<value_type> Mi2pp;
    std::vector<value_type> logMi2pp;
    std::vector<value_type> Mi2ppp;
    std::vector<value_type> logMi2ppp;
    this->computeMi2ppp(y,g,gp,gpp,gppp,Mi2,Mi2p,Mi2pp,Mi2ppp,logMi2,logMi2p,
                        logMi2pp,logMi2ppp);

    //Compute potential:
    return this->Vppp(y,g,gp,gpp,gppp,logMi2,Mi2,logMi2p,Mi2p,
                      logMi2pp,Mi2pp,logMi2ppp,Mi2ppp,include_xi,include_grav);
}
//4th derivative:
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::Vpppp
(value_type y,std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,std::vector<value_type>& gppp,
 std::vector<value_type>& gpppp,
 std::vector<value_type>& logMi2,std::vector<value_type>& Mi2v,
 std::vector<value_type>& logMi2p,std::vector<value_type>& Mi2vp,
 std::vector<value_type>& logMi2pp,std::vector<value_type>& Mi2vpp,
 std::vector<value_type>& logMi2ppp,std::vector<value_type>& Mi2vppp,
 std::vector<value_type>& logMi2pppp,std::vector<value_type>& Mi2vpppp,
 bool include_xi,bool include_grav)
{
    value_type lambda = this->eval.l(g);
    value_type lambdap = this->eval.l(gp);
    value_type lambdapp = this->eval.l(gpp);
    value_type lambdappp = this->eval.l(gppp);
    value_type lambdapppp = this->eval.l(gpppp);
    //value_type xi = g[6];
    //value_type xip = gp[6];
    //value_type xipp = gpp[6];
    //value_type xippp = gppp[6];
    //value_type xipppp = gpppp[6];
    value_type m2 = this->eval.m2(g);
    value_type m2p = this->eval.m2(gp);
    value_type m2pp = this->eval.m2(gpp);
    value_type m2ppp = this->eval.m2(gppp);
    value_type m2pppp = this->eval.m2(gpppp);
    value_type Z = this->eval.Z(g);
    value_type Zp = this->eval.Z(gp);
    value_type Zpp = this->eval.Z(gpp);
    value_type Zppp = this->eval.Z(gppp);
    value_type Zpppp = this->eval.Z(gpppp);
    value_type phi = Z*y;
    value_type phip = Zp*y + Z;
    value_type phipp = Zpp*y + this->_2p0*Zp;
    value_type phippp = Zppp*y + this->_3p0*Zpp;
    value_type phipppp = Zpppp*y + this->_4p0*Zppp;

    //Non-minimal coupling. If we are evaluating the Einstein-frame
    //potential, then we want to exclude the xi*R*phi^2/2 term, which is
    //conveniently done by setting xi = 0. Note that this doesn't apply to
    //the log terms! These are still xi dependent, but receive this via
    //the Mi2 functions, which are computed with the correct (non-zero in
    //general) xi.
    value_type xi = include_xi ? this->eval.xi(g) : this->_0p0;
    value_type xip = include_xi ? this->eval.xi(gp) : this->_0p0;
    value_type xipp = include_xi ? this->eval.xi(gpp) : this->_0p0;
    value_type xippp = include_xi ? this->eval.xi(gppp) : this->_0p0;
    value_type xipppp = include_xi ? this->eval.xi(gpppp) : this->_0p0;


    value_type vpppp = ((((this->_1p0)/(this->_4p0))
                        *(lambdapppp*phi*phi*phi*phi))
                        + ((this->_m6p0*(m2pp*phip*phip))
                        + ((this->_6p0*(lambda*phip*phip*phip*phip))
                        + ((this->_m12p0*(m2p*(phip*phipp)))
                        + ((this->_m3p0*(m2*phipp*phipp))
                        + ((this->_m4p0*(m2*(phip*phippp)))
                        + ((phi*phi*phi*((this->_4p0*(lambdappp*phip))
                        + ((this->_6p0*(lambdapp*phipp))
                        + ((this->_4p0*(lambdap*phippp))
                        + (lambda*phipppp)))))
                        + ((this->_3p0*(phipp*phipp*(this->R*xi)))
                        + ((this->_4p0*(phip*(phippp*(this->R*xi))))
                        + ((this->_12p0*(phip*(phipp*(this->R*xip))))
                        + ((this->_6p0*(phip*phip*(this->R*xipp)))
                        + ((phi*(((this->_m4p0*(m2ppp*phip))
                        + ((this->_24p0*(lambdap*phip*phip*phip))
                        + ((this->_m6p0*(m2pp*phipp))
                        + ((this->_36p0*(lambda*(phip*phip*phipp)))
                        + ((this->_m4p0*(m2p*phippp)) + ((phipppp*(this->R*xi))
                        + ((this->_4p0*(phippp*(this->R*xip)))
                        + ((this->_6p0*(phipp*(this->R*xipp)))
                        + (this->_4p0*(phip*(this->R*xippp)))))))))))
                        - (m2*phipppp))) + (((this->_1p0)/(this->_2p0))
                        *(phi*phi*(((this->_36p0*(lambdapp*phip*phip))
                        + ((this->_18p0*(lambda*phipp*phipp))
                        + ((this->_24p0*(phip*((this->_3p0*(lambdap*phipp))
                        + (lambda*phippp)))) + (this->R*xipppp))))
                        - m2pppp)))))))))))))));


    if(this->include_loops)
    {
        for(int i = 0;i < this->nLoops;i++)
        {
            value_type logFactor = logMi2[i];
            value_type logFactorp = logMi2p[i];
            value_type logFactorpp = logMi2pp[i];
            value_type logFactorppp = logMi2ppp[i];
            value_type logFactorpppp = logMi2pppp[i];
            value_type Mi2 = Mi2v[i];
            value_type Mi2p = Mi2vp[i];
            value_type Mi2pp = Mi2vpp[i];
            value_type Mi2ppp = Mi2vppp[i];
            value_type Mi2pppp = Mi2vpppp[i];
            vpppp += (((this->_1p0)/(this->_64p0))*(((logFactorpppp*Mi2*Mi2)
                    + ((this->_12p0*(logFactorpp*Mi2p*Mi2p))
                    + ((this->_6p0*((logFactor - this->ci[i])*Mi2pp*Mi2pp))
                    + ((this->_8p0*(Mi2p*((logFactorppp*Mi2)
                    + ((this->_3p0*(logFactorp*Mi2pp))
                    + ((logFactor - this->ci[i])*Mi2ppp)))))
                    + (this->_2p0*(Mi2*((this->_6p0*(logFactorpp*Mi2pp))
                    + ((this->_4p0*(logFactorp*Mi2ppp))
                    + ((logFactor - this->ci[i])*Mi2pppp)))))))))
                    *(this->ni[i]*((this->_1p0)/(this->PI*this->PI)))));
        }
    }

    if(include_grav)
    {
        value_type H2 = this->R/this->_12p0;
        value_type H4 = H2*H2;
        value_type V0pppp = this->eval.V0(gpppp);
        value_type kappapppp = this->eval.kappa(gpppp);
        value_type a1pppp = this->eval.alpha1(gpppp);
        value_type a2pppp = this->eval.alpha2(gpppp);
        value_type a3pppp = this->eval.alpha3(gpppp);
        value_type alphapppp = this->_144p0*a1pppp + this->_36p0*a2pppp
                            + this->_24p0*a3pppp;
        vpppp += V0pppp;
        if(this->include_R2)
        {
            vpppp += this->_12p0*kappapppp*H2 + alphapppp*H4;
        }
        if(this->include_loops)
        {
            for(int i = 0;i < this->nLoops;i++)
            {
                //log mass term for loop corrections:
                //value_type Mi2 = kiarray[i]*phi2 - kpiarray[i]
                //                    + thetaiarray[i]*this->R;
                //loop corrections:
                //These should sum to zero if the scale choice was a good one.
                //Numerical errors may break this, but this should account for that:
                vpppp += this->npi[i]*this->one_loop_factor*H4*logMi2pppp[i];
            }
        }
    }
    return vpppp;
}
template<class value_type>
value_type SMpotential_dS_1loop<value_type>::Vpppp
(value_type y,std::vector<value_type>& g,std::vector<value_type>& gp,
 std::vector<value_type>& gpp,std::vector<value_type>& gppp,
 std::vector<value_type>& gpppp,bool include_xi,bool include_grav)
{
    //Compute log terms:
    std::vector<value_type> Mi2;
    std::vector<value_type> logMi2;
    std::vector<value_type> Mi2p;
    std::vector<value_type> logMi2p;
    std::vector<value_type> Mi2pp;
    std::vector<value_type> logMi2pp;
    std::vector<value_type> Mi2ppp;
    std::vector<value_type> logMi2ppp;
    std::vector<value_type> Mi2pppp;
    std::vector<value_type> logMi2pppp;
    this->computeMi2pppp(y,g,gp,gpp,gppp,gpppp,Mi2,Mi2p,Mi2pp,Mi2ppp,Mi2pppp,
                        logMi2,logMi2p,logMi2pp,logMi2ppp,logMi2pppp);

    //Compute potential:
    return this->Vpppp(y,g,gp,gpp,gppp,gpppp,logMi2,Mi2,logMi2p,Mi2p,
                       logMi2pp,Mi2pp,logMi2ppp,Mi2ppp,logMi2pppp,Mi2pppp,
                       include_xi,include_grav);
}
template<class value_type>
DLL_EXPORT value_type SMpotential_dS_1loop<value_type>::operator()
(value_type y)
{
    if(splineCouplings)
    {
        value_type x = y*this->yscale;
        //std::cout << "\nx = " << x;
        //Choose scale to evaluate couplings at:
        int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
        //Use spline interpolation to figure out the couplings at this scale:
        std::vector<value_type> gdata (this->nCouplings,_0p0);

        //std::cout << "\ng = (";
        for(int i = 0;i < this->nCouplings;i++)
        {
            gdata[i] = this->Spline(y,this->ypdata[i],this->ydata[i],n);
            /*
            if(i != 0)
            {
                std::cout << ",";
            }*/
            //std::cout << gdata[i];
        }
        //std::cout << ")";
        //Return value of the potential at this scale:
        return this->V(x,gdata,this->include_xi,this->include_grav);
    }
    else
    {
        value_type x = y*this->yscale;
        //std::cout << "\nx = " << x;
        //Choose scale to evaluate couplings at:
        int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
        return this->Spline(y,this->vpdata.data(),this->vdata.data(),n);
    }
}
template<class value_type>
DLL_EXPORT void SMpotential_dS_1loop<value_type>::operator()
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i< numel;i++)
    {
        arrayOut[i] = this->operator()(y[i]);
    }
}
template<class value_type>
DLL_EXPORT value_type SMpotential_dS_1loop<value_type>::d
(value_type y)
{
    if(splineCouplings)
    {
        value_type x = y*this->yscale;
        //Choose scale to evaluate couplings at:
        int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
        //Use spline interpolation to figure out the couplings at this scale:
        std::vector<value_type> gdata (this->nCouplings,_0p0);//couplings
        std::vector<value_type> gpdata_phicl (this->nCouplings,_0p0);//couplings
        std::vector<value_type> gpdata_tcl (this->nCouplings,_0p0);//derivatives of couplings wrt
                //tcl = log(phicl^2/mt^2)
        value_type phicl = x;
        for(int i = 0;i < this->nCouplings;i++)
        {
            gdata[i] = this->Spline(y,this->ypdata[i],this->ydata[i],n);
            //gpdata_tcl[i] = this->Spline(y,this->yppdata[i],this->ypdata[i],n);
            gpdata_tcl[i] = this->SplineDerivative(y,this->ypdata[i],this->ydata[i],n);
            //We actually need derivatives of the couplings wrt phicl, not tcl:
            gpdata_phicl[i] = this->_2p0*gpdata_tcl[i]/phicl;
        }
        //Return value of the potential at this scale:
        value_type vp = this->Vp(x,gdata,gpdata_phicl,this->include_xi,this->include_grav);
        return vp;
    }
    else
    {
        value_type x = y*this->yscale;
        //std::cout << "\nx = " << x;
        //Choose scale to evaluate couplings at:
        int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
        return this->_2p0*this->SplineDerivative(y,this->vpdata.data(),
                                                 this->vdata.data(),n)/x;
    }
}
template<class value_type>
DLL_EXPORT void SMpotential_dS_1loop<value_type>::d
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i< numel;i++)
    {
        arrayOut[i] = this->d(y[i]);
    }
}
template<class value_type>
DLL_EXPORT value_type SMpotential_dS_1loop<value_type>::d2
(value_type y)
{
    value_type x = y*this->yscale;
    //Choose scale to evaluate couplings at:
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    //Use spline interpolation to figure out the couplings at this scale:
    std::vector<value_type> gdata (this->nCouplings,this->_0p0);//couplings
    std::vector<value_type> gpdata_phicl (this->nCouplings,this->_0p0);//couplings
    std::vector<value_type> gpdata_tcl (this->nCouplings,this->_0p0);//derivatives of couplings wrt
            //tcl = log(phicl^2/mt^2)
    std::vector<value_type> gppdata_phicl (this->nCouplings,this->_0p0);
    std::vector<value_type> gppdata_tcl (this->nCouplings,this->_0p0);
    value_type phicl = x;
    for(int i = 0;i < this->nCouplings;i++)
    {
        gdata[i] = this->Spline(y,this->ypdata[i],this->ydata[i],n);
        gpdata_tcl[i] = this->Spline(y,this->yppdata[i],this->ypdata[i],n);
        gppdata_tcl[i] = this->Spline(y,this->ypppdata[i],this->yppdata[i],n);
        //We actually need derivatives of the couplings wrt phicl, not tcl:
        gpdata_phicl[i] = this->_2p0*gpdata_tcl[i]/phicl;
        gppdata_phicl[i] = (-this->_2p0*gpdata_tcl[i]
                            + this->_4p0*gppdata_tcl[i])/(phicl*phicl);
    }
    //Return value of the potential at this scale:
    value_type vpp = this->Vpp(x,gdata,gpdata_phicl,gppdata_phicl,this->include_xi,
                     this->include_grav);
    return vpp;
}
template<class value_type>
DLL_EXPORT void SMpotential_dS_1loop<value_type>::d2
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i< numel;i++)
    {
        arrayOut[i] = this->d2(y[i]);
    }
}
template<class value_type>
DLL_EXPORT value_type SMpotential_dS_1loop<value_type>::d3
(value_type y)
{
    value_type x = y*this->yscale;
    //Choose scale to evaluate couplings at:
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    //Use spline interpolation to figure out the couplings at this scale:
    std::vector<value_type> gdata (this->nCouplings,this->_0p0);//couplings
    std::vector<value_type> gpdata_phicl (this->nCouplings,this->_0p0);//couplings
    std::vector<value_type> gpdata_tcl (this->nCouplings,this->_0p0);//derivatives of couplings wrt
            //tcl = log(phicl^2/mt^2)
    std::vector<value_type> gppdata_phicl (this->nCouplings,this->_0p0);
    std::vector<value_type> gppdata_tcl (this->nCouplings,this->_0p0);
    std::vector<value_type> gpppdata_phicl (this->nCouplings,this->_0p0);
    std::vector<value_type> gpppdata_tcl (this->nCouplings,this->_0p0);
    value_type phicl = x;
    for(int i = 0;i < this->nCouplings;i++)
    {
        gdata[i] = this->Spline(y,this->ypdata[i],this->ydata[i],n);
        gpdata_tcl[i] = this->Spline(y,this->yppdata[i],this->ypdata[i],n);
        gppdata_tcl[i] = this->Spline(y,this->ypppdata[i],this->yppdata[i],n);
        gpppdata_tcl[i] = this->Spline(y,this->yppppdata[i],
                                       this->ypppdata[i],n);
        //We actually need derivatives of the couplings wrt phicl, not tcl:
        gpdata_phicl[i] = this->_2p0*gpdata_tcl[i]/phicl;
        gppdata_phicl[i] = (-this->_2p0*gpdata_tcl[i]
                            + this->_4p0*gppdata_tcl[i])/(phicl*phicl);
        gpppdata_phicl[i] = (this->_4p0*gpdata_tcl[i]
                             -this->_12p0*gppdata_tcl[i]
                             + this->_8p0*gpppdata_tcl[i])/(phicl*phicl*phicl);
    }
    //Return value of the potential at this scale:
    value_type vppp =  this->Vppp(x,gdata,gpdata_phicl,gppdata_phicl,gpppdata_phicl,
                      this->include_xi,this->include_grav);
    return vppp;
}


template<class value_type>
DLL_EXPORT void SMpotential_dS_1loop<value_type>::d3
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i< numel;i++)
    {
        arrayOut[i] = this->d3(y[i]);
    }
}
template<class value_type>
DLL_EXPORT value_type SMpotential_dS_1loop<value_type>::d4
(value_type y)
{
    value_type x = y*this->yscale;
    //Choose scale to evaluate couplings at:
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    //Use spline interpolation to figure out the couplings at this scale:
    std::vector<value_type> gdata (this->nCouplings,this->_0p0);//couplings
    std::vector<value_type> gpdata_phicl (this->nCouplings,this->_0p0);//couplings
    std::vector<value_type> gpdata_tcl (this->nCouplings,this->_0p0);//derivatives of couplings wrt
            //tcl = log(phicl^2/mt^2)
    std::vector<value_type> gppdata_phicl (this->nCouplings,this->_0p0);
    std::vector<value_type> gppdata_tcl (this->nCouplings,this->_0p0);
    std::vector<value_type> gpppdata_phicl (this->nCouplings,this->_0p0);
    std::vector<value_type> gpppdata_tcl (this->nCouplings,this->_0p0);
    std::vector<value_type> gppppdata_phicl (this->nCouplings,this->_0p0);
    std::vector<value_type> gppppdata_tcl (this->nCouplings,this->_0p0);
    //Intermediates:
    value_type phicl = x;
    value_type phicl2 = phicl*phicl;
    value_type phicl3 = phicl2*phicl;
    value_type phicl4 = phicl3*phicl;
    //Loop to populate arrays of couplings and their derivatives
    //at the chosen scale:
    for(int i = 0;i < this->nCouplings;i++)
    {
        //derivatives wrt tcl (obtain by spline):
        gdata[i] = this->Spline(y,this->ypdata[i],this->ydata[i],n);
        gpdata_tcl[i] = this->Spline(y,this->yppdata[i],this->ypdata[i],n);
        gppdata_tcl[i] = this->Spline(y,this->ypppdata[i],this->yppdata[i],n);
        gpppdata_tcl[i] = this->Spline(y,this->yppppdata[i],
                                       this->ypppdata[i],n);
        gppppdata_tcl[i] = this->Spline(y,this->ypppppdata[i],
                                        this->yppppdata[i],n);
        //We actually need derivatives of the couplings wrt phicl, not tcl:
        gpdata_phicl[i] = this->_2p0*gpdata_tcl[i]/phicl;
        gppdata_phicl[i] = (-this->_2p0*gpdata_tcl[i]
                            + this->_4p0*gppdata_tcl[i])/(phicl2);
        gpppdata_phicl[i] = (this->_4p0*gpdata_tcl[i]
                             -this->_12p0*gppdata_tcl[i]
                             + this->_8p0*gpppdata_tcl[i])/(phicl3);
        gppppdata_phicl[i] = (-this->_12p0*gpdata_tcl[i]
                              + this->_44p0*gppdata_tcl[i]
                              -this->_48p0*gpppdata_tcl[i]
                              +this->_16p0*gppppdata_tcl[i])/(phicl4);

    }
    //Return value of the potential at this scale:
    value_type vpppp = this->Vpppp(x,gdata,gpdata_phicl,gppdata_phicl,gpppdata_phicl,
                       gpppdata_phicl,this->include_xi,this->include_grav);
    return vpppp;
}

template<class value_type>
DLL_EXPORT void SMpotential_dS_1loop<value_type>::d4
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i< numel;i++)
    {
        arrayOut[i] = this->d4(y[i]);
    }
}

template<class value_type>
void betaSMdS_1loop<value_type>::beta_basic(const std::vector<value_type>& g,
                                            std::vector<value_type>& beta_g,
                                            const value_type t)
{
    if(beta_g.size() != g.size())
    {
        throw("Beta function and coupling arrays must be the same size.");
    }


    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));


    //Beta functions:
    beta_g[0] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_27p0)/(this->_400p0))*g12*g12) + ((((this->_9p0)/(this->_40p0))*(g12*g22)) + ((((this->_9p0)/(this->_16p0))*g22*g22) + (l*((((this->_m9p0)/(this->_10p0))*g12) + ((((this->_m9p0)/(this->_2p0))*g22) + ((this->_12p0*l) + (this->_2p0*Y2S)))))))) - HS)));
    beta_g[1] = (((this->_1p0)/(this->_16p0))*(g12*g12*((((this->_1p0)/(this->_10p0)) + ((((this->_3p0)/(this->_5p0)) + (((this->_11p0)/(this->_45p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))));
    beta_g[2] = (((this->_1p0)/(this->_16p0))*(g22*g22*((((this->_m43p0)/(this->_6p0)) + ((((this->_1p0)/(this->_3p0)) + (((this->_1p0)/(this->_3p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))));
    beta_g[3] = (((this->_1p0)/(this->_16p0))*(g32*g32*(((this->_1p0)/(PI*PI))*((((this->_m11p0)/(this->_3p0))*Ca) + (((this->_8p0)/(this->_3p0))*(Ng*Tf))))));
    beta_g[4] = (((this->_1p0)/(this->_16p0))*(m2*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S))))));
    beta_g[5] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ye2*((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ye2))))));
    beta_g[6] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ymu2*((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ymu2))))));
    beta_g[7] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ytau2*((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ytau2))))));
    beta_g[8] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yu2*((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yu2 - yd2)))))))));;
    beta_g[9] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yd2 - yu2)))))))));
    beta_g[10] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2*((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yc2 - ys2)))))))));
    beta_g[11] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ys2*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(ys2 - yc2)))))))));
    beta_g[12] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yt2*((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yt2 - yb2)))))))));
    beta_g[13] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yb2 - yt2)))))))));
}

template<class value_type>
void betaSMdS_1loop<value_type>::betap_basic(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
                 const std::vector<value_type>& gp,const value_type t)
{
    if(beta_g.size() != g.size())
    {
        throw("Beta function and coupling arrays must be the same size.");
    }

    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));
    value_type HSp = ((this->_2p0*(ye2*ye2p)) + ((this->_2p0*(ymu2*ymu2p)) + ((this->_3p0*((this->_2p0*(yb2*yb2p)) + ((this->_2p0*(yd2*yd2p)) + (this->_2p0*(ys2*ys2p))))) + ((this->_2p0*(ytau2*ytau2p)) + (this->_3p0*((this->_2p0*(yc2*yc2p)) + ((this->_2p0*(yt2*yt2p)) + (this->_2p0*(yu2*yu2p)))))))));


    /*
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    this->beta_basic(g,gp,t);*/


    beta_g[0] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_27p0)/(this->_200p0))*(g12*g12p)) + ((((this->_9p0)/(this->_40p0))*(g12p*g22)) + ((((this->_9p0)/(this->_40p0))*(g12*g22p)) + ((((this->_9p0)/(this->_8p0))*(g22*g22p)) + ((lp*((((this->_m9p0)/(this->_10p0))*g12) + ((((this->_m9p0)/(this->_2p0))*g22) + ((this->_12p0*l) + (this->_2p0*Y2S))))) + (l*((((this->_m9p0)/(this->_10p0))*g12p) + ((((this->_m9p0)/(this->_2p0))*g22p) + ((this->_12p0*lp) + (this->_2p0*Y2Sp)))))))))) - HSp)));
    beta_g[1] = (((this->_1p0)/(this->_8p0))*(g12*(g12p*((((this->_1p0)/(this->_10p0)) + ((((this->_3p0)/(this->_5p0)) + (((this->_11p0)/(this->_45p0))*Nc))*Ng))*((this->_1p0)/(PI*PI))))));
    beta_g[2] = (((this->_1p0)/(this->_8p0))*(g22*(g22p*((((this->_m43p0)/(this->_6p0)) + ((((this->_1p0)/(this->_3p0)) + (((this->_1p0)/(this->_3p0))*Nc))*Ng))*((this->_1p0)/(PI*PI))))));
    beta_g[3] = (((this->_1p0)/(this->_8p0))*(g32*(g32p*(((this->_1p0)/(PI*PI))*((((this->_m11p0)/(this->_3p0))*Ca) + (((this->_8p0)/(this->_3p0))*(Ng*Tf)))))));
    beta_g[4] = ((((this->_1p0)/(this->_16p0))*(m2p*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S)))))) + (((this->_1p0)/(this->_16p0))*(m2*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_6p0*lp) + Y2Sp)))))));
    beta_g[5] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ye2)))*ye2p))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ye2*((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ye2p)))))));
    beta_g[6] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ymu2)))*ymu2p))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ymu2*((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ymu2p)))))));
    beta_g[7] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ytau2)))*ytau2p))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ytau2*((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ytau2p)))))));
    beta_g[8] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yu2 - yd2))))))*yu2p))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yu2*((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yu2p - yd2p))))))))));
    beta_g[9] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2p*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yd2 - yu2))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yd2p - yu2p))))))))));
    beta_g[10] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2p*((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yc2 - ys2))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2*((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yc2p - ys2p))))))))));
    beta_g[11] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(ys2 - yc2))))))*ys2p))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ys2*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(ys2p - yc2p))))))))));
    beta_g[12] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yt2 - yb2))))))*yt2p))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yt2*((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yt2p - yb2p))))))))));
    beta_g[13] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2p*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yb2 - yt2))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yb2p - yt2p))))))))));


}
template<class value_type>
void betaSMdS_1loop<value_type>::betap_basic(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
                 const value_type t)
{
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    this->beta_basic(g,gp,t);

    this->betap_basic(g,beta_g,gp,t);
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapp_basic(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
            const std::vector<value_type>& gp,
            const std::vector<value_type>& gpp,
            const value_type t)
{

    if(beta_g.size() != g.size())
    {
        throw("Beta function and coupling arrays must be the same size.");
    }
    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    //Coupling derivatives:
    value_type lpp = gpp[0];
    value_type g12pp = gpp[1];
    value_type g22pp = gpp[2];
    value_type g32pp = gpp[3];
    value_type m2pp = gpp[4];
    value_type ye2pp = gpp[5];
    value_type ymu2pp = gpp[6];
    value_type ytau2pp = gpp[7];
    value_type yu2pp = gpp[8];
    value_type yd2pp = gpp[9];
    value_type yc2pp = gpp[10];
    value_type ys2pp = gpp[11];
    value_type yt2pp = gpp[12];
    value_type yb2pp = gpp[13];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type Y2Spp = (ye2pp + (ymu2pp + ((this->_3p0*(yb2pp + (yd2pp + ys2pp))) + (ytau2pp + (this->_3p0*(yc2pp + (yt2pp + yu2pp)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));
    value_type HSp = ((this->_2p0*(ye2*ye2p)) + ((this->_2p0*(ymu2*ymu2p)) + ((this->_3p0*((this->_2p0*(yb2*yb2p)) + ((this->_2p0*(yd2*yd2p)) + (this->_2p0*(ys2*ys2p))))) + ((this->_2p0*(ytau2*ytau2p)) + (this->_3p0*((this->_2p0*(yc2*yc2p)) + ((this->_2p0*(yt2*yt2p)) + (this->_2p0*(yu2*yu2p)))))))));
    value_type HSpp = ((this->_2p0*ye2p*ye2p) + ((this->_2p0*(ye2*ye2pp)) + ((this->_2p0*ymu2p*ymu2p) + ((this->_2p0*(ymu2*ymu2pp)) + ((this->_3p0*((this->_2p0*yb2p*yb2p) + ((this->_2p0*(yb2*yb2pp)) + ((this->_2p0*yd2p*yd2p) + ((this->_2p0*(yd2*yd2pp)) + ((this->_2p0*ys2p*ys2p) + (this->_2p0*(ys2*ys2pp)))))))) + ((this->_2p0*ytau2p*ytau2p) + ((this->_2p0*(ytau2*ytau2pp)) + (this->_3p0*((this->_2p0*yc2p*yc2p) + ((this->_2p0*(yc2*yc2pp)) + ((this->_2p0*yt2p*yt2p) + ((this->_2p0*(yt2*yt2pp)) + ((this->_2p0*yu2p*yu2p) + (this->_2p0*(yu2*yu2pp)))))))))))))));

    beta_g[0] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_27p0)/(this->_200p0))*g12p*g12p) + ((((this->_27p0)/(this->_200p0))*(g12*g12pp)) + ((((this->_9p0)/(this->_40p0))*(g12pp*g22)) + ((((this->_9p0)/(this->_20p0))*(g12p*g22p)) + ((((this->_9p0)/(this->_8p0))*g22p*g22p) + ((((this->_9p0)/(this->_40p0))*(g12*g22pp)) + ((((this->_9p0)/(this->_8p0))*(g22*g22pp)) + ((lpp*((((this->_m9p0)/(this->_10p0))*g12) + ((((this->_m9p0)/(this->_2p0))*g22) + ((this->_12p0*l) + (this->_2p0*Y2S))))) + ((this->_2p0*(lp*((((this->_m9p0)/(this->_10p0))*g12p) + ((((this->_m9p0)/(this->_2p0))*g22p) + ((this->_12p0*lp) + (this->_2p0*Y2Sp)))))) + (l*((((this->_m9p0)/(this->_10p0))*g12pp) + ((((this->_m9p0)/(this->_2p0))*g22pp) + ((this->_12p0*lpp) + (this->_2p0*Y2Spp)))))))))))))) - HSpp)));
    beta_g[1] = ((((this->_1p0)/(this->_8p0))*(g12p*g12p*((((this->_1p0)/(this->_10p0)) + ((((this->_3p0)/(this->_5p0)) + (((this->_11p0)/(this->_45p0))*Nc))*Ng))*((this->_1p0)/(PI*PI))))) + (((this->_1p0)/(this->_8p0))*(g12*(g12pp*((((this->_1p0)/(this->_10p0)) + ((((this->_3p0)/(this->_5p0)) + (((this->_11p0)/(this->_45p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))))));
    beta_g[2] = ((((this->_1p0)/(this->_8p0))*(g22p*g22p*((((this->_m43p0)/(this->_6p0)) + ((((this->_1p0)/(this->_3p0)) + (((this->_1p0)/(this->_3p0))*Nc))*Ng))*((this->_1p0)/(PI*PI))))) + (((this->_1p0)/(this->_8p0))*(g22*(g22pp*((((this->_m43p0)/(this->_6p0)) + ((((this->_1p0)/(this->_3p0)) + (((this->_1p0)/(this->_3p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))))));
    beta_g[3] = ((((this->_1p0)/(this->_8p0))*(g32p*g32p*(((this->_1p0)/(PI*PI))*((((this->_m11p0)/(this->_3p0))*Ca) + (((this->_8p0)/(this->_3p0))*(Ng*Tf)))))) + (((this->_1p0)/(this->_8p0))*(g32*(g32pp*(((this->_1p0)/(PI*PI))*((((this->_m11p0)/(this->_3p0))*Ca) + (((this->_8p0)/(this->_3p0))*(Ng*Tf))))))));
    beta_g[4] = ((((this->_1p0)/(this->_16p0))*(m2pp*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S)))))) + ((((this->_1p0)/(this->_8p0))*(m2p*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_6p0*lp) + Y2Sp)))))) + (((this->_1p0)/(this->_16p0))*(m2*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_6p0*lpp) + Y2Spp))))))));
    beta_g[5] = ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(ye2p*((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ye2p)))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ye2)))*ye2pp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ye2*((((this->_m9p0)/(this->_4p0))*(g12pp + g22pp)) + (Y2Spp + (((this->_3p0)/(this->_2p0))*ye2pp))))))));
    beta_g[6] = ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(ymu2p*((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ymu2p)))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ymu2)))*ymu2pp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ymu2*((((this->_m9p0)/(this->_4p0))*(g12pp + g22pp)) + (Y2Spp + (((this->_3p0)/(this->_2p0))*ymu2pp))))))));
    beta_g[7] = ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(ytau2p*((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ytau2p)))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ytau2)))*ytau2pp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ytau2*((((this->_m9p0)/(this->_4p0))*(g12pp + g22pp)) + (Y2Spp + (((this->_3p0)/(this->_2p0))*ytau2pp))))))));
    beta_g[8] = ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yu2p*((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yu2p - yd2p))))))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yu2 - yd2))))))*yu2pp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yu2*((((this->_m17p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yu2pp - yd2pp)))))))))));
    beta_g[9] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2pp*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yd2 - yu2))))))))) + ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yd2p*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yd2p - yu2p))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2*((((this->_m1p0)/(this->_4p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yd2pp - yu2pp)))))))))));
    beta_g[10] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2pp*((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yc2 - ys2))))))))) + ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yc2p*((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yc2p - ys2p))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2*((((this->_m17p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yc2pp - ys2pp)))))))))));
    beta_g[11] = ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(ys2p*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(ys2p - yc2p))))))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(ys2 - yc2))))))*ys2pp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ys2*((((this->_m1p0)/(this->_4p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(ys2pp - yc2pp)))))))))));
    beta_g[12] = ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yt2p*((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yt2p - yb2p))))))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yt2 - yb2))))))*yt2pp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yt2*((((this->_m17p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yt2pp - yb2pp)))))))))));
    beta_g[13] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2pp*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yb2 - yt2))))))))) + ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yb2p*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yb2p - yt2p))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2*((((this->_m1p0)/(this->_4p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yb2pp - yt2pp)))))))))));


  }
template<class value_type>
void betaSMdS_1loop<value_type>::betapp_basic(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
            const value_type t)
{
    //Generate data about lower derivatives:
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    std::vector<value_type> gpp (int(g.size()),this->_0p0);
    this->beta_basic(g,gp,t);
    this->betap_basic(g,gpp,gp,t);

    //Generate second derivatives of the beta function:
    this->betapp_basic(g,beta_g,gp,gpp,t);
}
template<class value_type>
void betaSMdS_1loop<value_type>::betappp_basic(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
            const std::vector<value_type>& gp,
            const std::vector<value_type>& gpp,
            const std::vector<value_type>& gppp,
            const value_type t)
{

    if(beta_g.size() != g.size())
    {
        throw("Beta function and coupling arrays must be the same size.");
    }
    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    //Coupling derivatives:
    value_type lpp = gpp[0];
    value_type g12pp = gpp[1];
    value_type g22pp = gpp[2];
    value_type g32pp = gpp[3];
    value_type m2pp = gpp[4];
    value_type ye2pp = gpp[5];
    value_type ymu2pp = gpp[6];
    value_type ytau2pp = gpp[7];
    value_type yu2pp = gpp[8];
    value_type yd2pp = gpp[9];
    value_type yc2pp = gpp[10];
    value_type ys2pp = gpp[11];
    value_type yt2pp = gpp[12];
    value_type yb2pp = gpp[13];

    //Coupling derivatives:
    value_type lppp = gppp[0];
    value_type g12ppp = gppp[1];
    value_type g22ppp = gppp[2];
    value_type g32ppp = gppp[3];
    value_type m2ppp = gppp[4];
    value_type ye2ppp = gppp[5];
    value_type ymu2ppp = gppp[6];
    value_type ytau2ppp = gppp[7];
    value_type yu2ppp = gppp[8];
    value_type yd2ppp = gppp[9];
    value_type yc2ppp = gppp[10];
    value_type ys2ppp = gppp[11];
    value_type yt2ppp = gppp[12];
    value_type yb2ppp = gppp[13];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type Y2Spp = (ye2pp + (ymu2pp + ((this->_3p0*(yb2pp + (yd2pp + ys2pp))) + (ytau2pp + (this->_3p0*(yc2pp + (yt2pp + yu2pp)))))));
    value_type Y2Sppp = (ye2ppp + (ymu2ppp + ((this->_3p0*(yb2ppp + (yd2ppp + ys2ppp))) + (ytau2ppp + (this->_3p0*(yc2ppp + (yt2ppp + yu2ppp)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));
    value_type HSp = ((this->_2p0*(ye2*ye2p)) + ((this->_2p0*(ymu2*ymu2p)) + ((this->_3p0*((this->_2p0*(yb2*yb2p)) + ((this->_2p0*(yd2*yd2p)) + (this->_2p0*(ys2*ys2p))))) + ((this->_2p0*(ytau2*ytau2p)) + (this->_3p0*((this->_2p0*(yc2*yc2p)) + ((this->_2p0*(yt2*yt2p)) + (this->_2p0*(yu2*yu2p)))))))));
    value_type HSpp = ((this->_2p0*ye2p*ye2p) + ((this->_2p0*(ye2*ye2pp)) + ((this->_2p0*ymu2p*ymu2p) + ((this->_2p0*(ymu2*ymu2pp)) + ((this->_3p0*((this->_2p0*yb2p*yb2p) + ((this->_2p0*(yb2*yb2pp)) + ((this->_2p0*yd2p*yd2p) + ((this->_2p0*(yd2*yd2pp)) + ((this->_2p0*ys2p*ys2p) + (this->_2p0*(ys2*ys2pp)))))))) + ((this->_2p0*ytau2p*ytau2p) + ((this->_2p0*(ytau2*ytau2pp)) + (this->_3p0*((this->_2p0*yc2p*yc2p) + ((this->_2p0*(yc2*yc2pp)) + ((this->_2p0*yt2p*yt2p) + ((this->_2p0*(yt2*yt2pp)) + ((this->_2p0*yu2p*yu2p) + (this->_2p0*(yu2*yu2pp)))))))))))))));
    value_type HSppp = ((this->_6p0*(ye2p*ye2pp)) + ((this->_2p0*(ye2*ye2ppp)) + ((this->_6p0*(ymu2p*ymu2pp)) + ((this->_2p0*(ymu2*ymu2ppp)) + ((this->_3p0*((this->_6p0*(yb2p*yb2pp)) + ((this->_2p0*(yb2*yb2ppp)) + ((this->_6p0*(yd2p*yd2pp)) + ((this->_2p0*(yd2*yd2ppp)) + ((this->_6p0*(ys2p*ys2pp)) + (this->_2p0*(ys2*ys2ppp)))))))) + ((this->_6p0*(ytau2p*ytau2pp)) + ((this->_2p0*(ytau2*ytau2ppp)) + (this->_3p0*((this->_6p0*(yc2p*yc2pp)) + ((this->_2p0*(yc2*yc2ppp)) + ((this->_6p0*(yt2p*yt2pp)) + ((this->_2p0*(yt2*yt2ppp)) + ((this->_6p0*(yu2p*yu2pp)) + (this->_2p0*(yu2*yu2ppp)))))))))))))));

    beta_g[0] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_81p0)/(this->_200p0))*(g12p*g12pp)) + ((((this->_27p0)/(this->_200p0))*(g12*g12ppp)) + ((((this->_9p0)/(this->_40p0))*(g12ppp*g22)) + ((((this->_27p0)/(this->_40p0))*(g12pp*g22p)) + ((((this->_27p0)/(this->_40p0))*(g12p*g22pp)) + ((((this->_27p0)/(this->_8p0))*(g22p*g22pp)) + ((((this->_9p0)/(this->_40p0))*(g12*g22ppp)) + ((((this->_9p0)/(this->_8p0))*(g22*g22ppp)) + ((lppp*((((this->_m9p0)/(this->_10p0))*g12) + ((((this->_m9p0)/(this->_2p0))*g22) + ((this->_12p0*l) + (this->_2p0*Y2S))))) + ((this->_3p0*(lpp*((((this->_m9p0)/(this->_10p0))*g12p) + ((((this->_m9p0)/(this->_2p0))*g22p) + ((this->_12p0*lp) + (this->_2p0*Y2Sp)))))) + ((this->_3p0*(lp*((((this->_m9p0)/(this->_10p0))*g12pp) + ((((this->_m9p0)/(this->_2p0))*g22pp) + ((this->_12p0*lpp) + (this->_2p0*Y2Spp)))))) + (l*((((this->_m9p0)/(this->_10p0))*g12ppp) + ((((this->_m9p0)/(this->_2p0))*g22ppp) + ((this->_12p0*lppp) + (this->_2p0*Y2Sppp)))))))))))))))) - HSppp)));
    beta_g[1] = ((((this->_3p0)/(this->_8p0))*(g12p*(g12pp*((((this->_1p0)/(this->_10p0)) + ((((this->_3p0)/(this->_5p0)) + (((this->_11p0)/(this->_45p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))))) + (((this->_1p0)/(this->_8p0))*(g12*(g12ppp*((((this->_1p0)/(this->_10p0)) + ((((this->_3p0)/(this->_5p0)) + (((this->_11p0)/(this->_45p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))))));
    beta_g[2] = ((((this->_3p0)/(this->_8p0))*(g22p*(g22pp*((((this->_m43p0)/(this->_6p0)) + ((((this->_1p0)/(this->_3p0)) + (((this->_1p0)/(this->_3p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))))) + (((this->_1p0)/(this->_8p0))*(g22*(g22ppp*((((this->_m43p0)/(this->_6p0)) + ((((this->_1p0)/(this->_3p0)) + (((this->_1p0)/(this->_3p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))))));
    beta_g[3] = ((((this->_3p0)/(this->_8p0))*(g32p*(g32pp*(((this->_1p0)/(PI*PI))*((((this->_m11p0)/(this->_3p0))*Ca) + (((this->_8p0)/(this->_3p0))*(Ng*Tf))))))) + (((this->_1p0)/(this->_8p0))*(g32*(g32ppp*(((this->_1p0)/(PI*PI))*((((this->_m11p0)/(this->_3p0))*Ca) + (((this->_8p0)/(this->_3p0))*(Ng*Tf))))))));
    beta_g[4] = ((((this->_1p0)/(this->_16p0))*(m2ppp*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S)))))) + ((((this->_3p0)/(this->_16p0))*(m2pp*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_6p0*lp) + Y2Sp)))))) + ((((this->_3p0)/(this->_16p0))*(m2p*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_6p0*lpp) + Y2Spp)))))) + (((this->_1p0)/(this->_16p0))*(m2*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_6p0*lppp) + Y2Sppp)))))))));
    beta_g[5] = ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ye2p)))*ye2pp))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ye2p*((((this->_m9p0)/(this->_4p0))*(g12pp + g22pp)) + (Y2Spp + (((this->_3p0)/(this->_2p0))*ye2pp)))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ye2)))*ye2ppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ye2*((((this->_m9p0)/(this->_4p0))*(g12ppp + g22ppp)) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*ye2ppp)))))))));
    beta_g[6] = ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ymu2p)))*ymu2pp))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ymu2p*((((this->_m9p0)/(this->_4p0))*(g12pp + g22pp)) + (Y2Spp + (((this->_3p0)/(this->_2p0))*ymu2pp)))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ymu2)))*ymu2ppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ymu2*((((this->_m9p0)/(this->_4p0))*(g12ppp + g22ppp)) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*ymu2ppp)))))))));
    beta_g[7] = ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ytau2p)))*ytau2pp))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ytau2p*((((this->_m9p0)/(this->_4p0))*(g12pp + g22pp)) + (Y2Spp + (((this->_3p0)/(this->_2p0))*ytau2pp)))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ytau2)))*ytau2ppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ytau2*((((this->_m9p0)/(this->_4p0))*(g12ppp + g22ppp)) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*ytau2ppp)))))))));
    beta_g[8] = ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yu2p - yd2p))))))*yu2pp))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yu2p*((((this->_m17p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yu2pp - yd2pp))))))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yu2 - yd2))))))*yu2ppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yu2*((((this->_m17p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yu2ppp - yd2ppp))))))))))));
    beta_g[9] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2ppp*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yd2 - yu2))))))))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2pp*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yd2p - yu2p))))))))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2p*((((this->_m1p0)/(this->_4p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yd2pp - yu2pp))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2*((((this->_m1p0)/(this->_4p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yd2ppp - yu2ppp))))))))))));
    beta_g[10] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2ppp*((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yc2 - ys2))))))))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2pp*((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yc2p - ys2p))))))))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2p*((((this->_m17p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yc2pp - ys2pp))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2*((((this->_m17p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yc2ppp - ys2ppp))))))))))));
    beta_g[11] = ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(ys2p - yc2p))))))*ys2pp))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ys2p*((((this->_m1p0)/(this->_4p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(ys2pp - yc2pp))))))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(ys2 - yc2))))))*ys2ppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ys2*((((this->_m1p0)/(this->_4p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(ys2ppp - yc2ppp))))))))))));
    beta_g[12] = ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yt2p - yb2p))))))*yt2pp))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yt2p*((((this->_m17p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yt2pp - yb2pp))))))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yt2 - yb2))))))*yt2ppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yt2*((((this->_m17p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yt2ppp - yb2ppp))))))))))));
    beta_g[13] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2ppp*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yb2 - yt2))))))))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2pp*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yb2p - yt2p))))))))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2p*((((this->_m1p0)/(this->_4p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yb2pp - yt2pp))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2*((((this->_m1p0)/(this->_4p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yb2ppp - yt2ppp))))))))))));

}
template<class value_type>
void betaSMdS_1loop<value_type>::betappp_basic(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
            const value_type t)
{
    //Generate data about lower derivatives:
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    std::vector<value_type> gpp (int(g.size()),this->_0p0);
    std::vector<value_type> gppp (int(g.size()),this->_0p0);
    this->beta_basic(g,gp,t);
    this->betap_basic(g,gpp,gp,t);
    this->betapp_basic(g,gppp,gp,gpp,t);

    //Generate second derivatives of the beta function:
    this->betappp_basic(g,beta_g,gp,gpp,gppp,t);
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapppp_basic(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
            const std::vector<value_type>& gp,
            const std::vector<value_type>& gpp,
            const std::vector<value_type>& gppp,
            const std::vector<value_type>& gpppp,
            const value_type t)
{

    if(beta_g.size() != g.size())
    {
        throw("Beta function and coupling arrays must be the same size.");
    }
    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    //Coupling derivatives:
    value_type lpp = gpp[0];
    value_type g12pp = gpp[1];
    value_type g22pp = gpp[2];
    value_type g32pp = gpp[3];
    value_type m2pp = gpp[4];
    value_type ye2pp = gpp[5];
    value_type ymu2pp = gpp[6];
    value_type ytau2pp = gpp[7];
    value_type yu2pp = gpp[8];
    value_type yd2pp = gpp[9];
    value_type yc2pp = gpp[10];
    value_type ys2pp = gpp[11];
    value_type yt2pp = gpp[12];
    value_type yb2pp = gpp[13];

    //Coupling derivatives:
    value_type lppp = gppp[0];
    value_type g12ppp = gppp[1];
    value_type g22ppp = gppp[2];
    value_type g32ppp = gppp[3];
    value_type m2ppp = gppp[4];
    value_type ye2ppp = gppp[5];
    value_type ymu2ppp = gppp[6];
    value_type ytau2ppp = gppp[7];
    value_type yu2ppp = gppp[8];
    value_type yd2ppp = gppp[9];
    value_type yc2ppp = gppp[10];
    value_type ys2ppp = gppp[11];
    value_type yt2ppp = gppp[12];
    value_type yb2ppp = gppp[13];

    //Coupling derivatives:
    value_type lpppp = gpppp[0];
    value_type g12pppp = gpppp[1];
    value_type g22pppp = gpppp[2];
    value_type g32pppp = gpppp[3];
    value_type m2pppp = gpppp[4];
    value_type ye2pppp = gpppp[5];
    value_type ymu2pppp = gpppp[6];
    value_type ytau2pppp = gpppp[7];
    value_type yu2pppp = gpppp[8];
    value_type yd2pppp = gpppp[9];
    value_type yc2pppp = gpppp[10];
    value_type ys2pppp = gpppp[11];
    value_type yt2pppp = gpppp[12];
    value_type yb2pppp = gpppp[13];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type Y2Spp = (ye2pp + (ymu2pp + ((this->_3p0*(yb2pp + (yd2pp + ys2pp))) + (ytau2pp + (this->_3p0*(yc2pp + (yt2pp + yu2pp)))))));
    value_type Y2Sppp = (ye2ppp + (ymu2ppp + ((this->_3p0*(yb2ppp + (yd2ppp + ys2ppp))) + (ytau2ppp + (this->_3p0*(yc2ppp + (yt2ppp + yu2ppp)))))));
    value_type Y2Spppp = (ye2pppp + (ymu2pppp + ((this->_3p0*(yb2pppp + (yd2pppp + ys2pppp))) + (ytau2pppp + (this->_3p0*(yc2pppp + (yt2pppp + yu2pppp)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));
    value_type HSp = ((this->_2p0*(ye2*ye2p)) + ((this->_2p0*(ymu2*ymu2p)) + ((this->_3p0*((this->_2p0*(yb2*yb2p)) + ((this->_2p0*(yd2*yd2p)) + (this->_2p0*(ys2*ys2p))))) + ((this->_2p0*(ytau2*ytau2p)) + (this->_3p0*((this->_2p0*(yc2*yc2p)) + ((this->_2p0*(yt2*yt2p)) + (this->_2p0*(yu2*yu2p)))))))));
    value_type HSpp = ((this->_2p0*ye2p*ye2p) + ((this->_2p0*(ye2*ye2pp)) + ((this->_2p0*ymu2p*ymu2p) + ((this->_2p0*(ymu2*ymu2pp)) + ((this->_3p0*((this->_2p0*yb2p*yb2p) + ((this->_2p0*(yb2*yb2pp)) + ((this->_2p0*yd2p*yd2p) + ((this->_2p0*(yd2*yd2pp)) + ((this->_2p0*ys2p*ys2p) + (this->_2p0*(ys2*ys2pp)))))))) + ((this->_2p0*ytau2p*ytau2p) + ((this->_2p0*(ytau2*ytau2pp)) + (this->_3p0*((this->_2p0*yc2p*yc2p) + ((this->_2p0*(yc2*yc2pp)) + ((this->_2p0*yt2p*yt2p) + ((this->_2p0*(yt2*yt2pp)) + ((this->_2p0*yu2p*yu2p) + (this->_2p0*(yu2*yu2pp)))))))))))))));
    value_type HSppp = ((this->_6p0*(ye2p*ye2pp)) + ((this->_2p0*(ye2*ye2ppp)) + ((this->_6p0*(ymu2p*ymu2pp)) + ((this->_2p0*(ymu2*ymu2ppp)) + ((this->_3p0*((this->_6p0*(yb2p*yb2pp)) + ((this->_2p0*(yb2*yb2ppp)) + ((this->_6p0*(yd2p*yd2pp)) + ((this->_2p0*(yd2*yd2ppp)) + ((this->_6p0*(ys2p*ys2pp)) + (this->_2p0*(ys2*ys2ppp)))))))) + ((this->_6p0*(ytau2p*ytau2pp)) + ((this->_2p0*(ytau2*ytau2ppp)) + (this->_3p0*((this->_6p0*(yc2p*yc2pp)) + ((this->_2p0*(yc2*yc2ppp)) + ((this->_6p0*(yt2p*yt2pp)) + ((this->_2p0*(yt2*yt2ppp)) + ((this->_6p0*(yu2p*yu2pp)) + (this->_2p0*(yu2*yu2ppp)))))))))))))));
    value_type HSpppp = ((this->_6p0*ye2pp*ye2pp) + ((this->_8p0*(ye2p*ye2ppp)) + ((this->_2p0*(ye2*ye2pppp)) + ((this->_6p0*ymu2pp*ymu2pp) + ((this->_8p0*(ymu2p*ymu2ppp)) + ((this->_2p0*(ymu2*ymu2pppp)) + ((this->_3p0*((this->_6p0*yb2pp*yb2pp) + ((this->_8p0*(yb2p*yb2ppp)) + ((this->_2p0*(yb2*yb2pppp)) + ((this->_6p0*yd2pp*yd2pp) + ((this->_8p0*(yd2p*yd2ppp)) + ((this->_2p0*(yd2*yd2pppp)) + ((this->_6p0*ys2pp*ys2pp) + ((this->_8p0*(ys2p*ys2ppp)) + (this->_2p0*(ys2*ys2pppp))))))))))) + ((this->_6p0*ytau2pp*ytau2pp) + ((this->_8p0*(ytau2p*ytau2ppp)) + ((this->_2p0*(ytau2*ytau2pppp)) + (this->_3p0*((this->_6p0*yc2pp*yc2pp) + ((this->_8p0*(yc2p*yc2ppp)) + ((this->_2p0*(yc2*yc2pppp)) + ((this->_6p0*yt2pp*yt2pp) + ((this->_8p0*(yt2p*yt2ppp)) + ((this->_2p0*(yt2*yt2pppp)) + ((this->_6p0*yu2pp*yu2pp) + ((this->_8p0*(yu2p*yu2ppp)) + (this->_2p0*(yu2*yu2pppp)))))))))))))))))))));

    beta_g[0] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_81p0)/(this->_200p0))*g12pp*g12pp) + ((((this->_27p0)/(this->_50p0))*(g12p*g12ppp)) + ((((this->_27p0)/(this->_200p0))*(g12*g12pppp)) + ((((this->_9p0)/(this->_40p0))*(g12pppp*g22)) + ((((this->_9p0)/(this->_10p0))*(g12ppp*g22p)) + ((((this->_27p0)/(this->_20p0))*(g12pp*g22pp)) + ((((this->_27p0)/(this->_8p0))*g22pp*g22pp) + ((((this->_9p0)/(this->_10p0))*(g12p*g22ppp)) + ((((this->_9p0)/(this->_2p0))*(g22p*g22ppp)) + ((((this->_9p0)/(this->_40p0))*(g12*g22pppp)) + ((((this->_9p0)/(this->_8p0))*(g22*g22pppp)) + ((lpppp*((((this->_m9p0)/(this->_10p0))*g12) + ((((this->_m9p0)/(this->_2p0))*g22) + ((this->_12p0*l) + (this->_2p0*Y2S))))) + ((this->_4p0*(lppp*((((this->_m9p0)/(this->_10p0))*g12p) + ((((this->_m9p0)/(this->_2p0))*g22p) + ((this->_12p0*lp) + (this->_2p0*Y2Sp)))))) + ((this->_6p0*(lpp*((((this->_m9p0)/(this->_10p0))*g12pp) + ((((this->_m9p0)/(this->_2p0))*g22pp) + ((this->_12p0*lpp) + (this->_2p0*Y2Spp)))))) + ((this->_4p0*(lp*((((this->_m9p0)/(this->_10p0))*g12ppp) + ((((this->_m9p0)/(this->_2p0))*g22ppp) + ((this->_12p0*lppp) + (this->_2p0*Y2Sppp)))))) + (l*((((this->_m9p0)/(this->_10p0))*g12pppp) + ((((this->_m9p0)/(this->_2p0))*g22pppp) + ((this->_12p0*lpppp) + (this->_2p0*Y2Spppp)))))))))))))))))))) - HSpppp)));
    beta_g[1] = ((((this->_3p0)/(this->_8p0))*(g12pp*g12pp*((((this->_1p0)/(this->_10p0)) + ((((this->_3p0)/(this->_5p0)) + (((this->_11p0)/(this->_45p0))*Nc))*Ng))*((this->_1p0)/(PI*PI))))) + ((((this->_1p0)/(this->_2p0))*(g12p*(g12ppp*((((this->_1p0)/(this->_10p0)) + ((((this->_3p0)/(this->_5p0)) + (((this->_11p0)/(this->_45p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))))) + (((this->_1p0)/(this->_8p0))*(g12*(g12pppp*((((this->_1p0)/(this->_10p0)) + ((((this->_3p0)/(this->_5p0)) + (((this->_11p0)/(this->_45p0))*Nc))*Ng))*((this->_1p0)/(PI*PI))))))));
    beta_g[2] = ((((this->_3p0)/(this->_8p0))*(g22pp*g22pp*((((this->_m43p0)/(this->_6p0)) + ((((this->_1p0)/(this->_3p0)) + (((this->_1p0)/(this->_3p0))*Nc))*Ng))*((this->_1p0)/(PI*PI))))) + ((((this->_1p0)/(this->_2p0))*(g22p*(g22ppp*((((this->_m43p0)/(this->_6p0)) + ((((this->_1p0)/(this->_3p0)) + (((this->_1p0)/(this->_3p0))*Nc))*Ng))*((this->_1p0)/(PI*PI)))))) + (((this->_1p0)/(this->_8p0))*(g22*(g22pppp*((((this->_m43p0)/(this->_6p0)) + ((((this->_1p0)/(this->_3p0)) + (((this->_1p0)/(this->_3p0))*Nc))*Ng))*((this->_1p0)/(PI*PI))))))));
    beta_g[3] = ((((this->_3p0)/(this->_8p0))*(g32pp*g32pp*(((this->_1p0)/(PI*PI))*((((this->_m11p0)/(this->_3p0))*Ca) + (((this->_8p0)/(this->_3p0))*(Ng*Tf)))))) + ((((this->_1p0)/(this->_2p0))*(g32p*(g32ppp*(((this->_1p0)/(PI*PI))*((((this->_m11p0)/(this->_3p0))*Ca) + (((this->_8p0)/(this->_3p0))*(Ng*Tf))))))) + (((this->_1p0)/(this->_8p0))*(g32*(g32pppp*(((this->_1p0)/(PI*PI))*((((this->_m11p0)/(this->_3p0))*Ca) + (((this->_8p0)/(this->_3p0))*(Ng*Tf)))))))));
    beta_g[4] = ((((this->_1p0)/(this->_16p0))*(m2pppp*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S)))))) + ((((this->_1p0)/(this->_4p0))*(m2ppp*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_6p0*lp) + Y2Sp)))))) + ((((this->_3p0)/(this->_8p0))*(m2pp*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_6p0*lpp) + Y2Spp)))))) + ((((this->_1p0)/(this->_4p0))*(m2p*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_6p0*lppp) + Y2Sppp)))))) + (((this->_1p0)/(this->_16p0))*(m2*(((this->_1p0)/(PI*PI))*((((this->_m9p0)/(this->_20p0))*g12pppp) + ((((this->_m9p0)/(this->_4p0))*g22pppp) + ((this->_6p0*lpppp) + Y2Spppp))))))))));
    beta_g[5] = ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(ye2pp*((((this->_m9p0)/(this->_4p0))*(g12pp + g22pp)) + (Y2Spp + (((this->_3p0)/(this->_2p0))*ye2pp)))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ye2p)))*ye2ppp))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(ye2p*((((this->_m9p0)/(this->_4p0))*(g12ppp + g22ppp)) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*ye2ppp)))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ye2)))*ye2pppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ye2*((((this->_m9p0)/(this->_4p0))*(g12pppp + g22pppp)) + (Y2Spppp + (((this->_3p0)/(this->_2p0))*ye2pppp))))))))));
    beta_g[6] = ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(ymu2pp*((((this->_m9p0)/(this->_4p0))*(g12pp + g22pp)) + (Y2Spp + (((this->_3p0)/(this->_2p0))*ymu2pp)))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ymu2p)))*ymu2ppp))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(ymu2p*((((this->_m9p0)/(this->_4p0))*(g12ppp + g22ppp)) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*ymu2ppp)))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ymu2)))*ymu2pppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ymu2*((((this->_m9p0)/(this->_4p0))*(g12pppp + g22pppp)) + (Y2Spppp + (((this->_3p0)/(this->_2p0))*ymu2pppp))))))))));
    beta_g[7] = ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(ytau2pp*((((this->_m9p0)/(this->_4p0))*(g12pp + g22pp)) + (Y2Spp + (((this->_3p0)/(this->_2p0))*ytau2pp)))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12p + g22p)) + (Y2Sp + (((this->_3p0)/(this->_2p0))*ytau2p)))*ytau2ppp))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(ytau2p*((((this->_m9p0)/(this->_4p0))*(g12ppp + g22ppp)) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*ytau2ppp)))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_4p0))*(g12 + g22)) + (Y2S + (((this->_3p0)/(this->_2p0))*ytau2)))*ytau2pppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ytau2*((((this->_m9p0)/(this->_4p0))*(g12pppp + g22pppp)) + (Y2Spppp + (((this->_3p0)/(this->_2p0))*ytau2pppp))))))))));
    beta_g[8] = ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yu2pp*((((this->_m17p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yu2pp - yd2pp))))))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yu2p - yd2p))))))*yu2ppp))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(yu2p*((((this->_m17p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yu2ppp - yd2ppp))))))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yu2 - yd2))))))*yu2pppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yu2*((((this->_m17p0)/(this->_20p0))*g12pppp) + ((((this->_m9p0)/(this->_4p0))*g22pppp) + ((this->_m8p0*g32pppp) + (Y2Spppp + (((this->_3p0)/(this->_2p0))*(yu2pppp - yd2pppp)))))))))))));
    beta_g[9] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2pppp*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yd2 - yu2))))))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(yd2ppp*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yd2p - yu2p))))))))) + ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yd2pp*((((this->_m1p0)/(this->_4p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yd2pp - yu2pp))))))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(yd2p*((((this->_m1p0)/(this->_4p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yd2ppp - yu2ppp))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yd2*((((this->_m1p0)/(this->_4p0))*g12pppp) + ((((this->_m9p0)/(this->_4p0))*g22pppp) + ((this->_m8p0*g32pppp) + (Y2Spppp + (((this->_3p0)/(this->_2p0))*(yd2pppp - yu2pppp)))))))))))));
    beta_g[10] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2pppp*((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yc2 - ys2))))))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(yc2ppp*((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yc2p - ys2p))))))))) + ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yc2pp*((((this->_m17p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yc2pp - ys2pp))))))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(yc2p*((((this->_m17p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yc2ppp - ys2ppp))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yc2*((((this->_m17p0)/(this->_20p0))*g12pppp) + ((((this->_m9p0)/(this->_4p0))*g22pppp) + ((this->_m8p0*g32pppp) + (Y2Spppp + (((this->_3p0)/(this->_2p0))*(yc2pppp - ys2pppp)))))))))))));
    beta_g[11] = ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(ys2pp*((((this->_m1p0)/(this->_4p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(ys2pp - yc2pp))))))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(ys2p - yc2p))))))*ys2ppp))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(ys2p*((((this->_m1p0)/(this->_4p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(ys2ppp - yc2ppp))))))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(ys2 - yc2))))))*ys2pppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(ys2*((((this->_m1p0)/(this->_4p0))*g12pppp) + ((((this->_m9p0)/(this->_4p0))*g22pppp) + ((this->_m8p0*g32pppp) + (Y2Spppp + (((this->_3p0)/(this->_2p0))*(ys2pppp - yc2pppp)))))))))))));
    beta_g[12] = ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yt2pp*((((this->_m17p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yt2pp - yb2pp))))))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yt2p - yb2p))))))*yt2ppp))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(yt2p*((((this->_m17p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yt2ppp - yb2ppp))))))))) + ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((((this->_m17p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yt2 - yb2))))))*yt2pppp))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yt2*((((this->_m17p0)/(this->_20p0))*g12pppp) + ((((this->_m9p0)/(this->_4p0))*g22pppp) + ((this->_m8p0*g32pppp) + (Y2Spppp + (((this->_3p0)/(this->_2p0))*(yt2pppp - yb2pppp)))))))))))));
    beta_g[13] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2pppp*((((this->_m1p0)/(this->_4p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_m8p0*g32) + (Y2S + (((this->_3p0)/(this->_2p0))*(yb2 - yt2))))))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(yb2ppp*((((this->_m1p0)/(this->_4p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_m8p0*g32p) + (Y2Sp + (((this->_3p0)/(this->_2p0))*(yb2p - yt2p))))))))) + ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(yb2pp*((((this->_m1p0)/(this->_4p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_m8p0*g32pp) + (Y2Spp + (((this->_3p0)/(this->_2p0))*(yb2pp - yt2pp))))))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(yb2p*((((this->_m1p0)/(this->_4p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_m8p0*g32ppp) + (Y2Sppp + (((this->_3p0)/(this->_2p0))*(yb2ppp - yt2ppp))))))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(yb2*((((this->_m1p0)/(this->_4p0))*g12pppp) + ((((this->_m9p0)/(this->_4p0))*g22pppp) + ((this->_m8p0*g32pppp) + (Y2Spppp + (((this->_3p0)/(this->_2p0))*(yb2pppp - yt2pppp)))))))))))));
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapppp_basic(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
            const value_type t)
{
    //Generate data about lower derivatives:
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    std::vector<value_type> gpp (int(g.size()),this->_0p0);
    std::vector<value_type> gppp (int(g.size()),this->_0p0);
    std::vector<value_type> gpppp (int(g.size()),this->_0p0);
    this->beta_basic(g,gp,t);
    this->betap_basic(g,gpp,gp,t);
    this->betapp_basic(g,gppp,gp,gpp,t);
    this->betappp_basic(g,gpppp,gp,gpp,gppp,t);

    //Generate second derivatives of the beta function:
    this->betapppp_basic(g,beta_g,gp,gpp,gppp,gpppp,t);
}
//Beta functions including xi, Z, tcl, and t:
template<class value_type>
void betaSMdS_1loop<value_type>::beta_extended(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const value_type t)
{
    //Fill in the basics:
    this->beta_basic(g,beta_g,t);
    if(int(g.size()) < 17)
    {
        throw("Coefficient array too small.");
    }

    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];
    //value_type t = g[16];//Given by the function argument.

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    //value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));

    //Add the extended functions:
    beta_g[14] = (((this->_m1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))
                    *(Z*((((this->_m9p0)/(this->_40p0))*g12)
                    + ((((this->_m9p0)/(this->_8p0))*g22)
                    + ((((this->_1p0)/(this->_2p0))*Y2S)
                    + ((((this->_m1p0)/(this->_4p0))*(g22*this->zetaw))
                    + (((this->_m1p0)/(this->_4p0))
                    *((g22 + (((this->_3p0)/(this->_5p0))*g12))
                    *this->zetaz)))))))));
    beta_g[15] = this->_1p0;//PLACEHOLDER
    beta_g[16] = this->_1p0;
}
template<class value_type>
void betaSMdS_1loop<value_type>::betap_extended(const std::vector<value_type>& g,
                    std::vector<value_type>& beta_g,
                    const std::vector<value_type>& gp,const value_type t)
{
    if(int(g.size()) < 17)
    {
        throw("Coefficient array too small.");
    }

    //Fill in the basics:
    //std::vector<value_type> gp (int(g.size()),this->_0p0);
    //this->beta_extended(g,gp,t);
    this->betap_basic(g,beta_g,gp,t);


    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    value_type Zp = gp[14];
    value_type tclp = gp[15];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));



    //Add the extended functions:
    beta_g[14] = ((((this->_m1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))
                    *(Z*((((this->_m9p0)/(this->_40p0))*g12p)
                    + ((((this->_m9p0)/(this->_8p0))*g22p)
                    + ((((this->_1p0)/(this->_2p0))*Y2Sp)
                    + ((((this->_m1p0)/(this->_4p0))*(g22p*this->zetaw))
                    + (((this->_m1p0)/(this->_4p0))*((g22p
                    + (((this->_3p0)/(this->_5p0))*g12p))*this->zetaz)))))))))
                    + (((this->_m1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))
                    *(((((this->_m9p0)/(this->_40p0))*g12)
                    + ((((this->_m9p0)/(this->_8p0))*g22)
                    + ((((this->_1p0)/(this->_2p0))*Y2S)
                    + ((((this->_m1p0)/(this->_4p0))*(g22*this->zetaw))
                    + (((this->_m1p0)/(this->_4p0))
                    *((g22 + (((this->_3p0)/(this->_5p0))*g12))
                    *this->zetaz))))))*Zp))));
    beta_g[15] = this->_0p0;//PLACEHOLDER
    beta_g[16] = this->_0p0;
}
template<class value_type>
void betaSMdS_1loop<value_type>::betap_extended(const std::vector<value_type>& g,
                    std::vector<value_type>& beta_g,
                    const value_type t)
{
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    this->beta_extended(g,gp,t);

    //Compute beta functions:
    this->betap_extended(g,beta_g,gp,t);
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapp_extended(const std::vector<value_type>& g,
                     std::vector<value_type>& beta_g,
                     const std::vector<value_type>& gp,
                     const std::vector<value_type>& gpp,const value_type t)
{
    if(int(g.size()) < 17)
    {
        throw("Coefficient array too small.");
    }

    //Fill in the basics:
    this->betapp_basic(g,beta_g,gp,gpp,t);


    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    value_type Zp = gp[14];
    value_type tclp = gp[15];

    //Coupling derivatives:
    value_type lpp = gpp[0];
    value_type g12pp = gpp[1];
    value_type g22pp = gpp[2];
    value_type g32pp = gpp[3];
    value_type m2pp = gpp[4];
    value_type ye2pp = gpp[5];
    value_type ymu2pp = gpp[6];
    value_type ytau2pp = gpp[7];
    value_type yu2pp = gpp[8];
    value_type yd2pp = gpp[9];
    value_type yc2pp = gpp[10];
    value_type ys2pp = gpp[11];
    value_type yt2pp = gpp[12];
    value_type yb2pp = gpp[13];

    value_type Zpp = gpp[14];
    value_type tclpp = gpp[15];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type Y2Spp = (ye2pp + (ymu2pp + ((this->_3p0*(yb2pp + (yd2pp + ys2pp))) + (ytau2pp + (this->_3p0*(yc2pp + (yt2pp + yu2pp)))))));

    //Add the extended functions:
    beta_g[14] = ((((this->_m1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))
                *(Z*((((this->_m9p0)/(this->_40p0))*g12pp)
                + ((((this->_m9p0)/(this->_8p0))*g22pp)
                + ((((this->_1p0)/(this->_2p0))*Y2Spp)
                + ((((this->_m1p0)/(this->_4p0))*(g22pp*this->zetaw))
                + (((this->_m1p0)/(this->_4p0))
                *((g22pp + (((this->_3p0)/(this->_5p0))*g12pp))
                *this->zetaz))))))))) + ((((this->_m1p0)/(this->_8p0))
                *(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_40p0))*g12p)
                + ((((this->_m9p0)/(this->_8p0))*g22p)
                + ((((this->_1p0)/(this->_2p0))*Y2Sp)
                + ((((this->_m1p0)/(this->_4p0))*(g22p*this->zetaw))
                + (((this->_m1p0)/(this->_4p0))
                *((g22p + (((this->_3p0)/(this->_5p0))*g12p))*this->zetaz))))))
                *Zp))) + (((this->_m1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))
                *(((((this->_m9p0)/(this->_40p0))*g12)
                + ((((this->_m9p0)/(this->_8p0))*g22)
                + ((((this->_1p0)/(this->_2p0))*Y2S)
                + ((((this->_m1p0)/(this->_4p0))*(g22*this->zetaw))
                + (((this->_m1p0)/(this->_4p0))
                *((g22 + (((this->_3p0)/(this->_5p0))*g12))*this->zetaz))))))
                *Zpp)))));
    beta_g[15] = this->_0p0;//PLACEHOLDER
    beta_g[16] = this->_0p0;
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapp_extended(const std::vector<value_type>& g,
                     std::vector<value_type>& beta_g,
                     const value_type t)
{
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    std::vector<value_type> gpp (int(g.size()),this->_0p0);
    this->beta_extended(g,gp,t);
    this->betap_extended(g,gpp,gp,t);

    //Compute beta functions:
    this->betapp_extended(g,beta_g,gp,gpp,t);
}
template<class value_type>
void betaSMdS_1loop<value_type>::betappp_extended(const std::vector<value_type>& g,
                      std::vector<value_type>& beta_g,
                      const std::vector<value_type>& gp,
                      const std::vector<value_type>& gpp,
                      const std::vector<value_type>& gppp,
                      const value_type t)
{
    if(int(g.size()) < 17)
    {
        throw("Coefficient array too small.");
    }

    //Fill in the basics:
    this->betapp_basic(g,beta_g,gp,gpp,t);


    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    value_type Zp = gp[14];
    value_type tclp = gp[15];

    //Coupling derivatives:
    value_type lpp = gpp[0];
    value_type g12pp = gpp[1];
    value_type g22pp = gpp[2];
    value_type g32pp = gpp[3];
    value_type m2pp = gpp[4];
    value_type ye2pp = gpp[5];
    value_type ymu2pp = gpp[6];
    value_type ytau2pp = gpp[7];
    value_type yu2pp = gpp[8];
    value_type yd2pp = gpp[9];
    value_type yc2pp = gpp[10];
    value_type ys2pp = gpp[11];
    value_type yt2pp = gpp[12];
    value_type yb2pp = gpp[13];

    value_type Zpp = gpp[14];
    value_type tclpp = gpp[15];

    //Coupling derivatives:
    value_type lppp = gppp[0];
    value_type g12ppp = gppp[1];
    value_type g22ppp = gppp[2];
    value_type g32ppp = gppp[3];
    value_type m2ppp = gppp[4];
    value_type ye2ppp = gppp[5];
    value_type ymu2ppp = gppp[6];
    value_type ytau2ppp = gppp[7];
    value_type yu2ppp = gppp[8];
    value_type yd2ppp = gppp[9];
    value_type yc2ppp = gppp[10];
    value_type ys2ppp = gppp[11];
    value_type yt2ppp = gppp[12];
    value_type yb2ppp = gppp[13];

    value_type Zppp = gppp[14];
    value_type tclppp = gppp[15];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type Y2Spp = (ye2pp + (ymu2pp + ((this->_3p0*(yb2pp + (yd2pp + ys2pp))) + (ytau2pp + (this->_3p0*(yc2pp + (yt2pp + yu2pp)))))));
    value_type Y2Sppp = (ye2ppp + (ymu2ppp + ((this->_3p0*(yb2ppp + (yd2ppp + ys2ppp))) + (ytau2ppp + (this->_3p0*(yc2ppp + (yt2ppp + yu2ppp)))))));

    //Add the extended functions:
    beta_g[14] = ((((this->_m1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))
                *(Z*((((this->_m9p0)/(this->_40p0))*g12ppp)
                + ((((this->_m9p0)/(this->_8p0))*g22ppp)
                + ((((this->_1p0)/(this->_2p0))*Y2Sppp)
                + ((((this->_m1p0)/(this->_4p0))*(g22ppp*this->zetaw))
                + (((this->_m1p0)/(this->_4p0))
                *((g22ppp + (((this->_3p0)/(this->_5p0))*g12ppp))*this->zetaz)))))))))
                + ((((this->_m3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))
                *(((((this->_m9p0)/(this->_40p0))*g12pp)
                + ((((this->_m9p0)/(this->_8p0))*g22pp)
                + ((((this->_1p0)/(this->_2p0))*Y2Spp)
                + ((((this->_m1p0)/(this->_4p0))*(g22pp*this->zetaw))
                + (((this->_m1p0)/(this->_4p0))
                *((g22pp + (((this->_3p0)/(this->_5p0))*g12pp))
                *this->zetaz))))))*Zp))) + ((((this->_m3p0)/(this->_16p0))
                *(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_40p0))*g12p)
                + ((((this->_m9p0)/(this->_8p0))*g22p)
                + ((((this->_1p0)/(this->_2p0))*Y2Sp)
                + ((((this->_m1p0)/(this->_4p0))*(g22p*this->zetaw))
                + (((this->_m1p0)/(this->_4p0))
                *((g22p + (((this->_3p0)/(this->_5p0))*g12p))*this->zetaz))))))*Zpp)))
                + (((this->_m1p0)/(this->_16p0))
                *(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_40p0))*g12)
                + ((((this->_m9p0)/(this->_8p0))*g22)
                + ((((this->_1p0)/(this->_2p0))*Y2S)
                + ((((this->_m1p0)/(this->_4p0))*(g22*this->zetaw))
                + (((this->_m1p0)/(this->_4p0))*((g22
                + (((this->_3p0)/(this->_5p0))*g12))*this->zetaz))))))*Zppp))))));
    beta_g[15] = this->_0p0;//PLACEHOLDER
    beta_g[16] = this->_0p0;
}
template<class value_type>
void betaSMdS_1loop<value_type>::betappp_extended(const std::vector<value_type>& g,
                      std::vector<value_type>& beta_g,const value_type t)
{
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    std::vector<value_type> gpp (int(g.size()),this->_0p0);
    std::vector<value_type> gppp (int(g.size()),this->_0p0);
    this->beta_extended(g,gp,t);
    this->betap_extended(g,gpp,gp,t);
    this->betapp_extended(g,gppp,gp,gpp,t);

    //Compute beta functions:
    this->betappp_extended(g,beta_g,gp,gpp,gppp,t);
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapppp_extended(const std::vector<value_type>& g,
                       std::vector<value_type>& beta_g,
                       const std::vector<value_type>& gp,
                       const std::vector<value_type>& gpp,
                       const std::vector<value_type>& gppp,
                       const std::vector<value_type>& gpppp,
                       const value_type t)
{
    if(int(g.size()) < 17)
    {
        throw("Coefficient array too small.");
    }

    //Fill in the basics:
    this->betapp_basic(g,beta_g,gp,gpp,t);


    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    value_type Zp = gp[14];
    value_type tclp = gp[15];

    //Coupling derivatives:
    value_type lpp = gpp[0];
    value_type g12pp = gpp[1];
    value_type g22pp = gpp[2];
    value_type g32pp = gpp[3];
    value_type m2pp = gpp[4];
    value_type ye2pp = gpp[5];
    value_type ymu2pp = gpp[6];
    value_type ytau2pp = gpp[7];
    value_type yu2pp = gpp[8];
    value_type yd2pp = gpp[9];
    value_type yc2pp = gpp[10];
    value_type ys2pp = gpp[11];
    value_type yt2pp = gpp[12];
    value_type yb2pp = gpp[13];

    value_type Zpp = gpp[14];
    value_type tclpp = gpp[15];

    //Coupling derivatives:
    value_type lppp = gppp[0];
    value_type g12ppp = gppp[1];
    value_type g22ppp = gppp[2];
    value_type g32ppp = gppp[3];
    value_type m2ppp = gppp[4];
    value_type ye2ppp = gppp[5];
    value_type ymu2ppp = gppp[6];
    value_type ytau2ppp = gppp[7];
    value_type yu2ppp = gppp[8];
    value_type yd2ppp = gppp[9];
    value_type yc2ppp = gppp[10];
    value_type ys2ppp = gppp[11];
    value_type yt2ppp = gppp[12];
    value_type yb2ppp = gppp[13];

    value_type Zppp = gppp[14];
    value_type tclppp = gppp[15];

    //Coupling derivatives:
    value_type lpppp = gpppp[0];
    value_type g12pppp = gpppp[1];
    value_type g22pppp = gpppp[2];
    value_type g32pppp = gpppp[3];
    value_type m2pppp = gpppp[4];
    value_type ye2pppp = gpppp[5];
    value_type ymu2pppp = gpppp[6];
    value_type ytau2pppp = gpppp[7];
    value_type yu2pppp = gpppp[8];
    value_type yd2pppp = gpppp[9];
    value_type yc2pppp = gpppp[10];
    value_type ys2pppp = gpppp[11];
    value_type yt2pppp = gpppp[12];
    value_type yb2pppp = gpppp[13];

    value_type Zpppp = gpppp[14];
    value_type tclpppp = gpppp[15];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type Y2Spp = (ye2pp + (ymu2pp + ((this->_3p0*(yb2pp + (yd2pp + ys2pp))) + (ytau2pp + (this->_3p0*(yc2pp + (yt2pp + yu2pp)))))));
    value_type Y2Sppp = (ye2ppp + (ymu2ppp + ((this->_3p0*(yb2ppp + (yd2ppp + ys2ppp))) + (ytau2ppp + (this->_3p0*(yc2ppp + (yt2ppp + yu2ppp)))))));
    value_type Y2Spppp = (ye2pppp + (ymu2pppp + ((this->_3p0*(yb2pppp + (yd2pppp + ys2pppp))) + (ytau2pppp + (this->_3p0*(yc2pppp + (yt2pppp + yu2pppp)))))));

    //Add the extended functions:
    beta_g[14] = ((((this->_m1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))
                    *(Z*((((this->_m9p0)/(this->_40p0))*g12pppp)
                    + ((((this->_m9p0)/(this->_8p0))*g22pppp)
                    + ((((this->_1p0)/(this->_2p0))*Y2Spppp)
                    + ((((this->_m1p0)/(this->_4p0))*(g22pppp*this->zetaw))
                    + (((this->_m1p0)/(this->_4p0))
                    *((g22pppp + (((this->_3p0)/(this->_5p0))*g12pppp))
                    *this->zetaz))))))))) + ((((this->_m1p0)/(this->_4p0))
                    *(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_40p0))*g12ppp)
                    + ((((this->_m9p0)/(this->_8p0))*g22ppp)
                    + ((((this->_1p0)/(this->_2p0))*Y2Sppp)
                    + ((((this->_m1p0)/(this->_4p0))*(g22ppp*this->zetaw))
                    + (((this->_m1p0)/(this->_4p0))
                    *((g22ppp + (((this->_3p0)/(this->_5p0))*g12ppp))
                    *this->zetaz))))))*Zp))) + ((((this->_m3p0)/(this->_8p0))
                    *(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_40p0))*g12pp)
                    + ((((this->_m9p0)/(this->_8p0))*g22pp)
                    + ((((this->_1p0)/(this->_2p0))*Y2Spp)
                    + ((((this->_m1p0)/(this->_4p0))*(g22pp*this->zetaw))
                    + (((this->_m1p0)/(this->_4p0))
                    *((g22pp + (((this->_3p0)/(this->_5p0))*g12pp))
                    *this->zetaz))))))*Zpp))) + ((((this->_m1p0)/(this->_4p0))
                    *(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_40p0))*g12p)
                    + ((((this->_m9p0)/(this->_8p0))*g22p)
                    + ((((this->_1p0)/(this->_2p0))*Y2Sp)
                    + ((((this->_m1p0)/(this->_4p0))*(g22p*this->zetaw))
                    + (((this->_m1p0)/(this->_4p0))
                    *((g22p + (((this->_3p0)/(this->_5p0))*g12p))
                    *this->zetaz))))))*Zppp))) + (((this->_m1p0)/(this->_16p0))
                    *(((this->_1p0)/(PI*PI))*(((((this->_m9p0)/(this->_40p0))*g12)
                    + ((((this->_m9p0)/(this->_8p0))*g22)
                    + ((((this->_1p0)/(this->_2p0))*Y2S)
                    + ((((this->_m1p0)/(this->_4p0))*(g22*this->zetaw))
                    + (((this->_m1p0)/(this->_4p0))
                    *((g22 + (((this->_3p0)/(this->_5p0))*g12))
                    *this->zetaz))))))*Zpppp)))))));
    beta_g[15] = this->_0p0;//PLACEHOLDER
    beta_g[16] = this->_0p0;
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapppp_extended(const std::vector<value_type>& g,
                       std::vector<value_type>& beta_g,const value_type t)
{
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    std::vector<value_type> gpp (int(g.size()),this->_0p0);
    std::vector<value_type> gppp (int(g.size()),this->_0p0);
    std::vector<value_type> gpppp (int(g.size()),this->_0p0);
    this->beta_extended(g,gp,t);
    this->betap_extended(g,gpp,gp,t);
    this->betapp_extended(g,gppp,gp,gpp,t);
    this->betappp_extended(g,gpppp,gp,gpp,gppp,t);

    //Compute beta functions:
    this->betapppp_extended(g,beta_g,gp,gpp,gppp,gpppp,t);
}
template<class value_type>
void betaSMdS_1loop<value_type>::beta_grav(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const value_type t)
{
    //Fill in the basics:
    this->beta_extended(g,beta_g,t);
    if(int(g.size()) < 23)
    {
        throw("Coefficient array too small.");
    }

    //couplings
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];

    //value_type t = g[16];//Given by the function argument.

    value_type xi = g[17];
    value_type V0 = g[18];
    value_type kappa = g[19];
    value_type alpha1 = g[20];
    value_type alpha2 = g[21];
    value_type alpha3 = g[22];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));

    //Add the extended functions:
    beta_g[17] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*((((this->_m1p0)/(this->_6p0)) + xi)*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S))))));
    beta_g[18] = (((this->_1p0)/(this->_16p0))*(m2*m2*((this->_1p0)/(PI*PI))));
    beta_g[19] = (((this->_m1p0)/(this->_8p0))*(m2*(((this->_1p0)/(PI*PI))*(((this->_m1p0)/(this->_6p0)) + xi))));
    beta_g[20] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(((this->_m277p0)/(this->_288p0)) + ((((this->_m1p0)/(this->_3p0))*xi) + xi*xi))));
    beta_g[21] = (((this->_571p0)/(this->_2880p0))*((this->_1p0)/(PI*PI)));
    beta_g[22] = (((this->_m293p0)/(this->_23040p0))*((this->_1p0)/(PI*PI)));

}

template<class value_type>
void betaSMdS_1loop<value_type>::betap_grav(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const std::vector<value_type>& gp,
                   const value_type t)
{
    //Fill in the basics:
    if(int(g.size()) < 23)
    {
        throw("Coefficient array too small.");
    }

    //this->beta_grav(g,gp,t);

    this->betap_extended(g,beta_g,gp,t);

    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];

    value_type xi = g[17];
    value_type V0 = g[18];
    value_type kappa = g[19];
    value_type alpha1 = g[20];
    value_type alpha2 = g[21];
    value_type alpha3 = g[22];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    value_type Zp = gp[14];
    value_type tclp = gp[15];

    value_type xip = gp[17];
    value_type V0p = gp[18];
    value_type kappap = gp[19];
    value_type alpha1p = gp[20];
    value_type alpha2p = gp[21];
    value_type alpha3p = gp[22];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));
    value_type HSp = ((this->_2p0*(ye2*ye2p)) + ((this->_2p0*(ymu2*ymu2p)) + ((this->_3p0*((this->_2p0*(yb2*yb2p)) + ((this->_2p0*(yd2*yd2p)) + (this->_2p0*(ys2*ys2p))))) + ((this->_2p0*(ytau2*ytau2p)) + (this->_3p0*((this->_2p0*(yc2*yc2p)) + ((this->_2p0*(yt2*yt2p)) + (this->_2p0*(yu2*yu2p)))))))));


    //Add the extended functions:
    beta_g[17] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(xip*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S)))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*((((this->_m1p0)/(this->_6p0)) + xi)*((((this->_m9p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_6p0*lp) + Y2Sp)))))));
    beta_g[18] = (((this->_1p0)/(this->_8p0))*(m2*(m2p*((this->_1p0)/(PI*PI)))));
    beta_g[19] = ((((this->_m1p0)/(this->_8p0))*(m2p*(((this->_1p0)/(PI*PI))*(((this->_m1p0)/(this->_6p0)) + xi)))) + (((this->_m1p0)/(this->_8p0))*(m2*(((this->_1p0)/(PI*PI))*xip))));
    beta_g[20] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*((((this->_m1p0)/(this->_3p0))*xip) + (this->_2p0*(xi*xip)))));
    beta_g[21] = this->_0p0;
    beta_g[22] = this->_0p0;

}
template<class value_type>
void betaSMdS_1loop<value_type>::betap_grav(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const value_type t)
{
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    this->beta_grav(g,gp,t);
    this->betap_grav(g,beta_g,gp,t);
}

template<class value_type>
void betaSMdS_1loop<value_type>::betapp_grav(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const std::vector<value_type>& gp,
                   const std::vector<value_type>& gpp,
                   const value_type t)
{
    //Fill in the basics:
    if(int(g.size()) < 23)
    {
        throw("Coefficient array too small.");
    }

    //this->beta_grav(g,gp,t);

    this->betapp_extended(g,beta_g,gp,gpp,t);

    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];

    value_type xi = g[17];
    value_type V0 = g[18];
    value_type kappa = g[19];
    value_type alpha1 = g[20];
    value_type alpha2 = g[21];
    value_type alpha3 = g[22];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    value_type Zp = gp[14];
    value_type tclp = gp[15];

    value_type xip = gp[17];
    value_type V0p = gp[18];
    value_type kappap = gp[19];
    value_type alpha1p = gp[20];
    value_type alpha2p = gp[21];
    value_type alpha3p = gp[22];

    //Couppling derivatives:
    value_type lpp = gpp[0];
    value_type g12pp = gpp[1];
    value_type g22pp = gpp[2];
    value_type g32pp = gpp[3];
    value_type m2pp = gpp[4];
    value_type ye2pp = gpp[5];
    value_type ymu2pp = gpp[6];
    value_type ytau2pp = gpp[7];
    value_type yu2pp = gpp[8];
    value_type yd2pp = gpp[9];
    value_type yc2pp = gpp[10];
    value_type ys2pp = gpp[11];
    value_type yt2pp = gpp[12];
    value_type yb2pp = gpp[13];

    value_type Zpp = gpp[14];
    value_type tclpp = gpp[15];

    value_type xipp = gpp[17];
    value_type V0pp = gpp[18];
    value_type kappapp = gpp[19];
    value_type alpha1pp = gpp[20];
    value_type alpha2pp = gpp[21];
    value_type alpha3pp = gpp[22];



    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type Y2Spp = (ye2pp + (ymu2pp + ((this->_3p0*(yb2pp + (yd2pp + ys2pp))) + (ytau2pp + (this->_3p0*(yc2pp + (yt2pp + yu2pp)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));
    value_type HSp = ((this->_2p0*(ye2*ye2p)) + ((this->_2p0*(ymu2*ymu2p)) + ((this->_3p0*((this->_2p0*(yb2*yb2p)) + ((this->_2p0*(yd2*yd2p)) + (this->_2p0*(ys2*ys2p))))) + ((this->_2p0*(ytau2*ytau2p)) + (this->_3p0*((this->_2p0*(yc2*yc2p)) + ((this->_2p0*(yt2*yt2p)) + (this->_2p0*(yu2*yu2p)))))))));
    value_type HSpp = ((this->_2p0*ye2p*ye2p) + ((this->_2p0*(ye2*ye2pp)) + ((this->_2p0*ymu2p*ymu2p) + ((this->_2p0*(ymu2*ymu2pp)) + ((this->_3p0*((this->_2p0*yb2p*yb2p) + ((this->_2p0*(yb2*yb2pp)) + ((this->_2p0*yd2p*yd2p) + ((this->_2p0*(yd2*yd2pp)) + ((this->_2p0*ys2p*ys2p) + (this->_2p0*(ys2*ys2pp)))))))) + ((this->_2p0*ytau2p*ytau2p) + ((this->_2p0*(ytau2*ytau2pp)) + (this->_3p0*((this->_2p0*yc2p*yc2p) + ((this->_2p0*(yc2*yc2pp)) + ((this->_2p0*yt2p*yt2p) + ((this->_2p0*(yt2*yt2pp)) + ((this->_2p0*yu2p*yu2p) + (this->_2p0*(yu2*yu2pp)))))))))))))));

    //Add the extended functions:
    beta_g[17] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(xipp*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S)))))) + ((((this->_1p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(xip*((((this->_m9p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_6p0*lp) + Y2Sp)))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*((((this->_m1p0)/(this->_6p0)) + xi)*((((this->_m9p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_6p0*lpp) + Y2Spp))))))));
    beta_g[18] = ((((this->_1p0)/(this->_8p0))*(m2p*m2p*((this->_1p0)/(PI*PI)))) + (((this->_1p0)/(this->_8p0))*(m2*(m2pp*((this->_1p0)/(PI*PI))))));
    beta_g[19] = ((((this->_m1p0)/(this->_8p0))*(m2pp*(((this->_1p0)/(PI*PI))*(((this->_m1p0)/(this->_6p0)) + xi)))) + ((((this->_m1p0)/(this->_4p0))*(m2p*(((this->_1p0)/(PI*PI))*xip))) + (((this->_m1p0)/(this->_8p0))*(m2*(((this->_1p0)/(PI*PI))*xipp)))));
    beta_g[20] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*((this->_2p0*xip*xip) + ((((this->_m1p0)/(this->_3p0))*xipp) + (this->_2p0*(xi*xipp))))));
    beta_g[21] = this->_0p0;
    beta_g[22] = this->_0p0;

}
template<class value_type>
void betaSMdS_1loop<value_type>::betapp_grav(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const value_type t)
{
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    std::vector<value_type> gpp (int(g.size()),this->_0p0);
    this->beta_grav(g,gp,t);
    this->betap_grav(g,gpp,gp,t);

    this->betapp_grav(g,beta_g,gp,gpp,t);
}

template<class value_type>
void betaSMdS_1loop<value_type>::betappp_grav(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const std::vector<value_type>& gp,
                   const std::vector<value_type>& gpp,
                   const std::vector<value_type>& gppp,
                   const value_type t)
{
    //Fill in the basics:
    if(int(g.size()) < 23)
    {
        throw("Coefficient array too small.");
    }

    //this->beta_grav(g,gp,t);

    this->betappp_extended(g,beta_g,gp,gpp,gppp,t);

    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];

    value_type xi = g[17];
    value_type V0 = g[18];
    value_type kappa = g[19];
    value_type alpha1 = g[20];
    value_type alpha2 = g[21];
    value_type alpha3 = g[22];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    value_type Zp = gp[14];
    value_type tclp = gp[15];

    value_type xip = gp[17];
    value_type V0p = gp[18];
    value_type kappap = gp[19];
    value_type alpha1p = gp[20];
    value_type alpha2p = gp[21];
    value_type alpha3p = gp[22];

    //Couppling derivatives:
    value_type lpp = gpp[0];
    value_type g12pp = gpp[1];
    value_type g22pp = gpp[2];
    value_type g32pp = gpp[3];
    value_type m2pp = gpp[4];
    value_type ye2pp = gpp[5];
    value_type ymu2pp = gpp[6];
    value_type ytau2pp = gpp[7];
    value_type yu2pp = gpp[8];
    value_type yd2pp = gpp[9];
    value_type yc2pp = gpp[10];
    value_type ys2pp = gpp[11];
    value_type yt2pp = gpp[12];
    value_type yb2pp = gpp[13];

    value_type Zpp = gpp[14];
    value_type tclpp = gpp[15];

    value_type xipp = gpp[17];
    value_type V0pp = gpp[18];
    value_type kappapp = gpp[19];
    value_type alpha1pp = gpp[20];
    value_type alpha2pp = gpp[21];
    value_type alpha3pp = gpp[22];

    //Coupppling derivatives:
    value_type lppp = gppp[0];
    value_type g12ppp = gppp[1];
    value_type g22ppp = gppp[2];
    value_type g32ppp = gppp[3];
    value_type m2ppp = gppp[4];
    value_type ye2ppp = gppp[5];
    value_type ymu2ppp = gppp[6];
    value_type ytau2ppp = gppp[7];
    value_type yu2ppp = gppp[8];
    value_type yd2ppp = gppp[9];
    value_type yc2ppp = gppp[10];
    value_type ys2ppp = gppp[11];
    value_type yt2ppp = gppp[12];
    value_type yb2ppp = gppp[13];

    value_type Zppp = gppp[14];
    value_type tclppp = gppp[15];

    value_type xippp = gppp[17];
    value_type V0ppp = gppp[18];
    value_type kapppappp = gppp[19];
    value_type alpha1ppp = gppp[20];
    value_type alpha2ppp = gppp[21];
    value_type alpha3ppp = gppp[22];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type Y2Spp = (ye2pp + (ymu2pp + ((this->_3p0*(yb2pp + (yd2pp + ys2pp))) + (ytau2pp + (this->_3p0*(yc2pp + (yt2pp + yu2pp)))))));
    value_type Y2Sppp = (ye2ppp + (ymu2ppp + ((this->_3p0*(yb2ppp + (yd2ppp + ys2ppp))) + (ytau2ppp + (this->_3p0*(yc2ppp + (yt2ppp + yu2ppp)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));
    value_type HSp = ((this->_2p0*(ye2*ye2p)) + ((this->_2p0*(ymu2*ymu2p)) + ((this->_3p0*((this->_2p0*(yb2*yb2p)) + ((this->_2p0*(yd2*yd2p)) + (this->_2p0*(ys2*ys2p))))) + ((this->_2p0*(ytau2*ytau2p)) + (this->_3p0*((this->_2p0*(yc2*yc2p)) + ((this->_2p0*(yt2*yt2p)) + (this->_2p0*(yu2*yu2p)))))))));
    value_type HSpp = ((this->_2p0*ye2p*ye2p) + ((this->_2p0*(ye2*ye2pp)) + ((this->_2p0*ymu2p*ymu2p) + ((this->_2p0*(ymu2*ymu2pp)) + ((this->_3p0*((this->_2p0*yb2p*yb2p) + ((this->_2p0*(yb2*yb2pp)) + ((this->_2p0*yd2p*yd2p) + ((this->_2p0*(yd2*yd2pp)) + ((this->_2p0*ys2p*ys2p) + (this->_2p0*(ys2*ys2pp)))))))) + ((this->_2p0*ytau2p*ytau2p) + ((this->_2p0*(ytau2*ytau2pp)) + (this->_3p0*((this->_2p0*yc2p*yc2p) + ((this->_2p0*(yc2*yc2pp)) + ((this->_2p0*yt2p*yt2p) + ((this->_2p0*(yt2*yt2pp)) + ((this->_2p0*yu2p*yu2p) + (this->_2p0*(yu2*yu2pp)))))))))))))));
    value_type HSppp = ((this->_6p0*(ye2p*ye2pp)) + ((this->_2p0*(ye2*ye2ppp)) + ((this->_6p0*(ymu2p*ymu2pp)) + ((this->_2p0*(ymu2*ymu2ppp)) + ((this->_3p0*((this->_6p0*(yb2p*yb2pp)) + ((this->_2p0*(yb2*yb2ppp)) + ((this->_6p0*(yd2p*yd2pp)) + ((this->_2p0*(yd2*yd2ppp)) + ((this->_6p0*(ys2p*ys2pp)) + (this->_2p0*(ys2*ys2ppp)))))))) + ((this->_6p0*(ytau2p*ytau2pp)) + ((this->_2p0*(ytau2*ytau2ppp)) + (this->_3p0*((this->_6p0*(yc2p*yc2pp)) + ((this->_2p0*(yc2*yc2ppp)) + ((this->_6p0*(yt2p*yt2pp)) + ((this->_2p0*(yt2*yt2ppp)) + ((this->_6p0*(yu2p*yu2pp)) + (this->_2p0*(yu2*yu2ppp)))))))))))))));


    //Add the extended functions:
    beta_g[17] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(xippp*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S)))))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(xipp*((((this->_m9p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_6p0*lp) + Y2Sp)))))) + ((((this->_3p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(xip*((((this->_m9p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_6p0*lpp) + Y2Spp)))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*((((this->_m1p0)/(this->_6p0)) + xi)*((((this->_m9p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_6p0*lppp) + Y2Sppp)))))))));
    beta_g[18] = ((((this->_3p0)/(this->_8p0))*(m2p*(m2pp*((this->_1p0)/(PI*PI))))) + (((this->_1p0)/(this->_8p0))*(m2*(m2ppp*((this->_1p0)/(PI*PI))))));
    beta_g[19] = ((((this->_m1p0)/(this->_8p0))*(m2ppp*(((this->_1p0)/(PI*PI))*(((this->_m1p0)/(this->_6p0)) + xi)))) + ((((this->_m3p0)/(this->_8p0))*(m2pp*(((this->_1p0)/(PI*PI))*xip))) + ((((this->_m3p0)/(this->_8p0))*(m2p*(((this->_1p0)/(PI*PI))*xipp))) + (((this->_m1p0)/(this->_8p0))*(m2*(((this->_1p0)/(PI*PI))*xippp))))));
    beta_g[20] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*((this->_6p0*(xip*xipp)) + ((((this->_m1p0)/(this->_3p0))*xippp) + (this->_2p0*(xi*xippp))))));
    beta_g[21] = this->_0p0;
    beta_g[22] = this->_0p0;
}
template<class value_type>
void betaSMdS_1loop<value_type>::betappp_grav(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const value_type t)
{
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    std::vector<value_type> gpp (int(g.size()),this->_0p0);
    std::vector<value_type> gppp (int(g.size()),this->_0p0);
    this->beta_grav(g,gp,t);
    this->betap_grav(g,gpp,gp,t);
    this->betapp_grav(g,gppp,gp,gpp,t);

    this->betappp_grav(g,beta_g,gp,gpp,gppp,t);
}

template<class value_type>
void betaSMdS_1loop<value_type>::betapppp_grav(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const std::vector<value_type>& gp,
                   const std::vector<value_type>& gpp,
                   const std::vector<value_type>& gppp,
                   const std::vector<value_type>& gpppp,
                   const value_type t)
{
    //Fill in the basics:
    if(int(g.size()) < 23)
    {
        throw("Coefficient array too small.");
    }

    //this->beta_grav(g,gp,t);

    this->betapppp_extended(g,beta_g,gp,gpp,gppp,gpppp,t);

    //couplings:
    value_type l = g[0];
    value_type g12 = g[1];
    value_type g22 = g[2];
    value_type g32 = g[3];
    value_type m2 = g[4];
    value_type ye2 = g[5];
    value_type ymu2 = g[6];
    value_type ytau2 = g[7];
    value_type yu2 = g[8];
    value_type yd2 = g[9];
    value_type yc2 = g[10];
    value_type ys2 = g[11];
    value_type yt2 = g[12];
    value_type yb2 = g[13];

    value_type Z = g[14];
    value_type tcl = g[15];

    value_type xi = g[17];
    value_type V0 = g[18];
    value_type kappa = g[19];
    value_type alpha1 = g[20];
    value_type alpha2 = g[21];
    value_type alpha3 = g[22];

    //Coupling derivatives:
    value_type lp = gp[0];
    value_type g12p = gp[1];
    value_type g22p = gp[2];
    value_type g32p = gp[3];
    value_type m2p = gp[4];
    value_type ye2p = gp[5];
    value_type ymu2p = gp[6];
    value_type ytau2p = gp[7];
    value_type yu2p = gp[8];
    value_type yd2p = gp[9];
    value_type yc2p = gp[10];
    value_type ys2p = gp[11];
    value_type yt2p = gp[12];
    value_type yb2p = gp[13];

    value_type Zp = gp[14];
    value_type tclp = gp[15];

    value_type xip = gp[17];
    value_type V0p = gp[18];
    value_type kappap = gp[19];
    value_type alpha1p = gp[20];
    value_type alpha2p = gp[21];
    value_type alpha3p = gp[22];

    //Couppling derivatives:
    value_type lpp = gpp[0];
    value_type g12pp = gpp[1];
    value_type g22pp = gpp[2];
    value_type g32pp = gpp[3];
    value_type m2pp = gpp[4];
    value_type ye2pp = gpp[5];
    value_type ymu2pp = gpp[6];
    value_type ytau2pp = gpp[7];
    value_type yu2pp = gpp[8];
    value_type yd2pp = gpp[9];
    value_type yc2pp = gpp[10];
    value_type ys2pp = gpp[11];
    value_type yt2pp = gpp[12];
    value_type yb2pp = gpp[13];

    value_type Zpp = gpp[14];
    value_type tclpp = gpp[15];

    value_type xipp = gpp[17];
    value_type V0pp = gpp[18];
    value_type kappapp = gpp[19];
    value_type alpha1pp = gpp[20];
    value_type alpha2pp = gpp[21];
    value_type alpha3pp = gpp[22];

    //Coupppling derivatives:
    value_type lppp = gppp[0];
    value_type g12ppp = gppp[1];
    value_type g22ppp = gppp[2];
    value_type g32ppp = gppp[3];
    value_type m2ppp = gppp[4];
    value_type ye2ppp = gppp[5];
    value_type ymu2ppp = gppp[6];
    value_type ytau2ppp = gppp[7];
    value_type yu2ppp = gppp[8];
    value_type yd2ppp = gppp[9];
    value_type yc2ppp = gppp[10];
    value_type ys2ppp = gppp[11];
    value_type yt2ppp = gppp[12];
    value_type yb2ppp = gppp[13];

    value_type Zppp = gppp[14];
    value_type tclppp = gppp[15];

    value_type xippp = gppp[17];
    value_type V0ppp = gppp[18];
    value_type kapppappp = gppp[19];
    value_type alpha1ppp = gppp[20];
    value_type alpha2ppp = gppp[21];
    value_type alpha3ppp = gppp[22];

    //Couppppling derivatives:
    value_type lpppp = gpppp[0];
    value_type g12pppp = gpppp[1];
    value_type g22pppp = gpppp[2];
    value_type g32pppp = gpppp[3];
    value_type m2pppp = gpppp[4];
    value_type ye2pppp = gpppp[5];
    value_type ymu2pppp = gpppp[6];
    value_type ytau2pppp = gpppp[7];
    value_type yu2pppp = gpppp[8];
    value_type yd2pppp = gpppp[9];
    value_type yc2pppp = gpppp[10];
    value_type ys2pppp = gpppp[11];
    value_type yt2pppp = gpppp[12];
    value_type yb2pppp = gpppp[13];

    value_type Zpppp = gpppp[14];
    value_type tclpppp = gpppp[15];

    value_type xipppp = gpppp[17];
    value_type V0pppp = gpppp[18];
    value_type kappppapppp = gpppp[19];
    value_type alpha1pppp = gpppp[20];
    value_type alpha2pppp = gpppp[21];
    value_type alpha3pppp = gpppp[22];

    value_type Y2S = (ye2 + (ymu2 + ((this->_3p0*(yb2 + (yd2 + ys2))) + (ytau2 + (this->_3p0*(yc2 + (yt2 + yu2)))))));
    value_type Y2Sp = (ye2p + (ymu2p + ((this->_3p0*(yb2p + (yd2p + ys2p))) + (ytau2p + (this->_3p0*(yc2p + (yt2p + yu2p)))))));
    value_type Y2Spp = (ye2pp + (ymu2pp + ((this->_3p0*(yb2pp + (yd2pp + ys2pp))) + (ytau2pp + (this->_3p0*(yc2pp + (yt2pp + yu2pp)))))));
    value_type Y2Sppp = (ye2ppp + (ymu2ppp + ((this->_3p0*(yb2ppp + (yd2ppp + ys2ppp))) + (ytau2ppp + (this->_3p0*(yc2ppp + (yt2ppp + yu2ppp)))))));
    value_type Y2Spppp = (ye2pppp + (ymu2pppp + ((this->_3p0*(yb2pppp + (yd2pppp + ys2pppp))) + (ytau2pppp + (this->_3p0*(yc2pppp + (yt2pppp + yu2pppp)))))));
    value_type HS = (ye2*ye2 + (ymu2*ymu2 + ((this->_3p0*(yb2*yb2 + (yd2*yd2 + ys2*ys2))) + (ytau2*ytau2 + (this->_3p0*(yc2*yc2 + (yt2*yt2 + yu2*yu2)))))));
    value_type HSp = ((this->_2p0*(ye2*ye2p)) + ((this->_2p0*(ymu2*ymu2p)) + ((this->_3p0*((this->_2p0*(yb2*yb2p)) + ((this->_2p0*(yd2*yd2p)) + (this->_2p0*(ys2*ys2p))))) + ((this->_2p0*(ytau2*ytau2p)) + (this->_3p0*((this->_2p0*(yc2*yc2p)) + ((this->_2p0*(yt2*yt2p)) + (this->_2p0*(yu2*yu2p)))))))));
    value_type HSpp = ((this->_2p0*ye2p*ye2p) + ((this->_2p0*(ye2*ye2pp)) + ((this->_2p0*ymu2p*ymu2p) + ((this->_2p0*(ymu2*ymu2pp)) + ((this->_3p0*((this->_2p0*yb2p*yb2p) + ((this->_2p0*(yb2*yb2pp)) + ((this->_2p0*yd2p*yd2p) + ((this->_2p0*(yd2*yd2pp)) + ((this->_2p0*ys2p*ys2p) + (this->_2p0*(ys2*ys2pp)))))))) + ((this->_2p0*ytau2p*ytau2p) + ((this->_2p0*(ytau2*ytau2pp)) + (this->_3p0*((this->_2p0*yc2p*yc2p) + ((this->_2p0*(yc2*yc2pp)) + ((this->_2p0*yt2p*yt2p) + ((this->_2p0*(yt2*yt2pp)) + ((this->_2p0*yu2p*yu2p) + (this->_2p0*(yu2*yu2pp)))))))))))))));
    value_type HSppp = ((this->_6p0*(ye2p*ye2pp)) + ((this->_2p0*(ye2*ye2ppp)) + ((this->_6p0*(ymu2p*ymu2pp)) + ((this->_2p0*(ymu2*ymu2ppp)) + ((this->_3p0*((this->_6p0*(yb2p*yb2pp)) + ((this->_2p0*(yb2*yb2ppp)) + ((this->_6p0*(yd2p*yd2pp)) + ((this->_2p0*(yd2*yd2ppp)) + ((this->_6p0*(ys2p*ys2pp)) + (this->_2p0*(ys2*ys2ppp)))))))) + ((this->_6p0*(ytau2p*ytau2pp)) + ((this->_2p0*(ytau2*ytau2ppp)) + (this->_3p0*((this->_6p0*(yc2p*yc2pp)) + ((this->_2p0*(yc2*yc2ppp)) + ((this->_6p0*(yt2p*yt2pp)) + ((this->_2p0*(yt2*yt2ppp)) + ((this->_6p0*(yu2p*yu2pp)) + (this->_2p0*(yu2*yu2ppp)))))))))))))));
    value_type HSpppp = ((this->_6p0*ye2pp*ye2pp) + ((this->_8p0*(ye2p*ye2ppp)) + ((this->_2p0*(ye2*ye2pppp)) + ((this->_6p0*ymu2pp*ymu2pp) + ((this->_8p0*(ymu2p*ymu2ppp)) + ((this->_2p0*(ymu2*ymu2pppp)) + ((this->_3p0*((this->_6p0*yb2pp*yb2pp) + ((this->_8p0*(yb2p*yb2ppp)) + ((this->_2p0*(yb2*yb2pppp)) + ((this->_6p0*yd2pp*yd2pp) + ((this->_8p0*(yd2p*yd2ppp)) + ((this->_2p0*(yd2*yd2pppp)) + ((this->_6p0*ys2pp*ys2pp) + ((this->_8p0*(ys2p*ys2ppp)) + (this->_2p0*(ys2*ys2pppp))))))))))) + ((this->_6p0*ytau2pp*ytau2pp) + ((this->_8p0*(ytau2p*ytau2ppp)) + ((this->_2p0*(ytau2*ytau2pppp)) + (this->_3p0*((this->_6p0*yc2pp*yc2pp) + ((this->_8p0*(yc2p*yc2ppp)) + ((this->_2p0*(yc2*yc2pppp)) + ((this->_6p0*yt2pp*yt2pp) + ((this->_8p0*(yt2p*yt2ppp)) + ((this->_2p0*(yt2*yt2pppp)) + ((this->_6p0*yu2pp*yu2pp) + ((this->_8p0*(yu2p*yu2ppp)) + (this->_2p0*(yu2*yu2pppp)))))))))))))))))))));


    //Add the extended functions:
    beta_g[17] = ((((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*(xipppp*((((this->_m9p0)/(this->_20p0))*g12) + ((((this->_m9p0)/(this->_4p0))*g22) + ((this->_6p0*l) + Y2S)))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(xippp*((((this->_m9p0)/(this->_20p0))*g12p) + ((((this->_m9p0)/(this->_4p0))*g22p) + ((this->_6p0*lp) + Y2Sp)))))) + ((((this->_3p0)/(this->_8p0))*(((this->_1p0)/(PI*PI))*(xipp*((((this->_m9p0)/(this->_20p0))*g12pp) + ((((this->_m9p0)/(this->_4p0))*g22pp) + ((this->_6p0*lpp) + Y2Spp)))))) + ((((this->_1p0)/(this->_4p0))*(((this->_1p0)/(PI*PI))*(xip*((((this->_m9p0)/(this->_20p0))*g12ppp) + ((((this->_m9p0)/(this->_4p0))*g22ppp) + ((this->_6p0*lppp) + Y2Sppp)))))) + (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*((((this->_m1p0)/(this->_6p0)) + xi)*((((this->_m9p0)/(this->_20p0))*g12pppp) + ((((this->_m9p0)/(this->_4p0))*g22pppp) + ((this->_6p0*lpppp) + Y2Spppp))))))))));
    beta_g[18] = ((((this->_3p0)/(this->_8p0))*(m2pp*m2pp*((this->_1p0)/(PI*PI)))) + ((((this->_1p0)/(this->_2p0))*(m2p*(m2ppp*((this->_1p0)/(PI*PI))))) + (((this->_1p0)/(this->_8p0))*(m2*(m2pppp*((this->_1p0)/(PI*PI)))))));
    beta_g[19] = ((((this->_m1p0)/(this->_8p0))*(m2pppp*(((this->_1p0)/(PI*PI))*(((this->_m1p0)/(this->_6p0)) + xi)))) + ((((this->_m1p0)/(this->_2p0))*(m2ppp*(((this->_1p0)/(PI*PI))*xip))) + ((((this->_m3p0)/(this->_4p0))*(m2pp*(((this->_1p0)/(PI*PI))*xipp))) + ((((this->_m1p0)/(this->_2p0))*(m2p*(((this->_1p0)/(PI*PI))*xippp))) + (((this->_m1p0)/(this->_8p0))*(m2*(((this->_1p0)/(PI*PI))*xipppp)))))));
    beta_g[20] = (((this->_1p0)/(this->_16p0))*(((this->_1p0)/(PI*PI))*((this->_6p0*xipp*xipp) + ((this->_8p0*(xip*xippp)) + ((((this->_m1p0)/(this->_3p0))*xipppp) + (this->_2p0*(xi*xipppp)))))));
    beta_g[21] = this->_0p0;
    beta_g[22] = this->_0p0;
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapppp_grav(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,
                   const value_type t)
{
    std::vector<value_type> gp (int(g.size()),this->_0p0);
    std::vector<value_type> gpp (int(g.size()),this->_0p0);
    std::vector<value_type> gppp (int(g.size()),this->_0p0);
    std::vector<value_type> gpppp (int(g.size()),this->_0p0);
    this->beta_grav(g,gp,t);
    this->betap_grav(g,gpp,gp,t);
    this->betapp_grav(g,gppp,gp,gpp,t);
    this->betappp_grav(g,gpppp,gp,gpp,gppp,t);

    this->betapppp_grav(g,beta_g,gp,gpp,gppp,gpppp,t);
}

//Beta functions expressed in tcl co-ordinates:
template<class value_type>
void betaSMdS_1loop<value_type>::beta_tcl
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type tcl)
{
    //Compute original beta functions, in t co-ordinates:
    value_type t = this->eval.t(g);
    if(this->nType ==4)
    {
        this->beta_grav(g,beta_g,t);
    }
    else
    {
        this->beta_extended(g,beta_g,t);
    }
    //Multiply by an appropriate Jabocian factor:
    value_type dt_by_dtcl = this->jacobian.dt_dtcl(tcl,g,beta_g);
    for(int i = 0;i < 15;i++)
    {
        //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
        beta_g[i] *= dt_by_dtcl;
    }
    beta_g[15] = this->_1p0;//This is tcl itself. Obviously dtcl/dtcl = 1.
    beta_g[16] = dt_by_dtcl;//This tracks t(tcl).
    if(this->nType == 4)
    {
        for(int i = 17;i < 23;i++)
        {
            //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
            beta_g[i] *= dt_by_dtcl;
        }
    }
}
template<class value_type>
void betaSMdS_1loop<value_type>::betap_tcl
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type tcl)
{
    //Compute original beta functions, in t co-ordinates:
    //1st t derivative of g:
    std::vector<value_type> gp (int(g.size()),this->_0p0 );
    //2nd t derivative of g:
    std::vector<value_type> gpp (int(g.size()),this->_0p0 );

    value_type t = this->eval.t(g);
    if(this->nType == 4)
    {
        this->beta_grav(g,gp,t);
        this->betap_grav(g,gpp,gp,t);
    }
    else
    {
        this->beta_extended(g,gp,t);
        this->betap_extended(g,gpp,gp,t);
    }

    //2nd derivatives of g stored in beta_g
    //Multiply by an appropriate Jabocian factor:
    //Note on efficiency - there is some repetition of work in the two
    //functions below. If we wanted to improve this, we could try to restructure
    //the way they do the computation to share data between them, reducing the
    //workload slightly.
    value_type dt_by_dtcl = this->jacobian.dt_dtcl(tcl,g,gp);
    value_type d2t_by_dtcl2 = this->jacobian.d2t_dtcl2(tcl,g,gp,gpp,dt_by_dtcl);
    /*
    std::cout << "\ndt_by_dtcl = " << dt_by_dtcl;
    std::cout << "\nd2t_by_dtcl2 = " << d2t_by_dtcl2;
    std::cout << "\ng = (";
    for(int i = 0;i < 15;i++){ if(i != 0) {std::cout << ",";} std::cout << g[i];}
    std::cout << ")";
    std::cout << "\ngp = (";
    for(int i = 0;i < 15;i++){ if(i != 0) {std::cout << ",";} std::cout << gp[i];}
    std::cout << ")";
    std::cout << "\ngpp = (";
    for(int i = 0;i < 15;i++){ if(i != 0) {std::cout << ",";} std::cout << gpp[i];}
    std::cout << ")";
    */
    //Compute beta derivatives:
    for(int i = 0;i < 15;i++)
    {
        //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
        beta_g[i] = dt_by_dtcl*dt_by_dtcl*gpp[i] + d2t_by_dtcl2*gp[i];
    }
    beta_g[15] = this->_0p0;//This is tcl itself. Obviously dtcl/dtcl = 1
        //and higher derivatives vanish.
    beta_g[16] = d2t_by_dtcl2;//This tracks t(tcl).
    if(this->nType == 4)
    {
        for(int i = 17;i < 23;i++)
        {
            //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
            beta_g[i] = dt_by_dtcl*dt_by_dtcl*gpp[i] + d2t_by_dtcl2*gp[i];
        }
    }
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapp_tcl
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type tcl)
{
    //Compute original beta functions, in t co-ordinates:
    //1st t derivative of g:
    std::vector<value_type> gp (int(g.size()),this->_0p0 );
    //2nd t derivative of g:
    std::vector<value_type> gpp (int(g.size()),this->_0p0 );
    //3rd t derivative of g:
    std::vector<value_type> gppp (int(g.size()),this->_0p0 );

    value_type t = this->eval.t(g);
    if(this->nType == 4)
    {
        this->beta_grav(g,gp,t);
        this->betap_grav(g,gpp,gp,t);
        this->betapp_grav(g,gppp,gp,gpp,t);
    }
    else
    {
        this->beta_extended(g,gp,t);
        this->betap_extended(g,gpp,gp,t);
        this->betapp_extended(g,gppp,gp,gpp,t);
    }


    //Multiply by an appropriate Jabocian factor:
    //Note on efficiency - there is some repetition of work in the two
    //functions below. If we wanted to improve this, we could try to restructure
    //the way they do the computation to share data between them, reducing the
    //workload slightly.
    value_type dt_by_dtcl = this->jacobian.dt_dtcl(tcl,g,gp);
    value_type d2t_by_dtcl2 = this->jacobian.d2t_dtcl2(tcl,g,gp,gpp,dt_by_dtcl);
    value_type d3t_by_dtcl3 = this->jacobian.d3t_dtcl3(tcl,g,gp,gpp,gppp,
                                                       dt_by_dtcl,d2t_by_dtcl2);
    //Compute beta derivatives:
    for(int i = 0;i < 15;i++)
    {
        //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
        beta_g[i] = this->_3p0*dt_by_dtcl*d2t_by_dtcl2*gpp[i]
                    + dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*gppp[i]
                    + d3t_by_dtcl3*gp[i];
    }
    beta_g[15] = this->_0p0;//This is tcl itself. Obviously dtcl/dtcl = 1
        //and higher derivatives vanish.
    beta_g[16] = d3t_by_dtcl3;//This tracks t(tcl).
    if(this->nType == 4)
    {
        for(int i = 17;i < 23;i++)
        {
            //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
            beta_g[i] = this->_3p0*dt_by_dtcl*d2t_by_dtcl2*gpp[i]
                        + dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*gppp[i]
                        + d3t_by_dtcl3*gp[i];
        }
    }
}
template<class value_type>
void betaSMdS_1loop<value_type>::betappp_tcl
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type tcl)
{
    //Compute original beta functions, in t co-ordinates:
    //1st t derivative of g:
    std::vector<value_type> gp (int(g.size()),this->_0p0 );
    //2nd t derivative of g:
    std::vector<value_type> gpp (int(g.size()),this->_0p0 );
    //3rd t derivative of g:
    std::vector<value_type> gppp (int(g.size()),this->_0p0 );
    //4th t derivative of g:
    std::vector<value_type> gpppp (int(g.size()),this->_0p0 );

    value_type t = this->eval.t(g);
    if(this->nType == 4)
    {
        this->beta_grav(g,gp,t);
        this->betap_grav(g,gpp,gp,t);
        this->betapp_grav(g,gppp,gp,gpp,t);
        this->betappp_grav(g,gpppp,gp,gpp,gppp,t);
    }
    else
    {
        this->beta_extended(g,gp,t);
        this->betap_extended(g,gpp,gp,t);
        this->betapp_extended(g,gppp,gp,gpp,t);
        this->betappp_extended(g,gpppp,gp,gpp,gppp,t);
    }



    //Multiply by an appropriate Jabocian factor:
    //Note on efficiency - there is some repetition of work in the two
    //functions below. If we wanted to improve this, we could try to restructure
    //the way they do the computation to share data between them, reducing the
    //workload slightly.
    value_type dt_by_dtcl = this->jacobian.dt_dtcl(tcl,g,gp);
    value_type d2t_by_dtcl2 = this->jacobian.d2t_dtcl2(tcl,g,gp,gpp,dt_by_dtcl);
    value_type d3t_by_dtcl3 = this->jacobian.d3t_dtcl3(tcl,g,gp,gpp,gppp,
                                                       dt_by_dtcl,d2t_by_dtcl2);
    value_type d4t_by_dtcl4 = this->jacobian.d4t_dtcl4(tcl,g,gp,gpp,gppp,gpppp,
                                                       dt_by_dtcl,d2t_by_dtcl2,
                                                       d3t_by_dtcl3);

    //Compute beta derivatives:
    for(int i = 0;i < 15;i++)
    {
        //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
        beta_g[i] = this->_3p0*d2t_by_dtcl2*d2t_by_dtcl2*gpp[i]
                    + this->_6p0*dt_by_dtcl*dt_by_dtcl*d2t_by_dtcl2*gppp[i]
                    + this->_4p0*dt_by_dtcl*d3t_by_dtcl3*gpp[i]
                    + dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*gpppp[i]
                    + d4t_by_dtcl4*gp[i];
    }
    beta_g[15] = this->_0p0;//This is tcl itself. Obviously dtcl/dtcl = 1
        //and higher derivatives vanish.
    beta_g[16] = d4t_by_dtcl4;//This tracks t(tcl).
    if(this->nType == 4)
    {
        for(int i = 17;i < 23;i++)
        {
            //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
            beta_g[i] = this->_3p0*d2t_by_dtcl2*d2t_by_dtcl2*gpp[i]
                        + this->_6p0*dt_by_dtcl*dt_by_dtcl*d2t_by_dtcl2*gppp[i]
                        + this->_4p0*dt_by_dtcl*d3t_by_dtcl3*gpp[i]
                        + dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*gpppp[i]
                        + d4t_by_dtcl4*gp[i];
        }
    }
}
template<class value_type>
void betaSMdS_1loop<value_type>::betapppp_tcl
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type tcl)
{
    //Compute original beta functions, in t co-ordinates:
    //1st t derivative of g:
    std::vector<value_type> gp (int(g.size()),this->_0p0 );
    //2nd t derivative of g:
    std::vector<value_type> gpp (int(g.size()),this->_0p0 );
    //3rd t derivative of g:
    std::vector<value_type> gppp (int(g.size()),this->_0p0 );
    //4th t derivative of g:
    std::vector<value_type> gpppp (int(g.size()),this->_0p0 );
    //4th t derivative of g:
    std::vector<value_type> gppppp (int(g.size()),this->_0p0 );

    value_type t = this->eval.t(g);
    if(this->nType == 4)
    {
        this->beta_grav(g,gp,t);
        this->betap_grav(g,gpp,gp,t);
        this->betapp_grav(g,gppp,gp,gpp,t);
        this->betappp_grav(g,gpppp,gp,gpp,gppp,t);
        this->betapppp_grav(g,gppppp,gp,gpp,gppp,gpppp,t);
    }
    else
    {
        this->beta_extended(g,gp,t);
        this->betap_extended(g,gpp,gp,t);
        this->betapp_extended(g,gppp,gp,gpp,t);
        this->betappp_extended(g,gpppp,gp,gpp,gppp,t);
        this->betapppp_extended(g,gppppp,gp,gpp,gppp,gpppp,t);
    }


    //Multiply by an appropriate Jabocian factor:
    //Note on efficiency - there is some repetition of work in the two
    //functions below. If we wanted to improve this, we could try to restructure
    //the way they do the computation to share data between them, reducing the
    //workload slightly.
    value_type dt_by_dtcl = this->jacobian.dt_dtcl(tcl,g,gp);
    value_type d2t_by_dtcl2 = this->jacobian.d2t_dtcl2(tcl,g,gp,gpp,dt_by_dtcl);
    value_type d3t_by_dtcl3 = this->jacobian.d3t_dtcl3(tcl,g,gp,gpp,gppp,
                                                       dt_by_dtcl,d2t_by_dtcl2);
    value_type d4t_by_dtcl4 = this->jacobian.d4t_dtcl4(tcl,g,gp,gpp,gppp,gpppp,
                                                       dt_by_dtcl,d2t_by_dtcl2,
                                                       d3t_by_dtcl3);
    value_type d5t_by_dtcl5 = this->jacobian.d5t_dtcl5(tcl,g,gp,gpp,gppp,gpppp,
                                                          gppppp,dt_by_dtcl,
                                                          d2t_by_dtcl2,
                                                          d3t_by_dtcl3,
                                                          d4t_by_dtcl4);

    //Compute beta derivatives:
    for(int i = 0;i < 15;i++)
    {
        //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
        beta_g[i] = this->_15p0*dt_by_dtcl*d2t_by_dtcl2*d2t_by_dtcl2*gppp[i]
                    + this->_10p0*gpp[i]*d2t_by_dtcl2*d3t_by_dtcl3
                    + this->_10p0*dt_by_dtcl*dt_by_dtcl*gppp[i]*d3t_by_dtcl3
                    +this->_10p0*dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*d2t_by_dtcl2
                    *gpppp[i]
                    +this->_5p0*gpp[i]*d4t_by_dtcl4
                    + dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*dt_by_dtcl
                    *gppppp[i]
                    + gp[i]*d5t_by_dtcl5;
    }
    beta_g[15] = this->_0p0;//This is tcl itself. Obviously dtcl/dtcl = 1
        //and higher derivatives vanish.
    beta_g[16] = d5t_by_dtcl5;//This tracks t(tcl).
    if(this->nType == 4)
    {
        for(int i = 17;i < 23;i++)
        {
            //Use the fact that dg/dtcl = (dt/dtcl)*(dg/dt):
            beta_g[i] = this->_15p0*dt_by_dtcl*d2t_by_dtcl2*d2t_by_dtcl2*gppp[i]
                        + this->_10p0*gpp[i]*d2t_by_dtcl2*d3t_by_dtcl3
                        + this->_10p0*dt_by_dtcl*dt_by_dtcl*gppp[i]*d3t_by_dtcl3
                        +this->_10p0*dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*d2t_by_dtcl2
                        *gpppp[i]
                        +this->_5p0*gpp[i]*d4t_by_dtcl4
                        + dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*dt_by_dtcl*dt_by_dtcl
                        *gppppp[i]
                        + gp[i]*d5t_by_dtcl5;
        }
    }
}



template<class value_type>
void betaSMdS_1loop<value_type>::operator()
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type& t,int nDer,int NTYPE)
{
    switch(nDer)
    {
    case 0:
        this->beta_choice(g,beta_g,t,NTYPE,
                          &betaSMdS_1loop<value_type>::beta_basic,
                          &betaSMdS_1loop<value_type>::beta_extended,
                          &betaSMdS_1loop<value_type>::beta_tcl,
                          &betaSMdS_1loop<value_type>::beta_grav);
        break;
    case 1:
        this->beta_choice(g,beta_g,t,NTYPE,
                          &betaSMdS_1loop<value_type>::betap_basic,
                          &betaSMdS_1loop<value_type>::betap_extended,
                          &betaSMdS_1loop<value_type>::betap_tcl,
                          &betaSMdS_1loop<value_type>::betap_grav);
        break;
    case 2:
        this->beta_choice(g,beta_g,t,NTYPE,
                          &betaSMdS_1loop<value_type>::betapp_basic,
                          &betaSMdS_1loop<value_type>::betapp_extended,
                          &betaSMdS_1loop<value_type>::betapp_tcl,
                          &betaSMdS_1loop<value_type>::betapp_grav);
        break;
    case 3:
        this->beta_choice(g,beta_g,t,NTYPE,
                          &betaSMdS_1loop<value_type>::betappp_basic,
                          &betaSMdS_1loop<value_type>::betappp_extended,
                          &betaSMdS_1loop<value_type>::betappp_tcl,
                          &betaSMdS_1loop<value_type>::betappp_grav);
        break;
    case 4:
        this->beta_choice(g,beta_g,t,NTYPE,
                          &betaSMdS_1loop<value_type>::betapppp_basic,
                          &betaSMdS_1loop<value_type>::betapppp_extended,
                          &betaSMdS_1loop<value_type>::betapppp_tcl,
                          &betaSMdS_1loop<value_type>::betapppp_grav);
        break;
    default:
        throw("Error - invalid beta function derivative.");
    }
}
template<class value_type>
void betaSMdS_1loop<value_type>::operator()
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type& t,int nDer)
{
    this->operator()(g,beta_g,t,nDer,this->nType);
}
template<class value_type>
void betaSMdS_1loop<value_type>::betap
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type& t,int N)
{
    this->operator()(g,beta_g,t,N,this->nType);
}
template<class value_type>
void betaSMdS_1loop<value_type>::beta
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type& t)
{
    this->operator()(g,beta_g,t,0,this->nType);
}
template<class value_type>
void betaSMdS_1loop<value_type>::operator()
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type& t)
{
    this->operator()(g,beta_g,t,this->nDerivative,this->nType);
}
template<class value_type>
void betaSMdS_1loop<value_type>::beta_choice
(const std::vector<value_type>& g,std::vector<value_type>& beta_g,
 const value_type& t,int NTYPE,
 typename betaSMdS_1loop<value_type>::beta_pointer basic,
 typename betaSMdS_1loop<value_type>::beta_pointer extended,
 typename betaSMdS_1loop<value_type>::beta_pointer tcl,
 typename betaSMdS_1loop<value_type>::beta_pointer grav)
{
    switch(NTYPE)
    {
    case 0:
        //this->nType = 0;
        (this->*basic)(g,beta_g,t);
        break;
    case 1:
        //this->nType = 1;
        (this->*extended)(g,beta_g,t);
        break;
    case 2:
        //this->nType = 2;
        (this->*tcl)(g,beta_g,t);
        break;
    case 3:
        //this->nType = 3;
        (this->*grav)(g,beta_g,t);
        break;
    case 4:
        //this->nType = 4;
        (this->*tcl)(g,beta_g,t);
        break;
    default:
        throw("Invalid beta function type requested.");
    }
}

#endif // USING_SM_POTENTIAL
//Reference list of couplings:
    /*
    //Second argument, [i], refers to the element of the computed grid. First
    //argument determined the coupling.
    ydata[0][i] - lambda (Higgs self coupling)
    ydata[1][i] - g1^2 (U(1) coupling)
    ydata[2][i] - g2^2 (SU(2) coupling)
    ydata[3][i] - g3^2 (SU(3) coupling)
    ydata[4][i] - yt^2 (Top quark Yukawa coupling)
    ydata[5][i] - m^2 (Higgs tachyonic mass term)
    ydata[6][i] - xi (Higgs-curvature coupling)
    ydata[7][i] - Z (Higgs field renormalisation) - technically the RHS of this
                  is a gamma function, but we solve it with the beta functions)
    //Additionally, we have to solve for the non-beta functions:
    ydata[8][i] = t_cl = log(phicl^2/mt^2)
    ydata[9][i] = t = log(mu^2/mt^2)
    */
    /*
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::operator()
(value_type y)
{
    //y = phi_tilde (Einstein frame field)
    //Re-scale depending on units used:
    value_type x = y*this->yscale;
    //Choose scale to evaluate couplings at:
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    std::vector<value_type> gdata (10,_0p0);
    for(int i = 0;i < 10;i++)
    {
        gdata[i] = this->Spline(y,this->kdata,this->ydata,n);
    }
    //Compute phi_c(phi_tilde)
    value_type phicl = this->Spline(y,this->phiclp_data,this->phicl_data,n);

    //Intermediate variables:
    value_type h2 = this->h*this->h;
    value_type Z = gdata[7];//Field renormalisation
    value_type phi = Z*phicl;//renormalised field
    value_type xi = gdata[6];//non-minimal coupling.
    value_type xifactor = (this->_1p0 - h2*xi*phi*phi);//Einstein frame
        //potential denominator.

    //Compute potential:
    value_type A = this->V(phicl,gdata,false);
    value_type B = xifactor*xifactor;
    return A/B;
}
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::operator()
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel,i++)
    {
        arrayOut[i] = this->operator()(y[i]);
    }
}
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::d(value_type y)
{
    //y = phi_tilde (Einstein frame field)
    //Re-scale depending on units used:
    value_type x = y*this->yscale;
    //Choose scale to evaluate couplings at:
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    std::vector<value_type> gdata (10,_0p0);
    std::vector<value_type> gpdata_phixi(10,_0p0);//derivatives wrt to phi_tilde
    std::vector<value_type> gpdata_phicl(10,_0p0);//wrt to phicl
    for(int i = 0;i < 10;i++)
    {
        gdata[i] = this->Spline(y,this->kdata,this->ydata,n);
        gpdata_phixi[i] = this->Spline(y,this->kpdata_phixi,
                                       this->ypdata_phixi,n);
        gpdata_phicl[i] = this->Spline(y,this->kpdata_phicl,
                                       this->ypdata_phicl,n);
    }
    //Compute phi(phi_tilde) and derivatives (bear in mind that phi here
    // is not the classical field, phicl: phi = Z(t)*phicl).
    //value_type phicl = this->phicl_data[n];
    //value_type phiclp = this->phiclp_data[n];
    value_type phicl = this->Spline(y,this->dphicl_dphixi_data,
                                    this->phicl_data,n);
    value_type phiclp = this->Spline(y,this->d2phicl_dphixi2_data,
                                     this->dphicl_dphixi_data,n);

    //Temporary variables:
    value_type h2 = this->h*this->h;
    std::vector<value_type> Mi2;
    std::vector<value_type> logMi2;
    std::vector<value_type> Mi2p;
    std::vector<value_type> logMi2p;
    this->computeMi2p(phicl,gdata,gpdata_phicl,Mi2,Mi2p,logMi2,logMi2p);
    //Derivatives of t_xi wrt phi_tilde:
    value_type dt_xi_dphi_tilde = this->_2p0/x;
    //Derivatives of xi wrt to phi_xi:
    value_type xi = gdata[6];
    value_type xip = gpdata_phixi[6];
    value_type Z = gdata[7];//Field renormalisation
    value_type Zp = gpdata_phixi[7];
    value_type phi = Z*phicl;//renormalised field
    //Derivative of phi wrt phi_tilde:
    value_type phip = Zp*phicl + Z*phiclp;
    value_type xifactor = (this->_1p0 - h2*xi*phi*phi);


    //Parts of rational function:
    value_type A = this->V(phicl,gdata,logMi2,Mi2);
    value_type Ap = this->Vp(phicl,gdata,gpdata_phicl,logMi2,Mi2,logMi2p,
                             Mi2p)*phiclp;
    value_type B = xifactor*xifactor;
    value_type Bp = this->_2p0*xifactor*( -xip*h2*phi*phi - h2*xi
                                         *this->_2p0*phi*phip );
    return Ap/B - A*Bp/(B*B);
}
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::d
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel,i++)
    {
        arrayOut[i] = this->d(y[i]);
    }
}
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::d2(value_type y)
{
    //y = phi_tilde (Einstein frame field)
    //Re-scale depending on units used:
    value_type x = y*this->yscale;
    //Choose scale to evaluate couplings at:
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    std::vector<value_type> gdata (10,_0p0);
    std::vector<value_type> gpdata_phicl(10,_0p0);
    std::vector<value_type> gpdata_phixi(10,_0p0);
    std::vector<value_type> gppdata_phicl(10,_0p0);
    std::vector<value_type> gppdata_phixi(10,_0p0);
    for(int i = 0;i < 10;i++)
    {
        gdata[i] = this->Spline(y,this->kdata,this->ydata,n);
        gpdata_phixi[i] = this->Spline(y,this->kpdata_phixi,
                                       this->ypdata_phixi,n);
        gpdata_phicl[i] = this->Spline(y,this->kpdata_phicl,
                                       this->ypdata_phicl,n);
        gppdata_phixi[i] = this->Spline(y,this->kppdata_phixi,
                                       this->yppdata_phixi,n);
        gppdata_phicl[i] = this->Spline(y,this->kppdata_phicl,
                                       this->yppdata_phicl,n);
    }
    //Compute phi_c(phi_tilde)
    value_type phicl = this->Spline(y,this->dphicl_dphixi_data,
                                    this->phicl_data,n);
    value_type phiclp = this->Spline(y,this->d2phicl_dphixi2_data,
                                     this->dphicl_dphixi_data,n);
    value_type phiclpp = this->Spline(y,this->d3phicl_dphixi3_data,
                                     this->d2phicl_dphixi2_data,n);
    value_type phi = Z*phicl;//renormalised field
    //Derivatives of phi wrt phi_tilde:
    value_type Z = gdata[7];//Field renormalisation
    value_type Zp = gpdata_phixi[7];
    value_type Zpp = gppdata_phixi[7];
    value_type phip = Zp*phicl + Z*phiclp;
    value_type phipp = Zpp*phicl + this->_2p0*Zp*phiclp + Z*phiclpp;

    //Temporary variables:
    value_type h2 = this->h*this->h;
    value_type h4 = h2*h2;
    value_type xifactor = (this->_1p0 - h2*gdata[6]*phi*phi);

    std::vector<value_type> Mi2;
    std::vector<value_type> logMi2;
    std::vector<value_type> Mi2p;
    std::vector<value_type> logMi2p;
    std::vector<value_type> Mi2pp;
    std::vector<value_type> logMi2pp;
    this->computeMi2pp(phicl,gdata,gpdata,gppdata,
                       Mi2,Mi2p,Mi2pp,logMi2,logMi2p,logMi2pp);
    //Derivatives of t_xi wrt phi_tilde:
    value_type dt_xi_dphi_tilde = this->_2p0/x;
    value_type d2t_xi_dphi_tilde2 = -this->_2p0/(x*x);
    //Derivatives of xi wrt phi_tilde:
    value_type xi = gdata[6];
    value_type xip = gpdata[6]*dt_xi_dphi_tilde;
    value_type xipp = gppdata[6]*dt_xi_dphi_tilde*dt_xi_dphi_tilde
                        + gpdata[6]*d2t_xi_dphi_tilde2;
    value_type temp1 = phi*xip + this->_2p0*xi*phip;

    //Jordan frame potential and derivatives:
    value_type v = this->V(phicl,gdata,logMi2,Mi2);
    value_type vp = this->Vp(phicl,gdata,gpdata,logMi2,Mi2,logMi2p,
                             Mi2p);
    value_type vpp = this->Vpp(phicl,gdata,gpdata,gppdata,
                               logMi2,Mi2,logMi2p,Mi2p,logMi2pp,Mi2pp);

    //Parts of rational function:
    //Parts of rational function:
    value_type A = v;
    value_type Ap = vp*phiclp;
    value_type B = xifactor*xifactor;
    value_type Bp = this->_2p0*xifactor*( -xip*h2*phi*phi
                                         - h2*xi*this->_2p0*phi*phip );
    value_type App = vpp*phiclp*phiclp + vp*phiclpp;
    value_type Bpp = this->_2p0*(h4*phi*phi*(temp1*temp1) + h2*xifactor
                                 +(this->_2p0*xi*phip*phip
                                   + phi*phi*xipp
                                   + this->_2p0*phi
                                   *(this->_2p0*xip*phip + xi*phipp) ) );

    return -this->_2p0*Ap*Bp/(B*B) + this->_2p0*A*(Bp*Bp)/(B*B*B)
            + App/B - A*Bpp/(B*B);
}
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::d2
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel,i++)
    {
        arrayOut[i] = this->d2(y[i]);
    }
}
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::d3(value_type y)
{
    //y = phi_tilde (Einstein frame field)
    //Re-scale depending on units used:
    value_type x = y*this->yscale;
    //Choose scale to evaluate couplings at:
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    std::vector<value_type> gdata (10,_0p0);
    std::vector<value_type> gpdata_phicl(10,_0p0);
    std::vector<value_type> gpdata_phixi(10,_0p0);
    std::vector<value_type> gppdata_phicl(10,_0p0);
    std::vector<value_type> gppdata_phixi(10,_0p0);
    std::vector<value_type> gpppdata_phicl(10,_0p0);
    std::vector<value_type> gpppdata_phixi(10,_0p0);
    for(int i = 0;i < 10;i++)
    {
        gdata[i] = this->Spline(y,this->kdata,this->ydata,n);
        gpdata_phixi[i] = this->Spline(y,this->kpdata_phixi,
                                       this->ypdata_phixi,n);
        gpdata_phicl[i] = this->Spline(y,this->kpdata_phicl,
                                       this->ypdata_phicl,n);
        gppdata_phixi[i] = this->Spline(y,this->kppdata_phixi,
                                       this->yppdata_phixi,n);
        gppdata_phicl[i] = this->Spline(y,this->kppdata_phicl,
                                       this->yppdata_phicl,n);
        gpppdata_phixi[i] = this->Spline(y,this->kpppdata_phixi,
                                       this->ypppdata_phixi,n);
        gpppdata_phicl[i] = this->Spline(y,this->kpppdata_phicl,
                                       this->ypppdata_phicl,n);
    }
    //Compute phi_c(phi_tilde)
    value_type phicl = this->Spline(y,this->dphicl_dphixi_data,
                                    this->phicl_data,n);
    value_type phiclp = this->Spline(y,this->d2phicl_dphixi2_data,
                                     this->dphicl_dphixi_data,n);
    value_type phiclpp = this->Spline(y,this->d3phicl_dphixi3_data,
                                     this->d2phicl_dphixi2_data,n);
    value_type phiclppp = this->Spline(y,this->d4phicl_dphixi4_data,
                                     this->d3phicl_dphixi3_data,n);

    value_type Z = gdata[7];//Field renormalisation
    value_type Zp = gpdata_phixi[7];
    value_type Zpp = gppdata_phixi[7];
    value_type Zppp = gpppdata_phixi[7];
    value_type phi = Z*phicl;//renormalised field
    //Derivatives of phi wrt phi_tilde:
    value_type phip = Zp*phicl + Z*phiclp;
    value_type phipp = Zpp*phicl + this->_2p0*Zp*phiclp + Z*phiclpp;
    value_type phippp = Zppp*phicl + this->_3p0*Zpp*phiclp
                        + this->_3p0*Zp*phiclpp + Z*phiclppp;

    //Temporary variables:
    value_type h2 = this->h*this->h;
    value_type h4 = h2*h2;
    value_type xifactor = (this->_1p0 - h2*gdata[6]*phi*phi);

    std::vector<value_type> Mi2;
    std::vector<value_type> logMi2;
    std::vector<value_type> Mi2p;
    std::vector<value_type> logMi2p;
    std::vector<value_type> Mi2pp;
    std::vector<value_type> logMi2pp;
    std::vector<value_type> Mi2ppp;
    std::vector<value_type> logMi2ppp;
    this->computeMi2ppp(phicl,gdata,gpdata,gppdata,
                       gpppdata,
                       Mi2,Mi2p,Mi2pp,Mi2ppp,
                       logMi2,logMi2p,logMi2pp,logMi2ppp);
    //Derivatives of t_xi wrt phi_tilde:
    value_type dt_xi_dphi_tilde = this->_2p0/x;
    value_type d2t_xi_dphi_tilde2 = -this->_2p0/(x*x);
    value_type d3t_xi_dphi_tilde3 = this->_4p0/(x*x*x);
    //Derivatives of xi wrt phi_tilde:
    value_type xi = gdata[6];
    value_type xip = gpdata[6]*dt_xi_dphi_tilde;
    value_type xipp = gppdata[6]*dt_xi_dphi_tilde*dt_xi_dphi_tilde
                        + gpdata[6]*d2t_xi_dphi_tilde2;
    value_type xippp = this->_3p0*dt_xi_dphi_tilde*gppdata[6]
                        *d2t_xi_dphi_tilde2 + dt_xi_dphi_tilde*dt_xi_dphi_tilde
                        *dt_xi_dphi_tilde*gpppdata[6]
                        + gpdata[6]*d3t_xi_dphi_tilde3;
    value_type temp1 = phi*xip + this->_2p0*xi*phip;

    //Jordan frame potential and derivatives:
    value_type v = this->V(phicl,gdata,logMi2,Mi2);
    value_type vp = this->Vp(phicl,gdata,gpdata,logMi2,Mi2,logMi2p,
                             Mi2p);
    value_type vpp = this->Vpp(phicl,gdata,gpdata,gppdata,
                               logMi2,Mi2,logMi2p,Mi2p,logMi2pp,Mi2pp);
    value_type vppp = this->Vppp(phicl,gdata,gpdata,gppdata,
                                 gpppdata,
                                 logMi2,Mi2,logMi2p,Mi2p,logMi2pp,Mi2pp,
                                 logMi2ppp,Mi2ppp);

    //Parts of rational function:
    //Parts of rational function:
    value_type A = v;
    value_type Ap = vp*phiclp;
    value_type B = xifactor*xifactor;
    value_type Bp = this->_2p0*xifactor*( -xip*h2*phi*phi
                                         - h2*xi*this->_2p0*phi*phip );
    value_type App = vpp*phiclp*phiclp + vp*phiclpp;
    value_type Bpp = this->_2p0*(h4*phi*phi*(temp1*temp1) + h2*xifactor
                                 +(this->_2p0*xi*phip*phip
                                   + phi*phi*xipp
                                   + this->_2p0*phi
                                   *(this->_2p0*xip*phip + xi*phipp) ) );
    value_type Appp = this->_3p0*phiclp*vpp*phiclpp + phiclp*phiclp*phiclp*vppp
                        +vp*phiclppp;
    value_type Bppp = this->_6p0*h4*phi*temp1
                      *(this->_2p0*xi*phip*phip + phi*phi*xipp +
                        this->_2p0*phi*(this->_2p0*xip*phip + xi*phipp) )
                      +this->_2p0*h2*xifactor
                      *(this->_6p0*xi*phip*phipp + this->_6p0*xip
                        *(phip*phip + phi*phipp) + phi*phi*xippp
                        + this->_2p0*phi
                        *(this->_3p0*phip*xipp + xi*phippp) );

    return this->_6p0*Ap*Bp*Bp/(B*B) -this->_6p0*A*Bp*Bp*Bp/(B*B*B*B)
            -this->_3p0*Bp*App/(B*B) - this->_3p0*Ap*Bpp/(B*B)
            + this->_6p0*A*Bp*Bpp/(B*B*B) + Appp/B - A*Bppp/(B*B);
}
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::d3
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel,i++)
    {
        arrayOut[i] = this->d3(y[i]);
    }
}
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::d4(value_type y)
{
    //y = phi_tilde (Einstein frame field)
    //Re-scale depending on units used:
    value_type x = y*this->yscale;
    //Choose scale to evaluate couplings at:
    int n = this->ChooseSpline(x, this->xdata, this->Ndata - 1);
    std::vector<value_type> gdata (10,_0p0);
    std::vector<value_type> gpdata_tcl(10,_0p0);
    std::vector<value_type> gpdata_txi(10,_0p0);
    std::vector<value_type> gppdata_tcl(10,_0p0);
    std::vector<value_type> gppdata_txi(10,_0p0);
    std::vector<value_type> gpppdata_tcl(10,_0p0);
    std::vector<value_type> gpppdata_txi(10,_0p0);
    std::vector<value_type> gppppdata_tcl(10,_0p0);
    std::vector<value_type> gppppdata_txi(10,_0p0);
    for(int i = 0;i < 10;i++)
    {
        gdata[i] = this->Spline(y,this->kdata,this->ydata,n);
        gpdata_txi[i] = this->Spline(y,this->kpdata_txi,
                                       this->ypdata_txi,n);
        gpdata_tcl[i] = this->Spline(y,this->kpdata_tcl,
                                       this->ypdata_tcl,n);
        gppdata_txi[i] = this->Spline(y,this->kppdata_txi,
                                       this->yppdata_txi,n);
        gppdata_tcl[i] = this->Spline(y,this->kppdata_tcl,
                                       this->yppdata_tcl,n);
        gpppdata_txi[i] = this->Spline(y,this->kpppdata_txi,
                                       this->ypppdata_txi,n);
        gpppdata_tcl[i] = this->Spline(y,this->kpppdata_tcl,
                                       this->ypppdata_tcl,n);
        gppppdata_txi[i] = this->Spline(y,this->kppppdata_txi,
                                       this->yppppdata_txi,n);
        gppppdata_tcl[i] = this->Spline(y,this->kppppdata_tcl,
                                       this->yppppdata_tcl,n);
    }
    //Compute phi_c(phi_tilde)
    value_type tcl = this->Spline(y,this->dtcl_dtxi_data,
                                    this->tcl_data,n);
    value_type phicl = this->mt*exp(tcl/this->_2p0);
    value_type phi_tilde = x;//Einstein frame field.
    //Derivatives in terms of t_cl and t_xi. These are
    //in log space, so easier to construct spline approximations:
    value_type tclp = this->Spline(y,this->d2tcl_dtxi2_data,
                                     this->dtcl_dtxi_data,n);
    value_type tclpp = this->Spline(y,this->d3tcl_dtxi3_data,
                                     this->d2tcl_dtxi2_data,n);
    value_type tclppp = this->Spline(y,this->d4tcl_dtxi4_data,
                                     this->d3tcl_dtxi3_data,n);
    value_type tclpppp = this->Spline(y,this->d5tcl_dtxi5_data,
                                     this->d4tcl_dtxi4_data,n);
    //Derivatives of phi_cl (Jordan frame) wrt phi_tilde
    //(Einstein frame) field:
    value_type phiclp = (phicl/phi_tilde)*tclp;
    value_type phiclpp = (phicl/(phi_tilde*phi_tilde))*(-tclp + tclp*tclp
                                            + this->_2p0*tclpp);
    value_type phiclppp =(phicl/(phi_tilde*phi_tilde*phi_tilde))
                        *(-this->_3p0*tclp*tclp +
                                                tclp*tclp*tclp
                                                - this->_6p0*tclpp +
                                                tclp*(this->_2p0
                                                      + this->_6p0*tclpp)
                                                + this->_4p0*tclppp);
    value_type phiclpppp = (phicl/(phi_tilde*phi_tilde*phi_tilde*phi_tilde))
                            *(-this->_6p0*tclp*tclp*tclp + tclp*tclp*tclp*tclp
                              +this->_22p0*tclpp + this->_12p0*tclpp*tclpp
                              +tclp*tclp*(this->_11p0 + this->_12p0*tclpp)
                              - this->_2p0*tclp*(this->_3p0 + this->_18p0*tclpp
                                                  - this->_8p0&tclppp)
                              + this->_8p0*(-this->_3p0*tclppp + tclpppp));

    //Field renormalisation, Z
    value_type Z = gdata[7];
    //t_xi derivatives of Z:
    value_type Zd = gpdata_txi[7];
    value_type Zdd = gppdata_txi[7];
    value_type Zddd = gpppdata_txi[7];
    value_type Zdddd = gppppdata_txi[7];
    //phi_tilde derivatives of Z:
    value_type Zp = this->_2p0*Zd/phi_tilde;
    value_type Zpp = (this->_1p0/(phi_tilde*phi_tilde))
                    *(-this->_2p0*Zd + this->_4p0*Zdd);
    value_type Zppp = (this->_4p0*Zd - this->_12p0*Zdd + this->_8p0*Zddd)
                      /(phi_tilde*phi_tilde*phi_tilde);
    value_type Zpppp = (-this->_12p0*Zd + this->_44p0*Zdd - this->_48p0*Zddd
                        + this->_16p0*Zdddd)
                      /(phi_tilde*phi_tilde*phi_tilde*phi_tilde);
    value_type phi = Z*phicl;//renormalised field
    //Derivatives of phi wrt phi_tilde:
    value_type phip = Zp*phicl + Z*phiclp;
    value_type phipp = Zpp*phicl + this->_2p0*Zp*phiclp + Z*phiclpp;
    value_type phippp = Zppp*phicl + this->_3p0*Zpp*phiclp
                        + this->_3p0*Zp*phiclpp + Z*phiclppp;
    value_type phipppp = Zpppp*phicl + this->_4p0*Zppp*phiclp
                        + this->_6p0*Zpp*phiclpp + this->_4p0*Zp*phiclppp
                        + Z*phiclpppp;

    //Temporary variables:
    value_type h2 = this->h*this->h;
    value_type h4 = h2*h2;
    value_type xifactor = (this->_1p0 - h2*gdata[6]*phi*phi);

    std::vector<value_type> Mi2;
    std::vector<value_type> logMi2;
    std::vector<value_type> Mi2p;
    std::vector<value_type> logMi2p;
    std::vector<value_type> Mi2pp;
    std::vector<value_type> logMi2pp;
    std::vector<value_type> Mi2ppp;
    std::vector<value_type> logMi2ppp;
    std::vector<value_type> Mi2pppp;
    std::vector<value_type> logMi2pppp;
    this->computeMi2pppp(phicl,gdata,gpdata,gppdata,
                       gpppdata,gppppdata,
                       Mi2,Mi2p,Mi2pp,Mi2ppp,Mi2pppp,
                       logMi2,logMi2p,logMi2pp,logMi2ppp,logMi2pppp);
    //Derivatives of t_xi wrt phi_tilde:
    value_type dt_xi_dphi_tilde = this->_2p0/x;
    value_type d2t_xi_dphi_tilde2 = -this->_2p0/(x*x);
    value_type d3t_xi_dphi_tilde3 = this->_4p0/(x*x*x);
    value_type d4t_xi_dphi_tilde4 = -this->_12p0/(x*x*x*x);
    //Derivatives of xi wrt phi_tilde:
    value_type xi = gdata[6];
    value_type xip = gpdata[6]*dt_xi_dphi_tilde;
    value_type xipp = gppdata[6]*dt_xi_dphi_tilde*dt_xi_dphi_tilde
                        + gpdata[6]*d2t_xi_dphi_tilde2;
    value_type xippp = this->_3p0*dt_xi_dphi_tilde*gppdata[6]
                        *d2t_xi_dphi_tilde2 + dt_xi_dphi_tilde*dt_xi_dphi_tilde
                        *dt_xi_dphi_tilde*gpppdata[6]
                        + gpdata[6]*d3t_xi_dphi_tilde3;
    value_type xipppp = this->_3p0*gppdata[6]
                        *d2t_xi_dphi_tilde2*d2t_xi_dphi_tilde2
                        +this->_6p0*dt_xi_dphi_tilde*dt_xi_dphi_tilde
                        *d2t_xi_dphi_tilde2*gpppdata[6]
                        +this->_4p0*dt_xi_dphi_tilde*d3t_xi_dphi_tilde3
                        *gppdata[6]
                        + dt_xi_dphi_tilde*dt_xi_dphi_tilde*dt_xi_dphi_tilde
                        *dt_xi_dphi_tilde*gppppdata[6]
                        + gpdata[6]*d4t_xi_dphi_tilde4;

    //More temps:
    value_type temp1 = phi*xip + this->_2p0*xi*phip;
    value_type temp2 = this->_2p0*xi*phip*phip + phi*phi*xipp
                        +this->_2p0*phi*(this->_2p0*xip*phip + xi*phipp);

    //Jordan frame potential and derivatives:
    value_type v = this->V(phicl,gdata,logMi2,Mi2);
    value_type vp = this->Vp(phicl,gdata,gpdata,logMi2,Mi2,logMi2p,
                             Mi2p);
    value_type vpp = this->Vpp(phicl,gdata,gpdata,gppdata,
                               logMi2,Mi2,logMi2p,Mi2p,logMi2pp,Mi2pp);
    value_type vppp = this->Vppp(phicl,gdata,gpdata,gppdata,
                                 gpppdata,
                                 logMi2,Mi2,logMi2p,Mi2p,logMi2pp,Mi2pp,
                                 logMi2ppp,Mi2ppp);
    value_type vpppp = this->Vpppp(phicl,gdata,gpdata,gppdata,
                                   gpppdata,gppppdata,
                                   logMi2,Mi2,logMi2p,Mi2p,logMi2pp,Mi2pp,
                                   logMi2ppp,Mi2ppp,LogMi2pppp,Mi2pppp);

    //Parts of rational function:
    value_type A = v;
    value_type Ap = vp*phiclp;
    value_type B = xifactor*xifactor;
    value_type Bp = this->_2p0*xifactor*( -xip*h2*phi*phi
                                         - h2*xi*this->_2p0*phi*phip );
    value_type App = vpp*phiclp*phiclp + vp*phiclpp;
    value_type Bpp = this->_2p0*(h4*phi*phi*(temp1*temp1) + h2*xifactor
                                 +(this->_2p0*xi*phip*phip
                                   + phi*phi*xipp
                                   + this->_2p0*phi
                                   *(this->_2p0*xip*phip + xi*phipp) ) );
    value_type Appp = this->_3p0*phiclp*vpp*phiclpp + phiclp*phiclp*phiclp*vppp
                        +vp*phiclppp;
    value_type Bppp = this->_6p0*h4*phi*temp1
                      *(this->_2p0*xi*phip*phip + phi*phi*xipp +
                        this->_2p0*phi*(this->_2p0*xip*phip + xi*phipp) )
                      +this->_2p0*h2*xifactor
                      *(this->_6p0*xi*phip*phipp + this->_6p0*xip
                        *(phip*phip + phi*phipp) + phi*phi*xippp
                        + this->_2p0*phi
                        *(this->_3p0*phip*xipp + xi*phippp) );
    value_type Apppp = this->_6p0*phiclp*phiclp*phiclpp*vppp
                        +vpp*(this->_3p0*phiclpp*phiclpp
                              + this->_4p0*phiclp*phiclppp)
                        +phiclp*phiclp*phiclp*phiclp*vpppp + vp*phiclpppp;
    value_type Bpppp = this->_6p0*h4*temp2*temp2
                        + this->_8p0*h4*phi*(phi*xip + this->_2p0*xi*phip)
                        *(this->_6p0*xi*phip*phipp + this->_6p0*xip*
                          (phip*phip + phi*phipp)
                           + phi*phi*xippp + this->_2p0*phi
                           *(this->_3p0*phip*xipp + xi*phippp) )
                        + this->_2p0*h2*xifactor*
                        (this->_12p0*phip*phip*xipp
                         + this->_6p0*xi*phipp*phipp
                         + this->_8p0*phip*(this->_3p0*xip*phipp
                                              + phi*xippp + xi*phippp)
                         +phi*phi*xipppp
                         + this->_2p0*phi*(this->_6p0*xipp*phipp
                                             + this->_4p0*xip*phippp
                                             + xi*phipppp) );

    return -this->_24p0*Ap*Bp*Bp*Bp/(B*B*B*B)
           + this->_24p0*A*Bp*Bp*Bp*Bp/(B*B*B*B*B)
           + this->_12p0*Bp*Bp*App/(B*B*B)
           + this->_24p0*Ap*Bp*Bpp/(B*B*B)
           - this->_36p0*A*Bp*Bp*Bpp/(B*B*B*B)
           - this->_6p0*App*Bpp/(B*B)
           + this->_6p0*A*Bpp*Bpp/(B*B*B)
           - this->_4p0*Bp*Appp/(B*B)
           - this->_4p0*Ap*Bppp/(B*B)
           + this->_8p0*A*Bp*Bppp/(B*B*B)
           + Apppp/B - A*Bpppp/(B*B);
}
template<class value_type>
value_type SMHiggsPotential_dS_1loop_einstein<value_type>::d4
(value_type* y,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel,i++)
    {
        arrayOut[i] = this->d4(y[i]);
    }
}
//------------------------------------------------------------------------------
//Compute couplings and their derivatives for a given input scale



//------------------------------------------------------------------------------

/*
template< class value_type >
DLL_EXPORT void SMHiggsPotentialSplineDeSitter< value_type >
::SinglePointSpline(value_type& y, value_type* g , int m)
{
    int n = ChooseSpline(y, xdata, Ndata - 1);
    SinglePointSpline(y,g,m,n);
}
//------------------------------------------------------------------------------
//Same as the above, but assuming n is already known.
template< class value_type >
DLL_EXPORT void SMHiggsPotentialSplineDeSitter
< value_type >::SinglePointSpline(value_type& y, value_type* g, int m,int n)
{
	//Construct a cubic spline approximation of the running coupling:
	//NB - argument y MUST be in units of GeV.
	//int n = ChooseSpline(y, xdata, Ndata - 1);
	//const value_type _1p0 = value_type(1.0);
	value_type t = (log(y*y / (M*M)) - xdata[n]) / (xdata[n + 1] - xdata[n]);
	value_type a[8];
	value_type b[8];
	for(int i = 0;i < 8;i++)
    {
        value_type a[i] = kdata[i][m][n] * (xdata[n + 1] - xdata[n])
                          - (ydata[i][m][n +1] - ydata[i][m][n]);
        value_type b[i] = -kdata[i][m][n + 1] * (xdata[n + 1] - xdata[n])
                       + (ydata[i][m][n + 1] - ydata[i][m][n]);
        //Store result of interpolating each point:
        g[i] = (_1p0 - t)*ydata[i][m][n] + t*ydata[i][m][n + 1]
            + t*(_1p0 - t)*(a[i]*(_1p0 - t) + b[i]*t);
    }
}
//------------------------------------------------------------------------------
template< class value_type >
DLL_EXPORT void SMHiggsPotentialSplineDeSitter< value_type >::MultiPointSpline
(value_type* y, value_type* arrayOut, int numel, int nArrayToUse)
{
    switch(nArrayToUse)
    {
    case 0:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*yscale;//convert to GeV units
            int n = ChooseSpline(x, xdata, Ndata - 1);
            arrayOut[i] = SinglePointSpline(x, 0, n)*(y[i]*y[i]*y[i]*y[i])
                          /(value_type(4.0));
        }
        break;
    case 1:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*yscale;
            int n = ChooseSpline(x, xdata, Ndata - 1);
            value_type dldL = SinglePointSpline(x, 1, n);
            value_type l = SinglePointSpline(x, 0, n);
            value_type dldptp = value_type(2.0)*dldL;
            arrayOut[i] = (dldptp*(y[i]*y[i]*y[i])/(value_type(4.0)) +
            l*(y[i]*y[i]*y[i]));
        }
        break;
    case 2:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*yscale;
            int n = ChooseSpline(x, xdata, Ndata - 1);
            value_type d2ldL2 = SinglePointSpline(x, 2, n);
            value_type dldL = SinglePointSpline(x, 1, n);
            value_type l = SinglePointSpline(x, 0, n);
            value_type dldptp = value_type(2.0)*dldL;
            value_type d2ldp2tp2 = -value_type(2.0)*dldL
                                   + value_type(4.0)*d2ldL2;
            arrayOut[i] = (d2ldp2tp2*(y[i]*y[i])/(value_type(4.0)) +
            dldptp*(value_type(2.0)*y[i]*y[i]) +
            value_type(3.0)*y[i]*y[i]*l);
        }
        break;
    case 3:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*yscale;
            int n = ChooseSpline(x, xdata, Ndata - 1);
            value_type d3ldL3 = SinglePointSpline(x, 3, n);
            value_type d2ldL2 = SinglePointSpline(x, 2, n);
            value_type dldL = SinglePointSpline(x, 1, n);
            value_type l = SinglePointSpline(x, 0, n);
            value_type dldptp = value_type(2.0)*dldL;
            value_type d2ldp2tp2 = -value_type(2.0)*dldL
                                   + value_type(4.0)*d2ldL2;
            value_type d3ldp3tp3 = value_type(4.0)*dldL
                                   - value_type(12.0)*d2ldL2
                                   + value_type(8.0)*d3ldL3;
            arrayOut[i] = (d3ldp3tp3*(y[i])/(value_type(4.0)) +
            d2ldp2tp2*(value_type(3.0)*y[i]) +
            value_type(9.0)*y[i]*dldptp +
            value_type(6.0)*y[i]*l);
        }
        break;
    case 4:
        for (int i = 0; i < numel; i++)
        {
            value_type x = y[i]*yscale;
            int n = ChooseSpline(x, xdata, Ndata - 1);
            value_type d4ldL4 = SinglePointSpline(x, 4, n);
            value_type d3ldL3 = SinglePointSpline(x, 3, n);
            value_type d2ldL2 = SinglePointSpline(x, 2, n);
            value_type dldL = SinglePointSpline(x, 1, n);
            value_type l = SinglePointSpline(x, 0, n);
            value_type dldptp = value_type(2.0)*dldL;
            value_type d2ldp2tp2 = -value_type(2.0)*dldL
                                   + value_type(4.0)*d2ldL2;
            value_type d3ldp3tp3 = value_type(4.0)*dldL
                                   - value_type(12.0)*d2ldL2
                                   + value_type(8.0)*d3ldL3;
            value_type d4ldp4tp4 = -value_type(12.0)*dldL
                                   + value_type(44.0)*d2ldL2
                                   - value_type(48.0)*d3ldL3
                                   + value_type(16.0)*d4ldL4;
            arrayOut[i] = d4ldp4tp4/(value_type(4.0)) +
            d3ldp3tp3*(value_type(4.0)) +
            value_type(18.0)*d2ldp2tp2 +
            value_type(24.0)*dldptp +
            value_type(6.0)*l;
        }
        break;
    default:
        throw "Invalid derivative.\n";
    }
}
//------------------------------------------------------------------------------
template< class value_type >
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >::operator()
(value_type y)
{
    value_type x = y*yscale;
    int n = ChooseSpline(x, xdata, Ndata - 1);
    //Interpolate to get the couplings at this point:
    value_type g[8];
    //g = {g1^2, g2^2, gs^2, yt, l, m^2,xi, Z}
	SinglePointSpline(x, g , 0, n);
	//Evaluate the potential:
	value_type yz = g[7]*y;//Field with anomalous dimension.
	value_type yz2 = yz*yz;
	value_type V = -g[5]*yz2/_2p0 + g[6]*R*yz2/_2p0 + g[4]*(yz2*yz2)/_4p0;
	value_type Mi2;
	value_type ki[9] = {g[1]/_4p0,g[1]/_4p0,g[1]/_4p0,(g[1] + g[0])/_4p0,
                        (g[1] + g[0])/_4p0,(g[1] + g[0])/_4p0,g[3]*g[3]/_2p0,
                        _3p0*g[4],g[4]};
    value_type kip[9] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,g[5],g[5]};
	for(int i = 0;i < 9;i++)
    {
        Mi2 = abs(ki[i]*y2 - kip[i] + thetai[i]*R);
        V += ni[i]*Mi2*Mi2*(log(Mi2/(yz2 + R)) - ci[i])/one_loop_factor;
    }
    return V;
}
//------------------------------------------------------------------------------
template< class value_type >
DLL_EXPORT void SMHiggsPotentialSplineDeSitter< value_type >::operator()
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 0);
}
//------------------------------------------------------------------------------
//First derivatives:
template< class value_type >
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::d(value_type y)
{
    value_type x = y*yscale;
    int n = ChooseSpline(x, xdata, Ndata - 1);
    value_type dgdt[8];//Derivative with respect to t
    value_type g[8];
    //g = {g1^2, g2^2, gs^2, yt, l, m^2,xi, Z}
    SinglePointSpline(x, g , 0, n);
    SinglePointSpline(x, dgdt , 1, n);
    //Compute derivatives with respect to phi_c
    value_type dgdy[8];
    //Derivative of t with respect to classical field:
    //Follows from using the relation mu^2(t) = mt^2e^(2t)
    //And assigning mu^2(t) = Z(t)^2phi_c^2 + R. Differentiating
    //both sides with respect to phi_c and re-arranging gives:
    /*
    DLL_EXPORT value_type dtdy(value_type y,value_type R,value_type Z,
                               value_type Zd);
    DLL_EXPORT value_type d2tdy2(value_type y,value_type R,value_type Z,
                                 value_type Zd,value_type Zd2,value_type tp);
    DLL_EXPORT value_type d3tdy3(value_type y,value_type R,value_type Z,
                                 value_type Zd,value_type Zd2,value_type Zd3,
                                 value_type tp,value_type tp2);
    DLL_EXPORT value_type d4tdy4(value_type y,value_type R,value_type Z,
                                 value_type Zd,value_type Zd2,value_type Zd3,
                                 value_type,Zd4,value_type tp,value_type tp2,
                                 value_type tp3);
    */
    /*
    //value_type rc = g[7]*g[7]*y/(g[7]*g[7]*y*y + R - dgdt[7]*g[7]*y*y);
    value_type rc = this->dtdy(y,R,g[7],dgdt[7]);
    for(int i = 0;i < 8;i++)
    {
        dgdy[i] = dgdt[i]*rc;
    }
    //Compute the derivative of the potential:
    value_type yz = g[7]*y;//Field with anomalous dimension.
	value_type yz2 = yz*yz;
	value_type dyzdy = dgdt[7]*rc*y + g[7];
	//value_type V = -g[5]*yz2/_2p0 + g[6]*R*yz2/_2p0 + g[4]*(yz2*yz2)/_4p0;
    value_type Vp = -g[5]*yz*dyzdy -dgdy[5]*yz2/_2p0
                    + dgdy[6]*R*yz2/_2p0 + g[6]*R*yz*dyzdy
                     + dgdy[4]*(yz2*yz2)/_4p0 + g[4]*(yz2*yz*dyzdy);
    //Logarithmic part:
    value_type Mi2,Mi2p,mu2,mu2p;
	value_type ki[9] = {g[1]/_4p0,g[1]/_4p0,g[1]/_4p0,(g[1] + g[0])/_4p0,
                        (g[1] + g[0])/_4p0,(g[1] + g[0])/_4p0,g[3]*g[3]/_2p0,
                        _3p0*g[4],g[4]};
    value_type kip[9] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,g[5],g[5]};
    //Derivatives:
    value_type ki_d[9] = {dgdy[1]/_4p0,dgdy[1]/_4p0,dgdy[1]/_4p0,
                         (dgdy[1] + dgdy[0])/_4p0,(dgdy[1] + dgdy[0])/_4p0,
                         (dgdy[1] + dgdy[0])/_4p0,dgdy[3]*dgdy[3]/_2p0,
                        _3p0*dgdy[4],dgdy[4]};
    value_type kip_d[9] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,dgdy[5],dgdy[5]};
	for(int i = 0;i < 9;i++)
    {
        Mi2 = abs(ki[i]*yz2 - kip[i] + thetai[i]*R);
        Mi2p = (ki_d[i]*yz2 + ki[i]*yz*dyzdy - kip_d[i]);
        mu2 = abs(yz2 + R);
        mu2p = _2p0*yz*dyzdy;
        //V += ni[i]*Mi2*Mi2*(log(Mi2/(yz2 + R)) - ci[i])/one_loop_factor;
        Vp += (ni[i]/this->one_loop_factor)*(Mi2*Mi2p*(this->_1p0 -
                                             this->_2p0*ci[i] +
                                             this->_2p0*log(Mi2/mu2) )
                                             - Mi2*Mi2*mu2p/mu2);
    }
    return Vp;
}
//------------------------------------------------------------------------------
template< class value_type >
DLL_EXPORT void SMHiggsPotentialSplineDeSitter< value_type >::d
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 1);
}
//------------------------------------------------------------------------------
//Second derivatives:
template< class value_type >
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::d2(value_type y)
{
    value_type x = y*yscale;
    int n = ChooseSpline(x, xdata, Ndata - 1);
    value_type dgdt[8];//Derivative with respect to t
    value_type d2gdt2[8];//Second derivative with respect to t
    value_type g[8];
    //g = {g1^2, g2^2, gs^2, yt, l, m^2,xi, Z}
    SinglePointSpline(x, g , 0, n);
    SinglePointSpline(x, dgdt , 1, n);
    SinglePointSpline(x, d2gdt2 , 2, n);
    //Compute derivatives with respect to phi_c
    value_type dgdy[8];
    value_type d2gdy2[8];
    //Derivative of t with respect to classical field:
    //Follows from using the relation mu^2(t) = mt^2e^(2t)
    //And assigning mu^2(t) = Z(t)^2phi_c^2 + R. Differentiating
    //both sides with respect to phi_c and re-arranging gives:
    value_type denom = (g[7]*g[7]*y*y + R - dgdt[7]*g[7]*y*y);
    value_type Z2 = g[7]*g[7];
    value_type Z = g[7];
    value_type Zd = dgdt[7];//derivative of Z wrt t.
    value_type Zd2 = d2gdt2[7];//second derivative of Z wrt t.
    value_type dtdy = Z2*y/denom;
    value_type Zp = Zd*dtdy;//derivative of Z wrt y
    value_type d2tdy2 = (Z2/denom) +
                         this->_2p0*Zp*Z*y/denom -
                          (y*Z2*(this->_2p0*y*Z2 - this->_2p0*y*Z*Zd
                           + this->_2p0*y*y*Z*Zp - y*y*Zd*Zp - y*y*Z*Zd2*dtdy ))
                           /(denom*denom);
    for(int i = 0;i < 8;i++)
    {
        dgdy[i] = dgdt[i]*dtdy;
        d2gdy2[i] = d2gdt2[i]*dtdy*dtdy + dgdt[i]*d2tdy2;
    }
    //Compute the derivative of the potential:
    value_type yz = g[7]*y;//Field with anomalous dimension.
	value_type yz2 = yz*yz;
	value_type dyzdy = dgdt[7]*dtdy*y + g[7];
	value_type d2yzdy2 = d2gdt2[7]*dtdy*dtdy*y + dgdt[7]*d2tdy2*y
                        + this->_2p0*dgdt[7]*dtdy;
	//value_type V = -g[5]*yz2/_2p0 + g[6]*R*yz2/_2p0 + g[4]*(yz2*yz2)/_4p0;
    value_type Vpp = -this->_2p0*dgdy[5]*yz*dyzdy
                    -g[5]*dyzdy*dyzdy -g[5]*yz*d2yzdy2
                    -d2gdy2[5]*yz2/_2p0
                    + d2gdy2[6]*R*yz2/_2p0 + this->_2p0*dgdy[6]*R*yz*dyzdy
                    + g[6]*R*dyzdy*dyzdy
                    + g[6]*R*yz*d2yzdy2
                    + d2gdy2[4]*(yz2*yz2)/_4p0
                    + this->_2p0*dgdy[4]*(yz2*yz*dyzdy)
                    + g[4]*(this->_3p0*yz*dyzdy*yz*dyzdy)
                    + g[4]*(yz2*yz*d2yzdy2);
    //Logarithmic part:
    value_type Mi2,Mi2p,mu2,mu2p;
	value_type ki[9] = {g[1]/_4p0,g[1]/_4p0,g[1]/_4p0,(g[1] + g[0])/_4p0,
                        (g[1] + g[0])/_4p0,(g[1] + g[0])/_4p0,g[3]*g[3]/_2p0,
                        _3p0*g[4],g[4]};
    value_type kip[9] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,g[5],g[5]};
    //Derivatives:
    value_type ki_d[9] = {dgdy[1]/_4p0,dgdy[1]/_4p0,dgdy[1]/_4p0,
                         (dgdy[1] + dgdy[0])/_4p0,(dgdy[1] + dgdy[0])/_4p0,
                         (dgdy[1] + dgdy[0])/_4p0,dgdy[3]*dgdy[3]/_2p0,
                        _3p0*dgdy[4],dgdy[4]};
    value_type kip_d[9] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,dgdy[5],dgdy[5]};
    value_type ki_d2[9] = {d2gdy2[1]/_4p0,d2gdy2[1]/_4p0,d2gdy2[1]/_4p0,
                         (d2gdy2[1] + d2gdy2[0])/_4p0,
                         (d2gdy2[1] + d2gdy2[0])/_4p0,
                         (d2gdy2[1] + d2gdy2[0])/_4p0,d2gdy2[3]*d2gdy2[3]/_2p0,
                        _3p0*d2gdy2[4],d2gdy2[4]};
    value_type kip_d2[9] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,d2gdy2[5],
                            d2gdy2[5]};

	for(int i = 0;i < 9;i++)
    {
        Mi2 = abs(ki[i]*yz2 - kip[i] + thetai[i]*R);
        Mi2p = (ki_d[i]*yz2 + ki[i]*yz*dyzdy - kip_d[i]);
        Mi2pp = (ki_d2[i]*yz2 + ki_d[i]*this->_3p0*yz*dyzdy
                  + ki[i]*dyzdy*dyzdy + ki[i]*yz*d2yzdy2
                   - kip_d2[i]);
        mu2 = abs(yz2 + R);
        mu2p = _2p0*yz*dyzdy;
        mu2pp = _2p0*dyzdy*dyzdy + _2p0*yz*d2yzdy2;
        //V += ni[i]*Mi2*Mi2*(log(Mi2/(yz2 + R)) - ci[i])/one_loop_factor;
        value_type log_factor = log(Mi2/mu2);
        Vpp += (ni[i]/this->one_loop_factor)*(Mi2*Mi2*mu2p*mu2p/(mu2*mu2) +
                                             (this->_3p0 - this->_2p0*ci[i] +
                                              this->_2p0*log_factor)*Mi2p*Mi2p +
                                             (this->_1p0 - this->_2p0*ci[i]
                                               + this->_2p0*log_factor)*
                                             Mi2*Mi2pp
                                              - Mi2*(this->_4p0*Mi2p*mu2p
                                                      + Mi2*mu2pp)/mu2 );
    }
    return Vpp;
}
//------------------------------------------------------------------------------
template< class value_type >
DLL_EXPORT void SMHiggsPotentialSplineDeSitter< value_type >::d2
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 2);
}
//------------------------------------------------------------------------------
//Third derivatives:
template< class value_type >
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::d3(value_type y)
{
    value_type x = y*yscale;
    int n = ChooseSpline(x, xdata, Ndata - 1);
    value_type dgdt[8];//Derivative with respect to t
    value_type d2gdt2[8];//Second derivative with respect to t
    value_type d3gdt3[8];//Third derivative with respect to t
    value_type g[8];
    //g = {g1^2, g2^2, gs^2, yt, l, m^2,xi, Z}
    SinglePointSpline(x, g , 0, n);
    SinglePointSpline(x, dgdt , 1, n);
    SinglePointSpline(x, d2gdt2 , 2, n);
    SinglePointSpline(x, d3gdt3 , 3, n);

    //Compute derivatives with respect to phi_c
    value_type dgdy[8];
    value_type d2gdy2[8];
    value_type d3gdy3[8];
    //Derivative of t with respect to classical field:
    //Follows from using the relation mu^2(t) = mt^2e^(2t)
    //And assigning mu^2(t) = Z(t)^2phi_c^2 + R. Differentiating
    //both sides with respect to phi_c and re-arranging gives:
    value_type denom = (g[7]*g[7]*y*y + R - dgdt[7]*g[7]*y*y);
    value_type Z2 = g[7]*g[7];
    value_type Z = g[7];
    value_type Zd = dgdt[7];//derivative of Z wrt t.
    value_type Zd2 = d2gdt2[7];//second derivative of Z wrt t.
    value_type Zd3 = d3gdt3[7];
    value_type dtdy = Z2*y/denom;
    value_type Zp = Zd*dtdy;//derivative of Z wrt y
    value_type d2tdy2 = (Z2/denom) +
                         this->_2p0*Zp*Z*y/denom -
                          (y*Z2*(this->_2p0*y*Z2 - this->_2p0*y*Z*Zd
                           + this->_2p0*y*y*Z*Zp - y*y*Zd*Zp - y*y*Z*Zd2*dtdy ))
                           /(denom*denom);
    value_type num_factor_1 = ( this->_2p0*y*Z2 - this->_2p0*y*Z*Zd +
                                this->_2p0*y*y*Z*Zp - y*y*Zd*Zp
                                 - y*y*Z*dydt*Zd2 );
    value_type d3tdy3 = (this->_4p0**Z*Zp + this->_2p0*y*Zp*Zp +
                         this->_2p0*y*Z*Zd*d2tdy2
                         + this->_2p0*y*Z*dtdy*dtdy*Zd2)/denom
                         - (this->_2p0*Z2*num_factor_1)/(denom*denom)
                         - (this->_4p0*y*Z*Zp*num_factor_1)/(denom*denom)
                         + (this->_2p0*y*Z2*num_factor_1)/(denom*denom*denom)
                         - (y*Z2*( this->_2p0*Z2 - this->_2p0*Z*Zd
                                  + this->_8p0*y*Z*Zp - this->_4p0*y*Zd*Zp
                                  + this->_2p0*y*y*Zp*Zp
                                  + this->_2p0*y*yZ*Zd*d2tdy2 - y*y*Zd*d2tdy2
                                  - this->_4p0*y*Z*dtdy*Zd2
                                  + this->_2p0*y*y*Z*dtdy*dtdy*Zd2
                                  - this->_3p0*y*y*dtdy*Zp*Zd2
                                  - y*y*Z*d2tdy2*Zd2
                                  - y*y*Z*dtdy*dtdy*Zd3 ) )/(denom*denom);
    for(int i = 0;i < 8;i++)
    {
        dgdy[i] = dgdt[i]*dtdy;
        d2gdy2[i] = d2gdt2[i]*dtdy*dtdy + dgdt[i]*d2tdy2;
    }
    //Compute the derivative of the potential:
    value_type yz = g[7]*y;//Field with anomalous dimension.
	value_type yz2 = yz*yz;
	value_type dyzdy = dgdt[7]*dtdy*y + g[7];
	value_type d2yzdy2 = d2gdt2[7]*dtdy*dtdy*y + dgdt[7]*d2tdy2*y
                        + this->_2p0*dgdt[7]*dtdy;
	//value_type V = -g[5]*yz2/_2p0 + g[6]*R*yz2/_2p0 + g[4]*(yz2*yz2)/_4p0;
    value_type Vpp = -this->_2p0*dgdy[5]*yz*dyzdy
                    -g[5]*dyzdy*dyzdy -g[5]*yz*d2yzdy2
                    -d2gdy2[5]*yz2/_2p0
                    + d2gdy2[6]*R*yz2/_2p0 + this->_2p0*dgdy[6]*R*yz*dyzdy
                    + g[6]*R*dyzdy*dyzdy
                    + g[6]*R*yz*d2yzdy2
                    + d2gdy2[4]*(yz2*yz2)/_4p0
                    + this->_2p0*dgdy[4]*(yz2*yz*dyzdy)
                    + g[4]*(this->_3p0*yz*dyzdy*yz*dyzdy)
                    + g[4]*(yz2*yz*d2yzdy2);
    //Logarithmic part:
    value_type Mi2,Mi2p,mu2,mu2p;
	value_type ki[9] = {g[1]/_4p0,g[1]/_4p0,g[1]/_4p0,(g[1] + g[0])/_4p0,
                        (g[1] + g[0])/_4p0,(g[1] + g[0])/_4p0,g[3]*g[3]/_2p0,
                        _3p0*g[4],g[4]};
    value_type kip[9] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,g[5],g[5]};
    //Derivatives:
    value_type ki_d[9] = {dgdy[1]/_4p0,dgdy[1]/_4p0,dgdy[1]/_4p0,
                         (dgdy[1] + dgdy[0])/_4p0,(dgdy[1] + dgdy[0])/_4p0,
                         (dgdy[1] + dgdy[0])/_4p0,dgdy[3]*dgdy[3]/_2p0,
                        _3p0*dgdy[4],dgdy[4]};
    value_type kip_d[9] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,dgdy[5],dgdy[5]};
    value_type ki_d2[9] = {d2gdy2[1]/_4p0,d2gdy2[1]/_4p0,d2gdy2[1]/_4p0,
                         (d2gdy2[1] + d2gdy2[0])/_4p0,
                         (d2gdy2[1] + d2gdy2[0])/_4p0,
                         (d2gdy2[1] + d2gdy2[0])/_4p0,d2gdy2[3]*d2gdy2[3]/_2p0,
                        _3p0*d2gdy2[4],d2gdy2[4]};
    value_type kip_d2[9] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,d2gdy2[5],
                            d2gdy2[5]};

	for(int i = 0;i < 9;i++)
    {
        Mi2 = abs(ki[i]*yz2 - kip[i] + thetai[i]*R);
        Mi2p = (ki_d[i]*yz2 + ki[i]*yz*dyzdy - kip_d[i]);
        Mi2pp = (ki_d2[i]*yz2 + ki_d[i]*this->_3p0*yz*dyzdy
                  + ki[i]*dyzdy*dyzdy + ki[i]*yz*d2yzdy2
                   - kip_d2[i]);
        mu2 = abs(yz2 + R);
        mu2p = _2p0*yz*dyzdy;
        mu2pp = _2p0*dyzdy*dyzdy + _2p0*yz*d2yzdy2;
        //V += ni[i]*Mi2*Mi2*(log(Mi2/(yz2 + R)) - ci[i])/one_loop_factor;
        value_type log_factor = log(Mi2/mu2);
        Vpp += (ni[i]/this->one_loop_factor)*(Mi2*Mi2*mu2p*mu2p/(mu2*mu2) +
                                             (this->_3p0 - this->_2p0*ci[i] +
                                              this->_2p0*log_factor)*Mi2p*Mi2p +
                                             (this->_1p0 - this->_2p0*ci[i]
                                               + this->_2p0*log_factor)*
                                             Mi2*Mi2pp
                                              - Mi2*(this->_4p0*Mi2p*mu2p
                                                      + Mi2*mu2pp)/mu2 );
    }
    return Vpp;
}
//------------------------------------------------------------------------------
template< class value_type >
void SMHiggsPotentialSplineDeSitter< value_type >::d3
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 3);
}
//------------------------------------------------------------------------------
//Fourth derivatives:
template< class value_type >
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::d4(value_type y)
{
    value_type x = y*yscale;
    int n = ChooseSpline(x, xdata, Ndata - 1);
    value_type d4ldL4 = SinglePointSpline(x, 4, n);
    value_type d3ldL3 = SinglePointSpline(x, 3, n);
    value_type d2ldL2 = SinglePointSpline(x, 2, n);
    value_type dldL = SinglePointSpline(x, 1, n);
    value_type l = SinglePointSpline(x, 0, n);
    value_type dldptp = value_type(2.0)*dldL;
    value_type d2ldp2tp2 = -value_type(2.0)*dldL + value_type(4.0)*d2ldL2;
    value_type d3ldp3tp3 = value_type(4.0)*dldL - value_type(12.0)*d2ldL2
               + value_type(8.0)*d3ldL3;
    value_type d4ldp4tp4 = -value_type(12.0)*dldL + value_type(44.0)*d2ldL2
               - value_type(48.0)*d3ldL3 + value_type(16.0)*d4ldL4;
	return d4ldp4tp4/(value_type(4.0)) + d3ldp3tp3*(value_type(4.0))
               + value_type(18.0)*d2ldp2tp2 + value_type(24.0)*dldptp
               + value_type(6.0)*l;
}
//------------------------------------------------------------------------------
template< class value_type >
void SMHiggsPotentialSplineDeSitter< value_type >::d4
(value_type* y, value_type* arrayOut, int numel)
{
	MultiPointSpline(y, arrayOut, numel, 4);
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::a(value_type Z,value_type y)
{
    return this->alpha*Z*Z*y*y;
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::b(value_type Z,value_type Zd,value_type R,value_type y)
{
    return this->alpha*Z*Z*y*y +this->beta*R- this->_2p0*alpha*Z*Zd*y*y;
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::ap(value_type Z,value_type Zd,value_type tp,value_type y)
{
    return this->_2p0*this->alpha*Z*Zd*tp*y*y
                    + this->_2p0*this->alpha*Z*Z*y;
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::bp(value_type Z,value_type Zd,value_type Zd2,value_type tp,value_type R,
     value_type y)
{
    return this->_2p0*this->alpha*y*Z*Z - this->_4p0*this->alpha*y*Z*Zd
                     + this->_2p0*this->alpha*y*y*Z*tp*Zd
                     - this->_2p0*this->alpha*y*y*tp*Zd*Zd
                     - this->_2p0*this->_alpha*y*y*Z*tp*Zd2;
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::ap2(value_type Z,value_type Zd,value_type Zd2,value_type tp,value_type tp2,
      value_type y)
{
    return this->_2p0*this->alpha*Z*Z
                    + this->_8p0*this->alpha*y*Z*tp*Zd
                    + this->_2p0*this->alpha*y*y*tp*tp*Zd*Zd
                    + this->_2p0*this->alpha*y*y*Z*Zd*tp2
                    + this->_2p0*this->alpha*y*y*Z*tp*tp*Zd2;
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::bp2(value_type Z,value_type Zd,value_type Zd2,value_type Zd3,value_type tp,
      value_type tp2,value_type R,value_type y)
{
    return this->_2p0*this->alpha*Z*Z - this->_4p0*this->alpha*Z*Zd
                     + this->_8p0*this->alpha*this->y*this->Z*tp*Zd
                     - this->_8p0*this->alpha*y*tp*Zd*Zd
                     + this->_2p0*this->alpha*y*y*tp*tp*Zd*Zd
                     + this->_2p0*this->alpha*y*y*Z*Zd*tp2
                     - this->_2p0*this->alpha*y*y*Zd*Zd*tp2
                     - this->_8p0*this->alpha*y*Z*tp*Zd2
                     + this->_2p0*this->alpha*y*y*Z*tp*tp*Zd2
                     - this->_6p0*this->alpha*y*y*tp*tp*Zd*Zd2
                     - this->_2p0*this->alpha*y*y*Z*tp2*Zd2
                     - this->_2p0*this->alpha*y*y*Z*tp*tp*Zd3;
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::ap3(value_type Z,value_type Zd,value_type Zd2,value_type Zd3,value_type tp,
      value_type tp2,value_type tp3,value_type y)
{
    const value_type _12p0 = value_type(12.0);
    const value_type _6p0 = value_type(6.0);
    return this->_12p0*this->alpha*Z*tp*Zd
    + _12p0*this->alpha*y*tp*tp*Zd*Zd + _12p0*this->alpha*y*Z*Zd*tp2
    + _6p0*this->alpha*y*y*tp*Zd*Zd*tp2
    + _12p0*this->alpha*tp*tp*Zd2
    + _6p0*this->alpha*y*y*tp*tp*tp*Zd*Zd2
    + _6p0*this->alpha*y*y*Z*tp*tp2*Zd2
    + this->_2p0*this->alpha*Z*Zd*tp3
    + this->_2p0*this->alpha*y*y*Z*tp*tp*tp*Zd3;
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::bp3(value_type Z,value_type Zd,value_type Zd2,value_type Zd3,value_type tp,
      value_type tp2,value_type tp3,value_type R,value_type y)
{
    const value_type _12p0 = value_type(12.0);
    const value_type _6p0 = value_type(6.0);
    const value_type _36p0 = value_type(36.0);
    const value_type _18p0 = value_type(18.0);
    const value_type _8p0 = value_type(8.0);
    return _12p0*this->alpha*Z*tp*Zd - _12p0*this->alpha*tp*Zd*Zd
    +_12p0*this->alpha*y*tp*tp*Zd*Zd + _12p0*this->alpha*y*Z*Zd*tp2
    -_12p0*this->alpha*y*Zd*Zd*tp2 + _6p0*this->alpha*y*y*tp*Zd*Zd*tp2
    -_12p0*this->alpha*Z*tp*Zd2 + _12p0*this->alpha*y*Z*tp*tp*Zd2
    -_36p0*this->alpha*y*tp*tp*Zd*Zd2 + _6p0*this->alpha*y*y*tp*tp*tp*Zd*Zd2
    -_12p0*this->alpha*y*Z*tp2*Zd2 + _6p0*this->alpha*y*y*Z*tp*tp2*Zd2
    -_18p0*this->alpha*y*y*tp*Zd*tp2*Zd2 - _6p0*this->alpha*y*y*tp*tp*tp*Zd2*Zd2
    +this->_2p0*this->alpha*y*y*Z*Zd*tp3
    -this->_2p0*this->alpha*y*y*Zd*Zd*tp3
    -this->_2p0*this->alpha*y*y*Z*Zd2*tp3
    -_12p0*this->alpha*y*Z*tp*tp*Zd3
    +this->_2p0*this->alpha*y*y*Z*tp*tp*tp*Zd3
    - _8p0*this->alpha*y*y*tp*tp*tp*Zd*Zd3
    -_6p0*this->alpha*y*y*Z*tp*tp2*Zd3
    -this->_2p0*this->alpha*y*y*Z*tp*tp*tp*Zd4;
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::dtdy(value_type y,value_type R,value_type Z,value_type Zd)
{
    value_type A = this->a(Z,y);
    value_type B = this->b(Z,Zd,R,y);
    return A/B;
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::d2tdy2(value_type y,value_type R,value_type Z,value_type Zd,value_type Zd2,
         value_type tp)
{
    value_type A = this->a(Z,y);
    value_type B = this->b(Z,Zd,R,y);
    value_type Ap = this->ap(Z,Zd,tp,y);
    value_type Bp = this->bp(Z,Zd,Zd2,tp,R,y);
    return Ap/B - A*Bp/(B*B);
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::d3tdy3(value_type y,value_type R,value_type Z,value_type Zd,value_type Zd2,
         value_type Zd3,value_type tp,value_type tp2)
{
    value_type A = this->a(Z,y);
    value_type B = this->b(Z,Zd,R,y);
    value_type Ap = this->ap(Z,Zd,tp,y);
    value_type Bp = this->bp(Z,Zd,Zd2,tp,R,y);
    value_type Ap2 = this->ap2(Z,Zd,Zd2,tp,tp2,y);
    value_type Bp2 = this->bp2(Z,Zd,Zd2,Zd3,tp,tp2,R,y);
    return -this->_2p0*Ap*Bp/(B*B) + this->_2p0*A*Bp*Bp/(B*B*B)
            + Ap2/B - A*Bp2/(B*B);
}
//------------------------------------------------------------------------------
DLL_EXPORT value_type SMHiggsPotentialSplineDeSitter< value_type >
::d4tdy4(value_type y,value_type R,value_type Z,value_type Zd,value_type Zd2,
         value_type Zd3,value_type,Zd4,value_type tp,value_type tp2,
         value_type tp3)
{
    const value_type _6p0 = value_type(6.0);
    value_type A = this->a(Z,y);
    value_type B = this->b(Z,Zd,R,y);
    value_type Ap = this->ap(Z,Zd,tp,y);
    value_type Bp = this->bp(Z,Zd,Zd2,tp,R,y);
    value_type Ap2 = this->ap2(Z,Zd,Zd2,tp,tp2,y);
    value_type Bp2 = this->bp2(Z,Zd,Zd2,Zd3,tp,tp2,R,y);
    value_type Ap3 = this->ap3(Z,Zd,Zd2,Zd3,tp,tp2,tp3,y);
    value_type Bp3 = this->bp3(Z,Zd,Zd2,Zd3,tp,tp2,tp3,R,y);
    return _6p0*Ap*Bp*Bp/(B*B*B) - _6p0*A*Bp*Bp*Bp/(B*B*B*B)
            - this->_3p0*Bp*Ap2/(B*B) - this->_3p0*(Ap*Bp2)/(B*B)
            + _6p0*A*Bp*Bp2/(B*B*B) + Ap3/B - A*Bp3/(B*B);
}
//------------------------------------------------------------------------------
*/
//==============================================================================
//Log potential:
//------------------------------------------------------------------------------
template <class value_type>
DLL_EXPORT value_type log_potential<value_type>::operator()(value_type x)
{
    if(x != this->_0p0)
    {
        value_type x2 = x*x;
        value_type log_factor = log(x2/(this->M*this->M));
        return (this->lambda_0 + log_factor*(this->b*log_factor - this->a))*
                x2*x2/this->_4p0;
    }
    else
    {
        return this->_0p0;
    }
}
//------------------------------------------------------------------------------
template <class value_type>
DLL_EXPORT void log_potential<value_type>::operator()
(value_type* x,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel;i++)
    {
        arrayOut[i] = this->operator()(x[i]);
    }
}
//------------------------------------------------------------------------------
template <class value_type>
DLL_EXPORT value_type log_potential<value_type>::d(value_type x)
{
    if(x != this->_0p0)
    {
        value_type x2 = x*x;
        value_type log_factor = log(x2/(this->M*this->M));
        return -( a - this->_2p0*lambda_0 + log_factor*(this->_2p0*(a - b)
                    - this->_2p0*b*log_factor) )*
                x2*x/this->_2p0;
    }
    else
    {
        return this->_0p0;
    }
}
//------------------------------------------------------------------------------
template <class value_type>
DLL_EXPORT void log_potential<value_type>::d
(value_type* x,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel;i++)
    {
        arrayOut[i] = this->d(x[i]);
    }
}
//------------------------------------------------------------------------------

template <class value_type>
DLL_EXPORT value_type log_potential<value_type>::d2(value_type x)
{
    if(x != this->_0p0)
    {
        value_type x2 = x*x;
        value_type log_factor = log(x2/(this->M*this->M));
        return -(this->_7p0*a - this->_4p0*b - this->_6p0*lambda_0
                  + log_factor*(this->_2p0*(this->_3p0*a - this->_7p0*b) +
                                - this->_6p0*b*log_factor))*
                x2/this->_2p0;
    }
    else
    {
        return this->_0p0;
    }
}
//------------------------------------------------------------------------------
template <class value_type>
DLL_EXPORT void log_potential<value_type>::d2
(value_type* x,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel;i++)
    {
        arrayOut[i] = this->d2(x[i]);
    }
}
//------------------------------------------------------------------------------

template <class value_type>
DLL_EXPORT value_type log_potential<value_type>::d3(value_type x)
{
    if(x != this->_0p0)
    {
        value_type x2 = x*x;
        value_type log_factor = log(x2/(this->M*this->M));
        return (-this->_13p0*a + this->_6p0*(this->_3p0*b + this->_10p0) +
                log_factor*( (-this->_6p0*a + this->_26p0*b)
                             + this->_6p0*b*log_factor ) )*x;
    }
    else
    {
        return this->_0p0;
    }
}
//------------------------------------------------------------------------------
template <class value_type>
DLL_EXPORT void log_potential<value_type>::d3
(value_type* x,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel;i++)
    {
        arrayOut[i] = this->d3(x[i]);
    }
}
//------------------------------------------------------------------------------
template <class value_type>
DLL_EXPORT value_type log_potential<value_type>::d4(value_type x)
{
    value_type x2 = x*x;
    value_type log_factor = log(x2/(this->M*this->M));
    return (-this->_25p0*a + this->_70p0*b + this->_6p0*lambda_0 +
            log_factor*((-this->_6p0*a + this->_50p0*b)
                         + this->_6p0*b*log_factor)  );
}
//------------------------------------------------------------------------------
template <class value_type>
DLL_EXPORT void log_potential<value_type>::d4
(value_type* x,value_type* arrayOut,int numel)
{
    for(int i = 0;i < numel;i++)
    {
        arrayOut[i] = this->d4(x[i]);
    }
}
//------------------------------------------------------------------------------


//==============================================================================
#endif
