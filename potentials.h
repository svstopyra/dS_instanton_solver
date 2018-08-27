#ifndef POTENTIALS_FORWARD
#define POTENTIALS_FORWARD

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

#include "pchip.h"

//Forward declaration of the potentials. Does not contain any code.

//Selects the potential we want to use with the Instanton solver. Takes an
//integer and returns a function pointer.
//Format for a potential function:
/*
 *Should have an operator() for both value_type OR vector<value_type> giving the
 *potential value itself.
 *Should have derivatives (and vectorised derivatives) for each as well, up to
 *fourth order.
 */
//==============================================================================
//------------------------------------------------------------------------------
//Base class - basic potential with no derivatives specified.
template< class value_type >
class potential_basic
{
public:
    virtual value_type operator()(value_type x) =0;
};
//Extended potential - has first derivatives, and some other things we need
//for the instanton solver (first derivatives, and bools to record whether
//the potential has a max/min value).
template< class value_type >
class potential : public potential_basic<value_type>
{
public:

	//DLL_EXPORT virtual void operator()(value_type* x,value_type* arrayOut,
    //                                   int numel) =0;


    //1st derivative of the potential:
	virtual value_type d(value_type x) =0;

	//DLL_EXPORT virtual void d(value_type* x,value_type* arrayOut,int numel) =0;
	virtual value_type d2(value_type x) =0;
	/*
	DLL_EXPORT virtual void d2(value_type* x,value_type* arrayOut,int numel) =0;
	DLL_EXPORT virtual value_type d3(value_type x) =0;
	DLL_EXPORT virtual void d3(value_type* x,value_type* arrayOut,int numel) =0;
	DLL_EXPORT virtual value_type d4(value_type x) = 0;
	DLL_EXPORT virtual void d4(value_type* x,value_type* arrayOut,int numel) =0;
	*/
	//Parameters if the potential has a minimum and maximum.
	bool hasMaximum;
	bool hasMinimum;
	value_type Maximum;
	value_type Minimum;
	//Constructor:
	potential()
	{
	    //Default to false. Altered by derived class constructors.
        hasMaximum = false;
        hasMinimum = false;
        Maximum = value_type(0.0);
        Minimum = value_type(0.0);
    }
};
//Abstract classes for potentials with higher derivatives defined:
template<class value_type>
class potential_with_2nd_derivatives : public potential<value_type>
{
    //2nd derivative of the potential:
    DLL_EXPORT virtual value_type d2(value_type x) =0;
};
template<class value_type>
class potential_with_3nd_derivatives
: public potential_with_2nd_derivatives<value_type>
{
    //3rd derivative of the potential:
    DLL_EXPORT virtual value_type d3(value_type x) =0;
};
template<class value_type>
class potential_with_4th_derivatives
: public potential_with_3nd_derivatives<value_type>
{
    //3rd derivative of the potential:
    DLL_EXPORT virtual value_type d4(value_type x) =0;
};
//------------------------------------------------------------------------------
//==============================================================================
//Potential used to study new physics operators.
#ifdef USING_POTENTIAL_NEW_PHYSICS
template< class value_type >
class NewPhysicsPotential : public potential< value_type >
{
private:
    value_type m2;//Mass squared term, m^2\phi^2/2
    value_type l4; //Quartic coupling
    value_type l6; //phi^6 coupling
    value_type l8; //phi^8 coupling
    value_type M6; //phi^6 mass scale
    value_type M8; //phi^8 mass scale
    value_type yscale;
public:
    //Constructors:
    NewPhysicsPotential(value_type M2,value_type lambda_4,value_type lambda_6,
                        value_type lambda_8,value_type m6,value_type m8,
                        value_type YSCALE)
    {
        m2 = M2/(YSCALE*YSCALE);
        l4 = lambda_4;
        l6 = lambda_6;
        l8 = lambda_8;
        M6 = m6/YSCALE;
        M8 = m8/YSCALE;
        yscale = YSCALE;
        //Note - the potential expects these parameters in units of Mp. yscale
        //sets the scale that the potential is asked to give the answer V(y) in
        //and expects the argument, y, to be in, relative to the Planck scale.
        //Thus if we want the answers in terms of scale H, then we set
        //yscale = H/Mp; The constructor above re-scales the parameters
        //automatically to be in units of H (defined by yscale). Note that
        //the lambda parameters are all dimensionless.
    }
    //Potential and its derivatives:
    //Potential itself:
	DLL_EXPORT value_type operator()(value_type y);
	DLL_EXPORT void operator()(value_type* y,value_type* arrayOut,int numel);
    //First derivatives:
	DLL_EXPORT value_type d(value_type y);
	DLL_EXPORT void d(value_type* y,value_type* arrayOut,int numel);
    //Second derivatives:
	DLL_EXPORT value_type d2(value_type y);
	DLL_EXPORT void d2(value_type* y,value_type* arrayOut,int numel);
    //Third derivatives:
	DLL_EXPORT value_type d3(value_type y);
	DLL_EXPORT void d3(value_type* y,value_type* arrayOut,int numel);
    //Fourth derivatives:
	DLL_EXPORT value_type d4(value_type y);
	DLL_EXPORT void d4(value_type* y,value_type* arrayOut,int numel);
};
#endif
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//==============================================================================
//Derived classes:
#ifdef USING_POTENTIAL_TEST1
template< class value_type >
class testPotential : public potential< value_type >
{
private:
    value_type a;
    value_type b;
    value_type c;
    value_type m;
public:
    value_type yscale;//public, so that we can dynamically adjust the scale if
    //necessary.
    //Constructors
	testPotential(value_type& A, value_type& B, value_type& C, value_type& M,
                  value_type& YSCALE);
    //Potential itself:
	DLL_EXPORT value_type operator()(value_type y);
	DLL_EXPORT void operator()(value_type* y,value_type* arrayOut,int numel);
    //First derivatives:
	DLL_EXPORT value_type d(value_type y);
	DLL_EXPORT void d(value_type* y,value_type* arrayOut,int numel);
    //Second derivatives:
	DLL_EXPORT value_type d2(value_type y);
	DLL_EXPORT void d2(value_type* y,value_type* arrayOut,int numel);
    //Third derivatives:
	DLL_EXPORT value_type d3(value_type y);
	DLL_EXPORT void d3(value_type* y,value_type* arrayOut,int numel);
    //Fourth derivatives:
	DLL_EXPORT value_type d4(value_type y);
	DLL_EXPORT void d4(value_type* y,value_type* arrayOut,int numel);
};
//------------------------------------------------------------------------------
template< class value_type >
testPotential< value_type >::testPotential
(value_type& A, value_type& B, value_type& C, value_type& M, value_type& YSCALE)
{
	a = A;
	b = B;
	c = C;
	m = M;
	yscale = YSCALE;
}
//------------------------------------------------------------------------------
//==============================================================================
#endif
#ifdef USING_POTENTIAL_BRANCHINA
template< class value_type >
class branchinaPotential : public potential< value_type >
{
private:
    value_type a;
    value_type b;
    value_type g;
public:
    value_type yscale;//public, so that we can dynamically adjust the scale if
    //necessary.
    //Constructors
	branchinaPotential(value_type& A, value_type& B, value_type& G,
                       value_type& YSCALE);
    //Potential itself:
	DLL_EXPORT value_type operator()(value_type y);
	DLL_EXPORT void operator()(value_type* y,value_type* arrayOut,int numel);
    //First derivatives:
	DLL_EXPORT value_type d(value_type y);
	DLL_EXPORT void d(value_type* y,value_type* arrayOut,int numel);
    //Second derivatives:
	DLL_EXPORT value_type d2(value_type y);
	DLL_EXPORT void d2(value_type* y,value_type* arrayOut,int numel);
    //Third derivatives:
	DLL_EXPORT value_type d3(value_type y);
	DLL_EXPORT void d3(value_type* y,value_type* arrayOut,int numel);
    //Fourth derivatives:
	DLL_EXPORT value_type d4(value_type y);
	DLL_EXPORT void d4(value_type* y,value_type* arrayOut,int numel);
};
//------------------------------------------------------------------------------
template< class value_type >
branchinaPotential< value_type >::branchinaPotential
(value_type& A, value_type& B, value_type& G, value_type& YSCALE)
{
	a = A;
	b = B;
	g = G;
	yscale = YSCALE;
}
#endif
//------------------------------------------------------------------------------
//==============================================================================
#ifdef USING_SM_POTENTIAL
template< class value_type >
class SplinePotential : public potential< value_type >
{
protected:
    value_type* xdata;
     value_type M;
     int Ndata;
public:
    const value_type _3p0;// = value_type(3.0);
    const value_type _2p0;// = value_type(2.0);
    const value_type _1p0;// = value_type(1.0);
    const value_type _6p0;// = value_type(6.0);
    const value_type _4p0;// = value_type(4.0);
    const value_type _5p0;// = value_type(5.0);
    const value_type _12p0;// = value_type(12.0);
    const value_type _0p0;
     //Member functions:
     SplinePotential(value_type* XDATA,int NDATA,value_type m)
     : _3p0(3.0) , _1p0(1.0), _6p0(6.0), _5p0(5.0), _12p0(12.0), _2p0(2.0),
     _4p0(4.0), _0p0(0.0)
     {
         xdata = XDATA;
         Ndata = NDATA;
         M = m;
     }
     //Compute a spline for array (derivative) nArrayToUse, at point y.
     /*
	 DLL_EXPORT virtual value_type SinglePointSpline(value_type& y,
                                                     int nArrayToUse);
	 DLL_EXPORT virtual value_type SinglePointSpline(value_type& y,
                                                     int nArrayToUse,
                                                     int nSpline);
                                                     */
	 //Compute a spline of the  nArrayToUse'th derivative, for each of the set
	 //of points indicated by array y (with numel elements), putting output
	 //into array arrayOut (with numel elements):
	 /*DLL_EXPORT virtual void MultiPointSpline(value_type* y,
                                              value_type* arrayOut,
                                              int numel,int nArrayToUse);
                                              */

protected:
    DLL_EXPORT int ChooseSpline(value_type& x, value_type* tdata, int Nupper);
    DLL_EXPORT value_type SplineDerivative
    (value_type y,value_type* kdata,
     value_type* ydata,int n)
    {
        //Derivative of the cubic used for interpolation:
        //int n = this->ChooseSpline(y, this->xdata, this->Ndata - 1);
        //const value_type _2p0 = value_type(2.0);
        //const value_type _1p0 = value_type(1.0);
        value_type logfactor = (y != this->_0p0 ? log(y*y / (this->M*this->M))
                          : this->_0p0 );
        value_type t = (logfactor - this->xdata[n])
         / (this->xdata[n + 1] - this->xdata[n]);
        //value_type tp = (this->_1p0)/(this->xdata[n + 1] - this->xdata[n]);

        /*
        value_type a = kdata[n] * (this->xdata[n + 1] - this->xdata[n])
         - (ydata[n + 1] - ydata[n]);
        value_type b = -kdata[n + 1] *(this->xdata[n + 1] - this->xdata[n])
                       + (ydata[n + 1] - ydata[n]);
        return -tp*ydata[n] + tp*ydata[n + 1]
                + tp*(this->_1p0 - t)*(a*(this->_1p0 - t) + b*t)
                - t*tp*(a*(this->_1p0 - t) + b*t)
                + t*(this->_1p0 - t)*(-a*tp + b*tp);
                */
        return this->_6p0*t*(this->_1p0 - t)
        *(ydata[n+1] - ydata[n])/(xdata[n+1] - xdata[n])
        + kdata[n]*(this->_1p0 - t)*(this->_1p0 - this->_3p0*t)
        + kdata[n+1]*t*(this->_3p0*t*t - this->_2p0);

    }
    DLL_EXPORT value_type Spline
    (value_type y,value_type* kdata,
     value_type* ydata,int n)
    {
        //Derivative of the cubic used for interpolation:
        //int n = this->ChooseSpline(y, this->xdata, this->Ndata - 1);
        //const value_type _2p0 = value_type(2.0);
        //const value_type _1p0 = value_type(1.0);
        value_type logfactor = (y != this->_0p0 ? log(y*y / (this->M*this->M))
                          : this->_0p0 );
        value_type t = (logfactor - this->xdata[n])
         / (this->xdata[n + 1] - this->xdata[n]);
         /*
        value_type a = kdata[n] * (this->xdata[n + 1] - this->xdata[n])
         - (ydata[n + 1] - ydata[n]);
        value_type b = -kdata[n + 1] *(this->xdata[n + 1] - this->xdata[n])
                       + (ydata[n + 1] - ydata[n]);
        return (this->_1p0 - t)*ydata[n] + t*ydata[n + 1]
                +t*(this->_1p0 - t)*(a*(this->_1p0 - t) + b*t);
        */
        return (this->_1p0 + t*t*(this->_2p0*t - this->_3p0))*ydata[n]
                +t*t*(this->_3p0 - this->_2p0*t)*ydata[n+1]
                + (kdata[n]*t*(this->_1p0 - t)*(this->_1p0 - t)
                - t*t*(this->_1p0 - t)*kdata[n+1])*(xdata[n+1] - xdata[n]);
    }
    DLL_EXPORT value_type SinglePointSpline
    (value_type y,std::vector<value_type>& kdata,
     std::vector<value_type>& ydata,int n)
    {
        return this->Spline(y,xdata.data(),kdata.data(),
                                       ydata.data(),n);
    }
};




//==============================================================================
//Cubic Spline interpolation of the Standard Model potential, in
//logarithmic space.
template< class value_type >
class SMHiggsPotentialSpline : public SplinePotential< value_type >
{
private:
    //value_type* xdata;//Array of  x positions of
        //the known points for lambda and all its derivatives.
        //Note - use the same grid for ALL derivatives.
        /*
        xdata[i] = x_i grid for lambda
        //Note - x_i here is not phi, but log(phi^2/M^2)
        */
    value_type** ydata;//Array of 5 pointers to the y-values (or rather,
                       //lambda, lambda' etc...  values
        //in this case) at all the points xdata.
        /*
        ydata[0][i] = lambda(xdata[0][i])
        ydata[1][i] = lambda'(xdata[1][i])
        ydata[2][i] = lambda''(xdata[2][i])
        etc...
        */
    value_type** kdata;
        //Derivatives at each point (note that these are NOT the derivatives of
        //lambda directly, but those of the interpolating spline that best fits
        //the points and is C2 continuous). Each derivative of lamda
        //is separately represented as its own spline, to ensure continuity up
        //to fourth order.
        /*
        kdata[0][i] = y'[xdata[0][i]]
        etc...
        */
    //int Ndata;//No of points of data from which to interpolate.
    //value_type M;//Scale appearing in the log(phi/M) terms.
public:
    value_type yscale;//The potential as it stands is formulated in GeV,
    //usually. Set the yscale variable to tell the
    //potential what units we expect supplied arguments y to be in.
    //The potential will then attempt to return V(yscale*y)/yscale^4,
    //V'(yscale*y)/yscale^3
    //V''(yscale*y)/yscale^2, V'''(yscale*y)/yscale and V''''(yscale*y)
    //etc... yscale should be, in GeV, whatever y = 1 corresponds to.
    //The normalisation
    //is because W(y) = V(yscale*y)/yscale^4 is the effective potential
    //appearing in the equations of motion when y = phi/yscale, as opposed
    //to V(phi).
    //Constructors
    SMHiggsPotentialSpline(value_type* XDATA,value_type** YDATA,
                           value_type** KDATA,int NDATA,value_type m)
    : SplinePotential<value_type>(XDATA,NDATA,m)
    {
        //xdata = XDATA;
        ydata = YDATA;
        kdata = KDATA;
        //Ndata = NDATA;
        //M = m;
        yscale = value_type(1.0);
        //Setup maximum and minimum values:
        this->hasMaximum = true;
        this->hasMinimum = true;
        //const value_type _2p0 = value_type(2.0);
        //Maximum and minimum, but note that these
        //should be in the units of which we expect input,
        //Not GeV, as the potential expects:
        this->Maximum = (m/yscale)*exp(XDATA[NDATA - 1]/this->_2p0);
        this->Minimum = (m/yscale)*exp(XDATA[0]/this->_2p0);
    }
    SMHiggsPotentialSpline(value_type* XDATA,value_type** YDATA,
                           value_type** KDATA,int NDATA,value_type m,
                           value_type YSCALE)
    : SplinePotential<value_type>(XDATA,NDATA,m)
    {
        //xdata = XDATA;
        ydata = YDATA;
        kdata = KDATA;
        //Ndata = NDATA;
        //M = m;
        yscale = YSCALE;
        //Setup maximum and minimum values:
        this->hasMaximum = true;
        this->hasMinimum = true;
        //const value_type _2p0 = value_type(2.0);
        this->Maximum = (m/yscale)*exp(XDATA[NDATA - 1]/this->_2p0);
        this->Minimum = (m/yscale)*exp(XDATA[0]/this->_2p0);
    }
    //Potential and its derivatives:
    //Potential itself:
	DLL_EXPORT value_type operator()(value_type y);
	DLL_EXPORT void operator()(value_type* y,value_type* arrayOut,int numel);
    //First derivatives:
	DLL_EXPORT value_type d(value_type y);
	DLL_EXPORT void d(value_type* y,value_type* arrayOut,int numel);
    //Second derivatives:
	DLL_EXPORT value_type d2(value_type y);
	DLL_EXPORT void d2(value_type* y,value_type* arrayOut,int numel);
    //Third derivatives:
	DLL_EXPORT value_type d3(value_type y);
	DLL_EXPORT void d3(value_type* y,value_type* arrayOut,int numel);
    //Fourth derivatives:
	DLL_EXPORT value_type d4(value_type y);
	DLL_EXPORT void d4(value_type* y,value_type* arrayOut,int numel);
	//Function to identify which part of the spline a particular x lies in.
	//This tells us which spline to use, since the interpolating polynomial is
	//contructed piecewise between a given set of points. We identify which x
	//to use by means of a interpolative search, which is O(log(log(n)))
	//complexity if xdata is evenly spaced.
	//DLL_EXPORT int ChooseSpline(value_type& x,value_type* xdata,int Nupper);
	//The splines for each derivative are essentially the same, just applied,
	//to different sets of data, so we implement them all via a single
	//function call:

	//Compute a spline for array (derivative) nArrayToUse, at point y:
	DLL_EXPORT value_type SinglePointSpline(value_type& y,int nArrayToUse);
	DLL_EXPORT value_type SinglePointSpline(value_type& y,int nArrayToUse,
                                            int nSpline);
    //DLL_EXPORT value_type SplineDerivative(value_type& y,int m);
    //DLL_EXPORT value_type SplineDerivative(value_type& y,int m,int n);
	//Compute a spline of the  nArrayToUse'th derivative, for each of the set
	//of points indicated by array y (with numel elements), putting output
	//into array arrayOut (with numel elements):
	DLL_EXPORT void MultiPointSpline(value_type* y,value_type* arrayOut,
                                     int numel,int nArrayToUse);
};

//==============================================================================
#include "rational.h"
//Class which computes couplings as a function of an input
//Einstein-frame field. This only exists so that
//the potential can dual inherit from it, so that it counts
//as an instance of 'coupling_function'. This way dphit_dphi
//and others can refer to it without knowing what the potential is.
/*
template<class value_type>
class coupling_function
{
public:
    //g_n: n specifies which coupling to call,
    //and txi the value of txi at which to evaluate it. m is
    //the derivative that we desire (ie, m = 0 gives the coupling,
    // m = 1 gives its derivative, etc...).
    virtual value_type g_n(int n,int m,value_type txi) = 0;
};
//Function needed to compute the de Sitter potential:
template<class value_type>
class dphit_dphi : public rational_function<value_type>
{
private:
    const value_type _1p0;
    value_type h;//Scale, h = H/Mp where Mp is the planck mass
        //and H is the unit being used for energy.
    coupling_function<value_type>& couplings;
    //List of the relevant value of n needed to use g_n:
    const int[5] = {1,2,3,4,5};
    //{xi,tcl,}

public:
    //constructor:
    dphit_dphi(coupling_function<value_type>& COUPLINGS,value_type H)
    : rational_function<value_type>(), _1p0(1.0), couplings(COUPLINGS)
    {
        h = H;
    }
    value_type A(value_type)
    {
        return sqrt(this->_1p0 - )
    }
    value_type B(value_type)
    {
    }


};
*/
//------------------------------------------------------------------------------
/*
//struct to hold t and couplings, g:
template<class value_type>
struct t_g_pair
{
    value_type t;
    std::vector<value_type>& g;
    t_g_pair(value_type T,std::vector<value_type>& G)
    : t(T),g(G){}
};
//dtdtcl class:
template<class value_type>
class dtdtcl_1loop
: public rational_function_shared_data<value_type,t_g_pair&>
{

};*/
//------------------------------------------------------------------------------
//Abstract class used to transform beta functions from one variabe, t, to
//another variable, tcl, and computing the various derivatives thereof.
template<class value_type>
class dtdtcl_jacobian
{
public:
    virtual value_type dt_dtcl(const value_type& tcl,
                                 const std::vector<value_type>& g,
                                 const std::vector<value_type>& dgdt) = 0;
    virtual value_type d2t_dtcl2(const value_type& tcl,
                                 const std::vector<value_type>& g,
                                 const std::vector<value_type>& dgdt,
                                 const std::vector<value_type>& d2gdt2,
                                 const value_type& dtdtcl) = 0;
    virtual value_type d3t_dtcl3(const value_type& tcl,
                                 const std::vector<value_type>& g,
                                 const std::vector<value_type>& dgdt,
                                 const std::vector<value_type>& d2gdt2,
                                 const std::vector<value_type>& d3gdt3,
                                 const value_type& dtdtcl,
                                 const value_type& d2tdtcl2) = 0;
    virtual value_type d4t_dtcl4(const value_type& tcl,
                                 const std::vector<value_type>& g,
                                 const std::vector<value_type>& dgdt,
                                 const std::vector<value_type>& d2gdt2,
                                 const std::vector<value_type>& d3gdt3,
                                 const std::vector<value_type>& d4gdt4,
                                 const value_type& dtdtcl,
                                 const value_type& d2tdtcl2,
                                 const value_type& d3tdtcl3) = 0;
    virtual value_type d5t_dtcl5(const value_type& tcl,
                                 const std::vector<value_type>& g,
                                 const std::vector<value_type>& dgdt,
                                 const std::vector<value_type>& d2gdt2,
                                 const std::vector<value_type>& d3gdt3,
                                 const std::vector<value_type>& d4gdt4,
                                 const std::vector<value_type>& d5gdt5,
                                 const value_type& dtdtcl,
                                 const value_type& d2tdtcl2,
                                 const value_type& d3tdtcl3,
                                 const value_type& d4tdtcl4) = 0;
};
//------------------------------------------------------------------------------

//Function to evaluate a coupling given a coupling array.
template<class value_type>
class couplingEvaluator
{
public:

    value_type l(const std::vector<value_type>& g) {return g[0];}
    value_type g12(const std::vector<value_type>& g) {return g[1];}
    value_type g22(const std::vector<value_type>& g) {return g[2];}
    value_type g32(const std::vector<value_type>& g) {return g[3];}
    value_type m2(const std::vector<value_type>& g) {return g[4];}
    value_type ye2(const std::vector<value_type>& g) {return g[5];}
    value_type ymu2(const std::vector<value_type>& g) {return g[6];}
    value_type ytau2(const std::vector<value_type>& g) {return g[7];}
    value_type yu2(const std::vector<value_type>& g) {return g[8];}
    value_type yd2(const std::vector<value_type>& g) {return g[9];}
    value_type yc2(const std::vector<value_type>& g) {return g[10];}
    value_type ys2(const std::vector<value_type>& g) {return g[11];}
    value_type yt2(const std::vector<value_type>& g) {return g[12];}
    value_type yb2(const std::vector<value_type>& g) {return g[13];}
    value_type Z(const std::vector<value_type>& g) {return g[14];}
    value_type tcl(const std::vector<value_type>& g) {return g[15];}
    value_type t(const std::vector<value_type>& g) {return g[16];}
    value_type xi(const std::vector<value_type>& g) {return g[17];}
    value_type V0(const std::vector<value_type>& g) {return g[18];}
    value_type kappa(const std::vector<value_type>& g) {return g[19];}
    value_type alpha1(const std::vector<value_type>& g) {return g[20];}
    value_type alpha2(const std::vector<value_type>& g) {return g[21];}
    value_type alpha3(const std::vector<value_type>& g) {return g[22];}
};

//------------------------------------------------------------------------------
//Basic Jacobian transformation - simpler to code!
template<class value_type>
class dtdtcl_dS : public dtdtcl_jacobian<value_type>
{
public:
    //Constants we can vary to change the specific scale choice:
    value_type alpha;
    value_type beta;
    const value_type R;
    const value_type mt;

    couplingEvaluator<value_type> eval;

    //Numerical constants:
    const value_type _1p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _4p0;
    const value_type _5p0;
    const value_type _6p0;
    const value_type _24p0;
    const value_type _12p0;
    const value_type _36p0;
    const value_type _8p0;
    const value_type _9p0;
    const value_type _10p0;

    //Constructor
    dtdtcl_dS(value_type ALPHA,value_type BETA,value_type r,
              value_type MT)
    :_1p0(1), _2p0(2), _3p0(3), _4p0(4),_5p0(5.0), _6p0(6),
    _24p0(24.0),_12p0(12.0),_36p0(36.0),_8p0(8.0),_9p0(9.0),
    _10p0(10.0),R(r),mt(MT)
    {alpha = ALPHA;beta = BETA;}


    value_type dt_dtcl(const value_type& tcl,
                       const std::vector<value_type>& g,
                       const std::vector<value_type>& dgdt)
    {
        value_type phicl = this->mt*exp(tcl/this->_2p0);
        value_type a = this->A(tcl,g,dgdt,phicl);
        value_type b = this->B(tcl,g,phicl);
        return a/b;
    }
    value_type d2t_dtcl2(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const value_type& dtdtcl)
    {
        value_type phicl = this->mt*exp(tcl/this->_2p0);
        value_type a = this->A(tcl,g,dgdt,phicl);
        value_type b = this->B(tcl,g,phicl);
        value_type ap = this->Ap(tcl,g,dgdt,d2gdt2,dtdtcl,phicl);
        value_type bp = this->Bp(tcl,g,dgdt,dtdtcl,phicl);
        return ap/b - a*bp/(b*b);
    }
    value_type d3t_dtcl3(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const std::vector<value_type>& d3gdt3,
                         const value_type& dtdtcl,
                         const value_type& d2tdtcl2)
    {
        value_type phicl = this->mt*exp(tcl/this->_2p0);
        value_type a = this->A(tcl,g,dgdt,phicl);
        value_type b = this->B(tcl,g,phicl);
        value_type ap = this->Ap(tcl,g,dgdt,d2gdt2,dtdtcl,phicl);
        value_type bp = this->Bp(tcl,g,dgdt,dtdtcl,phicl);
        value_type app = this->App(tcl,g,dgdt,d2gdt2,d3gdt3,dtdtcl,d2tdtcl2,
                                   phicl);
        value_type bpp = this->Bpp(tcl,g,dgdt,d2gdt2,dtdtcl,d2tdtcl2,phicl);
        return app/b - this->_2p0*ap*bp/(b*b) - a*bpp/(b*b)
                + this->_2p0*a*bp*bp/(b*b*b);
    }
    value_type d4t_dtcl4(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const std::vector<value_type>& d3gdt3,
                         const std::vector<value_type>& d4gdt4,
                         const value_type& dtdtcl,
                         const value_type& d2tdtcl2,
                         const value_type& d3tdtcl3)
    {
        value_type phicl = this->mt*exp(tcl/this->_2p0);
        value_type a = this->A(tcl,g,dgdt,phicl);
        value_type b = this->B(tcl,g,phicl);
        value_type ap = this->Ap(tcl,g,dgdt,d2gdt2,dtdtcl,phicl);
        value_type bp = this->Bp(tcl,g,dgdt,dtdtcl,phicl);
        value_type app = this->App(tcl,g,dgdt,d2gdt2,d3gdt3,dtdtcl,d2tdtcl2,
                                   phicl);
        value_type bpp = this->Bpp(tcl,g,dgdt,d2gdt2,dtdtcl,d2tdtcl2,phicl);
        value_type appp = this->Appp(tcl,g,dgdt,d2gdt2,d3gdt3,d4gdt4,
                                     dtdtcl,d2tdtcl2,d3tdtcl3,phicl);
        value_type bppp = this->Bppp(tcl,g,dgdt,d2gdt2,d3gdt3,
                                     dtdtcl,d2tdtcl2,d3tdtcl3,phicl);
        return appp/b - this->_3p0*app*bp/(b*b)- this->_3p0*ap*bpp/(b*b)
                + this->_6p0*ap*bp*bp/(b*b*b) - a*bppp/(b*b)
                + this->_6p0*a*bpp*bp/(b*b*b)
                - this->_6p0*a*bp*bp*bp/(b*b*b*b);

    }
    value_type d5t_dtcl5(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const std::vector<value_type>& d3gdt3,
                         const std::vector<value_type>& d4gdt4,
                         const std::vector<value_type>& d5gdt5,
                         const value_type& dtdtcl,
                         const value_type& d2tdtcl2,
                         const value_type& d3tdtcl3,
                         const value_type& d4tdtcl4)
    {
        value_type phicl = this->mt*exp(tcl/this->_2p0);
        value_type a = this->A(tcl,g,dgdt,phicl);
        value_type b = this->B(tcl,g,phicl);
        value_type ap = this->Ap(tcl,g,dgdt,d2gdt2,dtdtcl,phicl);
        value_type bp = this->Bp(tcl,g,dgdt,dtdtcl,phicl);
        value_type app = this->App(tcl,g,dgdt,d2gdt2,d3gdt3,dtdtcl,d2tdtcl2,
                                   phicl);
        value_type bpp = this->Bpp(tcl,g,dgdt,d2gdt2,dtdtcl,d2tdtcl2,phicl);
        value_type appp = this->Appp(tcl,g,dgdt,d2gdt2,d3gdt3,d4gdt4,
                                     dtdtcl,d2tdtcl2,d3tdtcl3,phicl);
        value_type bppp = this->Bppp(tcl,g,dgdt,d2gdt2,d3gdt3,
                                     dtdtcl,d2tdtcl2,d3tdtcl3,phicl);
        value_type apppp = this->Apppp(tcl,g,dgdt,d2gdt2,d3gdt3,d4gdt4,d5gdt5,
                                       dtdtcl,d2tdtcl2,d3tdtcl3,d4tdtcl4,phicl);
        value_type bpppp = this->Bpppp(tcl,g,dgdt,d2gdt2,d3gdt3,d4gdt4,
                                       dtdtcl,d2tdtcl2,d3tdtcl3,d4tdtcl4,phicl);
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

private:
    value_type A(const value_type& tcl,const std::vector<value_type>& g,
                 const std::vector<value_type>& dgdt,const value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        value_type Zp = this->eval.Z(dgdt);
        return (this->alpha*(phicl*phicl*(Z*(Z + Zp))));
    }
    value_type B(const value_type& tcl,const std::vector<value_type>& g,
                 value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        return ((this->beta*this->R) + (this->alpha*(phicl*phicl*Z*Z)));
    }
    value_type Ap(const value_type& tcl,const std::vector<value_type>& g,
                 const std::vector<value_type>& dgdt,
                 const std::vector<value_type>& d2gdt2,
                 const value_type& dtdtcl,const value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        value_type Zp = this->eval.Z(dgdt);
        value_type Zpp = this->eval.Z(d2gdt2);

        return (this->alpha*(phicl*phicl*(Z*Z + ((dtdtcl*Zp*Zp)
                + (Z*(((this->_1p0 + (this->_2p0*dtdtcl))*Zp)
                + (dtdtcl*Zpp)))))));
    }
    value_type Bp(const value_type& tcl,const std::vector<value_type>& g,
                 const std::vector<value_type>& dgdt,
                 const value_type& dtdtcl,const value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        value_type Zp = this->eval.Z(dgdt);

        return (alpha*(phicl*phicl*(Z*(Z + (this->_2p0*(dtdtcl*Zp))))));
    }
    value_type App(const value_type& tcl,const std::vector<value_type>& g,
                   const std::vector<value_type>& dgdt,
                   const std::vector<value_type>& d2gdt2,
                   const std::vector<value_type>& d3gdt3,
                   const value_type& dtdtcl,const value_type& d2tdtcl2,
                   const value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        value_type Zp = this->eval.Z(dgdt);
        value_type Zpp = this->eval.Z(d2gdt2);
        value_type Zppp = this->eval.Z(d3gdt3);

        return (this->alpha*(phicl*phicl*(Z*Z + ((Zp*((d2tdtcl2*Zp)
                + ((this->_2p0*(dtdtcl*Zp)) + (dtdtcl*dtdtcl*((this->_2p0*Zp)
                + (this->_3p0*Zpp)))))) + (Z*(((this->_1p0
                + ((this->_2p0*d2tdtcl2) + (this->_4p0*dtdtcl)))*Zp)
                + ((d2tdtcl2*Zpp) + ((this->_2p0*(dtdtcl*Zpp))
                + (dtdtcl*dtdtcl*((this->_2p0*Zpp) + Zppp))))))))));
    }
    value_type Bpp(const value_type& tcl,const std::vector<value_type>& g,
                   const std::vector<value_type>& dgdt,
                   const std::vector<value_type>& d2gdt2,
                   const value_type& dtdtcl,const value_type& d2tdtcl2,
                   const value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        value_type Zp = this->eval.Z(dgdt);
        value_type Zpp = this->eval.Z(d2gdt2);

        return (alpha*(phicl*phicl*(Z*Z + ((this->_2p0*(dtdtcl*dtdtcl*Zp*Zp))
                + (this->_2p0*(Z*((d2tdtcl2*Zp) + ((this->_2p0*(dtdtcl*Zp))
                + (dtdtcl*dtdtcl*Zpp)))))))));
    }
    value_type Appp(const value_type& tcl,const std::vector<value_type>& g,
                   const std::vector<value_type>& dgdt,
                   const std::vector<value_type>& d2gdt2,
                   const std::vector<value_type>& d3gdt3,
                   const std::vector<value_type>& d4gdt4,
                   const value_type& dtdtcl,const value_type& d2tdtcl2,
                   const value_type& d3tdtcl3,const value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        value_type Zp = this->eval.Z(dgdt);
        value_type Zpp = this->eval.Z(d2gdt2);
        value_type Zppp = this->eval.Z(d3gdt3);
        value_type Zpppp = this->eval.Z(d4gdt4);

        return (this->alpha*(phicl*phicl*(Z*Z + ((((this->_3p0*d2tdtcl2)
                + d3tdtcl3)*Zp*Zp) + ((this->_3p0*(dtdtcl*dtdtcl*(Zp
                *((this->_2p0*Zp) + (this->_3p0*Zpp)))))
                + ((this->_3p0*(dtdtcl*(Zp*(((this->_1p0
                + (this->_2p0*d2tdtcl2))*Zp)
                + (this->_3p0*(d2tdtcl2*Zpp))))))
                + ((dtdtcl*dtdtcl*dtdtcl*((this->_3p0*Zpp*Zpp)
                + (Zp*((this->_6p0*Zpp) + (this->_4p0*Zppp)))))
                + (Z*(((this->_1p0 + ((this->_6p0*d2tdtcl2)
                + ((this->_2p0*d3tdtcl3) + (this->_6p0*dtdtcl))))*Zp)
                + ((((this->_3p0*d2tdtcl2) + d3tdtcl3)*Zpp)
                + ((this->_3p0*(dtdtcl*dtdtcl*((this->_2p0*Zpp) + Zppp)))
                + ((this->_3p0*(dtdtcl*(((this->_1p0 + (this->_2p0*d2tdtcl2))
                *Zpp) + (d2tdtcl2*Zppp)))) + (dtdtcl*dtdtcl*dtdtcl
                *((this->_2p0*Zppp) + Zpppp))))))))))))));
    }
    value_type Bppp(const value_type& tcl,const std::vector<value_type>& g,
                   const std::vector<value_type>& dgdt,
                   const std::vector<value_type>& d2gdt2,
                   const std::vector<value_type>& d3gdt3,
                   const value_type& dtdtcl,const value_type& d2tdtcl2,
                   const value_type& d3tdtcl3,const value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        value_type Zp = this->eval.Z(dgdt);
        value_type Zpp = this->eval.Z(d2gdt2);
        value_type Zppp = this->eval.Z(d3gdt3);

        return (this->alpha*(phicl*phicl*(Z*Z + ((this->_6p0*(dtdtcl*(Zp
                *((d2tdtcl2*Zp) + ((dtdtcl*Zp) + (dtdtcl*dtdtcl*Zpp))))))
                + (this->_2p0*(Z*((((this->_3p0*d2tdtcl2) + d3tdtcl3)*Zp)
                + ((this->_3p0*(dtdtcl*dtdtcl*Zpp))
                + ((this->_3p0*(dtdtcl*(Zp + (d2tdtcl2*Zpp))))
                   + (dtdtcl*dtdtcl*dtdtcl*Zppp))))))))));
    }
    value_type Apppp(const value_type& tcl,const std::vector<value_type>& g,
                   const std::vector<value_type>& dgdt,
                   const std::vector<value_type>& d2gdt2,
                   const std::vector<value_type>& d3gdt3,
                   const std::vector<value_type>& d4gdt4,
                   const std::vector<value_type>& d5gdt5,
                   const value_type& dtdtcl,const value_type& d2tdtcl2,
                   const value_type& d3tdtcl3,const value_type& d4tdtcl4,
                   const value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        value_type Zp = this->eval.Z(dgdt);
        value_type Zpp = this->eval.Z(d2gdt2);
        value_type Zppp = this->eval.Z(d3gdt3);
        value_type Zpppp = this->eval.Z(d4gdt4);
        value_type Zppppp = this->eval.Z(d5gdt5);

        return (this->alpha*(phicl*phicl*(Z*Z + ((Zp*((((this->_6p0*d2tdtcl2)
                + ((this->_6p0*d2tdtcl2*d2tdtcl2) + ((this->_4p0*d3tdtcl3)
                + d4tdtcl4)))*Zp) + (this->_9p0*(d2tdtcl2*d2tdtcl2*Zpp))))
                + ((this->_4p0*(dtdtcl*(Zp*(((this->_1p0
                + ((this->_6p0*d2tdtcl2) + (this->_2p0*d3tdtcl3)))*Zp)
                + (this->_3p0*(((this->_3p0*d2tdtcl2) + d3tdtcl3)*Zpp))))))
                + ((this->_4p0*(dtdtcl*dtdtcl*dtdtcl*((this->_3p0*Zpp*Zpp)
                + (Zp*((this->_6p0*Zpp) + (this->_4p0*Zppp))))))
                + ((this->_6p0*(dtdtcl*dtdtcl*((this->_2p0*Zp*Zp)
                + ((this->_3p0*(d2tdtcl2*Zpp*Zpp)) + (Zp*(((this->_3p0
                + (this->_6p0*d2tdtcl2))*Zpp)
                + (this->_4p0*(d2tdtcl2*Zppp))))))))
                + ((dtdtcl*dtdtcl*dtdtcl*dtdtcl*((this->_6p0*Zpp*Zpp)
                + ((this->_10p0*(Zpp*Zppp)) + (Zp*((this->_8p0*Zppp)
                + (this->_5p0*Zpppp)))))) + (Z*(((this->_1p0
                + ((this->_12p0*d2tdtcl2) + ((this->_8p0*d3tdtcl3)
                + ((this->_2p0*d4tdtcl4) + (this->_8p0*dtdtcl)))))*Zp)
                + ((this->_6p0*(d2tdtcl2*Zpp))
                + ((this->_6p0*(d2tdtcl2*d2tdtcl2*Zpp))
                + ((this->_4p0*(d3tdtcl3*Zpp)) + ((d4tdtcl4*Zpp)
                + ((this->_3p0*(d2tdtcl2*d2tdtcl2*Zppp))
                + ((this->_4p0*(dtdtcl*(((this->_1p0 + ((this->_6p0*d2tdtcl2)
                + (this->_2p0*d3tdtcl3)))*Zpp) + (((this->_3p0*d2tdtcl2)
                + d3tdtcl3)*Zppp))))
                + ((this->_4p0*(dtdtcl*dtdtcl*dtdtcl*((this->_2p0*Zppp)
                + Zpppp))) + ((this->_6p0*(dtdtcl*dtdtcl*((this->_2p0*Zpp)
                + (((this->_1p0 + (this->_2p0*d2tdtcl2))*Zppp)
                + (d2tdtcl2*Zpppp)))))
                + (dtdtcl*dtdtcl*dtdtcl*dtdtcl*((this->_2p0*Zpppp)
                + Zppppp))))))))))))))))))));
    }
    value_type Bpppp(const value_type& tcl,const std::vector<value_type>& g,
                     const std::vector<value_type>& dgdt,
                     const std::vector<value_type>& d2gdt2,
                     const std::vector<value_type>& d3gdt3,
                     const std::vector<value_type>& d4gdt4,
                     const value_type& dtdtcl,const value_type& d2tdtcl2,
                     const value_type& d3tdtcl3,const value_type& d4tdtcl4,
                     const value_type& phicl)
    {
        value_type Z = this->eval.Z(g);
        value_type Zp = this->eval.Z(dgdt);
        value_type Zpp = this->eval.Z(d2gdt2);
        value_type Zppp = this->eval.Z(d3gdt3);
        value_type Zpppp = this->eval.Z(d4gdt4);

        return (this->alpha*(phicl*phicl*(Z*Z
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
                + (dtdtcl*dtdtcl*dtdtcl*dtdtcl*Zpppp))))))))))))))));
    }
};
//------------------------------------------------------------------------------
template<class value_type>
void fill_constant_arrays(std::vector<value_type>& nilist,
                              std::vector<value_type>& cilist)
{
    const value_type _1p0 = value_type(1.0);
    const value_type _2p0 = value_type(2.0);
    const value_type _6p0 = value_type(6.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _12p0 = value_type(12.0);
    const value_type _5p0 = value_type(5.0);
    value_type nilist_array[] = {_2p0,_6p0,-_2p0,
                                 _1p0,_3p0,-_1p0,
                                 -_12p0,_1p0,_3p0};
    value_type cilist_array[] = {_3p0/_2p0,
                                 _5p0/_6p0,
                                _3p0/_2p0,_3p0/_2p0,
                                _5p0/_6p0,_3p0/_2p0,
                                _3p0/_2p0,_3p0/_2p0,
                                _3p0/_2p0};
    nilist.assign(nilist_array,nilist_array + 9);
    cilist.assign(cilist_array,cilist_array + 9);
}
template<class value_type>
void fill_constant_arrays(std::vector<value_type>& nilist,
                              std::vector<value_type>& cilist,
                              std::vector<value_type>& npilist)
{
    const value_type _1p0 = value_type(1.0);
    const value_type _2p0 = value_type(2.0);
    const value_type _6p0 = value_type(6.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _12p0 = value_type(12.0);
    const value_type _5p0 = value_type(5.0);
    const value_type _34p0 = value_type(34.0);
    const value_type _15p0 = value_type(15.0);
    const value_type _4p0 = value_type(4.0);
    const value_type _17p0 = value_type(17.0);
    const value_type _38p0 = value_type(38.0);
    const value_type _8p0 = value_type(8.0);
    const value_type _24p0 = value_type(24.0);
    const value_type _136p0 = value_type(136.0);
    const value_type _16p0 = value_type(16.0);
    const value_type _19p0 = value_type(19.0);
    value_type nilist_array[] = {_2p0,_6p0,-_2p0,
                                 _1p0,_3p0,-_1p0,
                                 -_12p0,-_12p0,-_12p0,
                                 -_12p0,-_12p0,-_12p0,
                                 -_4p0,-_4p0,-_4p0,
                                 _1p0,_2p0,_1p0,-_2p0,-_1p0,
                                 _1p0,_3p0,-_1p0,
                                 _8p0,_24p0,-_8p0,
                                 -_2p0,-_2p0,-_2p0,
                                 -_1p0,-_8p0};
    value_type cilist_array[] = {_3p0/_2p0,_5p0/_6p0,_3p0/_2p0,
                                 _3p0/_2p0,_5p0/_6p0,_3p0/_2p0,
                                _3p0/_2p0,_3p0/_2p0,_3p0/_2p0,
                                _3p0/_2p0,_3p0/_2p0,_3p0/_2p0,
                                _3p0/_2p0,_3p0/_2p0,_3p0/_2p0,
                                _3p0/_2p0,
                                _3p0/_2p0,_3p0/_2p0,
                                _3p0/_2p0,_3p0/_2p0,
                                _3p0/_2p0,_5p0/_6p0,_3p0/_2p0,
                                _3p0/_2p0,_5p0/_6p0,_3p0/_2p0,
                                _3p0/_2p0,_3p0/_2p0,_3p0/_2p0,
                                _3p0/_2p0,_3p0/_2p0};
    value_type npilist_array[] = {-_34p0/_15p0,-_34p0/_5p0,
                                  _4p0/_15p0,-_17p0/_15p0,
                                  -_17p0/_5p0,_2p0/_15p0,
                                  _38p0/_5p0,_38p0/_5p0,_38p0/_5p0,
                                  _38p0/_5p0,_38p0/_5p0,_38p0/_5p0,
                                  _38p0/_15p0,_38p0/_15p0,_38p0/_15p0,
                                  -_2p0/_15p0,
                                  -_4p0/_15p0,-_2p0/_15p0,
                                  _4p0/_15p0,_2p0/_15p0,
                                  -_17p0/_15p0,-_17p0/_5p0,_2p0/_15p0,
                                  -_136p0/_15p0,-_136p0/_5p0,_16p0/_15p0,
                                  _19p0/_15p0,_19p0/_15p0,_19p0/_15p0,
                                  _2p0/_15p0,_16p0/_15p0};
    nilist.assign(nilist_array,nilist_array + 31);
    cilist.assign(cilist_array,cilist_array + 31);
    npilist.assign(npilist_array,npilist_array + 31);
}
//Class for populating SM coefficient arrays
template<class value_type>
class coefficientArray
{
private:

public:

    const int nLoops;

    //Numerical constants:
    const value_type _3p0;
    const value_type _5p0;
    const value_type _4p0;
    const value_type _2p0;
    const value_type _0p0;
    const value_type _1p0;
    const value_type _12p0;
    const value_type _6p0;
    value_type su5;

    //Constructor:
    coefficientArray(int LOOPS)
    : nLoops(LOOPS), _3p0(3.0), _5p0(5.0),_4p0(4.0),_2p0(2.0),_0p0(0.0),
    _1p0(1.0),_12p0(12.0),_6p0(6.0)
    {
        this->su5 = this->_3p0/this->_5p0;
    }



    virtual void populateCoefficientArrays
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& g) = 0;

    virtual void populateCoefficientArraysDerivatives
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& gp) = 0;
};



template<class value_type>
class coefficientArrayNoGrav : public coefficientArray<value_type>
{
private:

public:
    coefficientArrayNoGrav() : coefficientArray(9){}

    void populateCoefficientArrays
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& g)
     {
         value_type kiarray[] = {g[2]/this->_4p0,g[2]/this->_4p0,
                                 g[2]/this->_4p0,
                                (g[2] + this->su5*g[1])/this->_4p0,
                                (g[2] + this->su5*g[1])/this->_4p0,
                                (g[2] + this->su5*g[1])/this->_4p0,
                                g[4]/this->_2p0,
                                this->_3p0*g[0],g[0]};
        value_type kpiarray[] = {this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,g[5],g[5]};
        value_type thetaiarray[] = {this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    -this->_1p0/this->_6p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    -this->_1p0/this->_6p0,
                                    this->_1p0/this->_12p0,
                                    g[6] - this->_1p0/this->_6p0,
                                    g[6] - this->_1p0/this->_6p0};


        k.assign(kiarray,kiarray + this->nLoops);
        kp.assign(kpiarray,kpiarray + this->nLoops);
        theta.assign(thetaiarray,thetaiarray + this->nLoops);
     }

     void populateCoefficientArraysDerivatives
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& gp)
     {
        value_type kiarray[] = {gp[2]/this->_4p0,gp[2]/this->_4p0,
                                    gp[2]/this->_4p0,
                                    (gp[2] + su5*gp[1])/this->_4p0,
                                    (gp[2] + this->su5*gp[1])/this->_4p0,
                                    (gp[2] + this->su5*gp[1])/this->_4p0,
                                    gp[4]/this->_2p0,
                                    this->_3p0*gp[0],gp[0]};
            value_type kpiarray[] = {this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0,gp[5],gp[5]};
            value_type thetaiarray[] = {this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,
                                    this->_0p0,gp[6],gp[6]};


            k.assign(kiarray,kiarray + this->nLoops);
            kp.assign(kpiarray,kpiarray + this->nLoops);
            theta.assign(thetaiarray,thetaiarray + this->nLoops);
     }
};



template<class value_type>
class coefficientArrayWithGrav : public coefficientArray<value_type>,
public couplingEvaluator<value_type>
{
public:
    const value_type zetaz;
    const value_type zetaw;
    coefficientArrayWithGrav(value_type ZETAZ,value_type ZETAW)
     : coefficientArray(31), zetaz(ZETAZ),zetaw(ZETAW){}

    void populateCoefficientArrays
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& g)
     {
        value_type Mz2 = (g22(g) + this->su5*g12(g))/this->_4p0;
        value_type Mw2 = g22(g)/this->_4p0;
        value_type kiarray[] = {Mw2,Mw2,Mw2,Mz2,Mz2,Mz2,
                                yu2(g)/this->_2p0,yd2(g)/this->_2p0,
                                yc2(g)/this->_2p0,ys2(g)/this->_2p0,
                                yt2(g)/this->_2p0,yb2(g)/this->_2p0,
                                ye2(g)/this->_2p0,ymu2(g)/this->_2p0,
                                ytau2(g)/this->_2p0,this->_3p0*l(g),
                                l(g) + this->zetaw*Mw2,
                                l(g) + this->zetaz*Mz2,
                                this->zetaw*Mw2,this->zetaz*Mz2,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0};
        value_type kpiarray[] = {this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                m2(g),m2(g),m2(g),this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0};
        value_type thetaiarray[] = {this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    -this->_1p0/this->_6p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    -this->_1p0/this->_6p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    xi(g) - this->_1p0/this->_6p0,
                                    xi(g) - this->_1p0/this->_6p0,
                                    xi(g) - this->_1p0/this->_6p0,
                                    -this->_1p0/this->_6p0,
                                    -this->_1p0/this->_6p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    -this->_1p0/this->_6p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    -this->_1p0/this->_6p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                    -this->_1p0/this->_6p0,
                                    -this->_1p0/this->_6p0};

        k.assign(kiarray,kiarray + this->nLoops);
        kp.assign(kpiarray,kpiarray + this->nLoops);
        theta.assign(thetaiarray,thetaiarray + this->nLoops);
     }

     void populateCoefficientArraysDerivatives
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& gp)
     {
        value_type Mz2p = (g22(gp) + this->su5*g12(gp))/this->_4p0;
        value_type Mw2p = g22(gp)/this->_4p0;
        value_type kiarray[] = {Mw2p,Mw2p,Mw2p,Mz2p,Mz2p,Mz2p,
                                yu2(gp)/this->_2p0,yd2(gp)/this->_2p0,
                                yc2(gp)/this->_2p0,ys2(gp)/this->_2p0,
                                yt2(gp)/this->_2p0,yb2(gp)/this->_2p0,
                                ye2(gp)/this->_2p0,ymu2(gp)/this->_2p0,
                                ytau2(gp)/this->_2p0,this->_3p0*l(gp),
                                l(gp) + this->zetaw*Mw2p,
                                l(gp) + this->zetaz*Mz2p,
                                this->zetaw*Mw2p,this->zetaz*Mz2p,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0};
        value_type kpiarray[] = {this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                m2(gp),m2(gp),m2(gp),this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0};
        value_type thetaiarray[] = {this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,this->_0p0,
                                    xi(gp),xi(gp),xi(gp),
                                    this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0,this->_0p0,this->_0p0,
                                    this->_0p0};

        k.assign(kiarray,kiarray + this->nLoops);
        kp.assign(kpiarray,kpiarray + this->nLoops);
        theta.assign(thetaiarray,thetaiarray + this->nLoops);
     }
};

//Jacobian for the transformation where we set the sum of the loop corrections
//to zero:
template<class value_type>
class dtdtcl_1loop : public dtdtcl_jacobian<value_type>
{
    //Reference list of the couplings vector and which
    //SM couplings they refer to:
    /*
    g[0] - lambda (Higgs self coupling)
    g[1] - g1^2 (U(1) coupling)
    g[2] - g2^2 (SU(2) coupling)
    g[3] - g3^2 (SU(3) coupling)
    g[4] - yt^2 (Top quark Yukawa coupling)
    g[5] - m^2 (Higgs tachyonic mass term)
    g[6] - xi (Higgs-curvature coupling)
    g[7] - Z (Higgs field renormalisation) - technically the RHS of this
                  is a gamma function, but we solve it with the beta functions)
    //Additionally, we have to solve for the non-beta functions:
    g[8] = t_cl = log(phicl^2/mt^2)
    g[9] = t = log(mu^2/mt^2)
    //Cosmological constant, if we include it:
    g[10] = Lambda*Mp^2 (cosmological constant, in units of the Planck
                                mass)
    g[11] = kappa (gravitational constant)
    g[12] = alpha1 (R^2 term)
    g[13] = alpha2 (R_{\mu\nu}R^{\mu\nu} term)
    g[14] = alpha3 (R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} term)

    */

    //Hard coded numerical constants:
    const value_type _0p0,_1p0,_16p0,_3p0,_8p0,_9p0,_m3p0,_2p0,_m9p0,_12p0;
    const value_type _6p0,_24p0,_m6p0,_4p0,_27p0,_m27p0,_72p0,_18p0,_m18p0;
    const value_type _96p0,_36p0,_m24p0,_25p0,_2304p0,_10p0,_5p0,_11p0,_45p0;
    const value_type _1152p0,_384p0,_288p0,_256p0,_m43p0,_128p0,_32p0,_m11p0;
    const value_type _m17p0,_m8p0;
    const value_type _m4p0,_20p0,_m1p0,_m2p0,_m12p0,_m16p0;

    //Physical constants:
    const value_type mt;//Top quark mass (used as the base scale wrt which mu
        //in the RG equations is defined, ie, t = log(mu^2/mt^2) is the
        //dimensionless running parameter.
    const value_type R;//Ricci scalar.
    //Coefficients appearing in the 1-loop corrections to the SM effective
        //potential:
    std::vector<value_type> nilist;
    std::vector<value_type> cilist;
    std::vector<value_type> npilist;
    //Function to fill these coefficients when the constructor is called
    //(separated into a function of its own right in case we want to add
    //multiple constructors at some point).

    //Tells the object how many loop correction terms are used in the
    //SM potential:
    static const int nLoops = 31;






public:

    coefficientArrayWithGrav<value_type> coeff;
    couplingEvaluator<value_type> eval;

    //Constructor:
    dtdtcl_1loop(value_type Mt,value_type r,value_type zetaz,value_type zetaw)
    : _1p0(1), _16p0(16), _3p0(3), _8p0(8), _9p0(9), _m3p0(-3), _2p0(2),
    _m9p0(-9),
    _12p0(12), _6p0(6), _24p0(24), _m6p0(-6), _4p0(4), _27p0(27), _m27p0(-27),
    _72p0(72), _18p0(18), _m18p0(-18), _96p0(96), _36p0(36), _m24p0(-24),
    _25p0(25), _2304p0(2304), _10p0(10), _5p0(5), _11p0(11), _45p0(45),
    _1152p0(1152), _384p0(384), _288p0(288), _256p0(256), _m43p0(-43),
    _128p0(128),_m4p0(-4.0),_20p0(20.0),_m1p0(-1.0),_m2p0(-2.0),_m12p0(-12.0),
    _32p0(32), _m11p0(-11), _m17p0(-17), _m8p0(-8), _0p0(0.0),_m16p0(-16.0),
    mt(Mt), R(r), coeff(zetaz,zetaw)
    {
        fill_constant_arrays(this->nilist,this->cilist,this->npilist);
    }
    //Parameters to the various functions used in this class:
    //
    //logfactor - log(Mi2/mu^2) (vector for different values of i)
    //phicl - classical field, ie, field evaluated at Electroweak scale.
    //          phicl mt*exp(tcl/2)
    //Mi2 - Mi2(phi) = k_i*phi^2(t) - k'_i + theta_i*R
    //      Where phi = Z(t)*phicl, R = Ricci scalar, and k_i, k'_i, theta_i
    //      are coupling dependent coefficients (depending on different
    //      couplings for each i).
    //kiarray - k_i coefficient array (see populateCoefficientArrays for details).
    //pMi2pt - partial derivative wrt t of Mi2, at constant phicl.
    //          pMi2pt = (dk_i/dt)*Z^2(t)*phicl^2 + k_i*2*Z(t)*phicl^2
    //                   -(dk'_i/dt) + (dtheta_i/dt)*R
    //kiparray - derivatives of kiarray wrt t (not tcl).
    //kipparray - higher t derivatives of k_i (# of ps = # of derivatives)
    //Mi2parray - derivatives of kiarray wrt tcl (not t!)
    //Mi2pparray - higher tcl derivatives of Mi2 (# of ps = # of derivatives)
    //dtdtcl - 1st derivative of t wrt tcl.
    //d2tdtcl2 - 2nd derivative of t wrt tcl.
    //...
    //dntdtcln - nth derivative of t wrt tcl.
    //Z - field renormalisation, defined as phi(t) = Z(t)*phicl, where phi(t)
    //      is the renormalised field.
    //Zd, Zdd... - t derivatives of Z (# of ds = # of derivatives).
    //pMi2ptparray, pMi2ptpparray... - Derivatives of pMi2pt wrt tcl (not t)
    //mu - renormalisation scale.
    //thetai - theta_i coefficient found in Mi2.
    //dthetaidt, d2thetaidt2 ... - derivatives (1st, 2nd...) of theta_i wrt t.
    //ki - k_i coefficient found in Mi2 (alias for kiarray - both are sometimes
    //      used in the same function because an array is used to define the
    //      coefficients, and a vector to store them, so two names are used)
    //kpi - k'_i coefficients found in Mi2
    //dkidt, d2kidt2,... - t derivatives (1st, 2nd...) of k_i
    //dkpidt, d2kpidt2,... - t derivatives (1st, 2nd...) of k'_i
    //logmu2 = log(mu^2/mt^2)
    //t - renormalisation running parameter, t = log(mu^2/mt^2)
    //g - couplings vector
    //dgdt, d2gdt2, d3gdt3... Derivatives wrt t (1st, 2nd, 3rd...) or the
    //      couplings.
    //tcl - log(phicl^2/mt^2). Parameter used to relate scalar field to
    //      renormalisation scale.


    //Functions which are used by clients of this class:
    //These functions compute the derivatives of t(tcl). They require
    //knowledge of the derivatives of the couplings with respect to t(not tcl!).
    //That information must be provided by (external) beta functions:
    value_type dt_dtcl(const value_type& t, const std::vector<value_type>& g,
                       const std::vector<value_type>& dgdt);
    value_type d2t_dtcl2(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const value_type& dtdtcl);
    value_type d3t_dtcl3(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const std::vector<value_type>& d3gdt3,
                         const value_type& dtdtcl,const value_type& d2tdtcl2);
    value_type d4t_dtcl4(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const std::vector<value_type>& d3gdt3,
                         const std::vector<value_type>& d4gdt4,
                         const value_type& dtdtcl,const value_type& d2tdtcl2,
                         const value_type& d3tdtcl3);
    value_type d5t_dtcl5(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const std::vector<value_type>& d3gdt3,
                         const std::vector<value_type>& d4gdt4,
                         const std::vector<value_type>& d5gdt5,
                         const value_type& dtdtcl,const value_type& d2tdtcl2,
                         const value_type& d3tdtcl3,const value_type& d4tdtcl4);

    //Versions of the above functions for use when additional information is
    //known (this can be useful for speeding up calculations by avoiding
    //repetition of work):
    value_type dt_dtcl(const value_type& t, const std::vector<value_type>& g,
                       const std::vector<value_type>& dgdt,
                       const std::vector<value_type>& logfactor,
                       const std::vector<value_type>& Mi2,
                       const std::vector<value_type>& pMi2pt,
                       const std::vector<value_type>& kiarray,
                       const value_type Z,const value_type phicl);
    value_type d2t_dtcl2(const value_type& t, const std::vector<value_type>& g,
                       const std::vector<value_type>& logfactor,
                       const std::vector<value_type>& Mi2,
                       const std::vector<value_type>& pMi2pt,
                       const std::vector<value_type>& Mi2p,
                       const std::vector<value_type>& pMi2ptp,
                       const std::vector<value_type>& kiarray,
                       const std::vector<value_type>& kiparray,
                       const value_type& Z,const value_type& Zd,
                       const value_type& Zdd,const value_type& mu,
                       const value_type phicl,const value_type& dtdtcl);
    value_type d3t_dtcl3(const value_type& t, const std::vector<value_type>& g,
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
                       const value_type& Z,const value_type& Zd,
                       const value_type& Zdd,const value_type& Zddd,
                       const value_type& mu,
                       const value_type phicl,const value_type& dtdtcl,
                       const value_type& d2tdtcl2);
    value_type d4t_dtcl4(const value_type& t, const std::vector<value_type>& g,
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
                         const value_type& Z,const value_type& Zd,
                         const value_type& Zdd,const value_type& Zddd,
                         const value_type& mu,
                         const value_type phicl,const value_type& dtdtcl,
                         const value_type& d2tdtcl2,const value_type& d3tdtcl3);
    value_type d5t_dtcl5(const value_type& t, const std::vector<value_type>& g,
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
                         const value_type& Z,const value_type& Zd,
                         const value_type& Zdd,const value_type& Zddd,
                         const value_type& Zdddd,
                         const value_type& mu,
                         const value_type phicl,const value_type& dtdtcl,
                         const value_type& d2tdtcl2,const value_type& d3tdtcl3,
                         const value_type& d4tdtcl4);


private:

    //Functions for internal use only:

    //Sets up arrays of coupling-dependent coefficients appearing in the
    //1-loop corrections. The storage vectors k,kp,& theta do not need to be
    //defined ahead of time - these functions will redefine them:


    void populateCoefficientArrays
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& g)
    {
        this->coeff.populateCoefficientArrays(k,kp,theta,g);
        /*
        value_type kiarray[] = {g[2]/this->_4p0,g[2]/this->_4p0,
                                g[2]/this->_4p0,
                                (g[2] + g[1])/this->_4p0,
                                (g[2] + g[1])/this->_4p0,
                                (g[2] + g[1])/this->_4p0,g[4]/this->_2p0,
                                this->_3p0*g[0],g[0]};
        value_type kpiarray[] = {this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,g[5],g[5]};
        value_type thetaiarray[] = {this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                -this->_1p0/this->_6p0,this->_1p0/this->_12p0,
                                this->_1p0/this->_12p0,-this->_1p0/this->_6p0,
                                this->_1p0/this->_12p0,
                                g[6] - this->_1p0/this->_6p0,
                                g[6] - this->_1p0/this->_6p0};


        k.assign(kiarray,kiarray + 9);
        kp.assign(kpiarray,kpiarray + 9);
        theta.assign(thetaiarray,thetaiarray + 9);
        */
    }
    //Derivatives of the coupling-dependent coefficients. If higher than 1st
    //derivatives are desired, then gp should be a vector of such
    //higher derivatives:
    void populateCoefficientArraysDerivatives
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& gp)
    {
        this->coeff.populateCoefficientArraysDerivatives(k,kp,theta,gp);
        /*
        value_type kiarray[] = {gp[2]/this->_4p0,gp[2]/this->_4p0,
                                gp[2]/this->_4p0,
                                (gp[2] + gp[1])/this->_4p0,
                                (gp[2] + gp[1])/this->_4p0,
                                (gp[2] + gp[1])/this->_4p0,gp[4]/this->_2p0,
                                this->_3p0*gp[0],gp[0]};
        value_type kpiarray[] = {this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,gp[5],gp[5]};
        value_type thetaiarray[] = {this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,
                                this->_0p0,gp[6],gp[6]};


        k.assign(kiarray,kiarray + 9);
        kp.assign(kpiarray,kpiarray + 9);
        theta.assign(thetaiarray,thetaiarray + 9);
        */
    }



    //Functions to compute Mi^2(phi) and its derivatives, which depend on the
    //above coupling dependent coefficients. Derivatives (wrt tcl)
    // are denoted with
    //the suffix p, second derivatives pp etc...

    void computeMi2
    (std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
     std::vector<value_type>& pMi2pt,
     const std::vector<value_type>& ki,const std::vector<value_type>& kpi,
     const std::vector<value_type>& thetai,
     const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
     const std::vector<value_type>& dthetaidt,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& logmu2);

     void computeMi2p
     (std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
     std::vector<value_type>& Mi2p,
     std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
     const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
     const std::vector<value_type>& thetai,
     const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
     const std::vector<value_type>& dthetaidt,
     const std::vector<value_type>& d2kidt2,
     const std::vector<value_type>& d2kpidt2,
     const std::vector<value_type>& d2thetaidt2,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& Zdd,const value_type& dtdtcl,const value_type& logmu2);

     void computeMi2pp
     (std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
     std::vector<value_type>& Mi2p,std::vector<value_type>& Mi2pp,
     std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
     std::vector<value_type>& pMi2ptpp,
     const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
     const std::vector<value_type>& thetai,
     const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
     const std::vector<value_type>& dthetaidt,
     const std::vector<value_type>& d2kidt2,
     const std::vector<value_type>& d2kpidt2,
     const std::vector<value_type>& d2thetaidt2,
     const std::vector<value_type>& d3kidt3,
     const std::vector<value_type>& d3kpidt3,
     const std::vector<value_type>& d3thetaidt3,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& Zdd,const value_type& Zddd,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& logmu2);

     void computeMi2ppp
     (std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
     std::vector<value_type>& Mi2p,std::vector<value_type>& Mi2pp,
     std::vector<value_type>& Mi2ppp,
     std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
     std::vector<value_type>& pMi2ptpp,std::vector<value_type>& pMi2ptppp,
     const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
     const std::vector<value_type>& thetai,
     const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
     const std::vector<value_type>& dthetaidt,
     const std::vector<value_type>& d2kidt2,
     const std::vector<value_type>& d2kpidt2,
     const std::vector<value_type>& d2thetaidt2,
     const std::vector<value_type>& d3kidt3,
     const std::vector<value_type>& d3kpidt3,
     const std::vector<value_type>& d3thetaidt3,
     const std::vector<value_type>& d4kidt4,
     const std::vector<value_type>& d4kpidt4,
     const std::vector<value_type>& d4thetaidt4,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& Zdd,const value_type& Zddd,const value_type& Zdddd,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,
     const value_type& logmu2);

     void computeMi2pppp
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
     const std::vector<value_type>& d2kidt2,
     const std::vector<value_type>& d2kpidt2,
     const std::vector<value_type>& d2thetaidt2,
     const std::vector<value_type>& d3kidt3,
     const std::vector<value_type>& d3kpidt3,
     const std::vector<value_type>& d3thetaidt3,
     const std::vector<value_type>& d4kidt4,
     const std::vector<value_type>& d4kpidt4,
     const std::vector<value_type>& d4thetaidt4,
     const std::vector<value_type>& d5kidt5,
     const std::vector<value_type>& d5kpidt5,
     const std::vector<value_type>& d5thetaidt5,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& Zdd,const value_type& Zddd,const value_type& Zdddd,
     const value_type& Zddddd,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,const value_type& d4tdtcl4,
     const value_type& logmu2);

    //dt/dtcl is actually defined as a quotient function, A(x)/B(x).
    //The derivatives are simplifed by first computing the derivatives of
    //A and B respectively, and combining the results. This also avoids a lot
    //of repetition because the results for A,A',A'' etc... can be saved,
    //making the calculation more efficient. We need up to 3rd derivaites
    //because d$4t/dtcl^4 is the highest derivative we compute.

    //A
    value_type A
    (const std::vector<value_type>& logfactor,const value_type& phicl,
     const std::vector<value_type>& Mi2,const std::vector<value_type>& kiarray,
     const value_type& Z);
    //B
    value_type B
    (const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2,
     const std::vector<value_type>& pMi2pt);
    //Ap
    value_type Ap
    (const value_type& Z,const value_type& Zd,
    const std::vector<value_type>& kiarray,
    const std::vector<value_type>& kiparray,
    const value_type& mu,const std::vector<value_type>& logfactor,
    const value_type& phicl,
    const std::vector<value_type>& Mi2array,
    const std::vector<value_type>& Mi2parray,
    const value_type& dtdtcl);
    //App
    value_type App
    (const value_type& Z,const value_type& Zd,const value_type& Zdd,
     const std::vector<value_type>& kiarray,
     const std::vector<value_type>& kiparray,
     const std::vector<value_type>& kipparray,
     const value_type& mu,const std::vector<value_type>& logfactor,
     const value_type& phicl,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const value_type& dtdtcl,const value_type& d2tdtcl2);
    //Appp
    value_type Appp
    (const value_type& Z,const value_type& Zd,const value_type& Zdd,
     const value_type& Zddd,
     const std::vector<value_type>& kiarray,
     const std::vector<value_type>& kiparray,
     const std::vector<value_type>& kipparray,
     const std::vector<value_type>& kippparray,
     const value_type& mu,const std::vector<value_type>& logfactor,
     const value_type& phicl,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& Mi2ppparray,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3);

     value_type Apppp
    (const value_type& Z,const value_type& Zd,const value_type& Zdd,
     const value_type& Zddd,
     const value_type& Zdddd,
     const std::vector<value_type>& kiarray,
     const std::vector<value_type>& kiparray,
     const std::vector<value_type>& kipparray,
     const std::vector<value_type>& kippparray,
     const std::vector<value_type>& kipppparray,
     const value_type& mu,const std::vector<value_type>& logfactor,
     const value_type& phicl,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& Mi2ppparray,
     const std::vector<value_type>& Mi2pppparray,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,
     const value_type& d4tdtcl4);
    //Bp
    value_type Bp
    (const value_type& dtdtcl,const value_type& mu,const value_type& phicl,
     const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& pMi2ptarray,
     const std::vector<value_type>& pMi2ptparray);
    //Bpp
    value_type Bpp
    (const value_type& mu,const value_type& phicl,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& pMi2ptarray,
     const std::vector<value_type>& pMi2ptparray,
     const std::vector<value_type>& pMi2ptpparray);
    //Bppp
    value_type Bppp
    (const value_type& mu,const value_type& phicl,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,
     const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& Mi2ppparray,
     const std::vector<value_type>& pMi2ptarray,
     const std::vector<value_type>& pMi2ptparray,
     const std::vector<value_type>& pMi2ptpparray,
     const std::vector<value_type>& pMi2ptppparray);

     value_type Bpppp
    (const value_type& mu,const value_type& phicl,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,
     const value_type& d4tdtcl4,
     const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& Mi2ppparray,
     const std::vector<value_type>& Mi2pppparray,
     const std::vector<value_type>& pMi2ptarray,
     const std::vector<value_type>& pMi2ptparray,
     const std::vector<value_type>& pMi2ptpparray,
     const std::vector<value_type>& pMi2ptppparray,
     const std::vector<value_type>& pMi2ptpppparray);

     //Functions used to compute various intermediaries that appear:
     //
     //Parameters:
     //See above for logfactor, Mi2, pMi2pt, Z etc...
     //t - renormalisation parameter (dimensionless): t = log(mu^2/mt^2)
     //g - vector of couplings at scale t.
     //dgdt - derivatives of g wrt t
    void compute_dt_dtcl_data
    (const value_type& t, const std::vector<value_type>& g,
     const std::vector<value_type>& dgdt,std::vector<value_type>& logfactor,
     std::vector<value_type>& Mi2,std::vector<value_type>& pMi2pt,
     std::vector<value_type>& kiarray,value_type& Z,value_type& phicl);

    void compute_d2t_dtcl2_data
    (const value_type& tcl, const std::vector<value_type>& g,
     const std::vector<value_type>& dgdt,const std::vector<value_type>& d2gdt2,
     const value_type& dtdtcl,
     std::vector<value_type>& logfactor,std::vector<value_type>& Mi2,
     std::vector<value_type>& pMi2pt,std::vector<value_type>& Mi2p,
     std::vector<value_type>& pMi2ptp,std::vector<value_type>& kiarray,
     std::vector<value_type>& kiparray,value_type& Z,value_type& Zp,
     value_type& Zpp,value_type& mu,value_type& phicl);

    void compute_d3t_dtcl3_data
    (const value_type& t, const std::vector<value_type>& g,
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
     const value_type& d2tdtcl2);

    void compute_d4t_dtcl4_data
    (const value_type& t, const std::vector<value_type>& g,
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
     const value_type& d2tdtcl2,const value_type& d3tdtcl3);

     void compute_d5t_dtcl5_data
    (const value_type& t, const std::vector<value_type>& g,
     const std::vector<value_type>& dgdt,const std::vector<value_type>& d2gdt2,
     const std::vector<value_type>& d3gdt3,
     const std::vector<value_type>& d4gdt4,
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
     std::vector<value_type>& kippparray,
     std::vector<value_type>& kipppparray,
     value_type& Z,
     value_type& Zp,value_type& Zpp,value_type& Zppp,value_type& Zpppp,
     value_type& Zppppp,
     value_type& mu,value_type& phicl,const value_type& dtdtcl,
     const value_type& d2tdtcl2,const value_type& d3tdtcl3,
     const value_type& d4tdtcl4);

};
/*
//Version of dtdtcl_1loop that includes the gravitational terms:
template<class value_type>
class dtdtcl_1loop_grav : public dtdtcl_jacobian<value_type>
{
    //Reference list of the couplings vector and which
    //SM couplings they refer to:
    /*
    g[0] - lambda (Higgs self coupling)
    g[1] - g1^2 (U(1) coupling)
    g[2] - g2^2 (SU(2) coupling)
    g[3] - g3^2 (SU(3) coupling)
    g[4] - yt^2 (Top quark Yukawa coupling)
    g[5] - m^2 (Higgs tachyonic mass term)
    g[6] - xi (Higgs-curvature coupling)
    g[7] - Z (Higgs field renormalisation) - technically the RHS of this
                  is a gamma function, but we solve it with the beta functions)
    //Additionally, we have to solve for the non-beta functions:
    g[8] = t_cl = log(phicl^2/mt^2)
    g[9] = t = log(mu^2/mt^2)
    //Cosmological constant, if we include it:
    g[10] = Lambda*Mp^2 (cosmological constant, in units of the Planck
                                mass)
    g[11] = kappa (gravitational constant)
    g[12] = alpha1 (R^2 term)
    g[13] = alpha2 (R_{\mu\nu}R^{\mu\nu} term)
    g[14] = alpha3 (R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma} term)

    */
    /*
    //Hard coded numerical constants:
    const value_type _0p0,_1p0,_16p0,_3p0,_8p0,_9p0,_m3p0,_2p0,_m9p0,_12p0;
    const value_type _6p0,_24p0,_m6p0,_4p0,_27p0,_m27p0,_72p0,_18p0,_m18p0;
    const value_type _96p0,_36p0,_m24p0,_25p0,_2304p0,_10p0,_5p0,_11p0,_45p0;
    const value_type _1152p0,_384p0,_288p0,_256p0,_m43p0,_128p0,_32p0,_m11p0;
    const value_type _m17p0,_m8p0;
    const value_type _m4p0,_20p0,_m1p0,_m2p0,_m12p0;

    //Physical constants:
    const value_type mt;//Top quark mass (used as the base scale wrt which mu
        //in the RG equations is defined, ie, t = log(mu^2/mt^2) is the
        //dimensionless running parameter.
    const value_type R;//Ricci scalar.
    //Coefficients appearing in the 1-loop corrections to the SM effective
        //potential:
    std::vector<value_type> nilist;
    std::vector<value_type> cilist;
    std::vector<value_type> npilist;
    //Function to fill these coefficients when the constructor is called
    //(separated into a function of its own right in case we want to add
    //multiple constructors at some point).

    //Tells the object how many loop correction terms are used in the
    //SM potential:
    static const int nLoops = 12;






public:

    //Constructor:
    dtdtcl_1loop(value_type Mt,value_type r)
    : _1p0(1), _16p0(16), _3p0(3), _8p0(8), _9p0(9), _m3p0(-3), _2p0(2),
    _m9p0(-9),
    _12p0(12), _6p0(6), _24p0(24), _m6p0(-6), _4p0(4), _27p0(27), _m27p0(-27),
    _72p0(72), _18p0(18), _m18p0(-18), _96p0(96), _36p0(36), _m24p0(-24),
    _25p0(25), _2304p0(2304), _10p0(10), _5p0(5), _11p0(11), _45p0(45),
    _1152p0(1152), _384p0(384), _288p0(288), _256p0(256), _m43p0(-43),
    _128p0(128),_m4p0(-4.0),_20p0(20.0),_m1p0(-1.0),_m2p0(-2.0),_m12p0(-12.0),
    _32p0(32), _m11p0(-11), _m17p0(-17), _m8p0(-8), _0p0(0.0),
    mt(Mt), R(r)
    {
        fill_constant_arrays(this->nilist,this->cilist,this->npilist);
    }
    //Parameters to the various functions used in this class:
    //
    //logfactor - log(Mi2/mu^2) (vector for different values of i)
    //phicl - classical field, ie, field evaluated at Electroweak scale.
    //          phicl mt*exp(tcl/2)
    //Mi2 - Mi2(phi) = k_i*phi^2(t) - k'_i + theta_i*R
    //      Where phi = Z(t)*phicl, R = Ricci scalar, and k_i, k'_i, theta_i
    //      are coupling dependent coefficients (depending on different
    //      couplings for each i).
    //kiarray - k_i coefficient array (see populateCoefficientArrays for details).
    //pMi2pt - partial derivative wrt t of Mi2, at constant phicl.
    //          pMi2pt = (dk_i/dt)*Z^2(t)*phicl^2 + k_i*2*Z(t)*phicl^2
    //                   -(dk'_i/dt) + (dtheta_i/dt)*R
    //kiparray - derivatives of kiarray wrt t (not tcl).
    //kipparray - higher t derivatives of k_i (# of ps = # of derivatives)
    //Mi2parray - derivatives of kiarray wrt tcl (not t!)
    //Mi2pparray - higher tcl derivatives of Mi2 (# of ps = # of derivatives)
    //dtdtcl - 1st derivative of t wrt tcl.
    //d2tdtcl2 - 2nd derivative of t wrt tcl.
    //...
    //dntdtcln - nth derivative of t wrt tcl.
    //Z - field renormalisation, defined as phi(t) = Z(t)*phicl, where phi(t)
    //      is the renormalised field.
    //Zd, Zdd... - t derivatives of Z (# of ds = # of derivatives).
    //pMi2ptparray, pMi2ptpparray... - Derivatives of pMi2pt wrt tcl (not t)
    //mu - renormalisation scale.
    //thetai - theta_i coefficient found in Mi2.
    //dthetaidt, d2thetaidt2 ... - derivatives (1st, 2nd...) of theta_i wrt t.
    //ki - k_i coefficient found in Mi2 (alias for kiarray - both are sometimes
    //      used in the same function because an array is used to define the
    //      coefficients, and a vector to store them, so two names are used)
    //kpi - k'_i coefficients found in Mi2
    //dkidt, d2kidt2,... - t derivatives (1st, 2nd...) of k_i
    //dkpidt, d2kpidt2,... - t derivatives (1st, 2nd...) of k'_i
    //logmu2 = log(mu^2/mt^2)
    //t - renormalisation running parameter, t = log(mu^2/mt^2)
    //g - couplings vector
    //dgdt, d2gdt2, d3gdt3... Derivatives wrt t (1st, 2nd, 3rd...) or the
    //      couplings.
    //tcl - log(phicl^2/mt^2). Parameter used to relate scalar field to
    //      renormalisation scale.


    //Functions which are used by clients of this class:
    //These functions compute the derivatives of t(tcl). They require
    //knowledge of the derivatives of the couplings with respect to t(not tcl!).
    //That information must be provided by (external) beta functions:
    value_type dt_dtcl(const value_type& t, const std::vector<value_type>& g,
                       const std::vector<value_type>& dgdt);
    value_type d2t_dtcl2(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const value_type& dtdtcl);
    value_type d3t_dtcl3(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const std::vector<value_type>& d3gdt3,
                         const value_type& dtdtcl,const value_type& d2tdtcl2);
    value_type d4t_dtcl4(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const std::vector<value_type>& d3gdt3,
                         const std::vector<value_type>& d4gdt4,
                         const value_type& dtdtcl,const value_type& d2tdtcl2,
                         const value_type& d3tdtcl3);
    value_type d5t_dtcl5(const value_type& tcl,
                         const std::vector<value_type>& g,
                         const std::vector<value_type>& dgdt,
                         const std::vector<value_type>& d2gdt2,
                         const std::vector<value_type>& d3gdt3,
                         const std::vector<value_type>& d4gdt4,
                         const std::vector<value_type>& d5gdt5,
                         const value_type& dtdtcl,const value_type& d2tdtcl2,
                         const value_type& d3tdtcl3,const value_type& d4tdtcl4);

    //Versions of the above functions for use when additional information is
    //known (this can be useful for speeding up calculations by avoiding
    //repetition of work):
    value_type dt_dtcl(const value_type& t, const std::vector<value_type>& g,
                       const std::vector<value_type>& dgdt,
                       const std::vector<value_type>& logfactor,
                       const std::vector<value_type>& Mi2,
                       const std::vector<value_type>& pMi2pt,
                       const std::vector<value_type>& kiarray,
                       const value_type Z,const value_type phicl);
    value_type d2t_dtcl2(const value_type& t, const std::vector<value_type>& g,
                       const std::vector<value_type>& logfactor,
                       const std::vector<value_type>& Mi2,
                       const std::vector<value_type>& pMi2pt,
                       const std::vector<value_type>& Mi2p,
                       const std::vector<value_type>& pMi2ptp,
                       const std::vector<value_type>& kiarray,
                       const std::vector<value_type>& kiparray,
                       const value_type& Z,const value_type& Zd,
                       const value_type& Zdd,const value_type& mu,
                       const value_type phicl,const value_type& dtdtcl);
    value_type d3t_dtcl3(const value_type& t, const std::vector<value_type>& g,
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
                       const value_type& Z,const value_type& Zd,
                       const value_type& Zdd,const value_type& Zddd,
                       const value_type& mu,
                       const value_type phicl,const value_type& dtdtcl,
                       const value_type& d2tdtcl2);
    value_type d4t_dtcl4(const value_type& t, const std::vector<value_type>& g,
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
                         const value_type& Z,const value_type& Zd,
                         const value_type& Zdd,const value_type& Zddd,
                         const value_type& mu,
                         const value_type phicl,const value_type& dtdtcl,
                         const value_type& d2tdtcl2,const value_type& d3tdtcl3);
    value_type d5t_dtcl5(const value_type& t, const std::vector<value_type>& g,
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
                         const value_type& Z,const value_type& Zd,
                         const value_type& Zdd,const value_type& Zddd,
                         const value_type& Zdddd,
                         const value_type& mu,
                         const value_type phicl,const value_type& dtdtcl,
                         const value_type& d2tdtcl2,const value_type& d3tdtcl3,
                         const value_type& d4tdtcl4);


private:

    //Functions for internal use only:

    //Sets up arrays of coupling-dependent coefficients appearing in the
    //1-loop corrections. The storage vectors k,kp,& theta do not need to be
    //defined ahead of time - these functions will redefine them:
    void populateCoefficientArrays
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& g)
    {
        value_type kiarray[] = {g[2]/this->_4p0,g[2]/this->_4p0,
                                g[2]/this->_4p0,
                                (g[2] + g[1])/this->_4p0,
                                (g[2] + g[1])/this->_4p0,
                                (g[2] + g[1])/this->_4p0,g[4]/this->_2p0,
                                this->_3p0*g[0],g[0]};
        value_type kpiarray[] = {this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,g[5],g[5]};
        value_type thetaiarray[] = {this->_1p0/this->_12p0,
                                    this->_1p0/this->_12p0,
                                -this->_1p0/this->_6p0,this->_1p0/this->_12p0,
                                this->_1p0/this->_12p0,-this->_1p0/this->_6p0,
                                this->_1p0/this->_12p0,
                                g[6] - this->_1p0/this->_6p0,
                                g[6] - this->_1p0/this->_6p0};


        k.assign(kiarray,kiarray + 9);
        kp.assign(kpiarray,kpiarray + 9);
        theta.assign(thetaiarray,thetaiarray + 9);
    }
    //Derivatives of the coupling-dependent coefficients. If higher than 1st
    //derivatives are desired, then gp should be a vector of such
    //higher derivatives:
    void populateCoefficientArraysDerivatives
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& gp)
    {
        value_type kiarray[] = {gp[2]/this->_4p0,gp[2]/this->_4p0,
                                gp[2]/this->_4p0,
                                (gp[2] + gp[1])/this->_4p0,
                                (gp[2] + gp[1])/this->_4p0,
                                (gp[2] + gp[1])/this->_4p0,gp[4]/this->_2p0,
                                this->_3p0*gp[0],gp[0]};
        value_type kpiarray[] = {this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,this->_0p0,
                                this->_0p0,gp[5],gp[5]};
        value_type thetaiarray[] = {this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,
                                this->_0p0,this->_0p0,
                                this->_0p0,gp[6],gp[6]};


        k.assign(kiarray,kiarray + 9);
        kp.assign(kpiarray,kpiarray + 9);
        theta.assign(thetaiarray,thetaiarray + 9);
    }


    //Functions to compute Mi^2(phi) and its derivatives, which depend on the
    //above coupling dependent coefficients. Derivatives (wrt tcl)
    // are denoted with
    //the suffix p, second derivatives pp etc...

    void computeMi2
    (std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
     std::vector<value_type>& pMi2pt,
     const std::vector<value_type>& ki,const std::vector<value_type>& kpi,
     const std::vector<value_type>& thetai,
     const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
     const std::vector<value_type>& dthetaidt,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& logmu2);

     void computeMi2p
     (std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
     std::vector<value_type>& Mi2p,
     std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
     const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
     const std::vector<value_type>& thetai,
     const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
     const std::vector<value_type>& dthetaidt,
     const std::vector<value_type>& d2kidt2,
     const std::vector<value_type>& d2kpidt2,
     const std::vector<value_type>& d2thetaidt2,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& Zdd,const value_type& dtdtcl,const value_type& logmu2);

     void computeMi2pp
     (std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
     std::vector<value_type>& Mi2p,std::vector<value_type>& Mi2pp,
     std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
     std::vector<value_type>& pMi2ptpp,
     const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
     const std::vector<value_type>& thetai,
     const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
     const std::vector<value_type>& dthetaidt,
     const std::vector<value_type>& d2kidt2,
     const std::vector<value_type>& d2kpidt2,
     const std::vector<value_type>& d2thetaidt2,
     const std::vector<value_type>& d3kidt3,
     const std::vector<value_type>& d3kpidt3,
     const std::vector<value_type>& d3thetaidt3,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& Zdd,const value_type& Zddd,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& logmu2);

     void computeMi2ppp
     (std::vector<value_type>& Mi2,std::vector<value_type>& logfactor,
     std::vector<value_type>& Mi2p,std::vector<value_type>& Mi2pp,
     std::vector<value_type>& Mi2ppp,
     std::vector<value_type>& pMi2pt,std::vector<value_type>& pMi2ptp,
     std::vector<value_type>& pMi2ptpp,std::vector<value_type>& pMi2ptppp,
     const std::vector<value_type>& kiarray,const std::vector<value_type>& kpi,
     const std::vector<value_type>& thetai,
     const std::vector<value_type>& dkidt,const std::vector<value_type>& dkpidt,
     const std::vector<value_type>& dthetaidt,
     const std::vector<value_type>& d2kidt2,
     const std::vector<value_type>& d2kpidt2,
     const std::vector<value_type>& d2thetaidt2,
     const std::vector<value_type>& d3kidt3,
     const std::vector<value_type>& d3kpidt3,
     const std::vector<value_type>& d3thetaidt3,
     const std::vector<value_type>& d4kidt4,
     const std::vector<value_type>& d4kpidt4,
     const std::vector<value_type>& d4thetaidt4,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& Zdd,const value_type& Zddd,const value_type& Zdddd,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,
     const value_type& logmu2);

     void computeMi2pppp
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
     const std::vector<value_type>& d2kidt2,
     const std::vector<value_type>& d2kpidt2,
     const std::vector<value_type>& d2thetaidt2,
     const std::vector<value_type>& d3kidt3,
     const std::vector<value_type>& d3kpidt3,
     const std::vector<value_type>& d3thetaidt3,
     const std::vector<value_type>& d4kidt4,
     const std::vector<value_type>& d4kpidt4,
     const std::vector<value_type>& d4thetaidt4,
     const std::vector<value_type>& d5kidt5,
     const std::vector<value_type>& d5kpidt5,
     const std::vector<value_type>& d5thetaidt5,
     const value_type& phicl,const value_type& Z,const value_type& Zd,
     const value_type& Zdd,const value_type& Zddd,const value_type& Zdddd,
     const value_type& Zddddd,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,const value_type& d4tdtcl4,
     const value_type& logmu2);

    //dt/dtcl is actually defined as a quotient function, A(x)/B(x).
    //The derivatives are simplifed by first computing the derivatives of
    //A and B respectively, and combining the results. This also avoids a lot
    //of repetition because the results for A,A',A'' etc... can be saved,
    //making the calculation more efficient. We need up to 3rd derivaites
    //because d$4t/dtcl^4 is the highest derivative we compute.

    //A
    value_type A
    (const std::vector<value_type>& logfactor,const value_type& phicl,
     const std::vector<value_type>& Mi2,const std::vector<value_type>& kiarray,
     const value_type& Z);
    //B
    value_type B
    (const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2,
     const std::vector<value_type>& pMi2pt);
    //Ap
    value_type Ap
    (const value_type& Z,const value_type& Zd,
    const std::vector<value_type>& kiarray,
    const std::vector<value_type>& kiparray,
    const value_type& mu,const std::vector<value_type>& logfactor,
    const value_type& phicl,
    const std::vector<value_type>& Mi2array,
    const std::vector<value_type>& Mi2parray,
    const value_type& dtdtcl);
    //App
    value_type App
    (const value_type& Z,const value_type& Zd,const value_type& Zdd,
     const std::vector<value_type>& kiarray,
     const std::vector<value_type>& kiparray,
     const std::vector<value_type>& kipparray,
     const value_type& mu,const std::vector<value_type>& logfactor,
     const value_type& phicl,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const value_type& dtdtcl,const value_type& d2tdtcl2);
    //Appp
    value_type Appp
    (const value_type& Z,const value_type& Zd,const value_type& Zdd,
     const value_type& Zddd,
     const std::vector<value_type>& kiarray,
     const std::vector<value_type>& kiparray,
     const std::vector<value_type>& kipparray,
     const std::vector<value_type>& kippparray,
     const value_type& mu,const std::vector<value_type>& logfactor,
     const value_type& phicl,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& Mi2ppparray,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3);

     value_type Apppp
    (const value_type& Z,const value_type& Zd,const value_type& Zdd,
     const value_type& Zddd,
     const value_type& Zdddd,
     const std::vector<value_type>& kiarray,
     const std::vector<value_type>& kiparray,
     const std::vector<value_type>& kipparray,
     const std::vector<value_type>& kippparray,
     const std::vector<value_type>& kipppparray,
     const value_type& mu,const std::vector<value_type>& logfactor,
     const value_type& phicl,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& Mi2ppparray,
     const std::vector<value_type>& Mi2pppparray,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,
     const value_type& d4tdtcl4);
    //Bp
    value_type Bp
    (const value_type& dtdtcl,const value_type& mu,const value_type& phicl,
     const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& pMi2ptarray,
     const std::vector<value_type>& pMi2ptparray);
    //Bpp
    value_type Bpp
    (const value_type& mu,const value_type& phicl,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& pMi2ptarray,
     const std::vector<value_type>& pMi2ptparray,
     const std::vector<value_type>& pMi2ptpparray);
    //Bppp
    value_type Bppp
    (const value_type& mu,const value_type& phicl,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,
     const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& Mi2ppparray,
     const std::vector<value_type>& pMi2ptarray,
     const std::vector<value_type>& pMi2ptparray,
     const std::vector<value_type>& pMi2ptpparray,
     const std::vector<value_type>& pMi2ptppparray);

     value_type Bpppp
    (const value_type& mu,const value_type& phicl,
     const value_type& dtdtcl,const value_type& d2tdtcl2,
     const value_type& d3tdtcl3,
     const value_type& d4tdtcl4,
     const std::vector<value_type>& logfactor,
     const std::vector<value_type>& Mi2array,
     const std::vector<value_type>& Mi2parray,
     const std::vector<value_type>& Mi2pparray,
     const std::vector<value_type>& Mi2ppparray,
     const std::vector<value_type>& Mi2pppparray,
     const std::vector<value_type>& pMi2ptarray,
     const std::vector<value_type>& pMi2ptparray,
     const std::vector<value_type>& pMi2ptpparray,
     const std::vector<value_type>& pMi2ptppparray,
     const std::vector<value_type>& pMi2ptpppparray);

     //Functions used to compute various intermediaries that appear:
     //
     //Parameters:
     //See above for logfactor, Mi2, pMi2pt, Z etc...
     //t - renormalisation parameter (dimensionless): t = log(mu^2/mt^2)
     //g - vector of couplings at scale t.
     //dgdt - derivatives of g wrt t
    void compute_dt_dtcl_data
    (const value_type& t, const std::vector<value_type>& g,
     const std::vector<value_type>& dgdt,std::vector<value_type>& logfactor,
     std::vector<value_type>& Mi2,std::vector<value_type>& pMi2pt,
     std::vector<value_type>& kiarray,value_type& Z,value_type& phicl);

    void compute_d2t_dtcl2_data
    (const value_type& tcl, const std::vector<value_type>& g,
     const std::vector<value_type>& dgdt,const std::vector<value_type>& d2gdt2,
     const value_type& dtdtcl,
     std::vector<value_type>& logfactor,std::vector<value_type>& Mi2,
     std::vector<value_type>& pMi2pt,std::vector<value_type>& Mi2p,
     std::vector<value_type>& pMi2ptp,std::vector<value_type>& kiarray,
     std::vector<value_type>& kiparray,value_type& Z,value_type& Zp,
     value_type& Zpp,value_type& mu,value_type& phicl);

    void compute_d3t_dtcl3_data
    (const value_type& t, const std::vector<value_type>& g,
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
     const value_type& d2tdtcl2);

    void compute_d4t_dtcl4_data
    (const value_type& t, const std::vector<value_type>& g,
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
     const value_type& d2tdtcl2,const value_type& d3tdtcl3);

     void compute_d5t_dtcl5_data
    (const value_type& t, const std::vector<value_type>& g,
     const std::vector<value_type>& dgdt,const std::vector<value_type>& d2gdt2,
     const std::vector<value_type>& d3gdt3,
     const std::vector<value_type>& d4gdt4,
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
     std::vector<value_type>& kippparray,
     std::vector<value_type>& kipppparray,
     value_type& Z,
     value_type& Zp,value_type& Zpp,value_type& Zppp,value_type& Zpppp,
     value_type& Zppppp,
     value_type& mu,value_type& phicl,const value_type& dtdtcl,
     const value_type& d2tdtcl2,const value_type& d3tdtcl3,
     const value_type& d4tdtcl4);
};
*/
//------------------------------------------------------------------------------
//SM beta functions in de Sitter space at one loop:
#include "beta_function_definition.h"
template<class value_type>
class betaSMdS_1loop : public betaFunctionWithDerivatives<value_type>
{

private:
    //We need to decide on an order for our beta functions:
        /*
    //Second argument, [i], refers to the element of the computed grid. First
    //argument determined the coupling.
    g[0] - lambda (Higgs self coupling)
    g[1] - g1^2 (U(1) coupling)
    g[2] - g2^2 (SU(2) coupling)
    g[3] - g3^2 (SU(3) coupling)
    g[4] - yt^2 (Top quark Yukawa coupling)
    g[5] - m^2 (Higgs tachyonic mass term)
    g[6] - xi (Higgs-curvature coupling)
    g[7] - Z (Higgs field renormalisation) - technically the RHS of this
                  is a gamma function, but we solve it with the beta functions)
    //Additionally, we have to solve for the non-beta functions:
    g[8] = t_cl = log(phicl^2/mt^2)
    g[9] = t = log(mu^2/mt^2)
    //Cosmological constant, if we include it:
    g[10] = Lambda*Mp^2 (cosmological constant, in units of the Planck
                                mass)
    */

    //Hard coded numerical constants:
    const value_type _0p0,_1p0,_16p0,_3p0,_8p0,_9p0,_m3p0,_2p0,_m9p0,_12p0;
    const value_type _6p0,_24p0,_m6p0,_4p0,_27p0,_m27p0,_72p0,_18p0,_m18p0;
    const value_type _96p0,_36p0,_m24p0,_25p0,_2304p0,_10p0,_5p0,_11p0,_45p0;
    const value_type _1152p0,_384p0,_288p0,_256p0,_m43p0,_128p0,_32p0,_m11p0;
    const value_type _m17p0,_m8p0,_15p0;
    const value_type _40p0,_20p0,_200p0,_81p0,_50p0,_400p0;

    const value_type _m29p0;
    const value_type _m2p0;
    const value_type _68p0;
    const value_type _m49p0;
    const value_type _360p0;
    const value_type _5760p0;
    const value_type _17p0;
    const value_type _180p0;

    const value_type _m277p0,_571p0,_2880p0,_293p0,_23040p0,_m293p0;

    const value_type _m1p0;

    value_type PI;
    value_type loop_factor;//1/(16pi^2)
    const value_type Ng;//Number of quark generations (3)
    const value_type Nc;//Number of colours (3)
    const value_type Tf;//Normalises SU(N) generators: Tr(t^a t^b) = Tf*delta_ab
    const value_type Ca;//Colour invariant associated to structure constants:
            //\sum_{a,b} f_{abc}f_{abd} = Ca*\delta_{cd}, gives Ca = N for SU(N)
            //so Ca = 3 for SU(3).
    const value_type Cf;//Colour factor, from trace over generators:
            //\sum_a t^a_{ik}t^a_{kj} = Cf*\delta_{ij}
            //Cf = Tf*(N^2 - 1)/N for SU(N) theories, which gives 4/3 for SU(3)
            //NB - Cf does not appear in the SM at one loop level, so isn't
            //necessarily used below. We include it for completeness.
    const value_type zetaw;
    const value_type zetaz;
    couplingEvaluator<value_type> eval;








    //Standard beta functions (in t variable, not t_cl = log(phicl^2/mt^2))
    //Basic beta functions include the couplings up to g[5] = mt^2
    void beta_basic(const std::vector<value_type>& g,
                    std::vector<value_type>& beta_g,
                    const value_type t);
    void betap_basic(const std::vector<value_type>& g,
                     std::vector<value_type>& beta_g,
                     const std::vector<value_type>& gp,const value_type t);
    void betap_basic(const std::vector<value_type>& g,
                     std::vector<value_type>& beta_g,
                     const value_type t);
    void betapp_basic(const std::vector<value_type>& g,
                      std::vector<value_type>& beta_g,
                      const std::vector<value_type>& gp,
                      const std::vector<value_type>& gpp,
                      const value_type t);
    void betapp_basic(const std::vector<value_type>& g,
                      std::vector<value_type>& beta_g,
                      const value_type t);
    void betappp_basic(const std::vector<value_type>& g,
                       std::vector<value_type>& beta_g,
                       const std::vector<value_type>& gp,
                       const std::vector<value_type>& gpp,
                       const std::vector<value_type>& gppp,
                       const value_type t);
    void betappp_basic(const std::vector<value_type>& g,
                       std::vector<value_type>& beta_g,
                       const value_type t);
    void betapppp_basic(const std::vector<value_type>& g,
                        std::vector<value_type>& beta_g,
                        const std::vector<value_type>& gp,
                        const std::vector<value_type>& gpp,
                        const std::vector<value_type>& gppp,
                        const std::vector<value_type>& gpppp,
                        const value_type t);
    void betapppp_basic(const std::vector<value_type>& g,
                        std::vector<value_type>& beta_g,
                        const value_type t);
    //Extended beta functions including xi, Z, tcl, and t:
    void beta_extended(const std::vector<value_type>& g,
                       std::vector<value_type>& beta_g,
                       const value_type t);
    void betap_extended(const std::vector<value_type>& g,
                        std::vector<value_type>& beta_g,
                        const std::vector<value_type>& gp,const value_type t);
    void betap_extended(const std::vector<value_type>& g,
                        std::vector<value_type>& beta_g,
                        const value_type t);
    void betapp_extended(const std::vector<value_type>& g,
                         std::vector<value_type>& beta_g,
                         const std::vector<value_type>& gp,
                         const std::vector<value_type>& gpp,const value_type t);
    void betapp_extended(const std::vector<value_type>& g,
                         std::vector<value_type>& beta_g,
                         const value_type t);
    void betappp_extended(const std::vector<value_type>& g,
                          std::vector<value_type>& beta_g,
                          const std::vector<value_type>& gp,
                          const std::vector<value_type>& gpp,
                          const std::vector<value_type>& gppp,
                          const value_type t);
    void betappp_extended(const std::vector<value_type>& g,
                          std::vector<value_type>& beta_g,const value_type t);
    void betapppp_extended(const std::vector<value_type>& g,
                           std::vector<value_type>& beta_g,
                           const std::vector<value_type>& gp,
                           const std::vector<value_type>& gpp,
                           const std::vector<value_type>& gppp,
                           const std::vector<value_type>& gpppp,
                           const value_type t);
    void betapppp_extended(const std::vector<value_type>& g,
                           std::vector<value_type>& beta_g,const value_type t);

    //Beta functions involving the gravitational couplings:
    void beta_grav(const std::vector<value_type>& g,
                       std::vector<value_type>& beta_g,
                       const value_type t);
    void betap_grav(const std::vector<value_type>& g,
                        std::vector<value_type>& beta_g,
                        const std::vector<value_type>& gp,const value_type t);
    void betap_grav(const std::vector<value_type>& g,
                        std::vector<value_type>& beta_g,
                        const value_type t);
    void betapp_grav(const std::vector<value_type>& g,
                         std::vector<value_type>& beta_g,
                         const std::vector<value_type>& gp,
                         const std::vector<value_type>& gpp,const value_type t);
    void betapp_grav(const std::vector<value_type>& g,
                         std::vector<value_type>& beta_g,
                         const value_type t);
    void betappp_grav(const std::vector<value_type>& g,
                          std::vector<value_type>& beta_g,
                          const std::vector<value_type>& gp,
                          const std::vector<value_type>& gpp,
                          const std::vector<value_type>& gppp,
                          const value_type t);
    void betappp_grav(const std::vector<value_type>& g,
                          std::vector<value_type>& beta_g,const value_type t);
    void betapppp_grav(const std::vector<value_type>& g,
                           std::vector<value_type>& beta_g,
                           const std::vector<value_type>& gp,
                           const std::vector<value_type>& gpp,
                           const std::vector<value_type>& gppp,
                           const std::vector<value_type>& gpppp,
                           const value_type t);
    void betapppp_grav(const std::vector<value_type>& g,
                           std::vector<value_type>& beta_g,const value_type t);

    //Mass of the top quark (used as the reference scale):
    const value_type mt;
    //Ricci scalar:
    const value_type R;

    //Chooses what version of the beta functions we want, according to NTYPE:
    typedef void (betaSMdS_1loop<value_type>::*beta_pointer)
                    (const std::vector<value_type>&,std::vector<value_type>&,
                     const value_type);
    void beta_choice(const std::vector<value_type>& g,
                    std::vector<value_type>& beta_g,const value_type& t,
                    int NTYPE,
                    beta_pointer basic,beta_pointer extended,beta_pointer tcl,
                    beta_pointer grav);

    //Reference to a dtdtcl_jacobian, which actually implements the dtdtcl
    //transformation:
    dtdtcl_jacobian<value_type>& jacobian;
    //Access the member functions of the jacobian:
    value_type dt_dtcl(const value_type& tcl,
                                 const std::vector<value_type>& g,
                                 const std::vector<value_type>& dgdt)
    {
        return jacobian.dt_dtcl(tcl,g,dgdt);
    }
    value_type d2t_dtcl2(const value_type& tcl,
                                 const std::vector<value_type>& g,
                                 const std::vector<value_type>& dgdt,
                                 const std::vector<value_type>& d2gdt2,
                                 const value_type& dtdtcl)
    {
        return jacobian.d2t_dtcl2(tcl,g,dgdt,d2gdt2,dtdtcl);
    }
    value_type d3t_dtcl3(const value_type& tcl,
                                 const std::vector<value_type>& g,
                                 const std::vector<value_type>& dgdt,
                                 const std::vector<value_type>& d2gdt2,
                                 const std::vector<value_type>& d3gdt3,
                                 const value_type& dtdtcl,
                                 const value_type& d2tdtcl2)
    {
        return jacobian.d3t_dtcl3(tcl,g,dgdt,d2gdt2,d3gdt3,dtdtcl,d2tdtcl2);
    }
    value_type d4t_dtcl4(const value_type& tcl,
                                 const std::vector<value_type>& g,
                                 const std::vector<value_type>& dgdt,
                                 const std::vector<value_type>& d2gdt2,
                                 const std::vector<value_type>& d3gdt3,
                                 const std::vector<value_type>& d4gdt4,
                                 const value_type& dtdtcl,
                                 const value_type& d2tdtcl2,
                                 const value_type& d3tdtcl3)
    {
        return jacobian.d4t_dtcl4(tcl,g,dgdt,d2gdt2,d3gdt3,d4gdt4,dtdtcl,
                                  d2tdtcl2,d3tdtcl3);
    }






public:

    //Type of beta function we want to be the default (that which we will
    //get if we call operator() without specifying the type).
    int nType;
    /*
    nType = 0 -> basic beta functions (lambda, g1^2, g2^2, g3^2, yt^2, m^2)
    nType = 1 -> extended (includes xi and Z)
    nType = 2 -> using tcl variable, not t
    nType = 3 -> gravity, t variable
    nType = 4 -> gravity, tcl variable.
    */


    //constructors
    //When we want to specify both NTYPE (which beta function we are interested
    //in) and nDer (the derivative we are interested in)
    betaSMdS_1loop(value_type Mt,value_type r,value_type ZETAZ,value_type ZETAW,
                   dtdtcl_jacobian<value_type>& JACOBIAN,
                   int NTYPE,int nDer)
    : betaFunctionWithDerivatives(nDer),jacobian(JACOBIAN),
    _1p0(1), _16p0(16), _3p0(3), _8p0(8), _9p0(9), _m3p0(-3), _2p0(2),
    _m9p0(-9),
    _12p0(12), _6p0(6), _24p0(24), _m6p0(-6), _4p0(4), _27p0(27), _m27p0(-27),
    _72p0(72), _18p0(18), _m18p0(-18), _96p0(96), _36p0(36), _m24p0(-24),
    _25p0(25), _2304p0(2304), _10p0(10), _5p0(5), _11p0(11), _45p0(45),
    _1152p0(1152), _384p0(384), _288p0(288), _256p0(256), _m43p0(-43),
    _128p0(128),_15p0(15.0),
    _32p0(32), _m11p0(-11), _m17p0(-17), _m8p0(-8), _0p0(0.0),_m1p0(-1.0),
    Ng(3),Nc(3),Ca(3),Tf(0.5),Cf(value_type(4.0)/value_type(3.0)),
    mt(Mt), R(r),
    _40p0(40.0),_20p0(20.0),_200p0(200.0),_81p0(81.0),_50p0(50.0),_400p0(400.0),
    _m29p0(-29.0),_m2p0(-2.0),_68p0(68.0),_m49p0(49.0),_360p0(360.0),
    _5760p0(5760.0),_17p0(17.0),_180p0(180.0),zetaz(ZETAZ),zetaw(ZETAW),
    _m277p0(-277.0),_571p0(571.0),_2880p0(2880.0),_293p0(293.0),
    _23040p0(23040.0),_m293p0(-293.0)
    {
        this->nType = NTYPE;
        this->PI = boost::math::constants::pi<value_type>();
        loop_factor = _1p0/(_16p0*PI*PI);
    }
    betaSMdS_1loop(value_type Mt,value_type r,value_type ZETAZ,value_type ZETAW,
                   dtdtcl_jacobian<value_type>& JACOBIAN,int nDer)
    : betaFunctionWithDerivatives(nDer),jacobian(JACOBIAN),
    _1p0(1), _16p0(16), _3p0(3), _8p0(8), _9p0(9), _m3p0(-3), _2p0(2),
    _m9p0(-9),_15p0(15.0),
    _12p0(12), _6p0(6), _24p0(24), _m6p0(-6), _4p0(4), _27p0(27), _m27p0(-27),
    _72p0(72), _18p0(18), _m18p0(-18), _96p0(96), _36p0(36), _m24p0(-24),
    _25p0(25), _2304p0(2304), _10p0(10), _5p0(5), _11p0(11), _45p0(45),
    _1152p0(1152), _384p0(384), _288p0(288), _256p0(256), _m43p0(-43),
    _128p0(128),
    _32p0(32), _m11p0(-11), _m17p0(-17), _m8p0(-8), _0p0(0.0),_m1p0(-1.0),
    Ng(3),Nc(3),Ca(3),Tf(0.5),Cf(value_type(4.0)/value_type(3.0)),
    mt(Mt),R(r),
    _40p0(40.0),_20p0(20.0),_200p0(200.0),_81p0(81.0),_50p0(50.0),_400p0(400.0),
    _m29p0(-29.0),_m2p0(-2.0),_68p0(68.0),_m49p0(49.0),_360p0(360.0),
    _5760p0(5760.0),_17p0(17.0),_180p0(180.0),zetaz(ZETAZ),zetaw(ZETAW),
    _m277p0(-277.0),_571p0(571.0),_2880p0(2880.0),_293p0(293.0),
    _23040p0(23040.0),_m293p0(-293.0)
    {
        this->nType = 0;
        this->PI = boost::math::constants::pi<value_type>();
        loop_factor = _1p0/(_16p0*PI*PI);
    }
    betaSMdS_1loop(value_type Mt,value_type r,value_type ZETAZ,value_type ZETAW,
                   dtdtcl_jacobian<value_type>& JACOBIAN)
                    : betaFunctionWithDerivatives(),jacobian(JACOBIAN),
    _1p0(1), _16p0(16), _3p0(3), _8p0(8), _9p0(9), _m3p0(-3), _2p0(2),
    _m9p0(-9),_15p0(15.0),
    _12p0(12), _6p0(6), _24p0(24), _m6p0(-6), _4p0(4), _27p0(27), _m27p0(-27),
    _72p0(72), _18p0(18), _m18p0(-18), _96p0(96), _36p0(36), _m24p0(-24),
    _25p0(25), _2304p0(2304), _10p0(10), _5p0(5), _11p0(11), _45p0(45),
    _1152p0(1152), _384p0(384), _288p0(288), _256p0(256), _m43p0(-43),
    _128p0(128),
    _32p0(32), _m11p0(-11), _m17p0(-17), _m8p0(-8), _0p0(0.0),_m1p0(-1.0),
    Ng(3),Nc(3),Ca(3),Tf(0.5),Cf(value_type(4.0)/value_type(3.0)),
    mt(Mt),R(r),
    _40p0(40.0),_20p0(20.0),_200p0(200.0),_81p0(81.0),_50p0(50.0),_400p0(400.0),
    _m29p0(-29.0),_m2p0(-2.0),_68p0(68.0),_m49p0(49.0),_360p0(360.0),
    _5760p0(5760.0),_17p0(17.0),_180p0(180.0),zetaz(ZETAZ),zetaw(ZETAW),
    _m277p0(-277.0),_571p0(571.0),_2880p0(2880.0),_293p0(293.0),
    _23040p0(23040.0),_m293p0(-293.0)
    {
        this->nType = 0;
        this->PI = boost::math::constants::pi<value_type>();
        loop_factor = _1p0/(_16p0*PI*PI);
    }


    //beta functions:

    //Default derivative and type:
    //By default, this is the function that will be called if the beta function
    //is used in an ode solver.
    void operator()(const std::vector<value_type>& g,
                    std::vector<value_type>& beta_g,const value_type& t);
    //The following functions are necessary for this to be a non-abstract
    //instance of the betaFunction class (with and without derivatives).
    //Returns the beta function:
    void beta(const std::vector<value_type>& g,
              std::vector<value_type>& beta_g,const value_type& t);
    //Necessary for betaFunctionWithDerivatives. Returns the
    //Nth derivative of the beta functions:
    void betap(const std::vector<value_type>& g,
              std::vector<value_type>& beta_g,const value_type& t,
              int N);
    //Custom derivative, default type:
    void operator()(const std::vector<value_type>& g,
                    std::vector<value_type>& beta_g,const value_type& t,
                    int nDer);
    //beta_choice.
    //Custom derivative and type:

    void operator()(const std::vector<value_type>& g,
                    std::vector<value_type>& beta_g,const value_type& t,
                    int nDer,int NTYPE);

    //Transformation to tcl coordinates:

    //Beta functions in tcl co-ordinates:
    void beta_tcl(const std::vector<value_type>& g,
                  std::vector<value_type>& beta_g,const value_type tcl);
    void betap_tcl(const std::vector<value_type>& g,
                   std::vector<value_type>& beta_g,const value_type tcl);
    void betapp_tcl(const std::vector<value_type>& g,
                    std::vector<value_type>& beta_g,const value_type tcl);
    void betappp_tcl(const std::vector<value_type>& g,
                     std::vector<value_type>& beta_g,const value_type tcl);
    void betapppp_tcl(const std::vector<value_type>& g,
                      std::vector<value_type>& beta_g,const value_type tcl);


};



//------------------------------------------------------------------------------
/*
//abstract class for structure containing data about the sm couplings
//at a particular scale.
template<class value_type>
class beta_data
{
private:

public:
    betaFunctionSet<value_type>& betaFunctions;
    value_type txi;
    value_type t;
    std::vector<value_type> g;
    std::vector<value_type> dgdt;
    std::vector<value_type> d2gdt2;
    std::vector<value_type> d3gdt3;
    std::vector<value_type> d4gdt4;
    //constructor:
    beta_data(betaFunctionSet<value_type>& beta_set)
    : betaFunctions(beta_set)
    {
    }
    //Virtual mmember function that sets g and its derivatives to a
    //particular value (among other setup, as required by the specific
    //implementation).
    virtual void initialise(value_type TXI,std::vector<value_type>& G) = 0;
};
*/
//------------------------------------------------------------------------------
/*
//Struct which stores the couplings and their derivatives at a particular
//value of t, so we don't have to recompute these. References to this
//can be passed between classes and functions.
template<class value_type>
class beta_data_1_loop : public beta_data<value_type>
{
private:
    //running parameters:
    value_type txi;//txi = log(phi_tilde^2/mt^2), phi_tilde = Einstein frame
        //scalar field.
    value_type t;//standard running parameter for SM beta functions.

    const value_type mt;//top quark mass
    const value_type R;//Ricci scalar of background.
    betaSMdS_1loop<value_type> betaFunctions;
public:
    //numerical constants:
    const value_type _4p0;
    const value_type _2p0;
    const value_type _3p0;
    const value_type _12p0;
    const value_type _6p0;
    const value_type _1p0;
    const value_type _0p0;
    //Coefficients in the 1-loop potential:
    std::vector<value_type> ki;
    std::vector<value_type> kpi;
    std::vector<value_type> thetai;
    //Mass terms in logs of 1-loop potential:
    std::vector<value_type> Mi2;// Mi2 = ki*Z^2(t)*phicl^2 - kpi + thetai*R
    std::vector<value_type> pMi2pt;//Partial derivative of Mi2 wrt t at constant
        //phicl, ie, pMi2pt = [(dki/dt)*Z^2 + 2*ki*Z*Z']*phicl^2 - dkpi/dt
        //+ R*dthetai/dt
    //1st Derivatives of coefficients:
    std::vector<value_type> dkidt;
    std::vector<value_type> dkpidt;
    std::vector<value_type> dthetaidt;
    //2nd derivatives:
    std::vector<value_type> d2kidt2;
    std::vector<value_type> d2kpidt2;
    std::vector<value_type> d2thetaidt2;
    //3rd derivatives:
    std::vector<value_type> d3kidt3;
    std::vector<value_type> d3kpidt3;
    std::vector<value_type> d3thetaidt3;
    //4th derivatives:
    std::vector<value_type> d4kidt4;
    std::vector<value_type> d4kpidt4;
    std::vector<value_type> d4thetaidt4;
    //Constructor - takes the top quark mass, Ricci scalar, and a set
        //of beta functions as arguments:
    beta_data_1_loop(betaFunctionSet<value_type>& beta_set,value_type MT,
                     value_type r): beta_data(beta_set), mt(MT),R(r),_0p0(0.0),
                     _1p0(1.0),_2p0(2.0),_3p0(3.0),_4p0(4.0),_6p0(6.0),
                     _12p0(12.0)
    {

    }
    //Code to set SM couplings to a particular scale, for
    //a given txi, and other associated coefficients
    //needed for the 1-loop potential.
    void initialise(value_type TXI,std::vector<value_type>& G)
    {
        txi = TXI;//t_xi = log(phi_tilde^2/mt^2), phi_tilde = Einstein frame
            //field, mt = top quark mass.
        t = G[9];//Usual beta function running parameter
            //NB - we use the convention where t = log(mu^2/mt^2),
            //NOT t = log(mu/t), thus mu(t) = mt*exp(t/2).
        //Copy across and/or compute relevant data:
        this->g = G;//couplings.
        //Compute derivatives of couplings:
        this->dgdt = this->betaFunctions.beta(t,G);
        this->d2gdt2 = this->betaFunctions.betap(t,G);
        this->d3gdt3 = this->betaFunctions.betapp(t,G);
        this->d4gdt4 = this->betaFunctions.betappp(t,G);
        //Setup coefficients for 1 loop potential:
        value_type ki_array[] = {g[2]/this->_4p0,g[2]/this->_4p0,
                                 g[2]/this->_4p0,(g[2] + g[1])/this->_4p0,
                                 (g[2] + g[1])/this->_4p0,
                                 (g[2] + g[1])/this->_4p0,g[4]/this->_2p0,
                                 this->_3p0*g[0],g[0]};
        value_type kpi_array[] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,g[5],g[5]};
        value_type thetai_array[] = {this->_1p0/this->_12p0,
                                     this->_1p0/this->_12p0,
                                     -this->_1p0/this->_6p0,
                                     this->_1p0/this->_12p0,
                                     this->_1p0/this->_12p0,
                                     -this->_1p0/this->_6p0,
                                     this->_1p0/this->_12p0,
                                     g[6] - this->_1p0/this->_6p0,
                                     g[6] - this->_1p0/this->_6p0};
        this->ki.assign(ki_array,ki_array + 9);
        this->kpi.assign(kpi_array,kpi_array + 9);
        this->thetai.assign(thetai_array,thetai_array + 9);
        //Setup 1st derivatives of coefficients:
        value_type kip_array[] = {dgdt[2]/this->_4p0,dgdt[2]/this->_4p0,
                                 dgdt[2]/this->_4p0,
                                 (dgdt[2] + dgdt[1])/this->_4p0,
                                 (dgdt[2] + dgdt[1])/this->_4p0,
                                 (dgdt[2] + dgdt[1])/this->_4p0,
                                 dgdt[4]/this->_2p0,
                                 this->_3p0*dgdt[0],dgdt[0]};
        value_type kpip_array[] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,dgdt[5],
                                   dgdt[5]};
        value_type thetaip_array[] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,
                                     dgdt[6],dgdt[6]};
        this->dkidt.assign(kip_array,kip_array + 9);
        this->dkpidt.assign(kpip_array,kpip_array + 9);
        this->dthetaidt.assign(thetaip_array,thetaip_array + 9);
        //2nd derivatives:
        value_type kipp_array[] = {d2gdt2[2]/this->_4p0,d2gdt2[2]/this->_4p0,
                                 d2gdt2[2]/this->_4p0,
                                 (d2gdt2[2] + d2gdt2[1])/this->_4p0,
                                 (d2gdt2[2] + d2gdt2[1])/this->_4p0,
                                 (d2gdt2[2] + d2gdt2[1])/this->_4p0,
                                 d2gdt2[4]/this->_2p0,
                                 this->_3p0*d2gdt2[0],d2gdt2[0]};
        value_type kpipp_array[] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,
                                    d2gdt2[5],d2gdt2[5]};
        value_type thetaipp_array[] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,
                                     d2gdt2[6],d2gdt2[6]};
        this->d2kidt2.assign(kipp_array,kipp_array + 9);
        this->d2kpidt2.assign(kpipp_array,kpipp_array + 9);
        this->d2thetaidt2.assign(thetaipp_array,thetaipp_array + 9);
        //3rd derivatives:
        value_type kippp_array[] = {d3gdt3[2]/this->_4p0,d3gdt3[2]/this->_4p0,
                                 d3gdt3[2]/this->_4p0,
                                 (d3gdt3[2] + d3gdt3[1])/this->_4p0,
                                 (d3gdt3[2] + d3gdt3[1])/this->_4p0,
                                 (d3gdt3[2] + d3gdt3[1])/this->_4p0,
                                 d3gdt3[4]/this->_2p0,
                                 this->_3p0*d3gdt3[0],d3gdt3[0]};
        value_type kpippp_array[] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,
                                    d3gdt3[5],d3gdt3[5]};
        value_type thetaippp_array[] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,
                                     d3gdt3[6],d3gdt3[6]};
        this->d3kidt3.assign(kippp_array,kippp_array + 9);
        this->d3kpidt3.assign(kpippp_array,kpippp_array + 9);
        this->d3thetaidt3.assign(thetaippp_array,thetaippp_array + 9);
        4th derivatives:
        value_type kipppp_array[] = {d4gdt4[2]/this->_4p0,d4gdt4[2]/this->_4p0,
                                 d4gdt4[2]/this->_4p0,
                                 (d4gdt4[2] + d4gdt4[1])/this->_4p0,
                                 (d4gdt4[2] + d4gdt4[1])/this->_4p0,
                                 (d4gdt4[2] + d4gdt4[1])/this->_4p0,
                                 d4gdt4[4]/this->_2p0,
                                 this->_3p0*d4gdt4[0],d4gdt4[0]};
        value_type kpipppp_array[] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,
                                    d4gdt4[5],d4gdt4[5]};
        value_type thetaipppp_array[] = {_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,_0p0,
                                     d4gdt4[6],d4gdt4[6]};
        this->d4kidt4.assign(kipppp_array,kipppp_array + 9);
        this->d4kpidt4.assign(kpipppp_array,kpipppp_array + 9);
        this->d4thetaidt4.assign(thetaipppp_array,thetaipppp_array + 9);

        //Field renormalisation:
        value_type Z = exp(g[7]);
        //gamma function (derivative of Z):
        value_type gamma = -dgdt[7];
        //Classical field:
        value_type phicl = this->mt*exp(g[9]);
        //Mass coefficients:
        this->Mi2.assign(9,_0p0);
        this->pMi2pt.assign(9,_0p0);
        for(int i = 0;i<9;i++)
        {
            this->Mi2[i] = ki_array[i]*Z*Z*phicl*phicl - kpi_array[i]
                            + thetai_array[i]*R;
            this->pMi2pt[i] = (dkidt[i]- this->_2p0*ki[i]*gamma)*Z*Z*phicl*phicl
                                - dkpidt[i] + this->R*dthetaidt[i];
        }
    }
};
*/
//------------------------------------------------------------------------------

//Class to implement the Einstein-Jordan frame field transformation
//for the scalar field. This computer dtcl/dtxi where:
/*
tcl = log(phicl^2/mt^2)
txi = log(phi_tilde^2/mt^2)
phicl = classical field (Jordan frame)
phi_tilde = classical field (Einstein frame)

Computing this is equivalent to computing the derivative
dphicl/dphi_tilde.
*/
/*
template<class value_type>
class dtcldtxi
: public rational_function_shared_data<value_type,beta_data<value_type>&>
{
private:

    //std::vector<value_type> ki;
    //std::vector<value_type> kpi;
    //std::vector<value_type> thetai;
    //std::vector<value_type> Mi2v;




    //Needed constants:
    value_type h;
    value_type mt;
    //Transformation chosen between t(running parameter) and t_cl.
    //This is the relationship that chooses the renormalisation scale,
    //usually so that the quantum corrections are minimised.
    rational_function<value_type, beta_data<value_type>&>& dtdtcl_fun;

    //Shared data used by all functions:
    value_type phicl_shared;
    value_type mu_shared;
    value_type dtdtcl_shared;
    value_type d2tdtcl2_shared;
    value_type d3tdtcl3_shared;
    value_type d4tdtcl4_shared;
    value_type dtcldtxi_shared;
    value_type d2tcldtxi2_shared;
    value_type d3tcldtxi3_shared;

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
    //Cosmological constant, if we include it:
    ydata[10][i] = Lambda*Mp^2 (cosmological constant, in units of the Planck
                                mass)
    */
    /*
public:
    //Numerical constants (initialised by constructor)
    const value_type _0p0;
    const value_type _3p0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _1p0;
    const value_type _2p0;
    const value_type _m1p0;
    const value_type _8p0;
    const value_type _12p0;
    const value_type _m3p0;
    const value_type _m2p0;
    const value_type _24p0;
    const value_type _m6p0;
    const value_type _m12p0;
    const value_type _5p0;
    const value_type _18p0;
    const value_type _36p0;


    //Shared data computation. Avoids having to repeatedly evaluate
    //these functions every time A,B,Ap etc... is called. Instead, these
    //functions are evaluated once and A,B,Ap... share the result of
    //the computation. Which functions are needed ('sharing_level')
    //depends on how many derivatives we want to compute.
    void computeSharedData(beta_data<value_type>& x,int sharing_level);


    //Private Member functions. These should only be accessed by the
    //object itself, because they use information sharing to avoid
    //repetition of work.

    //Code omitted here for brevity - found in "potentials_code.h"
    value_type A(beta_data<value_type>& x);
    value_type Ap(beta_data<value_type>& x);
    value_type App(beta_data<value_type>& x);
    value_type Appp(beta_data<value_type>& x);
    value_type B(beta_data<value_type>& x);
    value_type Bp(beta_data<value_type>& x);
    value_type Bpp(beta_data<value_type>& x);
    value_type Bppp(beta_data<value_type>& x);

public:
    //constructor:
    dtcldtxi(value_type H,value_type MT,
             rational_function<value_type, beta_data<value_type>&>& DTDTCL)
              :  dtdtcl(DTDTCL),_3p0(3), _4p0(4), _6p0(6), _1p0(1), _2p0(2),
                 _m1p0(-1), _8p0(8), _12p0(12), _m3p0(-3),
                 _m2p0(-2), _24p0(24), _m6p0(-6), _m12p0(-12),
                 _5p0(5), _18p0(18), _36p0(36),_0p0(0.0)
    {
        this->h = H;
        this->mt = MT;
    }

};
*/
//------------------------------------------------------------------------------
/*
template< class value_type >
class couplingsFunctionEinstein
{
private:
    //Einstein frame transformation
    dtcldtxi<value_type>& dtcldtxi_fun;

};
*/
//------------------------------------------------------------------------------
//Jordan frame potential, SM, 1-loop
template<class value_type>
class SMpotential_dS_1loop : public SplinePotential<value_type>
{
public:

    //Required couplings data:
    const static int nCouplings = 23;//Number of couplings we have data about.
    //Coupling data:
    int nData;
    //Number of loop corrections:
    //nLoops = 9 gives the SM case without extra gravitational corrections
        //from R^2 terms.
    //nLoops = 12 includes the R^2 corrections
    const static int nLoops = 31;
    value_type* tcldata;
    value_type** ydata;//Couplings:
    /*
    NEW LIST:

    ydata[0][i] = lambda (Higgs self couplings)
    ydata[1][i] = g_1^2 (U(1) coupling)
    ydata[2][i] = g_2^2 (SU(2) coupling)
    ydata[3][i] = g_3^3 (SU(3) coupling)
    ydata[4][i] = m^2 (Higgs mass parameter)
    ydata[5][i] = y_e^2 (electron Yukawa coupling)
    ydata[6][i] = y_mu^2 (muon Yukawa coupling)
    ydata[7][i] = y_tau^2 (tau Yukawa coupling)
    ydata[8][i] = y_u^2 (up quark Yukawa coupling)
    ydata[9][i] = y_d^2 (down quark Yukawa coupling)
    ydata[10][i] = y_c^2 (charm quark Yukawa coupling)
    ydata[11][i] = y_s^2 (strange quark Yukawa coupling)
    ydata[12][i] = y_t^2 (top quark Yukawa coupling)
    ydata[13][i] = y_b^2 (bottom quark Yukawa coupling)

    ydata[14][i] = Z (Higgs field renormalisation)
    ydata[15][i] = tcl = log(phicl^2/mt^2)
    ydata[16][i] = t = log(mu^2/mt^2)

    ydata[17][i] = xi (Higgs-curvature non-minimal coupling)
    ydata[18][i] = V0 (cosmological constant energy density)
    ydata[19][i] = kappa (Ricci coupling, ie, [Planck mass]^2/2)
    ydata[20][i] = alpha1 (R^2 coupling)
    ydata[21][i] = alpha2 (R_{\mu\nu})^2 coupling
    ydata[22][i] = alpha3 (R_{\mu\nu\rho\sigma})^2 coupling.

    OLD LIST:
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
    ydata[10][i] = V0/Mp^4 (cosmological constant, in units of the Planck
                                mass)
    ydata[11][i] = kappa/Mp^2 (gravitational constant, in planck units
                                as evaluated at EW scale)
    ydata[12][i] = alpha1 R^2 coupling.
    ydata[13][i] = alpha2 Ricci^2 coupling
    ydata[14][i] = alpha3 Riemanna^2 coupling.
    */
    value_type** ypdata;//1st Derivatives of couplings
    value_type** yppdata;//2nd derivatives
    value_type** ypppdata;//3rd derivatives
    value_type** yppppdata;//4th derivatives
    value_type** ypppppdata;//5th derivatives

    //Beta functions. Need to be expressed as derivatives
    // wrt t_cl = log(phicl^2/mt^2)
    //betaFunctionSet<value_type>& betaFunctions;

    //Constants:
    const value_type yscale;//Scale that inputs are
        //expected to be in, in GeV.
    const value_type R;
    const value_type mt;
    const value_type zetaw;//Gauge fixing parameter.
    const value_type zetaz;
    std::vector<value_type> ni;
    std::vector<value_type> ci;
    //For R^2 couplings:
    std::vector<value_type> npi;
    value_type one_loop_factor;
    //Numerical constants:
    const value_type _8p0;
    const value_type _44p0;
    const value_type _48p0;
    const value_type _16p0;
    const value_type _m4p0;
    const value_type _m2p0;
    const value_type _64p0;
    const value_type PI;
    const value_type _m3p0;
    const value_type _18p0;

    const value_type _m6p0;
    const value_type _m12p0;
    const value_type _24p0;
    const value_type _36p0;
    const value_type _1p0;
    const value_type _12p0;
    const value_type _6p0;
    const value_type _0p0;
    const value_type _144p0;
    const value_type _2p0;


    //Switches:
    bool include_xi;//Whether to include non-minimal couplings
        //in the potential.
    bool include_grav;//Whether to include the cosmological constant
        //in the potential, and other gravity terms like R^2.
    bool include_R2;
    bool include_loops;
    bool splineCouplings;//Determines whether we apply interpolation
        //to the individual couplings (true), or to the potential (false).
        //Use false is interpolation is problematic, because this results
        //in a smoother potential (otherwise interpolation problems may occur


    coefficientArrayWithGrav<value_type> coeff;
    couplingEvaluator<value_type> eval;

    //Internal storage for potential at known points:
    std::vector<value_type> vdata;
    std::vector<value_type> vpdata;

    //Constructor:
    SMpotential_dS_1loop(value_type MT,value_type r,int NDATA,
                         value_type YSCALE,bool INC_XI,bool INC_GRAV,
                         bool INC_R2,bool INC_LOOP,
                         value_type* TCLDATA,
                         value_type** YDATA,
                         value_type** YPDATA,
                         value_type** YPPDATA,
                         value_type** YPPPDATA,
                         value_type** YPPPPDATA,
                         value_type** YPPPPPDATA,
                         value_type ZETAW,value_type ZETAZ,
                         bool SPLINE_COUPLINGS)
    : SplinePotential(TCLDATA,NDATA,MT),
    _8p0(8.0),_44p0(44.0),_48p0(48.0),_16p0(16.0),_m4p0(-4.0),
    _m2p0(-2.0),_64p0(64.0),PI(boost::math::constants::pi<value_type>()),
    _m3p0(-3.0),_18p0(18.0),_m6p0(-6.0),_m12p0(-12.0),_24p0(24.0),
    _36p0(36.0),_1p0(1.0),_12p0(12.0),_6p0(6.0),_0p0(0.0),_144p0(144.0),
    mt(MT),R(r),nData(NDATA),yscale(YSCALE),zetaz(ZETAZ),_2p0(2.0),
    coeff(ZETAZ,ZETAW),zetaw(ZETAW)
    {
        this->one_loop_factor = _1p0/(_64p0*PI*PI);
        fill_constant_arrays(this->ni,this->ci,this->npi);
        include_xi = INC_XI;
        include_grav = INC_GRAV;
        include_R2 = INC_R2;
        include_loops = INC_LOOP;
        splineCouplings = SPLINE_COUPLINGS;
        ydata = YDATA;
        ypdata = YPDATA;
        yppdata = YPPDATA;
        ypppdata = YPPPDATA;
        yppppdata = YPPPPDATA;
        ypppppdata = YPPPPPDATA;
        if(!splineCouplings)
        {
            //Need to internally store potential computed at a few different
            //points:
            for(int i = 0;i < NDATA;i++)
            {
                value_type phi = MT*exp(TCLDATA[i]/this->_2p0);
                std::vector<value_type> g (this->nCouplings,this->_0p0);
                std::vector<value_type> gp = g;
                for(int j = 0;j < this->nCouplings;j++)
                {
                    g[j] = ydata[j][i];
                    gp[j] = ypdata[j][i];
                }

                std::vector<value_type> Mi2;
                std::vector<value_type> logMi2;
                std::vector<value_type> Mi2p;
                std::vector<value_type> logMi2p;
                this->computeMi2p(phi,g,gp,Mi2,Mi2p,logMi2,logMi2p);

                //Compute potential:
                vdata.push_back(this->V(phi,g,logMi2,Mi2,include_xi,include_grav));

                //Note, to avoid oscillation, we shouldn't use the actual
                //derivative of the potential here!
                /*vpdata.push_back(this->Vp(phi,g,gp,logMi2,Mi2,logMi2p,Mi2p,
                                          include_xi,include_grav));*/
            }
            //Determine k factors using the pchip method:
            std::vector<value_type> x;
            x.assign(xdata,xdata + NDATA);
            pchipk(vpdata,x,vdata);
        }
    }

    //Required definitions of virtual functions:
    DLL_EXPORT value_type operator()(value_type y);
	DLL_EXPORT void operator()(value_type* y,value_type* arrayOut,int numel);
    //First derivatives:
	DLL_EXPORT value_type d(value_type y);
	DLL_EXPORT void d(value_type* y,value_type* arrayOut,int numel);
    //Second derivatives:
	DLL_EXPORT value_type d2(value_type y);
	DLL_EXPORT void d2(value_type* y,value_type* arrayOut,int numel);
    //Third derivatives:
	DLL_EXPORT value_type d3(value_type y);
	DLL_EXPORT void d3(value_type* y,value_type* arrayOut,int numel);
    //Fourth derivatives:
	DLL_EXPORT value_type d4(value_type y);
	DLL_EXPORT void d4(value_type* y,value_type* arrayOut,int numel);


	//Potential functions, for computing V, V', V'' etc. These will be wrapped
    //by the functions above.
    value_type V(value_type y,std::vector<value_type>& g,bool include_xi,
                 bool include_grav);
    value_type V(value_type y,std::vector<value_type>& g,
                 std::vector<value_type>& logMi2,std::vector<value_type>& Mi2v,
                 bool include_xi,bool include_grav);
    value_type Vp(value_type y,std::vector<value_type>& g,
                  std::vector<value_type>& gp,std::vector<value_type>& logMi2,
                  std::vector<value_type>& Mi2v,
                  std::vector<value_type>& logMi2p,
                  std::vector<value_type>& Mi2vp,bool include_xi,
                  bool include_grav);
    value_type Vp(value_type y,std::vector<value_type>& g,
                  std::vector<value_type>& gp,bool include_xi,
                  bool include_grav);
    value_type Vpp(value_type y,std::vector<value_type>& g,
                   std::vector<value_type>& gp,std::vector<value_type>& gpp,
                   std::vector<value_type>& logMi2,
                   std::vector<value_type>& Mi2v,
                   std::vector<value_type>& logMi2p,
                   std::vector<value_type>& Mi2vp,
                   std::vector<value_type>& logMi2pp,
                   std::vector<value_type>& Mi2vpp,bool include_xi,
                   bool include_grav);
    value_type Vpp(value_type y,std::vector<value_type>& g,
                   std::vector<value_type>& gp,std::vector<value_type>& gpp,
                   bool include_xi,bool include_grav);
    value_type Vppp(value_type y,std::vector<value_type>& g,
                    std::vector<value_type>& gp, std::vector<value_type>& gpp,
                    std::vector<value_type>& gppp,
                    std::vector<value_type>& logMi2,
                    std::vector<value_type>& Mi2v,
                    std::vector<value_type>& logMi2p,
                    std::vector<value_type>& Mi2vp,
                    std::vector<value_type>& logMi2pp,
                    std::vector<value_type>& Mi2vpp,
                    std::vector<value_type>& logMi2ppp,
                    std::vector<value_type>& Mi2vppp,bool include_xi,
                    bool include_grav);
    value_type Vppp(value_type y,std::vector<value_type>& g,
                    std::vector<value_type>& gp,std::vector<value_type>& gpp,
                    std::vector<value_type>& gppp,
                    bool include_xi,bool include_grav);
    value_type Vpppp(value_type y,std::vector<value_type>& g,
                     std::vector<value_type>& gp,std::vector<value_type>& gpp,
                     std::vector<value_type>& gppp,
                     std::vector<value_type>& gpppp,
                     std::vector<value_type>& logMi2,
                     std::vector<value_type>& Mi2v,
                     std::vector<value_type>& logMi2p,
                     std::vector<value_type>& Mi2vp,
                     std::vector<value_type>& logMi2pp,
                     std::vector<value_type>& Mi2vpp,
                     std::vector<value_type>& logMi2ppp,
                     std::vector<value_type>& Mi2vppp,
                     std::vector<value_type>& logMi2pppp,
                     std::vector<value_type>& Mi2vpppp,bool include_xi,
                     bool include_grav);
    value_type Vpppp(value_type y,std::vector<value_type>& g,
                     std::vector<value_type>& gp,std::vector<value_type>& gpp,
                     std::vector<value_type>& gppp,
                     std::vector<value_type>& gpppp, bool include_xi,
                     bool include_grav);



private:
    //The Standard model loop corrections contain many coefficients which
    //depend on the couplings. The functions below populate arrays with
    //these coefficients once the couplings have been chosen, so that we
    //can quickly sum over the loop corrections in a loop (the contributions
    //from each particle have the same form, with different coefficients):

    //Class for computing standard  model coupling-dependent coefficients:

    void populateCoefficientArrays
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& gp);
    //Derivatives of these coefficients (applies to nth derivatives - gp
    //contains the actual derivatives).
    void populateCoefficientArraysDerivatives
    (std::vector<value_type>& k,std::vector<value_type>& kp,
     std::vector<value_type>& theta,const std::vector<value_type>& gp);
    //Functions for computing Mi^2(phi) and its derivatives, an important
    //intermediary in the loop corrections. These have many arguments as higher
    //derivatives depend on the result of the lower derivatives. When we want
    //to use them, we store the results and supply them as arguments to the
    //higher derivative functions to avoid having to recompute them.
    //The different overloaded versions differ in whether they assume some data
    //is already computed (in which case they require it as an argument) or
    //whether they are expected to compute it.
    void computeMi2(value_type y,std::vector<value_type>& g,
                    std::vector<value_type>& Mi2,
                    std::vector<value_type>& logMi2);
    void computeMi2(value_type y,std::vector<value_type>& g,
                    std::vector<value_type>& Mi2,
                    std::vector<value_type>& logMi2,
                    std::vector<value_type>& kiarray,
                    std::vector<value_type>& kpiarray,
                    std::vector<value_type>& thetaiarray);
    void computeMi2p(value_type y, std::vector<value_type>& g,
                     std::vector<value_type>& gp,
                     std::vector<value_type>& Mi2,
                     std::vector<value_type>& logMi2,
                     std::vector<value_type>& Mi2p,
                     std::vector<value_type>& logMi2p,
                     std::vector<value_type>& kiarray,
                     std::vector<value_type>& kpiarray,
                     std::vector<value_type>& thetaiarray,
                     std::vector<value_type>& kiparray,
                     std::vector<value_type>& kpiparray,
                     std::vector<value_type>& thetaiparray);
    void computeMi2p(value_type y, std::vector<value_type>& g,
                     std::vector<value_type>& gp,std::vector<value_type>& Mi2,
                     std::vector<value_type>& logMi2,
                     std::vector<value_type>& Mi2p,
                     std::vector<value_type>& logMi2p);
    void computeMi2pp(value_type y, std::vector<value_type>& g,
                      std::vector<value_type>& gp,std::vector<value_type>& gpp,
                      std::vector<value_type>& Mi2,
                      std::vector<value_type>& logMi2,
                      std::vector<value_type>& Mi2p,
                      std::vector<value_type>& logMi2p,
                      std::vector<value_type>& Mi2pp,
                      std::vector<value_type>& logMi2pp,
                      std::vector<value_type>& kiarray,
                      std::vector<value_type>& kpiarray,
                      std::vector<value_type>& thetaiarray,
                      std::vector<value_type>& kiparray,
                      std::vector<value_type>& kpiparray,
                      std::vector<value_type>& thetaiparray,
                      std::vector<value_type>& kipparray,
                      std::vector<value_type>& kpipparray,
                      std::vector<value_type>& thetaipparray);
    void computeMi2pp(value_type y, std::vector<value_type>& g,
                      std::vector<value_type>& gp,std::vector<value_type>& gpp,
                      std::vector<value_type>& Mi2,
                      std::vector<value_type>& logMi2,
                      std::vector<value_type>& Mi2p,
                      std::vector<value_type>& logMi2p,
                      std::vector<value_type>& Mi2pp,
                      std::vector<value_type>& logMi2pp);
    void computeMi2ppp(value_type y, std::vector<value_type>& g,
                       std::vector<value_type>& gp,std::vector<value_type>& gpp,
                       std::vector<value_type>& gppp,
                       std::vector<value_type>& Mi2,
                       std::vector<value_type>& logMi2,
                       std::vector<value_type>& Mi2p,
                       std::vector<value_type>& logMi2p,
                       std::vector<value_type>& Mi2pp,
                       std::vector<value_type>& logMi2pp,
                       std::vector<value_type>& Mi2ppp,
                       std::vector<value_type>& logMi2ppp,
                       std::vector<value_type>& kiarray,
                       std::vector<value_type>& kpiarray,
                       std::vector<value_type>& thetaiarray,
                       std::vector<value_type>& kiparray,
                       std::vector<value_type>& kpiparray,
                       std::vector<value_type>& thetaiparray,
                       std::vector<value_type>& kipparray,
                       std::vector<value_type>& kpipparray,
                       std::vector<value_type>& thetaipparray,
                       std::vector<value_type>& kippparray,
                       std::vector<value_type>& kpippparray,
                       std::vector<value_type>& thetaippparray);
    void computeMi2ppp(value_type y, std::vector<value_type>& g,
                       std::vector<value_type>& gp,
                       std::vector<value_type>& gpp,
                       std::vector<value_type>& gppp,
                       std::vector<value_type>& Mi2,
                       std::vector<value_type>& logMi2,
                       std::vector<value_type>& Mi2p,
                       std::vector<value_type>& logMi2p,
                       std::vector<value_type>& Mi2pp,
                       std::vector<value_type>& logMi2pp,
                       std::vector<value_type>& Mi2ppp,
                       std::vector<value_type>& logMi2ppp);
    void computeMi2pppp(value_type y, std::vector<value_type>& g,
                        std::vector<value_type>& gp,
                        std::vector<value_type>& gpp,
                        std::vector<value_type>& gppp,
                        std::vector<value_type>& gpppp,
                        std::vector<value_type>& Mi2,
                        std::vector<value_type>& logMi2,
                        std::vector<value_type>& Mi2p,
                        std::vector<value_type>& logMi2p,
                        std::vector<value_type>& Mi2pp,
                        std::vector<value_type>& logMi2pp,
                        std::vector<value_type>& Mi2ppp,
                        std::vector<value_type>& logMi2ppp,
                        std::vector<value_type>& Mi2pppp,
                        std::vector<value_type>& logMi2pppp,
                        std::vector<value_type>& kiarray,
                        std::vector<value_type>& kpiarray,
                        std::vector<value_type>& thetaiarray,
                        std::vector<value_type>& kiparray,
                        std::vector<value_type>& kpiparray,
                        std::vector<value_type>& thetaiparray,
                        std::vector<value_type>& kipparray,
                        std::vector<value_type>& kpipparray,
                        std::vector<value_type>& thetaipparray,
                        std::vector<value_type>& kippparray,
                        std::vector<value_type>& kpippparray,
                        std::vector<value_type>& thetaippparray,
                        std::vector<value_type>& kipppparray,
                        std::vector<value_type>& kpipppparray,
                        std::vector<value_type>& thetaipppparray);
    void computeMi2pppp(value_type y, std::vector<value_type>& g,
                        std::vector<value_type>& gp,
                        std::vector<value_type>& gpp,
                        std::vector<value_type> gppp,
                        std::vector<value_type> gpppp,
                        std::vector<value_type>& Mi2,
                        std::vector<value_type>& logMi2,
                        std::vector<value_type>& Mi2p,
                        std::vector<value_type>& logMi2p,
                        std::vector<value_type>& Mi2pp,
                        std::vector<value_type>& logMi2pp,
                        std::vector<value_type>& Mi2ppp,
                        std::vector<value_type>& logMi2ppp,
                        std::vector<value_type>& Mi2pppp,
                        std::vector<value_type>& logMi2pppp);
};


//------------------------------------------------------------------------------
//New version of the SM potential, designed for use in de-Sitter space,
//Including the effects of xi running.
/*
template< class value_type >
class SMHiggsPotential_dS_1loop_einstein : public SplinePotential< value_type >
{
private:

    //Stored data about couplings. These are used to compute
    //the spline interpolation:
    value_type** ydata;//Couplings
    value_type** kdata;//Derivatives of couplings.
    //value_type* xdata;//t_xi = log(phi_tilde^2/mt^2)
    //Reference list of couplings:
    /*
    ydata[0] - lambda (Higgs self coupling)
    ydata[1] - g1^2 (U(1) coupling)
    ydata[2] - g2^2 (SU(2) coupling)
    ydata[3] - g3^2 (SU(3) coupling)
    ydata[4] - yt^2 (Top quark Yukawa coupling)
    ydata[5] - m^2 (Higgs tachyonic mass term)
    ydata[6] - xi (Higgs-curvature coupling)
    ydata[7] - Z (Higgs field renormalisation) - technically the RHS of this
                  is a gamma function, but we solve it with the beta functions)
    //Additionally, we have to solve for the non-beta functions:
    ydata[8] = t_cl = log(phicl^2/mt^2)
    ydata[9] = t = log(mu^2/mt^2)
    */
    /*
    //Constants:
    std::vector<value_type> ni;
    std::vector<value_type> ci;
    const value_type PI;// = boost::math::constants::pi<value_type>();
    const value_type one_loop_factor;

    //Numerical constants:
    const value_type _3p0;
    const value_type _4p0;
    const value_type _6p0;
    const value_type _1p0;
    const value_type _2p0;
    const value_type _m1p0;
    const value_type _8p0;
    const value_type _12p0;
    const value_type _m3p0;
    const value_type _m2p0;
    const value_type _24p0;
    const value_type _m6p0;
    const value_type _m12p0;
    const value_type _5p0;
    const value_type _18p0;
    const value_type _36p0;
    const value_type _64p0;
    const value_type _m4p0;
    //More constants:
    const value_type _11p0;
    const value_type _22p0;
    const value_type _44p0;
    const value_type _48p0;
    const value_type _16p0;



public:
    value_type yscale;//The potential as it stands is formulated in GeV,
    //usually. Set the yscale variable to tell the
    //potential what units we expect supplied arguments y to be in.
    //The potential will then attempt to return V(yscale*y)/yscale^4,
    //V'(yscale*y)/yscale^3
    //V''(yscale*y)/yscale^2, V'''(yscale*y)/yscale and V''''(yscale*y)
    //etc... yscale should be, in GeV, whatever y = 1 corresponds to.
    //The normalisation
    //is because W(y) = V(yscale*y)/yscale^4 is the effective potential
    //appearing in the equations of motion when y = phi/yscale, as opposed
    //to V(phi).



    //Potential and its derivatives:


    //Einstein frame potential:
	DLL_EXPORT value_type operator()(value_type y);
	DLL_EXPORT void operator()(value_type* y,value_type* arrayOut,int numel);
    //First derivatives:
	DLL_EXPORT value_type d(value_type y);
	DLL_EXPORT void d(value_type* y,value_type* arrayOut,int numel);
    //Second derivatives:
	DLL_EXPORT value_type d2(value_type y);
	DLL_EXPORT void d2(value_type* y,value_type* arrayOut,int numel);
    //Third derivatives:
	DLL_EXPORT value_type d3(value_type y);
	DLL_EXPORT void d3(value_type* y,value_type* arrayOut,int numel);
    //Fourth derivatives:
	DLL_EXPORT value_type d4(value_type y);
	DLL_EXPORT void d4(value_type* y,value_type* arrayOut,int numel);


	//Jordan frame potential:




    /*
    //Compute a spline for array (derivative) nArrayToUse, at point y:
	DLL_EXPORT value_type SinglePointSpline(value_type& y,int nArrayToUse);
	DLL_EXPORT value_type SinglePointSpline(value_type& y,int nArrayToUse,
                                            int nSpline);
	//Compute a spline of the  nArrayToUse'th derivative, for each of the set
	//of points indicated by array y (with numel elements), putting output
	//into array arrayOut (with numel elements):
	DLL_EXPORT void MultiPointSpline(value_type* y,value_type* arrayOut,
                                     int numel,int nArrayToUse);
    */
    /*
    //Function to implement renormalisation scale choice. Computes:
    //dt/dt_cl where t = renormalisation parameter, and t_cl = log(phicl^2/mt^2)
    //where phicl is the classical field
    rational_function<value_type, beta_data<value_type>&>& dtdtcl_fun;

    //Function to implement transformation between Einstein and Jordan frame.
    //Computes dt_cl/dt_xi where t_cl is as above, and
    //t_xi = log(phi_tilde^2/mt^2), phi_tilde = Einstein frame field.
    dtcldtxi<value_type> dtcldtxi_fun;




    /*
    //Definitions:
    y - refers to einstein frame field. This is the input.
    yc - refers to Jordan frame (classical) field.
    txi = log(y*y/m*m);
    tc = log(yc*yc/m*m);
    R - Ricci scalar
    Z - field renormalisation
    gamma - anomalous dimension
    alpha,beta - constants involved in the scale choice:
        mu^2(t) = alpha*Z^2*yc^2 + beta*R
    xi - non-minimal coupling.
    H - unit choice for energy scales.
    Mp - reduced planck mass
    h = H/Mp
    t = original scale parameter. mu(t) = mt*exp(t)
    */
    /*
    //Constructors
    //Setup function, used in various constructors:
    void setup(value_type* XDATA,value_type** YDATA,
                           value_type** KDATA,int NDATA,value_type m)
    {
        //xdata = XDATA;
        ydata = YDATA;//couplings
        kdata = KDATA;//couplings derivatives wrt t_xi
        //Ndata = NDATA;
        //M = m;
        value_type ni_array = {_2p0 , _6p0 , -_2p0 , _1p0 , _3p0,-_1p0,-_12p0,
                               _1p0,_3p0};
        value_type ci_array = { _3p0/_2p0 , _5p0/_6p0 , _3p0/_2p0 , _3p0/_2p0 ,
                               _5p0/_6p0 , _3p0/_2p0 , _3p0/_2p0 , _3p0/_2p0 ,
                               _3p0/_2p0};
        this->ni.assign(ni_array,ni_array + 9);
        this->ci.assign(ci_array,ci_array + 9);
    }
    //Normal constructor:
    SMHiggsPotential_dS_1loop_einstein(value_type* XDATA,value_type** YDATA,
                           value_type** KDATA,int NDATA,value_type m,
                           rational_function<value_type,
                                            beta_data<value_type>&>& DTDTCL_fun)
    : SplinePotential(XDATA,NDATA,m),dtdtcl_fun(DTDTCL_fun),
    dtcldtxi_fun(yscale,m,DTDTCL_fun),
    PI(boost::math::constants::pi<value_type>()),
    one_loop_factor(value_type(64.0)*PI*PI),_3p0(3), _4p0(4), _6p0(6), _1p0(1),
    _2p0(2), _m1p0(-1), _8p0(8), _12p0(12),
    _m3p0(-3), _m2p0(-2), _24p0(24), _m6p0(-6), _m12p0(-12), _5p0(5), _18p0(18),
    _36p0(36), _64p0(64), _m4p0(-4),_11p0(11.0),_22p0(22.0),_44p0(44.0),
    _48p0(48.0),_16p0(16.0)
    {
        yscale = value_type(1.0);
        //Setup maximum and minimum values:
        this->hasMaximum = true;
        this->hasMinimum = true;
        //const value_type _2p0 = value_type(2.0);
        this->Maximum = (m/yscale)*exp(XDATA[NDATA - 1]/_2p0);
        this->Minimum = (m/yscale)*exp(XDATA[0]/_2p0);

        //Initialise arrays of constants:
        this->setup(XDATA,YDATA,KDATA,NDATA,m);
    }
    //Overloaded version where yscale is specified:
    SMHiggsPotential_dS_1loop_einstein(value_type* XDATA,value_type** YDATA,
                           value_type** KDATA,int NDATA,value_type m,
                           rational_function<value_type,
                                            beta_data<value_type>&>& DTDTCL_fun,
                           value_type YSCALE)
    : SplinePotential(XDATA,NDATA,m),dtdtcl_fun(DTDTCL_fun),
    dtcldtxi_fun(yscale,m,DTDTCL_fun),
    PI(boost::math::constants::pi<value_type>()),
    one_loop_factor(value_type(64.0)*PI*PI),_3p0(3), _4p0(4), _6p0(6), _1p0(1),
    _2p0(2), _m1p0(-1), _8p0(8), _12p0(12),
    _m3p0(-3), _m2p0(-2), _24p0(24), _m6p0(-6), _m12p0(-12), _5p0(5), _18p0(18),
    _36p0(36), _64p0(64), _m4p0(-4),_11p0(11.0),_22p0(22.0),_44p0(44.0),
    _48p0(48.0),_16p0(16.0)
    {
        this->yscale = YSCALE;
        //Setup maximum and minimum values:
        this->hasMaximum = true;
        this->hasMinimum = true;
        const value_type _2p0 = value_type(2.0);
        this->Maximum = (m/yscale)*exp(XDATA[NDATA - 1]/_2p0);
        this->Minimum = (m/yscale)*exp(XDATA[0]/_2p0);

        //Initialise arrays of constants:
        this->setup(XDATA,YDATA,KDATA,NDATA,m);
    }

};
*/
#endif
//------------------------------------------------------------------------------
//==============================================================================
template<class value_type>
class log_potential : public potential< value_type >
{
private:
    value_type lambda_0;
    value_type a;
    value_type b;
    value_type M;
    const value_type _4p0;
    const value_type _0p0;
    const value_type _2p0;
    const value_type _7p0;
    const value_type _6p0;
    const value_type _3p0;
    const value_type _13p0;
    const value_type _10p0;
    const value_type _26p0;
    const value_type _25p0;
    const value_type _70p0;
    const value_type _50p0;
public:
    //Constructors:
    log_potential(value_type L,value_type A,value_type B,value_type m)
    : _4p0(4.0), _0p0(0.0), _2p0(2.0), _7p0(7.0), _6p0(6.0), _3p0(3.0),
    _13p0(13.0), _10p0(10.0), _26p0(26.0), _25p0(25.0), _70p0(70.0), _50p0(50.0)
    {
        lambda_0 = L;
        a = A;
        b = B;
        M = m;
    }
    //Constructor which allows the use to specify the scale of expected inputs
    //The potential is then returned in units of that scale, rather than in the
    //original units.
    log_potential(value_type L,value_type A,value_type B,value_type m,
                  value_type scale)
    : _4p0(4.0), _0p0(0.0), _2p0(2.0), _7p0(7.0), _6p0(6.0), _3p0(3.0),
    _13p0(13.0), _10p0(10.0), _26p0(26.0), _25p0(25.0), _70p0(70.0), _50p0(50.0)
    {
        lambda_0 = L;
        a = A;
        b = B;
        M = m/scale;
    }
    DLL_EXPORT value_type operator()(value_type x);
	DLL_EXPORT void operator()(value_type* x,value_type* arrayOut,
                                       int numel);
	DLL_EXPORT value_type d(value_type x);
	DLL_EXPORT void d(value_type* x,value_type* arrayOut,int numel);
	DLL_EXPORT value_type d2(value_type x);
	DLL_EXPORT void d2(value_type* x,value_type* arrayOut,int numel);
	DLL_EXPORT value_type d3(value_type x);
	DLL_EXPORT void d3(value_type* x,value_type* arrayOut,int numel);
	DLL_EXPORT value_type d4(value_type x);
	DLL_EXPORT void d4(value_type* x,value_type* arrayOut,int numel);
};
//==============================================================================
#endif
