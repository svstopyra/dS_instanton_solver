#ifndef NEWTON_RAPHSON_INCLUDE
#define NEWTON_RAPHSON_INCLUDE

//DLL stuff - determines whether the header is being used by the dll itself (in which case we are exporting) or
//an application which references it (in which case we are importing):
#ifdef DLL_EXPORT
#define DLL_EXPORT __declspec(dllexport)
#else
#ifdef STATIC_EXPORT
#define DLL_EXPORT //If we specifically do not want to use either of the dll options, and just want to compile the functions normally.
#else
#define DLL_EXPORT __declspec(dllimport)
#endif
#endif


//c++ version of the globally convergent Newton-Raphson algorithm we are
//using for finding instantons.
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
//#include "errors.h"
#include <vector>
#include <boost\math\special_functions.hpp>

//Define abstract classes which we then derive functions and jacobians to use with the NR solver from:
template< class value_type >
class functionNR
{
public:
	virtual DLL_EXPORT void operator()(value_type* inPointer, value_type* outPointer) = 0;
};
template< class value_type >
class jacobianNR
{
public:
	void DLL_EXPORT operator()(std::vector< std::vector< value_type > >& J, std::vector< value_type >& x);
	virtual DLL_EXPORT void compute_jacobian(std::vector< std::vector< value_type > >& J, std::vector< value_type >& x) = 0;
};

//Abstract base class for Newton-Raphson solvers:
template<class value_type >
class gradientDescent
{
typedef Eigen::Matrix< value_type, Eigen::Dynamic , Eigen::Dynamic > mMatrix;
typedef Eigen::Matrix< value_type, Eigen::Dynamic, 1 > mVector;
public:
    //Broyden method flag:
    bool useApproxJacobian;
    bool usingBroyden;
    //Solver properties:
    value_type derivative_step;
    value_type precision;
    value_type stepMin;
    value_type stepMax;
    value_type alpha;//backtrack paramter.
	functionNR< value_type >* F_pointer;//Pointer to the function we are finding
        //a zero of.
    //bounded solver flag and bounds:
    bool bounded;
    std::vector< value_type > lowerBound;
    std::vector< value_type > upperBound;
    int nComponents;
    //Reason the solver stopped:
    int stopReason;
    //Global check for convergence:
    bool isConverging;
    //no. of backtracks we have performed:
    int backtracks;
    //Method to use for derivatives:
    int method;
    //Reciprocol condition number limit:
    value_type rconditionLimit;
    //Whether the user supplied a function for the exact Jacobian or not:
    //bool userSuppliedJacobian;
    //Pointer to user-supplied Jacobian:
	jacobianNR< value_type >* jacobian;
    //Constructors:
	gradientDescent(value_type DerivStep,value_type Prec, value_type step_min,value_type step_max,value_type Alpha = value_type(1e-4),value_type rcondLim = std::numeric_limits< value_type >::epsilon(),int METHOD = 4);
	gradientDescent(value_type DerivStep,value_type Prec, value_type step_min,value_type step_max,bool Broyden,value_type Alpha = value_type(1e-4),value_type rcondLim = std::numeric_limits< value_type >::epsilon(),int METHOD = 4);
	gradientDescent(std::vector< value_type >& lower,std::vector< value_type >& upper, value_type DerivStep,value_type Prec, value_type step_min,value_type step_max,value_type Alpha = value_type(1e-4),value_type rcondLim = std::numeric_limits< value_type >::epsilon(),int METHOD = 4);
	gradientDescent(std::vector< value_type >& lower,std::vector< value_type >& upper, value_type DerivStep,value_type Prec, value_type step_min,value_type step_max,bool Broyden,value_type Alpha = value_type(1e-4),value_type rcondLim = std::numeric_limits< value_type >::epsilon(),int METHOD = 4);
    //Member functions:
    //Accessing solver properties:
	DLL_EXPORT void bound(std::vector< value_type >& lower,std::vector< value_type >& upper);
	DLL_EXPORT void setPrecision(value_type newPrecision);
	DLL_EXPORT void setStepLimits(value_type newStepMin,value_type newStepMax);
	DLL_EXPORT void setDerivativeStep(value_type newDerivativeStep);
	DLL_EXPORT void useBroyden(bool useBroy);
	DLL_EXPORT void RefreshSolver();
    //NR solver itself:
    virtual DLL_EXPORT int operator()(value_type* x0,int n, functionNR< value_type >* f,value_type* output,value_type* Fout , jacobianNR< value_type >* jac = NULL) =0;
    //Subroutines used by NR solver:
	DLL_EXPORT bool checkWithinPrecision(mVector& Fnew);
	DLL_EXPORT void solveSystem(mMatrix& A,mVector& x,mVector& b);
	DLL_EXPORT value_type rcond(mMatrix& mat);
	DLL_EXPORT mVector F(mVector& x);
	DLL_EXPORT mMatrix Jacobian(mVector& x,int nMethod,mVector& F0);
	DLL_EXPORT mVector checkBounded(mVector& x, mVector& dx);
	DLL_EXPORT mVector imposeMaxStep(mVector& dx);
	DLL_EXPORT mVector ArrayToVector(std::vector< value_type >& in);
	DLL_EXPORT std::vector< value_type > VectorToArray(mVector& in);
};

template<class value_type >
gradientDescent< value_type >::gradientDescent(value_type DerivStep, value_type Prec, value_type step_min, value_type step_max, value_type Alpha, value_type rcondLim, int METHOD)
{
	derivative_step = DerivStep;
	precision = Prec;
	stepMin = step_min;
	stepMax = step_max;
	bounded = false;
	useApproxJacobian = false;
	alpha = Alpha;
	isConverging = false;
	method = METHOD;
	usingBroyden = useApproxJacobian;
	rconditionLimit = rcondLim;
	//userSuppliedJacobian = false;
}
template<class value_type >
gradientDescent< value_type >::gradientDescent(value_type DerivStep, value_type Prec, value_type step_min, value_type step_max, bool Broyden, value_type Alpha, value_type rcondLim, int METHOD)
{
	derivative_step = DerivStep;
	precision = Prec;
	stepMin = step_min;
	stepMax = step_max;
	bounded = false;
	useApproxJacobian = Broyden;
	alpha = Alpha;
	isConverging = false;
	method = METHOD;
	usingBroyden = useApproxJacobian;
	rconditionLimit = rcondLim;
	//userSuppliedJacobian = false;
}
template<class value_type >
gradientDescent< value_type >::gradientDescent(std::vector< value_type >& lower, std::vector< value_type >& upper, value_type DerivStep, value_type Prec, value_type step_min, value_type step_max, value_type Alpha, value_type rcondLim, int METHOD)
{
	derivative_step = DerivStep;
	precision = Prec;
	stepMin = step_min;
	stepMax = step_max;
	bounded = true;
	useApproxJacobian = false;
	lowerBound = lower;
	upperBound = upper;
	alpha = Alpha;
	isConverging = false;
	method = METHOD;
	usingBroyden = useApproxJacobian;
	rconditionLimit = rcondLim;
	//userSuppliedJacobian = false;
}
template<class value_type >
gradientDescent< value_type >::gradientDescent(std::vector< value_type >& lower, std::vector< value_type >& upper, value_type DerivStep, value_type Prec, value_type step_min, value_type step_max, bool Broyden, value_type Alpha, value_type rcondLim, int METHOD)
{
	derivative_step = DerivStep;
	precision = Prec;
	stepMin = step_min;
	stepMax = step_max;
	bounded = true;
	useApproxJacobian = Broyden;
	lowerBound = lower;
	upperBound = upper;
	alpha = Alpha;
	isConverging = false;
	method = METHOD;
	usingBroyden = useApproxJacobian;
	rconditionLimit = rcondLim;
	//userSuppliedJacobian = false;
}



//Basic newton-raphson class with backtracking and reflectng boundaries.
template<class value_type >
class NRSolver : public gradientDescent< value_type >
{
typedef Eigen::Matrix< value_type, Eigen::Dynamic , Eigen::Dynamic > mMatrix;
typedef Eigen::Matrix< value_type, Eigen::Dynamic, 1 > mVector;
private:
    //Whether we are using reflecting boundaries. Assumed true by default.
    bool reflecting;
    //Whether each component is sitting at its boundary. This needs to be set
    //by the Newton-Raphson solver when it is called.
    std::vector<bool> atUpperBoundary;
    std::vector<bool> atLowerBoundary;
public:
    //Constructors:
	NRSolver(value_type DerivStep,value_type Prec, value_type step_min,
            value_type step_max,value_type Alpha = value_type(1e-4),
            value_type rcondLim = std::numeric_limits< value_type >::epsilon(),
            bool use_reflection = true,int METHOD = 4);
	NRSolver(value_type DerivStep,value_type Prec, value_type step_min,
            value_type step_max,bool Broyden,
            value_type Alpha = value_type(1e-4),
            value_type rcondLim = std::numeric_limits< value_type >::epsilon(),
            bool use_reflection = true,int METHOD = 4);
	NRSolver(std::vector< value_type >& lower,
            std::vector< value_type >& upper,
            value_type DerivStep,value_type Prec,
            value_type step_min,value_type step_max,
            value_type Alpha = value_type(1e-4),
            value_type rcondLim = std::numeric_limits< value_type >::epsilon(),
            bool use_reflection = true,int METHOD = 4);
	NRSolver(std::vector< value_type >& lower,
            std::vector< value_type >& upper, value_type DerivStep,
            value_type Prec, value_type step_min,value_type step_max,
            bool Broyden,value_type Alpha = value_type(1e-4),
            value_type rcondLim = std::numeric_limits< value_type >::epsilon(),
            bool use_reflection = true,int METHOD = 4);
	DLL_EXPORT int operator()(value_type* x0,int n,functionNR< value_type >* f,value_type* output,value_type* Fout , jacobianNR< value_type >* jac = NULL);
	DLL_EXPORT bool checkBacktrackConvergence(mVector& Fnew,value_type l,value_type lmin,mMatrix& J,mVector& x);
	DLL_EXPORT mVector backtrack(mVector& x,mVector& dx,mVector& Fnew,mVector& Fold,mMatrix& J);
    //Bounds checking, with reflection:
	DLL_EXPORT mVector checkBoundedReflecting(mVector& x, mVector& dx);
};
template<class value_type >
NRSolver< value_type >::NRSolver(value_type DerivStep, value_type Prec, value_type step_min,
	value_type step_max, value_type Alpha,
	value_type rcondLim,
	bool use_reflection, int METHOD)
	: gradientDescent< value_type >(DerivStep, Prec, step_min, step_max, Alpha,
		rcondLim, METHOD)
{
	this->reflecting = use_reflection;
}
template<class value_type >
NRSolver< value_type >::NRSolver(value_type DerivStep, value_type Prec, value_type step_min,
	value_type step_max, bool Broyden,
	value_type Alpha,
	value_type rcondLim,
	bool use_reflection, int METHOD)
	: gradientDescent< value_type >(DerivStep, Prec, step_min, step_max, Broyden,
		Alpha, rcondLim, METHOD)
{
	this->reflecting = use_reflection;
}
template<class value_type >
NRSolver< value_type >::NRSolver(std::vector< value_type >& lower,
	std::vector< value_type >& upper,
	value_type DerivStep, value_type Prec,
	value_type step_min, value_type step_max,
	value_type Alpha,
	value_type rcondLim,
	bool use_reflection, int METHOD)
	: gradientDescent< value_type >(lower, upper, DerivStep, Prec, step_min,
		step_max, Alpha, rcondLim, METHOD)
{
	this->reflecting = use_reflection;
}
template<class value_type >
NRSolver< value_type >::NRSolver(std::vector< value_type >& lower,
	std::vector< value_type >& upper, value_type DerivStep,
	value_type Prec, value_type step_min, value_type step_max,
	bool Broyden, value_type Alpha,
	value_type rcondLim,
	bool use_reflection, int METHOD)
	: gradientDescent< value_type >(lower, upper, DerivStep, Prec, step_min,
		step_max, Broyden, Alpha, rcondLim, METHOD)
{
	this->reflecting = use_reflection;
}
/*
//For future implementation, if necessary:
//Version of gradient descent with simulated Annealing:
template<class function, class value_type>
class annealingSolver : public gradientDescent< function , value_type >
{
typedef Eigen::Matrix< value_type, Eigen::Dynamic , Eigen::Dynamic > mMatrix;
typedef Eigen::Matrix< value_type, Eigen::Dynamic, 1 > mVector;
public:
    //Constructors:
    annealingSolver(value_type DerivStep,value_type Prec, value_type step_min,
            value_type step_max,value_type Alpha = value_type(1e-4),
            value_type rcondLim = value_type(1e-15))
            : gradientDescent(DerivStep,Prec,step_min,step_max,Alpha,rcondLim){}
    annealingSolver(value_type DerivStep,value_type Prec, value_type step_min,
            value_type step_max,bool Broyden,
            value_type Alpha = value_type(1e-4),
            value_type rcondLim = value_type(1e-15))
            : gradientDescent(DerivStep,Prec,step_min,step_max,Broyden,
                    Alpha,rcondLim){}
    annealingSolver(std::vector< value_type >& lower,
            std::vector< value_type >& upper,
            value_type DerivStep,value_type Prec,
            value_type step_min,value_type step_max,
            value_type Alpha = value_type(1e-4),
            value_type rcondLim = value_type(1e-15))
            : gradientDescent(lower,upper,DerivStep,Prec,step_min,
                    step_max,Alpha,rcondLim){}
    annealingSolver(std::vector< value_type >& lower,
            std::vector< value_type >& upper, value_type DerivStep,
            value_type Prec, value_type step_min,value_type step_max,
            bool Broyden,value_type Alpha = value_type(1e-4),
            value_type rcondLim = value_type(1e-15))
            : gradientDescent(lower,upper,DerivStep,Prec,step_min,
                    step_max,Broyden,Alpha,rcondLim){}
    //Solver:
    int operator()(value_type* x0,int n,function& f,value_type* output,value_type* Fout)
    {
    }
};
 */



#endif
