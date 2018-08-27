#ifndef NEWTON_RAPHSON_CODE_H
#define NEWTON_RAPHSON_CODE_H
#include "NewtonRaphson.h"

//pre-processor macro to avoid typing (since typedef will only work within a give
//context of the term 'value_type':
#define mMatrix Eigen::Matrix< value_type, Eigen::Dynamic , Eigen::Dynamic >
#define mVector Eigen::Matrix< value_type, Eigen::Dynamic, 1 >

template< class value_type >
void jacobianNR< value_type >::operator()(std::vector< std::vector< value_type > >& J, std::vector< value_type >& x)
{
	compute_jacobian(J, x);
}

//Abstract base class for Newton-Raphson solvers:
//Member functions:
//Accessing solver properties:template<class value_type>
template<class value_type >
void gradientDescent< value_type >::bound(std::vector< value_type >& lower, std::vector< value_type >& upper)
{
	this->bounded = true;
	this->lowerbound = lower;
	this->upperbound = upper;
}
template<class value_type >
void gradientDescent< value_type >::setPrecision(value_type newPrecision)
{
	this->precision = newPrecision;
}
template<class value_type >
void gradientDescent< value_type >::setStepLimits(value_type newStepMin, value_type newStepMax)
{
	this->stepMin = newStepMin;
	this->stepMax = newStepMax;
}
template<class value_type >
void gradientDescent< value_type >::setDerivativeStep(value_type newDerivativeStep)
{
	this->derivative_step = newDerivativeStep;
}
template<class value_type >
void gradientDescent< value_type >::useBroyden(bool useBroy)
{
	this->useApproxJacobian = useBroy;
}
template<class value_type >
void gradientDescent< value_type >::RefreshSolver()
{
	//Sets the internal variables back to their defaults:
	this->StopReason = 0;
	this->backtracks = 0;
}
//NR solver itself:
//Subroutines used by NR solver:
template<class value_type >
bool gradientDescent< value_type >::checkWithinPrecision(mVector& Fnew)
{
	//Checks whether all components of F are sufficiently close to 0:
	bool withinPrecision = true;
	mVector Fabs = Fnew.cwiseAbs();
	for (int i = 0; i < int(Fnew.size()); i++)
	{
		withinPrecision = withinPrecision && (Fabs(i) < this->precision);
		if (!withinPrecision)
		{
			break;
		}
	}
	return withinPrecision;
}
template<class value_type >
void gradientDescent< value_type >::solveSystem(mMatrix& A, mVector& x, mVector& b)
{
	//Solves the system Ax = b for x, putting the result into x, via a
	//variety of methods.
	x = A.lu().solve(b);
}
template<class value_type >
value_type gradientDescent< value_type >::rcond(mMatrix& mat)
{
	//Returns the reciprocol of the condition number of the matrix supplied:

	//Singular value decomposition:
	Eigen::JacobiSVD<mMatrix> svd(mat);
	mVector sings = svd.singularValues();
	int N = int(sings.size());
	//Get the rcond number:
	return (sings(N - 1) / sings(0));
}
template<class value_type >
mVector gradientDescent< value_type >::F(mVector& x)
{
	//The supplied function needs only the following:
	/*  - Takes input data as an array of value_types.
	*  - Outputs result of function as an array of value_types.
	*The internals of how it deals with this information are up
	*to the user. This subroutine will allow us to interpret the result
	*as an mVector for use in the Newton-Raphson code.
	*
	*/
	mVector vToReturn(nComponents);
	//the function pointed to be F_pointer is assumed to have two
	//arguments:
	// (*F_pointer)(value_type* inPointer,value_type* outPointer)
	//It is crucial that this function does not overflow any of
	//the arrays it uses! It should take input data as an array of
	//value types, and output the result directly into the mVector:

	//value_type* outPointer = vToReturn.data();
	//value_type* inPointer = x.data();

	//Create a copy of x as a temporary. This guards against a badly written F
	//Accidentally corrupting x and messing up the solver routine.
	mVector x_copy = x;
	//Could also remove this step if x is very large and
	//optimisation is required.

	//F, otherwise, should not have a return type:
	(*F_pointer)(x_copy.data(), vToReturn.data());
	//Return the mVector:
	return vToReturn;
}
template<class value_type >
mMatrix gradientDescent< value_type >::Jacobian(mVector& x, int nMethod, mVector& F0)
{
	//PARAMETERS:
	/*
	*x - point at which to evaluate the Jacobian
	*F - function to evaluate the Jacobian of
	*method - what sort of derivative to use:
	*  method = 1 -> forwards difference
	*  method = 2 -> backwards difference
	*  method = 3 -> central difference (slower)
	*F0 - value of the function at x, so we don't have to recompute it.
	*
	*/
	//Setup Jacobian matrix.
	int N = int(x.size());
	mMatrix J(N, N);
	if (!jacobian)
	{
		//Use a finite difference scheme to try and estimate the Jacobian.
		//Get the number of variables. It is required that F
		//have as many output components as x has components, so
		//that the Jacobian is square and can be inverted.

		//Fill Jacobian by taking approximate derivatives:
		mMatrix dx(N, 1);
		value_type step;
		switch (method)
		{
		case 1:
			//Forwards difference:
			for (int i = 0; i < N; i++)
			{
				dx = mMatrix::Zero(N, 1);
				step = (this->derivative_step)*abs(x(i));
				dx(i) = step;
				mVector xpdx = x + dx;
				J.block(0, i, N, 1) = (F(xpdx) - F0) / step;
			}
			break;
		case 2:
			//Backwards difference:
			for (int i = 0; i < N; i++)
			{
				dx = mMatrix::Zero(N, 1);
				step = (this->derivative_step)*abs(x(i));
				dx(i) = step;
				mVector xmdx = x - dx;
				J.block(0, i, N, 1) = (F0 - F(xmdx)) / step;
			}
			break;
		case 3:
			//Central difference:
			for (int i = 0; i < N; i++)
			{
				dx = mMatrix::Zero(N, 1);
				step = (this->derivative_step)*abs(x(i));
				dx(i) = step;
				mVector xpdx = x + dx;
				mVector xmdx = x - dx;
				J.block(0, i, N, 1) = (F(xpdx) - F(xmdx)) / (value_type(2.0)*step);
			}
			break;
		case 4:
			//Adpative method. Attempts to choose the optimum value of the
			//step size assuming that the derivative_step given is the relative error in
			//F:
			for (int i = 0; i < N; i++)
			{
				dx = mMatrix::Zero(N, 1);
				//Initial guess at step:
				value_type d3Fdx = value_type(1.0);//Approximation of 3rd derivative
				value_type epsF = (this->derivative_step);//Scale of error in F:
				step = pow(abs(epsF*((value_type(5.0)*
					F0(i)) / (value_type(2.0)*d3Fdx))),
					value_type(1.0) / value_type(3.0));
				//Don't allow too large steps:
				if (abs(step) > value_type(0.1)*abs(x(i)))
				{
					step = value_type(0.1)*abs(x(i));
				}
				int nTimes = 3;//Maximum number of improvement loops.
				int counter = 0;//Counts number of improvement loops.
				std::cout << "Chosen step after " << counter
					<< " iterations: dx(" << i << ") = "
					<< step.str(50, std::ios_base::scientific)
					<< std::endl;
				//Now check whether this gives a consistent result:
				value_type ratio = value_type(0.1);
				value_type improved_step = step;
				//Displaced vectors:
				mVector xpdx(N, 1);
				mVector xmdx(N, 1);
				mVector xp2dx(N, 1);
				mVector xm2dx(N, 1);
				mVector Fxpdx(N, 1);
				mVector Fxmdx(N, 1);
				mVector Fxp2dx(N, 1);
				mVector Fxm2dx(N, 1);
				do
				{
					//Compute F with this step:
					step = improved_step;
					dx(i) = improved_step;
					xpdx = x + dx;
					xmdx = x - dx;
					xp2dx = x + value_type(2.0)*dx;
					xm2dx = x - value_type(2.0)*dx;
					Fxpdx = F(xpdx);
					Fxmdx = F(xmdx);
					Fxp2dx = F(xp2dx);
					Fxm2dx = F(xm2dx);
					//Estimate the third derivative, taking the average
					//over all components for a given dx:
					d3Fdx = value_type(0.0);
					for (int j = 0; j < N; j++)
					{
						d3Fdx += (Fxp2dx(j) - Fxm2dx(j) -
							value_type(2.0)*(Fxpdx(j) -
								Fxmdx(j))) / (value_type(2.0)*step*step*step);
					}
					d3Fdx /= value_type(N);
					//Optimum step:
					improved_step = pow(abs(epsF*((value_type(5.0)*
						F0(i)) / (value_type(2.0)*d3Fdx))),
						value_type(1.0) / value_type(3.0));
					counter++;
					std::cout << "Chosen step after " << counter
						<< " iterations: dx(" << i << ") = "
						<< improved_step.str(50, std::ios_base::scientific)
						<< std::endl;
				} while ((abs(improved_step / step) < ratio ||
					abs(improved_step / step) > value_type(1.0) / ratio)
					&& counter < nTimes);
				if (counter == nTimes)
				{
					std::cout << "Warning: optimum step-size not found.\n";
				}
				//Hopefully we have a good enough step size. Now use it
				//to take a central different approximation of the Jacobian:
				J.block(0, i, N, 1) = (Fxpdx - Fxmdx) / (value_type(2.0)*
					improved_step);
			}
		default:
			throw "Error - invalid method requested for computing Jacobian.";
		}
	}
	else
	{
		//Call the exact user supplied Jacobian. First define vectors to
		//hold all the elements:
		std::vector< std::vector< value_type > > Jexact;
		std::vector< value_type > x0;
		Jexact.reserve(nComponents);
		std::vector< value_type > zeros(nComponents, value_type(0.0));
		for (int i = 0; i < nComponents; i++)
		{
			x0.push_back(x(i));
			Jexact.push_back(zeros);
		}
		//Now call the user supplied Jacobian:
		(*(this->jacobian))(Jexact, x0);//will obviously fail if pointer is not of this form.
										//Take the resulting data and create a Jacobian matrix from it:
		for (int i = 0; i < this->nComponents; i++)
		{
			for (int j = 0; j < this->nComponents; j++)
			{
				J(i, j) = Jexact[i][j];
			}
		}
	}
	return J;
}
template<class value_type >
mVector gradientDescent< value_type >::checkBounded(mVector& x, mVector& dx)
{
	//Applies a battery of tests to the proposed step to ensure that it
	//does not cross a (finite) boundary, contains no NaN values (which would
	//indicate an error - this can happen sometimes) or contains infinite
	//inputs x or proposed steps dx, either of which indicates a serious error
	//with which the Newton-Step cannot continue.
	if (bounded)
	{
		int N = int(x.size());
		bool boundsBroken = false;
		//First ignore any bounds which are given as infinite:
		mVector xFiniteLowerBound;
		int Mlower = 0;
		mVector dxFiniteLowerBound;
		mVector xFiniteUpperBound;
		int Mupper = 0;
		mVector dxFiniteUpperBound;
		mVector lowerBoundFinite;
		mVector upperBoundFinite;
		int nOvershot = 0;
		mVector dxnew;
		mVector shortenFactorVector;

		for (int i = 0; i < N; i++)
		{
			//Filter out any non-numbers - these should generate errors.
			if (boost::math::isnan(x(i)))
			{
				this->stopReason = 7;
				return dx;
			}
			else if (boost::math::isnan(dx(i)))
			{
				this->stopReason = 8;
				return dx;
			}
			//Filter out directions with infinite bounds as we don't care about them:
			if (!boost::math::isinf(lowerBound[i]))
			{
				Mlower++;
				xFiniteLowerBound.conservativeResize(Mlower);
				xFiniteLowerBound(Mlower - 1) = x(i);
				dxFiniteLowerBound.conservativeResize(Mlower);
				dxFiniteLowerBound(Mlower - 1) = dx(i);
				lowerBoundFinite.conservativeResize(Mlower);
				lowerBoundFinite(Mlower - 1) = lowerBound[i];
			}
			if (!boost::math::isinf(upperBound[i]))
			{
				Mupper++;
				xFiniteUpperBound.conservativeResize(Mupper);
				xFiniteUpperBound(Mupper - 1) = x(i);
				dxFiniteUpperBound.conservativeResize(Mupper);
				dxFiniteUpperBound(Mupper - 1) = dx(i);
				upperBoundFinite.conservativeResize(Mupper);
				upperBoundFinite(Mupper - 1) = upperBound[i];
			}
			//Check for anything that breaks a bound.
			boundsBroken = boundsBroken || x(i) + dx(i) > upperBound[i] || x(i) + dx(i) < lowerBound[i];
			//Add an element to our shorten factor vector, to be used later:
			if ((x(i) + dx(i) > upperBound[i] && !boost::math::isinf(upperBound[i])) && dx(i) != value_type(0.0))
			{
				nOvershot++;
				shortenFactorVector.conservativeResize(nOvershot);
				shortenFactorVector(nOvershot - 1) = (upperBound[i] - x(i)) / (dx(i));
			}
			if ((x(i) + dx(i) < lowerBound[i] && !boost::math::isinf(lowerBound[i])) && dx(i) != value_type(0.0))
			{
				nOvershot++;
				shortenFactorVector.conservativeResize(nOvershot);
				shortenFactorVector(nOvershot - 1) = (lowerBound[i] - x(i)) / (dx(i));
			}
		}
		if (boundsBroken && nOvershot  > 0)
		{
			//Find the maximal overshooting of the boundary:
			value_type shortenFactor = shortenFactorVector.cwiseAbs().maxCoeff();
			//Now shorten dx so it remains within this:
			dxnew = (value_type(0.99)*shortenFactor)*dx;
			//Check whether this step is too short:
			value_type magNew = sqrt(dxnew.dot(dxnew));
			if (magNew < stepMin)
			{
				std::cout << "Warning: Newton search exited due to hitting boundary.\n";
				this->stopReason = 3;
			}
			return dxnew;
		}
		else
		{
			return dx;
		}
	}
	else
	{
		return dx;
	}
}
template<class value_type >
mVector gradientDescent< value_type >::imposeMaxStep(mVector& dx)
{
	value_type mag = sqrt(dx.dot(dx));
	mVector dxnew;
	if (mag > this->stepMax)
	{
		dxnew = (value_type(0.95)*(this->stepMax) / mag)*dx;
	}
	else
	{
		dxnew = dx;
	}
	return dxnew;
}
template<class value_type >
mVector gradientDescent< value_type >::ArrayToVector(std::vector< value_type >& in)
{
	//Converts an array (actually a std::vector) to an mVector.
	int N = int(in.size());
	mVector a(N);
	for (int i = 0; i < N; i++)
	{
		a(i) = in[i];
	}
	return a;
}
template<class value_type >
std::vector< value_type > gradientDescent< value_type >::VectorToArray(mVector& in)
{
	//Converts an mVector to an array (std::vector)
	int N = int(in.size());
	std::vector< value_type > a;
	a.reserve(N);
	for (int i = 0; i < N; i++)
	{
		a[i] = in(i);
	}
	return a;
}

//Basic newton-raphson class with backtracking and reflectng boundaries.
template<class value_type >
int NRSolver< value_type >::operator()(value_type* x0, int n, functionNR< value_type >* f, value_type* output, value_type* Fout, jacobianNR< value_type >* jac)
{
	/*
	*PARAMETERS:
	*x0 - pointer to first element of array containing the
	*      initial values. Could be an array or a vector.
	*      We pass this as a pointer because we don't want to
	*      pin the code down to one type of storage mechanism.
	*
	*n - number of components in the array described by x0.
	*
	*/

	//mexPrintf("1\n");
	//mexEvalString("pause(.001);");

	//Assign the function and the number of components:
	this->nComponents = n;
	this->F_pointer = f;
	this->jacobian = jac;

	//Allocate the boundary boolean vector:
	(this->atUpperBoundary).assign(this->nComponents, false);
	(this->atLowerBoundary).assign(this->nComponents, false);

	//Loop iteration counter:
	int loopCounter = 0;
	int loopMax = 200;

	//Rest backtrack counter and stopReason:
	this->backtracks = 0;
	this->stopReason = 0;

	//Logical test for loop exit:
	bool foundZero = false;

	//Convert the initial value into an mVector:
	mVector x(n);
	for (int i = 0; i < n; i++)
	{
		x(i) = x0[i];
	}
	//mexPrintf("2\n");
	//mexEvalString("pause(.001);");
	//Get the initial value of the function:
	mVector Fnew = this->F(x);
	//mexPrintf("3\n");
	//mexEvalString("pause(.001);");
	//Start off by numerically computing the Jacobian:
	mMatrix J = this->Jacobian(x, this->method, Fnew);
	//Output the Jacobian so we can see what it is doing:
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			std::cout << J(i, j).str(5, std::ios_base::scientific) << " ";
		}
		std::cout << std::endl;
	}
	this->usingBroyden = false;//Our approximation always starts off as
							   //exact (or exact as a finite difference can be - 'exact' here
							   // is relative to the Broyden method)
							   //mexPrintf("4\n");
							   //mexEvalString("pause(.001);");
							   //Pre-declaration of variables needed out of array scope:
	mVector Fold;

	while (!foundZero)
	{
		//Store last round of variables:
		Fold = Fnew;

		//Check for singular Jacobian:
		value_type rcondition = this->rcond(J);
		std::cout << "rcond = " << rcondition.str(50, std::ios_base::scientific) << std::endl;
		std::cout << "rcond = " << this->rconditionLimit.str(50, std::ios_base::scientific) << std::endl;
		if (rcondition < this->rconditionLimit || boost::math::isnan(rcondition))
		{
			std::cout << ".\nPossible convergence. Checking... " << std::endl;
			if (this->checkWithinPrecision(Fnew))
			{
				this->stopReason = 1;
				std::cout << "Converged to a zero within requested precision.";
			}
			else
			{
				//Try re-computing the Jacobian if we were using a Broyden
				//approximation. It may have got away from us!
				if (this->usingBroyden)
				{
					J = this->Jacobian(x, this->method, Fnew);
					this->usingBroyden = false;
					//Re-compute and check rcond:
					rcondition = this->rcond(J);
					if (rcondition < this->rconditionLimit || boost::math::isnan(rcondition))
					{
						//Still bad. Bisect and try again!
					}
				}
				else
				{
					//Even our exact Jacobian is badly conditioned!!
					this->stopReason = 4;
					std::cout << "Warning - spurious convergence detected. Try different initial conditions or adjust requested matchigng precision.\n";
				}
			}
			break;
		}

		//Check for non-finite inputs:
		bool nonFinite = false;
		for (int i = 0; i < x.size(); i++)
		{
			nonFinite = nonFinite || !boost::math::isfinite(x(i));
			if (nonFinite)
			{
				break;
			}
		}
		if (nonFinite)
		{
			this->stopReason = 7;
			break;
		}

		//Checks complete. Proceed with the Newton Step:
		//Aim to find dx = -J^{-1}Fnew, but inverting J is slow. We
		//simply solve J*dx = -Fnew instead:
		//mexPrintf("5\n");
		//mexEvalString("pause(.001);");
		mVector dx(x.size());
		mVector mFnew = -Fnew;
		this->solveSystem(J, dx, mFnew);

		//Check we remain within the requested bounds using this step:
		std::cout << "Step after solving Jacobian: " << std::endl;
		for (int i = 0; i < n; i++)
		{
			std::cout << "dx[" << i << "] = " << dx(i).str(50, std::ios_base::scientific) << std::endl;
		}
		dx = this->checkBoundedReflecting(x, dx);
		std::cout << "Step after reflecting off boundary: " << std::endl;
		for (int i = 0; i < n; i++)
		{
			std::cout << "dx[" << i << "] = " << dx(i).str(50, std::ios_base::scientific) << std::endl;
		}
		dx = this->imposeMaxStep(dx);
		std::cout << "Step after imposing max step size: " << std::endl;
		for (int i = 0; i < n; i++)
		{
			std::cout << "dx[" << i << "] = " << dx(i).str(50, std::ios_base::scientific) << std::endl;
		}
		//Stop if for whatever reason, solving cannot continue:
		if (this->stopReason != 0)
		{
			break;
		}
		//mexPrintf("6\n");
		//mexEvalString("pause(.001);");
		//Now evaluate F at the new point:
		mVector xnew = x + dx;
		Fnew = this->F(xnew);

		//Apply backtracking in case this step isn't good enough:
		dx = this->backtrack(x, dx, Fnew, Fold, J);
		std::cout << "Step after backtracking: " << std::endl;
		for (int i = 0; i < n; i++)
		{
			std::cout << "dx[" << i << "] = " << dx(i).str(50, std::ios_base::scientific) << std::endl;
		}
		//mexPrintf("7\n");
		//mexEvalString("pause(.001);");
		//Translate to the new point:
		x = x + dx;

		//Detect convergence:
		this->isConverging = this->checkWithinPrecision(Fnew);
		if (dx.cwiseAbs().maxCoeff() < this->stepMin)
		{
			if (this->isConverging)
			{
				this->stopReason = 1;
			}
			else
			{
				this->stopReason = 2;
			}
			break;
		}
		//mexPrintf("8\n");
		//mexEvalString("pause(.001);");
		//Recompute the Jacobian if necessary:
		if (!(this->useApproxJacobian))
		{
			J = this->Jacobian(x, this->method, Fnew);
		}
		else
		{
			//Perform a Broyden update:
			mVector dF = Fnew - Fold;

			J = J + ((dF - J*dx)*(dx.adjoint())) / (dx.dot(dx));
			this->usingBroyden = true;
		}
		//mexPrintf("9\n");
		//mexEvalString("pause(.001);");
		//Finish this loop iteration:
		loopCounter++;
		foundZero = loopCounter > loopMax || this->isConverging || this->stopReason != 0;
	}

	//Report why the loop exited:
	if (!(this->isConverging) && loopCounter > loopMax)
	{
		this->stopReason = 5;
	}
	else if (this->isConverging)
	{
		this->stopReason = 1;
	}
	//Transfer data to output arrays:
	for (int i = 0; i < n; i++)
	{
		output[i] = x(i);
		Fout[i] = Fnew(i);
	}
	return this->stopReason;
}
template<class value_type >
bool NRSolver< value_type >::checkBacktrackConvergence(mVector& Fnew, value_type l, value_type lmin, mMatrix& J, mVector& x)
{
	//Check for convergences during the backtrack subroutine.
	if (l < lmin)
	{
		//Check for convergence
		if (this->checkWithinPrecision(Fnew))
		{
			this->stopReason = 1;
			return true;
		}
		else
		{
			//Indicates spurious convergence (or too much requested precision
			//and we are at the limit - as close as we can get to zero
			//Try recomputing the Jacobian:
			if (this->usingBroyden)
			{
				this->backtracks = 0;
				J = this->Jacobian(x, this->method, Fnew);
				this->usingBroyden = false;
				return false;
			}
			else
			{
				//We have become trapped in a local minimum. Report this
				//and stop the search:
				this->stopReason = 2;
				return true;
			}
		}
	}
	else
	{
		return false;
	}
}
template<class value_type >
mVector NRSolver< value_type >::backtrack(mVector& x, mVector& dx, mVector& Fnew, mVector& Fold, mMatrix& J)
{
	//Backtracks the proposed newton step in order to ensure we always go
	//in a decreasing direction (guarantees global convergence, even if
	//it is to a spurious zero).
	/*PARAMETERS
	*x - current position
	*dx - proposed step
	*F - function we are searching for a zero of
	*Fnew - function evaluated at x + dx (avoids us having to recomputes it
	*      which is potentially expensive
	*Fold - function evaluated at x (provided for the same reason as Fnew)
	*J - Jacobian matrix.
	*/
	mVector dxnew;

	//Start backtrack cycle at zero:
	this->backtracks = 0;

	//Decide on a smallest allowed value of lambda, the rescaling parameter
	//for the step:
	value_type dxMax = dx.cwiseAbs().maxCoeff();
	value_type lmin = (this->stepMin) / dxMax;

	//Chosen step, lambda:
	value_type l = value_type(1.0);

	//Intermediates:
	value_type l1;
	value_type l2;
	value_type gp0;
	value_type g0;
	value_type g1;
	value_type g2;

	//Modulus square function:
	value_type fnew = value_type(0.5)*(Fnew.dot(Fnew));
	value_type fold = value_type(0.5)*(Fold.dot(Fold));

	//Store value of Fold.dot(J*dx) so we don't keep repeating ourselves:
	value_type FJdx = Fold.dot(J*dx);

	//Decide whether we need backtracking or not:
	bool backtrackNeeded = fnew > fold + l*(this->alpha)*(FJdx);
	if (!backtrackNeeded)
	{
		dxnew = dx;
	}
	else
	{
		while (backtrackNeeded)
		{
			if (this->backtracks == 0)
			{
				//Starting guess:
				l1 = value_type(1.0);

				//Parameters describing polynomial approximation of function
				//along the direction of the step:
				gp0 = -value_type(2.0)*fnew;
				g0 = fold;
				g1 = fnew;

				//Step to of this polynomial:
				l2 = -gp0 / (value_type(2.0)*(g1 - g0 - gp0));
				//Cutoff to ensure we don't take too large or too small steps:
				if (l2 > value_type(0.5)*l1)
				{
					l2 = value_type(0.5)*l1;
				}
				else if (l2 < value_type(0.1)*l1)
				{
					l2 = value_type(0.1)*l1;
				}
				std::cout << "lnew = " << l2.str(50, std::ios_base::scientific) << std::endl;

				//New step:
				dxnew = l2*dx;
				std::cout << "Step after " << this->backtracks + 1 << " backtracks : " << std::endl;
				for (int i = 0; i < this->nComponents; i++)
				{
					std::cout << "dx[" << i << "] = " << dx(i).str(50, std::ios_base::scientific) << std::endl;
				}

				//Try the function at this new point:
				mVector xnew = x + dxnew;
				//mexPrintf("x[0] = %s\nx[1] = %s\nx[2] = %s\n",x(0).str(50,std::ios_base::scientific).c_str(),x(1).str(50,std::ios_base::scientific).c_str(),x(2).str(50,std::ios_base::scientific).c_str());
				//mexPrintf("dxnew[0] = %s\ndxnew[1] = %s\ndxnew[2] = %s\n",dxnew(0).str(50,std::ios_base::scientific).c_str(),dxnew(1).str(50,std::ios_base::scientific).c_str(),dxnew(2).str(50,std::ios_base::scientific).c_str());
				Fnew = this->F(xnew);
				fnew = value_type(0.5)*(Fnew.dot(Fnew));
				g2 = fnew;
				this->backtracks = 1;
				l = l2;


				//Check for lambda being too small - possibly indicates convergence:
				this->isConverging = this->checkBacktrackConvergence(Fnew, l2, lmin, J, x);
				if (this->isConverging)
				{
					break;
				}
				//Re-assess whether to continue backtracking:
				backtrackNeeded = fnew > fold + l*(this->alpha)*(FJdx);
			}
			else
			{
				//After the first backtrack, we use cubic interpolation:
				//g(l) = al^3 + bl^2 + g'(0)l + g(0)

				//See Numerical recipes pg 479 for explanation. Find the
				//minimum of the cubic interpolant:
				mMatrix M(2, 2);
				M << value_type(1.0) / (l1*l1), -value_type(1.0) / (l2*l2),
					-l2 / (l1*l1), l1 / (l2*l2);
				mMatrix w(2, 1);
				w << g1 - gp0*l1 - g0,
					g2 - gp0*l2 - g0;
				mMatrix V = (value_type(1.0) / (l1 - l2))*(M*w);
				value_type a = V(0);
				value_type b = V(1);
				value_type discriminant = b*b - value_type(3.0)*a*gp0;
				value_type lnew = (-b + sqrt(b*b - value_type(3.0)*a*gp0)) / (value_type(3.0)*a);
				//mexPrintf("lnew = %s\n",lnew.str(50,std::ios_base::scientific).c_str());
				//Cutoff to ensure we don't take too large or too small steps:
				if (lnew > value_type(0.5)*l2)
				{
					//l must at least halve with every step
					lnew = value_type(0.5)*l2;
				}
				else if (l2 < value_type(0.1)*l2)
				{
					//Can't decrease faster than 10% of previous step each time.
					lnew = value_type(0.1)*l2;
				}
				//New step:
				dxnew = lnew*dx;
				for (int i = 0; i < this->nComponents; i++)
				{
					std::cout << "dx[" << i << "] = " << dx(i).str(50, std::ios_base::scientific) << std::endl;
				}
				//Try F again:
				mVector xnew = x + dxnew;
				//mexPrintf("x[0] = %s\nx[1] = %s\nx[2] = %s\n",x(0).str(50,std::ios_base::scientific).c_str(),x(1).str(50,std::ios_base::scientific).c_str(),x(2).str(50,std::ios_base::scientific).c_str());
				//mexPrintf("dxnew[0] = %s\ndxnew[1] = %s\ndxnew[2] = %s\n",dxnew(0).str(50,std::ios_base::scientific).c_str(),dxnew(1).str(50,std::ios_base::scientific).c_str(),dxnew(2).str(50,std::ios_base::scientific).c_str());
				Fnew = this->F(xnew);
				fnew = value_type(0.5)*(Fnew.dot(Fnew));
				this->backtracks++;

				//Reset lambda parameters for the next run:
				l1 = l2;
				l2 = lnew;
				g1 = g2;
				g2 = fnew;
				l = l2;

				//Check for lambda being too small - possibly indicates convergence:
				this->isConverging = this->checkBacktrackConvergence(Fnew, l2, lmin, J, x);
				if (this->isConverging)
				{
					break;
				}
				//Re-assess whether to continue backtracking:
				backtrackNeeded = fnew > fold + l*(this->alpha)*(FJdx);
			}
		}
	}
	//return the backtracked step:
	return dxnew;
}
template<class value_type >
mVector NRSolver< value_type >::checkBoundedReflecting(mVector& x, mVector& dx)
{
	//Applies a battery of tests to the proposed step to ensure that it
	//does not cross a (finite) boundary, contains no NaN values (which would
	//indicate an error - this can happen sometimes) or contains infinite
	//inputs x or proposed steps dx, either of which indicates a serious error
	//with which the Newton-Step cannot continue.
	if (this->bounded)
	{
		int N = int(x.size());
		bool boundsBroken = false;
		//First ignore any bounds which are given as infinite:
		mVector xFiniteLowerBound;
		int Mlower = 0;
		mVector dxFiniteLowerBound;
		mVector xFiniteUpperBound;
		int Mupper = 0;
		mVector dxFiniteUpperBound;
		mVector lowerBoundFinite;
		mVector upperBoundFinite;
		int nOvershot = 0;
		mVector dxnew;
		mVector shortenFactorVector;

		for (int i = 0; i < N; i++)
		{
			//Filter out any non-numbers - these should generate errors.
			if (boost::math::isnan(x(i)))
			{
				this->stopReason = 7;
				return dx;
			}
			else if (boost::math::isnan(dx(i)))
			{
				this->stopReason = 8;
				return dx;
			}
			//Filter out directions with infinite bounds as we don't care about them:
			if (!boost::math::isinf(this->lowerBound[i]))
			{
				Mlower++;
				xFiniteLowerBound.conservativeResize(Mlower);
				xFiniteLowerBound(Mlower - 1) = x(i);
				dxFiniteLowerBound.conservativeResize(Mlower);
				dxFiniteLowerBound(Mlower - 1) = dx(i);
				lowerBoundFinite.conservativeResize(Mlower);
				lowerBoundFinite(Mlower - 1) = this->lowerBound[i];
			}
			if (!boost::math::isinf(this->upperBound[i]))
			{
				Mupper++;
				xFiniteUpperBound.conservativeResize(Mupper);
				xFiniteUpperBound(Mupper - 1) = x(i);
				dxFiniteUpperBound.conservativeResize(Mupper);
				dxFiniteUpperBound(Mupper - 1) = dx(i);
				upperBoundFinite.conservativeResize(Mupper);
				upperBoundFinite(Mupper - 1) = this->upperBound[i];
			}
			//Check for anything that breaks a bound.
			boundsBroken = boundsBroken || x(i) + dx(i) > this->upperBound[i] || x(i) + dx(i) < this->lowerBound[i];
			//Add an element to our shorten factor vector, to be used later:
			if (this->reflecting)
			{
				if ((x(i) + dx(i) >= this->upperBound[i] && !boost::math::isinf(this->upperBound[i])) && dx(i) != value_type(0.0))
				{
					nOvershot++;
					//shortenFactorVector.conservativeResize(nOvershot);
					//shortenFactorVector(nOvershot - 1) = (upperBound[i] - x(i))/(dx(i));
					//Shorten dx in this direction so that it does not overshoot
					//the bound:
					if (int(atUpperBoundary.size()) != this->nComponents)
					{
						//Error - atBoundary vector hasn't been set!
						throw "Attempted to check reflecting bounds before problem is set.";
					}
					else if (atUpperBoundary[i])
					{
						//Reflect, since we started at the boundary and tried
						//to overshoot it:
						dx(i) *= -value_type(1.0);
						atUpperBoundary[i] = false;
						//Make sure we don't now overshoot the other boundary!!
						if (x(i) + dx(i) <= this->lowerBound[i])
						{
							dx(i) = this->lowerBound[i] - x(i);
							atLowerBoundary[i] = true;
						}
					}
					else
					{
						//We aren't sitting at the upper boundary, but attempted
						//to overshoot it. Position us at the upper boundary!
						dx(i) = this->upperBound[i] - x(i);
						atUpperBoundary[i] = true;
					}
				}
				if ((x(i) + dx(i) <= this->lowerBound[i] && !boost::math::isinf(this->lowerBound[i])) && dx(i) != value_type(0.0))
				{
					nOvershot++;
					//shortenFactorVector.conservativeResize(nOvershot);
					//shortenFactorVector(nOvershot - 1) = (upperBound[i] - x(i))/(dx(i));
					//Shorten dx in this direction so that it does not overshoot
					//the bound:
					if (int(atLowerBoundary.size()) != this->nComponents)
					{
						//Error - atBoundary vector hasn't been set!
						throw "Attempted to check reflecting bounds before problem is set.";
					}
					else if (atLowerBoundary[i])
					{
						//Reflect, since we started at the boundary and tried
						//to overshoot it:
						dx(i) *= -value_type(1.0);
						atLowerBoundary[i] = false;
						//Make sure we don't now overshoot the other boundary!!
						if (x(i) + dx(i) >= this->upperBound[i])
						{
							dx(i) = this->upperBound[i] - x(i);
							atUpperBoundary[i] = true;
						}
					}
					else
					{
						//We aren't sitting at the lower boundary, but attempted
						//to overshoot it. Position us at the lower boundary!
						dx(i) = this->lowerBound[i] - x(i);
						atLowerBoundary[i] = true;
					}
				}
			}
			else
			{
				if ((x(i) + dx(i) > this->upperBound[i] && !boost::math::isinf(this->upperBound[i])) && dx(i) != value_type(0.0))
				{
					nOvershot++;
					shortenFactorVector.conservativeResize(nOvershot);
					shortenFactorVector(nOvershot - 1) = (this->upperBound[i] - x(i)) / (dx(i));
				}
				if ((x(i) + dx(i) < this->lowerBound[i] && !boost::math::isinf(this->lowerBound[i])) && dx(i) != value_type(0.0))
				{
					nOvershot++;
					shortenFactorVector.conservativeResize(nOvershot);
					shortenFactorVector(nOvershot - 1) = (this->upperBound[i] - x(i)) / (dx(i));
				}
			}
		}
		//Apply dx shortening if not using reflecting boundaries. Otherwise
		//we are done:
		if ((!(this->reflecting)) && (boundsBroken && nOvershot > 0))
		{
			//Find the maximal overshooting of the boundary:
			value_type shortenFactor = shortenFactorVector.cwiseAbs().maxCoeff();
			//Now shorten dx so it remains within this:
			dxnew = (value_type(0.99)*shortenFactor)*dx;
			//Check whether this step is too short:
			value_type magNew = sqrt(dxnew.dot(dxnew));
			if (magNew < this->stepMin)
			{
				std::cout << "Warning: Newton search exited due to hitting boundary.\n";
				this->stopReason = 3;
			}
			return dxnew;
		}
		else
		{
			return dx;
		}
	}
	else
	{
		return dx;
	}
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
