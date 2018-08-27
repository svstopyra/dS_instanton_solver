#ifndef FIND_ZERO
#define FIND_ZERO

template<class value_type>
bool bounded_by(const value_type x,const value_type alpha,const value_type beta)
{
    return alpha < beta ? ( (x >= alpha) && (x <= beta) )
                : ( (x >= beta) && (x <= alpha) );
}
//Swaps the values of a and b:
template<class value_type>
void swap_vars(value_type& a,value_type& b)
{
    value_type temp = a;
    a = b;
    b = temp;
}


//Code to find the zero of an arbitrary function.
//Uses the Van Wijngaarden-Dekker-Brent method.
template<class value_type, class function>
value_type find_zero_brent
(function& f,const value_type& x0,const value_type& x1,
 const value_type tol = std::numeric_limits<value_type>::epsilon())
{
    //Method is based on using three points to form a quadratic interpolation
    //of the function, and using this to find the root. We store bounds on the
    //root, and if we go outside these with the quadratic interpolant, then we
    //use a bisection step instead.

    //Three points used for interpolation:
    value_type a = x0; //First bound on root.
    value_type b = x1; //Second bound on root, and best guess at the root.
    value_type c = x1; //guess_{n - 1} (previous iteration's guess).
    value_type d;//guess_{g-2} (guess before the previous iteration's guess)

    //Needed constants:
    const value_type _0p0 = value_type(0.0);
    const value_type _2p0 = value_type(2.0);
    const value_type _1p0 = value_type(1.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _4p0 = value_type(4.0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    const int nMax = 200;//Maximum number of iterations allowed by the method

    //Values of the function at each point:
    value_type fa = f(a);
    value_type fb = f(b);
    value_type fc;
    value_type fguess;
    //std::cout << "Fine up to here.";
    if( fa*fb > _0p0 )
    {
        throw "Error. Function f must differ either side of a root.";
    }
    //Ensure that b is always the point closest to 0. Swap a and b
    //if this is not the case:
    //std::cout << "Fine up to here.";
    if(abs(fa) < abs(fb))
    {
        swap_vars(a,b);
        swap_vars(fa,fb);
    }
    //Start with c at the same value as the first bound,
    //rather than at the best guess:
    c = a; fc = fa;


    //Intermediate variables:
    value_type int_1;
    value_type int_2;
    value_type int_3;
    value_type signPQ;// = sign(P*Q)
    bool cond1;
    bool cond2;
    bool cond3;
    value_type bound1;

    //Variables for inverse quadratic interpolation:
    value_type P;
    value_type Q;
    value_type R;
    value_type S;
    value_type T;

    //Other variables:
    value_type diff;//Size of bounding interval.
    value_type delta;//intermediate tolerance variable.
    value_type guess;//New guess for the zero.

    //Tells us whether we bisected on the last iteration:
    bool bisect_flag = true;


    //Zero finding iteration:
    for(int i = 0;i < nMax;i++)
    {
        //std::cout << "\n i = " << i;
        //Ensure that

        //Check for convergence:
        delta = _2p0*eps*abs(b) + tol/_2p0;//Sum of two tolerances ensures that
            //we cater for both small b (use tol/2) and large b (use c very
            //close to b).
        /*
        std::cout << "\nb - a = " << b - a;
        std::cout << "\nb = " << b;
        std::cout << "\nf(b) = " << fb;
        std::cout << "\ndelta = " << delta;
        */
        if( abs(b - a) <= delta || fb == _0p0 )
        {
            return b;
        }

        //Inverse quadratic interpolation step:
        S = fb/fa;
        int_2 = (c - b);
        if( a == c )
        {
            //Effectively use the Secant method instead
            //Since two of our points are co-incident:
            P = -int_2*S;// = (fb/fa)*(b - c)
            Q = _1p0 - S;// = (1 - fb/fa)
            //So:
            // P/Q = - fb*(b - c)/(fb - fa);
            // thus:
            // b + P/Q = b - fb*(b - c)/(fb - fa)
            // which is the same as the secant method.
        }
        else
        {
            //Use the zero of a quadratic interpolant:
            T = fa/fc;
            R = fb/fc;
            P = S*(int_2*T*(R - T) - (b - a)*(_1p0 - R));
            Q = (T - _1p0)*(R - _1p0)*(S - _1p0);
        }


        //Decide whether to accept or reject this. We
        //want to impose bounds before dividing by Q,
        //in case Q is very small:

        //We have three conditions, of which cond2 and cond3 differ
        //depending on whether we last used a bisection or interpolation
        //step. If any of these fail, we bisect:

        //cond1 makes sure that the new point is between
        // (3a + b)/4 and b, in other words, we are within the
        //[a,b] interval and not too close to a (the worse estimate).

        //cond2 ensures that consecutive interpolation
            //steps halve every two iterations.

        //cond3 ensures that the convergence is not too slow for
        //a given numerical tolerance.


        signPQ = P*Q > _0p0 ? _1p0 : -_1p0;
        int_3 = signPQ*abs(P);
        cond1 = !(bounded_by(int_3,_0p0,_3p0*(a - b)*abs(Q)/_4p0));

        if(bisect_flag)
        {
            bound1 = abs(b - c)*abs(Q)/_2p0;
            cond2 = (int_3 >= -bound1) && (int_3 <= bound1);
            cond3 = abs(b - c) < abs(delta);
        }
        else
        {
            bound1 = abs(c - d)*abs(Q)/_2p0;
            cond2 = (int_3 >= -bound1) && (int_3 <= bound1);
            cond3 = abs(c - d) < abs(delta);
        }
        if(cond1 || cond2 || cond3)
        {
            //Use bisection.
            guess = (a + b)/_2p0;
            bisect_flag = true;
        }
        else
        {
            //Use quadratic interpolation/secant method:
            guess = b + P/Q;
            bisect_flag = false;
        }

        /*
        std::cout << "\na = " << a << " f(a) = " << fa;
        std::cout << "\nb = " << b << " f(b) = " << fb;
        std::cout << "\nc = " << c << " f(c) = " << fc;
        std::cout << "\nguess = " << guess;
        */

        //Calculate f at the new point:
        fguess = f(guess);
        //std::cout << "\nf(guess) = " << fguess;
        //Store previous guesses:
        //guess = guess_n (current guess)
        d = c;//guess_{n - 2}
        c = b;//guess_{n - 1}
        fc = fb;

        //Re-assign the bounds depending on what side of the root
        //the guess is on:
        if(fa*fguess < _0p0)
        {
            b = guess;
            fb = fguess;
        }
        else
        {
            a = guess;
            fa = fguess;
        }

        /*
        std::cout << "\na = " << a << " f(a) = " << fa;
        std::cout << "\nb = " << b << " f(b) = " << fb;
        std::cout << "\nc = " << c << " f(c) = " << fc;
        */

        //Ensure that b is still always the better guess:
        if(abs(fa) < abs(fb))
        {
            swap_vars(a,b);
            swap_vars(fa,fb);
        }
    }
    return b;
}

#endif//FIND_ZERO
