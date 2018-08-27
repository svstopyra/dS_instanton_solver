#ifndef BCADJUST_HEADER
#define BCADJUST_HEADER

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

//Header file for boundary condition adjustments. Contains no code.

//Defines functions which adjust the boundary conditions for instanton
//calculations, dealing with the fact that a(0) = 0, leading to divergences
// in the equation.

//Needed core code:
#include "potentials.h"

//Needed libraries:
#include <vector>
#include <hypergeom.h>
#include <boost/math/special_functions.hpp>


//DLL_EXPORT int fastFactorial(int n);
//==============================================================================
//Bessel functions needed for flat space boundary conditions.
//------------------------------------------------------------------------------
//Functions needed for the flat space initial conditions (bessel function
//and their derivatives)
template<class value_type>
value_type I1ox(value_type x)
{
    const value_type thresh = value_type(0.1);
    const value_type _2p0 = value_type(2.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _1p0 = value_type(1.0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    if(abs(x) < thresh)
    {
        value_type factor = x/_2p0;
        value_type factor2 = factor*factor;
        value_type sum = _0p0;
        value_type counter = _0p0;
        value_type sum_max = value_type(10000);
        value_type an = _1p0;
        do
        {
            sum += an;
            counter += _1p0;
            an *= factor2/((counter + _1p0)*(counter + _2p0));
        }
        while(abs(an/sum) > eps && counter < sum_max && abs(an) > eps);
        return sum;
    }
    else
    {
        return boost::math::cyl_bessel_i(_1p0,x)/x;
    }
}
//------------------------------------------------------------------------------
//J1(x)/x
template<class value_type>
value_type J1ox(value_type x)
{
    const value_type thresh = value_type(0.1);
    const value_type _2p0 = value_type(2.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _1p0 = value_type(1.0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    if(abs(x) < thresh)
    {
        value_type factor = x/_2p0;
        value_type factor2 = -factor*factor;
        value_type sum = _0p0;
        value_type counter = _0p0;
        value_type sum_max = value_type(10000);
        value_type an = _1p0;
        do
        {
            sum += an;
            counter += _1p0;
            an *= factor2/((counter + _1p0)*(counter + _2p0));
        }
        while(abs(an/sum) > eps && counter < sum_max && abs(an) > eps);
        return sum;
    }
    else
    {
        return boost::math::cyl_bessel_j(_1p0,x)/x;
    }
}
//------------------------------------------------------------------------------
//d(J1(x)/x)/dx
template<class value_type>
value_type dJ1oxdx(value_type x)
{
    const value_type thresh = value_type(0.1);
    const value_type _2p0 = value_type(2.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _4p0 = value_type(4.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _1p0 = value_type(1.0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    if(abs(x) < thresh)
    {
        value_type factor2 = -x*x/_4p0;
        value_type sum = _0p0;
        value_type counter = _0p0;
        value_type sum_max = value_type(10000);
        value_type an = x/_4p0;
        do
        {
            sum += an;
            counter += _1p0;
            an *= factor2/((counter + _1p0)*(counter + _3p0));
        }
        while(abs(an/sum) > eps && counter < sum_max && abs(an) > eps);
        return sum;
    }
    else
    {
        return boost::math::cyl_bessel_j(_0p0,x)/x
                - _2p0*boost::math::cyl_bessel_j(_1p0,x)/(x*x);
    }
}
//------------------------------------------------------------------------------
//d(I1(x)/x)/dx
template<class value_type>
value_type dI1oxdx(value_type x)
{
    const value_type thresh = value_type(0.1);
    const value_type _2p0 = value_type(2.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _4p0 = value_type(4.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _1p0 = value_type(1.0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    if(abs(x) < thresh)
    {
        value_type factor2 = x*x/_4p0;
        value_type sum = _0p0;
        value_type counter = _0p0;
        value_type sum_max = value_type(10000);
        value_type an = x/_4p0;
        do
        {
            sum += an;
            counter += _1p0;
            an *= factor2/((counter + _1p0)*(counter + _3p0));
        }
        while(abs(an/sum) > eps && counter < sum_max && abs(an) > eps);
        return sum;
    }
    else
    {
        return boost::math::cyl_bessel_i(_0p0,x)/x
                - _2p0*boost::math::cyl_bessel_i(_1p0,x)/(x*x);
    }
}
//------------------------------------------------------------------------------
//==============================================================================


//==============================================================================
//------------------------------------------------------------------------------
//x_tilde:
//Returns x_tilde, the radial co-ordinate in the Einstein frame, as a function
//of x, the radial co-ordinate in the Jordan frame
template< class value_type >
DLL_EXPORT value_type x_tilde
    (value_type x,value_type y0,value_type xi,value_type h,value_type W0,
     potential< value_type >& W,value_type dy0);
//------------------------------------------------------------------------------
//Gi:
template< class value_type >
DLL_EXPORT value_type Gi(value_type x,value_type alpha);
//------------------------------------------------------------------------------
//Computes beta_x(a,b), the incomplete beta function, but uses the complement
//if 0.5 < x < 1
//to ensure the arguments are always small:
template<class value_type>
DLL_EXPORT value_type beta_efficient(value_type a,value_type b,value_type x)
{
    using boost::math::beta;
    const value_type _0p0 = value_type(0.0);
    const value_type _1p0 = value_type(1.0);
    const value_type _0p5 = value_type(0.5);
    //We have to use hypergeometric functions if a OR b < 0. Otherwise,
    //use the existing functions:
    if(a > _0p0 && b > _0p0)
    {
        if(x < _0p0 || x > _1p0)
        {
            throw "Incomplete beta function argument out of range.";
        }
        else if(x < _0p5)
        {
            return beta(a,b,x);
        }
        else
        {
            value_type xp = _1p0 - x;
            return _1p0 - beta(b,a,xp);
        }
    }//Using standard incomplete beta functions.
    else
    {
        return pow(x,a)*hypergeometric::hyperGeometric2F1(a,_1p0-b,a+_1p0,x)/a;
    }//Using hypergeometric function
}
//------------------------------------------------------------------------------
template<class value_type>
DLL_EXPORT value_type beta_diff
(value_type a,value_type b,value_type x2,value_type x1)
{
    //Computes \int_{x1}^{x2}t^{a - 1}(1 - t)^{b - 1}dt = B_x2(a,b) - B_x1(a,b)
    //but attempts to bypass problems with the incomplete beta function being
    //undefined at its
    //lower limit (0) which is an infinity which should cancel.

    //Need to check whether parameters are close to negative integers, as this
    //will cause problems:
    value_type ar = round(a);
    value_type br = round(b);
    //const value_type _0p1 = value_type(0.1);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    const value_type _1p0 = value_type(1.0);
    const value_type _0p0 = value_type(0.0);
    //std::cout << "Here 2.\n";
    //std::cout << "ar = " << ar <<std::endl;
    if(abs(ar - a) < eps && ar < _0p0)
    {
        //std::cout << "Here 3.\n";
        //Close to a negative integer. Check whether the same holds for b:
        if(abs(br - b) < eps && br < _0p0)
        {
            //std::cout << "Here 5.\n";
            //Have to get creative now...
            //Following the paper by Ozcag, Ege and Gurcay
            //at
            //http://www.sciencedirect.com/science/article/pii/S0022247X07007664
            //We need:
            int n = -ar.template convert_to<int>();
            int m = -br.template convert_to<int>();

            //Formula is complicated, but tractable:
            //Ignore phi(n) = 1 + 1/2 + 1/3 + ... + 1/n
            //terms as this cancels between the two beta functions.
            value_type sum1 = _0p0;
            value_type I = _1p0;
            value_type ai = _1p0/ar;
            value_type f_x1 = pow(x1,-ar)*pow(_1p0 - x1,br - _1p0);
            value_type f_x2 = pow(x2,-ar)*pow(_1p0 - x2,br - _1p0);
            value_type x1_factor = _1p0/(_1p0 - x1);
            value_type x2_factor = _1p0/(_1p0 - x2);
            value_type x1_factorx1 = x1*x1_factor;
            value_type x2_factorx2 = x2*x2_factor;
            value_type pre_factor = _1p0 + br;//will store (m + n)!/(n!m!)
            for(int i = 1;i < n;i++)
            {
                //perform necessary sums:
                sum1 += ai*f_x2 - ai*f_x1;
                //Recursively define next parameters for sum:
                ai*= (br + I)/(ar - I);
                f_x1 *= x1_factorx1;
                f_x2 *= x2_factorx2;
                I += _1p0;
                //Update pre-factor:
                pre_factor *= (_1p0 + br/I);
            }
            //Evaluate sum2:
            value_type sum2 = _0p0;
            value_type bi_x1 = x1_factor;
            value_type bi_x2 = x2_factor;
            I = _1p0;
            for(int i = 1;i < n + m;i++)
            {
                //Perform summation of terms in sum:
                sum2 += (bi_x2 - bi_x1)/I;
                //There is another term, but it is independent of
                //x1,x2, so cancels in the difference.

                //Update recursively:
                bi_x2 *= x2_factor;
                bi_x1 *= x2_factor;
            }
            //Finally, compute the actual difference between the two beta
            //functions:
            return -sum1 + pre_factor*( log(x2_factorx2)
                   - log(x1_factorx1) + sum2 );
        }
        else
        {
        //std::cout << "Here 6.\n";
            //Exploit:
            // B_x(a,b) = B(a,b) - B_{1-x}(b,a). For a equal to a
            //negative integer, B(a,b) diverges, but this cancels between the
            //two beta functions:
            return beta_efficient(b,a,_1p0 - x1) -beta_efficient(b,a,_1p0 - x2);
        }
    }
    else
    {
        //std::cout << "Here 4.\n";
        //We should be fine to go with things as they are:
        //Try to be efficient - if a < 0 we will likely encounter problems near
        //the t = 0 part of the integral.
        if(a <= _0p0)
        {
            if(b > _0p0)
            {
                return beta_efficient(b,a,_1p0 - x1)
                       - beta_efficient(b,a,_1p0 - x2);
            }
            else
            {
                //This will be extremely difficult, as both b and a are
                //negative... Try directly summing the differences in the
                //Taylor expansion of the hypergeometric function:
                value_type result = _0p0;
                value_type n = _0p0;
                value_type x2_power = pow(x2,a);
                value_type x1_power = pow(x1,a);
                value_type alpha_n = _1p0;
                value_type a_bar_n = a;
                value_type diff;
                value_type sum_max = value_type(10000);
                do
                {
                    //Perform sum:
                    diff = alpha_n*(x2_power - x1_power)/a_bar_n;
                    result += diff;

                    //Update:
                    alpha_n *= (_1p0 - b + n)/(n + _1p0);
                    a_bar_n += _1p0;
                    n += _1p0;
                    x2_power *= x2;
                    x1_power *= x1;
                }
                while(abs(diff/result) > eps && n < sum_max);
                return result;
            }
        }
        else
        {
            return beta_efficient(a,b,x2) - beta_efficient(a,b,x1);
        }
    }
}
//------------------------------------------------------------------------------
template<class value_type>
DLL_EXPORT value_type F2
(value_type a,value_type b,value_type c,value_type eps,value_type v)
{
    //Returns a numerical approximation of the derivative of the hypergeometric
    //function 2F1(a,b,c, v^2/(1 + v^2)), with respect to the parameter c.
    //To be precise, it returns: [ 2F1(a,b,c + eps, v^2/(1 + v^2))  -
    //2F1(a,b,c, v^2/(1 + v^2)) ]/eps,
    //in  manner which is numerically stable for small epsilon (by effectively
    //constructing a power series
    //in epsilon and summing it).

    //eps = 0 (which is safe to compute with this function numerically,.
    //by design) gives
    //the derivative of the gamma function with respect to the c argument.

    //Needed constants:
    const value_type _1p0 = value_type(1.0);
    const value_type _0p0 = value_type(0.0);
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();

    //Initial conditions for recursion:
    value_type a0 = _1p0;
    value_type C0 = _0p0;

    //Intermediate variables:
    value_type v2 = v*v;
    value_type arg = v2/(_1p0 + v2);
    value_type factor[4];
    factor[0] = a;//Stores a + n
    factor[1] = b;//Stores b + n
    factor[2] = c;//Stores c + n
    factor[3] = _1p0;//Stores n + 1


    //Loop to obtain result:
    value_type diff = _0p0;
    value_type sum = a0*C0;
    value_type n = _0p0;
    value_type sum_max = value_type(10000);
    value_type Cn = C0;
    value_type an = a0;
    do
    {
        //Recursion steps:
        //Accounts for the difference in the two sequences between 1/(c + eps)_n
        // and 1/(c)_n
        //Definition: Cn = [ (c)_n/(c + eps)_n -1 ]/eps, which satisfies the
        //following recursion relation:
        Cn = (_1p0/( _1p0 + eps/factor[2] ))*( Cn - _1p0/factor[2]);
        //Accounts for the rest of the factors in the usual expansion of the
        //hypergeometric function:
        an *= factor[0]*factor[1]*arg/(factor[2]*factor[3]);

        //Add to sum:
        diff = Cn*an;
        sum += diff;

        //Advance n:
        for(int i = 0;i < 4;i++)
        {
            factor[i] += _1p0;
        }
        n += _1p0;
    }
    while(abs(diff/sum) > epsilon && n < sum_max);
    return sum;
}
//------------------------------------------------------------------------------
template<class value_type>
value_type delta_nu_half
(value_type nu,value_type eps,value_type alpha,value_type v)
{
    const value_type _1p0 = value_type(1.0);
    const value_type _2p0 = value_type(2.0);
    value_type v2 = v*v;
    value_type arg = v2/(_1p0 + v2);
    value_type index = nu - alpha;
    value_type pre_factor = pow(v2,index)/(index*sqrt(_1p0 + v2));
    const value_type a = value_type(0.5);;
    const value_type b = value_type(1.0);
    value_type c = b + index;
    value_type _2F1 = hypergeometric::hyperGeometric2F1(a,b,c + eps,arg);
    value_type Eeps = hypergeometric::Eeps(eps,log(v2));

    return pre_factor*Eeps*_2F1 + pre_factor*( -_1p0/(index + eps))*_2F1
     + pre_factor*(_1p0/(_1p0 + eps/index) - _1p0)*Eeps*_2F1
     + pre_factor*F2(a,b,c,eps,v);
}
//------------------------------------------------------------------------------
//Gi definition:
template< class value_type >
DLL_EXPORT value_type Gi(value_type x,value_type alpha)
{
    //Evaluates the integral \int_{0}^{x}\frac{C_{\alpha}^{3/2}(u)}
    //{\sqrt{u^2 - 1}}du Achieves this by expressing C_{\alpha}^{3/2}(u) as a
    //power series in u. However, different
    //power series are used depending on the value of u, to ensure the series
    //always converges. This is achieved using two different expansions of the
    //hypergeometric function 2F1, which C_{alpha}^{3/2}(u) can be
    //expressed in termns of. The 2F1 power series converges only for |z| < 1,
    //but can be transformed by various formulae so that the parameter becomes
    //1/z or z/(z-1) (there are other transformations, but we do not need them
    //as we only consider real arguments x > 1, for which
    //the 1/z and z/(z-1) series are always the fastest converging).

    //Using the two expansions, we reduce the integral to a sum of simpler
    //integrals. To do this, we have to split the integral into two parts,
    //depending on the upper limit x, so that we use the expansion (1/z or
    //z/(z-1)) giving the smallest parameter. The transition can be determined
    //to occur at  x = 2 + sqrt(5) (note, this is where 1/z = z/(z-1), where
    //z = (1-x)/2, which is the parameter for 2F1 when
    //C_{\alpha}^{3/2} is expressed in terms of the hypergeometric function).

    //Thus, the first hald, \int_{1}^{2 + sqrt(5)}(...)du is done with
    //the z/(z-1) expansion.
    //The second half (if x is such that this is needed),
    //\int_{2 + \sqrt(5)}^{x}, is done with
    //the 1/z expansion.

    //Note that the 1/z expansion is complicated by the fact that the formula
    //has a cancellation between the poles of two gamma-functions when the
    //parameters conspire to supply a negative integer. This isn't a problem
    //mathematically, but poses a numerical challenge when the parameters are
    //a 'near integer'.
    //To handle this, we parameterise the distance from an integer by a
    //variable epsilon. When epsilon is small, (here we use epsilon < 0.1),
    //we have to re-write the formula in a form that is analytic in
    //epsilon. This is significantly more complicated, and involves different
    //integrals.

    //The end result of all this is that we transform the complicated integral
    //of the Gegenbauer function into a sum which should rapidly converge to
    //the correct value of the integral.

    //Locally needed constants:
    const value_type _1p0 = value_type(1.0);
    const value_type _2p0 = value_type(2.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _5p0 = value_type(5.0);
    const value_type _0p5 = value_type(0.5);
    const value_type _0p1 = value_type(0.1);
    const value_type _4p0 = value_type(4.0);

    //Special functions required:
    using boost::math::beta;
    using hypergeometric::sign_of;
    using hypergeometric::inverse_gamma;
    using boost::math::tgamma;
    using hypergeometric::sinc;
    using boost::math::pow;
    using hypergeometric::pochhammer;
    using hypergeometric::Geps;
    using hypergeometric::Peps;
    using hypergeometric::Eeps;

    //Boundary between regions where we use the z/(z-1) and 1/z expansions
    //respectively
    // of the hypergeometric function to do the integral:
    value_type boundary = _2p0 + sqrt(_5p0);
    //Minimum threshold for convergence - onece terms added to series get
    //smaller than this, stop adding.
    value_type threshold = std::numeric_limits<value_type>::epsilon();

    if(x < _1p0)
    {
        throw "Invalid argument to Gi Integral. x must exceed 1.";
    }
    else
    {
        //Perform the first part of the integral, up to x or 2 + sqrt(5),
        //whichever is smaller:
        value_type upper1 = x < boundary ? x : boundary;

        //Just use the first integral:
        value_type v = (upper1 - _1p0)/(upper1 + _1p0);

        //Terms for recursion series:
        //hypergeometric function parameters:
        value_type a = _0p5;
        value_type b = - alpha;
        //First terms in power series:
        value_type A0 = beta_efficient(a,b,v);
        //std::cout << "A0 = " << A0 << std::endl;
        value_type Atilde0 = (alpha + _2p0)*(alpha + _1p0)/_2p0;
        value_type A = A0;
        value_type Atilde = Atilde0;
        //pre-computed powers appearing in the formula (update by recursion, as
        //comnputing powers is expensive):
        value_type _1mv_power = pow(_1p0 - v, - alpha);
        value_type _v_power = pow(v,-_0p5);
        //Loop to sum the series for this integral:
        value_type diff;//Stores each term in the series.
        value_type n = _0p0;
        value_type sum_max = value_type(10000);
        value_type sum1 = A*Atilde;//Total sum of the series so far.
        //std::cout << "sum1 = " << sum1 << std::endl;
        value_type denom;//Check for when the denominator is near to zero
        //in the recursion relation.
        do
        {
            //Advance counter:
            n += _1p0;
            //Update A and Atilde via recursion relation:
            _v_power *= v;
            denom = (n - alpha - _0p5);
            //Special case when denom is near zero (recursion relation
            //inaccurate)
            if(abs(denom) < _0p1)
            {
                A = beta_efficient(n + _0p5,-alpha,v);
            }
            else
            {
                A = ((n - _0p5)*A - _1mv_power*_v_power)/denom;
            }
            //std::cout << "n = " << n << " A = " << A;
            Atilde = (-alpha + n - _1p0)*(-alpha + n - _2p0)
                     *Atilde/(n*(n + _1p0));
            //std::cout << " Atilde = " << Atilde << std::endl;
            //Add result to the sum:
            diff = A*Atilde;
            sum1 += diff;
        }
        while(abs(diff/sum1) > threshold && n < sum_max);
        //std::cout << "sum1 = " << sum1 << std::endl;

        //By this point, we have computed the first half of the integral.
        //We only proceed if
        //x > boundary, ie, we need to 1/z portion of the hypergeometric
        //function:
        //std::cout << "x = " << x << " boundary = " << boundary
                    //<< " x - boundary = " << x - boundary << std::endl;
        if(x > boundary)
        {
            //We now have to compute the second half of the integral
            //\int_{boundary}^{x} \frac{C_{\alpha}^{3/2}(u)}{\sqrt{Pu^2 - 1}}du
            //Using the 1/z expansion instead (which is MUCH more complicated
            //due to the need
            //to carefully perform a cancellation between two infinities
            //arising near certain values of the parameters - the naiive
            //approach would yield catastrophic rounding
            //errors otherwise).

            //Hypergeometric paramters are a = -alpha, b = alpha + 3, c = 2

            //Catastrophic rounding errors occur near M = integer, due to
            //cancellation of two infinities. Have to treat the system carefully
            //when eps < 0.1:
            value_type M = round(_2p0*alpha + _3p0);//Nearest integer to
                //parameter.
            int m = M.template convert_to<int>();
            value_type eps = _2p0*alpha + _3p0 - M;//parameter characterising
                //distance from
                //the nearest integer, which will cause problems with the gamma
                //functions. Variables to store two sums (sum21 a sum of a
                //finite number of terms,
            //sum22 a sum of an infinite number of terms which we sum until
            //convergence)
            value_type sum21;
            value_type sum22;

            //Change variables, so find the new boundary conditions for the
            //integral:
            value_type upper_re_arranged = (x - _1p0)/(x + _1p0);
            value_type boundary_re_arranged = (boundary - _1p0)/(boundary
                       + _1p0);

            //A few intermediates we will need:

            value_type gamma_c = _1p0; //Gamma(2)
            value_type inverse_gamma_malpha_pM_peps = inverse_gamma(-alpha
                                                                    + M + eps);
            value_type inverse_gamma_c_p_alpha = inverse_gamma(_2p0 + alpha);
            value_type g_0 = Atilde0*( abs(eps) > 0 ? gamma_c*
                                      inverse_gamma_malpha_pM_peps*
                                      inverse_gamma_c_p_alpha*
                                      inverse_gamma(_1p0 - M - eps)/(eps)
                                      : sign_of<value_type>(m)*tgamma(M)*
                                      (gamma_c*(inverse_gamma(-alpha + M)*
                                      inverse_gamma(alpha + _2p0))));

            //Over-all factor multiplying the sums:
            value_type overall_factor = sign_of<value_type>(m)/sinc(eps);

            //sum21, which is finite:
            //First terms:
            value_type beta_0 = beta_diff(alpha + _0p5,-alpha,upper_re_arranged,
                                          boundary_re_arranged);
            //std::cout << "\nbeta_0 = " << beta_0;
            value_type beta_n = beta_0;
            value_type g_n = g_0;
            value_type sign = _1p0;
            sum21 = g_n*beta_0;
            //Intermediates, update by recursion as theses are
            //expensive to do at every step:
            value_type vu = sqrt(_2p0/(x - _1p0));
            value_type vl = sqrt(_2p0/(boundary - _1p0));
            value_type sqrt_1_p_vu2 = sqrt(_1p0  + vu*vu);
            value_type sqrt_1_p_vl2 = sqrt(_1p0  + vl*vl);
            value_type vu_power = pow(vu,-_2p0*alpha - _2p0);
            value_type vl_power = pow(vl,-_2p0*alpha - _2p0);
            //std::cout << "sum21 = " << sum21 << std::endl;
            //std::cout << "beta_n = " << beta_n << std::endl;
            //std::cout << "g_n = " << g_n << std::endl;
            //Sum the finite series:
            n = _0p0;
            value_type denom_beta;
            //std::cout << "m = " << m << std::endl;
            while(n < M - _1p0)
            {
                //Recursion steps:
                sign *= -_1p0;
                g_n = (-alpha + n)*(-alpha + n - _1p0)*g_n/((n + _1p0)*
                      (n + _1p0 - M - eps));
                //std::cout << "g_n = " << g_n << " g_n error = " << g_n -
                //std::cout << "n = " << n << " g_n = " << g_n << std::endl;
                vu_power *= vu*vu;
                vl_power *= vl*vl;
                //If 1 - 2(n - alpha) is close to 0, then we need to use a
                //different method to evaluate it, as the recursion relation
                //fails: Fortunately, this can only happen for a single value
                // of n:
                n += _1p0;
                denom_beta = _1p0 - _2p0*(n - alpha);
                if( abs(denom_beta) < _0p1 )
                {
                    beta_n = beta_diff(alpha -n + _0p5,n-alpha,
                                       upper_re_arranged,boundary_re_arranged);
                }
                else
                {
                    beta_n = (_2p0*(vu_power*sqrt_1_p_vu2 -
                                    vl_power*sqrt_1_p_vl2) + (_2p0*(n - alpha)
                                    - _2p0)*beta_n)/(_1p0 - _2p0*(n - alpha));
                }
                //Add to sum:
                sum21 += g_n*sign*beta_n;
            }
            //std::cout << "n = " << n << " sum21 = " << sum21 << std::endl;

            //sum21 brought us up to beta_{m-1}. We first need to update this
            //into beta_m,
            //which appears in the first term of sum22:
            vu_power *= vu*vu;
            vl_power *= vl*vl;
            //std::cout << "denom = " << _1p0 - _2p0*(M - alpha) << std::endl;
            value_type beta_nu = abs(_1p0 - _2p0*(M - alpha)) < _0p1 ?
            (_2p0*(vu_power*sqrt_1_p_vu2 - vl_power*sqrt_1_p_vl2) +
            (_2p0*(M - alpha) - _2p0)*beta_n)/(_1p0 - _2p0*(M - alpha))
            : beta_diff(alpha - M + _0p5,M - alpha,upper_re_arranged,
                        boundary_re_arranged);
            //std::cout << "Here 1.\n";
            value_type beta_nu_p_eps = beta_diff(alpha - M - eps
                                                 + _0p5,M + eps - alpha,
                                                 upper_re_arranged,
                                                 boundary_re_arranged);
            sign = - sign;//(-1)^m
            value_type v2_upper = vu*vu;
            value_type v2_lower = vl*vl;
            value_type veps_power_u = pow(v2_upper,eps);
            value_type veps_power_l = pow(v2_lower,eps);

            //sum22 is vulnerable to round off errors when eps < 0.1.
            //Use a different
            //expansion when this is the case:
            if(abs(eps) < _0p1)
            {
                //Use the special method as we are near a pole of the gamma
                //function and need to express the result as an analytic
                //function of eps, while cancelling
                //the divergent gamma functions evaluated on negative integers.


                //Gamma function and pochhammer symbol intermediates:
                value_type inverse_gamma_alpha_p_2 = inverse_gamma(alpha+ _2p0);
                value_type inverse_gamma_1_m_eps = inverse_gamma(_1p0 - eps);
                value_type inverse_gamma_mapha_p_M_p_eps =
                    inverse_gamma(-alpha + M + eps);
                value_type inverse_gamma_M_p_1 = inverse_gamma(M + _1p0);
                value_type inverse_gamma_malpha_pM = inverse_gamma(-alpha + M);
                value_type inverse_gamma_M_p1_peps =
                    inverse_gamma(M + _1p0 + eps);
                value_type pochhammer_m1_malpha_peps =
                    pochhammer(-_1p0 - alpha + eps,m);
                value_type pochhammer_malpha = pochhammer(-alpha,m);

                //First terms in series:
                value_type b_tilde_0_1 = Atilde0*gamma_c*pochhammer_malpha
                                        *( (pochhammer_m1_malpha_peps
                                        *Geps(-eps,_1p0)
                                        - Peps(eps,-_1p0 - alpha,m)
                                        *inverse_gamma_1_m_eps )
                                        *(inverse_gamma_alpha_p_2
                                        *inverse_gamma_mapha_p_M_p_eps
                                        *inverse_gamma_M_p_1)
                                        + pochhammer_m1_malpha_peps
                                        *(Geps(eps,M + _1p0)
                                        *inverse_gamma_alpha_p_2
                                        *inverse_gamma_mapha_p_M_p_eps
                                        - Geps(eps,-alpha + M)
                                        *inverse_gamma_alpha_p_2
                                        *inverse_gamma_M_p1_peps
                                        - Geps(-eps,alpha + _2p0)
                                        *inverse_gamma_malpha_pM
                                        *inverse_gamma_M_p1_peps ) );
                value_type b_tilde_0_2 = -Atilde0*pochhammer_malpha*gamma_c
                                         *pochhammer_m1_malpha_peps
                                         *inverse_gamma_M_p1_peps
                                         *inverse_gamma_malpha_pM
                                         *inverse_gamma(alpha + _2p0 - eps);
                value_type g_tilde_0 = Atilde0*pochhammer_malpha
                                        *pochhammer( -_1p0 - alpha,m)*gamma_c
                                        *inverse_gamma_mapha_p_M_p_eps
                                        *inverse_gamma_alpha_p_2
                                        *inverse_gamma_1_m_eps
                                        *inverse_gamma_M_p_1;
                value_type b_tilde_n_1 = b_tilde_0_1;
                value_type b_tilde_n_2 = b_tilde_0_2;
                value_type g_tilde_n = g_tilde_0;

                //Determine initial condition for delta:
                value_type nu = M;
                value_type delta_0 = delta_nu_half(nu,eps,alpha,vl)
                                    - delta_nu_half(nu,eps,alpha,vu);
                value_type delta_nu = delta_0;



                //Evaluate sum22:
                n = _0p0;
                diff = _0p0;
                //Intermediate variables, to optimise number of operations:
                value_type fraction3;
                value_type factor1;
                value_type factor2;
                value_type factor3;
                value_type denom1;

                //Factors needed for b_tilde_n updates. This recusion
                // relation is of the form
                //b_{n+1} = f(n)b_{n} so we start from  n = 0;
                value_type _malpha_pM_pn = -alpha + M + n;
                value_type _malpha_pM_pn_peps = _malpha_pM_pn + eps;
                value_type _m1_malpha_pM_pn = -_1p0 - alpha + M + n;
                value_type _m1_malpha_pM_pn_peps = _m1_malpha_pM_pn + eps;
                value_type _M_pn_p1 = M + n + _1p0;
                value_type _n_p1 = n + _1p0;
                value_type _M_pn_p1_peps = _M_pn_p1 + eps;
                value_type _n_p1_meps = _n_p1 - eps;

                value_type _1_m2numalpha = _1p0 - _2p0*(nu - alpha);
                value_type Eeps_u = hypergeometric::Eeps(eps,log(v2_upper));
                value_type Eeps_l = hypergeometric::Eeps(eps,log(v2_lower));
                sum22 = sign*(b_tilde_0_1*beta_nu + b_tilde_0_2*delta_0);
                    //n = 0 term
                //std::cout << "sum22 = " << sum22 << std::endl;
                //std::cout << "sum21 = " << sum21 << std::endl;
                //std::cout << "sum1 = " << sum1 << std::endl;
                /*
                value_type beta_n_2;
                value_type beta_n_1;
                value_type a = -alpha;
                value_type b = alpha + _3p0;
                value_type c = _2p0;
                */
                value_type a_n;
                value_type b_n;
                value_type c_n;
                do
                {
                    a_n = _malpha_pM_pn_peps*_m1_malpha_pM_pn_peps
                          /(_M_pn_p1_peps*_n_p1);
                    b_n = ( (_malpha_pM_pn*_m1_malpha_pM_pn)
                           /(_M_pn_p1) -_malpha_pM_pn
                           -  _m1_malpha_pM_pn_peps
                           + (_malpha_pM_pn_peps*_m1_malpha_pM_pn_peps)
                           /(_n_p1) )/(_M_pn_p1_peps*_n_p1_meps);
                    c_n = _malpha_pM_pn*_m1_malpha_pM_pn/(_M_pn_p1*_n_p1_meps);

                    b_tilde_n_1 = a_n*b_tilde_n_1 + b_n*g_tilde_n;
                    g_tilde_n *= c_n;
                    b_tilde_n_2 *= a_n;

                    //delta updates:
                    vu_power *= vu*vu;
                    vl_power *= vl*vl;
                    nu = M + n + _1p0;
                    _1_m2numalpha = _1p0 - _2p0*(nu - alpha);
                    denom_beta = _1_m2numalpha;
                    if( abs(denom_beta) < _0p1
                       || abs(denom_beta - _2p0*eps) < _0p1 )
                    //if(true)
                    {
                        //Have to do this via the explicit formula instead:
                        delta_nu = delta_nu_half(nu,eps,alpha,vl)
                                    - delta_nu_half(nu,eps,alpha,vu);

                        //beta updates:
                        beta_nu = beta_diff(alpha - nu + _0p5,nu -
                                            alpha,upper_re_arranged,
                                            boundary_re_arranged);
                        beta_nu_p_eps = beta_diff(alpha - nu - eps + _0p5,
                                                  nu + eps - alpha,
                                                  upper_re_arranged,
                                                  boundary_re_arranged);
                    }
                    else
                    {
                        //Update via the recursion formula:
                        delta_nu = ( ( _4p0/(_1_m2numalpha - _2p0*eps) )
                                    *( vu_power*veps_power_u*sqrt_1_p_vu2
                                      - vl_power*veps_power_l*sqrt_1_p_vl2
                                      + (nu - alpha - _1p0)*beta_nu_p_eps
                                      + _2p0*eps*beta_nu_p_eps )
                                      + _2p0*( Eeps_u*vu_power*sqrt_1_p_vu2
                                      - Eeps_l*vl_power*sqrt_1_p_vl2
                                      + (nu - alpha - _1p0)
                                      *delta_nu + beta_nu_p_eps ) )
                                      /_1_m2numalpha;

                        //Beta updates
                        beta_nu = (_2p0*(vu_power*sqrt_1_p_vu2 - vl_power
                                  *sqrt_1_p_vl2) + (_2p0*(nu - alpha) - _2p0)
                                   *beta_nu)/(_1p0 - _2p0*(nu - alpha));
                        beta_nu_p_eps = (_2p0*(vu_power*veps_power_u
                                        *sqrt_1_p_vu2 - vl_power*veps_power_l
                                        *sqrt_1_p_vl2)
                                        + (_2p0*(nu + eps - alpha) - _2p0)
                                         *beta_nu_p_eps)
                                         /(_1p0 - _2p0*(nu + eps - alpha));
                    }

                    //Update sum:
                    sign *= -_1p0;
                    diff = sign*(b_tilde_n_1*beta_nu + b_tilde_n_2*delta_nu);
                    sum22 += diff;

                    //Update all factors for the next cycle (this avoids us
                    //having to
                    //recompute the factors each time, which can use more
                    //operations than
                    //simply adding 1 to them all).
                    _malpha_pM_pn += _1p0;
                    _malpha_pM_pn_peps += _1p0;
                    _m1_malpha_pM_pn += _1p0;
                    _m1_malpha_pM_pn_peps += _1p0;
                    _M_pn_p1 += _1p0;
                    _n_p1 += _1p0;
                    _M_pn_p1_peps += _1p0;
                    _n_p1_meps += _1p0;

                    n += _1p0;
                }
                while((abs(diff/sum22) > threshold )
                      && abs(diff) > threshold && n < sum_max);
                //std::cout << "n = " << n << " sum22 = " << sum22 << std::endl;
                //std::cout << "sum21 = " << sum21 << std::endl;
                //std::cout << "sum1 = " << sum1 << std::endl;
                return sum1 + overall_factor*(sum21 + sum22);
            }//if(epsilon < 0.1)
            else
            {
                //Use the normal method, as cancellation of the infinities is
                //not a
                //problem here (we aren't near the pole of a gamma function).

                //Intermediate gammas:
                value_type gamma_c_p_alpha = tgamma(_2p0 + alpha);
                value_type gamma_malpha_pM_peps = tgamma(-alpha + M + eps);



                //sum22, which is an infinite sum:
                //value_type strange_factor = pochhammer(alpha + _3p0 + eps,m);
                value_type b_0_1 = Atilde0*gamma_c*(pochhammer(-alpha,m)
                                   *pochhammer(-_1p0 - alpha,m)
                                   *inverse_gamma(_1p0 - eps)
                                   *inverse_gamma(M + _1p0)
                                   /(gamma_c_p_alpha*gamma_malpha_pM_peps))/eps;
                value_type b_0_2 = -Atilde0*gamma_c
                                   *(pochhammer(-_1p0 - alpha + eps,m)
                                   *inverse_gamma(-alpha)
                                   *inverse_gamma(_2p0 + alpha - eps)
                                   *inverse_gamma(M + _1p0 + eps))/eps;
                int N = 0;


                //To perform the sum, we need to track beta_{n + m} and
                //beta_{n + m + eps} separately:
                value_type nu = M + eps;
                //Expensive incomplete beta function evaluation:


                //Evaluate sum22:

                //Factors pre-generated to reduce number of operations inside
                //loop (NB offset from n -> n+ 1
                //due to formulation of the recursion relation):
                n = _0p0;
                diff = _0p0;
                //sign *= -_1p0;
                sum22 = sign*b_0_1*beta_nu + sign*b_0_2*beta_nu_p_eps;//Firat
                    //term in sum
                value_type b_n_1 = b_0_1;
                value_type b_n_2 = b_0_2;
                value_type _malpha_pn_pM_m1 = (-alpha + (n + _1p0) + M - _1p0);
                value_type _malpha_pn_pM_m2 = (-alpha + (n + _1p0) + M - _2p0);
                value_type _n_pM = ((n + _1p0) + M);
                value_type _n_meps = ((n + _1p0) - eps);
                value_type _malpha_peps_pn_pM_m1 = (-alpha + eps + (n + _1p0)
                                                    + M - _1p0);
                value_type _m2_malpha_peps_pn_pM = (-_2p0 - alpha + eps +
                                                    (n + _1p0) + M);
                value_type _n_pM_peps = ((n + _1p0) + M + eps);
                //std::cout << "sum22 = " << sum22 << std::endl;
                value_type _1_m2numalpha;

                //Summation loop:
                //std::cout << "m = " << m <<std::endl;
                do
                {
                    //Update recursion terms:
                    sign *= -_1p0;
                    b_n_1 = _malpha_pn_pM_m1*_malpha_pn_pM_m2*b_n_1
                             /(_n_pM*_n_meps);
                    b_n_2 = _malpha_peps_pn_pM_m1*_m2_malpha_peps_pn_pM*b_n_2
                             /((n + _1p0)*_n_pM_peps);

                    vu_power *= vu*vu;
                    vl_power *= vl*vl;
                    nu = M + n + _1p0;//+1 to account for b_n = f(n)b_{n-1}
                        //recursion relation
                    _1_m2numalpha = _1p0 - _2p0*(nu - alpha);
                    denom_beta = _1_m2numalpha;
                    if( abs(denom_beta) < _0p1
                       || abs(denom_beta - _2p0*eps) < _0p1 )
                    {
                        //beta updates:
                        beta_nu = beta_diff(alpha - nu + _0p5,nu - alpha,
                                            upper_re_arranged,
                                            boundary_re_arranged);
                        beta_nu_p_eps = beta_diff(alpha - nu - eps + _0p5,
                                                  nu + eps - alpha,
                                                  upper_re_arranged,
                                                  boundary_re_arranged);
                    }
                    else
                    {
                        //Beta updates
                        beta_nu = (_2p0*(vu_power*sqrt_1_p_vu2 - vl_power
                                  *sqrt_1_p_vl2) + (_2p0*(nu - alpha)
                                  - _2p0)*beta_nu)/(_1p0 - _2p0*(nu - alpha));
                        beta_nu_p_eps = (_2p0*(vu_power*veps_power_u
                                        *sqrt_1_p_vu2 - vl_power*veps_power_l
                                        *sqrt_1_p_vl2)
                                        + (_2p0*(nu + eps - alpha) - _2p0)
                                        *beta_nu_p_eps)/(_1p0 - _2p0
                                        *(nu + eps - alpha));
                    }
                    //Add to the sum:
                    diff = sign*(b_n_1*beta_nu + b_n_2*beta_nu_p_eps);
                    sum22 += diff;

                    //Update factors:
                    n += _1p0;
                    _malpha_pn_pM_m1 += _1p0;
                    _malpha_pn_pM_m2 += _1p0;
                    _n_pM += _1p0;
                    _n_meps += _1p0;
                    _malpha_peps_pn_pM_m1 += _1p0;
                    _m2_malpha_peps_pn_pM += _1p0;
                    _n_pM_peps += _1p0;
                }
                while((abs(diff/sum22) > threshold && abs(diff) > threshold )
                      && n < sum_max);
            }//if(epsilon > 0.1)
            return sum1 + overall_factor*(sum21 + sum22);
        }//if(x > boundary)
        else //x < 2 + sqrt(5), so we don't need the second integral using
            //the 1/z method
        {
            //What we have already computed is the entirety of the integral:
            return sum1;
        }//if(x < boundary)
    }//if(x > 1)
}
//==============================================================================

//==============================================================================
//------------------------------------------------------------------------------
//x_tilde definition:
//Returns x_tilde, the radial co-ordinate in the Einstein frame, as a function
//of x, the radial co-ordinate in the Jordan frame
template< class value_type >
DLL_EXPORT value_type x_tilde(value_type x,value_type y0,value_type xi,
                              value_type h,value_type W0,
                              potential< value_type >& W,value_type dy0)
{
    const value_type _1p0 = value_type(1.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _2p0 = value_type(2.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _6p0 = value_type(6.0);
    const value_type _1p5 = value_type(1.5);
    const value_type _4p0 = value_type(4.0);
    const value_type _8p0 = value_type(8.0);
    const value_type _24p0 = value_type(24.0);
    value_type factor = sqrt(_1p0 - xi*h*h*y0*y0);
    value_type conformal_factor = sqrt(_1p0 - xi*(_1p0 - _6p0*xi)*h*h*y0*y0);
    value_type factor2 = factor*factor;
    value_type Wy0 = W(y0);
    value_type Wy0pW0 = Wy0  + W0;
    value_type Wpy0 = W.d(y0);
    value_type Wppy0 = W.d2(y0);
    value_type alpha = -_1p5 + sqrt(_1p5*_1p5 + Wppy0);
    bool hyperbolic = Wy0pW0 < _0p0;
    value_type H_tilde;
    if(hyperbolic)
    {
        H_tilde = sqrt(-h*h*Wy0pW0/(_3p0*factor*factor));//H_tilde in units of H
    }
    else
    {
        //Positive starting position, which rather confuses everything, but we
        //simply convert to
        //trigonometric functions to get the right equations:
        H_tilde = sqrt(h*h*Wy0pW0/(_3p0*factor*factor));
        //Flipping the sign of Wy0pW0 has the effect of multiplying H_tilde by
        //i, so we exchange:
        //cosh(H_tilde*x) <-> cos(H_tilde*x)
        //sinh(H_tilde*x)/H_tilde <->  sin(H_tilde*x)/H_tilde
    }


    value_type A_num = Wpy0 + _4p0*xi*h*h*y0*(Wy0pW0)/factor2;
    value_type A_den1 = A_num*factor*(xi*xi*(_1p0 - _6p0*xi)*h*h*h*h*y0*y0*y0
                        - xi*h*h*y0*(_1p0  + _1p0*xi))/(conformal_factor
                        *conformal_factor*conformal_factor);
    value_type A_den2 = (Wppy0 + _8p0*xi*h*h*y0*Wpy0/factor2 + _4p0*xi*h*h
                         *(Wy0pW0)/factor2 + _24p0*xi*xi*h*h*h*h*y0*y0*Wy0pW0
                         /(factor2*factor2))*factor2/conformal_factor;
    value_type A = -(A_num)/(A_den1 + A_den2);
    value_type B = _2p0*(dy0 - A)/((alpha + _1p0)*(alpha + _2p0));

    value_type cosh_factor = hyperbolic ? cosh(H_tilde*factor*x)
                            : cos(H_tilde*factor*x);
    std::cout << "\nHtilde = " << H_tilde;
    std::cout << "\nfactor = " << factor;
    std::cout << "\nx = " << x;
    std::cout << "\ncosh_factor = " << cosh_factor;


    return factor*x - ((xi*h*h*y0)/(conformal_factor))*(A*factor*x
            + (B/H_tilde)*(Gi(cosh_factor,alpha)));
}
//------------------------------------------------------------------------------
//Struct to store data and compute intermediate variables, when calling the
//constructor pertaining to the computation of the boundary conditions by the
//linearisation method. Useful for passing all these arguments between
//different functions without having to give them all as separate arguments.
template<class value_type>
struct x_tilde_data
{
    value_type A;
    value_type B;
    value_type alpha;
    value_type& y0;
    bool hyperbolic;
    value_type& h;
    value_type& xi;
    value_type H_tilde;
    value_type factor;
    value_type factor2;
    value_type conformal_factor;
    //Potential and values calculated for y0:
    //Jordan Frame:
    value_type Wy0_J;
    value_type Wy0pW0_J;
    value_type Wpy0_J;
    value_type Wppy0_J;
    potential<value_type>& W;
    value_type factor4;
    value_type conformal_factor2;
    value_type conformal_factor4;
    value_type h2;
    value_type h4;
    value_type W0;
    value_type Wy0;
    value_type Wy0pW0;
    value_type Wpy0;
    value_type Wppy0;
    value_type dy0;
    value_type x_tilde;//chosen step in the Einstein frame

    //Useful constants:
    static const value_type _1p0;// = value_type(1.0);
    static const value_type _0p0;// = value_type(0.0);
    static const value_type _2p0;// = value_type(2.0);
    static const value_type _3p0;// = value_type(3.0);
    static const value_type _5p0;// = value_type(3.0);
    static const value_type _6p0;// = value_type(6.0);
    static const value_type _1p5;// = value_type(1.5);
    static const value_type _4p0;// = value_type(4.0);
    static const value_type _8p0;// = value_type(8.0);
    static const value_type _24p0;// = value_type(24.0);
    //Constructor:
    x_tilde_data(value_type& Y0,value_type& H,value_type& XI,
                 potential<value_type>& w,value_type& w0)
                    : y0(Y0),h(H),xi(XI),W(w),W0(w0)
    {
        //Setup the non-minimal coupling factors:
        factor = sqrt(_1p0 - xi*h*h*y0*y0);
        factor2 = factor*factor;
        conformal_factor = sqrt(_1p0 - xi*(_1p0 - _6p0*xi)*h*h*y0*y0);
        //Setup potential:
        //Jordan Frame:
        Wy0_J = W(y0);
        Wy0pW0_J = Wy0_J  + W0;
        Wpy0_J = W.d(y0);
        Wppy0_J = W.d2(y0);
        std::cout << "\nWppy0_J = " << Wppy0_J;
        //Various factors for efficiency:
        factor4 = factor2*factor2;
        conformal_factor2 = conformal_factor*conformal_factor;
        conformal_factor4 = conformal_factor2*conformal_factor2;
        h2 = h*h;
        h4 = h2*h2;
        //Einstein frame:
        Wy0 = Wy0_J/factor4;
        Wy0pW0 = Wy0pW0_J/factor4;
        Wpy0 = ( Wpy0_J/factor2  + _4p0*xi*h2*Wy0pW0/factor4 )/conformal_factor;
        Wppy0 = ( Wpy0_J/factor2 + _4p0*xi*h2*Wy0pW0/factor4 )
                *(xi*xi*(_1p0 - _6p0*xi)*h4*y0*y0*y0 - xi*h2*y0
                *(_1p0 + _6p0*xi))/conformal_factor4
                + ( Wppy0_J + _8p0*xi*h2*y0*Wpy0_J/factor2
                + _4p0*xi*h2*Wy0pW0/factor2 + _24p0*xi*xi*h4*y0*y0*Wy0pW0
                /factor4 )/conformal_factor2;
        std::cout << "\nWy0pW0_J = " << Wy0pW0_J;
        std::cout << "\nWppy0 = " << Wppy0;
        std::cout << "\nWpy0_J = " << Wpy0_J;
        std::cout << "\nfactor2 = " << factor2;
        std::cout << "\nWy0pW0 = " << Wy0pW0;
        std::cout << "\nfactor4 = " << factor4;
        std::cout << "\nconformal_factor = " << conformal_factor;
        std::cout << "\nconformal_factor4 = " << conformal_factor4;
        std::cout << "\nWy0pW0_J = " << Wy0pW0_J;
        //Compute A and B:
        alpha = -_1p5 + sqrt(_1p5*_1p5 + Wppy0);
        dy0 = _0p0;
        value_type A_num = Wpy0_J + _4p0*xi*h2*y0*(Wy0pW0_J)/factor2;
        value_type A_den1 = A_num*factor*(xi*xi*(_1p0 - _6p0*xi)*h4*y0*y0*y0
                            - xi*h2*y0*(_1p0  + _6p0*xi))/(conformal_factor
                            *conformal_factor*conformal_factor);
        value_type A_den2 = (Wppy0_J + _8p0*xi*h2*y0*Wpy0_J/factor2
                            + _4p0*xi*h2*(Wy0pW0_J)/factor2
                            + _24p0*xi*xi*h4*y0*y0*Wy0pW0_J/(factor2*factor2))
                            *factor2/conformal_factor;
        A = -(A_num)/(A_den1 + A_den2);
        B = _2p0*(dy0 - A)/((alpha + _1p0)*(alpha + _2p0));
        //Condition for hyperbolic or trigonometric scale factor function:
        hyperbolic = Wy0pW0_J < _0p0;
        if(hyperbolic)
        {
            H_tilde = sqrt(-h2*Wy0pW0_J/(_3p0*factor2));//H_tilde in units of H
        }
        else
        {
            //Positive starting position, which rather confuses everything,
            //but we simply convert to
            //trigonometric functions to get the right equations:
            H_tilde = sqrt(h2*Wy0pW0_J/(_3p0*factor2));
            //Flipping the sign of Wy0pW0 has the effect of multiplying
            //H_tilde by i, so we exchange:
            //cosh(H_tilde*x) <-> cos(H_tilde*x)
            //sinh(H_tilde*x)/H_tilde <->  sin(H_tilde*x)/H_tilde
        }
        x_tilde = _0p0;
    }
};
//------------------------------------------------------------------------------
template<class value_type> const value_type x_tilde_data<value_type>::_1p0
    = value_type(1.0);
template<class value_type> const value_type x_tilde_data<value_type>::_0p0
    = value_type(0.0);
template<class value_type> const value_type x_tilde_data<value_type>::_2p0
    = value_type(2.0);
template<class value_type> const value_type x_tilde_data<value_type>::_3p0
    = value_type(3.0);
template<class value_type> const value_type x_tilde_data<value_type>::_5p0
    = value_type(5.0);
template<class value_type> const value_type x_tilde_data<value_type>::_6p0
    = value_type(6.0);
template<class value_type> const value_type x_tilde_data<value_type>::_1p5
    = value_type(1.5);
template<class value_type> const value_type x_tilde_data<value_type>::_4p0
    = value_type(4.0);
template<class value_type> const value_type x_tilde_data<value_type>::_8p0
        = value_type(8.0);
template<class value_type> const value_type x_tilde_data<value_type>::_24p0
    = value_type(24.0);
//------------------------------------------------------------------------------
//Simple Newton-Raphson method:
template<class value_type,class function>
value_type NewtonRaphsonSimple(value_type guess,function f,function fp)
{
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    int counter = 0;
    value_type x = guess;
    value_type fx = f(x);
    value_type fpx = fp(x);
    while(abs(fx) > eps && counter < 10000)
    {
        //Solve using the Newton-Raphson solver until within the threshold:
        ++counter;
        x = x - fx/fpx;
        //Re-compute functions:
        fx = f(x);
        fpx = fp(x);
    }
    if(counter > 10000)
    {
        throw "Newton-Raphson loop took too long to converge.";
    }
    return x;
}
//------------------------------------------------------------------------------
//Class to store simple functions used by the linearisation of the bounce
//equations in the Einstein Frame. Used by several
//slightly different boundary conditions functions.
template<class value_type>
class instanton_equation_lineariser : public x_tilde_data<value_type>
{
public:
    //Constructor: automatically calls the x_tilde_data constructor which
    //initialises and constructs all the intermediate variables we will need.
    instanton_equation_lineariser(value_type& Y0,value_type& H,value_type& XI,
                 potential<value_type>& w,value_type& w0)
                    : x_tilde_data<value_type>(Y0,H,XI,w,w0)
    {

    }
    //Member functions:
    //Linearised deviation from the initial value (Einstein frame):
    value_type delta_y(value_type& x)
    {
        value_type lambda = this->_3p0/this->_2p0;
        std::cout << "\nlambda = " << lambda;
        std::cout << "\nH_tilde = " << this->H_tilde;
        std::cout << "\nx = " << x;
        value_type cosh_factor = this->hyperbolic ? cosh(this->H_tilde*x)
                                : cos(this->H_tilde*x);
        std::cout << "\ncosh_factor = " << cosh_factor;
        value_type vt_return = this->A + this->B*hypergeometric::GegenbauerC
            (this->alpha,lambda,cosh_factor);
        std::cout << "\nvt_return = " << vt_return;
        return vt_return;
    }
    //Linearised derivative of the deviation from the initial value
    //(Einstein frame):
    value_type delta_yp(value_type& x)
    {
        value_type lambda = this->_3p0/this->_2p0;
        value_type cosh_factor = this->hyperbolic ? cosh(this->H_tilde*x)
                                : cos(this->H_tilde*x);
        value_type sinh_factor = this->hyperbolic ? this->H_tilde
                                 *sinh(this->H_tilde*x) :
                                 - this->H_tilde*sin(this->H_tilde*x);
        return this->B*this->_2p0*lambda*hypergeometric::GegenbauerC
                                (this->alpha - this->_1p0,lambda + this->_1p0,
                                 cosh_factor)*sinh_factor;
    }
    //Flat space deltay:
    value_type delta_y_flat(value_type& x)
    {
        //Different function, depending on Vpp0:
        if(this->Wppy0_J < this->_0p0)
        {
            return (this->Wpy0_J/this->Wppy0_J)
                    *(this->_2p0*J1ox(sqrt(abs(this->Wppy0_J))*x) - this->_1p0);
        }
        else if(this->Wppy0_J > this->_0p0)
        {
            return (this->Wpy0_J/this->Wppy0_J)
                    *(this->_2p0*I1ox(sqrt(this->Wppy0_J)*x) - this->_1p0);
        }
        else
        {
            return this->Wpy0_J*(x*x)/this->_8p0;
        }
    }
    //Flat space deltayp:
    value_type delta_yp_flat(value_type& x)
    {
        //Different function, depending on Vpp0:
        if(this->Wppy0_J < this->_0p0)
        {
            return (this->Wpy0_J/this->Wppy0_J)
                    *(this->_2p0*dJ1oxdx(sqrt(abs(this->Wppy0_J))*x));
        }
        else if(this->Wppy0_J > this->_0p0)
        {
            return (this->Wpy0_J/this->Wppy0_J)
                    *(this->_2p0*dI1oxdx(sqrt(this->Wppy0_J)*x));
        }
        else
        {
            return this->Wpy0_J*x/this->_4p0;
        }
    }
    //Scale factor (Einstein frame):
    value_type a(value_type& x)
    {
        return this->hyperbolic ? sinh(this->H_tilde*x)/this->H_tilde
                     : sin(this->H_tilde*x)/this->H_tilde;
    }
    //Derivative of scale factor (Einstein frame)
    value_type ap(value_type& x)
    {
        return this->hyperbolic ? cosh(this->H_tilde*x) : cos(this->H_tilde*x);
    }
};
//==============================================================================


//==============================================================================
//DYNAMIC BOUNDARY CONDITIONS
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
class dynamicBCs
{
protected:
    value_type xi;//Non-minimal coupling
    value_type h;//Scale, h = H/Mp, as a function of the Planck mass
    value_type W0;// V(0)/H^4, ie, potential at 0 in units of H
    potential< value_type >& W;// W(y) = [V(Hy) - V(0)]/H^4, rescaled and
        //shifted potential
    value_type stepError;//Controls how large the solution is allowed to grow
public:
    //Constructors:
	dynamicBCs(value_type XI, value_type H, value_type w0,
               potential< value_type >& V,time_type sError);
    //Member functions:
	DLL_EXPORT virtual time_type operator()(solution_type& BCOld,
                                            solution_type& BCNew) = 0;
};
//------------------------------------------------------------------------------
//Constructor
template< class value_type, class solution_type, class time_type >
dynamicBCs< value_type, solution_type, time_type >::dynamicBCs
    (value_type XI, value_type H, value_type w0, potential< value_type >& V,
     time_type sError) : W(V)
{
	xi = XI;
	h = H;
	W0 = w0;
	stepError = sError;
}
//------------------------------------------------------------------------------
//==============================================================================
//Derived classes for each method:

//Basic BC function. Pure virtual, as it doesn't define the operator()
//Used to implement derived classes which utilise the same basic linearised form
//such as:
// - Basic Solution, no action.
// - Solution with action
// - Solution with derivatives with respect to parameters
template< class value_type , class solution_type , class time_type>
class dynamicBCs_basic : public dynamicBCs< value_type , solution_type ,
                                            time_type>
{
public:
    //Constructors:
	dynamicBCs_basic(value_type XI, value_type H, value_type w0,
                     potential< value_type >& V,time_type sError);
    //Member functions:
    time_type map_BCs(solution_type& BCOld,solution_type& BCNew);
    time_type map_BCs(solution_type& BCOld,solution_type& BCNew,
                      instanton_equation_lineariser<value_type>& d);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type>
dynamicBCs_basic< value_type, solution_type, time_type>::dynamicBCs_basic(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError)
	: dynamicBCs< value_type, solution_type, time_type >(XI,
		H, w0, V, sError) {}
//------------------------------------------------------------------------------
//Basic method, without action:
template< class value_type , class solution_type , class time_type>
class dynamicBCs_basicMethod : public dynamicBCs_basic< value_type ,
                                                        solution_type ,
                                                        time_type>
{
public:
    //Constructors:
	dynamicBCs_basicMethod(value_type XI, value_type H, value_type w0,
                           potential< value_type >& V,time_type sError);
    //Member functions:
    time_type operator()(solution_type& BCOld,solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type>
dynamicBCs_basicMethod< value_type, solution_type, time_type>
    ::dynamicBCs_basicMethod( value_type XI, value_type H, value_type w0,
                             potential< value_type >& V, time_type sError)
	: dynamicBCs_basic< value_type, solution_type, time_type >(XI,
		H, w0, V, sError) {}
//------------------------------------------------------------------------------
//Basic method with action included:
template< class value_type , class solution_type , class time_type>
class dynamicBCs_actionMethod : public dynamicBCs_basic< value_type ,
                                                         solution_type ,
                                                         time_type>
{
public:
    //Constructors:
	dynamicBCs_actionMethod(value_type XI, value_type H, value_type w0,
                            potential< value_type >& V,time_type sError);
    //Member functions:
    time_type operator()(solution_type& BCOld,solution_type& BCNew);
    //Return action as a function of x_tilde (i.e., evaluated in
    //Einstein frame):
    value_type S(value_type& x_tilde,
                 instanton_equation_lineariser<value_type>& d);
    value_type S(value_type& x_tilde,value_type& y0);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type>
dynamicBCs_actionMethod< value_type, solution_type, time_type>
    ::dynamicBCs_actionMethod
    (value_type XI, value_type H, value_type w0, potential< value_type >& V,
     time_type sError)
	: dynamicBCs_basic< value_type, solution_type, time_type >(XI,
		H, w0, V, sError) {}
//------------------------------------------------------------------------------
//==============================================================================


//==============================================================================
//BOUNDARY CONDITIONS FOR FLAT SPACE BOUNCES
//Boundary conditions for flat space:
template< class value_type , class solution_type , class time_type>
class dynamicBCs_flatSpace : public dynamicBCs< value_type, solution_type ,
                                                time_type >
{
public:
    //Constructor:
    dynamicBCs_flatSpace(potential< value_type >& V,time_type sError);
    //Member functions:
    time_type map_BCs(solution_type& BCOld,solution_type& BCNew);
    time_type map_BCs(solution_type& BCOld,solution_type& BCNew,
                      instanton_equation_lineariser<value_type>& d);
};
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
dynamicBCs_flatSpace< value_type, solution_type, time_type>
    ::dynamicBCs_flatSpace(potential< value_type >& V,time_type sError)
    : dynamicBCs< value_type, solution_type, time_type >
    (value_type(0.0),value_type(0.0), value_type(0.0), V, sError) {}
//------------------------------------------------------------------------------
//Boundary conditions for flat space, no action:
template< class value_type , class solution_type , class time_type>
class dynamicBCs_flatSpace_basic : public dynamicBCs_flatSpace< value_type,
                                                                solution_type ,
                                                                time_type >
{
public:
    //Constructor:
    dynamicBCs_flatSpace_basic(potential< value_type >& V,time_type sError);
    //Member functions:
    time_type operator()(solution_type& BCOld,solution_type& BCNew)
    {
        return this->map_BCs(BCOld,BCNew);
    }
};
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
dynamicBCs_flatSpace_basic< value_type, solution_type, time_type>
::dynamicBCs_flatSpace_basic(potential< value_type >& V,time_type sError)
    : dynamicBCs_flatSpace< value_type, solution_type, time_type >(V, sError) {}
//------------------------------------------------------------------------------
//Boundary conditions for flat space, including action:
template< class value_type , class solution_type , class time_type>
class dynamicBCs_flatSpace_action : public dynamicBCs_flatSpace< value_type,
                                                                 solution_type ,
                                                                 time_type >
{
public:
    //Constructor:
    dynamicBCs_flatSpace_action(potential< value_type >& V,time_type sError);
    //Member functions:
    time_type operator()(solution_type& BCOld,solution_type& BCNew)
    {
        //Perform the usual operation to get the first two components
        //(y and yp):
        value_type y0 = BCNew[0];
        value_type chi_step = this->map_BCs(BCOld,BCNew);
        //Now add the action into this mix:
        const value_type PI = boost::math::constants::pi<value_type>();
        const value_type _2p0 = value_type(2.0);
        const value_type _5p0 = value_type(5.0);
        const value_type _96p0 = value_type(96.0);

        value_type x2 = chi_step*chi_step;
        value_type x4 = x2*x2;
        value_type PI2 = PI*PI;
        value_type Vy0 = this->W(y0);
        value_type Vpy0 = this->W.d(y0);
        BCNew[2] = PI2*Vy0*x4/_2p0 + (_5p0/_96p0)*PI2*Vpy0*Vpy0*x4*x2;
    }
};
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
dynamicBCs_flatSpace_action< value_type, solution_type, time_type>
::dynamicBCs_flatSpace_action(potential< value_type >& V,time_type sError)
    : dynamicBCs_flatSpace< value_type, solution_type, time_type >(V, sError) {}
//------------------------------------------------------------------------------
//==============================================================================

//==============================================================================
//BOUNDARY CONDITIONS FOR AdS BOUNCES
#ifdef ANTI_DE_SITTER_FLAT_INSTANTONS
//Basically the same as ode1 so we don't need this:
template< class value_type , class solution_type , class time_type >
class dynamicBCs_taylor_basic : public dynamicBCs< value_type ,
                                                   solution_type , time_type >
{
public:
    int version;
    //Constructors:
	dynamicBCs_taylor_basic(value_type XI, value_type H, value_type w0,
                            potential< value_type >& V,time_type sError,
                            int VERSION);
    //Member functions:
    DLL_EXPORT time_type operator()(solution_type& BCOld, solution_type& BCNew)
    {
        switch(version)
        {
        case 1:
            return basic_bcs(BCOld,BCNew);
            break;
        case 2:
            return action_bcs(BCOld,BCNew);
            break;
        default:
            throw "Invalid version for dynamicBCs_taylor_basic.";
            break;
        }
    }
    //Virtual:
    DLL_EXPORT virtual time_type basic_bcs(solution_type& BCOld,
                                           solution_type& BCNew) = 0;
    DLL_EXPORT virtual time_type action_bcs(solution_type& BCOld,
                                            solution_type& BCNew) = 0;
    //Needed constants:
    const value_type _1p0;// = value_type(1.0);
    const value_type _2p0;// = value_type(2.0);
    const value_type _3p0;// = value_type(3.0);
    const value_type _4p0;// = value_type(4.0);
    const value_type _0p125;// = value_type(0.125);
    const value_type _0p25;// = value_type(0.25);
    const value_type _0p5;// = value_type(0.5);
    const value_type _10p0;// = value_type(10.0);
    const value_type _6p0;// = value_type(6.0);
    const value_type _0p0;// = value_type(0.0);
    const value_type PI;// = boost::math::constants::pi<value_type>();
    //Slots to store intermediates:
    value_type Wy0;
    value_type Wy0p;
    value_type y0pp;
    value_type a0pp;
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCs_taylor_basic< value_type, solution_type, time_type >
    ::dynamicBCs_taylor_basic(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs< value_type, solution_type, time_type >(XI,
		H, w0, V, sError),
		_1p0(1.0),
		_2p0(2.0),
		_3p0(3.0),
        _4p0(4.0),
        _0p125(0.125),
        _0p25(0.25),
        _0p5(0.5),
        _10p0(10.0),
        _6p0(6.0),
        _0p0(0.0),
        PI(boost::math::constants::pi<value_type>())
		 { version = VERSION; }
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type >
class dynamicBCsAdSFlat_taylor : public dynamicBCs_taylor_basic< value_type ,
                                                                 solution_type ,
                                                                 time_type >
{
public:
    //Constructors:
	dynamicBCsAdSFlat_taylor(value_type XI, value_type H, value_type w0,
                             potential< value_type >& V,time_type sError,
                             int VERSION);
    //Member functions:
    DLL_EXPORT time_type basic_bcs(solution_type& BCOld, solution_type& BCNew);
    DLL_EXPORT time_type action_bcs(solution_type& BCOld, solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCsAdSFlat_taylor< value_type, solution_type, time_type >
    ::dynamicBCsAdSFlat_taylor(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs_taylor_basic< value_type, solution_type, time_type >(XI,
		H, w0, V, sError, VERSION) {}
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type >
class dynamicBCs_flatfv_delta : public dynamicBCs_taylor_basic< value_type ,
                                                                 solution_type ,
                                                                 time_type >
{
public:
    //Constructors:
	dynamicBCs_flatfv_delta(value_type XI, value_type H, value_type w0,
                             potential< value_type >& V,time_type sError,
                             int VERSION);
    //Member functions:
    DLL_EXPORT time_type basic_bcs(solution_type& BCOld, solution_type& BCNew);
    DLL_EXPORT time_type action_bcs(solution_type& BCOld, solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCs_flatfv_delta< value_type, solution_type, time_type >
    ::dynamicBCs_flatfv_delta(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs_taylor_basic< value_type, solution_type, time_type >(XI,
		H, w0, V, sError, VERSION) {}
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type >
class dynamicBCsAdSFlat_taylor_delta_a : public dynamicBCs_taylor_basic
                                                < value_type , solution_type ,
                                                                 time_type >
{
public:
    //Constructors:
	dynamicBCsAdSFlat_taylor_delta_a(value_type XI, value_type H, value_type w0,
                             potential< value_type >& V,time_type sError,
                             int VERSION);
    //Member functions:
    DLL_EXPORT time_type basic_bcs(solution_type& BCOld, solution_type& BCNew);
    DLL_EXPORT time_type action_bcs(solution_type& BCOld, solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCsAdSFlat_taylor_delta_a< value_type, solution_type, time_type >
    ::dynamicBCsAdSFlat_taylor_delta_a(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs_taylor_basic< value_type, solution_type, time_type >(XI,
		H, w0, V, sError, VERSION) {}
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type >
class dynamicBCs_dS_fixed : public dynamicBCs_taylor_basic
                                                < value_type , solution_type ,
                                                                 time_type >
{
public:
    //Constructors:
	dynamicBCs_dS_fixed(value_type XI, value_type H, value_type w0,
                             potential< value_type >& V,time_type sError,
                             int VERSION);
    //Member functions:
    DLL_EXPORT time_type basic_bcs(solution_type& BCOld, solution_type& BCNew);
    DLL_EXPORT time_type action_bcs(solution_type& BCOld, solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCs_dS_fixed< value_type, solution_type, time_type >
    ::dynamicBCs_dS_fixed(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs_taylor_basic< value_type, solution_type, time_type >(XI,
		H, w0, V, sError, VERSION) {}
//------------------------------------------------------------------------------
//Flat space fixed background dynamic boundary conditions function:
template< class value_type, class solution_type, class time_type >
class dynamicBCsFlatFixed_taylor : public dynamicBCs_taylor_basic
    < value_type,solution_type, time_type >
{
public:
	//Constructors:
	dynamicBCsFlatFixed_taylor(potential< value_type >& V, time_type sError,
                                int VERSION);
	//Member functions:
	DLL_EXPORT time_type basic_bcs(solution_type& BCOld, solution_type& BCNew);
    DLL_EXPORT time_type action_bcs(solution_type& BCOld, solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCsFlatFixed_taylor< value_type, solution_type, time_type >
::dynamicBCsFlatFixed_taylor(potential< value_type >& V, time_type sError,
                             int VERSION)
	: dynamicBCs_taylor_basic< value_type, solution_type, time_type >
	(value_type(0.0),value_type(0.0),value_type(0.0), V, sError, VERSION) {}
//------------------------------------------------------------------------------
//Dynamic BCs for fixed AdS space:
/*
template< class value_type, class solution_type, class time_type >
class dynamicBCAdSFixed : public dynamicBCs
< value_type, solution_type, time_type >
{
public:
	//Constructors:
	dynamicBCAdSFixed(value_type XI, value_type H, value_type w0,
                   potential< value_type >& V, time_type sError, int VERSION);
	//Member functions:
	DLL_EXPORT time_type Lower(solution_type& BCOld, solution_type& BCNew);
	DLL_EXPORT time_type LowerAction(solution_type& BCOld,
                                     solution_type& BCNew);
	DLL_EXPORT time_type LowerLinear(solution_type& BCOld,
                                     solution_type& BCNew);
	DLL_EXPORT time_type UpperAction(solution_type& BCOld,
                                     solution_type& BCNew);
	DLL_EXPORT time_type Upper(solution_type& BCOld, solution_type& BCNew);
	DLL_EXPORT time_type UpperLinear(solution_type& BCOld,
                                     solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCAdSFixed< value_type, solution_type, time_type >::dynamicBCAdSFixed(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs< value_type, solution_type, time_type >(XI,
		H, w0, V, sError, VERSION) {}
//------------------------------------------------------------------------------
		*/
template< class value_type , class solution_type , class time_type >
class dynamicBCsGravFriction_taylor
 : public dynamicBCs_taylor_basic< value_type ,solution_type , time_type >
{
public:
    //Constructors:
	dynamicBCsGravFriction_taylor(value_type XI, value_type H, value_type w0,
                             potential< value_type >& V,time_type sError,
                             int VERSION);
    //Member functions:
    DLL_EXPORT time_type basic_bcs(solution_type& BCOld, solution_type& BCNew);
    DLL_EXPORT time_type action_bcs(solution_type& BCOld, solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCsGravFriction_taylor< value_type, solution_type, time_type >
    ::dynamicBCsGravFriction_taylor(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs_taylor_basic< value_type, solution_type, time_type >(XI,
		H, w0, V, sError, VERSION) {}
//------------------------------------------------------------------------------
#endif // ANTI_DE_SITTER_FLAT_INSTANTONS
//==============================================================================

//==============================================================================
//------------------------------------------------------------------------------
//For later implementation
#ifdef DE_SITTER_INSTANTONS
//Derived classes for each method:
template< class value_type , class solution_type , class time_type >
class dynamicBCsMethod2 : public dynamicBCs< value_type , solution_type ,
                                             time_type >
{
public:
    //Constructors:
	dynamicBCsMethod2(value_type XI, value_type H, value_type w0, potential
        < value_type >& V,time_type sError,int VERSION);
    //Member functions:
	DLL_EXPORT time_type Lower(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type LowerAction(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type LowerLinear(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type UpperAction(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type Upper(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type UpperLinear(solution_type& BCOld,solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type>
dynamicBCsMethod2< value_type, solution_type, time_type >::dynamicBCsMethod2(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs< value_type, solution_type, time_type >(XI, H, w0, V, sError,
    VERSION) {}
//------------------------------------------------------------------------------
//Derived classes for each method:
template< class value_type , class solution_type , class time_type >
class dynamicBCsMethod3 : public dynamicBCs< value_type , solution_type ,
                                             time_type >
{
public:
    //Constructors:
	dynamicBCsMethod3(value_type XI, value_type H, value_type w0, potential
        < value_type >& V,time_type sError,int VERSION);
    //Member functions:
	DLL_EXPORT time_type Lower(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type LowerAction(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type LowerLinear(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type UpperAction(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type Upper(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type UpperLinear(solution_type& BCOld,solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCsMethod3< value_type, solution_type, time_type >::dynamicBCsMethod3(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs< value_type, solution_type, time_type >(XI,
		H, w0, V, sError, VERSION) {}
//------------------------------------------------------------------------------
//Derived classes for each method:
template< class value_type , class solution_type , class time_type >
class dynamicBCsMethod4 : public dynamicBCs< value_type , solution_type ,
                                             time_type >
{
public:
    //Constructors:
	dynamicBCsMethod4(value_type XI, value_type H, value_type w0, potential
        < value_type >& V,time_type sError,int VERSION);
    //Member functions:
	DLL_EXPORT time_type Lower(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type LowerAction(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type LowerLinear(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type UpperAction(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type Upper(solution_type& BCOld,solution_type& BCNew);
	DLL_EXPORT time_type UpperLinear(solution_type& BCOld,solution_type& BCNew);
};
//------------------------------------------------------------------------------
//Constructor (essentially just calls the base class destructor).
template< class value_type, class solution_type, class time_type >
dynamicBCsMethod4< value_type, solution_type, time_type >::dynamicBCsMethod4(
	value_type XI, value_type H, value_type w0, potential< value_type >& V,
	time_type sError, int VERSION)
	: dynamicBCs< value_type, solution_type, time_type >(XI,
		H, w0, V, sError, VERSION) {}
//------------------------------------------------------------------------------
#endif // DE_SITTER_INSTANTONS
//==============================================================================



#endif
