//Function to evaluate the hypergeometric function
#include<boost/math/special_functions.hpp>
//For testing purposes only (outputs intermediate data to a file):
//#include <iostream>
//#include <fstream>

namespace hypergeometric{

template<class value_type>
value_type inverse_gamma(value_type z)
{
    using namespace boost::math;
    //Return 1/Gamma(z). which is actually a well behaved function because the gamma
    //function has no zeros.
    //NB - does not return the functional inverse of the gamma function, gamma^{-1}(z).
    //It returns the multiplicative inverse.

    //Poles are located at negative integers:
    value_type _0p0 = value_type(0.0);
    value_type _1p0 = value_type(1.0);
    value_type epsilon = std::numeric_limits<value_type>::epsilon();
    /*
    std::cout << "\nComputing inverse gamma: ";
    std::cout << "\nz = " << z;
    std::cout << "\nabs(round(z) - z) = " << abs(round(z) - z);
    std::cout << "\nepsilon = " << epsilon;
    std::cout << "\nz <= _0p0 = " << (z <= _0p0);
    std::cout << "\nabs(round(z) - z)< epsilon = " << (abs(round(z) - z)< epsilon);
    */
    if(z <= _0p0 && abs(round(z) - z)< epsilon)
    {
        //Pole of the gamma function = zero of the inverse gamma function
        return _0p0;
    }
    else
    {
        return _1p0/tgamma(z);
    }
}

template<class value_type>
value_type sinc(value_type x)
{
    using namespace boost::math;
    value_type epsilon = std::numeric_limits<value_type>::epsilon();
    value_type _1p0 = value_type(1.0);
    value_type PI = boost::math::constants::pi<value_type>();
    value_type PIx = PI*x;
    if(abs(x)  < epsilon)
    {
        return _1p0;
    }
    else
    {
        return sin(PIx)/(PIx);
    }
}

/*
int factorial(int n)
{
    if(n == 0)
    {
        return 1;
    }
    else if(n > 0)
    {
        int temp = 1;
        for(int i = 0;i < n;i++)
        {
            temp *= (n - i);
        }
        return temp;
    }
    else
    {
        throw "Error: argument to factorial(n) must be a non-negative integer.";
    }
}
*/

template<class value_type>
value_type sign_of(int m)
{
    if(m % 2 == 0)
    {
        return value_type(1.0);
    }
    else
    {
        return value_type(-1.0);
    }
}

//Most basic method possible: naiively sum the Taylor series until it converges:
template<class value_type>
value_type hyperGeometric2F1_basic_series(const value_type a,const value_type b,const value_type c,const value_type x)
{
    value_type error;
    value_type result;
    int iterations;
    hyperGeometric2F1_basic_series(a,b,c,x,error,iterations,result);
    return result;
}

//Overloaded version, accomplishing the same as the above, but recording error and number of iterations. Useful for comparisons.
template<class value_type>
void hyperGeometric2F1_basic_series(const value_type a,const value_type b,const value_type c,const value_type x,value_type& error,int& iterations, value_type& result)
{
    //std::ofstream testFile;
    //testFile.open("test_data.txt");


    //Sanity checks on input:
    value_type _0p0 = value_type(0.0);
    value_type _1p0 = value_type(1.0);//Evaluate this only once, to avoid repeated calls to the constructor.
    value_type EPS = std::numeric_limits< value_type >::epsilon();
    if(abs(c) < EPS)
    {
        throw "Argument 'c' to 2F1 is too small.";
    }
    if( c < _0p0 && abs(round(c) - c) < EPS)
    {
        throw "Argument 'c' to 2F1 is too close to a negative integer.";
    }
    if( abs(x) > _1p0)
    {
        throw "Argument 'x' to 2F1 is outside radius of convergence.";
    }

    //Sum the series until the error is less than epsilon:
    int Nmax = 10000;//Max number of iterations
    value_type k = _0p0;
    int i = 1;
    value_type sum_pos = _1p0;//Sum of positive terms
    value_type sum_neg = _0p0;//Sum of negative terms
        //Evaluate positive and negative terms in sequence separately, so that we
        //avoid losing precision
    value_type diff = _1p0;//Next term in series to be added.
    do
    {
        diff *= (a + k)*(b + k)*x/( (c + k)*( k + _1p0 ) );
        if(diff > _0p0)
        {
            sum_pos += diff;
        }
        else if(diff == _0p0)
        {
            //Terminate sequence (polynomial special case)
            break;
        }
        else
        {
            //diff < 0:
            sum_neg -= diff;
        }
        k += _1p0;
        //For testing:
        //testFile.precision(std::numeric_limits<value_type>::digits10);
        //testFile << std::endl << sum_pos - sum_neg;
    }
    while(i++ < Nmax && abs(diff) > EPS*(sum_pos - sum_neg));
    //testFile << "Completed " << i << " iterations. diff = " << diff << std::endl;
    if(i >= Nmax && abs(diff) > EPS*(sum_pos - sum_neg))
    {
        throw "Error: Convergence of 2F1 Taylor series too slow.";
    }
    //testFile.close();
    error = abs(diff);
    iterations = i;
    result = sum_pos - sum_neg;
}

//Sum Hypergeometric series using the z/(z-1) method
template<class value_type>
void hyperGeometric2F1_zo_zm1_series(const value_type a,const value_type b,const value_type c,const value_type x,value_type& error,int& iterations, value_type& result)
{
    //Sanity checks:
    value_type epsilon = std::numeric_limits<value_type>::epsilon();
    value_type _1p0 = value_type(1.0);
    if(abs(x - _1p0) < epsilon)
    {
        throw "Error: too close to z = 1 pole.";
    }

    //Setup a few needed variables:
    value_type base = _1p0 - x;
    value_type exponent = -a;
    value_type factor = pow(base,exponent);
    value_type bp = c - b;
    value_type new_param = x/(x - _1p0);
    //Call hypergeometric series, with the transformed paramters:
    hyperGeometric2F1_basic_series(a,bp,c,new_param,error,iterations,result);
    //Multiply by pre-factor resulting from transformation:
    result *= factor;
}
//Return value version:
template<class value_type>
value_type hyperGeometric2F1_zo_zm1_series(const value_type a,const value_type b,const value_type c,const value_type x)
{
    value_type error;
    value_type result;
    int iterations;
    hyperGeometric2F1_zo_zm1_series(a,b,c,x,error,iterations,result);
    return result;
}

template<class value_type>
value_type Eeps(value_type eps,value_type z)
{
    //Define needed functions:
    using namespace boost::math;

    const value_type _0p0 = value_type(0.0);
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();
    const value_type _1p0 = value_type(1.0);
    value_type arg = eps*z;
    const value_type sum_max = value_type(10000);
    value_type result = _0p0;
    value_type an = z;
    value_type np1 = value_type(2.0);
    do
    {
        result += an;
        an *= eps*z/np1;
        np1 += _1p0;
    }
    while(abs(an) < epsilon && np1 < sum_max);
    return result;
}

template<class value_type>
value_type Heps(value_type eps,value_type z)
{
    //Define needed functions:
    using namespace boost::math;

    value_type zarg = z + eps;
    value_type _0p5 = value_type(0.5);
    value_type _1p0 = value_type(1.0);
    value_type PI = boost::math::constants::pi<value_type>();
    value_type epsilon = std::numeric_limits<value_type>::epsilon();
    if(abs(z) < epsilon)
    {
        return -_1p0/eps;
    }
    else if(zarg < _0p5 && z < _0p5)
    {
        //Valid if the real parts of z or z + eps are less than 1/2:
        //Intermediate variables:
        value_type PIz = PI*z;
        value_type PIeps = PI*eps;
        value_type meps = -eps;
        value_type _1mzmeps = _1p0 - z - eps;
        value_type _1mz = _1p0 - z;
        value_type PIepso2 = PI*eps*_0p5;
        value_type epso2 = eps*_0p5;
        //Intermediate functions:
        value_type H_meps_1mz = expm1(lgamma(_1mzmeps) - lgamma(_1mz))/(meps);
        value_type H_meps_zpeps = (cos(PIeps) + sin(PIeps)/tan(PIz))*H_meps_1mz + PI*PIepso2*sinc(epso2)*sinc(epso2) - PI*sinc(eps)/(tan(PIz));
        //Final result:
        return H_meps_zpeps/(_1p0 - eps*H_meps_zpeps);
    }
    else
    {
        //Valid for all other cases:
        return expm1(lgamma(zarg) - lgamma(z))/eps;
    }

}

template<class value_type>
value_type Geps(value_type eps,value_type z)
{
    //Define needed functions:
    using namespace boost::math;

    value_type n = abs(round(z));//Technically defined in terms of the real
        //part of z, but we're only using this code for real z.
    //std::cout << "\neps = " << eps;
    //std::cout << "\nz = " << z;
    value_type zarg = z + eps;
    value_type m = abs(round(zarg));
    value_type _0p1 = value_type(0.1);
    value_type _1p0 = value_type(1.0);
    value_type _0p0 = value_type(0.0);
    value_type epsilon = std::numeric_limits<value_type>::epsilon();
    bool z_is_negative_int = (z <= _0p0) && (abs(z - round(z)) < epsilon);
    bool zarg_is_negative_int = (zarg <= _0p0) && (abs(zarg - round(zarg)) < epsilon);
    //Special case when eps = 0:
    if(abs(eps) < epsilon)
    {
        if(z_is_negative_int && abs(z) > _0p0)
        {
            //std::cout << "Geps = 2";
            value_type sign = _1p0;
            int n = -round(z).template convert_to<int>();
            if(n % 2 == 0)
            {
                sign = -_1p0;
                //Computes (-1)^{n+1}
            }
            return sign*tgamma(value_type(n) + _1p0);
        }
        else
        {
            //std::cout << "Geps = 1";
            return polygamma(0,z)*inverse_gamma(z);
        }
    }
    //Special case to deal with numerical instabilities for small
    //but non-zero eps:
    else if(zarg_is_negative_int)
    {
        //1/Gamma(z + eps)->0, so we are left with 1/(eps*Gamma(z)))
        return inverse_gamma(z)/eps;
    }
    else if(z_is_negative_int)
    {
        return -inverse_gamma(zarg)/eps;
    }
    else if(abs(eps) < _0p1)
    {
        //std::cout << "Geps = 3";
        //Two different (equivalent) versions, to avoid problems when z OR z + eps
        //are near to a negative integer:
        if(abs(z + n) < abs(z + eps + m))
        {
            //std::cout << "Geps = 4";
            return Heps(eps,z)*inverse_gamma(zarg);
        }
        else
        {
            //std::cout << "Geps = 5";
            value_type meps = -eps;
            return Heps(meps,zarg)*inverse_gamma(z);
        }
    }
    //General case:
    else
    {
        //std::cout << "Geps = 6";
        //Otherwise, return the general formula, which shouldn't be subject to
        //numerical instability:
        return (inverse_gamma(z) - inverse_gamma(zarg))/eps;
    }
}



template<class value_type>
value_type pochhammer(value_type q,int m)
{
    value_type _1p0 = value_type(1.0);
    value_type _0p0 = value_type(0.0);
    value_type M = value_type(m);
    if(m == 0)
    {
        return _1p0;
    }
    else if(m > 0)
    {
        value_type temp = _1p0;
        value_type I = _0p0;
        for(int i = 0;i < m;i++)
        {
            temp *= (q + I);
            I += _1p0;
        }
        return temp;
    }
    else
    {
        throw "Error: argument m to pochhammer symbol must be a non-negative integer.";
    }
}

//Horrible function appearing in the 1/z expansion:
template<class value_type>
value_type Peps(value_type eps,value_type z,int m)
{
    //Define needed functions:
    using namespace boost::math;


    value_type _0p0 = value_type(0.0);
    value_type _1p0 = value_type(1.0);
    value_type result;
    value_type N0 = -round(z);
    int n0 = N0.template convert_to<int>();
    bool delta_non_vanishing = n0 >= 0 && n0 < m;
    value_type first_term = _0p0;//Will contribute if delta_non_vanishing, irrelevant otherwise.
    value_type second_term = _0p0;
    value_type N = _0p0;//Loop iterator;
    if(eps == _0p0)
    {
        if(delta_non_vanishing)
        {
            first_term = _1p0;

            //Product (and sum for second term), skipping n0:
            N = _0p0;
            for(int n = 0;n<n0;n++)
            {
                first_term *= (z + N);
                second_term += _1p0/(z + N);
                N+= _1p0;
            }
            N = N0 + _1p0;
            for(int n = n0 + 1;n<m;n++)
            {
                first_term *= (z + N);
                second_term += _1p0/(z + N);
                N+= _1p0;
            }
            return first_term + pochhammer(z,m)*second_term;
        }
        else
        {
            //Sum the terms for the second term only. This time, delta vanishes,
            //which means n0 is outside [0,m-1]. Hence there is no need to skip it!
            N = _0p0;
            for(int n = 0;n<m;n++)
            {
                second_term += _1p0/(z + N);
                N+= _1p0;
            }
            return pochhammer(z,m)*second_term;
        }
    }
    else
    {
        if(delta_non_vanishing)
        {
            first_term = _1p0;

            //Product (and sum for second term), skipping n0:
            N = _0p0;
            for(int n = 0;n<n0;n++)
            {
                first_term *= (z + N + eps);
                value_type arg = eps/(z + N);
                second_term += log1p(arg);
                N+= _1p0;
            }
            N = N0 + _1p0;
            for(value_type n = n0 + 1;n<m;n++)
            {
                first_term *= (z + N + eps);
                value_type arg = eps/(z + N);
                second_term += log1p(arg);
                N+= _1p0;
            }
            return first_term + pochhammer(z,m)*expm1(second_term)/eps;
        }
        else
        {
            //Sum the terms for the second term only. This time, delta vanishes,
            //which means n0 is outside [0,m-1]. Hence there is no need to skip it!
            N = _0p0;
            for(int n = 0;n<m;n++)
            {
                value_type arg = eps/(z + N);
                second_term += log1p(arg);
                N+= _1p0;
            }
            return pochhammer(z,m)*expm1(second_term)/eps;
        }
    }
}

//Sum 1-z parameter version:
template<class value_type>
void hyperGeometric2F1_1mz_series(const value_type a,const value_type b,const value_type c,const value_type x,value_type& error,int& iterations, value_type& result)
{
    //Define needed functions:
    using namespace boost::math;

    //std::cout << "\nOk here. 1";
    //Transform if c - a - b is negative:
    value_type _1p0 = value_type(1.0);
    value_type _0p0 = value_type(0.0);
    value_type _0p1 = value_type(0.1);
    value_type epsilon = std::numeric_limits<value_type>::epsilon();
    value_type factor = _1p0;
    value_type R = c - a - b;
    value_type ap = a;
    value_type bp = b;
    value_type cp = c;
    value_type base = _1p0 - x;
    value_type I = _0p0;
    if(R < _0p0)
    {
        ap = c - a;
        bp = c - b;
        factor = pow(base,R);
        R = cp - ap - bp;
    }
    //std::cout << "\nfactor = " << factor;

    //std::cout << "\nOk here. 2";
    //Integer power appearing in the fomula:
    value_type M = round(R);
    int m = M.template convert_to<int>();
    //complex number appearing in the formula:
    value_type eps = R - M;//How much R differs from the nearest integer.

    //Actual result is:
    //F(z) = (-1)^m*(A(z) + B(z))/sinc(eps)
    //Need to compute A and B separately from power series.
    //Start from first coefficient in each series:
    value_type a0,b0,c0;

    //various intermediaries we will need:
    value_type _apm = ap + M;
    value_type _bpm = bp + M;
    value_type _1o_gamma_apm = inverse_gamma(_apm);
    value_type _1o_gamma_bpm = inverse_gamma(_bpm);
    value_type gamma_m = tgamma(M);
    value_type _1mMmeps = _1p0 - M - eps;
    value_type _apMpeps = ap + M + eps;
    value_type _bpMpeps = bp + M + eps;
    value_type _1o_gamma_apMpeps = inverse_gamma(_apMpeps);
    value_type _1o_gamma_bpMpeps = inverse_gamma(_bpMpeps);

    value_type _1meps = _1p0 - eps;
    value_type _mp1 = M + _1p0;
    value_type _1o_gamma_1meps = inverse_gamma(_1meps);
    value_type _1o_gamma_mp1 = inverse_gamma(_mp1);
    value_type _1o_gamma_a = inverse_gamma(ap);
    value_type _1o_gamma_b = inverse_gamma(bp);
    value_type _mp1peps = M + _1p0 + eps;
    value_type _1o_gamma_mp1peps = inverse_gamma(_mp1peps);
    value_type gamma_c = tgamma(cp);
    value_type meps = -eps;
    value_type log_base = log(base);


    //std::cout << "\nOk here. 3";
    //Compute a0:
    //std::cout << "\neps = " << eps;
    //std::cout << "\nm = " << m;
    if(eps == _0p0 && m!=0)//i.e., integer value
    {
        //std::cout << "\n1o_gamma_bpm = " << _1o_gamma_bpm;
        //std::cout << "\n1o_gamma_apm = " << _1o_gamma_apm;
        a0 = sign_of<value_type>(m)*gamma_m*(gamma_c*(_1o_gamma_bpm*_1o_gamma_apm));
    }
    else if(m != 0)
    {
        /*
        std::cout << "\n1o_gamma_apMpeps = " << _1o_gamma_apMpeps;
        std::cout << "\n1o_gamma_bpMpeps = " << _1o_gamma_bpMpeps;
        std::cout << "\n_1mMmeps = " << _1mMmeps;
        std::cout << "\n_1o_gamma_1mMmeps = " << inverse_gamma(_1mMmeps);
        std::cout << "\ngamma(c) = " << tgamma(cp);
        std::cout << "\ncp = " << cp;
        */
        a0 = (tgamma(cp)*inverse_gamma(_1mMmeps))*(_1o_gamma_apMpeps*_1o_gamma_bpMpeps)/(eps);
    }
    else
    {
        a0 = _0p0;
    }
    /*
    std::cout << "a0 = " << a0;

    std::cout << "\nOk here. 4";
    std::cout << "\na = " << ap;
    std::cout << "\nb = " << bp;
    std::cout << "\nc = " << cp;
    */
    //Compute b0:
    if(abs(eps) > _0p1)
    {
        /*
        std::cout << "\n (a)_m = " << pochhammer(ap,m);
        std::cout << "\n (b)_m = " << pochhammer(bp,m);
        std::cout << "\n_1o_gamma_1meps" << _1o_gamma_1meps;
        std::cout << "\n_1o_gamma_apMpeps" << _1o_gamma_apMpeps;
        std::cout << "\n_1o_gamma_bpMpeps" << _1o_gamma_bpMpeps;
        std::cout << "\n_1o_gamma_mp1" << _1o_gamma_mp1;
        std::cout << "\n(1-z)^eps" << pow(base,eps);
        std::cout << "\n_1o_gamma_a" << _1o_gamma_a;
        std::cout << "\n_1o_gamma_b" << _1o_gamma_b;
        std::cout << "\n_1o_gamma_mp1peps" << _1o_gamma_mp1peps;
        std::cout << "\ngamma_c" << gamma_c;
        */

        b0 = ( pochhammer(ap,m)*pochhammer(bp,m)*(_1o_gamma_1meps*_1o_gamma_apMpeps*_1o_gamma_bpMpeps*_1o_gamma_mp1)  - pow(base,eps)*(_1o_gamma_a*_1o_gamma_b*_1o_gamma_mp1peps) )*gamma_c/eps;
        for(int i = 0; i < m;i++)
        {
            b0 *= base;
        }
    }
    else
    {
        b0 = ( (Geps(meps,_1p0)*_1o_gamma_mp1 + Geps(eps,_mp1) )*(_1o_gamma_apMpeps*_1o_gamma_bpMpeps)
                 - (Geps(eps,_apm)*_1o_gamma_bpMpeps + Geps(eps,_bpm)*_1o_gamma_apm )*(_1o_gamma_mp1peps)
                    - Eeps(eps,log_base)*(_1o_gamma_apm*_1o_gamma_bpm*_1o_gamma_mp1peps))*gamma_c;
        I = _0p0;
        for(int i = 0; i < m;i++)
        {
            b0 *= base*(ap + I)*(bp + I);
            I += _1p0;
        }
    }
    //std::cout << "\nb0 = " << b0;

    //std::cout << "\nOk here. 5";
    //Compute g0:
    c0 = gamma_c*(_1o_gamma_apMpeps*_1o_gamma_bpMpeps*_1o_gamma_1meps*_1o_gamma_mp1);
    I = _0p0;
    for(int i = 0; i < m;i++)
    {
        c0 *= base*(ap + I)*(bp + I);
        I += _1p0;
    }
    //std::cout << "\nc0 = " << c0;

    //std::cout << "\nOk here. 6";
    //Compute A(z) series (finite):
    value_type A = a0;
    value_type alpha_n = a0;
    I = _0p0;
    for(int i = 0;i < m - 1;i++)
    {
        alpha_n *= base*(ap + I)*(bp + I)/((I + _1p0)*(_1p0 - M - eps + I));
        //std::cout << "\nalpha_n = " <<alpha_n;
        A += alpha_n;
        I += _1p0;
    }

    //std::cout << "\nOk here. 7";
    //Compute B(z) series, terminating when convergence is detected:
    value_type B = b0;
    value_type Nmax = value_type(10000);
    value_type n = _0p0;
    value_type beta_n = b0;
    value_type bn = b0;
    value_type cn = c0;
    value_type fn;
    value_type fn_data[] = {(ap + M + eps),(bp + M + eps),(M + _1p0 + eps),(_1p0)};
    value_type gn;
    value_type gn_data[] = {(ap + M),(bp + M),(M + _1p0),(ap + M),(bp + M),(ap + M + eps),(bp + M + eps),(_1p0),(M + _1p0 + eps),(_1p0 - eps)};
    value_type hn;
    value_type hn_data[] = {(ap + M),(bp + M),(M + _1p0),(_1p0 - eps)};
    do
    {
        //Intermediate functions:
        //Directly re-computing fn and hn involves a lot of additions/subtractions. We can speed this process up
        //by storing parts of it and updating them, reducing, for example, the three additions required to compute
        // (ap + M + n + eps) in each step to a single addition, fn1 += _1p0. This will shave off a fair amount of time since
        //we expect to repeat these steps hundreds of times.
        //fn = (ap + M + n + eps)*(bp + M + n + eps)/((M + n + _1p0 + eps)*(n + _1p0));
        //gn = ( (ap + M + n)*(bp + M + n)/(M + n + _1p0) - (ap + M + n) - (bp + M + n) - eps + (ap + M + n + eps)*(bp + M + n + eps)/(n + _1p0) )/((n + M + _1p0 + eps)*(n + _1p0 - eps));
        fn = fn_data[0]*fn_data[1]/(fn_data[2]*fn_data[3]);
        gn = ( gn_data[0]*gn_data[1]/gn_data[2] - gn_data[3] - gn_data[4] - eps + gn_data[5]*gn_data[6]/gn_data[7] )/(gn_data[8]*gn_data[9]);
        //Update beta_n = b_n (1-z)^n using a cunning trick:
        beta_n *= base*(fn + gn*cn/bn);
        //Add the next term to the Taylor series:
        B += beta_n;
        //Update bn and cn:
        bn = fn*bn + gn*cn;
        //cn *= (ap + M + n)*(bp + M + n)/((n + M + _1p0)*(n + _1p0 - eps));
        cn *= hn_data[0]*hn_data[1]/(hn_data[2]*hn_data[3]);


        //Advance n for the next step:
        n += _1p0;
        //Advance fn,hn,cn:
        for(int i = 0; i < 4;i++)
        {
            fn_data[i] += _1p0;
            hn_data[i] += _1p0;
        }
        for(int i = 0; i < 10;i++)
        {
            gn_data[i] += _1p0;
        }
        //std::cout << "\nbeta_n = " << beta_n << "\tB = " << B;
    }
    while(abs(beta_n) > abs(epsilon*B) && n < Nmax);

    //std::cout << "\nOk here. 8";

    //Compute result:
    value_type F;
    value_type PI = boost::math::constants::pi<value_type>();
    value_type PIeps = PI*eps;
    F = sign_of<value_type>(m)*(A + B)/sinc(eps);
    //F = (-1)^m(-z)^{-a}*(A(z) + B(z))/sinc(eps);



    //Output:
    result = factor*F;//To account for transformation, if R was negative.
    error = abs(beta_n);
    iterations = n.template convert_to<int>();
    //std::cout << "\nOk here. 9";
}
//Sum 1-z parameter version, without calling for error/no. of iterations.
template<class value_type>
value_type hyperGeometric2F1_1mz_series(const value_type a,const value_type b,const value_type c,const value_type x)
{
    value_type error;
    int iterations;
    value_type result;
    hyperGeometric2F1_1mz_series(a,b,c,x,error,iterations,result);
    return result;
}

//Sum 1/z method:
template<class value_type>
void hyperGeometric2F1_1oz_series(const value_type a,const value_type b,const value_type c,const value_type x,value_type& error,int& iterations, value_type& result)
{
    //Define needed functions:
    using namespace boost::math;
    //std::cout << "\nOk here. 1";
    //Transform if b - a is negative:
    value_type _1p0 = value_type(1.0);
    value_type _0p0 = value_type(0.0);
    value_type _0p1 = value_type(0.1);
    value_type I = _0p0;//Used as a value_type loop iterator.
    value_type epsilon = std::numeric_limits<value_type>::epsilon();
    //Sanity check:
    if(x == _0p0)
    {
        throw "Error: summing hypergeometric series using 1/z method not possible if z = 0.";
    }
    //std::cout << "\nOk here. 2";

    value_type factor = _1p0;
    value_type R = b - a;
    value_type ap = a;
    value_type bp = b;
    value_type cp = c;
    value_type base = _1p0/x;
    if(R < _0p0)
    {
        //Swap a and b:
        ap = b;
        bp = a;
        R = bp - ap;
    }
    //std::cout << "\nap = " << ap;
    //std::cout << "\nbp = " << bp;
    //std::cout << "\ncp = " << cp;

    //std::cout << "\nOk here. 3";
    //Integer power appearing in the formula:
    //std::cout << "\nR = " << R;
    value_type M = round(R);
    //std::cout << "\nM = " << M;
    int m = M.template convert_to<int>();
    //std::cout << "\nm = " << m;    //complex number appearing in the formula:
    value_type eps = R - M;//How much R differs from the nearest integer.
    //std::cout << "\neps = " << eps;

    //Actual result is:
    //F(z) = (-1)^m*(-z)^{-a}*(A(z) + B(z))/sinc(eps)
    //Need to compute A and B separately from power series.
    //Start from first coefficient in each series:
    value_type a0,b0,c0;

    //various intermediaries we will need, regardless of
    //the direction the if...else statements pick. Computing gamma
    //functions is expensive, so we try to keep them to a minimum by only computing
    //some of them within the if...else statements that require them.

    //Needed variables combinations:
    value_type _apm = ap + M;
    value_type _1mMmeps = _1p0 - M - eps;
    value_type _apMpeps = ap + M + eps;
    value_type _1meps = _1p0 - eps;
    value_type _mp1 = M + _1p0;
    value_type _mp1peps = M + _1p0 + eps;
    value_type meps = -eps;
    value_type cma = cp - ap;
    value_type cmameps = cp - ap - eps;
    value_type _1mcpa = _1p0 - cp + ap;
    value_type _1mcpapeps = _1p0 - cp + ap + eps;
    value_type mbase = -base;

    //Needed functions:
    value_type _1o_gamma_apMpeps = inverse_gamma(_apMpeps);
    value_type _1o_gamma_1meps = inverse_gamma(_1meps);
    value_type _1o_gamma_mp1 = inverse_gamma(_mp1);
    value_type gamma_c = tgamma(cp);
    value_type _1o_gamma_cma = inverse_gamma(cma);

    value_type pochhammer_1mppapeps = pochhammer(_1mcpapeps,m);
    //std::cout << "\nOk here. 4";

    //std::cout << "\nap = " << ap;
    value_type _1o_gamma_a = inverse_gamma(ap);
    //std::cout << "\ncmameps = " << cmameps;
    value_type _1o_gamma_cmameps = inverse_gamma(cmameps);
    //std::cout << "\n_mp1peps = " << _mp1peps;
    value_type _1o_gamma_mp1peps = inverse_gamma(_mp1peps);



    //std::cout << "\nOk here. 5";
    //Compute a0:
    if(eps == _0p0 && m != 0)//i.e., integer value
    {
        //Needed variables combinations:


        //Needed gamma functions:
        value_type gamma_m = tgamma(M);
        value_type _1o_gamma_apm = inverse_gamma(_apm);
        a0 = sign_of<value_type>(m)*gamma_m*gamma_c*_1o_gamma_apm*_1o_gamma_cma;
    }
    else if(m != 0)
    {
        //Needed variables combinations:


        //Needed functions:
        value_type _1o_gamma_1mMmeps = inverse_gamma(_1mMmeps);

        a0 = gamma_c*(_1o_gamma_1mMmeps*_1o_gamma_apMpeps*_1o_gamma_cma)/(eps);
    }
    else
    {
        a0 = _0p0;
    }
    //std::cout << "\nOk here. 6";

    //Compute b0:
    if(abs(eps) > _0p1)
    {
        /*
        std::cout << "\npochhammer(ap,m) = " << pochhammer(ap,m);
        std::cout << "\npochhammer(_1mcpa,m) = " << pochhammer(_1mcpa,m);
        std::cout << "\n_1o_gamma_cma = " << _1o_gamma_cma;
        std::cout << "\n_1o_gamma_apMpeps = " << _1o_gamma_apMpeps;
        std::cout << "\n_1o_gamma_1meps = " << _1o_gamma_1meps;
        std::cout << "\n_1o_gamma_mp1 = " << _1o_gamma_mp1;
        std::cout << "\npochhammer_1mppapeps = " << pochhammer_1mppapeps;
        std::cout << "\npow(mbase,meps) = " << pow(mbase,meps);
        std::cout << "\n_1o_gamma_a = " << _1o_gamma_a;
        std::cout << "\n_1o_gamma_cmameps = " << _1o_gamma_cmameps;
        std::cout << "\n_1o_gamma_mp1peps = " << _1o_gamma_mp1peps;
        std::cout << "\ngamma_c = " << gamma_c;
        std::cout << "\neps = " << eps;
        std::cout << "\nmbase = " << mbase;
        std::cout << "\nmeps = " << meps;
        */

        b0 = ( pochhammer(ap,m)*pochhammer(_1mcpa,m)*(_1o_gamma_cma*_1o_gamma_apMpeps*_1o_gamma_1meps*_1o_gamma_mp1)  - pochhammer_1mppapeps*pow(mbase,eps)*(_1o_gamma_a*_1o_gamma_cmameps*_1o_gamma_mp1peps) )*gamma_c/eps;
        for(int i = 0; i < m;i++)
        {
            b0 *= base;
        }
    }
    else
    {
        //Needed variables combinations:
        value_type mz = -x;

        //Needed functions:
        value_type logmz = log(mz);
        value_type _1o_gamma_apm = inverse_gamma(_apm);
        value_type _1o_gamma_mp1peps = inverse_gamma(_mp1peps);
        value_type _1o_gamma_cmameps = inverse_gamma(cmameps);

        /*
        std::cout << "\n_1o_gamma_apm = " << _1o_gamma_apm;
        std::cout << "\n_1o_gamma_mp1peps = " << _1o_gamma_mp1peps;
        std::cout << "\n_1o_gamma_cmameps = " << _1o_gamma_cmameps;
        std::cout << "\npochhammer_1mppapeps = " << pochhammer_1mppapeps;
        std::cout << "\nGeps(meps,_1p0) = " << Geps(meps,_1p0);
        std::cout << "\nPeps(eps,_1mcpa,m) = " << Peps(eps,_1mcpa,m);
        std::cout << "\n_1o_gamma_1meps = " << _1o_gamma_1meps;
        std::cout << "\n_1o_gamma_cma = " << _1o_gamma_cma;
        std::cout << "\n_1o_gamma_apMpeps = " << _1o_gamma_apMpeps;
        std::cout << "\n_1o_gamma_mp1 = " << _1o_gamma_mp1;
        std::cout << "\nGeps(eps,_apm) = " << Geps(eps,_apm);
        std::cout << "\nGeps(meps,cma) = " << Geps(meps,cma);
        std::cout << "\nEeps(meps,logmz) = " << Eeps(meps,logmz);
        std::cout << "\ngamma_c = " << gamma_c;
        */

        b0 = ( (pochhammer_1mppapeps*Geps(meps,_1p0) - Peps(eps,_1mcpa,m)*_1o_gamma_1meps )*(_1o_gamma_cma*_1o_gamma_apMpeps*_1o_gamma_mp1)
               +  pochhammer_1mppapeps*( (Geps(eps,_mp1))*(_1o_gamma_cma*_1o_gamma_apMpeps) - (Geps(eps,_apm))*(_1o_gamma_cma*_1o_gamma_mp1peps)
               - (Geps(meps,cma) - Eeps(meps,logmz)*_1o_gamma_cmameps )*(_1o_gamma_mp1peps*_1o_gamma_apm) )  )*gamma_c;
               I = _0p0;
        for(int i = 0; i < m;i++)
        {
            b0 *= base*(ap + I);
            I += _1p0;
        }
    }
    //std::cout << "\nb0 = " << b0;

    //std::cout << "\nOk here. 7";
    //Compute c0:
    c0 = gamma_c*(_1o_gamma_apMpeps*_1o_gamma_cma*_1o_gamma_1meps*_1o_gamma_mp1);
    //std::cout << "M = " << M;
    I = _0p0;
    for(int i = 0; i < m;i++)
    {
        //std::cout << "\nc0 = " << c0;
        c0 *= base*(ap + I)*(_1p0 - cp + ap + I);
        I += _1p0;
    }

    //std::cout << "\nOk here. 8";
    //Compute A(z) series (finite):
    value_type A = a0;
    value_type alpha_n = a0;
    I = _0p0;
    for(int i = 0;i < m - 1;i++)
    {
        alpha_n *= base*(ap + I)*(_1p0 - cp + ap + I)/((I + _1p0)*(_1p0 - M - eps + I));
        //std::cout << "\nalpha_n = " << alpha_n;
        A += alpha_n;
        I += _1p0;
    }
    //std::cout << "\nOk here. 9";

    //Compute B(z) series, terminating when convergence is detected:
    value_type B = b0;
    value_type Nmax = value_type(10000);
    value_type n = _0p0;
    value_type beta_n = b0;
    value_type bn = b0;
    value_type cn = c0;
    value_type fn;
    value_type fn_data [4];
    fn_data[0] = (ap + M + eps);
    fn_data[1] = (_1p0 - cp + ap + M + eps);
    fn_data[2] = (M + _1p0 + eps);
    fn_data[3] = _1p0;
    value_type gn;
    value_type gn_data[10];
    gn_data[0] = (ap + M);
    gn_data[1] = (_1p0 - cp + ap + M);
    gn_data[2] = (M + _1p0);
    gn_data[3] = (ap + M);
    gn_data[4] = (_1p0 - cp + ap + M);
    gn_data[5] = (ap + M + eps);
    gn_data[6] = (_1p0 - cp + ap + M + eps);
    gn_data[7] = (_1p0);
    gn_data[8] = (M + _1p0 + eps);
    gn_data[9] = (_1p0 - eps);
    value_type hn_data[4];
    hn_data[0] = (ap + M);
    hn_data[1] = (_1p0 - cp + ap + M);
    hn_data[2] = (M + _1p0);
    hn_data[3] = (_1p0 - eps);
    //std::cout << "\nOk here. 10";
    do
    {
        //Intermediate functions:
        //Directly re-computing fn and hn involves a lot of additions/subtractions. We can speed this process up
        //by storing parts of it and updating them, reducing, for example, the three additions required to compute
        // (ap + M + n + eps) in each step to a single addition, fn1 += _1p0. This will shave off a fair amount of time since
        //we expect to repeat these steps hundreds of times.
        /*fn = (ap + M + n + eps)*(_1p0 - cp + ap + M + n + eps)/((M + n + _1p0 + eps)*(n + _1p0));
        hn = ( (ap + M + n)*(_1p0 - cp + ap + M + n)/(M + n + _1p0) - (ap + M + n)
               - (_1p0 - cp + ap + M + n) - eps +
              (ap + M + n + eps)*(_1p0 - cp + ap + M + n + eps)/(n + _1p0) )/((n + M + _1p0 + eps)*(n + _1p0 - eps));
        */
        fn = fn_data[0]*fn_data[1]/(fn_data[2]*fn_data[3]);
        gn = ( gn_data[0]*gn_data[1]/gn_data[2] - gn_data[3] - gn_data[4] - eps + gn_data[5]*gn_data[6]/gn_data[7] )/(gn_data[8]*gn_data[9]);
        //Update beta_n = b_n (1-z)^n using a cunning trick:
        beta_n *= base*(fn + gn*cn/bn);
        //std::cout << "\nbase = " << base << " fn = " << fn << " gn = " << gn << " cn = " << cn;
        //Add the next term to the Taylor series:
        B += beta_n;
        //Update bn and cn:
        bn = fn*bn + gn*cn;
        cn *= hn_data[0]*hn_data[1]/(hn_data[2]*hn_data[3]);

        //Advance n for the next step:
        n += _1p0;
        //Advance fn,hn,cn:
        for(int i = 0; i < 4;i++)
        {
            fn_data[i] += _1p0;
            hn_data[i] += _1p0;
        }
        for(int i = 0; i < 10;i++)
        {
            gn_data[i] += _1p0;
        }
        //std::cout << "\nbeta_n = " << beta_n << "\tB = " << B;
    }
    while(abs(beta_n) > abs(epsilon*B) && n < Nmax);


    //std::cout << "\nOk here. 10";
    //Comput result:
    value_type PI = boost::math::constants::pi<value_type>();
    value_type PIeps = PI*eps;
    value_type mz = -x;
    value_type ma = -ap;
    result = sign_of<value_type>(m)*pow(mz,ma)*(A + B)/sinc(eps);
    //F = (-1)^m(-z)^{-a}*(A(z) + B(z))/sinc(eps);
    error = abs(beta_n);
    iterations = n.template convert_to<int>();
    //std::cout << "\nOk here. 11";
}

//Choose an appropriate series to use to evaluate the hypergeometric function.
template<class value_type>
void hyperGeometric2F1(value_type a,value_type b,value_type c,value_type x,value_type& error,int& iterations,value_type& result,int& method)
{
    //Define needed functions:
    using namespace boost::math;

    //Constants:
    const value_type _1p0 = value_type(1.0);
    const value_type _0p0 = value_type(0.0);
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();

    //First determine the radius for the different series:
    value_type R_norm = abs(x);//Expansion about x = 0 pole
    value_type R_zozm1_norm = abs(x - _1p0) > epsilon ? abs(x/(x - _1p0)) : _1p0; //Balanced expansion between 0 and 1
    value_type R_zm1_norm = abs(_1p0 - x);//Expansion about x = 1 pole
    value_type R_1oz_norm = abs(x) > epsilon ? abs(_1p0/x) : _1p0; //Expansion about x = \infty pole


    //Threshold to use no-trivial methods like 1-z and 1/z expansion.
    value_type R_thresh = value_type(0.7);
    //Basic idea: if R < R_thresh, use whichever is smaller of |x| or |x/(x-1)|.
    //If greater, then use whichever is smaller of |x|,|x/(x-1)|, |1-x| or |1/x|. We do this because
    //the |1-x| and |1/x| methods are generally slower (more operations required per step) despite their sometimes better convergence properties.
    //Thus is it sometimes more economical to use a slower converging, but faster evaluated method.



    //Special cases:
    //x = 1
    if(abs(x - _1p0) < epsilon)
    {
        value_type cma = c - a;
        value_type cmamb = cma - b;
        value_type cmb = c - b;
        if(cmamb > _0p0)
        {
            //result = tgamma_ratio(c,cma)*tgamma_ratio(cmamb,cmb);
            result = tgamma(c)*tgamma(cmamb)/(tgamma(cma)*tgamma(cmb));
            error = _0p0;
            iterations = 0;
            method = 0;
            //std::cout << "Special value at z = 1.\n";
        }
        else if(cmamb < _0p0)
        {
            //Transform first, then use the above:
            value_type cp = c;
            value_type ap = c - a;
            value_type bp = c - b;
            cma = cp - ap;
            cmamb = cma - bp;
            cmb = cp - bp;
            value_type base = _1p0 - x;
            value_type exponent = c - a - b;
            result = pow(base,exponent)*tgamma(c)*tgamma(cmamb)/(tgamma(cma)*tgamma(cmb));
            error = _0p0;
            iterations = 0;
            method = 0;
            //std::cout << "Special value at z = 1.\n";
        }
        else
        {
            //Special case for c = a + b, x = 1
            result = tgamma(c)/(tgamma(a)*tgamma(b));
            error = _0p0;
            iterations = 0;
            method = 0;
            //std::cout << "Special value at z = 1.\n";
        }
    }
    //x = 0:
    else if(abs(x) < epsilon)
    {
        result = _1p0;
        error = _0p0;
        iterations = 0;
        method = 0;
        //std::cout << "Special value at z = 0.\n";
    }
    else
    {
        if(R_norm < R_thresh || R_zozm1_norm < R_thresh)
        {
            if(R_norm < R_zozm1_norm)
            {
                //std::cout << "Using the basic series.\n";
                hyperGeometric2F1_basic_series(a,b,c,x,error,iterations,result);
                method = 1;

            }
            else
            {
                //std::cout << "Using the z/(z-1) series.\n";
                hyperGeometric2F1_zo_zm1_series(a,b,c,x,error,iterations,result);
                method = 2;

            }
        }
        else
        {
            if(R_norm < R_zozm1_norm)
            {
                if(R_norm < R_zm1_norm)
                {
                    if(R_norm < R_1oz_norm)
                    {
                        //Use the z version
                        //std::cout << "Using the basic series.\n";
                        hyperGeometric2F1_basic_series(a,b,c,x,error,iterations,result);
                        method = 1;

                    }
                    else
                    {
                        //Use the 1/z version:
                        //std::cout << "Using the 1/z series.\n";
                        hyperGeometric2F1_1oz_series(a,b,c,x,error,iterations,result);
                        method = 4;

                    }
                }
                else
                {
                    if(R_zm1_norm < R_1oz_norm)
                    {
                        //Use the z-1 version.
                        //std::cout << "Using the 1-z series.\n";
                        hyperGeometric2F1_1mz_series(a,b,c,x,error,iterations,result);
                        method = 3;

                    }
                    else
                    {
                        //Use the 1/z version:
                        //std::cout << "Using the 1/z series.\n";
                        hyperGeometric2F1_1oz_series(a,b,c,x,error,iterations,result);
                        method = 4;

                    }
                }
            }
            else
            {
                if(R_zozm1_norm < R_zm1_norm)
                {
                    if(R_zozm1_norm < R_1oz_norm)
                    {
                        //Use the z/(z-1) version
                        //std::cout << "Using the z/(z-1) series.\n";
                        hyperGeometric2F1_zo_zm1_series(a,b,c,x,error,iterations,result);
                        method = 2;

                    }
                    else
                    {
                        //Use the 1/z version:
                        //std::cout << "Using the 1/z series.\n";
                        hyperGeometric2F1_1oz_series(a,b,c,x,error,iterations,result);
                        method = 4;

                    }
                }
                else
                {
                    if(R_zm1_norm < R_1oz_norm)
                    {
                        //Use the z-1 version.
                        //std::cout << "Using the 1-z series.\n";
                        hyperGeometric2F1_1mz_series(a,b,c,x,error,iterations,result);
                        method = 3;

                    }
                    else
                    {
                        //Use the 1/z version:
                        //std::cout << "Using the 1/z series.\n";
                        hyperGeometric2F1_1oz_series(a,b,c,x,error,iterations,result);
                        method = 4;

                    }
                }
            }
        }
    }
}
//Direct result version:
template<class value_type>
value_type hyperGeometric2F1(value_type a,value_type b,value_type c,value_type x)
{
    value_type result;
    value_type error;
    int iterations;
    int method;
    hyperGeometric2F1(a,b,c,x,error,iterations,result,method);
    return result;
}


//Returns the Gegenbauer function, C_{alpha}^(lambda)(z), solution to the Gegenbauer differential equation.
template<class value_type>
value_type GegenbauerC(value_type alpha,value_type lambda,value_type z)
{
    //For using boost math functions:
    using namespace boost::math;
    //Needed constants:
    const value_type _0p5 = value_type(0.5);
    const value_type _0p0 = value_type(0.0);
    const value_type _2p0 = value_type(2.0);
    const value_type _1p0 = value_type(1.0);
    const value_type PI = boost::math::constants::pi<value_type>();
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();

    //Check for polynomial case:
    if(alpha - round(alpha) < epsilon)
    {
        //Polynomial!
        value_type r_alpha = round(alpha);
        int n = r_alpha.template convert_to<int>();
        if(n == 0)
        {
            return _1p0;
        }
        else if(n == 1)
        {
            return _2p0*lambda*z;
        }
        else
        {
            value_type Cn = _2p0*lambda*z;
            value_type Cnm1 = _1p0;
            value_type Cnp1;
            value_type I = _1p0;
            value_type factor1 = _2p0*lambda;
            value_type factor2 = factor1 - _1p0;
            for(int i = 1;i < n;i++)
            {
                //Update the factors:
                factor1 += _2p0;
                factor2 += _1p0;
                I += _1p0;
                //Compute the i+1 polynomial:
                Cnp1 = (factor1*z*Cn - factor2*Cnm1)/I;
                //Roll forward the polynomials for the next round:
                Cnm1 = Cn;
                Cn = Cnp1;
            }
            return Cnp1;
        }
    }
    // -lambda - 1/2 is not a positive integer:
    if(!(-lambda - _0p5 > _0p0 && abs(-lambda - _0p5 - round(-lambda - _0p5) < epsilon)))
    {
        return tgamma(alpha + _2p0*lambda)*inverse_gamma(_2p0*lambda)*inverse_gamma(alpha + _1p0)*hyperGeometric2F1(-alpha,_2p0*lambda+alpha,lambda + _0p5,_0p5*(_1p0 - z));
    }
    else
    {
        throw "Cannot represent this Gegenbauer function in terms of the hypergeometric function: (lambda + 1/2) cannot be a non-positive integer for this implementation.";
    }
}


//Returns a Legendre function, P_{\alpha}^{\lambda}(z)
template<class value_type>
value_type LegendreP(value_type alpha,value_type lambda,value_type z)
{
    using namespace boost::math;
    value_type _0p0 = value_type(0.0);
    value_type r_alpha = round(alpha);
    value_type r_lambda = round(lambda);
    if(abs(alpha - r_alpha) < _0p0 && abs(lambda - r_lambda) < _0p0)
    {
        //Use the existing boost Legendre polynomials:
        int l = r_alpha.template convert_to<int>();
        int m = r_lambda.template convert_to<int>();
        return legendre_p(l,m,z);
    }
    else
    {
        //Define these in terms of the hypergeometric function:
        value_type _1p0 = value_type(1.0);
        value_type _2p0 = value_type(2.0);
        value_type base = (_1p0 + z)/(_1p0 - z);
        value_type exponent = lambda/_2p0;
        return inverse_gamma(_1p0 - lambda)*pow(base,exponent)*hyperGeometric2F1(-alpha,alpha + _1p0,_1p0 - lambda,(_1p0 - z)/_2p0);
    }
}



}//hypergeometric
