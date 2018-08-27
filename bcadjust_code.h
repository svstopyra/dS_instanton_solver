#ifndef BCADJUST_CODE
#define BCADJUST_CODE
#include "bcadjust.h"

//------------------------------------------------------------------------------
//Finds the maximum of a set of values in an array:
template<class value_type>
value_type max_of_array(value_type* val_array,int n)
{
    if(n < 1)
    {
        throw "Invalid array for finding maximum: no elements.";
    }
    value_type maximum = val_array[0];
    for(int i = 0;i < n;i++)
    {
        if(val_array[i] > maximum)
        {
            maximum = val_array[i];
        }
    }
    return maximum;
}
//------------------------------------------------------------------------------
//Finds the minimum of a set of values in an array:
template<class value_type>
value_type min_of_array(value_type* val_array,int n)
{
    if(n < 1)
    {
        throw "Invalid array for finding minimum: no elements.";
    }
    value_type minimum = val_array[0];
    for(int i = 0;i < n;i++)
    {
        if(val_array[i] < minimum)
        {
            minimum = val_array[i];
        }
    }
    return minimum;
}
//------------------------------------------------------------------------------
//Computes x_tilde from the data in an x_tilde_data struct:
template<class value_type>
value_type x_tilde(value_type x,x_tilde_data<value_type>& d)
{
    value_type factor = sqrt(d._1p0 - d.xi*d.h*d.h*d.y0*d.y0);
    value_type conformal_factor = sqrt(d._1p0 - d.xi*(d._1p0 - d._6p0*d.xi)
                                       *d.h*d.h*d.y0*d.y0);
    value_type cosh_factor = d.hyperbolic ? cosh(d.H_tilde*factor*x)
                            : cos(d.H_tilde*factor*x);

    return factor*x - ((d.xi*d.h*d.h*d.y0)/(conformal_factor))*(d.A*factor*x
            + (d.B/d.H_tilde)*(Gi(cosh_factor,d.alpha)));
}
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
time_type dynamicBCs_basic< value_type ,solution_type, time_type >::map_BCs
    (solution_type& BCOld,solution_type& BCNew,
     instanton_equation_lineariser<value_type>& d)
{
    //Extract real initial conditions:
    value_type y0 = BCOld[0];
    value_type yp0 = BCOld[1];
    value_type a0 = BCOld[2];
    value_type ap0 = BCOld[3];



    //Quantities to compare delta_phi against:
    value_type y0_bound = abs(y0);
    //std::cout << "\nd.Wppy0 = " << d.Wppy0;
    value_type V0_bound = abs(d.Wy0/d.Wppy0);
    value_type V1_bound = abs(d.Wpy0/d.Wppy0);
    value_type phi_dot_bound = abs(d.Wy0);
    //std::cout << "\ny0_bound = " << y0_bound;
    //std::cout << "\nV0_bound = " << V0_bound;
    //std::cout << "\nV1_bound = " << V1_bound;
    //std::cout << "\nphi_dot_bound = " << phi_dot_bound;

    //We want A + BC_{\alpha}^{3/2}( cosh(x) ) < stepError

    //Determine threshold such that all of the required conditions are
    //satisfied:
    value_type threshold = std::min(y0_bound,std::min(V0_bound,V1_bound));
    std::cout << "\nthreshold = " << threshold;
    //Guess the initial value of chi to use. Try the expansion phi = phi_0 +
    //V'(phi_0)x^2/8 first, solving for phi - phi_0 = threshold.
    //If this doesn't work, because V' is too small (or zero) try
    //sqrt(threshold).
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    //value_type x_guess = abs(d.Wpy0) > eps ?
    //            sqrt(abs(d._8p0*(y0 - threshold)/(d.Wpy0))) : sqrt(threshold);
    value_type x_guess = sqrt(threshold);
    std::cout << "\nx_guess = " << x_guess;

    //Compute delta_y at this point:
    value_type delta_y = d.delta_y(x_guess);
    std::cout << "\ndelta_y = " << delta_y;
    value_type delta_yp = d.delta_yp(x_guess);
    std::cout << "\ndelta_yp = " << delta_yp;
    //Loop until we have delta_y smaller than the desired threshold:
    int counter = 0;
    while((abs(delta_y) > threshold || delta_yp*delta_yp > phi_dot_bound )
          && counter < 1000)
    {
        //Solve using the Newton-Raphson solver until within the threshold:
        ++counter;
        x_guess = x_guess - delta_y/delta_yp;
        //Re-compute functions:
        delta_y = d.delta_y(x_guess);
        delta_yp = d.delta_yp(x_guess);
    }
    if(counter > 1000)
    {
        throw "Newton-Raphson loop took too long.";
    }

    //Now compute a at this point:
    value_type a = d.a(x_guess);
    std::cout << "\na = " << a;
    value_type ap = d.ap(x_guess);
    std::cout << "\nap = " << ap;

    //Return the new boundary conditions. Remember that the above are all in
    //the Einstein
    //frame, NOT the Jordan frame, which we need to switch back to to solve
    //the ode:
    BCNew[0] = y0 + d.factor2*delta_y/d.conformal_factor;
    std::cout << "\ny0_new = " << BCNew[0];
    BCNew[1] = (d.factor2*delta_yp/d.conformal_factor)*d.factor;
        //extra factor since we need to convert from x_tilde to x units in the
        //derivative.
    std::cout << "\ny0p_new = " << BCNew[1];
    BCNew[2] = a*(d._1p0 + this->xi*d.h2*y0*BCNew[0]/d.factor2)/d.factor;
    std::cout << "\na0_new = " << BCNew[2];
    BCNew[3] = ap*(d._1p0 + this->xi*d.h2*y0*BCNew[0]/d.factor2)
                + a*( this->xi*d.h2*y0*delta_yp/d.conformal_factor );
    std::cout << "\na0p_new = " << BCNew[3];

    //Return the value of x, as a function of x_tilde:
    d.x_tilde = x_guess;//Store x_tilde value for later use.
    value_type cosh_factor = d.hyperbolic ? cosh(d.H_tilde*x_guess)
                            : cos(d.H_tilde*x_guess);
    std::cout << "\ncosh_factor = " << cosh_factor;
    std::cout << "\nx_guess = " << x_guess;
    std::cout << "\nd.hyperbolic = " << d.hyperbolic;
    return (x_guess + this->xi*d.h2*y0*(d.A*x_guess
            + d.B*Gi(cosh_factor,d.alpha)/d.H_tilde )/d.conformal_factor )
            /d.factor;
}
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
time_type dynamicBCs_basic< value_type ,solution_type, time_type >::map_BCs
    (solution_type& BCOld,solution_type& BCNew)
{
    //Setup storage class. This also computes several intermediate variables
    //that we will access for use below, and
    //defines the functions needed to compute the solutions of the linearised
    //equations.
    value_type y0 = BCOld[0];
    std::cout << "\ny0 = " << y0;
    instanton_equation_lineariser<value_type> d
        (y0,this->h,this->xi,this->W,this->W0);
    return this->map_BCs(BCOld,BCNew,d);
}
//------------------------------------------------------------------------------
//Generic case:
template< class value_type , class solution_type , class time_type>
time_type dynamicBCs_basicMethod< value_type ,solution_type, time_type >
    ::operator()(solution_type& BCOld,solution_type& BCNew)
{
    return this->map_BCs(BCOld,BCNew);
}
//------------------------------------------------------------------------------
//Member functions for the equation with action:
template< class value_type , class solution_type , class time_type>
time_type dynamicBCs_actionMethod< value_type ,solution_type, time_type >
    ::operator()(solution_type& BCOld,solution_type& BCNew)
{
    //Process the first 4 components as normal:
    value_type y0 = BCOld[0];
    instanton_equation_lineariser<value_type> d(y0,this->h,this->xi,
                                                this->W,this->W0);//Computes
                                                //most of the intermediates
    //as part of the constructor, so we don't have to repeat this word in
    //map_BCs and S.
    time_type x = this->map_BCs(BCOld,BCNew,d);//result is in the Jordan frame.

    //Now evaluate the action up to this point. Remember that we evaluate it
    //in the Einstein
    //frame, so use the co-ordinate x_tilde (computed as an intermediate by
    //map_BCs and
    //stored in d), rather than x:
    BCNew[4] = S(d.x_tilde,d);
    //Finally, return the time in the Jordan frame:
    return x;
}
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
value_type dynamicBCs_actionMethod< value_type ,solution_type, time_type >::S
    (value_type& x_tilde,instanton_equation_lineariser<value_type>& d)
{
    //cosh factor, as S = S(cosh(H*x)).
    value_type cosh_factor = d.hyperbolic ? cosh(d.H_tilde*x_tilde)
                            : cos(d.H_tilde*x_tilde);

    //Useful intermediates:
    value_type H_tilde2 = d.H_tilde*d.H_tilde;
    value_type H_tilde4 = H_tilde2*H_tilde2;
    const value_type PI = boost::math::constants::pi<value_type>();
    value_type _2PI2 = d._2p0*PI*PI;

    //Two constants appearing in the action formula:
    value_type A_tilde = -_2PI2*(d.Wy0 + d.Wpy0*d.A)/H_tilde4;
    value_type B_tilde = -_2PI2*d.Wpy0*d.B/H_tilde4;

    //Legendre functions needed to express the action;
    value_type P_n = hypergeometric::LegendreP(d.alpha,d._0p0,cosh_factor);
    value_type P_n_p1 = hypergeometric::LegendreP(d.alpha + d._1p0,d._0p0,
                                                  cosh_factor);
        //Use the recursion relation to fix the other two Legendre functions,
        //without having to
        //actually compute them, which is expensive:
    value_type P_n_p2 = ((d._2p0*d.alpha + d._3p0)*cosh_factor*P_n_p1
                         - (d.alpha + d._1p0)*P_n)/(d.alpha + d._2p0);
    value_type P_n_p3 = ((d._2p0*d.alpha + d._5p0)*cosh_factor*P_n_p2
                         - (d.alpha + d._2p0)*P_n_p1)/(d.alpha + d._3p0);

    //Compute S using the integratiion formula:
    return A_tilde*( cosh_factor*cosh_factor/d._3p0 - d._1p0 )*cosh_factor
     - B_tilde*(d.alpha + d._1p0)*( -( P_n_p3 - cosh_factor*P_n_p2 )
                                   /(d._2p0*d.alpha + d._3p0) +
        (d._1p0 - (d.alpha + d._1p0)/(d._2p0*d.alpha + d._3p0) )
                                   *( P_n_p1 - cosh_factor*P_n )/d.alpha );
}
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
value_type dynamicBCs_actionMethod< value_type ,solution_type, time_type >::S
    (value_type& x_tilde,value_type& y0)
{
    instanton_equation_lineariser<value_type> d
        (y0,this->h,this->xi,this->W,this->W0);
    return S(x_tilde,d);
}
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
time_type dynamicBCs_flatSpace< value_type, solution_type, time_type>::map_BCs
    (solution_type& BCOld,solution_type& BCNew,
     instanton_equation_lineariser<value_type>& d)
{
    //Fetch initial conditions:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];

    value_type chi_step;

    const value_type _0p0 = value_type(0.0);
    const value_type _8p0 = value_type(8.0);
    const value_type _10p0 = value_type(10.0);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    //instanton_equation_lineariser<value_type> d(y0,_0p0,_0p0,this->W,_0p0);

    //Initial guess at chi_step:
    chi_step = d.Wpy0_J == _0p0 ? eps : sqrt(_8p0*_10p0*eps/abs(d.Wpy0_J));
    if(d.Wppy0_J != _0p0)
    {
        //Need to use Newton-Raphson to improve the guess
        chi_step = NewtonRaphsonSimple(chi_step,&(d.delta_y_flat),
                                       &(d.delta_yp_flat));
    }
    //Return the function at these points:
    BCNew[0] = d.delta_y_flat(chi_step);
    BCNew[1] = d.delta_yp_flat(chi_step);
    return chi_step;
}
//------------------------------------------------------------------------------
template< class value_type , class solution_type , class time_type>
time_type dynamicBCs_flatSpace< value_type, solution_type, time_type>
    ::map_BCs(solution_type& BCOld,solution_type& BCNew)
{
    //Fetch initial conditions:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    const value_type _0p0 = value_type(0.0);
    instanton_equation_lineariser<value_type> d(y0,_0p0,_0p0,this->W,_0p0);
    return this->map_BCs(BCOld,BCNew,d);
}
//------------------------------------------------------------------------------
//Flat space decaying to AdS space:
#ifdef ANTI_DE_SITTER_FLAT_INSTANTONS
//Flat false vacuum to AdS true vacuum:
template< class value_type, class solution_type, class time_type >
time_type dynamicBCsAdSFlat_taylor< value_type, solution_type, time_type >
    ::basic_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    value_type a0 = BCOld[2];
    value_type a0p = BCOld[3];
    //std::cout << "\ny0 = " << y0;

    //Intermediates:
    this->Wy0 = this->W(y0);
    this->Wy0p = this->W.d(y0);
    std::cout << "\ny0p = " << y0p;
    std::cout << "\nWy0p = " << this->Wy0p;
    value_type y02 = y0*y0;
    value_type h2 = this->h*this->h;
    value_type conformal = (this->_1p0 - this->xi*h2*y02);
    value_type conf_factor = (this->_1p0 - this->xi*(this->_1p0
                - this->_6p0*this->xi)*h2*y02);
    this->y0pp = ( this->_0p25*this->Wy0p*conformal +h2*this->xi*y0*(this->Wy0
                + this->W0)  )/conf_factor;
    //this->a0pp = - h2*(this->W0 + this->Wy0
    //             - this->_6p0*this->xi*y0*this->y0pp)/(this->_3p0*conformal);
    //Third derivative of a at 0:
    this->a0pp = -a0p*h2*(this->Wy0 + this->W0 -
                          this->_3p0*this->xi*y0*this->Wy0p/this->_2p0)
                          /(this->_3p0*conf_factor);

    //Make phi - phi_0 as small as possible without underflowing:
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();
    const value_type _0p0 = value_type(0.0);
    const value_type x = this->Wy0p == this->_0p0 ? this->stepError/this->_2p0 :
        std::min(sqrt(abs(this->stepError*y0/this->y0pp)),
                 this->stepError/this->_2p0);
    std::cout << "\nx = " << x;

    //Return new BCs:
    BCNew[0] = y0 + this->_0p5*this->y0pp*x*x;
    BCNew[1] = y0p + this->y0pp*x;
    BCNew[2] = x*(a0p + this->a0pp*x*x/this->_6p0);
    BCNew[3] = a0p + this->_0p5*this->a0pp*x*x;
    std::cout << "\ny = " << BCNew[0];
    std::cout << "\nyp = " << BCNew[1];
    std::cout << "\na = " << BCNew[2];
    std::cout << "\nap = " << BCNew[3];
    return x;
}
//------------------------------------------------------------------------------
template< class value_type, class solution_type, class time_type >
time_type dynamicBCsAdSFlat_taylor< value_type, solution_type, time_type >
::action_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //First compute the ordinary boundary conditions:
    time_type x = basic_bcs(BCOld,BCNew);

    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    value_type a0 = BCOld[2];
    value_type a0p = BCOld[3];


    //Now add the action:

    BCNew[4] = this->_0p5*this->PI*this->PI*(this->_1p0 + x*this->a0pp)
               *(-(this->Wy0 + this->W0) + this->_3p0*this->xi*y0*this->Wy0p
               + this->_3p0*this->xi*this->xi*this->h*this->h*y0*y0*(this->_4p0
               *(this->Wy0 + this->W0))/(this->_1p0
                - this->h*this->h*y0*y0*this->xi
               *(this->_1p0 - this->_6p0*this->xi)) );
    return x;
}
//------------------------------------------------------------------------------
//Flat false vacuum to AdS true vacuum:
template< class value_type, class solution_type, class time_type >
time_type dynamicBCs_flatfv_delta< value_type, solution_type, time_type >
    ::basic_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    value_type a0 = BCOld[2];
    value_type a0p = BCOld[3];
    value_type da0 = BCOld[4];
    value_type da0p = BCOld[5];
    //std::cout << "\ny0 = " << y0;

    //Intermediates:
    this->Wy0 = this->W(y0);
    this->Wy0p = this->W.d(y0);
    std::cout << "\ny0p = " << y0p;
    std::cout << "\nWy0p = " << this->Wy0p;
    value_type y02 = y0*y0;
    value_type h2 = this->h*this->h;
    value_type conformal = (this->_1p0 - this->xi*h2*y02);
    value_type conf_factor = (this->_1p0 - this->xi*(this->_1p0
                - this->_6p0*this->xi)*h2*y02);
    this->y0pp = ( this->_0p25*this->Wy0p*conformal +h2*this->xi*y0*(this->Wy0
                + this->W0)  )/conf_factor;
    //this->a0pp = - h2*(this->W0 + this->Wy0
    //             - this->_6p0*this->xi*y0*this->y0pp)/(this->_3p0*conformal);
    //Third derivative of a at 0:
    this->a0pp = -a0p*h2*(this->Wy0 + this->W0 -
                          this->_3p0*this->xi*y0*this->Wy0p/this->_2p0)
                          /(this->_3p0*conf_factor);

    //Make phi - phi_0 as small as possible without underflowing:
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();
    const value_type _0p0 = value_type(0.0);
    const value_type x = this->Wy0p == this->_0p0 ? this->stepError/this->_2p0 :
        std::min(sqrt(abs(this->stepError*y0/this->y0pp)),
                 this->stepError/this->_2p0);
    std::cout << "\nx = " << x;

    //Return new BCs:
    BCNew[0] = y0 + this->_0p5*this->y0pp*x*x;
    BCNew[1] = y0p + this->y0pp*x;
    BCNew[2] = x*(a0p + this->a0pp*x*x/this->_6p0);
    BCNew[3] = a0p + this->_0p5*this->a0pp*x*x;
    BCNew[4] = BCNew[2] - x;
    BCNew[5] = BCNew[3] - this->_1p0;
    std::cout << "\ny = " << BCNew[0];
    std::cout << "\nyp = " << BCNew[1];
    std::cout << "\na = " << BCNew[2];
    std::cout << "\nap = " << BCNew[3];
    return x;
}
//------------------------------------------------------------------------------
template< class value_type, class solution_type, class time_type >
time_type dynamicBCs_flatfv_delta< value_type, solution_type, time_type >
::action_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //First compute the ordinary boundary conditions:
    time_type x = basic_bcs(BCOld,BCNew);

    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    value_type a0 = BCOld[2];
    value_type a0p = BCOld[3];


    //Now add the action:

    BCNew[6] = this->_0p5*this->PI*this->PI*(this->_1p0 + x*this->a0pp)
               *(-(this->Wy0 + this->W0) + this->_3p0*this->xi*y0*this->Wy0p
               + this->_3p0*this->xi*this->xi*this->h*this->h*y0*y0*(this->_4p0
               *(this->Wy0 + this->W0))/(this->_1p0
                - this->h*this->h*y0*y0*this->xi
               *(this->_1p0 - this->_6p0*this->xi)) );
    return x;
}
//------------------------------------------------------------------------------
template< class value_type, class solution_type, class time_type >
time_type dynamicBCsAdSFlat_taylor_delta_a< value_type, solution_type,time_type>
    ::basic_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    value_type da0 = BCOld[2];
    value_type da0p = BCOld[3];
    //std::cout << "\ny0 = " << y0;

    //Intermediates:
    this->Wy0 = this->W(y0);
    this->Wy0p = this->W.d(y0);
    const value_type _3p0 = value_type(3.0);
    const value_type _1p0 = value_type(1.0);
    const value_type _0p0 = value_type(0.0);
    value_type H0 = this->h*sqrt(this->W0/_3p0);
    value_type a0 = da0;
    value_type a0p = da0p + _1p0;
    std::cout << "\ny0p = " << y0p;
    std::cout << "\nWy0p = " << this->Wy0p;
    value_type y02 = y0*y0;
    value_type h2 = this->h*this->h;
    value_type conformal = (this->_1p0 - this->xi*h2*y02);
    value_type conf_factor = (this->_1p0 - this->xi*(this->_1p0
                - this->_6p0*this->xi)*h2*y02);
    this->y0pp = ( this->_0p25*this->Wy0p*conformal +h2*this->xi*y0*(this->Wy0
                + this->W0)  )/conf_factor;
    //this->a0pp = - h2*(this->W0 + this->Wy0
    //             - this->_6p0*this->xi*y0*this->y0pp)/(this->_3p0*conformal);
    //Third derivative of a at 0:
    this->a0pp = -a0p*h2*(this->Wy0 + this->W0 -
                          this->_3p0*this->xi*y0*this->Wy0p/this->_2p0)
                          /(this->_3p0*conf_factor);

    //Make phi - phi_0 as small as possible without underflowing:
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();
    const value_type x = this->Wy0p == _0p0 ? this->stepError/this->_2p0 :
        std::min(sqrt(abs(this->stepError*y0/this->y0pp)),
                 this->stepError/this->_2p0);
    std::cout << "\nx = " << x;

    //Return new BCs:
    value_type arg = H0*x;
    BCNew[0] = y0 + this->_0p5*this->y0pp*x*x;
    BCNew[1] = y0p + this->y0pp*x;
    BCNew[2] = x*(a0p + this->a0pp*x*x/this->_6p0) - sin(arg)/H0;
    BCNew[3] = a0p + this->_0p5*this->a0pp*x*x - cos(arg);
    std::cout << "\ny = " << BCNew[0];
    std::cout << "\nyp = " << BCNew[1];
    std::cout << "\na = " << BCNew[2];
    std::cout << "\nap = " << BCNew[3];
    return x;
}
//------------------------------------------------------------------------------
template< class value_type, class solution_type, class time_type >
time_type dynamicBCsAdSFlat_taylor_delta_a< value_type, solution_type,time_type>
::action_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //First compute the ordinary boundary conditions:
    time_type x = basic_bcs(BCOld,BCNew);

    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    value_type a0 = BCOld[2];
    value_type a0p = BCOld[3];
    //std::cout << "\nbasics_bcs works just fine.";


    //Now add the action:

    /*
    BCNew[4] = this->_0p5*this->PI*this->PI*(this->_1p0 + x*this->a0pp)
               *(-(this->Wy0 + this->W0) + this->_3p0*this->xi*y0*this->Wy0p
               + this->_3p0*this->xi*this->xi*this->h*this->h*y0*y0*(this->_4p0
               *(this->Wy0 + this->W0))/(this->_1p0
                - this->h*this->h*y0*y0*this->xi
               *(this->_1p0 - this->_6p0*this->xi)) );*/
    //Approximations to action after a small step:
    value_type H0 = this->h*sqrt(this->W0/this->_3p0);
    value_type arg = H0*x;
    value_type sinarg = sin(arg);
    value_type a = BCNew[2] + (sinarg/H0);
    /*std::cout << "\nGot here ok. H0 = " << H0 << " W0 = " << W0 << " h = "
              << h
              << " arg = " << arg << " sinarg = " << sinarg << " a = " << a;*/
    BCNew[4] = -this->_2p0*this->PI*this->PI*(a*a*a)
                *this->W(BCNew[0])*x;
    //std::cout << "\nGot BCNew[4] = " << BCNew[4];

    BCNew[5] = -this->_6p0*this->PI*this->PI*this->h*this->h*
                        (this->_3p0*sinarg*sinarg*BCNew[2] +
                        this->_3p0*sinarg*H0*BCNew[2]*BCNew[2] +
                        H0*H0*BCNew[2]*BCNew[2]*BCNew[2])*x;
    //std::cout << "\nGot BCNew[5] = " << BCNew[5];
    return x;
}
//------------------------------------------------------------------------------
template< class value_type, class solution_type, class time_type >
time_type dynamicBCs_dS_fixed< value_type, solution_type,time_type>
    ::basic_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    //std::cout << "\ny0 = " << y0;

    //Intermediates:
    this->Wy0 = this->W(y0);
    this->Wy0p = this->W.d(y0);
    const value_type _3p0 = value_type(3.0);
    const value_type _1p0 = value_type(1.0);
    const value_type _0p0 = value_type(0.0);
    value_type H0 = this->h*sqrt(this->W0/_3p0);
    value_type a0 = _0p0;
    value_type a0p = _1p0;
    std::cout << "\ny0p = " << y0p;
    std::cout << "\nWy0p = " << this->Wy0p;
    value_type y02 = y0*y0;
    value_type h2 = this->h*this->h;
    value_type conformal = (this->_1p0 - this->xi*h2*y02);
    value_type conf_factor = (this->_1p0 - this->xi*(this->_1p0
                - this->_6p0*this->xi)*h2*y02);
    this->y0pp = ( this->_0p25*this->Wy0p*conformal +h2*this->xi*y0*(this->Wy0
                + this->W0)  )/conf_factor;
    //this->a0pp = - h2*(this->W0 + this->Wy0
    //             - this->_6p0*this->xi*y0*this->y0pp)/(this->_3p0*conformal);
    //Third derivative of a at 0:
    this->a0pp = -a0p*h2*(this->Wy0 + this->W0 -
                          this->_3p0*this->xi*y0*this->Wy0p/this->_2p0)
                          /(this->_3p0*conf_factor);

    //Make phi - phi_0 as small as possible without underflowing:
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();
    const value_type x = this->Wy0p == _0p0 ? this->stepError/this->_2p0 :
        std::min(sqrt(abs(this->stepError*y0/this->y0pp)),
                 this->stepError/this->_2p0);
    std::cout << "\nx = " << x;

    //Return new BCs:
    value_type arg = H0*x;
    BCNew[0] = y0 + this->_0p5*this->y0pp*x*x;
    BCNew[1] = y0p + this->y0pp*x;
    std::cout << "\ny = " << BCNew[0];
    std::cout << "\nyp = " << BCNew[1];
    return x;
}
//------------------------------------------------------------------------------
template< class value_type, class solution_type, class time_type >
time_type dynamicBCs_dS_fixed< value_type, solution_type,time_type>
::action_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //First compute the ordinary boundary conditions:
    time_type x = basic_bcs(BCOld,BCNew);

    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    //std::cout << "\nbasics_bcs works just fine.";


    //Now add the action:

    //Approximations to action after a small step:
    value_type H0 = this->h*sqrt(this->W0/this->_3p0);
    value_type arg = H0*x;
    value_type sinarg = sin(arg);
    value_type a = (sinarg/H0);
    /*std::cout << "\nGot here ok. H0 = " << H0 << " W0 = " << W0 << " h = "
              << h
              << " arg = " << arg << " sinarg = " << sinarg << " a = " << a;*/
    BCNew[2] = -this->_2p0*this->PI*this->PI*(a*a*a)
                *this->Wy0*x;
    return x;
}
//------------------------------------------------------------------------------
//Gravitational case tracking friction rather than a and a'
template< class value_type, class solution_type, class time_type >
time_type dynamicBCsGravFriction_taylor< value_type, solution_type, time_type >
    ::basic_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    value_type theta = BCOld[2];
    value_type thetap = BCOld[3];
    //std::cout << "\ny0 = " << y0;

    //Intermediates:
    this->Wy0 = this->W(y0);
    this->Wy0p = this->W.d(y0);
    value_type y02 = y0*y0;
    value_type h2 = this->h*this->h;
    value_type conformal = (this->_1p0 - this->xi*h2*y02);
    this->y0pp = ( this->_0p25*this->Wy0p*conformal +h2*this->xi*y0*(this->Wy0
                + this->W0)  )/(this->_1p0 - this->xi*(this->_1p0
                - this->_6p0*this->xi)*h2*y02);
    this->a0pp = - h2*(this->W0 + this->Wy0
                 - this->_6p0*this->xi*y0*this->y0pp)/(this->_3p0*conformal);

    //Make phi - phi_0 as small as possible without underflowing:
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();
    const value_type _0p0 = value_type(0.0);
    const value_type x = this->Wy0p == this->_0p0 ? this->stepError/this->_2p0 :
        std::min(sqrt(abs(this->stepError*y0/this->y0pp)),
                 this->stepError/this->_2p0);

    //Return new BCs:
    BCNew[0] = y0 + this->_0p5*this->y0pp*x*x;
    BCNew[1] = y0p + this->y0pp*x;
    BCNew[2] = x - x*x/this->_2p0;
    BCNew[3] = this->_1p0 - x;
    return x;
}
//------------------------------------------------------------------------------
template< class value_type, class solution_type, class time_type >
time_type dynamicBCsGravFriction_taylor< value_type, solution_type, time_type >
::action_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //First compute the ordinary boundary conditions:
    time_type x = basic_bcs(BCOld,BCNew);

    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];
    value_type theta = BCOld[2];
    value_type thetap = BCOld[3];

    //Now add the action:
    BCNew[4] = this->_0p5*this->PI*this->PI*(this->_1p0 + x*this->a0pp)
               *(-this->Wy0 + this->_3p0*this->xi*y0*this->Wy0p
               + this->_3p0*this->xi*this->xi*this->h*this->h*y0*y0*(this->_4p0
               *(this->Wy0))/(this->_1p0 - this->h*this->h*y0*y0*this->xi
               *(this->_1p0 - this->_6p0*this->xi)) );
    return x;
}
//------------------------------------------------------------------------------
//Flat false vacuum to flat true vacuum (no gravity)
template< class value_type, class solution_type, class time_type >
time_type dynamicBCsFlatFixed_taylor< value_type, solution_type, time_type >
::basic_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];

    //Intermediates:
    this->Wy0 = this->W(y0);
    this->Wy0p = this->W.d(y0);
    value_type y02 = y0*y0;
    this->y0pp =this->_0p25*this->Wy0p;

    //Make phi - phi_0 as small as possible without underflowing:
    const value_type epsilon = std::numeric_limits<value_type>::epsilon();
    const value_type _0p0 = value_type(0.0);
    const value_type x = this->Wy0p == this->_0p0 ? this->stepError/this->_2p0 :
        std::min(sqrt(abs(this->stepError*y0/this->y0pp)),
                 this->stepError/this->_2p0);

    //Return new BCs:
    BCNew[0] = y0 + this->_0p5*this->y0pp*x*x;
    BCNew[1] = y0p + this->y0pp*x;
    return x;
}
//------------------------------------------------------------------------------
template< class value_type, class solution_type, class time_type >
time_type dynamicBCsFlatFixed_taylor< value_type, solution_type, time_type >
::action_bcs(solution_type& BCOld, solution_type& BCNew)
{
    //First compute the ordinary boundary conditions:
    time_type x = basic_bcs(BCOld,BCNew);

    //Get Old BCs:
    value_type y0 = BCOld[0];
    value_type y0p = BCOld[1];

    //Now add the action:
    BCNew[2] = this->_0p5*this->PI*this->PI*this->Wy0*x*x*x*x;
    return x;
}
//------------------------------------------------------------------------------
#endif // ANTI_DE_SITTER_FLAT_INSTANTONS

//For later implementation:
#ifdef DE_SITTER_INSTANTONS
#endif // DE_SITTER_INSTANTONS

#endif//BCADJUST_CODE
