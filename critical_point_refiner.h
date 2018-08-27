#ifndef CRITICAL_POINT_REFINER_H_INCLUDED
#define CRITICAL_POINT_REFINER_H_INCLUDED

#include "potentials.h"
#include "find_zero.h"

//Function to find 1-D zeros by the Secant method:
template< class value_type , class function >
value_type find_zero(function& f,const value_type& x0)
{
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    const value_type _10p0 = value_type(10.0);
    value_type diff;
    value_type xn = x0;
    value_type xnm1 = x0 + _10p0*eps;
    value_type fn;
    value_type fnm1;
    value_type xnold;
    int counter = 0;
    do
    {
        fn = f(xn);
        fnm1 = f(xnm1);
        xnold = xn;
        xn = xn - f(xn)*(xn - xnm1)/(fn - fnm1);
        ++counter;
    }
    while(abs(xn - xnold) > eps && counter < 1000);
    return xn;
}
//Secant, as above, but with bounds checking using bisection (i.e, false
//position method):
template< class value_type , class function >
value_type find_zero(function& f,const value_type& x0,const value_type& x1)
{
    //Use the secant method, unless this goes outside the bounds,
    //in which case, bisect. Guaranteed to converge.
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    const value_type _2p0 = value_type(2.0);
    const value_type _0p0 = value_type(0.0);

    value_type x_low = x0;
    value_type x_high = x1;
    value_type f0 = f(x0);
    value_type f1 = f(x1);
    std::cout << "\nf0 = " << f0;
    std::cout << "\nf1 = " << f1;
    value_type fn;
    value_type x_n1;
    value_type x_n0;
    value_type x_new;
    if(f0*f1 > _0p0)
    {
        throw "Bounds must bracket a root.";
    }
    int counter = 0;
    do
    {
        x_new = (x_high + x_low)/_2p0;
        fn = f(x_new);
        if(fn*f0 > _0p0)
        {
            x_low = x_new;
            f0 = fn;
        }
        else
        {
            x_high = x_new;
            f1 = fn;
        }
    }
    while(abs(fn) > eps && counter < 1000);
    return x_new;
}
template<class value_type>
class Vpxi
{
private:
    const value_type _1p0;// = value_type(1.0);
    const value_type _6p0;// = value_type(6.0);
    const value_type _4p0;// = value_type(4.0);
    const value_type h;
    const value_type xi;
    const value_type V0;
    potential<value_type>& V;
public:
    Vpxi(value_type H,potential<value_type>& W,value_type XI,value_type W0)
    : h(H),V(W), xi(XI), V0(W0), _1p0(1.0), _6p0(6.0), _4p0(4.0) {}
    value_type operator()(const value_type& x)
    {
        return V.d(x)*(this->_1p0 - this->xi*x*x*this->h*this->h)
            + this->_4p0*xi*this->h*this->h*x*(this->V(x) + this->V0);
    }
};

//Function to refine the locations of the true vacuum, barrier,
//and false vacuum:
template<class value_type>
void refine_critical_points(value_type& false_vacuum,value_type& barrier,
                            value_type& true_vacuum,potential<value_type>& V,
                            value_type xi,value_type h,value_type W0)
{
    const value_type _2p0 = value_type(2.0);
    const value_type _0p0 = value_type(0.0);
    const value_type _1p0 = value_type(1.0);
    const value_type _6p0 = value_type(6.0);
    const value_type _4p0 = value_type(4.0);
    const value_type _8p0 = value_type(8.0);
    const value_type _3p0 = value_type(3.0);
    const value_type _0p01 = value_type(0.01);
    const value_type _24p0 = value_type(24.0);
    const value_type _10p0 = value_type(10.0);
    const value_type PI = boost::math::constants::pi<value_type>();
    const value_type eps = std::numeric_limits<value_type>::epsilon();

    //std::cout << "\n1";
    //LHS_lower = false_vacuum;
    value_type bar_low = (barrier + false_vacuum)/_2p0;
    value_type bar_high = (barrier + true_vacuum)/_2p0;
    //std::cout << "\nbar_low = " << bar_low;
    //std::cout << "\nbar_high = " << bar_high;
    Vpxi<value_type> Vpxi_fun(h,V,xi,W0);
    //std::cout << "\n2";
    //std::cout << "\n " << V.d(bar_low);
    //std::cout << "\n " << V.d(bar_high);
    barrier = find_zero_brent(Vpxi_fun,bar_low,bar_high);
    value_type tv = true_vacuum;
    //LHS_upper = barrier;
    //RHS_lower = barrier;
    //std::cout << "\n3";
    if(xi*(_1p0 - _6p0*xi) > _0p0)
    {
        //Upper bound for phi if xi > 0, where the Einstein-Frame potential
        //diverges:
        value_type upperEnd = _1p0/(h*sqrt(xi*(_1p0 - _6p0*xi)));
        //outStream.precision(50);
        //std::cout << "\nupperEnd = " << upperEnd;
        //std::cout << "\ntv = " << tv;
        if(V.hasMaximum)
        {
            //We don't want to accidentally overflow the potential
            //if it is actually defined on a limited range: if we used the
            //above without checking for this, small xi would lead to the
            //code crashing if V is not defined for large phi.
            upperEnd = V.Maximum < upperEnd ? V.Maximum : upperEnd;
        }
        //Only refine the location if we are sure the zero actually
        //exists. We need to check whether there is a sign change around
        //the given value of tv:
        value_type range = abs(eps);
        bool zero = false;
        if(upperEnd > tv)
        {
            while(range < _1p0)
            {
                //Check for a change of sign:
                zero = Vpxi_fun(tv - range*(tv - barrier))
                        *Vpxi_fun(tv + range*(upperEnd - tv))
                            < _0p0;
                if(zero)
                {
                    //Find the zero in this range and exit:
                    //std::cout << "\nFound zero.";
                    true_vacuum = find_zero_brent(Vpxi_fun,
                                                  tv - range*(tv - barrier),
                                                  tv + range*(upperEnd - tv));
                    break;
                }
                //Else expand the range to keep searching:
                range *= _10p0;
            }
            if(!zero)
            {
                //We couldn't find a change of sign in this range,
                //so we assume there is no true vacuum. Use the upper bound
                //instead:
                //std::cout << "\nNo zero. Use upper bound.";
                true_vacuum = upperEnd - _10p0*eps;
            }
        }
        else
        {
            //Just use the upperEnd, since this appears to have subsumed the
            //true vacuum.
            true_vacuum = upperEnd - _10p0*eps;
        }
    }
    else
    {
        //xi(1-6xi) <= 0, so there is no divergence.
        //However, Einstein frame potential
        // will tend to zero, implying the existence
        //of a true vacuum somewhere. We will assume that it is somewhere
        //between the barrier and 2*tv, which is the user supplied guess.
        //If this is not the case, we will generate an error and ask the user
        //for better data.
        //Also generate an error if the user supplied tv does not describe a
        //negative V, as no bounce solutions exist for this.
        if(V(tv) > _0p0)
        {
            throw "Error: true vacuum must be lower energy than false vacuum.";
        }
        value_type range = abs(eps);
        bool zero = false;
        //std::cout << "\nV.Maximum = " << V.Maximum;
        //std::cout << "\ntv + range = " << tv + range;
        //std::cout << "\nV.hasMaximum = " << V.hasMaximum;
        if(!(V.hasMaximum && tv + range >= V.Maximum))
        {
            while(range < _1p0)
            {
                //std::cout << "\ntest 1";
                //Check for a change of sign:
                zero = Vpxi_fun(tv - range*(tv - barrier))
                        *Vpxi_fun(tv + range*tv) < _0p0;
                /*std::cout << "\nVpxi(lower) = "
                            << Vpxi_fun(tv -range*(tv - barrier));
                std::cout << "\nVpxi(upper) = " << Vpxi_fun(tv + range*tv);*/
                if(zero)
                {
                    //Find the zero in this range and exit:
                    true_vacuum = find_zero_brent(Vpxi_fun,
                                                  tv - range*(tv - barrier),
                                                  tv + range*tv);
                    break;
                }
                //Else expand the range to keep searching:
                range *= _10p0;
                /*
                std::cout.precision(50);
                std::cout << "\nV.Maximum = " << V.Maximum;
                std::cout << "\ntv + range = " << tv + range;
                */
                bool test = V.hasMaximum && tv + range >= V.Maximum;
                //std::cout << "\ntest = " << test;
                if(test)
                {
                    //std::cout << "Should break here!!!";
                    //We don't want to accidentally overflow the potential.
                    break;
                }
            }
        }
        if(!zero)
        {
            //No true vacuum available. Use the original guess:
            //throw "Error: no true vacuum found.";
            //true_vacuum = V.hasMaximum ? V.Maximum : tv;
            true_vacuum = tv;
        }
    }
}

#endif // CRITICAL_POINT_REFINER_H_INCLUDED
