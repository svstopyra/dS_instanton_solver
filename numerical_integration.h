#ifndef NUMERICAL_INTEGRATION_H
#define NUMERICAL_INTEGRATION_H
#include<vector>
//Various numerical integration techniques that may be useful
//------------------------------------------------------------------------------
//Numerical integration via the trapezium rule, for non-equally spaced points.
template< class value_type >
value_type trap(std::vector<value_type>& x,std::vector<value_type>& y)
{
    if(x.size() != y.size())
    {
        throw "Vectors supplied to trap must be the same size.";
    }
    value_type integral = value_type(0.0);
    const value_type _2p0 = value_type(2.0);
    for(int i = 0;i < int(x.size()) - 1;i++)
    {
        integral += (x[i + 1] - x[i])*(y[i + 1] + y[i])/_2p0;
    }
    return integral;
}
#endif//NUMERICAL_INTEGRATION_H
//------------------------------------------------------------------------------
