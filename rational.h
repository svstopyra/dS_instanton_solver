#ifndef RATIONAL_H
#define RATIONAL_H
//Code which defines rational functions and product functions
// and their derivatives, up to fourth.

//These are abstract classes - it is necessary to specify A and B to use them:


//------------------------------------------------------------------------------
//Rational functions, f(x) = A(x)/B(x).
//Separate arg_type and value_type because we may have a rational
//function for which the numerator and denominator are
//known in terms of other functions, which we may want to pass
//as arguments, rather than simply passing an argument like x, at
//which to evaluate them (for example - the rational function may
//appear on the RHS of an ode, dg/dt = A(t,g(t))/B(t,g(t)), and
//we may know t and g(t), but don't want to recompute g(t) - for
//this reason we may want to pass (t,g(t)) as arguments (or in some case
// just g(t) ) instead of t. Thus, we leave arg_type free to be different to
//value_type.
template<class value_type,class arg_type>
class rational_function
{
private:
    virtual value_type A(arg_type x) = 0;
    virtual value_type B(arg_type x) = 0;
    virtual value_type Ap(arg_type x) = 0;
    virtual value_type Bp(arg_type x) = 0;
    virtual value_type App(arg_type x) = 0;
    virtual value_type Bpp(arg_type x) = 0;
    virtual value_type Appp(arg_type x) = 0;
    virtual value_type Bppp(arg_type x) = 0;
public:
    const value_type _2p0;
    const value_type _6p0;
    const value_type _3p0;
    const value_type _24p0;
    const value_type _12p0;
    const value_type _36p0;
    const value_type _4p0;
    const value_type _8p0;
    //Constructor:
    rational_function() : _2p0(2.0),_6p0(6.0),_3p0(3.0),_24p0(24.0),
        _12p0(12.0),_36p0(36.0),_4p0(4.0),_8p0(8.0) {}
    //Virtual functions:
    //virtual void ComputeSharedData(arg_type x) = 0;
    //virtual value_type Apppp(arg_type x) = 0;
    //virtual value_type Bpppp(arg_type x) = 0;
    value_type f(arg_type x)
    {
        return this->A(x)/this->B(x);
    }
    value_type d(arg_type x)
    {
        value_type Bx = this->B(x);
        return this->Ap(x)/Bx - this->A(x)*this->Bp(x)/(Bx*Bx);
    }
    value_type d2(arg_type x)
    {
        value_type Bx = this->B(x);
        value_type Ax = this->A(x);
        value_type Bx2 = Bx*Bx;
        value_type Bpx = this->Bp(x);
        value_type Bpx2 = Bpx*Bpx;
        return -this->_2p0*this->Ap(x)*Bpx/Bx2 + this->_2p0*Ax*Bpx2/(Bx*Bx2)
            +this->App(x)/Bx - Ax*this->Bpp(x)/Bx2;
    }
    value_type d3(arg_type x)
    {
        value_type Bx = this->B(x);
        value_type Ax = this->A(x);
        value_type Bpx = this->Bp(x);
        value_type Apx = this->Ap(x);
        value_type Bppx = this->Bpp(x);
        value_type Appx = this->App(x);
        value_type Bx2 = Bx*Bx;
        value_type Bpx2 = Bpx*Bpx;
        return this->_6p0*Apx*Bpx2/(Bx2*Bx) - this->_6p0*A*Bpx2*Bpx/(Bx2*Bx2)
                - this->_3p0*Bpx*Appx/Bx2 - this->_3p0*Apx*Bppx/Bx2
                + this->_6p0*Ax*Bpx*Bppx/(Bx2*Bx) + this->Appp(x)/Bx
                - Ax*this->Bppp(x)/Bx2;
    }
    /*
    value_type d4(arg_type x)
    {
        value_type Bx = this->B(x);
        value_type Ax = this->A(x);
        value_type Bpx = this->Bp(x);
        value_type Apx = this->Ap(x);
        value_type Bppx = this->Bpp(x);
        value_type Appx = this->App(x);
        value_type Bpppx = this->Bppp(x);
        value_type Apppx = this->Appp(x);
        value_type Bx2 = Bx*Bx;
        value_type Bx3 = Bx2*Bx;
        value_type Bpx2 = Bpx*Bpx;
        value_type Bpx3 = Bpx*Bpx2;
        return -this->_24p0*Apx*Bpx3/(Bx2*Bx2)
               + this->_24p0*Ax*Bpx2*Bpx2/(Bx3*Bx2)
               + this->_12p0*Bpx2*Appx/Bx3 + this->_24p0*Apx*Bpx*Bppx/Bx3
               - this->_36p0*Ax*Bpx2*Bppx/(Bx2*Bx2) - this->_6p0*Appx*Bppx/(Bx2)
               + this->_6p0*Ax*Bppx*Bppx/Bx3 - this->_4p0*Bp*Apppx/Bx2
               - this->_4p0*Apx*Bpppx/Bx2 + this->_8p0*Ax*Bpx*Bpppx/Bx3
               + this->Apppp(x)/Bx - Ax*this->Bpppp(x)/Bx2;
    }*/
};
//------------------------------------------------------------------------------
//If we wish to define a simple rational or product function without
//any shared data, then use this.
template<class value_type,class arg_type>
class rational_function_shared_data
    : public rational_function<value_type,arg_type>
{
private:
    //Compute shared data, based on which derivative is requested:
    virtual void computeSharedData(arg_type x,int sharing_level) = 0;
    value_type operator()(arg_type x)
    {
        computeSharedData(x,0);
        return this->f(x);
    }
    value_type _1st_derivative(arg_type x)
    {
        computeSharedData(x,1);
        return this->d1(x);
    }
    value_type _2nd_derivative(arg_type x)
    {
        computeSharedData(x,2);
        return this->d2(x);
    }
    value_type _3rd_derivative(arg_type x)
    {
        computeSharedData(x,3);
        return this->d3(x);
    }
};

//------------------------------------------------------------------------------
//Product functions: f(x) = A(x)*B(x)
template<class value_type>
class product_function
{
private:
    const value_type _2p0;
    const value_type _6p0;
    const value_type _3p0;
    const value_type _4p0;
public:
    //Constructor:
    product_function() : _2p0(2.0),_6p0(6.0),_3p0(3.0),_4p0(4.0) {}
    //Virtual functions:
    virtual value_type A(value_type x) = 0;
    virtual value_type B(value_type x) = 0;
    virtual value_type Ap(value_type x) = 0;
    virtual value_type Bp(value_type x) = 0;
    virtual value_type App(value_type x) = 0;
    virtual value_type Bpp(value_type x) = 0;
    virtual value_type Appp(value_type x) = 0;
    virtual value_type Bppp(value_type x) = 0;
    virtual value_type Apppp(value_type x) = 0;
    virtual value_type Bpppp(value_type x) = 0;
    value_type operator()(value_type x)
    {
        return this->A(x)*this->B(x);
    }
    value_type d(value_type x)
    {
        return this->B(x)*this->Ap(x) + this->A(x)*this->Bp(x);
    }
    value_type d2(value_type x)
    {
        return this->_2p0*this->Ap(x)*this->Bp(x) + this->B(x)*this->App(x)
            + this->A(x)*this->Bpp(x);
    }
    value_type d3(value_type x)
    {
        return this->_3p0*this->Bp(x)*this->App(x)
        + this->_3p0*this->Ap(x)*this->Bpp(x)
        + this->A(x)*this->Bppp(x) + this->B(x)*this->Appp(x);
    }
    value_type d4(value_type x)
    {
        return this->A(x)*this->Bpppp(x) + this->_4p0*this->Ap(x)*this->Bppp(x)
            + this->_6p0*this->App(x)*this->Bpp(x)
            + this->_4p0*this->Appp(x)*this->Bp(x) + this->Apppp(x)*this->B(x);
    }
};
#endif //RATIONAL_H
