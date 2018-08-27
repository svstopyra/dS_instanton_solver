#ifndef MULTIPRECISION_DEFS
#define MULTIPRECISION_DEFS
//Header file defining the 'multi' and 'multi_etoff' types, so we can quickly change these.
//Shared among multiple parts of the project!

//#include <boost/multiprecision/cpp_dec_float.hpp>//Using mpfr, so don't need this.
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/constants/constants.hpp>//Needed by most files using multi.
#include <boost\math\special_functions.hpp>//Needed by most files using multi.

//Backend we are using to handle numbers (eg, cpp_dec_float from BOOST or
//mpfr_float_50 from MPFR library  etc ...)
typedef boost::multiprecision::mpfr_float_backend< 100 , boost::multiprecision::allocate_dynamic > mp_backend;
//Use a version with expression templates switched off, because the eigen
//library doesn't support them. This will slow down arithmetic, but it is
//the only way that is compatible with the eigen matrix library.
typedef boost::multiprecision::number< mp_backend , boost::multiprecision::et_off> multi_etoff;
//typedef boost::multiprecision::number< mp_backend , boost::multiprecision::et_on> multi;
typedef boost::multiprecision::number< mp_backend, boost::multiprecision::et_off> multi;//Test whether this compiles better (or at all, for that matter)


//Portable method for inter-converting between types:

//Overloaded version:
template<class type1,class type2>
type2 convert_type(type1 arg)
{
    return type2(arg);
}
//Special cases:
template<>
int convert_type<multi_etoff,int>(multi_etoff arg);
//Overloaded type for converting doubles to other types:
template<>
int convert_type<double,int>(double arg);

template<>
double convert_type<multi_etoff,double>(multi_etoff arg);

#endif
