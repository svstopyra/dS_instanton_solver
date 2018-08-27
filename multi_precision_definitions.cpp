#include "multi_precision_definitions.h"

//Explicit instantiation file for type conversion functions.
//We forward declare these so that they are available everywhere, and compile
//them here.


//Special cases:
template<>
int convert_type<multi_etoff,int>(multi_etoff arg)
{
    return arg.template convert_to<int>();
}
//Overloaded type for converting doubles to other types:
template<>
int convert_type<double,int>(double arg)
{
    return int(arg);
}

template<>
double convert_type<multi_etoff,double>(multi_etoff arg)
{
    return arg.template convert_to<double>();
}
