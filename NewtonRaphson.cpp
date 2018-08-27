//c++ version of the globally convergent Newton-Raphson algorithm we are
//using for finding instantons.

//Contains code for implementation and explicit instantiation of NewtonRaphson.h
#include <project_specific.h>//Options specific to the project (eg, DLL export
//etc...). Provide an empty file or comment out if not needed.
//#include "multi_precision_definitions.h"
#include "Newton_Raphson_Code.h"



//Explicit instantiation:
#ifdef DE_SITTER_INSTANTONS
template class NRSolver< multi_etoff >;
template class functionNR< multi >;
template class jacobianNR< multi >;
#endif // DE_SITTER_INSTANTONS

