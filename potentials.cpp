#include <project_specific.h>//Options specific to the project
//(eg, DLL export etc...). Provide an empty file or comment out if not needed.

//Code to be instantiated:
#include "potentials_code.h"

//Forward declarations for the versions of this that we need:
#ifdef USING_POTENTIAL_TEST1
template class testPotential< multi >;
#endif
#ifdef USING_SM_POTENTIAL
template class SMHiggsPotentialSpline< multi >;
template class SMHiggsPotentialSpline< double >;
//SM potential in dS space:
template class SMpotential_dS_1loop< multi >;
template class SMpotential_dS_1loop< double >;
template class betaSMdS_1loop< multi >;
template class betaSMdS_1loop< double >;
template class dtdtcl_1loop< multi >;
template class dtdtcl_1loop< double >;
template class dtdtcl_dS< multi >;
template class dtdtcl_dS< double >;
#endif
#ifdef USING_POLY_POTENTIAL
template class polyPotential< multi >;
#endif
template class log_potential< double >;
template class log_potential< multi >;
