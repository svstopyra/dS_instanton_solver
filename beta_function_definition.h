#ifndef BETA_FUNCTION_DEFINITION_H
#define BETA_FUNCTION_DEFINITION_H
//Definition of beta functions and their derivatives.
#include "ode_rhs.h"
template<class value_type>
class betaFunctionSet : public odeRHS<std::vector<value_type>,value_type,value_type>
{
private:
virtual void beta(const std::vector<value_type>& g,
                  std::vector<value_type>& beta_g,const value_type& t) = 0;
};
//Versions with higher derivatives:
template<class value_type>
class betaFunctionWithDerivatives
: public betaFunctionSet<value_type>
{
public:
    //Give the class the internal capacity to store a "default"
    //derivative to compute when operator() is called.
    int nDerivative;
    betaFunctionWithDerivatives(int N)
    {
        nDerivative = N;
    }
    betaFunctionWithDerivatives()
    {
        nDerivative = 0;
    }
    //Call an arbitrary derivative function:
    virtual void betap(const std::vector<value_type>& g,
                  std::vector<value_type>& beta_g,const value_type& t,
                  int nDer) = 0;
    //Call the default derivative function:
    void betap(const std::vector<value_type>& g,
                  std::vector<value_type>& beta_g,const value_type& t)
    {
        this->betap(g,beta_g,t,this->nDerivative);
    }
};
#endif // BETA_FUNCTION_DEFINITION_H
