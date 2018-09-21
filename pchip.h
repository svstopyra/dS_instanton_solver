#ifndef PCHIP_H
#define PCHIP_H

//Computes the coefficients requires for piecewise cubic hermite polynomial
//interpolation.

template<class value_type>
void pchipk(std::vector<value_type>& k,
            const std::vector<value_type>& x,const std::vector<value_type>& y)
{
    int nSize = int(x.size());
    //Compute a difference vector, h:
    const value_type _0p0 = value_type(0.0);
    const value_type _2p0 = value_type(2.0);
    const value_type _3p0 = value_type(3.0);
    //Vector of secants:
    std::vector<value_type> secants (nSize - 1,_0p0);
    std::vector<value_type> h (nSize - 1,_0p0);


    for(int i = 0;i < nSize - 1;i++)
    {
        h[i] = x[i+1] - x[i];
        secants[i] = (y[i+1] - y[i])/h[i];
    }

    //Initialise output vector:
    k.assign(nSize,_0p0);

    //Special case for n = 2:
    if(nSize == 2)
    {
        k[0] = secants[0];
        k[1] = secants[0];
    }
    else
    {
        k[0] = (_2p0*h[0] + h[1])*secants[0] - h[0]*secants[1]/(h[0] + h[1]);
        k[nSize - 1] = ( (_2p0*h[nSize - 2] + h[nSize - 3])*secants[nSize - 2]
                        - h[nSize - 2]*secants[nSize - 3])
                        /(h[nSize -2] + h[nSize - 3]);
        for(int i = 1;i < nSize - 1;i++)
        {
            if(secants[i-1]*secants[i] > _0p0)
            {
                //Regular point:
                value_type sum = h[i-1] + h[i];
                value_type kmin = std::min(abs(secants[i-1]),
                                           abs(secants[i]));
                value_type kmax = std::max(abs(secants[i-1]),
                                           abs(secants[i]));
                value_type weight1 = (h[i-1] + sum)/(_3p0*sum);
                value_type weight2 = (h[i] + sum)/(_3p0*sum);
                k[i] = kmin/( weight1*(secants[i-1]/kmax)
                             + weight2*(secants[i]/kmax) );
            }
            else
            {
                //Local extremum. Note that this
                //Includes the secants[i-1] = 0 or
                //secants[i] = 0 case:
                k[i] = _0p0;
            }
        }
    }
}


#endif // PCHIP_H
