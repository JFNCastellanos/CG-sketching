#ifndef METHODS_H
#define METHODS_H

#include "conjugate_gradient.h"

//Class for comparing different matrix-inversion methods.
class Methods {
public:
    Methods(const spinor& U, const spinor& rhs, const spinor& x0 ,const double m0, const double tol): 
    U(U), rhs(rhs), x0(x0), m0(m0),tol(tol){
        //Solution buffers
        xCG         = spinor(mpi::maxSizeH);
        xCG_S       = spinor(mpi::maxSizeH);

    }
    ~Methods(){}


    void CG(const bool print);
    void CG_S(const bool print);

    void check_solution(const spinor& x_sol);

    spinor xCG;
    spinor xCG_S;

private:

    const spinor U;
    const spinor rhs;
    const spinor x0;
    const double m0;
    const double tol;
    double start, end;

};


#endif