#include "methods.h"

void Methods::CG(const bool print){
    if (mpi::rank2d == 0)
        std::cout << "--------Inverting the normal equations with CG----------" << std::endl; 

    localFLOPS = 0;

    spinor D_dagg_psi(mpi::maxSizeH);
    D_dagger_phi(U, rhs, D_dagg_psi, m0);
    start = MPI_Wtime();
    iter_counters::CGIt = conjugate_gradient(U, D_dagg_psi, xCG, m0,print);
    end = MPI_Wtime();

    mpi_reduceFLOPS();
    printFLOPS(FLOPS);
    flop_counters::CGFlops=FLOPS;
    if (mpi::rank2d == 0)
        printf("[rank %d] time elapsed CG: %.4fs.\n\n", mpi::rank2d, end - start);  
}

void Methods::CG_S(const bool print){
    if (mpi::rank2d == 0)
        std::cout << "--------Inverting the normal equations with CG and sketching----------" << std::endl; 

    localFLOPS = 0;

    spinor D_dagg_psi(mpi::maxSizeH);
    D_dagger_phi(U, rhs, D_dagg_psi, m0);
    start = MPI_Wtime();
    iter_counters::CGIt = conjugate_gradient(U, D_dagg_psi, xCG, m0,print);
    end = MPI_Wtime();

    mpi_reduceFLOPS();
    printFLOPS(FLOPS);
    flop_counters::CGFlops=FLOPS;
    if (mpi::rank2d == 0)
        printf("[rank %d] time elapsed CG with sketching: %.4fs.\n\n", mpi::rank2d, end - start);  
}

void Methods::check_solution(const spinor& x_sol){
    spinor x(mpi::maxSizeH);
    D_phi(U,x_sol,x,m0);
    int index;
    for (int nx=1;nx<=mpi::width_x;nx++){
    for(int nt=1;nt<=mpi::width_t;nt++){
    for(int mu=0; mu<2; mu++){
        index = idx(nx,nt,mu);
        if (std::abs(x.val[index]-rhs.val[index]) > 1e-6 ){
            std::cout << "Solution differs at nx" << nx << " nt " << nt << " mu " << mu << " rank" << mpi::rank2d << std::endl;
            std::cout << x.val[index] << "   " << rhs.val[index] << std::endl;
        }
    }
    }
    }

    if (mpi::rank2d == 0)
        std::cout << "Solution verified" << std::endl;
}