#include <time.h> 
#include <ctime>
#include "methods.h"
#include "io.h"


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi::size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi::rank);
        
    srand(mpi::rank*time(0));
    //std::string confFile;
    //std::string rhsFile;
    
    
    double m0; //bare mass
    
    CG::max_iter = 10000; //Maximum number of iterations for the conjugate gradient method
    CG::tol = 1e-15; //Tolerance for convergence

    if (mpi::rank == 0){
         //---Input data---//
        std::cerr << " ------------------------------" << std::endl;
        std::cerr << "|      CG Schwinger Model      |" << std::endl;
        std::cerr << " ------------------------------" << std::endl;
        std::cerr << "Nx " << LV::Nx << " Nt " << LV::Nt << std::endl;
        std::cerr << "ranks_x: ";
        std::cin >> mpi::ranks_x;
        std::cerr << "ranks_t: ";
        std::cin >> mpi::ranks_t;
        std::cerr << "m0: ";
        std::cin >> m0;
        //std::cout << "Configuration file path: ";
        //std::cin >> confFile;
        //std::cout << "RHS file path: ";
        //std::cin >> rhsFile;
    }

    MPI_Bcast(&mpi::ranks_x, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&mpi::ranks_t, 1, MPI_INT,  0, MPI_COMM_WORLD);
    MPI_Bcast(&m0, 1, MPI_DOUBLE,  0, MPI_COMM_WORLD);
    //broadcast_file_name(confFile);
    //broadcast_file_name(rhsFile);
    mass::m0 = m0;


    allocate_lattice_arrays(); 
    initializeMPI(); //2D rank topology
    boundaries();
    //--------------------------------------//

    srand((mpi::rank2d+1));
    spinor U(mpi::maxSizeH);
    spinor rhs(mpi::maxSizeH);
    spinor sol(mpi::maxSizeH);   //Zero vector as initial solution
    spinor x0(mpi::maxSizeH);

    for(int x = 1; x<=mpi::width_x; x++){
		for(int t = 1; t<=mpi::width_t; t++){
			int n = x*(mpi::width_t+2)+t;
			//mu = 0
			U.val[2*n] = RandomU1();
            //mu = 1
            U.val[2*n+1] = RandomU1(); 

            rhs.val[2*n] = RandomU1();
            rhs.val[2*n+1] = RandomU1(); 

        }
    }
    
    //read_binary(confFile,U);
    //read_binary(rhsFile,rhs);
    exchange_halo(U.val);
    
    double tol = CG::tol;
    Methods Meth(U, rhs, x0 ,m0, tol);

    Meth.CG(true);
    Meth.CG_S(true);  

    //Free coordinate arrays
    free_lattice_arrays();
    MPI_Finalize();

	return 0;
}

//-0.18840579710144945
//-0.1868