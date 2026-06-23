#include "variables.h"

double pi=3.14159265359;
c_double I_number(0, 1); //imaginary number

long long int FLOPS=0;
long long int localFLOPS=0;

namespace iter_counters{
    long long int kCycleIt=-1;
    long long int vCycleIt=-1;
    long long int CGIt=-1;
    long long int BiCGIt=-1;
    long long int GMRESIt=-1;
    long long int SAPIt=-1;
    long long int FGMRES_SAPIt=-1;
}
namespace flop_counters{
    long long int CGFlops=-1;
}

MPI_Datatype global_conf_type;
MPI_Datatype global_conf_resized;
MPI_Datatype local_conf_type;
MPI_Datatype local_conf_resized;
MPI_Datatype column_type;

/*
	Vectorized lattice coords.*/
int Coords(const int& x, const int& t){
	return x*mpi::width_t + t;
}

namespace mass{double m0;}


namespace mpi{
    int rank = 0;
    int size = 1; 
    int maxSize = LV::Ntot; //Default value, will be updated in main
    int maxSizeH = LV::Ntot;
    int sitesH = LV::Ntot;
    int ranks_x = 1;
    int ranks_t = 1;
    int width_x = LV::Nx;
    int width_t = LV::Nt;
    int rank2d = 0; //linearize rank
    int coords[2] = {0,0}; //rank coords
    int top = 0; 
    int bot = 0; 
    int right = 0; 
    int left = 0;
    int bot_left = 0;
    int bot_right = 0;
    int top_left = 0;
    int top_right = 0;
    MPI_Comm cart_comm;
    
}

namespace CG{
	int max_iter = 10000; //Maximum number of iterations for the conjugate gradient method
	double tol = 1e-10; //Tolerance for convergence
}

int* lpb = nullptr;
int* rpb = nullptr;
c_double* lsign = nullptr;
c_double* rsign = nullptr;

void allocate_lattice_arrays() {
    lpb     = new int[mpi::maxSizeH]();
    rpb     = new int[mpi::maxSizeH]();
    lsign   = new c_double[mpi::maxSizeH]();
    rsign   = new c_double[mpi::maxSizeH]();
}

void free_lattice_arrays() {
    delete[] lpb;
    delete[] rpb;
    delete[] lsign;
    delete[] rsign;    
}