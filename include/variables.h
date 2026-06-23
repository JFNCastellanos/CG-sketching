#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED
#include "config.h"
#include <iostream>
#include <vector>
#include <complex>
#include <iomanip>
#include "mpi.h"

typedef std::complex<double> c_double;
extern double pi;
extern c_double I_number; //imaginary number

namespace mass{extern double m0;}

//FLOPs involved in different complex double operations
constexpr long long int ca = 2;     //  + complex addition = 2 FLOPs
constexpr long long int cm = 6;     //  * complex multiplication = 6 FLOPs
constexpr long long int cd = 11;    //  / complex division = 11 FLOPs (6 mul + 3 add + 2 div)
constexpr long long int da = 1;     //  double + double = 1 FLOP
constexpr long long int dm = 1;     //  double * double = 1 FLOP
constexpr long long int dd = 1;     //  double / double = 1 FLOP
constexpr long long int dcm = 2;    //  double * complex = 2 FLOPs
constexpr long long int dcd = 7;    //  double / complex_double = 7 FLOPs (4 mul + 1 add + 2 div)
constexpr long long int dsq = 1;    //  sqrt(double) = 1 FLOP

extern long long int FLOPS;
extern long long int localFLOPS;

namespace iter_counters{
    extern long long int CGIt;
}
namespace flop_counters{
    extern long long int CGFlops;
}


//Datatypes for reading/writing gauge confs and rhs
extern MPI_Datatype global_conf_type;
extern MPI_Datatype global_conf_resized;
extern MPI_Datatype local_conf_type;
extern MPI_Datatype local_conf_resized;
extern MPI_Datatype column_type;

//------------mpi settings-----------------//
namespace mpi{
    extern int rank;
    extern int size; 
    extern int maxSize;
    extern int maxSizeH; //maxSize with halos included
    extern int sitesH;
    extern int ranks_x;
    extern int ranks_t;
    extern int width_x;
    extern int width_t;
    extern int rank2d; //Rank id in the 2D communicator
    extern int coords[2];
    extern int top; 
    extern int bot; 
    extern int right; 
    extern int left;
    //Diagonal ranks necessary for staples
    extern int bot_left;
    extern int bot_right;
    extern int top_left;
    extern int top_right;
    extern MPI_Comm cart_comm;
}


//------------Lattice parameters--------------//
namespace LV {
    //Lattice dimensions//
    constexpr int Nx= NS; //We extract this value from config.h
    constexpr int Nt = NT; //We extract this value from config.h
    constexpr int Ntot = Nx*Nt; //Total number of lattice points
    constexpr int dof = 2;
}

//-----------CG solver------------//
namespace CG{
    extern int max_iter; //Maximum number of iterations for the conjugate gradient method
    extern double tol; //Tolerance for convergence
}

//Flattened spinor
struct spinor {
    c_double* val;  //Array with values
    int size;       //Size of array 
    //Constructor
    spinor(int N = LV::Ntot) : size(N) {
        val = new c_double[N]();
    }

    //Copy constructor (deep copy)
    spinor(const spinor& other) : size(other.size) {
        val = new c_double[size];
        std::copy(other.val, other.val + size, val);
    }

    //Assignment operator (deep copy)
    spinor& operator=(const spinor& other) {
        if (this != &other) {
            if (size != other.size) {
                delete[] val;
                size = other.size;
                val = new c_double[size];
            }
            std::copy(other.val, other.val + size, val);
        }
        return *this;
    }

    // Destructor
    ~spinor() {
        delete[] val;
    }

    inline void clearBuffer(){
        for(int n = 0; n<size; n++){
            val[n] = 0;
        }
    }
};


//Boundary conditions with padding.
extern int* lpb;
extern int* rpb;
extern c_double* lsign;
extern c_double* rsign;

//Allocation of many arrays used throughout the code 
void allocate_lattice_arrays();
void free_lattice_arrays();




#endif 