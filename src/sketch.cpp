#include "sketch.h"


/**
 * Generates a Gaussian (random) sketch matrix S ∈ ℝ^(m×n)
 * with entries drawn from N(0, 1/m).
 * 
 * n  Number of columns (original dimension)
 * m  Number of rows (sketch dimension, m << n)
 * 
 */
void gaussian_sketch(int n, int m, double* S){

    // Create random number generator with seed
    static std::random_device rd;
	static std::default_random_engine generator(rd());
	std::normal_distribution<double> distribution(0.0, 1.0); //mu, std
    
    // Scaling factor: 1/sqrt(m)
    double scale = 1.0 / std::sqrt(static_cast<double>(m));
    
    // Allocate and fill the sketch matrix

    std::vector<std::vector<double>> S(m, std::vector<double>(n));

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            S[i*(n+2)+j] = distribution(generator) * scale;
        }
    }
}