#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>

#ifndef _MODULE_H_
#define _MODULE_H_
#endif

#ifdef __CUDACC__
#define CUDA_GLOBAL __global__
#else
#define CUDA_GLOBAL
#endif

namespace withcuda {

struct cuda_matrix {
    double* data;
    unsigned int rows;
    unsigned int cols;
};

struct cuda_cmatrix {
    double* data_real;
    double* data_img;
    unsigned int rows;
    unsigned int cols;
};

CUDA_GLOBAL
void gain_compute_phase_kernel(const double* k_du_sin_cos, const double* k_du_sin_sin, double n, double m, double alpha, double PHI, double* phase_real, double* phase_img, size_t size);

/**
 * @param k_du_sin_cos the precomputed sin_cos matrix
 * @param k_du_sin_sin the precomputed sin_sin matrix
 * @param n the number part of the gain
 * @param m the number part of the gain
 * @param alpha alpha value
 * @param PHI PHI value
 * @param phase the output phase matrix
 */
void gain_compute_phase(const cuda_matrix& k_du_sin_cos, const cuda_matrix& k_du_sin_sin, double n, double m, double alpha, double PHI, cuda_cmatrix& phase);

/**
 * @param cuda output matrix
 * @param gsl input matrix
 */
void gsl_matrix_to_cuda_matrix(cuda_matrix& cuda, const gsl_matrix* gsl);

/**
 * @param gsl output matrix
 * @param cuda input matrix
 */
void cuda_matrix_to_gsl_matrix(gsl_matrix* gsl, const cuda_matrix& cuda);

/**
 * @param gsl output cmatrix
 * @param cuda input cmatrix
 */
void cuda_cmatrix_to_gsl_cmatrix(gsl_matrix_complex* gsl, const cuda_cmatrix& cuda);

/**
 * @param cuda output matrix
 * @param rows number of rows
 * @param cols number of cols
 */
void cuda_matrix_alloc(cuda_matrix& cuda, unsigned int rows, unsigned int cols);

/**
 * @param cuda output cmatrix
 * @param rows number of rows
 * @param cols number of cols
 */
void cuda_cmatrix_alloc(cuda_cmatrix& cuda, unsigned int rows, unsigned int cols);

/**
 * @param cuda martrix to free
 */
void cuda_matrix_free(cuda_matrix& cuda);

/**
 * @param cuda cmatrix to free
 */
void cuda_cmatrix_free(cuda_cmatrix& cuda);

} // namespace withcuda
