#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>

#ifndef _MODULE_H_
#define _MODULE_H_

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
void matrix_n_m_compute_kernel(const double* k_du_sin_cos, const double* k_du_sin_sin, double neg_n, double neg_m, double alpha, double PHI, double* dest_real, double* dest_img, size_t size);

void matrix_n_m_compute(const cuda_matrix& k_du_sin_cos, const cuda_matrix& k_du_sin_sin, double neg_n, double neg_m, double alpha, double PHI, cuda_cmatrix& dest);

void gsl_matrix_to_cuda_matrix(cuda_matrix& cuda, const gsl_matrix* gsl);

void cuda_matrix_to_gsl_matrix(gsl_matrix* gsl, const cuda_matrix& cuda);

void cuda_cmatrix_to_gsl_cmatrix(gsl_matrix_complex* gsl, const cuda_cmatrix& cuda);

void cuda_matrix_alloc(cuda_matrix& cuda, unsigned int rows, unsigned int cols);

void cuda_cmatrix_alloc(cuda_cmatrix& cuda, unsigned int rows, unsigned int cols);

void cuda_matrix_free(cuda_matrix& cuda);

void cuda_cmatrix_free(cuda_cmatrix& cuda);

} // namespace withcuda

#endif
