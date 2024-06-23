#include <driver_types.h>
#include <gsl/gsl_matrix.h>
#include <cuda_runtime.h>
#include <iostream>
#include "WithCuda.h"

namespace withcuda {

void cuda_assert(cudaError_t code, const char* file = __FILE__, int line = __LINE__)
{
    if (code != cudaSuccess) {
        std::cerr << "CUDA error: " << cudaGetErrorString(code) << " at " << file << ":" << line << std::endl;
        exit(code);
    }
}

__global__ void gain_compute_phase_kernel(const double* k_du_sin_cos, const double* k_du_sin_sin, double n, double m, double alpha, double PHI, double* phase_real, double* phase_img, size_t size)
{
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) {
        double tmp_img = -((-n * k_du_sin_cos[i]) + (-m * k_du_sin_sin[i]) + alpha + PHI);
        phase_real[i] += cos(tmp_img);
        phase_img[i] += sin(tmp_img);
    }
}

void gain_compute_phase(int max_threads_per_block, const cuda_matrix& k_du_sin_cos, const cuda_matrix& k_du_sin_sin, double n, double m, double alpha, double PHI, cuda_cmatrix& phase)
{
    int num_blocks = (k_du_sin_cos.rows * k_du_sin_cos.cols + max_threads_per_block - 1) / max_threads_per_block;
    gain_compute_phase_kernel<<<num_blocks, max_threads_per_block>>>(k_du_sin_cos.data, k_du_sin_sin.data, n, m, alpha, PHI, phase.data_real, phase.data_img, phase.rows * phase.cols);
}

void gsl_matrix_to_cuda_matrix(cuda_matrix& cuda, const gsl_matrix* gsl)
{
    cuda.rows = gsl->size1;
    cuda.cols = gsl->size2;
    cuda_assert(cudaMalloc(&cuda.data, cuda.rows * cuda.cols * sizeof(double)), __FILE__, __LINE__);
    cuda_assert(cudaMemcpy(cuda.data, gsl->data, cuda.rows * cuda.cols * sizeof(double), cudaMemcpyHostToDevice), __FILE__, __LINE__);
}

void cuda_matrix_to_gsl_matrix(gsl_matrix* gsl, const cuda_matrix& cuda)
{
    cuda_assert(cudaMemcpy(gsl->data, cuda.data, cuda.rows * cuda.cols * sizeof(double), cudaMemcpyDeviceToHost), __FILE__, __LINE__);
}

void cuda_cmatrix_to_gsl_cmatrix(gsl_matrix_complex* gsl, const cuda_cmatrix& cuda)
{
    double* data_real = (double*) malloc(cuda.rows * cuda.cols * sizeof(double));
    double* data_img = (double*) malloc(cuda.rows * cuda.cols * sizeof(double));
    cuda_assert(cudaMemcpy(data_real, cuda.data_real, cuda.rows * cuda.cols * sizeof(double), cudaMemcpyDeviceToHost), __FILE__, __LINE__);
    cuda_assert(cudaMemcpy(data_img, cuda.data_img, cuda.rows * cuda.cols * sizeof(double), cudaMemcpyDeviceToHost), __FILE__, __LINE__);
    for (unsigned int i = 0; i < cuda.rows; i++)
        for (unsigned int j = 0; j < cuda.cols; j++)
            gsl_matrix_complex_set(gsl, i, j, {data_real[i * cuda.cols + j], data_img[i * cuda.cols + j]});
    free(data_real);
    free(data_img);
}

void cuda_matrix_alloc(cuda_matrix& cuda, unsigned int rows, unsigned int cols)
{
    cuda.rows = rows;
    cuda.cols = cols;
    cuda_assert(cudaMalloc(&cuda.data, rows * cols * sizeof(double)), __FILE__, __LINE__);
    cuda_assert(cudaMemset(cuda.data, 0, rows * cols * sizeof(double)), __FILE__, __LINE__);
}

void cuda_cmatrix_alloc(cuda_cmatrix& cuda, unsigned int rows, unsigned int cols)
{
    cuda.rows = rows;
    cuda.cols = cols;
    cuda_assert(cudaMalloc(&cuda.data_real, rows * cols * sizeof(double)), __FILE__, __LINE__);
    cuda_assert(cudaMalloc(&cuda.data_img, rows * cols * sizeof(double)), __FILE__, __LINE__);
    cuda_assert(cudaMemset(cuda.data_real, 0, rows * cols * sizeof(double)), __FILE__, __LINE__);
    cuda_assert(cudaMemset(cuda.data_img, 0, rows * cols * sizeof(double)), __FILE__, __LINE__);
}

void cuda_matrix_free(cuda_matrix& cuda)
{
    cuda_assert(cudaFree(cuda.data), __FILE__, __LINE__);
    cuda.rows = 0;
    cuda.cols = 0;
    cuda.data = nullptr;
}

void cuda_cmatrix_free(cuda_cmatrix& cuda)
{
    cuda_assert(cudaFree(cuda.data_real), __FILE__, __LINE__);
    cuda_assert(cudaFree(cuda.data_img), __FILE__, __LINE__);
    cuda.rows = 0;
    cuda.cols = 0;
    cuda.data_real = nullptr;
    cuda.data_img = nullptr;
}

void set_cuda_device(int device)
{
    if (device < 0)
        throw std::runtime_error("Invalid cuda device id, must be >= 0");
    cuda_assert(cudaSetDevice(device), __FILE__, __LINE__);
}

int get_cuda_max_threads_per_block()
{
    int device;
#ifndef CUDA_DEVICE_ID
    device = cudaGetDevice(&device);
#else
    device = CUDA_DEVICE_ID;
#endif

    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device);

    return prop.maxThreadsPerBlock;
}

} // namespace withcuda
