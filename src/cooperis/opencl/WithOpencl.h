#include <gsl/gsl_matrix.h>
#include <iostream>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

class WithOpencl {
private:
    cl_context context;
    cl_device_id device;
    cl_command_queue queue;
    cl_program program;
    cl_kernel kernel;
    cl_platform_id platform;

    void cl_assert(cl_int, const char*, int, const char* = NULL);

public:
    struct cl_matrix {
        cl_mem data;
        unsigned int rows;
        unsigned int cols;
    };
    struct cl_cmatrix {
        cl_mem data_real;
        cl_mem data_img;
        unsigned int rows;
        unsigned int cols;
    };

    /**
     * @param platform_id the platform id
     * @param device_id the device id
     */
    WithOpencl(int platform_id, int device_id);
    virtual ~WithOpencl();

    /**
     * @param k_du_sin_cos the precomputed sin_cos matrix
     * @param k_du_sin_sin the precomputed sin_sin matrix
     * @param n the number part of the gain
     * @param m the number part of the gain
     * @param alpha alpha value
     * @param PHI PHI value
     * @param phase the output phase matrix
     */
    void gain_compute_phase(const cl_matrix& k_du_sin_cos, const cl_matrix& k_du_sin_sin, int& n, int& m, double& alpha, double& PHI, cl_cmatrix& phase);

    /**
     * @param gsl output gsl_matrix
     * @param cl input cl_matrix
     */
    void cl_cmatrix_to_gsl_cmatrix(gsl_matrix_complex* gsl, const cl_cmatrix& cl);

    /**
     * @param cl output cl_matrix
     * @param rows number of rows
     * @param cols number of columns
     * @param mem_flags opencl memory flags, if set to CL_MEM_COPY_HOST_PTR src_buff must be set
     * @param src_buff source buffer to copy from
     */
    void cl_matrix_alloc(cl_matrix& cl, unsigned int rows, unsigned int cols, cl_mem_flags mem_flags, double* src_buff = NULL);

    /**
     * @param cl output cl_cmatrix
     * @param rows number of rows
     * @param cols number of columns
     * @param mem_flags opencl memory flags, if set to CL_MEM_COPY_HOST_PTR src_real_buff and src_img_buff must be set
     * @param src_real_buff source buffer to copy from
     * @param src_img_buff source buffer to copy from
     */
    void cl_cmatrix_alloc(cl_cmatrix& cl, unsigned int rows, unsigned int cols, cl_mem_flags mem_flags, double* src_real_buff = NULL, double* src_img_buff = NULL);

    /**
     * @param cl input cl_matrix to free
     */
    void cl_matrix_free(cl_matrix&);

    /**
     * @param cl input cl_cmatrix to free
     */
    void cl_cmatrix_free(cl_cmatrix&);

    /**
     * @param cl input cl_cmatrix to zero
     */
    void zero_cl_cmatrix(cl_cmatrix&);
};
