#include <gsl/gsl_matrix.h>
#include <iostream>
#include <vector>
#include "WithOpencl.h"

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#ifndef CL_PLATFORM_ID
#define CL_PLATFORM_ID 0
#endif

#ifndef CL_DEVICE_ID
#define CL_DEVICE_ID 0
#endif

#define CL_MAX_PLATFORMS 10
#define CL_MAX_DEVICES 10

WithOpencl::WithOpencl()
{

    if (CL_PLATFORM_ID < 0 || CL_PLATFORM_ID >= CL_MAX_PLATFORMS)
        throw std::runtime_error("Invalid platform id, must be between 0 and " + std::to_string(CL_MAX_PLATFORMS));

    if (CL_DEVICE_ID < 0 || CL_DEVICE_ID >= CL_MAX_DEVICES)
        throw std::runtime_error("Invalid device id, must be between 0 and " + std::to_string(CL_MAX_DEVICES));

    std::vector<cl_platform_id> platforms(CL_MAX_PLATFORMS);
    std::vector<cl_device_id> devices(CL_MAX_DEVICES);
    cl_uint num_platforms;
    cl_uint num_devices;

    // get the platform
    cl_assert(clGetPlatformIDs(CL_MAX_PLATFORMS, &platforms[0], &num_platforms), __FILE__, __LINE__, "error getting the platforms");

    if (CL_PLATFORM_ID >= num_platforms)
        throw std::runtime_error("Invalid platform id, must be between 0 and " + std::to_string(num_platforms));

    this->platform = platforms[CL_PLATFORM_ID];

    // get the device
    cl_assert(clGetDeviceIDs(this->platform, CL_DEVICE_TYPE_GPU, CL_MAX_DEVICES, &devices[0], &num_devices), __FILE__, __LINE__, "error getting the device");

    if (CL_DEVICE_ID >= num_devices)
        throw std::runtime_error("Invalid device id, must be between 0 and " + std::to_string(num_devices));

    this->device = devices[CL_DEVICE_ID];

    // check if the device supports double precision
    cl_device_fp_config config;
    cl_assert(clGetDeviceInfo(this->device, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof(config), &config, NULL), __FILE__, __LINE__, "error getting the device info");
    if (config == 0) {
        throw std::runtime_error("OpenCL device does not support double precision");
    }

    // create the context
    cl_int err;
    this->context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    cl_assert(err, __FILE__, __LINE__, "error creating the context");

    // create the command queue
#ifdef __APPLE__
    this->queue = clCreateCommandQueueWithPropertiesAPPLE(this->context, this->device, NULL, &err);
#else
    this->queue = clCreateCommandQueueWithProperties(this->context, this->device, NULL, &err);
#endif
    cl_assert(err, __FILE__, __LINE__, "error creating the command queue");

    // load the kernel
    const char* gain_kernel =
#include "./kernels/gain_kernel.cl"
        ;

    // create the program
    program = clCreateProgramWithSource(this->context, 1, &gain_kernel, NULL, &err);
    cl_assert(err, __FILE__, __LINE__, "error creating the program");

    // build the program
    err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        std::cerr << "Error building the program: " << err << std::endl;
        size_t len;
        char buffer[2048];
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        std::cerr << buffer << std::endl;
        cl_assert(err, __FILE__, __LINE__, "error building the program");
    }

    // create the kernel
    this->kernel = clCreateKernel(program, "gain_kernel", &err);
    cl_assert(err, __FILE__, __LINE__, "error creating the kernel");
}
WithOpencl::~WithOpencl()
{
    cl_assert(clReleaseKernel(this->kernel), __FILE__, __LINE__);
    cl_assert(clReleaseProgram(this->program), __FILE__, __LINE__);
    cl_assert(clReleaseCommandQueue(this->queue), __FILE__, __LINE__);
    cl_assert(clReleaseContext(this->context), __FILE__, __LINE__);
    cl_assert(clReleaseDevice(this->device), __FILE__, __LINE__);
}

void WithOpencl::cl_assert(cl_int ret_code, const char* file, int line, const char* msg)
{
    if (ret_code != CL_SUCCESS) {
        std::cerr << "OpenCL error " << ret_code << ": " << msg << " at " << file << ":" << line << std::endl;
        throw std::runtime_error("OpenCL error");
    }
}

void WithOpencl::gain_compute_phase(const cl_matrix& k_du_sin_cos, const cl_matrix& k_du_sin_sin, int& n, int& m, double& alpha, double& PHI, cl_cmatrix& phase)
{
    // set the kernel arguments
    cl_assert(clSetKernelArg(this->kernel, 0, sizeof(cl_mem), &phase.data_real), __FILE__, __LINE__);
    cl_assert(clSetKernelArg(this->kernel, 1, sizeof(cl_mem), &phase.data_img), __FILE__, __LINE__);
    cl_assert(clSetKernelArg(this->kernel, 2, sizeof(cl_mem), &k_du_sin_cos.data), __FILE__, __LINE__);
    cl_assert(clSetKernelArg(this->kernel, 3, sizeof(cl_mem), &k_du_sin_sin.data), __FILE__, __LINE__);
    cl_assert(clSetKernelArg(this->kernel, 4, sizeof(double), &alpha), __FILE__, __LINE__);
    cl_assert(clSetKernelArg(this->kernel, 5, sizeof(double), &PHI), __FILE__, __LINE__);
    cl_assert(clSetKernelArg(this->kernel, 6, sizeof(int), &n), __FILE__, __LINE__);
    cl_assert(clSetKernelArg(this->kernel, 7, sizeof(int), &m), __FILE__, __LINE__);
    size_t global_work_size = phase.rows * phase.cols;
    cl_assert(clSetKernelArg(this->kernel, 8, sizeof(cl_ulong), &global_work_size), __FILE__, __LINE__);

    // extract global work size and max work group size
    size_t max_work_group_size = 0;
    cl_assert(clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_work_group_size, NULL), __FILE__, __LINE__, "error getting the device info");
    if (global_work_size > max_work_group_size)
        global_work_size = max_work_group_size * (global_work_size / max_work_group_size + 1);

    // enqueue the kernel
    cl_assert(clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global_work_size, &max_work_group_size, 0, NULL, NULL), __FILE__, __LINE__, "error while enqueuing the kernel");
    cl_assert(clFinish(queue), __FILE__, __LINE__, "error while finishing the queue");
}

void WithOpencl::cl_cmatrix_to_gsl_cmatrix(gsl_matrix_complex* gsl, const cl_cmatrix& cl)
{
    double* data_real = (double*) malloc(cl.rows * cl.cols * sizeof(double));
    double* data_img = (double*) malloc(cl.rows * cl.cols * sizeof(double));

    clEnqueueReadBuffer(this->queue, cl.data_real, CL_TRUE, 0, gsl->size1 * gsl->size2 * sizeof(double), data_real, 0, NULL, NULL);
    clEnqueueReadBuffer(this->queue, cl.data_img, CL_TRUE, 0, gsl->size1 * gsl->size2 * sizeof(double), data_img, 0, NULL, NULL);
    clFinish(this->queue);
    for (unsigned int i = 0; i < gsl->size1; i++)
        for (unsigned int j = 0; j < gsl->size2; j++)
            gsl_matrix_complex_set(gsl, i, j, {data_real[i * gsl->size2 + j], data_img[i * gsl->size2 + j]});
    free(data_real);
    free(data_img);
}

void WithOpencl::cl_matrix_alloc(cl_matrix& cl, unsigned int rows, unsigned int cols, cl_mem_flags mem_flags, double* src_buff)
{
    cl.rows = rows;
    cl.cols = cols;
    cl_int err;
    cl.data = clCreateBuffer(this->context, mem_flags, rows * cols * sizeof(double), src_buff, &err);
    cl_assert(err, __FILE__, __LINE__);
}

void WithOpencl::cl_cmatrix_alloc(cl_cmatrix& cl, unsigned int rows, unsigned int cols, cl_mem_flags mem_flags, double* src_real, double* src_img)
{
    cl.rows = rows;
    cl.cols = cols;
    cl_int err;
    cl.data_real = clCreateBuffer(this->context, mem_flags, rows * cols * sizeof(double), src_real, &err);
    cl_assert(err, __FILE__, __LINE__);

    cl.data_img = clCreateBuffer(this->context, mem_flags, rows * cols * sizeof(double), src_img, &err);
    cl_assert(err, __FILE__, __LINE__);
}

void WithOpencl::cl_matrix_free(cl_matrix& cl)
{
    cl_assert(clReleaseMemObject(cl.data), __FILE__, __LINE__);
    cl.rows = 0;
    cl.cols = 0;
    cl.data = nullptr;
}

void WithOpencl::cl_cmatrix_free(cl_cmatrix& cl)
{
    cl_assert(clReleaseMemObject(cl.data_real), __FILE__, __LINE__);
    cl_assert(clReleaseMemObject(cl.data_img), __FILE__, __LINE__);
    cl.rows = 0;
    cl.cols = 0;
    cl.data_real = nullptr;
    cl.data_img = nullptr;
}

void WithOpencl::zero_cl_cmatrix(cl_cmatrix& cl)
{
    double zero = 0;
    cl_assert(clEnqueueFillBuffer(this->queue, cl.data_real, &zero, sizeof(double), 0, cl.rows * cl.cols * sizeof(double), 0, NULL, NULL), __FILE__, __LINE__);
    cl_assert(clEnqueueFillBuffer(this->queue, cl.data_img, &zero, sizeof(double), 0, cl.rows * cl.cols * sizeof(double), 0, NULL, NULL), __FILE__, __LINE__);
}
