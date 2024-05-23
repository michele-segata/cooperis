R"(
  __kernel void gain_kernel (__global double *dest_real, __global double *dest_img, __global double *k_du_sin_cos, __global double *k_du_sin_sin, __private double alpha, __private double PHI, __private int n, __private int m, __private ulong size) {
    int i = get_global_id(0);

    if (i < size) {
        double tmp_img = -((-n * k_du_sin_cos[i]) + (-m * k_du_sin_sin[i]) + alpha + PHI);
        dest_real[i] += cos(tmp_img);
        dest_img[i] += sin(tmp_img);
    }
  }
)"
