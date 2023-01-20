//
// Copyright (C) 2022 Michele Segata <segata@ccs-labs.org>
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//

#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <veins-ris/utility/ReconfigurableIntelligentSurface.h>

ReconfigurableIntelligentSurface::ReconfigurableIntelligentSurface(double frequency, int n, int cellsPerLambda, int lambdaSize)
{

    printf("ReconfigurableIntelligentSurface: freq=%f Hz, n=%d, cellsPerLambda=%d, lambdaSize=%d, totalCells=%f\n", frequency, n, cellsPerLambda, lambdaSize, pow(cellsPerLambda*lambdaSize, 2));
    E = linspace(0, pow(2, n) - 1, pow(2, n));
    gsl_vector_scale(E, M_PI_X_2 / pow(2, n));
    c = 299792458;
    f = frequency;
    lambda = c / f;
    k = M_PI_X_2 / lambda;
    du = lambda / cellsPerLambda;
    DN = lambdaSize * lambda;
    DM = DN;
    M = round(DM / du);
    N = round(DN / du);
    du_k = du * k;
    dG = du;
    Ps = 361;
    Ts = 91;
    // phi = linspace(-1 * M_PI/180,M_PI,Ps);
    // theta = linspace(-1 * M_PI/180,92*M_PI/180,Ts);
    phi = linspace(0, M_PI_X_2, Ps);
    theta = linspace(0, M_PI_2, Ts);
    Dt = gsl_vector_get(theta, 1) - gsl_vector_get(theta, 0); // theta[1] - theta[0];
    Df = gsl_vector_get(phi, 1) - gsl_vector_get(phi, 0); // phi[1] - phi[0];
    cos_phi = matrix_from_vector(phi);
    apply_in_place(cos_phi, cos);
    sin_phi = matrix_from_vector(phi);
    apply_in_place(sin_phi, sin);
    sin_theta = matrix_from_vector(theta);
    apply_in_place(sin_theta, sin);
    k_dG_sin_cos = new_matrix(Ts, Ps);
    k_dG_sin_sin = new_matrix(Ts, Ps);
    // k * dG * sin(theta)*cos(phi)^T
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, sin_theta, cos_phi, 0, k_dG_sin_cos);
    gsl_matrix_scale(k_dG_sin_cos, k * dG);

    // k * dG * sin(theta)*sin(phi)^T
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, sin_theta, sin_phi, 0, k_dG_sin_sin);
    gsl_matrix_scale(k_dG_sin_sin, k * dG);

    GSL_REAL(C_0_M1j) = 0;
    GSL_IMAG(C_0_M1j) = -1;

    wt = 1;

    sin_abs_theta = matrix_from_vector(theta);
    abs(sin_abs_theta);
    apply_in_place(sin_abs_theta, sin);

    phi_Dt_Df = matrix_from_vector(phi);
    gsl_matrix_scale(phi_Dt_Df, Dt * Df);

    coding = new_matrix(M, N);
    coding_alpha = new_matrix(M, N);

    F2 = new_matrix(Ts, Ps);

    // startup configuration
    configureMetaSurface(0, 0, 0, 0);

}

ReconfigurableIntelligentSurface::~ReconfigurableIntelligentSurface()
{
    gsl_matrix_free(F2);
    gsl_matrix_free(coding_alpha);
    gsl_matrix_free(coding);
    gsl_matrix_free(phi_Dt_Df);
    gsl_matrix_free(sin_abs_theta);
    gsl_matrix_free(k_dG_sin_sin);
    gsl_matrix_free(k_dG_sin_cos);
    gsl_matrix_free(sin_theta);
    gsl_matrix_free(sin_phi);
    gsl_matrix_free(cos_phi);
    gsl_vector_free(theta);
    gsl_vector_free(phi);
    gsl_vector_free(E);
}

void ReconfigurableIntelligentSurface::configureMetaSurface(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad)
{

    // TODO: remove this when the model with incidence angle comes in
//    phiI_rad = 0;
//    thetaI_rad = 0;

    if (phiR_rad < KEEP_SAME_ANGLE)
        configPhiR = phiR_rad;
    else
        phiR_rad = configPhiR;
    if (phiI_rad < KEEP_SAME_ANGLE)
        configPhiI = phiI_rad;
    else
        phiI_rad = configPhiI;
    if (thetaR_rad < KEEP_SAME_ANGLE)
        configThetaR = thetaR_rad;
    else
        thetaR_rad = configThetaR;
    if (thetaI_rad < KEEP_SAME_ANGLE)
        configThetaI = thetaI_rad;
    else
        thetaI_rad = configThetaI;

//    printf("Configuring metasurface for phiR=%.2f, thetaR=%.2f, phiI=%.2f, thetaI=%.2f\n", RAD_TO_DEG(phiR_rad), RAD_TO_DEG(thetaR_rad), RAD_TO_DEG(phiI_rad), RAD_TO_DEG(thetaI_rad));

    phiR_rad = REF_TO_MATH_PHI(phiR_rad);
    thetaR_rad = REF_TO_MATH_THETA(thetaR_rad);
    phiI_rad = REF_TO_MATH_PHI(phiI_rad);
    thetaI_rad = REF_TO_MATH_THETA(thetaI_rad);

    double phiR = phiR_rad - M_PI;
    double thetaR = thetaR_rad;

    double dsx = du_k * (cos(phiR) * sin(thetaR) - sin(thetaI_rad) * cos(phiI_rad));
    double dsy = du_k * (sin(phiR) * sin(thetaR) - sin(thetaI_rad) * sin(phiI_rad));

    // B and coding are pointers, so we work on coding actually
    Matrix B = coding;
    // same for alpha and coding_alpha
    Matrix alpha = coding_alpha;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            gsl_matrix_set(B, i, j, (i+1) * dsx + (j+1) * dsy);
            gsl_matrix_set(alpha, i, j, M_PI_X_2 * ((i+1) * du * sin(thetaI_rad) * cos(phiI_rad) + (j + 1) * du * sin(thetaI_rad) * sin(phiI_rad)) / lambda);

            gsl_matrix_set(B, i, j, fmod(gsl_matrix_get(B, i, j), M_PI_X_2));

            Vector tmp = new_vector(E);
            gsl_vector_add_constant(tmp, -std::abs(gsl_matrix_get(B, i, j)));
            abs(tmp);
            int p = min_pos(tmp);
            gsl_vector_free(tmp);

            if (gsl_matrix_get(B, i, j) < 0)
                gsl_matrix_set(B, i, j, -gsl_vector_get(E, p));
            else
                gsl_matrix_set(B, i, j, gsl_vector_get(E, p));

        }
    }

    recomputeGainMap = true;

}

double ReconfigurableIntelligentSurface::cachedGain(double phiR_rad, double thetaR_rad) const
{

    // assume that no signal gets radiated behind the surface
    if (!(thetaR_rad >= 0 && thetaR_rad <= M_PI_2))
        return 0;

    // convert the ranges
    phiR_rad = REF_TO_MATH_PHI(phiR_rad);
    thetaR_rad = REF_TO_MATH_THETA(thetaR_rad);

    int phi_deg = RAD_TO_DEG_ROUND(phiR_rad);
    int theta_deg = RAD_TO_DEG_ROUND(thetaR_rad);
    // * 0.9 is an efficiency factor that was applied in the original matlab script
    return gsl_matrix_get(F2, theta_deg, phi_deg) * 4 * M_PI / v * 0.9;
}

bool ReconfigurableIntelligentSurface::canUseCache(double phiI_rad, double thetaI_rad) const
{
    if (recomputeGainMap)
        return false;
    if (RAD_TO_DEG_ROUND(phiI_rad) != cached_phiI_deg || RAD_TO_DEG_ROUND(thetaI_rad) != cached_thetaI_deg)
        return false;
    return true;
}

double ReconfigurableIntelligentSurface::gain(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad)
{

    // TODO: remove this when the model with incidence angle comes in
//    phiI_rad = 0;
//    thetaI_rad = 0;

    if (canUseCache(phiI_rad, thetaI_rad)) {
        return cachedGain(phiR_rad, thetaR_rad);
    }

    CMatrix phase = new_cmatrix(Ts, Ps);

    Matrix B = coding;
    Matrix alpha = coding_alpha;
    for (int a = 0; a < M; a++) {
        for (int b = 0; b < N; b++) {
            Matrix a_0_5 = new_matrix(k_dG_sin_cos);
            Matrix b_0_5 = new_matrix(k_dG_sin_sin);
            gsl_matrix_scale(a_0_5, (a+1) - 0.5);
            gsl_matrix_scale(b_0_5, (b+1) - 0.5);
            gsl_matrix_add(a_0_5, b_0_5);
            gsl_matrix_add_constant(a_0_5, gsl_matrix_get(alpha, a, b));
            gsl_matrix_add_constant(a_0_5, gsl_matrix_get(B, a, b));
            CMatrix complex_matrix = new_cmatrix(a_0_5);
            gsl_matrix_complex_scale(complex_matrix, C_0_M1j);
            //            double phase = exp(-1j*(alpha(a,b)+PH(a,b) + k*dG*(a-0.5).*sin(theta*M_PI/180)'*cos(phi*M_PI/180) + k*dG*(b-0.5).*sin(theta*M_PI/180)'*sin(phi*M_PI/180)));%%
            exp(complex_matrix);
            gsl_matrix_complex_add(phase, complex_matrix);

            gsl_matrix_complex_free(complex_matrix);
            gsl_matrix_free(a_0_5);
            gsl_matrix_free(b_0_5);
        }
    }

    Matrix Fa = abs(phase);
    gsl_matrix_complex_free(phase);

    // F2 = Fa ^ 2
    matrix_copy(F2, Fa);
    gsl_matrix_mul_elements(F2, F2);

    // F2^T * sin(abs(theta) * pi / 180)
    Matrix F2_sin_abs_theta = new_matrix(F2->size2, sin_abs_theta->size2);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, F2, sin_abs_theta, 0, F2_sin_abs_theta);

    Matrix vm = new_matrix(1, 1);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, F2_sin_abs_theta, phi_Dt_Df, 0, vm);
    v = gsl_matrix_get(vm, 0, 0);

    gsl_matrix_free(vm);
    gsl_matrix_free(F2_sin_abs_theta);

    gsl_matrix_free(Fa);

    recomputeGainMap = false;
    cached_phiI_deg = RAD_TO_DEG_ROUND(phiI_rad);
    cached_thetaI_deg = RAD_TO_DEG_ROUND(thetaI_rad);

    return cachedGain(phiR_rad, thetaR_rad);
}

void ReconfigurableIntelligentSurface::getConfiguration(double& phiR_rad, double& thetaR_rad, double& phiI_rad, double& thetaI_rad)
{
    phiR_rad = configPhiR;
    phiI_rad = configPhiI;
    thetaR_rad = configThetaR;
    thetaI_rad = configThetaI;
}

Vector ReconfigurableIntelligentSurface::linspace(double min, double max, int n)
{
    Vector v = new_vector(n);
    for (int i = 0; i < n; i++)
        gsl_vector_set(v, i, min + ((double)i / (n - 1)) * (max - min));
    return v;
}

Vector ReconfigurableIntelligentSurface::new_vector(int n, double value)
{
    Vector v = gsl_vector_alloc(n);
    gsl_vector_set_all(v, value);
    return v;
}

Vector ReconfigurableIntelligentSurface::new_vector(Vector src)
{
    Vector v = gsl_vector_alloc(src->size);
    for (int i = 0; i < src->size; i++)
        gsl_vector_set(v, i, gsl_vector_get(src, i));
    return v;
}

Matrix ReconfigurableIntelligentSurface::matrix_from_vector(Vector v)
{
    Matrix m = new_matrix(v->size, 1);
    for (int i = 0; i < v->size; i++)
        gsl_matrix_set(m, i, 0, gsl_vector_get(v, i));
    return m;
}

Matrix ReconfigurableIntelligentSurface::new_matrix(int rows, int columns, double value)
{
    Matrix m = gsl_matrix_alloc(rows, columns);
    gsl_matrix_set_all(m, value);
    return m;
}

Matrix ReconfigurableIntelligentSurface::new_matrix(Matrix src)
{
    Matrix m = gsl_matrix_alloc(src->size1, src->size2);
    for (int r = 0; r < src->size1; r++)
        for (int c = 0; c < src->size2; c++)
            gsl_matrix_set(m, r, c, gsl_matrix_get(src, r, c));


    return m;
}

CMatrix ReconfigurableIntelligentSurface::new_cmatrix(int rows, int columns, gsl_complex value)
{
    CMatrix m = gsl_matrix_complex_alloc(rows, columns);
    gsl_matrix_complex_set_all(m, value);
    return m;
}

CMatrix ReconfigurableIntelligentSurface::new_cmatrix(Matrix src)
{
    CMatrix m = gsl_matrix_complex_alloc(src->size1, src->size2);
    matrix_to_cmatrix(m, src);
    return m;
}

void ReconfigurableIntelligentSurface::matrix_to_cmatrix(CMatrix dst, const Matrix src)
{
    for (int r = 0; r < src->size1; r++)
        for (int c = 0; c < src->size2; c++) {
            gsl_complex v;
            GSL_REAL(v) = gsl_matrix_get(src, r, c);
            GSL_IMAG(v) = 0;
            gsl_matrix_complex_set(dst, r, c, v);
        }
}

void ReconfigurableIntelligentSurface::matrix_copy(Matrix dst, const Matrix src)
{
    for (int r = 0; r < src->size1; r++)
        for (int c = 0; c < src->size2; c++) {
            gsl_matrix_set(dst, r, c, gsl_matrix_get(src, r, c));
        }
}

void ReconfigurableIntelligentSurface::apply_in_place(Matrix m, double (*f)(double))
{
    for (int r = 0; r < m->size1; r++)
        for (int c = 0; c < m->size2; c++)
            gsl_matrix_set(m, r, c, f(gsl_matrix_get(m, r, c)));


}

void ReconfigurableIntelligentSurface::abs(Vector v)
{
    for (int i = 0; i < v->size; i++)
        gsl_vector_set(v, i, std::abs(gsl_vector_get(v, i)));
}

void ReconfigurableIntelligentSurface::abs(Matrix m)
{
    for (int r = 0; r < m->size1; r++)
        for (int c = 0; c < m->size2; c++)
            gsl_matrix_set(m, r, c, std::abs(gsl_matrix_get(m, r, c)));


}

Matrix ReconfigurableIntelligentSurface::abs(CMatrix m)
{
    Matrix n = gsl_matrix_alloc(m->size1, m->size2);
    for (int r = 0; r < m->size1; r++)
        for (int c = 0; c < m->size2; c++)
            gsl_matrix_set(n, r, c, gsl_complex_abs(gsl_matrix_complex_get(m, r, c)));


    return n;
}

int ReconfigurableIntelligentSurface::min_pos(Vector v)
{
    if (v->size == 0)
        return -1;
    int pos = gsl_vector_min_index(v);
    return pos;
}

void ReconfigurableIntelligentSurface::exp(CMatrix m)
{
    for (int r = 0; r < m->size1; r++)
        for (int c = 0; c < m->size2; c++)
            gsl_matrix_complex_set(m, r, c, gsl_complex_exp(gsl_matrix_complex_get(m, r, c)));


}

void ReconfigurableIntelligentSurface::print_matrix(Matrix m) {
    for (int r = 0; r < m->size1; r++) {
        for (int c = 0; c < m->size2; c++) {
            printf("%+.4f ", gsl_matrix_get(m, r, c));
        }
        printf("\n");
    }
}

void ReconfigurableIntelligentSurface::print_matrix(CMatrix m) {
    for (int r = 0; r < m->size1; r++) {
        for (int c = 0; c < m->size2; c++) {
            gsl_complex v = gsl_matrix_complex_get(m, r, c);
            printf("%+.4f %+.4fi ", GSL_REAL(v), GSL_IMAG(v));
        }
        printf("\n");
    }
}

void ReconfigurableIntelligentSurface::print_vector(Vector v) {
    for (int c = 0; c < v->size; c++)
        printf("%+.4f ", gsl_vector_get(v, c));
    printf("\n");
}
