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

#pragma once

#include <cmath>
#include <random>
using namespace std;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef gsl_vector* Vector;
typedef gsl_vector_complex* CVector;
typedef gsl_matrix* Matrix;
typedef gsl_matrix_complex* CMatrix;
typedef std::vector<Matrix> VMatrix;

#define M_PI_X_2 (2*M_PI)

#define RAD_TO_DEG(x) ((x)*180/M_PI)
#define RAD_TO_DEG_ROUND(x) (std::round((x)*180/M_PI))
#define DEG_TO_RAD(x) ((x)*M_PI/180)

#define KEEP_SAME_ANGLE 10000

#define ANGLE_EQUALS(x, y) (std::abs((x)-(y)) < 1e-9)

class ReconfigurableIntelligentSurface {

protected:
    double efficiency;

    double configPhiR = 0, configThetaR = 0, configPhiI = 0, configThetaI = 0;

    bool recomputeGainMap = true;

    bool canUseCache(double phiI_rad, double thetaI_rad) const;
    double cachedGain(double phiR_rad, double thetaR_rad) const;

private:
    // stored parameters
    // total number of states (number of available phases)
    int N_s;
    // size of the side of the RIS in multiples of lambda
    int N_lambda;
    // number of reflecting elements per lambda
    int rho_lambda;
    // set of pre-initialized variables that can be used multiple times
    Vector E;
    double c;
    double lambda;
    double k;
    double du;
    double DN;
    double DM;
    int M, N;
    double du_k;
    int Ps;
    int Ts;
    Vector phi;
    Vector theta;
    double d_theta;
    double d_phi;
    double min_phi, max_phi;
    double min_theta, max_theta;
    double phi_range, theta_range;
    Vector spherical_elements;
    Matrix cos_phi;
    Matrix sin_theta;
    Matrix sin_phi;
    Matrix k_du_sin_cos;
    Matrix k_du_sin_sin;
    // 0 - 1j
    gsl_complex C_0_M1j{};

    // matrices holding the coding of the metasurface (state)
    Matrix coding;
    // matrix and value caching the gains
    // gains are cached for a specific coding and a specific incidence angle
    Matrix P;
    double p_tot{};
    double cached_phiTX{};
    double cached_thetaTX{};

    std::mt19937_64 rng{};

    unsigned int n_max_threads;

    static void gain_compute_phase_CPU_routine(void* thread_args);

    void gain_compute_phase_CPU(CMatrix phase, double phiRX_rad, double thetaRX_rad, double phiTX_rad, double thetaTX_rad);

    struct thread_gain_args;

public:
    /**
     * @param frequency frequency in Hz at which the surface is operating
     * @param n number of states (discretization of phases)
     * @param cellsPerLambda number of unit cells per lambda (wavelength)
     * @param lambdaSize length of the side of the surface in units of lambda
     * @param efficiency reflective efficiency of the surface
     */
    ReconfigurableIntelligentSurface(int seed, double frequency, int n= 2, int cellsPerLambda= 3, int lambdaSize= 5, double efficiency= 1.0, double d_theta= M_PI/180, double d_phi= M_PI/180);
    virtual ~ReconfigurableIntelligentSurface();

    /**
     * Reconfigures the meta surface to optimize the reflection for a specific pair of incidence and reflection angles
     * This procedure enables caching
     *
     * @param phiR_rad the reflection angle phi (azimuth) in radians
     * @param thetaR_rad the reflection angle theta (elevation) in radians
     * @param phiI_rad the incidence angle phi (azimuth) in radians
     * @param thetaI_rad the incidence angle theta (elevation) in radians
     */
    void configureMetaSurface(double phiR_rad= KEEP_SAME_ANGLE, double thetaR_rad= KEEP_SAME_ANGLE, double phiI_rad= KEEP_SAME_ANGLE, double thetaI_rad= KEEP_SAME_ANGLE);

    /**
     * Computes and returns the phases to be applied to optimize the reflection for a specific pair of incidence and reflection angles
     * @param phases matrix where to store the computed phases. This must be pre-allocated
     * @param phiR_rad the reflection angle phi (azimuth) in radians
     * @param thetaR_rad the reflection angle theta (elevation) in radians
     * @param phiI_rad the incidence angle phi (azimuth) in radians
     * @param thetaI_rad the incidence angle theta (elevation) in radians
     */
    void computePhases(Matrix phases, double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad);

    /**
     * Computes and returns the phases to be applied to optimize the reflection for a specific pair of incidence and reflection angles
     * @param phiR_rad the reflection angle phi (azimuth) in radians
     * @param thetaR_rad the reflection angle theta (elevation) in radians
     * @param phiI_rad the incidence angle phi (azimuth) in radians
     * @param thetaI_rad the incidence angle theta (elevation) in radians
     * @return matrix of computed phases. This is allocated by the method and should be freed by the user with gsl_matrix_free
     */
    Matrix computePhases(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad);

    /**
     * Combines a set of configurations into a single one to perform, for example, beam splitting
     * @param configs vector of phase configurations to combine
     * @param combined matrix where to store the final combined configuration. This must be pre-allocated
     * @param random whether to use the random or the average strategy. Combining is either done using phases which occurs most,
     * choosing a random one in case of phases with the same number of occurrences, or simply averaging all the phases (baseline)
     */
    void combineConfigurations(VMatrix configs, Matrix combined, bool random);

    /**
     * Combines a set of configurations into a single one to perform, for example, beam splitting
     * @param configs vector of phase configurations to combine
     * @param random whether to use the random or the average strategy. Combining is either done using phases which occurs most,
     * choosing a random one in case of phases with the same number of occurrences, or simply averaging all the phases (baseline)
     * @return matrix with final combined configuration. This is allocated by the method and should be freed by the user with gsl_matrix_free
     */
    Matrix combineConfigurations(VMatrix configs, bool random);

    /**
     * Applies a configuration to the elements of the RIS. This could be used, for example, when combining multiple configurations
     * Applying this configuration disables caching exactly for that reason (caching works for a single incident-reflected beam config)
     * @param config matrix of phases to be applied
     */
    void applyConfiguration(Matrix config);

    /**
     * Computes the gain of the antenna given a specific pair of incidence and reflection angles
     * @param phiRX_rad the reflection angle phi (azimuth) in radians
     * @param thetaRX_rad the reflection angle theta (elevation) in radians
     * @param phiTX_rad the incidence angle phi (azimuth) in radians
     * @param thetaTX_rad the incidence angle theta (elevation) in radians
     * @return linear gain of the RIS
     */
    double gain(double phiRX_rad, double thetaRX_rad, double phiTX_rad, double thetaTX_rad);

    /**
     * Returns the incidence and reflection angles for which the metasurface has been configured
     */
    void getConfiguration(double& phiR_rad, double& thetaR_rad, double& phiI_rad, double& thetaI_rad) const;

    /**
     * Returns the matrix of phases assigned to the elements after running the configuration
     */
    const Matrix& getPhases() const;

    /**
     * Returns the matrix of gains after running the configuration
     * @param p_tot total power of the surface. The actual gain is P * 2 * pi / p_tot
     * @param phiTX_rad azimuth of the transmitter
     * @param thetaTX_rad elevation of the transmitter
     */
    const Matrix& getGains(double& p_tot, double phiTX_rad, double thetaTX_rad);

    /**
     * Returns the area of the spherical element at elevation theta, given the angular resolutions d_theta and d_phi
     * The element is considered to be centered in theta, +- d_theta/2
     * @param theta angle theta in radians
     * @param d_theta theta resolution in radians
     * @param d_phi phi resolution in radians
     */
    static double spherical_element(double theta, double d_theta, double d_phi);

    static double dot_product(Vector a, Vector b);

    void writeGains(string filename, double phiTX_rad, double thetaTX_rad);

    /**
     * Computes the modulo operation with negative dividend. The operations is the same as making the dividend positive
     * first and then computing the normal modulo operation. E.g., -7 mod 6 = (-7 + 6 + 6) mod 6 = 1
     * @param x dividend
     * @param y divisor
     * @return the modulo between x and y
     */
    static double nmod(double x, double y);

    /**
     * Performs a binary search for double values returning the position in the array if found or -1
     */
    static int binary_search(Vector v, double value);

    inline double fix_azimuth(double phi) const;

    /**
     * Converts an angle into an integer index for the matrix of gains. Sample parameters for azimuth:
     * if d_phi = 0.5 degrees then min = -179.5, max=180 (range = 359.5 degrees), n_angles = 720
     * @param angle_rad angle to convert
     * @param angle_range_rad range of angles in the matrix (max - min)
     * @param min_angle_rad minimum angle in the matrix
     * @param n_angles number of angles
     * @return the index for the given angle inside the matrix
     */
    static size_t angle_to_index(double angle_rad, double angle_range_rad, double min_angle_rad, size_t n_angles);

    static Vector linspace(double start, double end, int n);
    static Vector new_vector(size_t n, double value = 0);
    static Vector new_vector(Vector src);
    static Matrix matrix_from_vector(Vector v);
    static Matrix new_matrix(size_t rows, size_t columns, double value = 0);
    static Matrix new_matrix(Matrix src);
    static CMatrix new_cmatrix(size_t rows, size_t columns, gsl_complex value = {0, 0});
    static CMatrix new_cmatrix(Matrix src);
    static void matrix_to_cmatrix(CMatrix dst, const Matrix src);
    static void matrix_copy(Matrix dst, const Matrix src);
    static void apply_in_place(Matrix v, double (*f)(double));
    static void abs(Vector v);
    static void abs(Matrix m);
    static Matrix abs(CMatrix m);
    static size_t min_pos(Vector v);
    static double angle_distance(double a1, double a2);
    static size_t nearest_angle_pos(Vector phases, double phase);
    static void exp(CMatrix m);
    static double matrix_sum(Matrix m);
    static Vector matrix_sum(Matrix m, bool by_row);
    static void print_matrix(Matrix m);
    static void print_matrix(CMatrix m);
    static void print_vector(Vector v);

};
