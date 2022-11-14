#pragma once

#include <cmath>
using namespace std;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef gsl_vector* Vector;
typedef gsl_vector_complex* CVector;
typedef gsl_matrix* Matrix;
typedef gsl_matrix_complex* CMatrix;

#define M_PI_X_2 (2*M_PI)

//#define RIS_AZIMUTH_RIGHT 0
//#define RIS_AZIMUTH_CENTER (M_PI_2)
//#define RIS_AZIMUTH_LEFT (M_PI)
//#define RIS_ELEVATION_BOTTOM 0
//#define RIS_ELEVATION_MIDDLE (M_PI_2)
//#define RIS_ELEVATION_ABOVE (M_PI)

//#define RIS_AZIMUTH_RIGHT (-M_PI_2)
//#define RIS_AZIMUTH_CENTER 0
//#define RIS_AZIMUTH_LEFT M_PI_2
//#define RIS_ELEVATION_BOTTOM (-M_PI_2)
//#define RIS_ELEVATION_MIDDLE 0
//#define RIS_ELEVATION_ABOVE M_PI_2

//// in the math scripts the range goes from 0 to 180 both for azimuth and elevation
//// in the code we consider -90 (below and right) to 90 (above and left)
//#define REF_TO_MATH(x) (x+M_PI_2)

// theta goes from 0 to 90 both in the code and in the model
#define REF_TO_MATH_THETA(x) (x)
// phi goes from -180 to 180 in the code, but from 0 to 360 in the model
#define REF_TO_MATH_PHI(x) (x+M_PI)

#define RAD_TO_DEG(x) (x*180/M_PI)
#define RAD_TO_DEG_ROUND(x) (std::round(x*180/M_PI))
#define DEG_TO_RAD(x) (x*M_PI/180)

class ReconfigurableIntelligentSurface {

protected:
    double f;
    double n;

    double configPhiR = 0, configThetaR = 0, configPhiI = 0, configThetaI = 0;

    bool recomputeGainMap = true;

    bool canUseCache(double phiI_rad, double thetaI_rad) const;
    double cachedGain(double phiR_rad, double thetaR_rad) const;

private:
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
    double dG;
    int Ps;
    int Ts;
    Vector phi;
    Vector theta;
    double Dt;
    double Df;
    Matrix cos_phi;
    Matrix sin_theta;
    Matrix sin_phi;
    Matrix k_dG_sin_cos;
    Matrix k_dG_sin_sin;
    // 0 - 1j
    gsl_complex C_0_M1j;
    double wt;
    Matrix sin_abs_theta;
    Matrix phi_Dt_Df;

    // matrices holding the coding of the metasurface (state)
    Matrix coding;
    Matrix coding_alpha;
    // matrix and value caching the gains
    // gains are cached for a specific coding and a specific incidence angle
    Matrix F2;
    double v;
    int cached_phiI_deg;
    int cached_thetaI_deg;


public:
    ReconfigurableIntelligentSurface(double frequency, int n= 2);
    virtual ~ReconfigurableIntelligentSurface();

    /**
     * Reconfigures the meta surface to optimize the reflection for a specific pair of incidence and reflection angles
     *
     * @param phiR_rad the reflection angle phi (azimuth) in radians
     * @param thetaR_rad the reflection angle theta (elevation) in radians
     * @param phiI_rad the incidence angle phi (azimuth) in radians
     * @param thetaI_rad the incidence angle theta (elevation) in radians
     */
    void configureMetaSurface(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad);

    /**
     * Computes the gain of the antenna given a specific pair of incidence and reflection angles
     * @param phiR_rad the reflection angle phi (azimuth) in radians
     * @param thetaR_rad the reflection angle theta (elevation) in radians
     * @param phiI_rad the incidence angle phi (azimuth) in radians
     * @param thetaI_rad the incidence angle theta (elevation) in radians
     * @return linear gain of the RIS
     */
    double gain(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad);

    /**
     * Returns the incidence and reflection angles for which the metasurface has been configured
     */
    void getConfiguration(double& phiR_rad, double& thetaR_rad, double& phiI_rad, double& thetaI_rad);

    static Vector linspace(double start, double end, int n);
    static Vector new_vector(int n, double value = 0);
    static Vector new_vector(Vector src);
    static Matrix matrix_from_vector(Vector v);
    static Matrix new_matrix(int rows, int columns, double value = 0);
    static Matrix new_matrix(Matrix src);
    static CMatrix new_cmatrix(int rows, int columns, gsl_complex value = {0, 0});
    static CMatrix new_cmatrix(Matrix src);
    static void matrix_to_cmatrix(CMatrix dst, const Matrix src);
    static void matrix_copy(Matrix dst, const Matrix src);
    static void apply_in_place(Matrix v, double (*f)(double));
    static void abs(Vector v);
    static void abs(Matrix m);
    static Matrix abs(CMatrix m);
    static int min_pos(Vector v);
    static void exp(CMatrix m);
    static void print_matrix(Matrix m);
    static void print_matrix(CMatrix m);
    static void print_vector(Vector v);

};
