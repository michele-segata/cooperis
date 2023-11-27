//
// Copyright (C) 2018 Julien Jahneke <julien.jahneke@ccs-labs.org>
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

// #define STANDALONE 1

#include "catch_amalgamated.hpp"
#include <cmath>
#include "cooperis/utility/Utils.h"
#include "CsvReader.h"
#include "cooperis/utility/ReconfigurableIntelligentSurface.h"
#include <iostream>
#include <sstream>

#ifndef STANDALONE
using veins::Coord;
using veins::Angles;
#endif

#define EQUALS(x, y) (abs(x - y) < 1e-9)
#define EQUALS17(x, y) (abs(x - y) < 1e-17)
#define COORD_EQUALS(a, b) (EQUALS(a.x, b.x) && EQUALS(a.y, b.y) && EQUALS(a.z, b.z))
#define PI2 6.28

bool matrix_equals(const Matrix& m, const CsvReader& csv, double phiR) {
    if (m->size1 != csv.getRows())
        return false;
    if (m->size2 != csv.getColumns())
        return false;

    for (int r = 0; r < m->size1; r++)
        for (int c = 0; c < m->size2; c++)
            if (!EQUALS(gsl_matrix_get(m, r, c), csv.get(r, c))) {
                std::cout << "Test for phiR = " <<  phiR << ": Elements in position " << r << "," << c << " do not match: " << gsl_matrix_get(m, r, c) << " vs " << csv.get(r, c) << "\n";

                std::cout << "C++ Matrix:\n";
                ReconfigurableIntelligentSurface::print_matrix(m);

                std::cout << "\nMatlab Matrix:\n";
                csv.print();

                return false;
            }

    return true;
}

TEST_CASE("Modulo operations") {
    REQUIRE(EQUALS(ReconfigurableIntelligentSurface::nmod(7, PI2), 0.72));
    REQUIRE(EQUALS(ReconfigurableIntelligentSurface::nmod(14, PI2), 1.44));
    REQUIRE(EQUALS(ReconfigurableIntelligentSurface::nmod(-5, PI2), 1.28));
    REQUIRE(EQUALS(ReconfigurableIntelligentSurface::nmod(-7, PI2), 5.56));
    REQUIRE(EQUALS(ReconfigurableIntelligentSurface::nmod(-10.34, PI2), 2.22));
    REQUIRE(EQUALS(ReconfigurableIntelligentSurface::nmod(-6.2831853071795862, M_PI_X_2), 0));
    REQUIRE(EQUALS(ReconfigurableIntelligentSurface::nmod(-18.849555921538762, M_PI_X_2), 0));
}

TEST_CASE("Rotations") {
    Coord x(1, 0, 0);
    Coord y(0, 1, 0);
    Coord z(0, 0, 1);
    Coord nx(-1, 0, 0);
    Coord ny(0, -1, 0);
    Coord nz(0, 0, -1);
    REQUIRE(COORD_EQUALS(rotate(x, y, DEG_TO_RAD(90)), z));
    REQUIRE(COORD_EQUALS(rotate(y, z, DEG_TO_RAD(90)), x));
    REQUIRE(COORD_EQUALS(rotate(z, x, DEG_TO_RAD(90)), y));
    REQUIRE(COORD_EQUALS(rotate(x, z, DEG_TO_RAD(90)), ny));
    REQUIRE(COORD_EQUALS(rotate(y, x, DEG_TO_RAD(90)), nz));
    REQUIRE(COORD_EQUALS(rotate(z, y, DEG_TO_RAD(90)), nx));
}

TEST_CASE("Angle distance") {
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(0), DEG_TO_RAD(360))) == 0);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(360), DEG_TO_RAD(0))) == 0);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(0), DEG_TO_RAD(720))) == 0);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(720), DEG_TO_RAD(0))) == 0);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(-180), DEG_TO_RAD(179))) == -1);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(-181), DEG_TO_RAD(179))) == 0);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(179), DEG_TO_RAD(-180))) == 1);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(20), DEG_TO_RAD(-180))) == 160);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(-20), DEG_TO_RAD(-180))) == -160);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(-180), DEG_TO_RAD(20))) == -160);
    REQUIRE(RAD_TO_DEG_ROUND(ReconfigurableIntelligentSurface::angle_distance(DEG_TO_RAD(-180), DEG_TO_RAD(-20))) == 160);
}

TEST_CASE("Nearest phase") {
    Vector phases = gsl_vector_alloc(4);
    gsl_vector_set(phases, 0, 0);
    gsl_vector_set(phases, 1, M_PI_2);
    gsl_vector_set(phases, 2, M_PI);
    gsl_vector_set(phases, 3, M_PI / 2 * 3);
    REQUIRE(ReconfigurableIntelligentSurface::nearest_angle_pos(phases, 6.21) == 0);
    REQUIRE(ReconfigurableIntelligentSurface::nearest_angle_pos(phases, 6.29) == 0);
    gsl_vector_free(phases);
}

TEST_CASE("Metasurface phases") {
    ReconfigurableIntelligentSurface ris(0, 25e9, 4, 3, 5);
    CsvReader phases;

    int phiRs[] = {-180, -135, -90, -45, 0, 45, 90, 135, 180};
    for (auto &phiR : phiRs) {
        SECTION("Result for phiR = " + std::to_string(phiR))
        {
            std::stringstream filename;
            filename << "matlab/phases_phiR_";
            filename << phiR;
            filename << "_thetaR_45_phiI_0_thetaI_0_phiTX_0_thetaTX_0_n_4_pl_3_nl_5.csv";
            phases.read(filename.str());
            ris.configureMetaSurface(DEG_TO_RAD(phiR), DEG_TO_RAD(45), 0, 0);
            REQUIRE(matrix_equals(ris.getPhases(), phases, phiR));
        }
    }
}

TEST_CASE("Metasurface gains") {
    ReconfigurableIntelligentSurface ris(0, 25e9, 4, 3, 5);
    CsvReader gains;

    int phiRs[] = {-180, -135, -90, -45, 0, 45, 90, 135, 180};
    for (auto &phiR : phiRs) {
        SECTION("Result for phiR = " + std::to_string(phiR))
        {
            std::stringstream filename;
            filename << "matlab/gain_phiR_";
            filename << phiR;
            filename << "_thetaR_45_phiI_0_thetaI_0_phiTX_0_thetaTX_0_n_4_pl_3_nl_5.csv";
            gains.read(filename.str());
            ris.configureMetaSurface(DEG_TO_RAD(phiR), DEG_TO_RAD(45), 0, 0);
            double p_tot;
            Matrix gain = ReconfigurableIntelligentSurface::new_matrix(ris.getGains(p_tot, 0, 0));
            gsl_matrix_scale(gain, M_PI_X_2 / p_tot);
            REQUIRE(matrix_equals(gain, gains, phiR));
            gsl_matrix_free(gain);
        }
    }
}

TEST_CASE("Coordinates system") {

    Coord ris_v1(1, 0, 0);
    Coord ris_v2(0, 0, -1);
    Coord ris_pos(0, 0, 0);
    // if x axis is v1 and negative z axis is v2,
    // then vn = v1 x v2 = y axis
    Coord ris_vn = cross(ris_v1, ris_v2);

    SECTION("Plane normal") {
        REQUIRE((ris_vn.x == 0 && ris_vn.y == 1 && ris_vn.z == 0));
    }

    Angles a{};
    Coord node;

    SECTION("Coordinates test: surface normal") {
        // surface normal
        node = Coord(0, 1, 0);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        REQUIRE(EQUALS(a.theta, DEG_TO_RAD(0)));
    }

    SECTION("Coordinates test: north (phi = 0)") {
        // above should be north
        node = Coord(0, 0, 1);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        REQUIRE((EQUALS(a.theta, DEG_TO_RAD(90)) && EQUALS(a.phi, DEG_TO_RAD(0))));
    }

    SECTION("Coordinates test: east (phi = -90)") {
        node = Coord(1, 0, 0);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        REQUIRE((EQUALS(a.theta, DEG_TO_RAD(90)) && EQUALS(a.phi, DEG_TO_RAD(-90))));
    }

    SECTION("Coordinates test: west (phi = 90)") {
        node = Coord(-1, 0, 0);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        REQUIRE((EQUALS(a.theta, DEG_TO_RAD(90)) && EQUALS(a.phi, DEG_TO_RAD(90))));
    }

    SECTION("Coordinates test: south (phi = 180)") {
        node = Coord(-0.01, 0, -1);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        node = Coord(0.01, 0, -1);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        node = Coord(0, 0, -1);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        REQUIRE((EQUALS(a.theta, DEG_TO_RAD(90)) &&
                 (EQUALS(a.phi, DEG_TO_RAD(180)) || EQUALS(a.phi, DEG_TO_RAD(-180)))));
    }

    SECTION("Coordinate test: diagonal of a cube, direction north-west") {
        node = Coord(-1, 1, 1);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        REQUIRE((EQUALS(a.theta, DEG_TO_RAD(54.7356103)) && EQUALS(a.phi, DEG_TO_RAD(45))));
    }

    SECTION("Coordinate test: diagonal of a cube, direction north-east") {
        node = Coord(1, 1, 1);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        REQUIRE((EQUALS(a.theta, DEG_TO_RAD(54.7356103)) && EQUALS(a.phi, DEG_TO_RAD(-45))));
    }

    SECTION("Coordinate test: diagonal of a cube, direction south-west") {
        node = Coord(-1, 1, -1);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        REQUIRE((EQUALS(a.theta, DEG_TO_RAD(54.7356103)) && EQUALS(a.phi, DEG_TO_RAD(135))));
    }

    SECTION("Coordinate test: diagonal of a cube, direction south-east") {
        node = Coord(1, 1, -1);
        a = spherical_angles(ris_v1, ris_v2, ris_vn, ris_pos, node);
        REQUIRE((EQUALS(a.theta, DEG_TO_RAD(54.7356103)) && EQUALS(a.phi, DEG_TO_RAD(-135))));
    }

    SECTION("Beam direction: surface normal") {
        node = spherical_point_beam(ris_v1, ris_v2, ris_vn, ris_pos, DEG_TO_RAD(0), DEG_TO_RAD(0));
        REQUIRE((EQUALS(node.x, 0) && EQUALS(node.y, 1) && EQUALS(node.z, 0)));
    }

    SECTION("Beam direction: north (phi = 0)") {
        // above should be north
        node = spherical_point_beam(ris_v1, ris_v2, ris_vn, ris_pos, DEG_TO_RAD(0), DEG_TO_RAD(90));
        REQUIRE((EQUALS(node.x, 0) && EQUALS(node.y, 0) && EQUALS(node.z, 1)));
    }

    SECTION("Beam direction: east (phi = -90)") {
        node = spherical_point_beam(ris_v1, ris_v2, ris_vn, ris_pos, DEG_TO_RAD(-90), DEG_TO_RAD(90));
        REQUIRE((EQUALS(node.x, 1) && EQUALS(node.y, 0) && EQUALS(node.z, 0)));
    }

    SECTION("Beam direction: west (phi = 90)") {
        node = spherical_point_beam(ris_v1, ris_v2, ris_vn, ris_pos, DEG_TO_RAD(90), DEG_TO_RAD(90));
        REQUIRE((EQUALS(node.x, -1) && EQUALS(node.y, 0) && EQUALS(node.z, 0)));
    }

    SECTION("Beam direction: south (phi = 180)") {
        node = spherical_point_beam(ris_v1, ris_v2, ris_vn, ris_pos, DEG_TO_RAD(180), DEG_TO_RAD(90));
        REQUIRE((EQUALS(node.x, 0) && EQUALS(node.y, 0) && EQUALS(node.z, -1)));
    }

}

TEST_CASE("Spherical elements") {

    SECTION("Spherical element: theta = 0") {
        REQUIRE(EQUALS17(ReconfigurableIntelligentSurface::spherical_element(0, M_PI / 180, M_PI / 180),
                         0.000000664567899281));
    }

    SECTION("Spherical element: theta = 45") {
        REQUIRE(EQUALS17(ReconfigurableIntelligentSurface::spherical_element(M_PI_2 / 2, M_PI / 180, M_PI / 180),
                         0.000215394309305328));
    }

    SECTION("Spherical element: theta = 90") {
        REQUIRE(EQUALS17(ReconfigurableIntelligentSurface::spherical_element(M_PI_2, M_PI / 180, M_PI / 180),
                         0.000152306776738789));
    }

}

TEST_CASE("Azimuth angles to index conversion") {
    double d_phi[] = {1, 0.5, 2};
    vector<double> angles;
    for (double delta : d_phi) {
        SECTION("Azimuth with d_phi = " + std::to_string(delta) + " degree") {
            double min_angle = -180 + delta;
            double max_angle = 180;
            double angle_range = max_angle - min_angle;
            for (double a = min_angle; a <= max_angle; a += delta) {
                angles.push_back(DEG_TO_RAD(a));
            }
            vector<size_t> indexes;
            bool diff_one = true;
            for (double a : angles) {
                indexes.push_back(ReconfigurableIntelligentSurface::angle_to_index(a, DEG_TO_RAD(angle_range), DEG_TO_RAD(min_angle), angles.size()));
                if (indexes.size() > 1) {
                    // if we scan all angles, the index difference between successive angles must be one
                    if (indexes[indexes.size() - 1] - indexes[indexes.size() - 2] != 1) {
                        std::cout << "Non consecutive indexes for angle " << RAD_TO_DEG(a) << ": index for previous="
                                  << indexes[indexes.size() - 2] << " current index=" << indexes[indexes.size() - 1]
                                  << "\n";
                        diff_one = false;
                    }
                }
            }
            REQUIRE(diff_one);
            REQUIRE(indexes[0] == 0);
            REQUIRE(indexes[indexes.size() - 1] == indexes.size() - 1);
        }
    }
}

TEST_CASE("Elevation angles to index conversion") {
    double d_theta[] = {1, 0.5, 2};
    vector<double> angles;
    for (double delta : d_theta) {
        SECTION("Elevation with d_theta = " + std::to_string(delta) + " degree") {
            double min_angle = 0;
            double max_angle = 90;
            double angle_range = max_angle - min_angle;
            for (double a = min_angle; a <= max_angle; a += delta) {
                angles.push_back(DEG_TO_RAD(a));
            }
            vector<size_t> indexes;
            bool diff_one = true;
            for (double a : angles) {
                indexes.push_back(ReconfigurableIntelligentSurface::angle_to_index(a, DEG_TO_RAD(angle_range), DEG_TO_RAD(min_angle), angles.size()));
                if (indexes.size() > 1) {
                    // if we scan all angles, the index difference between successive angles must be one
                    if (indexes[indexes.size() - 1] - indexes[indexes.size() - 2] != 1) {
                        std::cout << "Non consecutive indexes for angle " << RAD_TO_DEG(a) << ": index for previous="
                                  << indexes[indexes.size() - 2] << " current index=" << indexes[indexes.size() - 1]
                                  << "\n";
                        diff_one = false;
                    }
                }
            }
            REQUIRE(diff_one);
            REQUIRE(indexes[0] == 0);
            REQUIRE(indexes[indexes.size() - 1] == indexes.size() - 1);
        }
    }
}

TEST_CASE("Azimuth angle to nearest index conversion") {

    double d_phi[] = {1, 0.5, 2};
    vector<double> angles;
    for (double delta : d_phi) {
        SECTION("Azimuth with d_phi = " + std::to_string(delta) + " degree") {
            double min_angle = -180 + delta;
            double max_angle = 180;
            double angle_range = max_angle - min_angle;
            for (double a = min_angle; a <= max_angle; a += delta) {
                angles.push_back(DEG_TO_RAD(a + delta/2 + delta/4));
            }
            vector<size_t> indexes;
            bool diff_one = true;
            // skip last angle which is going outside
            for (size_t i = 0; i < angles.size() - 1; i++) {
                double a = angles[i];
                indexes.push_back(ReconfigurableIntelligentSurface::angle_to_index(a, DEG_TO_RAD(angle_range), DEG_TO_RAD(min_angle), angles.size()));
                if (indexes.size() > 1) {
                    // if we scan all angles, the index difference between successive angles must be one
                    if (indexes[indexes.size() - 1] - indexes[indexes.size() - 2] != 1) {
                        std::cout << "Non consecutive indexes for angle " << RAD_TO_DEG(a) << ": index for previous="
                                  << indexes[indexes.size() - 2] << " angle=" << RAD_TO_DEG(angles[i-1]) << " current index=" << indexes[indexes.size() - 1]
                                  << "\n";
                        diff_one = false;
                    }
                }
            }
            REQUIRE(diff_one);
            // given that the first angle is the minimum angle plus more than half of the delta, the procedure should give 1 as first index
            REQUIRE(indexes[0] == 1);
            REQUIRE(indexes[indexes.size() - 1] == indexes.size());
        }
    }

}

TEST_CASE("Elevation angle to nearest index conversion") {
    double d_theta[] = {1, 0.5, 2};
    vector<double> angles;
    for (double delta : d_theta) {
        SECTION("Elevation with d_theta = " + std::to_string(delta) + " degree") {
            double min_angle = 0;
            double max_angle = 90;
            double angle_range = max_angle - min_angle;
            for (double a = min_angle; a <= max_angle; a += delta) {
                angles.push_back(DEG_TO_RAD(a + delta/2 + delta/4));
            }
            vector<size_t> indexes;
            bool diff_one = true;
            for (double a : angles) {
                indexes.push_back(ReconfigurableIntelligentSurface::angle_to_index(a, DEG_TO_RAD(angle_range), DEG_TO_RAD(min_angle), angles.size()));
                if (indexes.size() > 1) {
                    // if we scan all angles, the index difference between successive angles must be one
                    if (indexes[indexes.size() - 1] - indexes[indexes.size() - 2] != 1) {
                        std::cout << "Non consecutive indexes for angle " << RAD_TO_DEG(a) << ": index for previous="
                                  << indexes[indexes.size() - 2] << " current index=" << indexes[indexes.size() - 1]
                                  << "\n";
                        diff_one = false;
                    }
                }
            }
            REQUIRE(diff_one);
            REQUIRE(indexes[0] == 1);
            REQUIRE(indexes[indexes.size() - 1] == indexes.size());
        }
    }
}

TEST_CASE("Generation of gains") {
    SECTION("Output gains for plotting") {
        ReconfigurableIntelligentSurface ris(0, 25e9, 4, 3, 5);

        int phiRs[] = {-180, -135, -90, -45, 0, 45, 90, 135, 180};
        for (auto &phiR : phiRs) {
            SECTION("Result for phiR = " + std::to_string(phiR))
            {
                ris.configureMetaSurface(DEG_TO_RAD(phiR), DEG_TO_RAD(45), 0, 0);
                ris.writeGains("plotgains", 0, 0);
            }
        }
        REQUIRE(true);
    }
}

// https://stackoverflow.com/a/6098417
template <typename Iter>
std::string join(Iter begin, Iter end, std::string const &separator)
{
    std::ostringstream result;
    if (begin != end)
        result << *begin++;
    while (begin != end)
        result << separator << *begin++;
    return result.str();
}

TEST_CASE("Metasurface phases for combined configs (average strategy)") {
    ReconfigurableIntelligentSurface ris(0, 25e9, 4, 3, 5);
    CsvReader phases;

    vector<vector<int>> phiRs = {{-45, 45}, {-45, 45, 135}, {-45, 45, 130}};
    vector<vector<int>> thetaRs = {{45, 45}, {45, 45, 45}, {45, 45, 45}};
    for (int i = 0; i < phiRs.size(); i++) {
        auto phiR = phiRs[i];
        auto thetaR = thetaRs[i];
        SECTION("Result for phiR = " + join(phiR.begin(), phiR.end(), ", "))
        {
            std::stringstream filename;
            filename << "matlab/phases_phiR_";
            filename << join(phiR.begin(), phiR.end(), "_");
            filename << "_thetaR_";
            filename << join(thetaR.begin(), thetaR.end(), "_");
            filename << "_phiI_0_thetaI_0_phiTX_0_thetaTX_0_n_4_pl_3_nl_5.csv";
            phases.read(filename.str());

            VMatrix configs;
            for (int phi : phiR) {
                configs.push_back(ris.computePhases(DEG_TO_RAD(phi), DEG_TO_RAD(45), 0, 0));
            }
            Matrix combined = ris.combineConfigurations(configs, false);
            REQUIRE(matrix_equals(combined, phases, phiR[0]));
            gsl_matrix_free(combined);
            for (auto config : configs) {
                gsl_matrix_free(config);
            }
        }
    }
}

TEST_CASE("Combination of configs") {
    ReconfigurableIntelligentSurface ris(0, 25e9, 4, 3, 5);
    vector<vector<int>> phiRs = {{-45, 45}, {-45, 45, 135}, {-45, 45, 130}};
    bool random[] = {true, false};
    for (vector<int> phiR : phiRs) {
        for (bool rnd : random) {
            SECTION("Combining configs for phiRs = [" + join(phiR.begin(), phiR.end(), ", ") + "], strategy = " + (rnd ? "random" : "average")) {
                VMatrix configs;
                for (int phi : phiR) {
                    configs.push_back(ris.computePhases(DEG_TO_RAD(phi), DEG_TO_RAD(45), 0, 0));
                }
                Matrix combined = ris.combineConfigurations(configs, rnd);
                ris.applyConfiguration(combined);
                ris.writeGains("combogains_" + join(phiR.begin(), phiR.end(), "-") + "_" + (rnd ? "random" : "average"), 0, 0);
                gsl_matrix_free(combined);
                for (int i = 0; i < phiR.size(); i++) {
                    gsl_matrix_free(configs[i]);
                }
                REQUIRE(true);
            }
        }
    }
}
