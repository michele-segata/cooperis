// Copyright (C) 2018-2018 Julien Jahneke <julien.jahneke@ccs-labs.org>
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

#include "RisTest.h"
#include <fstream>

using namespace veins;

RisTest::RisTest()
{
    ris = new ReconfigurableMetaSurface(25e9);
}

void RisTest::writeGains(string filename, double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad)
{
    ris->configureMetaSurface(phiR_rad, thetaR_rad, phiI_rad, thetaI_rad);
    ofstream f;
    f.open(filename.c_str(), ios::out);
    f << "phi,theta,gain\n";
    for (int phi = -180; phi <= 180; phi++) {
        for (int theta = 0; theta <= 90; theta++) {
            if (theta == 10 && phi == -35) {
                printf("\n");
            }
            double phi_rad = DEG_TO_RAD(phi);
            double theta_rad = DEG_TO_RAD(theta);
            double gain = ris->gain(phi_rad, theta_rad, 0, 0);
            f << phi << "," << theta << "," << gain << "\n";
        }
    }
    f.close();
}

void RisTest::writeAngles(string filename)
{

    Coord risPos = {500, 20, 10};
    Coord ris_v1 = {0, 0, 1};
    Coord ris_v2 = {-1, 0, 0};
    Coord ris_vn = cross(ris_v2, ris_v1);
    Coord vPos = {0, 26.6, 1.895};

    double v = 50/3.6;
    Angles angles;

    ofstream f;
    f.open(filename.c_str(), ios::out);
    f << "t,x,phi,theta\n";
    for (double t = 0; t < 60; t += 0.1) {
        angles = spherical_angles(ris_v1, ris_v2, ris_vn, vPos, vPos);
        f << t << "," << vPos.x << "," << angles.phi << "," << angles.theta << "\n";
        vPos.x += v * 0.1;
    }
    f.close();

}
