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

#include "catch2/catch.hpp"
#include "veins/base/utils/Coord.h"
#include "RisTest.h"

#define PAIR(x,y) std::pair<int,int>(x,y)

using namespace veins;
SCENARIO("RisTestCases", "[risTestCases]")
{

    auto rist = new RisTest();

    GIVEN("Test gain")
    {
        WHEN("Configuring RIS")
        {
            rist->ris->gain(DEG_TO_RAD(-35), DEG_TO_RAD(10), 0, 0);
            THEN("True, False, False")
            {
                REQUIRE(true);
            }
        }
    }
    GIVEN("Computation of gains")
    {
        WHEN("Configuring RIS")
        {
            std::vector<std::pair<int, int>> angles;
            angles.push_back(PAIR(-80, 80));
            angles.push_back(PAIR(-30, 45));
            angles.push_back(PAIR(-45, 45));
            angles.push_back(PAIR(10, 10));
            angles.push_back(PAIR(75, 10));

            for (auto a : angles) {
                std::stringstream ss;
                ss << "gains_" << a.first << "_" << a.second << ".csv";
                rist->writeGains(ss.str(), DEG_TO_RAD(a.first), DEG_TO_RAD(a.second), 0, 0);
            }
            THEN("True, False, False")
            {
                REQUIRE(true);
            }
        }
    }
    GIVEN("Computation of angles")
    {
        WHEN("Car moving")
        {
            rist->writeAngles("angles.csv");
            THEN("True, False, False")
            {
                REQUIRE(true);
            }
        }

    }
}
