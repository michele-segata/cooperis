
//
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

#pragma once

#include "veins-rms/veins-rms.h"

#include "veins/base/utils/Coord.h"
#include "veins-rms/utility/ReconfigurableMetaSurface.h"
#include "veins-rms/utility/Utils.h"

namespace veins {


class VEINS_RMS_API RmsTest {
public:
    RmsTest();

    ~RmsTest(){};

    ReconfigurableMetaSurface* rms = nullptr;
    void writeGains(string filename, double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad);
    void writeAngles(string filename);
};
} // namespace veins
