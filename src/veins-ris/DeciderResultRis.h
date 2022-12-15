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
// Base files derived from Veins VLC by Agon Memedi and contributors

#pragma once

#include "veins/veins.h"

#include "veins/modules/phy/DeciderResult80211.h"

namespace veins {

class VEINS_API DeciderResultRis : public DeciderResult80211 {
protected:
    /** @brief Stores the gain applied by the ris TODO: transform this into an array */
    std::vector<double> gains;
    double phiR, thetaR, phiI, thetaI;
    bool reflected;

public:
    /**
     * @brief Initialises with the passed values.
     *
     * "bitrate" defines the bit-rate of the transmission of the packet.
     */
    DeciderResultRis(bool isCorrect, double bitrate, double snr, double recvPower_dBm = 0, bool collision = false, std::vector<double> gains = std::vector<double>(), double phiR= 0, double thetaR= 0, double phiI= 0, double thetaI= 0, bool reflected= false)
        : DeciderResult80211(isCorrect, bitrate, snr, recvPower_dBm, collision)
        , gains(gains)
        , phiR(phiR)
        , thetaR(thetaR)
        , phiI(phiI)
        , thetaI(thetaI)
        , reflected(reflected)
    {
    }

    std::vector<double> getGains() const
    {
        return gains;
    }

    double getPhiR() const
    {
        return phiR;
    }
    double getThetaR() const
    {
        return thetaR;
    }
    double getPhiI() const
    {
        return phiI;
    }
    double getThetaI() const
    {
        return thetaI;
    }
    double getReflected() const
    {
        return reflected;
    }
};

} // namespace veins
