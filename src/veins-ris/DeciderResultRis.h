//
// Copyright (C) 2007 Technische Universitaet Berlin (TUB), Germany, Telecommunication Networks Group
// Copyright (C) 2007 Technische Universiteit Delft (TUD), Netherlands
// Copyright (C) 2007 Universitaet Paderborn (UPB), Germany
// Copyright (C) 2014 Michele Segata <segata@ccs-labs.org>
//
// Documentation for these modules is at http://veins.car2x.org/
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

/*
 * DeciderResult80211.h
 *
 *  Created on: 04.02.2009
 *      Author: karl
 *
 *  Modified by Michele Segata (segata@ccs-labs.org)
 */

#pragma once

#include "veins/veins.h"

#include "veins/modules/phy/DeciderResult80211.h"

namespace veins {

class VEINS_API DeciderResultRis : public DeciderResult80211 {
protected:
    /** @brief Stores the gain applied by the ris TODO: transform this into an array */
    double gain;
    double phiR, thetaR, phiI, thetaI;
    bool reflected;

public:
    /**
     * @brief Initialises with the passed values.
     *
     * "bitrate" defines the bit-rate of the transmission of the packet.
     */
    DeciderResultRis(bool isCorrect, double bitrate, double snr, double recvPower_dBm = 0, bool collision = false, double gain=0, double phiR=0, double thetaR=0, double phiI=0, double thetaI=0, bool reflected=false)
        : DeciderResult80211(isCorrect, bitrate, snr, recvPower_dBm, collision)
        , gain(gain)
        , phiR(phiR)
        , thetaR(thetaR)
        , phiI(phiI)
        , thetaI(thetaI)
        , reflected(reflected)
    {
    }

    double getGain() const
    {
        return gain;
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
