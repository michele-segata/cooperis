//
// Copyright (C) 2022-2024 Michele Segata <segata@ccs-labs.org>
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
#include "cooperis/messages/AirFrameRis_m.h"

namespace veins {

class VEINS_API DeciderResultRis : public DeciderResult80211 {
protected:
    long int frameId;
    std::vector<double> gains_dB;
    std::vector<double> paths_m;
    std::vector<double> losses_dB;
    std::vector<double> actualLosses_dB;
    std::vector<double> phiRs;
    std::vector<double> thetaRs;
    std::vector<double> phiIs;
    std::vector<double> thetaIs;
    bool reflected;

public:
    /**
     * @brief Initialises with the passed values.
     *
     * "bitrate" defines the bit-rate of the transmission of the packet.
     */
    DeciderResultRis(bool isCorrect, double bitrate, double snr, double recvPower_dBm = 0, bool collision = false, const AirFrameRis* frame = nullptr)
        : DeciderResult80211(isCorrect, bitrate, snr, recvPower_dBm, collision)
    {
        frameId = frame->getId();

        // store list of gains
        gains_dB.reserve(frame->getRisGain_dBArraySize());
        for (int i = 0; i < frame->getRisGain_dBArraySize(); i++)
            gains_dB.push_back(frame->getRisGain_dB(i));

        // store the length of reflected paths
        paths_m.reserve(frame->getPathsArraySize());
        for (int i = 0; i < frame->getPathsArraySize(); i++)
            paths_m.push_back(frame->getPaths(i));

        // store list of path losses
        losses_dB.reserve(frame->getPathLoss_dBArraySize());
        for (int i = 0; i < frame->getPathLoss_dBArraySize(); i++)
            losses_dB.push_back(frame->getPathLoss_dB(i));

        // store list of actual losses (gain + loss)
        actualLosses_dB.reserve(frame->getActualLoss_dBArraySize());
        for (int i = 0; i < frame->getActualLoss_dBArraySize(); i++)
            actualLosses_dB.push_back(frame->getActualLoss_dB(i));

        // store list of incidence and reflection angles
        phiRs.reserve(frame->getReflectionPhiArraySize());
        for (int i = 0; i < frame->getReflectionPhiArraySize(); i++)
            phiRs.push_back(frame->getReflectionPhi(i));
        thetaRs.reserve(frame->getReflectionThetaArraySize());
        for (int i = 0; i < frame->getReflectionThetaArraySize(); i++)
            thetaRs.push_back(frame->getReflectionTheta(i));
        phiIs.reserve(frame->getIncidencePhiArraySize());
        for (int i = 0; i < frame->getIncidencePhiArraySize(); i++)
            phiIs.push_back(frame->getIncidencePhi(i));
        thetaIs.reserve(frame->getIncidenceThetaArraySize());
        for (int i = 0; i < frame->getIncidenceThetaArraySize(); i++)
            thetaIs.push_back(frame->getIncidenceTheta(i));

        // store whether frame has been reflected
        reflected = frame->getReflected();
    }

    DeciderResultRis(const DeciderResultRis& r)
        :DeciderResult80211(r)
        , frameId(r.frameId)
        , gains_dB(r.gains_dB)
        , paths_m(r.paths_m)
        , losses_dB(r.losses_dB)
        , actualLosses_dB(r.actualLosses_dB)
        , phiRs(r.phiRs)
        , thetaRs(r.thetaRs)
        , phiIs(r.phiIs)
        , thetaIs(r.thetaIs)
    {}

    long int getFrameId() const
    {
        return frameId;
    }

    std::vector<double> getGains_dB() const
    {
        return gains_dB;
    }

    std::vector<double> getLosses_dB() const
    {
        return losses_dB;
    }

    std::vector<double> getActualLosses_dB() const
    {
        return actualLosses_dB;
    }

    std::vector<double> getPaths_m() const
    {
        return paths_m;
    }

    std::vector<double> getPhiRs() const
    {
        return phiRs;
    }
    std::vector<double> getThetaRs() const
    {
        return thetaRs;
    }
    std::vector<double> getPhiIs() const
    {
        return phiIs;
    }
    std::vector<double> getThetaIs() const
    {
        return thetaIs;
    }

    bool getReflected() const
    {
        return reflected;
    }
};

} // namespace veins
