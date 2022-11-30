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

#include "veins/base/phyLayer/BaseDecider.h"

namespace veins {

using veins::AirFrame;

/**
 * @brief
 * Based on Decider80211.h from Karl Wessel
 * and modifications by Christopher Saloman
 *
 * @author Agon Memedi
 *
 * @ingroup decider
 *
 * @see PhyLayerRis
 * @see DeciderRis
 */
class DeciderRis : public BaseDecider {
public:
    enum PACKET_OK_RESULT {
        DECODED,
        NOT_DECODED,
        COLLISION
    };

protected:
    bool debug = true;
    double bitrate;

    bool ignoreNonReflectedSignals = false;
    bool ignoreNoiseAndInterference = false;

    double myBusyTime;
    double myStartTime;

    std::map<AirFrame*, int> signalStates;
    bool collectCollisionStats;
    unsigned int collisions;

protected:
    /**
     * @brief Checks a mapping against a specific threshold (element-wise).
     *
     * @return    true    , if every entry of the mapping is above threshold
     *             false    , otherwise
     *
     *
     */
    virtual DeciderResult* checkIfSignalOk(AirFrame* frame);

    virtual simtime_t processNewSignal(AirFrame* frame);

    /**
     * @brief Processes a received AirFrame.
     *
     * The SNR-mapping for the Signal is created and checked against the Deciders
     * SNR-threshold. Depending on that the received AirFrame is either sent up
     * to the MAC-Layer or dropped.
     *
     * @return    usually return a value for: 'do not pass it again'
     */
    virtual simtime_t processSignalEnd(AirFrame* frame);

    /** @brief computes if packet is ok or has errors*/
    enum DeciderRis::PACKET_OK_RESULT packetOk(double snirMin, double snrMin, int lengthMPDU, double bitrate);

public:
    /**
     * @brief Initializes the Decider with a pointer to its PhyLayer and
     * specific values for threshold and sensitivity
     */
    DeciderRis(cComponent* owner, DeciderToPhyInterface* phy, double sensitivity, double bRate, bool ignoreNonReflectedSignals, bool ignoreNoiseAndInterference, int myIndex = -1, bool collectCollisionStatistics = false)
        : BaseDecider(owner, phy, sensitivity, myIndex)
        , bitrate(bRate)
        , ignoreNonReflectedSignals(ignoreNonReflectedSignals)
        , ignoreNoiseAndInterference(ignoreNoiseAndInterference)
        , myBusyTime(0)
        , myStartTime(simTime().dbl())
        , collectCollisionStats(collectCollisionStatistics)
    {
    }

    int getSignalState(AirFrame* frame);
    virtual ~DeciderRis();
    /**
     * @brief invoke this method when the phy layer is also finalized,
     * so that statistics recorded by the decider can be written to
     * the output file
     */
    virtual void finish();
};

} // namespace veins
