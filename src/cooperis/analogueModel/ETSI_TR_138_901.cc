//
// Copyright (C) 2024 Michele Segata <segata@ccs-labs.org>
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

#include "ETSI_TR_138_901.h"
#include "cooperis/messages/AirFrameRis_m.h"
#include <omnetpp.h>

using namespace veins;

using veins::AirFrame;

ETSI_TR_138_901::ETSI_TR_138_901(cComponent* owner, ObstacleControl& obstacleControl, bool useTorus, const Coord& playgroundSize)
    : SimpleObstacleShadowing(owner, obstacleControl, useTorus, playgroundSize)
{
    if (useTorus) throw cRuntimeError("SimpleObstacleShadowing does not work on torus-shaped playgrounds");
}
void ETSI_TR_138_901::filterSignal(Signal* signal, AirFrame* frame)
{
    AirFrameRis* frameRis = check_and_cast<AirFrameRis*>(frame);

    auto senderPos = signal->getSenderPoa().pos.getPositionAt();
    auto receiverPos = signal->getReceiverPoa().pos.getPositionAt();

    // TODO: not the smartest way, but the quickest: use the base veins obstacle shadowing model to determine whether we are in LOS or NLOS
    double factor = obstacleControl.calculateAttenuation(senderPos, receiverPos);

    if (signal->getNumValues() > 1)
        EV_WARN << "the input signal has multiple frequencies. Applying path loss only on the first one" << endl;

    double freq = signal->getSpectrum().freqAt(0);

    EV_TRACE << "distance: " << senderPos.distance(receiverPos) << " m\n";

    double pl_db = 0;
    if (factor < 1) {
        // tag the frame as being shadowed
        frameRis->setShadowed(true);
        // we are in NLOS conditions
        pl_db = nonLineOfSightPathLoss(freq, senderPos, receiverPos);
        EV_TRACE << "NLOS pathloss: " << pl_db << " dB\n";
        EV_TRACE << "LOS pathloss would be " << lineOfSightPathLoss(freq, senderPos, receiverPos) << " dB\n";
    }
    else {
        // we are in LOS conditions. let other models take care of LOS pathloss
        pl_db = 0;
    }
    factor = 1 / pow(10, pl_db / 10);

    EV_TRACE << "value is: " << factor << endl;

    *signal *= factor;
}

void ETSI_TR_138_901::filterSignal(Signal* signal)
{
    throw cRuntimeError("SimpleObstacleShadowingRis::filterSignal(Signal* signal) must not be used as the model requires the AirFrame as well to perform computation");
}

double ETSI_TR_138_901::lineOfSightPathLoss(double f, const Coord &sender, const Coord &receiver)
{
    // shadow fading sigma
    static double sigma_sf = 4;
    static double h_e = 1.0;
    static double c = 299792458;
    double h_s = sender.z;
    double h_r = receiver.z;
    double d_bp = 4 * (h_s - h_e) * (h_r - h_e) * f / c;
    Coord sender2d = sender;
    Coord receiver2d = receiver;
    sender2d.z = 0;
    receiver2d.z = 0;
    double d2d = sender2d.distance(receiver2d);
    double d3d = sender.distance(receiver);
    double pl_db = 0;

    if (d2d <= d_bp) {
        // TODO: the formula in the TR is valid only for distances greater than 10 meters
        // Here we assume that below 10 meters we have the same attenuation
        pl_db = 32.4 + 21 * log10(d3d) + 20 * log10(f/1e9) + normal(getEnvir()->getRNG(0), 0, sigma_sf);
        if (d2d < 10)
            EV_WARN << "using ETSI TR 138.901 path loss model for a distance smaller than 10 m" << endl;
    }
    else
        // TODO: additional assumption. the base station is assumed to be 10 meters from ground
        // here we use the actual height of the transmitter
        pl_db = 32.4 + 40 * log10(d3d) + 20 * log10(f/1e9) - 9.5 * log10(pow(d_bp, 2) + pow(h_s - h_r, 2)) + normal(getEnvir()->getRNG(0), 0, sigma_sf);
    return pl_db;
}

double ETSI_TR_138_901::nonLineOfSightPathLoss(double f, const Coord &sender, const Coord &receiver)
{
    // shadow fading sigma
    static double sigma_sf = 7.82;
    double h_r = receiver.z;
    Coord sender2d = sender;
    Coord receiver2d = receiver;
    sender2d.z = 0;
    receiver2d.z = 0;
    double d2d = sender2d.distance(receiver2d);
    double d3d = sender.distance(receiver);
    if (d2d < 10)
        EV_WARN << "using ETSI TR 138.901 path loss model for a distance smaller than 10 m" << endl;
    else if (d2d > 5000)
        EV_WARN << "using ETSI TR 138.901 path loss model for a distance larger than 5 km" << endl;
    double pl_db = 35.3 * log10(d3d) + 22.4 + 21.3 * log10(f/1e9) - 0.3 * (h_r - 1.5) + normal(getEnvir()->getRNG(0), 0, sigma_sf);
    return std::max(lineOfSightPathLoss(f, sender, receiver), pl_db);
}
