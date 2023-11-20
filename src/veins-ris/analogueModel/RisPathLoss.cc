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

#include <cooperis/analogueModel/RisPathLoss.h>
#include <limits>

#include "veins/base/messages/AirFrame_m.h"
#include "cooperis/messages/AirFrameRis_m.h"

using namespace veins;

void RisPathLoss::filterSignal(Signal* signal)
{
    throw cRuntimeError("RisPathLoss::filterSignal(Signal* signal) must not be used as the model requires the AirFrame as well to perform computation");
}

void RisPathLoss::filterSignal(Signal* signal, AirFrame* frame)
{
    AirFrameRis* risMsg = check_and_cast<AirFrameRis*>(frame);

    auto sender = signal->getSenderPoa();
    auto receiver = signal->getReceiverPoa();

    auto senderPos = sender.pos.getPositionAt();
    auto receiverPos = receiver.pos.getPositionAt();

    double risGain = computeRisGain(signal, risMsg, senderPos, receiverPos);
    double pathLoss = computePathLoss(signal, risMsg, senderPos, receiverPos);

    double attenuation = std::min(1.0, risGain * pathLoss);
    double attenuation_dB = 10 * log10(attenuation);
    risMsg->appendActualLoss_dB(attenuation_dB);
    EV_TRACE << "RisPathLoss: applying actual attenuation of " << attenuation_dB << " dB\n";

    *signal *= attenuation;

}

double RisPathLoss::computePathLoss(Signal* signal, AirFrameRis* frame, const Coord& senderPos, const Coord& receiverPos)
{

    double attenuation = 1;

    /** Calculate the distance factor */
    double sqrDistance = useTorus ? receiverPos.sqrTorusDist(senderPos, playgroundSize) : receiverPos.sqrdist(senderPos);
    double distance = std::sqrt(sqrDistance);

    // save the length of the current path
    frame->appendPaths(distance);

    double totalDistance = frame->getTotalDistance() + distance;

    EV_TRACE << "sqrdistance is: " << sqrDistance << endl;

    if (sqrDistance <= 1.0) {
        // attenuation is negligible
        return attenuation;
    }

    if (useProductOfDistances || frame->getPathsArraySize() == 1) {
        // using far field model, simply compute the path loss independently
        // the part of the attenuation only depending on the distance
        double distFactor = pow(sqrDistance, -pathLossAlphaHalf) / (16.0 * M_PI * M_PI);
        EV_TRACE << "distance factor is: " << distFactor << endl;
        attenuation = wavelengthSquare * distFactor;
    }
    else {
        attenuation = pow(frame->getTotalDistance() / totalDistance, pathLossAlpha);
    }
    double attenuation_dB = 10*log10(attenuation);
    EV_TRACE << "Applying an attenuation of " << attenuation_dB << " dB\n";
    frame->appendPathLoss_dB(attenuation_dB);

    // save the total length travelled
    frame->setTotalDistance(totalDistance);
    return attenuation;

}

double RisPathLoss::computeRisGain(Signal* signal, AirFrameRis* frame, const Coord& senderPos, const Coord& receiverPos)
{

    PhyLayerRis* source = check_and_cast<PhyLayerRis*>(frame->getSenderModule());

    if (source->isReflectiveMetaSurface()) {
        Angles reflection = spherical_angles(frame->getRis_v1(), frame->getRis_v2(), frame->getRis_vn(), senderPos, receiverPos);
        double phiI = frame->getIncidencePhi(frame->getIncidencePhiArraySize()-1);
        double thetaI = frame->getIncidenceTheta(frame->getIncidenceThetaArraySize()-1);
        EV_TRACE << phyLayer->getFullPath() << " in position " << receiverPos.x << " " << receiverPos.y << " " << receiverPos.z << "\n";
        EV_TRACE << phyLayer->getFullPath() << " incoming signal from " << senderPos.x << " " << senderPos.y << " " << senderPos.z << "\n";
        EV_TRACE << phyLayer->getFullPath() << " Received signal from RIS. Incidence phi=" << RAD_TO_DEG(phiI) << " theta=" << RAD_TO_DEG(thetaI) << " Reflection phi=" << RAD_TO_DEG(reflection.phi) << " theta=" << RAD_TO_DEG(reflection.theta) << "\n";
        double gain = source->getMetasurfaceGain(reflection.phi, reflection.theta, phiI, thetaI);
        double gain_dB = 10 * log10(gain);
        frame->appendRisGain_dB(gain_dB);
        frame->appendReflectionPhi(reflection.phi);
        frame->appendReflectionTheta(reflection.theta);
        EV_TRACE << "Applying a gain of " << gain << "(" << gain_dB << " dB) to the incoming signal\n";

        return gain;
    }
    else {
        return 1;
    }

}
