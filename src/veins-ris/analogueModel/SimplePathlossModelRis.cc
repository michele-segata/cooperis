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

#include "SimplePathlossModelRis.h"

#include "veins-ris/messages/AirFrameRis_m.h"

using namespace veins;

using veins::AirFrame;

void SimplePathlossModelRis::filterSignal(Signal* signal)
{
    throw cRuntimeError("SimplePathlossModelRis::filterSignal(Signal* signal) must not be used as the model requires the AirFrame as well to perform computation");
}

void SimplePathlossModelRis::filterSignal(Signal* signal, AirFrame* frame)
{

    AirFrameRis* frameRis = check_and_cast<AirFrameRis*>(frame);

    auto senderPos = signal->getSenderPoa().pos.getPositionAt();
    auto receiverPos = signal->getReceiverPoa().pos.getPositionAt();

    /** Calculate the distance factor */
    double sqrDistance = useTorus ? receiverPos.sqrTorusDist(senderPos, playgroundSize) : receiverPos.sqrdist(senderPos);
    double distance = std::sqrt(sqrDistance);

    // save the length of the current path
    frameRis->appendPaths(distance);

    double totalDistance = frameRis->getTotalDistance() + distance;

    // we do not apply path loss at the metasurface, only at the final receiver
    //    if (receiver->isReflectiveMetaSurface())
    //        return;

    EV_TRACE << "sqrdistance is: " << sqrDistance << endl;

    if (sqrDistance <= 1.0) {
        // attenuation is negligible
        return;
    }

    Signal attenuation(signal->getSpectrum());

    if (useProductOfDistances || frameRis->getPathsArraySize() == 1) {
        // using far field model, simply compute the path loss independently
        // the part of the attenuation only depending on the distance
        double distFactor = pow(sqrDistance, -pathLossAlphaHalf) / (16.0 * M_PI * M_PI);
        EV_TRACE << "distance factor is: " << distFactor << endl;
        for (uint16_t i = 0; i < signal->getNumValues(); i++) {
            double wavelength = BaseWorldUtility::speedOfLight() / signal->getSpectrum().freqAt(i);
            double att = (wavelength * wavelength) * distFactor;
            attenuation.at(i) = att;
            if (i == 0)
                frameRis->appendPathLoss_dB(10*log10(att));
        }
    }
    else {
        double additionalPathFactor = pow(frameRis->getTotalDistance() / totalDistance, pathLossAlpha);
        for (uint16_t i = 0; i < signal->getNumValues(); i++) {
            attenuation.at(i) = additionalPathFactor;
            if (i == 0)
                frameRis->appendPathLoss_dB(10*log10(additionalPathFactor));
        }
    }

    // save the total length travelled
    frameRis->setTotalDistance(totalDistance);
    *signal *= attenuation;

}
