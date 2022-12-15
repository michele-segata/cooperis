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

#include <veins-ris/analogueModel/RisPathLoss.h>
#include <limits>

#include "veins/base/messages/AirFrame_m.h"
#include "veins-ris/messages/AirFrameRis_m.h"

using namespace veins;

void RisPathLoss::filterSignal(Signal* signal)
{
    throw cRuntimeError("RisPathLoss::filterSignal(Signal* signal) must not be used as the model requires the AirFrame as well to perform computation");
}

void RisPathLoss::filterSignal(Signal* signal, AirFrame* frame)
{
    AirFrameRis* risMsg = check_and_cast<AirFrameRis*>(frame);
    PhyLayerRis* source = check_and_cast<PhyLayerRis*>(frame->getSenderModule());

    auto sender = signal->getSenderPoa();
    auto receiver = signal->getReceiverPoa();

    auto senderPos = sender.pos.getPositionAt();
    auto receiverPos = receiver.pos.getPositionAt();

    if (source->isReflectiveMetaSurface()) {
        Angles reflection = spherical_angles(risMsg->getRis_v1(), risMsg->getRis_v2(), risMsg->getRis_vn(), senderPos, receiverPos);
        EV_TRACE << phyLayer->getFullPath() << " in position " << receiverPos.x << " " << receiverPos.y << " " << receiverPos.z << "\n";
        EV_TRACE << phyLayer->getFullPath() << " incoming signal from " << senderPos.x << " " << senderPos.y << " " << senderPos.z << "\n";
        EV_TRACE << phyLayer->getFullPath() << " Received signal from RIS. Incidence phi=" << RAD_TO_DEG(risMsg->getIncidencePhi()) << " theta=" << RAD_TO_DEG(risMsg->getIncidenceTheta()) << " Reflection phi=" << RAD_TO_DEG(reflection.phi) << " theta=" << RAD_TO_DEG(reflection.theta) << "\n";
        double gain = source->getMetasurfaceGain(reflection.phi, reflection.theta, risMsg->getIncidencePhi(), risMsg->getIncidenceTheta());
        risMsg->appendRisGains(gain);
        risMsg->setReflectionPhi(reflection.phi);
        risMsg->setReflectionTheta(reflection.theta);
        EV_TRACE << "Applying a gain of " << gain << "(" << (10*log10(gain)) << " dB) to the incoming signal\n";

        *signal *= gain;
    }

}
