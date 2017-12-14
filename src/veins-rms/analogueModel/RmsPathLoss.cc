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

#include <veins-rms/analogueModel/RmsPathLoss.h>
#include <limits>

#include "veins/base/messages/AirFrame_m.h"
#include "veins-rms/messages/AirFrameRms_m.h"

using namespace veins;

//#define EV_TRACE \
//   if (debug) EV_LOG(omnetpp::LOGLEVEL_TRACE, nullptr) << "[RsmPathLoss] "

#define EV_TRACE std::cout << "[RmsPathLoss] "

void RmsPathLoss::filterSignal(Signal* signal)
{
    throw cRuntimeError("RmsPathLoss::filterSignal(Signal* signal) must not be used as the model requires the AirFrame as well to perform computation");
}

void RmsPathLoss::filterSignal(Signal* signal, AirFrame* frame)
{
    AirFrameRms* rmsMsg = check_and_cast<AirFrameRms*>(frame);
    PhyLayerRms* source = check_and_cast<PhyLayerRms*>(frame->getSenderModule());

    auto sender = signal->getSenderPoa();
    auto receiver = signal->getReceiverPoa();

    auto senderPos = sender.pos.getPositionAt();
    auto receiverPos = receiver.pos.getPositionAt();

    // TODO: not yet considering shadowing. see later how to include it

    /* distance = computedistance
     * if (from RMS) {
     *   compute reflection angles at previous RMS
     *   compute and apply path loss given incidence, reflection, distance, previous RMS config
     * }
     * else { // from UE
     *   compute and apply classic path loss given distance
     * }
     *
     * // TODO: next part seems to belong to the phy layer, not the channel model
     * if (power < threshold)
     *   discard packet and stop here
     *
     * if (node type == RMS) {
     *   compute incidence angles
     *   add incidence angles and RMS config to AirFrame
     *   re-send frame
     * }
     *
     */

    if (source->isReflectiveMetaSurface()) {
        Angles reflection = spherical_angles(rmsMsg->getRms_v1(), rmsMsg->getRms_v2(), rmsMsg->getRms_vn(), senderPos, receiverPos);
//        EV_TRACE << phyLayer->getFullPath() << " in position " << receiverPos.x << " " << receiverPos.y << " " << receiverPos.z << "\n";
//        EV_TRACE << phyLayer->getFullPath() << " incoming signal from " << senderPos.x << " " << senderPos.y << " " << senderPos.z << "\n";
        EV_TRACE << phyLayer->getFullPath() << " Received signal from RMS. Incidence phi=" << RAD_TO_DEG(rmsMsg->getIncidencePhi()) << " theta=" << RAD_TO_DEG(rmsMsg->getIncidenceTheta()) << " Reflection phi=" << RAD_TO_DEG(reflection.phi) << " theta=" << RAD_TO_DEG(reflection.theta) << "\n";
        double gain = source->getMetasurfaceGain(reflection.phi, reflection.theta, rmsMsg->getIncidencePhi(), rmsMsg->getIncidenceTheta());
        rmsMsg->setRmsGain(gain);
        rmsMsg->setReflectionPhi(reflection.phi);
        rmsMsg->setReflectionTheta(reflection.theta);
        if (rmsMsg->getOriginalSenderModule() != phyLayer)
            EV_TRACE << "phi=" << reflection.phi << " gain=" << gain << "\n";
//            EV_TRACE << "Applying a gain of " << gain << "(" << (10*log10(gain)) << " dB) to the incoming signal\n";

        *signal *= gain;
    }
    // else branch (source == FROM_UE) handled by SimplePathlossModelRms
    // TODO: missing path loss on second path!!!! we have 2 freespace path loss. rms gain is in addition

    //    if (phyLayer->isReflectiveMetaSurface()) {
    //        Angles incident = get_angles(phyLayer->getRms_v1(), phyLayer->getRms_v2(), phyLayer->getRms_vn(), receiverPos, senderPos);
    //        EV_TRACE << "RMS: incoming signal. Incidence phi (azimuth): " << RAD_TO_DEG(incident.phi) << " incidence theta (elevation): " << RAD_TO_DEG(incident.theta) << "\n";
    //        // TODO: check that dup() properly works. when copying, the signal remains the pointer to the original one
    //        // if there are memory errors due to airframes and signals, this might be the cause
    //        RMS_INFO* info = new RMS_INFO();
    //        info->repropagate = true;
    //        info->incidence = incident;
    //        return info;
    //    }

}
