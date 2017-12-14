//
// Copyright (C) 2017 Agon Memedi <memedi@ccs-labs.org>
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

#include "veins-rms/Splitter.h"
#include "veins-rms/messages/RmsMessage_m.h"
#include "veins-rms/utility/Utils.h"

using namespace veins;

Define_Module(veins::Splitter);

Splitter::Splitter()
    : annotationManager(NULL)

{
}

Splitter::~Splitter()
{
}

void Splitter::initialize()
{
    // From upper layers --> lower layers
    fromApplication = findGate("applicationIn");
    toDsrcNic = findGate("nicOut");
    toRms = findGate("nicRmsOut");

    // From lower layers --> upper layers
    toApplication = findGate("applicationOut");
    fromDsrcNic = findGate("nicIn");
    fromRms = findGate("nicRmsIn");

    // Module parameters
    collectStatistics = par("collectStatistics").boolValue();
    debug = par("debug").boolValue();

    // Other simulation modules
    cModule* tmpMobility = getParentModule()->getSubmodule("veinsmobility");
    mobility = dynamic_cast<TraCIMobility*>(tmpMobility);
    ASSERT(mobility);

    rmsPhys = getSubmodulesOfType<PhyLayerRms>(getParentModule(), true);
    ASSERT(rmsPhys.size() > 0);

    annotationManager = AnnotationManagerAccess().getIfExists();
    ASSERT(annotationManager);
}

void Splitter::handleMessage(cMessage* msg)
{
    if (timerManager.handleMessage(msg)) return;

    if (msg->isSelfMessage()) {
        error("Self-message arrived!");
        delete msg;
        msg = NULL;
    }
    else {
        int arrivalGate = msg->getArrivalGateId();
        if (arrivalGate == fromApplication) {
            handleUpperMessage(msg);
        }
        // The arrival gate is not from the application, it'a from lower layers
        // TODO: add annotation drawings based on the destination technology
        else {
            EV_INFO << "Message from lower layers received!" << std::endl;
            handleLowerMessage(msg);
        }
    }
}

void Splitter::handleUpperMessage(cMessage* msg)
{
    // Cast the message to a subclass
    RmsMessage* rmsMsg = dynamic_cast<RmsMessage*>(msg);

    // Handle WSMs if the VLC has to be "retrofitted" to non-VLC app
    if (!rmsMsg) {
        BaseFrame1609_4* wsm = dynamic_cast<BaseFrame1609_4*>(msg);

        // If not a RmsMessage check whether it is a WSM to send directly
        if (wsm) {
            send(wsm, toDsrcNic);
            return;
        }
        else
            error("Not a RmsMessage, not BaseFrame1609_4");
    }

    // if (rmsMsg)...
    int networkType = rmsMsg->getAccessTechnology();
    if (networkType == DSRC) {
        EV_INFO << "DSRC message received from upper layer!" << std::endl;
        send(rmsMsg, toDsrcNic);
    }
    else if (networkType == RMS) {
        EV_INFO << "RMS message received from upper layer!" << std::endl;

        send(rmsMsg, toRms);
    }
    else {
        error("\tThe access technology has not been specified in the message!");
    }
}

void Splitter::handleLowerMessage(cMessage* msg)
{
    send(msg, toApplication);
}

void Splitter::drawRayLine(const AntennaPosition& ap, int length, double halfAngle, bool reverse)
{
    double heading = mobility->getHeading().getRad();
    // This is for the cone of the tail
    if (reverse) heading = reverseTraci(heading);

    annotationManager->scheduleErase(0.1,
        annotationManager->drawLine(ap.getPositionAt(),
        ap.getPositionAt() + Coord(length * cos(halfAngle + traci2myAngle(heading)), length * sin(halfAngle + traci2myAngle(heading))),
        "white"));
}

void Splitter::finish()
{
}
