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

#include <omnetpp.h>
#include <veins-rms/application/simpleRmsApp/SimpleRmsApp.h>

using namespace veins;

Define_Module(veins::SimpleRmsApp);

SimpleRmsApp::SimpleRmsApp()
    : traciManager(NULL)
    , mobility(NULL)
    , annotations(NULL)
{
}

SimpleRmsApp::~SimpleRmsApp()
{
}

void SimpleRmsApp::initialize(int stage)
{
    if (stage == 0) {
        toLower = findGate("lowerLayerOut");
        fromLower = findGate("lowerLayerIn");

        // Pointers to simulation modules
        traciManager = TraCIScenarioManagerAccess().get();
        ASSERT(traciManager);

        cModule* tmpMobility = getParentModule()->getSubmodule("veinsmobility");
        mobility = dynamic_cast<TraCIMobility*>(tmpMobility);
        ASSERT(mobility);

        annotations = AnnotationManagerAccess().getIfExists();
        ASSERT(annotations);

        sumoId = mobility->getExternalId();
        debug = par("debug").boolValue();
        byteLength = par("packetByteLength");
        transmissionPeriod = 1 / par("beaconingFrequency").doubleValue();
    }
    if (stage == 3) {
        auto dsrc = [this]() {
            RmsMessage* rmsMsg = new RmsMessage();
            rmsMsg->setAccessTechnology(DSRC);
            send(rmsMsg, toLower);
        };
        auto vlc = [this]() {
            RmsMessage* rmsMsg = generateRmsMessage(RMS);
            send(rmsMsg, toLower);
        };
        timerManager.create(veins::TimerSpecification(vlc).oneshotAt(SimTime(20, SIMTIME_S)));
    }
}

void SimpleRmsApp::handleMessage(cMessage* msg)
{
    // To handle the timer
    if (timerManager.handleMessage(msg)) return;

    if (msg->isSelfMessage()) {
        throw cRuntimeError("This module does not use custom self messages");
        return;
    }
    else {
        RmsMessage* rmsMsg = check_and_cast<RmsMessage*>(msg);
        int accessTech = rmsMsg->getAccessTechnology();
        switch (accessTech) {
        case DSRC: {
            EV_INFO << "DSRC message received!" << std::endl;
            delete rmsMsg;
            break;
        }
        case RMS: {
            EV_INFO << "RMS message received from: " << rmsMsg->getSourceNode() << std::endl;
            delete rmsMsg;
            break;
        }
        default:
            error("message neither from DSRC nor RMS");
            break;
        }
    }
}

RmsMessage* SimpleRmsApp::generateRmsMessage(int accessTechnology)
{
    RmsMessage* rmsMsg = new RmsMessage();

    // OMNeT-specific
    rmsMsg->setName("rmsMessage");

    // WSM fields

    // HeterogeneousMessage specific
    rmsMsg->setSourceNode(this->sumoId.c_str());
    rmsMsg->setDestinationNode("BROADCAST");
    rmsMsg->setAccessTechnology(accessTechnology);
    rmsMsg->setSentAt(simTime()); // There is timestamp field in WSM too

    // Set application layer packet length
    rmsMsg->setByteLength(byteLength);

    return rmsMsg;
}
