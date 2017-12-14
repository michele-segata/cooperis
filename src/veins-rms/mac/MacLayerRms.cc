//
// Copyright (C) 2019 Max Schettler <schettler@ccs-labs.org>
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

#include "veins-rms/mac/MacLayerRms.h"

#include "veins/base/phyLayer/PhyToMacControlInfo.h"
#include "veins/base/phyLayer/MacToPhyInterface.h"
#include "veins/base/messages/MacPkt_m.h"
#include "veins-rms/PhyLayerRms.h"
#include "veins-rms/DeciderResultRms.h"

using namespace veins;

using std::unique_ptr;

Define_Module(veins::MacLayerRms);

void MacLayerRms::initialize(int stage)
{
    BaseMacLayer::initialize(stage);
    if (stage == 0) {

        // TODO: port these parameters to MAC (?)
        // bit rate
        // int bitrate @unit(bps) = default(6 Mbps);
        // tx power [mW]
        // double txPower @unit(mW);

        queueSize = par("queueSize");
        transmitting = false;
    }
}

void MacLayerRms::handleUpperMsg(cMessage* msg)
{
    EV_TRACE << "Received a message from upper layer" << std::endl;
    ASSERT(dynamic_cast<cPacket*>(msg));

    if (transmitting) {
        enqueuePacket(check_and_cast<cPacket*>(msg));
    }
    else {
        if (!queue.isEmpty()) throw cRuntimeError("Radio not transmitting but packets in queue");
        sendDown(encapsMsg(check_and_cast<cPacket*>(msg)));
        transmitting = true;
    }
}

void MacLayerRms::handleLowerMsg(cMessage* msg)
{
    MacPkt* mac = static_cast<MacPkt*>(msg);
    LAddress::L2Type dest = mac->getDestAddr();
    LAddress::L2Type src = mac->getSrcAddr();

    // pass information about received frame to the upper layers
    DeciderResultRms* macRes = check_and_cast<DeciderResultRms*>(PhyToMacControlInfo::getDeciderResult(msg));
    DeciderResultRms* res = new DeciderResultRms(*macRes);

    // only foward to upper layer if message is for me or broadcast
    if ((dest == myMacAddr) || LAddress::isL2Broadcast(dest)) {
        EV_TRACE << "message with mac addr " << src << " for me (dest=" << dest << ") -> forward packet to upper layer\n";
        cPacket* encapsulatedMsg = mac->decapsulate();
        auto ctrlInfo = new PhyToMacControlInfo(res);
        ctrlInfo->setSourceAddress(mac->getSrcAddr());
        encapsulatedMsg->setControlInfo(ctrlInfo);
        sendUp(encapsulatedMsg);
    }
    else {
        EV_TRACE << "message with mac addr " << src << " not for me (dest=" << dest << ") -> delete (my MAC=" << myMacAddr << ")\n";
    }
    delete mac;
}

void MacLayerRms::handleLowerControl(cMessage* msg)
{
    switch (msg->getKind()) {
    case MacToPhyInterface::TX_OVER:
        transmitting = false;
        transmissionOpportunity();
        delete msg;
        break;
    default:
        EV << "MacLayerRms does not handle control messages of this type (name was " << msg->getName() << ")\n";
        delete msg;
        break;
    }
}

MacPkt* MacLayerRms::encapsMsg(cPacket* netwPkt)
{
    MacPkt* pkt = new MacPkt(netwPkt->getName(), netwPkt->getKind());
    pkt->addBitLength(headerLength);

    // TODO: setting up proper control info: according to the interfaces (?)

    pkt->setDestAddr(LAddress::L2BROADCAST());

    // set the src address to own mac address (nic module getId())
    pkt->setSrcAddr(myMacAddr);

    // encapsulate the network packet
    pkt->encapsulate(netwPkt);
    EV_TRACE << "pkt encapsulated\n";

    return pkt;
}

void MacLayerRms::enqueuePacket(cPacket* pkt)
{
    if (queueSize != -1 && queue.getLength() >= queueSize) {
        // emit(sigDropped, true);
        delete pkt;
        return; // TODO send FULL QUEUE message
    }
    queue.insert(pkt);
    // emit(sigQueueLength, queue.getLength());
}

void MacLayerRms::transmissionOpportunity()
{
    if (queue.isEmpty()) {
        return;
    }

    sendDown(encapsMsg(queue.pop()));
    // emit(sigQueueLength, queue.getLength());
    transmitting = true;
}
