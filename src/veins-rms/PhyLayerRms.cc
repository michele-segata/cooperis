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

/*
 * Based on PhyLayer80211p.cc from David Eckhoff
 */

#include "veins/veins.h"
#include "veins-rms/PhyLayerRms.h"
#include "veins-rms/utility/Utils.h"

#include "veins-rms/analogueModel/VehicleObstacleShadowingForVlc.h"
#include "veins/modules/analogueModel/SimpleObstacleShadowing.h"
#include "veins/base/connectionManager/BaseConnectionManager.h"
#include "veins/modules/messages/AirFrame11p_m.h"
#include "veins-rms/messages/AirFrameRms_m.h"
#include "veins-rms/utility/Utils.h"
#include "veins/base/toolbox/Spectrum.h"
#include "veins-rms/analogueModel/RmsPathLoss.h"
#include "veins-rms/analogueModel/SimplePathlossModelRms.h"

#include <memory>

using namespace veins;

using std::unique_ptr;

Define_Module(veins::PhyLayerRms);

// #define EV_TRACE std::cout << "[PhyLayerRms] "

void PhyLayerRms::initialize(int stage)
{
    if (stage == 0) {
        bitrate = par("bitrate");
        centerFrequency = par("centerFrequency");
        txPower = par("txPower");
        isRMS = par("isRMS");

        rms_v2.x = par("rms_v2_x");
        rms_v2.y = par("rms_v2_y");
        rms_v2.z = par("rms_v2_z");

        rms_v1.x = par("rms_v1_x");
        rms_v1.y = par("rms_v1_y");
        rms_v1.z = par("rms_v1_z");

        if (isRMS && ((rms_v2.squareLength() < 1e-6) || (rms_v1.squareLength() < 1e-6)))
            throw cRuntimeError("Cannot instantiate an RMS without setting the orientation vectors v2 and v1");

        if (isRMS) {
            // normalize vectors
            double v2Len = rms_v2.length();
            rms_v2.x /= v2Len;
            rms_v2.y /= v2Len;
            rms_v2.z /= v2Len;
            double v1Len = rms_v1.length();
            rms_v1.x /= v1Len;
            rms_v1.y /= v1Len;
            rms_v1.z /= v1Len;
            // compute facing direction
            rms_vn = cross(rms_v2, rms_v1);

            initMetasurface = new cMessage("initMetasurface");
            scheduleAt(simTime() + 1e-6, initMetasurface);

            rms = new ReconfigurableMetaSurface(centerFrequency);
            annotations = AnnotationManagerAccess().getIfExists();

            nodesAntennaHeight = par("nodesAntennaHeight");
        }

        ignoreNonReflectedSignals = par("ignoreNonReflectedSignals");
        ignoreNoiseAndInterference = par("ignoreNoiseAndInterference");

        // Create frequency mappings and initialize spectrum for signal representation
        overallSpectrum = Spectrum({centerFrequency});
    }
    BasePhyLayer::initialize(stage);
}

const Coord& PhyLayerRms::getRms_v2() const
{
    return rms_v2;
}
const Coord& PhyLayerRms::getRms_v1() const
{
    return rms_v1;
}
const Coord& PhyLayerRms::getRms_vn() const
{
    return rms_vn;
}

bool PhyLayerRms::isReflectiveMetaSurface() const
{
    return isRMS;
}

double PhyLayerRms::getMetasurfaceGain(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad) const
{
    if (!isRMS) throw cRuntimeError("Requesting metasurface gain on a PHY that is not a metasurface but a standard node");
    double configPhiR_rad, configPhiI_rad, configThetaR_rad, configThetaI_rad;
    rms->getConfiguration(configPhiR_rad, configThetaR_rad, configPhiI_rad, configThetaI_rad);
    double gain = rms->gain(phiR_rad, thetaR_rad, phiI_rad, thetaI_rad);

//    printf("RMS config for phiR=%.2f, thetaR=%.2f, phiI=%.2f, thetaI=%.2f\n", RAD_TO_DEG(configPhiR_rad), RAD_TO_DEG(configThetaR_rad), RAD_TO_DEG(configPhiI_rad), RAD_TO_DEG(configThetaI_rad));
//    printf("Req. gain  for phiR=%.2f, thetaR=%.2f, phiI=%.2f, thetaI=%.2f\n", RAD_TO_DEG(phiR_rad), RAD_TO_DEG(thetaR_rad), RAD_TO_DEG(phiI_rad), RAD_TO_DEG(thetaI_rad));
//    printf("Gain: %.2f (%.2f dB)\n", gain, 10 * log10(gain));

    return gain;
}

void PhyLayerRms::configureMetaSurface(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad)
{
    if (!isRMS) throw cRuntimeError("Requesting to configure a metasurface on a PHY that is not a metasurface but a standard node");
    rms->configureMetaSurface(phiR_rad, thetaR_rad, phiI_rad, thetaI_rad);
    // visualize using AnnotationManager
    if (configVisualizationR)
        annotations->erase(configVisualizationR);
    if (configVisualizationI)
        annotations->erase(configVisualizationI);

    if (annotations) {
        std::list<Coord> projectedBeam;
        projectedBeam = project_beam(rms_v1, rms_v2, rms_vn, antennaPosition.getPositionAt(), phiR_rad, thetaR_rad, nodesAntennaHeight);
        configVisualizationR = annotations->drawPolygon(projectedBeam, "red", nullptr);
        projectedBeam = project_beam(rms_v1, rms_v2, rms_vn, antennaPosition.getPositionAt(), phiI_rad, thetaI_rad, nodesAntennaHeight);
        configVisualizationI = annotations->drawPolygon(projectedBeam, "green", nullptr);
    }

}

unique_ptr<AnalogueModel> PhyLayerRms::getAnalogueModelFromName(std::string name, ParameterMap& params)
{

    if (name == "SimplePathlossModelRms") {
        return initializeSimplePathlossModelRms(params);
    }
    if (name == "SimpleObstacleShadowing") {
        return initializeSimpleObstacleShadowing(params);
    }
    if (name == "RmsPathLoss") {
        return initializeRmsPathLoss(params);
    }
    return BasePhyLayer::getAnalogueModelFromName(name, params);
}

unique_ptr<AnalogueModel> PhyLayerRms::initializeSimpleObstacleShadowing(ParameterMap& params)
{

    // init with default value
    bool useTorus = world->useTorus();
    const Coord& playgroundSize = *(world->getPgs());

    ParameterMap::iterator it;

    ObstacleControl* obstacleControlP = dynamic_cast<ObstacleControl*>(veins::findModuleByPath("obstaclesRms"));
    if (!obstacleControlP) throw cRuntimeError("initializeSimpleObstacleShadowing(): cannot find ObstacleControl module");
    return make_unique<SimpleObstacleShadowing>(this, *obstacleControlP, useTorus, playgroundSize);
}

unique_ptr<AnalogueModel> PhyLayerRms::initializeRmsPathLoss(ParameterMap& params)
{
    double centerFrequency = params["carrierFrequency"].doubleValue();
    int n = params["n"].longValue();
    return make_unique<RmsPathLoss>(this, centerFrequency, n);
}

unique_ptr<AnalogueModel> PhyLayerRms::initializeSimplePathlossModelRms(ParameterMap& params)
{

    // init with default value
    double alpha = 2.0;
    bool useTorus = world->useTorus();
    const Coord& playgroundSize = *(world->getPgs());

    // get alpha-coefficient from config
    ParameterMap::iterator it = params.find("alpha");

    if (it != params.end()) { // parameter alpha has been specified in config.xml
        // set alpha
        alpha = it->second.doubleValue();
        EV_TRACE << "createPathLossModel(): alpha set from config.xml to " << alpha << endl;

        // check whether alpha is not smaller than specified in ConnectionManager
        if (cc->hasPar("alpha") && alpha < cc->par("alpha").doubleValue()) {
            // throw error
            throw cRuntimeError("TestPhyLayer::createPathLossModel(): alpha can't be smaller than specified in \
                   ConnectionManager. Please adjust your config.xml file accordingly");
        }
    }
    else // alpha has not been specified in config.xml
    {
        if (cc->hasPar("alpha")) { // parameter alpha has been specified in ConnectionManager
            // set alpha according to ConnectionManager
            alpha = cc->par("alpha").doubleValue();
            EV_TRACE << "createPathLossModel(): alpha set from ConnectionManager to " << alpha << endl;
        }
        else // alpha has not been specified in ConnectionManager
        {
            // keep alpha at default value
            EV_TRACE << "createPathLossModel(): alpha set from default value to " << alpha << endl;
        }
    }

    return make_unique<SimplePathlossModelRms>(this, alpha, useTorus, playgroundSize);
}

unique_ptr<Decider> PhyLayerRms::getDeciderFromName(std::string name, ParameterMap& params)
{
    if (name == "DeciderRms") {
        protocolId = RMS;
        return initializeDeciderRms(params);
    }
    return BasePhyLayer::getDeciderFromName(name, params);
}

unique_ptr<AnalogueModel> PhyLayerRms::initializeObstacleShadowing(ParameterMap& params)
{
    // TODO: finish
    // init with default value
    bool useTorus = world->useTorus();
    const Coord& playgroundSize = *(world->getPgs());

    VehicleObstacleControl* vehicleObstacleControlP = VehicleObstacleControlAccess().getIfExists();
    if (!vehicleObstacleControlP) throw cRuntimeError("initializeVehicleObstacleShadowingForVlc(): cannot find VehicleObstacleControl module");
    return make_unique<VehicleObstacleShadowingForVlc>(this, *vehicleObstacleControlP, useTorus, playgroundSize);
}

unique_ptr<Decider> PhyLayerRms::initializeDeciderRms(ParameterMap& params)
{
    DeciderRms* dec = new DeciderRms(this, this, minPowerLevel, bitrate, ignoreNonReflectedSignals, ignoreNoiseAndInterference, -1, false);
    return std::unique_ptr<DeciderRms>(std::move(dec));
}

void PhyLayerRms::handleSelfMessage(cMessage* msg)
{

    switch (msg->getKind()) {
    // transmission over
    case TX_OVER:
        ASSERT(msg == txOverTimer);
        if (!isRMS)
            sendControlMsgToMac(new cMessage("Transmission over", TX_OVER));
        break;

    // radio switch over
    case RADIO_SWITCHING_OVER:
        ASSERT(msg == radioSwitchingOverTimer);
        finishRadioSwitching();
        break;

    // AirFrame
    case AIR_FRAME:
        handleAirFrame(static_cast<AirFrame*>(msg));
        break;

    default:
        break;
    }

    if (msg == initMetasurface) {
        double phiR = DEG_TO_RAD(-84.2719);
        double thetaR = DEG_TO_RAD(85.3535);
        double phiI = DEG_TO_RAD(11.1671);
        double thetaI = DEG_TO_RAD(4.74773);
        configureMetaSurface(phiR, thetaR, phiI, thetaI);
//        configureMetaSurface(RMS_AZIMUTH_CENTER, RMS_ELEVATION_MIDDLE, RMS_AZIMUTH_CENTER, RMS_ELEVATION_MIDDLE);
        cancelAndDelete(initMetasurface);
        initMetasurface = nullptr;
    }

}

void PhyLayerRms::handleMessage(cMessage* msg)
{
    // self messages
    if (msg->isSelfMessage()) {
        handleSelfMessage(msg);

        // MacPkts <- MacToPhyControlInfo
    }
    else if (msg->getArrivalGateId() == upperLayerIn) {
        setRadioState(veins::Radio::TX);
        BasePhyLayer::handleUpperMessage(msg);

        // controlmessages
    }
    else if (msg->getArrivalGateId() == upperControlIn) {
        BasePhyLayer::handleUpperControlMessage(msg);

        // AirFrames
        // msg received over air from other NICs
    }
    else if (msg->getKind() == AIR_FRAME) {
        AirFrameRms* rmsMsg = check_and_cast<AirFrameRms*>(msg);

        if (rmsMsg->getOriginalSenderModule() == this) {
            delete msg;
            return;
        }
        else {
            //            EV_TRACE << "AirFrameVlc id: " << rmsMsg->getId() << " handed to VLC PHY from: " << rmsMsg->getSenderModule()->getFullPath() << std::endl;
            bubble("Handing AirFrameVlc to lower layers to decide if it can be received");
            BasePhyLayer::handleAirFrame(static_cast<AirFrame*>(msg));
        }
    }
    else {
        //        EV_TRACE << "Unknown message received." << endl;
        delete msg;
    }
}

void PhyLayerRms::sendUp(AirFrame* frame, DeciderResult* result)
{
    if (!isRMS)
        BasePhyLayer::sendUp(frame, result);
}

void PhyLayerRms::handleAirFrameReceiving(AirFrame* msg)
{
    // if this is an RMS, do not attempt decoding
    if (!isRMS)
        BasePhyLayer::handleAirFrameReceiving(msg);
}

unique_ptr<AirFrame> PhyLayerRms::createAirFrame(cPacket* macPkt)
{
    return make_unique<AirFrameRms>(macPkt->getName(), AIR_FRAME);
}

unique_ptr<AirFrame> PhyLayerRms::encapsMsg(cPacket* macPkt, double txPower, double incidencePhi, double incidenceTheta, bool reflected, int originalId)
{
    auto airFrame = createAirFrame(macPkt);
    AirFrameRms* frame = dynamic_cast<AirFrameRms*>(airFrame.get());

    // set the members
    // set priority of AirFrames above the normal priority to ensure
    // channel consistency (before any thing else happens at a time
    // point t make sure that the channel has removed every AirFrame
    // ended at t and added every AirFrame started at t)
    frame->setSchedulingPriority(airFramePriority());
    frame->setProtocolId(myProtocolId());
    frame->setId(world->getUniqueAirFrameId());
    frame->setChannel(radio->getCurrentChannel());

    frame->setIncidencePhi(incidencePhi);
    frame->setIncidenceTheta(incidenceTheta);
    frame->setReflected(reflected);
    frame->setOriginalId(originalId);
    frame->setOriginalSenderModule(this);

    if (isRMS) {
        frame->setRms_v1(rms_v1);
        frame->setRms_v2(rms_v2);
        frame->setRms_vn(rms_vn);
    }
    else {
        frame->setTotalDistance(0);
        EV_TRACE << "\n\n\n\n";
    }

    // encapsulate the mac packet into the phy frame
    frame->encapsulate(macPkt);

    // attachSignal()
    // attach the spectrum-dependent Signal to the airFrame
    enum MCS mcs = getMCS(bitrate, Bandwidth::ofdm_400_mhz);
    const auto duration = getFrameDuration(frame->getEncapsulatedPacket()->getBitLength(), mcs);
    ASSERT(duration > 0);
    Signal signal(overallSpectrum, simTime(), duration);
    signal.at(0) = txPower;
    signal.setDataStart(0);
    signal.setDataEnd(0);
    signal.setCenterFrequencyIndex(0);
    // copy the signal into the AirFrame
    frame->setSignal(signal);
    frame->setDuration(signal.getDuration());
    frame->setMcs((int)mcs);
    // --- end of attachSignal() ---

    // --- from here on, the AirFrame is the owner of the MacPacket ---
    macPkt = nullptr;
    //    EV_TRACE << "AirFrame w/ id: " << frame->getId() << " encapsulated, bit length: " << frame->getBitLength() << "\n";

    return airFrame;
}

unique_ptr<AirFrame> PhyLayerRms::encapsMsg(cPacket* macPkt)
{
    return encapsMsg(macPkt, txPower);
}

simtime_t PhyLayerRms::setRadioState(int rs)
{
    if (rs == Radio::TX) decider->switchToTx();
    return BasePhyLayer::setRadioState(rs);
}

simtime_t PhyLayerRms::getFrameDuration(int payloadLengthBits, MCS mcs) const
{
    Enter_Method_Silent();
    ASSERT(mcs != MCS::undefined);
    auto ndbps = getNDBPS(mcs);
    // calculate frame duration according to Equation (17-29) of the IEEE 802.11-2007 standard
    // below values are set for 802.11p 10 MHz. We consider 400 MHz so we divide such values by 40
    return (PHY_HDR_PREAMBLE_DURATION + PHY_HDR_PLCPSIGNAL_DURATION)/40 + T_SYM_80211P / 40 * ceil(static_cast<double>(16 + payloadLengthBits + 6) / (ndbps));
}

void PhyLayerRms::filterSignal(AirFrame* frame)
{
    ASSERT(dynamic_cast<ChannelAccess* const>(frame->getArrivalModule()) == this);
    ASSERT(dynamic_cast<ChannelAccess* const>(frame->getSenderModule()));
    AirFrameRms* frameRms = dynamic_cast<AirFrameRms*>(frame);

    Signal& signal = frame->getSignal();

    // Extract position and orientation of sender and receiver (this module) first
    const AntennaPosition receiverPosition = antennaPosition;
    const Coord receiverOrientation = antennaHeading.toCoord();
    // get POA from frame with the sender's position, orientation and antenna
    POA& senderPOA = frame->getPoa();
    const AntennaPosition senderPosition = senderPOA.pos;
    const Coord senderOrientation = senderPOA.orientation;

    // add position information to signal
    signal.setSenderPoa(senderPOA);
    signal.setReceiverPoa({receiverPosition, receiverOrientation, antenna});

    // compute gains at sender and receiver antenna
    double receiverGain = antenna->getGain(receiverPosition.getPositionAt(), receiverOrientation, senderPosition.getPositionAt());
    double senderGain = senderPOA.antenna->getGain(senderPosition.getPositionAt(), senderOrientation, receiverPosition.getPositionAt());

    // add the resulting total gain to the attenuations list
    //    EV_TRACE << "Sender's antenna gain: " << senderGain << endl;
    //    EV_TRACE << "Own (receiver's) antenna gain: " << receiverGain << endl;
    signal *= receiverGain * senderGain;

    // go on with AnalogueModels
    // attach analogue models suitable for thresholding to signal (for later evaluation)
    signal.setAnalogueModelList(&analogueModelsThresholding);
    // TODO: shadowing properly works (visible effects) when thresholding is set to false
    // check whether effects are not visible when thresholding is enabled because the power level is already so low that the model is not even invoked

    // apply all analouge models that are *not* suitable for thresholding now
    for (auto& analogueModel : analogueModels) {
        FrameAnalogueModel* frameAnalogueModel = dynamic_cast<FrameAnalogueModel*>(analogueModel.get());
        if (frameAnalogueModel) {
            frameAnalogueModel->filterSignal(&signal, frame);
        }
        else
            analogueModel->filterSignal(&signal);
    }

    if (isRMS) {
        // if this is an RMS, we should repropagate, but we need to attach metadata first
        Angles incident = spherical_angles(getRms_v1(), getRms_v2(), getRms_vn(), receiverPosition.getPositionAt(), senderPosition.getPositionAt());
        double recvPower_dBm = 10 * log10(signal.getAtCenterFrequency());

        EV_TRACE << "RMS: AirFrame to repropagate. Power: " << recvPower_dBm << "\n";
        Coord senderCoord = senderPOA.pos.getPositionAt();
        EV_TRACE << "Frame coming from " << senderCoord.x << " " << senderCoord.y << " " << senderCoord.z << "\n";

        // TODO: set threshold parameter
        if (recvPower_dBm > -100) {

            unique_ptr<AirFrame> toRepropagate = encapsMsg(frame->getEncapsulatedPacket()->dup(), signal.getAtCenterFrequency(), incident.phi, incident.theta, true, frame->getId());

            AirFrameRms* frameToRepropagate = dynamic_cast<AirFrameRms*>(toRepropagate.get());
            frameToRepropagate->setOriginalSenderModule(frameRms->getOriginalSenderModule());
            frameToRepropagate->setRmsGain(frameRms->getRmsGain());
            frameToRepropagate->setReflected(true);
            frameToRepropagate->setTotalDistance(frameRms->getTotalDistance());

            // Prepare a POA object and attach it to the created Airframe
            AntennaPosition pos = antennaPosition;
            Coord orient = antennaHeading.toCoord();
            toRepropagate->setPoa({pos, orient, antenna});

            // TODO: might need to completely re-think the physical layer
            // TODO: RMS is passive so there is no notion of TX/RX status
            // make sure there is no self message of kind TX_OVER scheduled
            // and schedule the actual one
            ASSERT(!txOverTimer->isScheduled());
            sendSelfMessage(txOverTimer, simTime() + toRepropagate->getDuration());

            sendMessageDown(toRepropagate.release());

        }
        else {
            EV_TRACE << "Signal too low. RMS won't repropagate the signal\n";
        }
    }

}

void PhyLayerRms::updateVehiclePosition(int vehicleId, const Coord& position)
{
    vehicles[vehicleId] = position;
}


void PhyLayerRms::requestReconfiguration(int txId, int rxId)
{
    // if one of the positions is unknown, we cannot reconfigure
    if (vehicles.find(txId) == vehicles.end() || vehicles.find(rxId) == vehicles.end())
        return;

    Coord txPos = vehicles[txId];
    Coord rxPos = vehicles[rxId];

    Angles incident = spherical_angles(rms_v1, rms_v2, rms_vn, antennaPosition.getPositionAt(), txPos);
    Angles reflected = spherical_angles(rms_v1, rms_v2, rms_vn, antennaPosition.getPositionAt(), rxPos);
    configureMetaSurface(reflected.phi, reflected.theta, incident.phi, incident.theta);
}
