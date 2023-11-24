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

#include "veins/veins.h"
#include "cooperis/PhyLayerRis.h"

#include "cooperis/analogueModel/VehicleObstacleShadowingForVlc.h"
#include "veins/base/connectionManager/BaseConnectionManager.h"
#include "veins/base/toolbox/Spectrum.h"
#include "cooperis/analogueModel/RisPathLoss.h"
#include "cooperis/analogueModel/SimplePathlossModelRis.h"
#include "cooperis/analogueModel/SimpleObstacleShadowingRis.h"

#include <memory>

using namespace veins;

using std::unique_ptr;

Define_Module(veins::PhyLayerRis);

void PhyLayerRis::initialize(int stage)
{
    if (stage == 0) {
        bitrate = par("bitrate");
        centerFrequency = par("centerFrequency");
        txPower = par("txPower");
        isRIS = par("isRIS");

        ris_v1.x = par("ris_v1_x");
        ris_v1.y = par("ris_v1_y");
        ris_v1.z = par("ris_v1_z");

        ris_v2.x = par("ris_v2_x");
        ris_v2.y = par("ris_v2_y");
        ris_v2.z = par("ris_v2_z");

        if (isRIS && ((ris_v2.squareLength() < 1e-6) || (ris_v1.squareLength() < 1e-6)))
            throw cRuntimeError("Cannot instantiate an RIS without setting the orientation vectors v2 and v1");

        if (isRIS) {
            // normalize vectors
            double v1Len = ris_v1.length();
            ris_v1.x /= v1Len;
            ris_v1.y /= v1Len;
            ris_v1.z /= v1Len;
            double v2Len = ris_v2.length();
            ris_v2.x /= v2Len;
            ris_v2.y /= v2Len;
            ris_v2.z /= v2Len;
            // compute facing direction
            ris_vn = cross(ris_v1, ris_v2);

            initialConfigurationTime = par("initialConfigurationTime");
            focusBeamFrom = par("focusBeamFrom").stdstringValue();
            pointBeamTo = par("pointBeamTo").stdstringValue();
            initialIncidencePhi = par("initialIncidencePhi");
            initialIncidenceTheta = par("initialIncidenceTheta");
            initialReflectionPhi = par("initialReflectionPhi");
            initialReflectionTheta = par("initialReflectionTheta");
            destinationNodeToTrack = par("destinationNodeToTrack").stdstringValue();

            initMetasurface = new cMessage("initMetasurface");
            scheduleAt(initialConfigurationTime, initMetasurface);

            codingStates = par("codingStates");
            cellsPerLambda = par("cellsPerLambda");
            lambdaSize = par("lambdaSize");

            ris = new ReconfigurableIntelligentSurface(getSimulation()->getActiveEnvir()->getConfigEx()->getActiveRunNumber(), centerFrequency, codingStates, cellsPerLambda, lambdaSize);
            annotations = AnnotationManagerAccess().getIfExists();

            nodesAntennaHeight = par("nodesAntennaHeight");
        }

        ignoreNonReflectedSignals = par("ignoreNonReflectedSignals");
        ignoreShadowedSignals = par("ignoreShadowedSignals");
        ignoreNoiseAndInterference = par("ignoreNoiseAndInterference");

        useProductOfDistances = par("useProductOfDistances");

        repropagationThreshold_mW = par("repropagationThreshold");

        // Create frequency mappings and initialize spectrum for signal representation
        overallSpectrum = Spectrum({centerFrequency});
    }
    BasePhyLayer::initialize(stage);
}

const Coord& PhyLayerRis::getRis_v2() const
{
    return ris_v2;
}
const Coord& PhyLayerRis::getRis_v1() const
{
    return ris_v1;
}
const Coord& PhyLayerRis::getRis_vn() const
{
    return ris_vn;
}

bool PhyLayerRis::isReflectiveMetaSurface() const
{
    return isRIS;
}

double PhyLayerRis::getMetasurfaceGain(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad) const
{
    if (!isRIS) throw cRuntimeError("Requesting metasurface gain on a PHY that is not a metasurface but a standard node");
    double configPhiR_rad, configPhiI_rad, configThetaR_rad, configThetaI_rad;
    ris->getConfiguration(configPhiR_rad, configThetaR_rad, configPhiI_rad, configThetaI_rad);
    double gain = ris->gain(phiR_rad, thetaR_rad, phiI_rad, thetaI_rad);
    return gain;
}

void PhyLayerRis::configureMetaSurface(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad)
{
    if (!isRIS) throw cRuntimeError("Requesting to configure a metasurface on a PHY that is not a metasurface but a standard node");
    ris->configureMetaSurface(phiR_rad, thetaR_rad, phiI_rad, thetaI_rad);
    ris->getConfiguration(phiR_rad, thetaR_rad, phiI_rad, thetaI_rad);
    // visualize using AnnotationManager
    if (configVisualizationR)
        annotations->erase(configVisualizationR);
    if (configVisualizationI)
        annotations->erase(configVisualizationI);

    if (annotations) {
        std::list<Coord> projectedBeam;
        projectedBeam = project_beam(ris_v1, ris_v2, ris_vn, antennaPosition.getPositionAt(), phiR_rad, thetaR_rad, nodesAntennaHeight);
        configVisualizationR = annotations->drawPolygon(projectedBeam, "green", nullptr);
        projectedBeam = project_beam(ris_v1, ris_v2, ris_vn, antennaPosition.getPositionAt(), phiI_rad, thetaI_rad, nodesAntennaHeight);
        configVisualizationI = annotations->drawPolygon(projectedBeam, "red", nullptr);
    }

}

void PhyLayerRis::configureMetaSurfaceReflection(double phiR_rad, double thetaR_rad)
{
    configureMetaSurface(phiR_rad, thetaR_rad, KEEP_SAME_ANGLE, KEEP_SAME_ANGLE);
}
void PhyLayerRis::configureMetaSurfaceIncidence(double phiI_rad, double thetaI_rad)
{
    configureMetaSurface(KEEP_SAME_ANGLE, KEEP_SAME_ANGLE, phiI_rad, thetaI_rad);
}

unique_ptr<AnalogueModel> PhyLayerRis::getAnalogueModelFromName(std::string name, ParameterMap& params)
{

    if (name == "SimpleObstacleShadowingRis") {
        return initializeSimpleObstacleShadowingRis(params);
    }
    if (name == "RisPathLoss") {
        return initializeRisPathLoss(params);
    }
    return BasePhyLayer::getAnalogueModelFromName(name, params);
}

unique_ptr<AnalogueModel> PhyLayerRis::initializeSimpleObstacleShadowingRis(ParameterMap& params)
{

    // init with default value
    bool useTorus = world->useTorus();
    const Coord& playgroundSize = *(world->getPgs());

    ParameterMap::iterator it;

    ObstacleControl* obstacleControlP = dynamic_cast<ObstacleControl*>(veins::findModuleByPath("obstaclesRis"));
    if (!obstacleControlP) throw cRuntimeError("initializeSimpleObstacleShadowingRis(): cannot find ObstacleControl module");
    return make_unique<SimpleObstacleShadowingRis>(this, *obstacleControlP, useTorus, playgroundSize);
}

unique_ptr<AnalogueModel> PhyLayerRis::initializeRisPathLoss(ParameterMap& params)
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
        EV_TRACE << "initializeRisPathLoss(): alpha set from config.xml to " << alpha << endl;

        // check whether alpha is not smaller than specified in ConnectionManager
        if (cc->hasPar("alpha") && alpha < cc->par("alpha").doubleValue()) {
            // throw error
            throw cRuntimeError("initializeRisPathLoss(): alpha can't be smaller than specified in \
                   ConnectionManager. Please adjust your config.xml file accordingly");
        }
    }
    else // alpha has not been specified in config.xml
    {
        if (cc->hasPar("alpha")) { // parameter alpha has been specified in ConnectionManager
            // set alpha according to ConnectionManager
            alpha = cc->par("alpha").doubleValue();
            EV_TRACE << "initializeRisPathLoss(): alpha set from ConnectionManager to " << alpha << endl;
        }
        else // alpha has not been specified in ConnectionManager
        {
            // keep alpha at default value
            EV_TRACE << "initializeRisPathLoss(): alpha set from default value to " << alpha << endl;
        }
    }

    double centerFrequency = params["carrierFrequency"].doubleValue();
    int n = params["n"].longValue();
    return make_unique<RisPathLoss>(this, alpha, useTorus, useProductOfDistances, playgroundSize, centerFrequency, n);
}

unique_ptr<Decider> PhyLayerRis::getDeciderFromName(std::string name, ParameterMap& params)
{
    if (name == "DeciderRis") {
        protocolId = RIS;
        return initializeDeciderRis(params);
    }
    return BasePhyLayer::getDeciderFromName(name, params);
}

unique_ptr<AnalogueModel> PhyLayerRis::initializeObstacleShadowing(ParameterMap& params)
{
    // TODO: finish
    // init with default value
    bool useTorus = world->useTorus();
    const Coord& playgroundSize = *(world->getPgs());

    VehicleObstacleControl* vehicleObstacleControlP = VehicleObstacleControlAccess().getIfExists();
    if (!vehicleObstacleControlP) throw cRuntimeError("initializeVehicleObstacleShadowingForVlc(): cannot find VehicleObstacleControl module");
    return make_unique<VehicleObstacleShadowingForVlc>(this, *vehicleObstacleControlP, useTorus, playgroundSize);
}

unique_ptr<Decider> PhyLayerRis::initializeDeciderRis(ParameterMap& params)
{
    DeciderRis* dec = new DeciderRis(this, this, minPowerLevel, bitrate, ignoreNonReflectedSignals, ignoreShadowedSignals, ignoreNoiseAndInterference, -1, false);
    return std::unique_ptr<DeciderRis>(std::move(dec));
}

void PhyLayerRis::handleSelfMessage(cMessage* msg)
{

    switch (msg->getKind()) {
    // transmission over
    case TX_OVER:
        ASSERT(msg == txOverTimer);
        if (!isRIS)
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
        Coord myPos = antennaPosition.getPositionAt();
        if (focusBeamFrom.compare("") != 0) {
            cModule* module = findModuleByPath(focusBeamFrom.c_str());
            if (module) {
                BaseMobility* mobility = FindModule<BaseMobility*>::findSubModule(module);
                Coord otherPos = mobility->getPositionAt(simTime());
                Angles incident = spherical_angles(ris_v1, ris_v2, ris_vn, myPos, otherPos);
                initialIncidencePhi = incident.phi;
                initialIncidenceTheta = incident.theta;
            }
            else {
                std::cout << getFullPath() << ": Ignoring focusBeamFrom parameter as module " << focusBeamFrom << " cannot be found\n";
            }
        }
        if (pointBeamTo.compare("") != 0) {
            Angles reflected;
            if (pointBeamTowards(pointBeamTo, reflected)) {
                initialReflectionPhi = reflected.phi;
                initialReflectionTheta = reflected.theta;
            }
            else {
                std::cout << getFullPath() << ": Ignoring pointBeamTo parameter as module " << pointBeamTo << " cannot be found\n";
            }
        }
        configureMetaSurface(initialReflectionPhi, initialReflectionTheta, initialIncidencePhi, initialIncidenceTheta);
        cancelAndDelete(initMetasurface);
        initMetasurface = nullptr;
    }

}

void PhyLayerRis::handleMessage(cMessage* msg)
{
    // self messages
    if (msg->isSelfMessage()) {
        handleSelfMessage(msg);

        // MacPkts <- MacToPhyControlInfo
    }
    else if (msg->getArrivalGateId() == upperLayerIn) {
        if (getRadioState() == veins::Radio::TX) throw cRuntimeError("PhyLayerRIS: received a frame to transmit from MAC layer while already transmitting\n");
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
        AirFrameRis* risMsg = check_and_cast<AirFrameRis*>(msg);

        if (risMsg->getOriginalSenderModule() == this) {
            delete msg;
            return;
        }
        else {
            BasePhyLayer::handleAirFrame(static_cast<AirFrame*>(msg));
        }
    }
    else {
        delete msg;
    }
}

void PhyLayerRis::sendUp(AirFrame* frame, DeciderResult* result)
{
    if (!isRIS)
        BasePhyLayer::sendUp(frame, result);
}

void PhyLayerRis::handleAirFrameReceiving(AirFrame* msg)
{
    // if this is an RIS, do not attempt decoding
    if (!isRIS)
        BasePhyLayer::handleAirFrameReceiving(msg);
}

unique_ptr<AirFrame> PhyLayerRis::createAirFrame(cPacket* macPkt)
{
    return make_unique<AirFrameRis>(macPkt->getName(), AIR_FRAME);
}

unique_ptr<AirFrame> PhyLayerRis::encapsMsg(cPacket* macPkt, double txPower, bool reflected, int originalId)
{
    auto airFrame = createAirFrame(macPkt);
    AirFrameRis* frame = dynamic_cast<AirFrameRis*>(airFrame.get());

    // set the members
    // set priority of AirFrames above the normal priority to ensure
    // channel consistency (before any thing else happens at a time
    // point t make sure that the channel has removed every AirFrame
    // ended at t and added every AirFrame started at t)
    frame->setSchedulingPriority(airFramePriority());
    frame->setProtocolId(myProtocolId());
    frame->setId(world->getUniqueAirFrameId());
    frame->setChannel(radio->getCurrentChannel());

    frame->setReflected(reflected);
    if (originalId >= 0)
        frame->setOriginalId(originalId);
    else
        // set the original id of a new frame to its id, otherwise RIS will not repropagate it (ID 0 already seen!)
        frame->setOriginalId(frame->getId());
    frame->setOriginalSenderModule(this);

    if (isRIS) {
        frame->setRis_v1(ris_v1);
        frame->setRis_v2(ris_v2);
        frame->setRis_vn(ris_vn);
    }
    else {
        frame->setTotalDistance(0);
        frame->setPathsArraySize(0);
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

    return airFrame;
}

unique_ptr<AirFrame> PhyLayerRis::encapsMsg(cPacket* macPkt)
{
    return encapsMsg(macPkt, txPower);
}

//simtime_t PhyLayerRis::setRadioState(int rs)
//{
//    if (rs == Radio::TX) decider->switchToTx();
//    else if (rs == Radio::RX) decider->switchToRx();
//    return BasePhyLayer::setRadioState(rs);
//}

simtime_t PhyLayerRis::getFrameDuration(int payloadLengthBits, MCS mcs) const
{
    Enter_Method_Silent();
    ASSERT(mcs != MCS::undefined);
    auto ndbps = getNDBPS(mcs);
    // calculate frame duration according to Equation (17-29) of the IEEE 802.11-2007 standard
    // below values are set for 802.11p 10 MHz. We consider 400 MHz so we divide such values by 40
    return (PHY_HDR_PREAMBLE_DURATION + PHY_HDR_PLCPSIGNAL_DURATION)/40 + T_SYM_80211P / 40 * ceil(static_cast<double>(16 + payloadLengthBits + 6) / (ndbps));
}

bool PhyLayerRis::alreadyReflected(const AirFrameRis* frame)
{
    for (int i = 0; i < frame->getReflectingRISArraySize(); i++)
        if (frame->getReflectingRIS(i) == getId())
            return true;
    return false;
}

void PhyLayerRis::filterSignal(AirFrame* frame)
{
    ASSERT(dynamic_cast<ChannelAccess* const>(frame->getArrivalModule()) == this);
    ASSERT(dynamic_cast<ChannelAccess* const>(frame->getSenderModule()));
    AirFrameRis* frameRis = dynamic_cast<AirFrameRis*>(frame);

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

    if (isRIS) {

        EV_TRACE << "RIS: checking whether AirFrame " << frameRis->getOriginalId() << " from module " << frameRis->getSenderModuleId() << " should be reflected\n";

        if (frameRis->getShadowed() && ignoreShadowedSignals) {
            EV_TRACE << "Frame was crossing a building and won't be reflected\n";
        }
        else if (alreadyReflected(frameRis)) {
            EV_TRACE << "RIS: AirFrame " << frameRis->getOriginalId() << " already reflected\n";
        }
        else {
            EV_TRACE << "RIS: AirFrame " << frameRis->getOriginalId() << " should be reflected\n";

            // before reflecting, reconfigure RIS to point towards a certain destination
            if (destinationNodeToTrack.compare("") != 0) {
                Angles reflected;
                if (pointBeamTowards(destinationNodeToTrack, reflected)) {
                    configureMetaSurfaceReflection(reflected.phi, reflected.theta);
                }
                else {
                    std::cout << getFullPath() << ": Not reconfiguring RIS as module " << destinationNodeToTrack << " cannot be found\n";
                }
            }

            // if this is an RIS, we should reflect, but we need to attach metadata first
            Angles incident = spherical_angles(getRis_v1(), getRis_v2(), getRis_vn(), receiverPosition.getPositionAt(), senderPosition.getPositionAt());
            double recvPower_mW = signal.getAtCenterFrequency();

            EV_TRACE << "RIS: AirFrame to reflect. Power: " << 10*log10(recvPower_mW) << "\n";
            Coord senderCoord = senderPOA.pos.getPositionAt();
            EV_TRACE << "Frame coming from " << senderCoord.x << " " << senderCoord.y << " " << senderCoord.z << "\n";


            EV_TRACE << "recv mw = " << recvPower_mW << "\n";
            EV_TRACE << "threshold mw = " << repropagationThreshold_mW << "\n";
            // TODO: set threshold parameter
            if (recvPower_mW > repropagationThreshold_mW) {

                unique_ptr<AirFrame> toReflect = encapsMsg(frame->getEncapsulatedPacket()->dup(), signal.getAtCenterFrequency(), true, frameRis->getOriginalId());

                AirFrameRis* frametoReflect = dynamic_cast<AirFrameRis*>(toReflect.get());
                frametoReflect->setOriginalSenderModule(frameRis->getOriginalSenderModule());

                frametoReflect->setRisGain_dBArraySize(frameRis->getRisGain_dBArraySize());
                for (int i = 0; i < frameRis->getRisGain_dBArraySize(); i++)
                    frametoReflect->setRisGain_dB(i, frameRis->getRisGain_dB(i));

                frametoReflect->setPathsArraySize(frameRis->getPathsArraySize());
                for (int i = 0; i < frameRis->getPathsArraySize(); i++)
                    frametoReflect->setPaths(i, frameRis->getPaths(i));

                frametoReflect->setPathLoss_dBArraySize(frameRis->getPathLoss_dBArraySize());
                for (int i = 0; i < frameRis->getPathLoss_dBArraySize(); i++)
                    frametoReflect->setPathLoss_dB(i, frameRis->getPathLoss_dB(i));

                frametoReflect->setActualLoss_dBArraySize(frameRis->getActualLoss_dBArraySize());
                for (int i = 0; i < frameRis->getActualLoss_dBArraySize(); i++)
                    frametoReflect->setActualLoss_dB(i, frameRis->getActualLoss_dB(i));

                frametoReflect->setReflectingRISArraySize(frameRis->getReflectingRISArraySize());
                for (int i = 0; i < frameRis->getReflectingRISArraySize(); i++)
                    frametoReflect->setReflectingRIS(i, frameRis->getReflectingRIS(i));

                frametoReflect->setReflectionPhiArraySize(frameRis->getReflectionPhiArraySize());
                for (int i = 0; i < frameRis->getReflectionPhiArraySize(); i++)
                    frametoReflect->setReflectionPhi(i, frameRis->getReflectionPhi(i));

                frametoReflect->setReflectionThetaArraySize(frameRis->getReflectionThetaArraySize());
                for (int i = 0; i < frameRis->getReflectionThetaArraySize(); i++)
                    frametoReflect->setReflectionTheta(i, frameRis->getReflectionTheta(i));

                frametoReflect->setIncidencePhiArraySize(frameRis->getIncidencePhiArraySize());
                for (int i = 0; i < frameRis->getIncidencePhiArraySize(); i++)
                    frametoReflect->setIncidencePhi(i, frameRis->getIncidencePhi(i));

                frametoReflect->setIncidenceThetaArraySize(frameRis->getIncidenceThetaArraySize());
                for (int i = 0; i < frameRis->getIncidenceThetaArraySize(); i++)
                    frametoReflect->setIncidenceTheta(i, frameRis->getIncidenceTheta(i));

                frametoReflect->setTotalDistance(frameRis->getTotalDistance());

                // add this RIS to the list of RIS that have been reflecting the signal
                frametoReflect->appendReflectingRIS(getId());
                // append incidence angles
                frametoReflect->appendIncidencePhi(incident.phi);
                frametoReflect->appendIncidenceTheta(incident.theta);

                // Prepare a POA object and attach it to the created Airframe
                AntennaPosition pos = antennaPosition;
                Coord orient = antennaHeading.toCoord();
                toReflect->setPoa({pos, orient, antenna});

                sendMessageDown(toReflect.release());

            }
            else {
                EV_TRACE << "Signal too low. RIS won't reflect the signal\n";
            }

        }
    }

}

bool PhyLayerRis::pointBeamTowards(string nodeId, Angles &angles)
{
    cModule* module = findModuleByPath(nodeId.c_str());
    if (module) {
        Coord myPos = antennaPosition.getPositionAt();
        BaseMobility* mobility = FindModule<BaseMobility*>::findSubModule(module);
        Coord otherPos = mobility->getPositionAt(simTime());
        angles = spherical_angles(ris_v1, ris_v2, ris_vn, myPos, otherPos);
        return true;
    }
    else {
        return false;
    }
}

void PhyLayerRis::updateVehiclePosition(int vehicleId, const Coord& position)
{
    vehicles[vehicleId] = position;
}


void PhyLayerRis::requestReconfiguration(int txId, int rxId)
{
    // if one of the positions is unknown, we cannot reconfigure
    if (vehicles.find(txId) == vehicles.end() || vehicles.find(rxId) == vehicles.end())
        return;

    Coord txPos = vehicles[txId];
    Coord rxPos = vehicles[rxId];

    Angles incident = spherical_angles(ris_v1, ris_v2, ris_vn, antennaPosition.getPositionAt(), txPos);
    Angles reflected = spherical_angles(ris_v1, ris_v2, ris_vn, antennaPosition.getPositionAt(), rxPos);
    configureMetaSurface(reflected.phi, reflected.theta, incident.phi, incident.theta);
}

void PhyLayerRis::requestReconfiguration(int nodeId, bool incidence)
{
    if (vehicles.find(nodeId) == vehicles.end())
        return;

    Coord nodePos = vehicles[nodeId];
    Angles angles = spherical_angles(ris_v1, ris_v2, ris_vn, antennaPosition.getPositionAt(), nodePos);
    if (incidence)
        configureMetaSurfaceIncidence(angles.phi, angles.theta);
    else
        configureMetaSurfaceReflection(angles.phi, angles.theta);
}
