//
// Copyright (C) 2022-2024 Michele Segata <segata@ccs-labs.org>
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

#pragma once

#include "veins/base/phyLayer/BasePhyLayer.h"
#include "veins/base/toolbox/Spectrum.h"
#include "veins/modules/mac/ieee80211p/Mac80211pToPhy11pInterface.h"
#include "cooperis/DeciderRis.h"
#include "cooperis/utility/ReconfigurableIntelligentSurface.h"
#include "veins/modules/analogueModel/SimplePathlossModel.h"
#include "veins/base/connectionManager/BaseConnectionManager.h"
#include "veins/modules/phy/Decider80211pToPhy80211pInterface.h"
#include "veins/base/utils/Move.h"
#include "veins/modules/world/annotations/AnnotationManager.h"
#include "cooperis/messages/AirFrameRis_m.h"
#include "cooperis/utility/Utils.h"

namespace veins {

/**
 * @brief
 * Adaptation of the PhyLayer class for 802.11p.
 *
 * @ingroup phyLayer
 *
 * @see DemoBaseApplLayer
 * @see Mac1609_4
 * @see PhyLayer80211p
 * @see Decider80211p
 */

class PhyLayerRis : public BasePhyLayer {
public:
    void initialize(int stage) override;

    bool isReflectiveMetaSurface() const;

    const Coord& getRis_v2() const;
    const Coord& getRis_v1() const;
    const Coord& getRis_vn() const;

    double getMetasurfaceGain(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad) const;

    void configureMetaSurface(double phiR_rad, double thetaR_rad, double phiI_rad, double thetaI_rad);
    void configureMetaSurfaceReflection(double phiR_rad, double thetaR_rad);
    void configureMetaSurfaceIncidence(double phiI_rad, double thetaI_rad);

    virtual int numInitStages() const override { return 2; };

    /**
     * Updates the position of a vehicle
     */
    void updateVehiclePosition(int vehicleId, const Coord& position);

    /**
     * Request the reconfiguration of the RIS for a tx/rx pair of vehicles
     * If the vehicle positions are unknown reconfiguration is not performed
     */
    void requestReconfiguration(int txId, int rxId);
    void requestReconfiguration(int nodeId, bool incidence);

    bool pointBeamTowards(string nodeId, Angles& angles);

protected:

    enum ProtocolIds {
        RIS = 12125
    };

    int bitrate;
    // center freq in Hz
    double centerFrequency;
    // tx power in mW
    double txPower;

    // determine whether this is an RIS or an UE
    bool isRIS;
    // vector indicating the orientation of the RIS
    // see the .ned file for a detailed description
    Coord ris_v2;
    Coord ris_v1;
    // vector indicating the direction the RIS is facing
    // compute as the cross product between v1 and v2
    // cross product is between v1 and v2 as a convention, so that if x axis is v1 and negative z axis is v2,
    // then vn = v1 x v2 = y axis
    Coord ris_vn;

    bool randomCombinationStrategy = true;

    double initialConfigurationTime = 1e-6;
    string focusBeamFrom = "";
    string pointBeamTo = "";
    double initialIncidencePhi = 0;
    double initialIncidenceTheta = 0;
    double initialReflectionPhi = 0;
    double initialReflectionTheta = 0;
    vector<double> initialReflectionsPhi_rad = {};
    vector<double> initialReflectionsTheta_rad = {};
    vector<string> destinationNodesToTrack = {};

    int codingStates = 4;
    int cellsPerLambda = 3;
    int lambdaSize = 5;
    double d_phi = M_PI / 180;
    double d_theta = M_PI / 180;
    double txGain = 0;
    double rxGain = 0;

    ReconfigurableIntelligentSurface* ris = nullptr;

    AnnotationManager* annotations = nullptr;
    AnnotationManager::Annotation* configVisualizationR = nullptr;
    AnnotationManager::Annotation* configVisualizationI = nullptr;
    // just for graphically drawing beams
    double nodesAntennaHeight = 0;

    cMessage* initMetasurface = nullptr;

    bool ignoreNonReflectedSignals = false;
    bool ignoreShadowedSignals = false;
    bool ignoreNoiseAndInterference = false;

    bool useProductOfDistances = false;

    double repropagationThreshold_mW = 10e-10;

    // map from vehicle id to their coordinates, to enable reconfiguration of the RIS
    std::map<int, Coord> vehicles;

    bool alreadyReflected(const AirFrameRis* frame);

    /**
     * @brief Creates and returns an instance of the AnalogueModel with the
     * specified name.
     *
     * Is able to initialize the following AnalogueModels:
     */
    virtual std::unique_ptr<AnalogueModel> getAnalogueModelFromName(std::string name, ParameterMap& params) override;

    /**
     * @brief Creates and initializes ObstacleShadowing with the
     * passed parameter values.
     */
    std::unique_ptr<AnalogueModel> initializeObstacleShadowing(ParameterMap& params);

    unique_ptr<AnalogueModel> initializeSimpleObstacleShadowingRis(ParameterMap& params);

    /**
     * @brief Creates and initializes a RisPathLoss with the
     * passed parameter values.
     */
    std::unique_ptr<AnalogueModel> initializeRisPathLoss(ParameterMap& params);

    std::unique_ptr<AnalogueModel> initializeETSI_TR_138_901(ParameterMap& params);

    /**
     * @brief Creates and returns an instance of the Decider with the specified
     * name.
     *
     * Is able to initialize the following Deciders:
     *
     * - DeciderRis
     */
    virtual std::unique_ptr<Decider> getDeciderFromName(std::string name, ParameterMap& params) override;

    /**
     * @brief Initializes a new DeciderRis from the passed parameter map.
     */
    virtual std::unique_ptr<Decider> initializeDeciderRis(ParameterMap& params);

    /**
     * Create a protocol-specific AirFrame
     * Overloaded to create a specialize AirFrameRis.
     */
    std::unique_ptr<AirFrame> createAirFrame(cPacket* macPkt) override;

    /**
     * @brief This function encapsulates messages from the upper layer into an
     * AirFrame and sets all necessary attributes.
     */
    virtual std::unique_ptr<AirFrame> encapsMsg(cPacket* msg) override;
    virtual std::unique_ptr<AirFrame> encapsMsg(cPacket* macPkt, double txPower, bool reflected= false, int originalId= -1);

    virtual void handleAirFrameReceiving(AirFrame* msg) override;

    simtime_t getFrameDuration(int payloadLengthBits, MCS mcs) const;

    virtual void handleMessage(cMessage* msg) override;
    virtual void handleSelfMessage(cMessage* msg) override;
    virtual void sendUp(AirFrame* frame, DeciderResult* result) override;
    //    simtime_t setRadioState(int rs) override;

    /**
     * Override filter to be able to add metadata to the AirFrame and be able to repropagate it immediately (reflect)
     */
    virtual void filterSignal(AirFrame* frame) override;

};

} // namespace veins
