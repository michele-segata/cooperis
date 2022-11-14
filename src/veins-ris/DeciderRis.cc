//
// Copyright (C) 2016 Agon Memedi <memedi@ccs-labs.org>
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
 * Based on Decider80211p.cc from David Eckhoff
 * and modifications by Bastian Bloessl, Stefan Joerer, Michele Segata, Fabian Bronner
 */

#include "veins-ris/DeciderRis.h"
#include "veins-ris/DeciderResultRis.h"
#include "veins/modules/phy/NistErrorRate.h"
#include "veins/base/toolbox/Signal.h"
#include "veins-ris/messages/AirFrameRis_m.h"
#include "veins/modules/utility/ConstsPhy.h"

#include "veins/modules/utility/Consts80211p.h"
#include "veins-ris/utility/Utils.h"
#include "veins/base/utils/FWMath.h"

#include "veins/base/toolbox/SignalUtils.h"

using namespace veins;

simtime_t DeciderRis::processNewSignal(AirFrame* msg)
{

    AirFrameRis* frame = check_and_cast<AirFrameRis*>(msg);

    // get the receiving power of the Signal at start-time and center frequency
    Signal& signal = frame->getSignal();
    double recvPower = signal.getAtCenterFrequency();

    if (!frame->getReflected() && ignoreNonReflectedSignals) {
        EV_TRACE << "Ignoring non reflected signal\n";
        return signal.getReceptionEnd();
    }

    signalStates[frame] = EXPECT_END;

    if (signal.smallerAtCenterFrequency(minPowerLevel) && !ignoreNoiseAndInterference) {

        EV_TRACE << "AirFrame: " << frame->getId() << " with recvPower (" << recvPower << " < " << minPowerLevel << ") -> AirFrame can't be detected by the radio; discarded at its end." << std::endl;

        // annotate the frame, so that we won't try decoding it at its end
        return signal.getReceptionEnd();
    }
    else {

        // This value might be just an intermediate result (due to short circuiting)
        setChannelIdleStatus(false);

        myBusyTime += signal.getDuration().dbl();

        if (!currentSignal.first) {
            // NIC is not yet synced to any frame, so lock and try to decode this frame
            currentSignal.first = frame;
            EV_TRACE << "AirFrame: " << frame->getId() << " with (" << recvPower << " > " << minPowerLevel << ") -> Trying to receive AirFrame." << std::endl;
        }
        else {
            // NIC is currently trying to decode another frame. this frame will be simply treated as interference
            EV_TRACE << "AirFrame: " << frame->getId() << " with (" << recvPower << " > " << minPowerLevel << ") -> Already synced to another AirFrame. Treating AirFrame as interference." << std::endl;
        }
        return signal.getReceptionEnd();
    }
}

int DeciderRis::getSignalState(AirFrame* frame)
{

    if (signalStates.find(frame) == signalStates.end()) {
        return NEW;
    }
    else {
        return signalStates[frame];
    }
}

DeciderResult* DeciderRis::checkIfSignalOk(AirFrame* frame)
{
    auto frameRis = check_and_cast<AirFrameRis*>(frame);

    Signal& s = frame->getSignal();
    simtime_t start = s.getReceptionStart();
    simtime_t end = s.getReceptionEnd();

    // / 40 because of 400 mhz bandwidth instead of 10 mhz
    start = start + PHY_HDR_PREAMBLE_DURATION / 40; // its ok if something in the training phase is broken

    AirFrameVector airFrames;
    getChannelInfo(start, end, airFrames);

    double noise = phy->getNoiseFloorValue();

    // Make sure to use the adjusted starting-point (which ignores the preamble)
    double sinrMin = SignalUtils::getMinSINR(start, end, frame, airFrames, noise);
    double snrMin;
    if (collectCollisionStats) {
        // snrMin = SignalUtils::getMinSNR(start, end, frame, noise);
        snrMin = s.getDataMin() / noise;
    }
    else {
        // just set to any value. if collectCollisionStats != true
        // it will be ignored by packetOk
        snrMin = 1e200;
    }

    double payloadBitrate = getOfdmDatarate(static_cast<MCS>(frameRis->getMcs()), Bandwidth::ofdm_400_mhz);

    DeciderResultRis* result = nullptr;

    // compute receive power
    double recvPower_dBm = 10 * log10(s.getAtCenterFrequency());

    EV_TRACE << "Packet SINR: " << 10*log10(sinrMin) << "\n";

    if (ignoreNoiseAndInterference) {
        EV_TRACE << "Ignoring noise and interference. Packet is fine! We can decode it" << std::endl;
        result = new DeciderResultRis(true, payloadBitrate, sinrMin, recvPower_dBm, false, frameRis->getRisGain(), frameRis->getReflectionPhi(), frameRis->getReflectionTheta(), frameRis->getIncidencePhi(), frameRis->getIncidenceTheta(), frameRis->getReflected());
        return result;
    }

    switch (packetOk(sinrMin, snrMin, frame->getBitLength(), payloadBitrate)) {

    case DECODED:
        EV_TRACE << "Packet is fine! We can decode it" << std::endl;
        result = new DeciderResultRis(true, payloadBitrate, sinrMin, recvPower_dBm, false, frameRis->getRisGain(), frameRis->getReflectionPhi(), frameRis->getReflectionTheta(), frameRis->getIncidencePhi(), frameRis->getIncidenceTheta(), frameRis->getReflected());
        break;

    case NOT_DECODED:
        if (!collectCollisionStats) {
            EV_TRACE << "Packet has bit Errors. Lost " << std::endl;
        }
        else {
            EV_TRACE << "Packet has bit Errors due to low power. Lost " << std::endl;
        }
        result = new DeciderResultRis(false, payloadBitrate, sinrMin, recvPower_dBm, false, frameRis->getRisGain(), frameRis->getReflectionPhi(), frameRis->getReflectionTheta(), frameRis->getIncidencePhi(), frameRis->getIncidenceTheta(), frameRis->getReflected());
        break;

    case COLLISION:
        EV_TRACE << "Packet has bit Errors due to collision. Lost " << std::endl;
        collisions++;
        result = new DeciderResultRis(false, payloadBitrate, sinrMin, recvPower_dBm, true, frameRis->getRisGain(), frameRis->getReflectionPhi(), frameRis->getReflectionTheta(), frameRis->getIncidencePhi(), frameRis->getIncidenceTheta(), frameRis->getReflected());
        break;

    default:
        ASSERT2(false, "Impossible packet result returned by packetOk(). Check the code.");
        break;
    }

    return result;
}

enum DeciderRis::PACKET_OK_RESULT DeciderRis::packetOk(double sinrMin, double snrMin, int lengthMPDU, double bitrate)
{
    double packetOkSinr;
    double packetOkSnr;

    // compute success rate depending on mcs and bw
    packetOkSinr = NistErrorRate::getChunkSuccessRate(bitrate, Bandwidth::ofdm_400_mhz, sinrMin, PHY_HDR_SERVICE_LENGTH + lengthMPDU + PHY_TAIL_LENGTH);

    // check if header is broken
    double headerNoError = NistErrorRate::getChunkSuccessRate(PHY_HDR_BITRATE * 40, Bandwidth::ofdm_400_mhz, sinrMin, PHY_HDR_PLCPSIGNAL_LENGTH);

    double headerNoErrorSnr;
    // compute PER also for SNR only
    if (collectCollisionStats) {

        packetOkSnr = NistErrorRate::getChunkSuccessRate(bitrate, Bandwidth::ofdm_400_mhz, snrMin, PHY_HDR_SERVICE_LENGTH + lengthMPDU + PHY_TAIL_LENGTH);
        headerNoErrorSnr = NistErrorRate::getChunkSuccessRate(PHY_HDR_BITRATE * 40, Bandwidth::ofdm_400_mhz, snrMin, PHY_HDR_PLCPSIGNAL_LENGTH);

        // the probability of correct reception without considering the interference
        // MUST be greater or equal than when consider it
        ASSERT(packetOkSnr >= packetOkSinr);
        ASSERT(headerNoErrorSnr >= headerNoError);
    }

    // probability of no bit error in the PLCP header

    double rand = RNGCONTEXT dblrand();

    if (!collectCollisionStats) {
        if (rand > headerNoError) return NOT_DECODED;
    }
    else {

        if (rand > headerNoError) {
            // ups, we have a header error. is that due to interference?
            if (rand > headerNoErrorSnr) {
                // no. we would have not been able to receive that even
                // without interference
                return NOT_DECODED;
            }
            else {
                // yes. we would have decoded that without interference
                return COLLISION;
            }
        }
    }

    // probability of no bit error in the rest of the packet

    rand = RNGCONTEXT dblrand();

    if (!collectCollisionStats) {
        if (rand > packetOkSinr) {
            return NOT_DECODED;
        }
        else {
            return DECODED;
        }
    }
    else {

        if (rand > packetOkSinr) {
            // ups, we have an error in the payload. is that due to interference?
            if (rand > packetOkSnr) {
                // no. we would have not been able to receive that even
                // without interference
                return NOT_DECODED;
            }
            else {
                // yes. we would have decoded that without interference
                return COLLISION;
            }
        }
        else {
            return DECODED;
        }
    }
}

simtime_t DeciderRis::processSignalEnd(AirFrame* msg)
{

    AirFrameRis* frame = check_and_cast<AirFrameRis*>(msg);

    // remove this frame from our current signals
    signalStates.erase(frame);

    DeciderResult* result;

    if (false) {
        // this frame was not even detected by the radio card
        result = new DeciderResult(false);
    }
    else {

        // first check whether this is the frame NIC is currently synced on
        if (frame == currentSignal.first) {
            // check if the snr is above the Decider's specific threshold,
            // i.e. the Decider has received it correctly
            result = checkIfSignalOk(frame);

            // after having tried to decode the frame, the NIC is no more synced to the frame
            // and it is ready for syncing on a new one
            currentSignal.first = 0;
        }
        else {
            // if this is not the frame we are synced on, we cannot receive it
            result = new DeciderResult(false);
        }
    }

    if (result->isSignalCorrect()) {
        EV_TRACE << "packet was received correctly, it is now handed to upper layer...\n";
        // go on with processing this AirFrame, send it to the Mac-Layer
        phy->sendUp(frame, result);
    }
    else {
        // To investigate why the packet was not received correctly
        delete result;
    }

    return notAgain;
}

void DeciderRis::finish()
{
}

DeciderRis::~DeciderRis(){};
