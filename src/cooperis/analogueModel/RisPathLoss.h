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

#include <cassert>

#include "cooperis/cooperis.h"

#include "veins/base/phyLayer/AnalogueModel.h"
#include "veins/modules/mobility/traci/TraCIMobility.h"
#include "veins/modules/world/annotations/AnnotationManager.h"
#include "cooperis/utility/Utils.h"

#include "cooperis/PhyLayerRis.h"
#include "cooperis/FrameAnalogueModel.h"
#include "cooperis/messages/AirFrameRis_m.h"


using veins::AirFrame;
using veins::AnnotationManager;
using veins::AnnotationManagerAccess;
using veins::TraCIMobility;

namespace veins {

/**
 * @brief This class is used to compute the gain or the loss caused by RIS reflection
 *
 * @ ingroup analogueModels
 *
 */

class COOPERIS_API RisPathLoss : public SimplePathlossModel, public FrameAnalogueModel {
protected:
    AnnotationManager* annotations;

    bool debug = true;
    PhyLayerRis* phyLayer = nullptr;

    double frequency;
    double n = 2;

    bool useProductOfDistances;
    double pathLossAlpha = 2;

    double wavelength;
    double wavelengthSquare;

    void initialize();
    void finalize();

    double computeRisGain(Signal* signal, AirFrameRis* frame, const Coord& senderPos, const Coord& receiverPos);
    double computePathLoss(Signal* signal, AirFrameRis* frame, const Coord& senderPos, const Coord& receiverPos);

public:
    RisPathLoss(cComponent* owner, double alpha, bool useTorus, bool useProductOfDistances, const Coord& playgroundSize, double frequency, int n)
        : SimplePathlossModel(owner, alpha, useTorus, playgroundSize)
        , frequency(frequency)
        , n(n)
        , useProductOfDistances(useProductOfDistances)
        , pathLossAlpha(alpha)
    {
        phyLayer = check_and_cast<PhyLayerRis*>(owner);
        annotations = AnnotationManagerAccess().getIfExists();
        ASSERT(annotations);
        wavelength = BaseWorldUtility::speedOfLight() / frequency;
        wavelengthSquare = wavelength * wavelength;
    };

    virtual void filterSignal(Signal* signal) override;

    virtual void filterSignal(Signal* signal, AirFrame* frame) override;

    virtual bool neverIncreasesPower() override
    {
        return false;
    }

};
} // namespace veins
