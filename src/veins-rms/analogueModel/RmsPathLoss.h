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

#pragma once

#include <cassert>

#include "veins-rms/veins-rms.h"

#include "veins/base/phyLayer/AnalogueModel.h"
#include "veins-rms/utility/ConstsVlc.h"
#include "veins/modules/mobility/traci/TraCIMobility.h"
#include "veins/modules/world/annotations/AnnotationManager.h"
#include "veins-rms/utility/Utils.h"

#include "veins-rms/PhyLayerRms.h"
#include "veins-rms/FrameAnalogueModel.h"


using veins::AirFrame;
using veins::AnnotationManager;
using veins::AnnotationManagerAccess;
using veins::TraCIMobility;

namespace veins {

/**
 * @brief This class is used for two purposes:
 * 1) in the path from UE to RMS, it annotates the signal with incidence angles and computes the attenuation
 * 2) in the path from RMS to UE (or to RMS) (reflected) it computes the path loss considering both the incidence and the reflection angles
 *
 * @ ingroup analogueModels
 *
 */

class VEINS_RMS_API RmsPathLoss : public AnalogueModel, public FrameAnalogueModel {
protected:
    AnnotationManager* annotations;

    bool debug = true;
    PhyLayerRms* phyLayer;

    double frequency;
    double n = 2;

    void initialize();
    void finalize();

public:
    RmsPathLoss(cComponent* owner, double frequency, int n)
        : AnalogueModel(owner)
        , frequency(frequency)
        , n(n)
    {
        phyLayer = check_and_cast<PhyLayerRms*>(owner);
        annotations = AnnotationManagerAccess().getIfExists();
        ASSERT(annotations);
    };

    virtual void filterSignal(Signal* signal) override;

    virtual void filterSignal(Signal* signal, AirFrame* frame) override;

    virtual bool neverIncreasesPower() override
    {
        return false;
    }

};
} // namespace veins
