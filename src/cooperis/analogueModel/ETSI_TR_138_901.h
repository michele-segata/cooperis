//
// Copyright (C) 2024 Michele Segata <segata@ccs-labs.org>
//
// Documentation for these modules is at http://veins.car2x.org/
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

#include <cstdlib>

#include "veins/base/phyLayer/AnalogueModel.h"
#include "veins/base/modules/BaseWorldUtility.h"
#include "veins/modules/obstacle/ObstacleControl.h"
#include "veins/base/utils/Move.h"
#include "veins/base/messages/AirFrame_m.h"
#include "veins/modules/analogueModel/SimpleObstacleShadowing.h"
#include "cooperis/FrameAnalogueModel.h"

using veins::AirFrame;
using veins::ObstacleControl;

namespace veins {

class Signal;

/**
 * Implementation of path loss model in 3GPP TR 38.901 (ETSI TR 138.901).
 * The model currently only implements the Urban Micro model, which models the pathloss between
 * a base station and a UE, when the base station is located below the rooftop.
 * This is not the proper model for V2V modeling but currently this is mentioned as future
 * work in the report (ETSI TR 138 901 V16.1.0 (2020-11), page 14).
 */
class VEINS_API ETSI_TR_138_901 : public SimpleObstacleShadowing, public FrameAnalogueModel {
public:
    /**
     * @brief Initializes the analogue model. myMove and playgroundSize
     * need to be valid as long as this instance exists.
     *
     * The constructor needs some specific knowledge in order to create
     * its mapping properly:
     *
     * @param owner pointer to the cComponent that owns this AnalogueModel
     * @param obstacleControl the parent module
     * @param useTorus information about the playground the host is moving in
     * @param playgroundSize information about the playground the host is moving in
     */
    ETSI_TR_138_901(cComponent* owner, ObstacleControl& obstacleControl, bool useTorus, const Coord& playgroundSize);

    /**
     * @brief Filters a specified Signal by adding an attenuation
     * over time to the Signal.
     */
    void filterSignal(Signal*) override;
    virtual void filterSignal(Signal*, AirFrame*) override;

    bool neverIncreasesPower() override
    {
        return true;
    }

protected:
    double lineOfSightPathLoss(double f, const Coord &sender, const Coord &receiver);
    double nonLineOfSightPathLoss(double f, const Coord &sender, const Coord &receiver);
};

} // namespace veins
