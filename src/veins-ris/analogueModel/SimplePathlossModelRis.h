//
// Copyright (C) 2022 Michele Segata <segata@ccs-labs.org>
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

#include "veins/veins.h"
#include "veins/modules/analogueModel/SimplePathlossModel.h"
#include "veins-ris/FrameAnalogueModel.h"
#include "veins-ris/PhyLayerRis.h"

namespace veins {

using veins::AirFrame;

class SimplePathlossModelRis;

/**
 * @brief Basic implementation of a SimplePathlossModelRis
 * Derived from SimplePathlossModel, it is only applied in case the
 * source node of a signal is a UE and not an RIS
 *
 * An example config.xml for this AnalogueModel can be the following:
 * @verbatim
    <AnalogueModel type="SimplePathlossModelRis">
        <!-- Environment parameter of the pathloss formula
             If ommited default value is 3.5-->
        <parameter name="alpha" type="double" value="3.5"/>
    </AnalogueModel>
   @endverbatim
 *
 * @ingroup analogueModels
 */
class VEINS_API SimplePathlossModelRis : public SimplePathlossModel, public FrameAnalogueModel {

protected:

    PhyLayerRis* receiver = nullptr;

public:

    SimplePathlossModelRis(cComponent* owner, double alpha, bool useTorus, const Coord& playgroundSize)
        : SimplePathlossModel(owner, alpha, useTorus, playgroundSize)
    {
        receiver = check_and_cast<PhyLayerRis*>(owner);
    }

    /**
     * @brief Filters a specified AirFrame's Signal by adding an attenuation
     * over time to the Signal.
     */
    void filterSignal(Signal*) override;
    virtual void filterSignal(Signal*, AirFrame*) override;
};

} // namespace veins
