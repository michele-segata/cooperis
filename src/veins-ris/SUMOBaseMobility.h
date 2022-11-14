//
// Copyright (C) 2004 Telecommunication Networks Group (TKN) at Technische Universitaet Berlin, Germany.
// Copyright (C) 2005 Andras Varga
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

// author:      Daniel Willkomm, Andras Varga
// part of:     framework implementation developed by tkn

#pragma once

#include "veins/base/modules/BaseMobility.h"
#include "veins/modules/utility/SignalManager.h"
#include "veins/modules/mobility/traci/TraCIScenarioManager.h"
#include "veins/modules/mobility/traci/TraCICommandInterface.h"

namespace veins {

class VEINS_API SUMOBaseMobility : public BaseMobility {

public:
    SUMOBaseMobility();
    SUMOBaseMobility(unsigned stacksize);

    void initialize(int) override;
    int numInitStages() const override { return 1; }


    void initializePosition();
    //    using BaseMobility::receiveSignal;
    //    void receiveSignal(cComponent* source, simsignal_t signalID, bool v, cObject* details) override;
    //    void receiveSignal(cComponent* source, simsignal_t signalID, bool v)
    //    {
    //        receiveSignal(source, signalID, v, 0);
    //    }
protected:
    TraCIScenarioManager* manager;
    veins::SignalManager signalManager;

};

} // namespace veins
