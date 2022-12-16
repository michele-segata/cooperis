//
// Copyright (C) 2006-2018 Christoph Sommer <sommer@ccs-labs.org>
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

#include "SimpleObstacleShadowingRis.h"
#include "veins-ris/messages/AirFrameRis_m.h"

using namespace veins;

using veins::AirFrame;

SimpleObstacleShadowingRis::SimpleObstacleShadowingRis(cComponent* owner, ObstacleControl& obstacleControl, bool useTorus, const Coord& playgroundSize)
    : SimpleObstacleShadowing(owner, obstacleControl, useTorus, playgroundSize)
{
    if (useTorus) throw cRuntimeError("SimpleObstacleShadowing does not work on torus-shaped playgrounds");
}
void SimpleObstacleShadowingRis::filterSignal(Signal* signal, AirFrame* frame)
{
    AirFrameRis* frameRis = check_and_cast<AirFrameRis*>(frame);

    auto senderPos = signal->getSenderPoa().pos.getPositionAt();
    auto receiverPos = signal->getReceiverPoa().pos.getPositionAt();

    double factor = obstacleControl.calculateAttenuation(senderPos, receiverPos);

    // tag the frame as being shadowed
    if (factor < 1)
        frameRis->setShadowed(true);

    EV_TRACE << "value is: " << factor << endl;

    *signal *= factor;
}

void SimpleObstacleShadowingRis::filterSignal(Signal* signal)
{
    throw cRuntimeError("SimpleObstacleShadowingRis::filterSignal(Signal* signal) must not be used as the model requires the AirFrame as well to perform computation");
}
