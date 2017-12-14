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

#include "veins-rms/SUMOBaseMobility.h"
#include "veins/base/utils/FindModule.h"

using namespace veins;

Define_Module(veins::SUMOBaseMobility);

SUMOBaseMobility::SUMOBaseMobility()
    : BaseMobility()
{
}

SUMOBaseMobility::SUMOBaseMobility(unsigned stacksize)
    : BaseMobility(stacksize)
{
}

void SUMOBaseMobility::initializePosition()
{
    EV_TRACE << "initializing SUMOBaseMobility\n";

    // get utility pointers (world and host)
    world = FindModule<BaseWorldUtility*>::findGlobalModule();

    // initalize position with random values
    Coord pos = world->getRandomPosition();

    if (hasStartPosition) {
        TraCICoord sumoCoordinates;
        // override start position after converting from SUMO to omnet
        sumoCoordinates.x = startPosition.x;
        sumoCoordinates.y = startPosition.y;
        Coord omnetCoordinates = manager->getConnection()->traci2omnet(sumoCoordinates);
        startPosition.x = omnetCoordinates.x;
        startPosition.y = omnetCoordinates.y;
        // notice: z will remain as original, no translation
    }
    else {
        // read coordinates from parameters if available
        double x = hasPar("x") ? par("x").doubleValue() : pos.x;
        double y = hasPar("y") ? par("y").doubleValue() : pos.y;
        double z = hasPar("z") ? par("z").doubleValue() : pos.z;

        pos.x = x;
        pos.y = y;
        pos.z = z;

        TraCICoord sumoCoordinates;
        sumoCoordinates.x = pos.x;
        sumoCoordinates.y = pos.y;
        Coord omnetCoordinates = manager->getConnection()->traci2omnet(sumoCoordinates);
        startPosition.x = omnetCoordinates.x;
        startPosition.y = omnetCoordinates.y;
        startPosition.z = pos.z;
        // trick to bypass reading the parameters and tell BaseMobility we already have a position
        hasStartPosition = true;
    }
    std::list<veins::Coord> coords;
    coords.push_back(veins::Coord(startPosition.x-5, startPosition.y-5));
    coords.push_back(veins::Coord(startPosition.x-5, startPosition.y+5));
    coords.push_back(veins::Coord(startPosition.x+5, startPosition.y+5));
    coords.push_back(veins::Coord(startPosition.x+5, startPosition.y-5));
    manager->getCommandInterface()->addPolygon("rms", "rmstype", veins::TraCIColor(0, 0, 255, 255), true, 0, coords);

    BaseMobility::initialize(0);
    BaseMobility::initialize(1);

}

void SUMOBaseMobility::initialize(int stage)
{
    if (stage == 0) {
        manager = veins::TraCIScenarioManagerAccess().get();
        ASSERT(manager);
        auto init = [this](veins::SignalPayload<bool>) {
            initializePosition();
        };
        signalManager.subscribeCallback(manager, veins::TraCIScenarioManager::traciInitializedSignal, init);
    }
}
